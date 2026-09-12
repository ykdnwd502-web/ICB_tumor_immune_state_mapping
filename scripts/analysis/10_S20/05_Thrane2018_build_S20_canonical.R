# S20 Thrane 2018 legacy-ST external spatial recurrence audit
# CANONICAL v1.0.2 PARSER HOTFIX — analytical contract v1.0; source-locked to historical 18J v4 FINAL
#
# Purpose:
#   External melanoma spatial recurrence audit using Thrane et al. 2018
#   legacy spatial transcriptomics data.
#
# Key fixes:
#   1) Does NOT attach AnnotationDbi or org.Hs.eg.db with library().
#      Uses AnnotationDbi::mapIds() and org.Hs.eg.db::org.Hs.eg.db directly.
#      Therefore AnnotationDbi::select cannot mask dplyr::select.
#   2) Avoids naked select(); uses dplyr::select() or base R.
#   3) Supports Ensembl / Entrez / HGNC / mixed gene IDs.
#   4) Reads Thrane count matrices as genes x spots or spots x genes.
#   5) Selects orientation by overlap with locked scoring genes.
#   6) Outputs gene coverage, key readout, and inclusion decision.
#   7) v1.0.2 parser-only fix: Thrane first-column IDs are compound
#      identifiers of the form "SYMBOL ENSG...". Parse the leading HGNC
#      symbol explicitly before generic ID mapping. No analytical/statistical
#      definitions are changed.
#
# Boundary:
#   External spatial recurrence audit only.
#   Not ICB-response validation.
#   Not biomarker validation.
#   Not functional niche validation.
# ============================================================

rm(list = ls())
gc()

# -----------------------------
# 0. Configuration
# -----------------------------

PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)
PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)
DATASET_ID <- "Thrane2018_legacyST"

BASE_DIR <- file.path(
  PROJECT_DIR,
  "data_external_spatial_18J",
  DATASET_ID
)

RAW_DIR <- file.path(BASE_DIR, "raw_zenodo")
UNPACK_DIR <- file.path(BASE_DIR, "unpacked")
STD_DIR <- file.path(BASE_DIR, "standardized")

CANONICAL_ROOT <- file.path(
  PROJECT_DIR,
  "results", "tables",
  "S20_external_spatial_recurrence_CANONICAL_v1.0"
)
OUT_DIR <- file.path(CANONICAL_ROOT, DATASET_ID)

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(UNPACK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STD_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DOWNLOAD_ZENODO <- tolower(
  Sys.getenv("S20_DOWNLOAD_THRANE", unset = "false")
) %in% c("1", "true", "yes")

ZENODO_URL <- "https://zenodo.org/records/8215682/files/Thrane_et_al_2018_CAN_RES.zip?download=1"
ZENODO_ZIP <- file.path(RAW_DIR, "Thrane_et_al_2018_CAN_RES.zip")

EXPECTED_SIZE_MB <- 898.4
MIN_VALID_SIZE_MB <- 850

PRIMARY_X_STATE <- "Myeloid_Treg_Immunosuppressive"
PRIMARY_Y_STATE <- "Tumor_dedifferentiation_Stromal_remodeling"

PRIMARY_X_LABEL <- "Myeloid–Treg immunosuppressive"
PRIMARY_Y_LABEL <- "Tumor-dedifferentiation/stromal-remodeling"

PRIMARY_K <- 6
RUN_MORAN_DESCRIPTIVE <- TRUE

HISTORICAL_GENESET_GSE <- file.path(
  PROJECT_DIR, "results", "tables",
  "revision_external_spatial_recurrence_18J",
  "GSE250636",
  "18J_GSE250636_state_genes_used_long.csv"
)

HISTORICAL_GENESET_THRANE <- file.path(
  PROJECT_DIR, "results", "tables",
  "revision_external_spatial_recurrence_18J",
  "Thrane2018_legacyST",
  "18J_Thrane2018_state_genes_used_long.csv"
)

for (ff in c(HISTORICAL_GENESET_GSE, HISTORICAL_GENESET_THRANE)) {
  if (!file.exists(ff)) {
    stop("Locked historical state-gene file missing: ", ff, call. = FALSE)
  }
}

.gs1 <- utils::read.csv(HISTORICAL_GENESET_GSE, stringsAsFactors = FALSE, check.names = FALSE)
.gs2 <- utils::read.csv(HISTORICAL_GENESET_THRANE, stringsAsFactors = FALSE, check.names = FALSE)
.gs1 <- unique(.gs1[, c("state", "gene"), drop = FALSE])
.gs2 <- unique(.gs2[, c("state", "gene"), drop = FALSE])
.gs1 <- .gs1[order(.gs1$state, .gs1$gene), , drop = FALSE]
.gs2 <- .gs2[order(.gs2$state, .gs2$gene), , drop = FALSE]

if (!identical(.gs1, .gs2)) {
  stop("Historical GSE250636 and Thrane locked state-gene sets are not identical.", call. = FALSE)
}

STATE_GENESET_FILE <- HISTORICAL_GENESET_THRANE

GENESET_SEARCH_ROOTS <- c(
  dirname(STATE_GENESET_FILE)
)

MIN_STATE_GENES_REQUIRED <- 20
MIN_CRITICAL_DETECTED_GENES <- 10
MIN_SPOTS_PER_SAMPLE <- 30
MIN_SCORING_GENE_OVERLAP <- 10

# -----------------------------
# 1. Packages
# -----------------------------

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "stringr",
  "Matrix", "openxlsx", "AnnotationDbi", "org.Hs.eg.db"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Missing required package: ",
      pkg,
      "\nInstall with:\n",
      "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')\n",
      "BiocManager::install(c('AnnotationDbi', 'org.Hs.eg.db'), ask = FALSE, update = FALSE)",
      call. = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Matrix)
  library(openxlsx)
})

# IMPORTANT:
# Do not call library(AnnotationDbi)
# Do not call library(org.Hs.eg.db)
# We use AnnotationDbi::mapIds and org.Hs.eg.db::org.Hs.eg.db directly.

# -----------------------------
# 2. Helpers
# -----------------------------

message2 <- function(...) {
  cat(paste0("[S20-Thrane2018-CANONICAL-v1.0.2] ", paste0(..., collapse = ""), "\n"))
}

zscale <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

clean_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\..*$", "", x)
  x <- toupper(x)
  x <- gsub("^['\"]+|['\"]+$", "", x)
  x
}

is_gene_like <- function(x) {
  x <- clean_gene(x)
  !is.na(x) &
    x != "" &
    x != "NA" &
    !grepl("^[0-9.]+$", x) &
    grepl("^[A-Z0-9][A-Z0-9._-]{1,40}$", x)
}

write_csv <- function(df, filename) {
  data.table::fwrite(df, file.path(OUT_DIR, filename))
}

safe_sheet_name <- function(x) {
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  substr(x, 1, 31)
}

add_sheet <- function(wb, sheet, df) {
  sheet <- safe_sheet_name(sheet)
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df)
  if (ncol(df) > 0) {
    openxlsx::setColWidths(wb, sheet, cols = 1:ncol(df), widths = "auto")
  }
}

parse_gene_cell <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))

  x <- paste(x, collapse = ";")
  x <- gsub("[\r\n\t]+", ";", x)
  x <- unlist(strsplit(x, "[,;|]+|\\s+"))

  x <- clean_gene(x)
  x <- x[is_gene_like(x)]
  unique(x)
}

standardize_state_name <- function(x) {
  xx <- as.character(x)
  xx <- tolower(xx)
  xx <- gsub("[^a-z0-9]+", "_", xx)

  dplyr::case_when(
    grepl("immune", xx) & (grepl("cold", xx) | grepl("defective", xx)) ~
      "Immune_defective_Cold",

    grepl("myeloid", xx) & (grepl("treg", xx) | grepl("immunosuppress", xx)) ~
      "Myeloid_Treg_Immunosuppressive",

    (grepl("dediff", xx) | grepl("dedifferentiation", xx) | grepl("dedifferentiated", xx)) &
      (grepl("stromal", xx) | grepl("remodel", xx)) ~
      "Tumor_dedifferentiation_Stromal_remodeling",

    grepl("melanocytic", xx) |
      (grepl("mitf", xx) & !grepl("dediff|stromal|remodel", xx)) |
      (grepl("differentiation", xx) & !grepl("dediff|stromal|remodel", xx)) ~
      "Melanocytic_Differentiation",

    TRUE ~ NA_character_
  )
}

required_states <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

# -----------------------------
# 3. Robust download and unzip
# -----------------------------

file_size_mb <- function(x) {
  if (!file.exists(x)) return(0)
  round(file.info(x)$size / 1024^2, 2)
}

zip_ok <- function(x) {
  if (!file.exists(x)) return(FALSE)
  ok <- tryCatch({
    utils::unzip(x, list = TRUE)
    TRUE
  }, error = function(e) FALSE)
  ok
}

message2("Current local zip size: ", file_size_mb(ZENODO_ZIP), " MB")

if (DOWNLOAD_ZENODO) {
  if (file.exists(ZENODO_ZIP) &&
      file_size_mb(ZENODO_ZIP) >= MIN_VALID_SIZE_MB &&
      zip_ok(ZENODO_ZIP)) {
    message2("Existing Zenodo zip appears complete and readable; skipping download.")
  } else {
    message2("Downloading / resuming Thrane 2018 Zenodo zip using system curl...")

    curl_exe <- Sys.which("curl")

    if (curl_exe == "") {
      stop(
        "curl.exe not found. Manually download the file from Zenodo and save it as:\n",
        ZENODO_ZIP,
        call. = FALSE
      )
    }

    cmd_args <- c(
      "-L",
      "--retry", "20",
      "--retry-delay", "10",
      "--retry-all-errors",
      "-C", "-",
      "-o", shQuote(ZENODO_ZIP),
      shQuote(ZENODO_URL)
    )

    status <- system2(curl_exe, args = cmd_args)

    message2("Download finished with status: ", status)
    message2("Downloaded zip size: ", file_size_mb(ZENODO_ZIP), " MB")

    if (!file.exists(ZENODO_ZIP) || file_size_mb(ZENODO_ZIP) < MIN_VALID_SIZE_MB) {
      stop(
        "Downloaded file is incomplete. Current size = ",
        file_size_mb(ZENODO_ZIP),
        " MB; expected about ",
        EXPECTED_SIZE_MB,
        " MB. Re-run the script to resume.",
        call. = FALSE
      )
    }

    if (!zip_ok(ZENODO_ZIP)) {
      stop(
        "Downloaded zip exists but cannot be listed by unzip(). Delete and re-download:\n",
        ZENODO_ZIP,
        call. = FALSE
      )
    }
  }
}

existing_unpacked_files <- list.files(
  UNPACK_DIR,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

if (length(existing_unpacked_files) > 0) {
  message2(
    "Existing unpacked Thrane files detected (n=",
    length(existing_unpacked_files),
    "); reusing them without re-unzipping."
  )
} else if (file.exists(ZENODO_ZIP) && zip_ok(ZENODO_ZIP)) {
  message2("No unpacked files detected; unpacking local Zenodo zip...")
  utils::unzip(ZENODO_ZIP, exdir = UNPACK_DIR)
} else {
  stop(
    "No usable Thrane input found. Need either:\n",
    "1) existing files under: ", UNPACK_DIR, "\n",
    "or\n",
    "2) a valid local zip: ", ZENODO_ZIP, "\n",
    "or set S20_DOWNLOAD_THRANE=true to download.",
    call. = FALSE
  )
}

# -----------------------------
# 4. Build file inventory
# -----------------------------

message2("Building file inventory...")

all_files <- list.files(
  UNPACK_DIR,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

all_files <- unique(normalizePath(all_files, winslash = "/", mustWork = FALSE))

if (length(all_files) == 0) {
  stop("No files found after unzip. Check UNPACK_DIR: ", UNPACK_DIR, call. = FALSE)
}

file_inventory <- data.frame(
  file = basename(all_files),
  path = all_files,
  ext = tolower(tools::file_ext(all_files)),
  size_mb = round(file.info(all_files)$size / 1024^2, 3),
  stringsAsFactors = FALSE
)

file_inventory$ext[file_inventory$ext == ""] <- "no_extension"

write_csv(file_inventory, "S20_Thrane2018_file_inventory.csv")

ext_summary <- file_inventory %>%
  dplyr::group_by(ext) %>%
  dplyr::summarise(
    n_files = dplyr::n(),
    total_size_mb = round(sum(size_mb, na.rm = TRUE), 3),
    examples = paste(head(file, 5), collapse = " | "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_files))

write_csv(ext_summary, "S20_Thrane2018_file_extension_summary.csv")

message2("File extension summary:")
print(ext_summary)

# -----------------------------
# 5. Read locked state gene sets
# -----------------------------

read_table_any_sheet <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("csv", "txt", "tsv")) {
    sep <- ifelse(ext == "tsv", "\t", ",")
    return(list(
      sheet1 = as.data.frame(
        data.table::fread(path, sep = sep, check.names = FALSE),
        check.names = FALSE
      )
    ))
  }

  if (ext == "xlsx") {
    sheets <- openxlsx::getSheetNames(path)
    out <- list()
    for (s in sheets) {
      out[[s]] <- tryCatch(
        openxlsx::read.xlsx(path, sheet = s),
        error = function(e) NULL
      )
    }
    out <- out[!vapply(out, is.null, logical(1))]
    return(out)
  }

  stop("Unsupported file extension: ", path)
}

extract_genesets_from_df <- function(df) {
  cn <- colnames(df)

  state_candidates <- intersect(
    c(
      "state", "State", "STATE",
      "FinalState", "final_state", "Final_State",
      "state_name", "State_name", "signature", "Signature",
      "module", "Module", "gene_set", "gene_set_name",
      "state_label", "StateLabel", "state_short", "State_short",
      "StateName", "stateName"
    ),
    cn
  )

  gene_candidates <- intersect(
    c(
      "gene", "Gene", "GENE",
      "symbol", "Symbol", "gene_symbol",
      "GeneSymbol", "Gene.Symbol", "genes", "Genes",
      "gene_list", "Gene_list", "gene_pool", "Gene_pool",
      "retained_genes", "Retained_genes",
      "marker", "markers", "Marker", "Markers",
      "feature", "Feature"
    ),
    cn
  )

  out <- data.frame()

  if (length(state_candidates) > 0 && length(gene_candidates) > 0) {
    for (sc in state_candidates) {
      for (gc in gene_candidates) {
        for (i in seq_len(nrow(df))) {
          st_raw <- as.character(df[[sc]][i])
          st <- standardize_state_name(st_raw)
          genes <- parse_gene_cell(df[[gc]][i])

          if (!is.na(st) && length(genes) > 0) {
            out <- dplyr::bind_rows(
              out,
              data.frame(
                state_raw = st_raw,
                state = st,
                gene = genes,
                parse_mode = paste0("long:", sc, "+", gc),
                stringsAsFactors = FALSE
              )
            )
          }
        }
      }
    }
  }

  possible_state_cols <- cn[
    !is.na(standardize_state_name(cn))
  ]

  if (length(possible_state_cols) > 0) {
    for (cc in possible_state_cols) {
      st <- standardize_state_name(cc)
      genes <- parse_gene_cell(df[[cc]])

      if (!is.na(st) && length(genes) > 0) {
        out <- dplyr::bind_rows(
          out,
          data.frame(
            state_raw = cc,
            state = st,
            gene = genes,
            parse_mode = paste0("wide:", cc),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  out <- unique(out[, c("state_raw", "state", "gene", "parse_mode"), drop = FALSE])
  out
}

read_geneset_file_robust <- function(path) {
  tables <- read_table_any_sheet(path)
  out <- data.frame()

  for (nm in names(tables)) {
    df <- tables[[nm]]
    if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) next

    tmp <- tryCatch(
      extract_genesets_from_df(df),
      error = function(e) data.frame()
    )

    if (nrow(tmp) > 0) {
      tmp$sheet <- nm
      out <- dplyr::bind_rows(out, tmp)
    }
  }

  if (nrow(out) == 0) {
    return(out)
  }

  out <- out[out$state %in% required_states, , drop = FALSE]
  out <- unique(out[, c("state", "gene", "state_raw", "parse_mode", "sheet"), drop = FALSE])
  out
}

find_candidate_geneset_files <- function() {
  if (!is.na(STATE_GENESET_FILE) && file.exists(STATE_GENESET_FILE)) {
    return(STATE_GENESET_FILE)
  }

  roots <- GENESET_SEARCH_ROOTS[dir.exists(GENESET_SEARCH_ROOTS)]

  files <- unlist(lapply(roots, function(r) {
    list.files(
      r,
      recursive = TRUE,
      full.names = TRUE,
      pattern = "\\.(csv|tsv|txt|xlsx)$"
    )
  }))

  files <- unique(files)

  files[
    grepl(
      "state|signature|gene|geneset|gene_set|gene.pool|gene_pool|retained|final|locked|18S|Step14B",
      basename(files),
      ignore.case = TRUE
    )
  ]
}

candidate_files <- find_candidate_geneset_files()

if (length(candidate_files) == 0) {
  stop("No candidate gene-set file found.", call. = FALSE)
}

audit_rows <- data.frame()
parsed_list <- list()

for (f in candidate_files) {
  parsed <- tryCatch(
    read_geneset_file_robust(f),
    error = function(e) data.frame()
  )

  if (nrow(parsed) > 0) {
    cnt <- as.data.frame(table(parsed$state), stringsAsFactors = FALSE)
    colnames(cnt) <- c("state", "n_genes")

    for (st in required_states) {
      nst <- cnt$n_genes[cnt$state == st]
      if (length(nst) == 0) nst <- 0

      audit_rows <- dplyr::bind_rows(
        audit_rows,
        data.frame(
          file = f,
          basename = basename(f),
          state = st,
          n_genes = as.numeric(nst),
          stringsAsFactors = FALSE
        )
      )
    }

    parsed_list[[f]] <- parsed
  } else {
    for (st in required_states) {
      audit_rows <- dplyr::bind_rows(
        audit_rows,
        data.frame(
          file = f,
          basename = basename(f),
          state = st,
          n_genes = 0,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

audit_wide <- tidyr::pivot_wider(
  audit_rows,
  names_from = state,
  values_from = n_genes,
  values_fill = 0
)

for (st in required_states) {
  if (!st %in% colnames(audit_wide)) {
    audit_wide[[st]] <- 0
  }
}

audit_wide$n_states_ok <- rowSums(
  as.data.frame(audit_wide[, required_states, drop = FALSE]) >= MIN_STATE_GENES_REQUIRED
)

audit_wide$total_genes <- rowSums(
  as.data.frame(audit_wide[, required_states, drop = FALSE])
)

audit_wide <- audit_wide[order(-audit_wide$n_states_ok, -audit_wide$total_genes), , drop = FALSE]

write_csv(audit_rows, "S20_Thrane2018_state_geneset_candidate_audit.csv")
write_csv(audit_wide, "S20_Thrane2018_state_geneset_candidate_audit_wide.csv")

best_file <- audit_wide$file[
  audit_wide$n_states_ok == length(required_states)
][1]

if (is.na(best_file) || !best_file %in% names(parsed_list)) {
  stop(
    "No valid state gene-set file found. Check S20_Thrane2018_state_geneset_candidate_audit_wide.csv.",
    call. = FALSE
  )
}

tmp_gs <- parsed_list[[best_file]]

if (!all(c("state", "gene") %in% colnames(tmp_gs))) {
  stop(
    "Parsed gene-set file does not contain required columns: state and gene. Columns found: ",
    paste(colnames(tmp_gs), collapse = ", "),
    call. = FALSE
  )
}

gs_long <- unique(tmp_gs[, c("state", "gene"), drop = FALSE])

state_counts <- as.data.frame(table(gs_long$state), stringsAsFactors = FALSE)
colnames(state_counts) <- c("state", "n_genes")
state_counts <- state_counts[order(state_counts$state), , drop = FALSE]

message2("Selected state gene-set file: ", best_file)
message2("Parsed state gene counts:")
print(state_counts)

write_csv(gs_long, "S20_Thrane2018_state_genes_used_long.csv")
write_csv(state_counts, "S20_Thrane2018_state_genes_used_counts.csv")

state_gene_sets <- split(gs_long$gene, gs_long$state)

# -----------------------------
# 6. Composition modules
# -----------------------------

composition_gene_sets <- list(
  CAF_Stromal = clean_gene(c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2",
    "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "THY1", "FN1", "CXCL12",
    "PDGFRB", "RGS5"
  )),
  Endothelial = clean_gene(c(
    "PECAM1", "VWF", "KDR", "FLT1", "CLDN5", "ESAM", "RAMP2",
    "EMCN", "CDH5", "ACKR1", "PLVAP"
  )),
  Myeloid = clean_gene(c(
    "LST1", "LYZ", "C1QA", "C1QB", "C1QC", "CSF1R", "CD68",
    "FCGR3A", "TYROBP", "AIF1", "ITGAM", "CTSS"
  )),
  T_NK = clean_gene(c(
    "CD3D", "CD3E", "TRAC", "CD2", "IL7R", "CD8A", "CD4",
    "NKG7", "GNLY", "KLRD1", "GZMB", "PRF1"
  )),
  B_Plasma = clean_gene(c(
    "MS4A1", "CD79A", "CD79B", "CD74", "MZB1", "JCHAIN",
    "IGHG1", "IGKC", "XBP1", "SDC1", "BANK1"
  )),
  Melanoma_Lineage = clean_gene(c(
    "MITF", "MLANA", "PMEL", "TYR", "DCT", "SOX10",
    "TYRP1", "GPR143", "S100B", "MIA"
  ))
)

all_gene_sets <- c(state_gene_sets, composition_gene_sets)
all_scoring_genes <- unique(clean_gene(unlist(all_gene_sets)))

# -----------------------------
# 7. Gene ID mapping
# -----------------------------

map_gene_ids_to_symbols <- function(ids) {
  ids_raw <- as.character(ids)
  ids_raw <- trimws(ids_raw)

  ids_base <- sub("\\..*$", "", ids_raw)
  ids_clean <- clean_gene(ids_base)

  out <- rep(NA_character_, length(ids_clean))

  # ----------------------------------------------------------
  # v1.0.2 Thrane compound-ID parser hotfix
  # ----------------------------------------------------------
  # The legacy-ST count matrices store gene identifiers as:
  #   SYMBOL ENSG00000...
  # Example:
  #   ANXA2 ENSG00000182718
  #
  # v1.0.1 relied on a generic token-rescue branch. On the
  # affected Windows/R regex path, the whitespace token was not
  # split reliably, collapsing scoring-gene overlap to zero.
  # Parse this source-locked format explicitly before generic
  # SYMBOL / ENSEMBL / ENTREZ handling.
  compound_idx <- grepl(
    "^[[:space:]]*[^[:space:]]+[[:space:]]+ENSG[0-9]+(?:\\.[0-9]+)?[[:space:]]*$",
    ids_raw,
    perl = TRUE
  )

  if (any(compound_idx)) {
    compound_symbol <- sub(
      "^[[:space:]]*([^[:space:]]+)[[:space:]]+ENSG[0-9]+(?:\\.[0-9]+)?[[:space:]]*$",
      "\\1",
      ids_raw[compound_idx],
      perl = TRUE
    )

    compound_symbol <- clean_gene(compound_symbol)
    compound_symbol[!is_gene_like(compound_symbol)] <- NA_character_
    out[compound_idx] <- compound_symbol
  }

  symbol_like <- is_gene_like(ids_clean) &
    !grepl("^ENSG|^ENST|^ENSP", ids_clean) &
    !grepl("^[0-9]+$", ids_clean)

  direct_symbol_idx <- symbol_like & (is.na(out) | out == "")
  out[direct_symbol_idx] <- ids_clean[direct_symbol_idx]

  ens_idx <- grepl("^ENSG[0-9]+$", ids_clean)

  if (any(ens_idx)) {
    ens_keys <- unique(ids_clean[ens_idx])

    ens_map <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = ens_keys,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

    out[ens_idx] <- clean_gene(unname(ens_map[ids_clean[ens_idx]]))
  }

  entrez_idx <- grepl("^[0-9]+$", ids_clean)

  if (any(entrez_idx)) {
    entrez_keys <- unique(ids_clean[entrez_idx])

    entrez_map <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = entrez_keys,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )

    out[entrez_idx] <- clean_gene(unname(entrez_map[ids_clean[entrez_idx]]))
  }

  still_na <- is.na(out) | out == ""

  if (any(still_na)) {
    tokens <- strsplit(ids_raw[still_na], "[|;,\\s_]+", perl = TRUE)
    rescue <- vapply(tokens, function(v) {
      vv <- clean_gene(v)
      vv <- vv[
        is_gene_like(vv) &
          !grepl("^ENSG|^ENST|^ENSP", vv) &
          !grepl("^[0-9]+$", vv)
      ]

      if (length(vv) > 0) vv[1] else NA_character_
    }, character(1))

    out[still_na] <- rescue
  }

  out <- clean_gene(out)
  out[is.na(out) | out == "" | out == "NA"] <- NA_character_
  out
}

# -----------------------------
# 8. Candidate matrix detection
# -----------------------------

candidate_count_files <- file_inventory %>%
  dplyr::filter(
    ext %in% c("txt", "tsv", "csv", "gz"),
    !grepl("annotation|metadata|image|spot_image|aligned|he|h&e|json|readme", file, ignore.case = TRUE),
    grepl("count|counts|matrix|expression|expr|stdata|data", file, ignore.case = TRUE)
  ) %>%
  dplyr::arrange(dplyr::desc(size_mb))

if (nrow(candidate_count_files) == 0) {
  candidate_count_files <- file_inventory %>%
    dplyr::filter(ext %in% c("txt", "tsv", "csv", "gz")) %>%
    dplyr::arrange(dplyr::desc(size_mb))
}

write_csv(candidate_count_files, "S20_Thrane2018_candidate_count_files.csv")

message2("Candidate count files:")
print(head(candidate_count_files[, c("file", "ext", "size_mb", "path")], 20))

# -----------------------------
# 9. Matrix reader with orientation + ID mapping
# -----------------------------

read_expression_matrix_legacy <- function(path) {
  message2("Trying to read matrix: ", basename(path))

  sep <- ifelse(grepl("\\.csv(\\.gz)?$", basename(path), ignore.case = TRUE), ",", "\t")

  df <- tryCatch(
    as.data.frame(
      data.table::fread(
        path,
        sep = sep,
        check.names = FALSE,
        data.table = FALSE
      ),
      check.names = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(df) || nrow(df) < 20 || ncol(df) < 20) {
    return(list(expr = NULL, diag = data.frame()))
  }

  cn <- colnames(df)
  scoring_genes <- all_scoring_genes

  numeric_score <- function(x) {
    suppressWarnings({
      y <- as.numeric(x)
      mean(is.finite(y))
    })
  }

  numeric_frac_cols <- vapply(df, numeric_score, numeric(1))
  numeric_cols <- which(numeric_frac_cols > 0.80)

  candidate_mats <- list()
  diag_rows <- data.frame()

  # Candidate A: first column = gene IDs / gene symbols, columns = spots
  if (length(numeric_cols) >= 20) {
    id_col <- 1
    num_cols <- setdiff(numeric_cols, id_col)

    if (length(num_cols) >= 20) {
      ids_a <- as.character(df[[id_col]])
      symbols_a <- map_gene_ids_to_symbols(ids_a)

      compound_a <- grepl(
        "^[[:space:]]*[^[:space:]]+[[:space:]]+ENSG[0-9]+(?:\\.[0-9]+)?[[:space:]]*$",
        ids_a,
        perl = TRUE
      )

      if (any(compound_a)) {
        message2(
          "  compound SYMBOL+ENSG IDs detected=", sum(compound_a),
          " | parsed symbols=", sum(!is.na(symbols_a[compound_a]) & symbols_a[compound_a] != "")
        )
      }

      mat_a <- tryCatch(
        as.matrix(
          data.frame(
            lapply(df[, num_cols, drop = FALSE], function(x) suppressWarnings(as.numeric(x))),
            check.names = FALSE
          )
        ),
        error = function(e) NULL
      )

      if (!is.null(mat_a)) {
        keep_gene <- !is.na(symbols_a) &
          symbols_a != "" &
          rowSums(is.finite(mat_a)) > 0

        mat_a <- mat_a[keep_gene, , drop = FALSE]
        symbols_a <- symbols_a[keep_gene]

        if (nrow(mat_a) >= 20 && ncol(mat_a) >= 20) {
          rownames(mat_a) <- symbols_a
          colnames(mat_a) <- cn[num_cols]

          if (any(duplicated(rownames(mat_a)))) {
            mat_a <- rowsum(mat_a, group = rownames(mat_a), reorder = FALSE)
          }

          overlap_a <- length(intersect(unique(rownames(mat_a)), scoring_genes))

          candidate_mats[["genes_rows"]] <- list(
            mat = mat_a,
            orientation = "genes_in_rows_mapped",
            n_genes = nrow(mat_a),
            n_spots = ncol(mat_a),
            scoring_gene_overlap = overlap_a
          )

          diag_rows <- dplyr::bind_rows(
            diag_rows,
            data.frame(
              file = basename(path),
              orientation = "genes_in_rows_mapped",
              n_genes = nrow(mat_a),
              n_spots = ncol(mat_a),
              scoring_gene_overlap = overlap_a,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }

  # Candidate B: first column = spots, columns = gene IDs / gene symbols
  id_col <- 1
  possible_gene_cols <- setdiff(seq_len(ncol(df)), id_col)

  gene_ids_b <- cn[possible_gene_cols]
  symbols_b <- map_gene_ids_to_symbols(gene_ids_b)

  numeric_gene_cols <- possible_gene_cols[
    !is.na(symbols_b) &
      symbols_b != "" &
      numeric_frac_cols[possible_gene_cols] > 0.80
  ]

  if (length(numeric_gene_cols) >= 100) {
    spot_ids <- as.character(df[[id_col]])
    symbols_use <- map_gene_ids_to_symbols(cn[numeric_gene_cols])

    mat0 <- tryCatch(
      as.matrix(
        data.frame(
          lapply(df[, numeric_gene_cols, drop = FALSE], function(x) suppressWarnings(as.numeric(x))),
          check.names = FALSE
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(mat0)) {
      keep_spot <- !is.na(spot_ids) &
        spot_ids != "" &
        spot_ids != "NA" &
        rowSums(is.finite(mat0)) > 0

      keep_gene <- !is.na(symbols_use) &
        symbols_use != "" &
        symbols_use != "NA"

      mat0 <- mat0[keep_spot, keep_gene, drop = FALSE]
      spot_ids <- spot_ids[keep_spot]
      symbols_use <- symbols_use[keep_gene]

      if (nrow(mat0) >= 20 && ncol(mat0) >= 100) {
        mat_b <- t(mat0)
        rownames(mat_b) <- symbols_use
        colnames(mat_b) <- spot_ids

        if (any(duplicated(rownames(mat_b)))) {
          mat_b <- rowsum(mat_b, group = rownames(mat_b), reorder = FALSE)
        }

        overlap_b <- length(intersect(unique(rownames(mat_b)), scoring_genes))

        candidate_mats[["genes_columns_transposed"]] <- list(
          mat = mat_b,
          orientation = "genes_in_columns_transposed_mapped",
          n_genes = nrow(mat_b),
          n_spots = ncol(mat_b),
          scoring_gene_overlap = overlap_b
        )

        diag_rows <- dplyr::bind_rows(
          diag_rows,
          data.frame(
            file = basename(path),
            orientation = "genes_columns_transposed_mapped",
            n_genes = nrow(mat_b),
            n_spots = ncol(mat_b),
            scoring_gene_overlap = overlap_b,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  if (length(candidate_mats) == 0) {
    return(list(expr = NULL, diag = diag_rows))
  }

  overlaps <- vapply(candidate_mats, function(x) x$scoring_gene_overlap, numeric(1))
  best_name <- names(which.max(overlaps))
  best <- candidate_mats[[best_name]]

  message2(
    "  selected orientation: ", best$orientation,
    " | genes=", best$n_genes,
    " spots=", best$n_spots,
    " scoring_gene_overlap=", best$scoring_gene_overlap
  )

  if (best$scoring_gene_overlap < MIN_SCORING_GENE_OVERLAP) {
    message2("  WARNING: low scoring gene overlap after ID mapping; skipping this matrix.")
    return(list(expr = NULL, diag = diag_rows))
  }

  mat <- best$mat
  mat[!is.finite(mat)] <- 0

  list(expr = as(mat, "dgCMatrix"), diag = diag_rows)
}

extract_xy_from_spot <- function(spot_id) {
  s <- as.character(spot_id)
  nums <- stringr::str_extract_all(s, "-?[0-9]+\\.?[0-9]*")

  out <- lapply(nums, function(v) {
    if (length(v) >= 2) {
      vv <- suppressWarnings(as.numeric(v))
      c(x = vv[length(vv) - 1], y = vv[length(vv)])
    } else {
      c(x = NA_real_, y = NA_real_)
    }
  })

  mat <- do.call(rbind, out)

  data.frame(
    spot_barcode = spot_id,
    spatial_x = as.numeric(mat[, "x"]),
    spatial_y = as.numeric(mat[, "y"]),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 10. Module scoring
# -----------------------------

score_modules <- function(counts, gene_sets) {
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }

  row_clean <- map_gene_ids_to_symbols(rownames(counts))
  all_needed <- unique(clean_gene(unlist(gene_sets)))

  rows_needed <- which(!is.na(row_clean) & row_clean %in% all_needed)

  if (length(rows_needed) == 0) {
    stop("No scoring genes detected in count matrix.", call. = FALSE)
  }

  lib <- Matrix::colSums(counts)
  lib[lib <= 0] <- NA_real_

  sub <- counts[rows_needed, , drop = FALSE]
  sub_gene <- row_clean[rows_needed]

  sub_dense <- as.matrix(sub)

  if (any(duplicated(sub_gene))) {
    sub_dense <- rowsum(sub_dense, group = sub_gene, reorder = FALSE)
  } else {
    rownames(sub_dense) <- sub_gene
  }

  norm <- t(t(sub_dense) / lib * 1e4)
  log_norm <- log1p(norm)

  z <- t(scale(t(log_norm)))
  z[!is.finite(z)] <- NA_real_

  score_df <- data.frame(
    spot_barcode = colnames(counts),
    stringsAsFactors = FALSE
  )

  coverage <- data.frame()

  for (nm in names(gene_sets)) {
    genes <- unique(clean_gene(gene_sets[[nm]]))
    hit <- intersect(genes, rownames(z))

    coverage <- dplyr::bind_rows(
      coverage,
      data.frame(
        module = nm,
        n_genes_input = length(genes),
        n_genes_detected = length(hit),
        coverage_fraction = length(hit) / length(genes),
        detected_genes = paste(hit, collapse = ";"),
        missing_genes = paste(setdiff(genes, hit), collapse = ";"),
        stringsAsFactors = FALSE
      )
    )

    if (length(hit) < 3) {
      score_df[[nm]] <- NA_real_
    } else {
      score_df[[nm]] <- colMeans(z[hit, , drop = FALSE], na.rm = TRUE)
    }
  }

  list(scores = score_df, coverage = coverage)
}

# -----------------------------
# 11. Read, score, and collect samples
# -----------------------------

all_scores <- data.frame()
all_coverage <- data.frame()
sample_input_log <- data.frame()
matrix_read_log <- data.frame()
matrix_candidate_log <- data.frame()

used_sample_names <- character(0)

for (i in seq_len(nrow(candidate_count_files))) {
  f <- candidate_count_files$path[i]
  sample_raw <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(f)))
  sample_id <- gsub("[^A-Za-z0-9]+", "_", sample_raw)

  if (sample_id %in% used_sample_names) next

  rr <- read_expression_matrix_legacy(f)
  expr <- rr$expr

  if (!is.null(rr$diag) && nrow(rr$diag) > 0) {
    rr$diag$sample_id <- sample_id
    rr$diag$path <- f
    matrix_candidate_log <- dplyr::bind_rows(matrix_candidate_log, rr$diag)
  }

  if (is.null(expr)) {
    matrix_read_log <- dplyr::bind_rows(
      matrix_read_log,
      data.frame(
        file = basename(f),
        path = f,
        sample_id = sample_id,
        readable_expression = FALSE,
        n_genes = NA_integer_,
        n_spots = NA_integer_,
        reason = "not_detected_or_low_gene_overlap",
        stringsAsFactors = FALSE
      )
    )
    next
  }

  message2("Readable matrix: ", basename(f), " | genes=", nrow(expr), " spots=", ncol(expr))

  if (ncol(expr) < MIN_SPOTS_PER_SAMPLE) {
    matrix_read_log <- dplyr::bind_rows(
      matrix_read_log,
      data.frame(
        file = basename(f),
        path = f,
        sample_id = sample_id,
        readable_expression = TRUE,
        n_genes = nrow(expr),
        n_spots = ncol(expr),
        reason = "too_few_spots",
        stringsAsFactors = FALSE
      )
    )
    next
  }

  sc <- score_modules(expr, all_gene_sets)

  pos <- extract_xy_from_spot(colnames(expr))

  scores <- sc$scores %>%
    dplyr::mutate(
      dataset = DATASET_ID,
      sample_id = sample_id,
      source_file = basename(f)
    ) %>%
    dplyr::left_join(pos, by = "spot_barcode") %>%
    dplyr::mutate(
      spot_id = paste(sample_id, spot_barcode, sep = "__")
    )

  cov <- sc$coverage %>%
    dplyr::mutate(
      dataset = DATASET_ID,
      sample_id = sample_id,
      source_file = basename(f)
    )

  all_scores <- dplyr::bind_rows(all_scores, scores)
  all_coverage <- dplyr::bind_rows(all_coverage, cov)

  sample_input_log <- dplyr::bind_rows(
    sample_input_log,
    data.frame(
      dataset = DATASET_ID,
      sample_id = sample_id,
      source_file = basename(f),
      path = f,
      n_genes = nrow(expr),
      n_spots = ncol(expr),
      n_modules_scored = length(all_gene_sets),
      coordinate_inference = ifelse(
        mean(is.finite(pos$spatial_x) & is.finite(pos$spatial_y)) > 0.8,
        "from_spot_ids",
        "mostly_unavailable"
      ),
      coordinate_fraction = mean(is.finite(pos$spatial_x) & is.finite(pos$spatial_y)),
      stringsAsFactors = FALSE
    )
  )

  matrix_read_log <- dplyr::bind_rows(
    matrix_read_log,
    data.frame(
      file = basename(f),
      path = f,
      sample_id = sample_id,
      readable_expression = TRUE,
      n_genes = nrow(expr),
      n_spots = ncol(expr),
      reason = "used",
      stringsAsFactors = FALSE
    )
  )

  used_sample_names <- c(used_sample_names, sample_id)

  saveRDS(
    list(
      counts_dim = dim(expr),
      scores = scores,
      coverage = cov
    ),
    file.path(STD_DIR, paste0(sample_id, "_18J_Thrane_scores.rds"))
  )

  rm(expr, rr, sc, pos, scores, cov)
  gc()
}

write_csv(matrix_read_log, "S20_Thrane2018_matrix_read_log.csv")
write_csv(matrix_candidate_log, "S20_Thrane2018_matrix_candidate_orientation_log.csv")
write_csv(sample_input_log, "S20_Thrane2018_sample_input_log.csv")
write_csv(all_coverage, "S20_Thrane2018_gene_coverage.csv")
write_csv(all_scores, "S20_Thrane2018_spot_state_composition_scores.csv")

if (nrow(all_scores) == 0) {
  stop(
    "No usable expression matrix was detected. Check S20_Thrane2018_matrix_candidate_orientation_log.csv. ",
    "If all scoring_gene_overlap values are still 0, inspect gene IDs in the first column.",
    call. = FALSE
  )
}

coverage_summary <- all_coverage %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(
    min_detected = min(n_genes_detected, na.rm = TRUE),
    median_detected = median(n_genes_detected, na.rm = TRUE),
    max_detected = max(n_genes_detected, na.rm = TRUE),
    min_coverage = min(coverage_fraction, na.rm = TRUE),
    median_coverage = median(coverage_fraction, na.rm = TRUE),
    max_coverage = max(coverage_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(module)

write_csv(coverage_summary, "S20_Thrane2018_gene_coverage_summary.csv")

message2("Gene coverage summary:")
print(coverage_summary)

critical_bad <- coverage_summary %>%
  dplyr::filter(
    module %in% c(PRIMARY_X_STATE, PRIMARY_Y_STATE),
    median_detected < MIN_CRITICAL_DETECTED_GENES
  )

if (nrow(critical_bad) > 0) {
  stop(
    "Critical primary-pair signatures have low detected-gene coverage. Do not use Thrane unless gene symbols are fixed.",
    call. = FALSE
  )
}

# -----------------------------
# 12. Derived composition variables
# -----------------------------

df <- all_scores

df$Stromal_Vascular_Composite <- rowMeans(
  df[, c("CAF_Stromal", "Endothelial"), drop = FALSE],
  na.rm = TRUE
)

df$Tumor_Purity_Proxy <- zscale(df$Melanoma_Lineage) -
  rowMeans(
    cbind(
      zscale(df$CAF_Stromal),
      zscale(df$Endothelial),
      zscale(df$Myeloid),
      zscale(df$T_NK),
      zscale(df$B_Plasma)
    ),
    na.rm = TRUE
  )

comp_cols <- c(
  "CAF_Stromal", "Endothelial", "Myeloid",
  "T_NK", "B_Plasma", "Melanoma_Lineage"
)

pc_df <- df[, comp_cols, drop = FALSE]
pc_df[] <- lapply(pc_df, as.numeric)

pc_ok <- complete.cases(pc_df)

df$MarkerComp_PC1 <- NA_real_
df$MarkerComp_PC2 <- NA_real_

if (sum(pc_ok) > 50) {
  pc <- prcomp(scale(pc_df[pc_ok, , drop = FALSE]), center = TRUE, scale. = FALSE)
  df$MarkerComp_PC1[pc_ok] <- pc$x[, 1]
  df$MarkerComp_PC2[pc_ok] <- pc$x[, 2]
}

# -----------------------------
# 13. Residualization
# -----------------------------

residualize_score <- function(df, score_col, covariates, sample_col = "sample_id") {
  y <- as.numeric(df[[score_col]])
  out <- rep(NA_real_, length(y))

  covariates <- covariates[covariates %in% colnames(df)]

  covariates <- covariates[vapply(df[covariates], function(x) {
    x <- as.numeric(x)
    stats::sd(x, na.rm = TRUE) > 0 && sum(is.finite(x)) > 20
  }, logical(1))]

  use <- is.finite(y)

  for (cv in covariates) {
    use <- use & is.finite(as.numeric(df[[cv]]))
  }

  dat <- data.frame(y = y[use])

  for (cv in covariates) {
    dat[[cv]] <- as.numeric(df[[cv]][use])
  }

  include_sample <- FALSE

  if (sample_col %in% colnames(df)) {
    smp <- as.factor(df[[sample_col]][use])
    if (length(unique(smp)) > 1 && min(table(smp)) >= 20) {
      dat$sample_factor <- smp
      include_sample <- TRUE
    }
  }

  if (length(covariates) == 0 && !include_sample) {
    out[use] <- zscale(y[use])
    return(out)
  }

  rhs <- c(covariates, if (include_sample) "sample_factor" else NULL)

  if (length(rhs) == 0) {
    out[use] <- zscale(y[use])
    return(out)
  }

  form <- as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  fit <- tryCatch(lm(form, data = dat), error = function(e) NULL)

  if (is.null(fit)) {
    out[use] <- zscale(y[use])
  } else {
    out[use] <- zscale(residuals(fit))
  }

  out
}

adjustment_sets <- list(
  Raw_score = character(0),
  Stromal_vascular_residualized = c("CAF_Stromal", "Endothelial"),
  Broad_composition_residualized = c(
    "CAF_Stromal", "Endothelial", "Myeloid",
    "T_NK", "B_Plasma", "Melanoma_Lineage"
  ),
  Tumor_purity_proxy_residualized = c("Tumor_Purity_Proxy"),
  Composition_PC1_PC2_residualized = c("MarkerComp_PC1", "MarkerComp_PC2")
)

for (adj in names(adjustment_sets)) {
  for (st in c(PRIMARY_X_STATE, PRIMARY_Y_STATE)) {
    new_col <- paste(st, adj, sep = "__")

    if (adj == "Raw_score") {
      df[[new_col]] <- zscale(df[[st]])
    } else {
      df[[new_col]] <- residualize_score(
        df = df,
        score_col = st,
        covariates = adjustment_sets[[adj]],
        sample_col = "sample_id"
      )
    }
  }
}

write_csv(df, "S20_Thrane2018_spot_scores_with_residualized_primary_pair.csv")

# -----------------------------
# 14. State-composition correlations
# -----------------------------

state_cols <- required_states

composition_cols <- c(
  "CAF_Stromal", "Endothelial", "Myeloid",
  "T_NK", "B_Plasma", "Melanoma_Lineage",
  "Stromal_Vascular_Composite", "Tumor_Purity_Proxy"
)

cor_rows <- data.frame(
  dataset = character(),
  sample_id = character(),
  state = character(),
  composition_module = character(),
  spearman_rho = numeric(),
  p_value = numeric(),
  n_spots = integer(),
  stringsAsFactors = FALSE
)

for (smp in c(unique(df$sample_id), "pooled")) {
  dd <- if (smp == "pooled") df else df[df$sample_id == smp, , drop = FALSE]

  for (st in state_cols) {
    for (cp in composition_cols) {
      ok <- is.finite(dd[[st]]) & is.finite(dd[[cp]])

      if (sum(ok) >= 20) {
        ct <- suppressWarnings(cor.test(dd[[st]][ok], dd[[cp]][ok], method = "spearman"))

        cor_rows <- dplyr::bind_rows(
          cor_rows,
          data.frame(
            dataset = DATASET_ID,
            sample_id = smp,
            state = st,
            composition_module = cp,
            spearman_rho = unname(ct$estimate),
            p_value = ct$p.value,
            n_spots = sum(ok),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

write_csv(cor_rows, "S20_Thrane2018_state_composition_correlation.csv")

# -----------------------------
# 15. Primary pair correlation
# -----------------------------

primary_cor <- data.frame(
  dataset = character(),
  sample_id = character(),
  x_state = character(),
  y_state = character(),
  adjustment = character(),
  spearman_rho = numeric(),
  p_value = numeric(),
  n_spots = integer(),
  stringsAsFactors = FALSE
)

for (adj in names(adjustment_sets)) {
  x_col <- paste(PRIMARY_X_STATE, adj, sep = "__")
  y_col <- paste(PRIMARY_Y_STATE, adj, sep = "__")

  for (smp in c(unique(df$sample_id), "pooled")) {
    dd <- if (smp == "pooled") df else df[df$sample_id == smp, , drop = FALSE]
    ok <- is.finite(dd[[x_col]]) & is.finite(dd[[y_col]])

    if (sum(ok) >= 20) {
      ct <- suppressWarnings(cor.test(dd[[x_col]][ok], dd[[y_col]][ok], method = "spearman"))

      primary_cor <- dplyr::bind_rows(
        primary_cor,
        data.frame(
          dataset = DATASET_ID,
          sample_id = smp,
          x_state = PRIMARY_X_STATE,
          y_state = PRIMARY_Y_STATE,
          adjustment = adj,
          spearman_rho = unname(ct$estimate),
          p_value = ct$p.value,
          n_spots = sum(ok),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

write_csv(primary_cor, "S20_Thrane2018_primary_pair_raw_vs_residualized_correlation.csv")

# -----------------------------
# 16. Dual-high OR
# -----------------------------

fisher_dual_high <- function(x, y, group = NULL, q = 0.75) {
  ok <- is.finite(x) & is.finite(y)

  x <- x[ok]
  y <- y[ok]

  if (!is.null(group)) {
    group <- group[ok]
  }

  if (length(x) < 20) return(data.frame())

  if (is.null(group)) {
    x_cut <- quantile(x, q, na.rm = TRUE)
    y_cut <- quantile(y, q, na.rm = TRUE)

    xh <- x >= x_cut
    yh <- y >= y_cut
  } else {
    xh <- rep(FALSE, length(x))
    yh <- rep(FALSE, length(y))

    for (g in unique(group)) {
      idx <- which(group == g)

      if (length(idx) >= 20) {
        xh[idx] <- x[idx] >= quantile(x[idx], q, na.rm = TRUE)
        yh[idx] <- y[idx] >= quantile(y[idx], q, na.rm = TRUE)
      }
    }
  }

  tab <- matrix(
    c(
      sum(xh & yh),
      sum(xh & !yh),
      sum(!xh & yh),
      sum(!xh & !yh)
    ),
    nrow = 2,
    byrow = TRUE
  )

  ft <- suppressWarnings(fisher.test(tab))

  data.frame(
    n_spots = length(x),
    BothHigh = tab[1, 1],
    XHighOnly = tab[1, 2],
    YHighOnly = tab[2, 1],
    NeitherHigh = tab[2, 2],
    OR = unname(ft$estimate),
    CI_low = ft$conf.int[1],
    CI_high = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

dual_rows <- data.frame(
  dataset = character(),
  sample_id = character(),
  x_state = character(),
  y_state = character(),
  adjustment = character(),
  cutoff_definition = character(),
  n_spots = integer(),
  BothHigh = integer(),
  XHighOnly = integer(),
  YHighOnly = integer(),
  NeitherHigh = integer(),
  OR = numeric(),
  CI_low = numeric(),
  CI_high = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (adj in names(adjustment_sets)) {
  x_col <- paste(PRIMARY_X_STATE, adj, sep = "__")
  y_col <- paste(PRIMARY_Y_STATE, adj, sep = "__")

  for (smp in unique(df$sample_id)) {
    dd <- df[df$sample_id == smp, , drop = FALSE]

    res <- fisher_dual_high(dd[[x_col]], dd[[y_col]], group = NULL, q = 0.75)

    if (nrow(res) > 0) {
      dual_rows <- dplyr::bind_rows(
        dual_rows,
        data.frame(
          dataset = DATASET_ID,
          sample_id = smp,
          x_state = PRIMARY_X_STATE,
          y_state = PRIMARY_Y_STATE,
          adjustment = adj,
          cutoff_definition = "within_sample_top_quartile",
          res,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  res <- fisher_dual_high(df[[x_col]], df[[y_col]], group = df$sample_id, q = 0.75)

  if (nrow(res) > 0) {
    dual_rows <- dplyr::bind_rows(
      dual_rows,
      data.frame(
        dataset = DATASET_ID,
        sample_id = "pooled",
        x_state = PRIMARY_X_STATE,
        y_state = PRIMARY_Y_STATE,
        adjustment = adj,
        cutoff_definition = "within_sample_top_quartile_pooled",
        res,
        stringsAsFactors = FALSE
      )
    )
  }
}

write_csv(dual_rows, "S20_Thrane2018_primary_pair_dual_high_OR.csv")

# -----------------------------
# 17. Descriptive bivariate Moran-type statistic
# -----------------------------

make_knn_cache <- function(dd, k = 6) {
  coords <- as.matrix(dd[, c("spatial_x", "spatial_y")])
  ok <- is.finite(coords[, 1]) & is.finite(coords[, 2])
  idx <- which(ok)

  if (length(idx) < k + 5) return(NULL)

  d <- as.matrix(stats::dist(coords[idx, , drop = FALSE]))
  diag(d) <- Inf

  nn <- t(apply(d, 1, function(z) order(z)[seq_len(min(k, length(z) - 1))]))

  list(idx = idx, nn = nn)
}

bivar_moran <- function(x, y, cache) {
  if (is.null(cache)) return(NA_real_)

  idx <- cache$idx
  nn <- cache$nn

  xs <- zscale(x[idx])
  ys <- zscale(y[idx])

  xs[!is.finite(xs)] <- 0
  ys[!is.finite(ys)] <- 0

  lag_y <- rowMeans(matrix(ys[nn], nrow = nrow(nn)), na.rm = TRUE)
  lag_x <- rowMeans(matrix(xs[nn], nrow = nrow(nn)), na.rm = TRUE)

  I_xy <- mean(xs * lag_y, na.rm = TRUE)
  I_yx <- mean(ys * lag_x, na.rm = TRUE)

  mean(c(I_xy, I_yx), na.rm = TRUE)
}

moran_rows <- data.frame(
  dataset = character(),
  sample_id = character(),
  x_state = character(),
  y_state = character(),
  adjustment = character(),
  k = integer(),
  I_symmetric = numeric(),
  n_spots = integer(),
  stringsAsFactors = FALSE
)

if (RUN_MORAN_DESCRIPTIVE) {
  for (adj in names(adjustment_sets)) {
    x_col <- paste(PRIMARY_X_STATE, adj, sep = "__")
    y_col <- paste(PRIMARY_Y_STATE, adj, sep = "__")

    for (smp in unique(df$sample_id)) {
      dd <- df[df$sample_id == smp, , drop = FALSE]
      cache <- make_knn_cache(dd, k = PRIMARY_K)
      I <- bivar_moran(dd[[x_col]], dd[[y_col]], cache)

      moran_rows <- dplyr::bind_rows(
        moran_rows,
        data.frame(
          dataset = DATASET_ID,
          sample_id = smp,
          x_state = PRIMARY_X_STATE,
          y_state = PRIMARY_Y_STATE,
          adjustment = adj,
          k = PRIMARY_K,
          I_symmetric = I,
          n_spots = nrow(dd),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  if (nrow(moran_rows) > 0) {
    moran_pooled <- moran_rows %>%
      dplyr::group_by(adjustment) %>%
      dplyr::summarise(
        dataset = DATASET_ID,
        sample_id = "pooled_weighted",
        x_state = PRIMARY_X_STATE,
        y_state = PRIMARY_Y_STATE,
        k = PRIMARY_K,
        I_symmetric = weighted.mean(I_symmetric, w = n_spots, na.rm = TRUE),
        n_spots = sum(n_spots, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::select(
        dataset,
        sample_id,
        x_state,
        y_state,
        adjustment,
        k,
        I_symmetric,
        n_spots
      )

    moran_rows <- dplyr::bind_rows(moran_rows, moran_pooled)
  }
}

write_csv(moran_rows, "S20_Thrane2018_primary_pair_Moran_descriptive.csv")

# -----------------------------
# 18. Key readout
# -----------------------------

key_readout <- data.frame(
  readout = character(),
  dataset = character(),
  sample_id = character(),
  adjustment = character(),
  statistic_label = character(),
  statistic = numeric(),
  lower = numeric(),
  upper = numeric(),
  p_value = numeric(),
  n_spots = integer(),
  stringsAsFactors = FALSE
)

if (nrow(primary_cor) > 0) {
  key_readout <- dplyr::bind_rows(
    key_readout,
    primary_cor %>%
      dplyr::filter(sample_id == "pooled") %>%
      dplyr::mutate(
        readout = "primary_pair_spearman",
        statistic_label = "Spearman rho",
        statistic = spearman_rho,
        lower = NA_real_,
        upper = NA_real_
      ) %>%
      dplyr::select(
        readout,
        dataset,
        sample_id,
        adjustment,
        statistic_label,
        statistic,
        lower,
        upper,
        p_value,
        n_spots
      )
  )
}

if (nrow(dual_rows) > 0) {
  key_readout <- dplyr::bind_rows(
    key_readout,
    dual_rows %>%
      dplyr::filter(sample_id == "pooled") %>%
      dplyr::mutate(
        readout = "primary_pair_dual_high_OR",
        statistic_label = "Fisher OR",
        statistic = OR,
        lower = CI_low,
        upper = CI_high
      ) %>%
      dplyr::select(
        readout,
        dataset,
        sample_id,
        adjustment,
        statistic_label,
        statistic,
        lower,
        upper,
        p_value,
        n_spots
      )
  )
}

if (nrow(moran_rows) > 0) {
  key_readout <- dplyr::bind_rows(
    key_readout,
    moran_rows %>%
      dplyr::filter(sample_id == "pooled_weighted") %>%
      dplyr::mutate(
        readout = "primary_pair_bivariate_Moran_descriptive",
        statistic_label = "Bivariate Moran-type I",
        statistic = I_symmetric,
        lower = NA_real_,
        upper = NA_real_,
        p_value = NA_real_
      ) %>%
      dplyr::select(
        readout,
        dataset,
        sample_id,
        adjustment,
        statistic_label,
        statistic,
        lower,
        upper,
        p_value,
        n_spots
      )
  )
}

write_csv(key_readout, "S20_Thrane2018_key_readout.csv")

if (nrow(key_readout) == 0) {
  stop("Key readout is empty. Check matrix_read_log and gene coverage.", call. = FALSE)
}

# -----------------------------
# 19. Inclusion decision helper
# -----------------------------

raw_rho <- key_readout %>%
  dplyr::filter(readout == "primary_pair_spearman", adjustment == "Raw_score") %>%
  dplyr::pull(statistic)

broad_rho <- key_readout %>%
  dplyr::filter(readout == "primary_pair_spearman", adjustment == "Broad_composition_residualized") %>%
  dplyr::pull(statistic)

raw_or <- key_readout %>%
  dplyr::filter(readout == "primary_pair_dual_high_OR", adjustment == "Raw_score") %>%
  dplyr::pull(statistic)

broad_or <- key_readout %>%
  dplyr::filter(readout == "primary_pair_dual_high_OR", adjustment == "Broad_composition_residualized") %>%
  dplyr::pull(statistic)

include_decision <- data.frame(
  criterion = c(
    "critical_signature_coverage_ok",
    "raw_rho_positive",
    "raw_dual_high_OR_above_1",
    "composition_attenuation_rho",
    "composition_attenuation_OR",
    "recommended_use"
  ),
  value = NA_character_,
  stringsAsFactors = FALSE
)

include_decision$value[1] <- as.character(
  all(coverage_summary$median_detected[
    coverage_summary$module %in% c(PRIMARY_X_STATE, PRIMARY_Y_STATE)
  ] >= MIN_CRITICAL_DETECTED_GENES)
)

include_decision$value[2] <- as.character(
  length(raw_rho) == 1 && is.finite(raw_rho) && raw_rho > 0
)

include_decision$value[3] <- as.character(
  length(raw_or) == 1 && is.finite(raw_or) && raw_or > 1
)

include_decision$value[4] <- as.character(
  length(raw_rho) == 1 && length(broad_rho) == 1 &&
    is.finite(raw_rho) && is.finite(broad_rho) && abs(broad_rho) < abs(raw_rho)
)

include_decision$value[5] <- as.character(
  length(raw_or) == 1 && length(broad_or) == 1 &&
    is.finite(raw_or) && is.finite(broad_or) && broad_or < raw_or
)

recommended <- all(include_decision$value[1:5] == "TRUE")

include_decision$value[6] <- ifelse(
  recommended,
  "Candidate for manuscript as external legacy-ST recurrence audit",
  "Do not include as positive recurrence; keep internal or report as non-uniform only"
)

write_csv(include_decision, "S20_Thrane2018_inclusion_decision_helper.csv")

# -----------------------------
# 20. Workbook
# -----------------------------

wb <- openxlsx::createWorkbook()

readme <- data.frame(
  Item = c(
    "Dataset",
    "Purpose",
    "Primary pair",
    "Response annotation",
    "Interpretation boundary",
    "Input format",
    "Selected gene-set file",
    "Main outputs"
  ),
  Description = c(
    DATASET_ID,
    "External melanoma legacy-ST recurrence audit of composition-linked state-score co-variation.",
    paste(PRIMARY_X_LABEL, "vs", PRIMARY_Y_LABEL),
    "No matched ICB-response annotation was used or assumed.",
    "This analysis is not response validation, not biomarker validation, and not functional niche validation.",
    "Thrane et al. 2018 legacy spatial transcriptomics count matrices downloaded from Zenodo; gene IDs mapped to HGNC symbols before scoring.",
    best_file,
    "State-composition correlations, raw/residualized primary-pair correlations, dual-high odds ratios, and optional descriptive Moran-type statistics."
  ),
  stringsAsFactors = FALSE
)

add_sheet(wb, "README", readme)
add_sheet(wb, "Key_readout", key_readout)
add_sheet(wb, "Inclusion_decision", include_decision)
add_sheet(wb, "Sample_input_log", sample_input_log)
add_sheet(wb, "Matrix_read_log", matrix_read_log)
add_sheet(wb, "Matrix_candidate_log", matrix_candidate_log)
add_sheet(wb, "Gene_set_audit", audit_wide)
add_sheet(wb, "State_gene_counts", state_counts)
add_sheet(wb, "Gene_coverage_summary", coverage_summary)
add_sheet(wb, "Gene_coverage", all_coverage)
add_sheet(wb, "State_comp_cor", cor_rows)
add_sheet(wb, "Primary_pair_cor", primary_cor)
add_sheet(wb, "Dual_high_OR", dual_rows)
add_sheet(wb, "Moran_descriptive", moran_rows)

openxlsx::saveWorkbook(
  wb,
  file.path(OUT_DIR, "S20_Thrane2018_external_spatial_recurrence.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# 21. Console summary
# -----------------------------

message2("S20 Thrane 2018 CANONICAL v1.0 completed.")
message2("Output directory: ", OUT_DIR)

message2("Selected gene-set file:")
print(best_file)

message2("Sample input log:")
print(sample_input_log)

message2("Gene coverage summary:")
print(coverage_summary)

message2("Key readout:")
print(key_readout)

message2("Inclusion decision helper:")
print(include_decision)
