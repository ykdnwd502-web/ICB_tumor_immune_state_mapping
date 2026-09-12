############################################################
## 01_GSE250636_build_18J.R
##
## INPUT
##   data_external_spatial_18J/GSE250636/unpacked_recursive/GSE250636_RAW/
##
## OUTPUT
##   results/tables/revision_external_spatial_recurrence_18J/GSE250636/
##
## Role in public sequence:
##   raw GSE250636 -> historical 18J outputs
##
## Scientific/statistical definitions are retained from the archived
## "18J-light GSE250636 external spatial recurrence audit v3 FIXED".
## Public execution edits are limited to path configuration and explicit
## dplyr namespaces for package-compatibility.
############################################################

# 18J-light GSE250636 external spatial recurrence audit
# v3 FIXED
#
# Purpose:
#   External melanoma spatial recurrence audit using GSE250636
#   flat GEO files.
#
# Key fixes vs v2:
#   1) Robust state gene-set discovery and parsing.
#   2) Supports one-cell multi-gene lists separated by ; , | or spaces.
#   3) Fail-fast if primary-pair state signatures have poor coverage.
#   4) Memory-safe module scoring using only required scoring genes.
#   5) Empty result tables retain columns and do not crash silently.
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
GEO_ID <- "GSE250636"

FLAT_DIR <- file.path(
  PROJECT_DIR,
  "data_external_spatial_18J",
  GEO_ID,
  "unpacked_recursive",
  "GSE250636_RAW"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "revision_external_spatial_recurrence_18J",
  GEO_ID
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## Public sequential rerun: remove only this step's current core outputs
## so a failed rerun cannot be masked by stale files.
PUBLIC_CORE_OUTPUTS <- c(
  "18J_GSE250636_state_genes_used_long.csv",
  "18J_GSE250636_gene_coverage.csv",
  "18J_GSE250636_key_readout.csv",
  "18J_GSE250636_primary_pair_dual_high_OR.csv",
  "18J_GSE250636_primary_pair_Moran_descriptive.csv",
  "18J_GSE250636_primary_pair_raw_vs_residualized_correlation.csv",
  "18J_GSE250636_spot_scores_with_residualized_primary_pair.csv",
  "18J_GSE250636_spot_state_composition_scores.csv"
)
old_core <- file.path(OUT_DIR, PUBLIC_CORE_OUTPUTS)
old_core <- old_core[file.exists(old_core)]
if (length(old_core)) {
  unlink(old_core, force = TRUE)
}

PRIMARY_X_STATE <- "Myeloid_Treg_Immunosuppressive"
PRIMARY_Y_STATE <- "Tumor_dedifferentiation_Stromal_remodeling"

PRIMARY_X_LABEL <- "Myeloid–Treg immunosuppressive"
PRIMARY_Y_LABEL <- "Tumor-dedifferentiation/stromal-remodeling"

PRIMARY_K <- 6
RUN_MORAN_DESCRIPTIVE <- TRUE

# Strong recommendation:
# If auto-detection fails, manually set this to your final locked state-gene file.
# The file must contain the 4 final states and gene symbols.
# Examples of acceptable formats:
#   long format: state / gene
#   long format: FinalState / GeneSymbol
#   long format: state / retained_genes, where retained_genes = "GENE1;GENE2;..."
#   wide format: one column per state
STATE_GENESET_FILE <- file.path(
  PROJECT_DIR,
  "public_release",
  "S20_external_spatial_recurrence",
  "config",
  "S20_locked_state_genes.csv"
)

if (!file.exists(STATE_GENESET_FILE)) {
  stop(
    "Required public S20 state-gene lock is missing: ",
    STATE_GENESET_FILE,
    call. = FALSE
  )
}

GENESET_SEARCH_ROOTS <- c(
  file.path(PROJECT_DIR, "results"),
  file.path(PROJECT_DIR, "data_processed"),
  file.path(PROJECT_DIR, "data"),
  file.path(PROJECT_DIR, "tables")
)

MIN_STATE_GENES_REQUIRED <- 20
MIN_CRITICAL_DETECTED_GENES <- 10

# -----------------------------
# 1. Packages
# -----------------------------

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "stringr",
  "Matrix", "openxlsx", "Seurat"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required package: ", pkg, call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(Matrix)
  library(openxlsx)
  library(Seurat)
})

# -----------------------------
# 2. Helper functions
# -----------------------------

message2 <- function(...) {
  cat(paste0("[18J-GSE250636-v3] ", paste0(..., collapse = ""), "\n"))
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
    grepl("^[A-Z0-9][A-Z0-9._-]{1,30}$", x)
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
# 3. Detect flat GEO sample files
# -----------------------------

if (!dir.exists(FLAT_DIR)) {
  stop("FLAT_DIR does not exist: ", FLAT_DIR, call. = FALSE)
}

all_files <- list.files(FLAT_DIR, full.names = TRUE)

filtered_h5 <- all_files[
  grepl("_filtered_feature_bc_matrix\\.h5$", basename(all_files), ignore.case = TRUE)
]

if (length(filtered_h5) == 0) {
  stop("No filtered_feature_bc_matrix.h5 files found in: ", FLAT_DIR, call. = FALSE)
}

extract_sample_prefix <- function(path) {
  bn <- basename(path)
  sub("_filtered_feature_bc_matrix\\.h5$", "", bn, ignore.case = TRUE)
}

sample_prefixes <- extract_sample_prefix(filtered_h5)

sample_files <- data.frame(
  dataset = GEO_ID,
  sample_prefix = sample_prefixes,
  h5_file = filtered_h5,
  position_file = file.path(FLAT_DIR, paste0(sample_prefixes, "_tissue_positions_list.csv.gz")),
  scalefactor_file = file.path(FLAT_DIR, paste0(sample_prefixes, "_scalefactors_json.json.gz")),
  lowres_image = file.path(FLAT_DIR, paste0(sample_prefixes, "_tissue_lowres_image.png.gz")),
  hires_image = file.path(FLAT_DIR, paste0(sample_prefixes, "_tissue_hires_image.png.gz")),
  stringsAsFactors = FALSE
)

sample_files$position_file_exists <- file.exists(sample_files$position_file)
sample_files$h5_file_exists <- file.exists(sample_files$h5_file)

sample_files <- sample_files %>%
  mutate(
    sample_id = paste0(
      GEO_ID,
      "_",
      str_extract(sample_prefix, "sample[0-9]+")
    )
  ) %>%
  filter(h5_file_exists, position_file_exists) %>%
  arrange(sample_id)

write_csv(sample_files, "18J_GSE250636_flat_sample_file_map.csv")

message2("Detected complete filtered h5 + position-file samples: ", nrow(sample_files))
print(sample_files[, c("sample_id", "sample_prefix", "h5_file_exists", "position_file_exists")])

if (nrow(sample_files) == 0) {
  stop("No complete h5 + position-file sample pairs found.", call. = FALSE)
}

# -----------------------------
# 4. Robust state gene-set discovery
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

  # Candidate long-format columns
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

  # Long format: try all state/gene pairs
  if (length(state_candidates) > 0 && length(gene_candidates) > 0) {
    for (sc in state_candidates) {
      for (gc in gene_candidates) {
        tmp <- data.frame()

        for (i in seq_len(nrow(df))) {
          st_raw <- as.character(df[[sc]][i])
          st <- standardize_state_name(st_raw)
          genes <- parse_gene_cell(df[[gc]][i])

          if (!is.na(st) && length(genes) > 0) {
            tmp <- bind_rows(
              tmp,
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

        if (nrow(tmp) > 0) {
          out <- bind_rows(out, tmp)
        }
      }
    }
  }

  # Wide format: one column per state
  possible_state_cols <- cn[
    !is.na(standardize_state_name(cn))
  ]

  if (length(possible_state_cols) > 0) {
    for (cc in possible_state_cols) {
      st <- standardize_state_name(cc)
      genes <- parse_gene_cell(df[[cc]])

      if (!is.na(st) && length(genes) > 0) {
        out <- bind_rows(
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

  out %>%
    distinct(state, gene, .keep_all = TRUE)
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
      out <- bind_rows(out, tmp)
    }
  }

  out %>%
    filter(state %in% required_states) %>%
    distinct(state, gene, .keep_all = TRUE)
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

  if (length(files) == 0) return(character(0))

  files[
    grepl(
      "state|signature|gene|geneset|gene_set|gene.pool|gene_pool|retained|final|locked|18S",
      basename(files),
      ignore.case = TRUE
    )
  ]
}

candidate_files <- find_candidate_geneset_files()

if (length(candidate_files) == 0) {
  stop(
    "No candidate gene-set files found. Manually set STATE_GENESET_FILE near the top of this script.",
    call. = FALSE
  )
}

message2("Auditing candidate gene-set files: ", length(candidate_files))

audit_rows <- data.frame()
parsed_list <- list()

for (f in candidate_files) {
  parsed <- tryCatch(
    read_geneset_file_robust(f),
    error = function(e) data.frame()
  )

  if (nrow(parsed) > 0) {
    cnt <- parsed %>%
      dplyr::count(state, name = "n_genes")

    for (st in required_states) {
      nst <- cnt$n_genes[cnt$state == st]
      if (length(nst) == 0) nst <- 0

      audit_rows <- bind_rows(
        audit_rows,
        data.frame(
          file = f,
          basename = basename(f),
          state = st,
          n_genes = nst,
          stringsAsFactors = FALSE
        )
      )
    }

    parsed_list[[f]] <- parsed
  } else {
    for (st in required_states) {
      audit_rows <- bind_rows(
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

write_csv(audit_rows, "18J_GSE250636_state_geneset_candidate_audit.csv")

audit_wide <- audit_rows %>%
  tidyr::pivot_wider(names_from = state, values_from = n_genes, values_fill = 0) %>%
  mutate(
    n_states_ok = rowSums(across(all_of(required_states), ~ .x >= MIN_STATE_GENES_REQUIRED)),
    total_genes = rowSums(across(all_of(required_states)))
  ) %>%
  arrange(desc(n_states_ok), desc(total_genes))

write_csv(audit_wide, "18J_GSE250636_state_geneset_candidate_audit_wide.csv")

message2("Top gene-set candidates:")
print(head(audit_wide, 10))

best_file <- audit_wide$file[
  audit_wide$n_states_ok == length(required_states)
][1]

if (is.na(best_file) || !best_file %in% names(parsed_list)) {
  stop(
    "No valid state gene-set file found. ",
    "Check 18J_GSE250636_state_geneset_candidate_audit_wide.csv. ",
    "Manually set STATE_GENESET_FILE to the final locked state gene-set file.",
    call. = FALSE
  )
}

gs_long <- parsed_list[[best_file]] %>%
  dplyr::select(state, gene) %>%
  distinct()

message2("Selected state gene-set file: ", best_file)

state_counts <- gs_long %>%
  dplyr::count(state, name = "n_genes") %>%
  arrange(state)

message2("Parsed state gene counts:")
print(state_counts)

write_csv(gs_long, "18J_GSE250636_state_genes_used_long.csv")
write_csv(state_counts, "18J_GSE250636_state_genes_used_counts.csv")

missing_states <- setdiff(required_states, unique(gs_long$state))

if (length(missing_states) > 0) {
  stop(
    "Missing states after gene-set parsing: ",
    paste(missing_states, collapse = ", "),
    call. = FALSE
  )
}

too_small <- state_counts %>%
  filter(n_genes < MIN_STATE_GENES_REQUIRED)

if (nrow(too_small) > 0) {
  stop(
    "Some parsed state gene sets are too small. Check selected file: ",
    best_file,
    call. = FALSE
  )
}

state_gene_sets <- split(gs_long$gene, gs_long$state)

# -----------------------------
# 5. Composition marker modules
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

# -----------------------------
# 6. Read positions
# -----------------------------

read_positions_10x_old <- function(path) {
  pos <- data.table::fread(path, header = FALSE, check.names = FALSE)

  if (ncol(pos) < 6) {
    stop("Unexpected tissue_positions_list format: ", path, call. = FALSE)
  }

  colnames(pos)[1:6] <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )

  pos %>%
    transmute(
      barcode = as.character(barcode),
      in_tissue = as.integer(in_tissue),
      array_row = as.numeric(array_row),
      array_col = as.numeric(array_col),
      spatial_y = as.numeric(pxl_row_in_fullres),
      spatial_x = as.numeric(pxl_col_in_fullres)
    )
}

# -----------------------------
# 7. Memory-safe module scoring
# -----------------------------

score_modules <- function(counts, gene_sets) {
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }

  row_clean <- clean_gene(rownames(counts))
  all_needed <- unique(unlist(gene_sets))
  rows_needed <- which(row_clean %in% all_needed)

  if (length(rows_needed) == 0) {
    stop("No scoring genes detected in count matrix.", call. = FALSE)
  }

  lib <- Matrix::colSums(counts)
  lib[lib <= 0] <- NA_real_

  sub <- counts[rows_needed, , drop = FALSE]
  sub_gene <- row_clean[rows_needed]

  # Dense conversion only for scoring genes, not the full matrix.
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

    coverage <- bind_rows(
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
# 8. Read and score all samples
# -----------------------------

all_scores <- data.frame()
all_coverage <- data.frame()
sample_input_log <- data.frame()

for (i in seq_len(nrow(sample_files))) {
  sf <- sample_files[i, ]

  message2("Reading sample ", i, "/", nrow(sample_files), ": ", sf$sample_id)

  expr <- Seurat::Read10X_h5(sf$h5_file)

  if (is.list(expr)) {
    if ("Gene Expression" %in% names(expr)) {
      expr <- expr[["Gene Expression"]]
    } else {
      expr <- expr[[1]]
    }
  }

  if (!inherits(expr, "dgCMatrix")) {
    expr <- as(expr, "dgCMatrix")
  }

  pos <- read_positions_10x_old(sf$position_file)

  common <- intersect(colnames(expr), pos$barcode)

  if (length(common) == 0) {
    expr_bc_clean <- sub("-1$", "", colnames(expr))
    pos_bc_clean <- sub("-1$", "", pos$barcode)

    colnames(expr) <- expr_bc_clean
    pos$barcode <- pos_bc_clean

    common <- intersect(colnames(expr), pos$barcode)
  }

  message2("  barcode overlap: ", length(common))

  if (length(common) < 20) {
    warning("Very low barcode-position overlap for ", sf$sample_id)
    next
  }

  expr <- expr[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]

  keep <- which(pos$in_tissue == 1)

  expr <- expr[, keep, drop = FALSE]
  pos <- pos[keep, , drop = FALSE]

  message2("  in-tissue spots: ", ncol(expr))

  if (ncol(expr) < 20) {
    warning("Too few in-tissue spots for ", sf$sample_id)
    next
  }

  sc <- score_modules(expr, all_gene_sets)

  scores <- sc$scores %>%
    mutate(
      dataset = GEO_ID,
      sample_id = sf$sample_id,
      sample_prefix = sf$sample_prefix
    ) %>%
    left_join(
      pos %>%
        transmute(
          spot_barcode = barcode,
          spatial_x = spatial_x,
          spatial_y = spatial_y,
          array_row = array_row,
          array_col = array_col,
          in_tissue = in_tissue
        ),
      by = "spot_barcode"
    ) %>%
    mutate(
      spot_id = paste(sample_id, spot_barcode, sep = "__")
    )

  cov <- sc$coverage %>%
    mutate(
      dataset = GEO_ID,
      sample_id = sf$sample_id,
      sample_prefix = sf$sample_prefix
    )

  all_scores <- bind_rows(all_scores, scores)
  all_coverage <- bind_rows(all_coverage, cov)

  sample_input_log <- bind_rows(
    sample_input_log,
    data.frame(
      dataset = GEO_ID,
      sample_id = sf$sample_id,
      sample_prefix = sf$sample_prefix,
      h5_file = sf$h5_file,
      position_file = sf$position_file,
      n_genes = nrow(expr),
      n_in_tissue_spots = ncol(expr),
      n_modules_scored = length(all_gene_sets),
      stringsAsFactors = FALSE
    )
  )

  rm(expr, pos, sc, scores, cov)
  gc()
}

if (nrow(all_scores) == 0) {
  stop("No sample yielded usable scores.", call. = FALSE)
}

write_csv(sample_input_log, "18J_GSE250636_sample_input_log.csv")
write_csv(all_coverage, "18J_GSE250636_gene_coverage.csv")
write_csv(all_scores, "18J_GSE250636_spot_state_composition_scores.csv")

coverage_summary <- all_coverage %>%
  group_by(module) %>%
  summarise(
    min_detected = min(n_genes_detected, na.rm = TRUE),
    median_detected = median(n_genes_detected, na.rm = TRUE),
    max_detected = max(n_genes_detected, na.rm = TRUE),
    min_coverage = min(coverage_fraction, na.rm = TRUE),
    median_coverage = median(coverage_fraction, na.rm = TRUE),
    max_coverage = max(coverage_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(module)

write_csv(coverage_summary, "18J_GSE250636_gene_coverage_summary.csv")

message2("Gene coverage summary:")
print(coverage_summary)

critical_bad <- coverage_summary %>%
  filter(
    module %in% c(PRIMARY_X_STATE, PRIMARY_Y_STATE),
    median_detected < MIN_CRITICAL_DETECTED_GENES
  )

if (nrow(critical_bad) > 0) {
  stop(
    "Critical primary-pair state signatures have low detected-gene coverage. ",
    "Check 18J_GSE250636_gene_coverage_summary.csv and selected state gene-set file.",
    call. = FALSE
  )
}

# -----------------------------
# 9. Derived composition variables
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
# 10. Residualization
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

write_csv(df, "18J_GSE250636_spot_scores_with_residualized_primary_pair.csv")

# -----------------------------
# 11. State-composition correlations
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

        cor_rows <- bind_rows(
          cor_rows,
          data.frame(
            dataset = GEO_ID,
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

write_csv(cor_rows, "18J_GSE250636_state_composition_correlation.csv")

# -----------------------------
# 12. Primary pair correlation
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

      primary_cor <- bind_rows(
        primary_cor,
        data.frame(
          dataset = GEO_ID,
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

write_csv(primary_cor, "18J_GSE250636_primary_pair_raw_vs_residualized_correlation.csv")

# -----------------------------
# 13. Dual-high OR
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
      dual_rows <- bind_rows(
        dual_rows,
        data.frame(
          dataset = GEO_ID,
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
    dual_rows <- bind_rows(
      dual_rows,
      data.frame(
        dataset = GEO_ID,
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

write_csv(dual_rows, "18J_GSE250636_primary_pair_dual_high_OR.csv")

# -----------------------------
# 14. Descriptive bivariate Moran-type statistic
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

      moran_rows <- bind_rows(
        moran_rows,
        data.frame(
          dataset = GEO_ID,
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
      group_by(adjustment) %>%
      summarise(
        dataset = GEO_ID,
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

    moran_rows <- bind_rows(moran_rows, moran_pooled)
  }
}

write_csv(moran_rows, "18J_GSE250636_primary_pair_Moran_descriptive.csv")

# -----------------------------
# 15. Key readout
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
  key_readout <- bind_rows(
    key_readout,
    primary_cor %>%
      filter(sample_id == "pooled") %>%
      mutate(
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
  key_readout <- bind_rows(
    key_readout,
    dual_rows %>%
      filter(sample_id == "pooled") %>%
      mutate(
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
  key_readout <- bind_rows(
    key_readout,
    moran_rows %>%
      filter(sample_id == "pooled_weighted") %>%
      mutate(
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

write_csv(key_readout, "18J_GSE250636_key_readout.csv")

if (nrow(key_readout) == 0) {
  stop(
    "Key readout is empty. Check gene coverage and primary state score columns.",
    call. = FALSE
  )
}

# -----------------------------
# 16. Excel workbook
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
    GEO_ID,
    "External melanoma spatial recurrence audit of composition-linked state-score co-variation.",
    paste(PRIMARY_X_LABEL, "vs", PRIMARY_Y_LABEL),
    "No matched ICB-response annotation was used or assumed.",
    "This analysis is not response validation, not biomarker validation, and not functional niche validation.",
    "Flat GEO supplementary files: filtered_feature_bc_matrix.h5 plus tissue_positions_list.csv.gz for each sample.",
    best_file,
    "State-composition correlations, raw/residualized primary-pair correlations, dual-high odds ratios, and descriptive Moran-type statistics."
  ),
  stringsAsFactors = FALSE
)

add_sheet(wb, "README", readme)
add_sheet(wb, "Key_readout", key_readout)
add_sheet(wb, "Sample_file_map", sample_files)
add_sheet(wb, "Sample_input_log", sample_input_log)
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
  file.path(OUT_DIR, "historical_GSE250636_external_spatial_recurrence_audit.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# 17. Console summary
# -----------------------------

message2("18J-light GSE250636 v3 completed.")
message2("Output directory: ", OUT_DIR)

message2("Sample input log:")
print(sample_input_log)

message2("Selected gene-set file:")
print(best_file)

message2("Parsed state gene counts:")
print(state_counts)

message2("Gene coverage summary:")
print(coverage_summary)

message2("Key readout:")
print(key_readout)

# ------------------------------------------------------------
# Public sequential end gate
# ------------------------------------------------------------
expected_step01 <- PUBLIC_CORE_OUTPUTS

missing_step01 <- expected_step01[
  !file.exists(file.path(OUT_DIR, expected_step01))
]

if (length(missing_step01)) {
  stop(
    "STEP 01 finished but expected 18J file(s) are missing:\n",
    paste(missing_step01, collapse = "\n"),
    call. = FALSE
  )
}

cat(
  "\n============================================================\n",
  "STEP 01 PASS — GSE250636 raw -> 18J\n",
  "Expected core 18J outputs: 8/8 present\n",
  "Output: ", OUT_DIR, "\n",
  "============================================================\n",
  sep = ""
)
