############################################################
## 02_GSE78220_state_scoring_response_survival_analysis.R 
##
## Purpose:
##   Reproduce predefined tumor-immune state scoring in GSE78220 and perform
##   descriptive ICB-response association, baseline-only overall-survival
##   analysis, state-score correlation analysis, and manuscript Figure 7.
##
## Run after 01:
##   01_bulk_discovery_analysis.R
##   and after GSE78220 clinical metadata have been placed under
##   data_raw/external_validation/GSE78220/.
##
## Frozen analysis definitions:
##   - Response/AUC analysis uses Responder versus NonResponder bulk samples.
##   - Missing response labels are not imputed.
##   - Cox analysis is strictly restricted to pretreatment/baseline samples.
##   - Pt27A/Pt27B are collapsed to one patient for patient-level Cox analysis.
##   - Expected baseline-only Cox set: 25 patients and 11 death events.
##   - All outputs are written under results/tables/external_validation and
##     results/figures/external_validation.
############################################################


############################################################
## PART A. External GSE78220 state scoring using CLEAN GMT
############################################################

############################################################
## Part A: GSE78220 state scoring and response metadata matching
## Purpose:
##   Score external melanoma ICB cohort GSE78220 using the CLEAN
##   mechanism-defined tumor-immune state GMT generated in Step 01.
##
## Run after:
##   01_bulk_discovery_analysis.R
##
## Canonical state-score output used by downstream analyses:
##   results/tables/external_validation/GSE78220_external_state_scores.csv
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

############################################################
## 0. Project configuration
############################################################
project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") project_dir <- getwd()

message("Project directory: ", project_dir)

raw_dir       <- file.path(project_dir, "data_raw")
processed_dir <- file.path(project_dir, "data_processed")
table_dir     <- file.path(project_dir, "results/tables")
figure_dir    <- file.path(project_dir, "results/figures")
log_dir       <- file.path(project_dir, "logs")

out_table_dir <- file.path(table_dir, "external_validation")
out_fig_dir   <- file.path(figure_dir, "external_validation")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Packages
############################################################
required_pkgs <- c("readxl", "data.table", "dplyr", "tidyr", "stringr", "ggplot2", "ggrepel", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
    ". Restore the repository environment (for example with renv::restore()) before running this script."
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
  message("Saved table: ", file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
  message("Saved RDS: ", file)
}

safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = file, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved figure: ", file)
}

save_session_info <- function(script_name = "02_GSE78220_state_scoring_response_survival") {
  log_file <- file.path(log_dir, paste0(script_name, "_sessionInfo.txt"))
  sink(log_file)
  print(sessionInfo())
  sink()
  message("sessionInfo saved to: ", log_file)
}


############################################################
## 1A. Unified plotting theme and color palette
## Style anchor: Figure 1 / bulk discovery panel style
############################################################

theme_icb <- function(base_size = 12.5, base_family = "sans") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 2.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size - 1),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1.5),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1.5),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.6)
    )
}

state_colors <- c(
  "Immune_defective_Cold" = "#4DBBD5",
  "Myeloid_Treg_Immunosuppressive" = "#00A087",
  "Tumor_dedifferentiation_Stromal_remodeling" = "#E64B35",
  "Melanocytic_Differentiation" = "#3C5488"
)

heatmap_low  <- "#3B82F6"
heatmap_mid  <- "white"
heatmap_high <- "#EF4444"

response_group_colors <- c(
  "Responder" = "#E64B35",
  "StableDisease" = "#00A087",
  "NonResponder" = "#4DBBD5"
)

normalize_name <- function(x) {
  x <- as.character(x)
  x <- gsub("\\ufeff", "", x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- gsub("\\ufeff", "", x)
  x <- trimws(x)
  x <- gsub("^['\"]+|['\"]+$", "", x)
  toupper(x)
}

zvec <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

find_first_existing <- function(paths) {
  paths <- unique(paths)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

############################################################
## 2. Canonical state names and display labels
############################################################
state_cols_clean <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_label <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

############################################################
## 3. Locate GSE78220 expression file and CLEAN GMT
############################################################
gse78220_expr_candidates <- c(
  file.path(raw_dir, "external_validation/GSE78220/GSE78220_PatientFPKM.xlsx"),
  list.files(file.path(raw_dir, "external_validation/GSE78220"), pattern = "GSE78220.*(FPKM|fpkm|expr|expression).*\\.(xlsx|xls|csv|tsv|txt|gz)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
  list.files(project_dir, pattern = "GSE78220.*(PatientFPKM|FPKM|fpkm|expr|expression).*\\.(xlsx|xls|csv|tsv|txt|gz)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
)

gse78220_expr_file <- find_first_existing(gse78220_expr_candidates)
if (is.na(gse78220_expr_file)) {
  stop("No GSE78220 expression file found. Expected data_raw/external_validation/GSE78220/GSE78220_PatientFPKM.xlsx or equivalent.")
}
message("Using GSE78220 expression file: ", gse78220_expr_file)

gmt_candidates <- c(
  file.path(table_dir, "ICBcomb_final_input_gene_sets_CLEAN/02_final_ICBcomb_gene_sets/final_ICBcomb_gene_sets_CLEAN.gmt"),
  list.files(file.path(table_dir, "ICBcomb_final_input_gene_sets_CLEAN"), pattern = "final_ICBcomb_gene_sets_CLEAN\\.gmt$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
)

gmt_file <- find_first_existing(gmt_candidates)
if (is.na(gmt_file)) {
  stop("No CLEAN GMT found. Please generate the frozen GMT and run 01_GSE244982_bulk_state_scoring.R first.")
}
message("Using CLEAN GMT: ", gmt_file)

############################################################
## 4. Read GMT
############################################################
read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]
  gs <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 3) next
    nm <- parts[1]
    genes <- clean_gene_symbol(parts[-c(1, 2)])
    genes <- unique(genes[!is.na(genes) & genes != ""])
    gs[[nm]] <- genes
  }
  gs
}

gene_sets <- read_gmt(gmt_file)
message("Imported gene sets: ", length(gene_sets))

safe_write_csv(
  data.frame(GeneSet = names(gene_sets), n_genes = sapply(gene_sets, length), row.names = NULL),
  file.path(out_table_dir, "GSE78220_imported_GMT_gene_set_summary.csv")
)

############################################################
## 5. Read GSE78220 expression matrix
############################################################
read_table_any <- function(file) {
  if (grepl("\\.xlsx$|\\.xls$", file, ignore.case = TRUE)) {
    sheets <- readxl::excel_sheets(file)
    sheet_dfs <- lapply(sheets, function(sh) {
      df <- as.data.frame(readxl::read_excel(file, sheet = sh), check.names = FALSE)
      df$.sheet_name <- sh
      df
    })
    names(sheet_dfs) <- sheets
    return(sheet_dfs)
  }

  if (grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE)) {
    return(list(data = data.table::fread(file, data.table = FALSE, check.names = FALSE)))
  }

  if (grepl("\\.csv$|\\.csv\\.gz$", file, ignore.case = TRUE)) {
    return(list(data = data.table::fread(file, data.table = FALSE, check.names = FALSE)))
  }

  stop("Unsupported expression file type: ", file)
}

is_numeric_like <- function(x) {
  xx <- suppressWarnings(as.numeric(as.character(x)))
  mean(!is.na(xx)) >= 0.8
}

find_gene_col <- function(df) {
  cn <- colnames(df)
  cn_norm <- normalize_name(cn)
  preferred <- c("gene", "genes", "gene_symbol", "genesymbol", "hugo_symbol", "symbol", "external_gene_name")
  hit <- which(cn_norm %in% preferred)
  if (length(hit) > 0) return(cn[hit[1]])
  cn[1]
}

score_expression_sheet <- function(df) {
  if (!is.data.frame(df) || nrow(df) < 10 || ncol(df) < 5) return(NULL)
  df$.sheet_name <- NULL
  gene_col <- find_gene_col(df)
  numeric_cols <- setdiff(colnames(df), gene_col)
  numeric_cols <- numeric_cols[sapply(df[numeric_cols], is_numeric_like)]
  if (length(numeric_cols) < 5) return(NULL)
  n_gene_like <- sum(!is.na(df[[gene_col]]) & df[[gene_col]] != "")
  data.frame(
    gene_col = gene_col,
    n_rows = nrow(df),
    n_numeric_cols = length(numeric_cols),
    n_gene_like = n_gene_like,
    score = nrow(df) + 20 * length(numeric_cols) + n_gene_like,
    stringsAsFactors = FALSE
  )
}

all_tables <- read_table_any(gse78220_expr_file)
expr_candidates <- lapply(names(all_tables), function(nm) {
  s <- score_expression_sheet(all_tables[[nm]])
  if (is.null(s)) return(NULL)
  s$sheet <- nm
  s
})
expr_candidates <- dplyr::bind_rows(expr_candidates)

if (nrow(expr_candidates) == 0) {
  stop("Could not identify an expression-like table in GSE78220 file. Please inspect workbook sheets manually.")
}

expr_candidates <- expr_candidates %>% arrange(desc(score))
safe_write_csv(expr_candidates, file.path(out_table_dir, "GSE78220_expression_sheet_candidates_audit.csv"))

expr_sheet <- expr_candidates$sheet[1]
gene_col <- expr_candidates$gene_col[1]
message("Selected expression sheet/table: ", expr_sheet, " ; gene column: ", gene_col)

expr_df <- all_tables[[expr_sheet]]
expr_df$.sheet_name <- NULL

numeric_cols <- setdiff(colnames(expr_df), gene_col)
numeric_cols <- numeric_cols[sapply(expr_df[numeric_cols], is_numeric_like)]

expr_mat <- as.matrix(as.data.frame(lapply(expr_df[numeric_cols], function(x) as.numeric(as.character(x))), check.names = FALSE))
rownames(expr_mat) <- clean_gene_symbol(expr_df[[gene_col]])
colnames(expr_mat) <- numeric_cols

## Remove empty/duplicated genes and collapse duplicates by mean.
keep <- !is.na(rownames(expr_mat)) & rownames(expr_mat) != ""
expr_mat <- expr_mat[keep, , drop = FALSE]
expr_mat <- expr_mat[rowSums(!is.na(expr_mat)) > 0, , drop = FALSE]

if (any(duplicated(rownames(expr_mat)))) {
  message("Collapsing duplicated gene symbols by mean.")
  expr_mat <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE) / as.numeric(table(rownames(expr_mat)))[rownames(rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE))]
}

## Decide whether log2 transform is needed.
expr_max <- suppressWarnings(max(expr_mat, na.rm = TRUE))
expr_min <- suppressWarnings(min(expr_mat, na.rm = TRUE))
log2_transform_applied <- FALSE

if (is.finite(expr_max) && expr_max > 50 && expr_min >= 0) {
  expr_for_score <- log2(expr_mat + 1)
  log2_transform_applied <- TRUE
} else {
  expr_for_score <- expr_mat
}

input_audit <- data.frame(
  Dataset = "GSE78220",
  ExpressionFile = gse78220_expr_file,
  ExpressionSheet = expr_sheet,
  GeneColumn = gene_col,
  n_genes = nrow(expr_mat),
  n_samples = ncol(expr_mat),
  min_expr = expr_min,
  max_expr = expr_max,
  log2_transform_applied = log2_transform_applied,
  stringsAsFactors = FALSE
)

safe_write_csv(input_audit, file.path(out_table_dir, "GSE78220_external_expression_input_audit.csv"))
safe_save_rds(expr_mat, file.path(processed_dir, "GSE78220_external_expression_matrix_raw.rds"))
safe_save_rds(expr_for_score, file.path(processed_dir, "GSE78220_external_expression_matrix_for_scoring.rds"))

############################################################
## 6. Score GMT gene sets by same z-mean method as Step 01
############################################################
score_gene_set_zmean <- function(expr_mat, gene_sets) {
  expr_genes <- rownames(expr_mat)
  expr_z <- t(scale(t(expr_mat)))
  expr_z[is.na(expr_z)] <- 0

  score_list <- list()
  presence_list <- list()

  for (gs in names(gene_sets)) {
    genes <- clean_gene_symbol(gene_sets[[gs]])
    present <- intersect(genes, expr_genes)
    message(gs, ": ", length(present), "/", length(genes), " genes present")

    presence_list[[gs]] <- data.frame(
      GeneSet = gs,
      Gene = genes,
      Present = genes %in% expr_genes,
      stringsAsFactors = FALSE
    )

    if (length(present) == 0) {
      score <- rep(NA_real_, ncol(expr_mat))
      names(score) <- colnames(expr_mat)
    } else {
      score <- colMeans(expr_z[present, , drop = FALSE], na.rm = TRUE)
    }
    score_list[[gs]] <- score
  }

  score_mat <- do.call(rbind, score_list)
  rownames(score_mat) <- names(score_list)

  list(
    score_matrix = score_mat,
    raw_scores = as.data.frame(t(score_mat), check.names = FALSE),
    gene_presence = dplyr::bind_rows(presence_list)
  )
}

score_res <- score_gene_set_zmean(expr_for_score, gene_sets)
raw_scores <- score_res$raw_scores
raw_scores$Sample <- rownames(raw_scores)
raw_scores <- raw_scores[, c("Sample", setdiff(colnames(raw_scores), "Sample")), drop = FALSE]

safe_write_csv(raw_scores, file.path(out_table_dir, "GSE78220_external_raw_gene_set_scores.csv"))
safe_write_csv(score_res$gene_presence, file.path(out_table_dir, "GSE78220_external_gene_presence_all_gene_sets.csv"))

############################################################
## 7. Build external four-state scores using CLEAN GMT mapping
############################################################
get_signature_cols <- function(raw_scores, target_names, state_name, min_required = 1) {
  available <- setdiff(colnames(raw_scores), "Sample")
  available_norm <- normalize_name(available)
  names(available_norm) <- available
  target_norm <- normalize_name(target_names)

  matched <- available[available_norm %in% target_norm]
  if (length(matched) == 0) {
    for (tg in target_norm) {
      idx <- grep(tg, available_norm, fixed = TRUE)
      if (length(idx) > 0) matched <- c(matched, available[idx])
    }
  }
  matched <- unique(matched)
  message(state_name, " matched signature columns: ", paste(matched, collapse = ", "))
  if (length(matched) < min_required) {
    stop("Insufficient matched signature columns for ", state_name, ". Available columns: ", paste(available, collapse = ", "))
  }
  matched
}

immune_cols <- get_signature_cols(
  raw_scores,
  c("Immune_defective_Cold_RESTORE_hybrid", "Immune_defective_Cold_RESTORE_data_driven", "Immune_defective_Cold_RESTORE_curated"),
  "Immune_defective_Cold"
)
myeloid_cols <- get_signature_cols(
  raw_scores,
  c("Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid", "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven", "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated"),
  "Myeloid_Treg_Immunosuppressive"
)
tumor_cols <- get_signature_cols(
  raw_scores,
  c("Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid", "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven", "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated"),
  "Tumor_dedifferentiation_Stromal_remodeling"
)
melanocytic_cols <- get_signature_cols(
  raw_scores,
  c("Melanocytic_Differentiation_REFERENCE_hybrid", "Melanocytic_Differentiation_REFERENCE_data_driven", "Melanocytic_Differentiation_REFERENCE_curated"),
  "Melanocytic_Differentiation"
)

mapping_tbl <- data.frame(
  FinalState = c(rep("Immune_defective_Cold", length(immune_cols)),
                 rep("Myeloid_Treg_Immunosuppressive", length(myeloid_cols)),
                 rep("Tumor_dedifferentiation_Stromal_remodeling", length(tumor_cols)),
                 rep("Melanocytic_Differentiation", length(melanocytic_cols))),
  RawSignatureColumn = c(immune_cols, myeloid_cols, tumor_cols, melanocytic_cols),
  DirectionForFinalState = c(rep("inverse", length(immune_cols)),
                             rep("positive", length(myeloid_cols)),
                             rep("positive", length(tumor_cols)),
                             rep("positive", length(melanocytic_cols))),
  stringsAsFactors = FALSE
)
safe_write_csv(mapping_tbl, file.path(out_table_dir, "GSE78220_external_final_state_signature_mapping.csv"))

state_scores <- data.frame(
  Sample = raw_scores$Sample,
  Immune_defective_Cold = zvec(-rowMeans(raw_scores[, immune_cols, drop = FALSE], na.rm = TRUE)),
  Myeloid_Treg_Immunosuppressive = zvec(rowMeans(raw_scores[, myeloid_cols, drop = FALSE], na.rm = TRUE)),
  Tumor_dedifferentiation_Stromal_remodeling = zvec(rowMeans(raw_scores[, tumor_cols, drop = FALSE], na.rm = TRUE)),
  Melanocytic_Differentiation = zvec(rowMeans(raw_scores[, melanocytic_cols, drop = FALSE], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

state_scores$Dominant_state <- apply(
  state_scores[, state_cols_clean, drop = FALSE],
  1,
  function(x) state_cols_clean[which.max(x)]
)

############################################################
## 8. Try to attach response metadata
############################################################
make_key <- function(x) {
  ## Keep patient suffixes such as Pt27A/Pt27B when present.
  ## The previous suffix-stripping rule can collapse distinct samples.
  x <- toupper(as.character(x))
  x <- gsub("\\.baseline$|\\.pre$|\\.pretreatment$|\\.on$|\\.on_treatment$", "", x, ignore.case = TRUE)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

make_key_loose <- function(x) {
  ## Backup key only: collapse a trailing single letter after Pt number.
  y <- make_key(x)
  y <- gsub("^(PT[0-9]+)[A-Z]$", "\\1", y)
  y
}

map_response_group <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[_-]+", " ", y)
  out <- rep(NA_character_, length(y))

  out[grepl("complete response|partial response|\\bcr\\b|\\bpr\\b|\\br\\b|responder|responding|benefit|dcb", y)] <- "Responder"
  out[grepl("stable disease|\\bsd\\b|stable", y)] <- "StableDisease"
  out[grepl("progressive disease|\\bpd\\b|non responder|nonresponder|\\bnr\\b|no response|progression", y)] <- "NonResponder"

  out
}

find_metadata_candidates <- function(tables, selected_expr_sheet) {
  candidates <- list()
  for (nm in names(tables)) {
    if (identical(nm, selected_expr_sheet)) next
    df <- tables[[nm]]
    if (!is.data.frame(df) || nrow(df) < 2 || ncol(df) < 2) next
    df$.sheet_name <- NULL
    cn_norm <- normalize_name(colnames(df))
    sample_hits <- grep("sample|patient|title|geo|id", cn_norm, value = FALSE)
    response_hits <- grep("response|responder|benefit|recist|clinical|best", cn_norm, value = FALSE)
    if (length(sample_hits) > 0 && length(response_hits) > 0) {
      candidates[[length(candidates) + 1]] <- list(sheet = nm, df = df, sample_col = colnames(df)[sample_hits[1]], response_col = colnames(df)[response_hits[1]])
    }
  }
  candidates
}

state_scores$SampleKey <- make_key(state_scores$Sample)
state_scores$ResponseRaw <- NA_character_
state_scores$ResponseGroup <- NA_character_
state_scores$ResponseSource <- NA_character_

meta_candidates <- if (grepl("\\.xlsx$|\\.xls$", gse78220_expr_file, ignore.case = TRUE)) {
  find_metadata_candidates(all_tables, expr_sheet)
} else list()

response_audit <- data.frame(
  Status = "No response metadata matched automatically",
  CandidateSheet = NA_character_,
  SampleColumn = NA_character_,
  ResponseColumn = NA_character_,
  MatchedSamples = 0,
  stringsAsFactors = FALSE
)

if (length(meta_candidates) > 0) {
  best_join <- NULL
  best_n <- -1
  best_info <- NULL
  for (cand in meta_candidates) {
    md <- cand$df
    md$SampleKey <- make_key(md[[cand$sample_col]])
    md$ResponseRaw_tmp <- as.character(md[[cand$response_col]])
    md$ResponseGroup_tmp <- map_response_group(md$ResponseRaw_tmp)
    md2 <- md[, c("SampleKey", "ResponseRaw_tmp", "ResponseGroup_tmp"), drop = FALSE]
    md2 <- md2[!duplicated(md2$SampleKey), , drop = FALSE]
    joined <- dplyr::left_join(state_scores, md2, by = "SampleKey")
    n_match <- sum(!is.na(joined$ResponseGroup_tmp))
    if (n_match > best_n) {
      best_n <- n_match
      best_join <- joined
      best_info <- cand
    }
  }
  if (!is.null(best_join) && best_n > 0) {
    state_scores$ResponseRaw <- best_join$ResponseRaw_tmp
    state_scores$ResponseGroup <- best_join$ResponseGroup_tmp
    state_scores$ResponseSource <- paste0("workbook sheet: ", best_info$sheet)
    response_audit <- data.frame(
      Status = "Matched response metadata from workbook",
      CandidateSheet = best_info$sheet,
      SampleColumn = best_info$sample_col,
      ResponseColumn = best_info$response_col,
      MatchedSamples = best_n,
      stringsAsFactors = FALSE
    )
  }
}

## If workbook metadata did not work, try a separate local metadata file.
## Try all sample-column x response-column combinations and select
## the combination with the largest number of matched non-missing response groups.
if (all(is.na(state_scores$ResponseGroup))) {
  local_meta_candidates <- c(
    file.path(raw_dir, "external_validation/GSE78220/GSE78220_clinical_metadata_from_GEO.csv"),
    file.path(raw_dir, "external_validation/GSE78220/GSE78220_pData_raw_from_GEO.csv"),
    file.path(raw_dir, "external_validation/GSE78220/GSE78220_metadata.csv"),
    file.path(raw_dir, "external_validation/GSE78220/GSE78220_clinical_metadata.csv"),
    file.path(raw_dir, "external_validation/GSE78220/GSE78220_response_metadata.csv"),
    file.path(processed_dir, "GSE78220_pheno_from_GEOquery.csv"),
    file.path(processed_dir, "GSE78220_pheno.csv"),
    list.files(project_dir, pattern = "GSE78220.*(metadata|clinical|response|pheno|pData).*\\.(csv|tsv|txt)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  )
  local_meta_file <- find_first_existing(local_meta_candidates)

  if (!is.na(local_meta_file)) {
    message("Trying local response/phenotype metadata: ", local_meta_file)
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", local_meta_file, ignore.case = TRUE), "\t", ",")
    md <- tryCatch(
      read.table(local_meta_file, header = TRUE, sep = sep, quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )

    if (!is.null(md)) {
      cn <- colnames(md)
      cn_norm <- normalize_name(cn)

      ## Candidate sample columns. Try all because GEO_accession may appear before Sample,
      ## while expression columns are usually Pt1/Pt2/... rather than GSM IDs.
      sample_hits <- grep("sample|patient|title|geo|accession|id", cn_norm, value = FALSE)
      response_hits <- grep("anti.*pd.*response|response|responder|benefit|recist|clinical|best", cn_norm, value = FALSE)

      ## Prefer columns with clear names, but still test every candidate.
      sample_priority <- order(!grepl("^sample$|patient_id|patient|title", cn_norm[sample_hits]), cn_norm[sample_hits])
      response_priority <- order(!grepl("anti.*pd.*response|response", cn_norm[response_hits]), cn_norm[response_hits])
      sample_hits <- sample_hits[sample_priority]
      response_hits <- response_hits[response_priority]

      if (length(sample_hits) > 0 && length(response_hits) > 0) {
        candidate_audit <- list()
        best_join <- NULL
        best_n <- -1
        best_info <- NULL

        for (si in sample_hits) {
          for (ri in response_hits) {
            for (key_mode in c("strict", "loose")) {
              md_tmp <- md
              ss_tmp <- state_scores

              if (key_mode == "strict") {
                md_tmp$SampleKey <- make_key(md_tmp[[cn[si]]])
                ss_tmp$SampleKey <- make_key(ss_tmp$Sample)
              } else {
                md_tmp$SampleKey <- make_key_loose(md_tmp[[cn[si]]])
                ss_tmp$SampleKey <- make_key_loose(ss_tmp$Sample)
              }

              md_tmp$ResponseRaw_tmp <- as.character(md_tmp[[cn[ri]]])
              md_tmp$ResponseGroup_tmp <- map_response_group(md_tmp$ResponseRaw_tmp)

              ## Part A intentionally matches response metadata only.
              ## Overall-survival fields are parsed once, using stricter rules, in Part B.
              md2 <- md_tmp[, c("SampleKey", "ResponseRaw_tmp", "ResponseGroup_tmp"), drop = FALSE]
              md2 <- md2[!is.na(md2$SampleKey) & md2$SampleKey != "", , drop = FALSE]
              md2 <- md2[!duplicated(md2$SampleKey), , drop = FALSE]

              joined <- dplyr::left_join(ss_tmp, md2, by = "SampleKey")
              n_match_response <- sum(!is.na(joined$ResponseGroup_tmp))
              n_match_any <- sum(!is.na(joined$ResponseRaw_tmp))

              candidate_audit[[length(candidate_audit) + 1]] <- data.frame(
                MetadataFile = local_meta_file,
                SampleColumn = cn[si],
                ResponseColumn = cn[ri],
                KeyMode = key_mode,
                MatchedResponseGroups = n_match_response,
                MatchedRawResponses = n_match_any,
                stringsAsFactors = FALSE
              )

              if (n_match_response > best_n) {
                best_n <- n_match_response
                best_join <- joined
                best_info <- list(sample_col = cn[si], response_col = cn[ri], key_mode = key_mode)
              }
            }
          }
        }

        candidate_audit_df <- dplyr::bind_rows(candidate_audit) %>% arrange(desc(MatchedResponseGroups), desc(MatchedRawResponses))
        safe_write_csv(candidate_audit_df, file.path(out_table_dir, "GSE78220_external_local_metadata_matching_candidate_audit.csv"))

        if (!is.null(best_join) && best_n > 0) {
          state_scores$SampleKey <- best_join$SampleKey
          state_scores$ResponseRaw <- best_join$ResponseRaw_tmp
          state_scores$ResponseGroup <- best_join$ResponseGroup_tmp
          state_scores$ResponseSource <- local_meta_file
          response_audit <- data.frame(
            Status = "Matched response metadata from local file",
            CandidateSheet = basename(local_meta_file),
            SampleColumn = best_info$sample_col,
            ResponseColumn = best_info$response_col,
            KeyMode = best_info$key_mode,
            MatchedSamples = best_n,
            stringsAsFactors = FALSE
          )
        } else {
          warning("Local metadata file was found but no response groups matched. Check GSE78220_external_local_metadata_matching_candidate_audit.csv")
        }
      }
    }
  }
}

## Keep ResponseGroup column even if not found, because downstream 02 expects it.
response_audit$n_samples <- nrow(state_scores)
response_audit$n_response_nonmissing <- sum(!is.na(state_scores$ResponseGroup))
safe_write_csv(response_audit, file.path(out_table_dir, "GSE78220_external_response_matching_audit.csv"))
safe_write_csv(as.data.frame(table(ResponseGroup = state_scores$ResponseGroup, useNA = "ifany")), file.path(out_table_dir, "GSE78220_external_response_group_counts.csv"))

if (all(is.na(state_scores$ResponseGroup))) {
  warning(
    "No ResponseGroup could be matched automatically. The state-score table was created, but response/AUC analyses in 02 will require response metadata."
  )
}

############################################################
## 9. Export standardized table required by 02
############################################################
state_scores_out <- state_scores[, c(
  "Sample",
  state_cols_clean,
  "Dominant_state",
  "ResponseRaw",
  "ResponseGroup",
  "ResponseSource"
), drop = FALSE]

safe_write_csv(state_scores_out, file.path(out_table_dir, "GSE78220_external_state_scores.csv"))
safe_write_csv(state_scores_out, file.path(out_table_dir, "GSE78220_external_state_scores_display_labels.csv"))
safe_save_rds(state_scores_out, file.path(processed_dir, "GSE78220_external_state_scores.rds"))

############################################################
## 10. Quick figures for audit
############################################################
state_long <- state_scores_out %>%
  dplyr::select(Sample, all_of(state_cols_clean)) %>%
  pivot_longer(cols = all_of(state_cols_clean), names_to = "State", values_to = "Score")

sample_order <- state_scores_out %>%
  arrange(Dominant_state) %>%
  pull(Sample)

state_long$Sample <- factor(state_long$Sample, levels = sample_order)
state_long$State <- factor(state_long$State, levels = rev(state_cols_clean), labels = rev(state_label[state_cols_clean]))

p_heat <- ggplot(state_long, aes(x = Sample, y = State, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    name = "State\nscore"
  ) +
  theme_icb(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10.25), axis.title = element_blank(), panel.grid = element_blank()) +
  labs(title = "GSE78220 external state scores")

safe_ggsave(file.path(out_fig_dir, "GSE78220_external_state_score_heatmap.pdf"), p_heat, width = 12, height = 4.8)
safe_ggsave(file.path(out_fig_dir, "GSE78220_external_state_score_heatmap.png"), p_heat, width = 12, height = 4.8)

pca_mat <- as.matrix(state_scores_out[, state_cols_clean, drop = FALSE])
rownames(pca_mat) <- state_scores_out$Sample
if (nrow(pca_mat) >= 3 && ncol(pca_mat) >= 2) {
  pca <- prcomp(pca_mat, scale. = TRUE)
  pca_var <- pca$sdev^2 / sum(pca$sdev^2)
  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Dominant_state = state_scores_out$Dominant_state[match(rownames(pca$x), state_scores_out$Sample)],
    ResponseGroup = state_scores_out$ResponseGroup[match(rownames(pca$x), state_scores_out$Sample)],
    stringsAsFactors = FALSE
  )

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Dominant_state, label = Sample)) +
    geom_point(size = 2.8) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 100, show.legend = FALSE) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    theme_icb(base_size = 12.5) +
    labs(
      title = "PCA of GSE78220 external state scores",
      x = paste0("PC1 (", round(pca_var[1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(pca_var[2] * 100, 1), "%)"),
      color = "Dominant state"
    )

  safe_ggsave(file.path(out_fig_dir, "GSE78220_external_state_score_PCA.pdf"), p_pca, width = 8, height = 6)
  safe_ggsave(file.path(out_fig_dir, "GSE78220_external_state_score_PCA.png"), p_pca, width = 8, height = 6)
}

############################################################
## 11. Final output inventory
############################################################
output_inventory <- data.frame(
  File = c(
    file.path(out_table_dir, "GSE78220_external_state_scores.csv"),
    file.path(out_table_dir, "GSE78220_external_raw_gene_set_scores.csv"),
    file.path(out_table_dir, "GSE78220_external_gene_presence_all_gene_sets.csv"),
    file.path(out_table_dir, "GSE78220_external_response_matching_audit.csv"),
    file.path(out_table_dir, "GSE78220_external_response_group_counts.csv")
  ),
  Exists = file.exists(c(
    file.path(out_table_dir, "GSE78220_external_state_scores.csv"),
    file.path(out_table_dir, "GSE78220_external_raw_gene_set_scores.csv"),
    file.path(out_table_dir, "GSE78220_external_gene_presence_all_gene_sets.csv"),
    file.path(out_table_dir, "GSE78220_external_response_matching_audit.csv"),
    file.path(out_table_dir, "GSE78220_external_response_group_counts.csv")
  )),
  stringsAsFactors = FALSE
)
safe_write_csv(output_inventory, file.path(out_table_dir, "GSE78220_external_state_scoring_output_inventory.csv"))

save_session_info()

message("GSE78220 external state scoring completed.")
message("Main output for 02: ", file.path(out_table_dir, "GSE78220_external_state_scores.csv"))


############################################################
## Checkpoint after Part A
############################################################
partA_score_file <- file.path(project_dir, "results/tables/external_validation/GSE78220_external_state_scores.csv")
if (!file.exists(partA_score_file)) {
  stop("Part A finished but the expected GSE78220_external_state_scores.csv file was not found.")
}
partA_scores <- read.csv(partA_score_file, check.names = FALSE)
message("Part A response-group summary:")
print(table(partA_scores$ResponseGroup, useNA = "ifany"))
if (any(is.na(partA_scores$ResponseGroup))) {
  unmatched_response <- partA_scores[is.na(partA_scores$ResponseGroup), , drop = FALSE]
  write.csv(
    unmatched_response,
    file.path(project_dir, "results/tables/external_validation/GSE78220_samples_with_missing_ResponseGroup_AFTER_PART_A.csv"),
    row.names = FALSE
  )
  warning(
    "Some GSE78220 samples have missing ResponseGroup. This is allowed; they will be excluded from response/AUC analyses. See GSE78220_samples_with_missing_ResponseGroup_AFTER_PART_A.csv."
  )
}

############################################################
## PART B. Response/OS validation and Figure 7
############################################################

############################################################
## Part B: response association, baseline-only survival analysis, and Figure 7
##
## Purpose:
## 1) Load the standardized GSE78220 state scores generated in Part A.
## 2) Retrieve/parse phenotype metadata with overall-survival information.
## 3) Merge OS time/event fields using patient identifiers.
## 4) Perform baseline-only patient-level Cox analyses for the four states.
## 5) Generate descriptive response/AUC, Cox, correlation panels, and Figure 7.
##
## Canonical state names from Step 01:
##   Immune_defective_Cold
##   Myeloid_Treg_Immunosuppressive
##   Tumor_dedifferentiation_Stromal_remodeling
##   Melanocytic_Differentiation
##   Dominant_state / Dominant_final_state_CLEAN
##
## Plot font:
##   sans
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

############################################################
## 0. Project directory and folders
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- getwd()
}

message("Project directory: ", project_dir)

data_processed_dir <- file.path(project_dir, "data_processed")
out_table_dir <- file.path(project_dir, "results/tables/external_validation")
out_fig_dir   <- file.path(project_dir, "results/figures/external_validation")
log_dir       <- file.path(project_dir, "logs")

dir.create(data_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "stringr", "patchwork",
  "pROC", "survival", "data.table", "scales"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
    ". Restore the repository environment (for example with renv::restore()) before running this script."
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(pROC)
  library(survival)
  library(data.table)
  library(scales)
})

theme_icb <- function(base_size = 12.5, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2.5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 1.5),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1.5),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6)
    )
}

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
  message("Saved table: ", file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
  message("Saved RDS: ", file)
}

safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved figure: ", file)
}

save_session_info <- function(script_name = "02_GSE78220_state_scoring_response_survival_analysis") {
  log_file <- file.path(project_dir, "logs", paste0("sessionInfo_", script_name, ".txt"))
  sink(log_file)
  print(sessionInfo())
  sink()
  message("sessionInfo saved to: ", log_file)
}

find_files_recursive <- function(root, pattern) {
  if (!dir.exists(root)) return(character(0))
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}

find_first_existing <- function(paths) {
  paths <- unique(paths)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

filter_real_pheno_candidates <- function(paths) {
  paths <- unique(paths)
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(character(0))

  bn <- basename(paths)
  bn_norm <- tolower(bn)

  ## Exclude files generated by our own analysis pipeline that are not phenotype metadata.
  bad_pattern <- paste(
    c(
      "metadata_columns",
      "detected_columns",
      "merge_check",
      "outputs_inventory",
      "input_file_inventory",
      "parse_summary",
      "patient_level",
      "cox_input",
      "overall_survival_cox",
      "nonresponse_auc",
      "correlation",
      "state_scores",
      "external_state",
      "standardized",
      "with_os",
      "remerged",
      "fixed",
      "locked",
      "sessioninfo",
      "gene_presence",
      "dominant_state",
      "response_counts"
    ),
    collapse = "|"
  )

  paths[!grepl(bad_pattern, bn_norm)]
}

looks_like_pheno_table <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(FALSE)
  if (nrow(df) < 5 || ncol(df) < 3) return(FALSE)

  cn <- normalize_name(colnames(df))

  ## A column-inventory table often has only columns such as column/type/n_missing.
  inventory_like <- any(cn %in% c("column", "columns", "field")) &&
    ncol(df) <= 5 &&
    !any(grepl("title|geo_accession|patient|sample|survival|vital|status|os|response|characteristics", cn))

  if (inventory_like) return(FALSE)

  score <- sum(grepl("title|geo_accession|patient|sample|survival|overall|vital|status|response|os|characteristics", cn))

  score >= 2
}

normalize_name <- function(x) {
  x <- as.character(x)
  x <- gsub("\\ufeff", "", x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA)
  x[1]
}

coalesce_joined_column <- function(df, target, candidates) {
  hits <- candidates[candidates %in% colnames(df)]

  if (length(hits) == 0) {
    df[[target]] <- NA
    return(df)
  }

  out <- df[[hits[1]]]

  if (length(hits) >= 2) {
    for (h in hits[-1]) {
      out <- dplyr::coalesce(out, df[[h]])
    }
  }

  df[[target]] <- out
  df
}

############################################################
## 2. Unified canonical state naming
############################################################

bulk_state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

## Compatibility: rename legacy column names generated by earlier scripts to the
## current CLEAN names used after Step 01.
legacy_to_clean_state_cols <- c(
  "Immune_defective_cold" = "Immune_defective_Cold",
  "Myeloid_Treg_immunosuppressive" = "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiated" = "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_differentiated" = "Melanocytic_Differentiation"
)

bulk_state_display <- c(
  "Immune_defective_Cold" = "immune-defective/\ncold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic\ndifferentiation"
)

bulk_state_display_one_line <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

## Public Figure 7 display-label lock.
expected_public_state_labels <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation"
)
if (!identical(unname(bulk_state_display_one_line[bulk_state_cols]), expected_public_state_labels)) {
  stop("Public Figure 7 state-label lock failed: display labels do not match the frozen manuscript convention.")
}

bulk_state_order <- bulk_state_cols

state_colors <- c(
  "Immune_defective_Cold" = "#4DBBD5",
  "Myeloid_Treg_Immunosuppressive" = "#00A087",
  "Tumor_dedifferentiation_Stromal_remodeling" = "#E64B35",
  "Melanocytic_Differentiation" = "#3C5488"
)

heatmap_low  <- "#3B82F6"
heatmap_mid  <- "white"
heatmap_high <- "#EF4444"

response_group_colors <- c(
  "Responder" = "#E64B35",
  "StableDisease" = "#00A087",
  "NonResponder" = "#4DBBD5"
)

############################################################
## 3. Load standardized GSE78220 state-score table
############################################################

gse78220_score_file <- file.path(
  out_table_dir,
  "GSE78220_external_state_scores.csv"
)

if (!file.exists(gse78220_score_file)) {
  stop(
    "Canonical GSE78220 state-score table not found: ", gse78220_score_file,
    ". Part A must complete successfully before Part B."
  )
}

message("Using GSE78220 state-score table: ", gse78220_score_file)

final78220 <- read.csv(
  gse78220_score_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## Standardize legacy state-score column names to CLEAN names when needed.
for (old_nm in names(legacy_to_clean_state_cols)) {
  new_nm <- unname(legacy_to_clean_state_cols[[old_nm]])
  if (old_nm %in% colnames(final78220) && !new_nm %in% colnames(final78220)) {
    colnames(final78220)[colnames(final78220) == old_nm] <- new_nm
  }
}

if ("Dominant_final_state_CLEAN" %in% colnames(final78220) && !"Dominant_state" %in% colnames(final78220)) {
  final78220$Dominant_state <- final78220$Dominant_final_state_CLEAN
}

missing_state_cols <- setdiff(bulk_state_cols, colnames(final78220))
if (length(missing_state_cols) > 0) {
  stop(
    "GSE78220 table lacks canonical state columns: ",
    paste(missing_state_cols, collapse = ", ")
  )
}

if (!"Sample" %in% colnames(final78220)) {
  final78220$Sample <- rownames(final78220)
}

if (!"ResponseGroup" %in% colnames(final78220)) {
  stop("GSE78220 state-score table lacks ResponseGroup. Please rerun 02 first.")
}

############################################################
## 4. Patient-ID normalization
############################################################

extract_patient_id <- function(x, collapse_AB_suffix = TRUE) {
  x0 <- as.character(x)
  x0 <- gsub("\\.baseline$|\\.pre$|\\.pretreatment$|\\.on$|\\.on_treatment$", "", x0, ignore.case = TRUE)
  x0 <- gsub("[[:space:]_]+", "", x0)

  ## Extract Pt-number with optional trailing letter
  m <- stringr::str_extract(x0, regex("pt\\d+[a-z]?", ignore_case = TRUE))
  m[is.na(m)] <- x0[is.na(m)]
  m <- toupper(m)

  if (collapse_AB_suffix) {
    m <- gsub("^(PT[0-9]+)[A-Z]$", "\\1", m)
  }

  m <- gsub("^PT", "Pt", m)
  m
}

final78220$PatientID_from_score <- extract_patient_id(final78220$Sample, collapse_AB_suffix = TRUE)

############################################################
## 5. Locate or retrieve GSE78220 phenotype metadata
############################################################

## This script first searches locally.
## If no local phenotype metadata is found, it tries GEOquery download.
## You may set this to an exact local phenotype file path if automatic search is wrong.
manual_pheno_file <- NULL
allow_geoquery_download <- TRUE

parse_characteristics_columns <- function(pheno) {

  pheno <- as.data.frame(pheno, stringsAsFactors = FALSE)
  char_cols <- grep("^characteristics|characteristics_ch1", colnames(pheno), value = TRUE, ignore.case = TRUE)

  if (length(char_cols) == 0) return(pheno)

  for (cc in char_cols) {
    vals <- as.character(pheno[[cc]])
    has_colon <- grepl(":", vals)
    if (!any(has_colon, na.rm = TRUE)) next

    keys <- trimws(sub(":.*$", "", vals))
    vals2 <- trimws(sub("^[^:]+:\\s*", "", vals))

    key_norm <- normalize_name(keys)

    for (ky in unique(key_norm[!is.na(key_norm) & key_norm != ""])) {
      new_col <- paste0("char_", ky)
      if (!new_col %in% colnames(pheno)) {
        pheno[[new_col]] <- NA_character_
      }
      idx <- which(key_norm == ky)
      pheno[[new_col]][idx] <- vals2[idx]
    }
  }

  pheno
}

load_any_pheno_from_rds <- function(file) {
  obj <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(obj)) return(NULL)

  candidates <- list()

  scan_obj <- function(x, nm = "root") {
    if (is.data.frame(x)) {
      cn <- normalize_name(colnames(x))
      score <- sum(grepl("title|geo_accession|patient|survival|overall|vital|status|response|os", cn))
      if (score >= 2) {
        candidates[[length(candidates) + 1]] <<- list(name = nm, data = x, score = score)
      }
    } else if (is.list(x)) {
      for (n in names(x)) {
        scan_obj(x[[n]], paste0(nm, "$", n))
      }
    }
  }

  scan_obj(obj, basename(file))

  if (length(candidates) == 0) return(NULL)

  scores <- sapply(candidates, function(z) z$score)
  best <- candidates[[which.max(scores)]]
  message("Found phenotype-like data frame in RDS: ", file, " -> ", best$name)
  as.data.frame(best$data, stringsAsFactors = FALSE)
}

local_pheno_candidates_raw <- c(
  file.path(project_dir, "data_processed/GSE78220_pheno.rds"),
  file.path(project_dir, "data_processed/GSE78220_pheno.csv"),
  file.path(project_dir, "data_processed/GSE78220_metadata.csv"),
  file.path(project_dir, "data_processed/GSE78220_clinical_metadata.csv"),
  file.path(project_dir, "data_raw/GSE78220/GSE78220_pheno.csv"),
  file.path(project_dir, "data_raw/GSE78220/GSE78220_metadata.csv"),
  file.path(project_dir, "data_raw/GSE78220/GSE78220_clinical_metadata.csv"),
  find_files_recursive(project_dir, "GSE78220.*(pheno|metadata|clinical|pData).*\\.(csv|tsv|txt|rds)$"),
  find_files_recursive(project_dir, "external_validation_results.*\\.rds$")
)

if (!is.null(manual_pheno_file) && file.exists(manual_pheno_file)) {
  local_pheno_candidates <- manual_pheno_file
} else {
  local_pheno_candidates <- filter_real_pheno_candidates(local_pheno_candidates_raw)
}

candidate_inventory <- data.frame(
  Path = unique(local_pheno_candidates_raw),
  Exists = file.exists(unique(local_pheno_candidates_raw)),
  UsedAsCandidate = unique(local_pheno_candidates_raw) %in% local_pheno_candidates,
  stringsAsFactors = FALSE
)

safe_write_csv(
  candidate_inventory,
  file.path(out_table_dir, "GSE78220_local_pheno_candidate_inventory.csv"),
  row.names = FALSE
)

local_pheno_file <- find_first_existing(local_pheno_candidates)

gse78220_pheno <- NULL

if (!is.na(local_pheno_file)) {
  message("Local GSE78220 phenotype candidate found: ", local_pheno_file)

  if (grepl("\\.rds$", local_pheno_file, ignore.case = TRUE)) {
    gse78220_pheno <- load_any_pheno_from_rds(local_pheno_file)
  } else {
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", local_pheno_file, ignore.case = TRUE), "\t", ",")
    tmp_pheno <- tryCatch(
      read.table(
        local_pheno_file,
        header = TRUE,
        sep = sep,
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      error = function(e) {
        message("Failed to read local candidate: ", e$message)
        NULL
      }
    )

    if (!is.null(tmp_pheno) && looks_like_pheno_table(tmp_pheno)) {
      gse78220_pheno <- tmp_pheno
    } else {
      message("Local candidate does not look like a phenotype metadata table. Skipping: ", local_pheno_file)
      gse78220_pheno <- NULL
    }
  }

  if (!looks_like_pheno_table(gse78220_pheno)) {
    message("Selected local candidate failed phenotype-table check. Will try GEOquery if allowed.")
    gse78220_pheno <- NULL
  }
}

## If still not found, try an in-memory object from current R session.
if (is.null(gse78220_pheno) && exists("gse78220_pheno", envir = .GlobalEnv)) {
  message("Using in-memory object: gse78220_pheno")
  gse78220_pheno <- get("gse78220_pheno", envir = .GlobalEnv)
}

## If no local data, use GEOquery.
if (is.null(gse78220_pheno) && allow_geoquery_download) {

  message("No local GSE78220 phenotype table found. Trying GEOquery download.")

  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop(
      "Package 'GEOquery' is required when local phenotype metadata are unavailable. ",
      "Restore the repository environment (for example with renv::restore()) and rerun."
    )
  }

  suppressPackageStartupMessages(library(GEOquery))

  gset <- tryCatch(
    GEOquery::getGEO("GSE78220", GSEMatrix = TRUE),
    error = function(e) {
      message("GEOquery failed: ", e$message)
      NULL
    }
  )

  if (!is.null(gset)) {
    if (is.list(gset)) {
      idx <- which.max(sapply(gset, ncol))
      gse78220_pheno <- Biobase::pData(gset[[idx]])
    } else {
      gse78220_pheno <- Biobase::pData(gset)
    }

    safe_save_rds(
      gse78220_pheno,
      file.path(data_processed_dir, "GSE78220_pheno_from_GEOquery.rds")
    )
    safe_write_csv(
      gse78220_pheno,
      file.path(data_processed_dir, "GSE78220_pheno_from_GEOquery.csv"),
      row.names = FALSE
    )
  }
}

if (is.null(gse78220_pheno)) {
  stop(
    "Unable to find or download GSE78220 phenotype metadata. ",
    "Please place a phenotype file containing patient/sample ID, OS time, and OS status into data_processed/ ",
    "with a name such as GSE78220_pheno.csv, then rerun this script."
  )
}

gse78220_pheno <- parse_characteristics_columns(gse78220_pheno)

safe_write_csv(
  gse78220_pheno,
  file.path(out_table_dir, "GSE78220_pheno_metadata_raw_or_parsed.csv"),
  row.names = FALSE
)

############################################################
## 6. Detect patient, OS-time, and OS-event columns
############################################################

detect_patient_col <- function(pheno) {
  cn <- colnames(pheno)
  cn_norm <- normalize_name(cn)

  priority_patterns <- c(
    "^title$",
    "patient.*id",
    "^patient$",
    "char_patient",
    "char_patient_id",
    "sample.*id",
    "geo_accession"
  )

  for (pat in priority_patterns) {
    hit <- cn[grepl(pat, cn_norm)]
    if (length(hit) > 0) return(hit[1])
  }

  cn[1]
}

detect_os_time_col <- function(pheno) {
  cn <- colnames(pheno)
  cn_norm <- normalize_name(cn)

  ## Prefer OS/overall survival columns, avoid status/event/response.
  preferred <- cn[
    grepl("overall.*survival|survival.*time|os.*time|os_month|os_day|char_overall.*survival|char_survival", cn_norm) &
      !grepl("status|event|death|dead|vital|response|group", cn_norm)
  ]

  if (length(preferred) > 0) return(preferred[1])

  fallback <- cn[
    grepl("survival|os", cn_norm) &
      !grepl("status|event|death|dead|vital|response|group", cn_norm)
  ]

  if (length(fallback) > 0) return(fallback[1])

  NA_character_
}

detect_os_event_col <- function(pheno) {
  cn <- colnames(pheno)
  cn_norm <- normalize_name(cn)

  preferred <- cn[
    grepl("vital.*status|survival.*status|os.*event|death|dead|event|char_vital|char_status", cn_norm) &
      !grepl("time|month|day|response|group", cn_norm)
  ]

  if (length(preferred) > 0) return(preferred[1])

  fallback <- cn[
    grepl("status|event|vital", cn_norm) &
      !grepl("time|month|day|response|group", cn_norm)
  ]

  if (length(fallback) > 0) return(fallback[1])

  NA_character_
}

## Manual override if automatic detection is wrong.
## Example:
## manual_patient_col <- "title"
## manual_os_time_col <- "char_overall_survival_months"
## manual_os_event_col <- "char_vital_status"
manual_patient_col <- NULL
manual_os_time_col <- NULL
manual_os_event_col <- NULL

patient_col <- if (!is.null(manual_patient_col)) manual_patient_col else detect_patient_col(gse78220_pheno)
os_time_col <- if (!is.null(manual_os_time_col)) manual_os_time_col else detect_os_time_col(gse78220_pheno)
os_event_col <- if (!is.null(manual_os_event_col)) manual_os_event_col else detect_os_event_col(gse78220_pheno)

message("Detected patient column: ", patient_col)
message("Detected OS time column: ", os_time_col)
message("Detected OS event/status column: ", os_event_col)

detected_cols <- data.frame(
  Field = c("PatientID", "OS_time", "OS_event"),
  Column = c(patient_col, os_time_col, os_event_col),
  stringsAsFactors = FALSE
)

safe_write_csv(
  detected_cols,
  file.path(out_table_dir, "GSE78220_OS_metadata_detected_columns.csv"),
  row.names = FALSE
)

if (is.na(os_time_col) || is.na(os_event_col)) {
  stop(
    "Unable to automatically detect OS time or OS event/status column. ",
    "Please open GSE78220_pheno_metadata_raw_or_parsed.csv, identify the columns, ",
    "then set manual_os_time_col and manual_os_event_col near the top of section 6."
  )
}

parse_numeric_time <- function(x) {
  x <- as.character(x)
  x <- gsub(",", "", x)
  x <- gsub("[^0-9\\.\\-]+", "", x)
  suppressWarnings(as.numeric(x))
}

parse_event_status <- function(x) {
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA_integer_, ifelse(x > 0, 1L, 0L)))
  }

  x0 <- tolower(trimws(as.character(x)))

  out <- rep(NA_integer_, length(x0))

  out[x0 %in% c("1", "dead", "deceased", "death", "event", "yes", "true", "died")] <- 1L
  out[x0 %in% c("0", "alive", "living", "censored", "no", "false")] <- 0L

  ## Handles strings like "dead: 1", "vital status: deceased"
  out[is.na(out) & grepl("dead|deceased|death|died", x0)] <- 1L
  out[is.na(out) & grepl("alive|living|censor", x0)] <- 0L

  ## Last fallback numeric conversion
  numeric_x <- suppressWarnings(as.numeric(gsub("[^0-9\\.\\-]+", "", x0)))
  out[is.na(out) & !is.na(numeric_x)] <- ifelse(numeric_x[is.na(out) & !is.na(numeric_x)] > 0, 1L, 0L)

  out
}

pheno_os <- gse78220_pheno %>%
  mutate(
    PatientID_from_pheno_raw = as.character(.data[[patient_col]]),
    PatientID_from_pheno = extract_patient_id(PatientID_from_pheno_raw, collapse_AB_suffix = TRUE),
    OS_time_raw = as.character(.data[[os_time_col]]),
    OS_event_raw = as.character(.data[[os_event_col]]),
    OS_time = parse_numeric_time(OS_time_raw),
    OS_event = parse_event_status(OS_event_raw)
  ) %>%
  dplyr::select(
    PatientID_from_pheno,
    PatientID_from_pheno_raw,
    OS_time,
    OS_event,
    OS_time_raw,
    OS_event_raw,
    everything()
  )

pheno_os_summary <- pheno_os %>%
  summarise(
    n_rows = n(),
    n_unique_patient_ids = n_distinct(PatientID_from_pheno),
    n_with_OS_time = sum(!is.na(OS_time)),
    n_with_OS_event = sum(!is.na(OS_event)),
    n_events = sum(OS_event == 1, na.rm = TRUE),
    n_censored = sum(OS_event == 0, na.rm = TRUE)
  )

safe_write_csv(
  pheno_os_summary,
  file.path(out_table_dir, "GSE78220_OS_metadata_parse_summary.csv"),
  row.names = FALSE
)

safe_write_csv(
  pheno_os,
  file.path(out_table_dir, "GSE78220_OS_metadata_parsed.csv"),
  row.names = FALSE
)

############################################################
## 7. Collapse metadata to one row per patient and merge
############################################################

pheno_os_patient <- pheno_os %>%
  group_by(PatientID_from_pheno) %>%
  summarise(
    PatientID_from_pheno_raw = first_non_missing(PatientID_from_pheno_raw),
    OS_time = first_non_missing(OS_time),
    OS_event = first_non_missing(OS_event),
    OS_time_raw = first_non_missing(OS_time_raw),
    OS_event_raw = first_non_missing(OS_event_raw),
    n_metadata_rows = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(PatientID_from_pheno) & PatientID_from_pheno != "")

safe_write_csv(
  pheno_os_patient,
  file.path(out_table_dir, "GSE78220_OS_metadata_patient_level.csv"),
  row.names = FALSE
)

merged78220 <- final78220 %>%
  left_join(
    pheno_os_patient,
    by = c("PatientID_from_score" = "PatientID_from_pheno")
  )

## Resolve possible .x/.y suffixes after left_join.
## This happens when the standardized score table already contains OS_time/OS_event
## columns, while the phenotype table also provides these columns.
merged78220 <- coalesce_joined_column(
  merged78220,
  "OS_time",
  c("OS_time", "OS_time.y", "OS_time_pheno", "OS_time.x")
)

merged78220 <- coalesce_joined_column(
  merged78220,
  "OS_event",
  c("OS_event", "OS_event.y", "OS_event_pheno", "OS_event.x")
)

merged78220 <- coalesce_joined_column(
  merged78220,
  "OS_time_raw",
  c("OS_time_raw", "OS_time_raw.y", "OS_time_raw_pheno", "OS_time_raw.x")
)

merged78220 <- coalesce_joined_column(
  merged78220,
  "OS_event_raw",
  c("OS_event_raw", "OS_event_raw.y", "OS_event_raw_pheno", "OS_event_raw.x")
)

merged78220$OS_time <- suppressWarnings(as.numeric(merged78220$OS_time))
merged78220$OS_event <- suppressWarnings(as.numeric(merged78220$OS_event))

merge_check <- merged78220 %>%
  transmute(
    Sample,
    PatientID_from_score,
    matched_OS = !is.na(OS_time) & !is.na(OS_event),
    OS_time,
    OS_event,
    ResponseGroup
  )

safe_write_csv(
  merge_check,
  file.path(out_table_dir, "GSE78220_OS_merge_check_by_sample.csv"),
  row.names = FALSE
)

safe_write_csv(
  merged78220,
  file.path(out_table_dir, "GSE78220_external_state_scores_WITH_OS.csv"),
  row.names = FALSE
)

safe_save_rds(
  merged78220,
  file.path(data_processed_dir, "GSE78220_external_state_scores_WITH_OS.rds")
)

message("GSE78220 OS merge summary:")
print(table(merge_check$matched_OS, useNA = "ifany"))

############################################################
## 8. Patient-level Cox input: strict baseline-only analysis
############################################################

## The survival association is interpreted as a baseline prognostic analysis.
## Therefore, on-treatment expression samples are excluded before patient-level
## aggregation. Current GSE78220 expression sample names use the suffix
## ".baseline" for pretreatment samples. Additional common pretreatment suffixes
## are accepted for portability, but on-treatment samples are never included.
is_baseline_sample <- function(x) {
  x <- as.character(x)
  grepl("(\\.baseline|\\.pre|\\.pretreatment)$", x, ignore.case = TRUE)
}

cox_nonbaseline_excluded <- merged78220 %>%
  filter(!is_baseline_sample(Sample)) %>%
  select(Sample, PatientID_from_score, OS_time, OS_event, ResponseGroup)

safe_write_csv(
  cox_nonbaseline_excluded,
  file.path(out_table_dir, "GSE78220_Cox_excluded_nonbaseline_samples.csv"),
  row.names = FALSE
)

cox_input <- merged78220 %>%
  filter(is_baseline_sample(Sample)) %>%
  filter(!is.na(OS_time), !is.na(OS_event), OS_time > 0) %>%
  group_by(PatientID_from_score) %>%
  summarise(
    Sample = first_non_missing(Sample),
    ResponseGroup = first_non_missing(ResponseGroup),
    OS_time = first_non_missing(OS_time),
    OS_event = first_non_missing(OS_event),
    across(all_of(bulk_state_cols), ~ mean(.x, na.rm = TRUE)),
    n_score_rows = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(OS_time), !is.na(OS_event))

safe_write_csv(
  cox_input,
  file.path(out_table_dir, "GSE78220_Cox_input_patient_level_baseline_only.csv"),
  row.names = FALSE
)

cox_input_summary <- cox_input %>%
  summarise(
    n_patients = n(),
    n_events = sum(OS_event == 1, na.rm = TRUE),
    n_censored = sum(OS_event == 0, na.rm = TRUE),
    median_OS_time = median(OS_time, na.rm = TRUE)
  )

safe_write_csv(
  cox_input_summary,
  file.path(out_table_dir, "GSE78220_Cox_input_summary_baseline_only.csv"),
  row.names = FALSE
)

## Frozen reproducibility gate for the public release.
expected_cox_patients <- 25L
expected_cox_events <- 11L
if (nrow(cox_input) != expected_cox_patients ||
    sum(cox_input$OS_event == 1, na.rm = TRUE) != expected_cox_events) {
  stop(
    "Baseline-only Cox cohort does not match the frozen analysis set. Expected ",
    expected_cox_patients, " patients and ", expected_cox_events, " events; observed ",
    nrow(cox_input), " patients and ", sum(cox_input$OS_event == 1, na.rm = TRUE), " events. ",
    "Inspect GSE78220_Cox_excluded_nonbaseline_samples.csv and the OS merge audit."
  )
}

if (length(unique(cox_input$OS_event)) < 2) {
  stop("Baseline-only OS data do not contain both event and censored observations.")
}

############################################################
## 9. Cox regression for each canonical state
############################################################

cox_list <- list()

for (st in bulk_state_cols) {

  tmp <- cox_input %>%
    transmute(
      OS_time = suppressWarnings(as.numeric(OS_time)),
      OS_event = suppressWarnings(as.numeric(OS_event)),
      score = suppressWarnings(as.numeric(.data[[st]]))
    ) %>%
    filter(!is.na(OS_time), !is.na(OS_event), !is.na(score), OS_time > 0)

  if (nrow(tmp) >= 8 && length(unique(tmp$OS_event)) == 2) {
    fit <- survival::coxph(survival::Surv(OS_time, OS_event) ~ score, data = tmp)
    sm <- summary(fit)

    cox_list[[st]] <- data.frame(
      Dataset = "GSE78220",
      State = st,
      StateLabel = bulk_state_display_one_line[st],
      HR = sm$coefficients[1, "exp(coef)"],
      CI_lower = sm$conf.int[1, "lower .95"],
      CI_upper = sm$conf.int[1, "upper .95"],
      P_value = sm$coefficients[1, "Pr(>|z|)"],
      n_patients = nrow(tmp),
      n_events = sum(tmp$OS_event == 1),
      stringsAsFactors = FALSE
    )
  }
}

cox_table <- bind_rows(cox_list)

if (nrow(cox_table) == 0) {
  stop("Cox model failed for all state scores.")
}

cox_table <- cox_table %>%
  mutate(
    State = factor(State, levels = bulk_state_order),
    StateLabelPlot = factor(
      bulk_state_display[as.character(State)],
      levels = rev(bulk_state_display[bulk_state_order])
    ),
    Significant = !is.na(P_value) & P_value < 0.05,
    LabelText = paste0(
      "HR=", sprintf("%.2f", HR),
      "\nP=", format_p(P_value),
      ifelse(Significant, " *", "")
    )
  )

safe_write_csv(
  cox_table,
  file.path(out_table_dir, "GSE78220_overall_survival_cox_baseline_only.csv"),
  row.names = FALSE
)

############################################################
## 10. Build Figure 7 panels
############################################################

## Panel A: boxplot
box_df <- merged78220 %>%
  filter(ResponseGroup %in% c("Responder", "NonResponder")) %>%
  mutate(
    ResponseGroup = factor(ResponseGroup, levels = c("Responder", "NonResponder"))
  )

box_long <- box_df %>%
  dplyr::select(Sample, ResponseGroup, all_of(bulk_state_cols)) %>%
  pivot_longer(
    cols = all_of(bulk_state_cols),
    names_to = "State",
    values_to = "Score"
  ) %>%
  mutate(
    StateLabel = factor(bulk_state_display[State], levels = bulk_state_display[bulk_state_order])
  )

pA <- ggplot(box_long, aes(x = ResponseGroup, y = Score, fill = ResponseGroup)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.45, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.75, color = "black") +
  facet_wrap(~ StateLabel, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = response_group_colors, drop = FALSE) +
  labs(
    title = "GSE78220: state scores by response",
    x = NULL,
    y = "State score",
    fill = "Response group"
  ) +
  theme_icb(base_size = 12.5) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

safe_ggsave(
  file.path(out_fig_dir, "Figure7A_GSE78220_state_scores_by_response.pdf"),
  pA,
  width = 9.5,
  height = 7.2
)

safe_ggsave(
  file.path(out_fig_dir, "Figure7A_GSE78220_state_scores_by_response.png"),
  pA,
  width = 9.5,
  height = 7.2,
  dpi = 300
)

## Panel B: AUC with stratified bootstrap 95% CI
## GSE78220 binary response analysis uses only Responder vs NonResponder.
## StableDisease or missing response labels are excluded from this descriptive AUC analysis.
auc_df <- merged78220 %>%
  dplyr::filter(ResponseGroup %in% c("Responder", "NonResponder")) %>%
  dplyr::mutate(NonResponse = ifelse(ResponseGroup == "NonResponder", 1, 0))

auc_list <- list()

## Fixed seed for reproducible bootstrap confidence intervals.
auc_bootstrap_seed <- 20260506
auc_bootstrap_n <- 10000
set.seed(auc_bootstrap_seed)

for (st in bulk_state_cols) {
  x <- auc_df[[st]]
  keep <- is.finite(x) & !is.na(auc_df$NonResponse)

  if (sum(keep) >= 5 && length(unique(auc_df$NonResponse[keep])) == 2) {
    roc_obj <- pROC::roc(
      response = auc_df$NonResponse[keep],
      predictor = x[keep],
      levels = c(0, 1),   # 0 = Responder, 1 = NonResponder
      direction = "<",    # higher score indicates higher non-response tendency
      quiet = TRUE
    )

    auc_raw <- as.numeric(pROC::auc(roc_obj))

    ## pROC::ci.auc returns c(lower, median, upper).
    ## stratified bootstrap is used to preserve the response/non-response balance.
    ci_auc <- tryCatch(
      {
        as.numeric(pROC::ci.auc(
          roc_obj,
          conf.level = 0.95,
          method = "bootstrap",
          boot.n = auc_bootstrap_n,
          boot.stratified = TRUE,
          progress = "none"
        ))
      },
      error = function(e) {
        message("AUC CI failed for ", st, ": ", conditionMessage(e))
        c(NA_real_, NA_real_, NA_real_)
      }
    )

    auc_list[[st]] <- data.frame(
      Dataset = "GSE78220",
      State = st,
      StateLabel = bulk_state_display_one_line[st],
      AUC_raw_score_for_nonresponse = auc_raw,
      AUC_95CI_low = ci_auc[1],
      AUC_95CI_mid = ci_auc[2],
      AUC_95CI_high = ci_auc[3],
      n_samples = sum(keep),
      n_nonresponders = sum(auc_df$NonResponse[keep] == 1),
      n_responders = sum(auc_df$NonResponse[keep] == 0),
      AUC_CI_method = "stratified bootstrap",
      AUC_bootstrap_n = auc_bootstrap_n,
      AUC_seed = auc_bootstrap_seed,
      stringsAsFactors = FALSE
    )
  }
}

auc_table <- dplyr::bind_rows(auc_list)

auc_table <- auc_table %>%
  dplyr::mutate(
    State = factor(State, levels = bulk_state_order),
    StateLabel2 = factor(
      bulk_state_display[as.character(State)],
      levels = bulk_state_display[bulk_state_order]
    ),
    AUC_highlight = ifelse(
      as.character(State) == "Tumor_dedifferentiation_Stromal_remodeling",
      "Primary signal",
      "Other states"
    ),
    AUC_label = ifelse(
      is.finite(AUC_95CI_low) & is.finite(AUC_95CI_high),
      sprintf("%.3f\n95%% CI %.3f–%.3f",
              AUC_raw_score_for_nonresponse,
              AUC_95CI_low,
              AUC_95CI_high),
      sprintf("%.3f", AUC_raw_score_for_nonresponse)
    )
  )

## Updated table containing AUC and bootstrap 95% CI.
safe_write_csv(
  auc_table,
  file.path(out_table_dir, "GSE78220_nonresponse_AUC_with_95CI.csv"),
  row.names = FALSE
)

pB <- ggplot(
  auc_table,
  aes(x = StateLabel2, y = AUC_raw_score_for_nonresponse, fill = State)
) +
  geom_col(width = 0.65) +
  geom_errorbar(
    aes(ymin = AUC_95CI_low, ymax = AUC_95CI_high),
    width = 0.18,
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_text(
    aes(label = sprintf("%.3f", AUC_raw_score_for_nonresponse)),
    vjust = -0.6,
    size = 3.3,
    family = "sans"
  ) +
  scale_fill_manual(values = state_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(
    title = "GSE78220: state-score AUC for non-response",
    subtitle = "AUC < 0.5 indicates an inverse association with non-response",
    x = NULL,
    y = "AUC"
  ) +
  theme_icb(base_size = 12.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10.75),
    plot.subtitle = element_text(hjust = 0.5, size = 11.5)
  )

safe_ggsave(
  file.path(out_fig_dir, "Figure7B_GSE78220_nonresponse_AUC_with_95CI.pdf"),
  pB,
  width = 8.8,
  height = 5.5
)

safe_ggsave(
  file.path(out_fig_dir, "Figure7B_GSE78220_nonresponse_AUC_with_95CI.png"),
  pB,
  width = 8.8,
  height = 5.5,
  dpi = 300
)

## Panel C: Cox forest
## Move HR/P annotation farther to the right to avoid overlap with CI bars.
x_min <- min(cox_table$CI_lower, na.rm = TRUE) * 0.90
x_text <- max(cox_table$CI_upper, na.rm = TRUE) * 1.35
x_max <- x_text * 1.30

pC <- ggplot(cox_table, aes(y = StateLabelPlot, x = HR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_segment(
    aes(
      x = CI_lower,
      xend = CI_upper,
      y = StateLabelPlot,
      yend = StateLabelPlot,
      color = Significant
    ),
    linewidth = 0.7
  ) +
  geom_point(aes(color = Significant), size = 2.8) +
  geom_text(
    aes(x = x_text, label = LabelText),
    hjust = 0,
    lineheight = 0.95,
    size = 3.1,
    family = "sans"
  ) +
  scale_x_log10(
    limits = c(x_min, x_max),
    breaks = c(0.3, 0.5, 1, 3, 5),
    labels = c("0.3", "0.5", "1.0", "3.0", "5.0")
  ) +
  scale_color_manual(
    values = c("FALSE" = "#3C5488", "TRUE" = "#E64B35"),
    guide = "none"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "GSE78220 overall survival analysis (baseline samples only)",
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +
  theme_icb(base_size = 12.5) +
  theme(plot.margin = margin(8, 75, 8, 8))

safe_ggsave(
  file.path(out_fig_dir, "Figure7C_GSE78220_overall_survival_cox_baseline_only.pdf"),
  pC,
  width = 8.5,
  height = 4.8
)

safe_ggsave(
  file.path(out_fig_dir, "Figure7C_GSE78220_overall_survival_cox_baseline_only.png"),
  pC,
  width = 8.5,
  height = 4.8,
  dpi = 300
)

## Panel D: state correlation
cor_mat <- cor(
  merged78220[, bulk_state_cols, drop = FALSE],
  method = "spearman",
  use = "pairwise.complete.obs"
)

safe_write_csv(
  as.data.frame(cor_mat),
  file.path(out_table_dir, "GSE78220_state_spearman_correlation.csv"),
  row.names = TRUE
)

cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("State1", "State2", "Rho")

cor_long$State1Label <- factor(
  bulk_state_display[as.character(cor_long$State1)],
  levels = bulk_state_display[bulk_state_order]
)
cor_long$State2Label <- factor(
  bulk_state_display[as.character(cor_long$State2)],
  levels = rev(bulk_state_display[bulk_state_order])
)

pD <- ggplot(cor_long, aes(x = State1Label, y = State2Label, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Rho)), size = 3.1, family = "sans") +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  coord_equal() +
  labs(
    title = "GSE78220: state correlation",
    x = NULL,
    y = NULL
  ) +
  theme_icb(base_size = 12.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10.75),
    panel.grid = element_blank()
  )

safe_ggsave(
  file.path(out_fig_dir, "Figure7D_GSE78220_state_correlation.pdf"),
  pD,
  width = 6.8,
  height = 5.8
)

safe_ggsave(
  file.path(out_fig_dir, "Figure7D_GSE78220_state_correlation.png"),
  pD,
  width = 6.8,
  height = 5.8,
  dpi = 300
)

############################################################
## 11. Final manuscript Figure 7
############################################################

fig7 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    tag_levels = "A",
    title = "Descriptive association of predefined tumor–immune state scores with ICB response in GSE78220",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16.5, family = "sans"),
      plot.tag = element_text(face = "bold", size = 16, family = "sans")
    )
  )

safe_ggsave(
  file.path(out_fig_dir, "Figure7_GSE78220_ICB_response_association.pdf"),
  fig7,
  width = 13,
  height = 9.2
)

safe_ggsave(
  file.path(out_fig_dir, "Figure7_GSE78220_ICB_response_association.png"),
  fig7,
  width = 13,
  height = 9.2,
  dpi = 300
)

############################################################
## 12. Output inventory and final RDS
############################################################

safe_save_rds(
  list(
    final_state_scores_with_OS = merged78220,
    phenotype_OS_patient_level = pheno_os_patient,
    cox_input_baseline_only = cox_input,
    cox_table_baseline_only = cox_table,
    auc_table = auc_table,
    state_spearman_correlation = cor_mat
  ),
  file.path(data_processed_dir, "GSE78220_state_response_survival_analysis.rds")
)

canonical_outputs <- c(
  file.path(out_table_dir, "GSE78220_external_state_scores.csv"),
  file.path(out_table_dir, "GSE78220_Cox_input_patient_level_baseline_only.csv"),
  file.path(out_table_dir, "GSE78220_Cox_input_summary_baseline_only.csv"),
  file.path(out_table_dir, "GSE78220_Cox_excluded_nonbaseline_samples.csv"),
  file.path(out_table_dir, "GSE78220_overall_survival_cox_baseline_only.csv"),
  file.path(out_table_dir, "GSE78220_nonresponse_AUC_with_95CI.csv"),
  file.path(out_table_dir, "GSE78220_state_spearman_correlation.csv"),
  file.path(out_fig_dir, "Figure7A_GSE78220_state_scores_by_response.pdf"),
  file.path(out_fig_dir, "Figure7B_GSE78220_nonresponse_AUC_with_95CI.pdf"),
  file.path(out_fig_dir, "Figure7C_GSE78220_overall_survival_cox_baseline_only.pdf"),
  file.path(out_fig_dir, "Figure7D_GSE78220_state_correlation.pdf"),
  file.path(out_fig_dir, "Figure7_GSE78220_ICB_response_association.pdf"),
  file.path(data_processed_dir, "GSE78220_state_response_survival_analysis.rds")
)

output_inventory <- data.frame(
  File = basename(canonical_outputs),
  Path = canonical_outputs,
  Exists = file.exists(canonical_outputs),
  Size_MB = ifelse(file.exists(canonical_outputs),
                   round(file.info(canonical_outputs)$size / 1024^2, 4),
                   NA_real_),
  stringsAsFactors = FALSE
)

safe_write_csv(
  output_inventory,
  file.path(out_table_dir, "GSE78220_public_output_inventory.csv"),
  row.names = FALSE
)

save_session_info("02_GSE78220_state_scoring_response_survival_analysis")

message("02_GSE78220_state_scoring_response_survival_analysis.R finished successfully.")

