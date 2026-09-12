############################################################
## 04_GSE91061_state_scoring_response_analysis.R
##
## Purpose:
##   Reproduce the supplementary external response-association assessment
##   in the GSE91061 anti-PD-1 melanoma cohort.
##
## Frozen analytical definitions:
##   1) Use the hg19KnownGene RLD expression matrix.
##   2) Map Entrez Gene IDs to HGNC gene symbols with org.Hs.eg.db.
##   3) Do NOT apply an additional log transformation to RLD values.
##   4) Score the same frozen tumor–immune gene sets used by scripts 01–03.
##   5) Preserve pretreatment and on-treatment samples as separate samples.
##   6) Primary response assessment: pretreatment CR/PR vs PD only.
##      SD and unknown/missing response are excluded.
##   7) Clinical-benefit sensitivity analysis: pretreatment CR/PR/SD vs PD.
##   8) Paired Pre-to-On change is exploratory; if multiple On samples exist
##      for one patient, On-treatment state scores are averaged by patient.
##
## Interpretation boundary:
##   GSE91061 is a supplementary external response-association assessment,
##   not a direct validation cohort and not a predictive model validation.
##
## Reproducibility:
##   Run from the repository root, or set ICB_PROJECT_DIR to the repository
##   root before running. Restore package versions from renv.lock.
############################################################

options(stringsAsFactors = FALSE)

############################################################
## 0. Project configuration
############################################################

resolve_project_dir <- function() {
  env_dir <- Sys.getenv("ICB_PROJECT_DIR")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, winslash = "/", mustWork = FALSE))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  wd_has_project_anchors <- (
    dir.exists(file.path(wd, "data_raw")) &&
    dir.exists(file.path(wd, "results"))
  ) || file.exists(file.path(wd, "renv.lock"))

  if (wd_has_project_anchors) {
    return(wd)
  }

  stop(
    "Cannot determine the repository/project root. ",
    "Either run R from the repository root or set:\n",
    'Sys.setenv(ICB_PROJECT_DIR = "D:/ICB_resistance_project")\n',
    "Current working directory is: ", wd
  )
}

project_dir <- resolve_project_dir()
message("Project directory: ", project_dir)

raw_dir          <- file.path(project_dir, "data_raw")
processed_dir    <- file.path(project_dir, "data_processed")
table_dir        <- file.path(project_dir, "results", "tables", "GSE91061")
intermediate_dir <- file.path(project_dir, "results", "intermediate", "GSE91061")
log_dir          <- file.path(project_dir, "logs")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "data.table",
  "dplyr",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Restore the repository environment with renv::restore() before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

############################################################
## 2. Helpers
############################################################

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
  message("Saved table: ", file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = file)
  message("Saved RDS: ", file)
}

save_session_info <- function() {
  log_file <- file.path(
    log_dir,
    "sessionInfo_04_GSE91061_state_scoring_response_analysis.txt"
  )
  sink(log_file)
  print(sessionInfo())
  sink()
  message("sessionInfo saved to: ", log_file)
}

path_relative_to_project <- function(path) {
  x <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  ifelse(startsWith(x, prefix), substring(x, nchar(prefix) + 1L), x)
}

safe_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  gs <- lapply(lines, function(x) {
    parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
    unique(parts[-c(1, 2)])
  })
  names(gs) <- vapply(
    lines,
    function(x) strsplit(x, "\t", fixed = TRUE)[[1]][1],
    character(1)
  )
  gs
}

row_z <- function(mat) {
  m <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, stats::sd, na.rm = TRUE)
  keep <- is.finite(s) & s > 0

  out <- mat[keep, , drop = FALSE]
  out <- sweep(out, 1, m[keep], "-")
  out <- sweep(out, 1, s[keep], "/")
  out
}

normalize_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", x))
}

zvec <- function(x) {
  out <- as.numeric(scale(as.numeric(x)))
  if (any(!is.finite(out))) {
    stop("Non-finite values generated during final state-score standardization.")
  }
  out
}

get_signature_cols <- function(df, target_names, state_name, min_required = 1L) {
  available <- setdiff(colnames(df), "Sample")
  available_norm <- normalize_name(available)
  names(available_norm) <- available
  target_norm <- normalize_name(target_names)

  matched <- available[available_norm %in% target_norm]

  if (length(matched) == 0) {
    for (tg in target_norm) {
      idx <- grep(tg, available_norm, fixed = TRUE)
      if (length(idx) > 0) {
        matched <- c(matched, available[idx])
      }
    }
  }

  matched <- unique(matched)

  message(
    state_name,
    " matched signature columns: ",
    paste(matched, collapse = ", ")
  )

  if (length(matched) < min_required) {
    stop("Insufficient matched signature columns for ", state_name)
  }

  matched
}

calc_auc_rank <- function(score, label01) {
  ok <- is.finite(score) & !is.na(label01)
  score <- score[ok]
  label01 <- label01[ok]

  if (length(unique(label01)) != 2) return(NA_real_)

  pos <- score[label01 == 1]
  neg <- score[label01 == 0]

  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)

  comp <- outer(pos, neg, "-")
  mean(comp > 0) + 0.5 * mean(comp == 0)
}

count_level <- function(x, level) {
  sum(x == level, na.rm = TRUE)
}

gate_numeric <- function(name, observed, expected, tolerance = 0) {
  pass <- is.finite(observed) &&
    is.finite(expected) &&
    abs(observed - expected) <= tolerance

  data.frame(
    Gate = name,
    Observed = format(observed, digits = 17, scientific = FALSE, trim = TRUE),
    Expected = format(expected, digits = 17, scientific = FALSE, trim = TRUE),
    Tolerance = format(tolerance, digits = 17, scientific = FALSE, trim = TRUE),
    Pass = pass,
    stringsAsFactors = FALSE
  )
}

gate_text <- function(name, observed, expected) {
  data.frame(
    Gate = name,
    Observed = as.character(observed),
    Expected = as.character(expected),
    Tolerance = "",
    Pass = identical(as.character(observed), as.character(expected)),
    stringsAsFactors = FALSE
  )
}

############################################################
## 3. Canonical inputs
############################################################

gse_dir <- file.path(
  raw_dir,
  "external_validation",
  "GSE91061"
)

rld_file <- file.path(
  gse_dir,
  "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz"
)

pdat_file <- file.path(
  gse_dir,
  "GSE91061_pData_raw_from_GEO.csv"
)

gmt_file <- file.path(
  project_dir,
  "results",
  "tables",
  "ICBcomb_final_input_gene_sets_CLEAN",
  "02_final_ICBcomb_gene_sets",
  "final_ICBcomb_gene_sets_CLEAN.gmt"
)

input_files <- c(
  RLD_expression = rld_file,
  GEO_phenotype = pdat_file,
  Frozen_gene_sets = gmt_file
)

input_audit <- data.frame(
  Input = names(input_files),
  RelativePath = vapply(input_files, path_relative_to_project, character(1)),
  Exists = file.exists(input_files),
  MD5 = vapply(input_files, safe_md5, character(1)),
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_audit,
  file.path(table_dir, "GSE91061_input_file_audit.csv")
)

if (any(!input_audit$Exists)) {
  missing_inputs <- input_audit$RelativePath[!input_audit$Exists]
  stop(
    "Missing canonical GSE91061 input(s):\n",
    paste(missing_inputs, collapse = "\n")
  )
}

############################################################
## 4. RLD expression matrix and Entrez-to-HGNC mapping
############################################################

message("Reading GSE91061 RLD matrix: ", rld_file)

expr_dt <- data.table::fread(rld_file)

if (!"V1" %in% colnames(expr_dt)) {
  stop("Expected first RLD column named V1 containing Entrez Gene IDs.")
}

entrez_ids <- as.character(expr_dt[["V1"]])

expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- entrez_ids

sample_names <- colnames(expr_mat)

message("Mapping Entrez IDs to HGNC symbols using org.Hs.eg.db...")

map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(entrez_ids),
  keytype = "ENTREZID",
  columns = "SYMBOL"
) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(ENTREZID, SYMBOL)

safe_write_csv(
  map_df,
  file.path(table_dir, "GSE91061_entrez_to_symbol_mapping.csv")
)

symbol_vec <- map_df$SYMBOL[match(entrez_ids, map_df$ENTREZID)]

keep_mapped <- !is.na(symbol_vec) & symbol_vec != ""
expr_mat <- expr_mat[keep_mapped, , drop = FALSE]
symbol_vec <- symbol_vec[keep_mapped]

## Preserve the frozen duplicate-symbol handling: mean expression by HGNC symbol.
expr_symbol_dt <- data.table::as.data.table(expr_mat)
expr_symbol_dt[, gene_symbol := symbol_vec]

expr_symbol_dt <- expr_symbol_dt[
  ,
  lapply(.SD, mean, na.rm = TRUE),
  by = gene_symbol
]

expr_for_scoring <- as.matrix(
  expr_symbol_dt[, -1, with = FALSE]
)

rownames(expr_for_scoring) <- expr_symbol_dt$gene_symbol
mode(expr_for_scoring) <- "numeric"

## RLD is already log-like: no further log transformation.
transform_used <- "RLD values, no additional log transform"

expression_audit <- data.frame(
  Source = "rld",
  TransformUsed = transform_used,
  n_input_rows = nrow(expr_dt),
  n_mapped_gene_symbols = nrow(expr_for_scoring),
  n_samples = ncol(expr_for_scoring),
  first_sample = sample_names[1],
  last_sample = sample_names[length(sample_names)],
  stringsAsFactors = FALSE
)

safe_write_csv(
  expression_audit,
  file.path(table_dir, "GSE91061_expression_matrix_audit.csv")
)

safe_save_rds(
  expr_for_scoring,
  file.path(
    intermediate_dir,
    "GSE91061_expression_symbol_matrix_for_state_scoring.rds"
  )
)

############################################################
## 5. Frozen gene sets and gene-wise z-mean scoring
############################################################

gene_sets <- read_gmt(gmt_file)

safe_write_csv(
  data.frame(
    GeneSet = names(gene_sets),
    n_genes = vapply(gene_sets, length, integer(1)),
    stringsAsFactors = FALSE
  ),
  file.path(table_dir, "GSE91061_gene_set_summary.csv")
)

zmat <- row_z(expr_for_scoring)

score_gene_set <- function(gs_name, genes) {
  present <- intersect(genes, rownames(zmat))

  gene_presence <- data.frame(
    GeneSet = gs_name,
    Gene = genes,
    Present = genes %in% rownames(zmat),
    stringsAsFactors = FALSE
  )

  if (length(present) < 3) {
    warning(
      "Too few genes present for ",
      gs_name,
      ": ",
      length(present),
      "/",
      length(genes)
    )
    score <- rep(NA_real_, ncol(zmat))
  } else {
    score <- colMeans(
      zmat[present, , drop = FALSE],
      na.rm = TRUE
    )
  }

  list(
    score = score,
    presence = gene_presence,
    n_present = length(present),
    n_total = length(genes)
  )
}

raw_scores <- data.frame(
  Sample = colnames(zmat),
  stringsAsFactors = FALSE
)

presence_all <- data.frame()
presence_summary <- data.frame()

for (nm in names(gene_sets)) {
  tmp <- score_gene_set(nm, gene_sets[[nm]])

  raw_scores[[nm]] <- tmp$score

  presence_all <- rbind(
    presence_all,
    tmp$presence
  )

  presence_summary <- rbind(
    presence_summary,
    data.frame(
      GeneSet = nm,
      n_total = tmp$n_total,
      n_present = tmp$n_present,
      stringsAsFactors = FALSE
    )
  )
}

presence_summary$coverage_fraction <-
  presence_summary$n_present / presence_summary$n_total

safe_write_csv(
  raw_scores,
  file.path(table_dir, "GSE91061_raw_gene_set_scores.csv")
)

safe_write_csv(
  presence_all,
  file.path(table_dir, "GSE91061_gene_presence_all_gene_sets.csv")
)

safe_write_csv(
  presence_summary,
  file.path(table_dir, "GSE91061_gene_presence_summary.csv")
)

############################################################
## 6. Four predefined tumor–immune state scores
############################################################

state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

display_labels_one_line <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" =
    "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

immune_cols <- get_signature_cols(
  raw_scores,
  c(
    "Immune_defective_Cold_RESTORE_hybrid",
    "Immune_defective_Cold_RESTORE_data_driven",
    "Immune_defective_Cold_RESTORE_curated"
  ),
  "Immune_defective_Cold"
)

myeloid_cols <- get_signature_cols(
  raw_scores,
  c(
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated"
  ),
  "Myeloid_Treg_Immunosuppressive"
)

tumor_cols <- get_signature_cols(
  raw_scores,
  c(
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated"
  ),
  "Tumor_dedifferentiation_Stromal_remodeling"
)

melanocytic_cols <- get_signature_cols(
  raw_scores,
  c(
    "Melanocytic_Differentiation_REFERENCE_hybrid",
    "Melanocytic_Differentiation_REFERENCE_data_driven",
    "Melanocytic_Differentiation_REFERENCE_curated"
  ),
  "Melanocytic_Differentiation"
)

state_mapping <- data.frame(
  FinalState = c(
    rep("Immune_defective_Cold", length(immune_cols)),
    rep("Myeloid_Treg_Immunosuppressive", length(myeloid_cols)),
    rep(
      "Tumor_dedifferentiation_Stromal_remodeling",
      length(tumor_cols)
    ),
    rep("Melanocytic_Differentiation", length(melanocytic_cols))
  ),
  RawSignatureColumn = c(
    immune_cols,
    myeloid_cols,
    tumor_cols,
    melanocytic_cols
  ),
  DirectionForFinalState = c(
    rep("inverse", length(immune_cols)),
    rep("positive", length(myeloid_cols)),
    rep("positive", length(tumor_cols)),
    rep("positive", length(melanocytic_cols))
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  state_mapping,
  file.path(table_dir, "GSE91061_final_state_signature_mapping.csv")
)

state_scores <- data.frame(
  Sample = raw_scores$Sample,
  Immune_defective_Cold = zvec(
    -rowMeans(
      raw_scores[, immune_cols, drop = FALSE],
      na.rm = TRUE
    )
  ),
  Myeloid_Treg_Immunosuppressive = zvec(
    rowMeans(
      raw_scores[, myeloid_cols, drop = FALSE],
      na.rm = TRUE
    )
  ),
  Tumor_dedifferentiation_Stromal_remodeling = zvec(
    rowMeans(
      raw_scores[, tumor_cols, drop = FALSE],
      na.rm = TRUE
    )
  ),
  Melanocytic_Differentiation = zvec(
    rowMeans(
      raw_scores[, melanocytic_cols, drop = FALSE],
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

state_scores$Dominant_state <- apply(
  state_scores[, state_cols, drop = FALSE],
  1,
  function(x) state_cols[which.max(x)]
)

if (anyNA(state_scores[, state_cols, drop = FALSE])) {
  stop("Missing values detected in final GSE91061 state scores.")
}

############################################################
## 7. Clinical metadata and frozen response definitions
############################################################

pd <- read.csv(
  pdat_file,
  check.names = FALSE,
  row.names = 1,
  stringsAsFactors = FALSE
)

if (!"title" %in% colnames(pd)) {
  stop("GSE91061 phenotype metadata must contain column 'title'.")
}

response_col <- grep(
  "^response",
  colnames(pd),
  value = TRUE,
  ignore.case = TRUE
)[1]

visit_col <- grep(
  "visit.*pre.*on|visit",
  colnames(pd),
  value = TRUE,
  ignore.case = TRUE
)[1]

tissue_col <- grep(
  "^tissue",
  colnames(pd),
  value = TRUE,
  ignore.case = TRUE
)[1]

if (is.na(response_col) || is.na(visit_col)) {
  stop(
    "Could not identify the frozen response/visit metadata columns in GSE91061 pData."
  )
}

metadata_column_audit <- data.frame(
  Field = c("sample_title", "response", "visit", "tissue"),
  Column = c(
    "title",
    response_col,
    visit_col,
    ifelse(is.na(tissue_col), "", tissue_col)
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  metadata_column_audit,
  file.path(table_dir, "GSE91061_metadata_column_audit.csv")
)

meta <- data.frame(
  Sample = as.character(pd$title),
  GEO_accession = if ("geo_accession" %in% colnames(pd)) {
    as.character(pd$geo_accession)
  } else {
    rownames(pd)
  },
  ResponseRaw = as.character(pd[[response_col]]),
  VisitRaw = as.character(pd[[visit_col]]),
  Tissue = if (!is.na(tissue_col)) {
    as.character(pd[[tissue_col]])
  } else {
    NA_character_
  },
  stringsAsFactors = FALSE
)

meta$PatientID <- sub(
  "^(Pt[^_]+)_.*$",
  "\\1",
  meta$Sample
)

meta$VisitGroup <- NA_character_

meta$VisitGroup[
  grepl("Pre", meta$VisitRaw, ignore.case = TRUE) |
    grepl("_Pre_", meta$Sample)
] <- "Pre"

meta$VisitGroup[
  grepl("On", meta$VisitRaw, ignore.case = TRUE) |
    grepl("_On_", meta$Sample)
] <- "On"

meta$ResponseGroup <- NA_character_

meta$ResponseGroup[
  grepl(
    "CR|Complete|PR|Partial",
    meta$ResponseRaw,
    ignore.case = TRUE
  )
] <- "Responder"

meta$ResponseGroup[
  grepl(
    "PD|Progressive",
    meta$ResponseRaw,
    ignore.case = TRUE
  )
] <- "NonResponder"

meta$ResponseGroup[
  grepl(
    "SD|Stable",
    meta$ResponseRaw,
    ignore.case = TRUE
  )
] <- "StableDisease"

meta$ResponseGroup[
  grepl(
    "UNK|Unknown",
    meta$ResponseRaw,
    ignore.case = TRUE
  )
] <- NA_character_

meta$BinaryResponseGroup <- NA_character_
meta$BinaryResponseGroup[
  meta$ResponseGroup == "Responder"
] <- "Responder"
meta$BinaryResponseGroup[
  meta$ResponseGroup == "NonResponder"
] <- "NonResponder"

meta$ClinicalBenefitGroup <- NA_character_
meta$ClinicalBenefitGroup[
  meta$ResponseGroup %in% c("Responder", "StableDisease")
] <- "ClinicalBenefit"
meta$ClinicalBenefitGroup[
  meta$ResponseGroup == "NonResponder"
] <- "NoClinicalBenefit"

safe_write_csv(
  meta,
  file.path(table_dir, "GSE91061_clinical_metadata_formatted.csv")
)

matching_audit <- data.frame(
  n_expression_samples = nrow(state_scores),
  n_metadata_rows = nrow(meta),
  n_matched_by_title = length(
    intersect(state_scores$Sample, meta$Sample)
  ),
  n_unmatched_expression_samples = sum(
    !state_scores$Sample %in% meta$Sample
  ),
  n_unmatched_metadata_samples = sum(
    !meta$Sample %in% state_scores$Sample
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  matching_audit,
  file.path(table_dir, "GSE91061_metadata_matching_audit.csv")
)

merged_scores <- state_scores %>%
  dplyr::left_join(meta, by = "Sample")

pre_scores <- merged_scores %>%
  dplyr::filter(VisitGroup == "Pre")

on_scores <- merged_scores %>%
  dplyr::filter(VisitGroup == "On")

safe_write_csv(
  merged_scores,
  file.path(table_dir, "GSE91061_state_scores_all_samples.csv")
)

safe_write_csv(
  pre_scores,
  file.path(table_dir, "GSE91061_state_scores_pretreatment.csv")
)

safe_write_csv(
  on_scores,
  file.path(table_dir, "GSE91061_state_scores_on_treatment.csv")
)

safe_save_rds(
  merged_scores,
  file.path(
    intermediate_dir,
    "GSE91061_state_scores_all_samples.rds"
  )
)

visit_counts <- as.data.frame(
  table(merged_scores$VisitGroup, useNA = "ifany"),
  stringsAsFactors = FALSE
)
colnames(visit_counts) <- c("VisitGroup", "n_samples")

response_counts <- as.data.frame(
  table(pre_scores$ResponseGroup, useNA = "ifany"),
  stringsAsFactors = FALSE
)
colnames(response_counts) <- c("ResponseGroup", "n_samples")

binary_counts <- as.data.frame(
  table(pre_scores$BinaryResponseGroup, useNA = "ifany"),
  stringsAsFactors = FALSE
)
colnames(binary_counts) <- c("BinaryResponseGroup", "n_samples")

benefit_counts <- as.data.frame(
  table(pre_scores$ClinicalBenefitGroup, useNA = "ifany"),
  stringsAsFactors = FALSE
)
colnames(benefit_counts) <- c("ClinicalBenefitGroup", "n_samples")

safe_write_csv(
  visit_counts,
  file.path(table_dir, "GSE91061_visit_group_counts.csv")
)

safe_write_csv(
  response_counts,
  file.path(table_dir, "GSE91061_pretreatment_response_group_counts.csv")
)

safe_write_csv(
  binary_counts,
  file.path(
    table_dir,
    "GSE91061_pretreatment_binary_response_group_counts.csv"
  )
)

safe_write_csv(
  benefit_counts,
  file.path(
    table_dir,
    "GSE91061_pretreatment_clinical_benefit_group_counts.csv"
  )
)

############################################################
## 8. Primary pretreatment response association
##    Strict binary comparison: CR/PR vs PD
############################################################

binary_df <- pre_scores %>%
  dplyr::filter(
    BinaryResponseGroup %in% c("Responder", "NonResponder")
  )

binary_label <- ifelse(
  binary_df$BinaryResponseGroup == "NonResponder",
  1,
  0
)

auc_binary <- lapply(
  state_cols,
  function(st) {
    data.frame(
      Analysis = "CRPR_vs_PD_binary_nonresponse",
      State = st,
      StateLabel = display_labels_one_line[st],
      n_nonresponder = sum(binary_label == 1, na.rm = TRUE),
      n_responder = sum(binary_label == 0, na.rm = TRUE),
      AUC_for_nonresponse = calc_auc_rank(
        binary_df[[st]],
        binary_label
      ),
      stringsAsFactors = FALSE
    )
  }
) %>%
  dplyr::bind_rows()

safe_write_csv(
  auc_binary,
  file.path(
    table_dir,
    "GSE91061_pretreatment_nonresponse_AUC_CRPR_vs_PD.csv"
  )
)

############################################################
## 9. Clinical-benefit sensitivity analysis
##    CR/PR/SD vs PD
############################################################

benefit_df <- pre_scores %>%
  dplyr::filter(
    ClinicalBenefitGroup %in% c(
      "ClinicalBenefit",
      "NoClinicalBenefit"
    )
  )

benefit_label <- ifelse(
  benefit_df$ClinicalBenefitGroup == "NoClinicalBenefit",
  1,
  0
)

auc_benefit <- lapply(
  state_cols,
  function(st) {
    data.frame(
      Analysis = "ClinicalBenefit_vs_PD_nonbenefit",
      State = st,
      StateLabel = display_labels_one_line[st],
      n_no_benefit = sum(benefit_label == 1, na.rm = TRUE),
      n_benefit = sum(benefit_label == 0, na.rm = TRUE),
      AUC_for_nonbenefit = calc_auc_rank(
        benefit_df[[st]],
        benefit_label
      ),
      stringsAsFactors = FALSE
    )
  }
) %>%
  dplyr::bind_rows()

safe_write_csv(
  auc_benefit,
  file.path(
    table_dir,
    "GSE91061_pretreatment_nonbenefit_AUC_clinical_benefit.csv"
  )
)

############################################################
## 10. Pretreatment state correlation
############################################################

cor_mat <- stats::cor(
  pre_scores[, state_cols, drop = FALSE],
  method = "spearman",
  use = "pairwise.complete.obs"
)

cor_out <- data.frame(
  State = rownames(cor_mat),
  cor_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

safe_write_csv(
  cor_out,
  file.path(
    table_dir,
    "GSE91061_pretreatment_state_spearman_correlation.csv"
  )
)

############################################################
## 11. Exploratory paired Pre-to-On state-score changes
############################################################

paired_ids <- intersect(
  pre_scores$PatientID,
  on_scores$PatientID
)

paired_pre <- pre_scores %>%
  dplyr::filter(PatientID %in% paired_ids) %>%
  dplyr::select(
    PatientID,
    ResponseRaw,
    ResponseGroup,
    BinaryResponseGroup,
    ClinicalBenefitGroup,
    dplyr::all_of(state_cols)
  ) %>%
  dplyr::group_by(PatientID) %>%
  dplyr::summarise(
    ResponseRaw = dplyr::first(ResponseRaw),
    ResponseGroup = dplyr::first(ResponseGroup),
    BinaryResponseGroup = dplyr::first(BinaryResponseGroup),
    ClinicalBenefitGroup = dplyr::first(ClinicalBenefitGroup),
    dplyr::across(
      dplyr::all_of(state_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

paired_on <- on_scores %>%
  dplyr::filter(PatientID %in% paired_ids) %>%
  dplyr::select(
    PatientID,
    dplyr::all_of(state_cols)
  ) %>%
  dplyr::group_by(PatientID) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(state_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

paired_change <- paired_pre %>%
  dplyr::inner_join(
    paired_on,
    by = "PatientID",
    suffix = c("_Pre", "_On")
  )

for (st in state_cols) {
  paired_change[[paste0(st, "_Delta_On_minus_Pre")]] <-
    paired_change[[paste0(st, "_On")]] -
    paired_change[[paste0(st, "_Pre")]]
}

safe_write_csv(
  paired_change,
  file.path(
    table_dir,
    "GSE91061_paired_pre_on_state_change.csv"
  )
)

safe_save_rds(
  paired_change,
  file.path(
    intermediate_dir,
    "GSE91061_paired_pre_on_state_change.rds"
  )
)

############################################################
## 12. Frozen public-release reproducibility gates
##
## These benchmark values were established from the validated mother
## workflow. They are audit checks only and are not used to calculate
## or select any analytical result.
############################################################

expected_primary_auc <- c(
  "Immune_defective_Cold" = 0.660869565217391,
  "Myeloid_Treg_Immunosuppressive" = 0.365217391304348,
  "Tumor_dedifferentiation_Stromal_remodeling" = 0.434782608695652,
  "Melanocytic_Differentiation" = 0.621739130434783
)

expected_benefit_auc <- c(
  "Immune_defective_Cold" = 0.586956521739130,
  "Myeloid_Treg_Immunosuppressive" = 0.446488294314381,
  "Tumor_dedifferentiation_Stromal_remodeling" = 0.471571906354515,
  "Melanocytic_Differentiation" = 0.525083612040134
)

gate_list <- list(
  gate_text(
    "expression_source",
    expression_audit$Source[1],
    "rld"
  ),
  gate_text(
    "expression_transform",
    expression_audit$TransformUsed[1],
    "RLD values, no additional log transform"
  ),
  gate_numeric(
    "expression_input_rows",
    nrow(expr_dt),
    22187
  ),
  gate_numeric(
    "mapped_gene_symbols",
    nrow(expr_for_scoring),
    22068
  ),
  gate_numeric(
    "expression_samples",
    nrow(merged_scores),
    109
  ),
  gate_numeric(
    "unique_patients_all",
    dplyr::n_distinct(merged_scores$PatientID),
    65
  ),
  gate_numeric(
    "metadata_rows",
    nrow(meta),
    109
  ),
  gate_numeric(
    "metadata_matched_by_title",
    matching_audit$n_matched_by_title[1],
    109
  ),
  gate_numeric(
    "metadata_unmatched_expression",
    matching_audit$n_unmatched_expression_samples[1],
    0
  ),
  gate_numeric(
    "metadata_unmatched_metadata",
    matching_audit$n_unmatched_metadata_samples[1],
    0
  ),
  gate_numeric(
    "pretreatment_samples",
    nrow(pre_scores),
    51
  ),
  gate_numeric(
    "on_treatment_samples",
    nrow(on_scores),
    58
  ),
  gate_numeric(
    "pretreatment_unique_patients",
    dplyr::n_distinct(pre_scores$PatientID),
    51
  ),
  gate_numeric(
    "on_treatment_unique_patients",
    dplyr::n_distinct(on_scores$PatientID),
    57
  ),
  gate_numeric(
    "pretreatment_responder",
    count_level(pre_scores$ResponseGroup, "Responder"),
    10
  ),
  gate_numeric(
    "pretreatment_nonresponder",
    count_level(pre_scores$ResponseGroup, "NonResponder"),
    23
  ),
  gate_numeric(
    "pretreatment_stable_disease",
    count_level(pre_scores$ResponseGroup, "StableDisease"),
    16
  ),
  gate_numeric(
    "pretreatment_response_missing",
    sum(is.na(pre_scores$ResponseGroup)),
    2
  ),
  gate_numeric(
    "strict_binary_responder",
    count_level(pre_scores$BinaryResponseGroup, "Responder"),
    10
  ),
  gate_numeric(
    "strict_binary_nonresponder",
    count_level(pre_scores$BinaryResponseGroup, "NonResponder"),
    23
  ),
  gate_numeric(
    "strict_binary_missing_or_excluded",
    sum(is.na(pre_scores$BinaryResponseGroup)),
    18
  ),
  gate_numeric(
    "clinical_benefit",
    count_level(
      pre_scores$ClinicalBenefitGroup,
      "ClinicalBenefit"
    ),
    26
  ),
  gate_numeric(
    "no_clinical_benefit",
    count_level(
      pre_scores$ClinicalBenefitGroup,
      "NoClinicalBenefit"
    ),
    23
  ),
  gate_numeric(
    "clinical_benefit_missing",
    sum(is.na(pre_scores$ClinicalBenefitGroup)),
    2
  ),
  gate_numeric(
    "paired_pre_on_patients",
    nrow(paired_change),
    43
  ),
  gate_numeric(
    "final_state_score_NA_count",
    sum(is.na(state_scores[, state_cols, drop = FALSE])),
    0
  )
)

for (st in state_cols) {
  observed_auc <- auc_binary$AUC_for_nonresponse[
    match(st, auc_binary$State)
  ]

  gate_list[[length(gate_list) + 1L]] <- gate_numeric(
    paste0("primary_AUC_", st),
    observed_auc,
    expected_primary_auc[[st]],
    tolerance = 1e-12
  )
}

for (st in state_cols) {
  observed_auc <- auc_benefit$AUC_for_nonbenefit[
    match(st, auc_benefit$State)
  ]

  gate_list[[length(gate_list) + 1L]] <- gate_numeric(
    paste0("clinical_benefit_AUC_", st),
    observed_auc,
    expected_benefit_auc[[st]],
    tolerance = 1e-12
  )
}

reproducibility_gates <- dplyr::bind_rows(gate_list)

safe_write_csv(
  reproducibility_gates,
  file.path(
    table_dir,
    "GSE91061_reproducibility_gates.csv"
  )
)

if (any(!reproducibility_gates$Pass)) {
  failed <- reproducibility_gates[
    !reproducibility_gates$Pass,
    ,
    drop = FALSE
  ]

  stop(
    "GSE91061 public-release reproducibility gate failed:\n",
    paste(
      paste0(
        failed$Gate,
        " [observed=",
        failed$Observed,
        "; expected=",
        failed$Expected,
        "]"
      ),
      collapse = "\n"
    )
  )
}

############################################################
## 13. Canonical output inventory and session information
############################################################

save_session_info()

canonical_outputs <- c(
  file.path(table_dir, "GSE91061_input_file_audit.csv"),
  file.path(table_dir, "GSE91061_expression_matrix_audit.csv"),
  file.path(table_dir, "GSE91061_entrez_to_symbol_mapping.csv"),
  file.path(table_dir, "GSE91061_gene_set_summary.csv"),
  file.path(table_dir, "GSE91061_gene_presence_summary.csv"),
  file.path(table_dir, "GSE91061_raw_gene_set_scores.csv"),
  file.path(table_dir, "GSE91061_final_state_signature_mapping.csv"),
  file.path(table_dir, "GSE91061_clinical_metadata_formatted.csv"),
  file.path(table_dir, "GSE91061_metadata_matching_audit.csv"),
  file.path(table_dir, "GSE91061_state_scores_all_samples.csv"),
  file.path(table_dir, "GSE91061_state_scores_pretreatment.csv"),
  file.path(table_dir, "GSE91061_state_scores_on_treatment.csv"),
  file.path(table_dir, "GSE91061_pretreatment_nonresponse_AUC_CRPR_vs_PD.csv"),
  file.path(table_dir, "GSE91061_pretreatment_nonbenefit_AUC_clinical_benefit.csv"),
  file.path(table_dir, "GSE91061_pretreatment_state_spearman_correlation.csv"),
  file.path(table_dir, "GSE91061_paired_pre_on_state_change.csv"),
  file.path(table_dir, "GSE91061_reproducibility_gates.csv"),
  file.path(
    intermediate_dir,
    "GSE91061_expression_symbol_matrix_for_state_scoring.rds"
  ),
  file.path(
    intermediate_dir,
    "GSE91061_state_scores_all_samples.rds"
  ),
  file.path(
    intermediate_dir,
    "GSE91061_paired_pre_on_state_change.rds"
  ),
  file.path(
    log_dir,
    "sessionInfo_04_GSE91061_state_scoring_response_analysis.txt"
  )
)

output_inventory <- data.frame(
  File = basename(canonical_outputs),
  RelativePath = vapply(
    canonical_outputs,
    path_relative_to_project,
    character(1)
  ),
  Exists = file.exists(canonical_outputs),
  SizeBytes = ifelse(
    file.exists(canonical_outputs),
    file.info(canonical_outputs)$size,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  output_inventory,
  file.path(table_dir, "GSE91061_output_inventory.csv")
)

message("04_GSE91061_state_scoring_response_analysis.R finished successfully.")
message("Primary analysis: pretreatment CR/PR vs PD (SD/unknown excluded).")
message("Clinical-benefit sensitivity: pretreatment CR/PR/SD vs PD.")
message("Exploratory paired analysis: On-treatment minus pretreatment.")
message("All frozen reproducibility gates passed.")
