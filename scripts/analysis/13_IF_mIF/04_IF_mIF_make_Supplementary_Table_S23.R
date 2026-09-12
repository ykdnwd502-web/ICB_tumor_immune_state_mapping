############################################################
## 04_IF_mIF_make_Supplementary_Table_S23_v1.2_PURE_XLSX.R
##
## FINAL PUBLIC RELEASE PATCH
## Supplementary Table S23 pure XLSX export.
##
## Scope:
## - Export layer only.
## - No analytical definition changed.
## - No input data changed.
## - No statistical calculation changed.
##
## Strategy:
## Use writexl::write_xlsx() to create a minimal,
## clean workbook without drawing/VML/comment relationships.
############################################################

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR")

if (!nzchar(PROJECT_DIR)) {
  PROJECT_DIR <- "D:/ICB_resistance_project"
}

PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)

if (!requireNamespace("writexl", quietly = TRUE)) {
  stop(
    "Package 'writexl' is required.",
    call. = FALSE
  )
}

TABLE_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "IF_mIF_validation"
)

FIG_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "figures",
  "IF_mIF_validation"
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

corr_file <- file.path(
  TABLE_DIR,
  "IF_mIF_protein_state_spearman.csv"
)

harm_file <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S23_IF_mIF_marker_signal_state_correlation_source.csv"
)

if (!file.exists(corr_file)) {
  stop(
    "Missing canonical Spearman source table.",
    call. = FALSE
  )
}

if (!file.exists(harm_file)) {
  stop(
    "Missing canonical S23 harmonized source table.",
    call. = FALSE
  )
}


corr <- read.csv(
  corr_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

harm <- read.csv(
  harm_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


############################################################
## HARD VALIDATION BEFORE EXPORT
############################################################

if (nrow(corr) != 32L) {
  stop(
    paste0(
      "IF_mIF_Spearman rows incorrect: ",
      nrow(corr)
    ),
    call. = FALSE
  )
}

if (nrow(harm) != 32L) {
  stop(
    paste0(
      "S23_harmonized rows incorrect: ",
      nrow(harm)
    ),
    call. = FALSE
  )
}


readme <- data.frame(
  Item = c(
    "Table",
    "Purpose",
    "Analysis unit",
    "Association",
    "Multiple testing",
    "Matched cells",
    "Interpretation boundary",
    "Source"
  ),
  Value = c(
    "Supplementary Table S23. IF/mIF protein–state Spearman correlations",
    "Source-locked protein/marker-signal association with the four RNA-defined tumor–immune state scores.",
    "Matched CosMx cell after FOV–cell matching",
    "Spearman correlation; cor.test(exact = FALSE)",
    "Benjamini–Hochberg FDR across all 32 protein-feature x state tests",
    "86,572",
    "Association support only; not functional protein validation, causal mechanism validation, or predictive biomarker validation.",
    "Public sequential rebuild from raw CosMx expression through matched-cell IF/mIF analysis"
  ),
  stringsAsFactors = FALSE
)


source_manifest <- data.frame(
  role = c(
    "CANONICAL_SPEARMAN_TABLE",
    "CANONICAL_S23_HARMONIZED_SOURCE"
  ),
  path = c(
    normalizePath(
      corr_file,
      winslash = "/",
      mustWork = TRUE
    ),
    normalizePath(
      harm_file,
      winslash = "/",
      mustWork = TRUE
    )
  ),
  stringsAsFactors = FALSE
)


############################################################
## PURE XLSX WRITE
############################################################

out_file <- file.path(
  TABLE_DIR,
  "Supplementary_Table_S23_IF_mIF.xlsx"
)

if (file.exists(out_file)) {
  unlink(out_file)
}


writexl::write_xlsx(
  x = list(
    README = readme,
    IF_mIF_Spearman = corr,
    S23_harmonized = harm,
    Source_manifest = source_manifest
  ),
  path = out_file
)


############################################################
## POST EXPORT CHECK
############################################################

if (!file.exists(out_file)) {
  stop(
    "Output xlsx was not created.",
    call. = FALSE
  )
}


cat("\n============================================================\n")
cat("IF/mIF SUPPLEMENTARY TABLE S23 PURE XLSX EXPORT COMPLETED\n")
cat("Output: ", out_file, "\n", sep = "")
cat("Sheets:\n")
cat(" - README\n")
cat(" - IF_mIF_Spearman\n")
cat(" - S23_harmonized\n")
cat(" - Source_manifest\n")
cat("Export engine: writexl::write_xlsx\n")
cat("Mode: PURE_NO_DRAWING_NO_VML\n")
cat("============================================================\n")
