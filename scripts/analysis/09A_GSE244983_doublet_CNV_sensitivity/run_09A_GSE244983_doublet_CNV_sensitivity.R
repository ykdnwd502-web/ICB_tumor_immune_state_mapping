############################################################
## Run 09A GSE244983 doublet/CNV sensitivity module
##
## Intended repository location:
## scripts/analysis/09A_GSE244983_doublet_CNV_sensitivity/
##
## Run after analysis 09 and before S20/CosMx/Visium modules.
############################################################



options(stringsAsFactors = FALSE)

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
Sys.setenv(ICB_PROJECT_DIR = project_dir)

module_dir <- file.path(
  project_dir, "scripts", "analysis",
  "09A_GSE244983_doublet_CNV_sensitivity"
)
if (!dir.exists(module_dir)) {
  stop("09A module directory not found at expected path: ", module_dir, call. = FALSE)
}

steps <- c(
  "01_input_lineage_audit.R",
  "02_pANN_doublet_detection.R",
  "03_finalize_doublet_filter.R",
  "04_primary_doublet_robustness.R",
  "05_doublet_rate_sensitivity.R",
  "06_generate_inferCNV_gene_order.R",
  "07_run_inferCNV_support.R",
  "08_final_module_audit.R"
)
step_paths <- file.path(module_dir, steps)
if (!all(file.exists(step_paths))) {
  stop("Missing 09A script(s):\n", paste(step_paths[!file.exists(step_paths)], collapse = "\n"), call. = FALSE)
}

## Full preflight before any expensive pANN/inferCNV work.
required_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "RANN", "readr", "dplyr", "tidyr",
  "AnnotationDbi", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "GenomicFeatures", "GenomeInfoDb", "BiocGenerics",
  "infercnv", "matrixStats", "digest"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop(
    "09A preflight failed. Missing package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Restore/install these before starting the long full clean run.",
    call. = FALSE
  )
}

input_probe <- file.path(
  project_dir, "results", "intermediate", "GSE244983",
  "GSE244983_seurat_state_localized.rds"
)
if (!file.exists(input_probe)) {
  stop(
    "09A preflight failed: current 26,053-cell state-localized object is missing. Run analysis 06/07 first.",
    call. = FALSE
  )
}

message("09A preflight PASS. Starting 8-step doublet/CNV sensitivity module.")
for (i in seq_along(step_paths)) {
  message("\n============================================================")
  message(sprintf("09A step %d/%d: %s", i, length(step_paths), basename(step_paths[[i]])))
  message("============================================================")
  source(step_paths[[i]], local = new.env(parent = globalenv()), chdir = FALSE)
}
message("\n09A GSE244983 doublet/CNV sensitivity module completed successfully.")
