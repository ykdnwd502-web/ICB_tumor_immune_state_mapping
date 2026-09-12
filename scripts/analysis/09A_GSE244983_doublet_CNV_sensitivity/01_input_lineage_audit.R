############################################################
## 09A-01. GSE244983 26,053-cell doublet/CNV sensitivity
## Input-lineage audit
##
## Scientific boundary
## -------------------
## This module is a revision-stage technical sensitivity branch.
## It starts from the current 26,053-cell state-localized object.
## It does NOT reconstruct or target the historical 25,972-cell branch,
## and it does NOT replace the primary 26,053-cell analysis lineage.
############################################################

options(stringsAsFactors = FALSE)

SCRIPT_ID <- "09A_01_input_lineage_audit"

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  probe <- file.path(
    cwd, "results", "intermediate", "GSE244983",
    "GSE244983_seurat_state_localized.rds"
  )
  if (file.exists(probe)) {
    project_dir <- cwd
  } else {
    stop(
      "Set ICB_PROJECT_DIR to the repository root, or run from the project root.",
      call. = FALSE
    )
  }
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

required_pkgs <- c("Seurat", "SeuratObject", "readr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Restore the frozen environment and rerun.",
    call. = FALSE
  )
}

input_rds <- file.path(
  project_dir, "results", "intermediate", "GSE244983",
  "GSE244983_seurat_state_localized.rds"
)
out_table_dir <- file.path(
  project_dir, "results", "tables", "GSE244983",
  "doublet_CNV_sensitivity"
)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) {
  stop("Required 26,053-cell state-localized object not found: ", input_rds, call. = FALSE)
}

md5_file <- function(path) {
  unname(tools::md5sum(path))
}

cellset_md5 <- function(cells) {
  tf <- tempfile(fileext = ".txt")
  on.exit(unlink(tf), add = TRUE)
  writeLines(sort(as.character(cells)), tf, useBytes = TRUE)
  unname(tools::md5sum(tf))
}

write_csv <- function(x, name) {
  readr::write_csv(x, file.path(out_table_dir, name), na = "")
}

expected_n_cells <- 26053L
expected_n_features <- 11616L
expected_cellset_md5 <- "a574b6a23054be747b74e2e09f5ddea9"

expected_samples <- c(
  "Pat_ICBnaive1" = 7414L,
  "Pat_ICBnaive2" = 9081L,
  "Pat42" = 7188L,
  "Pat5" = 2370L
)

expected_major_counts <- c(
  "Malignant" = 17708L,
  "Cycling malignant" = 1888L,
  "CAF/stromal-like cells" = 333L,
  "Endothelial" = 1586L,
  "Myeloid cells" = 231L,
  "T/NK cells" = 2973L,
  "T/NK/Treg-like cells" = 382L,
  "B/Plasma cells" = 952L,
  "Cycling cells" = 0L
)

expected_cluster_counts <- c(
  "0" = 4547L, "1" = 3896L, "2" = 2705L, "3" = 2343L,
  "4" = 2205L, "5" = 2012L, "6" = 1888L, "7" = 1565L,
  "8" = 1541L, "9" = 1103L, "10" = 616L, "11" = 382L,
  "12" = 333L, "13" = 305L, "14" = 231L, "15" = 193L,
  "16" = 143L, "17" = 45L
)

state_cols <- c(
  "Immune_defective_Cold_score",
  "Myeloid_Treg_Immunosuppressive_score",
  "Tumor_dedifferentiation_Stromal_remodeling_score",
  "Melanocytic_Differentiation_score"
)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input RDS is not a Seurat object.", call. = FALSE)
}

required_meta <- c("Sample", "MajorCellType", "SeuratCluster", state_cols)
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0L) {
  stop(
    "Required metadata columns missing from current state-localized object: ",
    paste(missing_meta, collapse = ", "),
    call. = FALSE
  )
}

observed_sample <- table(factor(as.character(obj$Sample), levels = names(expected_samples)))
observed_major <- table(factor(as.character(obj$MajorCellType), levels = names(expected_major_counts)))
observed_cluster <- table(factor(as.character(obj$SeuratCluster), levels = names(expected_cluster_counts)))

lineage_gates <- data.frame(
  gate = c(
    "n_cells_26053",
    "n_features_11616",
    "cellset_md5_exact",
    "sample_counts_exact",
    "major_celltype_counts_exact",
    "cluster_counts_exact",
    "four_state_scores_present",
    "four_state_scores_finite"
  ),
  pass = c(
    ncol(obj) == expected_n_cells,
    nrow(obj) == expected_n_features,
    identical(cellset_md5(colnames(obj)), expected_cellset_md5),
    identical(as.integer(observed_sample), as.integer(expected_samples)),
    identical(as.integer(observed_major), as.integer(expected_major_counts)),
    identical(as.integer(observed_cluster), as.integer(expected_cluster_counts)),
    all(state_cols %in% colnames(obj@meta.data)),
    all(vapply(state_cols, function(x) all(is.finite(obj@meta.data[[x]])), logical(1)))
  ),
  stringsAsFactors = FALSE
)
write_csv(lineage_gates, "D01_input_lineage_gates.csv")

if (!all(lineage_gates$pass)) {
  stop(
    "09A input lineage gate failed. Inspect D01_input_lineage_gates.csv before continuing.",
    call. = FALSE
  )
}

input_audit <- data.frame(
  input_rds = input_rds,
  input_md5 = md5_file(input_rds),
  n_cells = ncol(obj),
  n_features = nrow(obj),
  cellset_md5 = cellset_md5(colnames(obj)),
  default_assay = Seurat::DefaultAssay(obj),
  n_samples = length(unique(as.character(obj$Sample))),
  n_major_celltypes_observed = length(unique(as.character(obj$MajorCellType[!is.na(obj$MajorCellType)]))),
  n_clusters = length(unique(as.character(obj$SeuratCluster))),
  state_score_columns = paste(state_cols, collapse = ";"),
  stringsAsFactors = FALSE
)
write_csv(input_audit, "D01_input_object_audit.csv")

write_csv(
  data.frame(Sample = names(expected_samples), N = as.integer(observed_sample)),
  "D01_sample_counts.csv"
)
write_csv(
  data.frame(MajorCellType = names(expected_major_counts), N = as.integer(observed_major)),
  "D01_major_celltype_counts.csv"
)
write_csv(
  data.frame(SeuratCluster = names(expected_cluster_counts), N = as.integer(observed_cluster)),
  "D01_cluster_counts.csv"
)

parameter_contract <- data.frame(
  parameter = c(
    "analysis_role",
    "primary_cell_universe",
    "historical_25972_targeted",
    "primary_pANN_expected_rate",
    "pANN_rate_sensitivity",
    "pN_artificial",
    "pK_fixed",
    "pANN_PCs",
    "pANN_variable_features",
    "pANN_max_neighbors",
    "inferCNV_cutoff",
    "inferCNV_HMM",
    "inferCNV_denoise",
    "inferCNV_reference_major_types"
  ),
  value = c(
    "technical_sensitivity_only",
    "26053_current_revised_cells",
    "FALSE",
    "0.075",
    "0.05;0.075;0.10",
    "0.25",
    "0.09",
    "1:20",
    "2000",
    "300",
    "0.1",
    "FALSE",
    "TRUE",
    "T/NK cells;T/NK/Treg-like cells;B/Plasma cells;Myeloid cells"
  ),
  stringsAsFactors = FALSE
)
write_csv(parameter_contract, "D01_method_parameter_contract.csv")

writeLines(capture.output(sessionInfo()), file.path(out_table_dir, paste0(SCRIPT_ID, "_sessionInfo.txt")))
message("09A-01 input-lineage audit PASS: 26,053 cells / 11,616 features.")
