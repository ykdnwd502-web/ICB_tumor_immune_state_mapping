############################################################
## 09A-03. Finalize primary pANN doublet sensitivity objects
##
## Primary pANN call = homotypic-adjusted 7.5% expected rate.
## v1.0.1 patch: use Seurat object column subsetting rather than
## Seurat::subset(), because subset is a base generic and is not exported
## from the Seurat namespace. No scientific logic or filtering rule changed.
## The singlet object is a strict subset of the current 26,053-cell
## state-localized object. It is NOT renormalized, reclustered, or
## substituted into the primary biological analysis.
############################################################

options(stringsAsFactors = FALSE)
SCRIPT_ID <- "09A_03_finalize_doublet_filter"

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

required_pkgs <- c("Seurat", "SeuratObject", "readr", "dplyr", "tidyr")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)

out_table_dir <- file.path(project_dir, "results", "tables", "GSE244983", "doublet_CNV_sensitivity")
out_obj_dir <- file.path(project_dir, "results", "objects", "GSE244983", "doublet_CNV_sensitivity")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_obj_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(out_obj_dir, "GSE244983_all_cells_with_pANN_doublet_calls.rds")
singlet_rds <- file.path(out_obj_dir, "GSE244983_singlets_after_pANN_doublet_filter.rds")

write_csv <- function(x, name) readr::write_csv(x, file.path(out_table_dir, name), na = "")

if (!file.exists(input_rds)) stop("09A-02 labeled object not found: ", input_rds, call. = FALSE)
obj <- readRDS(input_rds)
if (ncol(obj) != 26053L) stop("Expected 26,053 labeled cells.", call. = FALSE)

required_meta <- c(
  "Sample", "MajorCellType", "SeuratCluster", "pANN_manual",
  "pANN_doublet_5pct", "pANN_doublet_primary", "pANN_doublet_10pct"
)
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0L) stop("Missing required pANN metadata: ", paste(missing_meta, collapse = ", "), call. = FALSE)

primary_label <- as.character(obj$pANN_doublet_primary)
if (!all(primary_label %in% c("Singlet", "Doublet"))) stop("Invalid primary doublet labels.", call. = FALSE)

rate_audit_file <- file.path(out_table_dir, "D02_pANN_per_sample_rate_audit.csv")
if (!file.exists(rate_audit_file)) stop("Missing D02 parameter audit.", call. = FALSE)
rate_audit <- readr::read_csv(rate_audit_file, show_col_types = FALSE)
expected_primary_n <- sum(rate_audit$nExp_adjusted[abs(rate_audit$expected_rate - 0.075) < 1e-12])
observed_primary_n <- sum(primary_label == "Doublet")
if (!identical(as.integer(observed_primary_n), as.integer(expected_primary_n))) {
  stop("Primary predicted-doublet count does not match D02 homotypic-adjusted expectation.", call. = FALSE)
}

singlet_cells <- colnames(obj)[primary_label == "Singlet"]
obj_singlet <- obj[, singlet_cells]
if (!identical(colnames(obj_singlet), singlet_cells)) {
  stop("Singlet subset cell order/content mismatch.", call. = FALSE)
}
if (!identical(rownames(obj_singlet), rownames(obj))) {
  stop("Singlet subset unexpectedly changed the feature universe.", call. = FALSE)
}
if (ncol(obj_singlet) + observed_primary_n != ncol(obj)) stop("Singlet/doublet cardinality mismatch.", call. = FALSE)
saveRDS(obj_singlet, singlet_rds)

meta <- obj@meta.data
meta$Cell <- rownames(meta)

write_csv(
  data.frame(
    doublet_status = c("Singlet", "Doublet"),
    N = c(sum(primary_label == "Singlet"), sum(primary_label == "Doublet")),
    Fraction = c(mean(primary_label == "Singlet"), mean(primary_label == "Doublet")),
    stringsAsFactors = FALSE
  ),
  "D03_primary_doublet_global_summary.csv"
)

by_sample <- as.data.frame(table(
  Sample = as.character(meta$Sample),
  doublet_status = primary_label
), stringsAsFactors = FALSE)
by_sample <- by_sample[by_sample$Freq > 0, , drop = FALSE]
by_sample <- dplyr::group_by(by_sample, Sample)
by_sample <- dplyr::mutate(by_sample, Fraction = Freq / sum(Freq))
by_sample <- dplyr::ungroup(by_sample)
write_csv(by_sample, "D03_primary_doublet_by_sample.csv")

by_major <- as.data.frame(table(
  MajorCellType = as.character(meta$MajorCellType),
  doublet_status = primary_label
), stringsAsFactors = FALSE)
by_major <- by_major[by_major$Freq > 0, , drop = FALSE]
by_major <- dplyr::group_by(by_major, MajorCellType)
by_major <- dplyr::mutate(by_major, TotalMajor = sum(Freq), FractionWithinMajor = Freq / TotalMajor)
by_major <- dplyr::ungroup(by_major)
write_csv(by_major, "D03_primary_doublet_by_major_celltype.csv")

by_cluster <- as.data.frame(table(
  SeuratCluster = as.character(meta$SeuratCluster),
  doublet_status = primary_label
), stringsAsFactors = FALSE)
by_cluster <- by_cluster[by_cluster$Freq > 0, , drop = FALSE]
by_cluster <- dplyr::group_by(by_cluster, SeuratCluster)
by_cluster <- dplyr::mutate(by_cluster, TotalCluster = sum(Freq), FractionWithinCluster = Freq / TotalCluster)
by_cluster <- dplyr::ungroup(by_cluster)
write_csv(by_cluster, "D03_primary_doublet_by_cluster.csv")

qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt"), colnames(meta))
if (length(qc_cols) > 0L) {
  qcd <- meta[, c("Cell", qc_cols), drop = FALSE]
  qcd$doublet_status <- primary_label
  qcl <- tidyr::pivot_longer(qcd, cols = dplyr::all_of(qc_cols), names_to = "metric", values_to = "value")
  qcs <- dplyr::summarise(
    dplyr::group_by(qcl, doublet_status, metric),
    N = sum(!is.na(value)),
    Mean = mean(value, na.rm = TRUE),
    Median = stats::median(value, na.rm = TRUE),
    Q25 = stats::quantile(value, 0.25, na.rm = TRUE),
    Q75 = stats::quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
  write_csv(qcs, "D03_singlet_vs_doublet_QC_metrics.csv")
}

composition <- do.call(rbind, lapply(c("All cells", "Retained singlets", "Predicted doublets"), function(group) {
  keep <- switch(
    group,
    "All cells" = rep(TRUE, nrow(meta)),
    "Retained singlets" = primary_label == "Singlet",
    "Predicted doublets" = primary_label == "Doublet"
  )
  tt <- table(as.character(meta$MajorCellType[keep]))
  data.frame(
    Group = group,
    MajorCellType = names(tt),
    N = as.integer(tt),
    Fraction = as.integer(tt) / sum(tt),
    stringsAsFactors = FALSE
  )
}))
write_csv(composition, "D03_major_celltype_composition_all_singlet_doublet.csv")

final_summary <- data.frame(
  n_cells_all = ncol(obj),
  n_predicted_doublets_primary = observed_primary_n,
  n_retained_singlets = ncol(obj_singlet),
  predicted_doublet_fraction = observed_primary_n / ncol(obj),
  primary_expected_rate = 0.075,
  singlet_object_is_subset_only = TRUE,
  singlet_object_reprocessed = FALSE,
  all_cell_rds = input_rds,
  singlet_rds = singlet_rds,
  stringsAsFactors = FALSE
)
write_csv(final_summary, "D03_final_doublet_filter_summary.csv")
writeLines(capture.output(sessionInfo()), file.path(out_table_dir, paste0(SCRIPT_ID, "_sessionInfo.txt")))
message("09A-03 complete: ", ncol(obj_singlet), " singlets retained from 26,053 cells.")
