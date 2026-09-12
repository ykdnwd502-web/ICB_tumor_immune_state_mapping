############################################################
## 07_S20_aggregate.R
##
## Purpose:
##   Combine the source-locked canonical reruns for:
##     - GSE250636 Visium
##     - Thrane 2018 legacy-ST
##
##   into the exact five input families required by current
##   Supplementary Figure S20.
##
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)
PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

CANONICAL_ROOT <- file.path(
  PROJECT_DIR,
  "results", "tables",
  "S20_external_spatial_recurrence_CANONICAL_v1.0"
)

GSE_DIR <- file.path(CANONICAL_ROOT, "GSE250636")
THRANE_DIR <- file.path(CANONICAL_ROOT, "Thrane2018_legacyST")
OUT_DIR <- file.path(CANONICAL_ROOT, "combined")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  file.path(GSE_DIR, "S20_GSE250636_key_readout.csv"),
  file.path(GSE_DIR, "S20_GSE250636_gene_coverage_summary.csv"),
  file.path(GSE_DIR, "S20_GSE250636_sample_input_log.csv"),
  file.path(THRANE_DIR, "S20_Thrane2018_key_readout.csv"),
  file.path(THRANE_DIR, "S20_Thrane2018_gene_coverage_summary.csv"),
  file.path(THRANE_DIR, "S20_Thrane2018_sample_input_log.csv")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required canonical S20 dataset outputs:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

read_csv0 <- function(path) {
  data.table::fread(path, data.table = FALSE)
}
write_csv0 <- function(df, name) {
  data.table::fwrite(df, file.path(OUT_DIR, name))
}

gse_key <- read_csv0(file.path(GSE_DIR, "S20_GSE250636_key_readout.csv"))
gse_cov <- read_csv0(file.path(GSE_DIR, "S20_GSE250636_gene_coverage_summary.csv"))
gse_sample <- read_csv0(file.path(GSE_DIR, "S20_GSE250636_sample_input_log.csv"))

thr_key <- read_csv0(file.path(THRANE_DIR, "S20_Thrane2018_key_readout.csv"))
thr_cov <- read_csv0(file.path(THRANE_DIR, "S20_Thrane2018_gene_coverage_summary.csv"))
thr_sample <- read_csv0(file.path(THRANE_DIR, "S20_Thrane2018_sample_input_log.csv"))

dataset_labels <- c(
  GSE250636 = "GSE250636 Visium",
  Thrane2018_legacyST = "Thrane 2018 legacy-ST"
)

adjustment_labels <- c(
  Raw_score = "Raw",
  Stromal_vascular_residualized = "Stromal/vascular residualized",
  Broad_composition_residualized = "Broad composition residualized",
  Tumor_purity_proxy_residualized = "Tumor-purity proxy residualized",
  Composition_PC1_PC2_residualized = "Composition PC1/PC2 residualized"
)

module_labels <- c(
  Myeloid_Treg_Immunosuppressive = "Myeloid–Treg immunosuppressive",
  Tumor_dedifferentiation_Stromal_remodeling =
    "Tumor-dedifferentiation/stromal-remodeling"
)

## ---- Dataset overview ----
gse_nspots_col <- if ("n_in_tissue_spots" %in% names(gse_sample)) {
  "n_in_tissue_spots"
} else {
  "n_spots"
}

gse_ngenes_col <- if ("n_genes" %in% names(gse_sample)) {
  "n_genes"
} else {
  NA_character_
}

thr_ngenes_col <- if ("n_genes" %in% names(thr_sample)) {
  "n_genes"
} else {
  NA_character_
}

dataset_overview <- dplyr::bind_rows(
  data.frame(
    dataset_id = "GSE250636",
    dataset_label = dataset_labels[["GSE250636"]],
    n_samples = length(unique(gse_sample$sample_id)),
    total_spots = sum(gse_sample[[gse_nspots_col]], na.rm = TRUE),
    min_genes = if (!is.na(gse_ngenes_col)) min(gse_sample[[gse_ngenes_col]], na.rm = TRUE) else NA_real_,
    median_genes = if (!is.na(gse_ngenes_col)) median(gse_sample[[gse_ngenes_col]], na.rm = TRUE) else NA_real_,
    max_genes = if (!is.na(gse_ngenes_col)) max(gse_sample[[gse_ngenes_col]], na.rm = TRUE) else NA_real_,
    coordinate_fraction_median = NA_real_,
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset_id = "Thrane2018_legacyST",
    dataset_label = dataset_labels[["Thrane2018_legacyST"]],
    n_samples = length(unique(thr_sample$sample_id)),
    total_spots = sum(thr_sample$n_spots, na.rm = TRUE),
    min_genes = if (!is.na(thr_ngenes_col)) min(thr_sample[[thr_ngenes_col]], na.rm = TRUE) else NA_real_,
    median_genes = if (!is.na(thr_ngenes_col)) median(thr_sample[[thr_ngenes_col]], na.rm = TRUE) else NA_real_,
    max_genes = if (!is.na(thr_ngenes_col)) max(thr_sample[[thr_ngenes_col]], na.rm = TRUE) else NA_real_,
    coordinate_fraction_median = if ("coordinate_fraction" %in% names(thr_sample)) {
      median(thr_sample$coordinate_fraction, na.rm = TRUE)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
)

coverage_combined <- dplyr::bind_rows(
  transform(
    gse_cov,
    dataset_id = "GSE250636",
    dataset_label = dataset_labels[["GSE250636"]]
  ),
  transform(
    thr_cov,
    dataset_id = "Thrane2018_legacyST",
    dataset_label = dataset_labels[["Thrane2018_legacyST"]]
  )
)

decorate_key <- function(df, dataset_id) {
  df$dataset_id <- dataset_id
  df$dataset_label <- unname(dataset_labels[[dataset_id]])
  df$adjustment_label <- unname(adjustment_labels[df$adjustment])
  df
}

key_combined <- dplyr::bind_rows(
  decorate_key(gse_key, "GSE250636"),
  decorate_key(thr_key, "Thrane2018_legacyST")
)

plot_coverage <- coverage_combined %>%
  dplyr::filter(
    module %in% c(
      "Myeloid_Treg_Immunosuppressive",
      "Tumor_dedifferentiation_Stromal_remodeling"
    )
  ) %>%
  dplyr::mutate(
    module_label = unname(module_labels[module])
  )

plot_spearman <- key_combined %>%
  dplyr::filter(readout == "primary_pair_spearman")

plot_dual_high <- key_combined %>%
  dplyr::filter(readout == "primary_pair_dual_high_OR")

plot_moran <- key_combined %>%
  dplyr::filter(readout == "primary_pair_bivariate_Moran_descriptive")

focused <- key_combined %>%
  dplyr::mutate(
    statistic_rounded = round(statistic, 3),
    lower_rounded = round(lower, 3),
    upper_rounded = round(upper, 3),
    p_value_formatted = ifelse(
      is.finite(p_value),
      format.pval(p_value, digits = 3, eps = 1e-300),
      NA_character_
    )
  )

write_csv0(dataset_overview, "S20_dataset_overview.csv")
write_csv0(dataset_overview, "S20_plot_input_overview.csv")
write_csv0(coverage_combined, "S20_combined_gene_coverage_summary.csv")
write_csv0(coverage_combined, "S20_primary_and_composition_gene_coverage.csv")
write_csv0(plot_coverage, "S20_plot_input_coverage.csv")
write_csv0(key_combined, "S20_combined_key_readout.csv")
write_csv0(focused, "S20_primary_pair_focused_summary.csv")
write_csv0(plot_spearman, "S20_plot_input_spearman.csv")
write_csv0(plot_dual_high, "S20_plot_input_dual_high_OR.csv")
write_csv0(plot_moran, "S20_plot_input_moran.csv")

cat("\n============================================================\n")
cat("S20 CANONICAL AGGREGATE COMPLETED\n")
cat("============================================================\n")
cat(
  "GSE250636: ",
  dataset_overview$total_spots[dataset_overview$dataset_id == "GSE250636"],
  " spots / ",
  dataset_overview$n_samples[dataset_overview$dataset_id == "GSE250636"],
  " samples\n",
  sep = ""
)
cat(
  "Thrane2018: ",
  dataset_overview$total_spots[dataset_overview$dataset_id == "Thrane2018_legacyST"],
  " spots / ",
  dataset_overview$n_samples[dataset_overview$dataset_id == "Thrane2018_legacyST"],
  " samples\n",
  sep = ""
)
cat("Key readouts: ", nrow(key_combined), "\n", sep = "")
cat("Output: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("============================================================\n")
