############################################################
## 09A-08. Final structural audit and S7-source manifest
##
## Structural reproducibility gates are blocking.
## Biological/robustness outcomes are reported but are NOT forced to
## match historical values or a prespecified favorable conclusion.
############################################################

options(stringsAsFactors = FALSE)
SCRIPT_ID <- "09A_08_final_module_audit"

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

required_pkgs <- c("readr", "digest")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)

out_table_dir <- file.path(project_dir, "results", "tables", "GSE244983", "doublet_CNV_sensitivity")
out_obj_dir <- file.path(project_dir, "results", "objects", "GSE244983", "doublet_CNV_sensitivity")

f <- function(name) file.path(out_table_dir, name)
required_files <- c(
  f("D01_input_object_audit.csv"),
  f("D02_pANN_global_rate_summary.csv"),
  f("D03_final_doublet_filter_summary.csv"),
  f("D04_state_localization_concordance.csv"),
  f("D04_high_state_doublet_enrichment.csv"),
  f("D05_doublet_rate_nestedness_gates.csv"),
  f("D05_doublet_rate_sensitivity_summary.csv"),
  f("D06_inferCNV_gene_order_audit.csv"),
  f("D06_inferCNV_gene_order_hg38_noheader.txt"),
  f("D07_inferCNV_run_status.csv"),
  f("D07_inferCNV_final_summary.csv"),
  f("D07_inferCNV_CNV_like_signal_by_major_celltype.csv"),
  file.path(out_obj_dir, "GSE244983_all_cells_with_pANN_doublet_calls.rds"),
  file.path(out_obj_dir, "GSE244983_singlets_after_pANN_doublet_filter.rds"),
  file.path(out_obj_dir, "GSE244983_inferCNV_result.rds")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("09A required outputs missing:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
}

read_csv <- function(name) readr::read_csv(f(name), show_col_types = FALSE)
d01 <- read_csv("D01_input_object_audit.csv")
d02 <- read_csv("D02_pANN_global_rate_summary.csv")
d03 <- read_csv("D03_final_doublet_filter_summary.csv")
d04 <- read_csv("D04_state_localization_concordance.csv")
d05g <- read_csv("D05_doublet_rate_nestedness_gates.csv")
d05 <- read_csv("D05_doublet_rate_sensitivity_summary.csv")
d06 <- read_csv("D06_inferCNV_gene_order_audit.csv")
d07s <- read_csv("D07_inferCNV_run_status.csv")
d07 <- read_csv("D07_inferCNV_final_summary.csv")

state_names <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation"
)

gates <- data.frame(
  gate = c(
    "input_26053_cells",
    "input_11616_features",
    "three_pANN_rates_present",
    "primary_rate_0.075_present",
    "final_cardinality_closes",
    "four_major_state_robustness_rows",
    "doublet_rate_calls_nested",
    "three_rate_sensitivity_rows",
    "gene_order_at_least_3000",
    "inferCNV_completed",
    "inferCNV_gene_overlap_at_least_1000",
    "inferCNV_reference_cells_at_least_100"
  ),
  pass = c(
    as.integer(d01$n_cells[[1]]) == 26053L,
    as.integer(d01$n_features[[1]]) == 11616L,
    setequal(round(as.numeric(d02$expected_rate), 3), c(0.05, 0.075, 0.10)),
    any(abs(as.numeric(d02$expected_rate) - 0.075) < 1e-12 & as.logical(d02$is_primary)),
    as.integer(d03$n_predicted_doublets_primary[[1]]) + as.integer(d03$n_retained_singlets[[1]]) == 26053L,
    sum(d04$Grouping == "MajorCellType" & d04$State %in% state_names) == 4L,
    all(as.logical(d05g$pass)),
    nrow(d05) == 3L,
    as.integer(d06$n_gene_order_genes[[1]]) >= 3000L,
    identical(as.character(d07s$status[[1]]), "completed"),
    as.integer(d07$n_genes_used[[1]]) >= 1000L,
    as.integer(d07$n_reference_cells_used[[1]]) >= 100L
  ),
  blocking = TRUE,
  stringsAsFactors = FALSE
)
readr::write_csv(gates, f("D08_09A_final_gate_summary.csv"), na = "")
if (!all(gates$pass)) {
  stop("09A final structural gate failed. Inspect D08_09A_final_gate_summary.csv.", call. = FALSE)
}

major_rob <- d04[d04$Grouping == "MajorCellType", , drop = FALSE]
key_results <- data.frame(
  item = c(
    "primary_predicted_doublets",
    "primary_retained_singlets",
    "primary_predicted_doublet_fraction",
    "major_top_group_unchanged_states",
    "major_state_spearman_min",
    "major_state_spearman_max",
    "major_state_max_abs_delta_max",
    "inferCNV_cells_used",
    "inferCNV_genes_used",
    "inferCNV_reference_cells_used",
    "analysis_role"
  ),
  value = c(
    d03$n_predicted_doublets_primary[[1]],
    d03$n_retained_singlets[[1]],
    d03$predicted_doublet_fraction[[1]],
    sum(as.logical(major_rob$TopGroupUnchanged)),
    min(as.numeric(major_rob$Spearman), na.rm = TRUE),
    max(as.numeric(major_rob$Spearman), na.rm = TRUE),
    max(as.numeric(major_rob$MaxAbsDelta), na.rm = TRUE),
    d07$n_cells_used[[1]],
    d07$n_genes_used[[1]],
    d07$n_reference_cells_used[[1]],
    "revision_stage_technical_sensitivity_not_primary_lineage"
  ),
  stringsAsFactors = FALSE
)
readr::write_csv(key_results, f("D08_09A_key_result_summary.csv"), na = "")

s7_sources <- data.frame(
  file = c(
    "D02_pANN_per_sample_rate_audit.csv",
    "D02_pANN_global_rate_summary.csv",
    "D03_primary_doublet_global_summary.csv",
    "D03_primary_doublet_by_sample.csv",
    "D03_singlet_vs_doublet_QC_metrics.csv",
    "D03_major_celltype_composition_all_singlet_doublet.csv",
    "D04_state_localization_all_vs_singlets.csv",
    "D04_state_localization_concordance.csv",
    "D04_high_state_doublet_enrichment.csv",
    "D05_doublet_rate_sensitivity_summary.csv",
    "D07_inferCNV_input_parameter_audit.csv",
    "D07_inferCNV_reference_and_downsampling_audit.csv",
    "D07_inferCNV_CNV_like_signal_by_cell.csv",
    "D07_inferCNV_CNV_like_signal_by_major_celltype.csv",
    "D07_inferCNV_CNV_like_signal_by_cluster.csv",
    "D07_inferCNV_final_summary.csv"
  ),
  role = c(
    "pANN_per_sample_parameters",
    "pANN_rate_summary",
    "primary_doublet_global_summary",
    "primary_doublet_by_sample",
    "QC_metric_sensitivity",
    "composition_sensitivity",
    "state_localization_before_after",
    "state_localization_concordance",
    "high_state_doublet_enrichment",
    "doublet_rate_sensitivity",
    "inferCNV_parameters",
    "inferCNV_reference_downsampling",
    "inferCNV_cell_level_support",
    "inferCNV_major_celltype_support",
    "inferCNV_cluster_support",
    "inferCNV_final_summary"
  ),
  stringsAsFactors = FALSE
)
s7_sources$path <- file.path(out_table_dir, s7_sources$file)
s7_sources$exists <- file.exists(s7_sources$path)
readr::write_csv(s7_sources, f("D08_S7_source_manifest.csv"), na = "")

manifest_paths <- unique(c(
  list.files(out_table_dir, full.names = TRUE, recursive = FALSE),
  required_files[grepl("[.]rds$", required_files, ignore.case = TRUE)]
))
manifest_paths <- manifest_paths[file.exists(manifest_paths)]
manifest <- data.frame(
  path = manifest_paths,
  bytes = file.info(manifest_paths)$size,
  sha256 = vapply(
    manifest_paths,
    function(x) digest::digest(file = x, algo = "sha256", serialize = FALSE),
    character(1)
  ),
  stringsAsFactors = FALSE
)
readr::write_csv(manifest, f("D08_09A_output_SHA256_manifest.csv"), na = "")
writeLines(capture.output(sessionInfo()), f(paste0(SCRIPT_ID, "_sessionInfo.txt")))
message("09A module PASS: 26,053-cell pANN -> final doublet -> robustness -> inferCNV completed.")
