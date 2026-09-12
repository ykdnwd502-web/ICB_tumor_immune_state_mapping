############################################################
## 03_build_historical_S28_from_18J.R
##
## INPUT
##   results/tables/revision_external_spatial_recurrence_18J/GSE250636/
##   results/tables/revision_external_spatial_recurrence_18J/Thrane2018_legacyST/
##
## OUTPUT
##   results/tables/revision_external_spatial_recurrence_18J/
##   S28_external_spatial_recurrence_audits/
##
## v4.1 compatibility patch:
##   Map GSE250636 historical `n_in_tissue_spots` -> `n_spots`
##   for the S28 dataset-overview interface only.
##   No scientific/statistical calculation is changed.
##
## Role in sequential reproduction:
##   Step 01 + Step 02 historical 18J outputs
##      -> historical S28 combined tables / plot inputs / figure
##
## Historical scientific/statistical aggregation logic is retained.
## Only PROJECT_DIR is made configurable through ICB_PROJECT_DIR.
############################################################

# ============================================================
# S28 External melanoma spatial recurrence audits
# Table + Supplementary Figure S28 generator
#
# Purpose:
#   Combine GSE250636 Visium and Thrane 2018 legacy-ST external
#   melanoma spatial recurrence audits into one Supplementary Table S28
#   and one Supplementary Figure S28.
#
# Outputs:
#   1) historical_external_spatial_recurrence_summary.xlsx
#   2) historical_external_spatial_recurrence_audit.png
#   3) historical_external_spatial_recurrence_audit.pdf
#   4) Combined CSV audit inputs
#
# Interpretation boundary:
#   External spatial recurrence audit only.
#   Not ICB-response validation.
#   Not validated predictive biomarker.
#   Not functional resistance niche validation.
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

BASE_OUT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "revision_external_spatial_recurrence_18J"
)

S28_OUT_DIR <- file.path(BASE_OUT_DIR, "S28_external_spatial_recurrence_audits")

dir.create(S28_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Expected dataset output directories
GSE_DIR <- file.path(BASE_OUT_DIR, "GSE250636")
THRANE_DIR <- file.path(BASE_OUT_DIR, "Thrane2018_legacyST")

# If your GSE250636 folder has a slightly different name, the script will search automatically.
DATASET_LABELS <- c(
  GSE250636 = "GSE250636 Visium",
  Thrane2018_legacyST = "Thrane 2018 legacy-ST"
)

PRIMARY_X_STATE <- "Myeloid_Treg_Immunosuppressive"
PRIMARY_Y_STATE <- "Tumor_dedifferentiation_Stromal_remodeling"

PRIMARY_X_LABEL <- "Myeloid–Treg immunosuppressive"
PRIMARY_Y_LABEL <- "Tumor-dedifferentiation/stromal-remodeling"

PRIMARY_PAIR_LABEL <- paste(PRIMARY_X_LABEL, "vs", PRIMARY_Y_LABEL)

# Adjustment order for table and figure
ADJ_ORDER <- c(
  "Raw_score",
  "Stromal_vascular_residualized",
  "Broad_composition_residualized",
  "Tumor_purity_proxy_residualized",
  "Composition_PC1_PC2_residualized"
)

ADJ_LABELS <- c(
  Raw_score = "Raw",
  Stromal_vascular_residualized = "Stromal/vascular residualized",
  Broad_composition_residualized = "Broad composition residualized",
  Tumor_purity_proxy_residualized = "Tumor-purity proxy residualized",
  Composition_PC1_PC2_residualized = "Composition PC1/PC2 residualized"
)

# -----------------------------
# 1. Packages
# -----------------------------

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "stringr",
  "openxlsx", "ggplot2", "patchwork", "scales"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# -----------------------------
# 2. Helper functions
# -----------------------------

message2 <- function(...) {
  cat(paste0("[S28] ", paste0(..., collapse = ""), "\n"))
}

safe_read_csv <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    warning("File not found: ", path)
    return(data.frame())
  }

  as.data.frame(
    data.table::fread(path, check.names = FALSE, data.table = FALSE),
    check.names = FALSE
  )
}

write_csv <- function(df, filename) {
  data.table::fwrite(df, file.path(S28_OUT_DIR, filename))
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

find_first_file <- function(root_dirs, pattern) {
  root_dirs <- root_dirs[dir.exists(root_dirs)]

  if (length(root_dirs) == 0) {
    return(NA_character_)
  }

  files <- unlist(lapply(root_dirs, function(r) {
    list.files(
      r,
      recursive = TRUE,
      full.names = TRUE,
      pattern = pattern
    )
  }))

  files <- unique(files)

  if (length(files) == 0) {
    return(NA_character_)
  }

  files[1]
}

find_dataset_dir <- function(base_dir, keyword) {
  candidates <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

  hit <- candidates[
    grepl(keyword, basename(candidates), ignore.case = TRUE)
  ]

  if (length(hit) > 0) {
    return(hit[1])
  }

  return(file.path(base_dir, keyword))
}

normalize_key_readout <- function(df, dataset_id, dataset_label) {
  if (nrow(df) == 0) {
    return(data.frame())
  }

  required_cols <- c(
    "readout", "dataset", "sample_id", "adjustment",
    "statistic_label", "statistic", "lower", "upper",
    "p_value", "n_spots"
  )

  for (cc in required_cols) {
    if (!cc %in% colnames(df)) {
      df[[cc]] <- NA
    }
  }

  out <- df[, required_cols, drop = FALSE]

  out$dataset_id <- dataset_id
  out$dataset_label <- dataset_label

  out$adjustment <- as.character(out$adjustment)
  out$adjustment_label <- ifelse(
    out$adjustment %in% names(ADJ_LABELS),
    ADJ_LABELS[out$adjustment],
    out$adjustment
  )

  out$adjustment <- factor(out$adjustment, levels = ADJ_ORDER)
  out$adjustment_label <- factor(
    out$adjustment_label,
    levels = unname(ADJ_LABELS[ADJ_ORDER])
  )

  out$readout <- as.character(out$readout)
  out$statistic <- as.numeric(out$statistic)
  out$lower <- as.numeric(out$lower)
  out$upper <- as.numeric(out$upper)
  out$p_value <- as.numeric(out$p_value)
  out$n_spots <- as.integer(out$n_spots)

  out
}

normalize_coverage <- function(df, dataset_id, dataset_label) {
  if (nrow(df) == 0) {
    return(data.frame())
  }

  # Allow both coverage summary and full coverage formats
  if (!"module" %in% colnames(df)) {
    return(data.frame())
  }

  if (!"median_detected" %in% colnames(df)) {
    if ("n_genes_detected" %in% colnames(df)) {
      df <- df %>%
        dplyr::group_by(module) %>%
        dplyr::summarise(
          min_detected = min(n_genes_detected, na.rm = TRUE),
          median_detected = median(n_genes_detected, na.rm = TRUE),
          max_detected = max(n_genes_detected, na.rm = TRUE),
          min_coverage = min(coverage_fraction, na.rm = TRUE),
          median_coverage = median(coverage_fraction, na.rm = TRUE),
          max_coverage = max(coverage_fraction, na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      return(data.frame())
    }
  }

  req <- c(
    "module", "min_detected", "median_detected", "max_detected",
    "min_coverage", "median_coverage", "max_coverage"
  )

  for (cc in req) {
    if (!cc %in% colnames(df)) {
      df[[cc]] <- NA
    }
  }

  out <- df[, req, drop = FALSE]

  out$dataset_id <- dataset_id
  out$dataset_label <- dataset_label

  out$module <- as.character(out$module)
  out$median_detected <- as.numeric(out$median_detected)
  out$median_coverage <- as.numeric(out$median_coverage)

  out
}

normalize_sample_log <- function(df, dataset_id, dataset_label) {
  if (nrow(df) == 0) {
    return(data.frame())
  }

  ## Sequential-public compatibility:
  ## historical GSE250636 v3 names the spot-count column
  ## `n_in_tissue_spots`, whereas this historical S28 aggregator expects
  ## `n_spots`. They represent the same in-tissue spot count.
  if (
    !"n_spots" %in% colnames(df) &&
    "n_in_tissue_spots" %in% colnames(df)
  ) {
    df$n_spots <- df$n_in_tissue_spots
  }

  req <- c(
    "dataset", "sample_id", "source_file", "path",
    "n_genes", "n_spots", "n_modules_scored",
    "coordinate_inference", "coordinate_fraction"
  )

  for (cc in req) {
    if (!cc %in% colnames(df)) {
      df[[cc]] <- NA
    }
  }

  out <- df[, req, drop = FALSE]

  out$dataset_id <- dataset_id
  out$dataset_label <- dataset_label
  out$n_genes <- as.integer(out$n_genes)
  out$n_spots <- as.integer(out$n_spots)

  out
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 1e-300, "<1e-300", formatC(p, format = "e", digits = 2))
  )
}

# -----------------------------
# 3. Locate input files
# -----------------------------

# Try exact directories first; otherwise search by keyword
if (!dir.exists(GSE_DIR)) {
  GSE_DIR <- find_dataset_dir(BASE_OUT_DIR, "GSE250636")
}

if (!dir.exists(THRANE_DIR)) {
  THRANE_DIR <- find_dataset_dir(BASE_OUT_DIR, "Thrane2018")
}

message2("GSE250636 directory: ", GSE_DIR)
message2("Thrane directory: ", THRANE_DIR)

gse_key_file <- find_first_file(
  c(GSE_DIR, BASE_OUT_DIR),
  "18J_GSE250636.*key_readout.*\\.csv$"
)

gse_cov_file <- find_first_file(
  c(GSE_DIR, BASE_OUT_DIR),
  "18J_GSE250636.*gene_coverage_summary.*\\.csv$"
)

if (is.na(gse_cov_file)) {
  gse_cov_file <- find_first_file(
    c(GSE_DIR, BASE_OUT_DIR),
    "18J_GSE250636.*gene_coverage.*\\.csv$"
  )
}

gse_sample_file <- find_first_file(
  c(GSE_DIR, BASE_OUT_DIR),
  "18J_GSE250636.*sample_input_log.*\\.csv$"
)

thrane_key_file <- find_first_file(
  c(THRANE_DIR, BASE_OUT_DIR),
  "18J_Thrane2018.*key_readout.*\\.csv$"
)

thrane_cov_file <- find_first_file(
  c(THRANE_DIR, BASE_OUT_DIR),
  "18J_Thrane2018.*gene_coverage_summary.*\\.csv$"
)

thrane_sample_file <- find_first_file(
  c(THRANE_DIR, BASE_OUT_DIR),
  "18J_Thrane2018.*sample_input_log.*\\.csv$"
)

thrane_inclusion_file <- find_first_file(
  c(THRANE_DIR, BASE_OUT_DIR),
  "18J_Thrane2018.*inclusion_decision_helper.*\\.csv$"
)

thrane_matrix_candidate_file <- find_first_file(
  c(THRANE_DIR, BASE_OUT_DIR),
  "18J_Thrane2018.*matrix_candidate_orientation_log.*\\.csv$"
)

message2("Input files:")
print(data.frame(
  item = c(
    "gse_key_file", "gse_cov_file", "gse_sample_file",
    "thrane_key_file", "thrane_cov_file", "thrane_sample_file",
    "thrane_inclusion_file", "thrane_matrix_candidate_file"
  ),
  path = c(
    gse_key_file, gse_cov_file, gse_sample_file,
    thrane_key_file, thrane_cov_file, thrane_sample_file,
    thrane_inclusion_file, thrane_matrix_candidate_file
  ),
  stringsAsFactors = FALSE
))

# -----------------------------
# 4. Read and combine inputs
# -----------------------------

gse_key <- safe_read_csv(gse_key_file)
gse_cov <- safe_read_csv(gse_cov_file)
gse_sample <- safe_read_csv(gse_sample_file)

thrane_key <- safe_read_csv(thrane_key_file)
thrane_cov <- safe_read_csv(thrane_cov_file)
thrane_sample <- safe_read_csv(thrane_sample_file)
thrane_inclusion <- safe_read_csv(thrane_inclusion_file)
thrane_matrix_candidate <- safe_read_csv(thrane_matrix_candidate_file)

key_combined <- dplyr::bind_rows(
  normalize_key_readout(gse_key, "GSE250636", DATASET_LABELS[["GSE250636"]]),
  normalize_key_readout(thrane_key, "Thrane2018_legacyST", DATASET_LABELS[["Thrane2018_legacyST"]])
)

coverage_combined <- dplyr::bind_rows(
  normalize_coverage(gse_cov, "GSE250636", DATASET_LABELS[["GSE250636"]]),
  normalize_coverage(thrane_cov, "Thrane2018_legacyST", DATASET_LABELS[["Thrane2018_legacyST"]])
)

sample_combined <- dplyr::bind_rows(
  normalize_sample_log(gse_sample, "GSE250636", DATASET_LABELS[["GSE250636"]]),
  normalize_sample_log(thrane_sample, "Thrane2018_legacyST", DATASET_LABELS[["Thrane2018_legacyST"]])
)

write_csv(key_combined, "S28_combined_key_readout.csv")
write_csv(coverage_combined, "S28_combined_gene_coverage_summary.csv")
write_csv(sample_combined, "S28_combined_sample_input_log.csv")

# -----------------------------
# 5. Dataset overview
# -----------------------------

dataset_overview <- sample_combined %>%
  dplyr::group_by(dataset_id, dataset_label) %>%
  dplyr::summarise(
    n_samples = dplyr::n_distinct(sample_id),
    total_spots = sum(n_spots, na.rm = TRUE),
    min_genes = min(n_genes, na.rm = TRUE),
    median_genes = median(n_genes, na.rm = TRUE),
    max_genes = max(n_genes, na.rm = TRUE),
    coordinate_fraction_median = median(coordinate_fraction, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(dataset_overview, "S28_dataset_overview.csv")

# -----------------------------
# 6. Primary-pair focused summary
# -----------------------------

primary_key_summary <- key_combined %>%
  dplyr::filter(
    readout %in% c(
      "primary_pair_spearman",
      "primary_pair_dual_high_OR",
      "primary_pair_bivariate_Moran_descriptive"
    )
  ) %>%
  dplyr::mutate(
    p_value_formatted = format_p(p_value),
    statistic_rounded = round(statistic, 3),
    lower_rounded = round(lower, 3),
    upper_rounded = round(upper, 3)
  ) %>%
  dplyr::arrange(dataset_label, readout, adjustment)

write_csv(primary_key_summary, "S28_primary_pair_focused_summary.csv")

primary_coverage_summary <- coverage_combined %>%
  dplyr::filter(
    module %in% c(
      PRIMARY_X_STATE,
      PRIMARY_Y_STATE,
      "CAF_Stromal",
      "Endothelial",
      "Myeloid",
      "Melanoma_Lineage",
      "T_NK",
      "B_Plasma"
    )
  ) %>%
  dplyr::mutate(
    module_label = dplyr::case_when(
      module == PRIMARY_X_STATE ~ PRIMARY_X_LABEL,
      module == PRIMARY_Y_STATE ~ PRIMARY_Y_LABEL,
      module == "CAF_Stromal" ~ "CAF/stromal",
      module == "Endothelial" ~ "Endothelial",
      module == "Myeloid" ~ "Myeloid",
      module == "Melanoma_Lineage" ~ "Melanoma lineage",
      module == "T_NK" ~ "T/NK",
      module == "B_Plasma" ~ "B/plasma",
      TRUE ~ module
    )
  )

write_csv(primary_coverage_summary, "S28_primary_and_composition_gene_coverage.csv")

# -----------------------------
# 7. Create Supplementary Table S28
# -----------------------------

readme <- data.frame(
  Item = c(
    "Supplementary table",
    "Purpose",
    "Primary state pair",
    "Datasets",
    "Response annotation",
    "Interpretation boundary",
    "Main conclusion",
    "Generated on"
  ),
  Description = c(
    "Supplementary Table S28. External melanoma spatial recurrence audits of composition-linked state-score co-variation.",
    "To evaluate whether the primary state-score spatial co-variation pattern recurs in external melanoma spatial transcriptomics resources.",
    PRIMARY_PAIR_LABEL,
    "GSE250636 Visium and Thrane et al. 2018 legacy spatial transcriptomics.",
    "No matched ICB-response annotation was used or assumed.",
    "External spatial recurrence audit only; not ICB-response validation, not predictive biomarker validation, and not functional niche validation.",
    "Both external spatial resources support recurrent raw co-variation and dual-high enrichment for the primary pair, with substantial attenuation after composition-aware adjustment.",
    as.character(Sys.time())
  ),
  stringsAsFactors = FALSE
)

input_inventory <- data.frame(
  item = c(
    "GSE250636 key readout",
    "GSE250636 gene coverage",
    "GSE250636 sample log",
    "Thrane 2018 key readout",
    "Thrane 2018 gene coverage",
    "Thrane 2018 sample log",
    "Thrane 2018 inclusion decision",
    "Thrane 2018 matrix candidate log"
  ),
  path = c(
    gse_key_file,
    gse_cov_file,
    gse_sample_file,
    thrane_key_file,
    thrane_cov_file,
    thrane_sample_file,
    thrane_inclusion_file,
    thrane_matrix_candidate_file
  ),
  exists = file.exists(c(
    gse_key_file,
    gse_cov_file,
    gse_sample_file,
    thrane_key_file,
    thrane_cov_file,
    thrane_sample_file,
    thrane_inclusion_file,
    thrane_matrix_candidate_file
  )),
  stringsAsFactors = FALSE
)

wb <- openxlsx::createWorkbook()

add_sheet(wb, "README", readme)
add_sheet(wb, "Input_inventory", input_inventory)
add_sheet(wb, "Dataset_overview", dataset_overview)
add_sheet(wb, "Primary_key_summary", primary_key_summary)
add_sheet(wb, "Primary_gene_coverage", primary_coverage_summary)
add_sheet(wb, "Combined_key_readout", key_combined)
add_sheet(wb, "Combined_gene_coverage", coverage_combined)
add_sheet(wb, "Combined_sample_log", sample_combined)
add_sheet(wb, "GSE250636_key_readout", gse_key)
add_sheet(wb, "GSE250636_gene_coverage", gse_cov)
add_sheet(wb, "GSE250636_sample_log", gse_sample)
add_sheet(wb, "Thrane_key_readout", thrane_key)
add_sheet(wb, "Thrane_gene_coverage", thrane_cov)
add_sheet(wb, "Thrane_sample_log", thrane_sample)
add_sheet(wb, "Thrane_inclusion", thrane_inclusion)
add_sheet(wb, "Thrane_matrix_log", thrane_matrix_candidate)

table_path <- file.path(
  S28_OUT_DIR,
  "historical_external_spatial_recurrence_summary.xlsx"
)

openxlsx::saveWorkbook(
  wb,
  table_path,
  overwrite = TRUE
)

message2("Saved table: ", table_path)

# -----------------------------
# 8. Prepare plotting data
# -----------------------------

# Panel A: dataset overview
plot_overview <- dataset_overview %>%
  dplyr::mutate(
    dataset_label = factor(
      dataset_label,
      levels = c("GSE250636 Visium", "Thrane 2018 legacy-ST")
    )
  )

# Panel B: primary-pair coverage
plot_coverage <- primary_coverage_summary %>%
  dplyr::filter(module %in% c(PRIMARY_X_STATE, PRIMARY_Y_STATE)) %>%
  dplyr::mutate(
    dataset_label = factor(
      dataset_label,
      levels = c("GSE250636 Visium", "Thrane 2018 legacy-ST")
    ),
    module_label = factor(
      module_label,
      levels = c(PRIMARY_X_LABEL, PRIMARY_Y_LABEL)
    )
  )

# Panel C: Spearman rho
plot_rho <- key_combined %>%
  dplyr::filter(readout == "primary_pair_spearman") %>%
  dplyr::filter(adjustment %in% ADJ_ORDER) %>%
  dplyr::mutate(
    dataset_label = factor(
      dataset_label,
      levels = c("GSE250636 Visium", "Thrane 2018 legacy-ST")
    ),
    adjustment_label = factor(
      adjustment_label,
      levels = unname(ADJ_LABELS[ADJ_ORDER])
    )
  )

# Panel D: dual-high OR
plot_or <- key_combined %>%
  dplyr::filter(readout == "primary_pair_dual_high_OR") %>%
  dplyr::filter(adjustment %in% ADJ_ORDER) %>%
  dplyr::mutate(
    dataset_label = factor(
      dataset_label,
      levels = c("GSE250636 Visium", "Thrane 2018 legacy-ST")
    ),
    adjustment_label = factor(
      adjustment_label,
      levels = unname(ADJ_LABELS[ADJ_ORDER])
    )
  )

# Panel E optional: Moran-type I
plot_moran <- key_combined %>%
  dplyr::filter(readout == "primary_pair_bivariate_Moran_descriptive") %>%
  dplyr::filter(adjustment %in% ADJ_ORDER) %>%
  dplyr::mutate(
    dataset_label = factor(
      dataset_label,
      levels = c("GSE250636 Visium", "Thrane 2018 legacy-ST")
    ),
    adjustment_label = factor(
      adjustment_label,
      levels = unname(ADJ_LABELS[ADJ_ORDER])
    )
  )

write_csv(plot_overview, "S28_plot_input_overview.csv")
write_csv(plot_coverage, "S28_plot_input_coverage.csv")
write_csv(plot_rho, "S28_plot_input_spearman.csv")
write_csv(plot_or, "S28_plot_input_dual_high_OR.csv")
write_csv(plot_moran, "S28_plot_input_moran.csv")

# -----------------------------
# 9. Plot style
# -----------------------------

theme_s28 <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, color = "#E8E8E8"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey60"),
      strip.text = ggplot2::element_text(face = "bold", size = 11.5, color = "black"),
      axis.title = ggplot2::element_text(face = "bold", size = 12, color = "black"),
      axis.text = ggplot2::element_text(size = 10.5, color = "black"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1, size = 10.25),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 14.5, color = "black"),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 11, color = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 10.5)
    )
}

# -----------------------------
# 10. Build panels
# -----------------------------

pA <- ggplot(plot_overview, aes(x = dataset_label, y = total_spots)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0("n=", total_spots, "\n", n_samples, " samples")),
    vjust = -0.25,
    size = 2.8
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "A. External spatial datasets",
    x = NULL,
    y = "Total spatial spots"
  ) +
  theme_s28()

pB <- ggplot(
  plot_coverage,
  aes(x = module_label, y = median_coverage, group = dataset_label)
) +
  geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
  geom_errorbar(
    aes(ymin = min_coverage, ymax = max_coverage),
    width = 0.15,
    position = position_dodge(width = 0.45)
  ) +
  facet_wrap(~ dataset_label, nrow = 1) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05)
  ) +
  labs(
    title = "B. Primary signature gene coverage",
    x = NULL,
    y = "Detected gene coverage"
  ) +
  theme_s28()

pC <- ggplot(
  plot_rho,
  aes(x = adjustment_label, y = statistic, group = dataset_label)
) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(size = 2.2) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "C. Primary-pair Spearman correlation",
    x = NULL,
    y = "Spearman rho"
  ) +
  theme_s28()

pD <- ggplot(
  plot_or,
  aes(x = adjustment_label, y = statistic, group = dataset_label)
) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dashed") +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.15,
    linewidth = 0.35
  ) +
  facet_wrap(~ dataset_label, nrow = 1) +
  scale_y_continuous(trans = "log10") +
  labs(
    title = "D. Primary-pair dual-high enrichment",
    x = NULL,
    y = "Fisher OR, log10 scale"
  ) +
  theme_s28()

pE <- ggplot(
  plot_moran,
  aes(x = adjustment_label, y = statistic, group = dataset_label)
) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(size = 2.2) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ dataset_label, nrow = 1) +
  labs(
    title = "E. Descriptive bivariate Moran-type statistic",
    x = NULL,
    y = "Moran-type I"
  ) +
  theme_s28()

# If Moran is unavailable for one dataset, still plot available values.
# Final layout: A+B on top, C+D+E below.
fig_s28 <- (pA | pB) / (pC | pD) / pE +
  patchwork::plot_annotation(
    title = "Supplementary Figure S28. External spatial recurrence audits in melanoma spatial transcriptomics datasets",
    subtitle = paste0(
      "Primary state pair: ",
      PRIMARY_PAIR_LABEL,
      ". Analyses are spatial recurrence audits and do not use matched ICB-response annotation."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15.5),
      plot.subtitle = ggplot2::element_text(size = 11.5)
    )
  )

# -----------------------------
# 11. Save figure
# -----------------------------

png_path <- file.path(
  S28_OUT_DIR,
  "historical_external_spatial_recurrence_audit.png"
)

pdf_path <- file.path(
  S28_OUT_DIR,
  "historical_external_spatial_recurrence_audit.pdf"
)

tiff_path <- file.path(
  S28_OUT_DIR,
  "Supplementary_Figure_S28_external_spatial_recurrence_audits.tiff"
)

ggplot2::ggsave(
  filename = png_path,
  plot = fig_s28,
  width = 12,
  height = 11,
  dpi = 450,
  bg = "white"
)

ggplot2::ggsave(
  filename = pdf_path,
  plot = fig_s28,
  width = 12,
  height = 11,
  bg = "white"
)

ggplot2::ggsave(
  filename = tiff_path,
  plot = fig_s28,
  width = 12,
  height = 11,
  dpi = 450,
  compression = "lzw",
  bg = "white"
)

message2("Saved figure PNG: ", png_path)
message2("Saved figure PDF: ", pdf_path)
message2("Saved figure TIFF: ", tiff_path)

# -----------------------------
# 12. Console summary
# -----------------------------

message2("S28 generation completed.")
message2("Output directory: ", S28_OUT_DIR)

message2("Dataset overview:")
print(dataset_overview)

message2("Primary key summary:")
print(primary_key_summary)

message2("Primary coverage summary:")
print(primary_coverage_summary)

message2("Supplementary Table S28:")
print(table_path)

message2("Supplementary Figure S28 PNG:")
print(png_path)






# ------------------------------------------------------------
# Public sequential end gate
# ------------------------------------------------------------
EXPECTED_S28_CORE <- c(
  "S28_plot_input_overview.csv",
  "S28_dataset_overview.csv",
  "S28_combined_gene_coverage_summary.csv",
  "S28_primary_and_composition_gene_coverage.csv",
  "S28_plot_input_coverage.csv",
  "S28_plot_input_spearman.csv",
  "S28_primary_pair_focused_summary.csv",
  "S28_combined_key_readout.csv",
  "S28_plot_input_dual_high_OR.csv",
  "S28_plot_input_moran.csv"
)

missing_s28 <- EXPECTED_S28_CORE[
  !file.exists(
    file.path(
      S28_OUT_DIR,
      EXPECTED_S28_CORE
    )
  )
]

if (length(missing_s28)) {
  stop(
    "STEP 03 finished but expected historical S28 file(s) are missing:\n",
    paste(missing_s28, collapse = "\n"),
    call. = FALSE
  )
}

cat(
  "\n============================================================\n",
  "STEP 03 PASS — regenerated 18J -> historical S28 layer\n",
  "Historical S28 core CSV outputs: 10/10 present\n",
  "Output: ", S28_OUT_DIR, "\n",
  "============================================================\n",
  sep = ""
)
