############################################################
## Supplementary Figure S3. Bulk state-score heatmap and raw gene-set z-mean visualization in GSE244982
##
## Purpose:
##   Generate Supplementary Figure S3:
##   A) Predefined bulk tumor–immune state scores heatmap
##   B) Bulk gene-set z-mean scores heatmap
##
## Output:
##   D:/ICB_resistance_project/results/figures/bulk_discovery/
##     - Supplementary Figure S3. Bulk state-score heatmap and raw gene-set z-mean visualization in GSE244982.png
##     - Supplementary Figure S3. Bulk state-score heatmap and raw gene-set z-mean visualization in GSE244982.jpg
##     - Supplementary Figure S3. Bulk state-score heatmap and raw gene-set z-mean visualization in GSE244982.pdf
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
  library(patchwork)
})

############################################################
## 0. User settings
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

bulk_table_dir  <- file.path(project_dir, "results", "tables", "bulk_discovery")
bulk_figure_dir <- file.path(project_dir, "results", "figures", "bulk_discovery")
dir.create(bulk_figure_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

HEATMAP_LOW  <- "#3B82F6"
HEATMAP_MID  <- "#FFFFFF"
HEATMAP_HIGH <- "#EF4444"

show_panel_A_numbers <- FALSE

FIGS3_PANEL_A_TEXT_BOOST <- 1.10
FIGS3_PANEL_B_TEXT_BOOST <- 1.20

############################################################
## 1. Helpers
############################################################

theme_icb_supp <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      plot.tag = element_text(face = "bold", size = base_size * 1.6, color = "black"),
      plot.margin = margin(5, 8, 5, 8)
    )
}

safe_ggsave <- function(file, plot, width, height, dpi = 300, device = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE,
    device = device
  )
  message("Saved: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

############################################################
## 2. Input files (对齐您文件夹中的真实文件名)
############################################################

final_state_file <- file.path(bulk_table_dir, "GSE244982_final_tumor_immune_state_scores.csv")
raw_score_file   <- file.path(bulk_table_dir, "GSE244982_bulk_raw_signature_scores.csv")

if (!file.exists(final_state_file)) {
  stop("找不到状态评分文件: ", final_state_file)
}
if (!file.exists(raw_score_file)) {
  stop("找不到原始基因集得分文件: ", raw_score_file)
}

message("Using final state-score table: ", final_state_file)
message("Using raw gene-set score table: ", raw_score_file)

final_state <- readr::read_csv(final_state_file, show_col_types = FALSE)
raw_scores  <- readr::read_csv(raw_score_file, show_col_types = FALSE)

if (!"Sample" %in% colnames(final_state)) {
  colnames(final_state)[1] <- "Sample"
}
if (!"Sample" %in% colnames(raw_scores)) {
  colnames(raw_scores)[1] <- "Sample"
}

############################################################
## 3. Define standardized labels and sample order
############################################################

state_cols_current <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_cols_legacy <- c(
  "Immune_defective_cold",
  "Myeloid_Treg_immunosuppressive",
  "Tumor_dedifferentiated",
  "Melanocytic_differentiated"
)

if (all(state_cols_current %in% colnames(final_state))) {
  state_cols <- unname(state_cols_current)
} else if (all(state_cols_legacy %in% colnames(final_state))) {
  state_cols <- unname(state_cols_legacy)
} else {
  colnames(final_state) <- gsub("–", "_", colnames(final_state))
  state_cols <- intersect(colnames(final_state), c(
    "Immune_defective_Cold", "Immune_defective_cold",
    "Myeloid_Treg_Immunosuppressive", "Myeloid_Treg_immunosuppressive",
    "Tumor_dedifferentiation_Stromal_remodeling", "Tumor_dedifferentiated",
    "Melanocytic_Differentiation", "Melanocytic_differentiated"
  ))
}

state_display <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation",
  "Immune_defective_cold" = "immune-defective/cold",
  "Myeloid_Treg_immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiated" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_differentiated" = "melanocytic differentiation"
)

sample_order <- final_state %>%
  mutate(.immune_for_order = .data[[state_cols[1]]]) %>%
  arrange(desc(.immune_for_order)) %>%
  pull(Sample)

dominant_col <- intersect(colnames(final_state), c("Dominant_final_state_CLEAN", "Domin_final_state", "Dominant_final_state"))
if (length(dominant_col) > 0) {
  dominant_col <- dominant_col[1]
  dominant_order_map <- c(
    "Immune_defective_Cold" = 1,
    "Immune_defective_cold" = 1,
    "Myeloid_Treg_Immunosuppressive" = 2,
    "Myeloid_Treg_immunosuppressive" = 2,
    "Tumor_dedifferentiation_Stromal_remodeling" = 3,
    "Tumor_dedifferentiated" = 3,
    "Melanocytic_Differentiation" = 4,
    "Melanocytic_differentiated" = 4
  )
  sample_order <- final_state %>%
    mutate(
      .dominant_order = unname(dominant_order_map[as.character(.data[[dominant_col]])]),
      .dominant_order = ifelse(is.na(.dominant_order), 999, .dominant_order),
      .immune_for_order = .data[[state_cols[1]]]
    ) %>%
    arrange(.dominant_order, desc(.immune_for_order)) %>%
    pull(Sample)
}

############################################################
## 4. Panel A: predefined final state-score heatmap
############################################################

state_long <- final_state %>%
  select(Sample, all_of(state_cols)) %>%
  pivot_longer(
    cols = all_of(state_cols),
    names_to = "State",
    values_to = "Score"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    StateLabel = factor(
      unname(state_display[State]),
      levels = rev(unname(state_display[state_cols]))
    )
  )

pA <- ggplot(state_long, aes(x = Sample, y = StateLabel, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = HEATMAP_LOW,
    mid = HEATMAP_MID,
    high = HEATMAP_HIGH,
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish,
    name = "State\nscore"
  ) +
  labs(
    title = "Predefined bulk tumor–immune state scores",
    x = NULL,
    y = NULL
  ) +
  theme_icb_supp(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.2 * FIGS3_PANEL_A_TEXT_BOOST),
    axis.text.y = element_text(size = 9.4 * FIGS3_PANEL_A_TEXT_BOOST, lineheight = 0.95),
    legend.title = element_text(size = 9.2 * FIGS3_PANEL_A_TEXT_BOOST, face = "bold"),
    legend.text = element_text(size = 8.4 * FIGS3_PANEL_A_TEXT_BOOST)
  )

if (isTRUE(show_panel_A_numbers)) {
  pA <- pA + geom_text(aes(label = sprintf("%.2f", Score)), size = 2.0 * FIGS3_PANEL_A_TEXT_BOOST, family = "sans")
}

############################################################
## 5. Panel B: raw gene-set z-mean heatmap with display labels
############################################################

raw_order <- c(
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated",
  "Melanocytic_Differentiation_REFERENCE_hybrid",
  "Melanocytic_Differentiation_REFERENCE_data_driven",
  "Melanocytic_Differentiation_REFERENCE_curated",
  "Immune_defective_Cold_RESTORE_hybrid",
  "Immune_defective_Cold_RESTORE_data_driven",
  "Immune_defective_Cold_RESTORE_curated",
  "Differentiated_inflammatory_REFERENCE",
  "Differentiated_inflammatory_REFERENCE_curated",
  "Differentiated_inflammatory_REFERENCE_hybrid",
  "Differentiated_inflammatory_REFERENCE_data_driven"
)

raw_cols_resolved <- intersect(raw_order, colnames(raw_scores))
if (length(raw_cols_resolved) == 0) {
  raw_cols_resolved <- setdiff(colnames(raw_scores), "Sample")
}

raw_display_map <- c(
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid" = "tumor-dedifferentiation/stromal-remodeling, hybrid",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven" = "tumor-dedifferentiation/stromal-remodeling, data-driven",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated" = "tumor-dedifferentiation/stromal-remodeling, curated",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid" = "myeloid–Treg immunosuppressive, hybrid",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven" = "myeloid–Treg immunosuppressive, data-driven",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated" = "myeloid–Treg immunosuppressive, curated",
  "Melanocytic_Differentiation_REFERENCE_hybrid" = "melanocytic differentiation, hybrid",
  "Melanocytic_Differentiation_REFERENCE_data_driven" = "melanocytic differentiation, data-driven",
  "Melanocytic_Differentiation_REFERENCE_curated" = "melanocytic differentiation, curated",
  "Immune_defective_Cold_RESTORE_hybrid" = "immune-defective/cold, hybrid",
  "Immune_defective_Cold_RESTORE_data_driven" = "immune-defective/cold, data-driven",
  "Immune_defective_Cold_RESTORE_curated" = "immune-defective/cold, curated",
  "Differentiated_inflammatory_REFERENCE" = "differentiated inflammatory reference",
  "Differentiated_inflammatory_REFERENCE_curated" = "differentiated inflammatory reference",
  "Differentiated_inflammatory_REFERENCE_hybrid" = "differentiated inflammatory reference, hybrid",
  "Differentiated_inflammatory_REFERENCE_data_driven" = "differentiated inflammatory reference, data-driven"
)

raw_long <- raw_scores %>%
  select(Sample, all_of(unname(raw_cols_resolved))) %>%
  pivot_longer(
    cols = all_of(unname(raw_cols_resolved)),
    names_to = "GeneSet",
    values_to = "Score"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    GeneSetLabel = raw_display_map[GeneSet],
    GeneSetLabel = ifelse(is.na(GeneSetLabel), str_replace_all(GeneSet, "_", " "), GeneSetLabel)
  )

raw_label_levels <- raw_display_map[raw_cols_resolved]
raw_label_levels <- ifelse(is.na(raw_label_levels), str_replace_all(raw_cols_resolved, "_", " "), raw_label_levels)
raw_label_levels <- unique(as.character(raw_label_levels))

raw_long <- raw_long %>%
  mutate(
    GeneSetLabel = factor(GeneSetLabel, levels = rev(raw_label_levels))
  )

pB <- ggplot(raw_long, aes(x = Sample, y = GeneSetLabel, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.28) +
  scale_fill_gradient2(
    low = HEATMAP_LOW,
    mid = HEATMAP_MID,
    high = HEATMAP_HIGH,
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish,
    name = "Mean\nz-score"
  ) +
  labs(
    title = "Bulk gene-set z-mean scores in GSE244982",
    x = NULL,
    y = NULL
  ) +
  theme_icb_supp(base_size = 9.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.2 * FIGS3_PANEL_B_TEXT_BOOST),
    axis.text.y = element_text(size = 7.1 * FIGS3_PANEL_B_TEXT_BOOST, lineheight = 0.90),
    legend.title = element_text(size = 8.4 * FIGS3_PANEL_B_TEXT_BOOST, face = "bold"),
    legend.text = element_text(size = 7.8 * FIGS3_PANEL_B_TEXT_BOOST)
  )

############################################################
## 6. Combined Supplementary Figure S3
############################################################

supp_s3 <- (pA / pB) +
  plot_layout(heights = c(0.86, 1.20)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 22, family = "sans"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

base_filename <- "Supplementary Figure S3. Bulk state-score heatmap and raw gene-set z-mean visualization in GSE244982"
out_png <- file.path(bulk_figure_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(bulk_figure_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(bulk_figure_dir, paste0(base_filename, ".pdf"))

safe_ggsave(out_png, supp_s3, width = 15.2, height = 9.4, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s3, width = 15.2, height = 9.4, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s3, width = 15.2, height = 9.4, dpi = dpi_out, device = cairo_pdf)

cat("\nSupplementary Figure S3 generated successfully in PNG, JPG, and PDF formats (300 DPI)!\n")