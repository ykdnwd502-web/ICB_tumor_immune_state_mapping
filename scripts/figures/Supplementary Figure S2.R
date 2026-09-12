############################################################
## Supplementary Figure S2
## Bulk PCA and detected-gene assessment in GSE244982
##
## Output:
##   D:/ICB_resistance_project/results/figures/bulk_discovery/
##     - Supplementary Figure S2. Bulk PCA and detected-gene assessment in GSE244982.png
##     - Supplementary Figure S2. Bulk PCA and detected-gene assessment in GSE244982.jpg
##     - Supplementary Figure S2. Bulk PCA and detected-gene assessment in GSE244982.pdf
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
  library(patchwork)
  library(grid)
})

############################################################
## 0. User settings
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

bulk_table_dir  <- file.path(project_dir, "results", "tables", "bulk_discovery")
bulk_figure_dir <- file.path(project_dir, "results", "figures", "bulk_discovery")
out_dir         <- bulk_figure_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

label_all_pca_samples <- TRUE
n_extreme_labels <- 16

fig_width  <- 9.2
fig_height <- 10.0

############################################################
## 1. Helpers
############################################################

theme_icb_supp <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(hjust = 0.5, color = "black"),
      plot.tag = element_text(face = "bold", size = base_size * 1.55, color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      plot.margin = margin(5.5, 8, 5.5, 8)
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
## 2. Input files
############################################################

state_file <- file.path(bulk_table_dir, "GSE244982_final_tumor_immune_state_scores.csv")
qc_file    <- file.path(bulk_table_dir, "GSE244982_bulk_sample_QC.csv")

if (!file.exists(state_file)) {
  stop("找不到状态评分文件: ", state_file)
}
if (!file.exists(qc_file)) {
  stop("找不到QC质控文件: ", qc_file)
}

state_df <- readr::read_csv(state_file, show_col_types = FALSE)
qc_df    <- readr::read_csv(qc_file, show_col_types = FALSE)

############################################################
## 3. PCA from final state-score table
############################################################

candidate_state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

fallback_state_cols <- c(
  "Immune_defective_cold",
  "Myeloid_Treg_immunosuppressive",
  "Tumor_dedifferentiated",
  "Melanocytic_differentiated"
)

if (all(candidate_state_cols %in% colnames(state_df))) {
  state_cols <- candidate_state_cols
} else if (all(fallback_state_cols %in% colnames(state_df))) {
  state_cols <- fallback_state_cols
} else {
  colnames(state_df) <- gsub("–", "_", colnames(state_df))
  state_cols <- intersect(colnames(state_df), c(
    "Immune_defective_Cold", "Immune_defective_cold",
    "Myeloid_Treg_Immunosuppressive", "Myeloid_Treg_immunosuppressive",
    "Tumor_dedifferentiation_Stromal_remodeling", "Tumor_dedifferentiated",
    "Melanocytic_Differentiation", "Melanocytic_differentiated"
  ))
}

if (!"Sample" %in% colnames(state_df)) {
  colnames(state_df)[1] <- "Sample"
}

pca_mat <- as.matrix(state_df[, state_cols, drop = FALSE])
rownames(pca_mat) <- state_df$Sample
pca_mat[is.na(pca_mat)] <- 0

pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
var_exp <- summary(pca)$importance[2, 1:2] * 100

pca_df <- data.frame(
  Sample = rownames(pca_mat),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

pca_df <- pca_df %>%
  mutate(
    pca_radius = sqrt(PC1^2 + PC2^2),
    extreme_rank = rank(-pca_radius, ties.method = "first"),
    pca_label = if (isTRUE(label_all_pca_samples)) {
      Sample
    } else {
      ifelse(extreme_rank <= n_extreme_labels, Sample, "")
    }
  )

############################################################
## 4. Panel A: PCA with ggrepel labels
############################################################

pA <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, color = "#F0F0F0", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "#F0F0F0", linewidth = 0.35) +
  geom_point(size = 2.7, color = "#4C72B0", alpha = 0.92) +
  ggrepel::geom_text_repel(
    aes(label = pca_label),
    size = 3.1,
    family = "sans",
    color = "black",
    min.segment.length = 0,
    segment.color = "grey70",
    segment.linewidth = 0.25,
    box.padding = 0.35,
    point.padding = 0.18,
    max.overlaps = Inf,
    seed = 123
  ) +
  labs(
    title = "PCA of bulk tumor–immune state scores",
    x = paste0("PC1 (", sprintf("%.1f", var_exp[1]), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", var_exp[2]), "%)")
  ) +
  theme_icb_supp(base_size = 11) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10.5)
  )

############################################################
## 5. Panel B: Detected-gene lollipop plot
############################################################

if (!"Sample" %in% colnames(qc_df)) {
  colnames(qc_df)[1] <- "Sample"
}

detected_col <- intersect(
  colnames(qc_df),
  c("n_detected_genes", "Detected_genes", "detected_genes", "nGenes", "n_genes", "n_detected")
)
if (length(detected_col) == 0) {
  detected_col <- colnames(qc_df)[2]
} else {
  detected_col <- detected_col[1]
}

qc_df$n_detected_genes <- as.numeric(qc_df[[detected_col]])

qc_plot_df <- qc_df %>%
  mutate(
    Sample = as.character(Sample),
    Sample_ordered = factor(Sample, levels = Sample[order(n_detected_genes, decreasing = FALSE)])
  )

x_upper <- ceiling(max(qc_plot_df$n_detected_genes, na.rm = TRUE) / 500) * 500

pB <- ggplot(qc_plot_df, aes(x = n_detected_genes, y = Sample_ordered)) +
  geom_segment(
    aes(x = 0, xend = n_detected_genes, y = Sample_ordered, yend = Sample_ordered),
    color = "#B8C7DD",
    linewidth = 0.45
  ) +
  geom_point(
    size = 2.4,
    color = "#4C72B0",
    alpha = 0.95
  ) +
  scale_x_continuous(
    limits = c(0, x_upper),
    breaks = pretty(c(0, x_upper), n = 4),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Detected genes per bulk sample",
    x = "Number of detected genes",
    y = "Bulk sample"
  ) +
  theme_icb_supp(base_size = 10) +
  theme(
    axis.title.x = element_text(size = 11.5, face = "bold"),
    axis.title.y = element_text(size = 11.5, face = "bold"),
    axis.text.x = element_text(size = 9.8),
    axis.text.y = element_text(size = 7.1),
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.25)
  )

############################################################
## 6. Combined Supplementary Figure S2
############################################################

supp_s2 <- (pA / pB) +
  plot_layout(heights = c(1.0, 1.15)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 22, family = "sans"),
      plot.margin = margin(5, 5, 5, 5)
    )
  )

base_filename <- "Supplementary Figure S2. Bulk PCA and detected-gene assessment in GSE244982"
out_png <- file.path(out_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(out_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(out_dir, paste0(base_filename, ".pdf"))

safe_ggsave(out_png, supp_s2, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s2, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s2, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

cat("\nSupplementary Figure S2 generated successfully in PNG, JPG, and PDF formats!\n")