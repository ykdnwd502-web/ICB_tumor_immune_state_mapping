############################################################
## Figure_2_Bulk_State_Plotting_Only.R
##
## Purpose:
## 1) LOAD the already processed bulk state scores from CSV.
## 2) Generate Figure 2 panels matching the strict established style:
##    - Integrated tags (A, B, C inside titles)
##    - Left-aligned panel titles
##    - Full black panel borders
##    - Fix title truncation with explicit wide right-margins
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

############################################################
## 0. Paths, Data Loading & Helpers
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
bulk_table_dir <- file.path(project_dir, "results/tables/bulk_discovery")
bulk_figure_dir <- file.path(project_dir, "results/figures/bulk_discovery")
dir.create(bulk_figure_dir, recursive = TRUE, showWarnings = FALSE)

# 安全绘图保存函数
safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  message("Saved figure: ", file)
}

# 强制读取已固化的结果 CSV[cite: 10]
score_file <- file.path(bulk_table_dir, "GSE244982_final_tumor_immune_state_scores.csv")

if (!file.exists(score_file)) {
  stop("找不到已固化的打分文件: ", score_file, "\n请确认 01 脚本是否已经生成过该文件。")
}

message("Loading processed scores from: ", score_file)
final_state_public <- read.csv(score_file, check.names = FALSE)

############################################################
## 1. Unified Theme & Palettes (Strict Style)
############################################################

# 全包围边框、无网格线，Panel 标题内嵌标签并严格左对齐[cite: 10]
theme_icb <- function(base_size = 12.5, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2), 
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6), 
      axis.line = element_blank(), 
      panel.background = element_rect(fill = "white", color = NA)
    )
}

numeric_state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_order <- numeric_state_cols

state_label_map <- c(
  "Immune_defective_Cold" = "immune-defective/\ncold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic\ndifferentiation"
)

state_label_one_line <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

state_colors <- c(
  "Immune_defective_Cold" = "#4DBBD5",
  "Myeloid_Treg_Immunosuppressive" = "#00A087",
  "Tumor_dedifferentiation_Stromal_remodeling" = "#E64B35",
  "Melanocytic_Differentiation" = "#3C5488"
)

############################################################
## 2. Generate Panels
############################################################

## ---------- Panel A: Heatmap ----------
sample_order <- final_state_public %>%
  arrange(
    factor(Dominant_state, levels = state_order),
    desc(Immune_defective_Cold)
  ) %>%
  pull(Sample)

state_long <- final_state_public %>%
  select(Sample, all_of(numeric_state_cols)) %>%
  pivot_longer(
    cols = all_of(numeric_state_cols),
    names_to = "State",
    values_to = "Score"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    StateLabel = factor(state_label_map[State], levels = rev(unname(state_label_map)))
  )

p_state_heat <- ggplot(state_long, aes(x = Sample, y = StateLabel, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_gradient2(
    low = "#3B82F6", mid = "white", high = "#EF4444", midpoint = 0, name = "State\nscore"
  ) +
  labs(
    title = "A  Predefined bulk tumor–immune state scores", 
    x = NULL, y = NULL
  ) +
  theme_icb() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10.75),
    axis.text.y = element_text(size = 11.5)
  )

## ---------- Panel B: PCA ----------
pca_mat <- as.matrix(final_state_public[, numeric_state_cols, drop = FALSE])
pca_mat[is.na(pca_mat)] <- 0
pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
var_exp <- summary(pca)$importance[2, 1:2] * 100

pca_df <- data.frame(
  Sample = final_state_public$Sample,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Dominant_state = final_state_public$Dominant_state,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Dominant_state = factor(Dominant_state, levels = state_order),
    DominantStateLabel = factor(state_label_one_line[as.character(Dominant_state)], levels = state_label_one_line[state_order])
  )

p_pca_clean <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Dominant_state)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = state_colors, breaks = state_order, labels = state_label_one_line[state_order], name = "Dominant state"
  ) +
  labs(
    title = "B  PCA of bulk tumor–immune state scores", 
    x = paste0("PC1 (", sprintf("%.1f", var_exp[1]), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", var_exp[2]), "%)")
  ) +
  theme_icb() +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.75, "cm"),
    legend.text = element_text(lineheight = 0.85),
    plot.margin = margin(t = 5, r = 0, b = 5, l = 5, unit = "pt")
  )

## ---------- Panel C: Bar Plot ----------
dominant_count_plot <- final_state_public %>%
  count(Dominant_state, name = "n_samples") %>%
  mutate(
    Dominant_state = factor(Dominant_state, levels = state_order),
    StateLabel = factor(state_label_map[as.character(Dominant_state)], levels = rev(state_label_map[state_order]))
  )

p_dom <- ggplot(dominant_count_plot, aes(x = StateLabel, y = n_samples, fill = Dominant_state)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = state_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "C  Distribution of dominant tumor–immune state", 
    x = NULL, y = "Number of samples"
  ) +
  theme_icb() +
  theme(
    axis.text.y = element_text(size = 11.5),
    # 强制增加右侧 60pt 保险边距，彻底杜绝长标题被右侧画布截断
    plot.margin = margin(t = 5, r = 60, b = 5, l = 0, unit = "pt") 
  )

############################################################
## 3. Assemble and Save Figure
############################################################

bulk_state_overview <- p_state_heat / (p_pca_clean | p_dom) +
  plot_layout(heights = c(0.92, 1.15), widths = c(0.85, 1.15)) +
  plot_annotation(
    title = "Discovery of predefined bulk tumor–immune states in GSE244982",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16.5, family = "Arial"),
      # 为整张海报的右侧也加上足够的呼吸空间
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt")
    )
  )

fig2_name <- "Figure 2. Predefined tumor–immune state scores in GSE244982"

safe_ggsave(
  file.path(bulk_figure_dir, paste0(fig2_name, ".pdf")),
  bulk_state_overview,
  width = 13.2, height = 9.2
)

safe_ggsave(
  file.path(bulk_figure_dir, paste0(fig2_name, ".jpg")),
  bulk_state_overview,
  width = 13.2, height = 9.2, dpi = 300
)

message("Figure 2 Plotting finished successfully.")