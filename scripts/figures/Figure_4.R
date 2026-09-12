############################################################
## Figure_4_CAF_Plotting_Only.R
##
## Purpose:
## 1) LOAD the already processed and frozen Seurat object (NO recalculation)
## 2) Calculate the CAF/stromal core score exactly from frozen data
## 3) Generate Figure 4 panels matching the requested integrated-label style
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

############################################################
## 0. Paths & Input
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
out_fig_dir <- file.path(project_dir, "results/figures/scRNA_major_annotation")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# 严格锁死输入路径，强制读取保留了旧图精确 UMAP 坐标和聚类的对象
rds_file <- file.path(project_dir, "results", "intermediate", "GSE244983", "GSE244983_seurat_major_annotated.rds")

if (!file.exists(rds_file)) {
  stop("找不到指定的 .rds 文件，请确认路径是否正确: ", rds_file)
}

message("Loading frozen Seurat object to preserve exact UMAP and clusters: ", rds_file)
obj <- readRDS(rds_file)
DefaultAssay(obj) <- "RNA"

############################################################
## 1. Unified Theme & Palettes (Integrated Label Style)
############################################################

# 恢复全包围边框主题，匹配附图中的 integrated label 风格
theme_icb <- function(base_size = 11, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 4), # Panel title 左对齐, 包含 A/B/C/D
      axis.title = element_text(face = "bold", size = base_size + 1.5),
      axis.text = element_text(color = "black", size = base_size),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6), # 恢复全包围黑框
      axis.line = element_blank(), # 移除 L 型单线
      panel.background = element_rect(fill = "white", color = NA)
    )
}

heatmap_low  <- "#3B82F6"
heatmap_mid  <- "white"
heatmap_high <- "#EF4444"

major_celltype_levels <- c(
  "Malignant", "Cycling malignant", "CAF/stromal-like cells",
  "Endothelial", "Myeloid cells", "T/NK cells",
  "T/NK/Treg-like cells", "B/Plasma cells", "Cycling cells"
)

major_celltype_colors_unified <- c(
  "Malignant" = "#3C5488", "Cycling malignant" = "#4DBBD5", "CAF/stromal-like cells" = "#E64B35",
  "Endothelial" = "#00A087", "Myeloid cells" = "#00A087", "T/NK cells" = "#8491B4",
  "T/NK/Treg-like cells" = "#91D1C2", "B/Plasma cells" = "#F39B7F", "Cycling cells" = "#B09C85"
)
major_celltype_colors <- major_celltype_colors_unified[major_celltype_levels]

############################################################
## 2. Calculate CAF/stromal core score exactly
############################################################

caf_markers <- c(
  "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3",
  "DCN", "LUM", "FAP", "ACTA2", "PDPN", "THY1", "VCAN", "TNC", "TGFB1",
  "FBN1", "MRC2", "ADAM12", "TIMP1", "INHBA"
)
caf_markers_present <- intersect(caf_markers, rownames(obj))

# 取 normalized data 进行严格的均值打分
expr_data <- tryCatch(GetAssayData(obj, layer = "data"), error = function(e) GetAssayData(obj, slot = "data"))

# 计算均值并 z-score 标准化
caf_core_raw <- Matrix::colMeans(expr_data[caf_markers_present, , drop = FALSE])
zscore_safe <- function(x) {
  sx <- stats::sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}
obj$CAF_stromal_like_score_z <- zscore_safe(as.numeric(caf_core_raw))

############################################################
## 3. Generate Figure 4 Panels
############################################################

message("Generating Figure panels from frozen data...")

## 准备绘图用的因子水平
obj$MajorCellType <- factor(obj$MajorCellType, levels = major_celltype_levels)
obj$SeuratCluster <- factor(obj$SeuratCluster, levels = sort(unique(as.numeric(as.character(obj$SeuratCluster)))))

## Panel A: Dotplot (标签内嵌)
p_fig4a <- DotPlot(
  obj,
  features = caf_markers_present,
  group.by = "SeuratCluster",
  dot.scale = 6
) +
  scale_color_gradient2(
    low = heatmap_low, mid = heatmap_mid, high = heatmap_high, midpoint = 0, name = "Average\nExpression"
  ) +
  labs(
    title = "A  CAF/stromal markers across Seurat clusters",
    x = NULL, y = "Seurat cluster"
  ) +
  theme_icb() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right"
  )

## Panel B: UMAP (标签内嵌)
p_fig4b <- FeaturePlot(
  obj,
  features = "CAF_stromal_like_score_z",
  reduction = "umap",
  pt.size = 0.12,
  order = TRUE,
  min.cutoff = "q02",
  max.cutoff = "q98",
  cols = c("lightgrey", "#E64B35"),
  raster = FALSE
) +
  coord_equal() + 
  labs(
    title = "B  CAF/stromal core score",
    x = "UMAP 1", y = "UMAP 2",
    color = "CAF/stromal\ncore score"
  ) +
  theme_icb()

## 提取数据用于小提琴图
caf_plot_df <- obj@meta.data

## Panel C: Violin by Major Cell Type (带有内部箱线图，标签内嵌)
p_fig4c <- ggplot(
  caf_plot_df,
  aes(x = MajorCellType, y = CAF_stromal_like_score_z, fill = MajorCellType)
) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.9) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.8, linewidth = 0.25) +
  scale_fill_manual(values = major_celltype_colors, drop = TRUE) +
  labs(
    title = "C  CAF/stromal core score by major cell type",
    x = NULL, y = "CAF/stromal\ncore score"
  ) +
  theme_icb() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

## Panel D: Violin by Seurat Cluster (带有内部箱线图，标签内嵌)
p_fig4d <- ggplot(
  caf_plot_df,
  aes(x = SeuratCluster, y = CAF_stromal_like_score_z, fill = MajorCellType)
) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.9) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.8, linewidth = 0.25) +
  scale_fill_manual(values = major_celltype_colors, drop = TRUE) +
  labs(
    title = "D  CAF/stromal core score by Seurat cluster",
    x = "Seurat cluster", y = "CAF/stromal\ncore score", fill = "Cell type"
  ) +
  theme_icb() +
  theme(legend.position = "right")

## Assemble Full Figure
plot_title_text <- "Single-cell evidence for a CAF/stromal-associated component"
output_file_name <- "Figure 4. Single-cell evidence for a CAF-stromal-associated component"

figure4_caf_reannotation <- (
  p_fig4a /
    (p_fig4b | p_fig4c) /
    p_fig4d
) +
  patchwork::plot_layout(heights = c(1.05, 1.0, 0.95)) +
  patchwork::plot_annotation(
    title = plot_title_text,
    # 移除了自动的 tag_levels，因为已经内嵌到标题中了
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16, family = "Arial")
    )
  )

safe_ggsave(
  file.path(out_fig_dir, paste0(output_file_name, ".pdf")),
  figure4_caf_reannotation,
  width = 15.0, height = 14.6
)
safe_ggsave(
  file.path(out_fig_dir, paste0(output_file_name, ".jpg")),
  figure4_caf_reannotation,
  width = 15.0, height = 14.6, dpi = 300
)

message("Figure 4 Plotting finished successfully.")