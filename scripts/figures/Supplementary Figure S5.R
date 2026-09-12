############################################################
## Supplementary Figure S5. Single-cell quality-control and annotation overview in GSE244983
##
## Purpose:
##   Generate Supplementary Figure S5:
##   A) nCount_RNA vs nFeature_RNA scatter plot
##   B) nFeature_RNA and nCount_RNA violin plots by sample
##   C) UMAP visualization colored by Seurat cluster
##   D) UMAP visualization colored by sample
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S5. Single-cell quality-control and annotation overview in GSE244983.png
##     - Supplementary Figure S5. Single-cell quality-control and annotation overview in GSE244983.jpg
##     - Supplementary Figure S5. Single-cell quality-control and annotation overview in GSE244983.pdf
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(cowplot)
  library(scales)
  library(tibble)
})

############################################################
## 0. Paths and multi-path fallback setup
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

results_dir <- file.path(project_dir, "results")
fig_dir     <- file.path(results_dir, "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## 多路径候选列表：按优先级检查各个可能的子文件夹
candidate_rds_paths <- c(
  file.path(project_dir, "results", "intermediate", "GSE244983", "GSE244983_seurat_major_annotated.rds"),
  file.path(project_dir, "data_processed", "GSE244983_seurat_major_annotated_LOCKED.rds"),
  file.path(project_dir, "data_processed", "GSE244983_seurat_major_annotated.rds"),
  file.path(project_dir, "data_processed", "GSE244983_seurat_QC_integrated_clustered.rds")
)

############################################################
## 1. Locate Seurat object (命中即用)
############################################################

rds_path <- NA_character_
for (path in candidate_rds_paths) {
  norm_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (file.exists(norm_path)) {
    rds_path <- norm_path
    message("Successfully matched Seurat object at: ", rds_path)
    break
  }
}

if (is.na(rds_path)) {
  all_rds <- list.files(project_dir, pattern = "GSE244983.*\\.rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(all_rds) > 0) {
    rds_path <- normalizePath(all_rds[1], winslash = "/", mustWork = FALSE)
    message("Fallback: Found Seurat object via recursive search: ", rds_path)
  } else {
    stop("错误：在项目目录中未找到任何 GSE244983 的 Seurat 对象文件，请检查路径。")
  }
}

message("Using Seurat object: ", rds_path)
obj <- readRDS(rds_path)

if (!inherits(obj, "Seurat")) {
  stop("加载的文件不是有效的 Seurat 对象: ", rds_path)
}

############################################################
## 2. Metadata harmonization
############################################################

meta <- obj@meta.data

cluster_candidates <- c("seurat_clusters", "SeuratCluster", "seurat_cluster", "cluster", "clusters", "RNA_snn_res.0.8")
cluster_col <- cluster_candidates[cluster_candidates %in% colnames(meta)][1]
if (is.na(cluster_col)) {
  stop("未在元数据中找到 Seurat cluster 列。可用列名: ", paste(colnames(meta), collapse = ", "))
}

sample_candidates <- c("Sample", "sample", "patient", "Patient", "patient_id", "donor", "Donor", "orig.ident")
sample_col <- sample_candidates[sample_candidates %in% colnames(meta)][1]
if (is.na(sample_col)) {
  stop("未在元数据中找到 sample 列。可用列名: ", paste(colnames(meta), collapse = ", "))
}

if (!"umap" %in% names(obj@reductions)) {
  stop("Seurat 对象中未找到 umap 降维结果。")
}

umap_embed <- Embeddings(obj, reduction = "umap")
if (ncol(umap_embed) < 2) {
  stop("UMAP 嵌入维度不足两维。")
}

plot_df <- meta %>%
  rownames_to_column("cell_id") %>%
  mutate(
    SeuratCluster = as.character(.data[[cluster_col]]),
    Sample = as.character(.data[[sample_col]])
  )

plot_df$UMAP_1 <- umap_embed[, 1]
plot_df$UMAP_2 <- umap_embed[, 2]

if (!"nCount_RNA" %in% colnames(plot_df)) stop("元数据中缺少 nCount_RNA.")
if (!"nFeature_RNA" %in% colnames(plot_df)) stop("元数据中缺少 nFeature_RNA.")

############################################################
## 3. Factor ordering
############################################################

cluster_num <- suppressWarnings(as.numeric(plot_df$SeuratCluster))
if (all(!is.na(cluster_num))) {
  cluster_levels <- as.character(sort(unique(cluster_num)))
} else {
  cluster_levels <- sort(unique(plot_df$SeuratCluster))
}
plot_df$SeuratCluster <- factor(plot_df$SeuratCluster, levels = cluster_levels)

preferred_sample_order <- c("Pat_ICBnaive1", "Pat_ICBnaive2", "Pat42", "Pat5")
sample_levels <- unique(plot_df$Sample)
if (all(preferred_sample_order %in% sample_levels)) {
  sample_levels <- preferred_sample_order
}
plot_df$Sample <- factor(plot_df$Sample, levels = sample_levels)

############################################################
## 4. Colors and theme (恢复标准字体大小 base_size = 12)
############################################################

n_clusters <- length(levels(plot_df$SeuratCluster))
cluster_cols <- hue_pal()(n_clusters)
names(cluster_cols) <- levels(plot_df$SeuratCluster)

n_samples <- length(levels(plot_df$Sample))
sample_cols <- hue_pal()(n_samples)
names(sample_cols) <- levels(plot_df$Sample)

theme_s5 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 11, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

############################################################
## 5. Panel A: QC scatter
############################################################

pA <- ggplot(plot_df, aes(x = nCount_RNA, y = nFeature_RNA, color = SeuratCluster)) +
  geom_point(size = 0.75, alpha = 0.80) +
  scale_color_manual(values = cluster_cols, name = "Seurat cluster") +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "nCount_RNA vs nFeature_RNA",
    subtitle = "Points colored by Seurat cluster",
    x = "nCount_RNA",
    y = "nFeature_RNA"
  ) +
  theme_s5 +
  theme(legend.position = "right")

############################################################
## 6. Panel B: Violin plots by sample
############################################################

pB1 <- ggplot(plot_df, aes(x = Sample, y = nFeature_RNA, fill = Sample)) +
  geom_violin(scale = "width", color = "grey30", linewidth = 0.25) +
  scale_fill_manual(values = sample_cols) +
  labs(title = "nFeature_RNA", x = "Sample", y = NULL) +
  theme_s5 +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)
  )

pB2 <- ggplot(plot_df, aes(x = Sample, y = nCount_RNA, fill = Sample)) +
  geom_violin(scale = "width", color = "grey30", linewidth = 0.25) +
  scale_fill_manual(values = sample_cols) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "nCount_RNA", x = "Sample", y = NULL) +
  theme_s5 +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)
  )

pB <- plot_grid(pB1, pB2, nrow = 1, align = "hv", labels = NULL)

############################################################
## 7. Panel C: UMAP by Seurat cluster
############################################################

cluster_centers <- plot_df %>%
  group_by(SeuratCluster) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

pC <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = SeuratCluster)) +
  geom_point(size = 0.55, alpha = 0.90) +
  scale_color_manual(values = cluster_cols, name = "Seurat cluster", breaks = cluster_levels) +
  geom_text(
    data = cluster_centers,
    aes(x = UMAP_1, y = UMAP_2, label = SeuratCluster),
    inherit.aes = FALSE,
    size = 4.0,
    color = "black"
  ) +
  labs(
    title = "GSE244983 UMAP by Seurat cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_s5 +
  theme(legend.position = "right")

############################################################
## 8. Panel D: UMAP by sample
############################################################

pD <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Sample)) +
  geom_point(size = 0.55, alpha = 0.90) +
  scale_color_manual(values = sample_cols, name = "Sample") +
  labs(
    title = "GSE244983 UMAP by sample",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_s5 +
  theme(legend.position = "right")

############################################################
## 9. Assemble and save in PNG, JPG, PDF (300 DPI)
############################################################

top_row <- plot_grid(
  pA, pB,
  ncol = 2,
  rel_widths = c(1.15, 1.0),
  labels = c("A", "B"),
  label_size = 22,
  label_fontface = "bold"
)

bottom_row <- plot_grid(
  pC, pD,
  ncol = 2,
  rel_widths = c(1.15, 1.0),
  labels = c("C", "D"),
  label_size = 22,
  label_fontface = "bold"
)

suppfig_s5 <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "v"
)

base_filename <- "Supplementary Figure S5. Single-cell quality-control and annotation overview in GSE244983"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

ggsave(filename = out_png, plot = suppfig_s5, width = 16, height = 12, units = "in", dpi = 300, bg = "white")
ggsave(filename = out_jpg, plot = suppfig_s5, width = 16, height = 12, units = "in", dpi = 300, bg = "white")
ggsave(filename = out_pdf, plot = suppfig_s5, width = 16, height = 12, units = "in", dpi = 300, bg = "white", device = cairo_pdf)

cat("\nSupplementary Figure S5 generated successfully in PNG, JPG, and PDF formats (300 DPI)!\n")