############################################################
## Supplementary Figure S9. Extended CAF-stromal marker dot plots in GSE244983
##
## Purpose:
##   1) Load Seurat object from: <project>/results/intermediate/GSE244983/GSE244983_seurat_state_localized.rds
##   2) Harmonize major cell-type labels and enforce natural numeric Seurat cluster ordering (0 to max)
##   3) Generate and export Supplementary Figure S9 with unified heatmap color palette and collected legends
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S9. Extended CAF-stromal marker dot plots in GSE244983.png
##     - Supplementary Figure S9. Extended CAF-stromal marker dot plots in GSE244983.jpg
##     - Supplementary Figure S9. Extended CAF-stromal marker dot plots in GSE244983.pdf
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(scales)
  library(patchwork)
})

############################################################
## 0. Paths and settings
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

## 优先指定用户要求的 Seurat 对象路径
manual_rds_path <- file.path(project_dir, "results", "intermediate", "GSE244983", "GSE244983_seurat_state_localized.rds")

HEATMAP_LOW  <- "#3B82F6"
HEATMAP_MID  <- "#FFFFFF"
HEATMAP_HIGH <- "#EF4444"

############################################################
## 1. Helper functions
############################################################

safe_ggsave <- function(file, plot, width, height, dpi = 300, device = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE,
    device = device
  )
  message("Saved: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

find_first_col <- function(meta, patterns, label) {
  meta_cols <- colnames(meta)
  for (pat in patterns) {
    hit <- grep(pat, meta_cols, value = TRUE, ignore.case = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  stop(
    "Could not find column for ", label, ". Tried patterns:\n",
    paste(patterns, collapse = "\n"),
    "\nAvailable metadata columns:\n",
    paste(meta_cols, collapse = ", ")
  )
}

clean_celltype <- function(x) {
  x <- as.character(x)
  dplyr::recode(
    x,
    "T/NK" = "T/NK cells",
    "T_NK" = "T/NK cells",
    "T NK" = "T/NK cells",
    "Myeloid" = "Myeloid cells",
    "Myeloid cells" = "Myeloid cells",
    "B/plasma" = "B/Plasma cells",
    "B/Plasma" = "B/Plasma cells",
    "B_Plasma" = "B/Plasma cells",
    "B/Plasma cells" = "B/Plasma cells",
    "CAF/stromal-like" = "CAF/stromal-like cells",
    "CAF_stromal_like" = "CAF/stromal-like cells",
    "CAF/stromal-like cells" = "CAF/stromal-like cells",
    "Endothelial" = "Endothelial",
    "Malignant" = "Malignant",
    "Cycling malignant" = "Cycling malignant",
    "Cycling_malignant" = "Cycling malignant",
    "T/NK/Treg-like" = "T/NK/Treg-like cells",
    "T_NK_Treg_like" = "T/NK/Treg-like cells",
    "T/NK/Treg-like cells" = "T/NK/Treg-like cells",
    .default = x
  )
}

theme_s9 <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 11.5, face = "bold", color = "black"),
    axis.text = element_text(size = 9.4, color = "black"),
    legend.title = element_text(size = 10.2, face = "bold", color = "black"),
    legend.text = element_text(size = 9.0, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.32),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 6, 5, 6)
  )

match_features_case_insensitive <- function(obj, genes) {
  rn <- rownames(obj)
  rn_upper <- toupper(rn)
  names(rn_upper) <- rn
  
  out <- character(0)
  for (g in genes) {
    hit <- names(rn_upper)[rn_upper == toupper(g)]
    if (length(hit) > 0) out <- c(out, hit[1])
  }
  unique(out)
}

############################################################
## 2. Load Seurat object
############################################################

if (!file.exists(manual_rds_path)) {
  stop("manual_rds_path does not exist: ", manual_rds_path)
}

message("Using Seurat object: ", manual_rds_path)
obj <- readRDS(manual_rds_path)

if (!inherits(obj, "Seurat")) {
  stop("Selected RDS is not a Seurat object: ", manual_rds_path)
}

message("Object cells: ", ncol(obj))
message("Object features: ", nrow(obj))

if ("RNA" %in% names(obj@assays)) {
  DefaultAssay(obj) <- "RNA"
}

meta <- obj@meta.data

############################################################
## 3. Metadata harmonization & cluster sorting
############################################################

celltype_col <- find_first_col(
  meta,
  patterns = c(
    "^MajorCellType_LOCKED$",
    "^MajorCellType$",
    "^major_cell_type$",
    "^major_celltype$",
    "^CellType$",
    "^cell_type$"
  ),
  label = "major cell-type annotation"
)

cluster_col <- find_first_col(
  meta,
  patterns = c(
    "^SeuratCluster$",
    "^seurat_clusters$",
    "^seurat_cluster$",
    "^RNA_snn_res\\.0\\.5$",
    "^cluster$",
    "^clusters$"
  ),
  label = "Seurat cluster annotation"
)

obj$MajorCellType_S9 <- clean_celltype(obj@meta.data[[celltype_col]])

major_celltype_levels <- c(
  "Malignant",
  "Cycling malignant",
  "CAF/stromal-like cells",
  "Endothelial",
  "Myeloid cells",
  "T/NK cells",
  "T/NK/Treg-like cells",
  "B/Plasma cells"
)

observed_major <- unique(as.character(obj$MajorCellType_S9))
major_celltype_levels_final <- c(
  major_celltype_levels[major_celltype_levels %in% observed_major],
  setdiff(sort(observed_major), major_celltype_levels)
)

obj$MajorCellType_S9 <- factor(
  as.character(obj$MajorCellType_S9),
  levels = major_celltype_levels_final
)

## 强制对 Seurat cluster 进行绝对干净的自然数字排序（0, 1, 2, ..., n）
cluster_raw <- as.character(obj@meta.data[[cluster_col]])
cluster_num <- suppressWarnings(as.numeric(cluster_raw))

if (all(!is.na(cluster_num))) {
  cluster_levels <- as.character(sort(unique(cluster_num)))
} else {
  cluster_levels <- sort(unique(cluster_raw))
}

obj$SeuratCluster_S9 <- factor(cluster_raw, levels = cluster_levels)

############################################################
## 4. CAF/stromal marker set
############################################################

caf_markers <- c(
  "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3",
  "DCN", "LUM", "FAP", "ACTA2", "PDPN", "THY1", "VCAN", "TNC", "TGFB1",
  "FBN1", "MRC2", "ADAM12", "TIMP1", "INHBA"
)

caf_markers_present <- match_features_case_insensitive(obj, caf_markers)
caf_markers_missing <- setdiff(caf_markers, toupper(caf_markers_present))

if (length(caf_markers_present) < 2) {
  stop("Too few CAF/stromal markers detected for DotPlot.")
}

############################################################
## 5. Dot plots
############################################################

shared_dot_scale <- scale_color_gradient2(
  low = HEATMAP_LOW,
  mid = HEATMAP_MID,
  high = HEATMAP_HIGH,
  midpoint = 0,
  limits = c(-2, 2),
  oob = scales::squish,
  name = "Average\nExpression"
)

shared_guides <- guides(
  color = guide_colorbar(
    title = "Average\nExpression",
    barheight = unit(2.4, "cm"),
    barwidth = unit(0.35, "cm")
  ),
  size = guide_legend(
    title = "Percent\nExpressed",
    override.aes = list(color = "black")
  )
)

pA <- DotPlot(
  obj,
  features = caf_markers_present,
  group.by = "MajorCellType_S9",
  dot.scale = 5
) +
  shared_dot_scale +
  shared_guides +
  labs(
    title = "CAF/stromal markers across major cell types in GSE244983",
    x = NULL,
    y = "Major cell type"
  ) +
  theme_s9 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.2),
    axis.text.y = element_text(size = 9.6),
    legend.position = "right"
  )

pB <- DotPlot(
  obj,
  features = caf_markers_present,
  group.by = "SeuratCluster_S9",
  dot.scale = 5
) +
  shared_dot_scale +
  shared_guides +
  labs(
    title = "CAF/stromal markers across Seurat clusters in GSE244983",
    x = NULL,
    y = "Seurat cluster"
  ) +
  theme_s9 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.2),
    axis.text.y = element_text(size = 9.6),
    legend.position = "right"
  )

############################################################
## 6. Assemble with collected legends
############################################################

supp_s9 <- (pA / pB) +
  plot_layout(heights = c(0.82, 1.00), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18, family = "sans"),
      plot.margin = margin(4, 4, 4, 4)
    )
  ) &
  theme(legend.position = "right")

############################################################
## 7. Save outputs in PNG, JPG, and PDF (300 DPI)
############################################################

base_filename <- "Supplementary Figure S9. Extended CAF-stromal marker dot plots in GSE244983"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- max(12.0, 0.36 * length(caf_markers_present) + 4.6)
fig_height <- 10.2

safe_ggsave(out_png, supp_s9, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s9, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s9, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

############################################################
## 8. Audit outputs
############################################################

audit <- data.frame(
  seurat_rds = manual_rds_path,
  major_celltype_col = celltype_col,
  cluster_col = cluster_col,
  n_cells = ncol(obj),
  n_features = nrow(obj),
  n_markers_requested = length(caf_markers),
  n_markers_detected = length(caf_markers_present),
  markers_detected = paste(caf_markers_present, collapse = "; "),
  markers_missing = paste(caf_markers_missing, collapse = "; "),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S9_CAF_stromal_marker_dotplots_REVIEWER_FIXED_300dpi_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

cat("\nSupplementary Figure S9 generated successfully in PNG, JPG, and PDF formats (300 DPI)!\n")