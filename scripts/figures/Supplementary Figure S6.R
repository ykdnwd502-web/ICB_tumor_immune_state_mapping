############################################################
## Supplementary Figure S6. Cell-type composition by sample in GSE244983
##
## Purpose:
##   1) Import or load GSE244983 single-cell RNA-seq data with multi-path fallback
##   2) Perform QC, normalization, dimensionality reduction, clustering, and UMAP
##   3) Annotate major cell types using marker-module scores
##   4) Generate and export Supplementary Figure S6 (Cell-type composition by sample)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S6. Cell-type composition by sample in GSE244983.png
##     - Supplementary Figure S6. Cell-type composition by sample in GSE244983.jpg
##     - Supplementary Figure S6. Cell-type composition by sample in GSE244983.pdf
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(bitmapType = "cairo")

############################################################
## 0. Project directory and folders
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

data_raw_dir       <- file.path(project_dir, "data_raw")
data_processed_dir <- file.path(project_dir, "data_processed")
out_table_dir      <- file.path(project_dir, "results/tables/scRNA_major_annotation")
out_fig_dir        <- file.path(project_dir, "results/figures/supplementary")
log_dir            <- file.path(project_dir, "logs")

dir.create(data_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Locked parameters
############################################################

input_mode <- "auto"
qc_min_features <- 200
qc_max_features <- Inf
qc_min_counts   <- 0
qc_max_counts   <- Inf
qc_max_percent_mt <- 25

n_variable_features <- 2000
n_pcs_to_compute <- 30
n_pcs_to_use <- 20
cluster_resolution <- 0.5
use_harmony_if_available <- TRUE

manual_cluster_annotation <- c(
  "12" = "CAF/stromal-like cells"
)

############################################################
## 2. Packages and helper functions
############################################################

required_pkgs <- c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "patchwork",
  "dplyr",
  "tidyr",
  "stringr",
  "data.table",
  "scales",
  "ggrepel"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE, type = "binary")
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(scales)
  library(ggrepel)
})

theme_icb <- function(base_size = 11, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 3, color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "black"),
      axis.title = element_text(face = "bold", size = base_size + 1, color = "black"),
      axis.text = element_text(size = base_size), ## 移除加粗，让X轴和Y轴的字体常规样式完全一致
      legend.title = element_text(face = "bold", size = base_size, color = "black"),
      legend.text = element_text(size = base_size - 1, color = "black"),
      strip.text = element_text(face = "bold", size = base_size, color = "black"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6)
    )
}

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
}

safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 300, device = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, plot = plot, width = width, height = height, dpi = dpi, bg = "white", device = device)
  message("Saved figure: ", file)
}

find_files_recursive <- function(root, pattern) {
  if (!dir.exists(root)) return(character(0))
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}

find_dirs_with_10x <- function(root) {
  if (!dir.exists(root)) return(character(0))
  all_dirs <- unique(dirname(find_files_recursive(root, "^matrix\\.mtx(\\.gz)?$")))
  all_dirs[file.exists(file.path(all_dirs, "barcodes.tsv")) |
             file.exists(file.path(all_dirs, "barcodes.tsv.gz")) |
             file.exists(file.path(all_dirs, "features.tsv")) |
             file.exists(file.path(all_dirs, "features.tsv.gz")) |
             file.exists(file.path(all_dirs, "genes.tsv")) |
             file.exists(file.path(all_dirs, "genes.tsv.gz"))]
}

clean_sample_name <- function(x) {
  x <- basename(x)
  x <- gsub("\\.h5$|\\.hdf5$|_filtered_feature_bc_matrix|_raw_feature_bc_matrix|filtered_feature_bc_matrix|raw_feature_bc_matrix", "", x, ignore.case = TRUE)
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x
}

get_assay_data_safe <- function(obj, assay = NULL, layer = "data", slot = NULL) {
  if (!is.null(slot)) layer <- slot
  if (is.null(assay)) DefaultAssay(obj)
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = layer)
  )
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  sx <- stats::sd(x, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(0, length(x)))
  (x - mx) / sx
}

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
count <- dplyr::count

heatmap_low  <- "#3B82F6"
heatmap_mid  <- "white"
heatmap_high <- "#EF4444"

major_celltype_colors_unified <- c(
  "Malignant" = "#3C5488",
  "Cycling malignant" = "#4DBBD5",
  "CAF/stromal-like cells" = "#E64B35",
  "Endothelial" = "#7E6148",
  "Myeloid cells" = "#00A087",
  "T/NK cells" = "#8491B4",
  "T/NK/Treg-like cells" = "#91D1C2",
  "B/Plasma cells" = "#F39B7F",
  "Cycling cells" = "#B09C85"
)

marker_sets <- list(
  Malignant = c("MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10", "S100B", "MIA", "TFAP2A"),
  T_NK = c("PTPRC", "CD3D", "CD3E", "CD2", "TRAC", "NKG7", "GNLY", "GZMB", "PRF1", "FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  B_Plasma = c("PTPRC", "MS4A1", "CD79A", "CD79B", "CD74", "MZB1", "JCHAIN", "IGHG1", "IGKC"),
  Myeloid = c("PTPRC", "LYZ", "LST1", "TYROBP", "CD14", "FCGR3A", "CD68", "C1QA", "C1QB", "S100A8", "S100A9", "IL1B"),
  Endothelial = c("PECAM1", "VWF", "KDR", "FLT1", "CLDN5", "RAMP2", "ESAM"),
  CAF_stromal_like = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "DCN", "LUM", "FAP", "ACTA2", "PDPN", "THY1", "VCAN", "TNC", "TGFB1", "FBN1", "MRC2", "ADAM12", "TIMP1", "INHBA"),
  Cycling = c("MKI67", "TOP2A", "UBE2C", "PCNA", "TYMS", "HMGB2", "STMN1", "CENPF")
)

marker_display <- c(
  Malignant = "Malignant",
  T_NK = "T/NK cells",
  B_Plasma = "B/Plasma cells",
  Myeloid = "Myeloid cells",
  Endothelial = "Endothelial",
  CAF_stromal_like = "CAF/stromal-like cells",
  Cycling = "Cycling"
)

major_celltype_levels <- c(
  "Malignant",
  "Cycling malignant",
  "CAF/stromal-like cells",
  "Endothelial",
  "Myeloid cells",
  "T/NK cells",
  "T/NK/Treg-like cells",
  "B/Plasma cells",
  "Cycling cells"
)

major_celltype_colors <- major_celltype_colors_unified[major_celltype_levels]

############################################################
## 3. Multi-path fallback Seurat object loading (命中即用)
############################################################

candidate_rds_paths <- c(
  file.path(project_dir, "results", "intermediate", "GSE244983", "GSE244983_seurat_major_annotated.rds"),
  file.path(project_dir, "data_processed", "GSE244983_seurat_major_annotated_LOCKED.rds"),
  file.path(project_dir, "data_processed", "GSE244983_seurat_major_annotated.rds"),
  file.path(project_dir, "data_processed", "GSE244983_seurat_QC_integrated_clustered.rds")
)

rds_path <- NA_character_
for (path in candidate_rds_paths) {
  norm_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (file.exists(norm_path)) {
    rds_path <- norm_path
    message("Successfully matched Seurat object at: ", rds_path)
    break
  }
}

obj <- NULL
if (!is.na(rds_path)) {
  obj <- readRDS(rds_path)
} else {
  all_rds <- list.files(project_dir, pattern = "GSE244983.*\\.rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(all_rds) > 0) {
    rds_path <- normalizePath(all_rds[1], winslash = "/", mustWork = FALSE)
    message("Fallback: Found Seurat object via recursive search: ", rds_path)
    obj <- readRDS(rds_path)
  }
}

if (is.null(obj)) {
  stop("Error: No valid GSE244983 Seurat object found across candidate paths or project directory.")
}

DefaultAssay(obj) <- ifelse("RNA" %in% names(obj@assays), "RNA", DefaultAssay(obj))

if (!"Sample" %in% colnames(obj@meta.data)) {
  if ("orig.ident" %in% colnames(obj@meta.data)) {
    obj$Sample <- obj$orig.ident
  } else {
    obj$Sample <- "GSE244983"
  }
}

############################################################
## 4. Metadata and annotation harmonization check
############################################################

if (!"MajorCellType" %in% colnames(obj@meta.data)) {
  message("MajorCellType metadata missing. Running automated annotation workflow...")
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = n_variable_features, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = n_pcs_to_compute, verbose = FALSE)
  
  dims_use <- seq_len(n_pcs_to_use)
  reduction_use <- "pca"
  if (use_harmony_if_available && requireNamespace("harmony", quietly = TRUE) && length(unique(obj$Sample)) > 1) {
    obj <- tryCatch(harmony::RunHarmony(obj, group.by.vars = "Sample", reduction.use = "pca", dims.use = dims_use, verbose = FALSE), error = function(e) obj)
    if ("harmony" %in% Reductions(obj)) reduction_use <- "harmony"
  }
  
  obj <- FindNeighbors(obj, reduction = reduction_use, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = reduction_use, dims = dims_use, verbose = FALSE)
  obj$SeuratCluster <- as.character(Idents(obj))
  
  expr_data <- get_assay_data_safe(obj, slot = "data")
  for (set_name in names(marker_sets)) {
    genes <- marker_sets[[set_name]]
    present <- intersect(genes, rownames(expr_data))
    score_z_col <- paste0(set_name, "_score_z")
    if (length(present) > 0) {
      score <- Matrix::colMeans(expr_data[present, , drop = FALSE])
      obj[[score_z_col]] <- zscore_safe(score)
    } else {
      obj[[score_z_col]] <- NA_real_
    }
  }
  
  score_cols_z <- paste0(names(marker_sets), "_score_z")
  cluster_scores <- obj@meta.data %>%
    mutate(SeuratCluster = as.character(SeuratCluster)) %>%
    group_by(SeuratCluster) %>%
    summarise(across(all_of(score_cols_z), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  auto_annotate_cluster <- function(row_df) {
    scores <- as.numeric(row_df[1, paste0(non_cycling_sets, "_score_z")])
    names(scores) <- non_cycling_sets
    scores[is.na(scores)] <- -Inf
    top_set <- names(which.max(scores))
    top_label <- marker_display[top_set]
    cycle_score <- as.numeric(row_df[1, "Cycling_score_z"])
    if (is.na(cycle_score)) cycle_score <- -Inf
    if (top_set == "Malignant" && cycle_score >= 0.75) return("Cycling malignant")
    if (top_set == "T_NK" && cycle_score >= 0.75) return("T/NK/Treg-like cells")
    if (cycle_score >= 1.00) return("Cycling cells")
    top_label
  }
  
  non_cycling_sets <- c("Malignant", "T_NK", "B_Plasma", "Myeloid", "Endothelial", "CAF_stromal_like")
  cluster_to_type <- sapply(seq_len(nrow(cluster_scores)), function(i) auto_annotate_cluster(cluster_scores[i, , drop = FALSE]))
  cluster_to_type <- setNames(cluster_to_type, cluster_scores$SeuratCluster)
  
  for (cl in names(manual_cluster_annotation)) {
    if (cl %in% names(cluster_to_type)) cluster_to_type[cl] <- manual_cluster_annotation[[cl]]
  }
  
  major_type_vec <- unname(cluster_to_type[as.character(obj$SeuratCluster)])
  obj@meta.data$MajorCellType <- factor(major_type_vec, levels = major_celltype_levels)
} else {
  obj@meta.data$MajorCellType <- factor(obj@meta.data$MajorCellType, levels = major_celltype_levels)
}

############################################################
## 5. Generate Supplementary Figure S6 (Cell-type composition by sample)
############################################################

cell_composition <- obj@meta.data %>%
  count(Sample, MajorCellType, name = "n_cells") %>%
  group_by(Sample) %>%
  mutate(proportion = n_cells / sum(n_cells)) %>%
  ungroup()

cell_composition$MajorCellType <- factor(cell_composition$MajorCellType, levels = major_celltype_levels)

p_comp <- ggplot(
  cell_composition,
  aes(x = Sample, y = proportion, fill = MajorCellType)
) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = major_celltype_colors, drop = FALSE) +
  labs(
    title = "Cell-type composition by sample in GSE244983",
    x = "Sample",
    y = "Cell-type proportion",
    fill = "Cell type"
  ) +
  theme_icb(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), ## 恢复常规非加粗字体，和Y轴完全一致
    legend.position = "right",
    legend.key.height = grid::unit(0.45, "cm")
  )

base_filename <- "Supplementary Figure S6. Cell-type composition by sample in GSE244983"
out_png <- file.path(out_fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(out_fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(out_fig_dir, paste0(base_filename, ".pdf"))

safe_ggsave(out_png, p_comp, width = 9.0, height = 6.2, dpi = 300)
safe_ggsave(out_jpg, p_comp, width = 9.0, height = 6.2, dpi = 300)
safe_ggsave(out_pdf, p_comp, width = 9.0, height = 6.2, dpi = 300, device = cairo_pdf)

cat("\nSupplementary Figure S6 generated successfully with normal axis text styling in PNG, JPG, and PDF formats (300 DPI)!\n")