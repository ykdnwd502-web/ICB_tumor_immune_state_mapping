############################################################
## 05_scRNA_malignant_subclustering.R
## Project: ICB_resistance_project
## Purpose:
##   Recluster malignant-lineage cells from the GSE244983 single-cell
##   melanoma ICB dataset and summarize malignant subcluster programs.
##
## Locked output logic:
##   - Input: <project>/results/intermediate/GSE244983/GSE244983_seurat_state_localized.rds
##   - Subset: Malignant + Cycling malignant cells
##   - Reprocess malignant-lineage cells
##   - Generate malignant subclusters named Mal_0, Mal_1, ...
##   - Save malignant-only Seurat object, tables, and Figure 4 panels (PDF/JPG)
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(Matrix)
})

## Avoid function masking
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
count <- dplyr::count
across <- dplyr::across
all_of <- dplyr::all_of
any_of <- dplyr::any_of
case_when <- dplyr::case_when
recode <- dplyr::recode
bind_rows <- dplyr::bind_rows

############################################################
## 0. Paths and output folders
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
if (!dir.exists(project_dir)) {
  stop("Project directory not found: ", project_dir)
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
message("Project directory: ", project_dir)

dir_data_processed <- file.path(project_dir, "data_processed")
dir_intermediate <- file.path(project_dir, "results", "intermediate", "GSE244983")
dir_tables <- file.path(project_dir, "results", "tables", "scRNA_malignant_subclustering")
dir_figs <- file.path(project_dir, "results", "figures", "scRNA_malignant_subclustering")
dir_logs <- file.path(project_dir, "logs")

dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_logs, recursive = TRUE, showWarnings = FALSE)

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = file, row.names = row.names, fileEncoding = "UTF-8")
  message("Saved table: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

safe_ggsave <- function(file, plot, width, height, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )
  message("Saved figure: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

############################################################
## 0.1 Typography and Aesthetics (Frozen Rules)
############################################################

theme_icb <- function(base_family = "Arial") {
  ggplot2::theme_bw(base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 15), # Left aligned panel titles
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 11),
      axis.title = ggplot2::element_text(face = "bold", size = 12.5),
      axis.text = ggplot2::element_text(color = "black", size = 11),
      legend.title = ggplot2::element_text(face = "bold", size = 11),
      legend.text = ggplot2::element_text(size = 10),
      strip.background = ggplot2::element_rect(fill = "#D9D9D9", color = "black", linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "bold", size = 12.5),
      panel.grid.major = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.60),
      legend.key = ggplot2::element_rect(fill = "white", color = NA)
    )
}

heatmap_low  <- "#3B82F6"
heatmap_mid  <- "white"
heatmap_high <- "#EF4444"

malignant_program_colors <- c(
  "Tumor_plasticity_dedifferentiation_score" = "#E64B35",
  "Stromal_ECM_remodeling_score" = "#E64B35",
  "Melanocytic_differentiation_score" = "#3C5488",
  "Antigen_presentation_score" = "#00A087",
  "IFN_response_score" = "#4DBBD5"
)

############################################################
## 1. Locate input Seurat object (Fixed Path)
############################################################

input_rds <- file.path(dir_intermediate, "GSE244983_seurat_state_localized.rds")

if (!file.exists(input_rds)) {
  stop("Input Seurat object not found at expected path: ", input_rds)
}

message("Loading input Seurat object: ", input_rds)
obj <- readRDS(input_rds)

############################################################
## 2. Helper functions compatible with Seurat v4/v5
############################################################

get_assay_matrix <- function(seu, assay = NULL, slot = "data") {
  assay <- assay %||% DefaultAssay(seu)
  aa <- seu[[assay]]
  
  if ("Assay5" %in% class(aa)) {
    layer_try <- slot
    if (slot == "data" && !"data" %in% Layers(aa)) {
      layer_try <- if ("counts" %in% Layers(aa)) "counts" else Layers(aa)[1]
    }
    return(GetAssayData(seu, assay = assay, layer = layer_try))
  } else {
    return(GetAssayData(seu, assay = assay, slot = slot))
  }
}

score_module_simple <- function(seu, genes, score_name, assay = NULL) {
  assay <- assay %||% DefaultAssay(seu)
  mat <- get_assay_matrix(seu, assay = assay, slot = "data")
  genes_present <- intersect(genes, rownames(mat))
  
  if (length(genes_present) < 2) {
    warning("Too few genes present for score: ", score_name)
    seu@meta.data[[score_name]] <- NA_real_
    return(seu)
  }
  
  score <- Matrix::colMeans(mat[genes_present, , drop = FALSE])
  score <- score[rownames(seu@meta.data)]
  seu@meta.data[[score_name]] <- as.numeric(score)
  return(seu)
}

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

############################################################
## 3. Identify malignant-lineage cells
############################################################

meta_cols <- colnames(obj@meta.data)

major_col_candidates <- c(
  "MajorCellType_LOCKED",
  "MajorCellType",
  "major_celltype",
  "celltype_major",
  "CellType",
  "cell_type",
  "annotation"
)

major_col <- major_col_candidates[major_col_candidates %in% meta_cols][1]
if (is.na(major_col) || length(major_col) == 0) {
  stop("No major cell-type annotation column found.")
}

message("Using major cell-type column: ", major_col)

major_values <- as.character(obj@meta.data[[major_col]])
malignant_labels <- c("Malignant", "Cycling malignant")

mal_cells <- rownames(obj@meta.data)[major_values %in% malignant_labels]
if (length(mal_cells) < 200) {
  stop("Too few malignant-lineage cells detected: ", length(mal_cells))
}

message("Malignant-lineage cells selected: ", length(mal_cells))

mal_obj <- subset(obj, cells = mal_cells)
DefaultAssay(mal_obj) <- DefaultAssay(obj)

try({
  assay_now <- DefaultAssay(mal_obj)
  if (!"Assay5" %in% class(mal_obj[[assay_now]])) {
    mal_obj[[assay_now]]@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
  }
}, silent = TRUE)

############################################################
## 4. Reprocess malignant-lineage cells
############################################################

set.seed(20260430)

if ("RNA" %in% Assays(mal_obj)) {
  DefaultAssay(mal_obj) <- "RNA"
}

rna_data <- tryCatch(get_assay_matrix(mal_obj, assay = DefaultAssay(mal_obj), slot = "data"), error = function(e) NULL)
if (is.null(rna_data) || nrow(rna_data) == 0) {
  mal_obj <- NormalizeData(mal_obj, verbose = FALSE)
}

mal_obj <- FindVariableFeatures(mal_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
mal_obj <- ScaleData(mal_obj, features = VariableFeatures(mal_obj), verbose = FALSE)
mal_obj <- RunPCA(mal_obj, features = VariableFeatures(mal_obj), npcs = 40, verbose = FALSE)

sample_col_candidates <- c("Sample", "sample", "patient", "Patient", "orig.ident")
sample_col <- sample_col_candidates[sample_col_candidates %in% colnames(mal_obj@meta.data)][1]
use_reduction <- "pca"

run_harmony_compatible <- function(seu, batch_col, dims_use = 1:30) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    message("Harmony package not installed; proceeding with PCA reduction.")
    return(list(object = seu, reduction = "pca", used = FALSE))
  }
  try_new <- tryCatch({
    harmony::RunHarmony(object = seu, group.by.vars = batch_col, reduction.use = "pca", dims.use = dims_use, reduction.save = "harmony", verbose = FALSE)
  }, error = function(e) e)
  
  if (!inherits(try_new, "error")) return(list(object = try_new, reduction = "harmony", used = TRUE))
  
  try_old <- tryCatch({
    harmony::RunHarmony(object = seu, group.by.vars = batch_col, dims.use = dims_use, reduction.save = "harmony", verbose = FALSE)
  }, error = function(e) e)
  
  if (!inherits(try_old, "error")) return(list(object = try_old, reduction = "harmony", used = TRUE))
  
  return(list(object = seu, reduction = "pca", used = FALSE))
}

if (!is.na(sample_col) && length(unique(mal_obj@meta.data[[sample_col]])) > 1) {
  harmony_res <- run_harmony_compatible(mal_obj, sample_col, dims_use = 1:30)
  mal_obj <- harmony_res$object
  use_reduction <- harmony_res$reduction
}

mal_obj <- FindNeighbors(mal_obj, reduction = use_reduction, dims = 1:30, verbose = FALSE)
mal_resolution <- 0.6
mal_obj <- FindClusters(mal_obj, resolution = mal_resolution, verbose = FALSE)
mal_obj <- RunUMAP(mal_obj, reduction = use_reduction, dims = 1:30, verbose = FALSE)

cluster_col <- "seurat_clusters"
mal_clusters_raw <- as.character(mal_obj@meta.data[[cluster_col]])
cluster_levels <- sort(as.numeric(unique(mal_clusters_raw)))
cluster_levels_chr <- as.character(cluster_levels)
mal_cluster_labels <- paste0("Mal_", cluster_levels_chr)
names(mal_cluster_labels) <- cluster_levels_chr

mal_obj@meta.data$MalignantSubcluster <- factor(
  mal_cluster_labels[as.character(mal_obj@meta.data[[cluster_col]])],
  levels = paste0("Mal_", cluster_levels_chr)
)

############################################################
## 5. Malignant marker/program gene sets
############################################################

malignant_marker_sets <- list(
  "Tumor plasticity/dedifferentiation" = c("AXL", "NGFR", "WNT5A", "EGFR", "SOX9", "ZEB1"),
  "Melanocytic differentiation" = c("MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"),
  "Antigen presentation" = c("B2M", "HLA-A", "HLA-B", "TAP1", "PSMB9", "STAT1", "IRF1"),
  "IFN response" = c("IFIT1", "ISG15", "IFIT3", "CXCL10", "GBP1", "OAS1")
)

score_gene_sets <- list(
  Tumor_plasticity_dedifferentiation_score = malignant_marker_sets[["Tumor plasticity/dedifferentiation"]],
  Melanocytic_differentiation_score = malignant_marker_sets[["Melanocytic differentiation"]],
  Antigen_presentation_score = malignant_marker_sets[["Antigen presentation"]],
  IFN_response_score = malignant_marker_sets[["IFN response"]]
)

for (nm in names(score_gene_sets)) {
  mal_obj <- score_module_simple(mal_obj, score_gene_sets[[nm]], nm, assay = DefaultAssay(mal_obj))
}

marker_order <- unique(unlist(malignant_marker_sets, use.names = FALSE))
marker_present <- marker_order[marker_order %in% rownames(mal_obj)]
program_score_cols <- names(score_gene_sets)

############################################################
## 6. Tables: subcluster composition and program summary
############################################################

program_summary <- mal_obj@meta.data %>%
  as.data.frame() %>%
  group_by(MalignantSubcluster) %>%
  summarise(
    N_cells = dplyr::n(),
    across(all_of(program_score_cols), list(mean = ~mean(.x, na.rm = TRUE), median = ~median(.x, na.rm = TRUE))),
    .groups = "drop"
  )

mean_cols <- paste0(program_score_cols, "_mean")
program_matrix <- program_summary %>% select(MalignantSubcluster, all_of(mean_cols))
mat <- as.matrix(program_matrix[, mean_cols, drop = FALSE])
rownames(mat) <- as.character(program_matrix$MalignantSubcluster)
mat_z <- scale(mat)
mat_z[is.na(mat_z)] <- 0

# Convert internal scoring column names to lowercase proper display labels
program_heatmap_df <- as.data.frame(mat_z) %>%
  tibble::rownames_to_column("MalignantSubcluster") %>%
  pivot_longer(-MalignantSubcluster, names_to = "Program", values_to = "Average_z_score") %>%
  mutate(
    Program = gsub("_score_mean$", "", Program),
    ProgramLabel = case_when(
      Program == "Tumor_plasticity_dedifferentiation" ~ "tumor plasticity/\ndedifferentiation",
      Program == "Melanocytic_differentiation" ~ "melanocytic\ndifferentiation",
      Program == "Antigen_presentation" ~ "antigen\npresentation",
      Program == "IFN_response" ~ "IFN\nresponse"
    ),
    ProgramLabel = factor(
      ProgramLabel,
      levels = c(
        "tumor plasticity/\ndedifferentiation",
        "melanocytic\ndifferentiation",
        "antigen\npresentation",
        "IFN\nresponse"
      )
    ),
    MalignantSubcluster = factor(MalignantSubcluster, levels = rownames(mat_z))
  )

############################################################
## 7. Figures
############################################################

mal_subcluster_palette <- c(
  "#EF4636", "#D99000", "#80B000", "#00A98B",
  "#00A6C8", "#55BFD1", "#3F5E9A", "#7C4DFF",
  "#E83E9F", "#F15A99", "#9E77ED", "#64748B",
  "#94A3B8", "#B45309", "#0F766E", "#BE123C",
  "#0369A1", "#9333EA", "#16A34A", "#DB2777"
)
n_sub <- length(levels(mal_obj$MalignantSubcluster))
mal_subcluster_colors <- mal_subcluster_palette[seq_len(min(n_sub, length(mal_subcluster_palette)))]
if (n_sub > length(mal_subcluster_palette)) mal_subcluster_colors <- grDevices::rainbow(n_sub)
names(mal_subcluster_colors) <- levels(mal_obj$MalignantSubcluster)

## Panel A
p_umap_mal <- DimPlot(
  mal_obj,
  reduction = "umap",
  group.by = "MalignantSubcluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.15,
  raster = FALSE,
  cols = mal_subcluster_colors
) +
  labs(
    title = "A  Malignant-cell subclusters",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Subcluster"
  ) +
  theme_icb(base_family = "Arial")

## Panel B
if (length(marker_present) >= 5) {
  marker_program_df <- bind_rows(lapply(names(malignant_marker_sets), function(pg) {
    data.frame(Program = pg, Feature = malignant_marker_sets[[pg]], stringsAsFactors = FALSE)
  })) %>%
    filter(Feature %in% marker_present) %>%
    mutate(
      FeatureIndex = match(Feature, marker_present),
      ProgramLabel = case_when(
        Program == "Tumor plasticity/dedifferentiation" ~ "tumor plasticity/\ndediff.",
        Program == "Melanocytic differentiation" ~ "melanocytic",
        Program == "Antigen presentation" ~ "antigen\npresentation",
        Program == "IFN response" ~ "IFN\nresponse",
        TRUE ~ Program
      )
    )
  
  program_bounds <- marker_program_df %>%
    group_by(Program, ProgramLabel) %>%
    summarise(
      xmin = min(FeatureIndex) - 0.5,
      xmax = max(FeatureIndex) + 0.5,
      xmid = mean(range(FeatureIndex)),
      .groups = "drop"
    ) %>%
    arrange(xmin)
  
  program_separators <- program_bounds %>%
    filter(xmax < length(marker_present) + 0.5) %>%
    transmute(xintercept = xmax)
  
  dot_base <- suppressWarnings(
    DotPlot(
      mal_obj,
      features = marker_present,
      group.by = "MalignantSubcluster",
      dot.scale = 5
    )
  )
  
  dot_df <- dot_base$data %>%
    mutate(
      Feature = as.character(features.plot),
      FeatureIndex = match(Feature, marker_present),
      Subcluster = as.character(id),
      SubclusterIndex = match(Subcluster, levels(mal_obj@meta.data$MalignantSubcluster))
    ) %>%
    filter(!is.na(FeatureIndex), !is.na(SubclusterIndex))
  
  n_subclusters <- length(levels(mal_obj@meta.data$MalignantSubcluster))
  
  p_dot_mal <- ggplot(dot_df, aes(x = FeatureIndex, y = SubclusterIndex)) +
    geom_vline(
      data = program_separators,
      aes(xintercept = xintercept),
      inherit.aes = FALSE,
      color = "#E5E7EB",
      linewidth = 0.55
    ) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled), alpha = 0.95) +
    geom_text(
      data = program_bounds,
      aes(x = xmid, y = n_subclusters + 0.92, label = ProgramLabel),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = 3.5, # approx 10pt
      lineheight = 0.90
    ) +
    scale_x_continuous(
      breaks = seq_along(marker_present),
      labels = marker_present,
      expand = expansion(add = c(0.45, 0.45))
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(mal_obj@meta.data$MalignantSubcluster)),
      labels = levels(mal_obj@meta.data$MalignantSubcluster),
      expand = expansion(add = c(0.45, 1.45))
    ) +
    scale_color_gradient2(
      low = heatmap_low, mid = heatmap_mid, high = heatmap_high, midpoint = 0,
      name = "Average\nExpression"
    ) +
    scale_size(
      range = c(0.15, 4.8),
      breaks = c(0, 25, 50, 75, 100), limits = c(0, 100),
      name = "Percent\nExpressed"
    ) +
    labs(
      title = "B  Marker programs across malignant subclusters",
      x = NULL, y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_icb(base_family = "Arial") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.5), # Standard rule size
      axis.text.y = element_text(size = 11),
      legend.position = "right",
      plot.margin = margin(t = 12, r = 8, b = 6, l = 6)
    )
} else {
  p_dot_mal <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Too few malignant markers available") +
    theme_void(base_family = "Arial")
}

## Panel C
p_heat_mal <- ggplot(program_heatmap_df, aes(x = ProgramLabel, y = MalignantSubcluster, fill = Average_z_score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Average_z_score)), size = 3.8, family = "Arial") +
  scale_fill_gradient2(
    low = heatmap_low, mid = heatmap_mid, high = heatmap_high, midpoint = 0,
    name = "Average\nz-score"
  ) +
  labs(
    title = "C  Malignant subcluster programs",
    x = NULL, y = NULL
  ) +
  theme_icb(base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 11),
    axis.text.y = element_text(size = 11)
  )

## Panel D
feature_titles <- c(
  Tumor_plasticity_dedifferentiation_score = "tumor plasticity/\ndedifferentiation",
  Melanocytic_differentiation_score = "melanocytic\ndifferentiation",
  Antigen_presentation_score = "antigen\npresentation",
  IFN_response_score = "IFN\nresponse"
)

feature_plots <- lapply(program_score_cols, function(sc) {
  FeaturePlot(
    mal_obj,
    features = sc,
    reduction = "umap",
    pt.size = 0.12,
    raster = FALSE,
    order = TRUE
  ) +
    scale_color_gradient(
      low = "#F7F7F7",
      high = malignant_program_colors[[sc]] %||% "#E64B35",
      name = "Score"
    ) +
    labs(
      title = feature_titles[[sc]],
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_icb(base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12.5), # Specific sub-title sizing
      legend.position = "right"
    )
})
names(feature_plots) <- program_score_cols

p_panel_d_title <- ggplot() +
  annotate(
    "text",
    x = 0, y = 0.5,
    label = "D  Malignant-cell distribution of tumor-intrinsic programs",
    family = "Arial",
    fontface = "bold",
    size = 5.0, # Approx 15pt mapping
    hjust = 0
  ) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void(base_family = "Arial")

p_fig3d_combined <- (feature_plots[[1]] | feature_plots[[2]]) / (feature_plots[[3]] | feature_plots[[4]])
p_fig3d_wrapped <- p_panel_d_title / p_fig3d_combined + plot_layout(heights = c(0.07, 0.93))

############################################################
## 8. Combine & Save Figure 4
############################################################

p_figure4 <- (
  (p_umap_mal | p_dot_mal) /
    (p_heat_mal | p_fig3d_wrapped)
) +
  plot_layout(heights = c(0.95, 1.05), widths = c(0.42, 0.58)) +
  plot_annotation(
    title = "Malignant-cell subclustering and tumor-intrinsic program heterogeneity in GSE244983",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15.5, family = "Arial"))
  )

output_file_name <- "Figure 5. Malignant-cell subclustering and tumor-intrinsic program heterogeneity in GSE244983"

# Save standardized outputs
safe_ggsave(file.path(dir_figs, paste0(output_file_name, ".pdf")), p_figure4, width = 16.8, height = 12.2)
safe_ggsave(file.path(dir_figs, paste0(output_file_name, ".jpg")), p_figure4, width = 16.8, height = 12.2, dpi = 300)

############################################################
## 9. Save Object and Finish
############################################################

saveRDS(mal_obj, file.path(dir_data_processed, "GSE244983_malignant_subclustered.rds"))
message("Saved RDS: ", file.path(dir_data_processed, "GSE244983_malignant_subclustered.rds"))

message("05_scRNA_malignant_subclustering.R finished successfully.")

