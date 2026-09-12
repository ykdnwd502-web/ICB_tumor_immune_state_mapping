############################################################
## Figure_3_scRNA_localization.R
##
## Purpose:
## 1) Load the major-cell-type annotated GSE244983 Seurat object
## 2) Reuse the bulk-derived GMT signatures and final-state mapping
## 3) Score tumor–immune state programs at single-cell level
## 4) Localize these programs across major cell types
## 5) Generate Figure 3 panels with formal typography and layout.
##
## Plot font:
##   sans
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

############################################################
## 0. Project directory and folders
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}

message("Project directory: ", project_dir)

data_processed_dir <- file.path(project_dir, "data_processed")
table_dir          <- file.path(project_dir, "results/tables")
figure_dir         <- file.path(project_dir, "results/figures")
out_table_dir      <- file.path(project_dir, "results/tables/scRNA_major_program_localization")
out_fig_dir        <- file.path(project_dir, "results/figures/scRNA_major_program_localization")
log_dir            <- file.path(project_dir, "logs")

dir.create(data_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. User-adjustable parameters
############################################################

seurat_candidates <- c(
  file.path(project_dir, "results/intermediate/GSE244983/GSE244983_seurat_major_annotated.rds"),
  file.path(project_dir, "results/intermediate/GSE244983/GSE244983_seurat_state_localized.rds")
)

state_mapping_candidates <- c(
  file.path(project_dir, "results/tables/bulk_discovery/GSE244982_final_state_signature_mapping.csv"),
  file.path(project_dir, "results/tables/GSE244982_final_state_signature_mapping_FIXED.csv")
)

gmt_candidates <- c(
  file.path(project_dir, "results/tables/ICBcomb_final_input_gene_sets_CLEAN/02_final_ICBcomb_gene_sets/final_ICBcomb_gene_sets_CLEAN.gmt"),
  file.path(project_dir, "results/tables/ICBcomb_final_input_gene_sets_v2/final_ICBcomb_gene_sets_v2.gmt"),
  file.path(project_dir, "results/tables/ICBcomb_gene_sets/GSE244982_state_gene_sets_for_ICBcomb.gmt")
)

score_suffix <- "_score"

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

heatmap_abs_limit <- NULL

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
  "scales"
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
})

theme_icb <- function(base_size = 10, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0, size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6)
    )
}

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
  message("Saved table: ", file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
  message("Saved RDS: ", file)
}

safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved figure: ", file)
}

find_first_existing <- function(paths) {
  paths <- unique(paths)
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) return(NA_character_)
  paths[1]
}

find_files_recursive <- function(root, pattern) {
  if (!dir.exists(root)) return(character(0))
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}

clean_gene_symbol <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\s+", "", x)
  x <- gsub("^HUMAN_", "", x)
  x <- gsub("^MOUSE_", "", x)
  x
}

normalize_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

zvec <- function(x) {
  x <- as.numeric(x)
  sx <- stats::sd(x, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(0, length(x)))
  (x - mx) / sx
}

get_assay_data_safe <- function(obj, assay = NULL, layer = "data", slot = NULL) {
  if (!is.null(slot)) layer <- slot
  if (is.null(assay)) assay <- DefaultAssay(obj)
  out <- tryCatch(
    { Seurat::GetAssayData(obj, assay = assay, layer = layer) },
    error = function(e_layer) {
      tryCatch(
        { Seurat::GetAssayData(obj, assay = assay, slot = layer) },
        error = function(e_slot) {
          if (layer != "counts") {
            return(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
          }
          stop(e_slot)
        }
      )
    }
  )
  out
}

read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file, warn = FALSE)
  lines <- lines[nchar(lines) > 0]
  gene_sets <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 3) next
    set_name <- parts[1]
    genes <- unique(clean_gene_symbol(parts[-c(1, 2)]))
    genes <- genes[!is.na(genes) & genes != ""]
    gene_sets[[set_name]] <- genes
  }
  gene_sets
}

match_gene_set_name <- function(target, available) {
  if (target %in% available) return(target)
  target_norm <- normalize_name(target)
  available_norm <- normalize_name(available)
  idx <- which(available_norm == target_norm)
  if (length(idx) >= 1) return(available[idx[1]])
  idx2 <- grep(target_norm, available_norm, fixed = TRUE)
  if (length(idx2) >= 1) return(available[idx2[1]])
  idx3 <- grep(available_norm, target_norm, fixed = TRUE)
  if (length(idx3) >= 1) return(available[idx3[1]])
  NA_character_
}

############################################################
## 2A. Unified plotting style and palettes
############################################################

if (!exists("heatmap_low"))  heatmap_low  <- "#3B82F6"
if (!exists("heatmap_mid"))  heatmap_mid  <- "white"
if (!exists("heatmap_high")) heatmap_high <- "#EF4444"

major_celltype_colors_unified <- c(
  "Malignant" = "#3C5488",
  "Cycling malignant" = "#4DBBD5",
  "CAF/stromal-like cells" = "#E64B35",
  "Endothelial" = "#00A087",
  "Myeloid cells" = "#00A087",
  "T/NK cells" = "#8491B4",
  "T/NK/Treg-like cells" = "#91D1C2",
  "B/Plasma cells" = "#F39B7F",
  "Cycling cells" = "#7E6148"
)

############################################################
## 3. Load annotated Seurat object
############################################################

seurat_file <- find_first_existing(seurat_candidates)
if (is.na(seurat_file)) {
  stop("No annotated Seurat object found.")
}

message("Loading Seurat object: ", seurat_file)
obj <- readRDS(seurat_file)
DefaultAssay(obj) <- ifelse("RNA" %in% names(obj@assays), "RNA", DefaultAssay(obj))

# Standardize column name for plotting
if ("MajorCellType_LOCKED" %in% colnames(obj@meta.data)) {
  obj$Plot_MajorCellType <- factor(as.character(obj@meta.data$MajorCellType_LOCKED), levels = major_celltype_levels)
} else {
  obj$Plot_MajorCellType <- factor(as.character(obj@meta.data$MajorCellType), levels = major_celltype_levels)
}

# Invert the Y-axis of the UMAP embedding to match manuscript orientation
if ("umap" %in% names(obj@reductions)) {
  obj@reductions$umap@cell.embeddings[, 2] <- -obj@reductions$umap@cell.embeddings[, 2]
}

############################################################
## 4. Load bulk-derived gene sets and final-state mapping
############################################################

gmt_recursive <- find_files_recursive(project_dir, ".*gene.*set.*\\.gmt$")
gmt_file <- find_first_existing(unique(c(gmt_candidates, gmt_recursive)))

message("Using GMT file: ", gmt_file)
gene_sets_all <- read_gmt(gmt_file)

mapping_recursive <- find_files_recursive(project_dir, "GSE244982_final_state_signature_mapping.*\\.csv$")
state_mapping_file <- find_first_existing(unique(c(state_mapping_candidates, mapping_recursive)))

message("Using final-state mapping table: ", state_mapping_file)
state_mapping <- read.csv(state_mapping_file, check.names = FALSE, stringsAsFactors = FALSE)

clean_to_legacy_state <- c(
  "Immune_defective_Cold" = "Immune_defective_cold",
  "Myeloid_Treg_Immunosuppressive" = "Myeloid_Treg_immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "Tumor_dedifferentiated",
  "Melanocytic_Differentiation" = "Melanocytic_differentiated"
)
if ("FinalState" %in% colnames(state_mapping)) {
  state_mapping$FinalState <- ifelse(
    state_mapping$FinalState %in% names(clean_to_legacy_state),
    unname(clean_to_legacy_state[state_mapping$FinalState]),
    state_mapping$FinalState
  )
}

available_gs <- names(gene_sets_all)
state_mapping$MatchedGeneSet <- vapply(
  state_mapping$RawSignatureColumn,
  match_gene_set_name,
  FUN.VALUE = character(1),
  available = available_gs
)

state_mapping_use <- state_mapping %>%
  filter(!is.na(MatchedGeneSet)) %>%
  distinct(FinalState, RawSignatureColumn, DirectionForFinalState, MatchedGeneSet)

state_order <- c(
  "Immune_defective_cold",
  "Myeloid_Treg_immunosuppressive",
  "Tumor_dedifferentiated",
  "Melanocytic_differentiated"
)

state_label_map <- c(
  "Immune_defective_cold" = "immune-defective/\ncold",
  "Myeloid_Treg_immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiated" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_differentiated" = "melanocytic\ndifferentiation"
)

state_label_one_line <- c(
  "Immune_defective_cold" = "immune-defective/cold",
  "Myeloid_Treg_immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiated" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_differentiated" = "melanocytic differentiation"
)

score_col_map <- paste0(state_order, score_suffix)
names(score_col_map) <- state_order

state_colors <- c(
  "Immune_defective_cold" = "#4DBBD5",
  "Myeloid_Treg_immunosuppressive" = "#00A087",
  "Tumor_dedifferentiated" = "#E64B35",
  "Melanocytic_differentiated" = "#3C5488"
)

major_celltype_colors <- major_celltype_colors_unified[major_celltype_levels]
major_celltype_colors <- major_celltype_colors[!is.na(names(major_celltype_colors))]

major_marker_sets <- list(
  "T/NK cells" = c("PTPRC", "CD3D", "CD3E", "CD2", "NKG7", "GNLY", "GZMB", "PRF1", "FOXP3"),
  "B/Plasma cells" = c("MS4A1", "CD79A", "MZB1"),
  "Myeloid cells" = c("LYZ", "S100A8", "S100A9", "CD68", "C1QA"),
  "CAF/stromal-like cells" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "PDPN", "THY1", "VCAN", "TNC"),
  "Endothelial" = c("PECAM1", "VWF", "KDR", "FLT1"),
  "Malignant" = c("MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"),
  "Cycling" = c("MKI67", "TOP2A")
)
major_marker_order <- unique(unlist(major_marker_sets, use.names = FALSE))

############################################################
## 5. Score raw signatures at single-cell level
############################################################

expr_data <- get_assay_data_safe(obj, layer = "data")
expr_genes_original <- rownames(expr_data)
expr_genes_clean <- clean_gene_symbol(expr_genes_original)
names(expr_genes_clean) <- expr_genes_original

score_one_gene_set <- function(genes_clean, set_name) {
  matched_rows <- expr_genes_original[expr_genes_clean %in% genes_clean]
  matched_rows <- unique(matched_rows)
  if (length(matched_rows) == 0) return(rep(NA_real_, ncol(expr_data)))
  small_mat <- as.matrix(expr_data[matched_rows, , drop = FALSE])
  small_z <- t(scale(t(small_mat)))
  small_z[is.na(small_z)] <- 0
  colMeans(small_z, na.rm = TRUE)
}

raw_score_list <- list()
for (gs in unique(state_mapping_use$MatchedGeneSet)) {
  genes <- gene_sets_all[[gs]]
  genes_clean <- unique(clean_gene_symbol(genes))
  raw_score_list[[gs]] <- score_one_gene_set(genes_clean, gs)
}

raw_score_mat <- as.data.frame(do.call(cbind, raw_score_list), check.names = FALSE)
rownames(raw_score_mat) <- colnames(obj)
raw_score_mat$Cell <- rownames(raw_score_mat)

############################################################
## 6. Build final state scores at single-cell level
############################################################

raw_score_numeric <- raw_score_mat[, setdiff(colnames(raw_score_mat), "Cell"), drop = FALSE]
final_score_df <- data.frame(Cell = colnames(obj), stringsAsFactors = FALSE)
rownames(final_score_df) <- colnames(obj)

for (st in state_order) {
  map_st <- state_mapping_use %>% filter(FinalState == st)
  matched_cols <- intersect(map_st$MatchedGeneSet, colnames(raw_score_numeric))
  if (length(matched_cols) == 0) {
    final_score_df[[score_col_map[[st]]]] <- NA_real_
    next
  }
  raw_state_score <- rowMeans(raw_score_numeric[, matched_cols, drop = FALSE], na.rm = TRUE)
  direction <- unique(map_st$DirectionForFinalState)
  if (any(tolower(direction[!is.na(direction)]) == "inverse")) {
    raw_state_score <- -raw_state_score
  }
  final_score_df[[score_col_map[[st]]]] <- zvec(raw_state_score)
}

score_cols <- unname(score_col_map)
for (cc in score_cols) {
  if (cc %in% colnames(final_score_df)) {
    obj@meta.data[[cc]] <- final_score_df[colnames(obj), cc]
  }
}

############################################################
## 7. Summarize state scores by major cell type
############################################################

meta_score <- obj@meta.data %>%
  mutate(
    Cell = rownames(obj@meta.data),
    Plot_MajorCellType = obj$Plot_MajorCellType
  ) %>%
  select(Cell, Plot_MajorCellType, all_of(score_cols))

score_long <- meta_score %>%
  pivot_longer(cols = all_of(score_cols), names_to = "ScoreColumn", values_to = "Score") %>%
  mutate(
    State = names(score_col_map)[match(ScoreColumn, score_col_map)],
    State = factor(State, levels = state_order),
    StateLabel = factor(state_label_map[as.character(State)], levels = state_label_map[state_order]),
    StateLabelOneLine = factor(state_label_one_line[as.character(State)], levels = state_label_one_line[state_order])
  )

summary_by_celltype <- score_long %>%
  group_by(Plot_MajorCellType, State, StateLabel, StateLabelOneLine) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  arrange(State, Plot_MajorCellType)

############################################################
## 8. Figures: Figure 3 panels A-D
############################################################

## Panel A: UMAP
p_fig3a <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "Plot_MajorCellType",
  label = FALSE,
  pt.size = 0.12,
  raster = FALSE,
  cols = major_celltype_colors
) +
  labs(
    title = "A  Major cell-type annotation",
    x = "UMAP 1", 
    y = "UMAP 2", 
    color = "Cell type"
  ) +
  theme_icb(base_size = 10) +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.42, "cm"),
    plot.title = element_text(face = "bold", hjust = 0, size = 12) 
  )

## Panel B: DotPlot
marker_features_present <- major_marker_order[major_marker_order %in% rownames(obj)]

p_fig3b <- suppressWarnings(
  DotPlot(
    obj,
    features = marker_features_present, 
    group.by = "Plot_MajorCellType",
    dot.scale = 5
  )
) +
  scale_color_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    name = "Average\nExpression"
  ) +
  labs(
    title = "B  Representative markers across major cell types",
    x = NULL,
    y = NULL,
    size = "Percent\nExpressed"
  ) +
  theme_icb(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.2), 
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.3),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    plot.title = element_text(face = "bold", hjust = 0, size = 12) 
  )

## Panel C: Heatmap
heat_df <- summary_by_celltype %>%
  mutate(
    Plot_MajorCellType = factor(as.character(Plot_MajorCellType), levels = rev(major_celltype_levels)),
    StateLabel = factor(state_label_map[as.character(State)], levels = state_label_map[state_order])
  ) %>%
  filter(!is.na(Plot_MajorCellType), !is.na(StateLabel))

if (is.null(heatmap_abs_limit)) {
  max_abs <- as.numeric(stats::quantile(abs(heat_df$mean_score), probs = 0.98, na.rm = TRUE))
  max_abs <- max(1.5, max_abs)
} else {
  max_abs <- heatmap_abs_limit
}

p_heat <- ggplot(
  heat_df,
  aes(x = StateLabel, y = Plot_MajorCellType, fill = mean_score)
) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", mean_score)), size = 3.5, family = "sans") +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    oob = scales::squish,
    name = "Average\nz-score\n(clipped)"
  ) +
  labs(
    title = "C  Mean tumor–immune state scores across major cell types", 
    subtitle = paste0("Color scale clipped at ±", sprintf("%.2f", max_abs), "; numbers show original average z-scores"),
    x = NULL,
    y = NULL
  ) +
  theme_icb(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11), 
    axis.text.y = element_text(size = 11),                       
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0, size = 12)
  )

## Panel D: UMAP Feature Plots
umap_df <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
umap_df$Cell <- rownames(umap_df)
umap_score_df <- cbind(
  umap_df,
  obj@meta.data[rownames(umap_df), c("Plot_MajorCellType", score_cols), drop = FALSE]
)

plot_score_umap <- function(df, score_col, title_text, high_col = "#4DBBD5") {
  dd <- df %>% arrange(.data[[score_col]])
  ggplot(dd, aes(x = umap_1, y = umap_2, color = .data[[score_col]])) +
    geom_point(size = 0.08, alpha = 0.85) +
    scale_color_gradient2(
      low = heatmap_low, mid = heatmap_mid, high = high_col, midpoint = 0, name = "Score"
    ) +
    labs(
      title = title_text, 
      x = "UMAP 1", 
      y = "UMAP 2"  
    ) +
    theme_icb(base_size = 10) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 10)
    )
}

umap_plots <- list()
for (st in state_order) {
  umap_plots[[st]] <- plot_score_umap(
    umap_score_df,
    score_col_map[[st]],
    state_label_one_line[[st]],
    high_col = state_colors[[st]]
  )
}

p_panel_d_title <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 0.5,
    label = "D  Single-cell distribution of tumor–immune state scores",
    fontface = "bold",
    family = "sans",
    size = 4.25,
    hjust = 0 
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void(base_family = "sans")

p_umap_grid_combined <- (
  umap_plots[["Immune_defective_cold"]] | umap_plots[["Myeloid_Treg_immunosuppressive"]]
) / (
  umap_plots[["Tumor_dedifferentiated"]] | umap_plots[["Melanocytic_differentiated"]]
)

p_umap_scores_combined <- p_panel_d_title / p_umap_grid_combined +
  plot_layout(heights = c(0.07, 0.93))

############################################################
## 8.5 Combine Figure 3 A-D
############################################################

plot_formal_title <- "Single-cell localization of predefined tumor-immune state scores in GSE2449833"

p_figure3 <- (
  (p_fig3a | p_fig3b) /
    (p_heat | p_umap_scores_combined)
) +
  plot_layout(
    heights = c(0.95, 1.05),
    widths = c(0.38, 0.62)
  ) +
  plot_annotation(
    title = plot_formal_title, 
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16, family = "sans")
    )
  )

############################################################
## 9. Save Figures
############################################################

output_file_name <- "Figure 3. Single-cell localization of predefined tumor-immune state scores in GSE244983"

safe_ggsave(
  file.path(out_fig_dir, paste0(output_file_name, ".pdf")),
  p_figure3,
  width = 17,
  height = 12
)

safe_ggsave(
  file.path(out_fig_dir, paste0(output_file_name, ".jpg")),
  p_figure3,
  width = 17,
  height = 12,
  dpi = 300
)

message("Figure_3_scRNA_localization.R finished successfully.")

