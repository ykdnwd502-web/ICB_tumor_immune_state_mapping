############################################################
## 06_scRNA_malignant_program_split.R
## Project: ICB_resistance_project
## Purpose:
##   Split malignant-cell tumor-intrinsic programs into five interpretable
##   modules: tumor plasticity/dedifferentiation, stromal/ECM remodeling,
##   melanocytic differentiation, antigen presentation, and IFN response.
##
## Locked output logic:
##   - Input: D:/ICB_resistance_project/data_processed/GSE244983_malignant_subclustered.rds
##   - Scores: simple mean expression of genes present in each module
##   - Outputs: tables, split-program figures matching rigorous style rules
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
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
## 0. Paths and helper functions
############################################################

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

project_dir <- Sys.getenv("ICB_PROJECT_DIR", unset = "")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
if (!dir.exists(project_dir)) {
  stop("Project directory not found: ", project_dir)
}

dir_data_processed <- file.path(project_dir, "data_processed")
dir_tables <- file.path(project_dir, "results", "tables", "scRNA_malignant_program_split")
dir_figs <- file.path(project_dir, "results", "figures", "scRNA_malignant_program_split")

dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = file, row.names = row.names, fileEncoding = "UTF-8")
  message("Saved table: ", file)
}

safe_ggsave <- function(file, plot, width, height, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved figure: ", file)
}

## 字体放大系数
FIG5_TEXT_SCALE <- 1.40
nt <- function(x) x * FIG5_TEXT_SCALE
FIG5_PROGRAM_HEADER_SIZE <- 3.55
FIG5_FEATURE_TITLE_SIZE <- 7.4

# 全局绘图规则：Arial 字体，左对齐标题，全包围黑框，无网格线
theme_icb <- function(base_size = 11, base_family = "Arial") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 2), # Panel 标题左对齐
      plot.subtitle = ggplot2::element_text(hjust = 0, size = base_size),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1),
      strip.background = ggplot2::element_rect(fill = "#D9D9D9", color = "black", linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      panel.grid.major = ggplot2::element_blank(), # 移除主网格线
      panel.grid.minor = ggplot2::element_blank(), # 移除副网格线
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.60), # 全包围边框
      axis.line = ggplot2::element_blank(),
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
## 1. Load malignant-only Seurat object
############################################################

# Clean-run dependency: Figure 5 writes the first candidate.
# The analysis-09 intermediate is retained only as a standalone fallback.
input_candidates <- c(
  file.path(project_dir, "data_processed", "GSE244983_malignant_subclustered.rds"),
  file.path(project_dir, "results", "intermediate", "GSE244983", "GSE244983_malignant_subclustered.rds")
)
input_candidates <- input_candidates[file.exists(input_candidates)]
if (length(input_candidates) == 0L) {
  stop("找不到 malignant-subclustered Seurat 对象；请先运行 Figure 5 或上游 malignant analysis。")
}
input_rds <- input_candidates[[1]]

message("Loading malignant-only Seurat object: ", input_rds)
mal_obj <- readRDS(input_rds)

if (!"MalignantSubcluster_LOCKED" %in% colnames(mal_obj@meta.data)) {
  if("MalignantSubcluster" %in% colnames(mal_obj@meta.data)) {
    mal_obj$MalignantSubcluster_LOCKED <- mal_obj$MalignantSubcluster
  } else {
    stop("找不到恶性亚群聚类列，请检查输入对象。")
  }
}

subcluster_levels <- unique(as.character(mal_obj@meta.data$MalignantSubcluster_LOCKED))
subcluster_levels <- subcluster_levels[order(as.integer(gsub("^Mal_", "", subcluster_levels)))]
mal_obj@meta.data$MalignantSubcluster_LOCKED <- factor(
  as.character(mal_obj@meta.data$MalignantSubcluster_LOCKED),
  levels = subcluster_levels
)

############################################################
## 2. Matrix and scoring functions
############################################################

get_assay_matrix <- function(seu, assay = NULL, slot = "data") {
  assay <- assay %||% DefaultAssay(seu)
  aa <- seu[[assay]]
  if ("Assay5" %in% class(aa)) {
    available_layers <- Layers(aa)
    layer_try <- slot
    if (!layer_try %in% available_layers) {
      layer_try <- if ("data" %in% available_layers) "data" else if ("counts" %in% available_layers) "counts" else available_layers[1]
    }
    return(GetAssayData(seu, assay = assay, layer = layer_try))
  } else {
    return(GetAssayData(seu, assay = assay, slot = slot))
  }
}

match_genes_to_object <- function(genes, object_genes) {
  genes <- unique(genes)
  direct <- intersect(genes, object_genes)
  missing <- setdiff(genes, direct)
  if (length(missing) > 0) {
    obj_upper <- toupper(object_genes)
    names(object_genes) <- obj_upper
    mapped <- object_genes[toupper(missing)]
    mapped <- mapped[!is.na(mapped)]
    direct <- unique(c(direct, unname(mapped)))
  }
  direct <- direct[direct %in% object_genes]
  unique(direct)
}

score_module_simple <- function(seu, genes, score_name, assay = NULL) {
  assay <- assay %||% DefaultAssay(seu)
  mat <- get_assay_matrix(seu, assay = assay, slot = "data")
  genes_present <- match_genes_to_object(genes, rownames(mat))
  
  if (length(genes_present) < 2) {
    seu@meta.data[[score_name]] <- NA_real_
    seu@meta.data[[paste0(score_name, "_z")]] <- NA_real_
    return(seu)
  }
  
  score <- Matrix::colMeans(mat[genes_present, , drop = FALSE])
  score <- score[rownames(seu@meta.data)]
  score <- as.numeric(score)
  score_z <- as.numeric(scale(score))
  
  seu@meta.data[[score_name]] <- score
  seu@meta.data[[paste0(score_name, "_z")]] <- score_z
  return(seu)
}

############################################################
## 3. Define split tumor-intrinsic programs
############################################################

program_gene_sets <- list(
  Tumor_plasticity_dedifferentiation = c("AXL", "NGFR", "WNT5A", "EGFR", "SOX9", "ZEB1", "FN1", "VIM"),
  Stromal_ECM_remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3",
                             "DCN", "LUM", "FAP", "ACTA2", "PDPN", "THY1", "VCAN", "TNC", "TGFBI",
                             "TIMP1", "INHBA", "MRC2", "ADAM12"),
  Melanocytic_differentiation = c("MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"),
  Antigen_presentation = c("B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9", "STAT1", "IRF1"),
  IFN_response = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "CXCL10", "GBP1", "OAS1", "IFI6", "MX1")
)

program_order <- names(program_gene_sets)
program_display <- c(
  Tumor_plasticity_dedifferentiation = "Tumor plasticity/\nDedifferentiation",
  Stromal_ECM_remodeling = "Stromal/ECM\nremodeling",
  Melanocytic_differentiation = "Melanocytic\nDifferentiation",
  Antigen_presentation = "Antigen\nPresentation",
  IFN_response = "IFN\nResponse"
)
program_short <- c(
  Tumor_plasticity_dedifferentiation = "Plasticity/\ndediff.",
  Stromal_ECM_remodeling = "Stromal/ECM",
  Melanocytic_differentiation = "Melanocytic",
  Antigen_presentation = "Antigen\npresentation",
  IFN_response = "IFN\nresponse"
)

############################################################
## 4. Compute split-program module scores
############################################################

score_cols <- paste0(program_order, "_score")
for (i in seq_along(program_order)) {
  pg <- program_order[i]
  score_nm <- score_cols[i]
  mal_obj <- score_module_simple(mal_obj, program_gene_sets[[pg]], score_nm, assay = DefaultAssay(mal_obj))
}
score_z_cols <- paste0(score_cols, "_z")

############################################################
## 5. Subcluster-level summary and heatmap matrix
############################################################

program_summary <- mal_obj@meta.data %>%
  as.data.frame() %>%
  group_by(MalignantSubcluster_LOCKED) %>%
  summarise(
    across(all_of(score_cols), list(mean = ~mean(.x, na.rm = TRUE))),
    .groups = "drop"
  )

mean_cols <- paste0(score_cols, "_mean")
program_matrix <- program_summary %>% select(MalignantSubcluster_LOCKED, all_of(mean_cols))
mat <- as.matrix(program_matrix[, mean_cols, drop = FALSE])
rownames(mat) <- as.character(program_matrix$MalignantSubcluster_LOCKED)

mat_z <- scale(mat)
mat_z[is.na(mat_z)] <- 0
colnames(mat_z) <- program_order

program_heatmap_df <- as.data.frame(mat_z) %>%
  tibble::rownames_to_column("MalignantSubcluster") %>%
  pivot_longer(-MalignantSubcluster, names_to = "Program", values_to = "Average_z_score") %>%
  mutate(
    Program = factor(Program, levels = program_order),
    ProgramLabel = factor(program_display[as.character(Program)], levels = program_display[program_order]),
    MalignantSubcluster = factor(MalignantSubcluster, levels = rev(subcluster_levels))
  )

############################################################
## 6. Correlation among split programs
############################################################

cor_mat <- suppressWarnings(cor(
  mal_obj@meta.data[, score_cols, drop = FALSE],
  method = "spearman",
  use = "pairwise.complete.obs"
))
rownames(cor_mat) <- program_order
colnames(cor_mat) <- program_order

cor_df <- as.data.frame(cor_mat) %>%
  tibble::rownames_to_column("Program_row") %>%
  pivot_longer(-Program_row, names_to = "Program_col", values_to = "Spearman_rho") %>%
  mutate(
    Program_row = factor(Program_row, levels = rev(program_order)),
    Program_col = factor(Program_col, levels = program_order),
    Program_row_label = factor(program_display[as.character(Program_row)], levels = rev(program_display[program_order])),
    Program_col_label = factor(program_display[as.character(Program_col)], levels = program_display[program_order])
  )

############################################################
## 7. UMAP data
############################################################

umap_df <- as.data.frame(Embeddings(mal_obj, "umap"))
colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
umap_df$Cell <- rownames(umap_df)
umap_df <- bind_cols(umap_df, mal_obj@meta.data[rownames(umap_df), , drop = FALSE])

subcluster_centers <- umap_df %>%
  group_by(MalignantSubcluster_LOCKED) %>%
  summarise(
    umap_1 = median(umap_1, na.rm = TRUE),
    umap_2 = median(umap_2, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 8. Generate Figures
############################################################

subcluster_palette <- c(
  "Mal_0" = "#EF4438", "Mal_1" = "#D99500", "Mal_2" = "#7CAE00", "Mal_3" = "#00A88A",
  "Mal_4" = "#00A6C8", "Mal_5" = "#4DBBD5", "Mal_6" = "#3C5488", "Mal_7" = "#6C43FF",
  "Mal_8" = "#E83E8C", "Mal_9" = "#F45B9A", "Mal_10" = "#8B5CF6", "Mal_11" = "#64748B"
)
subcluster_palette <- subcluster_palette[subcluster_levels]

## Panel A
p_umap_subclusters <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = MalignantSubcluster_LOCKED)) +
  geom_point(size = 0.18, alpha = 0.9) +
  geom_text(
    data = subcluster_centers,
    aes(x = umap_1, y = umap_2, label = MalignantSubcluster_LOCKED),
    inherit.aes = FALSE,
    family = "Arial",
    size = nt(3.3),
    color = "black"
  ) +
  scale_color_manual(values = subcluster_palette, drop = FALSE, name = "Subcluster") +
  labs(title = "A  Malignant-cell subclusters", x = "UMAP 1", y = "UMAP 2") +
  coord_equal() +
  theme_icb(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 15, hjust = 0) # 左对齐
  )

## Panel B
marker_gene_sets_reduced <- list(
  Tumor_plasticity_dedifferentiation = c("AXL", "NGFR", "WNT5A", "EGFR", "SOX9", "ZEB1", "FN1", "VIM"),
  Stromal_ECM_remodeling = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "THY1", "VCAN", "TGFBI"),
  Melanocytic_differentiation = c("MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"),
  Antigen_presentation = c("B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2"),
  IFN_response = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "CXCL10", "GBP1")
)

make_split_marker_dotplot <- function(seu, gene_sets, title) {
  marker_order <- unique(unlist(gene_sets[program_order], use.names = FALSE))
  marker_present <- intersect(marker_order, rownames(seu))
  
  marker_program_df <- bind_rows(lapply(program_order, function(pg) {
    data.frame(Program = pg, Feature = intersect(gene_sets[[pg]], marker_present), stringsAsFactors = FALSE)
  })) %>%
    mutate(
      FeatureIndex = match(Feature, marker_present),
      Program = factor(Program, levels = program_order),
      ProgramLabel = program_short[as.character(Program)]
    )
  
  program_bounds <- marker_program_df %>%
    group_by(Program, ProgramLabel) %>%
    summarise(
      xmin = min(FeatureIndex) - 0.5,
      xmax = max(FeatureIndex) + 0.5,
      xmid = mean(range(FeatureIndex)),
      .groups = "drop"
    )
  
  program_separators <- program_bounds %>%
    filter(xmax < length(marker_present) + 0.5) %>%
    transmute(xintercept = xmax)
  
  dot_base <- suppressWarnings(DotPlot(seu, features = marker_present, group.by = "MalignantSubcluster_LOCKED", dot.scale = 5.2))
  dot_df <- dot_base$data %>%
    mutate(
      FeatureIndex = match(as.character(features.plot), marker_present),
      SubclusterIndex = match(as.character(id), subcluster_levels)
    )
  
  p <- ggplot(dot_df, aes(x = FeatureIndex, y = SubclusterIndex)) +
    geom_vline(data = program_separators, aes(xintercept = xintercept), inherit.aes = FALSE, color = "#E5E7EB", linewidth = 0.55) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled), alpha = 0.95) +
    geom_text(
      data = program_bounds,
      aes(x = xmid, y = length(subcluster_levels) + 0.95, label = ProgramLabel),
      inherit.aes = FALSE, family = "Arial", fontface = "bold", size = FIG5_PROGRAM_HEADER_SIZE, lineheight = 0.88
    ) +
    scale_x_continuous(breaks = seq_along(marker_present), labels = marker_present, expand = expansion(add = c(0.45, 0.45))) +
    scale_y_continuous(breaks = seq_along(subcluster_levels), labels = subcluster_levels, expand = expansion(add = c(0.45, 2.05))) +
    scale_color_gradient2(low = heatmap_low, mid = heatmap_mid, high = heatmap_high, midpoint = 0, name = "Average\nExpression") +
    scale_size(range = c(0.15, 4.8), breaks = c(0, 25, 50, 75, 100), limits = c(0, 100), name = "Percent\nExpressed") +
    labs(title = title, x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_icb(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = nt(8.2)),
      plot.title = element_text(face = "bold", hjust = 0, size = 15), # 左对齐
      plot.margin = margin(t = 18, r = 8, b = 6, l = 6)
    )
  p
}

p_marker_dot_reduced <- make_split_marker_dotplot(mal_obj, marker_gene_sets_reduced, "B  Split marker programs across malignant subclusters")

## Panel C
p_heat_split <- ggplot(program_heatmap_df, aes(x = ProgramLabel, y = MalignantSubcluster, fill = Average_z_score)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", Average_z_score)), family = "Arial", size = nt(2.65)) +
  scale_fill_gradient2(low = heatmap_low, mid = heatmap_mid, high = heatmap_high, midpoint = 0, name = "Average\nz-score") +
  labs(title = "C  Split tumor-intrinsic programs across malignant subclusters", x = NULL, y = NULL) +
  theme_icb(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = nt(9)),
    plot.title = element_text(face = "bold", hjust = 0, size = 15) # 左对齐
  )

## Panel D
plot_feature_umap <- function(df, score_col, title, q_low = 0.01, q_high = 0.99) {
  vals <- df[[score_col]]
  lo <- quantile(vals, q_low, na.rm = TRUE)
  hi <- quantile(vals, q_high, na.rm = TRUE)
  df$Score_plot <- pmin(pmax(vals, lo), hi)
  
  ggplot(df, aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = Score_plot), size = 0.16, alpha = 0.9) +
    scale_color_gradient(low = "#F7F7F7", high = malignant_program_colors[[score_col]] %||% "#E64B35", name = "Score") +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    coord_equal() +
    theme_icb(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = FIG5_FEATURE_TITLE_SIZE, lineheight = 0.82),
      plot.margin = margin(t = 2, r = 3, b = 3, l = 3)
    )
}

feature_titles <- c(
  Tumor_plasticity_dedifferentiation_score = "Tumor plasticity/\nDedifferentiation",
  Stromal_ECM_remodeling_score = "Stromal/ECM\nremodeling",
  Melanocytic_differentiation_score = "Melanocytic\nDifferentiation",
  Antigen_presentation_score = "Antigen\nPresentation",
  IFN_response_score = "IFN\nResponse"
)

feature_plots <- lapply(names(feature_titles), function(sc) { plot_feature_umap(umap_df, sc, feature_titles[[sc]]) })
p_feature_grid <- wrap_plots(feature_plots, ncol = 3) +
  plot_annotation(title = "D  Malignant-cell distribution of tumor-intrinsic programs") &
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 15, family = "Arial")) # 左对齐嵌套标签

## Panel E
p_cor_split <- ggplot(cor_df, aes(x = Program_col_label, y = Program_row_label, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Spearman_rho)), family = "Arial", size = nt(3.4)) +
  scale_fill_gradient2(low = heatmap_low, mid = heatmap_mid, high = heatmap_high, midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
  labs(title = "E  Correlation among split malignant-cell programs", x = NULL, y = NULL) +
  theme_icb(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = nt(9.5)),
    plot.title = element_text(face = "bold", hjust = 0, size = 15, margin = margin(b = 6)) # 左对齐
  )

############################################################
## 9. Assemble and Save
############################################################

p_figure6_top <- p_umap_subclusters | p_marker_dot_reduced
p_figure6_mid <- p_heat_split | p_feature_grid

p_figure6 <- p_figure6_top /
  p_figure6_mid /
  p_cor_split +
  plot_layout(heights = c(1.0, 1.18, 0.82), widths = c(1.0, 1.15)) +
  plot_annotation(
    title = "Split tumor-intrinsic programs within malignant-cell subclusters in GSE244983",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18, family = "Arial")) # 大图名居中
  )

output_file_name <- "Figure 6. Separable tumor-intrinsic programs within malignant-cell subclusters"

safe_ggsave(file.path(dir_figs, paste0(output_file_name, ".pdf")), p_figure6, width = 17.2, height = 16.2)
safe_ggsave(file.path(dir_figs, paste0(output_file_name, ".jpg")), p_figure6, width = 17.2, height = 16.2, dpi = 300)

message("Figure 6 Plotting finished successfully.")