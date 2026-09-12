############################################################
## STYLE-ONLY v1.1 S20-REFERENCE LOCK
## Figure typography is harmonized to the accepted S20 visual ratio.
## Long subtitles are locally reduced and/or wrapped when needed.
## Scientific inputs, thresholds, statistics, state definitions,
## analytical logic, numerical outputs, and figure canvas sizes are unchanged.
############################################################

## PUBLIC FIGURE STYLE CONTRACT
##   Internal analytical state IDs are unchanged.
##   Display labels are lowercase-first:
##     immune-defective/cold
##     myeloid–Treg immunosuppressive
##     tumor-dedifferentiation/stromal-remodeling
##     melanocytic differentiation
##   Fixed state palette:
##     #4DBBD5 / #00A087 / #E64B35 / #3C5488
##   Manuscript-facing figures use sans, base size 10, white background,
##   light major grid, no minor grid, and black panel border where applicable.
##   No numerical/statistical definition is changed by the style patch.

############################################################
## REPORTING-SUBSET LOCK v1.0
## Figure 8C / Figure 8D display the top 20 ranked positive
## candidate genes, matching the frozen reporting subset.
## IMPORTANT: this is a reporting/display subset only.
## The >=5% detection filter, ranking rule, candidate table,
## dual-high definition, statistics, and all scientific logic
## are unchanged.
############################################################

# ============================================================
# Step 12. Dual-high niche mechanism nomination (Figure 8B-D)
# Clean-room reproducibility version
# v5 layout-fix: refine standalone Figure 8B/C/D only;
# 8B compact forest-style panel; 8C unified low-saturation gene-group colors;
# 8D four-group candidate-gene dotplot; no Figure 8B-D combined export
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir,winslash="/",mustWork=TRUE)
manual_spatial_rds <- file.path(
  project_dir,
  "data_processed",
  "spatial_10x_human_melanoma_IF_FFPE_with_ICBcomb_state_scores.rds"
)

step11_or_table <- file.path(
  project_dir,
  "results", "tables", "spatial_melanoma_validation",
  "Step11_spatial_dual_high_spot_level_OR_top25_pooled.csv"
)

out_table_dir <- file.path(project_dir, "results", "tables", "spatial_melanoma_validation")
out_fig_dir   <- file.path(project_dir, "results", "figures", "spatial_melanoma_validation")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

high_cutoff_quantile <- 0.75
min_detect_frac <- 0.05
n_top_genes_plot <- 20
n_top_genes_heat <- 20

state_cols <- c(
  immune = "Immune_defective_Cold",
  myeloid_treg = "Myeloid_Treg_Immunosuppressive",
  dediff_stromal = "Tumor_dedifferentiation_Stromal_remodeling",
  melanocytic = "Melanocytic_Differentiation"
)

dual_group_colors <- c(
  "Neither high" = "#CFCFCF",
  "Myeloid–Treg high only" = "#10A38B",
  "Dediff/Stromal high only" = "#EB4D34",
  "Both high" = "#8E44AD"
)

gene_group_colors <- c(
  "ECM/CAF-like remodeling" = "#D28A4A",
  "B-cell/plasma-cell/humoral immune" = "#7E65A5",
  "Complement/inflammatory remodeling" = "#7A9B35",
  "Antigen presentation/MHC" = "#4AA3BF",
  "Myeloid/macrophage inflammatory remodeling" = "#D65C6A",
  "Dedifferentiation/invasive phenotype" = "#6F9B86",
  "Other" = "#8A8A8A"
)

message("Project directory: ", project_dir)
message("Output table directory: ", out_table_dir)
message("Output figure directory: ", out_fig_dir)

save_csv_clean <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
  message("Saved table: ", path)
}

save_plot_both <- function(p, stem, width = 10, height = 8, dpi = 300) {
  png_file <- file.path(out_fig_dir, paste0(stem, ".png"))
  pdf_file <- file.path(out_fig_dir, paste0(stem, ".pdf"))
  ggplot2::ggsave(png_file, plot = p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(pdf_file, plot = p, width = width, height = height, device = cairo_pdf)
  message("Saved figure: ", png_file)
  message("Saved figure: ", pdf_file)
}

theme_clean_pub <- function(base_size = 12.5, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2.5, hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "black"),
      axis.title = element_text(face = "bold", size = base_size, color = "black"),
      axis.title.y = element_text(face = "bold", size = base_size - 1, color = "black"),
      axis.title.x = element_text(face = "bold", size = base_size, color = "black"),
      axis.text = element_text(size = base_size - 1.5, color = "black"),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.60),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1.5),
      strip.background = element_rect(fill = "grey92", color = "grey50", linewidth = 0.30),
      strip.text = element_text(face = "bold", color = "black", size = base_size - 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

first_existing <- function(paths) {
  idx <- which(file.exists(paths))
  if (length(idx) == 0) return(NA_character_)
  paths[idx[1]]
}

calc_or_from_binary <- function(a, b) {
  tab <- table(A = a, B = b)
  ft <- fisher.test(tab)
  list(
    table = as.data.frame.matrix(tab),
    estimate = unname(ft$estimate),
    conf.low = unname(ft$conf.int[1]),
    conf.high = unname(ft$conf.int[2]),
    p.value = unname(ft$p.value)
  )
}

safe_get_expr_matrix <- function(obj, assay = NULL) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)

  assay_obj <- obj[[assay]]

  normalize_mat_names <- function(mat) {
    if (is.null(mat)) return(NULL)
    if (nrow(mat) == 0 || ncol(mat) == 0) return(NULL)

    if (is.null(rownames(mat))) {
      feature_names <- tryCatch(rownames(assay_obj), error = function(e) NULL)
      if (is.null(feature_names)) feature_names <- tryCatch(Seurat::Features(obj, assay = assay), error = function(e) NULL)
      if (!is.null(feature_names) && length(feature_names) == nrow(mat)) rownames(mat) <- feature_names
    }
    if (is.null(colnames(mat))) {
      cell_names <- tryCatch(colnames(obj), error = function(e) NULL)
      if (!is.null(cell_names) && length(cell_names) == ncol(mat)) colnames(mat) <- cell_names
    }
    mat
  }

  get_layers_safe <- function() {
    lyr <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) NULL)
    if (is.null(lyr)) lyr <- tryCatch(Seurat::Layers(assay_obj), error = function(e) NULL)
    if (is.null(lyr)) lyr <- character(0)
    as.character(lyr)
  }

  get_layer_safe <- function(layer_name) {
    mat <- tryCatch(SeuratObject::LayerData(assay_obj, layer = layer_name), error = function(e) NULL)
    if (is.null(mat)) {
      mat <- tryCatch(Seurat::LayerData(obj, assay = assay, layer = layer_name), error = function(e) NULL)
    }
    normalize_mat_names(mat)
  }

  layers <- get_layers_safe()
  preferred_layers <- unique(c(
    "data",
    grep("^data", layers, value = TRUE),
    "counts",
    grep("^counts", layers, value = TRUE),
    layers
  ))
  preferred_layers <- preferred_layers[preferred_layers %in% layers]

  for (lyr in preferred_layers) {
    mat_try <- get_layer_safe(lyr)
    if (!is.null(mat_try) && nrow(mat_try) > 0 && ncol(mat_try) > 0) {
      if (grepl("^counts", lyr, ignore.case = TRUE)) {
        message("Using expression layer '", lyr, "' from assay '", assay, "' with log1p(counts) fallback.")
        mat_try <- log1p(mat_try)
      } else {
        message("Using expression layer '", lyr, "' from assay '", assay, "'.")
      }
      return(mat_try)
    }
  }

  mat_try <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "data"), error = function(e) NULL)
  mat_try <- normalize_mat_names(mat_try)
  if (!is.null(mat_try) && nrow(mat_try) > 0 && ncol(mat_try) > 0) {
    message("Using GetAssayData(slot='data') from assay '", assay, "'.")
    return(mat_try)
  }

  mat_try <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  mat_try <- normalize_mat_names(mat_try)
  if (!is.null(mat_try) && nrow(mat_try) > 0 && ncol(mat_try) > 0) {
    message("Using GetAssayData(slot='counts') from assay '", assay, "' with log1p(counts) fallback.")
    return(log1p(mat_try))
  }

  for (sl in c("data", "counts")) {
    if (sl %in% methods::slotNames(assay_obj)) {
      mat_try <- tryCatch(methods::slot(assay_obj, sl), error = function(e) NULL)
      mat_try <- normalize_mat_names(mat_try)
      if (!is.null(mat_try) && nrow(mat_try) > 0 && ncol(mat_try) > 0) {
        if (sl == "counts") {
          message("Using direct assay slot 'counts' with log1p(counts) fallback.")
          return(log1p(mat_try))
        } else {
          message("Using direct assay slot 'data'.")
          return(mat_try)
        }
      }
    }
  }

  stop(
    "Could not retrieve expression matrix from assay: ", assay,
    ". Available layers: ", paste(layers, collapse = ", "),
    ". This usually means the Spatial assay has no accessible data/counts layer; please inspect Layers(obj[['", assay, "']])."
  )
}

get_assay_layer_audit <- function(obj, assay = NULL) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  if (!assay %in% names(obj@assays)) {
    return(data.frame(Assay = assay, Layer = NA_character_, Status = "ASSAY_NOT_FOUND"))
  }
  assay_obj <- obj[[assay]]
  layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) NULL)
  if (is.null(layers)) layers <- tryCatch(Seurat::Layers(assay_obj), error = function(e) NULL)
  if (is.null(layers) || length(layers) == 0) {
    return(data.frame(Assay = assay, Layer = NA_character_, Status = "NO_LAYER_LISTED"))
  }
  do.call(rbind, lapply(as.character(layers), function(lyr) {
    mat <- tryCatch(SeuratObject::LayerData(assay_obj, layer = lyr), error = function(e) NULL)
    if (is.null(mat)) mat <- tryCatch(Seurat::LayerData(obj, assay = assay, layer = lyr), error = function(e) NULL)
    data.frame(
      Assay = assay,
      Layer = lyr,
      n_features = ifelse(is.null(mat), NA_integer_, nrow(mat)),
      n_spots = ifelse(is.null(mat), NA_integer_, ncol(mat)),
      HasRowNames = ifelse(is.null(mat), NA, !is.null(rownames(mat))),
      HasColNames = ifelse(is.null(mat), NA, !is.null(colnames(mat))),
      Status = ifelse(is.null(mat), "UNREADABLE", "READABLE"),
      stringsAsFactors = FALSE
    )
  }))
}

classify_gene_group <- function(gene) {
  gene <- toupper(gene)
  if (gene %in% c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","COL6A3","DCN","LUM","COMP","SPARC","THBS1","FN1","MMP2","IGFBP4")) return("ECM/CAF-like remodeling")
  if (gene %in% c("IGKC","IGLC1","IGHG1","IGHA1","JCHAIN","MZB1","CD79A","MS4A1")) return("B-cell/plasma-cell/humoral immune")
  if (gene %in% c("C1R","C1S","C3","CXCL9","CXCL10","CXCL11","IFITM3")) return("Complement/inflammatory remodeling")
  if (gene %in% c("CD74","B2M","HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","TAP1","TAP2","PSMB8","PSMB9")) return("Antigen presentation/MHC")
  if (gene %in% c("LYZ","S100A6","FCER1G","TYROBP","AIF1","LST1","CTSB","CTSD")) return("Myeloid/macrophage inflammatory remodeling")
  if (gene %in% c("KRT14","VIM","ZEB1","SOX9","AXL","NGFR","EGFR")) return("Dedifferentiation/invasive phenotype")
  "Other"
}

spatial_rds <- manual_spatial_rds
if (!file.exists(spatial_rds)) {
  candidates <- list.files(
    file.path(project_dir, "data_processed"),
    pattern = "spatial.*ICBcomb.*state.*scores.*\\.rds$|spatial.*state.*scores.*\\.rds$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  spatial_rds <- first_existing(candidates)
}
if (is.na(spatial_rds) || !file.exists(spatial_rds)) {
  stop("Could not locate the state-scored Visium RDS. Please check manual_spatial_rds.")
}

obj <- readRDS(spatial_rds)
if (!inherits(obj, "Seurat")) stop("Input object is not a Seurat object: ", class(obj)[1])
message("Loaded spatial object: ", spatial_rds)

meta <- obj@meta.data
meta$Spot <- rownames(meta)
missing_states <- setdiff(unname(state_cols), colnames(meta))
if (length(missing_states) > 0) {
  stop("Missing expected state columns in metadata: ", paste(missing_states, collapse = ", "))
}

object_audit <- data.frame(
  InputRDS = spatial_rds,
  ObjectClass = class(obj)[1],
  Assay = Seurat::DefaultAssay(obj),
  n_spots = ncol(obj),
  n_genes = nrow(obj),
  stringsAsFactors = FALSE
)
save_csv_clean(object_audit, file.path(out_table_dir, "Step12_dual_high_input_object_audit.csv"))

assay_layer_audit <- get_assay_layer_audit(obj, assay = Seurat::DefaultAssay(obj))
save_csv_clean(assay_layer_audit, file.path(out_table_dir, "Step12_expression_assay_layer_audit.csv"))

my_cut <- as.numeric(stats::quantile(meta[[state_cols[["myeloid_treg"]]]], probs = high_cutoff_quantile, na.rm = TRUE))
dd_cut <- as.numeric(stats::quantile(meta[[state_cols[["dediff_stromal"]]]], probs = high_cutoff_quantile, na.rm = TRUE))

meta$myeloid_treg_high <- meta[[state_cols[["myeloid_treg"]]]] >= my_cut
meta$dediff_stromal_high <- meta[[state_cols[["dediff_stromal"]]]] >= dd_cut
meta$dual_high_group <- dplyr::case_when(
  meta$myeloid_treg_high & meta$dediff_stromal_high ~ "Both high",
  meta$myeloid_treg_high & !meta$dediff_stromal_high ~ "Myeloid–Treg high only",
  !meta$myeloid_treg_high & meta$dediff_stromal_high ~ "Dediff/Stromal high only",
  TRUE ~ "Neither high"
)
meta$dual_high_group <- factor(meta$dual_high_group, levels = c("Neither high", "Myeloid–Treg high only", "Dediff/Stromal high only", "Both high"))

save_csv_clean(
  data.frame(CutoffMode = "pooled_top25", HighQuantile = high_cutoff_quantile, MyeloidTregCutoff = my_cut, DediffStromalCutoff = dd_cut),
  file.path(out_table_dir, "Step12_dual_high_cutoff_audit.csv")
)

spot_group_counts <- meta %>% dplyr::count(dual_high_group, name = "n_spots") %>% dplyr::rename(Group = dual_high_group)
save_csv_clean(spot_group_counts, file.path(out_table_dir, "Step12_dual_high_group_counts.csv"))

or_res <- calc_or_from_binary(meta$myeloid_treg_high, meta$dediff_stromal_high)
or_table <- data.frame(
  AnalysisUnit = "spot",
  CutoffMode = "pooled_top25",
  MyeloidTregHigh_n = sum(meta$myeloid_treg_high, na.rm = TRUE),
  DediffStromalHigh_n = sum(meta$dediff_stromal_high, na.rm = TRUE),
  BothHigh_n = sum(meta$myeloid_treg_high & meta$dediff_stromal_high, na.rm = TRUE),
  NeitherHigh_n = sum(!meta$myeloid_treg_high & !meta$dediff_stromal_high, na.rm = TRUE),
  OR = or_res$estimate,
  CI_lower = or_res$conf.low,
  CI_upper = or_res$conf.high,
  P_value = or_res$p.value,
  stringsAsFactors = FALSE
)
save_csv_clean(or_table, file.path(out_table_dir, "Step12_dual_high_spot_level_OR_top25_pooled.csv"))
save_csv_clean(or_res$table %>% tibble::rownames_to_column("MyeloidTregHigh"), file.path(out_table_dir, "Step12_dual_high_2x2_table_top25_pooled.csv"))

if (file.exists(step11_or_table)) {
  step11_or <- tryCatch(read.csv(step11_or_table, check.names = FALSE), error = function(e) NULL)
  if (!is.null(step11_or)) {
    step12_or_audit <- data.frame(
      Step11File = step11_or_table,
      Step12OR = or_table$OR,
      Step12Lower = or_table$CI_lower,
      Step12Upper = or_table$CI_upper,
      Step12P = or_table$P_value,
      stringsAsFactors = FALSE
    )
    save_csv_clean(step12_or_audit, file.path(out_table_dir, "Step12_dual_high_OR_step11_comparison_audit.csv"))
  }
}

# ------------------------------
# Figure 8B (standalone)
# ------------------------------
fig8b_data <- data.frame(
  Label = "Dual-high\nco-occurrence",
  OR = or_table$OR,
  CI_lower = or_table$CI_lower,
  CI_upper = or_table$CI_upper,
  P_value = or_table$P_value,
  stringsAsFactors = FALSE
)

fig8b_annot <- sprintf(
  "OR = %.2f\n95%% CI %.2f–%.2f\nP = %.2e",
  fig8b_data$OR,
  fig8b_data$CI_lower,
  fig8b_data$CI_upper,
  fig8b_data$P_value
)

x_max_8b <- max(11.5, fig8b_data$CI_upper * 1.23)
annot_y_8b <- min(fig8b_data$CI_upper * 1.05, x_max_8b * 0.86)

p8b <- ggplot(fig8b_data, aes(x = Label, y = OR)) +
  geom_hline(
    yintercept = 1,
    linetype = 2,
    color = "grey45",
    linewidth = 0.7
  ) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.09,
    linewidth = 1.0,
    color = "black"
  ) +
  geom_point(
    size = 4.8,
    color = dual_group_colors[["Both high"]]
  ) +
  annotate(
    "text",
    x = 1,
    y = annot_y_8b,
    label = fig8b_annot,
    hjust = 0,
    vjust = 0.5,
    size = 4.7
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(0, x_max_8b),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Figure 8B. Spot-level dual-high enrichment",
    subtitle = "myeloid–Treg high and tumor-dedifferentiation/stromal-remodeling high;\nunit of analysis: Visium spot",
    x = NULL,
    y = "Odds ratio"
  ) +
  theme_clean_pub(base_size = 12.5) +
  theme(
    panel.border = element_blank(),
    axis.line.x = element_line(color = "grey25", linewidth = 0.5),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12.5),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5, lineheight = 1.0),
    plot.margin = margin(8, 70, 8, 8)
  )

save_plot_both(
  p8b,
  "Figure8B_dual_high_enrichment_OR",
  width = 8.6,
  height = 4.2
)

# Expression-based candidate nomination
# ------------------------------
expr_mat <- safe_get_expr_matrix(obj, assay = Seurat::DefaultAssay(obj))
common_spots <- intersect(colnames(expr_mat), meta$Spot)
if (length(common_spots) == 0) stop("No overlapping spots between expression matrix columns and metadata rownames.")
expr_mat <- expr_mat[, common_spots, drop = FALSE]
meta_sub <- meta[match(common_spots, meta$Spot), , drop = FALSE]

keep <- meta_sub$dual_high_group %in% c("Neither high", "Both high")
expr_sub <- expr_mat[, keep, drop = FALSE]
meta_grp <- meta_sub[keep, , drop = FALSE]

both_idx <- which(meta_grp$dual_high_group == "Both high")
neither_idx <- which(meta_grp$dual_high_group == "Neither high")

frac_both <- Matrix::rowMeans(expr_sub[, both_idx, drop = FALSE] > 0)
frac_neither <- Matrix::rowMeans(expr_sub[, neither_idx, drop = FALSE] > 0)
mean_both <- Matrix::rowMeans(expr_sub[, both_idx, drop = FALSE])
mean_neither <- Matrix::rowMeans(expr_sub[, neither_idx, drop = FALSE])

candidate_tbl <- data.frame(
  Gene = rownames(expr_sub),
  Mean_BothHigh = as.numeric(mean_both),
  Mean_NeitherHigh = as.numeric(mean_neither),
  Diff_Both_minus_Neither = as.numeric(mean_both - mean_neither),
  FracExpr_BothHigh = as.numeric(frac_both),
  FracExpr_NeitherHigh = as.numeric(frac_neither),
  MaxFracExpr = pmax(as.numeric(frac_both), as.numeric(frac_neither)),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(GeneGroup = vapply(Gene, classify_gene_group, FUN.VALUE = character(1))) %>%
  dplyr::arrange(dplyr::desc(Diff_Both_minus_Neither), dplyr::desc(MaxFracExpr), Gene)

save_csv_clean(candidate_tbl, file.path(out_table_dir, "Step12_dual_high_candidate_gene_screen_Both_vs_Neither.csv"))
candidate_tbl_filt <- candidate_tbl %>% dplyr::filter(MaxFracExpr >= min_detect_frac)
save_csv_clean(candidate_tbl_filt, file.path(out_table_dir, "Step12_dual_high_candidate_gene_screen_Both_vs_Neither_filtered.csv"))

# ------------------------------
# Figure 8C (standalone)
# ------------------------------
top_genes <- candidate_tbl_filt %>%
  dplyr::filter(Diff_Both_minus_Neither > 0) %>%
  dplyr::slice_head(n = n_top_genes_plot)
save_csv_clean(top_genes, file.path(out_table_dir, "Step12_dual_high_top_candidate_genes_Figure8C.csv"))

top_genes_plot <- top_genes %>%
  dplyr::mutate(
    GeneGroup = factor(GeneGroup, levels = names(gene_group_colors)),
    Gene = forcats::fct_reorder(Gene, Diff_Both_minus_Neither)
  )

used_gene_groups_8c <- names(gene_group_colors)[
  names(gene_group_colors) %in% as.character(unique(top_genes_plot$GeneGroup))
]

p8c <- ggplot(
  top_genes_plot,
  aes(x = Gene, y = Diff_Both_minus_Neither, fill = GeneGroup)
) +
  geom_col(
    width = 0.82,
    color = "grey25",
    linewidth = 0.25
  ) +
  coord_flip() +
  scale_fill_manual(
    values = gene_group_colors[used_gene_groups_8c],
    breaks = used_gene_groups_8c,
    drop = TRUE
  ) +
  labs(
    title = "Figure 8C. Candidate mechanism genes in the dual-high niche",
    subtitle = "Ranked by mean log-normalized expression difference: Both high minus Neither high",
    x = NULL,
    y = "Mean expression difference (Both high - Neither high)",
    fill = "Gene group"
  ) +
  theme_clean_pub(base_size = 12.5) +
  theme(
    axis.text.y = element_text(size = 10.75),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12.5),
    legend.position = "right",
    legend.title = element_text(size = 11.5),
    legend.text = element_text(size = 10.75),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width = unit(0.45, "cm"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 8)
  )

save_plot_both(
  p8c,
  "Figure8C_candidate_mechanism_genes_Both_vs_Neither",
  width = 12.0,
  height = 8.2
)

# ------------------------------
# Figure 8D (standalone; four-group candidate-gene dotplot)
# ------------------------------
# Rationale:
# The previous two-group row-z heatmap was mathematically correct but visually redundant:
# for genes selected as Both-high > Neither-high, a two-column row z-score almost always
# becomes blue in Neither-high and red in Both-high. This v6 panel summarizes nominated
# genes across all four spot-level groups and uses compact wrapped group labels.

fig8d_genes <- candidate_tbl_filt %>%
  dplyr::filter(Diff_Both_minus_Neither > 0) %>%
  dplyr::slice_head(n = n_top_genes_heat) %>%
  dplyr::mutate(
    GeneGroup = factor(GeneGroup, levels = names(gene_group_colors))
  ) %>%
  dplyr::arrange(GeneGroup, dplyr::desc(Diff_Both_minus_Neither))

fig8d_genes$y_pos <- rev(seq_len(nrow(fig8d_genes)))

group_levels_8d <- c(
  "Neither high",
  "Myeloid–Treg high only",
  "Dediff/Stromal high only",
  "Both high"
)

group_labels_8d_x <- c(
  "Neither\nhigh",
  "Myeloid–Treg\nhigh only",
  "Dediff/Stromal\nhigh only",
  "Both\nhigh"
)

gene_group_label_short <- c(
  "ECM/CAF-like remodeling" = "ECM/CAF-like\nremodeling",
  "B-cell/plasma-cell/humoral immune" = "B-cell/plasma-cell\nhumoral immune",
  "Complement/inflammatory remodeling" = "Complement/\ninflammatory remodeling",
  "Antigen presentation/MHC" = "Antigen presentation/\nMHC",
  "Myeloid/macrophage inflammatory remodeling" = "Myeloid/macrophage\ninflammatory remodeling",
  "Dedifferentiation/invasive phenotype" = "Dedifferentiation/\ninvasive phenotype",
  "Other" = "Other"
)

# Calculate mean log-normalized expression and percent-expressing spots in each group.
# Dot color: row-scaled mean expression across the four groups for each gene.
# Dot size: percent of spots with expression > 0 in each group.
fig8d_summary_list <- lapply(fig8d_genes$Gene, function(g) {
  vals_gene <- as.numeric(expr_mat[g, meta_sub$Spot, drop = TRUE])
  names(vals_gene) <- meta_sub$Spot

  do.call(rbind, lapply(group_levels_8d, function(grp) {
    idx <- which(as.character(meta_sub$dual_high_group) == grp)
    vals <- vals_gene[idx]
    data.frame(
      Gene = g,
      Group = grp,
      MeanExpr = mean(vals, na.rm = TRUE),
      PctExpr = mean(vals > 0, na.rm = TRUE) * 100,
      stringsAsFactors = FALSE
    )
  }))
})

fig8d_summary <- do.call(rbind, fig8d_summary_list) %>%
  dplyr::left_join(
    fig8d_genes[, c("Gene", "GeneGroup", "Diff_Both_minus_Neither", "y_pos")],
    by = "Gene"
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(
    RowScaledMeanExpr = {
      s <- stats::sd(MeanExpr, na.rm = TRUE)
      if (is.na(s) || s == 0) rep(0, dplyr::n()) else (MeanExpr - mean(MeanExpr, na.rm = TRUE)) / s
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Group = factor(Group, levels = group_levels_8d),
    GroupX = as.numeric(Group),
    GeneGroup = factor(GeneGroup, levels = names(gene_group_colors))
  )

save_csv_clean(
  fig8d_summary,
  file.path(out_table_dir, "Step12_dual_high_Figure8D_four_group_candidate_gene_expression_summary.csv")
)

fig8d_gene_label_df <- fig8d_genes[, c("Gene", "GeneGroup", "y_pos")]

fig8d_group_label_df <- fig8d_genes %>%
  dplyr::group_by(GeneGroup) %>%
  dplyr::summarise(
    y_mid = mean(y_pos),
    y_min = min(y_pos) - 0.5,
    y_max = max(y_pos) + 0.5,
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(GeneGroup)) %>%
  dplyr::mutate(
    GeneGroupShort = gene_group_label_short[as.character(GeneGroup)]
  )

p8d <- ggplot() +
  geom_segment(
    data = fig8d_gene_label_df,
    aes(
      x = 0.52,
      xend = 0.52,
      y = y_pos - 0.42,
      yend = y_pos + 0.42,
      color = GeneGroup
    ),
    linewidth = 4.0,
    lineend = "butt",
    show.legend = FALSE
  ) +
  geom_segment(
    data = fig8d_group_label_df,
    aes(
      x = 0.72,
      xend = 4.35,
      y = y_min,
      yend = y_min
    ),
    color = "#DADADA",
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = fig8d_summary,
    aes(x = GroupX, y = y_pos, size = PctExpr, fill = RowScaledMeanExpr),
    shape = 21,
    color = "grey30",
    stroke = 0.28,
    alpha = 0.95
  ) +
  geom_text(
    data = fig8d_group_label_df,
    aes(
      x = 4.62,
      y = y_mid,
      label = GeneGroupShort,
      color = GeneGroup
    ),
    hjust = 0,
    vjust = 0.5,
    size = 3.45,
    lineheight = 0.93,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    low = "#6F90DA",
    mid = "#F7F7F7",
    high = "#EE705D",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    oob = scales::squish
  ) +
  scale_color_manual(values = gene_group_colors, drop = TRUE) +
  scale_size_continuous(
    range = c(0.8, 5.2),
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100)
  ) +
  scale_x_continuous(
    limits = c(0.38, 6.75),
    breaks = seq_along(group_levels_8d),
    labels = group_labels_8d_x,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = fig8d_genes$y_pos,
    labels = fig8d_genes$Gene,
    expand = expansion(add = c(0.35, 0.35))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Figure 8D. Four-group expression pattern of nominated genes",
    subtitle = "Dot color shows row-scaled mean log expression; dot size shows percent-expressing spots",
    x = NULL,
    y = NULL,
    fill = "Row-scaled\nmean expression",
    size = "Percent\nexpressed"
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = unit(3.2, "cm"), barwidth = unit(0.45, "cm")),
    size = guide_legend(order = 2, override.aes = list(fill = "grey75", color = "grey30"))
  ) +
  theme_clean_pub(base_size = 12.5) +
  theme(
    panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.35),
    panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.28),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10.75),
    axis.text.x = element_text(size = 10.75),
    legend.position = "right",
    legend.title = element_text(size = 10.75, face = "bold"),
    legend.text = element_text(size = 10.5),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.75),
    plot.margin = margin(8, 70, 8, 8)
  )

save_plot_both(
  p8d,
  "Figure8D_four_group_candidate_gene_dotplot",
  width = 13.4,
  height = 9.6
)

message("Step 12 finished successfully.")
message("Standalone Figure 8B, Figure 8C and Figure 8D were exported.")
message("Figure 8D is a four-group candidate-gene dotplot with compact wrapped group labels in this v6 version.")
message("No Figure 8B-D combined panel was generated in this v6 layout-fix version.")

## Public sequential end gate
.step03_out <- file.path(out_table_dir, "Step12_dual_high_Figure8D_four_group_candidate_gene_expression_summary.csv")
if (!file.exists(.step03_out)) stop("STEP 03 output missing: ", .step03_out, call. = FALSE)
cat("\nSTEP 03 PASS — dual-high niche analysis / Figure8B-D completed\n")
