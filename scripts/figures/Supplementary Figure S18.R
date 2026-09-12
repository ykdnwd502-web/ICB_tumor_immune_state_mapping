################--------------------------------------------
## Supplementary Figure S18. Composition-aware and permutation-supported spatial sensitivity analysis of Visium state-score co-variation
##
## Purpose:
##   1) Load spatial spot-level composition-residualized scores and audit tables from revision_visium_composition
##   2) Generate and export Supplementary Figure S18 (Panel A state-composition correlation, Panel B residualized spatial correlation, Panel C residualized dual-high OR, Panel D bivariate Moran's I with permutation support)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S18. Composition-aware and permutation-supported spatial sensitivity analysis of Visium state-score co-variation.png
##     - Supplementary Figure S18. Composition-aware and permutation-supported spatial sensitivity analysis of Visium state-score co-variation.jpg
##     - Supplementary Figure S18. Composition-aware and permutation-supported spatial sensitivity analysis of Visium state-score co-variation.pdf
################--------------------------------------------

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(scales)
  library(patchwork)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
transmute <- dplyr::transmute

################--------------------------------------------
## 0. Paths and fixed settings
################--------------------------------------------

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

spot_file <- file.path(
  project_dir,
  "results", "tables", "revision_visium_composition",
  "18E_Visium_spot_state_composition_residualized_scores.csv"
)

dpi_out <- 300
high_quantile <- 0.75
knn_k <- 6
n_perm <- 999
perm_seed <- 20260629

state_order <- c(
  "melanocytic differentiation",
  "tumor-dedifferentiation/stromal-remodeling",
  "myeloid–Treg immunosuppressive",
  "immune-defective/cold"
)

module_order <- c(
  "CAF/stromal",
  "myeloid",
  "T/NK",
  "B/plasma",
  "melanoma lineage",
  "endothelial"
)

method_order <- c(
  "Raw",
  "QC-residualized",
  "Composition-residualized",
  "QC/composition-residualized"
)

score_pair_map <- list(
  "Raw" = c(
    myeloid = "Myeloid_Treg_Immunosuppressive",
    dediff = "Tumor_dedifferentiation_Stromal_remodeling"
  ),
  "QC-residualized" = c(
    myeloid = "MyeloidTreg_resid_QC",
    dediff = "DediffStromal_resid_QC"
  ),
  "Composition-residualized" = c(
    myeloid = "MyeloidTreg_resid_Composition",
    dediff = "DediffStromal_resid_Composition"
  ),
  "QC/composition-residualized" = c(
    myeloid = "MyeloidTreg_resid_QC_Composition",
    dediff = "DediffStromal_resid_QC_Composition"
  )
)

################--------------------------------------------
## 1. Helper functions
################--------------------------------------------

read_any_table <- function(file) {
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  
  ext <- tolower(tools::file_ext(file))
  
  if (ext %in% c("csv", "txt")) {
    out <- suppressMessages(readr::read_csv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "tsv") {
    out <- suppressMessages(readr::read_tsv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "rds") {
    out <- readRDS(file)
    if (!is.data.frame(out)) {
      stop("RDS object is not a data.frame: ", file)
    }
  } else {
    stop("Unsupported file extension: ", file)
  }
  
  as.data.frame(out, stringsAsFactors = FALSE)
}

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

require_cols <- function(df, cols, label) {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(miss, collapse = ", "),
      "\nAvailable columns:\n",
      paste(colnames(df), collapse = ", ")
    )
  }
  invisible(TRUE)
}

wrap_state <- function(x) {
  dplyr::recode(
    as.character(x),
    "immune-defective/cold" = "immune-defective/\ncold",
    "myeloid–Treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
    "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
    "melanocytic differentiation" = "melanocytic\ndifferentiation",
    .default = as.character(x)
  )
}

wrap_method <- function(x) {
  dplyr::recode(
    as.character(x),
    "Raw" = "Raw",
    "QC-residualized" = "QC-residualized",
    "Composition-residualized" = "Composition-\nresidualized",
    "QC/composition-residualized" = "QC/composition-\nresidualized",
    .default = as.character(x)
  )
}

theme_s18 <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 12.8, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 9.2, hjust = 0, color = "black"),
    axis.title = element_text(size = 10.4, face = "bold", color = "black"),
    axis.text = element_text(size = 8.8, color = "black"),
    legend.title = element_text(size = 9.2, face = "bold", color = "black"),
    legend.text = element_text(size = 8.0, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

################--------------------------------------------
## 2. Load spot-level data
################--------------------------------------------

spot_df <- read_any_table(spot_file)

state_col_map <- c(
  "melanocytic differentiation" = "Melanocytic_Differentiation",
  "tumor-dedifferentiation/stromal-remodeling" = "Tumor_dedifferentiation_Stromal_remodeling",
  "myeloid–Treg immunosuppressive" = "Myeloid_Treg_Immunosuppressive",
  "immune-defective/cold" = "Immune_defective_Cold"
)

module_col_map <- c(
  "CAF/stromal" = "CAF_Stromal",
  "myeloid" = "Myeloid",
  "T/NK" = "T_NK",
  "B/plasma" = "B_Plasma",
  "melanoma lineage" = "Melanoma_Lineage",
  "endothelial" = "Endothelial"
)

required_cols <- c(
  "spatial_x",
  "spatial_y",
  unname(state_col_map),
  unname(module_col_map),
  unlist(score_pair_map)
)

require_cols(spot_df, required_cols, "Spot-level table")

score_pair_audit <- data.frame(
  Method = method_order,
  myeloid_col = vapply(score_pair_map, function(x) x[["myeloid"]], character(1)),
  dediff_col = vapply(score_pair_map, function(x) x[["dediff"]], character(1)),
  usable = vapply(
    score_pair_map,
    function(x) all(c(x[["myeloid"]], x[["dediff"]]) %in% colnames(spot_df)),
    logical(1)
  ),
  stringsAsFactors = FALSE
)

usable_methods <- score_pair_audit$Method[score_pair_audit$usable]

################--------------------------------------------
## 3. Panel A: state-composition correlation heatmap
################--------------------------------------------

panelA_df <- expand.grid(
  State = names(state_col_map),
  Module = names(module_col_map),
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    rho = suppressWarnings(cor(
      as.numeric(spot_df[[state_col_map[[State]]]]),
      as.numeric(spot_df[[module_col_map[[Module]]]]),
      method = "spearman",
      use = "pairwise.complete.obs"
    ))
  ) %>%
  ungroup() %>%
  mutate(
    State = factor(State, levels = state_order),
    Module = factor(Module, levels = module_order),
    label = sprintf("%.2f", rho)
  ) %>%
  filter(is.finite(rho), !is.na(State), !is.na(Module))

################--------------------------------------------
## 4. Panel B: raw versus residualized spatial correlation
################--------------------------------------------

panelB_df <- lapply(usable_methods, function(m) {
  pair <- score_pair_map[[m]]
  x <- as.numeric(spot_df[[pair[["myeloid"]]]])
  y <- as.numeric(spot_df[[pair[["dediff"]]]])
  rho <- suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
  data.frame(Method = m, rho = rho, stringsAsFactors = FALSE)
}) %>%
  bind_rows() %>%
  filter(is.finite(rho)) %>%
  mutate(
    Method = factor(Method, levels = method_order),
    label = sprintf("%.2f", rho)
  ) %>%
  arrange(Method)

################--------------------------------------------
## 5. Panel C: Fisher OR from top-quartile dual-high definition
################--------------------------------------------

compute_fisher_or <- function(x, y, q = 0.75) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  x_cut <- as.numeric(stats::quantile(x, probs = q, na.rm = TRUE, names = FALSE))
  y_cut <- as.numeric(stats::quantile(y, probs = q, na.rm = TRUE, names = FALSE))
  
  x_high <- x >= x_cut
  y_high <- y >= y_cut
  
  tab <- table(
    y_high = factor(y_high, levels = c(FALSE, TRUE)),
    x_high = factor(x_high, levels = c(FALSE, TRUE))
  )
  
  ft <- suppressWarnings(stats::fisher.test(tab))
  
  data.frame(
    OR = as.numeric(ft$estimate),
    CI_low = as.numeric(ft$conf.int[1]),
    CI_high = as.numeric(ft$conf.int[2]),
    x_cutoff = x_cut,
    y_cutoff = y_cut,
    n = length(x),
    n_x_high = sum(x_high),
    n_y_high = sum(y_high),
    n_both_high = sum(x_high & y_high),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

panelC_df <- lapply(usable_methods, function(m) {
  pair <- score_pair_map[[m]]
  out <- compute_fisher_or(
    x = spot_df[[pair[["myeloid"]]]],
    y = spot_df[[pair[["dediff"]]]],
    q = high_quantile
  )
  out$Method <- m
  out
}) %>%
  bind_rows() %>%
  filter(is.finite(OR), OR > 0) %>%
  mutate(
    Method = factor(Method, levels = method_order),
    label = sprintf("%.2f", OR)
  ) %>%
  arrange(Method)

################--------------------------------------------
## 6. Panel D: symmetric bivariate Moran's I and permutation
################--------------------------------------------

compute_knn_indices <- function(coords, k) {
  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knn(coords, k = k)
    return(nn$nn.index)
  }
  
  d <- as.matrix(stats::dist(coords))
  diag(d) <- Inf
  idx <- t(apply(d, 1, function(v) order(v)[seq_len(k)]))
  idx
}

symmetric_bivar_I_from_index <- function(x, y, idx) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  
  xz <- as.numeric(scale(x[keep]))
  yz <- as.numeric(scale(y[keep]))
  
  if (!all(keep)) {
    idx <- idx[keep, , drop = FALSE]
  }
  
  I_xy <- mean(vapply(
    seq_along(xz),
    function(i) mean(xz[i] * yz[idx[i, ]]),
    numeric(1)
  ))
  
  I_yx <- mean(vapply(
    seq_along(yz),
    function(i) mean(yz[i] * xz[idx[i, ]]),
    numeric(1)
  ))
  
  (I_xy + I_yx) / 2
}

complete_rows <- which(
  is.finite(as.numeric(spot_df$spatial_x)) &
    is.finite(as.numeric(spot_df$spatial_y))
)

for (m in usable_methods) {
  pair <- score_pair_map[[m]]
  xm <- as.numeric(spot_df[[pair[["myeloid"]]]])
  ym <- as.numeric(spot_df[[pair[["dediff"]]]])
  complete_rows <- intersect(complete_rows, which(is.finite(xm) & is.finite(ym)))
}

coords <- as.matrix(spot_df[complete_rows, c("spatial_x", "spatial_y")])
coords <- apply(coords, 2, as.numeric)

knn_idx <- compute_knn_indices(coords, k = knn_k)

panelD_df <- lapply(usable_methods, function(m) {
  pair <- score_pair_map[[m]]
  x <- as.numeric(spot_df[[pair[["myeloid"]]]])[complete_rows]
  y <- as.numeric(spot_df[[pair[["dediff"]]]])[complete_rows]
  
  data.frame(
    Method = m,
    I = symmetric_bivar_I_from_index(x, y, knn_idx),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  filter(is.finite(I)) %>%
  mutate(
    Method = factor(Method, levels = method_order),
    label = sprintf("%.3f", I)
  ) %>%
  arrange(Method)

raw_pair <- score_pair_map[["Raw"]]
raw_x <- as.numeric(spot_df[[raw_pair[["myeloid"]]]])[complete_rows]
raw_y <- as.numeric(spot_df[[raw_pair[["dediff"]]]])[complete_rows]
raw_obs_I <- symmetric_bivar_I_from_index(raw_x, raw_y, knn_idx)

set.seed(perm_seed)
raw_perm_I <- vapply(seq_len(n_perm), function(i) {
  symmetric_bivar_I_from_index(raw_x, sample(raw_y), knn_idx)
}, numeric(1))

perm_p <- (1 + sum(raw_perm_I >= raw_obs_I, na.rm = TRUE)) / (n_perm + 1)
perm_p_label <- if (perm_p <= 0.001) "P <= 0.001" else paste0("P = ", sprintf("%.3f", perm_p))

################--------------------------------------------
## 7. Plots (Panel C numeric labels returned to original height, but nudged right to avoid error bar line)
################--------------------------------------------

pA <- ggplot(panelA_df, aes(x = Module, y = State, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 2.8, color = "black") +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(labels = wrap_state, drop = FALSE) +
  labs(
    title = "Visium state–composition correlation",
    subtitle = "Spearman correlation across spots",
    x = "Composition marker module",
    y = NULL
  ) +
  theme_s18 +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  )

pB <- ggplot(panelB_df, aes(x = Method, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_col(width = 0.62, fill = "grey35") +
  geom_text(aes(label = label), vjust = -0.28, size = 3.1) +
  scale_x_discrete(labels = wrap_method) +
  scale_y_continuous(
    limits = c(0, max(panelB_df$rho, na.rm = TRUE) * 1.20),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Raw versus residualized spatial correlation",
    subtitle = "Myeloid–Treg immunosuppressive versus tumor-dedifferentiation/stromal-remodeling",
    x = NULL,
    y = "Spearman rho"
  ) +
  theme_s18 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

# 核心修改：恢复原先垂直高度（vjust = -0.35），并使用 nudge_x 将数值向右水平错开，非加粗字体
pC <- ggplot(panelC_df, aes(x = Method, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_col(width = 0.62, fill = "grey35") +
  geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high),
    width = 0.18,
    linewidth = 0.50,
    color = "black"
  ) +
  geom_text(
    aes(label = label), 
    vjust = -0.35,  # 恢复原高度
    nudge_x = 0.15, # 向右平移错开误差条竖线
    size = 3.0,
    fontface = "plain"  # 非加粗
  ) +
  scale_x_discrete(labels = wrap_method) +
  scale_y_continuous(
    limits = c(0, max(panelC_df$CI_high, na.rm = TRUE) * 1.15),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Raw versus residualized dual-high enrichment",
    subtitle = "Top-quartile high-score spot definition; dashed line indicates odds ratio = 1",
    x = NULL,
    y = "Fisher odds ratio"
  ) +
  theme_s18 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

pD <- ggplot(panelD_df, aes(x = Method, y = I)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_col(width = 0.62, fill = "grey35") +
  geom_text(aes(label = label), vjust = -0.35, size = 3.0) +
  scale_x_discrete(labels = wrap_method) +
  scale_y_continuous(
    limits = c(0, max(panelD_df$I, na.rm = TRUE) * 1.20),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Permutation-supported bivariate Moran-type statistic",
    subtitle = paste0("kNN k = ", knn_k, "; within-section permutation ", perm_p_label),
    x = NULL,
    y = "Symmetric bivariate Moran's I"
  ) +
  theme_s18 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

################--------------------------------------------
## 8. Assemble and save
################--------------------------------------------

supp_s18 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 19, face = "bold", family = "sans"),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

base_filename <- "Supplementary Figure S18. Composition-aware and permutation-supported spatial sensitivity analysis of Visium state-score co-variation"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 14.2
fig_height <- 10.0

safe_ggsave(out_png, supp_s18, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s18, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s18, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 9. Audit outputs
################--------------------------------------------

audit <- data.frame(
  spot_file = spot_file,
  n_spots = nrow(spot_df),
  high_quantile = high_quantile,
  knn_k = knn_k,
  n_perm = n_perm,
  perm_seed = perm_seed,
  raw_observed_I = raw_obs_I,
  raw_permutation_p = perm_p,
  usable_methods = paste(usable_methods, collapse = "; "),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S18_composition_permutation_spatial_sensitivity_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

write.csv(score_pair_audit, file.path(fig_dir, "Supplementary_Figure_S18_score_pair_column_audit.csv"), row.names = FALSE)
write.csv(panelA_df, file.path(fig_dir, "Supplementary_Figure_S18_panelA_state_composition_correlation.csv"), row.names = FALSE)
write.csv(panelB_df, file.path(fig_dir, "Supplementary_Figure_S18_panelB_residualized_spatial_correlation.csv"), row.names = FALSE)
write.csv(panelC_df, file.path(fig_dir, "Supplementary_Figure_S18_panelC_residualized_dual_high_OR.csv"), row.names = FALSE)
write.csv(panelD_df, file.path(fig_dir, "Supplementary_Figure_S18_panelD_bivariate_Moran.csv"), row.names = FALSE)
write.csv(data.frame(null_I = raw_perm_I), file.path(fig_dir, "Supplementary_Figure_S18_raw_bivariate_Moran_permutation_null.csv"), row.names = FALSE)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S18 script finished successfully.")