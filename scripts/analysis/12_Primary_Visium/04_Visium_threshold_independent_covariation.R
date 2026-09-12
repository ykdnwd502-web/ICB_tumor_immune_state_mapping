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
## 18D_threshold_independent_spatial_covariation_revision_FIXED_INPUT.R
##
## Purpose:
##   Threshold-independent spatial co-variation analysis using the
##   confirmed Visium spot-level file.
##
## Confirmed input:
##   results/tables/spatial_melanoma_validation/
##   Step11_spatial_analysis_metadata_with_coordinates_and_scores.csv
##
## Main goal:
##   Address Reviewer 3 concern that top-quartile thresholding may
##   artificially generate the "dual-high" category.
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

project_dir <- Sys.getenv("ICB_PROJECT_DIR", unset = "D:/ICB_resistance_project")
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

table_dir  <- file.path(project_dir, "results/tables")
figure_dir <- file.path(project_dir, "results/figures")
log_dir    <- file.path(project_dir, "logs")

out_table_dir <- file.path(table_dir, "revision_spatial_covariation")
out_fig_dir <- file.path(project_dir, "results", "diagnostics", "12_Primary_Visium", "04_threshold_independent_covariation")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
reviewer_style_file <- file.path(project_dir, "scripts", "00_reviewer_figure_style_helpers.R")
if (file.exists(reviewer_style_file)) source(reviewer_style_file)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "pROC",
  "tibble",
  "purrr",
  "FNN"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE, type = "binary")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(FNN)
})

theme_revision <- function(base_size = 12.5, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2.5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      axis.title.y = element_text(face = "bold", size = base_size - 1),
      axis.title.x = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 1.5),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1.5),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.60)
    )
}

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
  message("Wrote: ", file)
}

safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 320) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file, plot = plot, width = width, height = height, dpi = dpi)
  message("Wrote: ", file)
}

zvec <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

cor_test_pair <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  
  if (n < 4) {
    return(data.frame(
      n = n,
      estimate = NA_real_,
      p_value = NA_real_,
      method = method,
      stringsAsFactors = FALSE
    ))
  }
  
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
  
  data.frame(
    n = n,
    estimate = unname(ct$estimate),
    p_value = ct$p.value,
    method = method,
    stringsAsFactors = FALSE
  )
}

partial_spearman_residual <- function(df, x_col, y_col, covariates) {
  covariates <- covariates[covariates %in% colnames(df)]
  covariates <- covariates[sapply(df[covariates], function(v) any(is.finite(as.numeric(v))))]
  
  vars <- c(x_col, y_col, covariates)
  dat <- df[, vars, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  
  if (nrow(dat) < 10 || length(covariates) == 0) {
    return(data.frame(
      n = nrow(dat),
      estimate = NA_real_,
      p_value = NA_real_,
      method = "partial Spearman by rank residualization",
      covariates = paste(covariates, collapse = " + "),
      status = ifelse(length(covariates) == 0, "Skipped: no covariates", "Skipped: insufficient complete cases"),
      stringsAsFactors = FALSE
    ))
  }
  
  rank_df <- as.data.frame(lapply(dat, function(v) rank(as.numeric(v), ties.method = "average", na.last = "keep")))
  
  fml_x <- as.formula(paste(x_col, "~", paste(covariates, collapse = " + ")))
  fml_y <- as.formula(paste(y_col, "~", paste(covariates, collapse = " + ")))
  
  rx <- residuals(lm(fml_x, data = rank_df))
  ry <- residuals(lm(fml_y, data = rank_df))
  
  ct <- suppressWarnings(stats::cor.test(rx, ry, method = "pearson"))
  
  data.frame(
    n = nrow(dat),
    estimate = unname(ct$estimate),
    p_value = ct$p.value,
    method = "partial Spearman by rank residualization",
    covariates = paste(covariates, collapse = " + "),
    status = "Completed",
    stringsAsFactors = FALSE
  )
}

make_knn_index <- function(coords, k = 6) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"
  nn <- FNN::get.knn(coords, k = k)
  nn$nn.index
}

spatial_lag <- function(v, knn_idx) {
  v <- as.numeric(v)
  apply(knn_idx, 1, function(idx) mean(v[idx], na.rm = TRUE))
}

global_bivariate_moran <- function(x, y, coords, k = 6, n_perm = 999, seed = 20260608) {
  set.seed(seed + k)
  
  ok <- is.finite(x) & is.finite(y) & is.finite(coords[, 1]) & is.finite(coords[, 2])
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  coords <- as.matrix(coords[ok, , drop = FALSE])
  
  n <- length(x)
  
  if (n < k + 5) {
    return(list(
      summary = data.frame(
        k = k,
        n = n,
        I_xy = NA_real_,
        I_yx = NA_real_,
        I_symmetric = NA_real_,
        perm_mean = NA_real_,
        perm_sd = NA_real_,
        z_perm = NA_real_,
        p_greater = NA_real_,
        p_two_sided = NA_real_,
        n_perm = n_perm
      ),
      local = data.frame(),
      perm = data.frame()
    ))
  }
  
  knn_idx <- make_knn_index(coords, k = k)
  
  zx <- zvec(x)
  zy <- zvec(y)
  
  lag_y <- spatial_lag(zy, knn_idx)
  lag_x <- spatial_lag(zx, knn_idx)
  
  I_xy <- mean(zx * lag_y, na.rm = TRUE)
  I_yx <- mean(zy * lag_x, na.rm = TRUE)
  I_sym <- mean(c(I_xy, I_yx), na.rm = TRUE)
  
  local_score <- 0.5 * (zx * lag_y + zy * lag_x)
  
  perm_vals <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    zy_perm <- sample(zy, size = length(zy), replace = FALSE)
    lag_y_perm <- spatial_lag(zy_perm, knn_idx)
    I_xy_perm <- mean(zx * lag_y_perm, na.rm = TRUE)
    
    zx_perm <- sample(zx, size = length(zx), replace = FALSE)
    lag_x_perm <- spatial_lag(zx_perm, knn_idx)
    I_yx_perm <- mean(zy * lag_x_perm, na.rm = TRUE)
    
    perm_vals[b] <- mean(c(I_xy_perm, I_yx_perm), na.rm = TRUE)
  }
  
  p_greater <- (sum(perm_vals >= I_sym, na.rm = TRUE) + 1) / (n_perm + 1)
  p_two <- (sum(abs(perm_vals - mean(perm_vals, na.rm = TRUE)) >= abs(I_sym - mean(perm_vals, na.rm = TRUE)), na.rm = TRUE) + 1) / (n_perm + 1)
  
  perm_mean <- mean(perm_vals, na.rm = TRUE)
  perm_sd <- stats::sd(perm_vals, na.rm = TRUE)
  z_perm <- (I_sym - perm_mean) / perm_sd
  
  summary <- data.frame(
    k = k,
    n = n,
    I_xy = I_xy,
    I_yx = I_yx,
    I_symmetric = I_sym,
    perm_mean = perm_mean,
    perm_sd = perm_sd,
    z_perm = z_perm,
    p_greater = p_greater,
    p_two_sided = p_two,
    n_perm = n_perm,
    stringsAsFactors = FALSE
  )
  
  local <- data.frame(
    spatial_x = coords[, 1],
    spatial_y = coords[, 2],
    z_MyeloidTreg = zx,
    z_DediffStromal = zy,
    lag_z_MyeloidTreg = lag_x,
    lag_z_DediffStromal = lag_y,
    local_bivariate_colocalization = local_score,
    stringsAsFactors = FALSE
  )
  
  perm <- data.frame(
    k = k,
    permutation = seq_len(n_perm),
    I_symmetric_perm = perm_vals,
    stringsAsFactors = FALSE
  )
  
  list(summary = summary, local = local, perm = perm)
}

safe_fisher_or <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
  
  if (is.null(ft)) {
    return(data.frame(
      OR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p_value = NA_real_
    ))
  }
  
  data.frame(
    OR = unname(ft$estimate),
    CI_low = unname(ft$conf.int[1]),
    CI_high = unname(ft$conf.int[2]),
    p_value = ft$p.value
  )
}

threshold_grid_or <- function(df, x_col, y_col) {
  cut_defs <- list(
    top15 = 0.85,
    top20 = 0.80,
    top25 = 0.75,
    top30 = 0.70,
    top35 = 0.65,
    median = 0.50
  )
  
  out <- list()
  
  for (nm in names(cut_defs)) {
    q <- cut_defs[[nm]]
    
    x_cut <- stats::quantile(df[[x_col]], probs = q, na.rm = TRUE, names = FALSE)
    y_cut <- stats::quantile(df[[y_col]], probs = q, na.rm = TRUE, names = FALSE)
    
    x_high <- df[[x_col]] >= x_cut
    y_high <- df[[y_col]] >= y_cut
    
    both_high <- sum(x_high & y_high, na.rm = TRUE)
    x_only <- sum(x_high & !y_high, na.rm = TRUE)
    y_only <- sum(!x_high & y_high, na.rm = TRUE)
    neither <- sum(!x_high & !y_high, na.rm = TRUE)
    
    or_tbl <- safe_fisher_or(
      a = both_high,
      b = x_only,
      c = y_only,
      d = neither
    )
    
    out[[length(out) + 1]] <- cbind(
      data.frame(
        CutoffName = nm,
        QuantileHighDefinition = q,
        MyeloidTreg_cutoff = x_cut,
        DediffStromal_cutoff = y_cut,
        BothHigh = both_high,
        MyeloidTregHighOnly = x_only,
        DediffStromalHighOnly = y_only,
        NeitherHigh = neither,
        stringsAsFactors = FALSE
      ),
      or_tbl
    )
  }
  
  x_high <- df[[x_col]] > 0
  y_high <- df[[y_col]] > 0
  
  both_high <- sum(x_high & y_high, na.rm = TRUE)
  x_only <- sum(x_high & !y_high, na.rm = TRUE)
  y_only <- sum(!x_high & y_high, na.rm = TRUE)
  neither <- sum(!x_high & !y_high, na.rm = TRUE)
  
  or_tbl <- safe_fisher_or(
    a = both_high,
    b = x_only,
    c = y_only,
    d = neither
  )
  
  out[[length(out) + 1]] <- cbind(
    data.frame(
      CutoffName = "z_gt_0",
      QuantileHighDefinition = NA_real_,
      MyeloidTreg_cutoff = 0,
      DediffStromal_cutoff = 0,
      BothHigh = both_high,
      MyeloidTregHighOnly = x_only,
      DediffStromalHighOnly = y_only,
      NeitherHigh = neither,
      stringsAsFactors = FALSE
    ),
    or_tbl
  )
  
  dplyr::bind_rows(out)
}

############################################################
## 1. Load fixed Visium table
############################################################

visium_file <- file.path(
  project_dir,
  "results/tables/spatial_melanoma_validation/Step11_spatial_analysis_metadata_with_coordinates_and_scores.csv"
)

if (!file.exists(visium_file)) {
  stop("Fixed Visium file not found: ", visium_file)
}

message("Using fixed Visium spot-level table: ", visium_file)

visium_df <- data.table::fread(
  visium_file,
  data.table = FALSE,
  check.names = FALSE
)

colnames(visium_df) <- gsub("\\s+", "_", colnames(visium_df))
colnames(visium_df) <- gsub("-", "_", colnames(visium_df))
colnames(visium_df) <- gsub("/", "_", colnames(visium_df))
colnames(visium_df) <- gsub("\\.+", "_", colnames(visium_df))
colnames(visium_df) <- gsub("__+", "_", colnames(visium_df))

if ("SpotID" %in% colnames(visium_df) && !"Spot" %in% colnames(visium_df)) {
  colnames(visium_df)[colnames(visium_df) == "SpotID"] <- "Spot"
}

if ("nFeature_Spatial_audit" %in% colnames(visium_df) && !"nFeature_Spatial" %in% colnames(visium_df)) {
  colnames(visium_df)[colnames(visium_df) == "nFeature_Spatial_audit"] <- "nFeature_Spatial"
}

if ("nCount_Spatial_audit" %in% colnames(visium_df) && !"nCount_Spatial" %in% colnames(visium_df)) {
  colnames(visium_df)[colnames(visium_df) == "nCount_Spatial_audit"] <- "nCount_Spatial"
}

if ("percent_mt_audit" %in% colnames(visium_df) && !"percent_mt" %in% colnames(visium_df)) {
  colnames(visium_df)[colnames(visium_df) == "percent_mt_audit"] <- "percent_mt"
}

required_numeric_cols <- c(
  "spatial_x",
  "spatial_y",
  "nFeature_Spatial",
  "nCount_Spatial",
  "percent_mt",
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

for (cc in intersect(required_numeric_cols, colnames(visium_df))) {
  visium_df[[cc]] <- suppressWarnings(as.numeric(visium_df[[cc]]))
}

required_cols <- c(
  "Spot",
  "spatial_x",
  "spatial_y",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling"
)

missing_cols <- setdiff(required_cols, colnames(visium_df))

if (length(missing_cols) > 0) {
  stop("Missing required columns in fixed Visium file: ", paste(missing_cols, collapse = ", "))
}

x_col <- "Myeloid_Treg_Immunosuppressive"
y_col <- "Tumor_dedifferentiation_Stromal_remodeling"

visium_df <- visium_df %>%
  dplyr::filter(
    is.finite(.data[[x_col]]),
    is.finite(.data[[y_col]]),
    is.finite(spatial_x),
    is.finite(spatial_y)
  )

message("Loaded fixed Visium table: ", nrow(visium_df), " spots x ", ncol(visium_df), " columns")

safe_write_csv(
  data.frame(
    SelectedInputFile = visium_file,
    n_spots_retained = nrow(visium_df),
    n_cols = ncol(visium_df),
    has_nCount_Spatial = "nCount_Spatial" %in% colnames(visium_df),
    has_nFeature_Spatial = "nFeature_Spatial" %in% colnames(visium_df),
    has_percent_mt = "percent_mt" %in% colnames(visium_df),
    stringsAsFactors = FALSE
  ),
  file.path(out_table_dir, "18D_selected_visium_input_summary.csv")
)

safe_write_csv(
  visium_df,
  file.path(out_table_dir, "18D_visium_spot_state_coordinate_table_used.csv")
)

############################################################
## 2. Continuous correlation
############################################################

raw_spearman <- cor_test_pair(visium_df[[x_col]], visium_df[[y_col]], method = "spearman") %>%
  dplyr::mutate(
    Analysis = "Raw continuous Spearman correlation",
    Covariates = "None"
  )

raw_pearson <- cor_test_pair(visium_df[[x_col]], visium_df[[y_col]], method = "pearson") %>%
  dplyr::mutate(
    Analysis = "Raw continuous Pearson correlation",
    Covariates = "None"
  )

qc_covariates <- intersect(c("nCount_Spatial", "nFeature_Spatial", "percent_mt"), colnames(visium_df))

partial_qc <- partial_spearman_residual(
  df = visium_df,
  x_col = x_col,
  y_col = y_col,
  covariates = qc_covariates
) %>%
  dplyr::mutate(
    Analysis = "QC-adjusted partial Spearman correlation",
    Covariates = paste(qc_covariates, collapse = " + ")
  )

cor_summary <- dplyr::bind_rows(
  raw_spearman,
  raw_pearson,
  partial_qc
) %>%
  dplyr::select(
    Analysis,
    method,
    n,
    estimate,
    p_value,
    Covariates,
    dplyr::everything()
  )

safe_write_csv(
  cor_summary,
  file.path(out_table_dir, "Visium_continuous_state_correlation_and_partial_correlation.csv")
)

############################################################
## 3. kNN bivariate spatial Moran
############################################################

coords <- as.matrix(visium_df[, c("spatial_x", "spatial_y")])

k_values <- c(4, 6, 8, 12)
n_perm <- 999

moran_results <- list()
local_results <- list()
perm_results <- list()

for (kk in k_values) {
  message("Running bivariate spatial Moran analysis: k = ", kk)
  
  res <- global_bivariate_moran(
    x = visium_df[[x_col]],
    y = visium_df[[y_col]],
    coords = coords,
    k = kk,
    n_perm = n_perm,
    seed = 20260608
  )
  
  moran_results[[as.character(kk)]] <- res$summary
  
  local_tmp <- res$local
  if (nrow(local_tmp) > 0) {
    local_tmp$k <- kk
    local_tmp$Spot <- visium_df$Spot[seq_len(nrow(local_tmp))]
    local_results[[as.character(kk)]] <- local_tmp
  }
  
  perm_results[[as.character(kk)]] <- res$perm
}

moran_summary <- dplyr::bind_rows(moran_results) %>%
  dplyr::mutate(
    Interpretation = "Positive I_symmetric indicates that high values of one state tend to be spatially adjacent to high values of the other state under row-standardized kNN weights."
  )

safe_write_csv(
  moran_summary,
  file.path(out_table_dir, "Visium_bivariate_spatial_Moran_kNN_sensitivity.csv")
)

local_all <- dplyr::bind_rows(local_results)
perm_all <- dplyr::bind_rows(perm_results)

safe_write_csv(
  local_all,
  file.path(out_table_dir, "Visium_local_bivariate_colocalization_scores_all_k.csv")
)

safe_write_csv(
  perm_all,
  file.path(out_table_dir, "Visium_bivariate_spatial_Moran_permutation_nulls.csv")
)

local_k6 <- local_all %>% dplyr::filter(k == 6)

safe_write_csv(
  local_k6,
  file.path(out_table_dir, "Visium_local_bivariate_colocalization_scores_k6_primary.csv")
)

############################################################
## 4. Threshold-grid OR sensitivity
############################################################

threshold_or <- threshold_grid_or(visium_df, x_col, y_col)

safe_write_csv(
  threshold_or,
  file.path(out_table_dir, "Visium_threshold_grid_dual_high_Fisher_OR_sensitivity.csv")
)

############################################################
## 5. Figures
############################################################

p_scatter <- ggplot(
  visium_df,
  aes(
    x = .data[[x_col]],
    y = .data[[y_col]]
  )
) +
  geom_point(alpha = 0.55, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  labs(
    title = "Continuous co-variation between spatial state scores",
    subtitle = "Visium spots; threshold-independent state-score relationship",
    x = "myeloid–Treg immunosuppressive state score",
    y = "tumor-dedifferentiation/stromal-remodeling state score"
  ) +
  theme_revision(base_size = 12.0)

safe_ggsave(
  file.path(out_fig_dir, "18D_Visium_continuous_state_score_scatter.pdf"),
  p_scatter,
  width = 6.5,
  height = 5.5
)

safe_ggsave(
  file.path(out_fig_dir, "18D_Visium_continuous_state_score_scatter.png"),
  p_scatter,
  width = 6.5,
  height = 5.5
)

plot_df_long <- visium_df %>%
  dplyr::select(Spot, spatial_x, spatial_y, all_of(c(x_col, y_col))) %>%
  tidyr::pivot_longer(
    cols = all_of(c(x_col, y_col)),
    names_to = "State",
    values_to = "Score"
  ) %>%
  dplyr::mutate(
    State = dplyr::recode(
      State,
      Myeloid_Treg_Immunosuppressive = "myeloid–Treg immunosuppressive",
      Tumor_dedifferentiation_Stromal_remodeling = "tumor-dedifferentiation/stromal-remodeling"
    )
  )

p_map_states <- ggplot(
  plot_df_long,
  aes(x = spatial_x, y = spatial_y, color = Score)
) +
  geom_point(size = 1.1, alpha = 0.9) +
  scale_y_reverse() +
  facet_wrap(~ State, ncol = 2) +
  scale_color_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    na.value = "grey85"
  ) +
  coord_equal() +
  labs(
    title = "Continuous spatial distribution of the two state scores",
    subtitle = "Visium spot-level state scores",
    x = "Spatial x",
    y = "Spatial y",
    color = "Score"
  ) +
  theme_revision(base_size = 12.0)

safe_ggsave(
  file.path(out_fig_dir, "18D_Visium_continuous_state_score_spatial_maps.pdf"),
  p_map_states,
  width = 9,
  height = 4.8
)

safe_ggsave(
  file.path(out_fig_dir, "18D_Visium_continuous_state_score_spatial_maps.png"),
  p_map_states,
  width = 9,
  height = 4.8
)

p_moran <- moran_summary %>%
  dplyr::mutate(k = factor(k, levels = k_values)) %>%
  ggplot(aes(x = k, y = I_symmetric)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.6) +
  geom_line(aes(group = 1), linewidth = 0.5) +
  geom_errorbar(
    aes(
      ymin = perm_mean - 1.96 * perm_sd,
      ymax = perm_mean + 1.96 * perm_sd
    ),
    width = 0.15,
    linewidth = 0.4
  ) +
  labs(
    title = "Threshold-independent bivariate spatial co-variation",
    subtitle = "Observed symmetric bivariate Moran-type statistic across kNN definitions;\nerror bars show permutation-null 95% range",
    x = "k nearest neighbors",
    y = "Symmetric bivariate spatial Moran statistic"
  ) +
  theme_revision(base_size = 12.0) +
  theme(plot.subtitle = element_text(size = 10.5, lineheight = 1.0))

safe_ggsave(
  file.path(out_fig_dir, "18D_bivariate_spatial_Moran_kNN_sensitivity.pdf"),
  p_moran,
  width = 6.5,
  height = 5
)

safe_ggsave(
  file.path(out_fig_dir, "18D_bivariate_spatial_Moran_kNN_sensitivity.png"),
  p_moran,
  width = 6.5,
  height = 5
)

p_local <- ggplot(
  local_k6,
  aes(x = spatial_x, y = spatial_y, color = local_bivariate_colocalization)
) +
  geom_point(size = 1.1, alpha = 0.9) +
  scale_y_reverse() +
  scale_color_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    na.value = "grey85"
  ) +
  coord_equal() +
  labs(
    title = "Local bivariate co-localization score",
    subtitle = "Primary k = 6 kNN definition;\npositive values indicate local co-enrichment of the two states",
    x = "Spatial x",
    y = "Spatial y",
    color = "Local score"
  ) +
  theme_revision(base_size = 12.0) +
  theme(plot.subtitle = element_text(size = 10.5, lineheight = 1.0))

safe_ggsave(
  file.path(out_fig_dir, "18D_local_bivariate_colocalization_map_k6.pdf"),
  p_local,
  width = 6.5,
  height = 5.8
)

safe_ggsave(
  file.path(out_fig_dir, "18D_local_bivariate_colocalization_map_k6.png"),
  p_local,
  width = 6.5,
  height = 5.8
)

p_or <- threshold_or %>%
  dplyr::mutate(CutoffName = factor(CutoffName, levels = CutoffName)) %>%
  ggplot(aes(x = CutoffName, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.4) +
  scale_y_log10() +
  labs(
    title = "Threshold-grid sensitivity of dual-high co-enrichment",
    subtitle = "Fisher odds ratios across operational high-state definitions;\nsecondary to continuous spatial analyses",
    x = "High-state definition",
    y = "Fisher OR, log10 scale"
  ) +
  theme_revision(base_size = 12.0) +
  theme(plot.subtitle = element_text(size = 10.5, lineheight = 1.0), axis.text.x = element_text(angle = 35, hjust = 1, size = 10.75))

safe_ggsave(
  file.path(out_fig_dir, "18D_threshold_grid_dual_high_OR_sensitivity.pdf"),
  p_or,
  width = 7,
  height = 5
)

safe_ggsave(
  file.path(out_fig_dir, "18D_threshold_grid_dual_high_OR_sensitivity.png"),
  p_or,
  width = 7,
  height = 5
)

############################################################
## 6. Interpretation helper
############################################################

interpretation_summary <- data.frame(
  Item = c(
    "Purpose",
    "Primary threshold-independent test",
    "Continuous correlation",
    "Spatial co-variation",
    "Threshold-grid OR",
    "Interpretation boundary"
  ),
  RecommendedInterpretation = c(
    "Evaluate whether spatial co-enrichment of myeloid–Treg and tumor-dedifferentiation/stromal-remodeling states depends on top-quartile thresholding.",
    "Use continuous state-score correlation and kNN-based bivariate spatial Moran-type statistics with permutation testing.",
    "A positive continuous correlation supports co-variation independent of categorical high-state definitions.",
    "A positive bivariate spatial Moran statistic indicates that high values of one state tend to be spatially adjacent to high values of the other state within the analyzed section.",
    "Threshold-grid Fisher OR is retained as an operational sensitivity analysis rather than the sole evidence for dual-high organization.",
    "The analysis supports within-section spatial co-variation but does not establish patient-level generalizability, causality, or functional interaction."
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  interpretation_summary,
  file.path(out_table_dir, "18D_interpretation_summary.csv")
)

output_inventory <- data.frame(
  Type = c(
    "table",
    "table",
    "table",
    "table",
    "table",
    "figure",
    "figure",
    "figure",
    "figure",
    "figure"
  ),
  File = c(
    file.path(out_table_dir, "18D_selected_visium_input_summary.csv"),
    file.path(out_table_dir, "Visium_continuous_state_correlation_and_partial_correlation.csv"),
    file.path(out_table_dir, "Visium_bivariate_spatial_Moran_kNN_sensitivity.csv"),
    file.path(out_table_dir, "Visium_threshold_grid_dual_high_Fisher_OR_sensitivity.csv"),
    file.path(out_table_dir, "Visium_local_bivariate_colocalization_scores_k6_primary.csv"),
    file.path(out_fig_dir, "18D_Visium_continuous_state_score_scatter.pdf"),
    file.path(out_fig_dir, "18D_Visium_continuous_state_score_spatial_maps.pdf"),
    file.path(out_fig_dir, "18D_bivariate_spatial_Moran_kNN_sensitivity.pdf"),
    file.path(out_fig_dir, "18D_local_bivariate_colocalization_map_k6.pdf"),
    file.path(out_fig_dir, "18D_threshold_grid_dual_high_OR_sensitivity.pdf")
  ),
  stringsAsFactors = FALSE
)

output_inventory$Exists <- file.exists(output_inventory$File)

safe_write_csv(
  output_inventory,
  file.path(out_table_dir, "18D_output_inventory.csv")
)

sink(file.path(log_dir, "sessionInfo_18D_threshold_independent_spatial_covariation.txt"))
print(sessionInfo())
sink()

message("============================================================")
message("18D fixed-input threshold-independent spatial co-variation completed.")
message("Main correlation table: ", file.path(out_table_dir, "Visium_continuous_state_correlation_and_partial_correlation.csv"))
message("Main spatial table: ", file.path(out_table_dir, "Visium_bivariate_spatial_Moran_kNN_sensitivity.csv"))
message("Threshold-grid OR table: ", file.path(out_table_dir, "Visium_threshold_grid_dual_high_Fisher_OR_sensitivity.csv"))
message("============================================================")
## Public sequential end gate
.step04_out <- file.path(out_table_dir, "Visium_bivariate_spatial_Moran_kNN_sensitivity.csv")
if (!file.exists(.step04_out)) stop("STEP 04 output missing: ", .step04_out, call. = FALSE)
cat("\nSTEP 04 PASS — threshold-independent spatial covariation completed\n")
