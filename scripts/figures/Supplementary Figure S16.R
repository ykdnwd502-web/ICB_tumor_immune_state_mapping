################--------------------------------------------
## Supplementary Figure S16. Threshold-independent spatial co-variation between myeloid–Treg and dedifferentiation-stromal scores in Visium
##
## Purpose:
##   1) Load spatial state scores and threshold sensitivity OR tables from spatial_melanoma_validation
##   2) Generate and export Supplementary Figure S16 (Panel A continuous score scatter, Panel B bivariate Moran's I across kNN, Panel C threshold grid sensitivity OR)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S16. Threshold-independent spatial co-variation between myeloid–Treg and dedifferentiation-stromal scores in Visium.png
##     - Supplementary Figure S16. Threshold-independent spatial co-variation between myeloid–Treg and dedifferentiation-stromal scores in Visium.jpg
##     - Supplementary Figure S16. Threshold-independent spatial co-variation between myeloid–Treg and dedifferentiation-stromal scores in Visium.pdf
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

## Avoid common namespace conflicts in long sessions.
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
transmute <- dplyr::transmute

################--------------------------------------------
## 0. Paths and settings
################--------------------------------------------

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

table_root <- file.path(project_dir, "results", "tables", "spatial_melanoma_validation")
fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

# Updated path definition to point to the new filename
panelA_file <- file.path(
  table_root,
  "Step11A_spatial_final_four_ICB_state_scores.csv"
)

panelC_file <- file.path(
  table_root,
  "Step11_spatial_dual_high_threshold_sensitivity_OR.csv"
)

knn_values <- c(4, 6, 8, 12)

myeloid_label <- "myeloid–Treg immunosuppressive"
dediff_label <- "tumor-dedifferentiation/stromal-remodeling"

################--------------------------------------------
## 1. Helper functions
################--------------------------------------------

read_any_table <- function(file) {
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  
  ext <- tolower(tools::file_ext(file))
  
  if (ext %in% c("csv", "txt")) {
    out <- suppressMessages(
      readr::read_csv(file, show_col_types = FALSE, guess_max = 100000)
    )
  } else if (ext == "tsv") {
    out <- suppressMessages(
      readr::read_tsv(file, show_col_types = FALSE, guess_max = 100000)
    )
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

find_col <- function(df, patterns, label, required = TRUE, exclude_patterns = character(0)) {
  cols <- colnames(df)
  
  for (pat in patterns) {
    hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE)
    
    if (length(exclude_patterns) > 0 && length(hit) > 0) {
      for (ep in exclude_patterns) {
        hit <- hit[!grepl(ep, hit, ignore.case = TRUE)]
      }
    }
    
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  
  if (required) {
    stop(
      "Could not find column for ", label, ". Tried patterns:\n",
      paste(patterns, collapse = "\n"),
      "\nAvailable columns:\n",
      paste(cols, collapse = ", ")
    )
  }
  
  NA_character_
}

top_fraction_label <- function(x) {
  xf <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    is.finite(xf) & abs(xf - 0.15) < 1e-8 ~ "top 15%",
    is.finite(xf) & abs(xf - 0.20) < 1e-8 ~ "top 20%",
    is.finite(xf) & abs(xf - 0.25) < 1e-8 ~ "top 25%",
    is.finite(xf) & abs(xf - 0.30) < 1e-8 ~ "top 30%",
    is.finite(xf) & abs(xf - 0.35) < 1e-8 ~ "top 35%",
    is.finite(xf) & xf > 0 & xf < 1 ~ paste0("top ", round(100 * xf), "%"),
    is.finite(xf) & xf >= 1 ~ paste0("top ", round(xf), "%"),
    TRUE ~ NA_character_
  )
}

clean_mode_label <- function(x) {
  x0 <- as.character(x)
  xl <- tolower(x0)
  dplyr::case_when(
    grepl("pooled", xl) ~ "pooled",
    grepl("per.*sample|sample", xl) ~ "per-sample",
    TRUE ~ x0
  )
}

theme_s16 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 13.0, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 9.5, hjust = 0, color = "black"),
    axis.title = element_text(size = 11.0, face = "bold", color = "black"),
    axis.text = element_text(size = 9.4, color = "black"),
    legend.title = element_text(size = 9.7, face = "bold", color = "black"),
    legend.text = element_text(size = 8.6, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

################--------------------------------------------
## 2. Load and harmonize Panel A spot-level continuous scores
################--------------------------------------------

panelA_raw <- read_any_table(panelA_file)

message("Panel A/Panel B source table: ", panelA_file)
message("Columns in spot-level state-score table:")
message(paste(colnames(panelA_raw), collapse = ", "))

myeloid_col_A <- find_col(
  panelA_raw,
  patterns = c(
    "myeloid.*treg.*immunosuppressive.*score",
    "myeloid.*treg.*state.*score",
    "myeloid.*treg.*score",
    "myeloid.*treg",
    "Myeloid_Treg",
    "Myeloid.*Treg",
    "Treg.*myeloid"
  ),
  label = "myeloid–Treg immunosuppressive state score",
  exclude_patterns = c("dominant", "label", "high", "cutoff", "threshold", "rank")
)

dediff_col_A <- find_col(
  panelA_raw,
  patterns = c(
    "tumor.*dedifferentiation.*stromal.*remodeling.*score",
    "tumor.*dediff.*stromal.*score",
    "dediff.*stromal.*score",
    "dediff.*stromal",
    "Tumor_dedifferentiation_Stromal_remodeling",
    "Stromal.*remodel"
  ),
  label = "tumor-dedifferentiation/stromal-remodeling state score",
  exclude_patterns = c("dominant", "label", "high", "cutoff", "threshold", "rank")
)

x_col_A <- find_col(
  panelA_raw,
  patterns = c("^imagecol$", "^pxl_col_in_fullres$", "^x$", "pixel.*x", "full.*x", "coord.*x", "array_col", "col"),
  label = "Visium x coordinate",
  required = FALSE
)

y_col_A <- find_col(
  panelA_raw,
  patterns = c("^imagerow$", "^pxl_row_in_fullres$", "^y$", "pixel.*y", "full.*y", "coord.*y", "array_row", "row"),
  label = "Visium y coordinate",
  required = FALSE
)

panelA_df <- panelA_raw %>%
  transmute(
    MyeloidTreg = as.numeric(.data[[myeloid_col_A]]),
    DediffStromal = as.numeric(.data[[dediff_col_A]]),
    coord_x = if (!is.na(x_col_A)) as.numeric(.data[[x_col_A]]) else NA_real_,
    coord_y = if (!is.na(y_col_A)) as.numeric(.data[[y_col_A]]) else NA_real_
  ) %>%
  filter(is.finite(MyeloidTreg), is.finite(DediffStromal))

if (nrow(panelA_df) == 0) {
  stop("Panel A harmonized table has zero rows.")
}

if (is.na(x_col_A) || is.na(y_col_A)) {
  stop(
    "Panel B is computed from kNN spatial coordinates, but coordinate columns were not found.\n",
    "Please confirm Step11A_spatial_final_four_ICB_state_scores.csv contains imagecol and imagerow."
  )
}

panelB_input <- panelA_df %>%
  filter(is.finite(coord_x), is.finite(coord_y))

if (nrow(panelB_input) < 100) {
  stop("Too few spots with finite coordinates for Panel B kNN computation: n = ", nrow(panelB_input))
}

message("Panel A score columns:")
message("  Myeloid–Treg: ", myeloid_col_A)
message("  Tumor-dedifferentiation/stromal-remodeling: ", dediff_col_A)
message("Coordinate columns:")
message("  x: ", x_col_A)
message("  y: ", y_col_A)

################--------------------------------------------
## 3. Compute Panel B from continuous scores and kNN graph
################--------------------------------------------

compute_knn_indices <- function(coords, k) {
  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knn(coords, k = k)
    return(nn$nn.index)
  }
  
  ## Fallback: base R distance matrix.
  ## For ~3,500 Visium spots this is acceptable.
  d <- as.matrix(stats::dist(coords))
  diag(d) <- Inf
  idx <- t(apply(d, 1, function(v) order(v)[seq_len(k)]))
  idx
}

symmetric_bivariate_moran <- function(x, y, coords, k) {
  xz <- as.numeric(scale(x))
  yz <- as.numeric(scale(y))
  keep <- is.finite(xz) & is.finite(yz) & is.finite(coords[, 1]) & is.finite(coords[, 2])
  
  xz <- xz[keep]
  yz <- yz[keep]
  coords <- coords[keep, , drop = FALSE]
  
  idx <- compute_knn_indices(coords, k = k)
  
  n <- length(xz)
  I_xy <- mean(vapply(seq_len(n), function(i) mean(xz[i] * yz[idx[i, ]]), numeric(1)))
  I_yx <- mean(vapply(seq_len(n), function(i) mean(yz[i] * xz[idx[i, ]]), numeric(1)))
  
  (I_xy + I_yx) / 2
}

coords_mat <- as.matrix(panelB_input[, c("coord_x", "coord_y")])

panelB_df <- data.frame(
  k = knn_values,
  I_symmetric = vapply(
    knn_values,
    function(k) symmetric_bivariate_moran(
      x = panelB_input$MyeloidTreg,
      y = panelB_input$DediffStromal,
      coords = coords_mat,
      k = k
    ),
    numeric(1)
  )
)

message("Panel B computed kNN symmetric bivariate Moran values:")
print(panelB_df)

################--------------------------------------------
## 4. Load and harmonize Panel C threshold-grid OR
################--------------------------------------------

panelC_raw <- read_any_table(panelC_file)

message("Panel C source table: ", panelC_file)
message("Panel C columns:")
message(paste(colnames(panelC_raw), collapse = ", "))

mode_col_C <- find_col(
  panelC_raw,
  patterns = c("^Mode$", "mode"),
  label = "Panel C mode"
)

top_fraction_col_C <- find_col(
  panelC_raw,
  patterns = c("^TopFraction$", "top.*fraction", "fraction", "quantile", "HighQuantile"),
  label = "Panel C top fraction"
)

or_col_C <- find_col(
  panelC_raw,
  patterns = c("^OR$", "^or$", "fisher.*or", "odds.*ratio", "odds_ratio", "ratio"),
  label = "Panel C Fisher odds ratio"
)

low_col_C <- find_col(
  panelC_raw,
  patterns = c("^CI_low$", "ci.*low", "lower", "or_low", "OR_low", "conf.low", "lcl", "lo95", "q025", "p025"),
  label = "Panel C lower CI",
  required = FALSE
)

high_col_C <- find_col(
  panelC_raw,
  patterns = c("^CI_high$", "ci.*high", "upper", "or_high", "OR_high", "conf.high", "ucl", "hi95", "q975", "p975"),
  label = "Panel C upper CI",
  required = FALSE
)

panelC_df <- panelC_raw %>%
  transmute(
    Mode_raw = as.character(.data[[mode_col_C]]),
    Mode = clean_mode_label(.data[[mode_col_C]]),
    TopFraction = as.numeric(.data[[top_fraction_col_C]]),
    Definition = top_fraction_label(.data[[top_fraction_col_C]]),
    OR = as.numeric(.data[[or_col_C]]),
    CI_low = if (!is.na(low_col_C)) as.numeric(.data[[low_col_C]]) else NA_real_,
    CI_high = if (!is.na(high_col_C)) as.numeric(.data[[high_col_C]]) else NA_real_
  ) %>%
  filter(is.finite(TopFraction), is.finite(OR), OR > 0) %>%
  mutate(
    CI_low = ifelse(is.finite(CI_low) & CI_low > 0, CI_low, OR),
    CI_high = ifelse(is.finite(CI_high) & CI_high > 0, CI_high, OR)
  )

if (any(is.na(panelC_df$Definition))) {
  stop("Some TopFraction values could not be converted into high-state definition labels.")
}

definition_levels <- panelC_df %>%
  distinct(TopFraction, Definition) %>%
  arrange(TopFraction) %>%
  pull(Definition)

mode_levels <- c("pooled", "per-sample")
mode_levels <- c(mode_levels[mode_levels %in% unique(panelC_df$Mode)],
                 setdiff(unique(panelC_df$Mode), mode_levels))

panelC_df <- panelC_df %>%
  mutate(
    Definition = factor(Definition, levels = definition_levels),
    Mode = factor(Mode, levels = mode_levels)
  ) %>%
  arrange(Definition, Mode)

if (nrow(panelC_df) == 0) {
  stop("Panel C harmonized table has zero rows.")
}

message("Panel C threshold labels after cleaning:")
print(as.data.frame(panelC_df[, c("Mode_raw", "Mode", "TopFraction", "Definition", "OR", "CI_low", "CI_high")]))

################--------------------------------------------
## 5. Plot Panel A
################--------------------------------------------

pA <- ggplot(panelA_df, aes(x = MyeloidTreg, y = DediffStromal)) +
  geom_point(size = 1.15, alpha = 0.55, color = "black") +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.75,
    color = "#3B82F6",
    fill = "grey80"
  ) +
  labs(
    title = "Continuous co-variation between spatial state scores",
    subtitle = "Visium spots; threshold-independent state-score relationship",
    x = "Myeloid–Treg immunosuppressive state score",
    y = "Tumor-dedifferentiation/stromal-remodeling state score"
  ) +
  theme_s16 +
  theme(
    plot.title = element_text(size = 13.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.3, hjust = 0)
  )

################--------------------------------------------
## 6. Plot Panel B
################--------------------------------------------

pB <- ggplot(panelB_df, aes(x = k, y = I_symmetric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_line(linewidth = 0.55, color = "black") +
  geom_point(size = 2.3, color = "black") +
  scale_x_continuous(
    breaks = sort(unique(panelB_df$k)),
    labels = function(x) as.character(as.integer(x))
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  labs(
    title = "Bivariate spatial co-variation across kNN definitions",
    subtitle = "Symmetric Moran-type statistic for continuous state scores",
    x = "k nearest neighbors",
    y = "Symmetric bivariate Moran's I"
  ) +
  theme_s16 +
  theme(
    plot.title = element_text(size = 13.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.3, hjust = 0)
  )

################--------------------------------------------
## 7. Plot Panel C
################--------------------------------------------

dodge <- position_dodge(width = 0.35)

pC <- ggplot(panelC_df, aes(x = Definition, y = OR, group = Mode)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high, linetype = Mode),
    width = 0.16,
    linewidth = 0.50,
    color = "black",
    position = dodge
  ) +
  geom_point(
    aes(shape = Mode),
    size = 2.7,
    color = "black",
    position = dodge
  ) +
  scale_shape_manual(
    values = c("pooled" = 16, "per-sample" = 1),
    name = "Sensitivity stratum"
  ) +
  scale_linetype_manual(
    values = c("pooled" = "solid", "per-sample" = "solid"),
    name = "Sensitivity stratum"
  ) +
  scale_y_log10(
    breaks = c(1, 3, 10),
    minor_breaks = NULL,
    labels = c("1", "3", "10"),
    expand = expansion(mult = c(0.04, 0.10))
  ) +
  labs(
    title = "Threshold-grid sensitivity of dual-high co-enrichment",
    subtitle = "Fisher odds ratios across TopFraction definitions; dashed line indicates odds ratio = 1",
    x = "High-state definition",
    y = "Fisher odds ratio (log10 scale)"
  ) +
  theme_s16 +
  theme(
    plot.title = element_text(size = 13.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9.3, hjust = 0),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.justification = "right"
  )

################--------------------------------------------
## 8. Assemble and save
################--------------------------------------------

top_row <- pA | pB

supp_s16 <- top_row / pC +
  plot_layout(heights = c(1.0, 0.96)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 19, face = "bold", family = "sans"),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

base_filename <- "Supplementary Figure S16. Threshold-independent spatial co-variation between myeloid–Treg and dedifferentiation-stromal scores in Visium"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 13.2
fig_height <- 9.4

safe_ggsave(out_png, supp_s16, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s16, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s16, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 9. Audit outputs
################--------------------------------------------

audit <- data.frame(
  panelA_file = panelA_file,
  panelC_file = panelC_file,
  panelB_source = "computed_from_panelA_spot_scores_and_coordinates",
  panelA_myeloid_col = myeloid_col_A,
  panelA_dediff_col = dediff_col_A,
  panelA_x_col = x_col_A,
  panelA_y_col = y_col_A,
  panelC_mode_col = mode_col_C,
  panelC_top_fraction_col = top_fraction_col_C,
  panelC_or_col = or_col_C,
  panelC_low_col = low_col_C,
  panelC_high_col = high_col_C,
  n_panelA_rows = nrow(panelA_df),
  n_panelB_spots_with_coordinates = nrow(panelB_input),
  n_panelB_rows = nrow(panelB_df),
  n_panelC_rows = nrow(panelC_df),
  k_values = paste(panelB_df$k, collapse = "; "),
  I_symmetric_values = paste(sprintf("%.6f", panelB_df$I_symmetric), collapse = "; "),
  threshold_labels = paste(as.character(unique(panelC_df$Definition)), collapse = "; "),
  sensitivity_strata = paste(as.character(unique(panelC_df$Mode)), collapse = "; "),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S16_threshold_independent_spatial_covariation_REVIEWER_FIXED_v3_300dpi_audit.csv")

write.csv(audit, audit_file, row.names = FALSE)
write.csv(panelA_df, file.path(fig_dir, "Supplementary_Figure_S16_panelA_continuous_score_scatter_harmonized_v3.csv"), row.names = FALSE)
write.csv(panelB_df, file.path(fig_dir, "Supplementary_Figure_S16_panelB_bivariate_Moran_computed_v3.csv"), row.names = FALSE)
write.csv(panelC_df, file.path(fig_dir, "Supplementary_Figure_S16_panelC_threshold_grid_OR_harmonized_v3.csv"), row.names = FALSE)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S16 script finished successfully.")