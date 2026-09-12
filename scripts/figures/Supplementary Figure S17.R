################--------------------------------------------
## Supplementary Figure S17. Composition-aware analysis of Visium dual-high tissue regions
##
## Purpose:
##   1) Load spatial spot-level composition-residualized scores and Fisher OR tables from revision_visium_composition
##   2) Generate and export Supplementary Figure S17 (Panel A marker modules by category, Panel B state-module correlations, Panel C residualized OR, Panel D spatial categories)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S17. Composition-aware analysis of Visium dual-high tissue regions.png
##     - Supplementary Figure S17. Composition-aware analysis of Visium dual-high tissue regions.jpg
##     - Supplementary Figure S17. Composition-aware analysis of Visium dual-high tissue regions.pdf
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
## 0. Paths
################--------------------------------------------

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

spot_file <- file.path(
  project_dir,
  "results", "tables", "revision_visium_composition",
  "18E_Visium_spot_state_composition_residualized_scores.csv"
)

or_file <- file.path(
  project_dir,
  "results", "tables", "revision_visium_composition",
  "18E_raw_vs_residualized_dual_high_Fisher_OR.csv"
)

################--------------------------------------------
## 1. Helpers
################--------------------------------------------

read_any_table <- function(file) {
  if (!file.exists(file)) stop("File does not exist: ", file)
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("csv", "txt")) {
    out <- suppressMessages(readr::read_csv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "tsv") {
    out <- suppressMessages(readr::read_tsv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "rds") {
    out <- readRDS(file)
    if (!is.data.frame(out)) stop("RDS object is not a data.frame: ", file)
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

standardize_state <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("/", " ", x1)
  x1 <- gsub("\\blike\\b", " ", x1, ignore.case = TRUE)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("immune.*defective|defective.*cold|immune.*cold|cold", xl) ~ "immune-defective/cold",
    grepl("myeloid.*treg|treg.*immunosuppress|myeloid treg|myeloid.treg", xl) ~ "myeloid–Treg immunosuppressive",
    grepl("tumor.*dediff|dediff.*stromal|stromal.*remodel|dediff.*stroma", xl) ~ "tumor-dedifferentiation/stromal-remodeling",
    grepl("melanocytic|melanocyte|melanoma.*differentiation|lineage", xl) ~ "melanocytic differentiation",
    TRUE ~ x0
  )
}

standardize_category <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("/", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("both", xl) ~ "both high",
    grepl("neither|none|no high|low", xl) ~ "neither high",
    grepl("myeloid.*treg", xl) & !grepl("dediff|stromal|both", xl) ~ "myeloid–Treg high only",
    grepl("dediff|stromal", xl) & !grepl("myeloid|treg|both", xl) ~ "tumor-dedifferentiation/stromal-remodeling high only",
    grepl("myeloid.*only|treg.*only", xl) ~ "myeloid–Treg high only",
    grepl("dediff.*only|stromal.*only", xl) ~ "tumor-dedifferentiation/stromal-remodeling high only",
    TRUE ~ x0
  )
}

derive_category <- function(myeloid_high, dediff_high) {
  dplyr::case_when(
    myeloid_high & dediff_high ~ "both high",
    myeloid_high & !dediff_high ~ "myeloid–Treg high only",
    !myeloid_high & dediff_high ~ "tumor-dedifferentiation/stromal-remodeling high only",
    TRUE ~ "neither high"
  )
}

standardize_method <- function(x) {
  x0 <- as.character(x)
  xl <- tolower(gsub("[_\\-]+", " ", x0))
  dplyr::case_when(
    grepl("qc", xl) & grepl("composition|comp", xl) ~ "QC/composition-residualized state scores",
    grepl("composition|comp", xl) & grepl("resid|residual", xl) ~ "composition-residualized state scores",
    grepl("qc", xl) & grepl("resid|residual", xl) ~ "QC-residualized state scores",
    grepl("raw", xl) ~ "raw state scores",
    TRUE ~ x0
  )
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

wrap_category <- function(x) {
  dplyr::recode(
    as.character(x),
    "neither high" = "neither high",
    "myeloid–Treg high only" = "myeloid–Treg\nhigh only",
    "tumor-dedifferentiation/stromal-remodeling high only" = "dediff/stromal\nhigh only",
    "both high" = "both high",
    .default = as.character(x)
  )
}

wrap_method <- function(x) {
  dplyr::recode(
    as.character(x),
    "raw state scores" = "raw state scores",
    "QC-residualized state scores" = "QC-residualized\nstate scores",
    "composition-residualized state scores" = "composition-residualized\nstate scores",
    "QC/composition-residualized state scores" = "QC/composition-residualized\nstate scores",
    .default = as.character(x)
  )
}

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

category_order <- c(
  "neither high",
  "myeloid–Treg high only",
  "tumor-dedifferentiation/stromal-remodeling high only",
  "both high"
)

category_colors <- c(
  "neither high" = "#E64B35",
  "myeloid–Treg high only" = "#00A087",
  "tumor-dedifferentiation/stromal-remodeling high only" = "#4DBBD5",
  "both high" = "#B26BFF"
)

theme_s17 <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 12.8, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 9.2, hjust = 0, color = "black"),
    axis.title = element_text(size = 10.4, face = "bold", color = "black"),
    axis.text = element_text(size = 8.8, color = "black"),
    strip.text = element_text(size = 8.8, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey85", color = "grey45", linewidth = 0.30),
    legend.title = element_text(size = 9.2, face = "bold", color = "black"),
    legend.text = element_text(size = 8.0, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

################--------------------------------------------
## 2. Load inputs
################--------------------------------------------

spot_df <- read_any_table(spot_file)
or_df <- read_any_table(or_file)

message("Spot-level table: ", spot_file)
message("Spot-level columns: ", paste(colnames(spot_df), collapse = ", "))
message("OR table: ", or_file)
message("OR columns: ", paste(colnames(or_df), collapse = ", "))

required_spot_cols <- c(
  "spatial_x",
  "spatial_y",
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation",
  "CAF_Stromal",
  "Myeloid",
  "T_NK",
  "B_Plasma",
  "Melanoma_Lineage",
  "Endothelial",
  "Raw_DualHighCategory",
  "MyeloidTreg_resid_QC_Composition",
  "DediffStromal_resid_QC_Composition"
)

require_cols(spot_df, required_spot_cols, "Spot-level table")

required_or_cols <- c("Analysis", "OR", "CI_low", "CI_high")
require_cols(or_df, required_or_cols, "OR table")

################--------------------------------------------
## 3. Panel A: composition module boxplots by raw category
################--------------------------------------------

module_col_map <- c(
  "B/plasma" = "B_Plasma",
  "CAF/stromal" = "CAF_Stromal",
  "endothelial" = "Endothelial",
  "melanoma lineage" = "Melanoma_Lineage",
  "myeloid" = "Myeloid",
  "T/NK" = "T_NK"
)

panelA_df <- spot_df %>%
  transmute(
    Category = standardize_category(Raw_DualHighCategory),
    !!!setNames(lapply(unname(module_col_map), function(cc) as.numeric(spot_df[[cc]])), names(module_col_map))
  ) %>%
  filter(Category %in% category_order) %>%
  pivot_longer(
    cols = all_of(names(module_col_map)),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(
    Module = factor(Module, levels = c("B/plasma", "CAF/stromal", "endothelial", "melanoma lineage", "myeloid", "T/NK")),
    Category = factor(Category, levels = category_order)
  ) %>%
  filter(is.finite(Score), !is.na(Module), !is.na(Category))

if (nrow(panelA_df) == 0) stop("Panel A harmonized table has zero rows.")

################--------------------------------------------
## 4. Panel B: state-score/module Spearman heatmap
################--------------------------------------------

state_col_map <- c(
  "melanocytic differentiation" = "Melanocytic_Differentiation",
  "tumor-dedifferentiation/stromal-remodeling" = "Tumor_dedifferentiation_Stromal_remodeling",
  "myeloid–Treg immunosuppressive" = "Myeloid_Treg_Immunosuppressive",
  "immune-defective/cold" = "Immune_defective_Cold"
)

panelB_long <- expand.grid(
  State = names(state_col_map),
  Module = c("CAF/stromal", "myeloid", "T/NK", "B/plasma", "melanoma lineage", "endothelial"),
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

if (nrow(panelB_long) == 0) stop("Panel B correlation table has zero rows.")

################--------------------------------------------
## 5. Panel C: OR table
################--------------------------------------------

panelC_df <- or_df %>%
  transmute(
    Method_raw = as.character(Analysis),
    Method = standardize_method(Analysis),
    OR = as.numeric(OR),
    CI_low = as.numeric(CI_low),
    CI_high = as.numeric(CI_high)
  ) %>%
  filter(is.finite(OR), OR > 0) %>%
  mutate(
    CI_low = ifelse(is.finite(CI_low) & CI_low > 0, CI_low, OR),
    CI_high = ifelse(is.finite(CI_high) & CI_high > 0, CI_high, OR),
    Method = factor(
      Method,
      levels = c(
        "raw state scores",
        "QC-residualized state scores",
        "composition-residualized state scores",
        "QC/composition-residualized state scores"
      )
    )
  ) %>%
  filter(!is.na(Method)) %>%
  arrange(Method)

if (nrow(panelC_df) == 0) {
  stop(
    "Panel C harmonized table has zero rows. Inspect unique(or_df$Analysis):\n",
    paste(unique(or_df$Analysis), collapse = "; ")
  )
}

################--------------------------------------------
## 6. Panel D: raw and derived residualized spatial categories
################--------------------------------------------

cutoff_q <- NA_real_
if ("cutoff_quantile" %in% colnames(or_df)) {
  qc_row <- grep("qc.*comp|comp.*qc|qc.*composition|composition.*qc", or_df$Analysis, ignore.case = TRUE)
  if (length(qc_row) > 0) {
    cutoff_q <- suppressWarnings(as.numeric(or_df$cutoff_quantile[qc_row[1]]))
  }
  if (!is.finite(cutoff_q)) {
    cutoff_q <- suppressWarnings(as.numeric(or_df$cutoff_quantile[1]))
  }
}
if (!is.finite(cutoff_q) || cutoff_q <= 0 || cutoff_q >= 1) {
  cutoff_q <- 0.75
}

mt_resid <- as.numeric(spot_df$MyeloidTreg_resid_QC_Composition)
ds_resid <- as.numeric(spot_df$DediffStromal_resid_QC_Composition)
mt_cut <- as.numeric(stats::quantile(mt_resid, probs = cutoff_q, na.rm = TRUE, names = FALSE))
ds_cut <- as.numeric(stats::quantile(ds_resid, probs = cutoff_q, na.rm = TRUE, names = FALSE))

panelD_df <- spot_df %>%
  transmute(
    x = as.numeric(spatial_x),
    y = as.numeric(spatial_y),
    `QC/composition-residualized state scores` = derive_category(
      as.numeric(MyeloidTreg_resid_QC_Composition) >= mt_cut,
      as.numeric(DediffStromal_resid_QC_Composition) >= ds_cut
    ),
    `raw state scores` = standardize_category(Raw_DualHighCategory)
  ) %>%
  filter(is.finite(x), is.finite(y)) %>%
  pivot_longer(
    cols = c("QC/composition-residualized state scores", "raw state scores"),
    names_to = "Score type",
    values_to = "Category"
  ) %>%
  filter(Category %in% category_order) %>%
  mutate(
    Category = factor(Category, levels = category_order),
    `Score type` = factor(
      `Score type`,
      levels = c("QC/composition-residualized state scores", "raw state scores")
    )
  )

if (nrow(panelD_df) == 0) stop("Panel D harmonized table has zero rows.")

################--------------------------------------------
## 7. Plots (With compact labels & safe margins)
################--------------------------------------------

pA <- ggplot(panelA_df, aes(x = Category, y = Score)) +
  geom_boxplot(width = 0.62, outlier.size = 0.30, outlier.alpha = 0.45, linewidth = 0.35) +
  facet_wrap(~ Module, ncol = 3, scales = "free_y") +
  scale_x_discrete(labels = wrap_category) +
  labs(
    title = "Composition marker modules across Visium high-state categories",
    subtitle = "Marker-module scores provide composition-aware context for mixed-cell Visium spots",
    x = "Raw high-state category",
    y = "Marker-module score"
  ) +
  theme_s17 +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.5),
    strip.text = element_text(size = 8.1, face = "bold"),
    plot.margin = margin(5, 5, 12, 5)
  )

pB <- ggplot(panelB_long, aes(x = Module, y = State, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 2.6, color = "black") +
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
    title = "Visium state-score associations with composition marker modules",
    subtitle = "Spearman correlations across spatial spots",
    x = "Composition marker module",
    y = "Tumor–immune state"
  ) +
  theme_s17 +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  )

pC <- ggplot(panelC_df, aes(y = Method, x = OR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.18,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_point(size = 2.2, color = "black") +
  scale_x_log10(
    breaks = c(1, 3, 10),
    labels = c("1", "3", "10"),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  scale_y_discrete(labels = wrap_method) +
  labs(
    title = "Raw and residualized dual-high co-enrichment in Visium",
    subtitle = "Composition-aware sensitivity analysis; dashed line indicates odds ratio = 1",
    x = "Fisher odds ratio (log10 scale)",
    y = NULL
  ) +
  theme_s17 +
  theme(axis.text.y = element_text(size = 8.4))

pD <- ggplot(panelD_df, aes(x = x, y = y, color = Category)) +
  geom_point(size = 0.42, alpha = 0.90) +
  coord_equal() +
  scale_y_reverse() +
  facet_wrap(~ `Score type`, nrow = 1) +
  scale_color_manual(
    values = category_colors,
    breaks = category_order,
    labels = wrap_category,
    name = "Category",
    drop = FALSE
  ) +
  labs(
    title = "Raw and composition-residualized high-state categories",
    subtitle = paste0(
      "Top-quantile operational categories before and after QC/composition residualization; cutoff quantile = ",
      sprintf("%.2f", cutoff_q)
    ),
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_s17 +
  theme(
    strip.text = element_text(size = 8.2, face = "bold"),
    axis.title = element_text(size = 9.2, face = "bold"),
    axis.text = element_text(size = 7.4),
    legend.position = "right",
    legend.text = element_text(size = 6.8),
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.35, "cm")
  )

################--------------------------------------------
## 8. Assemble and save (Expanded width to prevent crowding)
################--------------------------------------------

supp_s17 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.02, 0.98)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 19, face = "bold", family = "sans"),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

base_filename <- "Supplementary Figure S17. Composition-aware analysis of Visium dual-high tissue regions"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 17.2  # 适当放大画幅宽度
fig_height <- 10.8

safe_ggsave(out_png, supp_s17, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s17, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s17, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 9. Audit outputs
################--------------------------------------------

audit <- data.frame(
  spot_file = spot_file,
  or_file = or_file,
  n_spots = nrow(spot_df),
  n_panelA_rows = nrow(panelA_df),
  n_panelB_rows = nrow(panelB_long),
  n_panelC_rows = nrow(panelC_df),
  n_panelD_rows = nrow(panelD_df),
  cutoff_quantile_used_for_panelD = cutoff_q,
  MyeloidTreg_resid_QC_Composition_cutoff = mt_cut,
  DediffStromal_resid_QC_Composition_cutoff = ds_cut,
  panelC_methods = paste(as.character(panelC_df$Method), collapse = "; "),
  category_levels = paste(category_order, collapse = "; "),
  module_levels = paste(module_order, collapse = "; "),
  state_levels = paste(state_order, collapse = "; "),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S17_composition_aware_Visium_dual_high_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

write.csv(panelA_df, file.path(fig_dir, "Supplementary_Figure_S17_panelA_composition_modules_by_category_harmonized.csv"), row.names = FALSE)
write.csv(panelB_long, file.path(fig_dir, "Supplementary_Figure_S17_panelB_state_module_correlations_harmonized.csv"), row.names = FALSE)
write.csv(panelC_df, file.path(fig_dir, "Supplementary_Figure_S17_panelC_residualized_OR_harmonized.csv"), row.names = FALSE)
write.csv(panelD_df, file.path(fig_dir, "Supplementary_Figure_S17_panelD_spatial_categories_harmonized.csv"), row.names = FALSE)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S17 script finished successfully.")