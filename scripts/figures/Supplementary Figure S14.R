################--------------------------------------------
## Supplementary Figure S14. Spatial QC-state correlation and QC-residualization sensitivity analysis
##
## Purpose:
##   1) Load spatial QC Spearman correlation and QC-residualized concordance tables from spatial_melanoma_validation
##   2) Generate and export Supplementary Figure S14 (Panel A QC correlation heatmap, Panel B raw vs QC-residualized concordance barplot)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S14. Spatial QC-state correlation and QC-residualization sensitivity analysis.png
##     - Supplementary Figure S14. Spatial QC-state correlation and QC-residualization sensitivity analysis.jpg
##     - Supplementary Figure S14. Spatial QC-state correlation and QC-residualization sensitivity analysis.pdf
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

################--------------------------------------------
## 0. Paths and settings
################--------------------------------------------

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

manual_qc_corr_file <- file.path(
  project_dir,
  "results", "tables", "spatial_melanoma_validation",
  "Step11_spatial_QC_state_Spearman_correlation.csv"
)

manual_qc_concordance_file <- file.path(
  project_dir,
  "results", "tables", "spatial_melanoma_validation",
  "Step11_spatial_raw_vs_QCresidual_state_correlation.csv"
)

HEATMAP_LOW  <- "#3B82F6"
HEATMAP_MID  <- "#FFFFFF"
HEATMAP_HIGH <- "#EF4444"

state_order <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation"
)

qc_metric_order <- c("nFeature_Spatial", "nCount_Spatial", "percent.mt")

################--------------------------------------------
## 1. Helper functions
################--------------------------------------------

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

theme_s14 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 15.0, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 10.8, hjust = 0, color = "black"),
    axis.title = element_text(size = 12.0, face = "bold", color = "black"),
    axis.text = element_text(size = 10.6, color = "black"),
    legend.title = element_text(size = 10.8, face = "bold", color = "black"),
    legend.text = element_text(size = 9.4, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.34),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 6, 6, 6)
  )

read_any_table <- function(file) {
  if (!file.exists(file)) {
    stop(
      "Required Supplementary Figure S14 input file does not exist:\n",
      file,
      "\nPlease confirm Step11 spatial QC outputs were generated."
    )
  }
  suppressMessages(readr::read_csv(file, show_col_types = FALSE, guess_max = 100000)) %>%
    as.data.frame(stringsAsFactors = FALSE)
}

find_col <- function(df, patterns, label, required = TRUE) {
  cols <- colnames(df)
  for (pat in patterns) {
    hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE)
    if (length(hit) > 0) return(hit[1])
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

standardize_state <- function(x) {
  x0 <- as.character(x)
  x1 <- x0
  x1 <- gsub("_", " ", x1)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  x1_low <- tolower(x1)
  
  dplyr::case_when(
    grepl("immune.*defective|defective.*cold|immune.*cold|cold", x1_low) ~ "immune-defective/cold",
    grepl("myeloid.*treg|treg.*immunosuppress|myeloid treg|myeloid.treg", x1_low) ~ "myeloid–Treg immunosuppressive",
    grepl("tumor.*dediff|dediff.*stromal|stromal.*remodel|dediff/stromal|dediff.*stroma", x1_low) ~
      "tumor-dedifferentiation/stromal-remodeling",
    grepl("melanocytic|melanocyte|melanoma.*differentiation|lineage", x1_low) ~ "melanocytic differentiation",
    TRUE ~ x0
  )
}

standardize_qc_metric <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", ".", x0)
  x1_low <- tolower(x1)
  
  dplyr::case_when(
    grepl("nfeature|feature", x1_low) ~ "nFeature_Spatial",
    grepl("ncount|count|umi", x1_low) ~ "nCount_Spatial",
    grepl("percent.*mt|percent.mt|mito|mitochond", x1_low) ~ "percent.mt",
    TRUE ~ x0
  )
}

wrap_state <- function(x) {
  dplyr::recode(
    as.character(x),
    "immune-defective/cold" = "immune-defective/cold",
    "myeloid–Treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
    "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
    "melanocytic differentiation" = "melanocytic\ndifferentiation",
    .default = as.character(x)
  )
}

################--------------------------------------------
## 2. Read correct inputs
################--------------------------------------------

qc_corr_file <- manual_qc_corr_file
qc_concordance_file <- manual_qc_concordance_file

qc_corr_raw <- read_any_table(qc_corr_file)
qc_conc_raw <- read_any_table(qc_concordance_file)

message("Panel A table: ", qc_corr_file)
message("Panel A columns: ", paste(colnames(qc_corr_raw), collapse = ", "))
message("Panel B table: ", qc_concordance_file)
message("Panel B columns: ", paste(colnames(qc_conc_raw), collapse = ", "))

################--------------------------------------------
## 3. Harmonize Panel A
################--------------------------------------------

state_col_A <- find_col(
  qc_corr_raw,
  patterns = c("^State$", "^state$", "state_name", "StateLabel", "state_label", "final.*state", "tumor.*immune.*state", "signature", "score"),
  label = "Panel A state"
)

metric_col_A <- find_col(
  qc_corr_raw,
  patterns = c("^QC_metric$", "^qc_metric$", "^metric$", "qc.*metric", "spatial.*metric", "variable", "feature", "covariate"),
  label = "Panel A QC metric"
)

rho_col_A <- find_col(
  qc_corr_raw,
  patterns = c("^rho$", "spearman.*rho", "spearman", "correlation", "^cor$", "estimate"),
  label = "Panel A Spearman rho"
)

panelA_df <- qc_corr_raw %>%
  transmute(
    State = standardize_state(.data[[state_col_A]]),
    QCMetric = standardize_qc_metric(.data[[metric_col_A]]),
    SpearmanRho = as.numeric(.data[[rho_col_A]])
  ) %>%
  filter(State %in% state_order, QCMetric %in% qc_metric_order, is.finite(SpearmanRho)) %>%
  group_by(State, QCMetric) %>%
  summarise(SpearmanRho = SpearmanRho[1], .groups = "drop") %>%
  mutate(
    State = factor(State, levels = rev(state_order)),
    QCMetric = factor(QCMetric, levels = qc_metric_order),
    label = sprintf("%.2f", SpearmanRho)
  )

if (nrow(panelA_df) == 0) {
  stop("Panel A harmonized table has zero rows.")
}

################--------------------------------------------
## 4. Harmonize Panel B
################--------------------------------------------

state_col_B <- find_col(
  qc_conc_raw,
  patterns = c("^State$", "^state$", "state_name", "StateLabel", "state_label", "final.*state", "tumor.*immune.*state", "signature", "score"),
  label = "Panel B state"
)

rho_col_B <- find_col(
  qc_conc_raw,
  patterns = c("^rho$", "raw.*resid.*rho", "raw.*qc.*rho", "concordance", "spearman.*rho", "spearman", "correlation", "^cor$", "estimate"),
  label = "Panel B raw vs QC-residualized Spearman rho"
)

panelB_df <- qc_conc_raw %>%
  transmute(
    State = standardize_state(.data[[state_col_B]]),
    SpearmanRho = as.numeric(.data[[rho_col_B]])
  ) %>%
  filter(State %in% state_order, is.finite(SpearmanRho)) %>%
  group_by(State) %>%
  summarise(SpearmanRho = SpearmanRho[1], .groups = "drop") %>%
  mutate(
    State = factor(State, levels = rev(state_order)),
    label = sprintf("%.2f", SpearmanRho)
  )

if (nrow(panelB_df) == 0) {
  stop("Panel B harmonized table has zero rows.")
}

message("Panel B harmonized values:")
print(as.data.frame(panelB_df), row.names = FALSE)

################--------------------------------------------
## 5. Plot Panel A
################--------------------------------------------

pA <- ggplot(panelA_df, aes(x = QCMetric, y = State, fill = SpearmanRho)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 4.1, color = "black") +
  scale_fill_gradient2(
    low = HEATMAP_LOW,
    mid = HEATMAP_MID,
    high = HEATMAP_HIGH,
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  scale_x_discrete(labels = c(
    "nFeature_Spatial" = "nFeature_Spatial",
    "nCount_Spatial" = "nCount_Spatial",
    "percent.mt" = "percent.mt"
  )) +
  scale_y_discrete(labels = wrap_state) +
  labs(
    title = "State-score correlations with spatial QC metrics",
    x = NULL,
    y = NULL
  ) +
  theme_s14 +
  theme(
    plot.title = element_text(size = 15.0, face = "bold", hjust = 0),
    axis.text.x = element_text(size = 10.8, angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10.8),
    legend.title = element_text(size = 10.8, face = "bold"),
    legend.text = element_text(size = 9.6),
    panel.grid = element_blank()
  )

################--------------------------------------------
## 6. Plot Panel B
################--------------------------------------------

pB <- ggplot(panelB_df, aes(x = SpearmanRho, y = State)) +
  geom_col(width = 0.58, fill = "grey35", color = "grey20", linewidth = 0.25) +
  geom_text(aes(label = label), hjust = -0.15, size = 3.9, color = "black") +
  scale_y_discrete(labels = wrap_state) +
  scale_x_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1.0, 0.25),
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Raw versus QC-residualized\nstate-score concordance",
    x = "Spearman rho",
    y = NULL
  ) +
  theme_s14 +
  theme(
    plot.title = element_text(size = 14.8, face = "bold", hjust = 0),
    axis.text.x = element_text(size = 10.8),
    axis.text.y = element_text(size = 10.8),
    axis.title.x = element_text(size = 12.0, face = "bold")
  )

################--------------------------------------------
## 7. Assemble and save
################--------------------------------------------

supp_s14 <- pA | pB
supp_s14 <- supp_s14 +
  plot_layout(widths = c(1.05, 1.12)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 21, face = "bold", family = "sans"),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

base_filename <- "Supplementary Figure S14. Spatial QC-state correlation and QC-residualization sensitivity analysis"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 14.2
fig_height <- 5.9

safe_ggsave(out_png, supp_s14, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s14, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s14, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

cat("\nSupplementary Figure S14 generated successfully with specified file paths!\n")