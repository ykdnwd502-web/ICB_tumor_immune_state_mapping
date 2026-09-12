################################----------------------------
## Supplementary Figure S13. Formal composition-adjusted response analysis and
## leave-one-sample-out stability assessment in bulk ICB cohorts
##
## Purpose:
##   1) Load composition-adjusted AUC summary and LOSO stability tables from bulk_response_sensitivity
##   2) Generate and export Supplementary Figure S13 with adjusted font sizing
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S13. Formal composition-adjusted response analysis and leave-one-sample-out stability assessment in bulk ICB cohorts.png
##     - Supplementary Figure S13. Formal composition-adjusted response analysis and leave-one-sample-out stability assessment in bulk ICB cohorts.jpg
##     - Supplementary Figure S13. Formal composition-adjusted response analysis and leave-one-sample-out stability assessment in bulk ICB cohorts.pdf
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

table_root <- file.path(project_dir, "results", "tables")
fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

preferred_auc_file  <- file.path(project_dir, "results", "tables", "bulk_response_sensitivity", "S13_formal_composition_adjusted_auc_summary.csv")
preferred_loso_file <- file.path(project_dir, "results", "tables", "bulk_response_sensitivity", "S13_loso_auc_values.csv")

manual_auc_file  <- ""
manual_loso_file <- ""

target_state_label <- "tumor-dedifferentiation/stromal-remodeling"
require_composition_pc <- TRUE

################--------------------------------------------
## 1. Helper functions & Updated Theme
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

## 除了标题和亚标题所有文字放大15%
theme_s13 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14.4, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 11.0, hjust = 0, color = "black"),
    axis.title = element_text(size = 12.0 * 1.15, face = "bold", color = "black"),
    axis.text = element_text(size = 10.6 * 1.15, color = "black"),
    strip.background = element_rect(fill = "#D9D9D9", color = "grey45", linewidth = 0.35),
    strip.text = element_text(size = 9.8 * 1.15, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.34),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(6, 6, 6, 6)
  )

read_any_table <- function(file) {
  if (!file.exists(file)) stop("File does not exist: ", file)
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("csv", "txt")) {
    out <- suppressMessages(readr::read_csv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "tsv") {
    out <- suppressMessages(readr::read_tsv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "rds") {
    out <- readRDS(file)
    if (!is.data.frame(out)) stop("RDS is not a data.frame: ", file)
  } else {
    stop("Unsupported file extension: ", file)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

find_col <- function(df, patterns, label, required = TRUE) {
  cols <- colnames(df)
  for (pat in patterns) {
    hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  if (required) {
    stop("Could not find column for ", label, ".\nAvailable columns:\n", paste(cols, collapse = ", "))
  }
  NA_character_
}

first_finite_or_na <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

standardize_cohort <- function(x) {
  x0 <- as.character(x)
  x1 <- tolower(x0)
  dplyr::case_when(
    grepl("78220", x1) ~ "GSE78220",
    grepl("91061", x1) ~ "GSE91061",
    TRUE ~ x0
  )
}

standardize_state <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  x1_low <- tolower(x1)
  
  dplyr::case_when(
    grepl("tumor.*dediff|dediff.*stromal|stromal.*remodel", x1_low) ~ "tumor-dedifferentiation/stromal-remodeling",
    grepl("immune.*defective|defective.*cold|immune.*cold", x1_low) ~ "immune-defective/cold",
    grepl("myeloid.*treg|treg.*immunosuppress", x1_low) ~ "myeloid–Treg immunosuppressive",
    grepl("melanocytic|melanocyte", x1_low) ~ "melanocytic differentiation",
    TRUE ~ x0
  )
}

standardize_method <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  x1_low <- tolower(x1)
  
  dplyr::case_when(
    grepl("^raw$|raw score|unadjusted|full sample raw|^none$|no adjustment", x1_low) ~ "Raw",
    grepl("caf.*stromal|stromal.*caf|^caf$|caf/stromal", x1_low) ~ "CAF/stromal",
    grepl("caf.*melanoma|caf.*lineage", x1_low) ~ "CAF + melanoma-lineage",
    grepl("composition.*pc|pc.*composition|principal|pca|pc1|pc2|marker.*pc|all.*pc|pc.*adjust|pc.*residual", x1_low) ~ "Composition PC",
    grepl("all.*marker|marker.*comp|all.*composition|all.*module", x1_low) ~ "All marker composition",
    grepl("purity", x1_low) ~ "Purity proxy",
    grepl("mcp", x1_low) ~ "MCP-counter",
    TRUE ~ x0
  )
}

method_order <- c(
  "Raw",
  "CAF/stromal",
  "CAF + melanoma-lineage",
  "All marker composition",
  "Purity proxy",
  "Composition PC",
  "MCP-counter"
)

################--------------------------------------------
## 2. Load tables
################--------------------------------------------

auc_file <- if (nzchar(manual_auc_file)) manual_auc_file else preferred_auc_file
loso_file <- if (nzchar(manual_loso_file)) manual_loso_file else preferred_loso_file

auc_raw <- read_any_table(auc_file)

################--------------------------------------------
## 3. Harmonize AUC summary table
################--------------------------------------------

cohort_col <- find_col(auc_raw, c("^cohort$", "^dataset$", "^Dataset$", "^Cohort$", "source", "GSE", "Study"), "cohort")
state_col <- find_col(auc_raw, c("^state$", "^State$", "state_name", "final.*state", "signature"), "state")
outcome_col <- find_col(auc_raw, c("^outcome_mode$", "outcome", "endpoint", "response_mode"), "outcome mode", required = FALSE)
method_col <- find_col(auc_raw, c("^adjustment$", "adjustment", "AdjustmentSet", "Model", "model", "Method", "method", "Predictor", "predictor", "Analysis", "analysis"), "adjustment/model")
auc_col <- find_col(auc_raw, c("^auc$", "^AUC$", "AUC_raw", "AUC_adjusted", "AUC_full", "full.*auc", "observed.*auc"), "auc")
low_col <- find_col(auc_raw, c("^bootstrap_ci_low$", "bootstrap.*low", "ci.*low", "lower", "auc_low"), "lower CI", required = FALSE)
high_col <- find_col(auc_raw, c("^bootstrap_ci_high$", "bootstrap.*high", "ci.*high", "upper", "auc_high"), "upper CI", required = FALSE)
loso_med_col <- find_col(auc_raw, c("^loso_auc_median$", "loso.*median", "median.*auc"), "LOSO median", required = FALSE)
loso_min_col <- find_col(auc_raw, c("^loso_auc_min$", "loso.*min", "min.*auc"), "LOSO min", required = FALSE)
loso_max_col <- find_col(auc_raw, c("^loso_auc_max$", "loso.*max", "max.*auc"), "LOSO max", required = FALSE)

auc_long <- auc_raw %>%
  transmute(
    Cohort = standardize_cohort(.data[[cohort_col]]),
    OutcomeMode = if (!is.na(outcome_col)) as.character(.data[[outcome_col]]) else NA_character_,
    State = standardize_state(.data[[state_col]]),
    Method = standardize_method(.data[[method_col]]),
    AUC = as.numeric(.data[[auc_col]]),
    CI_low = if (!is.na(low_col)) as.numeric(.data[[low_col]]) else NA_real_,
    CI_high = if (!is.na(high_col)) as.numeric(.data[[high_col]]) else NA_real_,
    LOSO_median = if (!is.na(loso_med_col)) as.numeric(.data[[loso_med_col]]) else NA_real_,
    LOSO_min = if (!is.na(loso_min_col)) as.numeric(.data[[loso_min_col]]) else NA_real_,
    LOSO_max = if (!is.na(loso_max_col)) as.numeric(.data[[loso_max_col]]) else NA_real_
  )

target_rows <- auc_long %>%
  filter(
    Cohort %in% c("GSE78220", "GSE91061"),
    State == target_state_label,
    Method %in% method_order,
    is.finite(AUC)
  )

if (nrow(target_rows) == 0) {
  stop("No rows remained after filtering to target state: ", target_state_label)
}

auc_df <- target_rows %>%
  group_by(Cohort, Method) %>%
  summarise(
    AUC = first_finite_or_na(AUC),
    CI_low = first_finite_or_na(CI_low),
    CI_high = first_finite_or_na(CI_high),
    LOSO_median = first_finite_or_na(LOSO_median),
    LOSO_min = first_finite_or_na(LOSO_min),
    LOSO_max = first_finite_or_na(LOSO_max),
    .groups = "drop"
  ) %>%
  mutate(
    CI_low = ifelse(is.finite(CI_low), CI_low, AUC),
    CI_high = ifelse(is.finite(CI_high), CI_high, AUC)
  )

method_order_display <- method_order[method_order %in% as.character(unique(auc_df$Method))]
auc_df <- auc_df %>%
  mutate(Method = factor(as.character(Method), levels = rev(method_order_display)))

################--------------------------------------------
## 4. LOSO data for Panel C
################--------------------------------------------

if (all(c(!is.na(loso_med_col), !is.na(loso_min_col), !is.na(loso_max_col)))) {
  loso_df <- auc_df %>%
    transmute(
      Cohort,
      Method,
      LOSO_min = LOSO_min,
      LOSO_median = LOSO_median,
      LOSO_max = LOSO_max,
      Full_AUC = AUC
    ) %>%
    filter(is.finite(LOSO_min), is.finite(LOSO_median), is.finite(LOSO_max), is.finite(Full_AUC))
} else {
  loso_raw <- read_any_table(loso_file)
  cohort_col_loso <- find_col(loso_raw, c("^cohort$", "^dataset$", "^Dataset$", "^Cohort$", "source", "GSE", "Study"), "LOSO cohort")
  state_col_loso <- find_col(loso_raw, c("^state$", "^State$", "state_name"), "LOSO state")
  method_col_loso <- find_col(loso_raw, c("^adjustment$", "adjustment", "Model", "model", "Method", "method"), "LOSO adjustment")
  loso_auc_col <- find_col(loso_raw, c("LeaveOneOut_AUC", "LOSO_AUC", "leave.*auc", "loso.*auc"), "LOSO AUC")
  
  loso_df <- loso_raw %>%
    transmute(
      Cohort = standardize_cohort(.data[[cohort_col_loso]]),
      State = standardize_state(.data[[state_col_loso]]),
      Method = standardize_method(.data[[method_col_loso]]),
      LOSO_AUC = as.numeric(.data[[loso_auc_col]])
    ) %>%
    filter(Cohort %in% c("GSE78220", "GSE91061"), State == target_state_label, Method %in% method_order, is.finite(LOSO_AUC)) %>%
    group_by(Cohort, Method) %>%
    summarise(
      LOSO_min = min(LOSO_AUC, na.rm = TRUE),
      LOSO_median = median(LOSO_AUC, na.rm = TRUE),
      LOSO_max = max(LOSO_AUC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(auc_df %>% select(Cohort, Method, Full_AUC = AUC), by = c("Cohort", "Method")) %>%
    filter(is.finite(LOSO_min), is.finite(LOSO_median), is.finite(LOSO_max), is.finite(Full_AUC))
}

loso_df <- loso_df %>%
  mutate(Method = factor(as.character(Method), levels = rev(method_order_display)))

################--------------------------------------------
## 5. Plotting panels (Panel C strip.text 放大25%)
################--------------------------------------------

plot_auc_forest <- function(df, title, subtitle = NULL) {
  df <- df %>% filter(!is.na(Method), is.finite(AUC))
  ggplot(df, aes(x = AUC, y = Method)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.42, color = "black") +
    geom_segment(aes(x = CI_low, xend = CI_high, y = Method, yend = Method), linewidth = 0.62, color = "black") +
    geom_point(size = 2.4, color = "black") +
    scale_x_continuous(limits = c(0.16, 1.02), breaks = seq(0.2, 1.0, 0.2), labels = number_format(accuracy = 0.01)) +
    labs(title = title, subtitle = subtitle, x = "AUC with bootstrap 95% CI", y = NULL) +
    theme_s13
}

pA <- plot_auc_forest(auc_df %>% filter(Cohort == "GSE78220"), "GSE78220: raw and composition-adjusted AUC", "Tumor-dedifferentiation/stromal-remodeling score")
pB <- plot_auc_forest(auc_df %>% filter(Cohort == "GSE91061"), "GSE91061: non-replication and adjustment sensitivity", "Tumor-dedifferentiation/stromal-remodeling score")

pC <- ggplot(loso_df %>% filter(!is.na(Method)), aes(y = Method)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.42, color = "black") +
  geom_segment(aes(x = LOSO_min, xend = LOSO_max, yend = Method), linewidth = 0.62, color = "black") +
  geom_point(aes(x = LOSO_median), size = 2.3, shape = 16, color = "black") +
  geom_point(aes(x = Full_AUC), size = 2.4, shape = 21, stroke = 0.7, fill = "white", color = "black") +
  facet_grid(. ~ Cohort) +
  scale_x_continuous(limits = c(0.16, 1.02), breaks = seq(0.2, 1.0, 0.2), labels = number_format(accuracy = 0.01)) +
  labs(title = "Leave-one-sample-out AUC stability", subtitle = "Line: LOSO min–max; filled point: LOSO median; open circle: full-sample AUC", x = "AUC", y = NULL) +
  theme_s13 +
  theme(
    strip.text = element_text(size = 9.8 * 1.15 * 1.25, face = "bold", color = "black")
  )

supp_s13 <- (pA | pB) / pC +
  plot_layout(heights = c(0.88, 1.12)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 21, face = "bold", family = "sans"),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

base_filename <- "Supplementary Figure S13. Formal composition-adjusted response analysis and leave-one-sample-out stability assessment in bulk ICB cohorts"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

safe_ggsave(out_png, supp_s13, width = 14.2, height = 10.6, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s13, width = 14.2, height = 10.6, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s13, width = 14.2, height = 10.6, dpi = dpi_out, device = cairo_pdf)

cat("\nSupplementary Figure S13 generated successfully with requested font scaling!\n")