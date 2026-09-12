################--------------------------------------------
## Supplementary Figure S15. Visium spatial organization and state distribution
##
## Purpose:
##   1) Load spatial dominant state map, Moran's I raw vs QC-residualized, and neighborhood enrichment tables
##   2) Generate and export Supplementary Figure S15 (Panel A spatial dominant state map, Panel B Moran's I, Panel C neighborhood odds ratio)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S15. Visium spatial organization and state distribution.png
##     - Supplementary Figure S15. Visium spatial organization and state distribution.jpg
##     - Supplementary Figure S15. Visium spatial organization and state distribution.pdf
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

## 硬锁定空间验证表格路径
manual_panelA_file <- file.path(
  project_dir,
  "results", "tables", "spatial_melanoma_validation",
  "Step11_spatial_analysis_metadata_with_coordinates_and_scores.csv"
)

manual_panelB_file <- file.path(
  project_dir,
  "results", "tables", "spatial_melanoma_validation",
  "Step11_spatial_Morans_I_raw_and_QCresidual.csv"
)

manual_panelC_file <- file.path(
  project_dir,
  "results", "tables", "spatial_melanoma_validation",
  "Step11_spatial_high_state_neighbor_enrichment_top25_pooled.csv"
)

state_order <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation"
)

state_colors <- c(
  "immune-defective/cold" = "#4DBBD5",
  "myeloid–Treg immunosuppressive" = "#00A087",
  "tumor-dedifferentiation/stromal-remodeling" = "#E64B35",
  "melanocytic differentiation" = "#3C5488"
)

score_type_order <- c("QC-residualized", "Raw")
score_type_colors <- c("QC-residualized" = "#4DBBD5", "Raw" = "#E64B35")

################--------------------------------------------
## 1. Helper functions (支持多状态评分自动计算主导状态)
################--------------------------------------------

read_any_table <- function(file) {
  if (!file.exists(file)) stop("File does not exist: ", file)
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("csv", "txt")) {
    as.data.frame(suppressMessages(readr::read_csv(file, show_col_types = FALSE, guess_max = 100000)))
  } else if (ext == "tsv") {
    as.data.frame(suppressMessages(readr::read_tsv(file, show_col_types = FALSE, guess_max = 100000)))
  } else if (ext == "rds") {
    x <- readRDS(file)
    if (!is.data.frame(x)) stop("RDS is not a data.frame: ", file)
    as.data.frame(x)
  } else {
    stop("Unsupported file extension: ", file)
  }
}

find_col <- function(df, patterns, label, required = TRUE) {
  for (pat in patterns) {
    hit <- grep(pat, colnames(df), value = TRUE, ignore.case = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  if (required) {
    stop("Could not find column for ", label, ". Available columns: ", paste(colnames(df), collapse = ", "))
  }
  NA_character_
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

standardize_score_type <- function(x) {
  xl <- tolower(as.character(x))
  dplyr::case_when(
    grepl("resid|qc", xl) & !grepl("^raw", xl) ~ "QC-residualized",
    grepl("raw|unadjust|original", xl) ~ "Raw",
    TRUE ~ as.character(x)
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

harmonize_panelA_candidate <- function(dat) {
  message("Scanning Panel A columns: ", paste(colnames(dat), collapse = ", "))
  
  # 1. 检查是否存在显式的 DominantState 或状态分类列
  state_col <- find_col(dat, c(
    "dominant.*state", "DominantState", "dominant_state", "Raw_DominantState", "Dominant_state",
    "dominant.*category", "dominant.*label", "^state$", "^State$", "StateLabel", "category"
  ), "Panel A state", required = FALSE)
  
  x_col <- find_col(dat, c("^spatial_x$", "^x$", "spatial.*x", "pixel.*x", "full.*x", "coord.*x", "array_col", "^col$"), "Panel A x", required = FALSE)
  y_col <- find_col(dat, c("^spatial_y$", "^y$", "spatial.*y", "pixel.*y", "full.*y", "coord.*y", "array_row", "^row$"), "Panel A y", required = FALSE)
  
  if (is.na(x_col)) x_col <- find_col(dat, c("pxl_col", "imagecol"), "Panel A x fallback")
  if (is.na(y_col)) y_col <- find_col(dat, c("pxl_row", "imagerow"), "Panel A y fallback")
  
  # 如果没有显式分类列，但包含四个经典评分列，则通过计算最大值动态推导主导状态
  candidate_state_cols <- grep("Immune_defective_Cold|Myeloid_Treg_Immunosuppressive|Tumor_dedifferentiation_Stromal_remodeling|Melanocytic_Differentiation", colnames(dat), value = TRUE, ignore.case = TRUE)
  
  if (length(candidate_state_cols) >= 4) {
    message("Found 4 state score columns. Dynamically computing dominant state per spot.")
    score_mat <- as.matrix(dat[, candidate_state_cols])
    storage.mode(score_mat) <- "numeric"
    max_idx <- max.col(score_mat, ties.method = "first")
    dat$Computed_DominantState <- standardize_state(colnames(score_mat)[max_idx])
    state_col <- "Computed_DominantState"
  }
  
  if (is.na(state_col) || is.na(x_col) || is.na(y_col)) return(NULL)
  
  message(sprintf("Matched Panel A columns -> State: %s, X: %s, Y: %s", state_col, x_col, y_col))
  
  out <- dat %>%
    transmute(
      x = suppressWarnings(as.numeric(.data[[x_col]])),
      y = suppressWarnings(as.numeric(.data[[y_col]])),
      State = standardize_state(.data[[state_col]])
    ) %>%
    filter(is.finite(x), is.finite(y), State %in% state_order)
  
  if (nrow(out) == 0) return(NULL)
  attr(out, "state_col") <- state_col
  attr(out, "x_col") <- x_col
  attr(out, "y_col") <- y_col
  out
}

################--------------------------------------------
## 2. Read tables
################--------------------------------------------

panelA_file <- manual_panelA_file
panelB_file <- manual_panelB_file
panelC_file <- manual_panelC_file

panelA_raw <- read_any_table(panelA_file)
panelB_raw <- read_any_table(panelB_file)
panelC_raw <- read_any_table(panelC_file)

message("Panel A table: ", panelA_file)
message("Panel B table: ", panelB_file)
message("Panel C table: ", panelC_file)

################--------------------------------------------
## 3. Harmonize data panels
################--------------------------------------------

panelA_df <- harmonize_panelA_candidate(panelA_raw)
if (is.null(panelA_df) || nrow(panelA_df) == 0) stop("Panel A harmonization failed. Please check table structure.")
state_col_A <- attr(panelA_df, "state_col")
x_col_A <- attr(panelA_df, "x_col")
y_col_A <- attr(panelA_df, "y_col")
panelA_df <- panelA_df %>% mutate(State = factor(State, levels = state_order))

state_col_B <- find_col(panelB_raw, c("^state$", "^State$", "StateLabel", "state_label", "final.*state", "signature", "score"), "Panel B state")
type_col_B <- find_col(panelB_raw, c("score.*type", "score_type", "type", "analysis", "mode", "raw.*resid", "residual"), "Panel B score type", required = FALSE)
moran_col_B <- find_col(panelB_raw, c("Moran", "Morans", "Moran_I", "morans_i", "moran_i", "moran", "^I$", "spatial.*auto"), "Panel B Moran's I", required = FALSE)

if (!is.na(type_col_B) && !is.na(moran_col_B)) {
  panelB_df <- panelB_raw %>%
    transmute(State = standardize_state(.data[[state_col_B]]),
              ScoreType = standardize_score_type(.data[[type_col_B]]),
              MoransI = as.numeric(.data[[moran_col_B]]))
} else {
  raw_cols <- grep("raw.*moran|moran.*raw|raw.*I$|^raw$", colnames(panelB_raw), value = TRUE, ignore.case = TRUE)
  resid_cols <- grep("resid.*moran|moran.*resid|qc.*moran|moran.*qc|residual", colnames(panelB_raw), value = TRUE, ignore.case = TRUE)
  if (length(raw_cols) == 0 || length(resid_cols) == 0) stop("Could not detect Panel B Moran's I columns.")
  panelB_df <- panelB_raw %>%
    transmute(State = standardize_state(.data[[state_col_B]]),
              Raw = as.numeric(.data[[raw_cols[1]]]),
              `QC-residualized` = as.numeric(.data[[resid_cols[1]]])) %>%
    pivot_longer(cols = c("Raw", "QC-residualized"), names_to = "ScoreType", values_to = "MoransI")
}

panelB_df <- panelB_df %>%
  filter(State %in% state_order, ScoreType %in% score_type_order, is.finite(MoransI)) %>%
  group_by(State, ScoreType) %>%
  summarise(MoransI = MoransI[1], .groups = "drop") %>%
  mutate(State = factor(State, levels = state_order),
         ScoreType = factor(ScoreType, levels = score_type_order))

state_col_C <- find_col(panelC_raw, c("^state$", "^State$", "StateLabel", "state_label", "final.*state", "signature", "score"), "Panel C state")
or_col_C <- find_col(panelC_raw, c("^OR$", "^or$", "odds.*ratio", "odds_ratio", "edge.*odds", "ratio"), "Panel C odds ratio")
low_col_C <- find_col(panelC_raw, c("ci.*low", "lower", "or_low", "OR_low", "conf.low", "lcl", "lo95", "q025", "p025"), "Panel C lower CI", required = FALSE)
high_col_C <- find_col(panelC_raw, c("ci.*high", "upper", "or_high", "OR_high", "conf.high", "ucl", "hi95", "q975", "p975"), "Panel C upper CI", required = FALSE)

panelC_df <- panelC_raw %>%
  transmute(State = standardize_state(.data[[state_col_C]]),
            OR = as.numeric(.data[[or_col_C]]),
            CI_low = if (!is.na(low_col_C)) as.numeric(.data[[low_col_C]]) else NA_real_,
            CI_high = if (!is.na(high_col_C)) as.numeric(.data[[high_col_C]]) else NA_real_) %>%
  filter(State %in% state_order, is.finite(OR)) %>%
  group_by(State) %>%
  summarise(OR = OR[1], CI_low = CI_low[1], CI_high = CI_high[1], .groups = "drop") %>%
  mutate(CI_low = ifelse(is.finite(CI_low), CI_low, OR),
         CI_high = ifelse(is.finite(CI_high), CI_high, OR),
         State = factor(State, levels = state_order))

################--------------------------------------------
## 4. Plotting panels & Assembly
################--------------------------------------------

theme_s15 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 13.2, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 9.2, hjust = 0, color = "black"),
    axis.title = element_text(size = 11.0, face = "bold", color = "black"),
    axis.text = element_text(size = 9.6, color = "black"),
    legend.title = element_text(size = 10.0, face = "bold", color = "black"),
    legend.text = element_text(size = 8.8, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

pA <- ggplot(panelA_df, aes(x = x, y = -y, color = State)) +
  geom_point(size = 0.72, alpha = 0.92) +
  coord_equal() +
  scale_color_manual(values = state_colors, breaks = state_order, labels = wrap_state, name = "Dominant state", drop = FALSE) +
  labs(title = "Spatial distribution of dominant tumor–immune states",
       subtitle = "Coordinates from Visium tissue positions; y-axis is inverted for display",
       x = "Full-resolution pixel x", y = "Full-resolution pixel y") +
  theme_s15 +
  theme(legend.position = "right", legend.text = element_text(size = 8.2), legend.key.height = unit(0.45, "cm"))

pB <- ggplot(panelB_df, aes(x = State, y = MoransI, fill = ScoreType)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  scale_fill_manual(values = score_type_colors, breaks = score_type_order, name = "Score type") +
  scale_x_discrete(labels = wrap_state) +
  scale_y_continuous(labels = number_format(accuracy = 0.001), expand = expansion(mult = c(0.10, 0.12))) +
  labs(title = "Moran's I for raw and QC-residualized state scores",
       subtitle = "Global spatial autocorrelation; values close to zero indicate weak global autocorrelation",
       x = NULL, y = "Moran's I") +
  theme_s15 +
  theme(axis.text.x = element_text(size = 8.5, angle = 35, hjust = 1, vjust = 1),
        legend.position = "top", legend.justification = "right",
        legend.title = element_text(size = 8.5, face = "bold"), legend.text = element_text(size = 8.2))

pC <- ggplot(panelC_df, aes(x = State, y = OR, fill = State)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_col(width = 0.62) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.20, linewidth = 0.50, color = "black") +
  scale_fill_manual(values = state_colors, guide = "none") +
  scale_x_discrete(labels = wrap_state) +
  scale_y_continuous(limits = c(0, max(1.15, max(panelC_df$CI_high, na.rm = TRUE) * 1.10)),
                     breaks = pretty_breaks(n = 5), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "High-state neighborhood assessment",
       subtitle = "Dashed line indicates odds ratio = 1",
       x = NULL, y = "Edge-level odds ratio") +
  theme_s15 +
  theme(axis.text.x = element_text(size = 8.5, angle = 35, hjust = 1, vjust = 1))

right_col <- pB / pC + plot_layout(heights = c(0.92, 1.08))
supp_s15 <- pA | right_col
supp_s15 <- supp_s15 +
  plot_layout(widths = c(1.03, 1.55)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(size = 17, face = "bold", family = "sans"),
                                plot.margin = margin(4, 6, 4, 6)))

out_png <- file.path(fig_dir, "Supplementary Figure S15. Visium spatial organization and state distribution.png")
out_jpg <- file.path(fig_dir, "Supplementary Figure S15. Visium spatial organization and state distribution.jpg")
out_pdf <- file.path(fig_dir, "Supplementary Figure S15. Visium spatial organization and state distribution.pdf")

ggsave(out_png, supp_s15, width = 15.5, height = 7.8, units = "in", dpi = dpi_out, bg = "white", limitsize = FALSE)
ggsave(out_jpg, supp_s15, width = 15.5, height = 7.8, units = "in", dpi = dpi_out, bg = "white", limitsize = FALSE)
ggsave(out_pdf, supp_s15, width = 15.5, height = 7.8, units = "in", dpi = dpi_out, bg = "white", limitsize = FALSE, device = cairo_pdf)

cat("\nSupplementary Figure S15 generated successfully with dynamic dominant state derivation!\n")