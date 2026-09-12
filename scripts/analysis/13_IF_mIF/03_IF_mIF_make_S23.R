############################################################
## 03_IF_mIF_make_S23.R
##
## Public/frozen Supplementary Figure S23.
##
## Fixed source only; NO fuzzy auto-search.
## Scientific data/statistics are unchanged.
## Figure typography follows the frozen project-wide S20 reference.
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(scales)
})

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(PROJECT_DIR)) PROJECT_DIR <- "D:/ICB_resistance_project"
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)

TABLE_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "IF_mIF_validation"
)

FIG_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "figures",
  "IF_mIF_validation"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

CORR_FILE <- file.path(
  TABLE_DIR,
  "IF_mIF_protein_state_spearman.csv"
)

if (!file.exists(CORR_FILE)) {
  stop("Canonical IF/mIF correlation table not found. Run 02_IF_mIF_match_and_correlation.R first.", call. = FALSE)
}

corr_raw <- suppressMessages(
  readr::read_csv(
    CORR_FILE,
    show_col_types = FALSE
  )
)
corr_raw <- as.data.frame(corr_raw, stringsAsFactors = FALSE)

standardize_state <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("/", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)

  dplyr::case_when(
    grepl("immune.*defective|cold", xl) ~ "immune-defective/cold",
    grepl("myeloid.*treg|treg.*immunosuppress", xl) ~ "myeloid–Treg immunosuppressive",
    grepl("dediff|stromal.*remodel|tumor.*stromal|stromal remodeling", xl) ~ "tumor-dedifferentiation/stromal-remodeling",
    grepl("melanocytic.*differentiation|melanoma.*differentiation|^melanocytic$", xl) ~ "melanocytic differentiation",
    TRUE ~ x0
  )
}

standardize_marker <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)

  dplyr::case_when(
    grepl("s100b.*pmel17|pmel17.*s100b|s100b/pmel17", xl) & grepl("mean", xl) ~ "S100B/PMEL17 (mean)",
    grepl("s100b.*pmel17|pmel17.*s100b|s100b/pmel17", xl) & grepl("max", xl) ~ "S100B/PMEL17 (max)",
    grepl("\\bcd3\\b", xl) & grepl("mean", xl) ~ "CD3 (mean)",
    grepl("\\bcd3\\b", xl) & grepl("max", xl) ~ "CD3 (max)",
    grepl("\\bcd45\\b", xl) & grepl("mean", xl) ~ "CD45 (mean)",
    grepl("\\bcd45\\b", xl) & grepl("max", xl) ~ "CD45 (max)",
    grepl("\\bdapi\\b", xl) & grepl("mean", xl) ~ "DAPI (mean)",
    grepl("\\bdapi\\b", xl) & grepl("max", xl) ~ "DAPI (max)",
    TRUE ~ x0
  )
}

state_order <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation"
)

state_labels <- c(
  "immune-defective/cold" = "immune-defective/\ncold",
  "myeloid–Treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "melanocytic differentiation" = "melanocytic\ndifferentiation"
)

marker_order <- c(
  "CD3 (max)",
  "CD45 (max)",
  "DAPI (max)",
  "S100B/PMEL17 (max)",
  "CD3 (mean)",
  "CD45 (mean)",
  "DAPI (mean)",
  "S100B/PMEL17 (mean)"
)

corr_df <- corr_raw %>%
  dplyr::transmute(
    marker = standardize_marker(.data$protein_feature),
    state = standardize_state(.data$state_score),
    rho = as.numeric(.data$rho),
    n = as.numeric(.data$n)
  ) %>%
  dplyr::filter(
    marker %in% marker_order,
    state %in% state_order,
    is.finite(rho)
  ) %>%
  dplyr::group_by(marker, state) %>%
  dplyr::summarise(
    rho = median(rho, na.rm = TRUE),
    n = max(n, na.rm = TRUE),
    .groups = "drop"
  )

if (
  nrow(corr_df) != 32L ||
  length(unique(corr_df$marker)) != 8L ||
  length(unique(corr_df$state)) != 4L
) {
  stop(
    "S23 display contract failed: expected 32 rows, 8 markers, 4 states.",
    call. = FALSE
  )
}

n_cells <- unique(corr_df$n)
if (
  length(n_cells) != 1L ||
  !is.finite(n_cells) ||
  n_cells != 86572
) {
  stop("S23 n-cell contract failed.", call. = FALSE)
}

corr_df <- corr_df %>%
  dplyr::mutate(
    marker = factor(marker, levels = rev(marker_order)),
    state = factor(state, levels = state_order)
  )

subtitle_s23 <- paste0(
  "CosMx cell-level Spearman correlations after FOV–cell matching; ",
  "n = ", scales::comma(n_cells),
  "; values are rho."
)

theme_s23 <- theme_bw(
  base_size = 12.5,
  base_family = "sans"
) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0,
      color = "black"
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0,
      color = "black"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(
      size = 10.5,
      color = "black",
      angle = 38,
      hjust = 1,
      vjust = 1,
      lineheight = 0.92
    ),
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    legend.title = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    legend.text = element_text(
      size = 10.5,
      color = "black"
    ),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "grey35",
      fill = NA,
      linewidth = 0.45
    ),
    plot.margin = margin(8, 16, 8, 8)
  )

p <- ggplot(corr_df, aes(x = state, y = marker, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.75) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 4.1, color = "black") +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient2(
    low = "#3B7DDD",
    mid = "white",
    high = "#EF3B2C",
    midpoint = 0,
    limits = c(-0.35, 0.35),
    oob = scales::squish,
    breaks = c(-0.3, 0, 0.3),
    labels = c("-0.3", "0", "0.3"),
    name = "Spearman\nrho"
  ) +
  labs(
    title = "IF/mIF marker-signal–state correlations",
    subtitle = subtitle_s23
  ) +
  coord_fixed(ratio = 0.34) +
  theme_s23

prefix <-
  "Supplementary_Figure_S23_IF_mIF_marker_signal_state_correlation"

out_png <- file.path(FIG_DIR, paste0(prefix, ".png"))
out_jpg <- file.path(FIG_DIR, paste0(prefix, ".jpg"))
out_pdf <- file.path(FIG_DIR, paste0(prefix, ".pdf"))

ggplot2::ggsave(
  out_png, p,
  width = 10.6, height = 5.6,
  units = "in", dpi = 300,
  bg = "white", limitsize = FALSE
)
ggplot2::ggsave(
  out_jpg, p,
  width = 10.6, height = 5.6,
  units = "in", dpi = 300,
  bg = "white", limitsize = FALSE
)
ggplot2::ggsave(
  out_pdf, p,
  width = 10.6, height = 5.6,
  units = "in", dpi = 300,
  bg = "white", limitsize = FALSE
)

harmonized_out <- corr_df %>%
  dplyr::mutate(
    marker = as.character(marker),
    state = as.character(state)
  ) %>%
  dplyr::arrange(
    match(marker, marker_order),
    match(state, state_order)
  )

harmonized_file <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S23_IF_mIF_marker_signal_state_correlation_source.csv"
)
utils::write.csv(
  harmonized_out,
  harmonized_file,
  row.names = FALSE
)

audit <- data.frame(
  input_mode = "fixed_public_sequential_correlation_table",
  corr_file = normalizePath(CORR_FILE, winslash = "/", mustWork = TRUE),
  n_cells = as.integer(n_cells),
  n_harmonized_rows = nrow(harmonized_out),
  n_markers = length(unique(harmonized_out$marker)),
  n_states = length(unique(harmonized_out$state)),
  rho_min = min(harmonized_out$rho),
  rho_max = max(harmonized_out$rho),
  rho_limits = "-0.35, 0.35",
  coord_fixed_ratio = 0.34,
  fuzzy_autosearch = FALSE,
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

utils::write.csv(
  audit,
  file.path(FIG_DIR, paste0(prefix, "_audit.csv")),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("IF/mIF SUPPLEMENTARY FIGURE S23 COMPLETED\n")
cat("Rows: ", nrow(harmonized_out), " | markers: 8 | states: 4\n", sep = "")
cat("n cells: ", n_cells, "\n", sep = "")
cat("============================================================\n")
