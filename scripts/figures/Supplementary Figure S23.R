################--------------------------------------------
## Supplementary Figure S23. IF-mIF protein-state correlation analysis
##
## Purpose:
##   1) Load IF/mIF marker-signal–state Spearman correlation table from IF_mIF_validation
##   2) Generate and export Supplementary Figure S23 rectangular heatmap correlation analysis
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S23. IF-mIF protein-state correlation analysis.png
##     - Supplementary Figure S23. IF-mIF protein-state correlation analysis.jpg
##     - Supplementary Figure S23. IF-mIF protein-state correlation analysis.pdf
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
## 0. Paths and manual overrides
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
base_filename <- "Supplementary Figure S23. IF-mIF protein-state correlation analysis"

manual_corr_file <- file.path(
  project_dir,
  "results", "tables",
  "IF_mIF_validation",
  "IF_mIF_protein_state_spearman.csv"
)
manual_cell_file <- ""

rho_limits <- c(-0.35, 0.35)

################--------------------------------------------
## 1. Helper functions
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

find_col <- function(df, patterns, label, required = TRUE, exclude_patterns = character(0)) {
  cols <- colnames(df)
  
  for (pat in patterns) {
    hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE, perl = TRUE)
    
    if (length(hit) > 0 && length(exclude_patterns) > 0) {
      for (ep in exclude_patterns) {
        hit <- hit[!grepl(ep, hit, ignore.case = TRUE, perl = TRUE)]
      }
    }
    
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
  x1 <- gsub("/", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("immune.*defective|cold", xl) ~ "immune-defective/cold",
    grepl("myeloid.*treg|treg.*immunosuppress", xl) ~ "myeloid–Treg immunosuppressive",
    grepl("dediff|stromal.*remodel|tumor.*stromal|stromal remodeling", xl) ~ "tumor-dedifferentiation/stromal-remodeling",
    grepl("melanocytic|melanoma.*differentiation|differentiation", xl) ~ "melanocytic differentiation",
    TRUE ~ x0
  )
}

standardize_marker <- function(x) {
  x0 <- as.character(x)
  x1 <- x0
  x1 <- gsub("_", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("s100b.*pmel17|pmel17.*s100b|s100b/pmel17", xl) & grepl("mean", xl) ~ "S100B/PMEL17 (mean)",
    grepl("s100b.*pmel17|pmel17.*s100b|s100b/pmel17", xl) & grepl("max", xl) ~ "S100B/PMEL17 (max)",
    grepl("s100b.*pmel17|pmel17.*s100b|s100b/pmel17", xl) ~ "S100B/PMEL17",
    grepl("\\bcd3\\b", xl) & grepl("mean", xl) ~ "CD3 (mean)",
    grepl("\\bcd3\\b", xl) & grepl("max", xl) ~ "CD3 (max)",
    grepl("\\bcd3\\b", xl) ~ "CD3",
    grepl("\\bcd45\\b", xl) & grepl("mean", xl) ~ "CD45 (mean)",
    grepl("\\bcd45\\b", xl) & grepl("max", xl) ~ "CD45 (max)",
    grepl("\\bcd45\\b", xl) ~ "CD45",
    grepl("\\bdapi\\b", xl) & grepl("mean", xl) ~ "DAPI (mean)",
    grepl("\\bdapi\\b", xl) & grepl("max", xl) ~ "DAPI (max)",
    grepl("\\bdapi\\b", xl) ~ "DAPI",
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

theme_s23 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 16.2, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 10.2, hjust = 0, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10.0, color = "black", angle = 38, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10.8, color = "black"),
    legend.title = element_text(size = 11.0, face = "bold", color = "black"),
    legend.text = element_text(size = 9.5, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.45),
    plot.margin = margin(8, 16, 8, 8)
  )

################--------------------------------------------
## 2. Load and harmonize correlation data
################--------------------------------------------

input_mode <- "precomputed_correlation_table"
corr_file <- manual_corr_file
if (!nzchar(corr_file) || !file.exists(corr_file)) {
  stop("Correlation file does not exist: ", manual_corr_file)
}
corr_file <- normalizePath(corr_file, winslash = "/", mustWork = TRUE)
corr_raw <- read_any_table(corr_file)

marker_col <- find_col(
  corr_raw,
  patterns = c("^marker$", "^Marker$", "protein", "IF_marker", "mIF_marker", "channel", "feature", "signal"),
  label = "marker/signal column"
)

state_col <- find_col(
  corr_raw,
  patterns = c("^state$", "^State$", "tumor.*immune.*state", "signature", "program", "state_score"),
  label = "state column",
  required = FALSE
)

rho_col <- find_col(
  corr_raw,
  patterns = c("^rho$", "spearman.*rho", "spearman", "correlation", "^corr$", "^r$"),
  label = "Spearman rho column",
  required = FALSE
)

if (!is.na(state_col) && !is.na(rho_col)) {
  corr_df <- corr_raw %>%
    transmute(
      marker = standardize_marker(.data[[marker_col]]),
      state = standardize_state(.data[[state_col]]),
      rho = as.numeric(.data[[rho_col]])
    )
} else {
  wide_state_cols <- grep(
    "immune|myeloid|treg|dediff|stromal|melanocytic|differentiation|cold",
    colnames(corr_raw),
    value = TRUE,
    ignore.case = TRUE,
    perl = TRUE
  )
  
  corr_df <- corr_raw %>%
    select(all_of(c(marker_col, wide_state_cols))) %>%
    tidyr::pivot_longer(
      cols = all_of(wide_state_cols),
      names_to = "state",
      values_to = "rho"
    ) %>%
    transmute(
      marker = standardize_marker(.data[[marker_col]]),
      state = standardize_state(state),
      rho = as.numeric(rho)
    )
}

n_cells <- NA_integer_
n_col <- find_col(
  corr_raw,
  patterns = c("^n$", "n_cells", "n_cell", "n_matched", "n_used", "sample_size"),
  label = "n cells",
  required = FALSE
)
if (!is.na(n_col)) {
  n_vals <- suppressWarnings(as.numeric(corr_raw[[n_col]]))
  n_cells <- suppressWarnings(max(n_vals, na.rm = TRUE))
  if (!is.finite(n_cells)) n_cells <- NA_integer_
}

corr_df <- corr_df %>%
  mutate(
    marker = standardize_marker(marker),
    state = standardize_state(state),
    rho = as.numeric(rho)
  ) %>%
  filter(
    marker %in% marker_order,
    state %in% state_order,
    is.finite(rho)
  ) %>%
  group_by(marker, state) %>%
  summarise(rho = median(rho, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    marker = factor(marker, levels = rev(marker_order)),
    state = factor(state, levels = state_order)
  )

if (nrow(corr_df) == 0) {
  stop("S23 harmonized correlation table has zero rows after filtering.")
}

if (!is.finite(n_cells) || is.na(n_cells)) {
  n_cells_label <- "n not available"
} else {
  n_cells_label <- paste0("n = ", scales::comma(round(n_cells)))
}

################--------------------------------------------
## 3. Plot
################--------------------------------------------

subtitle_s23 <- paste0(
  "CosMx cell-level Spearman correlations after FOV–cell matching; ",
  n_cells_label,
  "; values are rho."
)

p_s23 <- ggplot(corr_df, aes(x = state, y = marker, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.75) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 4.1, color = "black") +
  scale_x_discrete(labels = state_labels, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient2(
    low = "#3B7DDD",
    mid = "white",
    high = "#EF3B2C",
    midpoint = 0,
    limits = rho_limits,
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

################--------------------------------------------
## 4. Save outputs
################--------------------------------------------

out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 10.6
fig_height <- 5.6

safe_ggsave(out_png, p_s23, width = fig_width, height = fig_height, dpi = dpi_out, device = NULL)
safe_ggsave(out_jpg, p_s23, width = fig_width, height = fig_height, dpi = dpi_out, device = NULL)
safe_ggsave(out_pdf, p_s23, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 5. Audit outputs
################--------------------------------------------

audit <- data.frame(
  input_mode = input_mode,
  corr_file = corr_file,
  n_cells_label = n_cells_label,
  n_harmonized_rows = nrow(corr_df),
  rho_min = min(corr_df$rho, na.rm = TRUE),
  rho_max = max(corr_df$rho, na.rm = TRUE),
  rho_limits = paste(rho_limits, collapse = ", "),
  title = "IF/mIF marker-signal–state correlations",
  subtitle = subtitle_s23,
  rectangular_heatmap_cells = TRUE,
  coord_fixed_ratio = 0.34,
  state_order = paste(state_order, collapse = "; "),
  marker_order = paste(marker_order, collapse = "; "),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S23_IF_mIF_protein_state_correlation_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

write.csv(
  corr_df %>%
    mutate(marker = as.character(marker), state = as.character(state)),
  file.path(fig_dir, "Supplementary_Figure_S23_IF_mIF_marker_signal_state_correlation_harmonized.csv"),
  row.names = FALSE
)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S23 script finished successfully.")