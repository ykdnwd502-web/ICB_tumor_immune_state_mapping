############################################################
## 02_CosMx_cell_level_spatial_validation.R
## v1.1 intermediate plots -> results/diagnostics/11_CosMx.
##
## PUBLIC SEQUENTIAL ROLE
##   Step01 state scores + raw Run5611_MK3 metadata coordinates
##     -> cell-level spatial/state/dual-high table
##
## The frozen canonical audit identifies Run5611_MK3_metadata_file.csv,
## not the polygon file, as the actual Step15 coordinate source.
## Scientific/statistical logic is unchanged from Step15 v2 FIRST_FIX.
############################################################

# ============================================================
# Step 15. CosMx cell-level state validation + dual-high spatial map
# Clean-room reproducibility version v2
# ------------------------------------------------------------
# Purpose:
#   1) Use Step14B CosMx cell-level state scores as input.
#   2) Validate state-score distributions and dominant-state composition.
#   3) Define Myeloid–Treg high and Dediff/Stromal high cells.
#   4) Test cell-level dual-high co-occurrence enrichment.
#   5) Merge CosMx global coordinates and draw spatial maps:
#      - dominant state spatial map
#      - dual-high co-occurrence spatial map
#
# Important:
#   - This script uses fov + cell_ID as the cell-level composite key.
#   - It automatically screens candidate coordinate files and selects the
#     file with the highest overlap and usable x/y coordinates.
#   - If a polygon/vertex-level coordinate file is selected, coordinates
#     are aggregated to cell-level centroids.
#
# Run:
#   source("D:/ICB_resistance_project/scripts/15_CosMx_cell_level_spatial_validation_CLEAN_INPUT_v2_FIRST_FIX.R")
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ------------------------------
# 0. Config
# ------------------------------
project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir,winslash="/",mustWork=TRUE)

out_table_dir <- file.path(project_dir, "results", "tables", "CosMx_spatial_validation")
out_fig_dir <- file.path(
  project_dir,
  "results", "diagnostics", "11_CosMx",
  "02_spatial_validation"
)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Step14B v4 output. Do not recalculate state scores here.
manual_state_score_file <- file.path(
  project_dir,
  "results", "tables", "IF_mIF_validation",
  "CosMx_cell_state_scores.csv"
)

# Optional manual coordinate file.
# Leave NA to auto-detect. If auto-detection fails, set to the full path of
# Run5611_MK3-polygons.csv or another cell-coordinate table.
manual_coordinate_file <- file.path(
  project_dir,
  "data_raw", "CosMx_melanocytic_tumors_Dryad",
  "Slide_4", "Run5611_MK3",
  "Run5611_MK3_metadata_file.csv"
)

# Main dual-high cutoff.
main_top_fraction <- 0.25

# Sensitivity cutoffs.
sensitivity_top_fractions <- c(0.20, 0.25, 0.30)

# For publication display; pixel-like coordinates are usually shown with y reversed.
invert_y_for_display <- TRUE

# Downsample only the grey background if needed. Set Inf to keep all cells.
max_background_points <- Inf

state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_display <- c(
  Immune_defective_Cold = "Immune-defective/Cold",
  Myeloid_Treg_Immunosuppressive = "Myeloid–Treg\nImmunosuppressive",
  Tumor_dedifferentiation_Stromal_remodeling = "Tumor dediff./\nStromal remodeling",
  Melanocytic_Differentiation = "Melanocytic\nDifferentiation"
)

## Public-figure display lock; analytical/table labels above remain untouched.
state_display_public <- c(
  "Immune_defective_Cold" = "immune-defective/\ncold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic\ndifferentiation"
)

# Colors aligned with earlier Figure 8/dual-high style.
dual_high_colors <- c(
  "Neither high" = "#D9D9D9",
  "Myeloid–Treg high only" = "#00A087",
  "Dediff/Stromal high only" = "#E64B35",
  "Both high" = "#8E44AD"
)

dominant_state_colors <- c(
  "Immune_defective_Cold" = "#4DBBD5",
  "Myeloid_Treg_Immunosuppressive" = "#00A087",
  "Tumor_dedifferentiation_Stromal_remodeling" = "#E64B35",
  "Melanocytic_Differentiation" = "#3C5488"
)

# ------------------------------
# 1. Utilities
# ------------------------------
safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE)
  message("Saved table: ", path)
}

safe_num <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

sanitize_names <- function(df) {
  nm <- colnames(df)
  nm[is.na(nm) | nm == ""] <- paste0("blank_col_", which(is.na(nm) | nm == ""))
  colnames(df) <- make.unique(nm, sep = "_")
  df
}

read_csv_safely <- function(path, nrows = -1, col_classes = NA) {
  if (!file.exists(path)) return(NULL)
  if (requireNamespace("data.table", quietly = TRUE) && nrows < 0) {
    df <- tryCatch(
      as.data.frame(data.table::fread(path, data.table = FALSE, showProgress = FALSE)),
      error = function(e) NULL
    )
  } else {
    df <- tryCatch(
      read.csv(path, check.names = FALSE, nrows = nrows, colClasses = col_classes),
      error = function(e) NULL
    )
  }
  if (is.null(df)) return(NULL)
  sanitize_names(df)
}

clean_text <- function(x) {
  x <- as.character(x)
  trimws(x)
}

make_key_raw <- function(fov, cell_id) {
  paste0(clean_text(fov), "__", clean_text(cell_id))
}

norm_numeric_like <- function(x) {
  x <- clean_text(x)
  x <- gsub("\\.0$", "", x)
  has_digit <- grepl("[0-9]", x)
  out <- x
  if (any(has_digit)) {
    first_num <- regmatches(x[has_digit], regexpr("[0-9]+", x[has_digit]))
    first_num <- gsub("^0+", "", first_num)
    first_num[first_num == ""] <- "0"
    out[has_digit] <- first_num
  }
  tolower(out)
}

make_key_norm <- function(fov, cell_id) {
  paste0("fov", norm_numeric_like(fov), "__cell", norm_numeric_like(cell_id))
}

guess_col <- function(cols, patterns) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

candidate_xy_pairs <- function(cols) {
  cols_low <- tolower(cols)

  # Prefer global/centroid coordinate fields.
  x_patterns <- c(
    "global.*x", "x.*global",
    "center.*x", "centroid.*x",
    "x_center", "x_centroid",
    "^x$", "x_local", "local.*x",
    "pixel.*x", "x.*pixel",
    "cell.*x"
  )
  y_patterns <- c(
    "global.*y", "y.*global",
    "center.*y", "centroid.*y",
    "y_center", "y_centroid",
    "^y$", "y_local", "local.*y",
    "pixel.*y", "y.*pixel",
    "cell.*y"
  )

  x_hits <- unique(unlist(lapply(x_patterns, function(p) cols[grepl(p, cols_low)])))
  y_hits <- unique(unlist(lapply(y_patterns, function(p) cols[grepl(p, cols_low)])))

  if (length(x_hits) == 0 || length(y_hits) == 0) {
    return(data.frame(x_col = character(0), y_col = character(0), xy_priority = numeric(0)))
  }

  pairs <- expand.grid(x_col = x_hits, y_col = y_hits, stringsAsFactors = FALSE)

  # Score likely pairs. Higher is better.
  score_one <- function(x, y) {
    xl <- tolower(x); yl <- tolower(y)
    s <- 0

    if (grepl("global", xl) && grepl("global", yl)) s <- s + 100
    if (grepl("center|centroid", xl) && grepl("center|centroid", yl)) s <- s + 50
    if (grepl("pixel|px", xl) && grepl("pixel|px", yl)) s <- s + 20
    if (grepl("local", xl) && grepl("local", yl)) s <- s + 10
    if (grepl("^x$", xl) && grepl("^y$", yl)) s <- s + 5

    # Penalize mismatched names that are unlikely to be paired.
    x_base <- gsub("x", "", xl)
    y_base <- gsub("y", "", yl)
    if (nchar(x_base) > 0 && nchar(y_base) > 0 && x_base == y_base) s <- s + 5

    s
  }

  pairs$xy_priority <- mapply(score_one, pairs$x_col, pairs$y_col)
  pairs <- pairs[order(-pairs$xy_priority), ]
  pairs
}

check_xy_numeric <- function(df, x_col, y_col) {
  if (!all(c(x_col, y_col) %in% colnames(df))) return(FALSE)
  x <- safe_num(df[[x_col]])
  y <- safe_num(df[[y_col]])
  mean(is.finite(x) & is.finite(y), na.rm = TRUE) > 0.30
}

safe_fisher_or <- function(myeloid_high, dediff_high) {
  tab <- table(
    Myeloid_Treg_high = factor(myeloid_high, levels = c(FALSE, TRUE)),
    Dediff_Stromal_high = factor(dediff_high, levels = c(FALSE, TRUE))
  )

  ft <- fisher.test(tab)
  data.frame(
    neither = as.integer(tab["FALSE", "FALSE"]),
    myeloid_only = as.integer(tab["TRUE", "FALSE"]),
    dediff_only = as.integer(tab["FALSE", "TRUE"]),
    both = as.integer(tab["TRUE", "TRUE"]),
    odds_ratio = unname(ft$estimate),
    conf_low = unname(ft$conf.int[1]),
    conf_high = unname(ft$conf.int[2]),
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

cutoff_by_mode <- function(df, score_col, top_fraction = 0.25, mode = c("pooled", "per_fov")) {
  mode <- match.arg(mode)
  prob <- 1 - top_fraction
  if (mode == "pooled") {
    q <- stats::quantile(df[[score_col]], probs = prob, na.rm = TRUE, names = FALSE)
    return(df[[score_col]] >= q)
  }

  df %>%
    group_by(fov) %>%
    mutate(.tmp_cutoff = stats::quantile(.data[[score_col]], probs = prob, na.rm = TRUE, names = FALSE)) %>%
    ungroup() %>%
    transmute(.high = .data[[score_col]] >= .tmp_cutoff) %>%
    pull(.high)
}

# ------------------------------
# 2. Read CosMx state-score table
# ------------------------------
if (!file.exists(manual_state_score_file)) {
  stop("State-score file not found. Please run Step14B v4 first or set manual_state_score_file.")
}

state_df <- read_csv_safely(manual_state_score_file)
if (is.null(state_df)) stop("Could not read state-score file: ", manual_state_score_file)

missing_states <- setdiff(state_cols, colnames(state_df))
if (length(missing_states) > 0) {
  stop("Missing state-score columns: ", paste(missing_states, collapse = ", "))
}

fov_col <- guess_col(colnames(state_df), c("^fov$", "FOV", "field"))
cell_col <- guess_col(colnames(state_df), c("^cell_ID$", "cellid", "cell_id", "CellID", "Cell_ID", "^cell$"))
if (is.na(fov_col) || is.na(cell_col)) {
  stop("State-score table must contain fov and cell_ID columns.")
}

state_df <- state_df %>%
  mutate(
    fov = clean_text(.data[[fov_col]]),
    cell_ID = clean_text(.data[[cell_col]]),
    .key_raw = make_key_raw(fov, cell_ID),
    .key_norm = make_key_norm(fov, cell_ID)
  )

for (cc in state_cols) {
  state_df[[cc]] <- safe_num(state_df[[cc]])
}

# Remove duplicated keys if any, but preserve audit.
dup_n <- sum(duplicated(state_df$.key_raw))
state_keep <- state_df %>%
  select(fov, cell_ID, .key_raw, .key_norm, all_of(state_cols)) %>%
  group_by(.key_raw) %>%
  summarise(
    # Use dplyr::first() explicitly because other packages may mask first().
    fov = dplyr::first(fov),
    cell_ID = dplyr::first(cell_ID),
    .key_norm = dplyr::first(.key_norm),
    across(all_of(state_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

for (cc in state_cols) state_keep[[cc]][is.nan(state_keep[[cc]])] <- NA_real_

safe_write_csv(
  data.frame(
    StateScoreFile = manual_state_score_file,
    n_rows = nrow(state_df),
    n_unique_raw_key = length(unique(state_df$.key_raw)),
    n_duplicate_raw_key = dup_n,
    fov_col = fov_col,
    cell_col = cell_col,
    state_columns = paste(state_cols, collapse = ";"),
    stringsAsFactors = FALSE
  ),
  file.path(out_table_dir, "CosMx_state_score_input_audit.csv")
)

# ------------------------------
# 3. State score summary and dominant state
# ------------------------------
state_summary <- lapply(state_cols, function(st) {
  x <- state_keep[[st]]
  data.frame(
    state = st,
    display_state = state_display[[st]],
    n_finite = sum(is.finite(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    q75 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

safe_write_csv(state_summary, file.path(out_table_dir, "CosMx_state_score_summary.csv"))

state_mat <- as.matrix(state_keep[, state_cols])
dominant_idx <- max.col(state_mat, ties.method = "first")
state_keep$Dominant_state <- state_cols[dominant_idx]
state_keep$Dominant_state_label <- unname(state_display[state_keep$Dominant_state])

dominant_counts <- state_keep %>%
  count(Dominant_state, Dominant_state_label, name = "n_cells") %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  arrange(desc(n_cells))

safe_write_csv(dominant_counts, file.path(out_table_dir, "CosMx_dominant_state_counts.csv"))

# ------------------------------
# 4. Dual-high classification and OR enrichment
# ------------------------------
myeloid_col <- "Myeloid_Treg_Immunosuppressive"
dediff_col  <- "Tumor_dedifferentiation_Stromal_remodeling"

myeloid_cut <- quantile(state_keep[[myeloid_col]], probs = 1 - main_top_fraction, na.rm = TRUE, names = FALSE)
dediff_cut  <- quantile(state_keep[[dediff_col]],  probs = 1 - main_top_fraction, na.rm = TRUE, names = FALSE)

state_keep <- state_keep %>%
  mutate(
    Myeloid_Treg_high = .data[[myeloid_col]] >= myeloid_cut,
    Dediff_Stromal_high = .data[[dediff_col]] >= dediff_cut,
    dual_high_group = case_when(
      Myeloid_Treg_high & Dediff_Stromal_high ~ "Both high",
      Myeloid_Treg_high & !Dediff_Stromal_high ~ "Myeloid–Treg high only",
      !Myeloid_Treg_high & Dediff_Stromal_high ~ "Dediff/Stromal high only",
      TRUE ~ "Neither high"
    ),
    dual_high_group = factor(
      dual_high_group,
      levels = c("Neither high", "Myeloid–Treg high only", "Dediff/Stromal high only", "Both high")
    )
  )

cutoff_audit <- data.frame(
  cutoff_mode = "pooled",
  top_fraction = main_top_fraction,
  Myeloid_Treg_cutoff = myeloid_cut,
  Dediff_Stromal_cutoff = dediff_cut,
  n_cells = nrow(state_keep),
  stringsAsFactors = FALSE
)
safe_write_csv(cutoff_audit, file.path(out_table_dir, "CosMx_dual_high_cutoff_audit.csv"))

group_counts <- state_keep %>%
  count(dual_high_group, name = "n_cells") %>%
  mutate(prop = n_cells / sum(n_cells))

safe_write_csv(group_counts, file.path(out_table_dir, "CosMx_dual_high_group_counts.csv"))

or_main <- safe_fisher_or(state_keep$Myeloid_Treg_high, state_keep$Dediff_Stromal_high) %>%
  mutate(
    cutoff_mode = "pooled",
    top_fraction = main_top_fraction,
    analysis_unit = "CosMx cell",
    n_cells = nrow(state_keep),
    .before = 1
  )

safe_write_csv(or_main, file.path(out_table_dir, "CosMx_dual_high_OR_top25_pooled.csv"))

# Sensitivity: pooled vs per_fov.
sens_list <- list()
k <- 1
for (frac in sensitivity_top_fractions) {
  for (mode in c("pooled", "per_fov")) {
    mh <- cutoff_by_mode(state_keep, myeloid_col, top_fraction = frac, mode = mode)
    dh <- cutoff_by_mode(state_keep, dediff_col, top_fraction = frac, mode = mode)
    sens_list[[k]] <- safe_fisher_or(mh, dh) %>%
      mutate(
        cutoff_mode = mode,
        top_fraction = frac,
        analysis_unit = "CosMx cell",
        n_cells = nrow(state_keep),
        .before = 1
      )
    k <- k + 1
  }
}

sens_tbl <- bind_rows(sens_list)
safe_write_csv(sens_tbl, file.path(out_table_dir, "CosMx_dual_high_OR_threshold_sensitivity.csv"))

# Save state table with dual-high labels before coordinate merge.
safe_write_csv(
  state_keep,
  file.path(out_table_dir, "CosMx_cell_state_scores_with_dual_high_groups.csv")
)

# ------------------------------
# 5. Coordinate file screening and merge
# ------------------------------
if (!is.na(manual_coordinate_file) && file.exists(manual_coordinate_file)) {
  coord_candidates <- manual_coordinate_file
} else {
  coord_candidates <- list.files(
    file.path(project_dir, "data_raw", "CosMx_melanocytic_tumors_Dryad"),
    pattern = "polygon|coordinate|coords|cell.*meta|metadata_file|fov_positions",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  coord_candidates <- coord_candidates[grepl("\\.csv$", coord_candidates, ignore.case = TRUE)]

  # Keep Slide 4 / Run5611 candidates first, but audit all relevant candidates.
  coord_candidates <- unique(c(
    coord_candidates[grepl("Slide_4|Slide4|Run5611|MK3", coord_candidates, ignore.case = TRUE)],
    coord_candidates
  ))
}

if (length(coord_candidates) == 0) {
  stop("No coordinate candidate files found. Set manual_coordinate_file.")
}

coord_audit_list <- list()
coord_cache <- list()

for (f in unique(coord_candidates)) {
  # Read only header + small number first to detect columns.
  preview <- read_csv_safely(f, nrows = 2000, col_classes = "character")
  if (is.null(preview)) {
    coord_audit_list[[length(coord_audit_list) + 1]] <- data.frame(
      file = f, basename = basename(f), n_preview_rows = NA_integer_,
      fov_col = NA_character_, cell_col = NA_character_,
      x_col = NA_character_, y_col = NA_character_,
      xy_priority = NA_real_, raw_overlap = 0, norm_overlap = 0,
      selected_candidate = FALSE, stringsAsFactors = FALSE
    )
    next
  }

  fov_c <- guess_col(colnames(preview), c("^fov$", "FOV", "field"))
  cell_c <- guess_col(colnames(preview), c("^cell_ID$", "cellid", "cell_id", "CellID", "Cell_ID", "^cell$"))

  xy_pairs <- candidate_xy_pairs(colnames(preview))

  if (is.na(fov_c) || is.na(cell_c) || nrow(xy_pairs) == 0) {
    coord_audit_list[[length(coord_audit_list) + 1]] <- data.frame(
      file = f, basename = basename(f), n_preview_rows = nrow(preview),
      fov_col = fov_c, cell_col = cell_c,
      x_col = NA_character_, y_col = NA_character_,
      xy_priority = NA_real_, raw_overlap = 0, norm_overlap = 0,
      selected_candidate = FALSE, stringsAsFactors = FALSE
    )
    next
  }

  # Test candidate xy pairs on preview; keep best numeric pair.
  xy_pairs$numeric_ok <- mapply(function(xc, yc) check_xy_numeric(preview, xc, yc), xy_pairs$x_col, xy_pairs$y_col)
  xy_pairs <- xy_pairs[xy_pairs$numeric_ok, , drop = FALSE]
  if (nrow(xy_pairs) == 0) {
    coord_audit_list[[length(coord_audit_list) + 1]] <- data.frame(
      file = f, basename = basename(f), n_preview_rows = nrow(preview),
      fov_col = fov_c, cell_col = cell_c,
      x_col = NA_character_, y_col = NA_character_,
      xy_priority = NA_real_, raw_overlap = 0, norm_overlap = 0,
      selected_candidate = FALSE, stringsAsFactors = FALSE
    )
    next
  }

  x_c <- xy_pairs$x_col[1]
  y_c <- xy_pairs$y_col[1]

  # Read full file only when it has promising columns.
  coord_full <- read_csv_safely(f)
  if (is.null(coord_full)) next

  # Some coordinate files are huge polygon-vertex files. Keep only required columns.
  coord_full <- coord_full[, unique(c(fov_c, cell_c, x_c, y_c)), drop = FALSE]
  coord_full$.key_raw <- make_key_raw(coord_full[[fov_c]], coord_full[[cell_c]])
  coord_full$.key_norm <- make_key_norm(coord_full[[fov_c]], coord_full[[cell_c]])
  coord_full$.x <- safe_num(coord_full[[x_c]])
  coord_full$.y <- safe_num(coord_full[[y_c]])

  coord_cell <- coord_full %>%
    filter(is.finite(.x), is.finite(.y)) %>%
    group_by(.key_raw, .key_norm) %>%
    summarise(
      spatial_x = mean(.x, na.rm = TRUE),
      spatial_y = mean(.y, na.rm = TRUE),
      n_coordinate_rows = n(),
      .groups = "drop"
    )

  raw_overlap <- length(intersect(coord_cell$.key_raw, state_keep$.key_raw))
  norm_overlap <- length(intersect(coord_cell$.key_norm, state_keep$.key_norm))

  selected_candidate <- max(raw_overlap, norm_overlap) >= 1000

  coord_audit_list[[length(coord_audit_list) + 1]] <- data.frame(
    file = f,
    basename = basename(f),
    n_preview_rows = nrow(preview),
    n_coordinate_cells = nrow(coord_cell),
    fov_col = fov_c,
    cell_col = cell_c,
    x_col = x_c,
    y_col = y_c,
    xy_priority = xy_pairs$xy_priority[1],
    raw_overlap = raw_overlap,
    norm_overlap = norm_overlap,
    best_strategy = ifelse(raw_overlap >= norm_overlap, "raw", "norm"),
    best_overlap = max(raw_overlap, norm_overlap),
    selected_candidate = selected_candidate,
    stringsAsFactors = FALSE
  )

  coord_cache[[f]] <- list(
    coord_cell = coord_cell,
    fov_col = fov_c,
    cell_col = cell_c,
    x_col = x_c,
    y_col = y_c,
    raw_overlap = raw_overlap,
    norm_overlap = norm_overlap,
    best_strategy = ifelse(raw_overlap >= norm_overlap, "raw", "norm"),
    best_overlap = max(raw_overlap, norm_overlap)
  )
}

coord_audit <- bind_rows(coord_audit_list) %>%
  arrange(desc(best_overlap), desc(xy_priority))

safe_write_csv(coord_audit, file.path(out_table_dir, "CosMx_coordinate_candidate_overlap_audit.csv"))

if (nrow(coord_audit) == 0 || max(coord_audit$best_overlap, na.rm = TRUE) < 1000) {
  stop("No usable coordinate file found. Inspect CosMx_coordinate_candidate_overlap_audit.csv.")
}

selected_coord_file <- coord_audit$file[1]
coord_info <- coord_cache[[selected_coord_file]]
coord_cell <- coord_info$coord_cell
coord_strategy <- coord_info$best_strategy

message("Selected coordinate file: ", selected_coord_file)
message("Coordinate merge strategy: ", coord_strategy)
message("Coordinate best overlap: ", coord_info$best_overlap)

if (coord_strategy == "raw") {
  plot_df <- state_keep %>%
    inner_join(coord_cell %>% select(.key_raw, spatial_x, spatial_y, n_coordinate_rows), by = ".key_raw")
} else {
  plot_df <- state_keep %>%
    inner_join(coord_cell %>% select(.key_norm, spatial_x, spatial_y, n_coordinate_rows), by = ".key_norm")
}

coord_merge_audit <- data.frame(
  SelectedCoordinateFile = selected_coord_file,
  CoordinateStrategy = coord_strategy,
  CoordinateFovCol = coord_info$fov_col,
  CoordinateCellCol = coord_info$cell_col,
  CoordinateXCol = coord_info$x_col,
  CoordinateYCol = coord_info$y_col,
  StateCells = nrow(state_keep),
  CoordinateCells = nrow(coord_cell),
  MergedCells = nrow(plot_df),
  CoordinateRowsMedianPerCell = median(plot_df$n_coordinate_rows, na.rm = TRUE),
  CoordinateRowsMaxPerCell = max(plot_df$n_coordinate_rows, na.rm = TRUE),
  stringsAsFactors = FALSE
)

safe_write_csv(coord_merge_audit, file.path(out_table_dir, "CosMx_coordinate_merge_audit.csv"))
safe_write_csv(plot_df, file.path(out_table_dir, "CosMx_state_dual_high_with_coordinates.csv"))

# ------------------------------
# 6. Spatial maps
# ------------------------------
theme_spatial <- theme_bw(base_size = 12.5, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5, color = "black"),
    axis.title = element_text(face = "bold", size = 12.5, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(face = "bold", size = 12.5, color = "black"),
    legend.text = element_text(size = 11, color = "black"),
    panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.60),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(6, 6, 6, 6)
  )

# Dominant-state spatial map.
p_dom <- ggplot(plot_df, aes(x = spatial_x, y = spatial_y, color = Dominant_state)) +
  geom_point(size = 0.18, alpha = 0.85) +
  scale_color_manual(
    values = dominant_state_colors,
    breaks = state_cols,
    labels = state_display_public[state_cols],
    name = "Dominant state",
    drop = FALSE
  ) +
  coord_fixed() +
  labs(
    title = "CosMx remapping of RNA-defined tumor–immune states",
    subtitle = paste0("Unit: CosMx cell; coordinate source: ", basename(selected_coord_file)),
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_spatial

if (invert_y_for_display) p_dom <- p_dom + scale_y_reverse()

ggsave(
  file.path(out_fig_dir, "CosMx_dominant_state_spatial_map.png"),
  p_dom, width = 11.5, height = 9, dpi = 500, bg = "white"
)
ggsave(
  file.path(out_fig_dir, "CosMx_dominant_state_spatial_map.pdf"),
  p_dom, width = 11.5, height = 9, device = cairo_pdf, bg = "white"
)

# Dual-high map: grey background plus colored high-state groups.
set.seed(20260504)
bg_df <- plot_df %>% filter(dual_high_group == "Neither high")
if (is.finite(max_background_points) && nrow(bg_df) > max_background_points) {
  bg_df <- bg_df %>% slice_sample(n = max_background_points)
}
fg_df <- plot_df %>% filter(dual_high_group != "Neither high")

p_dual <- ggplot() +
  geom_point(
    data = bg_df,
    aes(x = spatial_x, y = spatial_y),
    color = dual_high_colors[["Neither high"]],
    size = 0.12,
    alpha = 0.35
  ) +
  geom_point(
    data = fg_df,
    aes(x = spatial_x, y = spatial_y, color = dual_high_group),
    size = 0.22,
    alpha = 0.90
  ) +
  scale_color_manual(
    values = dual_high_colors[names(dual_high_colors) != "Neither high"],
    labels = c(
      "Myeloid–Treg high only" = "myeloid–Treg high only",
      "Dediff/Stromal high only" = "tumor-dedifferentiation/\nstromal-remodeling high only",
      "Both high" = "both high"
    ),
    name = "Dual-high group",
    drop = FALSE
  ) +
  coord_fixed() +
  labs(
    title = "CosMx cell-level dual-high co-occurrence map",
    subtitle = paste0(
      "High-state cutoff: top ", main_top_fraction * 100,
      "% pooled; unit: CosMx cell; y-axis",
      ifelse(invert_y_for_display, " inverted for display", " not inverted")
    ),
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_spatial

if (invert_y_for_display) p_dual <- p_dual + scale_y_reverse()

ggsave(
  file.path(out_fig_dir, "CosMx_dual_high_spatial_map.png"),
  p_dual, width = 11.5, height = 9, dpi = 500, bg = "white"
)
ggsave(
  file.path(out_fig_dir, "CosMx_dual_high_spatial_map.pdf"),
  p_dual, width = 11.5, height = 9, device = cairo_pdf, bg = "white"
)

# OR sensitivity plot.
sens_plot_tbl <- sens_tbl %>%
  mutate(
    cutoff_mode = factor(cutoff_mode, levels = c("pooled", "per_fov")),
    top_fraction = as.numeric(top_fraction)
  )

p_sens <- ggplot(sens_plot_tbl, aes(x = top_fraction, y = odds_ratio, color = cutoff_mode, group = cutoff_mode)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.01, linewidth = 0.5) +
  geom_point(size = 2.6) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = sensitivity_top_fractions) +
  labs(
    title = "CosMx dual-high co-occurrence threshold sensitivity",
    subtitle = "myeloid–Treg high and tumor-dedifferentiation/stromal-remodeling high; unit of analysis: CosMx cell",
    x = "Top fraction used to define high-state",
    y = "Cell-level odds ratio",
    color = "Cutoff mode"
  ) +
  theme_spatial

ggsave(
  file.path(out_fig_dir, "CosMx_dual_high_OR_threshold_sensitivity.png"),
  p_sens, width = 8.2, height = 5.8, dpi = 400, bg = "white"
)
ggsave(
  file.path(out_fig_dir, "CosMx_dual_high_OR_threshold_sensitivity.pdf"),
  p_sens, width = 8.2, height = 5.8, device = cairo_pdf, bg = "white"
)

# ------------------------------
# 7. Console summary
# ------------------------------
message("\nStep 15 finished.")
message("State cells: ", nrow(state_keep))
message("Coordinate-merged cells: ", nrow(plot_df))
message("Dual-high OR: ", sprintf("%.3f", or_main$odds_ratio),
        " (95% CI ", sprintf("%.3f", or_main$conf_low), "–",
        sprintf("%.3f", or_main$conf_high), "), P = ",
        format(or_main$p_value, scientific = TRUE, digits = 3))
message("Tables: ", out_table_dir)
message("Figures: ", out_fig_dir)


## Public sequential end gate
.step02_out <- file.path(
  out_table_dir,
  "CosMx_state_dual_high_with_coordinates.csv"
)
if (!file.exists(.step02_out)) stop("STEP 02 output missing: ", .step02_out, call. = FALSE)
.step02_n <- nrow(utils::read.csv(.step02_out, check.names = FALSE))
cat("\n============================================================\n")
cat("STEP 02 PASS — state scores + raw metadata -> spatial table\n")
cat("Spatial-table rows: ", .step02_n, "\n", sep = "")
cat("Output: ", .step02_out, "\n", sep = "")
cat("============================================================\n")
