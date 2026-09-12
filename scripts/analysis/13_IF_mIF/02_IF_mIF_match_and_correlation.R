############################################################
## 02_IF_mIF_match_and_correlation.R
##
## PUBLIC SEQUENTIAL IF/mIF — MATCHED-CELL MERGE + SPEARMAN
##
## Scientific logic is inherited unchanged from frozen Step14A v8:
##   - candidate IF metadata audited by FOV-cell overlap
##   - selected metadata = largest overlap
##   - raw versus normalized composite-key strategy chosen by overlap
##   - minimum matched cells = 1000
##   - protein-marker identification logic unchanged
##   - Spearman cor.test(exact = FALSE)
##   - BH FDR across all protein-feature x state pairs
##
## Historical heatmap production is intentionally omitted here because the
## manuscript-facing reporting figure is rebuilt formally by Step 03 (S23).
##
## Only formal output names/reporting text are changed.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir,winslash="/",mustWork=TRUE)

out_table_dir <- file.path(project_dir, "results", "tables", "IF_mIF_validation")
out_fig_dir   <- file.path(project_dir, "results", "figures", "IF_mIF_validation")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Usually keep this empty. If auto-selection fails, set this manually to the
# metadata file with the highest overlap shown in the v8 audit table.
manual_if_metadata_files <- character(0)

manual_state_score_file <- file.path(
  project_dir,
  "results", "tables", "IF_mIF_validation",
  "IF_mIF_CosMx_cell_state_scores.csv"
)

min_overlap_cells <- 1000

state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_display_labels <- c(
  "Immune-defective/\nCold",
  "Myeloid–Treg\nImmunosuppressive",
  "Tumor dedifferentiation/\nStromal remodeling",
  "Melanocytic\nDifferentiation"
)
names(state_display_labels) <- state_cols

# Same blue-white-red heatmap style used in earlier main figures.
heatmap_low  <- "#3B82F6"
heatmap_mid  <- "white"
heatmap_high <- "#EF4444"

protein_marker_keywords <- c(
  "CD3", "CD45", "S100", "PMEL", "PMEL17", "S100B", "S100b",
  "DAPI", "Vimentin", "VIM", "PCNA", "CD8", "CD8a"
)

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE)
  message("Saved table: ", path)
}

sanitize_names <- function(df) {
  nm <- colnames(df)
  nm[is.na(nm) | nm == ""] <- paste0("blank_col_", which(is.na(nm) | nm == ""))
  colnames(df) <- make.unique(nm, sep = "_")
  df
}

safe_num <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

numeric_fraction <- function(x) {
  xx <- safe_num(x)
  mean(is.finite(xx), na.rm = TRUE)
}

guess_col <- function(cols, patterns) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

clean_text <- function(x) {
  x <- as.character(x)
  trimws(x)
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

make_key_raw <- function(fov, cell_id) {
  paste0(clean_text(fov), "__", clean_text(cell_id))
}

make_key_norm <- function(fov, cell_id) {
  paste0("fov", norm_numeric_like(fov), "__cell", norm_numeric_like(cell_id))
}

identify_protein_features <- function(df, min_numeric_fraction = 0.30) {
  cols <- colnames(df)
  pats <- paste(protein_marker_keywords, collapse = "|")
  marker_cols <- cols[grepl(pats, cols, ignore.case = TRUE)]
  marker_cols <- marker_cols[vapply(marker_cols, function(cc) {
    numeric_fraction(df[[cc]]) >= min_numeric_fraction
  }, logical(1))]
  marker_cols
}

read_metadata_as_character <- function(f) {
  df <- tryCatch(
    utils::read.csv(
      f,
      check.names = FALSE,
      colClasses = "character",
      na.strings = c("", "NA", "NaN", "NULL", "null")
    ),
    error = function(e) NULL
  )
  if (is.null(df)) return(NULL)
  df <- sanitize_names(df)
  df$IF_source_file <- basename(f)
  df
}

sig_label <- function(fdr) {
  ifelse(is.na(fdr), "",
         ifelse(fdr < 0.001, "***",
                ifelse(fdr < 0.01, "**",
                       ifelse(fdr < 0.05, "*", ""))))
}

interpret_feature <- function(feature, state, rho) {
  f <- tolower(feature)
  s <- tolower(state)
  if (is.na(rho)) return("insufficient data")
  if (grepl("dapi", f)) return("QC/image-related feature; interpret cautiously")
  if (grepl("cd3", f) && grepl("myeloid_treg", s) && rho > 0) return("supportive immune-related association")
  if (grepl("cd45", f) && grepl("tumor_dedifferentiation", s) && rho < 0) return("negative pan-leukocyte association; not merely immune-high")
  if ((grepl("s100", f) || grepl("pmel", f)) && grepl("melanocytic", s) && rho > 0) {
    return("lineage-marker association; interpret with targeted-panel overlap")
  }
  if (abs(rho) < 0.10) return("weak association")
  if (abs(rho) < 0.30) return("modest association")
  "moderate association"
}

# ------------------------------
# 1. Read state-score table and construct keys
# ------------------------------
if (!file.exists(manual_state_score_file)) {
  stop("State-score file not found. Run 01_IF_mIF_build_CosMx_state_scores.R first or check manual_state_score_file: ", manual_state_score_file)
}

state_df <- read.csv(manual_state_score_file, check.names = FALSE)
state_df <- sanitize_names(state_df)

st_fov_col <- guess_col(colnames(state_df), c("^fov$", "FOV", "field"))
st_cell_col <- guess_col(colnames(state_df), c("^cell_ID$", "cellid", "cell_id", "CellID", "Cell_ID", "^cell$"))

if (is.na(st_fov_col) || is.na(st_cell_col)) {
  stop("State-score table must contain both fov and cell_ID.")
}

missing_states <- setdiff(state_cols, colnames(state_df))
if (length(missing_states) > 0) {
  stop("Missing state columns in state-score file: ", paste(missing_states, collapse = ", "))
}

state_df$.key_raw  <- make_key_raw(state_df[[st_fov_col]], state_df[[st_cell_col]])
state_df$.key_norm <- make_key_norm(state_df[[st_fov_col]], state_df[[st_cell_col]])

state_keep_raw <- state_df %>%
  dplyr::select(.key_raw, dplyr::all_of(state_cols)) %>%
  dplyr::group_by(.key_raw) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(state_cols), ~ mean(safe_num(.x), na.rm = TRUE)), .groups = "drop")

state_keep_norm <- state_df %>%
  dplyr::select(.key_norm, dplyr::all_of(state_cols)) %>%
  dplyr::group_by(.key_norm) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(state_cols), ~ mean(safe_num(.x), na.rm = TRUE)), .groups = "drop")

for (cc in state_cols) {
  state_keep_raw[[cc]][is.nan(state_keep_raw[[cc]])] <- NA_real_
  state_keep_norm[[cc]][is.nan(state_keep_norm[[cc]])] <- NA_real_
}

# ------------------------------
# 2. Candidate metadata files
# ------------------------------
if (length(manual_if_metadata_files) > 0) {
  if_candidates <- manual_if_metadata_files[file.exists(manual_if_metadata_files)]
} else {
  raw_root <- file.path(project_dir, "data_raw")
  if_candidates <- list.files(
    raw_root,
    pattern = "processed_metadata.*\\.csv$|metadata_file\\.csv$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  # Prefer Slide 4 / Run5611-related files, but keep all matching files for audit.
  if_candidates <- if_candidates[grepl("CosMx|Run5611|Slide_4|Slide4|Dryad|processed_metadata", if_candidates, ignore.case = TRUE)]
}

if (length(if_candidates) == 0) stop("No candidate IF metadata files found.")

audit_list <- list()
meta_cache <- list()

for (f in unique(if_candidates)) {
  md <- read_metadata_as_character(f)

  if (is.null(md)) {
    audit_list[[length(audit_list) + 1]] <- data.frame(
      file = f, basename = basename(f),
      n_rows = NA_integer_, fov_col = NA_character_, cell_col = NA_character_,
      n_protein_features = 0,
      raw_overlap = 0, norm_overlap = 0,
      best_strategy = NA_character_, best_overlap = 0,
      first_5_if_raw = "", first_5_if_norm = "",
      stringsAsFactors = FALSE
    )
    next
  }

  fov_col <- guess_col(colnames(md), c("^fov$", "FOV", "field"))
  cell_col <- guess_col(colnames(md), c("^cell_ID$", "cellid", "cell_id", "CellID", "Cell_ID", "^cell$"))
  features <- identify_protein_features(md)

  if (is.na(fov_col) || is.na(cell_col) || length(features) < 2) {
    audit_list[[length(audit_list) + 1]] <- data.frame(
      file = f, basename = basename(f),
      n_rows = nrow(md), fov_col = fov_col, cell_col = cell_col,
      n_protein_features = length(features),
      raw_overlap = 0, norm_overlap = 0,
      best_strategy = NA_character_, best_overlap = 0,
      first_5_if_raw = "", first_5_if_norm = "",
      stringsAsFactors = FALSE
    )
    next
  }

  md$.key_raw  <- make_key_raw(md[[fov_col]], md[[cell_col]])
  md$.key_norm <- make_key_norm(md[[fov_col]], md[[cell_col]])

  raw_overlap <- length(intersect(unique(md$.key_raw), state_keep_raw$.key_raw))
  norm_overlap <- length(intersect(unique(md$.key_norm), state_keep_norm$.key_norm))

  if (raw_overlap >= norm_overlap) {
    best_strategy <- "raw"
    best_overlap <- raw_overlap
  } else {
    best_strategy <- "norm"
    best_overlap <- norm_overlap
  }

  audit_list[[length(audit_list) + 1]] <- data.frame(
    file = f,
    basename = basename(f),
    n_rows = nrow(md),
    fov_col = fov_col,
    cell_col = cell_col,
    n_protein_features = length(features),
    protein_features = paste(features, collapse = ";"),
    raw_overlap = raw_overlap,
    norm_overlap = norm_overlap,
    best_strategy = best_strategy,
    best_overlap = best_overlap,
    first_5_if_raw = paste(head(unique(md$.key_raw), 5), collapse = ";"),
    first_5_if_norm = paste(head(unique(md$.key_norm), 5), collapse = ";"),
    stringsAsFactors = FALSE
  )

  meta_cache[[f]] <- list(
    data = md,
    fov_col = fov_col,
    cell_col = cell_col,
    features = features,
    raw_overlap = raw_overlap,
    norm_overlap = norm_overlap,
    best_strategy = best_strategy,
    best_overlap = best_overlap
  )
}

audit_df <- bind_rows(audit_list) %>% arrange(desc(best_overlap))
safe_write_csv(audit_df, file.path(out_table_dir, "IF_mIF_metadata_candidate_overlap_audit.csv"))

if (nrow(audit_df) == 0 || max(audit_df$best_overlap, na.rm = TRUE) < min_overlap_cells) {
  stop(
    "No IF metadata candidate had sufficient overlap with the state-score table. ",
    "Inspect IF_mIF_metadata_candidate_overlap_audit.csv."
  )
}

selected_file <- audit_df$file[1]
selected_strategy <- audit_df$best_strategy[1]
message("Selected IF metadata file: ", selected_file)
message("Selected merge-key strategy: ", selected_strategy)
message("Best overlap: ", audit_df$best_overlap[1])

sel <- meta_cache[[selected_file]]
if_meta <- sel$data
if_features <- sel$features
if_fov_col <- sel$fov_col
if_cell_col <- sel$cell_col

# ------------------------------
# 3. Prepare selected metadata and merge
# ------------------------------
if (selected_strategy == "raw") {
  if_meta_reduced <- if_meta %>%
    dplyr::select(.key_raw, dplyr::all_of(if_features)) %>%
    dplyr::group_by(.key_raw) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(if_features), ~ mean(safe_num(.x), na.rm = TRUE)), .groups = "drop")

  for (cc in if_features) if_meta_reduced[[cc]][is.nan(if_meta_reduced[[cc]])] <- NA_real_

  merged <- dplyr::inner_join(if_meta_reduced, state_keep_raw, by = ".key_raw")
  key_col_used <- ".key_raw"
  if_unique_n <- nrow(if_meta_reduced)
  state_unique_n <- nrow(state_keep_raw)
} else {
  if_meta_reduced <- if_meta %>%
    dplyr::select(.key_norm, dplyr::all_of(if_features)) %>%
    dplyr::group_by(.key_norm) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(if_features), ~ mean(safe_num(.x), na.rm = TRUE)), .groups = "drop")

  for (cc in if_features) if_meta_reduced[[cc]][is.nan(if_meta_reduced[[cc]])] <- NA_real_

  merged <- dplyr::inner_join(if_meta_reduced, state_keep_norm, by = ".key_norm")
  key_col_used <- ".key_norm"
  if_unique_n <- nrow(if_meta_reduced)
  state_unique_n <- nrow(state_keep_norm)
}

merge_audit <- data.frame(
  selected_IF_metadata_file = selected_file,
  selected_strategy = selected_strategy,
  key_col_used = key_col_used,
  IF_unique_key_n = if_unique_n,
  State_unique_key_n = state_unique_n,
  Merged_key_n = nrow(merged),
  IF_fov_col = if_fov_col,
  IF_cell_col = if_cell_col,
  State_fov_col = st_fov_col,
  State_cell_col = st_cell_col,
  State_columns = paste(state_cols, collapse = ";"),
  Protein_features = paste(if_features, collapse = ";"),
  stringsAsFactors = FALSE
)
safe_write_csv(merge_audit, file.path(out_table_dir, "IF_mIF_state_merge_audit.csv"))

if (nrow(merged) < min_overlap_cells) {
  stop("Merged IF/state table has too few overlapping records: ", nrow(merged))
}

merged_out <- file.path(out_table_dir, "IF_mIF_state_protein_merged_table.csv")
safe_write_csv(merged, merged_out)

# ------------------------------
# 4. Spearman correlation
# ------------------------------
cor_list <- list()
k <- 1

for (pf in if_features) {
  for (st in state_cols) {
    x <- safe_num(merged[[pf]])
    y <- safe_num(merged[[st]])
    ok <- is.finite(x) & is.finite(y)
    n_ok <- sum(ok)

    if (n_ok >= 10 && length(unique(x[ok])) > 2 && length(unique(y[ok])) > 2) {
      ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      rho <- unname(ct$estimate)
      pval <- ct$p.value
    } else {
      rho <- NA_real_
      pval <- NA_real_
    }

    cor_list[[k]] <- data.frame(
      dataset_source = basename(selected_file),
      analysis_unit = paste0("CosMx_cell_level_", selected_strategy, "_FOV_cell_ID_key"),
      protein_feature = pf,
      state_score = st,
      n = n_ok,
      rho = rho,
      p_value = pval,
      stringsAsFactors = FALSE
    )
    k <- k + 1
  }
}

cor_tbl <- bind_rows(cor_list) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    significance = sig_label(FDR),
    interpretation_note = mapply(interpret_feature, protein_feature, state_score, rho)
  ) %>%
  arrange(state_score, protein_feature)

safe_write_csv(cor_tbl, file.path(out_table_dir, "IF_mIF_protein_state_spearman.csv"))

############################################################
## 5. End gate
############################################################

expected_corr_rows <- length(if_features) * length(state_cols)

if (nrow(cor_tbl) != expected_corr_rows) {
  stop(
    "Unexpected IF/mIF Spearman row count: ",
    nrow(cor_tbl),
    "; expected ",
    expected_corr_rows,
    ".",
    call. = FALSE
  )
}

message("IF/mIF Step 02 matched-cell merge + correlation finished.")
message("Selected IF metadata file: ", selected_file)
message("Selected strategy: ", selected_strategy)
message("Merged n: ", nrow(merged))
message("Merged table: ", merged_out)
message(
  "Correlation table: ",
  file.path(
    out_table_dir,
    "IF_mIF_protein_state_spearman.csv"
  )
)
