############################################################
## 03_CosMx_metadata_QC_and_technical_sensitivity.R
##
## PUBLIC SEQUENTIAL ROLE
##   Step02 spatial/state table + raw Run5611_MK3 metadata QC
##     -> metadata-QC v2 enhanced cell table + technical sensitivity
##
## IMPORTANT
##   This v2 enhanced table is the locked primary source authority used by
##   the final canonical S21 and S22 reporting scripts.
##
## Public execution edits only:
##   * deterministic Step02 state table
##   * deterministic raw metadata file
## Scientific/statistical logic is unchanged from frozen 19B authority.
############################################################

# ============================================================
# 19B_CosMx_find_best_metadata_QC_join_and_rerun_sensitivity.R
# Purpose:
#   Robustly identify the correct CosMx metadata/QC table, merge it with
#   CosMx state-score + coordinate table, and rerun spatial-neighborhood
#   dual-high QC sensitivity analysis.
#
# Key improvements vs earlier script:
#   1) Audits ALL candidate metadata files and multiple key-normalization strategies.
#   2) Does NOT accept a join with zero overlap.
#   3) Requires usable, non-missing QC columns before calling it QC-residualized analysis.
#   4) Preserves original spatial_x/spatial_y columns from the state-coordinate file.
#   5) Outputs detailed join audit tables for manual troubleshooting.
# ============================================================

############################################################
## v1.1 PUBLIC AUTHORITY NOTE
##
## Historical 19B metadata-v2 is the final S21/S22 source authority.
## Public output names are neutralized to avoid reader confusion.
##
## Historical 19C metadata-v3 is not regenerated: it is a later
## sensitivity/standardization branch and is not consumed by final S21/S22.
## The required technical sensitivity is already produced by this script.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
})

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir,winslash="/",mustWork=TRUE)

# ---- Optional manual overrides ----
state_file <- file.path(
  project_dir,
  "results", "tables", "CosMx_spatial_validation",
  "CosMx_state_dual_high_with_coordinates.csv"
)
metadata_file <- file.path(
  project_dir,
  "data_raw", "CosMx_melanocytic_tumors_Dryad",
  "Slide_4", "Run5611_MK3",
  "Run5611_MK3_metadata_file.csv"
)

manual_state_cell_id_col <- NA_character_
manual_state_fov_col     <- NA_character_
manual_meta_cell_id_col  <- NA_character_
manual_meta_fov_col      <- NA_character_

manual_x_col <- NA_character_
manual_y_col <- NA_character_
manual_myeloid_state_col <- NA_character_
manual_tumor_state_col   <- NA_character_

# If automatic QC detection misses known columns, add them here.
manual_qc_cols <- character(0)

# Sensitivity settings
threshold_grid <- c(0.20, 0.25, 0.30)
k_grid <- c(10, 20, 30)
n_perm <- 999
set.seed(20260507)

out_dir <- file.path(project_dir, "results", "tables", "CosMx_metadata_QC")
fig_dir <- file.path(project_dir, "results", "diagnostics", "11_CosMx", "03_metadata_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv2 <- function(x, filename) readr::write_csv(x, file.path(out_dir, filename), na = "")

norm_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    tolower()
}

norm_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "\\.0$", "")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  tolower(x)
}

numeric_id <- function(x) {
  x <- as.character(x)
  # Prefer last number in an ID such as c_12_345 or cell_345.
  nums <- stringr::str_extract_all(x, "[0-9]+")
  vapply(nums, function(z) if (length(z) == 0) NA_character_ else tail(z, 1), character(1))
}

first_hit <- function(nms, patterns, require_numeric = FALSE, df = NULL) {
  nn <- norm_name(nms)
  for (pat in patterns) {
    idx <- which(stringr::str_detect(nn, pat))
    if (length(idx) > 0) {
      if (require_numeric && !is.null(df)) {
        idx <- idx[vapply(df[idx], function(z) suppressWarnings(mean(!is.na(as.numeric(as.character(z)))) > 0.5), logical(1))]
      }
      if (length(idx) > 0) return(nms[idx[1]])
    }
  }
  NA_character_
}

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(readr::read_csv(path, show_col_types = FALSE, guess_max = 100000))
  if (ext %in% c("tsv", "txt")) return(readr::read_tsv(path, show_col_types = FALSE, guess_max = 100000))
  stop("Unsupported file type: ", path)
}

# ---- 1) Load state-coordinate file ----
if (is.na(state_file)) {
  cand <- list.files(project_dir, pattern = "CosMx.*state.*dual.*high.*coordinate.*CLEAN\\.csv$|Step15_CosMx_state_dual_high_with_coordinates_CLEAN\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(cand) == 0) stop("No CosMx state-coordinate file found. Set state_file manually.")
  # Prefer Step15 cleaned file.
  cand <- cand[order(!grepl("Step15_CosMx_state_dual_high_with_coordinates_CLEAN", cand), nchar(cand))]
  state_file <- cand[1]
}
message("Loading CosMx state-coordinate file: ", state_file)
state_df <- read_any(state_file)
state_nms <- names(state_df)

state_cell_col <- if (!is.na(manual_state_cell_id_col)) manual_state_cell_id_col else first_hit(state_nms, c("^cell_id$", "cell.*id", "cellid", "barcode"))
state_fov_col  <- if (!is.na(manual_state_fov_col)) manual_state_fov_col else first_hit(state_nms, c("^fov$", "field.*view", "fov_id"))
x_col <- if (!is.na(manual_x_col)) manual_x_col else first_hit(state_nms, c("spatial_x", "global_x", "x_global", "center_x", "centroid_x", "^x$"), TRUE, state_df)
y_col <- if (!is.na(manual_y_col)) manual_y_col else first_hit(state_nms, c("spatial_y", "global_y", "y_global", "center_y", "centroid_y", "^y$"), TRUE, state_df)
myeloid_col <- if (!is.na(manual_myeloid_state_col)) manual_myeloid_state_col else first_hit(state_nms, c("myeloid.*treg.*immunosuppressive", "myeloid.*treg"), TRUE, state_df)
tumor_col   <- if (!is.na(manual_tumor_state_col)) manual_tumor_state_col else first_hit(state_nms, c("tumor.*dedifferentiation.*stromal.*remodeling", "dediff.*stromal", "tumor.*stromal"), TRUE, state_df)

state_col_audit <- tibble(
  role = c("state_cell_id", "state_fov", "x", "y", "myeloid_state", "tumor_state"),
  selected_col = c(state_cell_col, state_fov_col, x_col, y_col, myeloid_col, tumor_col)
)
write_csv2(state_col_audit, "QC01_state_file_selected_columns.csv")
write_csv2(tibble(column_index = seq_along(state_nms), column_name = state_nms, normalized = norm_name(state_nms)), "QC01A_state_file_column_names.csv")

if (any(is.na(c(state_cell_col, x_col, y_col, myeloid_col, tumor_col)))) {
  stop("Required state-coordinate columns not detected. Inspect QC01_state_file_selected_columns.csv and set manual overrides.")
}

# ---- 2) Search candidate metadata files ----
if (is.na(metadata_file)) {
  meta_cand <- list.files(project_dir, pattern = "metadata.*\\.(csv|tsv|txt)$|processed_metadata.*\\.(csv|tsv|txt)$|.*metadata_file.*\\.(csv|tsv|txt)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  meta_cand <- meta_cand[!grepl("CosMx_state_dual_high|QC|summary|audit|candidate", meta_cand, ignore.case = TRUE)]
} else {
  meta_cand <- metadata_file
}
if (length(meta_cand) == 0) stop("No candidate CosMx metadata files found. Set metadata_file manually.")

qc_patterns <- c(
  "ncount", "n_count", "total.*count", "total.*transcript", "transcript.*count", "num.*transcript", "ntranscript", "counts$",
  "nfeature", "n_feature", "detected.*gene", "gene.*detected", "num.*gene", "ngene",
  "cell.*area", "area", "nucleus.*area", "nuc.*area",
  "transcript.*area", "count.*area", "density"
)

get_qc_cols <- function(df) {
  nms <- names(df); nn <- norm_name(nms)
  hits <- unique(unlist(lapply(qc_patterns, function(p) nms[stringr::str_detect(nn, p)])))
  hits <- unique(c(hits, manual_qc_cols))
  hits <- hits[hits %in% names(df)]
  # Keep numeric-like columns with at least some non-missing values.
  hits <- hits[vapply(df[hits], function(z) {
    suppressWarnings({
      zz <- as.numeric(as.character(z))
      sum(!is.na(zz)) > 0
    })
  }, logical(1))]
  hits
}

make_keys <- function(df, cell_col, fov_col = NA_character_) {
  cell_raw <- if (!is.na(cell_col) && cell_col %in% names(df)) df[[cell_col]] else rep(NA_character_, nrow(df))
  fov_raw  <- if (!is.na(fov_col) && fov_col %in% names(df)) df[[fov_col]] else rep(NA_character_, nrow(df))
  tibble(
    key_cell_norm = norm_id(cell_raw),
    key_cell_numeric = numeric_id(cell_raw),
    key_fov_cell_norm = if (!is.na(fov_col) && fov_col %in% names(df)) paste0(norm_id(fov_raw), "__", norm_id(cell_raw)) else NA_character_,
    key_fov_cell_numeric = if (!is.na(fov_col) && fov_col %in% names(df)) paste0(numeric_id(fov_raw), "__", numeric_id(cell_raw)) else NA_character_
  )
}

state_keys <- make_keys(state_df, state_cell_col, state_fov_col)

candidate_audits <- list()
metadata_objects <- list()

for (p in meta_cand) {
  message("Auditing metadata candidate: ", p)
  md <- tryCatch(read_any(p), error = function(e) NULL)
  if (is.null(md)) next
  md_nms <- names(md)
  md_cell_col <- first_hit(md_nms, c("^cell_id$", "cell.*id", "cellid", "barcode"))
  md_fov_col  <- first_hit(md_nms, c("^fov$", "field.*view", "fov_id"))
  md_qc_cols <- get_qc_cols(md)
  md_keys <- make_keys(md, md_cell_col, md_fov_col)
  overlap_tbl <- tibble(
    file = p,
    metadata_n = nrow(md),
    meta_cell_col = md_cell_col,
    meta_fov_col = md_fov_col,
    qc_cols_detected = paste(md_qc_cols, collapse = ";"),
    n_qc_cols_detected = length(md_qc_cols),
    key_strategy = c("cell_norm", "cell_numeric", "fov_cell_norm", "fov_cell_numeric"),
    overlap_n = c(
      length(intersect(state_keys$key_cell_norm, md_keys$key_cell_norm)),
      length(intersect(state_keys$key_cell_numeric, md_keys$key_cell_numeric)),
      length(intersect(state_keys$key_fov_cell_norm, md_keys$key_fov_cell_norm)),
      length(intersect(state_keys$key_fov_cell_numeric, md_keys$key_fov_cell_numeric))
    )
  )
  # Estimate nonmissing QC after possible best key merge later; for candidate ranking count in metadata only.
  if (length(md_qc_cols) > 0) {
    qc_nonmiss <- sum(vapply(md[md_qc_cols], function(z) suppressWarnings(sum(!is.na(as.numeric(as.character(z))))), numeric(1)))
  } else qc_nonmiss <- 0
  overlap_tbl$total_qc_nonmissing_in_metadata <- qc_nonmiss
  candidate_audits[[p]] <- overlap_tbl
  metadata_objects[[p]] <- list(df = md, cell_col = md_cell_col, fov_col = md_fov_col, qc_cols = md_qc_cols, keys = md_keys)
}

candidate_audit <- bind_rows(candidate_audits) %>% arrange(desc(overlap_n), desc(n_qc_cols_detected), desc(total_qc_nonmissing_in_metadata))
write_csv2(candidate_audit, "QC02A_metadata_candidate_join_audit_all.csv")

best <- candidate_audit %>% filter(overlap_n > 0, n_qc_cols_detected > 0) %>% slice_max(overlap_n, n = 1, with_ties = FALSE)
if (nrow(best) == 0) {
  write_csv2(tibble(message = "No metadata candidate had both positive overlap and numeric QC columns. Inspect QC02A_metadata_candidate_join_audit_all.csv."), "QC02B_metadata_join_failure_reason.csv")
  stop("No valid metadata join found. Inspect QC02A_metadata_candidate_join_audit_all.csv and set metadata_file / ID columns manually.")
}

best_file <- best$file[1]
best_strategy <- best$key_strategy[1]
md_info <- metadata_objects[[best_file]]
metadata_df <- md_info$df
metadata_keys <- md_info$keys
metadata_qc_cols <- md_info$qc_cols

message("Selected metadata file: ", best_file)
message("Selected join strategy: ", best_strategy)

state_join <- state_df %>% mutate(.state_row_id = row_number()) %>% bind_cols(state_keys)
metadata_join <- metadata_df %>% mutate(.metadata_row_id = row_number()) %>% bind_cols(metadata_keys)

key_col <- switch(best_strategy,
                  cell_norm = "key_cell_norm",
                  cell_numeric = "key_cell_numeric",
                  fov_cell_norm = "key_fov_cell_norm",
                  fov_cell_numeric = "key_fov_cell_numeric")

# Avoid duplicated metadata keys by keeping the first. Audit duplicates.
md_dup <- metadata_join %>% count(.data[[key_col]], name = "n") %>% filter(!is.na(.data[[key_col]]), n > 1)
write_csv2(md_dup, "QC02B_metadata_duplicate_selected_keys.csv")
metadata_join_dedup <- metadata_join %>% filter(!is.na(.data[[key_col]])) %>% distinct(.data[[key_col]], .keep_all = TRUE)

metadata_cols_keep <- unique(c(key_col, md_info$cell_col, md_info$fov_col, metadata_qc_cols))
metadata_cols_keep <- metadata_cols_keep[metadata_cols_keep %in% names(metadata_join_dedup)]

merged <- state_join %>%
  left_join(metadata_join_dedup %>% select(all_of(metadata_cols_keep)), by = key_col, suffix = c("", ".metadata"))

# Standardize selected QC columns into safe names.
qc_renames <- tibble(original_qc_col = metadata_qc_cols) %>%
  mutate(standard_qc_col = paste0("CosMx_QC_", make.names(original_qc_col)))

for (i in seq_len(nrow(qc_renames))) {
  old <- qc_renames$original_qc_col[i]
  new <- qc_renames$standard_qc_col[i]
  if (old %in% names(merged)) merged[[new]] <- suppressWarnings(as.numeric(as.character(merged[[old]])))
}
standard_qc_cols <- qc_renames$standard_qc_col[qc_renames$standard_qc_col %in% names(merged)]
usable_qc_cols <- standard_qc_cols[vapply(merged[standard_qc_cols], function(z) sum(!is.na(z)) > 100 && stats::sd(z, na.rm = TRUE) > 0, logical(1))]

join_audit <- tibble(
  state_n = nrow(state_df),
  metadata_file = best_file,
  metadata_n = nrow(metadata_df),
  selected_join_strategy = best_strategy,
  selected_key_col = key_col,
  overlap_n = best$overlap_n[1],
  n_rows_after_merge = nrow(merged),
  n_rows_with_any_QC = if (length(usable_qc_cols) > 0) sum(rowSums(!is.na(merged[usable_qc_cols])) > 0) else 0,
  usable_qc_cols = paste(usable_qc_cols, collapse = ";"),
  n_usable_qc_cols = length(usable_qc_cols)
)
write_csv2(join_audit, "QC02C_selected_metadata_join_strategy_audit.csv")
write_csv2(qc_renames, "QC02D_QC_column_standardization_map.csv")

if (length(usable_qc_cols) == 0) {
  stop("Metadata merged but no usable QC columns remained after merge. Inspect QC02C/QC02D and candidate audit.")
}

# Required columns nonmissing audit
req_audit <- tibble(
  n_total = nrow(merged),
  n_nonmissing_x = sum(!is.na(merged[[x_col]])),
  n_nonmissing_y = sum(!is.na(merged[[y_col]])),
  n_nonmissing_myeloid = sum(!is.na(merged[[myeloid_col]])),
  n_nonmissing_tumor = sum(!is.na(merged[[tumor_col]]))
)
for (qc in usable_qc_cols) req_audit[[paste0("n_nonmissing_", qc)]] <- sum(!is.na(merged[[qc]]))
write_csv2(req_audit, "QC03_required_and_QC_nonmissing_audit.csv")

# Correlation state vs QC
corr_tbl <- tidyr::expand_grid(state = c(myeloid_col, tumor_col), qc_metric = usable_qc_cols) %>%
  mutate(
    spearman_rho = map2_dbl(state, qc_metric, ~ suppressWarnings(cor(merged[[.x]], merged[[.y]], method = "spearman", use = "complete.obs"))),
    n_complete = map2_int(state, qc_metric, ~ sum(complete.cases(merged[, c(.x, .y)])))
  )
write_csv2(corr_tbl, "QC04_state_QC_spearman_correlation.csv")

# Residualize state scores against QC metrics.
complete_for_resid <- complete.cases(merged[, c(myeloid_col, tumor_col, usable_qc_cols, x_col, y_col)])
merged$Myeloid_QCresid <- NA_real_
merged$Tumor_QCresid <- NA_real_
if (sum(complete_for_resid) > 100) {
  f_my <- as.formula(paste(myeloid_col, "~", paste(usable_qc_cols, collapse = " + ")))
  f_tu <- as.formula(paste(tumor_col, "~", paste(usable_qc_cols, collapse = " + ")))
  merged$Myeloid_QCresid[complete_for_resid] <- resid(lm(f_my, data = merged[complete_for_resid, ]))
  merged$Tumor_QCresid[complete_for_resid] <- resid(lm(f_tu, data = merged[complete_for_resid, ]))
  resid_status <- "QC_residualization_completed"
} else {
  resid_status <- "QC_residualization_not_enough_complete_cases"
}
write_csv2(tibble(residualization_status = resid_status, n_complete = sum(complete_for_resid), usable_qc_cols = paste(usable_qc_cols, collapse = ";")), "QC05_QC_residualization_audit.csv")

# QC-extreme flag: outside 1st-99th percentile for any usable QC metric.
merged$QC_extreme_flag <- FALSE
for (qc in usable_qc_cols) {
  qs <- quantile(merged[[qc]], probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  if (all(is.finite(qs))) {
    merged$QC_extreme_flag <- merged$QC_extreme_flag | (!is.na(merged[[qc]]) & (merged[[qc]] < qs[1] | merged[[qc]] > qs[2]))
  }
}
write_csv2(tibble(n_total = nrow(merged), n_QC_extreme_flagged = sum(merged$QC_extreme_flag, na.rm = TRUE), fraction_QC_extreme = mean(merged$QC_extreme_flag, na.rm = TRUE)), "QC06_QC_extreme_filter_summary.csv")

# Save enhanced table.
readr::write_csv(merged, file.path(out_dir, "CosMx_state_dual_high_with_coordinates_AND_metadata_QC.csv"), na = "")

# ---- 3) Spatial neighborhood enrichment ----
calc_neighbor_enrichment <- function(df, my_score, tu_score, top_prop = 0.25, k = 20, n_perm = 999) {
  dd <- df %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]), !is.na(.data[[my_score]]), !is.na(.data[[tu_score]]))
  n <- nrow(dd)
  if (n <= k + 10) {
    return(tibble(top_prop = top_prop, k = k, n_cells = n, error = "too_few_complete_cells"))
  }
  my_cut <- quantile(dd[[my_score]], probs = 1 - top_prop, na.rm = TRUE, names = FALSE)
  tu_cut <- quantile(dd[[tu_score]], probs = 1 - top_prop, na.rm = TRUE, names = FALSE)
  my_high <- dd[[my_score]] >= my_cut
  tu_high <- dd[[tu_score]] >= tu_cut
  if (sum(my_high) == 0 || sum(tu_high) == 0) {
    return(tibble(top_prop = top_prop, k = k, n_cells = n, error = "no_high_cells"))
  }
  coords <- as.matrix(dd[, c(x_col, y_col)])
  # Use RANN if available; otherwise FNN.
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(coords, coords, k = k + 1)$nn.idx[, -1, drop = FALSE]
  } else if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knn(coords, k = k)$nn.index
  } else {
    stop("Please install either RANN or FNN for kNN computation.")
  }
  idx_my <- which(my_high)
  obs <- mean(tu_high[as.vector(nn[idx_my, , drop = FALSE])])
  perm_vals <- replicate(n_perm, {
    perm_tu <- sample(tu_high)
    mean(perm_tu[as.vector(nn[idx_my, , drop = FALSE])])
  })
  tibble(
    top_prop = top_prop,
    k = k,
    n_cells = n,
    n_myeloid_high = sum(my_high),
    n_tumor_high = sum(tu_high),
    observed_neighbor_fraction = obs,
    permuted_expected_fraction = mean(perm_vals),
    log2_enrichment = log2((obs + 1e-9) / (mean(perm_vals) + 1e-9)),
    empirical_p_high = (sum(perm_vals >= obs) + 1) / (n_perm + 1),
    n_perm = n_perm,
    error = NA_character_
  )
}

analysis_sets <- list(
  raw_all_cells = merged,
  raw_exclude_QC_extremes = merged %>% filter(!QC_extreme_flag),
  QCresid_all_cells = merged %>% filter(!is.na(Myeloid_QCresid), !is.na(Tumor_QCresid)),
  QCresid_exclude_QC_extremes = merged %>% filter(!QC_extreme_flag, !is.na(Myeloid_QCresid), !is.na(Tumor_QCresid))
)

neigh_results <- list()
for (aset in names(analysis_sets)) {
  for (top in threshold_grid) {
    for (kk in k_grid) {
      message("Running neighborhood sensitivity: ", aset, ", top=", top, ", k=", kk)
      if (grepl("QCresid", aset)) {
        res <- calc_neighbor_enrichment(analysis_sets[[aset]], "Myeloid_QCresid", "Tumor_QCresid", top, kk, n_perm)
        score_type <- "QC_residualized"
      } else {
        res <- calc_neighbor_enrichment(analysis_sets[[aset]], myeloid_col, tumor_col, top, kk, n_perm)
        score_type <- "raw"
      }
      neigh_results[[paste(aset, top, kk, sep = "_")]] <- res %>% mutate(analysis_set = aset, score_type = score_type)
    }
  }
}
neigh_tbl <- bind_rows(neigh_results) %>% relocate(any_of(c("analysis_set", "score_type", "top_prop", "k")))
write_csv2(neigh_tbl, "QC07_CosMx_spatial_neighborhood_enrichment_sensitivity_with_metadataQC.csv")

# Same-cell overlap descriptive for raw top thresholds, not identity proof.
same_overlap <- map_dfr(threshold_grid, function(top) {
  dd <- merged %>% filter(!is.na(.data[[myeloid_col]]), !is.na(.data[[tumor_col]]))
  my_cut <- quantile(dd[[myeloid_col]], probs = 1 - top, na.rm = TRUE, names = FALSE)
  tu_cut <- quantile(dd[[tumor_col]], probs = 1 - top, na.rm = TRUE, names = FALSE)
  tab <- table(my_high = dd[[myeloid_col]] >= my_cut, tumor_high = dd[[tumor_col]] >= tu_cut)
  ft <- fisher.test(tab)
  tibble(
    top_prop = top,
    n_total = nrow(dd),
    n_both_high = sum(dd[[myeloid_col]] >= my_cut & dd[[tumor_col]] >= tu_cut),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    interpretation = "descriptive_same_cell_overlap_only_not_single_cell_identity_proof"
  )
})
write_csv2(same_overlap, "QC08_CosMx_same_cell_overlap_descriptive.csv")

summary_tbl <- tibble(
  item = c(
    "state_file", "metadata_file", "n_cells_state", "n_cells_after_merge", "selected_join_strategy", "metadata_overlap_n",
    "n_usable_qc_cols", "usable_qc_cols", "n_rows_with_any_QC", "residualization_status", "n_QC_extreme_flagged",
    "primary_interpretation", "same_cell_overlap_interpretation"
  ),
  value = c(
    state_file, best_file, nrow(state_df), nrow(merged), best_strategy, best$overlap_n[1],
    length(usable_qc_cols), paste(usable_qc_cols, collapse = ";"), join_audit$n_rows_with_any_QC[1], resid_status, sum(merged$QC_extreme_flag, na.rm = TRUE),
    "dual-high interpreted as spatial-neighborhood co-localization",
    "same-cell overlap is descriptive only, not proof of hybrid single-cell identity"
  )
)
write_csv2(summary_tbl, "QC99_CosMx_metadata_QC_sensitivity_final_summary.csv")

message("Done. Enhanced CosMx metadata-QC sensitivity analysis saved to: ", out_dir)


## Public sequential end gate
.step03_out <- file.path(
  out_dir,
  "CosMx_state_dual_high_with_coordinates_AND_metadata_QC.csv"
)
.step03_tech <- file.path(
  out_dir,
  "QC07_CosMx_spatial_neighborhood_enrichment_sensitivity_with_metadataQC.csv"
)
if (!file.exists(.step03_out)) stop("STEP 03 enhanced table missing: ", .step03_out, call. = FALSE)
if (!file.exists(.step03_tech)) stop("STEP 03 technical sensitivity missing: ", .step03_tech, call. = FALSE)
.step03_n <- nrow(utils::read.csv(.step03_out, check.names = FALSE))
cat("\n============================================================\n")
cat("STEP 03 PASS — primary metadata-QC + technical sensitivity rebuilt\n")
cat("Matched metadata-QC cells: ", .step03_n, "\n", sep = "")
cat("Expected frozen S21/S22 universe: 86572\n")
cat("============================================================\n")
