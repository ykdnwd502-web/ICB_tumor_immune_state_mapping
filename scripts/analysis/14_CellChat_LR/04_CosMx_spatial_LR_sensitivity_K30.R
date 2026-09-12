############################################################
## 04_CosMx_spatial_LR_sensitivity_K30.R
##
## Public sequential CosMx spatial ligand-receptor proximity
## analysis; k = 30 (sensitivity).
##
## SCIENTIFIC CONTRACT — UNCHANGED
##   seed = 20260504
##   k = 30
##   n_perm = 999
##   high = top 25% among expressing cells
##   if <50 positive cells: all positive cells are high
##   min eligible edges = 100
##   max_cells_for_knn = Inf
##   BH FDR within this k-specific 60-test universe
##
## Legacy Step15 input filename is retained READ ONLY.
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260504)

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

out_table_dir <- file.path(project_dir, "results", "tables", "CosMx_spatial_LR_validation", "KNN30")
out_fig_dir   <- file.path(project_dir, "results", "figures", "CosMx_spatial_LR_validation", "KNN30")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Inputs from previous steps.
step15_state_coord_file <- file.path(
  project_dir,
  "results", "tables", "CosMx_spatial_validation",
  "Step15_CosMx_state_dual_high_with_coordinates_CLEAN.csv"
)

step16_candidate_pair_file <- file.path(
  project_dir,
  "results", "tables", "CellChat_candidate_LR",
  "CellChat_candidate_pairs_for_CosMx_validation.csv"
)

cosmx_expression_file <- file.path(
  project_dir,
  "data_raw", "CosMx_melanocytic_tumors_Dryad",
  "Slide_4", "Run5611_MK3",
  "Run5611_MK3_exprMat_file.csv"
)

# Spatial-neighbor settings.
# K30 sensitivity version: all output files include K30 to avoid overwriting K10/K20 results.
k_neighbors <- 30

# Permutation settings.
# Use 199 for faster exploratory runs; set to 999 for final.
n_perm <- 999

# Expression-high threshold settings.
# High expression is defined as top 25% among expressing cells.
# If too few cells express a gene, all positive cells are treated as high.
expr_high_top_fraction <- 0.25
min_positive_cells_for_quantile <- 50

# Minimum edge count to report a reliable contrast.
min_eligible_edges <- 100

# To speed debugging, set max_cells_for_knn to e.g. 30000.
# For final analysis, keep Inf.
max_cells_for_knn <- Inf

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

clean_text <- function(x) {
  x <- as.character(x)
  trimws(x)
}

as_logical_robust <- function(x) {
  if (is.logical(x)) return(x)
  xx <- tolower(clean_text(x))
  xx %in% c("true", "t", "1", "yes", "y")
}

make_key_raw <- function(fov, cell_id) {
  paste0(clean_text(fov), "__", clean_text(cell_id))
}

read_csv_safely <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(path, data.table = FALSE, showProgress = FALSE))
  } else {
    read.csv(path, check.names = FALSE)
  }
}

read_selected_expression <- function(path, select_cols) {
  if (!file.exists(path)) stop("CosMx expression matrix not found: ", path)

  # First read header only.
  header <- names(read.csv(path, check.names = FALSE, nrows = 0))
  select_present <- intersect(select_cols, header)
  missing <- setdiff(select_cols, header)

  if (length(missing) > 0) {
    message("Missing expression columns: ", paste(missing, collapse = ", "))
  }

  if (requireNamespace("data.table", quietly = TRUE)) {
    df <- as.data.frame(
      data.table::fread(path, select = select_present, data.table = FALSE, showProgress = FALSE)
    )
  } else {
    message("data.table not installed; falling back to read.csv() for the full expression file.")
    df_all <- read.csv(path, check.names = FALSE)
    df <- df_all[, select_present, drop = FALSE]
    rm(df_all)
    gc()
  }

  list(data = df, header = header, selected = select_present, missing = missing)
}

get_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Could not find any column among: ", paste(candidates, collapse = ", "))
  NA_character_
}

parse_pair_def <- function(pair_family) {
  # Explicit mapping to avoid ambiguous regex parsing.
  map <- list(
    "MIF-CD74" = list(ligand = "MIF", receptor_genes = c("CD74")),
    "MIF-CD44" = list(ligand = "MIF", receptor_genes = c("CD44")),
    "MIF-CXCR4" = list(ligand = "MIF", receptor_genes = c("CXCR4")),
    "SPP1-CD44" = list(ligand = "SPP1", receptor_genes = c("CD44")),
    "SPP1-ITGAV/ITGB1" = list(ligand = "SPP1", receptor_genes = c("ITGAV", "ITGB1")),
    "SPP1-ITGAV/ITGB5" = list(ligand = "SPP1", receptor_genes = c("ITGAV", "ITGB5")),
    "FN1-CD44" = list(ligand = "FN1", receptor_genes = c("CD44")),
    "APP-CD74" = list(ligand = "APP", receptor_genes = c("CD74")),
    "GDF15-TGFBR2" = list(ligand = "GDF15", receptor_genes = c("TGFBR2")),
    "COL1A1-CD44" = list(ligand = "COL1A1", receptor_genes = c("CD44")),
    "COL1A2-CD44" = list(ligand = "COL1A2", receptor_genes = c("CD44")),
    "COL3A1-CD44" = list(ligand = "COL3A1", receptor_genes = c("CD44"))
  )

  if (!pair_family %in% names(map)) {
    parts <- strsplit(pair_family, "-", fixed = TRUE)[[1]]
    if (length(parts) < 2) stop("Cannot parse pair_family: ", pair_family)
    ligand <- parts[1]
    receptor <- paste(parts[-1], collapse = "-")
    receptor_genes <- unlist(strsplit(receptor, "/", fixed = TRUE))
    return(list(ligand = ligand, receptor_genes = receptor_genes))
  }

  map[[pair_family]]
}

calc_complex_expr <- function(df, genes) {
  genes <- unique(genes)
  if (!all(genes %in% colnames(df))) {
    return(rep(NA_real_, nrow(df)))
  }

  mat <- as.data.frame(lapply(genes, function(g) log1p(safe_num(df[[g]]))))
  colnames(mat) <- genes

  if (length(genes) == 1) {
    return(mat[[1]])
  }

  # For receptor complexes, use the minimum expression among subunits.
  # This is conservative and requires all receptor subunits to be present.
  apply(as.matrix(mat), 1, min, na.rm = TRUE)
}

define_high_expr <- function(expr, top_fraction = 0.25, min_positive = 50) {
  expr <- safe_num(expr)
  positive <- is.finite(expr) & expr > 0
  n_pos <- sum(positive)

  if (n_pos == 0) {
    return(list(
      high = rep(FALSE, length(expr)),
      cutoff = NA_real_,
      mode = "no_positive_expression",
      n_positive = 0,
      pct_positive = 0,
      pct_high = 0
    ))
  }

  if (n_pos < min_positive) {
    high <- positive
    return(list(
      high = high,
      cutoff = 0,
      mode = "all_positive_as_high_due_to_low_positive_n",
      n_positive = n_pos,
      pct_positive = mean(positive),
      pct_high = mean(high)
    ))
  }

  q <- stats::quantile(expr[positive], probs = 1 - top_fraction, na.rm = TRUE, names = FALSE)
  high <- positive & expr >= q

  list(
    high = high,
    cutoff = q,
    mode = paste0("top_", top_fraction * 100, "pct_among_positive"),
    n_positive = n_pos,
    pct_positive = mean(positive),
    pct_high = mean(high)
  )
}

build_knn_edges <- function(df, k = 20) {
  coords <- as.matrix(df[, c("spatial_x", "spatial_y")])
  storage.mode(coords) <- "double"

  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knnx(data = coords, query = coords, k = k + 1)
    idx <- nn$nn.index[, -1, drop = FALSE]
    dist <- nn$nn.dist[, -1, drop = FALSE]
  } else if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = coords, query = coords, k = k + 1)
    idx <- nn$nn.idx[, -1, drop = FALSE]
    dist <- nn$nn.dists[, -1, drop = FALSE]
  } else {
    stop("Please install either FNN or RANN for kNN spatial graph construction.")
  }

  data.frame(
    source = rep(seq_len(nrow(df)), each = k),
    target = as.vector(t(idx)),
    distance = as.vector(t(dist))
  )
}

run_one_lr_contrast <- function(ligand_high, receptor_high,
                                source_mask, target_mask,
                                edges, n_perm = 999,
                                min_edges = 100) {
  edge_keep <- source_mask[edges$source] & target_mask[edges$target]
  n_edges <- sum(edge_keep)

  n_source <- sum(source_mask)
  n_target <- sum(target_mask)

  if (n_edges < min_edges || n_source < 5 || n_target < 5) {
    return(data.frame(
      n_source_cells = n_source,
      n_target_cells = n_target,
      n_eligible_edges = n_edges,
      observed_fraction = NA_real_,
      expected_fraction_mean = NA_real_,
      expected_fraction_sd = NA_real_,
      enrichment_ratio = NA_real_,
      log2_spatial_enrichment = NA_real_,
      empirical_p = NA_real_,
      ligand_high_source_pct = mean(ligand_high[source_mask], na.rm = TRUE),
      receptor_high_target_pct = mean(receptor_high[target_mask], na.rm = TRUE),
      status = "too_few_edges"
    ))
  }

  src_edge <- edges$source[edge_keep]
  tgt_edge <- edges$target[edge_keep]

  obs_vec <- ligand_high[src_edge] & receptor_high[tgt_edge]
  obs_frac <- mean(obs_vec, na.rm = TRUE)

  source_pool <- which(source_mask)
  target_pool <- which(target_mask)

  ligand_source <- ligand_high[source_pool]
  receptor_target <- receptor_high[target_pool]

  src_pos <- match(src_edge, source_pool)
  tgt_pos <- match(tgt_edge, target_pool)

  perm_frac <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    lig_perm <- sample(ligand_source, length(ligand_source), replace = FALSE)
    rec_perm <- sample(receptor_target, length(receptor_target), replace = FALSE)
    perm_frac[i] <- mean(lig_perm[src_pos] & rec_perm[tgt_pos], na.rm = TRUE)
  }

  eps <- 1e-10
  expected_mean <- mean(perm_frac, na.rm = TRUE)
  empirical_p <- (1 + sum(perm_frac >= obs_frac, na.rm = TRUE)) / (n_perm + 1)

  data.frame(
    n_source_cells = n_source,
    n_target_cells = n_target,
    n_eligible_edges = n_edges,
    observed_fraction = obs_frac,
    expected_fraction_mean = expected_mean,
    expected_fraction_sd = sd(perm_frac, na.rm = TRUE),
    enrichment_ratio = (obs_frac + eps) / (expected_mean + eps),
    log2_spatial_enrichment = log2((obs_frac + eps) / (expected_mean + eps)),
    empirical_p = empirical_p,
    ligand_high_source_pct = mean(ligand_source, na.rm = TRUE),
    receptor_high_target_pct = mean(receptor_target, na.rm = TRUE),
    status = "ok"
  )
}

# ------------------------------
# 2. Read inputs
# ------------------------------
if (!file.exists(step15_state_coord_file)) {
  stop("Step15 state-coordinate file not found: ", step15_state_coord_file)
}
if (!file.exists(step16_candidate_pair_file)) {
  stop("Step16 candidate pair file not found: ", step16_candidate_pair_file)
}
if (!file.exists(cosmx_expression_file)) {
  stop("CosMx expression matrix not found: ", cosmx_expression_file)
}

state_df <- read_csv_safely(step15_state_coord_file)
pair_df <- read_csv_safely(step16_candidate_pair_file)

required_state_cols <- c("fov", "cell_ID", "spatial_x", "spatial_y", "dual_high_group")
missing_state_cols <- setdiff(required_state_cols, colnames(state_df))
if (length(missing_state_cols) > 0) {
  stop("Missing required columns in Step15 table: ", paste(missing_state_cols, collapse = ", "))
}

# Standardize columns.
state_df <- state_df %>%
  mutate(
    fov = clean_text(fov),
    cell_ID = clean_text(cell_ID),
    .key_raw = make_key_raw(fov, cell_ID),
    spatial_x = safe_num(spatial_x),
    spatial_y = safe_num(spatial_y),
    dual_high_group = as.character(dual_high_group)
  ) %>%
  filter(is.finite(spatial_x), is.finite(spatial_y))

# Robust logical flags.
if ("Myeloid_Treg_high" %in% colnames(state_df)) {
  state_df$Myeloid_Treg_high <- as_logical_robust(state_df$Myeloid_Treg_high)
} else {
  state_df$Myeloid_Treg_high <- state_df$dual_high_group %in% c("Myeloid–Treg high only", "Both high")
}

if ("Dediff_Stromal_high" %in% colnames(state_df)) {
  state_df$Dediff_Stromal_high <- as_logical_robust(state_df$Dediff_Stromal_high)
} else {
  state_df$Dediff_Stromal_high <- state_df$dual_high_group %in% c("Dediff/Stromal high only", "Both high")
}

# Optional debugging downsample for kNN.
if (is.finite(max_cells_for_knn) && nrow(state_df) > max_cells_for_knn) {
  set.seed(20260504)
  state_df <- state_df %>% slice_sample(n = max_cells_for_knn)
}

# Pair definitions.
if (!"pair_family" %in% colnames(pair_df)) {
  stop("Step16 candidate file must contain pair_family.")
}

pair_df <- pair_df %>%
  mutate(
    pair_family = as.character(pair_family),
    nominated_by_cellchat = if ("nominated_by_cellchat" %in% colnames(.)) as_logical_robust(nominated_by_cellchat) else NA,
    validation_priority = if ("validation_priority" %in% colnames(.)) as.character(validation_priority) else NA_character_
  )

pair_defs <- lapply(pair_df$pair_family, parse_pair_def)
pair_df$ligand_gene <- vapply(pair_defs, function(x) x$ligand, character(1))
pair_df$receptor_genes <- vapply(pair_defs, function(x) paste(x$receptor_genes, collapse = ";"), character(1))

genes_needed <- unique(c(
  "fov", "cell_ID",
  pair_df$ligand_gene,
  unlist(strsplit(pair_df$receptor_genes, ";", fixed = TRUE))
))
genes_needed <- genes_needed[!is.na(genes_needed) & genes_needed != ""]

# ------------------------------
# 3. Read expression matrix and merge
# ------------------------------
expr_read <- read_selected_expression(cosmx_expression_file, genes_needed)
expr_df <- expr_read$data

if (!all(c("fov", "cell_ID") %in% colnames(expr_df))) {
  stop("Expression matrix must contain fov and cell_ID.")
}

expr_df <- expr_df %>%
  mutate(
    fov = clean_text(fov),
    cell_ID = clean_text(cell_ID),
    .key_raw = make_key_raw(fov, cell_ID)
  )

gene_overlap <- data.frame(
  requested_gene_or_col = genes_needed,
  found_in_expression_matrix = genes_needed %in% expr_read$header,
  stringsAsFactors = FALSE
)

safe_write_csv(gene_overlap, file.path(out_table_dir, "CosMx_LR_gene_overlap_audit.csv"))

merged <- state_df %>%
  inner_join(expr_df, by = ".key_raw", suffix = c("", ".expr"))

# Use state fov/cell_ID if duplicated after join.
if ("fov.expr" %in% colnames(merged)) merged$fov.expr <- NULL
if ("cell_ID.expr" %in% colnames(merged)) merged$cell_ID.expr <- NULL

input_audit <- data.frame(
  Step15_state_coordinate_file = step15_state_coord_file,
  Step16_candidate_pair_file = step16_candidate_pair_file,
  CosMx_expression_file = cosmx_expression_file,
  State_coordinate_cells = nrow(state_df),
  Expression_cells_read = nrow(expr_df),
  Merged_cells = nrow(merged),
  k_neighbors = k_neighbors,
  n_perm = n_perm,
  expr_high_top_fraction = expr_high_top_fraction,
  min_positive_cells_for_quantile = min_positive_cells_for_quantile,
  max_cells_for_knn = as.character(max_cells_for_knn),
  stringsAsFactors = FALSE
)

safe_write_csv(input_audit, file.path(out_table_dir, "CosMx_LR_input_audit.csv"))

if (nrow(merged) < 1000) {
  stop("Too few merged cells for Step17: ", nrow(merged), ". Inspect Step17 input audit.")
}

# ------------------------------
# 4. Build expression-high flags for ligands and receptors
# ------------------------------
expr_threshold_list <- list()
pair_ready <- pair_df

ligand_high_map <- list()
receptor_high_map <- list()

for (i in seq_len(nrow(pair_df))) {
  pf <- pair_df$pair_family[i]
  ligand <- pair_df$ligand_gene[i]
  receptors <- strsplit(pair_df$receptor_genes[i], ";", fixed = TRUE)[[1]]

  ligand_available <- ligand %in% colnames(merged)
  receptor_available <- all(receptors %in% colnames(merged))

  if (ligand_available) {
    lig_expr <- calc_complex_expr(merged, ligand)
    lig_high <- define_high_expr(
      lig_expr,
      top_fraction = expr_high_top_fraction,
      min_positive = min_positive_cells_for_quantile
    )
  } else {
    lig_expr <- rep(NA_real_, nrow(merged))
    lig_high <- list(
      high = rep(FALSE, nrow(merged)),
      cutoff = NA_real_,
      mode = "ligand_missing",
      n_positive = 0,
      pct_positive = NA_real_,
      pct_high = NA_real_
    )
  }

  if (receptor_available) {
    rec_expr <- calc_complex_expr(merged, receptors)
    rec_high <- define_high_expr(
      rec_expr,
      top_fraction = expr_high_top_fraction,
      min_positive = min_positive_cells_for_quantile
    )
  } else {
    rec_expr <- rep(NA_real_, nrow(merged))
    rec_high <- list(
      high = rep(FALSE, nrow(merged)),
      cutoff = NA_real_,
      mode = "receptor_missing",
      n_positive = 0,
      pct_positive = NA_real_,
      pct_high = NA_real_
    )
  }

  ligand_high_map[[pf]] <- lig_high$high
  receptor_high_map[[pf]] <- rec_high$high

  expr_threshold_list[[i]] <- data.frame(
    pair_family = pf,
    ligand_gene = ligand,
    receptor_genes = paste(receptors, collapse = ";"),
    ligand_available = ligand_available,
    receptor_available = receptor_available,
    ligand_cutoff_log1p = lig_high$cutoff,
    receptor_cutoff_log1p = rec_high$cutoff,
    ligand_threshold_mode = lig_high$mode,
    receptor_threshold_mode = rec_high$mode,
    ligand_n_positive = lig_high$n_positive,
    receptor_n_positive = rec_high$n_positive,
    ligand_pct_positive = lig_high$pct_positive,
    receptor_pct_positive = rec_high$pct_positive,
    ligand_pct_high = lig_high$pct_high,
    receptor_pct_high = rec_high$pct_high,
    pair_usable = ligand_available && receptor_available,
    stringsAsFactors = FALSE
  )
}

expr_threshold_tbl <- bind_rows(expr_threshold_list)
safe_write_csv(expr_threshold_tbl, file.path(out_table_dir, "CosMx_LR_expression_threshold_audit.csv"))

usable_pairs <- expr_threshold_tbl %>%
  filter(pair_usable) %>%
  pull(pair_family)

if (length(usable_pairs) == 0) {
  stop("No usable L-R pairs after gene overlap. Inspect CosMx_LR_expression_threshold_audit.csv.")
}

# ------------------------------
# 5. Build spatial kNN graph
# ------------------------------
message("Building kNN spatial graph; cells = ", nrow(merged), ", k = ", k_neighbors)
edges <- build_knn_edges(merged, k = k_neighbors)

knn_audit <- data.frame(
  n_cells = nrow(merged),
  k_neighbors = k_neighbors,
  n_directed_edges = nrow(edges),
  median_neighbor_distance = median(edges$distance, na.rm = TRUE),
  mean_neighbor_distance = mean(edges$distance, na.rm = TRUE),
  max_neighbor_distance = max(edges$distance, na.rm = TRUE),
  stringsAsFactors = FALSE
)

safe_write_csv(knn_audit, file.path(out_table_dir, "CosMx_spatial_knn_audit.csv"))

# ------------------------------
# 6. Define focused source-target contrasts
# ------------------------------
contrast_defs <- list(
  list(
    contrast = "Dediff/Stromal-high → Myeloid–Treg-high",
    source_mask = merged$Dediff_Stromal_high,
    target_mask = merged$Myeloid_Treg_high
  ),
  list(
    contrast = "Myeloid–Treg-high → Dediff/Stromal-high",
    source_mask = merged$Myeloid_Treg_high,
    target_mask = merged$Dediff_Stromal_high
  ),
  list(
    contrast = "Strict Dediff-only → Myeloid-only",
    source_mask = merged$dual_high_group == "Dediff/Stromal high only",
    target_mask = merged$dual_high_group == "Myeloid–Treg high only"
  ),
  list(
    contrast = "Strict Myeloid-only → Dediff-only",
    source_mask = merged$dual_high_group == "Myeloid–Treg high only",
    target_mask = merged$dual_high_group == "Dediff/Stromal high only"
  ),
  list(
    contrast = "Within Both-high niche",
    source_mask = merged$dual_high_group == "Both high",
    target_mask = merged$dual_high_group == "Both high"
  )
)

contrast_audit <- bind_rows(lapply(contrast_defs, function(x) {
  edge_keep <- x$source_mask[edges$source] & x$target_mask[edges$target]
  data.frame(
    contrast = x$contrast,
    n_source_cells = sum(x$source_mask),
    n_target_cells = sum(x$target_mask),
    n_eligible_edges = sum(edge_keep),
    stringsAsFactors = FALSE
  )
}))

safe_write_csv(contrast_audit, file.path(out_table_dir, "CosMx_LR_contrast_edge_audit.csv"))

# ------------------------------
# 7. Spatial L-R validation
# ------------------------------
message("Running spatial L-R validation...")
res_list <- list()
idx <- 1

for (pf in usable_pairs) {
  lig_high <- ligand_high_map[[pf]]
  rec_high <- receptor_high_map[[pf]]

  pair_meta <- pair_df %>% filter(pair_family == pf) %>% slice(1)
  thresh_meta <- expr_threshold_tbl %>% filter(pair_family == pf) %>% slice(1)

  for (ct in contrast_defs) {
    message("  Pair: ", pf, " | contrast: ", ct$contrast)

    one <- run_one_lr_contrast(
      ligand_high = lig_high,
      receptor_high = rec_high,
      source_mask = ct$source_mask,
      target_mask = ct$target_mask,
      edges = edges,
      n_perm = n_perm,
      min_edges = min_eligible_edges
    )

    res_list[[idx]] <- cbind(
      data.frame(
        pair_family = pf,
        ligand_gene = pair_meta$ligand_gene,
        receptor_genes = pair_meta$receptor_genes,
        nominated_by_cellchat = pair_meta$nominated_by_cellchat,
        validation_priority = pair_meta$validation_priority,
        cellchat_max_prob = if ("max_prob" %in% colnames(pair_meta)) pair_meta$max_prob else NA_real_,
        cellchat_top_source_target = if ("top_source_target" %in% colnames(pair_meta)) pair_meta$top_source_target else NA_character_,
        contrast = ct$contrast,
        k_neighbors = k_neighbors,
        n_perm = n_perm,
        stringsAsFactors = FALSE
      ),
      one,
      data.frame(
        ligand_pct_positive_all = thresh_meta$ligand_pct_positive,
        receptor_pct_positive_all = thresh_meta$receptor_pct_positive,
        ligand_pct_high_all = thresh_meta$ligand_pct_high,
        receptor_pct_high_all = thresh_meta$receptor_pct_high,
        stringsAsFactors = FALSE
      )
    )

    idx <- idx + 1
  }
}

res <- bind_rows(res_list)

res <- res %>%
  mutate(
    empirical_p = as.numeric(empirical_p),
    FDR = p.adjust(empirical_p, method = "BH"),
    neighbor_fraction_percent = observed_fraction * 100,
    expected_fraction_percent = expected_fraction_mean * 100,
    interpretation = case_when(
      status != "ok" ~ "not evaluable",
      is.finite(FDR) & FDR < 0.05 & log2_spatial_enrichment > 0 ~ "spatially enriched",
      is.finite(FDR) & FDR < 0.05 & log2_spatial_enrichment < 0 ~ "spatially depleted",
      TRUE ~ "not significant"
    )
  ) %>%
  arrange(FDR, desc(log2_spatial_enrichment), pair_family, contrast)

safe_write_csv(res, file.path(out_table_dir, "CosMx_spatial_LR_validation_all.csv"))

summary_tbl <- res %>%
  group_by(pair_family, ligand_gene, receptor_genes, nominated_by_cellchat, validation_priority) %>%
  summarise(
    n_evaluable_contrasts = sum(status == "ok"),
    n_enriched_contrasts = sum(interpretation == "spatially enriched", na.rm = TRUE),
    max_log2_spatial_enrichment = max(log2_spatial_enrichment, na.rm = TRUE),
    max_neighbor_fraction_percent = max(neighbor_fraction_percent, na.rm = TRUE),
    min_FDR = min(FDR, na.rm = TRUE),
    top_contrast = {
      tmp_score <- ifelse(is.finite(log2_spatial_enrichment), log2_spatial_enrichment, -Inf)
      if (all(is.infinite(tmp_score))) NA_character_ else contrast[which.max(tmp_score)][1]
    },
    .groups = "drop"
  ) %>%
  mutate(
    max_log2_spatial_enrichment = ifelse(is.infinite(max_log2_spatial_enrichment), NA_real_, max_log2_spatial_enrichment),
    max_neighbor_fraction_percent = ifelse(is.infinite(max_neighbor_fraction_percent), NA_real_, max_neighbor_fraction_percent),
    min_FDR = ifelse(is.infinite(min_FDR), NA_real_, min_FDR),
    validation_priority = as.character(validation_priority),
    overall_support = dplyr::case_when(
      n_enriched_contrasts > 0 & validation_priority %in% c("High", "Medium") ~ "CellChat-nominated and spatially enriched",
      n_enriched_contrasts > 0 ~ "Spatially enriched exploratory pair",
      TRUE ~ "No positive spatial enrichment"
    )
  ) %>%
  arrange(desc(n_enriched_contrasts), min_FDR, desc(max_log2_spatial_enrichment))

safe_write_csv(summary_tbl, file.path(out_table_dir, "CosMx_spatial_LR_validation_summary.csv"))

# Positive subset for manuscript inspection.
positive_tbl <- res %>%
  filter(interpretation == "spatially enriched") %>%
  arrange(FDR, desc(log2_spatial_enrichment))

safe_write_csv(positive_tbl, file.path(out_table_dir, "CosMx_spatial_LR_validation_positive_hits.csv"))

# ------------------------------
# 8. Figures
# ------------------------------
plot_res <- res %>%
  filter(status == "ok") %>%
  mutate(
    pair_family = factor(pair_family, levels = rev(unique(pair_df$pair_family))),
    contrast = factor(
      contrast,
      levels = c(
        "Dediff/Stromal-high → Myeloid–Treg-high",
        "Myeloid–Treg-high → Dediff/Stromal-high",
        "Strict Dediff-only → Myeloid-only",
        "Strict Myeloid-only → Dediff-only",
        "Within Both-high niche"
      )
    ),
    significant = ifelse(!is.na(FDR) & FDR < 0.05, "FDR < 0.05", "NS")
  )

p_dot <- ggplot(plot_res, aes(x = contrast, y = pair_family)) +
  geom_point(
    aes(
      size = neighbor_fraction_percent,
      color = log2_spatial_enrichment,
      alpha = significant
    )
  ) +
  scale_color_gradient2(
    low = "#3C5488",
    mid = "white",
    high = "#E64B35",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "log2 spatial\nenrichment"
  ) +
  scale_size_continuous(
    range = c(1.5, 7),
    name = "L-high/R-high\nneighbor fraction (%)"
  ) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "NS" = 0.45), name = NULL) +
  labs(
    title = "Step17 CosMx spatial proximity-constrained L-R validation",
    subtitle = paste0("kNN = ", k_neighbors, "; permutations = ", n_perm,
                      "; CellChat candidates tested in dual-high-related spatial neighborhoods"),
    x = "Spatial source → target contrast",
    y = "Candidate L-R axis"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_fig_dir, "CosMx_spatial_LR_validation_dotplot.png"),
  p_dot, width = 12.5, height = 7.8, dpi = 500, bg = "white"
)
ggsave(
  file.path(out_fig_dir, "CosMx_spatial_LR_validation_dotplot.pdf"),
  p_dot, width = 12.5, height = 7.8, device = cairo_pdf, bg = "white"
)

# Summary bar plot for top positive hits.
top_hits <- positive_tbl %>%
  mutate(hit_label = paste0(pair_family, "\n", contrast)) %>%
  slice_head(n = 20) %>%
  mutate(hit_label = factor(hit_label, levels = rev(hit_label)))

if (nrow(top_hits) > 0) {
  p_top <- ggplot(top_hits, aes(x = log2_spatial_enrichment, y = hit_label, fill = validation_priority)) +
    geom_col(width = 0.72, color = "grey30", linewidth = 0.2) +
    labs(
      title = "Top spatially enriched candidate L-R axes in CosMx",
      subtitle = "Ranked by FDR and log2 spatial enrichment",
      x = "log2 spatial enrichment",
      y = NULL,
      fill = "Priority"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(out_fig_dir, "CosMx_spatial_LR_top_hits.png"),
    p_top, width = 9, height = 7.2, dpi = 450, bg = "white"
  )
  ggsave(
    file.path(out_fig_dir, "CosMx_spatial_LR_top_hits.pdf"),
    p_top, width = 9, height = 7.2, device = cairo_pdf, bg = "white"
  )
}

# ------------------------------
# 9. Console summary
# ------------------------------
message("\nCosMx spatial L-R proximity analysis finished.")
message("Merged CosMx cells: ", nrow(merged))
message("Usable L-R pairs: ", length(usable_pairs), " / ", nrow(pair_df))
message("Directed kNN edges: ", nrow(edges))
message("Positive spatially enriched tests: ", nrow(positive_tbl))
message("Tables: ", out_table_dir)
message("Figures: ", out_fig_dir)
