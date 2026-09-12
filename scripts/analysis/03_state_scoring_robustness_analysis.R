############################################################
## 03_state_scoring_robustness_analysis.R
##
## Purpose:
##   Assess robustness of the predefined tumor–immune state scoring framework
##   across GSE244982 and GSE78220.
##
## Main analyses:
##   1) Gene-set retention and final-state union-gene-pool audit.
##   2) Pairwise final-state gene-pool overlap (Jaccard index).
##   3) Concordance of alternative scoring methods with the frozen primary
##      constituent z-mean score.
##   4) Leave-one-signature-out robustness.
##   5) Random gene-subsampling robustness (70%, 80%, 90% retained; 100
##      iterations per retained fraction).
##   6) Descriptive GSE78220 non-response AUC robustness across scoring methods.
##      Confidence intervals use 10,000 stratified bootstrap replicates
##      with seed 20260506, matching script 02.
##
## Important definition:
##   In the gene-subsampling analysis, each subsampled union-gene score is
##   correlated with the corresponding FULL PRIMARY CONSTITUENT Z-MEAN SCORE.
##   The full union-gene score is not used as the reference.
##
## Reproducibility:
##   Run after scripts 01 and 02. Run from the repository root, or set
##   ICB_PROJECT_DIR. Package versions should be restored from renv.lock.
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

############################################################
## 0. Project directory
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- getwd()
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

message("Project directory: ", project_dir)

data_processed_dir <- file.path(project_dir, "data_processed")
table_dir          <- file.path(project_dir, "results/tables")
figure_dir         <- file.path(project_dir, "results/figures")
log_dir            <- file.path(project_dir, "logs")

out_table_dir <- file.path(table_dir, "state_scoring_robustness")
out_fig_dir   <- file.path(figure_dir, "state_scoring_robustness")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "tibble",
  "purrr",
  "pROC"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
    ". Restore the repository environment with renv::restore() before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(purrr)
  library(pROC)
})

theme_robustness <- function(base_size = 12.5, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2.5, color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "black"),
      axis.title = element_text(face = "bold", size = base_size, color = "black"),
      axis.text = element_text(color = "black", size = base_size - 1.5),
      legend.title = element_text(face = "bold", size = base_size, color = "black"),
      legend.text = element_text(size = base_size - 1.5, color = "black"),
      strip.text = element_text(face = "bold", size = base_size - 0.5, color = "black"),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.60),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

############################################################
## 2. Utility functions
############################################################

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
  message("Wrote: ", file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = file)
  message("Wrote: ", file)
}

safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 320) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file, plot = plot, width = width, height = height, dpi = dpi)
  message("Wrote: ", file)
}

find_first_existing <- function(candidates) {
  candidates <- unique(unlist(candidates))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}


clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\|.*$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

normalize_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

zvec <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

make_sample_key <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("\\.BASELINE$|\\.PRE$|\\.PRETREATMENT$|\\.ON$|\\.ON_TREATMENT$", "", x, ignore.case = TRUE)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

collapse_duplicate_genes <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- clean_gene_symbol(rownames(mat))
  
  keep <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep, , drop = FALSE]
  
  dt <- data.table(GeneSymbol = rownames(mat), mat, check.names = FALSE)
  dt2 <- dt[, lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)), by = GeneSymbol]
  genes <- dt2$GeneSymbol
  mat2 <- as.matrix(dt2[, !"GeneSymbol"])
  rownames(mat2) <- genes
  storage.mode(mat2) <- "numeric"
  mat2
}

read_expression_matrix <- function(file) {
  message("Reading expression matrix: ", file)
  
  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    obj <- readRDS(file)
    
    if (is.matrix(obj) || is.data.frame(obj)) {
      mat <- as.matrix(obj)
      if (is.null(rownames(mat))) {
        stop("RDS object has no rownames; cannot identify gene symbols: ", file)
      }
      storage.mode(mat) <- "numeric"
      return(collapse_duplicate_genes(mat))
    }
    
    stop("Unsupported RDS expression object type: ", paste(class(obj), collapse = ", "))
  }
  
  if (grepl("\\.csv$|\\.tsv$|\\.txt$", file, ignore.case = TRUE)) {
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE), "\t", ",")
    dt <- data.table::fread(file, sep = sep, data.table = FALSE, check.names = FALSE)
    
    cn <- colnames(dt)
    cn_norm <- tolower(gsub("[^a-z0-9]", "", cn))
    gene_col_idx <- which(cn_norm %in% c("genesymbol", "gene", "genes", "symbol", "hgncsymbol"))
    
    if (length(gene_col_idx) == 0) {
      gene_col_idx <- 1
      warning("No explicit gene column detected. Using first column as gene symbol: ", cn[1])
    } else {
      gene_col_idx <- gene_col_idx[1]
    }
    
    genes <- dt[[gene_col_idx]]
    expr_df <- dt[, -gene_col_idx, drop = FALSE]
    
    numeric_ok <- sapply(expr_df, function(x) {
      xx <- suppressWarnings(as.numeric(as.character(x)))
      sum(!is.na(xx)) > 0
    })
    
    expr_df <- expr_df[, numeric_ok, drop = FALSE]
    expr_mat <- as.matrix(sapply(expr_df, function(x) as.numeric(as.character(x))))
    colnames(expr_mat) <- colnames(expr_df)
    rownames(expr_mat) <- genes
    
    return(collapse_duplicate_genes(expr_mat))
  }
  
  stop("Unsupported expression file format: ", file)
}

read_table_any <- function(file) {
  message("Reading table: ", file)
  
  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    obj <- readRDS(file)
    if (is.data.frame(obj)) return(as.data.frame(obj))
    if (is.matrix(obj)) return(as.data.frame(obj))
    stop("Unsupported RDS table object: ", paste(class(obj), collapse = ", "))
  }
  
  if (grepl("\\.csv$|\\.tsv$|\\.txt$", file, ignore.case = TRUE)) {
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE), "\t", ",")
    return(data.table::fread(file, sep = sep, data.table = FALSE, check.names = FALSE))
  }
  
  stop("Unsupported table file format: ", file)
}

read_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file, warn = FALSE)
  lines <- lines[nchar(lines) > 0]
  
  gene_sets <- list()
  
  for (ln in lines) {
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 3) next
    
    set_name <- parts[1]
    genes <- parts[-c(1, 2)]
    genes <- unique(clean_gene_symbol(genes))
    genes <- genes[!is.na(genes) & genes != ""]
    
    gene_sets[[set_name]] <- genes
  }
  
  gene_sets
}

standardize_response_columns <- function(df) {
  if (!"Sample" %in% colnames(df)) {
    df$Sample <- rownames(df)
  }
  
  if (!"ResponseGroup" %in% colnames(df)) {
    possible <- c("response_group", "BinaryResponse", "binary_response", "Response", "response", "BestResponse", "best_response")
    hit <- intersect(possible, colnames(df))
    if (length(hit) > 0) {
      colnames(df)[match(hit[1], colnames(df))] <- "ResponseGroup"
    }
  }
  
  df
}

############################################################
## 3. Canonical input files
############################################################

gmt_candidates <- c(
  Sys.getenv("ICB_CLEAN_GMT"),
  file.path(
    project_dir,
    "results/tables/ICBcomb_final_input_gene_sets_CLEAN/02_final_ICBcomb_gene_sets/final_ICBcomb_gene_sets_CLEAN.gmt"
  )
)
gmt_candidates <- unique(gmt_candidates[nzchar(gmt_candidates)])

gse244982_expr_candidates <- c(
  file.path(data_processed_dir, "GSE244982_bulk_matrix.rds"),
  file.path(table_dir, "bulk_discovery/GSE244982_bulk_expression_matrix.csv")
)

gse78220_expr_candidates <- c(
  file.path(data_processed_dir, "GSE78220_external_expression_matrix_for_scoring.rds")
)

gse78220_state_candidates <- c(
  file.path(table_dir, "external_validation/GSE78220_external_state_scores.csv"),
  file.path(data_processed_dir, "GSE78220_external_state_scores.rds")
)

gmt_file <- find_first_existing(gmt_candidates)
gse244982_expr_file <- find_first_existing(gse244982_expr_candidates)
gse78220_expr_file  <- find_first_existing(gse78220_expr_candidates)
gse78220_state_file <- find_first_existing(gse78220_state_candidates)

input_inventory <- data.frame(
  Input = c(
    "Frozen tumor-immune GMT",
    "GSE244982 expression",
    "GSE78220 expression",
    "GSE78220 state/response table"
  ),
  File = c(
    gmt_file,
    gse244982_expr_file,
    gse78220_expr_file,
    gse78220_state_file
  ),
  Exists = !is.na(c(gmt_file, gse244982_expr_file, gse78220_expr_file, gse78220_state_file)) &
    file.exists(c(gmt_file, gse244982_expr_file, gse78220_expr_file, gse78220_state_file)),
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_inventory,
  file.path(out_table_dir, "input_file_inventory.csv")
)

if (any(!input_inventory$Exists)) {
  stop(
    "One or more canonical inputs were not found. Run scripts 01 and 02 first, ",
    "and verify the paths listed in input_file_inventory.csv."
  )
}

############################################################
## 4. Load gene sets and expression matrices
############################################################

gene_sets <- read_gmt(gmt_file)

gene_set_summary <- data.frame(
  RawSignature = names(gene_sets),
  RawSignature_nGenes = sapply(gene_sets, length),
  stringsAsFactors = FALSE
)

safe_write_csv(
  gene_set_summary,
  file.path(out_table_dir, "raw_GMT_gene_set_summary.csv")
)

expr244 <- read_expression_matrix(gse244982_expr_file)
expr782 <- read_expression_matrix(gse78220_expr_file)

message("GSE244982 expression: ", nrow(expr244), " genes x ", ncol(expr244), " samples")
message("GSE78220 expression: ", nrow(expr782), " genes x ", ncol(expr782), " samples")

############################################################
## 5. Define final-state mapping from raw signatures
############################################################

get_signature_names <- function(gene_sets, target_names, state_name, min_required = 1) {
  available <- names(gene_sets)
  available_norm <- normalize_name(available)
  names(available_norm) <- available
  
  target_norm <- normalize_name(target_names)
  
  matched <- available[available_norm %in% target_norm]
  
  if (length(matched) == 0) {
    for (tg in target_norm) {
      idx <- grep(tg, available_norm, fixed = TRUE)
      if (length(idx) > 0) matched <- c(matched, available[idx])
    }
  }
  
  matched <- unique(matched)
  
  message(state_name, " matched raw signatures: ", paste(matched, collapse = ", "))
  
  if (length(matched) < min_required) {
    stop(
      "Insufficient matched signature columns for ", state_name,
      ". Required at least ", min_required,
      ", found ", length(matched),
      ". Available signatures:\n",
      paste(available, collapse = "\n")
    )
  }
  
  matched
}

immune_sigs <- get_signature_names(
  gene_sets,
  c(
    "Immune_defective_Cold_RESTORE_hybrid",
    "Immune_defective_Cold_RESTORE_data_driven",
    "Immune_defective_Cold_RESTORE_curated",
    "Immune_cold_RESTORE_hybrid",
    "Immune_cold_RESTORE_data_driven",
    "Immune_cold_RESTORE_curated"
  ),
  "Immune_defective_Cold"
)

myeloid_sigs <- get_signature_names(
  gene_sets,
  c(
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated",
    "Myeloid_Treg_SUPPRESS_hybrid",
    "Myeloid_Treg_SUPPRESS_data_driven",
    "Myeloid_Treg_SUPPRESS_curated"
  ),
  "Myeloid_Treg_Immunosuppressive"
)

tumor_sigs <- get_signature_names(
  gene_sets,
  c(
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated",
    "Tumor_dedifferentiation_SUPPRESS_hybrid",
    "Tumor_dedifferentiation_SUPPRESS_data_driven",
    "Tumor_dedifferentiation_SUPPRESS_curated"
  ),
  "Tumor_dedifferentiation_Stromal_remodeling"
)

melanocytic_sigs <- get_signature_names(
  gene_sets,
  c(
    "Melanocytic_Differentiation_REFERENCE_hybrid",
    "Melanocytic_Differentiation_REFERENCE_data_driven",
    "Melanocytic_Differentiation_REFERENCE_curated",
    "Melanocytic_differentiation_REFERENCE_hybrid",
    "Melanocytic_differentiation_REFERENCE_data_driven",
    "Melanocytic_differentiation_REFERENCE_curated"
  ),
  "Melanocytic_Differentiation"
)

state_mapping <- data.frame(
  FinalState = c(
    rep("Immune_defective_Cold", length(immune_sigs)),
    rep("Myeloid_Treg_Immunosuppressive", length(myeloid_sigs)),
    rep("Tumor_dedifferentiation_Stromal_remodeling", length(tumor_sigs)),
    rep("Melanocytic_Differentiation", length(melanocytic_sigs))
  ),
  RawSignature = c(
    immune_sigs,
    myeloid_sigs,
    tumor_sigs,
    melanocytic_sigs
  ),
  DirectionForFinalState = c(
    rep("inverse", length(immune_sigs)),
    rep("positive", length(myeloid_sigs)),
    rep("positive", length(tumor_sigs)),
    rep("positive", length(melanocytic_sigs))
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  state_mapping,
  file.path(out_table_dir, "final_state_raw_signature_mapping_used.csv")
)

state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

## Public-figure display lock. Internal analytical state IDs remain unchanged.
state_display_wrap <- c(
  "Immune_defective_Cold" = "immune-defective/\ncold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic\ndifferentiation"
)

state_direction <- c(
  Immune_defective_Cold = -1,
  Myeloid_Treg_Immunosuppressive = 1,
  Tumor_dedifferentiation_Stromal_remodeling = 1,
  Melanocytic_Differentiation = 1
)

############################################################
## 6. Build final-state gene pools
############################################################

build_state_gene_pools <- function(gene_sets, mapping) {
  pools <- list()
  
  for (st in unique(mapping$FinalState)) {
    sigs <- mapping$RawSignature[mapping$FinalState == st]
    genes <- unique(unlist(gene_sets[sigs]))
    genes <- clean_gene_symbol(genes)
    genes <- genes[!is.na(genes) & genes != ""]
    pools[[st]] <- unique(genes)
  }
  
  pools[state_cols]
}

state_gene_pools <- build_state_gene_pools(gene_sets, state_mapping)

state_gene_pool_summary <- data.frame(
  FinalState = names(state_gene_pools),
  n_union_genes = sapply(state_gene_pools, length),
  DirectionForFinalState = ifelse(names(state_gene_pools) == "Immune_defective_Cold", "inverse", "positive"),
  stringsAsFactors = FALSE
)

safe_write_csv(
  state_gene_pool_summary,
  file.path(out_table_dir, "final_state_union_gene_pool_summary.csv")
)

## Gene retention audit in each dataset.
gene_retention_audit <- function(expr_mat, state_gene_pools, dataset_name) {
  expr_genes <- rownames(expr_mat)
  
  out <- list()
  
  for (st in names(state_gene_pools)) {
    genes <- state_gene_pools[[st]]
    present <- genes %in% expr_genes
    
    out[[length(out) + 1]] <- data.frame(
      Dataset = dataset_name,
      FinalState = st,
      RequestedGene = genes,
      Present = present,
      stringsAsFactors = FALSE
    )
  }
  
  dplyr::bind_rows(out) %>%
    group_by(Dataset, FinalState) %>%
    mutate(
      n_requested = n(),
      n_present = sum(Present),
      n_missing = sum(!Present),
      retention_fraction = n_present / n_requested
    ) %>%
    ungroup()
}

ret244 <- gene_retention_audit(expr244, state_gene_pools, "GSE244982")
ret782 <- gene_retention_audit(expr782, state_gene_pools, "GSE78220")

ret_all <- bind_rows(ret244, ret782)

safe_write_csv(
  ret_all,
  file.path(out_table_dir, "gene_retention_audit_long.csv")
)

ret_summary <- ret_all %>%
  distinct(Dataset, FinalState, n_requested, n_present, n_missing, retention_fraction)

safe_write_csv(
  ret_summary,
  file.path(out_table_dir, "gene_retention_summary_by_dataset_state.csv")
)

############################################################
## 7. Jaccard overlap among final-state gene pools
############################################################

jaccard_pair <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(NA_real_)
  inter / union
}

jaccard_mat <- matrix(
  NA_real_,
  nrow = length(state_cols),
  ncol = length(state_cols),
  dimnames = list(state_cols, state_cols)
)

overlap_long <- list()

for (s1 in state_cols) {
  for (s2 in state_cols) {
    g1 <- state_gene_pools[[s1]]
    g2 <- state_gene_pools[[s2]]
    
    inter <- length(intersect(g1, g2))
    union_n <- length(union(g1, g2))
    jac <- jaccard_pair(g1, g2)
    
    jaccard_mat[s1, s2] <- jac
    
    overlap_long[[length(overlap_long) + 1]] <- data.frame(
      State1 = s1,
      State2 = s2,
      n_State1 = length(g1),
      n_State2 = length(g2),
      n_intersection = inter,
      n_union = union_n,
      Jaccard = jac,
      stringsAsFactors = FALSE
    )
  }
}

overlap_long <- bind_rows(overlap_long)

safe_write_csv(
  as.data.frame(jaccard_mat),
  file.path(out_table_dir, "final_state_gene_pool_Jaccard_matrix.csv"),
  row.names = TRUE
)

safe_write_csv(
  overlap_long,
  file.path(out_table_dir, "final_state_gene_pool_overlap_long.csv")
)

p_jaccard <- overlap_long %>%
  mutate(
    State1 = factor(State1, levels = state_cols),
    State2 = factor(State2, levels = state_cols),
    label = sprintf("%.2f", Jaccard)
  ) %>%
  ggplot(aes(x = State2, y = State1, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = "white", high = "#2166AC", na.value = "grey90") +
  scale_x_discrete(labels = state_display_wrap) +
  scale_y_discrete(labels = state_display_wrap) +
  labs(
    title = "Gene-pool overlap among predefined tumor–immune states",
    subtitle = "Jaccard index based on final-state union gene pools",
    x = NULL,
    y = NULL,
    fill = "Jaccard"
  ) +
  theme_robustness(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10.75)
  )

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_final_state_gene_pool_Jaccard_heatmap.pdf"),
  p_jaccard,
  width = 7,
  height = 5.5
)

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_final_state_gene_pool_Jaccard_heatmap.png"),
  p_jaccard,
  width = 7,
  height = 5.5
)

############################################################
## 8. Scoring functions
############################################################

score_raw_signatures_zmean <- function(expr_mat, gene_sets) {
  expr_mat <- collapse_duplicate_genes(expr_mat)
  expr_genes <- rownames(expr_mat)
  
  expr_z <- t(scale(t(expr_mat)))
  expr_z[!is.finite(expr_z)] <- 0
  
  score_list <- list()
  presence_list <- list()
  
  for (gs in names(gene_sets)) {
    genes <- clean_gene_symbol(gene_sets[[gs]])
    present <- intersect(genes, expr_genes)
    
    presence_list[[gs]] <- data.frame(
      RawSignature = gs,
      Gene = genes,
      Present = genes %in% expr_genes,
      stringsAsFactors = FALSE
    )
    
    if (length(present) == 0) {
      score_list[[gs]] <- rep(NA_real_, ncol(expr_z))
    } else {
      score_list[[gs]] <- colMeans(expr_z[present, , drop = FALSE], na.rm = TRUE)
    }
  }
  
  raw_scores <- as.data.frame(score_list, check.names = FALSE)
  raw_scores$Sample <- colnames(expr_mat)
  raw_scores <- raw_scores[, c("Sample", setdiff(colnames(raw_scores), "Sample")), drop = FALSE]
  
  list(
    raw_scores = raw_scores,
    gene_presence = bind_rows(presence_list)
  )
}

build_final_scores_from_raw <- function(raw_scores, mapping, method_label) {
  out <- data.frame(Sample = raw_scores$Sample, stringsAsFactors = FALSE)
  
  for (st in state_cols) {
    sigs <- mapping$RawSignature[mapping$FinalState == st]
    direction <- unique(mapping$DirectionForFinalState[mapping$FinalState == st])
    
    if (length(sigs) == 0) {
      out[[st]] <- NA_real_
      next
    }
    
    raw_mat <- as.matrix(raw_scores[, sigs, drop = FALSE])
    storage.mode(raw_mat) <- "numeric"
    
    sign <- ifelse(direction[1] == "inverse", -1, 1)
    out[[st]] <- zvec(sign * rowMeans(raw_mat, na.rm = TRUE))
  }
  
  out$ScoringMethod <- method_label
  out
}

score_state_pools_zmean <- function(expr_mat, state_gene_pools, method_label = "zmean_union_gene_pool") {
  expr_mat <- collapse_duplicate_genes(expr_mat)
  expr_z <- t(scale(t(expr_mat)))
  expr_z[!is.finite(expr_z)] <- 0
  
  out <- data.frame(Sample = colnames(expr_mat), stringsAsFactors = FALSE)
  
  for (st in state_cols) {
    genes <- intersect(state_gene_pools[[st]], rownames(expr_z))
    sign <- state_direction[[st]]
    
    if (length(genes) == 0) {
      out[[st]] <- NA_real_
    } else {
      out[[st]] <- zvec(sign * colMeans(expr_z[genes, , drop = FALSE], na.rm = TRUE))
    }
  }
  
  out$ScoringMethod <- method_label
  out
}

score_state_pools_rank_mean <- function(expr_mat, state_gene_pools, method_label = "rank_mean_singscore_like") {
  expr_mat <- collapse_duplicate_genes(expr_mat)
  
  rank_mat <- apply(expr_mat, 2, function(v) {
    rank(v, ties.method = "average", na.last = "keep")
  })
  
  rank_mat <- as.matrix(rank_mat)
  rownames(rank_mat) <- rownames(expr_mat)
  colnames(rank_mat) <- colnames(expr_mat)
  
  ## Scale ranks to approximately 0-1 within each sample.
  rank_scaled <- sweep(rank_mat, 2, colSums(!is.na(rank_mat)), "/")
  
  out <- data.frame(Sample = colnames(expr_mat), stringsAsFactors = FALSE)
  
  for (st in state_cols) {
    genes <- intersect(state_gene_pools[[st]], rownames(rank_scaled))
    sign <- state_direction[[st]]
    
    if (length(genes) == 0) {
      out[[st]] <- NA_real_
    } else {
      out[[st]] <- zvec(sign * colMeans(rank_scaled[genes, , drop = FALSE], na.rm = TRUE))
    }
  }
  
  out$ScoringMethod <- method_label
  out
}

score_state_pools_gsva_ssgsea_optional <- function(expr_mat, state_gene_pools, method_label = "GSVA_ssGSEA_optional") {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    message("Package GSVA is not installed. Skipping optional ssGSEA scoring.")
    return(NULL)
  }
  
  expr_mat <- collapse_duplicate_genes(expr_mat)
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  
  gene_pools_present <- lapply(state_gene_pools, function(g) intersect(g, rownames(expr_mat)))
  gene_pools_present <- gene_pools_present[sapply(gene_pools_present, length) > 1]
  
  if (length(gene_pools_present) == 0) {
    message("No gene pools with enough genes for GSVA ssGSEA.")
    return(NULL)
  }
  
  gsva_res <- tryCatch(
    {
      GSVA::gsva(
        expr = expr_mat,
        gset.idx.list = gene_pools_present,
        method = "ssgsea",
        kcdf = "Gaussian",
        verbose = FALSE
      )
    },
    error = function(e1) {
      message("Old GSVA::gsva interface failed: ", e1$message)
      tryCatch(
        {
          param <- GSVA::ssgseaParam(
            exprData = expr_mat,
            geneSets = gene_pools_present
          )
          GSVA::gsva(param, verbose = FALSE)
        },
        error = function(e2) {
          message("New GSVA interface also failed: ", e2$message)
          NULL
        }
      )
    }
  )
  
  if (is.null(gsva_res)) return(NULL)
  
  gsva_res <- as.matrix(gsva_res)
  
  out <- data.frame(Sample = colnames(expr_mat), stringsAsFactors = FALSE)
  
  for (st in state_cols) {
    if (!st %in% rownames(gsva_res)) {
      out[[st]] <- NA_real_
    } else {
      sign <- state_direction[[st]]
      out[[st]] <- zvec(sign * as.numeric(gsva_res[st, ]))
    }
  }
  
  out$ScoringMethod <- method_label
  out
}

############################################################
## 9. Compute primary and alternative scores
############################################################

raw244 <- score_raw_signatures_zmean(expr244, gene_sets)
raw782 <- score_raw_signatures_zmean(expr782, gene_sets)

safe_write_csv(
  raw244$raw_scores,
  file.path(out_table_dir, "GSE244982_raw_signature_scores_zmean.csv")
)

safe_write_csv(
  raw782$raw_scores,
  file.path(out_table_dir, "GSE78220_raw_signature_scores_zmean.csv")
)

primary244 <- build_final_scores_from_raw(raw244$raw_scores, state_mapping, "primary_constituent_zmean")
primary782 <- build_final_scores_from_raw(raw782$raw_scores, state_mapping, "primary_constituent_zmean")

union_z244 <- score_state_pools_zmean(expr244, state_gene_pools, "union_gene_zmean")
union_z782 <- score_state_pools_zmean(expr782, state_gene_pools, "union_gene_zmean")

rank244 <- score_state_pools_rank_mean(expr244, state_gene_pools, "rank_mean_singscore_like")
rank782 <- score_state_pools_rank_mean(expr782, state_gene_pools, "rank_mean_singscore_like")

gsva244 <- score_state_pools_gsva_ssgsea_optional(expr244, state_gene_pools, "GSVA_ssGSEA_optional")
gsva782 <- score_state_pools_gsva_ssgsea_optional(expr782, state_gene_pools, "GSVA_ssGSEA_optional")

scores244_all <- bind_rows(
  primary244,
  union_z244,
  rank244,
  gsva244
) %>%
  mutate(Dataset = "GSE244982")

scores782_all <- bind_rows(
  primary782,
  union_z782,
  rank782,
  gsva782
) %>%
  mutate(Dataset = "GSE78220")

safe_write_csv(
  scores244_all,
  file.path(out_table_dir, "GSE244982_scores_by_scoring_method.csv")
)

safe_write_csv(
  scores782_all,
  file.path(out_table_dir, "GSE78220_scores_by_scoring_method.csv")
)

############################################################
## 10. Scoring-method concordance
############################################################

compute_method_concordance <- function(score_df, dataset_name) {
  out <- list()
  
  methods <- unique(score_df$ScoringMethod)
  
  for (st in state_cols) {
    for (m in methods) {
      if (m == "primary_constituent_zmean") next
      
      dat_primary <- score_df %>%
        filter(ScoringMethod == "primary_constituent_zmean") %>%
        select(Sample, primary_score = all_of(st))
      
      dat_alt <- score_df %>%
        filter(ScoringMethod == m) %>%
        select(Sample, alt_score = all_of(st))
      
      dat <- left_join(dat_primary, dat_alt, by = "Sample")
      
      ok <- is.finite(dat$primary_score) & is.finite(dat$alt_score)
      
      if (sum(ok) < 4) {
        rho <- NA_real_
        pval <- NA_real_
        pear <- NA_real_
        pvalp <- NA_real_
      } else {
        ct <- suppressWarnings(cor.test(dat$primary_score[ok], dat$alt_score[ok], method = "spearman", exact = FALSE))
        ct2 <- suppressWarnings(cor.test(dat$primary_score[ok], dat$alt_score[ok], method = "pearson"))
        rho <- unname(ct$estimate)
        pval <- ct$p.value
        pear <- unname(ct2$estimate)
        pvalp <- ct2$p.value
      }
      
      out[[length(out) + 1]] <- data.frame(
        Dataset = dataset_name,
        FinalState = st,
        AlternativeScoringMethod = m,
        n = sum(ok),
        Spearman_rho_vs_primary = rho,
        Spearman_p = pval,
        Pearson_r_vs_primary = pear,
        Pearson_p = pvalp,
        stringsAsFactors = FALSE
      )
    }
  }
  
  bind_rows(out)
}

method_conc <- bind_rows(
  compute_method_concordance(scores244_all, "GSE244982"),
  compute_method_concordance(scores782_all, "GSE78220")
)

safe_write_csv(
  method_conc,
  file.path(out_table_dir, "scoring_method_concordance_vs_primary.csv")
)

p_method_conc <- method_conc %>%
  mutate(
    Dataset = factor(Dataset, levels = c("GSE244982", "GSE78220")),
    FinalState = factor(FinalState, levels = state_cols),
    label = ifelse(is.na(Spearman_rho_vs_primary), "", sprintf("%.2f", Spearman_rho_vs_primary))
  ) %>%
  ggplot(aes(x = AlternativeScoringMethod, y = FinalState, fill = Spearman_rho_vs_primary)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3) +
  facet_wrap(~ Dataset, ncol = 1) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  scale_y_discrete(labels = state_display_wrap) +
  labs(
    title = "Alternative scoring-method concordance",
    subtitle = "Spearman correlation versus the primary constituent z-mean score",
    x = "Alternative scoring method",
    y = "Final tumor–immune state",
    fill = "Spearman\nrho"
  ) +
  theme_robustness(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.75))

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_scoring_method_concordance_heatmap.pdf"),
  p_method_conc,
  width = 8.5,
  height = 6
)

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_scoring_method_concordance_heatmap.png"),
  p_method_conc,
  width = 8.5,
  height = 6
)

############################################################
## 11. Leave-one-signature-out robustness
############################################################

leave_one_signature_out <- function(raw_scores, mapping, dataset_name) {
  full_scores <- build_final_scores_from_raw(raw_scores, mapping, "primary_constituent_zmean")
  out_conc <- list()
  out_scores <- list()
  
  for (st in state_cols) {
    sigs <- mapping$RawSignature[mapping$FinalState == st]
    
    if (length(sigs) < 2) {
      out_conc[[length(out_conc) + 1]] <- data.frame(
        Dataset = dataset_name,
        FinalState = st,
        OmittedRawSignature = NA_character_,
        n_remaining_signatures = length(sigs),
        n = NA_integer_,
        Spearman_rho_vs_full = NA_real_,
        Spearman_p = NA_real_,
        Status = "Skipped: only one constituent signature",
        stringsAsFactors = FALSE
      )
      next
    }
    
    for (omit in sigs) {
      mapping_sub <- mapping %>%
        filter(!(FinalState == st & RawSignature == omit))
      
      sub_scores <- build_final_scores_from_raw(
        raw_scores = raw_scores,
        mapping = mapping_sub,
        method_label = paste0("leave_one_out_", omit)
      )
      
      dat <- full_scores %>%
        select(Sample, full_score = all_of(st)) %>%
        left_join(
          sub_scores %>% select(Sample, loo_score = all_of(st)),
          by = "Sample"
        )
      
      ok <- is.finite(dat$full_score) & is.finite(dat$loo_score)
      
      if (sum(ok) < 4) {
        rho <- NA_real_
        pval <- NA_real_
      } else {
        ct <- suppressWarnings(cor.test(dat$full_score[ok], dat$loo_score[ok], method = "spearman", exact = FALSE))
        rho <- unname(ct$estimate)
        pval <- ct$p.value
      }
      
      out_conc[[length(out_conc) + 1]] <- data.frame(
        Dataset = dataset_name,
        FinalState = st,
        OmittedRawSignature = omit,
        n_remaining_signatures = length(sigs) - 1,
        n = sum(ok),
        Spearman_rho_vs_full = rho,
        Spearman_p = pval,
        Status = "Completed",
        stringsAsFactors = FALSE
      )
      
      tmp_scores <- sub_scores %>%
        select(Sample, all_of(st)) %>%
        rename(LOO_score = all_of(st)) %>%
        mutate(
          Dataset = dataset_name,
          FinalState = st,
          OmittedRawSignature = omit
        )
      
      out_scores[[length(out_scores) + 1]] <- tmp_scores
    }
  }
  
  list(
    concordance = bind_rows(out_conc),
    scores = bind_rows(out_scores)
  )
}

loo244 <- leave_one_signature_out(raw244$raw_scores, state_mapping, "GSE244982")
loo782 <- leave_one_signature_out(raw782$raw_scores, state_mapping, "GSE78220")

loo_conc <- bind_rows(loo244$concordance, loo782$concordance)

safe_write_csv(
  loo_conc,
  file.path(out_table_dir, "leave_one_signature_out_concordance.csv")
)

p_loo <- loo_conc %>%
  filter(Status == "Completed") %>%
  mutate(
    Dataset = factor(Dataset, levels = c("GSE244982", "GSE78220")),
    FinalState = factor(FinalState, levels = state_cols)
  ) %>%
  ggplot(aes(x = FinalState, y = Spearman_rho_vs_full)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8) +
  facet_wrap(~ Dataset, ncol = 1) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(labels = state_display_wrap) +
  labs(
    title = "Leave-one-signature-out robustness",
    subtitle = "Concordance of leave-one-signature-out scores with full composite scores",
    x = "Final tumor–immune state",
    y = "Spearman rho vs full primary score"
  ) +
  theme_robustness(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.75))

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_leave_one_signature_out_robustness.pdf"),
  p_loo,
  width = 8,
  height = 6
)

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_leave_one_signature_out_robustness.png"),
  p_loo,
  width = 8,
  height = 6
)

############################################################
## 12. Gene subsampling robustness
############################################################

score_single_state_zmean_from_genes <- function(expr_mat, genes, st) {
  expr_mat <- collapse_duplicate_genes(expr_mat)
  expr_z <- t(scale(t(expr_mat)))
  expr_z[!is.finite(expr_z)] <- 0
  
  present <- intersect(clean_gene_symbol(genes), rownames(expr_z))
  sign <- state_direction[[st]]
  
  if (length(present) == 0) {
    return(rep(NA_real_, ncol(expr_z)))
  }
  
  zvec(sign * colMeans(expr_z[present, , drop = FALSE], na.rm = TRUE))
}

run_gene_subsampling <- function(expr_mat, state_gene_pools, full_scores, dataset_name,
                                 fractions = c(0.70, 0.80, 0.90),
                                 n_iter = 100,
                                 seed = 20260608) {
  set.seed(seed)
  
  out <- list()
  
  for (st in state_cols) {
    genes_all <- intersect(state_gene_pools[[st]], rownames(expr_mat))
    n_all <- length(genes_all)
    
    full_primary_dat <- full_scores %>%
      filter(ScoringMethod == "primary_constituent_zmean") %>%
      select(Sample, full_score = all_of(st))

    if (nrow(full_primary_dat) != ncol(expr_mat)) {
      stop(
        "Gene-subsampling reference mismatch for ", dataset_name, " / ", st,
        ": expected one full primary constituent score per expression sample."
      )
    }
    
    for (frac in fractions) {
      n_keep <- max(2, floor(n_all * frac))
      
      for (iter in seq_len(n_iter)) {
        if (n_all < 3) {
          out[[length(out) + 1]] <- data.frame(
            Dataset = dataset_name,
            FinalState = st,
            FractionRetained = frac,
            Iteration = iter,
            n_available_genes = n_all,
            n_sampled_genes = NA_integer_,
            Spearman_rho_vs_full = NA_real_,
            Spearman_p = NA_real_,
            Status = "Skipped: too few available genes",
            stringsAsFactors = FALSE
          )
          next
        }
        
        sampled <- sample(genes_all, size = n_keep, replace = FALSE)
        sub_score <- score_single_state_zmean_from_genes(expr_mat, sampled, st)
        
        sub_df <- data.frame(
          Sample = colnames(expr_mat),
          sub_score = sub_score,
          stringsAsFactors = FALSE
        )
        
        dat <- full_primary_dat %>%
          left_join(sub_df, by = "Sample")
        
        ok <- is.finite(dat$full_score) & is.finite(dat$sub_score)
        
        if (sum(ok) < 4) {
          rho <- NA_real_
          pval <- NA_real_
        } else {
          ct <- suppressWarnings(cor.test(dat$full_score[ok], dat$sub_score[ok], method = "spearman", exact = FALSE))
          rho <- unname(ct$estimate)
          pval <- ct$p.value
        }
        
        out[[length(out) + 1]] <- data.frame(
          Dataset = dataset_name,
          FinalState = st,
          FractionRetained = frac,
          Iteration = iter,
          n_available_genes = n_all,
          n_sampled_genes = n_keep,
          Spearman_rho_vs_full = rho,
          Spearman_p = pval,
          Status = "Completed",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  bind_rows(out)
}

sub244 <- run_gene_subsampling(
  expr_mat = expr244,
  state_gene_pools = state_gene_pools,
  full_scores = scores244_all,
  dataset_name = "GSE244982",
  fractions = c(0.70, 0.80, 0.90),
  n_iter = 100,
  seed = 20260608
)

sub782 <- run_gene_subsampling(
  expr_mat = expr782,
  state_gene_pools = state_gene_pools,
  full_scores = scores782_all,
  dataset_name = "GSE78220",
  fractions = c(0.70, 0.80, 0.90),
  n_iter = 100,
  seed = 20260609
)

sub_all <- bind_rows(sub244, sub782)

safe_write_csv(
  sub_all,
  file.path(out_table_dir, "gene_subsampling_all_iterations.csv")
)

sub_summary <- sub_all %>%
  filter(Status == "Completed") %>%
  group_by(Dataset, FinalState, FractionRetained) %>%
  summarise(
    n_iterations = n(),
    median_rho = median(Spearman_rho_vs_full, na.rm = TRUE),
    IQR_low = quantile(Spearman_rho_vs_full, 0.25, na.rm = TRUE),
    IQR_high = quantile(Spearman_rho_vs_full, 0.75, na.rm = TRUE),
    min_rho = min(Spearman_rho_vs_full, na.rm = TRUE),
    max_rho = max(Spearman_rho_vs_full, na.rm = TRUE),
    .groups = "drop"
  )

safe_write_csv(
  sub_summary,
  file.path(out_table_dir, "gene_subsampling_stability_summary.csv")
)

p_sub <- sub_all %>%
  filter(Status == "Completed") %>%
  mutate(
    Dataset = factor(Dataset, levels = c("GSE244982", "GSE78220")),
    FinalState = factor(FinalState, levels = state_cols),
    FractionRetained = factor(paste0(round(FractionRetained * 100), "% retained"),
                              levels = c("70% retained", "80% retained", "90% retained"))
  ) %>%
  ggplot(aes(x = FractionRetained, y = Spearman_rho_vs_full)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  facet_grid(
    Dataset ~ FinalState,
    labeller = labeller(FinalState = state_display_wrap)
  ) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Gene-subsampling robustness of predefined state scores",
    subtitle = "100 random subsampling iterations per retained-gene fraction",
    x = "Gene-retention fraction",
    y = "Spearman rho vs full primary score"
  ) +
  theme_robustness(base_size = 11.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25))

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_gene_subsampling_robustness_boxplots.pdf"),
  p_sub,
  width = 12,
  height = 6.5
)

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_gene_subsampling_robustness_boxplots.png"),
  p_sub,
  width = 12,
  height = 6.5
)

############################################################
## 13. Optional GSE78220 AUC robustness across scoring methods
############################################################

## Match the frozen GSE78220 Figure 7 bootstrap protocol exactly.
## The primary constituent z-mean rows therefore reproduce the AUC/95% CI
## estimates from 02_GSE78220_state_scoring_response_survival.R.
auc_bootstrap_seed <- 20260506L
auc_bootstrap_n <- 10000L

calc_auc <- function(df, predictor_col, label) {
  if (!"ResponseGroup" %in% colnames(df)) {
    return(data.frame(
      Analysis = label,
      n_total = NA_integer_,
      n_responder = NA_integer_,
      n_nonresponder = NA_integer_,
      AUC = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      Status = "Skipped: ResponseGroup column unavailable",
      stringsAsFactors = FALSE
    ))
  }
  
  dat <- df %>%
    filter(ResponseGroup %in% c("Responder", "NonResponder")) %>%
    mutate(ResponseGroup = factor(ResponseGroup, levels = c("Responder", "NonResponder"))) %>%
    filter(is.finite(.data[[predictor_col]]))
  
  n_total <- nrow(dat)
  n_resp <- sum(dat$ResponseGroup == "Responder")
  n_non <- sum(dat$ResponseGroup == "NonResponder")
  
  if (n_total < 6 || n_resp < 2 || n_non < 2) {
    return(data.frame(
      Analysis = label,
      n_total = n_total,
      n_responder = n_resp,
      n_nonresponder = n_non,
      AUC = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      Status = "Skipped: insufficient binary samples",
      stringsAsFactors = FALSE
    ))
  }
  
  roc_obj <- pROC::roc(
    response = dat$ResponseGroup,
    predictor = dat[[predictor_col]],
    levels = c("Responder", "NonResponder"),
    direction = "<",
    quiet = TRUE
  )
  
  ci_obj <- tryCatch(
    pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = auc_bootstrap_n, stratified = TRUE),
    error = function(e) c(NA_real_, NA_real_, NA_real_)
  )
  
  data.frame(
    Analysis = label,
    n_total = n_total,
    n_responder = n_resp,
    n_nonresponder = n_non,
    AUC = as.numeric(pROC::auc(roc_obj)),
    CI_low = as.numeric(ci_obj[1]),
    CI_high = as.numeric(ci_obj[3]),
    Status = "Completed",
    stringsAsFactors = FALSE
  )
}

gse78220_response_tbl <- read_table_any(gse78220_state_file) %>%
  standardize_response_columns()

response_cols_keep <- intersect(
  c("Sample", "ResponseGroup", "BestResponse", "best_recist_response", "response", "ClinicalGroup"),
  colnames(gse78220_response_tbl)
)

response_tbl <- gse78220_response_tbl[, response_cols_keep, drop = FALSE]

if (!"ResponseGroup" %in% colnames(response_tbl)) {
  warning("ResponseGroup column not found in GSE78220 state table. AUC robustness will be skipped.")
}

auc_method_out <- list()

## Set the seed once before the ordered method/state loop.
## primary_constituent_zmean is evaluated first, and state_cols uses the same
## frozen state order as script 02; therefore the four primary-score bootstrap
## confidence intervals use the same random-number stream as script 02.
set.seed(auc_bootstrap_seed)

for (m in unique(scores782_all$ScoringMethod)) {
  tmp <- scores782_all %>%
    filter(ScoringMethod == m) %>%
    select(Sample, all_of(state_cols)) %>%
    left_join(response_tbl, by = "Sample")
  
  for (st in state_cols) {
    auc_method_out[[length(auc_method_out) + 1]] <- calc_auc(
      tmp,
      predictor_col = st,
      label = paste0(m, " | ", st)
    ) %>%
      mutate(
        ScoringMethod = m,
        FinalState = st
      )
  }
}

auc_method_tbl <- bind_rows(auc_method_out) %>%
  mutate(
    BootstrapSeed = auc_bootstrap_seed,
    BootstrapN = auc_bootstrap_n
  )

safe_write_csv(
  auc_method_tbl,
  file.path(out_table_dir, "GSE78220_AUC_by_scoring_method_and_state.csv")
)

p_auc <- auc_method_tbl %>%
  filter(Status == "Completed") %>%
  mutate(
    FinalState = factor(FinalState, levels = state_cols)
  ) %>%
  ggplot(aes(x = ScoringMethod, y = AUC)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.4) +
  facet_wrap(
    ~ FinalState,
    ncol = 2,
    labeller = labeller(FinalState = state_display_wrap)
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "GSE78220 non-response AUC across scoring methods",
    subtitle = "Exploratory robustness audit; not a predictive model",
    x = "Scoring method",
    y = "AUC for non-response"
  ) +
  theme_robustness(base_size = 11.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25))

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_GSE78220_AUC_by_scoring_method.pdf"),
  p_auc,
  width = 10,
  height = 7
)

safe_ggsave(
  file.path(out_fig_dir, "diagnostic_GSE78220_AUC_by_scoring_method.png"),
  p_auc,
  width = 10,
  height = 7
)

############################################################
## 14. Interpretation helper table
############################################################

interpretation_summary <- data.frame(
  Item = c(
    "Purpose",
    "What this analysis addresses",
    "What this analysis cannot address",
    "Signature construction audit",
    "Gene overlap audit",
    "Alternative scoring audit",
    "Leave-one-signature-out audit",
    "Gene subsampling audit",
    "GSE78220 AUC robustness"
  ),
  RecommendedInterpretation = c(
    "Assess whether predefined tumor–immune state scores are robust to gene-set construction, gene overlap, and scoring strategy.",
    "This analysis directly addresses concerns that state interpretation may be biased by signature selection, gene overlap, retained genes, or scoring method.",
    "This analysis does not convert predefined signatures into de novo-discovered states and does not provide functional validation.",
    "Report retained genes and missing genes per state and per dataset.",
    "Report Jaccard overlap among final-state gene pools to show whether the states are strongly overlapping or largely separable at the gene-set level.",
    "Report concordance between primary z-mean scoring and alternative rank-based or optional ssGSEA scoring.",
    "Report whether removing one constituent signature materially changes the final composite state score.",
    "Report whether random retention of 70%, 80%, or 90% of genes preserves sample-level state ranking.",
    "Report GSE78220 AUC across scoring methods only as descriptive robustness, not as biomarker validation."
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  interpretation_summary,
  file.path(out_table_dir, "interpretation_summary.csv")
)

############################################################
## 15. Output inventory and session info
############################################################

output_inventory <- data.frame(
  Type = c(
    "table",
    "table",
    "table",
    "table",
    "table",
    "table",
    "table",
    "table",
    "figure",
    "figure",
    "figure",
    "figure"
  ),
  File = c(
    file.path(out_table_dir, "input_file_inventory.csv"),
    file.path(out_table_dir, "final_state_raw_signature_mapping_used.csv"),
    file.path(out_table_dir, "gene_retention_summary_by_dataset_state.csv"),
    file.path(out_table_dir, "final_state_gene_pool_Jaccard_matrix.csv"),
    file.path(out_table_dir, "scoring_method_concordance_vs_primary.csv"),
    file.path(out_table_dir, "leave_one_signature_out_concordance.csv"),
    file.path(out_table_dir, "gene_subsampling_stability_summary.csv"),
    file.path(out_table_dir, "GSE78220_AUC_by_scoring_method_and_state.csv"),
    file.path(out_fig_dir, "diagnostic_final_state_gene_pool_Jaccard_heatmap.pdf"),
    file.path(out_fig_dir, "diagnostic_scoring_method_concordance_heatmap.pdf"),
    file.path(out_fig_dir, "diagnostic_leave_one_signature_out_robustness.pdf"),
    file.path(out_fig_dir, "diagnostic_gene_subsampling_robustness_boxplots.pdf")
  ),
  stringsAsFactors = FALSE
)

output_inventory$Exists <- file.exists(output_inventory$File)

safe_write_csv(
  output_inventory,
  file.path(out_table_dir, "output_inventory.csv")
)

sink(file.path(log_dir, "sessionInfo_03_state_scoring_robustness_analysis.txt"))
print(sessionInfo())
sink()

message("============================================================")
message("03_state_scoring_robustness_analysis.R completed successfully.")
message("Main table 1: ", file.path(out_table_dir, "gene_retention_summary_by_dataset_state.csv"))
message("Main table 2: ", file.path(out_table_dir, "final_state_gene_pool_overlap_long.csv"))
message("Main table 3: ", file.path(out_table_dir, "scoring_method_concordance_vs_primary.csv"))
message("Main table 4: ", file.path(out_table_dir, "leave_one_signature_out_concordance.csv"))
message("Main table 5: ", file.path(out_table_dir, "gene_subsampling_stability_summary.csv"))
message("============================================================")