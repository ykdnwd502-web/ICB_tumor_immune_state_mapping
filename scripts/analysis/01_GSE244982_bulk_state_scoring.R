############################################################
## 01_GSE244982_bulk_state_scoring.R
##
## Purpose:
##   Reproduce predefined tumor–immune state scoring in the GSE244982
##   bulk RNA-seq discovery cohort.
##
## Main steps:
##   1) Import the GSE244982 bulk expression matrix.
##   2) Read the frozen tumor–immune gene-set GMT.
##   3) Compute constituent gene-set scores by gene-wise z-score mean.
##   4) Build four predefined final tumor–immune state scores.
##   5) Export canonical state-score tables, QC summaries, diagnostic plots,
##      and session information for downstream analyses.
##
## Canonical final-state names:
##   Immune_defective_Cold
##   Myeloid_Treg_Immunosuppressive
##   Tumor_dedifferentiation_Stromal_remodeling
##   Melanocytic_Differentiation
##
## Reproducibility:
##   Run from the repository root, or set ICB_PROJECT_DIR to the repository
##   root before running. Package versions should be restored from renv.lock.
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Clean-run path contract: set ICB_PROJECT_DIR, or run from the repository root.
# Do not force a machine-specific working directory here.

############################################################
## 0. Project directory and output folders
############################################################

## Prefer explicit environment variable, otherwise use the current clean-room project path.
## Old C:/Users/ykdnw/Downloads/... paths are intentionally not used.
project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- getwd()
}

project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

message("Project directory: ", project_dir)

dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)

data_raw_dir <- file.path(project_dir, "data_raw")
data_processed_dir <- file.path(project_dir, "data_processed")
table_dir <- file.path(project_dir, "results/tables")
figure_dir <- file.path(project_dir, "results/figures")
bulk_table_dir <- file.path(table_dir, "bulk_discovery")
bulk_figure_dir <- file.path(figure_dir, "bulk_discovery")
log_dir <- file.path(project_dir, "logs")

dir.create(data_raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bulk_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bulk_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Load required packages
############################################################

required_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "scales",
  "patchwork"
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
  library(scales)
  library(patchwork)
})

theme_icb <- function(base_size = 12.5, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2.5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 1.5),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1.5),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6)
    )
}

############################################################
## 2. Utility functions
############################################################

save_session_info <- function(script_name = "01_GSE244982_bulk_state_scoring") {
  log_file <- file.path(
    log_dir,
    paste0("sessionInfo_", script_name, ".txt")
  )
  
  sink(log_file)
  print(sessionInfo())
  sink()
  
  message("sessionInfo saved to: ", log_file)
}

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
  message("Saved table: ", file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = file)
  message("Saved RDS: ", file)
}

safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  message("Saved figure: ", file)
}

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s+", "", x)
  x <- toupper(x)
  x
}

normalize_name <- function(x) {
  x <- as.character(x)
  x <- gsub("\\ufeff", "", x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

find_first_existing <- function(paths) {
  paths <- unique(paths)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

zvec <- function(x) {
  x <- as.numeric(x)
  out <- as.numeric(scale(x))
  out[is.na(out)] <- 0
  out
}

############################################################
## 3. Locate GSE244982 bulk expression matrix
############################################################

bulk_candidates <- c(
  file.path(project_dir, "data_raw/GSE244982_ProcessedData_bulkRNAseq.txt.gz"),
  file.path(project_dir, "data_raw/GSE244982_ProcessedData_bulkRNAseq.txt"),
  file.path(project_dir, "data_raw/GSE244982_ProcessedData_bulkRNAseq.csv"),
  file.path(project_dir, "data_raw/GSE244982_ProcessedData_bulkRNAseq.tsv"),
  file.path(project_dir, "data_raw/GSE244982/GSE244982_ProcessedData_bulkRNAseq.txt.gz"),
  file.path(project_dir, "data_raw/GSE244982/GSE244982_ProcessedData_bulkRNAseq.txt"),
  file.path(project_dir, "data_raw/GSE244982/GSE244982_ProcessedData_bulkRNAseq.csv"),
  file.path(project_dir, "data_raw/GSE244982/GSE244982_ProcessedData_bulkRNAseq.tsv"),
  file.path(project_dir, "data_processed/GSE244982_ProcessedData_bulkRNAseq.txt.gz"),
  file.path(project_dir, "data_processed/GSE244982_ProcessedData_bulkRNAseq.txt"),
  file.path(project_dir, "data_processed/GSE244982_bulk_expression_matrix.csv")
)

bulk_rds_candidates <- c(
  file.path(project_dir, "data_processed/GSE244982_bulk_matrix.rds"),
  file.path(project_dir, "data_processed/GSE244982_bulk_expr.rds"),
  file.path(project_dir, "data_processed/GSE244982_bulk_expression_matrix.rds")
)

bulk_file <- find_first_existing(bulk_candidates)
bulk_rds_file <- find_first_existing(bulk_rds_candidates)

message("Bulk text/csv candidate selected: ", ifelse(is.na(bulk_file), "none", bulk_file))
message("Bulk RDS candidate selected: ", ifelse(is.na(bulk_rds_file), "none", bulk_rds_file))

############################################################
## 4. Import or load bulk matrix
############################################################

expr_mat <- NULL

if (!is.na(bulk_file)) {
  message("Reading bulk expression file: ", bulk_file)
  
  expr_dt <- data.table::fread(
    bulk_file,
    data.table = FALSE,
    check.names = FALSE
  )
  
  if (nrow(expr_dt) == 0 || ncol(expr_dt) < 2) {
    stop("Bulk expression file was read but appears empty or has fewer than 2 columns: ", bulk_file)
  }
  
  col_names <- colnames(expr_dt)
  gene_col_candidates <- c(
    "Gene", "gene", "GeneSymbol", "gene_symbol", "Symbol", "symbol",
    "genes", "Genes", "GeneName", "gene_name", "ID", "id"
  )
  
  gene_col <- gene_col_candidates[gene_col_candidates %in% col_names][1]
  
  if (is.na(gene_col)) {
    gene_col <- col_names[1]
    message("No explicit gene column detected. Using first column as gene column: ", gene_col)
  } else {
    message("Detected gene column: ", gene_col)
  }
  
  gene_symbols <- clean_gene_symbol(expr_dt[[gene_col]])
  sample_cols <- setdiff(col_names, gene_col)
  
  numeric_test <- sapply(expr_dt[, sample_cols, drop = FALSE], function(x) {
    suppressWarnings(mean(!is.na(as.numeric(x)))) > 0.5
  })
  
  sample_cols_numeric <- sample_cols[numeric_test]
  
  if (length(sample_cols_numeric) == 0) {
    stop("No numeric expression columns detected in bulk expression file.")
  }
  
  expr_numeric <- expr_dt[, sample_cols_numeric, drop = FALSE]
  expr_numeric[] <- lapply(expr_numeric, function(x) suppressWarnings(as.numeric(x)))
  
  expr_mat0 <- as.matrix(expr_numeric)
  rownames(expr_mat0) <- gene_symbols
  
  keep_gene <- !is.na(rownames(expr_mat0)) & rownames(expr_mat0) != ""
  expr_mat0 <- expr_mat0[keep_gene, , drop = FALSE]
  
  keep_nonempty <- rowSums(!is.na(expr_mat0)) > 0
  expr_mat0 <- expr_mat0[keep_nonempty, , drop = FALSE]
  
  if (any(duplicated(rownames(expr_mat0)))) {
    message("Duplicated gene symbols detected. Aggregating duplicated genes by mean.")
    
    expr_df <- as.data.frame(expr_mat0)
    expr_df$GeneSymbol <- rownames(expr_mat0)
    
    expr_df_agg <- expr_df %>%
      group_by(GeneSymbol) %>%
      summarise(
        across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )
    
    expr_mat <- as.matrix(expr_df_agg[, setdiff(colnames(expr_df_agg), "GeneSymbol"), drop = FALSE])
    rownames(expr_mat) <- expr_df_agg$GeneSymbol
  } else {
    expr_mat <- expr_mat0
  }
  
} else if (!is.na(bulk_rds_file)) {
  message("No text/csv bulk matrix found. Loading existing RDS: ", bulk_rds_file)
  expr_mat <- readRDS(bulk_rds_file)
  
  if (is.data.frame(expr_mat)) {
    expr_mat <- as.matrix(expr_mat)
  }
} else {
  stop(
    "No GSE244982 bulk expression file found. Expected one of:\n",
    paste(unique(c(bulk_candidates, bulk_rds_candidates)), collapse = "\n")
  )
}

expr_mat <- as.matrix(expr_mat)
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- clean_gene_symbol(rownames(expr_mat))
expr_mat <- expr_mat[rowSums(!is.na(expr_mat)) > 0, colSums(!is.na(expr_mat)) > 0, drop = FALSE]

message("Bulk matrix dimensions: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples")

safe_save_rds(
  expr_mat,
  file.path(data_processed_dir, "GSE244982_bulk_matrix.rds")
)

safe_write_csv(
  data.frame(GeneSymbol = rownames(expr_mat), expr_mat, check.names = FALSE),
  file.path(bulk_table_dir, "GSE244982_bulk_expression_matrix.csv"),
  row.names = FALSE
)

qc_df <- data.frame(
  Sample = colnames(expr_mat),
  n_detected_genes = colSums(!is.na(expr_mat) & expr_mat != 0),
  total_expression = colSums(expr_mat, na.rm = TRUE),
  mean_expression = colMeans(expr_mat, na.rm = TRUE),
  median_expression = apply(expr_mat, 2, median, na.rm = TRUE),
  stringsAsFactors = FALSE
)

safe_write_csv(
  qc_df,
  file.path(bulk_table_dir, "GSE244982_bulk_sample_QC.csv"),
  row.names = FALSE
)

p_detected <- ggplot(qc_df, aes(x = reorder(Sample, n_detected_genes), y = n_detected_genes)) +
  geom_col(width = 0.75, fill = "#4C72B0") +
  coord_flip() +
  labs(
    title = "Detected genes per bulk sample",
    x = NULL,
    y = "Detected genes"
  ) +
  theme_icb()

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_sample_detected_genes.pdf"),
  p_detected,
  width = 7.5,
  height = max(4, 0.25 * nrow(qc_df) + 2)
)

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_sample_detected_genes.png"),
  p_detected,
  width = 7.5,
  height = max(4, 0.25 * nrow(qc_df) + 2),
  dpi = 300
)

############################################################
## 5. Locate and read GMT gene sets
############################################################

gmt_candidates <- c(
  Sys.getenv("ICB_CLEAN_GMT"),
  file.path(
    project_dir,
    "results/tables/ICBcomb_final_input_gene_sets_CLEAN/02_final_ICBcomb_gene_sets/final_ICBcomb_gene_sets_CLEAN.gmt"
  )
)

gmt_candidates <- unique(gmt_candidates[nzchar(gmt_candidates)])
gmt_file <- find_first_existing(gmt_candidates)

message("CLEAN GMT candidate selected: ", ifelse(is.na(gmt_file), "none", gmt_file))

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

if (is.na(gmt_file)) {
  stop(
    "No CLEAN gene-set GMT found. Expected one of:\n",
    paste(unique(gmt_candidates), collapse = "\n")
  )
}

gene_sets <- read_gmt(gmt_file)
message("Number of gene sets imported from GMT: ", length(gene_sets))

gene_set_summary <- data.frame(
  GeneSet = names(gene_sets),
  n_genes = sapply(gene_sets, length),
  stringsAsFactors = FALSE
)

safe_write_csv(
  gene_set_summary,
  file.path(bulk_table_dir, "GSE244982_imported_gene_set_summary.csv"),
  row.names = FALSE
)

############################################################
## 6. Bulk gene-set scoring by simple z-mean method
############################################################

score_gene_set_zmean <- function(expr_mat, gene_sets) {
  expr_genes <- rownames(expr_mat)
  
  expr_z <- t(scale(t(expr_mat)))
  expr_z[is.na(expr_z)] <- 0
  
  score_list <- list()
  presence_list <- list()
  
  for (gs in names(gene_sets)) {
    genes <- clean_gene_symbol(gene_sets[[gs]])
    present <- intersect(genes, expr_genes)
    
    message(gs, ": ", length(present), "/", length(genes), " genes present")
    
    presence_list[[gs]] <- data.frame(
      GeneSet = gs,
      Gene = genes,
      Present = genes %in% expr_genes,
      stringsAsFactors = FALSE
    )
    
    if (length(present) == 0) {
      score <- rep(NA_real_, ncol(expr_mat))
      names(score) <- colnames(expr_mat)
    } else {
      score <- colMeans(expr_z[present, , drop = FALSE], na.rm = TRUE)
    }
    
    score_list[[gs]] <- score
  }
  
  score_mat <- do.call(rbind, score_list)
  rownames(score_mat) <- names(score_list)
  
  list(
    score_matrix = score_mat,
    raw_scores = as.data.frame(t(score_mat), check.names = FALSE),
    gene_presence = dplyr::bind_rows(presence_list)
  )
}

bulk_score_result <- score_gene_set_zmean(expr_mat, gene_sets)

raw_scores <- bulk_score_result$raw_scores
raw_scores$Sample <- rownames(raw_scores)
raw_scores <- raw_scores[, c("Sample", setdiff(colnames(raw_scores), "Sample")), drop = FALSE]

safe_write_csv(
  raw_scores,
  file.path(bulk_table_dir, "GSE244982_bulk_raw_signature_scores.csv"),
  row.names = FALSE
)

safe_write_csv(
  raw_scores,
  file.path(table_dir, "GSE244982_bulk_raw_signature_scores.csv"),
  row.names = FALSE
)

safe_save_rds(
  raw_scores,
  file.path(data_processed_dir, "GSE244982_bulk_raw_signature_scores.rds")
)

safe_write_csv(
  bulk_score_result$gene_presence,
  file.path(bulk_table_dir, "GSE244982_bulk_gene_presence_all_gene_sets.csv"),
  row.names = FALSE
)

for (gs in unique(bulk_score_result$gene_presence$GeneSet)) {
  gs_df <- bulk_score_result$gene_presence %>% filter(GeneSet == gs)
  safe_write_csv(
    gs_df,
    file.path(
      bulk_table_dir,
      paste0("Bulk_gene_presence_", make.names(gs), ".csv")
    ),
    row.names = FALSE
  )
}

############################################################
## 7. Build final mechanism-defined state scores
############################################################

get_signature_cols <- function(raw_scores, target_names, state_name, min_required = 1) {
  available <- setdiff(colnames(raw_scores), "Sample")
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
  message(state_name, " matched signature columns: ", paste(matched, collapse = ", "))
  
  if (length(matched) < min_required) {
    stop(
      "Insufficient matched signature columns for ", state_name,
      ". Required at least ", min_required,
      ", found ", length(matched),
      ". Available columns:\n",
      paste(available, collapse = "\n")
    )
  }
  
  matched
}

immune_cols <- get_signature_cols(
  raw_scores,
  c(
    "Immune_defective_Cold_RESTORE_hybrid",
    "Immune_defective_Cold_RESTORE_data_driven",
    "Immune_defective_Cold_RESTORE_curated",
    "Immune_cold_RESTORE_hybrid",
    "Immune_cold_RESTORE_data_driven",
    "Immune_cold_RESTORE_curated"
  ),
  "Immune_defective_cold",
  min_required = 1
)

myeloid_cols <- get_signature_cols(
  raw_scores,
  c(
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated",
    "Myeloid_Treg_SUPPRESS_hybrid",
    "Myeloid_Treg_SUPPRESS_data_driven",
    "Myeloid_Treg_SUPPRESS_curated"
  ),
  "Myeloid_Treg_immunosuppressive",
  min_required = 1
)

tumor_cols <- get_signature_cols(
  raw_scores,
  c(
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated",
    "Tumor_dedifferentiation_SUPPRESS_hybrid",
    "Tumor_dedifferentiation_SUPPRESS_data_driven",
    "Tumor_dedifferentiation_SUPPRESS_curated"
  ),
  "Tumor_dedifferentiated",
  min_required = 1
)

melanocytic_cols <- get_signature_cols(
  raw_scores,
  c(
    "Melanocytic_Differentiation_REFERENCE_hybrid",
    "Melanocytic_Differentiation_REFERENCE_data_driven",
    "Melanocytic_Differentiation_REFERENCE_curated",
    "Melanocytic_differentiation_REFERENCE_hybrid",
    "Melanocytic_differentiation_REFERENCE_data_driven",
    "Melanocytic_differentiation_REFERENCE_curated"
  ),
  "Melanocytic_differentiated",
  min_required = 1
)

state_mapping_table <- data.frame(
  FinalState = c(
    rep("Immune_defective_cold", length(immune_cols)),
    rep("Myeloid_Treg_immunosuppressive", length(myeloid_cols)),
    rep("Tumor_dedifferentiated", length(tumor_cols)),
    rep("Melanocytic_differentiated", length(melanocytic_cols))
  ),
  RawSignatureColumn = c(
    immune_cols,
    myeloid_cols,
    tumor_cols,
    melanocytic_cols
  ),
  DirectionForFinalState = c(
    rep("inverse", length(immune_cols)),
    rep("positive", length(myeloid_cols)),
    rep("positive", length(tumor_cols)),
    rep("positive", length(melanocytic_cols))
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  state_mapping_table,
  file.path(bulk_table_dir, "GSE244982_final_state_signature_mapping.csv"),
  row.names = FALSE
)

final_state <- data.frame(
  Immune_defective_cold = zvec(-rowMeans(raw_scores[, immune_cols, drop = FALSE], na.rm = TRUE)),
  Myeloid_Treg_immunosuppressive = zvec(rowMeans(raw_scores[, myeloid_cols, drop = FALSE], na.rm = TRUE)),
  Tumor_dedifferentiated = zvec(rowMeans(raw_scores[, tumor_cols, drop = FALSE], na.rm = TRUE)),
  Melanocytic_differentiated = zvec(rowMeans(raw_scores[, melanocytic_cols, drop = FALSE], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

rownames(final_state) <- raw_scores$Sample

numeric_state_cols <- c(
  "Immune_defective_cold",
  "Myeloid_Treg_immunosuppressive",
  "Tumor_dedifferentiated",
  "Melanocytic_differentiated"
)

na_check <- sapply(final_state[, numeric_state_cols, drop = FALSE], function(x) sum(is.na(x)))
print(na_check)

safe_write_csv(
  data.frame(State = names(na_check), NA_count = as.integer(na_check)),
  file.path(bulk_table_dir, "GSE244982_final_state_NA_check.csv"),
  row.names = FALSE
)

dominant_state <- apply(
  final_state[, numeric_state_cols, drop = FALSE],
  1,
  function(x) numeric_state_cols[which.max(x)]
)

final_state$Dominant_final_state <- dominant_state

final_state_export <- data.frame(
  Sample = rownames(final_state),
  final_state,
  check.names = FALSE,
  row.names = NULL
)

final_state_export_CLEAN_NAMES <- final_state_export
final_state_export_CLEAN_NAMES <- final_state_export_CLEAN_NAMES[, c("Sample", numeric_state_cols, "Dominant_final_state"), drop = FALSE]
colnames(final_state_export_CLEAN_NAMES) <- c(
  "Sample",
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation",
  "Dominant_final_state"
)
final_state_export_CLEAN_NAMES$Dominant_final_state_CLEAN <- dplyr::recode(
  final_state_export_CLEAN_NAMES$Dominant_final_state,
  "Immune_defective_cold" = "Immune_defective_Cold",
  "Myeloid_Treg_immunosuppressive" = "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiated" = "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_differentiated" = "Melanocytic_Differentiation",
  .default = final_state_export_CLEAN_NAMES$Dominant_final_state
)

final_state_public <- final_state_export_CLEAN_NAMES %>%
  dplyr::select(
    Sample,
    Immune_defective_Cold,
    Myeloid_Treg_Immunosuppressive,
    Tumor_dedifferentiation_Stromal_remodeling,
    Melanocytic_Differentiation,
    Dominant_final_state_CLEAN
  ) %>%
  dplyr::rename(Dominant_state = Dominant_final_state_CLEAN)

safe_write_csv(
  final_state_public,
  file.path(bulk_table_dir, "GSE244982_final_tumor_immune_state_scores.csv"),
  row.names = FALSE
)

safe_save_rds(
  final_state_public,
  file.path(data_processed_dir, "GSE244982_final_tumor_immune_state_scores.rds")
)

dominant_count <- final_state_export %>%
  count(Dominant_final_state, name = "n_samples") %>%
  arrange(desc(n_samples))

dominant_count_public <- final_state_public %>%
  dplyr::count(Dominant_state, name = "n_samples") %>%
  dplyr::arrange(desc(n_samples))

safe_write_csv(
  dominant_count_public,
  file.path(bulk_table_dir, "GSE244982_dominant_state_counts.csv"),
  row.names = FALSE
)

cor_mat <- cor(
  final_state_public[, c(
    "Immune_defective_Cold",
    "Myeloid_Treg_Immunosuppressive",
    "Tumor_dedifferentiation_Stromal_remodeling",
    "Melanocytic_Differentiation"
  ), drop = FALSE],
  method = "spearman",
  use = "pairwise.complete.obs"
)

safe_write_csv(
  as.data.frame(cor_mat),
  file.path(bulk_table_dir, "GSE244982_final_state_correlation_matrix.csv"),
  row.names = TRUE
)

safe_write_csv(
  as.data.frame(cor_mat),
  file.path(table_dir, "GSE244982_final_state_correlation_matrix.csv"),
  row.names = TRUE
)

############################################################
## 8. Raw signature heatmap
############################################################

raw_long <- raw_scores %>%
  pivot_longer(
    cols = -Sample,
    names_to = "GeneSet",
    values_to = "Score"
  )

p_raw_heat <- ggplot(raw_long, aes(x = Sample, y = GeneSet, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    name = "z-mean\nscore"
  ) +
  labs(
    title = "Bulk gene-set scores in GSE244982",
    x = NULL,
    y = NULL
  ) +
  theme_icb(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5),
    panel.grid = element_blank()
  )

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_gene_set_scores_simple_zmean.pdf"),
  p_raw_heat,
  width = max(8, 0.35 * length(unique(raw_long$Sample)) + 4),
  height = max(5, 0.25 * length(unique(raw_long$GeneSet)) + 2)
)

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_gene_set_scores_simple_zmean.png"),
  p_raw_heat,
  width = max(8, 0.35 * length(unique(raw_long$Sample)) + 4),
  height = max(5, 0.25 * length(unique(raw_long$GeneSet)) + 2),
  dpi = 300
)

############################################################
## 9. Corrected final-state heatmap and PCA
############################################################

state_label_map <- c(
  "Immune_defective_cold" = "immune-defective/\ncold",
  "Myeloid_Treg_immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiated" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_differentiated" = "melanocytic\ndifferentiation"
)

# 【修改 1】将第三个标签拆为两行，避免图例过长被向右挤压
state_label_one_line <- c(
  "Immune_defective_cold" = "immune-defective/cold",
  "Myeloid_Treg_immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiated" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_differentiated" = "melanocytic differentiation"
)

state_colors <- c(
  "Immune_defective_cold" = "#4DBBD5",
  "Myeloid_Treg_immunosuppressive" = "#00A087",
  "Tumor_dedifferentiated" = "#E64B35",
  "Melanocytic_differentiated" = "#3C5488"
)

state_order <- numeric_state_cols

sample_order <- final_state_export %>%
  arrange(
    factor(Dominant_final_state, levels = state_order),
    desc(Immune_defective_cold)
  ) %>%
  pull(Sample)

state_long <- final_state_export %>%
  select(Sample, all_of(numeric_state_cols)) %>%
  pivot_longer(
    cols = all_of(numeric_state_cols),
    names_to = "State",
    values_to = "Score"
  ) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    StateLabel = factor(
      state_label_map[State],
      levels = rev(unname(state_label_map))
    )
  )

p_state_heat <- ggplot(
  state_long,
  aes(x = Sample, y = StateLabel, fill = Score)
) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    name = "State\nscore"
  ) +
  labs(
    title = "Predefined bulk tumor–immune state scores",
    x = NULL,
    y = NULL
  ) +
  theme_icb(base_size = 12.5) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10.75),
    axis.text.y = element_text(size = 11.5),
    panel.grid = element_blank()
  )

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_resistance_state_scores_heatmap.pdf"),
  p_state_heat,
  width = 13,
  height = 4.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_resistance_state_scores_heatmap.png"),
  p_state_heat,
  width = 13,
  height = 4.8,
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, "GSE244982_state_score_heatmap_diagnostic.pdf"),
  p_state_heat,
  width = 13,
  height = 4.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "GSE244982_state_score_heatmap_diagnostic.png"),
  p_state_heat,
  width = 13,
  height = 4.8,
  dpi = 300
)

p_state_heat_annotated <- p_state_heat +
  geom_text(aes(label = sprintf("%.2f", Score)), size = 2.5, family = "sans")

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_resistance_state_scores_heatmap_annotated_audit.pdf"),
  p_state_heat_annotated,
  width = 13,
  height = 4.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_resistance_state_scores_heatmap_annotated_audit.png"),
  p_state_heat_annotated,
  width = 13,
  height = 4.8,
  dpi = 300
)

dominant_count_plot <- dominant_count %>%
  mutate(
    Dominant_final_state = factor(Dominant_final_state, levels = state_order),
    StateLabel = factor(
      state_label_map[as.character(Dominant_final_state)],
      levels = rev(state_label_map[state_order])
    )
  )

# 【修改 2】添加适当右边距，防止 panel c 标题右侧被裁剪
p_dom <- ggplot(
  dominant_count_plot,
  aes(
    x = StateLabel,
    y = n_samples,
    fill = Dominant_final_state
  )
) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = state_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Distribution of dominant tumor–immune states",
    x = NULL,
    y = "Number of samples"
  ) +
  theme_icb(base_size = 12.5) +
  theme(
    axis.text.y = element_text(size = 11.5),
    plot.margin = margin(t = 5, r = 18, b = 5, l = 5, unit = "pt")
  )

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_dominant_state_counts.pdf"),
  p_dom,
  width = 6.5,
  height = 4.5
)

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_dominant_state_counts.png"),
  p_dom,
  width = 6.5,
  height = 4.5,
  dpi = 300
)

pca_mat <- as.matrix(final_state[, numeric_state_cols, drop = FALSE])
pca_mat[is.na(pca_mat)] <- 0

pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  Sample = rownames(final_state),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Dominant_final_state = final_state$Dominant_final_state,
  stringsAsFactors = FALSE
)

safe_write_csv(
  pca_df,
  file.path(bulk_table_dir, "GSE244982_final_state_PCA_coordinates.csv"),
  row.names = FALSE
)

safe_write_csv(
  pca_df,
  file.path(table_dir, "GSE244982_PCA_coordinates.csv"),
  row.names = FALSE
)

var_exp <- summary(pca)$importance[2, 1:2] * 100

pca_df <- pca_df %>%
  mutate(
    Dominant_final_state = factor(Dominant_final_state, levels = state_order),
    DominantStateLabel = factor(
      state_label_one_line[as.character(Dominant_final_state)],
      levels = state_label_one_line[state_order]
    )
  )

p_pca_clean <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Dominant_final_state)
) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(
    values = state_colors,
    breaks = state_order,
    labels = state_label_one_line[state_order],
    name = "Dominant state"
  ) +
  labs(
    title = "PCA of bulk tumor–immune state scores",
    x = paste0("PC1 (", sprintf("%.1f", var_exp[1]), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", var_exp[2]), "%)")
  ) +
  theme_icb(base_size = 12.5) +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(0.75, "cm"),
    legend.text = element_text(lineheight = 0.85)
  )

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_resistance_state_PCA.pdf"),
  p_pca_clean,
  width = 7,
  height = 5.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_resistance_state_PCA.png"),
  p_pca_clean,
  width = 7,
  height = 5.8,
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, "GSE244982_state_heatmap.pdf"),
  p_state_heat,
  width = 13,
  height = 4.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "GSE244982_state_heatmap.png"),
  p_state_heat,
  width = 13,
  height = 4.8,
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, "GSE244982_state_PCA.pdf"),
  p_pca_clean,
  width = 7,
  height = 5.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "GSE244982_state_PCA.png"),
  p_pca_clean,
  width = 7,
  height = 5.8,
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, "GSE244982_dominant_state_counts.pdf"),
  p_dom,
  width = 6.5,
  height = 4.5
)

safe_ggsave(
  file.path(bulk_figure_dir, "GSE244982_dominant_state_counts.png"),
  p_dom,
  width = 6.5,
  height = 4.5,
  dpi = 300
)

# 【修改 3】拼图宽度权重与边距协调，让 panel C 适当往左移并完整展示
bulk_state_overview <- p_state_heat / (p_pca_clean | p_dom) +
  plot_layout(heights = c(0.92, 1.15), widths = c(1.02, 0.98)) +
  plot_annotation(
    tag_levels = "A",
    title = "Discovery of mechanism-defined bulk tumor–immune states in GSE244982",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16.5, family = "sans"),
      plot.tag = element_text(face = "bold", size = 16, family = "sans")
    )
  )

# 【修改 4】修改文件名为 "Figure 2. Predefined tumor–immune state scores in GSE244982"，并新增 300 dpi JPG 保存
fig2_name <- "Figure 2. Predefined tumor–immune state scores in GSE244982"

safe_ggsave(
  file.path(bulk_figure_dir, paste0(fig2_name, ".pdf")),
  bulk_state_overview,
  width = 13.2,
  height = 9.2
)

safe_ggsave(
  file.path(bulk_figure_dir, paste0(fig2_name, ".png")),
  bulk_state_overview,
  width = 13.2,
  height = 9.2,
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, paste0(fig2_name, ".jpg")),
  bulk_state_overview,
  width = 13.2,
  height = 9.2,
  dpi = 300
)

# 审计用带样本名散点图
p_pca_label <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, label = Sample)
) +
  geom_point(size = 3, color = "#4C72B0") +
  geom_text(vjust = -0.8, size = 3, family = "sans") +
  labs(
    title = "PCA of bulk tumor–immune state scores",
    x = paste0("PC1 (", sprintf("%.1f", var_exp[1]), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", var_exp[2]), "%)")
  ) +
  theme_icb(base_size = 12.5)

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_resistance_state_PCA_labeled_audit.pdf"),
  p_pca_label,
  width = 7,
  height = 5.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "Bulk_resistance_state_PCA_labeled_audit.png"),
  p_pca_label,
  width = 7,
  height = 5.8,
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, "diagnostic_GSE244982_PCA_labeled.pdf"),
  p_pca_label,
  width = 7,
  height = 5.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "diagnostic_GSE244982_PCA_labeled.png"),
  p_pca_label,
  width = 7,
  height = 5.8,
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, "diagnostic_GSE244982_state_heatmap_annotated.pdf"),
  p_state_heat_annotated,
  width = 13,
  height = 4.8
)

safe_ggsave(
  file.path(bulk_figure_dir, "diagnostic_GSE244982_state_heatmap_annotated.png"),
  p_state_heat_annotated,
  width = 13,
  height = 4.8,
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, "diagnostic_GSE244982_gene_set_scores_zmean.pdf"),
  p_raw_heat,
  width = max(8, 0.35 * length(unique(raw_long$Sample)) + 4),
  height = max(5, 0.25 * length(unique(raw_long$GeneSet)) + 2)
)

safe_ggsave(
  file.path(bulk_figure_dir, "diagnostic_GSE244982_gene_set_scores_zmean.png"),
  p_raw_heat,
  width = max(8, 0.35 * length(unique(raw_long$Sample)) + 4),
  height = max(5, 0.25 * length(unique(raw_long$GeneSet)) + 2),
  dpi = 300
)

safe_ggsave(
  file.path(bulk_figure_dir, "diagnostic_GSE244982_detected_genes.pdf"),
  p_detected,
  width = 7.5,
  height = max(4, 0.25 * nrow(qc_df) + 2)
)

safe_ggsave(
  file.path(bulk_figure_dir, "diagnostic_GSE244982_detected_genes.png"),
  p_detected,
  width = 7.5,
  height = max(4, 0.25 * nrow(qc_df) + 2),
  dpi = 300
)

############################################################
## 10. Key output inventory and session information
############################################################

key_outputs <- data.frame(
  Type = c("table", "table", "table", "table", "rds"),
  File = c(
    file.path(bulk_table_dir, "GSE244982_final_tumor_immune_state_scores.csv"),
    file.path(bulk_table_dir, "GSE244982_bulk_sample_QC.csv"),
    file.path(bulk_table_dir, "GSE244982_bulk_raw_signature_scores.csv"),
    file.path(bulk_table_dir, "GSE244982_final_state_correlation_matrix.csv"),
    file.path(data_processed_dir, "GSE244982_bulk_matrix.rds")
  ),
  stringsAsFactors = FALSE
)
key_outputs$Exists <- file.exists(key_outputs$File)

safe_write_csv(
  key_outputs,
  file.path(bulk_table_dir, "GSE244982_bulk_state_scoring_output_inventory.csv"),
  row.names = FALSE
)

save_session_info("01_GSE244982_bulk_state_scoring")

message("01_GSE244982_bulk_state_scoring.R finished successfully.")

