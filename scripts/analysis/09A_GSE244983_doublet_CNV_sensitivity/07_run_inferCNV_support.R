############################################################
## 09A-07. inferCNV support analysis on pANN-retained singlets
##
## Exploratory CNV-like support only; not a malignant-cell classifier.
## Reference cells are defined from the current major-cell annotation:
##   T/NK cells, T/NK/Treg-like cells, B/Plasma cells, Myeloid cells.
## CAF/stromal, endothelial, malignant, and cycling malignant cells are
## not used as primary inferCNV references.
############################################################

options(stringsAsFactors = FALSE)
SCRIPT_ID <- "09A_07_run_inferCNV_support"
set.seed(12345)

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

required_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "readr", "dplyr", "tibble",
  "infercnv", "matrixStats"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required inferCNV package(s): ", paste(missing_pkgs, collapse = ", "),
    ". Restore the frozen environment before the full clean run.", call. = FALSE
  )
}

out_table_dir <- file.path(project_dir, "results", "tables", "GSE244983", "doublet_CNV_sensitivity")
out_obj_dir <- file.path(project_dir, "results", "objects", "GSE244983", "doublet_CNV_sensitivity")
infercnv_out_dir <- file.path(project_dir, "results", "inferCNV", "GSE244983_doublet_CNV_sensitivity")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(infercnv_out_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(out_obj_dir, "GSE244983_singlets_after_pANN_doublet_filter.rds")
gene_order_file <- file.path(out_table_dir, "D06_inferCNV_gene_order_hg38_noheader.txt")
result_rds <- file.path(out_obj_dir, "GSE244983_inferCNV_result.rds")
write_csv <- function(x, name) readr::write_csv(x, file.path(out_table_dir, name), na = "")

reference_major_types <- c(
  "T/NK cells",
  "T/NK/Treg-like cells",
  "B/Plasma cells",
  "Myeloid cells"
)
max_reference_cells_per_type <- 750L
max_nonreference_cells_per_cluster <- 600L
infercnv_cutoff <- 0.1
infercnv_denoise <- TRUE
infercnv_HMM <- FALSE
infercnv_cluster_by_groups <- TRUE
thread_env <- suppressWarnings(as.integer(Sys.getenv("ICB_INFERCNV_THREADS", unset = NA_character_)))
if (is.na(thread_env) || thread_env < 1L) {
  detected <- parallel::detectCores(logical = TRUE)
  if (is.na(detected)) detected <- 2L
  infercnv_num_threads <- max(1L, min(8L, detected - 1L))
} else {
  infercnv_num_threads <- thread_env
}

get_counts_matrix <- function(obj, assay = "RNA") {
  if (!assay %in% names(obj@assays)) assay <- Seurat::DefaultAssay(obj)
  mat <- tryCatch(
    SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(mat)) stop("Could not extract RNA counts matrix.", call. = FALSE)
  mat
}

cnv_signal_chunked <- function(expr, ref_median, chunk_size = 500L) {
  n <- ncol(expr)
  out <- numeric(n)
  names(out) <- colnames(expr)
  starts <- seq.int(1L, n, by = chunk_size)
  for (st in starts) {
    en <- min(n, st + chunk_size - 1L)
    m <- as.matrix(expr[, st:en, drop = FALSE])
    m <- sweep(m, 1L, ref_median, FUN = "-")
    out[st:en] <- colMeans(abs(m), na.rm = TRUE)
  }
  out
}

if (!file.exists(input_rds)) stop("09A-03 singlet object missing.", call. = FALSE)
if (!file.exists(gene_order_file)) stop("09A-06 gene-order file missing.", call. = FALSE)
obj <- readRDS(input_rds)
if (!all(c("Sample", "MajorCellType", "SeuratCluster") %in% colnames(obj@meta.data))) {
  stop("Required current annotations missing from singlet object.", call. = FALSE)
}

major_present <- unique(as.character(obj$MajorCellType))
missing_reference_types <- setdiff(reference_major_types, major_present)
if (length(missing_reference_types) > 0L) {
  stop("Required immune reference types missing: ", paste(missing_reference_types, collapse = ", "), call. = FALSE)
}

cell_df <- data.frame(
  cell_id = colnames(obj),
  Sample = as.character(obj$Sample),
  MajorCellType = as.character(obj$MajorCellType),
  SeuratCluster = as.character(obj$SeuratCluster),
  stringsAsFactors = FALSE
)
cell_df$is_reference <- cell_df$MajorCellType %in% reference_major_types
cell_df$infercnv_group <- ifelse(cell_df$is_reference, "Reference", paste0("Cluster_", cell_df$SeuratCluster))

## Balanced reference sampling by major immune type.
ref_df <- cell_df[cell_df$is_reference, , drop = FALSE]
ref_parts <- lapply(reference_major_types, function(ct) {
  x <- ref_df[ref_df$MajorCellType == ct, , drop = FALSE]
  n_take <- min(max_reference_cells_per_type, nrow(x))
  if (n_take < nrow(x)) x <- x[sample.int(nrow(x), n_take), , drop = FALSE]
  x
})
ref_used <- do.call(rbind, ref_parts)

## Non-reference sampling is stratified by the current SeuratCluster.
nonref_df <- cell_df[!cell_df$is_reference, , drop = FALSE]
nonref_split <- split(nonref_df, nonref_df$SeuratCluster)
nonref_parts <- lapply(nonref_split, function(x) {
  n_take <- min(max_nonreference_cells_per_cluster, nrow(x))
  if (n_take < nrow(x)) x <- x[sample.int(nrow(x), n_take), , drop = FALSE]
  x
})
nonref_used <- do.call(rbind, nonref_parts)
infer_cells_df <- rbind(ref_used, nonref_used)
infer_cells_df <- infer_cells_df[!duplicated(infer_cells_df$cell_id), , drop = FALSE]

if (sum(infer_cells_df$is_reference) < 100L) stop("Too few reference cells after downsampling.", call. = FALSE)

reference_audit <- aggregate(
  cell_id ~ MajorCellType + is_reference,
  data = cell_df,
  FUN = length
)
names(reference_audit)[3] <- "N_before"
used_ref_audit <- aggregate(
  cell_id ~ MajorCellType + is_reference,
  data = infer_cells_df,
  FUN = length
)
names(used_ref_audit)[3] <- "N_used"
reference_audit <- merge(reference_audit, used_ref_audit, by = c("MajorCellType", "is_reference"), all.x = TRUE)
reference_audit$N_used[is.na(reference_audit$N_used)] <- 0L
write_csv(reference_audit, "D07_inferCNV_reference_and_downsampling_audit.csv")

annotation_file <- file.path(out_table_dir, "D07_inferCNV_cell_annotations.txt")
readr::write_tsv(
  infer_cells_df[, c("cell_id", "infercnv_group"), drop = FALSE],
  annotation_file,
  col_names = FALSE
)

counts <- get_counts_matrix(obj, assay = "RNA")
counts <- counts[, infer_cells_df$cell_id, drop = FALSE]
gene_order <- readr::read_tsv(
  gene_order_file,
  col_names = c("gene", "chr", "start", "stop"),
  show_col_types = FALSE
)
gene_order <- as.data.frame(gene_order, stringsAsFactors = FALSE)
gene_order <- gene_order[gene_order$gene %in% rownames(counts), , drop = FALSE]
gene_order <- gene_order[!duplicated(gene_order$gene), , drop = FALSE]
if (nrow(gene_order) < 1000L) stop("Too few count genes overlap the inferCNV gene order.", call. = FALSE)
counts <- counts[gene_order$gene, , drop = FALSE]
subset_gene_order_file <- file.path(out_table_dir, "D07_inferCNV_gene_order_subset_used.txt")
readr::write_tsv(gene_order, subset_gene_order_file, col_names = FALSE)

input_audit <- data.frame(
  n_singlet_cells_input = ncol(obj),
  n_cells_used = ncol(counts),
  n_reference_cells_used = sum(infer_cells_df$is_reference),
  n_nonreference_cells_used = sum(!infer_cells_df$is_reference),
  n_genes_used = nrow(counts),
  reference_major_types = paste(reference_major_types, collapse = ";"),
  max_reference_cells_per_type = max_reference_cells_per_type,
  max_nonreference_cells_per_cluster = max_nonreference_cells_per_cluster,
  cutoff = infercnv_cutoff,
  denoise = infercnv_denoise,
  HMM = infercnv_HMM,
  cluster_by_groups = infercnv_cluster_by_groups,
  num_threads = infercnv_num_threads,
  interpretation = "exploratory_CNV_like_support_not_classifier",
  stringsAsFactors = FALSE
)
write_csv(input_audit, "D07_inferCNV_input_parameter_audit.csv")

run_status <- data.frame(status = "not_started", error_message = NA_character_, stringsAsFactors = FALSE)
write_csv(run_status, "D07_inferCNV_run_status.csv")

infer_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = annotation_file,
  delim = "\t",
  gene_order_file = subset_gene_order_file,
  ref_group_names = c("Reference")
)

result <- tryCatch(
  infercnv::run(
    infer_obj,
    cutoff = infercnv_cutoff,
    out_dir = infercnv_out_dir,
    cluster_by_groups = infercnv_cluster_by_groups,
    denoise = infercnv_denoise,
    HMM = infercnv_HMM,
    num_threads = infercnv_num_threads,
    no_prelim_plot = TRUE
  ),
  error = function(e) {
    run_status <<- data.frame(status = "failed", error_message = conditionMessage(e), stringsAsFactors = FALSE)
    write_csv(run_status, "D07_inferCNV_run_status.csv")
    NULL
  }
)
if (is.null(result)) stop("inferCNV failed; inspect D07_inferCNV_run_status.csv.", call. = FALSE)

saveRDS(result, result_rds)
run_status <- data.frame(status = "completed", error_message = NA_character_, stringsAsFactors = FALSE)
write_csv(run_status, "D07_inferCNV_run_status.csv")

expr <- result@expr.data
ann <- infer_cells_df[infer_cells_df$cell_id %in% colnames(expr), , drop = FALSE]
ref_cells <- ann$cell_id[ann$is_reference]
ref_cells <- intersect(ref_cells, colnames(expr))
if (length(ref_cells) < 100L) stop("Too few reference cells in inferCNV output.", call. = FALSE)

ref_mat <- as.matrix(expr[, ref_cells, drop = FALSE])
ref_median <- matrixStats::rowMedians(ref_mat, na.rm = TRUE)
rm(ref_mat)
gc(verbose = FALSE)

cell_signal_vec <- cnv_signal_chunked(expr, ref_median, chunk_size = 500L)
cell_signal <- data.frame(
  cell_id = names(cell_signal_vec),
  infercnv_CNV_like_signal = as.numeric(cell_signal_vec),
  stringsAsFactors = FALSE
)
cell_signal <- merge(cell_signal, ann, by = "cell_id", all.x = TRUE, sort = FALSE)
write_csv(cell_signal, "D07_inferCNV_CNV_like_signal_by_cell.csv")

major_signal <- dplyr::summarise(
  dplyr::group_by(cell_signal, MajorCellType),
  N = dplyr::n(),
  MeanSignal = mean(infercnv_CNV_like_signal, na.rm = TRUE),
  MedianSignal = stats::median(infercnv_CNV_like_signal, na.rm = TRUE),
  Q25 = stats::quantile(infercnv_CNV_like_signal, 0.25, na.rm = TRUE),
  Q75 = stats::quantile(infercnv_CNV_like_signal, 0.75, na.rm = TRUE),
  .groups = "drop"
)
write_csv(major_signal, "D07_inferCNV_CNV_like_signal_by_major_celltype.csv")

cluster_signal <- dplyr::summarise(
  dplyr::group_by(cell_signal, SeuratCluster),
  N = dplyr::n(),
  MeanSignal = mean(infercnv_CNV_like_signal, na.rm = TRUE),
  MedianSignal = stats::median(infercnv_CNV_like_signal, na.rm = TRUE),
  Q25 = stats::quantile(infercnv_CNV_like_signal, 0.25, na.rm = TRUE),
  Q75 = stats::quantile(infercnv_CNV_like_signal, 0.75, na.rm = TRUE),
  .groups = "drop"
)
write_csv(cluster_signal, "D07_inferCNV_CNV_like_signal_by_cluster.csv")

reference_summary <- dplyr::summarise(
  cell_signal[cell_signal$is_reference, , drop = FALSE],
  N_reference = dplyr::n(),
  MeanSignal = mean(infercnv_CNV_like_signal, na.rm = TRUE),
  MedianSignal = stats::median(infercnv_CNV_like_signal, na.rm = TRUE),
  Q75 = stats::quantile(infercnv_CNV_like_signal, 0.75, na.rm = TRUE)
)
write_csv(reference_summary, "D07_inferCNV_reference_signal_summary.csv")

final_summary <- data.frame(
  infercnv_status = "completed",
  n_singlet_cells_input = ncol(obj),
  n_cells_used = ncol(counts),
  n_genes_used = nrow(counts),
  n_reference_cells_used = length(ref_cells),
  reference_median_CNV_like_signal = reference_summary$MedianSignal[[1]],
  result_rds = result_rds,
  runtime_output_dir = infercnv_out_dir,
  interpretation = "exploratory_CNV_like_support_not_classifier",
  stringsAsFactors = FALSE
)
write_csv(final_summary, "D07_inferCNV_final_summary.csv")
writeLines(capture.output(sessionInfo()), file.path(out_table_dir, paste0(SCRIPT_ID, "_sessionInfo.txt")))
message("09A-07 inferCNV support analysis completed.")
