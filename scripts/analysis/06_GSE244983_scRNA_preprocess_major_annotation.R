############################################################
## 06_GSE244983_scRNA_preprocess_major_annotation.R
##
## Purpose
## -------
## Rebuild the primary GSE244983 single-cell preprocessing and
## major-cell-type annotation lineage directly from the canonical
## raw count matrix and author annotation table.
##
## Primary cell universe
## ---------------------
## The public primary Figure 3--6 / Supplementary Figure S5--S10
## lineage is frozen to the 26,053-cell branch established by the
## 07A/07B provenance audits. The historical 25,972-cell branch is
## retained for later sensitivity analyses and is not constructed here.
##
## This script:
##   1) uses explicit canonical raw inputs only;
##   2) does not recursively discover files;
##   3) does not auto-install packages;
##   4) requires Harmony rather than silently falling back to PCA;
##   5) reproduces the historical preprocessing/annotation algorithm;
##   6) writes analysis/source tables only (no manuscript figures);
##   7) applies hard runtime gates before saving the canonical object;
##   8) records input MD5s, cell-set MD5s, parameters, and sessionInfo().
############################################################

options(stringsAsFactors = FALSE)

SCRIPT_ID <- "07_GSE244983_scRNA_preprocess_major_annotation"

############################################################
## 0. Project root
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")

if (!nzchar(project_dir)) {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  required_probe <- c(
    file.path(cwd, "data_raw", "GSE244983", "GSE244983_RawCounts_scRNAseq.txt.gz"),
    file.path(cwd, "data_raw", "GSE244983", "GSE244983_SingleCellAnnotations.txt.gz")
  )

  if (all(file.exists(required_probe))) {
    project_dir <- cwd
  } else {
    stop(
      "Set ICB_PROJECT_DIR to the repository root, or run from the root containing ",
      "data_raw/GSE244983/GSE244983_RawCounts_scRNAseq.txt.gz and ",
      "data_raw/GSE244983/GSE244983_SingleCellAnnotations.txt.gz."
    )
  }
}

project_dir <- normalizePath(
  project_dir,
  winslash = "/",
  mustWork = TRUE
)

message("Project root: ", project_dir)

############################################################
## 1. Packages -- NO auto-install
############################################################

required_pkgs <- c(
  "Seurat",
  "Matrix",
  "data.table",
  "dplyr",
  "harmony"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Restore the repository environment (for example, renv::restore()) and rerun."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(harmony)
})

############################################################
## 2. Canonical paths and helpers
############################################################

raw_count_file <- file.path(
  project_dir,
  "data_raw",
  "GSE244983",
  "GSE244983_RawCounts_scRNAseq.txt.gz"
)

annotation_file <- file.path(
  project_dir,
  "data_raw",
  "GSE244983",
  "GSE244983_SingleCellAnnotations.txt.gz"
)

out_table_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "GSE244983",
  "preprocess_major_annotation"
)

out_intermediate_dir <- file.path(
  project_dir,
  "results",
  "intermediate",
  "GSE244983"
)

log_dir <- file.path(project_dir, "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_count_file)) {
  stop("Missing canonical raw count file: ", raw_count_file)
}
if (!file.exists(annotation_file)) {
  stop("Missing canonical annotation file: ", annotation_file)
}

safe_write_csv <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    x,
    file = path,
    row.names = row.names,
    na = ""
  )
  message("Saved: ", path)
  invisible(path)
}

safe_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
  message("Saved: ", path)
  invisible(path)
}

md5_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

cellset_md5 <- function(ids) {
  ids <- sort(unique(as.character(ids)))
  tf <- tempfile(fileext = ".txt")
  on.exit(unlink(tf), add = TRUE)
  writeLines(ids, tf, useBytes = TRUE)
  unname(tools::md5sum(tf))
}

relative_to_project <- function(path) {
  p <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- paste0(normalizePath(project_dir, winslash = "/", mustWork = TRUE), "/")
  if (startsWith(p, root)) {
    substring(p, nchar(root) + 1L)
  } else {
    p
  }
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  sx <- stats::sd(x, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(0, length(x)))
  (x - mx) / sx
}

get_assay_data_safe <- function(obj, assay = NULL, layer = "data") {
  if (is.null(assay)) assay <- DefaultAssay(obj)

  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e_layer) {
      Seurat::GetAssayData(obj, assay = assay, slot = layer)
    }
  )
}

session_file <- file.path(
  log_dir,
  paste0("sessionInfo_", SCRIPT_ID, ".txt")
)

save_session_info <- function() {
  writeLines(capture.output(sessionInfo()), session_file)
  message("Saved: ", session_file)
}

############################################################
## 3. Frozen parameters and primary-lineage gates
############################################################

## Historical preprocessing parameters preserved from the final mother script.
qc_min_features <- 200L
qc_max_features <- Inf
qc_min_counts <- 0
qc_max_counts <- Inf
qc_max_percent_mt <- 25

normalization_method <- "LogNormalize"
normalization_scale_factor <- 10000
variable_feature_method <- "vst"
n_variable_features <- 2000L
n_pcs_to_compute <- 30L
n_pcs_to_use <- 20L
cluster_resolution <- 0.5
findclusters_random_seed <- 0L
umap_seed <- 42L

## Manual cluster-level override preserved from the historical accepted lineage.
manual_cluster_annotation <- c(
  "12" = "CAF/stromal-like cells"
)

## 07A/07B primary 26,053-cell universe gate.
## 07C then demonstrated exact old-vs-new reproduction for all 26,053 cells:
## identical cell set, exact numeric cluster labels, exact cluster partition,
## exact major-cell annotations, and cluster-12 Jaccard = 1.
expected_n_cells <- 26053L
expected_n_features <- 11616L
expected_cellset_md5 <- "a574b6a23054be747b74e2e09f5ddea9"

expected_sample_counts <- c(
  "Pat_ICBnaive1" = 7414L,
  "Pat_ICBnaive2" = 9081L,
  "Pat42" = 7188L,
  "Pat5" = 2370L
)

## Frozen major-cell-type counts verified by 07C against
## data_processed/GSE244983_seurat_major_annotated_LOCKED.rds.
expected_major_counts <- c(
  "Malignant" = 17708L,
  "Cycling malignant" = 1888L,
  "CAF/stromal-like cells" = 333L,
  "Endothelial" = 1586L,
  "Myeloid cells" = 231L,
  "T/NK cells" = 2973L,
  "T/NK/Treg-like cells" = 382L,
  "B/Plasma cells" = 952L,
  "Cycling cells" = 0L
)

## 17,708 malignant + 1,888 cycling malignant = 19,596.
expected_malignant_lineage_cells <- 19596L

parameter_table <- data.frame(
  Parameter = c(
    "qc_min_features",
    "qc_max_features",
    "qc_min_counts",
    "qc_max_counts",
    "qc_max_percent_mt",
    "normalization_method",
    "normalization_scale_factor",
    "variable_feature_method",
    "n_variable_features",
    "n_pcs_to_compute",
    "n_pcs_to_use",
    "harmony_group_by",
    "cluster_resolution",
    "findclusters_random_seed",
    "umap_seed",
    "manual_cluster_12_annotation"
  ),
  Value = c(
    as.character(qc_min_features),
    as.character(qc_max_features),
    as.character(qc_min_counts),
    as.character(qc_max_counts),
    as.character(qc_max_percent_mt),
    normalization_method,
    as.character(normalization_scale_factor),
    variable_feature_method,
    as.character(n_variable_features),
    as.character(n_pcs_to_compute),
    as.character(n_pcs_to_use),
    "Sample",
    as.character(cluster_resolution),
    as.character(findclusters_random_seed),
    as.character(umap_seed),
    unname(manual_cluster_annotation[["12"]])
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  parameter_table,
  file.path(out_table_dir, "GSE244983_preprocessing_parameters.csv")
)

input_audit <- data.frame(
  Input = c("Raw count matrix", "Author single-cell annotation"),
  RelativePath = c(
    relative_to_project(raw_count_file),
    relative_to_project(annotation_file)
  ),
  Exists = c(file.exists(raw_count_file), file.exists(annotation_file)),
  MD5 = c(md5_file(raw_count_file), md5_file(annotation_file)),
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_audit,
  file.path(out_table_dir, "GSE244983_input_file_audit.csv")
)

############################################################
## 4. Raw-count header and author-annotation barcode audit
############################################################

message("Reading raw-count header...")
raw_header <- data.table::fread(
  raw_count_file,
  nrows = 0,
  check.names = FALSE,
  showProgress = FALSE
)

raw_columns <- colnames(raw_header)
if (length(raw_columns) < 2L) {
  stop("Raw count file contains fewer than two columns.")
}

raw_gene_column <- raw_columns[[1]]
raw_cell_ids <- raw_columns[-1]

raw_header_audit <- data.frame(
  GeneColumn = raw_gene_column,
  N_cell_columns = length(raw_cell_ids),
  N_unique_cell_columns = length(unique(raw_cell_ids)),
  N_duplicate_cell_columns = sum(duplicated(raw_cell_ids)),
  CellSetMD5 = cellset_md5(raw_cell_ids),
  stringsAsFactors = FALSE
)

safe_write_csv(
  raw_header_audit,
  file.path(out_table_dir, "GSE244983_raw_count_header_audit.csv")
)

message("Reading author annotation table...")
anno_dt <- data.table::fread(
  annotation_file,
  check.names = FALSE,
  showProgress = FALSE
)
anno_df <- as.data.frame(anno_dt, check.names = FALSE)

barcode_scores <- lapply(colnames(anno_df), function(cn) {
  vals <- as.character(anno_df[[cn]])
  vals_valid <- vals[!is.na(vals) & vals != ""]

  data.frame(
    Column = cn,
    N_nonmissing = length(vals_valid),
    N_unique = length(unique(vals_valid)),
    DirectOverlapWithRaw = length(intersect(raw_cell_ids, vals_valid)),
    DotDashNormalizedOverlapWithRaw = length(
      intersect(gsub("\\.", "-", raw_cell_ids), gsub("\\.", "-", vals_valid))
    ),
    stringsAsFactors = FALSE
  )
})

barcode_scores <- do.call(rbind, barcode_scores)
barcode_scores$BestOverlap <- pmax(
  barcode_scores$DirectOverlapWithRaw,
  barcode_scores$DotDashNormalizedOverlapWithRaw
)
barcode_scores <- barcode_scores[
  order(-barcode_scores$DirectOverlapWithRaw, -barcode_scores$BestOverlap),
  ,
  drop = FALSE
]
rownames(barcode_scores) <- NULL

safe_write_csv(
  barcode_scores,
  file.path(out_table_dir, "GSE244983_annotation_barcode_column_audit.csv")
)

best_barcode_col <- barcode_scores$Column[[1]]
annotation_barcodes <- as.character(anno_df[[best_barcode_col]])

direct_overlap <- intersect(raw_cell_ids, annotation_barcodes)
normalized_overlap <- intersect(
  gsub("\\.", "-", raw_cell_ids),
  gsub("\\.", "-", annotation_barcodes)
)

if (length(direct_overlap) != expected_n_cells) {
  stop(
    "Canonical author annotation did not reproduce the expected 26,053 direct barcode matches. ",
    "Best barcode column: ", best_barcode_col,
    "; direct overlap: ", length(direct_overlap), "."
  )
}

if (anyDuplicated(annotation_barcodes[annotation_barcodes %in% raw_cell_ids]) > 0L) {
  stop("Author annotation contains duplicated barcodes among raw-count cells.")
}

annotation_index <- match(raw_cell_ids, annotation_barcodes)
if (anyNA(annotation_index)) {
  stop("At least one raw-count cell lacks a directly matched author annotation row.")
}

annotation_aligned <- anno_df[annotation_index, , drop = FALSE]
rownames(annotation_aligned) <- raw_cell_ids

annotation_overlap_audit <- data.frame(
  AnnotationRows = nrow(anno_df),
  BestBarcodeColumn = best_barcode_col,
  DirectOverlapWithRaw = length(direct_overlap),
  DotDashNormalizedOverlapWithRaw = length(normalized_overlap),
  MatchModeUsed = "direct",
  MatchedRate = length(direct_overlap) / length(raw_cell_ids),
  stringsAsFactors = FALSE
)

safe_write_csv(
  annotation_overlap_audit,
  file.path(out_table_dir, "GSE244983_annotation_overlap_audit.csv")
)

############################################################
## 5. Read count matrix and create the primary Seurat object
############################################################

message("Reading full GSE244983 raw count matrix. This step is memory intensive...")
counts_dt <- data.table::fread(
  raw_count_file,
  check.names = FALSE,
  showProgress = TRUE
)

if (ncol(counts_dt) != length(raw_cell_ids) + 1L) {
  stop("Full raw count matrix column count differs from the audited header.")
}

if (!identical(colnames(counts_dt)[-1], raw_cell_ids)) {
  stop("Raw count cell-column order changed between header audit and full read.")
}

genes_raw <- as.character(counts_dt[[1]])
if (anyNA(genes_raw) || any(genes_raw == "")) {
  stop("Raw count gene column contains missing or empty feature names.")
}

genes_unique <- make.unique(genes_raw)

count_mat <- as.matrix(counts_dt[, -1, with = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- genes_unique
colnames(count_mat) <- raw_cell_ids

count_sparse <- Matrix::Matrix(count_mat, sparse = TRUE)

rm(count_mat, counts_dt)
gc()

obj <- Seurat::CreateSeuratObject(
  counts = count_sparse,
  project = "GSE244983",
  min.cells = 3,
  min.features = 0
)
rm(count_sparse)
gc()

DefaultAssay(obj) <- "RNA"
obj$Dataset <- "GSE244983"
obj$InputSource <- "GSE244983_RawCounts_scRNAseq.txt.gz"

## Add all author metadata in exact Seurat-cell order.
annotation_for_object <- annotation_aligned[colnames(obj), , drop = FALSE]
if (!identical(rownames(annotation_for_object), colnames(obj))) {
  stop("Aligned annotation rownames do not match the Seurat cell order.")
}

obj <- Seurat::AddMetaData(obj, metadata = annotation_for_object)

if (!"Sample" %in% colnames(obj@meta.data)) {
  stop("Canonical author annotation did not provide the required 'Sample' metadata column.")
}
obj$Sample <- as.character(obj$Sample)

post_create_audit <- data.frame(
  N_cells = ncol(obj),
  N_features = nrow(obj),
  CellSetMD5 = cellset_md5(colnames(obj)),
  N_samples = length(unique(obj$Sample)),
  stringsAsFactors = FALSE
)

safe_write_csv(
  post_create_audit,
  file.path(out_table_dir, "GSE244983_post_creation_audit.csv")
)

############################################################
## 6. Mitochondrial field and frozen QC rule
############################################################

all_genes <- rownames(obj)
mt_genes_upper <- grep("^MT-", all_genes, value = TRUE)
mt_genes_lower <- grep("^mt-", all_genes, value = TRUE)

if (length(mt_genes_upper) >= length(mt_genes_lower)) {
  mt_pattern_used <- "^MT-"
  mt_genes_used <- mt_genes_upper
} else {
  mt_pattern_used <- "^mt-"
  mt_genes_used <- mt_genes_lower
}

if (length(mt_genes_used) > 0L) {
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern_used)
  mt_available <- TRUE
} else {
  ## Historical final lineage: mitochondrial percentage is unavailable/non-informative.
  obj$percent.mt <- 0
  mt_available <- FALSE
}

mitochondrial_audit <- data.frame(
  Pattern = c("^MT-", "^mt-", "used"),
  N_genes = c(length(mt_genes_upper), length(mt_genes_lower), length(mt_genes_used)),
  stringsAsFactors = FALSE
)

safe_write_csv(
  mitochondrial_audit,
  file.path(out_table_dir, "GSE244983_mitochondrial_feature_audit.csv")
)

qc_before <- data.frame(
  N_cells = ncol(obj),
  Median_nFeature_RNA = stats::median(obj$nFeature_RNA),
  Median_nCount_RNA = stats::median(obj$nCount_RNA),
  Median_percent_mt = stats::median(obj$percent.mt),
  Min_nFeature_RNA = min(obj$nFeature_RNA),
  Max_nFeature_RNA = max(obj$nFeature_RNA),
  Min_nCount_RNA = min(obj$nCount_RNA),
  Max_nCount_RNA = max(obj$nCount_RNA),
  Min_percent_mt = min(obj$percent.mt),
  Max_percent_mt = max(obj$percent.mt),
  stringsAsFactors = FALSE
)

safe_write_csv(
  qc_before,
  file.path(out_table_dir, "GSE244983_QC_summary_before_filtering.csv")
)

keep_cells <-
  obj$nFeature_RNA >= qc_min_features &
  obj$nFeature_RNA <= qc_max_features &
  obj$nCount_RNA >= qc_min_counts &
  obj$nCount_RNA <= qc_max_counts &
  obj$percent.mt <= qc_max_percent_mt

qc_cell_status <- data.frame(
  Cell = colnames(obj),
  nFeature_RNA = obj$nFeature_RNA,
  nCount_RNA = obj$nCount_RNA,
  percent_mt = obj$percent.mt,
  Pass_QC = keep_cells,
  stringsAsFactors = FALSE
)

safe_write_csv(
  qc_cell_status,
  file.path(out_table_dir, "GSE244983_cell_level_QC_status.csv")
)

obj <- subset(obj, cells = colnames(obj)[keep_cells])

qc_after <- data.frame(
  N_cells = ncol(obj),
  N_features = nrow(obj),
  CellSetMD5 = cellset_md5(colnames(obj)),
  Median_nFeature_RNA = stats::median(obj$nFeature_RNA),
  Median_nCount_RNA = stats::median(obj$nCount_RNA),
  Median_percent_mt = stats::median(obj$percent.mt),
  Min_nFeature_RNA = min(obj$nFeature_RNA),
  Max_nFeature_RNA = max(obj$nFeature_RNA),
  Min_nCount_RNA = min(obj$nCount_RNA),
  Max_nCount_RNA = max(obj$nCount_RNA),
  Min_percent_mt = min(obj$percent.mt),
  Max_percent_mt = max(obj$percent.mt),
  stringsAsFactors = FALSE
)

safe_write_csv(
  qc_after,
  file.path(out_table_dir, "GSE244983_QC_summary_after_filtering.csv")
)

############################################################
## 7. Normalization, PCA, Harmony, clustering, and UMAP
############################################################

message("Normalizing and identifying variable features...")
obj <- Seurat::NormalizeData(
  obj,
  normalization.method = normalization_method,
  scale.factor = normalization_scale_factor,
  verbose = FALSE
)

obj <- Seurat::FindVariableFeatures(
  obj,
  selection.method = variable_feature_method,
  nfeatures = n_variable_features,
  verbose = FALSE
)

## Historical accepted lineage scaled all retained features.
obj <- Seurat::ScaleData(
  obj,
  features = rownames(obj),
  verbose = FALSE
)

obj <- Seurat::RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = n_pcs_to_compute,
  verbose = FALSE
)

available_pcs <- ncol(Embeddings(obj, "pca"))
if (available_pcs < n_pcs_to_use) {
  stop(
    "PCA returned only ", available_pcs,
    " components; expected at least ", n_pcs_to_use, "."
  )
}
dims_use <- seq_len(n_pcs_to_use)

message("Running Harmony using Sample...")
obj <- harmony::RunHarmony(
  object = obj,
  group.by.vars = "Sample",
  reduction.use = "pca",
  dims.use = dims_use,
  verbose = FALSE
)

if (!"harmony" %in% Reductions(obj)) {
  stop("Harmony reduction was not created.")
}

obj <- Seurat::FindNeighbors(
  obj,
  reduction = "harmony",
  dims = dims_use,
  verbose = FALSE
)

obj <- Seurat::FindClusters(
  obj,
  resolution = cluster_resolution,
  random.seed = findclusters_random_seed,
  verbose = FALSE
)

## Preserve the historical explicit pre-UMAP seed call while also fixing the
## Seurat UMAP seed argument to its historical default for version robustness.
set.seed(1234)
obj <- Seurat::RunUMAP(
  obj,
  reduction = "harmony",
  dims = dims_use,
  seed.use = umap_seed,
  verbose = FALSE
)

obj$SeuratCluster <- as.character(Idents(obj))

cluster_counts <- as.data.frame(
  table(obj$SeuratCluster),
  stringsAsFactors = FALSE
)
colnames(cluster_counts) <- c("SeuratCluster", "N_cells")
cluster_counts$SeuratCluster <- as.character(cluster_counts$SeuratCluster)
cluster_counts <- cluster_counts[
  order(suppressWarnings(as.integer(cluster_counts$SeuratCluster))),
  ,
  drop = FALSE
]

safe_write_csv(
  cluster_counts,
  file.path(out_table_dir, "GSE244983_cluster_cell_counts.csv")
)

umap_coords <- as.data.frame(Seurat::Embeddings(obj, "umap"))
umap_coords$Cell <- rownames(umap_coords)
umap_coords$Sample <- as.character(obj$Sample)
umap_coords$SeuratCluster <- as.character(obj$SeuratCluster)
umap_coords <- umap_coords[
  ,
  c("Cell", "Sample", "SeuratCluster", setdiff(colnames(umap_coords), c("Cell", "Sample", "SeuratCluster"))),
  drop = FALSE
]

safe_write_csv(
  umap_coords,
  file.path(out_table_dir, "GSE244983_UMAP_coordinates.csv")
)

############################################################
## 8. Frozen marker-module major-cell-type annotation
############################################################

marker_sets <- list(
  Malignant = c(
    "MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10", "S100B", "MIA", "TFAP2A"
  ),
  T_NK = c(
    "PTPRC", "CD3D", "CD3E", "CD2", "TRAC", "NKG7", "GNLY", "GZMB", "PRF1",
    "FOXP3", "IL2RA", "CTLA4", "TIGIT"
  ),
  B_Plasma = c(
    "PTPRC", "MS4A1", "CD79A", "CD79B", "CD74", "MZB1", "JCHAIN", "IGHG1", "IGKC"
  ),
  Myeloid = c(
    "PTPRC", "LYZ", "LST1", "TYROBP", "CD14", "FCGR3A", "CD68", "C1QA", "C1QB",
    "S100A8", "S100A9", "IL1B"
  ),
  Endothelial = c(
    "PECAM1", "VWF", "KDR", "FLT1", "CLDN5", "RAMP2", "ESAM"
  ),
  CAF_stromal_like = c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3",
    "DCN", "LUM", "FAP", "ACTA2", "PDPN", "THY1", "VCAN", "TNC", "TGFB1",
    "FBN1", "MRC2", "ADAM12", "TIMP1", "INHBA"
  ),
  Cycling = c(
    "MKI67", "TOP2A", "UBE2C", "PCNA", "TYMS", "HMGB2", "STMN1", "CENPF"
  )
)

marker_display <- c(
  Malignant = "Malignant",
  T_NK = "T/NK cells",
  B_Plasma = "B/Plasma cells",
  Myeloid = "Myeloid cells",
  Endothelial = "Endothelial",
  CAF_stromal_like = "CAF/stromal-like cells",
  Cycling = "Cycling"
)

major_celltype_levels <- names(expected_major_counts)

expr_data <- get_assay_data_safe(obj, assay = "RNA", layer = "data")

gene_presence_list <- list()

for (set_name in names(marker_sets)) {
  genes <- marker_sets[[set_name]]
  present <- intersect(genes, rownames(expr_data))

  gene_presence_list[[set_name]] <- data.frame(
    MarkerSet = set_name,
    DisplayName = unname(marker_display[[set_name]]),
    Gene = genes,
    Present = genes %in% present,
    stringsAsFactors = FALSE
  )

  score_col <- paste0(set_name, "_score")
  score_z_col <- paste0(set_name, "_score_z")

  if (length(present) == 0L) {
    obj[[score_col]] <- NA_real_
    obj[[score_z_col]] <- NA_real_
  } else {
    score <- Matrix::colMeans(expr_data[present, , drop = FALSE])
    obj[[score_col]] <- as.numeric(score)
    obj[[score_z_col]] <- zscore_safe(score)
  }

  message(
    "Marker set ", set_name, ": ",
    length(present), "/", length(genes), " genes present"
  )
}

gene_presence <- dplyr::bind_rows(gene_presence_list)
safe_write_csv(
  gene_presence,
  file.path(out_table_dir, "GSE244983_major_annotation_marker_gene_presence.csv")
)

score_cols_z <- paste0(names(marker_sets), "_score_z")

cluster_scores <- obj@meta.data |>
  dplyr::mutate(SeuratCluster = as.character(SeuratCluster)) |>
  dplyr::group_by(SeuratCluster) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    dplyr::across(dplyr::all_of(score_cols_z), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

non_cycling_sets <- c(
  "Malignant",
  "T_NK",
  "B_Plasma",
  "Myeloid",
  "Endothelial",
  "CAF_stromal_like"
)
non_cycling_score_cols <- paste0(non_cycling_sets, "_score_z")

auto_annotate_cluster <- function(row_df) {
  scores <- vapply(
    non_cycling_score_cols,
    function(cn) as.numeric(row_df[[cn]][[1]]),
    numeric(1)
  )
  names(scores) <- non_cycling_sets
  scores[is.na(scores)] <- -Inf

  top_set <- names(which.max(scores))
  top_label <- unname(marker_display[[top_set]])

  cycle_score <- as.numeric(row_df[["Cycling_score_z"]][[1]])
  malignant_score <- as.numeric(row_df[["Malignant_score_z"]][[1]])
  tnk_score <- as.numeric(row_df[["T_NK_score_z"]][[1]])

  if (is.na(cycle_score)) cycle_score <- -Inf
  if (is.na(malignant_score)) malignant_score <- -Inf
  if (is.na(tnk_score)) tnk_score <- -Inf

  if (top_set == "Malignant" && cycle_score >= 0.75) {
    return("Cycling malignant")
  }

  if (top_set == "T_NK" && cycle_score >= 0.75 && tnk_score >= 0.25) {
    return("T/NK/Treg-like cells")
  }

  if (cycle_score >= 1.00 && malignant_score < 0.25) {
    return("Cycling cells")
  }

  top_label
}

auto_labels <- vapply(
  seq_len(nrow(cluster_scores)),
  function(i) auto_annotate_cluster(cluster_scores[i, , drop = FALSE]),
  character(1)
)

cluster_annotation <- cluster_scores
cluster_annotation$AutoMajorCellType <- auto_labels
cluster_annotation$MajorCellType <- cluster_annotation$AutoMajorCellType
cluster_annotation$ManualOverride <- FALSE

for (cl in names(manual_cluster_annotation)) {
  hit <- cluster_annotation$SeuratCluster == cl
  if (any(hit)) {
    cluster_annotation$MajorCellType[hit] <- unname(manual_cluster_annotation[[cl]])
    cluster_annotation$ManualOverride[hit] <- TRUE
  }
}

if (any(!cluster_annotation$MajorCellType %in% major_celltype_levels)) {
  stop(
    "Unexpected major-cell-type label(s): ",
    paste(
      unique(cluster_annotation$MajorCellType[
        !cluster_annotation$MajorCellType %in% major_celltype_levels
      ]),
      collapse = ", "
    )
  )
}

cluster_annotation$MajorCellType <- factor(
  cluster_annotation$MajorCellType,
  levels = major_celltype_levels
)

safe_write_csv(
  cluster_annotation,
  file.path(out_table_dir, "GSE244983_cluster_to_major_celltype.csv")
)

cluster_to_type <- setNames(
  as.character(cluster_annotation$MajorCellType),
  cluster_annotation$SeuratCluster
)

major_type_vec <- unname(cluster_to_type[as.character(obj$SeuratCluster)])
if (length(major_type_vec) != ncol(obj) || anyNA(major_type_vec)) {
  stop("Failed to assign a major cell type to every cell.")
}

obj$MajorCellType <- factor(
  major_type_vec,
  levels = major_celltype_levels
)

############################################################
## 9. Canonical source tables
############################################################

sample_counts <- as.data.frame(
  table(as.character(obj$Sample)),
  stringsAsFactors = FALSE
)
colnames(sample_counts) <- c("Sample", "N_cells")
sample_counts$Sample <- as.character(sample_counts$Sample)
sample_counts <- sample_counts[order(sample_counts$Sample), , drop = FALSE]
rownames(sample_counts) <- NULL

safe_write_csv(
  sample_counts,
  file.path(out_table_dir, "GSE244983_sample_cell_counts.csv")
)

major_counts_observed <- table(
  factor(as.character(obj$MajorCellType), levels = major_celltype_levels)
)
major_counts <- data.frame(
  MajorCellType = major_celltype_levels,
  N_cells = as.integer(major_counts_observed[major_celltype_levels]),
  Percentage = as.integer(major_counts_observed[major_celltype_levels]) / ncol(obj) * 100,
  stringsAsFactors = FALSE
)

safe_write_csv(
  major_counts,
  file.path(out_table_dir, "GSE244983_major_celltype_counts.csv")
)

cell_composition <- obj@meta.data |>
  dplyr::mutate(
    Sample = as.character(Sample),
    MajorCellType = as.character(MajorCellType)
  ) |>
  dplyr::count(Sample, MajorCellType, name = "N_cells") |>
  dplyr::group_by(Sample) |>
  dplyr::mutate(Proportion = N_cells / sum(N_cells)) |>
  dplyr::ungroup()

safe_write_csv(
  cell_composition,
  file.path(out_table_dir, "GSE244983_cell_composition_by_sample.csv")
)

cell_metadata_source <- data.frame(
  Cell = colnames(obj),
  Sample = as.character(obj$Sample),
  SeuratCluster = as.character(obj$SeuratCluster),
  MajorCellType = as.character(obj$MajorCellType),
  nFeature_RNA = obj$nFeature_RNA,
  nCount_RNA = obj$nCount_RNA,
  percent_mt = obj$percent.mt,
  stringsAsFactors = FALSE
)

safe_write_csv(
  cell_metadata_source,
  file.path(out_table_dir, "GSE244983_cell_metadata_source.csv")
)

############################################################
## 10. Hard reproducibility gates
############################################################

observed_sample_counts <- setNames(sample_counts$N_cells, sample_counts$Sample)

observed_major_counts_raw <- setNames(
  major_counts$N_cells,
  major_counts$MajorCellType
)

## Reindex against the complete frozen category universe and fill absent
## zero-count categories explicitly with 0 rather than NA.
observed_major_counts <- setNames(
  integer(length(expected_major_counts)),
  names(expected_major_counts)
)

common_major_names <- intersect(
  names(observed_major_counts_raw),
  names(expected_major_counts)
)

observed_major_counts[common_major_names] <-
  as.integer(observed_major_counts_raw[common_major_names])

malignant_lineage_n <- sum(
  as.character(obj$MajorCellType) %in% c("Malignant", "Cycling malignant")
)

gates <- data.frame(
  Gate = c(
    "Raw count cell columns = 26,053",
    "Raw count cell columns unique",
    "Raw cell-set MD5 matches primary lineage",
    "Author annotation direct overlap = 26,053",
    "Post-creation cells = 26,053",
    "Post-creation features = 11,616",
    "Post-creation cell-set MD5 matches primary lineage",
    "QC removes zero cells",
    "Post-QC cells = 26,053",
    "Post-QC cell-set MD5 matches primary lineage",
    "Harmony reduction exists",
    "UMAP reduction exists",
    "No missing Sample values",
    "Sample counts match frozen 26,053 lineage",
    "Every marker set has at least one present gene",
    "Every cell has a major-cell-type annotation",
    "Major-cell-type counts match frozen 26,053 lineage",
    "Malignant + cycling malignant = 19,596"
  ),
  Observed = c(
    as.character(length(raw_cell_ids)),
    as.character(sum(duplicated(raw_cell_ids)) == 0L),
    cellset_md5(raw_cell_ids),
    as.character(length(direct_overlap)),
    as.character(post_create_audit$N_cells[[1]]),
    as.character(post_create_audit$N_features[[1]]),
    post_create_audit$CellSetMD5[[1]],
    as.character(sum(!keep_cells)),
    as.character(ncol(obj)),
    cellset_md5(colnames(obj)),
    as.character("harmony" %in% Reductions(obj)),
    as.character("umap" %in% Reductions(obj)),
    as.character(sum(is.na(obj$Sample) | obj$Sample == "")),
    paste(paste0(names(observed_sample_counts), "=", observed_sample_counts), collapse = ";"),
    paste(
      vapply(
        names(marker_sets),
        function(s) sum(gene_presence$MarkerSet == s & gene_presence$Present),
        integer(1)
      ),
      collapse = ";"
    ),
    as.character(sum(is.na(obj$MajorCellType))),
    paste(paste0(names(observed_major_counts), "=", observed_major_counts), collapse = ";"),
    as.character(malignant_lineage_n)
  ),
  Expected = c(
    as.character(expected_n_cells),
    "TRUE",
    expected_cellset_md5,
    as.character(expected_n_cells),
    as.character(expected_n_cells),
    as.character(expected_n_features),
    expected_cellset_md5,
    "0",
    as.character(expected_n_cells),
    expected_cellset_md5,
    "TRUE",
    "TRUE",
    "0",
    paste(paste0(names(expected_sample_counts), "=", expected_sample_counts), collapse = ";"),
    ">0 genes in each marker set",
    "0",
    paste(paste0(names(expected_major_counts), "=", expected_major_counts), collapse = ";"),
    as.character(expected_malignant_lineage_cells)
  ),
  Pass = c(
    length(raw_cell_ids) == expected_n_cells,
    sum(duplicated(raw_cell_ids)) == 0L,
    identical(cellset_md5(raw_cell_ids), expected_cellset_md5),
    length(direct_overlap) == expected_n_cells,
    post_create_audit$N_cells[[1]] == expected_n_cells,
    post_create_audit$N_features[[1]] == expected_n_features,
    identical(post_create_audit$CellSetMD5[[1]], expected_cellset_md5),
    sum(!keep_cells) == 0L,
    ncol(obj) == expected_n_cells,
    identical(cellset_md5(colnames(obj)), expected_cellset_md5),
    "harmony" %in% Reductions(obj),
    "umap" %in% Reductions(obj),
    sum(is.na(obj$Sample) | obj$Sample == "") == 0L,
    setequal(names(observed_sample_counts), names(expected_sample_counts)) &&
      identical(
        as.integer(observed_sample_counts[names(expected_sample_counts)]),
        as.integer(expected_sample_counts)
      ),
    all(
      vapply(
        names(marker_sets),
        function(s) sum(gene_presence$MarkerSet == s & gene_presence$Present) > 0L,
        logical(1)
      )
    ),
    sum(is.na(obj$MajorCellType)) == 0L,
    identical(
      as.integer(observed_major_counts[names(expected_major_counts)]),
      as.integer(expected_major_counts)
    ),
    malignant_lineage_n == expected_malignant_lineage_cells
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  gates,
  file.path(out_table_dir, "GSE244983_reproducibility_gates.csv")
)

if (!all(gates$Pass)) {
  failed <- gates$Gate[!gates$Pass]
  save_session_info()
  stop(
    "07 reproducibility gate failure: ",
    paste(failed, collapse = " | "),
    ". Canonical Seurat object was NOT saved."
  )
}

############################################################
## 11. Save canonical intermediate object only after all gates pass
############################################################

canonical_object_file <- file.path(
  out_intermediate_dir,
  "GSE244983_seurat_major_annotated.rds"
)

safe_save_rds(obj, canonical_object_file)

canonical_object_audit <- data.frame(
  RelativePath = relative_to_project(canonical_object_file),
  MD5 = md5_file(canonical_object_file),
  N_cells = ncol(obj),
  N_features = nrow(obj),
  CellSetMD5 = cellset_md5(colnames(obj)),
  N_clusters = length(unique(obj$SeuratCluster)),
  N_major_cell_types_observed = length(unique(as.character(obj$MajorCellType))),
  MalignantLineageCells = malignant_lineage_n,
  stringsAsFactors = FALSE
)

safe_write_csv(
  canonical_object_audit,
  file.path(out_table_dir, "GSE244983_canonical_object_audit.csv")
)

############################################################
## 12. Session information and explicit output inventory
############################################################

save_session_info()

generated_files <- c(
  file.path(out_table_dir, "GSE244983_preprocessing_parameters.csv"),
  file.path(out_table_dir, "GSE244983_input_file_audit.csv"),
  file.path(out_table_dir, "GSE244983_raw_count_header_audit.csv"),
  file.path(out_table_dir, "GSE244983_annotation_barcode_column_audit.csv"),
  file.path(out_table_dir, "GSE244983_annotation_overlap_audit.csv"),
  file.path(out_table_dir, "GSE244983_post_creation_audit.csv"),
  file.path(out_table_dir, "GSE244983_mitochondrial_feature_audit.csv"),
  file.path(out_table_dir, "GSE244983_QC_summary_before_filtering.csv"),
  file.path(out_table_dir, "GSE244983_cell_level_QC_status.csv"),
  file.path(out_table_dir, "GSE244983_QC_summary_after_filtering.csv"),
  file.path(out_table_dir, "GSE244983_cluster_cell_counts.csv"),
  file.path(out_table_dir, "GSE244983_UMAP_coordinates.csv"),
  file.path(out_table_dir, "GSE244983_major_annotation_marker_gene_presence.csv"),
  file.path(out_table_dir, "GSE244983_cluster_to_major_celltype.csv"),
  file.path(out_table_dir, "GSE244983_sample_cell_counts.csv"),
  file.path(out_table_dir, "GSE244983_major_celltype_counts.csv"),
  file.path(out_table_dir, "GSE244983_cell_composition_by_sample.csv"),
  file.path(out_table_dir, "GSE244983_cell_metadata_source.csv"),
  file.path(out_table_dir, "GSE244983_reproducibility_gates.csv"),
  file.path(out_table_dir, "GSE244983_canonical_object_audit.csv"),
  canonical_object_file,
  session_file
)

generated_files <- generated_files[file.exists(generated_files)]

output_inventory <- data.frame(
  Type = ifelse(
    grepl("\\.rds$", generated_files, ignore.case = TRUE),
    "RDS",
    ifelse(grepl("\\.txt$", generated_files, ignore.case = TRUE), "log", "table")
  ),
  RelativePath = vapply(generated_files, relative_to_project, character(1)),
  SizeBytes = file.info(generated_files)$size,
  MD5 = vapply(generated_files, md5_file, character(1)),
  stringsAsFactors = FALSE
)

safe_write_csv(
  output_inventory,
  file.path(out_table_dir, "GSE244983_output_inventory.csv")
)

cat("\n================ 07 SUMMARY ================\n")
cat("Primary cells: ", ncol(obj), "\n", sep = "")
cat("Features: ", nrow(obj), "\n", sep = "")
cat("Cell-set MD5: ", cellset_md5(colnames(obj)), "\n", sep = "")
cat("Clusters: ", length(unique(obj$SeuratCluster)), "\n", sep = "")
cat("Malignant-lineage cells: ", malignant_lineage_n, "\n", sep = "")
cat("All reproducibility gates: PASS\n")
cat("Canonical object: ", relative_to_project(canonical_object_file), "\n", sep = "")
cat("============================================\n")
