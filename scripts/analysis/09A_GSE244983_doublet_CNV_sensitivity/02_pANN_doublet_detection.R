############################################################
## 09A-02. Per-sample manual pANN doublet detection
##
## Method is a DoubletFinder-like manual pANN workflow without
## calling DoubletFinder package functions. Artificial doublets are
## generated per sample, PCA is recomputed on real + artificial cells,
## and pANN is the artificial-neighbor fraction in PCA space.
##
## Primary call: expected rate 7.5%, homotypic-adjusted.
## Prespecified rate sensitivity: 5%, 7.5%, 10%.
############################################################

options(stringsAsFactors = FALSE)

SCRIPT_ID <- "09A_02_pANN_doublet_detection"
BASE_SEED <- 1234L
PRIMARY_RATE <- 0.075
RATE_SET <- c(0.05, 0.075, 0.10)
PN_ARTIFICIAL <- 0.25
PK_FIXED <- 0.09
PCS_USE <- 1:20
N_VARIABLE_FEATURES <- 2000L
MAX_K_NEIGHBORS <- 300L

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

required_pkgs <- c(
  "Seurat", "SeuratObject", "Matrix", "RANN", "readr", "dplyr"
)
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required R package(s): ", paste(missing_pkgs, collapse = ", "),
    ". No package auto-install is performed.", call. = FALSE
  )
}

input_rds <- file.path(
  project_dir, "results", "intermediate", "GSE244983",
  "GSE244983_seurat_state_localized.rds"
)
out_table_dir <- file.path(
  project_dir, "results", "tables", "GSE244983", "doublet_CNV_sensitivity"
)
out_obj_dir <- file.path(
  project_dir, "results", "objects", "GSE244983", "doublet_CNV_sensitivity"
)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_obj_dir, recursive = TRUE, showWarnings = FALSE)

output_rds <- file.path(out_obj_dir, "GSE244983_all_cells_with_pANN_doublet_calls.rds")

write_csv <- function(x, name) {
  readr::write_csv(x, file.path(out_table_dir, name), na = "")
}

get_counts_matrix <- function(seu, assay = "RNA") {
  if (!assay %in% names(seu@assays)) assay <- Seurat::DefaultAssay(seu)
  mat <- tryCatch(
    SeuratObject::LayerData(seu, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    mat <- tryCatch(
      Seurat::GetAssayData(seu, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(mat)) stop("Could not extract RNA counts matrix.", call. = FALSE)
  mat
}

homotypic_prop <- function(labels) {
  labels <- as.character(labels)
  labels <- labels[!is.na(labels) & labels != ""]
  if (length(labels) == 0L) return(0)
  p <- prop.table(table(labels))
  as.numeric(sum(p^2))
}

make_artificial_doublets <- function(counts, n_artificial, sample_id) {
  n_real <- ncol(counts)
  pair1 <- sample.int(n_real, n_artificial, replace = TRUE)
  pair2 <- sample.int(n_real, n_artificial, replace = TRUE)
  artificial <- counts[, pair1, drop = FALSE] + counts[, pair2, drop = FALSE]
  safe_id <- gsub("[^A-Za-z0-9]+", "_", sample_id)
  colnames(artificial) <- sprintf("ArtificialDoublet_%s_%06d", safe_id, seq_len(n_artificial))
  artificial
}

classify_exact_top <- function(pann, n_expected) {
  cells <- names(pann)
  ord <- order(-as.numeric(pann), cells)
  out <- rep("Singlet", length(pann))
  names(out) <- cells
  take <- min(max(0L, as.integer(n_expected)), length(ord))
  if (take > 0L) out[ord[seq_len(take)]] <- "Doublet"
  out
}

calc_pann_one_sample <- function(counts, sample_clusters, sample_id, sample_index) {
  n_real <- ncol(counts)
  keep_genes <- Matrix::rowSums(counts) > 0
  counts <- counts[keep_genes, , drop = FALSE]

  n_artificial <- max(50L, as.integer(round(n_real * PN_ARTIFICIAL / (1 - PN_ARTIFICIAL))))
  set.seed(BASE_SEED + as.integer(sample_index) * 1000L)
  artificial <- make_artificial_doublets(counts, n_artificial, sample_id)
  combined_counts <- Matrix::cbind2(counts, artificial)
  combined_counts <- as(combined_counts, "dgCMatrix")

  combined <- Seurat::CreateSeuratObject(
    counts = combined_counts,
    project = paste0("pANN_", sample_id)
  )
  Seurat::DefaultAssay(combined) <- "RNA"
  combined <- Seurat::NormalizeData(combined, verbose = FALSE)
  combined <- Seurat::FindVariableFeatures(
    combined, selection.method = "vst", nfeatures = N_VARIABLE_FEATURES,
    verbose = FALSE
  )
  combined <- Seurat::ScaleData(combined, features = Seurat::VariableFeatures(combined), verbose = FALSE)
  combined <- Seurat::RunPCA(
    combined, features = Seurat::VariableFeatures(combined),
    npcs = max(PCS_USE), verbose = FALSE
  )

  emb <- Seurat::Embeddings(combined, reduction = "pca")[, PCS_USE, drop = FALSE]
  k_neighbors <- as.integer(round(PK_FIXED * nrow(emb)))
  k_neighbors <- max(10L, min(k_neighbors, MAX_K_NEIGHBORS, nrow(emb) - 1L))

  real_idx <- seq_len(n_real)
  nn <- RANN::nn2(
    data = emb,
    query = emb[real_idx, , drop = FALSE],
    k = k_neighbors + 1L
  )
  artificial_flag <- c(rep(FALSE, n_real), rep(TRUE, n_artificial))
  pann <- numeric(n_real)
  for (i in seq_len(n_real)) {
    idx <- nn$nn.idx[i, ]
    idx <- idx[idx != i]
    if (length(idx) > k_neighbors) idx <- idx[seq_len(k_neighbors)]
    pann[i] <- mean(artificial_flag[idx])
  }
  names(pann) <- colnames(counts)

  hp <- homotypic_prop(sample_clusters[colnames(counts)])
  rate_audit <- do.call(rbind, lapply(RATE_SET, function(rate) {
    n_raw <- max(1L, as.integer(round(n_real * rate)))
    n_adj <- max(1L, as.integer(round(n_raw * (1 - hp))))
    cls <- classify_exact_top(pann, n_adj)
    data.frame(
      Sample = sample_id,
      expected_rate = rate,
      n_real_cells = n_real,
      n_artificial_doublets = n_artificial,
      pN_artificial = PN_ARTIFICIAL,
      pK_fixed = PK_FIXED,
      k_neighbors = k_neighbors,
      homotypic_prop = hp,
      nExp_raw = n_raw,
      nExp_adjusted = n_adj,
      n_predicted_doublets = sum(cls == "Doublet"),
      observed_fraction = mean(cls == "Doublet"),
      stringsAsFactors = FALSE
    )
  }))

  classifications <- lapply(RATE_SET, function(rate) {
    row <- rate_audit[abs(rate_audit$expected_rate - rate) < 1e-12, , drop = FALSE]
    classify_exact_top(pann, row$nExp_adjusted[[1]])
  })
  names(classifications) <- sprintf("rate_%03d", as.integer(round(RATE_SET * 1000)))

  list(pANN = pann, audit = rate_audit, classifications = classifications)
}

if (!file.exists(input_rds)) stop("Input state-localized RDS not found: ", input_rds, call. = FALSE)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.", call. = FALSE)
if (ncol(obj) != 26053L || nrow(obj) != 11616L) {
  stop("09A-02 requires the current 26,053-cell / 11,616-feature object.", call. = FALSE)
}
required_meta <- c("Sample", "SeuratCluster")
if (!all(required_meta %in% colnames(obj@meta.data))) {
  stop("Sample and SeuratCluster metadata are required for pANN.", call. = FALSE)
}

counts_all <- get_counts_matrix(obj, assay = "RNA")
if (ncol(counts_all) != ncol(obj)) stop("Counts/object cell-count mismatch.", call. = FALSE)

sample_levels <- unique(as.character(obj$Sample))
all_pann <- setNames(rep(NA_real_, ncol(obj)), colnames(obj))
all_cls <- lapply(RATE_SET, function(x) setNames(rep(NA_character_, ncol(obj)), colnames(obj)))
names(all_cls) <- sprintf("rate_%03d", as.integer(round(RATE_SET * 1000)))
audit_list <- vector("list", length(sample_levels))

cluster_vec <- as.character(obj$SeuratCluster)
names(cluster_vec) <- colnames(obj)

for (i in seq_along(sample_levels)) {
  sid <- sample_levels[[i]]
  cells <- colnames(obj)[as.character(obj$Sample) == sid]
  message("09A-02 pANN sample ", sid, " (", length(cells), " cells)")
  res <- calc_pann_one_sample(
    counts = counts_all[, cells, drop = FALSE],
    sample_clusters = cluster_vec,
    sample_id = sid,
    sample_index = i
  )
  all_pann[names(res$pANN)] <- res$pANN
  for (nm in names(all_cls)) {
    all_cls[[nm]][names(res$classifications[[nm]])] <- res$classifications[[nm]]
  }
  audit_list[[i]] <- res$audit
  rm(res)
  gc(verbose = FALSE)
}

if (any(!is.finite(all_pann))) stop("Non-finite or missing pANN values detected.", call. = FALSE)
if (any(vapply(all_cls, function(x) any(is.na(x)), logical(1)))) {
  stop("Missing pANN classification labels detected.", call. = FALSE)
}

audit <- do.call(rbind, audit_list)
write_csv(audit, "D02_pANN_per_sample_rate_audit.csv")

obj@meta.data$pANN_manual <- unname(all_pann[colnames(obj)])
obj@meta.data$pANN_doublet_5pct <- unname(all_cls[["rate_050"]][colnames(obj)])
obj@meta.data$pANN_doublet_primary <- unname(all_cls[["rate_075"]][colnames(obj)])
obj@meta.data$pANN_doublet_10pct <- unname(all_cls[["rate_100"]][colnames(obj)])

call_table <- data.frame(
  Cell = colnames(obj),
  Sample = as.character(obj$Sample),
  SeuratCluster = as.character(obj$SeuratCluster),
  pANN_manual = obj$pANN_manual,
  pANN_doublet_5pct = obj$pANN_doublet_5pct,
  pANN_doublet_primary = obj$pANN_doublet_primary,
  pANN_doublet_10pct = obj$pANN_doublet_10pct,
  stringsAsFactors = FALSE
)
write_csv(call_table, "D02_pANN_cell_calls.csv")

rate_summary <- do.call(rbind, lapply(seq_along(RATE_SET), function(i) {
  nm <- names(all_cls)[[i]]
  data.frame(
    expected_rate = RATE_SET[[i]],
    label_column = c("pANN_doublet_5pct", "pANN_doublet_primary", "pANN_doublet_10pct")[[i]],
    n_predicted_doublets = sum(all_cls[[nm]] == "Doublet"),
    n_retained_singlets = sum(all_cls[[nm]] == "Singlet"),
    predicted_doublet_fraction = mean(all_cls[[nm]] == "Doublet"),
    is_primary = abs(RATE_SET[[i]] - PRIMARY_RATE) < 1e-12,
    stringsAsFactors = FALSE
  )
}))
write_csv(rate_summary, "D02_pANN_global_rate_summary.csv")

pann_quantiles <- do.call(rbind, lapply(sample_levels, function(sid) {
  x <- obj$pANN_manual[as.character(obj$Sample) == sid]
  qs <- stats::quantile(x, probs = c(0, .25, .5, .75, .9, .95, .99, 1), na.rm = TRUE)
  data.frame(Sample = sid, quantile = names(qs), pANN = as.numeric(qs), stringsAsFactors = FALSE)
}))
write_csv(pann_quantiles, "D02_pANN_quantiles_by_sample.csv")

saveRDS(obj, output_rds)
writeLines(capture.output(sessionInfo()), file.path(out_table_dir, paste0(SCRIPT_ID, "_sessionInfo.txt")))
message("09A-02 complete. Primary pANN labels are stored in pANN_doublet_primary.")
