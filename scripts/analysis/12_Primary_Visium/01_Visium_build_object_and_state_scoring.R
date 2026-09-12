############################################################
## STYLE-ONLY v1.1 S20-REFERENCE LOCK
## Figure typography is harmonized to the accepted S20 visual ratio.
## Long subtitles are locally reduced and/or wrapped when needed.
## Scientific inputs, thresholds, statistics, state definitions,
## analytical logic, numerical outputs, and figure canvas sizes are unchanged.
############################################################

############################################################
## 01_Visium_build_object_and_state_scoring.R
## PUBLIC FIGURE STYLE CONTRACT
##   Internal analytical state IDs are unchanged.
##   Display labels are lowercase-first:
##     immune-defective/cold
##     myeloid–Treg immunosuppressive
##     tumor-dedifferentiation/stromal-remodeling
##     melanocytic differentiation
##   Fixed state palette:
##     #4DBBD5 / #00A087 / #E64B35 / #3C5488
##   Manuscript-facing figures use sans, base size 10, white background,
##   light major grid, no minor grid, and black panel border where applicable.
##   No numerical/statistical definition is changed by the style patch.

## Project: ICB resistance / melanoma tumor–immune states
## Purpose:
##   Clean-room construction of a Visium spatial Seurat object from
##   raw 10x spatial files, followed by final ICBcomb state scoring.
##
## This script should be run BEFORE Step 11 Moran's I / colocalization
## if no processed spatial RDS containing both spatial coordinates and
## four final state-score columns is available.
############################################################

options(stringsAsFactors = FALSE)
set.seed(1234)

############################################################
## 0. Project paths
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir,winslash="/",mustWork=TRUE)

raw_spatial_dir <- file.path(
  project_dir,
  "data_raw",
  "spatial_10x_human_melanoma_IF_FFPE"
)

gmt_file <- file.path(
  project_dir,
  "results",
  "tables",
  "ICBcomb_final_input_gene_sets_CLEAN",
  "02_final_ICBcomb_gene_sets",
  "final_ICBcomb_gene_sets_CLEAN.gmt"
)

out_tab_dir <- file.path(project_dir, "results", "tables", "spatial_melanoma_validation")
out_fig_dir <- file.path(project_dir, "results", "diagnostics", "12_Primary_Visium", "01_build_object")
out_rds_dir <- file.path(project_dir, "data_processed")
log_dir <- file.path(project_dir, "logs")

dir.create(out_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

message("Project directory: ", project_dir)
message("Raw Visium directory: ", raw_spatial_dir)
message("GMT file: ", gmt_file)

if (!dir.exists(raw_spatial_dir)) {
  stop("Raw Visium directory not found: ", raw_spatial_dir)
}
if (!file.exists(gmt_file)) {
  stop("GMT file not found: ", gmt_file)
}

############################################################
## 1. Package setup
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(patchwork)
})

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file)
  message("Saved table: ", file)
}

safe_ggsave <- function(filename, plot, width, height, dpi = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, dpi = dpi, limitsize = FALSE)
  message("Saved figure: ", filename)
}

theme_clean <- function(base_size = 12.5, base_family = "sans") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 2.5, color = "black"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size - 1, color = "black"),
      axis.title = ggplot2::element_text(face = "bold", size = base_size, color = "black"),
      axis.title.y = ggplot2::element_text(face = "bold", size = base_size - 1, color = "black"),
      axis.title.x = ggplot2::element_text(face = "bold", size = base_size, color = "black"),
      axis.text = ggplot2::element_text(size = base_size - 1.5, color = "black"),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1.5),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.5),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey50", linewidth = 0.30),
      panel.grid.major = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.30),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.60),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

state_palette <- c(
  "Immune_defective_Cold" = "#4DBBD5",
  "Myeloid_Treg_Immunosuppressive" = "#00A087",
  "Tumor_dedifferentiation_Stromal_remodeling" = "#E64B35",
  "Melanocytic_Differentiation" = "#3C5488"
)

state_label_map <- c(
  "Immune_defective_Cold" = "immune-defective/\ncold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic\ndifferentiation"
)

############################################################
## 2. Input inventory and locate 10x matrix directory
############################################################

raw_files <- list.files(raw_spatial_dir, recursive = TRUE, full.names = TRUE)
input_inventory <- tibble::tibble(
  File = normalizePath(raw_files, winslash = "/", mustWork = FALSE),
  Basename = basename(raw_files),
  IsDir = dir.exists(raw_files),
  SizeMB = ifelse(file.exists(raw_files), round(file.info(raw_files)$size / 1024^2, 3), NA_real_)
)

safe_write_csv(
  input_inventory,
  file.path(out_tab_dir, "Step11A_spatial_raw_input_file_inventory.csv")
)

tar_candidates <- list.files(
  raw_spatial_dir,
  pattern = "filtered_feature_bc_matrix.*\\.tar\\.gz$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(tar_candidates) > 0) {
  extract_dir <- file.path(raw_spatial_dir, "filtered_feature_bc_matrix_extracted")
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
    message("Extracting 10x matrix tar.gz: ", tar_candidates[1])
    utils::untar(tar_candidates[1], exdir = extract_dir)
  } else {
    message("Using existing extracted matrix directory: ", extract_dir)
  }
}

all_dirs <- unique(c(
  raw_spatial_dir,
  dirname(list.files(raw_spatial_dir, recursive = TRUE, full.names = TRUE))
))

has_matrix_file <- function(d) {
  any(file.exists(file.path(d, c("matrix.mtx", "matrix.mtx.gz")))) &&
    any(file.exists(file.path(d, c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")))) &&
    any(file.exists(file.path(d, c("barcodes.tsv", "barcodes.tsv.gz"))))
}

matrix_dirs <- all_dirs[vapply(all_dirs, has_matrix_file, logical(1))]

matrix_audit <- tibble::tibble(
  MatrixDir = normalizePath(matrix_dirs, winslash = "/", mustWork = FALSE),
  HasMatrix = vapply(matrix_dirs, function(d) any(file.exists(file.path(d, c("matrix.mtx", "matrix.mtx.gz")))), logical(1)),
  HasFeatures = vapply(matrix_dirs, function(d) any(file.exists(file.path(d, c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")))), logical(1)),
  HasBarcodes = vapply(matrix_dirs, function(d) any(file.exists(file.path(d, c("barcodes.tsv", "barcodes.tsv.gz")))), logical(1))
)

safe_write_csv(
  matrix_audit,
  file.path(out_tab_dir, "Step11A_spatial_matrix_directory_audit.csv")
)

if (length(matrix_dirs) < 1) {
  stop(
    "No valid 10x matrix directory found. Expected matrix.mtx(.gz), features/genes.tsv(.gz), and barcodes.tsv(.gz). ",
    "Inspect Step11A_spatial_raw_input_file_inventory.csv."
  )
}

## Prefer extracted filtered_feature_bc_matrix if present
matrix_dir <- matrix_dirs[order(
  !grepl("filtered_feature_bc_matrix", matrix_dirs, ignore.case = TRUE),
  nchar(matrix_dirs)
)][1]

message("Using 10x matrix directory: ", matrix_dir)

############################################################
## 3. Read matrix and create Seurat object
############################################################

counts0 <- Seurat::Read10X(data.dir = matrix_dir)

if (is.list(counts0)) {
  if ("Gene Expression" %in% names(counts0)) {
    counts <- counts0[["Gene Expression"]]
  } else {
    counts <- counts0[[1]]
  }
} else {
  counts <- counts0
}

if (is.null(rownames(counts)) || is.null(colnames(counts))) {
  stop("Counts matrix has no rownames or colnames.")
}

## Remove duplicated gene symbols if present, preserving sparse matrix structure
rownames(counts) <- make.unique(rownames(counts))

seu <- Seurat::CreateSeuratObject(
  counts = counts,
  assay = "Spatial",
  project = "spatial_10x_human_melanoma_IF_FFPE",
  min.cells = 0,
  min.features = 0
)

DefaultAssay(seu) <- "Spatial"

message("Created Seurat object with ", ncol(seu), " spots and ", nrow(seu), " genes.")

############################################################
## 4. Read and attach spatial coordinates
############################################################

coord_candidates <- list.files(
  file.path(raw_spatial_dir, "spatial"),
  pattern = "tissue_positions.*\\.csv$|tissue_positions_list.*\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(coord_candidates) < 1) {
  coord_candidates <- list.files(
    raw_spatial_dir,
    pattern = "tissue_positions.*\\.csv$|tissue_positions_list.*\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
}

coord_audit0 <- tibble::tibble(
  CoordinateFile = normalizePath(coord_candidates, winslash = "/", mustWork = FALSE),
  Basename = basename(coord_candidates),
  SizeMB = ifelse(file.exists(coord_candidates), round(file.info(coord_candidates)$size / 1024^2, 3), NA_real_)
)

safe_write_csv(
  coord_audit0,
  file.path(out_tab_dir, "Step11A_spatial_coordinate_file_candidate_audit.csv")
)

if (length(coord_candidates) < 1) {
  stop("No tissue_positions.csv or tissue_positions_list.csv found under raw spatial directory.")
}

coord_file <- coord_candidates[1]
message("Using coordinate file: ", coord_file)

coord_raw <- readr::read_csv(coord_file, show_col_types = FALSE, col_names = FALSE)

## Detect whether first row is header
first_cell <- as.character(coord_raw[[1]][1])
has_header <- grepl("barcode|Barcode", first_cell)

if (has_header) {
  coord <- readr::read_csv(coord_file, show_col_types = FALSE)
  cn <- colnames(coord)
  ## Harmonize 10x v2 names
  cn_low <- tolower(cn)
  names(cn_low) <- cn
  barcode_col <- cn[cn_low %in% c("barcode", "barcodes")][1]
  in_tissue_col <- cn[cn_low %in% c("in_tissue", "intissue")][1]
  array_row_col <- cn[cn_low %in% c("array_row", "arrayrow")][1]
  array_col_col <- cn[cn_low %in% c("array_col", "arraycol")][1]
  pxl_row_col <- cn[cn_low %in% c("pxl_row_in_fullres", "pxl_row", "imagerow")][1]
  pxl_col_col <- cn[cn_low %in% c("pxl_col_in_fullres", "pxl_col", "imagecol")][1]
  coord <- coord %>%
    dplyr::rename(
      barcode = dplyr::all_of(barcode_col),
      in_tissue = dplyr::all_of(in_tissue_col),
      array_row = dplyr::all_of(array_row_col),
      array_col = dplyr::all_of(array_col_col),
      pxl_row_in_fullres = dplyr::all_of(pxl_row_col),
      pxl_col_in_fullres = dplyr::all_of(pxl_col_col)
    )
} else {
  coord <- coord_raw
  colnames(coord) <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )[seq_len(ncol(coord))]
}

coord <- coord %>%
  dplyr::mutate(
    barcode = as.character(.data$barcode),
    in_tissue = suppressWarnings(as.integer(.data$in_tissue)),
    array_row = suppressWarnings(as.numeric(.data$array_row)),
    array_col = suppressWarnings(as.numeric(.data$array_col)),
    pxl_row_in_fullres = suppressWarnings(as.numeric(.data$pxl_row_in_fullres)),
    pxl_col_in_fullres = suppressWarnings(as.numeric(.data$pxl_col_in_fullres))
  )

spot_barcodes <- colnames(seu)
overlap <- intersect(spot_barcodes, coord$barcode)

coord_overlap_audit <- tibble::tibble(
  MatrixSpots = length(spot_barcodes),
  CoordinateRows = nrow(coord),
  OverlapSpots = length(overlap),
  MatrixOnly = length(setdiff(spot_barcodes, coord$barcode)),
  CoordinateOnly = length(setdiff(coord$barcode, spot_barcodes)),
  CoordinateSource = "10x_tissue_positions_csv",
  CoordinateFile = normalizePath(coord_file, winslash = "/", mustWork = FALSE)
)

safe_write_csv(
  coord_overlap_audit,
  file.path(out_tab_dir, "Step11A_spatial_coordinate_overlap_audit.csv")
)

if (length(overlap) < 10) {
  stop("Too few overlapping barcodes between matrix and tissue_positions file. Please inspect coordinate overlap audit.")
}

## Subset to overlap and attach coordinates to metadata
seu <- subset(seu, cells = overlap)

coord_meta <- coord %>%
  dplyr::filter(.data$barcode %in% colnames(seu)) %>%
  dplyr::distinct(.data$barcode, .keep_all = TRUE)

rownames(coord_meta) <- coord_meta$barcode
coord_meta <- coord_meta[colnames(seu), , drop = FALSE]

seu$barcode <- coord_meta$barcode
seu$in_tissue <- coord_meta$in_tissue
seu$array_row <- coord_meta$array_row
seu$array_col <- coord_meta$array_col
seu$pxl_row_in_fullres <- coord_meta$pxl_row_in_fullres
seu$pxl_col_in_fullres <- coord_meta$pxl_col_in_fullres

## Add coordinate aliases for downstream Step 11
seu$spatial_x <- seu$pxl_col_in_fullres
seu$spatial_y <- seu$pxl_row_in_fullres
seu$imagecol <- seu$pxl_col_in_fullres
seu$imagerow <- seu$pxl_row_in_fullres
seu$coordinate_source <- "10x_tissue_positions_csv_fullres_pixel"
seu$Sample <- "spatial_10x_human_melanoma_IF_FFPE"

############################################################
## 5. QC metrics and normalization
############################################################

## Mitochondrial percentage
mt_genes <- grep("^MT-", rownames(seu), value = TRUE, ignore.case = FALSE)
if (length(mt_genes) > 0) {
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
} else {
  seu$percent.mt <- 0
}

qc_summary <- seu@meta.data %>%
  tibble::rownames_to_column("Spot") %>%
  dplyr::summarise(
    n_spots = dplyr::n(),
    n_genes = nrow(seu),
    median_nFeature_Spatial = median(.data$nFeature_Spatial, na.rm = TRUE),
    median_nCount_Spatial = median(.data$nCount_Spatial, na.rm = TRUE),
    median_percent_mt = median(.data$percent.mt, na.rm = TRUE),
    min_spatial_x = min(.data$spatial_x, na.rm = TRUE),
    max_spatial_x = max(.data$spatial_x, na.rm = TRUE),
    min_spatial_y = min(.data$spatial_y, na.rm = TRUE),
    max_spatial_y = max(.data$spatial_y, na.rm = TRUE)
  )

safe_write_csv(
  qc_summary,
  file.path(out_tab_dir, "Step11A_spatial_object_QC_summary.csv")
)

seu <- Seurat::NormalizeData(
  seu,
  assay = "Spatial",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

############################################################
## 6. Read GMT and compute gene-set scores
############################################################

read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  out <- lapply(lines, function(z) {
    sp <- strsplit(z, "\t", fixed = TRUE)[[1]]
    list(
      name = sp[1],
      desc = ifelse(length(sp) >= 2, sp[2], ""),
      genes = unique(sp[-c(1, 2)])
    )
  })
  names(out) <- vapply(out, `[[`, character(1), "name")
  out
}

gmt <- read_gmt(gmt_file)
message("Imported gene sets: ", length(gmt))

get_assay_matrix <- function(obj, assay = "Spatial", layer = "data") {
  ## Seurat v5 layer-compatible; falls back to slot for older versions
  mat <- tryCatch(
    {
      Seurat::GetAssayData(obj, assay = assay, layer = layer)
    },
    error = function(e) {
      Seurat::GetAssayData(obj, assay = assay, slot = layer)
    }
  )
  mat
}

data_mat <- get_assay_matrix(seu, assay = "Spatial", layer = "data")
gene_universe <- rownames(data_mat)

gene_presence <- lapply(names(gmt), function(gs) {
  genes <- unique(gmt[[gs]]$genes)
  present <- intersect(genes, gene_universe)
  tibble::tibble(
    GeneSet = gs,
    n_total = length(genes),
    n_present = length(present),
    present_fraction = ifelse(length(genes) > 0, length(present) / length(genes), NA_real_),
    Present_genes = paste(present, collapse = ";"),
    Missing_genes = paste(setdiff(genes, gene_universe), collapse = ";")
  )
}) %>%
  dplyr::bind_rows()

safe_write_csv(
  gene_presence,
  file.path(out_tab_dir, "Step11A_spatial_gene_set_presence_summary.csv")
)

score_one_gene_set <- function(mat, genes) {
  present <- intersect(unique(genes), rownames(mat))
  if (length(present) < 3) {
    return(rep(NA_real_, ncol(mat)))
  }
  sub <- as.matrix(mat[present, , drop = FALSE])
  ## Gene-wise z-score across spots; robust to constant genes
  z <- t(scale(t(sub)))
  z[!is.finite(z)] <- 0
  colMeans(z, na.rm = TRUE)
}

gene_set_scores <- lapply(names(gmt), function(gs) {
  score_one_gene_set(data_mat, gmt[[gs]]$genes)
})
gene_set_scores <- as.data.frame(gene_set_scores, check.names = FALSE)
colnames(gene_set_scores) <- names(gmt)
rownames(gene_set_scores) <- colnames(seu)

## Attach individual GMT score columns with a prefix to avoid ambiguity
for (nm in colnames(gene_set_scores)) {
  seu[[paste0("GMT_", nm)]] <- gene_set_scores[[nm]]
}

safe_write_csv(
  tibble::rownames_to_column(gene_set_scores, "Spot"),
  file.path(out_tab_dir, "Step11A_spatial_raw_GMT_gene_set_scores.csv")
)

############################################################
## 7. Collapse 13 GMT signatures into four final state scores
############################################################

final_state_prefixes <- list(
  Immune_defective_Cold = "^Immune_defective_Cold_",
  Myeloid_Treg_Immunosuppressive = "^Myeloid_Treg_Immunosuppressive_",
  Tumor_dedifferentiation_Stromal_remodeling = "^Tumor_dedifferentiation_Stromal_remodeling_",
  Melanocytic_Differentiation = "^Melanocytic_Differentiation_"
)

mapping <- lapply(names(final_state_prefixes), function(state) {
  matched <- grep(final_state_prefixes[[state]], colnames(gene_set_scores), value = TRUE)
  tibble::tibble(
    FinalState = state,
    MatchedGeneSet = matched,
    n_matched = length(matched)
  )
}) %>% dplyr::bind_rows()

safe_write_csv(
  mapping,
  file.path(out_tab_dir, "Step11A_spatial_final_state_signature_mapping.csv")
)

for (state in names(final_state_prefixes)) {
  matched <- grep(final_state_prefixes[[state]], colnames(gene_set_scores), value = TRUE)
  if (length(matched) < 1) {
    stop("No matched GMT signature columns for final state: ", state)
  }
  tmp <- as.matrix(gene_set_scores[, matched, drop = FALSE])
  if (ncol(tmp) == 1) {
    final_score <- as.numeric(scale(tmp[, 1]))
  } else {
    tmp_z <- scale(tmp)
    tmp_z[!is.finite(tmp_z)] <- 0
    final_score <- rowMeans(tmp_z, na.rm = TRUE)
  }
  seu[[state]] <- final_score
}

final_state_cols <- names(final_state_prefixes)

final_scores <- seu@meta.data %>%
  tibble::rownames_to_column("Spot") %>%
  dplyr::select(
    .data$Spot,
    .data$Sample,
    .data$nFeature_Spatial,
    .data$nCount_Spatial,
    .data$percent.mt,
    .data$spatial_x,
    .data$spatial_y,
    .data$imagecol,
    .data$imagerow,
    .data$coordinate_source,
    dplyr::all_of(final_state_cols)
  )

final_scores$Dominant_state <- final_state_cols[
  max.col(as.matrix(final_scores[, final_state_cols, drop = FALSE]), ties.method = "first")
]

seu$Dominant_state <- final_scores$Dominant_state

safe_write_csv(
  final_scores,
  file.path(out_tab_dir, "Step11A_spatial_final_four_ICB_state_scores.csv")
)

dominant_counts <- final_scores %>%
  dplyr::count(.data$Dominant_state, name = "n_spots") %>%
  dplyr::arrange(dplyr::desc(.data$n_spots))

safe_write_csv(
  dominant_counts,
  file.path(out_tab_dir, "Step11A_spatial_dominant_state_counts.csv")
)

############################################################
## 8. Spatial coordinate preview
############################################################

plot_df <- final_scores %>%
  dplyr::mutate(
    spatial_y_plot = -.data$spatial_y,
    Dominant_state = factor(.data$Dominant_state, levels = final_state_cols)
  )

p_coord <- ggplot(plot_df, aes(x = spatial_x, y = spatial_y_plot)) +
  geom_point(aes(color = Dominant_state), size = 1.0, alpha = 0.9) +
  scale_color_manual(
    values = state_palette,
    breaks = final_state_cols,
    labels = state_label_map[final_state_cols],
    drop = FALSE
  ) +
  coord_equal() +
  labs(
    title = "Visium dominant-state spatial map",
    subtitle = "Coordinates from 10x tissue_positions.csv; y-axis inverted for display",
    x = "Full-resolution pixel x",
    y = "Full-resolution pixel y",
    color = "Dominant state"
  ) +
  theme_clean(base_size = 12.5)

safe_ggsave(
  file.path(out_fig_dir, "Step11A_spatial_coordinate_preview_state_scored.png"),
  p_coord,
  width = 8.5,
  height = 6.5,
  dpi = 300
)
safe_ggsave(
  file.path(out_fig_dir, "Step11A_spatial_coordinate_preview_state_scored.pdf"),
  p_coord,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

p_qc <- ggplot(seu@meta.data, aes(x = nCount_Spatial, y = nFeature_Spatial)) +
  geom_point(alpha = 0.75, size = 1.0, color = "#4C72B0") +
  labs(
    title = "Visium spatial QC",
    subtitle = "RNA depth versus detected spatial features",
    x = "Spatial RNA counts (nCount_Spatial)",
    y = "Detected spatial features (nFeature_Spatial)"
  ) +
  theme_clean(base_size = 12.5)

safe_ggsave(
  file.path(out_fig_dir, "Step11A_spatial_QC_nCount_vs_nFeature.png"),
  p_qc,
  width = 6.0,
  height = 4.8,
  dpi = 300
)
safe_ggsave(
  file.path(out_fig_dir, "Step11A_spatial_QC_nCount_vs_nFeature.pdf"),
  p_qc,
  width = 6.0,
  height = 4.8,
  dpi = 300
)

############################################################
## 9. Save state-scored spatial object
############################################################

out_rds <- file.path(
  out_rds_dir,
  "spatial_10x_human_melanoma_IF_FFPE_with_ICBcomb_state_scores.rds"
)

saveRDS(seu, out_rds)
message("Saved RDS: ", out_rds)

object_audit <- tibble::tibble(
  OutputRDS = normalizePath(out_rds, winslash = "/", mustWork = FALSE),
  ObjectClass = paste(class(seu), collapse = ";"),
  Assay = DefaultAssay(seu),
  n_spots = ncol(seu),
  n_genes = nrow(seu),
  HasSpatialX = "spatial_x" %in% colnames(seu@meta.data),
  HasSpatialY = "spatial_y" %in% colnames(seu@meta.data),
  CoordinateSource = unique(seu$coordinate_source)[1],
  HasAllFinalStateColumns = all(final_state_cols %in% colnames(seu@meta.data)),
  FinalStateColumns = paste(final_state_cols, collapse = ";"),
  OutputTime = as.character(Sys.time())
)

safe_write_csv(
  object_audit,
  file.path(out_tab_dir, "Step11A_spatial_state_scored_object_audit.csv")
)

output_inventory <- tibble::tibble(
  OutputType = c("RDS", "Table", "Table", "Table", "Figure", "Figure"),
  Description = c(
    "State-scored Visium Seurat object",
    "Raw input file inventory",
    "Coordinate overlap audit",
    "Final four state scores",
    "Coordinate preview",
    "QC scatter preview"
  ),
  Path = c(
    out_rds,
    file.path(out_tab_dir, "Step11A_spatial_raw_input_file_inventory.csv"),
    file.path(out_tab_dir, "Step11A_spatial_coordinate_overlap_audit.csv"),
    file.path(out_tab_dir, "Step11A_spatial_final_four_ICB_state_scores.csv"),
    file.path(out_fig_dir, "Step11A_spatial_coordinate_preview_state_scored.png"),
    file.path(out_fig_dir, "Step11A_spatial_QC_nCount_vs_nFeature.png")
  )
) %>%
  dplyr::mutate(
    Exists = file.exists(.data$Path),
    Path = normalizePath(.data$Path, winslash = "/", mustWork = FALSE)
  )

safe_write_csv(
  output_inventory,
  file.path(out_tab_dir, "Step11A_spatial_build_and_scoring_output_inventory.csv")
)

sink(file.path(log_dir, "Step11A_spatial_build_and_state_scoring_sessionInfo.txt"))
print(sessionInfo())
sink()

message("\nStep 11A completed successfully.")
message("Main output RDS: ", out_rds)
message("Next step:")
message("  Run Step 11 v3 script again, or set manual_spatial_rds to this RDS path.")

## Public sequential end gate
.step01_out <- file.path(out_tab_dir, "Step11A_spatial_final_four_ICB_state_scores.csv")
if (!file.exists(.step01_out)) stop("STEP 01 output missing: ", .step01_out, call. = FALSE)
.step01_n <- nrow(utils::read.csv(.step01_out, check.names = FALSE))
cat("\nSTEP 01 PASS — raw Visium -> state-scored object; rows = ", .step01_n, "\n", sep = "")
