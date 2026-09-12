############################################################
## STYLE-ONLY v1.1 S20-REFERENCE LOCK
## Figure typography is harmonized to the accepted S20 visual ratio.
## Long subtitles are locally reduced and/or wrapped when needed.
## Scientific inputs, thresholds, statistics, state definitions,
## analytical logic, numerical outputs, and figure canvas sizes are unchanged.
############################################################

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

############################################################
## 11_spatial_QC_Moran_colocalization_CLEAN_INPUT_v7_FIG8A_SUBTITLEFIX.R
## Project: ICB_resistance_project
## Purpose:
##   Clean-room Visium spatial validation of mechanism-defined
##   tumor-immune states:
##   1) input/object/coordinate audit
##   2) state-QC Spearman correlation
##   3) QC residualization and raw-vs-residual concordance
##   4) kNN spatial graph construction
##   5) Moran's I for raw and QC-residualized state scores
##   6) high-state neighbor enrichment
##   7) Myeloid-Treg + Dediff/Stromal dual-high colocalization
##   8) top 20/25/30% and pooled/per-sample threshold sensitivity
##
## Recommended run:
##   source("D:/ICB_resistance_project/scripts/11_spatial_QC_Moran_colocalization_CLEAN_INPUT_v7_FIG8A_SUBTITLEFIX.R")
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

############################################################
## 0. Project configuration
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir,winslash="/",mustWork=TRUE)
if (!dir.exists(project_dir)) {
  stop("Project directory not found: ", project_dir)
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

## Parameters
KNN_K <- 6
MORAN_N_PERM <- 999
MAIN_HIGH_TOP_FRACTION <- 0.25
SENSITIVITY_TOP_FRACTIONS <- c(0.20, 0.25, 0.30)
RANDOM_SEED <- 20260503
set.seed(RANDOM_SEED)

## Optional manual override. Leave as NA_character_ for automatic discovery.
## If automatic discovery fails, set this to the true Visium/spatial Seurat RDS path.
default_spatial_rds <- file.path(project_dir, "data_processed", "spatial_10x_human_melanoma_IF_FFPE_with_ICBcomb_state_scores.rds")
manual_spatial_rds <- if (file.exists(default_spatial_rds)) default_spatial_rds else NA_character_

## Output directories
out_table_dir <- file.path(project_dir, "results", "tables", "spatial_melanoma_validation")
out_fig_dir   <- file.path(project_dir, "results", "figures", "spatial_melanoma_validation")
out_int_dir   <- file.path(project_dir, "results", "intermediate", "spatial_melanoma_validation")
out_log_dir   <- file.path(project_dir, "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
diagnostic_fig_dir <- file.path(
  project_dir, "results", "diagnostics", "12_Primary_Visium", "02_spatial_QC_Moran"
)
dir.create(diagnostic_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_int_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

message("Project directory: ", project_dir)
message("Output table directory: ", out_table_dir)
message("Output figure directory: ", out_fig_dir)

############################################################
## 1. Packages and helpers
############################################################

cran_pkgs <- c("dplyr", "tidyr", "ggplot2", "patchwork", "FNN", "Matrix", "tibble", "stringr")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(FNN)
  library(Matrix)
  library(tibble)
  library(stringr)
})

## Load Seurat only if available. The script can inspect Seurat objects if package is present.
if (!requireNamespace("Seurat", quietly = TRUE)) {
  warning("Seurat is not installed. The script can still read plain data frames, but Seurat RDS extraction may fail.")
} else {
  suppressPackageStartupMessages(library(Seurat))
}

## Avoid function masking.
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
left_join <- dplyr::left_join
bind_rows <- dplyr::bind_rows
count <- dplyr::count
all_of <- dplyr::all_of
any_of <- dplyr::any_of

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file = file, row.names = row.names, fileEncoding = "UTF-8")
  message("Saved table: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

safe_ggsave <- function(file, plot, width, height, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )
  message("Saved figure: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

theme_icb <- function(base_size = 12.5, base_family = "sans") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 2.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base_size - 1),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.title.y = ggplot2::element_text(face = "bold", size = base_size - 1),
      axis.title.x = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1.5),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1.5),
      strip.background = ggplot2::element_rect(fill = "#D9D9D9", color = "black", linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.60),
      legend.key = ggplot2::element_rect(fill = "white", color = NA)
    )
}

state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_labels <- c(
  Immune_defective_Cold = "Immune-defective/Cold",
  Myeloid_Treg_Immunosuppressive = "Myeloid-Treg Immunosuppressive",
  Tumor_dedifferentiation_Stromal_remodeling = "Tumor-dedifferentiation/Stromal-remodeling",
  Melanocytic_Differentiation = "Melanocytic Differentiation"
)

state_colors <- c(
  "Immune_defective_Cold" = "#4DBBD5",
  "Myeloid_Treg_Immunosuppressive" = "#00A087",
  "Tumor_dedifferentiation_Stromal_remodeling" = "#E64B35",
  "Melanocytic_Differentiation" = "#3C5488"
)

state_display_public <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

state_label_to_public <- setNames(
  unname(state_display_public[state_cols]),
  unname(state_labels[state_cols])
)
plot_state_label <- function(x) {
  y <- as.character(x)
  hit <- y %in% names(state_label_to_public)
  y[hit] <- unname(state_label_to_public[y[hit]])
  y
}

heatmap_low  <- "#3B82F6"
heatmap_mid  <- "white"
heatmap_high <- "#EF4444"

focused_colors <- c(
  "Neither high" = "#D9D9D9",
  "Myeloid-Treg high only" = "#00A087",
  "Dediff/Stromal high only" = "#E64B35",
  "Both high" = "#7E2F8E"
)

############################################################
## 2. Input discovery and object loading
############################################################

## The script prioritizes locked/clean spatial RDS objects with state-score outputs.
rds_candidates <- list.files(
  project_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

## Prefer spatial melanoma validation objects; avoid final figure collectors and old archive if possible.
score_pattern <- paste(c(
  "spatial", "mel", "melanoma", "Visium", "10x", "state", "ICBcomb", "with_spatial_state_outputs", "with_ICBcomb_state_scores"
), collapse = "|")

rds_rank <- data.frame(
  CandidateFile = rds_candidates,
  FileName = basename(rds_candidates),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    lower_path = tolower(CandidateFile),
    score = 0,
    score = score + ifelse(grepl("spatial|visium|10x", lower_path), 6, 0),
    score = score + ifelse(grepl("melanoma|mel_sp", lower_path), 3, 0),
    score = score + ifelse(grepl("with_spatial_state_outputs|with_icbcomb_state_scores|spatial_state", lower_path), 6, 0),
    score = score + ifelse(grepl("state|icbcomb", lower_path), 2, 0),
    score = score + ifelse(grepl("final|locked|clean", lower_path), 2, 0),
    score = score - ifelse(grepl("bulk_discovery|external_validation|final_tumor_immune_state_scores|gse244982_final_tumor_immune_state_scores|gse78220|gse91061|final_main_figures|final_supplementary", lower_path), 30, 0),
    score = score - ifelse(grepl("archive|old|backup", lower_path), 10, 0)
  ) %>%
  dplyr::arrange(desc(score), CandidateFile)

if (!is.na(manual_spatial_rds)) {
  if (!file.exists(manual_spatial_rds)) stop("manual_spatial_rds does not exist: ", manual_spatial_rds)
  manual_spatial_rds_norm <- normalizePath(manual_spatial_rds, winslash = "/", mustWork = TRUE)
  rds_rank <- rbind(
    data.frame(
      CandidateFile = manual_spatial_rds_norm,
      FileName = basename(manual_spatial_rds_norm),
      lower_path = tolower(manual_spatial_rds_norm),
      score = 999,
      stringsAsFactors = FALSE
    ),
    rds_rank
  )
}

safe_write_csv(rds_rank, file.path(out_table_dir, "Step11_spatial_input_rds_candidate_inventory.csv"))

## Try to select the first RDS that contains required state columns or aliases.
normalize_col <- function(x) tolower(gsub("[^a-zA-Z0-9]", "", x))
state_aliases <- list(
  Immune_defective_Cold = c("Immune_defective_Cold", "Immune_defective_cold", "Immune_defective_Cold_score", "Immune_defective_cold_score"),
  Myeloid_Treg_Immunosuppressive = c("Myeloid_Treg_Immunosuppressive", "Myeloid_Treg_immunosuppressive", "Myeloid_Treg_Immunosuppressive_score", "Myeloid_Treg_immunosuppressive_score"),
  Tumor_dedifferentiation_Stromal_remodeling = c("Tumor_dedifferentiation_Stromal_remodeling", "Tumor_dedifferentiated", "Tumor_dedifferentiation_Stromal_remodeling_score", "Tumor_dedifferentiated_score"),
  Melanocytic_Differentiation = c("Melanocytic_Differentiation", "Melanocytic_differentiated", "Melanocytic_Differentiation_score", "Melanocytic_differentiated_score")
)

find_state_aliases <- function(cols) {
  ncols <- normalize_col(cols)
  out <- setNames(rep(NA_character_, length(state_cols)), state_cols)
  for (st in state_cols) {
    aliases <- normalize_col(state_aliases[[st]])
    idx <- match(aliases, ncols)
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) out[[st]] <- cols[idx[1]]
  }
  out
}

read_spatial_candidate <- function(file) {
  obj <- readRDS(file)
  meta <- NULL
  object_class <- paste(class(obj), collapse = ";")
  if (inherits(obj, "Seurat")) {
    meta <- obj@meta.data
  } else if (is.data.frame(obj)) {
    meta <- obj
  } else if (is.list(obj)) {
    if (!is.null(obj$metadata) && is.data.frame(obj$metadata)) meta <- obj$metadata
    if (is.null(meta) && !is.null(obj$meta.data) && is.data.frame(obj$meta.data)) meta <- obj$meta.data
  }
  list(obj = obj, meta = meta, object_class = object_class)
}

has_candidate_coordinates <- function(obj, meta) {
  if (is.null(meta) || !is.data.frame(meta)) return(FALSE)
  meta_cols <- colnames(meta)
  pairs <- list(
    c("global_x", "global_y"),
    c("Global.X", "Global.Y"),
    c("x_global", "y_global"),
    c("X_global", "Y_global"),
    c("x", "y"),
    c("X", "Y"),
    c("imagecol", "imagerow"),
    c("array_col", "array_row"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("spatial_1", "spatial_2")
  )
  for (pair in pairs) {
    if (all(pair %in% meta_cols)) {
      xx <- suppressWarnings(as.numeric(meta[[pair[1]]]))
      yy <- suppressWarnings(as.numeric(meta[[pair[2]]]))
      if (sum(is.finite(xx) & is.finite(yy)) > 10) return(TRUE)
    }
  }
  if (inherits(obj, "Seurat")) {
    if (length(obj@images) > 0) return(TRUE)
  }
  FALSE
}

## Candidate screening is stricter than state-column matching.
## A bulk state-score table may contain all four state columns but has no spatial coordinates.
## Therefore a candidate must have all four state columns AND spatial coordinates/Seurat images.
selected <- NULL
selected_info <- list()
candidate_screen <- list()
if (nrow(rds_rank) > 0) {
  for (f in rds_rank$CandidateFile) {
    tmp <- tryCatch(read_spatial_candidate(f), error = function(e) NULL)
    if (is.null(tmp) || is.null(tmp$meta)) {
      candidate_screen[[length(candidate_screen) + 1]] <- data.frame(
        CandidateFile = f,
        ObjectClass = NA_character_,
        N_meta_rows = NA_integer_,
        N_meta_cols = NA_integer_,
        StateColumnsFoundN = 0,
        HasCandidateCoordinates = FALSE,
        NonSpatialBulkLikeExcluded = grepl("final_tumor_immune_state_scores|bulk_discovery|external_validation|gse78220|gse91061", tolower(f)),
        Selected = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }
    aliases <- find_state_aliases(colnames(tmp$meta))
    n_found <- sum(!is.na(aliases))
    has_coord <- has_candidate_coordinates(tmp$obj, tmp$meta)
    is_nonspatial_bulk_like <- nrow(tmp$meta) < 100 || grepl("final_tumor_immune_state_scores|bulk_discovery|external_validation|gse78220|gse91061", tolower(f))
    selectable <- (n_found >= 4) && has_coord && !is_nonspatial_bulk_like
    candidate_screen[[length(candidate_screen) + 1]] <- data.frame(
      CandidateFile = f,
      ObjectClass = tmp$object_class,
      N_meta_rows = nrow(tmp$meta),
      N_meta_cols = ncol(tmp$meta),
      StateColumnsFoundN = n_found,
      HasCandidateCoordinates = has_coord,
      NonSpatialBulkLikeExcluded = is_nonspatial_bulk_like,
      Selected = selectable,
      stringsAsFactors = FALSE
    )
    if (selectable) {
      selected <- f
      selected_info <- tmp
      selected_info$aliases <- aliases
      break
    }
  }
}

candidate_screen_df <- dplyr::bind_rows(candidate_screen)
safe_write_csv(candidate_screen_df, file.path(out_table_dir, "Step11_spatial_rds_candidate_screening_audit.csv"))

if (is.null(selected)) {
  stop(
    "No usable SPATIAL RDS containing all four state-score columns AND spatial coordinates was found.
",
    "The previous failure happened because a bulk state-score table can contain all four state columns but no coordinates.
",
    "Please inspect: ", file.path(out_table_dir, "Step11_spatial_rds_candidate_screening_audit.csv"), "
",
    "If the true Visium RDS exists but is not selected, set manual_spatial_rds near the top of this script to its full path."
  )
}

sp_obj <- selected_info$obj
meta <- selected_info$meta
state_alias_map <- selected_info$aliases

input_audit <- data.frame(
  SelectedFile = selected,
  ObjectClass = selected_info$object_class,
  N_meta_rows = nrow(meta),
  N_meta_cols = ncol(meta),
  UnitOfAnalysis = "spot",
  KNN_k = KNN_K,
  Moran_n_perm = MORAN_N_PERM,
  MainHighTopFraction = MAIN_HIGH_TOP_FRACTION,
  stringsAsFactors = FALSE
)
safe_write_csv(input_audit, file.path(out_table_dir, "Step11_spatial_selected_input_object_audit.csv"))

safe_write_csv(
  data.frame(StandardState = names(state_alias_map), MatchedColumn = unname(state_alias_map), stringsAsFactors = FALSE),
  file.path(out_table_dir, "Step11_spatial_state_column_alias_audit.csv")
)

message("Selected spatial object: ", selected)
message("Object class: ", selected_info$object_class)
message("Meta rows: ", nrow(meta), "; columns: ", ncol(meta))
message("State aliases: ", paste(names(state_alias_map), unname(state_alias_map), sep = "=", collapse = "; "))

############################################################
## 3. Standardize metadata, QC metrics and state columns
############################################################

meta <- as.data.frame(meta)
if (is.null(rownames(meta)) || any(rownames(meta) == "")) {
  meta$SpotID <- paste0("spot_", seq_len(nrow(meta)))
} else {
  meta$SpotID <- rownames(meta)
}

## Create standard state-score columns.
for (st in state_cols) {
  meta[[st]] <- suppressWarnings(as.numeric(meta[[state_alias_map[[st]]]]))
}

## QC metrics. Prefer Spatial names, otherwise detect common alternatives.
find_first_col <- function(cols, patterns) {
  for (pat in patterns) {
    hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

qc_alias <- c(
  nFeature = find_first_col(colnames(meta), c("^nFeature_Spatial$", "nFeature", "detected", "nGenes", "n_gene")),
  nCount   = find_first_col(colnames(meta), c("^nCount_Spatial$", "nCount", "UMI", "total_counts", "n_umi")),
  percent_mt = find_first_col(colnames(meta), c("^percent\\.mt$", "percent_mt", "pct.*mt", "mito"))
)

## Use explicit if/else blocks here to avoid parsing issues in some R sessions.
if (is.na(qc_alias[["nFeature"]])) {
  meta$nFeature_Spatial_audit <- NA_real_
} else {
  meta$nFeature_Spatial_audit <- suppressWarnings(
    as.numeric(meta[[qc_alias[["nFeature"]]]])
  )
}

if (is.na(qc_alias[["nCount"]])) {
  meta$nCount_Spatial_audit <- NA_real_
} else {
  meta$nCount_Spatial_audit <- suppressWarnings(
    as.numeric(meta[[qc_alias[["nCount"]]]])
  )
}

if (is.na(qc_alias[["percent_mt"]])) {
  meta$percent.mt_audit <- NA_real_
} else {
  meta$percent.mt_audit <- suppressWarnings(
    as.numeric(meta[[qc_alias[["percent_mt"]]]])
  )
}

## If percent.mt is unavailable and a Seurat assay is accessible, try to compute it.
if (all(is.na(meta$percent.mt_audit)) && inherits(sp_obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
  try({
    gene_names <- rownames(sp_obj)
    mt_genes <- grep("^MT-|^mt-", gene_names, value = TRUE)
    if (length(mt_genes) > 0) {
      pct <- Seurat::PercentageFeatureSet(sp_obj, features = mt_genes)
      pct <- pct[rownames(meta)]
      meta$percent.mt_audit <- as.numeric(pct)
    }
  }, silent = TRUE)
}

qc_cols <- c("nFeature_Spatial_audit", "nCount_Spatial_audit", "percent.mt_audit")
qc_audit <- data.frame(
  StandardQC = c("nFeature", "nCount", "percent.mt"),
  MatchedColumn = unname(qc_alias),
  AuditColumn = qc_cols,
  NonMissingN = vapply(qc_cols, function(x) sum(!is.na(meta[[x]])), numeric(1)),
  UniqueN = vapply(qc_cols, function(x) length(unique(na.omit(meta[[x]]))), numeric(1)),
  stringsAsFactors = FALSE
)
safe_write_csv(qc_audit, file.path(out_table_dir, "Step11_spatial_QC_metric_alias_audit.csv"))

############################################################
## 4. Coordinate extraction and audit
############################################################

extract_coordinates <- function(obj, meta) {
  ## 4A. Common metadata coordinate columns.
  meta_cols <- colnames(meta)
  candidates <- list(
    c("spatial_x", "spatial_y"),
    c("global_x", "global_y"),
    c("Global.X", "Global.Y"),
    c("x_global", "y_global"),
    c("X_global", "Y_global"),
    c("x", "y"),
    c("X", "Y"),
    c("imagecol", "imagerow"),
    c("array_col", "array_row"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("spatial_1", "spatial_2"),
    c("umap_1", "umap_2")
  )
  for (pair in candidates) {
    if (all(pair %in% meta_cols)) {
      coor <- data.frame(
        SpotID = meta$SpotID,
        spatial_x = suppressWarnings(as.numeric(meta[[pair[1]]])),
        spatial_y = suppressWarnings(as.numeric(meta[[pair[2]]])),
        coordinate_source = paste0("metadata:", pair[1], "+", pair[2]),
        stringsAsFactors = FALSE
      )
      if (sum(is.finite(coor$spatial_x) & is.finite(coor$spatial_y)) > 10) return(coor)
    }
  }

  ## 4B. Seurat image coordinates.
  if (inherits(obj, "Seurat") && requireNamespace("Seurat", quietly = TRUE)) {
    image_names <- names(obj@images)
    if (length(image_names) > 0) {
      for (img in image_names) {
        coor <- tryCatch({
          x <- Seurat::GetTissueCoordinates(obj, image = img)
          x <- as.data.frame(x)
          x$SpotID <- rownames(x)
          cn <- colnames(x)
          ## Seurat may return x/y, imagecol/imagerow, or col/row.
          xcol <- find_first_col(cn, c("^x$", "imagecol", "col$", "pxl_col"))
          ycol <- find_first_col(cn, c("^y$", "imagerow", "row$", "pxl_row"))
          if (is.na(xcol) || is.na(ycol)) stop("No x/y coordinate columns returned by GetTissueCoordinates")
          data.frame(
            SpotID = x$SpotID,
            spatial_x = suppressWarnings(as.numeric(x[[xcol]])),
            spatial_y = suppressWarnings(as.numeric(x[[ycol]])),
            coordinate_source = paste0("SeuratImage:", img, ":", xcol, "+", ycol),
            stringsAsFactors = FALSE
          )
        }, error = function(e) NULL)
        if (!is.null(coor) && sum(is.finite(coor$spatial_x) & is.finite(coor$spatial_y)) > 10) return(coor)
      }
    }
  }

  NULL
}

coord_df <- extract_coordinates(sp_obj, meta)
if (is.null(coord_df)) {
  stop("No usable spatial coordinates found. Please inspect object metadata and image coordinates before continuing.")
}

## Join coordinates to metadata.
## Deduplicate coordinate rows by spot ID using base R.
## This avoids conflicts with Bioconductor/IRanges::slice().
coord_df <- coord_df[!duplicated(coord_df$SpotID), , drop = FALSE]

## Avoid duplicated spatial_x/spatial_y columns after joining.
## The Step11A object already stores spatial_x/spatial_y in metadata; joining another
## coordinate table would otherwise create spatial_x.x / spatial_x.y and remove the
## plain spatial_x name, which breaks downstream select/filter calls.
pre_existing_coord_cols <- intersect(c("spatial_x", "spatial_y", "coordinate_source"), colnames(meta))
meta_no_coord <- meta[, setdiff(colnames(meta), pre_existing_coord_cols), drop = FALSE]

meta2 <- dplyr::left_join(meta_no_coord, coord_df, by = "SpotID")

## Defensive standardization in case a future object still creates suffixes.
if (!("spatial_x" %in% colnames(meta2))) {
  sx_candidates <- intersect(c("spatial_x.y", "spatial_x.x", "imagecol", "pxl_col_in_fullres"), colnames(meta2))
  if (length(sx_candidates) > 0) meta2$spatial_x <- suppressWarnings(as.numeric(meta2[[sx_candidates[1]]]))
}
if (!("spatial_y" %in% colnames(meta2))) {
  sy_candidates <- intersect(c("spatial_y.y", "spatial_y.x", "imagerow", "pxl_row_in_fullres"), colnames(meta2))
  if (length(sy_candidates) > 0) meta2$spatial_y <- suppressWarnings(as.numeric(meta2[[sy_candidates[1]]]))
}
if (!("coordinate_source" %in% colnames(meta2))) {
  cs_candidates <- intersect(c("coordinate_source.y", "coordinate_source.x"), colnames(meta2))
  if (length(cs_candidates) > 0) meta2$coordinate_source <- as.character(meta2[[cs_candidates[1]]])
}

if (!("spatial_x" %in% colnames(meta2)) || !("spatial_y" %in% colnames(meta2))) {
  stop("Coordinate merge failed: spatial_x/spatial_y were not available after joining. Please inspect Step11_spatial_coordinates_used.csv and object metadata.")
}

coord_audit <- data.frame(
  n_meta_spots = nrow(meta),
  n_coordinate_rows = nrow(coord_df),
  n_matched_spots = sum(is.finite(meta2$spatial_x) & is.finite(meta2$spatial_y)),
  matched_rate = mean(is.finite(meta2$spatial_x) & is.finite(meta2$spatial_y)),
  coordinate_source = paste(unique(na.omit(coord_df$coordinate_source)), collapse = ";"),
  pre_existing_coord_cols_removed_before_join = paste(pre_existing_coord_cols, collapse = ";"),
  unit_of_analysis = "spot",
  stringsAsFactors = FALSE
)
safe_write_csv(coord_audit, file.path(out_table_dir, "Step11_spatial_coordinate_source_and_merge_audit.csv"))
safe_write_csv(
  meta2[, intersect(c("SpotID", "spatial_x", "spatial_y", "coordinate_source"), colnames(meta2)), drop = FALSE],
  file.path(out_table_dir, "Step11_spatial_coordinates_used.csv")
)

## Drop spots without coordinates or state scores.
analysis_df <- meta2 %>%
  dplyr::filter(is.finite(spatial_x), is.finite(spatial_y)) %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(state_cols), ~ is.finite(.x)))

if (nrow(analysis_df) < 50) {
  stop("Too few spots with both spatial coordinates and state scores: ", nrow(analysis_df))
}

## Sample/slide grouping for graph construction. Avoid cross-sample edges if multiple samples exist.
sample_col <- find_first_col(colnames(analysis_df), c("^orig.ident$", "sample", "slide", "library", "section", "slice", "patient"))
if (!is.na(sample_col)) {
  analysis_df$SpatialSampleID <- as.character(analysis_df[[sample_col]])
} else {
  analysis_df$SpatialSampleID <- "Sample_1"
}
analysis_df$SpatialSampleID[is.na(analysis_df$SpatialSampleID) | analysis_df$SpatialSampleID == ""] <- "Sample_1"

sample_audit <- analysis_df %>%
  dplyr::count(SpatialSampleID, name = "n_spots") %>%
  dplyr::arrange(desc(n_spots))
safe_write_csv(sample_audit, file.path(out_table_dir, "Step11_spatial_sample_group_audit.csv"))

saveRDS(analysis_df, file.path(out_int_dir, "Step11_spatial_analysis_metadata_with_coordinates_and_scores.rds"))
safe_write_csv(
  analysis_df %>% dplyr::select(SpotID, SpatialSampleID, spatial_x, spatial_y, dplyr::all_of(qc_cols), dplyr::all_of(state_cols)),
  file.path(out_table_dir, "Step11_spatial_analysis_metadata_with_coordinates_and_scores.csv")
)

############################################################
## 5. QC-state Spearman correlation and residualization
############################################################

spearman_one <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10 || length(unique(x[ok])) < 3 || length(unique(y[ok])) < 3) {
    return(data.frame(rho = NA_real_, p_value = NA_real_, n = sum(ok)))
  }
  tt <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  data.frame(rho = unname(tt$estimate), p_value = tt$p.value, n = sum(ok))
}

qc_cor <- dplyr::bind_rows(lapply(state_cols, function(st) {
  dplyr::bind_rows(lapply(qc_cols, function(qc) {
    z <- spearman_one(analysis_df[[st]], analysis_df[[qc]])
    data.frame(State = st, StateLabel = state_labels[[st]], QC_metric = qc, z, stringsAsFactors = FALSE)
  }))
})) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH"))
safe_write_csv(qc_cor, file.path(out_table_dir, "Step11_spatial_QC_state_Spearman_correlation.csv"))

## Residualize state scores against available QC covariates.
qc_for_lm <- qc_cols[vapply(qc_cols, function(qc) {
  x <- analysis_df[[qc]]
  sum(is.finite(x)) >= 10 && length(unique(na.omit(x))) >= 3
}, logical(1))]

resid_audit <- data.frame(
  QC_metric = qc_cols,
  UsedForResidualization = qc_cols %in% qc_for_lm,
  NonMissingN = vapply(qc_cols, function(qc) sum(is.finite(analysis_df[[qc]])), numeric(1)),
  UniqueN = vapply(qc_cols, function(qc) length(unique(na.omit(analysis_df[[qc]]))), numeric(1)),
  stringsAsFactors = FALSE
)
safe_write_csv(resid_audit, file.path(out_table_dir, "Step11_spatial_QC_residualization_covariate_audit.csv"))

analysis_resid <- analysis_df
for (st in state_cols) {
  resid_col <- paste0(st, "_QCresid")
  if (length(qc_for_lm) == 0) {
    analysis_resid[[resid_col]] <- scale(analysis_resid[[st]])[, 1]
  } else {
    lm_df <- analysis_resid[, c(st, qc_for_lm), drop = FALSE]
    colnames(lm_df)[1] <- "state_score"
    ok <- complete.cases(lm_df)
    rr <- rep(NA_real_, nrow(analysis_resid))
    if (sum(ok) >= 20) {
      fit <- lm(as.formula(paste("state_score ~", paste(qc_for_lm, collapse = " + "))), data = lm_df[ok, , drop = FALSE])
      rr[ok] <- residuals(fit)
      rr[ok] <- as.numeric(scale(rr[ok]))
    }
    analysis_resid[[resid_col]] <- rr
  }
}
resid_cols <- paste0(state_cols, "_QCresid")

raw_resid_cor <- dplyr::bind_rows(lapply(seq_along(state_cols), function(i) {
  st <- state_cols[i]
  rc <- resid_cols[i]
  z <- spearman_one(analysis_resid[[st]], analysis_resid[[rc]])
  data.frame(State = st, StateLabel = state_labels[[st]], RawColumn = st, ResidualColumn = rc, z, stringsAsFactors = FALSE)
})) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH"))
safe_write_csv(raw_resid_cor, file.path(out_table_dir, "Step11_spatial_raw_vs_QCresidual_state_correlation.csv"))

saveRDS(analysis_resid, file.path(out_int_dir, "Step11_spatial_analysis_metadata_with_QCresidual_scores.rds"))

## QC correlation heatmap and raw-residual concordance plot.
qc_cor_plot_df <- qc_cor %>%
  dplyr::mutate(
    StateLabel = factor(StateLabel, levels = state_labels[state_cols]),
    QC_metric = factor(QC_metric, levels = qc_cols)
  )
p_qc_cor <- ggplot(qc_cor_plot_df, aes(x = QC_metric, y = StateLabel, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(rho), "NA", sprintf("%.2f", rho))), size = 3.4) +
  scale_fill_gradient2(low = heatmap_low, mid = heatmap_mid, high = heatmap_high, midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
  labs(title = "A  State-score correlations with spatial QC metrics", x = NULL, y = NULL) +
  theme_icb(base_size = 11.5) +
  scale_y_discrete(labels = plot_state_label) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25), panel.grid = element_blank())

p_raw_resid <- ggplot(raw_resid_cor, aes(x = StateLabel, y = rho, fill = State)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(is.na(rho), "NA", sprintf("%.2f", rho))), vjust = -0.25, size = 3.5) +
  scale_fill_manual(values = state_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "B  Raw vs QC-residualized state-score concordance", x = NULL, y = "Spearman rho") +
  theme_icb(base_size = 11.5) +
  scale_x_discrete(labels = plot_state_label) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25))

p_qc_combined <- p_qc_cor | p_raw_resid
safe_ggsave(file.path(diagnostic_fig_dir, "Supplementary_Figure_S11_spatial_QC_state_correlation_and_residualization.png"), p_qc_combined, 12.5, 5.5)
safe_ggsave(file.path(diagnostic_fig_dir, "Supplementary_Figure_S11_spatial_QC_state_correlation_and_residualization.pdf"), p_qc_combined, 12.5, 5.5)

############################################################
## 6. Spatial graph construction
############################################################

make_knn_edges <- function(df, k = 6) {
  out <- list()
  idx_global <- seq_len(nrow(df))
  for (sid in unique(df$SpatialSampleID)) {
    ii <- idx_global[df$SpatialSampleID == sid]
    sub <- df[ii, , drop = FALSE]
    if (nrow(sub) <= k + 1) next
    coor <- as.matrix(sub[, c("spatial_x", "spatial_y")])
    knn <- FNN::get.knn(coor, k = k)
    edges <- data.frame(
      from = rep(ii, each = k),
      to = as.vector(t(ii[knn$nn.index])),
      distance = as.vector(t(knn$nn.dist)),
      SpatialSampleID = sid,
      stringsAsFactors = FALSE
    )
    out[[sid]] <- edges
  }
  dplyr::bind_rows(out)
}

edges <- make_knn_edges(analysis_resid, k = KNN_K)
if (nrow(edges) == 0) stop("No kNN edges created. Check coordinates and sample group sizes.")
safe_write_csv(edges, file.path(out_table_dir, "Step11_spatial_knn6_edge_list.csv"))
saveRDS(edges, file.path(out_int_dir, "Step11_spatial_knn6_edge_list.rds"))

graph_audit <- edges %>%
  dplyr::group_by(SpatialSampleID) %>%
  dplyr::summarise(
    n_edges = dplyr::n(),
    median_distance = median(distance, na.rm = TRUE),
    mean_distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  )
safe_write_csv(graph_audit, file.path(out_table_dir, "Step11_spatial_knn6_graph_audit.csv"))

############################################################
## 7. Moran's I for raw and QC-residualized state scores
############################################################

moran_I <- function(x, edges) {
  ok <- is.finite(x)
  valid_edge <- ok[edges$from] & ok[edges$to]
  ee <- edges[valid_edge, , drop = FALSE]
  xx <- x
  xm <- mean(xx[ok], na.rm = TRUE)
  z <- xx - xm
  denom <- sum(z[ok]^2, na.rm = TRUE)
  if (!is.finite(denom) || denom == 0 || nrow(ee) == 0) return(NA_real_)
  n <- sum(ok)
  S0 <- nrow(ee)
  (n / S0) * sum(z[ee$from] * z[ee$to], na.rm = TRUE) / denom
}

moran_perm <- function(x, edges, n_perm = 999) {
  obs <- moran_I(x, edges)
  if (!is.finite(obs)) return(data.frame(Moran_I = NA_real_, p_empirical = NA_real_, n_perm = n_perm))
  ok <- is.finite(x)
  perms <- numeric(n_perm)
  x_ok <- x[ok]
  idx_ok <- which(ok)
  for (b in seq_len(n_perm)) {
    xp <- x
    xp[idx_ok] <- sample(x_ok, length(x_ok), replace = FALSE)
    perms[b] <- moran_I(xp, edges)
  }
  p <- (sum(abs(perms) >= abs(obs), na.rm = TRUE) + 1) / (sum(is.finite(perms)) + 1)
  data.frame(Moran_I = obs, p_empirical = p, n_perm = n_perm)
}

moran_results <- dplyr::bind_rows(lapply(state_cols, function(st) {
  raw_res <- moran_perm(analysis_resid[[st]], edges, n_perm = MORAN_N_PERM)
  res_col <- paste0(st, "_QCresid")
  qc_res <- moran_perm(analysis_resid[[res_col]], edges, n_perm = MORAN_N_PERM)
  dplyr::bind_rows(
    data.frame(State = st, StateLabel = state_labels[[st]], ScoreType = "Raw", raw_res, stringsAsFactors = FALSE),
    data.frame(State = st, StateLabel = state_labels[[st]], ScoreType = "QC-residualized", qc_res, stringsAsFactors = FALSE)
  )
})) %>%
  dplyr::mutate(FDR = p.adjust(p_empirical, method = "BH"))
safe_write_csv(moran_results, file.path(out_table_dir, "Step11_spatial_Morans_I_raw_and_QCresidual.csv"))

p_moran <- ggplot(moran_results, aes(x = StateLabel, y = Moran_I, fill = ScoreType)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("Raw" = "#E64B35", "QC-residualized" = "#4DBBD5")) +
  labs(title = "A  Moran's I for raw and QC-residualized state scores", x = NULL, y = "Moran's I", fill = "Score type") +
  theme_icb(base_size = 11.5) +
  scale_x_discrete(labels = plot_state_label) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25), legend.position = "top")

############################################################
## 8. High-state calls, neighbor enrichment and dual-high OR
############################################################

call_high <- function(df, state_col, top_fraction = 0.25, mode = c("pooled", "per_sample")) {
  mode <- match.arg(mode)
  cutoff_prob <- 1 - top_fraction
  out <- rep(FALSE, nrow(df))
  cutoff_vec <- rep(NA_real_, nrow(df))
  if (mode == "pooled") {
    cutoff <- as.numeric(quantile(df[[state_col]], probs = cutoff_prob, na.rm = TRUE, names = FALSE))
    out <- df[[state_col]] >= cutoff
    cutoff_vec[] <- cutoff
  } else {
    for (sid in unique(df$SpatialSampleID)) {
      ii <- which(df$SpatialSampleID == sid)
      cutoff <- as.numeric(quantile(df[[state_col]][ii], probs = cutoff_prob, na.rm = TRUE, names = FALSE))
      out[ii] <- df[[state_col]][ii] >= cutoff
      cutoff_vec[ii] <- cutoff
    }
  }
  list(high = out, cutoff = cutoff_vec)
}

compute_high_tables <- function(df, top_fraction = 0.25, mode = "pooled") {
  res <- df
  cutoff_rows <- list()
  for (st in state_cols) {
    h <- call_high(res, st, top_fraction = top_fraction, mode = mode)
    high_col <- paste0(st, "_high")
    cutoff_col <- paste0(st, "_cutoff")
    res[[high_col]] <- h$high
    res[[cutoff_col]] <- h$cutoff
    if (mode == "pooled") {
      cutoff_rows[[st]] <- data.frame(State = st, SpatialSampleID = "pooled", TopFraction = top_fraction, Mode = mode, Cutoff = unique(h$cutoff)[1])
    } else {
      cutoff_rows[[st]] <- res %>%
        dplyr::group_by(SpatialSampleID) %>%
        dplyr::summarise(Cutoff = unique(.data[[cutoff_col]])[1], .groups = "drop") %>%
        dplyr::mutate(State = st, TopFraction = top_fraction, Mode = mode) %>%
        dplyr::select(State, SpatialSampleID, TopFraction, Mode, Cutoff)
    }
  }
  list(df = res, cutoffs = dplyr::bind_rows(cutoff_rows))
}

main_high <- compute_high_tables(analysis_resid, top_fraction = MAIN_HIGH_TOP_FRACTION, mode = "pooled")
high_df <- main_high$df
cutoff_table <- main_high$cutoffs
safe_write_csv(cutoff_table, file.path(out_table_dir, "Step11_spatial_main_top25_pooled_state_score_cutoffs.csv"))

## Neighbor enrichment per state: source high vs neighbor high edge-level OR.
edge_or <- function(high_vec, edges) {
  a <- high_vec[edges$from]
  b <- high_vec[edges$to]
  tab <- table(factor(a, levels = c(TRUE, FALSE)), factor(b, levels = c(TRUE, FALSE)))
  ## Rows source high TRUE/FALSE, cols neighbor high TRUE/FALSE.
  m <- matrix(as.numeric(tab), nrow = 2)
  rownames(m) <- c("source_high", "source_not_high")
  colnames(m) <- c("neighbor_high", "neighbor_not_high")
  ft <- suppressWarnings(fisher.test(m))
  data.frame(
    n_source_high_neighbor_high = m[1, 1],
    n_source_high_neighbor_not_high = m[1, 2],
    n_source_not_high_neighbor_high = m[2, 1],
    n_source_not_high_neighbor_not_high = m[2, 2],
    OR = unname(ft$estimate),
    CI_low = ft$conf.int[1],
    CI_high = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}

neighbor_enrichment <- dplyr::bind_rows(lapply(state_cols, function(st) {
  high_col <- paste0(st, "_high")
  data.frame(State = st, StateLabel = state_labels[[st]], edge_or(high_df[[high_col]], edges), stringsAsFactors = FALSE)
})) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH"))
safe_write_csv(neighbor_enrichment, file.path(out_table_dir, "Step11_spatial_high_state_neighbor_enrichment_top25_pooled.csv"))

p_neighbor <- ggplot(neighbor_enrichment, aes(x = StateLabel, y = OR, fill = State)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.18) +
  scale_fill_manual(values = state_colors, guide = "none") +
  labs(title = "B  High-state neighbor enrichment", x = NULL, y = "Edge-level odds ratio") +
  theme_icb(base_size = 11.5) +
  scale_x_discrete(labels = plot_state_label) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25))

p_moran_neighbor <- p_moran / p_neighbor +
  plot_annotation(title = "Supplementary Figure. Spatial autocorrelation and high-state neighborhood enrichment") &
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
safe_ggsave(file.path(diagnostic_fig_dir, "Supplementary_Figure_S12_spatial_Morans_I_and_neighbor_enrichment.png"), p_moran_neighbor, 12.5, 9.5)
safe_ggsave(file.path(diagnostic_fig_dir, "Supplementary_Figure_S12_spatial_Morans_I_and_neighbor_enrichment.pdf"), p_moran_neighbor, 12.5, 9.5)

## Focused Myeloid-Treg and Tumor/Stromal dual-high classification.
myeloid_col <- "Myeloid_Treg_Immunosuppressive_high"
tumor_col <- "Tumor_dedifferentiation_Stromal_remodeling_high"
high_df$DualHighGroup <- dplyr::case_when(
  high_df[[myeloid_col]] & high_df[[tumor_col]] ~ "Both high",
  high_df[[myeloid_col]] & !high_df[[tumor_col]] ~ "Myeloid-Treg high only",
  !high_df[[myeloid_col]] & high_df[[tumor_col]] ~ "Dediff/Stromal high only",
  TRUE ~ "Neither high"
)
high_df$DualHighGroup <- factor(high_df$DualHighGroup, levels = c("Neither high", "Myeloid-Treg high only", "Dediff/Stromal high only", "Both high"))

dual_counts <- high_df %>%
  dplyr::count(DualHighGroup, name = "n_spots") %>%
  dplyr::mutate(fraction = n_spots / sum(n_spots), unit_of_analysis = "spot", high_definition = "top25 pooled")
safe_write_csv(dual_counts, file.path(out_table_dir, "Step11_spatial_focused_dual_high_group_counts_top25_pooled.csv"))

## Spot-level co-occurrence OR.
spot_tab <- table(
  MyeloidTreg_high = factor(high_df[[myeloid_col]], levels = c(TRUE, FALSE)),
  DediffStromal_high = factor(high_df[[tumor_col]], levels = c(TRUE, FALSE))
)
spot_mat <- matrix(as.numeric(spot_tab), nrow = 2)
rownames(spot_mat) <- c("MyeloidTreg_high", "MyeloidTreg_not_high")
colnames(spot_mat) <- c("DediffStromal_high", "DediffStromal_not_high")
spot_ft <- suppressWarnings(fisher.test(spot_mat))
dual_or <- data.frame(
  Analysis = "MyeloidTreg_high_vs_DediffStromal_high_spot_level",
  UnitOfAnalysis = "spot",
  HighDefinition = "top25 pooled",
  BothHigh = spot_mat[1, 1],
  MyeloidOnly = spot_mat[1, 2],
  DediffStromalOnly = spot_mat[2, 1],
  NeitherHigh = spot_mat[2, 2],
  OR = unname(spot_ft$estimate),
  CI_low = spot_ft$conf.int[1],
  CI_high = spot_ft$conf.int[2],
  p_value = spot_ft$p.value,
  stringsAsFactors = FALSE
)
safe_write_csv(dual_or, file.path(out_table_dir, "Step11_spatial_dual_high_spot_level_OR_top25_pooled.csv"))

## Spatial map.
## Display-only change in v6: reverse the y-axis to match the 10x Visium
## tissue-image orientation used in Step 11A preview. This does not alter
## any downstream statistics, high-state labels, OR, Moran's I or tables.
p_dual_map <- ggplot(high_df, aes(x = spatial_x, y = spatial_y)) +
  geom_point(aes(color = DualHighGroup), size = 0.65, alpha = 0.88) +
  scale_color_manual(
    values = focused_colors,
    drop = FALSE,
    name = "Dual-high group",
    labels = c(
      "Neither high" = "neither high",
      "Myeloid-Treg high only" = "myeloid–Treg high only",
      "Dediff/Stromal high only" = "tumor-dedifferentiation/stromal-remodeling high only",
      "Both high" = "both high"
    )
  ) +
  scale_y_reverse() +
  coord_equal() +
  labs(
    title = "Spot-level dual-high co-occurrence map",
    subtitle = "Visium spots; high-state cutoff: top 25% pooled; y-axis inverted for display",
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_icb(base_size = 12.5) +
  theme(
    legend.position = "right",
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11, margin = ggplot2::margin(b = 6)),
    plot.margin = ggplot2::margin(t = 10, r = 14, b = 8, l = 14)
  )

safe_write_csv(
  data.frame(
    Figure = "Figure8A_dual_high_spatial_colocalization",
    DisplayOnlyChange = TRUE,
    YAxisInvertedForDisplay = TRUE,
    CoordinateSource = coord_audit$coordinate_source[1],
    SubtitleUsed = "Visium spots; high-state cutoff: top 25% pooled; y-axis inverted for display",
    stringsAsFactors = FALSE
  ),
  file.path(out_table_dir, "Step11_Figure8A_display_note.csv")
)

safe_ggsave(file.path(out_fig_dir, "Figure8A_dual_high_spatial_colocalization.png"), p_dual_map, 8.6, 6.8)
safe_ggsave(file.path(out_fig_dir, "Figure8A_dual_high_spatial_colocalization.pdf"), p_dual_map, 8.6, 6.8)

############################################################
## 9. Threshold sensitivity for dual-high co-occurrence
############################################################

sensitivity_rows <- list()
for (mode in c("pooled", "per_sample")) {
  for (tf in SENSITIVITY_TOP_FRACTIONS) {
    tmp <- compute_high_tables(analysis_resid, top_fraction = tf, mode = mode)$df
    mh <- tmp[["Myeloid_Treg_Immunosuppressive_high"]]
    th <- tmp[["Tumor_dedifferentiation_Stromal_remodeling_high"]]
    tab <- table(
      MyeloidTreg_high = factor(mh, levels = c(TRUE, FALSE)),
      DediffStromal_high = factor(th, levels = c(TRUE, FALSE))
    )
    mm <- matrix(as.numeric(tab), nrow = 2)
    rownames(mm) <- c("MyeloidTreg_high", "MyeloidTreg_not_high")
    colnames(mm) <- c("DediffStromal_high", "DediffStromal_not_high")
    ft <- suppressWarnings(fisher.test(mm))
    sensitivity_rows[[paste(mode, tf, sep = "_")]] <- data.frame(
      Mode = mode,
      TopFraction = tf,
      UnitOfAnalysis = "spot",
      BothHigh = mm[1, 1],
      MyeloidOnly = mm[1, 2],
      DediffStromalOnly = mm[2, 1],
      NeitherHigh = mm[2, 2],
      OR = unname(ft$estimate),
      CI_low = ft$conf.int[1],
      CI_high = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
}
threshold_sensitivity <- dplyr::bind_rows(sensitivity_rows) %>%
  dplyr::mutate(FDR = p.adjust(p_value, method = "BH"))
safe_write_csv(threshold_sensitivity, file.path(out_table_dir, "Step11_spatial_dual_high_threshold_sensitivity_OR.csv"))

p_sens <- ggplot(threshold_sensitivity, aes(x = factor(TopFraction), y = OR, color = Mode, group = Mode)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_point(size = 2.6) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.06) +
  scale_color_manual(values = c("pooled" = "#E64B35", "per_sample" = "#4DBBD5")) +
  labs(
    title = "Dual-high co-occurrence threshold sensitivity",
    subtitle = "myeloid–Treg high and tumor-dedifferentiation/stromal-remodeling high;\nunit of analysis: spot",
    x = "Top fraction used to define high-state",
    y = "Spot-level odds ratio",
    color = "Cutoff mode"
  ) +
  theme_icb(base_size = 12.0) +
  theme(plot.subtitle = element_text(size = 10.5), axis.text.x = element_text(angle = 35, hjust = 1, size = 10.75)) +
  theme(legend.position = "top")
safe_ggsave(file.path(diagnostic_fig_dir, "Supplementary_Figure_S13_spatial_dual_high_threshold_sensitivity.png"), p_sens, 7.5, 5.8)
safe_ggsave(file.path(diagnostic_fig_dir, "Supplementary_Figure_S13_spatial_dual_high_threshold_sensitivity.pdf"), p_sens, 7.5, 5.8)

############################################################
## 10. Output inventory and session info
############################################################

output_inventory <- data.frame(
  File = c(
    list.files(out_table_dir, pattern = "Step11_.*\\.csv$|Supplementary_.*\\.csv$|Figure.*\\.csv$", full.names = TRUE),
    list.files(out_fig_dir, pattern = "(Step11|Supplementary|Figure).*\\.(png|pdf)$", full.names = TRUE),
    list.files(out_int_dir, pattern = "Step11_.*\\.rds$", full.names = TRUE)
  ),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    File = normalizePath(File, winslash = "/", mustWork = FALSE),
    FileSizeKB = round(file.info(File)$size / 1024, 1),
    Exists = file.exists(File)
  )
safe_write_csv(output_inventory, file.path(out_table_dir, "Step11_spatial_QC_Moran_colocalization_output_inventory.csv"))

sink(file.path(out_log_dir, "sessionInfo_11_spatial_QC_Moran_colocalization.txt"))
print(sessionInfo())
sink()

message("Step 11 completed successfully.")
message("Selected input object: ", selected)
message("Coordinate source: ", coord_audit$coordinate_source[1])
message("Main dual-high OR saved to: ", file.path(out_table_dir, "Step11_spatial_dual_high_spot_level_OR_top25_pooled.csv"))
message("Main dual-high map saved to: ", file.path(out_fig_dir, "Figure8A_dual_high_spatial_colocalization.png"))

## Public sequential end gate
.step02_out <- file.path(out_table_dir, "Step11_spatial_dual_high_spot_level_OR_top25_pooled.csv")
if (!file.exists(.step02_out)) stop("STEP 02 output missing: ", .step02_out, call. = FALSE)
cat("\nSTEP 02 PASS — spatial QC/Moran/dual-high analysis completed\n")
