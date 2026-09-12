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
## 18F_humoral_TLS_like_spatial_interpretation_revision.R
##
## Purpose:
##   Revision analysis for Reviewer 3 humoral / plasma-cell signal concern.
##
## Reviewer concern:
##   Humoral or plasma-cell-related signals can be associated with either
##   response or resistance depending on tertiary lymphoid structure (TLS)
##   architecture. Visium data alone cannot prove mature TLS.
##
## Main questions:
##   1) Are humoral/plasma/TLS-like transcriptional modules enriched
##      in raw dual-high spots?
##   2) Are these signals spatially organized or mainly diffuse?
##   3) Do TLS-like chemokine/B-cell/T-cell modules co-localize?
##   4) How should IGKC/IGHG1 / plasma-cell signals be interpreted?
##
## Interpretation:
##   This analysis does NOT perform mature TLS calling.
##   It provides TLS-like / humoral spatial context only.
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

############################################################
## 0. Project directory
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

message("Project directory: ", project_dir)

data_processed_dir <- file.path(project_dir, "data_processed")
table_dir          <- file.path(project_dir, "results/tables")
figure_dir         <- file.path(project_dir, "results/figures")
log_dir            <- file.path(project_dir, "logs")

out_table_dir <- file.path(table_dir, "revision_humoral_TLS_like")
out_fig_dir <- file.path(project_dir, "results", "diagnostics", "12_Primary_Visium", "06_humoral_TLS_context")

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
  "Matrix",
  "FNN"
)

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE, type = "binary")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(purrr)
  library(Matrix)
  library(FNN)
})

theme_revision <- function(base_size = 12.5, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2.5),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      axis.title.y = element_text(face = "bold", size = base_size - 1),
      axis.title.x = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "black", size = base_size - 1.5),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1.5),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.60)
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

find_files_recursive <- function(root, pattern) {
  if (!dir.exists(root)) return(character(0))
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
}

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\|.*$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

normalize_colnames <- function(x) {
  x2 <- x
  x2 <- gsub("\\s+", "_", x2)
  x2 <- gsub("-", "_", x2)
  x2 <- gsub("/", "_", x2)
  x2 <- gsub("\\.+", "_", x2)
  x2 <- gsub("__+", "_", x2)
  x2
}

spot_key <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("-1$", "", x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

zvec <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
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

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  
  row_mean <- rowMeans(mat, na.rm = TRUE)
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  
  z <- sweep(mat, 1, row_mean, "-")
  z <- sweep(z, 1, row_sd, "/")
  z[!is.finite(z)] <- NA_real_
  
  z
}

############################################################
## 3. Load 18E Visium state/QC/composition table
############################################################

############################################################
## 3. Load 18E Visium state/QC/composition table
##    Prefer the final 18E table containing Raw_DualHighCategory.
############################################################

visium_comp_candidates <- c(
  file.path(
    table_dir,
    "revision_visium_composition/18E_Visium_spot_state_composition_residualized_scores.csv"
  ),
  file.path(
    table_dir,
    "revision_visium_composition/18E_Visium_spot_state_QC_composition_table.csv"
  ),
  file.path(
    data_processed_dir,
    "18E_Visium_spot_state_QC_composition_table.rds"
  )
)

visium_comp_file <- visium_comp_candidates[file.exists(visium_comp_candidates)][1]

if (is.na(visium_comp_file) || !file.exists(visium_comp_file)) {
  stop(
    "18E Visium state/QC/composition table not found. ",
    "Please run 18E first."
  )
}

message("Reading 18E Visium table: ", visium_comp_file)

if (grepl("\\.rds$", visium_comp_file, ignore.case = TRUE)) {
  visium_df <- readRDS(visium_comp_file)
} else {
  visium_df <- data.table::fread(
    visium_comp_file,
    data.table = FALSE,
    check.names = FALSE
  )
}

colnames(visium_df) <- normalize_colnames(colnames(visium_df))

## If the table does not contain Raw_DualHighCategory, reconstruct it
## from the raw Myeloid–Treg and Dediff/Stromal state scores.
if (!"Raw_DualHighCategory" %in% colnames(visium_df)) {
  message("Raw_DualHighCategory not found. Reconstructing from raw state-score top-quartile cutoffs.")
  
  required_state_for_category <- c(
    "Myeloid_Treg_Immunosuppressive",
    "Tumor_dedifferentiation_Stromal_remodeling"
  )
  
  missing_for_category <- setdiff(required_state_for_category, colnames(visium_df))
  
  if (length(missing_for_category) > 0) {
    stop(
      "Cannot reconstruct Raw_DualHighCategory. Missing columns: ",
      paste(missing_for_category, collapse = ", ")
    )
  }
  
  visium_df$Myeloid_Treg_Immunosuppressive <- suppressWarnings(
    as.numeric(visium_df$Myeloid_Treg_Immunosuppressive)
  )
  
  visium_df$Tumor_dedifferentiation_Stromal_remodeling <- suppressWarnings(
    as.numeric(visium_df$Tumor_dedifferentiation_Stromal_remodeling)
  )
  
  myeloid_cut <- quantile(
    visium_df$Myeloid_Treg_Immunosuppressive,
    probs = 0.75,
    na.rm = TRUE,
    names = FALSE
  )
  
  dediff_cut <- quantile(
    visium_df$Tumor_dedifferentiation_Stromal_remodeling,
    probs = 0.75,
    na.rm = TRUE,
    names = FALSE
  )
  
  myeloid_high <- visium_df$Myeloid_Treg_Immunosuppressive >= myeloid_cut
  dediff_high  <- visium_df$Tumor_dedifferentiation_Stromal_remodeling >= dediff_cut
  
  visium_df$Raw_DualHighCategory <- ifelse(
    myeloid_high & dediff_high,
    "Both-High",
    ifelse(
      myeloid_high & !dediff_high,
      "Myeloid-Treg-High Only",
      ifelse(
        !myeloid_high & dediff_high,
        "Dediff/Stromal-High Only",
        "Neither-High"
      )
    )
  )
  
  message(
    "Reconstructed Raw_DualHighCategory using top-quartile cutoffs: ",
    "Myeloid-Treg cutoff = ", signif(myeloid_cut, 4),
    "; Dediff/Stromal cutoff = ", signif(dediff_cut, 4)
  )
}

required_cols <- c(
  "Spot",
  "spatial_x",
  "spatial_y",
  "Raw_DualHighCategory",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling"
)

missing_cols <- setdiff(required_cols, colnames(visium_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns from 18E table: ", paste(missing_cols, collapse = ", "))
}

visium_df$Raw_DualHighCategory <- factor(
  visium_df$Raw_DualHighCategory,
  levels = c(
    "Neither-High",
    "Myeloid-Treg-High Only",
    "Dediff/Stromal-High Only",
    "Both-High"
  )
)

for (cc in c("spatial_x", "spatial_y",
             "Myeloid_Treg_Immunosuppressive",
             "Tumor_dedifferentiation_Stromal_remodeling")) {
  visium_df[[cc]] <- suppressWarnings(as.numeric(visium_df[[cc]]))
}

message("Visium table retained: ", nrow(visium_df), " spots")

safe_write_csv(
  data.frame(
    InputVisiumTable = visium_comp_file,
    n_spots = nrow(visium_df),
    n_cols = ncol(visium_df),
    stringsAsFactors = FALSE
  ),
  file.path(out_table_dir, "18F_visium_input_summary.csv")
)

############################################################
## 4. Find expression source from 18E summary or manual env
############################################################

extract_seurat_expression_matrix <- function(obj) {
  if (!("Seurat" %in% class(obj))) {
    stop("Object is not a Seurat object.")
  }
  
  if (!requireNamespace("SeuratObject", quietly = TRUE) &&
      !requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat object found, but neither SeuratObject nor Seurat package is available.")
  }
  
  assay_names <- tryCatch(names(obj@assays), error = function(e) character(0))
  
  if (length(assay_names) == 0) {
    stop("No assay found in Seurat object.")
  }
  
  preferred_assays <- intersect(c("Spatial", "SCT", "RNA"), assay_names)
  if (length(preferred_assays) == 0) preferred_assays <- assay_names
  
  assay_use <- preferred_assays[1]
  message("Using Seurat assay: ", assay_use)
  
  mat <- NULL
  
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    mat <- tryCatch(
      {
        assay_obj <- obj[[assay_use]]
        available_layers <- tryCatch(
          SeuratObject::Layers(assay_obj),
          error = function(e) character(0)
        )
        
        message(
          "Available layers in assay ",
          assay_use,
          ": ",
          paste(available_layers, collapse = ", ")
        )
        
        preferred_layers <- c(
          "data",
          "counts",
          "scale.data",
          grep("^data", available_layers, value = TRUE),
          grep("^counts", available_layers, value = TRUE)
        )
        
        preferred_layers <- unique(preferred_layers)
        layer_use <- intersect(preferred_layers, available_layers)
        
        if (length(layer_use) == 0 && length(available_layers) > 0) {
          layer_use <- available_layers[1]
        }
        
        if (length(layer_use) == 0) {
          stop("No layer found in Seurat assay: ", assay_use)
        }
        
        layer_use <- layer_use[1]
        message("Using Seurat layer: ", layer_use)
        
        m1 <- tryCatch(
          SeuratObject::LayerData(object = assay_obj, layer = layer_use),
          error = function(e1) {
            message("LayerData(assay_obj) failed: ", paste(as.character(e1$message), collapse = " | "))
            
            tryCatch(
              SeuratObject::LayerData(object = obj, assay = assay_use, layer = layer_use),
              error = function(e2) {
                message("LayerData(Seurat object) failed: ", paste(as.character(e2$message), collapse = " | "))
                NULL
              }
            )
          }
        )
        
        if (is.null(m1)) stop("LayerData extraction returned NULL.")
        m1
      },
      error = function(e_layer) {
        message("LayerData route failed: ", paste(as.character(e_layer$message), collapse = " | "))
        
        tryCatch(
          SeuratObject::GetAssayData(
            object = obj,
            assay = assay_use,
            layer = "data"
          ),
          error = function(e_data) {
            message("GetAssayData(layer='data') failed: ", paste(as.character(e_data$message), collapse = " | "))
            
            SeuratObject::GetAssayData(
              object = obj,
              assay = assay_use,
              layer = "counts"
            )
          }
        )
      }
    )
  }
  
  if (is.null(mat) && requireNamespace("Seurat", quietly = TRUE)) {
    mat <- tryCatch(
      Seurat::GetAssayData(
        object = obj,
        assay = assay_use,
        layer = "data"
      ),
      error = function(e_data) {
        message("Seurat::GetAssayData(layer='data') failed: ", paste(as.character(e_data$message), collapse = " | "))
        
        Seurat::GetAssayData(
          object = obj,
          assay = assay_use,
          layer = "counts"
        )
      }
    )
  }
  
  if (is.null(mat)) stop("Failed to extract expression matrix from Seurat object.")
  
  mat <- as.matrix(mat)
  collapse_duplicate_genes(mat)
}

read_expression_table <- function(file) {
  message("Trying expression source: ", file)
  
  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    obj <- readRDS(file)
    
    if ("Seurat" %in% class(obj)) {
      message("Detected Seurat object.")
      return(extract_seurat_expression_matrix(obj))
    }
    
    if (is.matrix(obj) || is.data.frame(obj)) {
      mat <- as.matrix(obj)
      if (is.null(rownames(mat))) stop("RDS matrix/data.frame has no rownames.")
      storage.mode(mat) <- "numeric"
      return(collapse_duplicate_genes(mat))
    }
    
    if (is.list(obj)) {
      candidates <- list()
      
      for (nm in names(obj)) {
        item <- obj[[nm]]
        if (is.matrix(item) || is.data.frame(item)) {
          mat <- as.matrix(item)
          if (!is.null(rownames(mat)) && !is.null(colnames(mat))) {
            candidates[[nm]] <- mat
          }
        }
      }
      
      if (length(candidates) == 0) {
        stop("RDS list contains no matrix/data.frame with rownames and colnames.")
      }
      
      spot_keys <- spot_key(visium_df$Spot)
      
      score_one <- function(mat) {
        ckeys <- spot_key(colnames(mat))
        rkeys <- spot_key(rownames(mat))
        overlap_cols <- sum(spot_keys %in% ckeys)
        overlap_rows <- sum(spot_keys %in% rkeys)
        
        if (overlap_cols >= overlap_rows) {
          return(c(overlap = overlap_cols, orientation = "genes_rows"))
        } else {
          return(c(overlap = overlap_rows, orientation = "genes_cols"))
        }
      }
      
      sc <- lapply(candidates, score_one)
      overlap <- sapply(sc, function(x) as.numeric(x[["overlap"]]))
      best_nm <- names(which.max(overlap))
      best_mat <- candidates[[best_nm]]
      
      if (as.numeric(sc[[best_nm]][["overlap"]]) < 10) {
        stop("No list matrix has meaningful spot overlap.")
      }
      
      if (sc[[best_nm]][["orientation"]] == "genes_cols") {
        best_mat <- t(best_mat)
      }
      
      storage.mode(best_mat) <- "numeric"
      message("Using matrix from list element: ", best_nm)
      return(collapse_duplicate_genes(best_mat))
    }
    
    stop("Unsupported RDS object class: ", paste(class(obj), collapse = ", "))
  }
  
  if (grepl("\\.csv$|\\.tsv$|\\.txt$", file, ignore.case = TRUE)) {
    sep <- ifelse(grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE), "\t", ",")
    
    dt <- data.table::fread(
      file,
      sep = sep,
      data.table = FALSE,
      check.names = FALSE
    )
    
    cn <- colnames(dt)
    cn_norm <- tolower(gsub("[^a-z0-9]", "", cn))
    
    gene_col_idx <- which(
      cn_norm %in% c("genesymbol", "gene", "genes", "symbol", "hgncsymbol", "features")
    )
    
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
    
    expr_mat <- as.matrix(
      sapply(expr_df, function(x) suppressWarnings(as.numeric(as.character(x))))
    )
    
    colnames(expr_mat) <- colnames(expr_df)
    rownames(expr_mat) <- genes
    
    return(collapse_duplicate_genes(expr_mat))
  }
  
  if (grepl("\\.h5$", file, ignore.case = TRUE)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("H5 file found, but Seurat package is not installed. Cannot use Read10X_h5.")
    }
    
    mat <- Seurat::Read10X_h5(file)
    
    if (is.list(mat)) {
      if ("Gene Expression" %in% names(mat)) {
        mat <- mat[["Gene Expression"]]
      } else {
        mat <- mat[[1]]
      }
    }
    
    mat <- as.matrix(mat)
    return(collapse_duplicate_genes(mat))
  }
  
  stop("Unsupported expression file format: ", file)
}

read_10x_mtx_folder <- function(folder) {
  message("Trying 10x folder: ", folder)
  
  matrix_file <- find_first_existing(c(
    file.path(folder, "matrix.mtx"),
    file.path(folder, "matrix.mtx.gz")
  ))
  
  features_file <- find_first_existing(c(
    file.path(folder, "features.tsv"),
    file.path(folder, "features.tsv.gz"),
    file.path(folder, "genes.tsv"),
    file.path(folder, "genes.tsv.gz")
  ))
  
  barcodes_file <- find_first_existing(c(
    file.path(folder, "barcodes.tsv"),
    file.path(folder, "barcodes.tsv.gz")
  ))
  
  if (!file.exists(matrix_file) || !file.exists(features_file) || !file.exists(barcodes_file)) {
    stop("Incomplete 10x matrix folder.")
  }
  
  mat <- Matrix::readMM(matrix_file)
  features <- data.table::fread(features_file, header = FALSE, data.table = FALSE)
  barcodes <- data.table::fread(barcodes_file, header = FALSE, data.table = FALSE)
  
  genes <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  rownames(mat) <- clean_gene_symbol(genes)
  colnames(mat) <- as.character(barcodes[[1]])
  
  mat <- as.matrix(mat)
  collapse_duplicate_genes(mat)
}

get_expression_source_from_18E <- function() {
  manual <- Sys.getenv("VISIUM_EXPR_FILE")
  if (manual != "" && file.exists(manual)) return(manual)
  
  summary_file <- file.path(
    table_dir,
    "revision_visium_composition/18E_selected_visium_expression_source_summary.csv"
  )
  
  if (file.exists(summary_file)) {
    sm <- data.table::fread(summary_file, data.table = FALSE)
    if ("SelectedExpressionSource" %in% colnames(sm)) {
      src <- sm$SelectedExpressionSource[1]
      if (!is.na(src) && file.exists(src)) return(src)
    }
  }
  
  stop(
    "No expression source found. Please run 18E successfully first, ",
    "or set Sys.setenv(VISIUM_EXPR_FILE='path/to/expression_or_seurat.rds')."
  )
}

expr_source <- get_expression_source_from_18E()
message("Using expression source: ", expr_source)

if (dir.exists(expr_source)) {
  visium_expr <- read_10x_mtx_folder(expr_source)
} else {
  visium_expr <- read_expression_table(expr_source)
}

safe_write_csv(
  data.frame(
    SelectedExpressionSource = expr_source,
    n_genes = nrow(visium_expr),
    n_expression_spots = ncol(visium_expr),
    n_meta_spots = nrow(visium_df),
    exact_spot_overlap = sum(visium_df$Spot %in% colnames(visium_expr)),
    key_spot_overlap = sum(spot_key(visium_df$Spot) %in% spot_key(colnames(visium_expr))),
    stringsAsFactors = FALSE
  ),
  file.path(out_table_dir, "18F_selected_expression_source_summary.csv")
)

############################################################
## 5. Align expression to spots
############################################################

align_expr_to_meta <- function(expr_mat, meta_df) {
  expr_cols <- colnames(expr_mat)
  meta_spots <- as.character(meta_df$Spot)
  
  idx_exact <- match(meta_spots, expr_cols)
  exact_matched <- !is.na(idx_exact)
  
  if (sum(exact_matched) >= 0.8 * nrow(meta_df)) {
    message("Using exact spot matching: ", sum(exact_matched), " spots")
    idx <- idx_exact
    method <- "exact"
  } else {
    expr_key <- spot_key(expr_cols)
    meta_key <- spot_key(meta_spots)
    idx_key <- match(meta_key, expr_key)
    key_matched <- !is.na(idx_key)
    message("Using key-based spot matching: ", sum(key_matched), " spots")
    idx <- idx_key
    method <- "key"
  }
  
  matched <- !is.na(idx)
  
  audit <- data.frame(
    Spot = meta_spots,
    SpotKey = spot_key(meta_spots),
    MatchedExpressionColumn = ifelse(matched, expr_cols[idx], NA_character_),
    Matched = matched,
    MatchingMethod = method,
    stringsAsFactors = FALSE
  )
  
  if (sum(matched) < 1000) {
    stop("Too few matched Visium spots between metadata and expression matrix: ", sum(matched))
  }
  
  expr2 <- expr_mat[, idx[matched], drop = FALSE]
  colnames(expr2) <- meta_spots[matched]
  
  meta2 <- meta_df[matched, , drop = FALSE]
  
  list(expr = expr2, meta = meta2, audit = audit)
}

aligned <- align_expr_to_meta(visium_expr, visium_df)

visium_expr_aligned <- aligned$expr
visium_df_aligned <- aligned$meta

safe_write_csv(
  aligned$audit,
  file.path(out_table_dir, "18F_expression_metadata_spot_matching_audit.csv")
)

message("Aligned expression: ", nrow(visium_expr_aligned), " genes x ", ncol(visium_expr_aligned), " spots")

############################################################
## 6. Define humoral / TLS-like modules
############################################################

humoral_tls_modules <- list(
  TLS_Chemokine_Core = c(
    "CXCL13", "CCL19", "CCL21", "LTA", "LTB", "CXCR5", "CCR7"
  ),
  B_Cell = c(
    "MS4A1", "CD79A", "CD79B", "CD19", "CD22", "BANK1", "PAX5", "CD74"
  ),
  Plasma_Cell = c(
    "MZB1", "JCHAIN", "IGKC", "IGHG1", "IGHG3", "XBP1", "SDC1", "PRDM1"
  ),
  T_Cell_Zone = c(
    "CD3D", "CD3E", "CD3G", "TRAC", "IL7R", "CCR7", "LTB", "CXCR5"
  ),
  Tfh_Like = c(
    "CXCR5", "PDCD1", "ICOS", "BCL6", "IL21", "TOX", "SH2D1A"
  ),
  FDC_HEV_Like = c(
    "CR2", "FDCSP", "CXCL13", "CCL19", "CCL21", "PDPN", "ACKR1", "SELE", "ICAM1"
  ),
  IG_Humoral = c(
    "IGKC", "IGHG1", "IGHG3", "IGHG4", "IGHM", "IGHA1", "JCHAIN", "MZB1"
  )
)

module_cols <- names(humoral_tls_modules)

compute_module_scores <- function(expr_mat, modules) {
  expr_mat <- collapse_duplicate_genes(expr_mat)
  gene_z <- zscore_rows(expr_mat)
  
  score_list <- list()
  gene_audit <- list()
  
  for (module_name in names(modules)) {
    genes <- clean_gene_symbol(modules[[module_name]])
    matched <- intersect(genes, rownames(gene_z))
    
    gene_audit[[module_name]] <- data.frame(
      Module = module_name,
      RequestedGene = genes,
      Present = genes %in% rownames(gene_z),
      stringsAsFactors = FALSE
    )
    
    if (length(matched) == 0) {
      warning("No matched genes for module: ", module_name)
      score <- rep(NA_real_, ncol(gene_z))
    } else {
      score <- colMeans(gene_z[matched, , drop = FALSE], na.rm = TRUE)
    }
    
    score_list[[module_name]] <- score
  }
  
  score_df <- as.data.frame(score_list, check.names = FALSE)
  score_df$Spot <- colnames(expr_mat)
  score_df <- score_df[, c("Spot", names(modules)), drop = FALSE]
  
  gene_audit_df <- bind_rows(gene_audit) %>%
    group_by(Module) %>%
    mutate(
      n_requested = n(),
      n_present = sum(Present),
      n_missing = sum(!Present)
    ) %>%
    ungroup()
  
  list(scores = score_df, gene_audit = gene_audit_df)
}

tls_scores <- compute_module_scores(visium_expr_aligned, humoral_tls_modules)

safe_write_csv(
  tls_scores$scores,
  file.path(out_table_dir, "18F_humoral_TLS_like_module_scores.csv")
)

safe_write_csv(
  tls_scores$gene_audit,
  file.path(out_table_dir, "18F_humoral_TLS_like_marker_gene_presence_audit.csv")
)

tls_df <- visium_df_aligned %>%
  left_join(tls_scores$scores, by = "Spot")

safe_write_csv(
  tls_df,
  file.path(out_table_dir, "18F_Visium_spot_state_composition_TLS_like_table.csv")
)

safe_save_rds(
  tls_df,
  file.path(data_processed_dir, "18F_Visium_spot_state_composition_TLS_like_table.rds")
)

############################################################
## 7. Compare TLS-like modules across raw dual-high categories
############################################################

tls_long <- tls_df %>%
  select(Spot, Raw_DualHighCategory, all_of(module_cols)) %>%
  pivot_longer(
    cols = all_of(module_cols),
    names_to = "TLSModule",
    values_to = "ModuleScore"
  )

tls_group_summary <- tls_long %>%
  group_by(TLSModule, Raw_DualHighCategory) %>%
  summarise(
    n = sum(is.finite(ModuleScore)),
    mean = mean(ModuleScore, na.rm = TRUE),
    median = median(ModuleScore, na.rm = TRUE),
    IQR_low = quantile(ModuleScore, 0.25, na.rm = TRUE),
    IQR_high = quantile(ModuleScore, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

safe_write_csv(
  tls_group_summary,
  file.path(out_table_dir, "18F_TLS_like_module_summary_by_raw_dual_high_category.csv")
)

tls_both_vs_neither <- tls_long %>%
  filter(Raw_DualHighCategory %in% c("Both-High", "Neither-High")) %>%
  group_by(TLSModule) %>%
  summarise(
    n_both = sum(Raw_DualHighCategory == "Both-High" & is.finite(ModuleScore)),
    n_neither = sum(Raw_DualHighCategory == "Neither-High" & is.finite(ModuleScore)),
    median_both = median(ModuleScore[Raw_DualHighCategory == "Both-High"], na.rm = TRUE),
    median_neither = median(ModuleScore[Raw_DualHighCategory == "Neither-High"], na.rm = TRUE),
    median_difference_both_minus_neither = median_both - median_neither,
    wilcox_p = tryCatch(
      wilcox.test(ModuleScore ~ Raw_DualHighCategory)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(wilcox_p_adj_BH = p.adjust(wilcox_p, method = "BH"))

safe_write_csv(
  tls_both_vs_neither,
  file.path(out_table_dir, "18F_TLS_like_BothHigh_vs_NeitherHigh_Wilcoxon.csv")
)

p_tls_box <- tls_long %>%
  mutate(
    Raw_DualHighCategory = factor(
      Raw_DualHighCategory,
      levels = c(
        "Neither-High",
        "Myeloid-Treg-High Only",
        "Dediff/Stromal-High Only",
        "Both-High"
      )
    )
  ) %>%
  ggplot(aes(x = Raw_DualHighCategory, y = ModuleScore)) +
  geom_boxplot(outlier.size = 0.4, width = 0.65) +
  facet_wrap(~ TLSModule, scales = "free_y", ncol = 3) +
  labs(
    title = "Humoral and TLS-like transcriptional modules across Visium high-state categories",
    subtitle = "Used for cautious TLS-like / humoral spatial interpretation, not mature TLS calling",
    x = "Raw high-state category",
    y = "Module score"
  ) +
  theme_revision(base_size = 11.5) +
  theme(plot.subtitle = element_text(size = 10.5), axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25))

safe_ggsave(
  file.path(out_fig_dir, "18F_TLS_like_modules_by_raw_dual_high_category.pdf"),
  p_tls_box,
  width = 11,
  height = 7
)

safe_ggsave(
  file.path(out_fig_dir, "18F_TLS_like_modules_by_raw_dual_high_category.png"),
  p_tls_box,
  width = 11,
  height = 7
)

############################################################
## 8. Correlations with states and composition modules
############################################################

cor_test_pair <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  
  if (n < 4) {
    return(data.frame(n = n, rho = NA_real_, p_value = NA_real_))
  }
  
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  
  data.frame(
    n = n,
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

state_cols <- intersect(
  c(
    "Immune_defective_Cold",
    "Myeloid_Treg_Immunosuppressive",
    "Tumor_dedifferentiation_Stromal_remodeling",
    "Melanocytic_Differentiation"
  ),
  colnames(tls_df)
)

state_display_public <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

composition_reference_cols <- intersect(
  c(
    "CAF_Stromal",
    "Myeloid",
    "T_NK",
    "B_Plasma",
    "Melanoma_Lineage",
    "Endothelial"
  ),
  colnames(tls_df)
)

cor_out <- list()

for (m in module_cols) {
  for (target in c(state_cols, composition_reference_cols)) {
    tmp <- cor_test_pair(tls_df[[m]], tls_df[[target]], method = "spearman")
    tmp$TLSModule <- m
    tmp$Target <- target
    tmp$TargetType <- ifelse(target %in% state_cols, "State", "Composition")
    cor_out[[length(cor_out) + 1]] <- tmp
  }
}

tls_cor_tbl <- bind_rows(cor_out) %>%
  select(TLSModule, TargetType, Target, n, rho, p_value) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

safe_write_csv(
  tls_cor_tbl,
  file.path(out_table_dir, "18F_TLS_like_module_correlations_with_states_and_composition.csv")
)

p_tls_cor <- tls_cor_tbl %>%
  mutate(
    TLSModule = factor(TLSModule, levels = module_cols),
    label = sprintf("%.2f", rho)
  ) %>%
  ggplot(aes(x = Target, y = TLSModule, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 2.8) +
  facet_wrap(~ TargetType, scales = "free_x") +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  labs(
    title = "Humoral/TLS-like module associations with states and composition modules",
    subtitle = "Spearman correlations across Visium spots",
    x = NULL,
    y = "Humoral/TLS-like module",
    fill = "Spearman\nrho"
  ) +
  scale_x_discrete(labels = function(x) { y <- x; hit <- y %in% names(state_display_public); y[hit] <- unname(state_display_public[y[hit]]); y }) +
  theme_revision(base_size = 11.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25))

safe_ggsave(
  file.path(out_fig_dir, "18F_TLS_like_module_correlation_heatmap.pdf"),
  p_tls_cor,
  width = 11,
  height = 6
)

safe_ggsave(
  file.path(out_fig_dir, "18F_TLS_like_module_correlation_heatmap.png"),
  p_tls_cor,
  width = 11,
  height = 6
)

############################################################
## 9. Spatial maps of humoral/TLS-like modules
############################################################

tls_map_long <- tls_df %>%
  select(Spot, spatial_x, spatial_y, all_of(module_cols)) %>%
  pivot_longer(
    cols = all_of(module_cols),
    names_to = "TLSModule",
    values_to = "ModuleScore"
  )

p_tls_map <- ggplot(
  tls_map_long,
  aes(x = spatial_x, y = spatial_y, color = ModuleScore)
) +
  geom_point(size = 0.85, alpha = 0.9) +
  scale_y_reverse() +
  coord_equal() +
  facet_wrap(~ TLSModule, ncol = 3) +
  scale_color_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    na.value = "grey85"
  ) +
  labs(
    title = "Spatial distribution of humoral/TLS-like modules",
    subtitle = "Transcriptional spatial context; not mature TLS identification",
    x = "Spatial x",
    y = "Spatial y",
    color = "Module\nscore"
  ) +
  theme_revision(base_size = 11.5)

safe_ggsave(
  file.path(out_fig_dir, "18F_spatial_maps_of_TLS_like_modules.pdf"),
  p_tls_map,
  width = 10.5,
  height = 7.5
)

safe_ggsave(
  file.path(out_fig_dir, "18F_spatial_maps_of_TLS_like_modules.png"),
  p_tls_map,
  width = 10.5,
  height = 7.5
)

############################################################
## 10. Spatial autocorrelation and co-localization
############################################################

make_knn_index <- function(coords, k = 6) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"
  FNN::get.knn(coords, k = k)$nn.index
}

spatial_lag <- function(v, knn_idx) {
  v <- as.numeric(v)
  apply(knn_idx, 1, function(idx) mean(v[idx], na.rm = TRUE))
}

global_univariate_moran <- function(v, coords, k = 6, n_perm = 999, seed = 20260608) {
  set.seed(seed + k)
  
  ok <- is.finite(v) & is.finite(coords[, 1]) & is.finite(coords[, 2])
  v <- as.numeric(v[ok])
  coords <- as.matrix(coords[ok, , drop = FALSE])
  
  if (length(v) < k + 5) {
    return(data.frame(
      k = k,
      n = length(v),
      Moran_I = NA_real_,
      perm_mean = NA_real_,
      perm_sd = NA_real_,
      z_perm = NA_real_,
      p_greater = NA_real_,
      n_perm = n_perm
    ))
  }
  
  knn_idx <- make_knn_index(coords, k = k)
  zv <- zvec(v)
  lag_v <- spatial_lag(zv, knn_idx)
  
  I_obs <- mean(zv * lag_v, na.rm = TRUE)
  
  perm_vals <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    v_perm <- sample(zv, length(zv), replace = FALSE)
    lag_perm <- spatial_lag(v_perm, knn_idx)
    perm_vals[b] <- mean(zv * lag_perm, na.rm = TRUE)
  }
  
  perm_mean <- mean(perm_vals, na.rm = TRUE)
  perm_sd <- sd(perm_vals, na.rm = TRUE)
  
  data.frame(
    k = k,
    n = length(v),
    Moran_I = I_obs,
    perm_mean = perm_mean,
    perm_sd = perm_sd,
    z_perm = (I_obs - perm_mean) / perm_sd,
    p_greater = (sum(perm_vals >= I_obs, na.rm = TRUE) + 1) / (n_perm + 1),
    n_perm = n_perm
  )
}

global_bivariate_moran <- function(x, y, coords, k = 6, n_perm = 999, seed = 20260608) {
  set.seed(seed + k)
  
  ok <- is.finite(x) & is.finite(y) & is.finite(coords[, 1]) & is.finite(coords[, 2])
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  coords <- as.matrix(coords[ok, , drop = FALSE])
  
  if (length(x) < k + 5) {
    return(data.frame(
      k = k,
      n = length(x),
      I_symmetric = NA_real_,
      perm_mean = NA_real_,
      perm_sd = NA_real_,
      z_perm = NA_real_,
      p_greater = NA_real_,
      n_perm = n_perm
    ))
  }
  
  knn_idx <- make_knn_index(coords, k = k)
  
  zx <- zvec(x)
  zy <- zvec(y)
  
  lag_y <- spatial_lag(zy, knn_idx)
  lag_x <- spatial_lag(zx, knn_idx)
  
  I_xy <- mean(zx * lag_y, na.rm = TRUE)
  I_yx <- mean(zy * lag_x, na.rm = TRUE)
  I_sym <- mean(c(I_xy, I_yx), na.rm = TRUE)
  
  perm_vals <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    zy_perm <- sample(zy, length(zy), replace = FALSE)
    lag_y_perm <- spatial_lag(zy_perm, knn_idx)
    I_xy_perm <- mean(zx * lag_y_perm, na.rm = TRUE)
    
    zx_perm <- sample(zx, length(zx), replace = FALSE)
    lag_x_perm <- spatial_lag(zx_perm, knn_idx)
    I_yx_perm <- mean(zy * lag_x_perm, na.rm = TRUE)
    
    perm_vals[b] <- mean(c(I_xy_perm, I_yx_perm), na.rm = TRUE)
  }
  
  perm_mean <- mean(perm_vals, na.rm = TRUE)
  perm_sd <- sd(perm_vals, na.rm = TRUE)
  
  data.frame(
    k = k,
    n = length(x),
    I_symmetric = I_sym,
    perm_mean = perm_mean,
    perm_sd = perm_sd,
    z_perm = (I_sym - perm_mean) / perm_sd,
    p_greater = (sum(perm_vals >= I_sym, na.rm = TRUE) + 1) / (n_perm + 1),
    n_perm = n_perm
  )
}

coords <- as.matrix(tls_df[, c("spatial_x", "spatial_y")])
k_values <- c(4, 6, 8, 12)

moran_uni <- list()

for (m in module_cols) {
  for (kk in k_values) {
    tmp <- global_univariate_moran(
      v = tls_df[[m]],
      coords = coords,
      k = kk,
      n_perm = 999,
      seed = 20260608
    )
    tmp$TLSModule <- m
    moran_uni[[length(moran_uni) + 1]] <- tmp
  }
}

moran_uni_tbl <- bind_rows(moran_uni) %>%
  select(TLSModule, k, n, Moran_I, perm_mean, perm_sd, z_perm, p_greater, n_perm)

safe_write_csv(
  moran_uni_tbl,
  file.path(out_table_dir, "18F_TLS_like_module_univariate_spatial_Moran.csv")
)

bivar_pairs <- list(
  TLSCore_vs_TCellZone = c("TLS_Chemokine_Core", "T_Cell_Zone"),
  TLSCore_vs_BCell = c("TLS_Chemokine_Core", "B_Cell"),
  TLSCore_vs_Plasma = c("TLS_Chemokine_Core", "Plasma_Cell"),
  BCell_vs_TCellZone = c("B_Cell", "T_Cell_Zone"),
  Plasma_vs_TCellZone = c("Plasma_Cell", "T_Cell_Zone"),
  IGHumoral_vs_TCellZone = c("IG_Humoral", "T_Cell_Zone")
)

moran_bi <- list()

for (pair_name in names(bivar_pairs)) {
  pair <- bivar_pairs[[pair_name]]
  
  if (!all(pair %in% colnames(tls_df))) next
  
  for (kk in k_values) {
    tmp <- global_bivariate_moran(
      x = tls_df[[pair[1]]],
      y = tls_df[[pair[2]]],
      coords = coords,
      k = kk,
      n_perm = 999,
      seed = 20260608
    )
    tmp$PairName <- pair_name
    tmp$Module1 <- pair[1]
    tmp$Module2 <- pair[2]
    moran_bi[[length(moran_bi) + 1]] <- tmp
  }
}

moran_bi_tbl <- bind_rows(moran_bi) %>%
  select(PairName, Module1, Module2, k, n, I_symmetric, perm_mean, perm_sd, z_perm, p_greater, n_perm)

safe_write_csv(
  moran_bi_tbl,
  file.path(out_table_dir, "18F_TLS_like_bivariate_spatial_colocalization_Moran.csv")
)

p_moran_uni <- moran_uni_tbl %>%
  mutate(
    TLSModule = factor(TLSModule, levels = module_cols),
    k = factor(k, levels = k_values)
  ) %>%
  ggplot(aes(x = k, y = Moran_I, group = TLSModule)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line() +
  geom_point(size = 1.8) +
  facet_wrap(~ TLSModule, ncol = 3) +
  labs(
    title = "Spatial autocorrelation of humoral/TLS-like modules",
    subtitle = "Positive values indicate spatial aggregation within the Visium section",
    x = "k nearest neighbors",
    y = "Moran-type statistic"
  ) +
  theme_revision(base_size = 11.5)

safe_ggsave(
  file.path(out_fig_dir, "18F_TLS_like_univariate_spatial_Moran.pdf"),
  p_moran_uni,
  width = 10.5,
  height = 7
)

safe_ggsave(
  file.path(out_fig_dir, "18F_TLS_like_univariate_spatial_Moran.png"),
  p_moran_uni,
  width = 10.5,
  height = 7
)

############################################################
## 11. Operational TLS-like co-high spatial categories
############################################################

assign_high <- function(v, q = 0.75) {
  cut <- quantile(v, q, na.rm = TRUE, names = FALSE)
  v >= cut
}

tls_df <- tls_df %>%
  mutate(
    TLSCore_high = assign_high(TLS_Chemokine_Core, 0.75),
    BCell_high = assign_high(B_Cell, 0.75),
    Plasma_high = assign_high(Plasma_Cell, 0.75),
    TCellZone_high = assign_high(T_Cell_Zone, 0.75),
    IGHumoral_high = assign_high(IG_Humoral, 0.75),
    Operational_TLS_like_cohigh = TLSCore_high & BCell_high & TCellZone_high,
    Humoral_Tcell_cohigh = IGHumoral_high & TCellZone_high,
    Plasma_Tcell_cohigh = Plasma_high & TCellZone_high
  )

safe_write_csv(
  tls_df,
  file.path(out_table_dir, "18F_Visium_TLS_like_operational_cohigh_table.csv")
)

category_tls_counts <- tls_df %>%
  group_by(Raw_DualHighCategory) %>%
  summarise(
    n_spots = n(),
    n_TLSCore_high = sum(TLSCore_high, na.rm = TRUE),
    frac_TLSCore_high = n_TLSCore_high / n_spots,
    n_BCell_high = sum(BCell_high, na.rm = TRUE),
    frac_BCell_high = n_BCell_high / n_spots,
    n_TCellZone_high = sum(TCellZone_high, na.rm = TRUE),
    frac_TCellZone_high = n_TCellZone_high / n_spots,
    n_Operational_TLS_like_cohigh = sum(Operational_TLS_like_cohigh, na.rm = TRUE),
    frac_Operational_TLS_like_cohigh = n_Operational_TLS_like_cohigh / n_spots,
    n_Humoral_Tcell_cohigh = sum(Humoral_Tcell_cohigh, na.rm = TRUE),
    frac_Humoral_Tcell_cohigh = n_Humoral_Tcell_cohigh / n_spots,
    n_Plasma_Tcell_cohigh = sum(Plasma_Tcell_cohigh, na.rm = TRUE),
    frac_Plasma_Tcell_cohigh = n_Plasma_Tcell_cohigh / n_spots,
    .groups = "drop"
  )

safe_write_csv(
  category_tls_counts,
  file.path(out_table_dir, "18F_TLS_like_operational_cohigh_counts_by_dual_high_category.csv")
)

safe_fisher_or <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
  
  if (is.null(ft)) {
    return(data.frame(OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_))
  }
  
  data.frame(
    OR = unname(ft$estimate),
    CI_low = unname(ft$conf.int[1]),
    CI_high = unname(ft$conf.int[2]),
    p_value = ft$p.value
  )
}

cohigh_vars <- c(
  "TLSCore_high",
  "BCell_high",
  "TCellZone_high",
  "Operational_TLS_like_cohigh",
  "Humoral_Tcell_cohigh",
  "Plasma_Tcell_cohigh"
)

cohigh_or <- list()

for (vv in cohigh_vars) {
  dat <- tls_df %>%
    filter(Raw_DualHighCategory %in% c("Both-High", "Neither-High"))
  
  both_group <- dat$Raw_DualHighCategory == "Both-High"
  positive <- dat[[vv]]
  
  a <- sum(both_group & positive, na.rm = TRUE)
  b <- sum(both_group & !positive, na.rm = TRUE)
  c <- sum(!both_group & positive, na.rm = TRUE)
  d <- sum(!both_group & !positive, na.rm = TRUE)
  
  tmp <- safe_fisher_or(a, b, c, d)
  tmp$OperationalFeature <- vv
  tmp$n_BothHigh_positive <- a
  tmp$n_BothHigh_negative <- b
  tmp$n_NeitherHigh_positive <- c
  tmp$n_NeitherHigh_negative <- d
  cohigh_or[[length(cohigh_or) + 1]] <- tmp
}

cohigh_or_tbl <- bind_rows(cohigh_or) %>%
  select(OperationalFeature, everything()) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

safe_write_csv(
  cohigh_or_tbl,
  file.path(out_table_dir, "18F_TLS_like_operational_cohigh_enrichment_BothHigh_vs_NeitherHigh.csv")
)

p_tls_cohigh_map <- tls_df %>%
  mutate(
    TLS_like_category = case_when(
      Operational_TLS_like_cohigh ~ "TLS-core/B-cell/T-cell co-high",
      Humoral_Tcell_cohigh ~ "Humoral/T-cell co-high",
      TLSCore_high ~ "TLS-core high only/other",
      IGHumoral_high ~ "Humoral high only/other",
      TRUE ~ "Other"
    ),
    TLS_like_category = factor(
      TLS_like_category,
      levels = c(
        "Other",
        "Humoral high only/other",
        "TLS-core high only/other",
        "Humoral/T-cell co-high",
        "TLS-core/B-cell/T-cell co-high"
      )
    )
  ) %>%
  ggplot(aes(x = spatial_x, y = spatial_y, color = TLS_like_category)) +
  geom_point(size = 0.95, alpha = 0.9) +
  scale_y_reverse() +
  coord_equal() +
  labs(
    title = "Operational humoral/TLS-like co-high spatial features",
    subtitle = "Top-quartile co-high categories; not mature TLS calling",
    x = "Spatial x",
    y = "Spatial y",
    color = "Operational feature"
  ) +
  theme_revision(base_size = 12.0)

safe_ggsave(
  file.path(out_fig_dir, "18F_operational_TLS_like_cohigh_spatial_map.pdf"),
  p_tls_cohigh_map,
  width = 8.5,
  height = 6
)

safe_ggsave(
  file.path(out_fig_dir, "18F_operational_TLS_like_cohigh_spatial_map.png"),
  p_tls_cohigh_map,
  width = 8.5,
  height = 6
)

############################################################
## 12. Interpretation helper
############################################################

interpretation_summary <- data.frame(
  Item = c(
    "Purpose",
    "What can be concluded",
    "What cannot be concluded",
    "If TLS-like modules are enriched in Both-high spots",
    "If TLS-core/B-cell/T-cell co-high features are focal",
    "If humoral/plasma signals are diffuse",
    "Recommended wording"
  ),
  RecommendedInterpretation = c(
    "Evaluate whether humoral/plasma-cell signals in dual-high regions show TLS-like transcriptional and spatial context.",
    "Humoral/plasma-cell and TLS-like transcriptional features can be described as enriched, depleted, diffuse, or partially co-localized within the analyzed Visium section.",
    "These analyses cannot identify mature TLS architecture, germinal centers, follicular dendritic-cell networks, or HEV structures without histology/protein validation.",
    "Describe humoral/TLS-like enrichment as a feature of the mixed-cell dual-high region, not as proof of favorable TLS-associated immunity.",
    "Use cautious language such as partial TLS-like spatial organization or operational TLS-like co-high features.",
    "Interpret IGKC/IGHG1/plasma-cell signals as unresolved or incompletely organized humoral features.",
    "The humoral signal should be interpreted as a TLS-like/humoral transcriptional context within the dual-high region rather than as mature TLS or as a definitive favorable/resistance-promoting feature."
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  interpretation_summary,
  file.path(out_table_dir, "18F_interpretation_summary_for_manuscript_revision.csv")
)

############################################################
## 13. Output inventory and session info
############################################################

output_inventory <- data.frame(
  Type = c(
    "table", "table", "table", "table", "table",
    "table", "table", "table", "table", "table",
    "figure", "figure", "figure", "figure", "figure"
  ),
  File = c(
    file.path(out_table_dir, "18F_selected_expression_source_summary.csv"),
    file.path(out_table_dir, "18F_humoral_TLS_like_marker_gene_presence_audit.csv"),
    file.path(out_table_dir, "18F_humoral_TLS_like_module_scores.csv"),
    file.path(out_table_dir, "18F_TLS_like_module_summary_by_raw_dual_high_category.csv"),
    file.path(out_table_dir, "18F_TLS_like_BothHigh_vs_NeitherHigh_Wilcoxon.csv"),
    file.path(out_table_dir, "18F_TLS_like_module_correlations_with_states_and_composition.csv"),
    file.path(out_table_dir, "18F_TLS_like_module_univariate_spatial_Moran.csv"),
    file.path(out_table_dir, "18F_TLS_like_bivariate_spatial_colocalization_Moran.csv"),
    file.path(out_table_dir, "18F_TLS_like_operational_cohigh_counts_by_dual_high_category.csv"),
    file.path(out_table_dir, "18F_TLS_like_operational_cohigh_enrichment_BothHigh_vs_NeitherHigh.csv"),
    file.path(out_fig_dir, "18F_TLS_like_modules_by_raw_dual_high_category.pdf"),
    file.path(out_fig_dir, "18F_TLS_like_module_correlation_heatmap.pdf"),
    file.path(out_fig_dir, "18F_spatial_maps_of_TLS_like_modules.pdf"),
    file.path(out_fig_dir, "18F_TLS_like_univariate_spatial_Moran.pdf"),
    file.path(out_fig_dir, "18F_operational_TLS_like_cohigh_spatial_map.pdf")
  ),
  stringsAsFactors = FALSE
)

output_inventory$Exists <- file.exists(output_inventory$File)

safe_write_csv(
  output_inventory,
  file.path(out_table_dir, "18F_output_inventory.csv")
)

sink(file.path(log_dir, "sessionInfo_18F_humoral_TLS_like_spatial_interpretation_revision.txt"))
print(sessionInfo())
sink()

message("============================================================")
message("18F humoral/TLS-like spatial interpretation completed.")
message("Main table 1: ", file.path(out_table_dir, "18F_TLS_like_BothHigh_vs_NeitherHigh_Wilcoxon.csv"))
message("Main table 2: ", file.path(out_table_dir, "18F_TLS_like_module_correlations_with_states_and_composition.csv"))
message("Main table 3: ", file.path(out_table_dir, "18F_TLS_like_module_univariate_spatial_Moran.csv"))
message("Main table 4: ", file.path(out_table_dir, "18F_TLS_like_bivariate_spatial_colocalization_Moran.csv"))
message("Main table 5: ", file.path(out_table_dir, "18F_TLS_like_operational_cohigh_enrichment_BothHigh_vs_NeitherHigh.csv"))
message("============================================================")
## Public sequential end gate
.step06_out <- file.path(out_table_dir, "18F_humoral_TLS_like_module_scores.csv")
if (!file.exists(.step06_out)) stop("STEP 06 output missing: ", .step06_out, call. = FALSE)
cat("\nSTEP 06 PASS — humoral/TLS-like context rebuilt for S19\n")
