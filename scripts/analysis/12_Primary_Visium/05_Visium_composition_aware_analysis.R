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
## 18E_Visium_composition_aware_spatial_analysis_revision.R
##
## Purpose:
##   Revision analysis for Reviewer 3 mixed-spot / composition concern.
##
## Reviewer concern:
##   Visium spots contain mixed cell populations. Co-enrichment of
##   myeloid–Treg and stromal states within spots may reflect local
##   tissue composition rather than functional interaction.
##
## Main questions:
##   1) What composition marker modules characterize Both-high spots?
##   2) Is the dual-high enrichment attenuated after QC/composition residualization?
##   3) Should the dual-high result be interpreted as a mixed-cell tissue-level niche?
##
## Confirmed spot-level metadata/state input:
##   D:/ICB_resistance_project/results/tables/spatial_melanoma_validation/
##   Step11_spatial_analysis_metadata_with_coordinates_and_scores.csv
##
## Additional required input:
##   A Visium expression matrix, 10x matrix folder, 10x h5, or Seurat object
##   containing gene expression for the same 3,458 Visium spots.
##
## Optional manual input:
##   Sys.setenv(VISIUM_EXPR_FILE = "path/to/your/visium_expression_or_seurat.rds")
##
## Outputs:
##   results/tables/revision_visium_composition/
##   results/figures/revision_visium_composition/
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
intermediate_dir   <- file.path(project_dir, "results/intermediate")
log_dir            <- file.path(project_dir, "logs")

out_table_dir <- file.path(table_dir, "revision_visium_composition")
out_fig_dir <- file.path(project_dir, "results", "diagnostics", "12_Primary_Visium", "05_composition_aware")

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
  "Matrix"
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
## 3. Marker modules for composition-aware interpretation
############################################################

composition_modules <- list(
  CAF_Stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "PDGFRA"),
  Myeloid = c("LYZ", "CD68", "CSF1R", "LST1", "AIF1"),
  T_NK = c("CD3D", "CD3E", "CD8A", "NKG7", "GZMB"),
  B_Plasma = c("MS4A1", "CD79A", "MZB1", "JCHAIN", "IGKC"),
  Melanoma_Lineage = c("MLANA", "PMEL", "TYR", "MITF", "SOX10"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ENG")
)

composition_cols <- names(composition_modules)

############################################################
## 4. Load confirmed Visium spot-level state-score table
############################################################

visium_meta_file <- file.path(
  project_dir,
  "results/tables/spatial_melanoma_validation/Step11_spatial_analysis_metadata_with_coordinates_and_scores.csv"
)

if (!file.exists(visium_meta_file)) {
  stop("Confirmed Visium metadata/state-score file not found: ", visium_meta_file)
}

message("Reading Visium metadata/state-score table: ", visium_meta_file)

visium_df <- data.table::fread(
  visium_meta_file,
  data.table = FALSE,
  check.names = FALSE
)

colnames(visium_df) <- normalize_colnames(colnames(visium_df))

if ("SpotID" %in% colnames(visium_df) && !"Spot" %in% colnames(visium_df)) {
  colnames(visium_df)[colnames(visium_df) == "SpotID"] <- "Spot"
}

if (!"Spot" %in% colnames(visium_df)) {
  visium_df$Spot <- paste0("spot_", seq_len(nrow(visium_df)))
}

if ("nFeature_Spatial_audit" %in% colnames(visium_df) && !"nFeature_Spatial" %in% colnames(visium_df)) {
  colnames(visium_df)[colnames(visium_df) == "nFeature_Spatial_audit"] <- "nFeature_Spatial"
}

if ("nCount_Spatial_audit" %in% colnames(visium_df) && !"nCount_Spatial" %in% colnames(visium_df)) {
  colnames(visium_df)[colnames(visium_df) == "nCount_Spatial_audit"] <- "nCount_Spatial"
}

if ("percent_mt_audit" %in% colnames(visium_df) && !"percent_mt" %in% colnames(visium_df)) {
  colnames(visium_df)[colnames(visium_df) == "percent_mt_audit"] <- "percent_mt"
}

required_numeric_cols <- c(
  "spatial_x",
  "spatial_y",
  "nFeature_Spatial",
  "nCount_Spatial",
  "percent_mt",
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

for (cc in intersect(required_numeric_cols, colnames(visium_df))) {
  visium_df[[cc]] <- suppressWarnings(as.numeric(visium_df[[cc]]))
}

required_cols <- c(
  "Spot",
  "spatial_x",
  "spatial_y",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling"
)

missing_cols <- setdiff(required_cols, colnames(visium_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in Visium metadata table: ", paste(missing_cols, collapse = ", "))
}

visium_df <- visium_df %>%
  filter(
    is.finite(spatial_x),
    is.finite(spatial_y),
    is.finite(Myeloid_Treg_Immunosuppressive),
    is.finite(Tumor_dedifferentiation_Stromal_remodeling)
  )

message("Visium metadata retained: ", nrow(visium_df), " spots")

safe_write_csv(
  data.frame(
    InputFile = visium_meta_file,
    n_spots = nrow(visium_df),
    n_cols = ncol(visium_df),
    stringsAsFactors = FALSE
  ),
  file.path(out_table_dir, "18E_visium_metadata_input_summary.csv")
)

############################################################
## 5. Seurat v4/v5 expression extraction helpers
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
  if (length(preferred_assays) == 0) {
    preferred_assays <- assay_names
  }
  
  assay_use <- preferred_assays[1]
  message("Using Seurat assay: ", assay_use)
  
  mat <- NULL
  
  ## Seurat v5 / Assay5 preferred route.
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
        
        ## First try LayerData on assay object.
        m1 <- tryCatch(
          SeuratObject::LayerData(object = assay_obj, layer = layer_use),
          error = function(e1) {
            message("LayerData(assay_obj) failed: ", paste(as.character(e1$message), collapse = " | "))
            
            ## Then try LayerData on Seurat object with assay.
            tryCatch(
              SeuratObject::LayerData(object = obj, assay = assay_use, layer = layer_use),
              error = function(e2) {
                message("LayerData(Seurat object) failed: ", paste(as.character(e2$message), collapse = " | "))
                NULL
              }
            )
          }
        )
        
        if (is.null(m1)) {
          stop("LayerData extraction returned NULL.")
        }
        
        m1
      },
      error = function(e_layer) {
        message("LayerData route failed: ", paste(as.character(e_layer$message), collapse = " | "))
        
        ## Fallback: GetAssayData with layer argument, not slot.
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
  
  ## Seurat package fallback.
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
  
  if (is.null(mat)) {
    stop("Failed to extract expression matrix from Seurat object.")
  }
  
  mat <- as.matrix(mat)
  collapse_duplicate_genes(mat)
}

############################################################
## 6. Locate and read Visium expression matrix
############################################################

read_expression_table <- function(file) {
  message("Trying expression table: ", file)
  
  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    obj <- readRDS(file)
    
    ## Seurat object.
    if ("Seurat" %in% class(obj)) {
      message("Detected Seurat object.")
      return(extract_seurat_expression_matrix(obj))
    }
    
    ## Matrix or data.frame.
    if (is.matrix(obj) || is.data.frame(obj)) {
      mat <- as.matrix(obj)
      if (is.null(rownames(mat))) {
        stop("RDS matrix/data.frame has no rownames.")
      }
      storage.mode(mat) <- "numeric"
      return(collapse_duplicate_genes(mat))
    }
    
    ## List object: choose matrix/data.frame with best spot overlap.
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
  
  ## CSV / TSV / TXT expression matrix.
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
  
  ## 10x H5.
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

score_expression_candidate <- function(file_or_folder, is_folder = FALSE) {
  out <- data.frame(
    Candidate = file_or_folder,
    IsFolder = is_folder,
    Status = NA_character_,
    n_genes = NA_integer_,
    n_spots = NA_integer_,
    n_spot_overlap_exact = NA_integer_,
    n_spot_overlap_key = NA_integer_,
    Has_CAF_genes = NA_integer_,
    Has_Myeloid_genes = NA_integer_,
    Has_TNK_genes = NA_integer_,
    Has_BPlasma_genes = NA_integer_,
    Has_Melanoma_genes = NA_integer_,
    Has_Endothelial_genes = NA_integer_,
    Score = NA_real_,
    Error = NA_character_,
    stringsAsFactors = FALSE
  )
  
  expr <- tryCatch(
    {
      if (is_folder) {
        read_10x_mtx_folder(file_or_folder)
      } else {
        read_expression_table(file_or_folder)
      }
    },
    error = function(e) {
      out$Status <<- "read_error"
      out$Error <<- paste(as.character(e$message), collapse = " | ")
      NULL
    }
  )
  
  if (is.null(expr)) {
    return(out)
  }
  
  spot_overlap_exact <- sum(visium_df$Spot %in% colnames(expr))
  spot_overlap_key <- sum(spot_key(visium_df$Spot) %in% spot_key(colnames(expr)))
  
  gene_counts <- sapply(
    composition_modules,
    function(g) sum(clean_gene_symbol(g) %in% rownames(expr))
  )
  
  out$Status <- "read_ok"
  out$n_genes <- nrow(expr)
  out$n_spots <- ncol(expr)
  out$n_spot_overlap_exact <- spot_overlap_exact
  out$n_spot_overlap_key <- spot_overlap_key
  out$Has_CAF_genes <- gene_counts[["CAF_Stromal"]]
  out$Has_Myeloid_genes <- gene_counts[["Myeloid"]]
  out$Has_TNK_genes <- gene_counts[["T_NK"]]
  out$Has_BPlasma_genes <- gene_counts[["B_Plasma"]]
  out$Has_Melanoma_genes <- gene_counts[["Melanoma_Lineage"]]
  out$Has_Endothelial_genes <- gene_counts[["Endothelial"]]
  
  out$Score <-
    min(spot_overlap_key, nrow(visium_df)) * 10 +
    sum(gene_counts) * 5 +
    ifelse(nrow(expr) > 1000, 100, 0) +
    ifelse(ncol(expr) >= 3000 && ncol(expr) <= 5000, 100, 0)
  
  out
}

locate_visium_expression <- function() {
  manual <- Sys.getenv("VISIUM_EXPR_FILE")
  if (manual != "" && file.exists(manual)) {
    message("Using manual VISIUM_EXPR_FILE: ", manual)
    expr <- read_expression_table(manual)
    return(list(source = manual, expr = expr, inventory = data.frame()))
  }
  
  candidate_files <- unique(c(
    find_files_recursive(data_processed_dir, "Visium.*\\.(rds|csv|tsv|txt|h5)$"),
    find_files_recursive(data_processed_dir, "spatial.*\\.(rds|csv|tsv|txt|h5)$"),
    find_files_recursive(intermediate_dir, "Visium.*\\.(rds|csv|tsv|txt|h5)$"),
    find_files_recursive(intermediate_dir, "spatial.*\\.(rds|csv|tsv|txt|h5)$"),
    find_files_recursive(file.path(project_dir, "results"), "Visium.*\\.(rds|csv|tsv|txt|h5)$"),
    find_files_recursive(file.path(project_dir, "results"), "spatial.*\\.(rds|csv|tsv|txt|h5)$"),
    find_files_recursive(project_dir, "filtered_feature_bc_matrix\\.h5$"),
    find_files_recursive(project_dir, "raw_feature_bc_matrix\\.h5$")
  ))
  
  candidate_files <- candidate_files[
    !grepl(
      "18E_|18D_|18S_|18A_|candidate|inventory|summary|correlation|AUC|Moran|threshold|output|sessionInfo",
      basename(candidate_files),
      ignore.case = TRUE
    )
  ]
  
  candidate_files <- candidate_files[
    !normalizePath(candidate_files, winslash = "/", mustWork = FALSE) %in%
      normalizePath(visium_meta_file, winslash = "/", mustWork = FALSE)
  ]
  
  matrix_files <- find_files_recursive(project_dir, "^matrix\\.mtx(\\.gz)?$")
  candidate_folders <- unique(dirname(matrix_files))
  
  inventory_list <- list()
  
  for (ff in candidate_files) {
    message("Scoring candidate file: ", ff)
    inventory_list[[length(inventory_list) + 1]] <- score_expression_candidate(ff, is_folder = FALSE)
  }
  
  for (dd in candidate_folders) {
    message("Scoring candidate 10x folder: ", dd)
    inventory_list[[length(inventory_list) + 1]] <- score_expression_candidate(dd, is_folder = TRUE)
  }
  
  if (length(inventory_list) == 0) {
    inventory <- data.frame()
    safe_write_csv(inventory, file.path(out_table_dir, "18E_visium_expression_candidate_inventory_EMPTY.csv"))
    stop(
      "No candidate Visium expression source was found. ",
      "Set Sys.setenv(VISIUM_EXPR_FILE='path/to/expression_or_seurat.rds') and rerun."
    )
  }
  
  inventory <- bind_rows(inventory_list) %>%
    arrange(desc(Score))
  
  safe_write_csv(
    inventory,
    file.path(out_table_dir, "18E_visium_expression_candidate_inventory.csv")
  )
  
  best <- inventory %>%
    filter(Status == "read_ok") %>%
    arrange(desc(Score)) %>%
    slice(1)
  
  if (nrow(best) == 0 || is.na(best$Score) || best$n_spot_overlap_key < 100) {
    stop(
      "No usable Visium expression source found. ",
      "Check 18E_visium_expression_candidate_inventory.csv, or set VISIUM_EXPR_FILE manually."
    )
  }
  
  source <- best$Candidate[1]
  is_folder <- best$IsFolder[1]
  
  message("Selected Visium expression source: ", source)
  
  expr <- if (is_folder) {
    read_10x_mtx_folder(source)
  } else {
    read_expression_table(source)
  }
  
  list(source = source, expr = expr, inventory = inventory)
}

expr_loc <- locate_visium_expression()
visium_expr <- expr_loc$expr
visium_expr_source <- expr_loc$source

safe_write_csv(
  data.frame(
    SelectedExpressionSource = visium_expr_source,
    n_genes = nrow(visium_expr),
    n_expression_spots = ncol(visium_expr),
    n_meta_spots = nrow(visium_df),
    exact_spot_overlap = sum(visium_df$Spot %in% colnames(visium_expr)),
    key_spot_overlap = sum(spot_key(visium_df$Spot) %in% spot_key(colnames(visium_expr))),
    stringsAsFactors = FALSE
  ),
  file.path(out_table_dir, "18E_selected_visium_expression_source_summary.csv")
)

############################################################
## 7. Align expression to metadata spots
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
  file.path(out_table_dir, "18E_visium_expression_metadata_spot_matching_audit.csv")
)

message("Aligned Visium expression: ", nrow(visium_expr_aligned), " genes x ", ncol(visium_expr_aligned), " spots")

############################################################
## 8. Compute Visium marker-module composition scores
############################################################

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

comp <- compute_module_scores(visium_expr_aligned, composition_modules)

safe_write_csv(
  comp$scores,
  file.path(out_table_dir, "18E_Visium_composition_marker_module_scores.csv")
)

safe_write_csv(
  comp$gene_audit,
  file.path(out_table_dir, "18E_Visium_composition_marker_gene_presence_audit.csv")
)

visium_comp_df <- visium_df_aligned %>%
  left_join(comp$scores, by = "Spot")

safe_write_csv(
  visium_comp_df,
  file.path(out_table_dir, "18E_Visium_spot_state_QC_composition_table.csv")
)

safe_save_rds(
  visium_comp_df,
  file.path(data_processed_dir, "18E_Visium_spot_state_QC_composition_table.rds")
)

############################################################
## 9. Define raw high-state categories
############################################################

x_col <- "Myeloid_Treg_Immunosuppressive"
y_col <- "Tumor_dedifferentiation_Stromal_remodeling"

assign_dual_category <- function(df, x_col, y_col, cutoff_quantile = 0.75,
                                 prefix = "") {
  x_cut <- quantile(df[[x_col]], probs = cutoff_quantile, na.rm = TRUE, names = FALSE)
  y_cut <- quantile(df[[y_col]], probs = cutoff_quantile, na.rm = TRUE, names = FALSE)
  
  x_high <- df[[x_col]] >= x_cut
  y_high <- df[[y_col]] >= y_cut
  
  cat <- ifelse(
    x_high & y_high, "Both-High",
    ifelse(
      x_high & !y_high, "Myeloid-Treg-High Only",
      ifelse(!x_high & y_high, "Dediff/Stromal-High Only", "Neither-High")
    )
  )
  
  df[[paste0(prefix, "MyeloidTreg_high")]] <- x_high
  df[[paste0(prefix, "DediffStromal_high")]] <- y_high
  df[[paste0(prefix, "DualHighCategory")]] <- factor(
    cat,
    levels = c(
      "Neither-High",
      "Myeloid-Treg-High Only",
      "Dediff/Stromal-High Only",
      "Both-High"
    )
  )
  
  attr(df, paste0(prefix, "cutoffs")) <- c(MyeloidTreg = x_cut, DediffStromal = y_cut)
  df
}

visium_comp_df <- assign_dual_category(
  visium_comp_df,
  x_col = x_col,
  y_col = y_col,
  cutoff_quantile = 0.75,
  prefix = "Raw_"
)

raw_category_counts <- as.data.frame(table(visium_comp_df$Raw_DualHighCategory))
colnames(raw_category_counts) <- c("Raw_DualHighCategory", "n_spots")

safe_write_csv(
  raw_category_counts,
  file.path(out_table_dir, "18E_raw_dual_high_category_counts.csv")
)

############################################################
## 10. Composition profile across raw dual-high categories
############################################################

composition_long <- visium_comp_df %>%
  select(Spot, Raw_DualHighCategory, all_of(composition_cols)) %>%
  pivot_longer(
    cols = all_of(composition_cols),
    names_to = "CompositionModule",
    values_to = "ModuleScore"
  )

composition_group_summary <- composition_long %>%
  group_by(CompositionModule, Raw_DualHighCategory) %>%
  summarise(
    n = sum(is.finite(ModuleScore)),
    mean = mean(ModuleScore, na.rm = TRUE),
    median = median(ModuleScore, na.rm = TRUE),
    IQR_low = quantile(ModuleScore, 0.25, na.rm = TRUE),
    IQR_high = quantile(ModuleScore, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

safe_write_csv(
  composition_group_summary,
  file.path(out_table_dir, "18E_composition_module_summary_by_raw_dual_high_category.csv")
)

kw_tests <- composition_long %>%
  group_by(CompositionModule) %>%
  summarise(
    n = sum(is.finite(ModuleScore)),
    kruskal_p = tryCatch(
      kruskal.test(ModuleScore ~ Raw_DualHighCategory)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(kruskal_p_adj_BH = p.adjust(kruskal_p, method = "BH"))

safe_write_csv(
  kw_tests,
  file.path(out_table_dir, "18E_composition_module_Kruskal_tests_across_categories.csv")
)

both_vs_neither <- composition_long %>%
  filter(Raw_DualHighCategory %in% c("Both-High", "Neither-High")) %>%
  group_by(CompositionModule) %>%
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
  both_vs_neither,
  file.path(out_table_dir, "18E_composition_BothHigh_vs_NeitherHigh_Wilcoxon.csv")
)

p_comp_box <- composition_long %>%
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
  facet_wrap(~ CompositionModule, scales = "free_y", ncol = 3) +
  labs(
    title = "Composition marker modules across Visium high-state categories",
    subtitle = "Marker-module scores provide composition-aware context for mixed-cell Visium spots",
    x = "Raw high-state category",
    y = "Marker-module score"
  ) +
  theme_revision(base_size = 11.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25))

safe_ggsave(
  file.path(out_fig_dir, "18E_composition_modules_by_raw_dual_high_category.pdf"),
  p_comp_box,
  width = 11,
  height = 7
)

safe_ggsave(
  file.path(out_fig_dir, "18E_composition_modules_by_raw_dual_high_category.png"),
  p_comp_box,
  width = 11,
  height = 7
)

############################################################
## 11. Spatial maps of composition modules
############################################################

comp_map_long <- visium_comp_df %>%
  select(Spot, spatial_x, spatial_y, all_of(composition_cols)) %>%
  pivot_longer(
    cols = all_of(composition_cols),
    names_to = "CompositionModule",
    values_to = "ModuleScore"
  )

p_comp_map <- ggplot(
  comp_map_long,
  aes(x = spatial_x, y = spatial_y, color = ModuleScore)
) +
  geom_point(size = 0.9, alpha = 0.9) +
  scale_y_reverse() +
  coord_equal() +
  facet_wrap(~ CompositionModule, ncol = 3) +
  scale_color_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    na.value = "grey85"
  ) +
  labs(
    title = "Spatial distribution of Visium composition marker modules",
    subtitle = "Used for composition-aware interpretation of mixed-cell spots",
    x = "Spatial x",
    y = "Spatial y",
    color = "Module\nscore"
  ) +
  theme_revision(base_size = 11.5)

safe_ggsave(
  file.path(out_fig_dir, "18E_spatial_maps_of_composition_modules.pdf"),
  p_comp_map,
  width = 10.5,
  height = 7
)

safe_ggsave(
  file.path(out_fig_dir, "18E_spatial_maps_of_composition_modules.png"),
  p_comp_map,
  width = 10.5,
  height = 7
)

############################################################
## 12. State-composition correlations in Visium
############################################################

state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_display_public <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

cor_test_pair <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  
  if (n < 4) {
    return(data.frame(n = n, estimate = NA_real_, p_value = NA_real_))
  }
  
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  
  data.frame(
    n = n,
    estimate = unname(ct$estimate),
    p_value = ct$p.value
  )
}

visium_state_comp_cor <- list()

for (st in intersect(state_cols, colnames(visium_comp_df))) {
  for (cc in composition_cols) {
    tmp <- cor_test_pair(visium_comp_df[[st]], visium_comp_df[[cc]], method = "spearman")
    tmp$State <- st
    tmp$CompositionModule <- cc
    visium_state_comp_cor[[length(visium_state_comp_cor) + 1]] <- tmp
  }
}

visium_state_comp_cor <- bind_rows(visium_state_comp_cor) %>%
  select(State, CompositionModule, n, rho = estimate, p_value) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

safe_write_csv(
  visium_state_comp_cor,
  file.path(out_table_dir, "18E_Visium_state_composition_Spearman_correlation.csv")
)

p_state_comp_heat <- visium_state_comp_cor %>%
  mutate(
    State = factor(State, levels = state_cols),
    CompositionModule = factor(CompositionModule, levels = composition_cols),
    label = sprintf("%.2f", rho)
  ) %>%
  ggplot(aes(x = CompositionModule, y = State, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  labs(
    title = "Visium state-score associations with composition marker modules",
    subtitle = "Spearman correlations across spatial spots",
    x = "Composition marker module",
    y = "Tumor–immune state",
    fill = "Spearman\nrho"
  ) +
  scale_y_discrete(labels = state_display_public) +
  theme_revision(base_size = 11.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 10.25))

safe_ggsave(
  file.path(out_fig_dir, "18E_Visium_state_composition_correlation_heatmap.pdf"),
  p_state_comp_heat,
  width = 8.5,
  height = 5.5
)

safe_ggsave(
  file.path(out_fig_dir, "18E_Visium_state_composition_correlation_heatmap.png"),
  p_state_comp_heat,
  width = 8.5,
  height = 5.5
)

############################################################
## 13. Composition-aware residualization of state scores
############################################################

add_residualized_score <- function(df, outcome_col, covariates, new_col) {
  out <- df
  vars <- c(outcome_col, covariates)
  vars <- vars[vars %in% colnames(out)]
  
  if (!outcome_col %in% vars) {
    out[[new_col]] <- NA_real_
    return(out)
  }
  
  covariates2 <- setdiff(vars, outcome_col)
  
  complete <- complete.cases(out[, vars, drop = FALSE])
  n_complete <- sum(complete)
  
  resid_vec <- rep(NA_real_, nrow(out))
  
  if (length(covariates2) == 0 || n_complete < length(covariates2) + 5) {
    warning("Residualization skipped/unstable for ", new_col)
    out[[new_col]] <- resid_vec
    return(out)
  }
  
  fml <- as.formula(paste(outcome_col, "~", paste(covariates2, collapse = " + ")))
  fit <- lm(fml, data = out[complete, , drop = FALSE])
  resid_vec[complete] <- zvec(residuals(fit))
  
  out[[new_col]] <- resid_vec
  out
}

qc_covariates <- intersect(c("nCount_Spatial", "nFeature_Spatial", "percent_mt"), colnames(visium_comp_df))
comp_covariates <- composition_cols
qc_comp_covariates <- c(qc_covariates, comp_covariates)

visium_comp_df <- visium_comp_df %>%
  add_residualized_score(
    outcome_col = x_col,
    covariates = qc_covariates,
    new_col = "MyeloidTreg_resid_QC"
  ) %>%
  add_residualized_score(
    outcome_col = y_col,
    covariates = qc_covariates,
    new_col = "DediffStromal_resid_QC"
  ) %>%
  add_residualized_score(
    outcome_col = x_col,
    covariates = comp_covariates,
    new_col = "MyeloidTreg_resid_Composition"
  ) %>%
  add_residualized_score(
    outcome_col = y_col,
    covariates = comp_covariates,
    new_col = "DediffStromal_resid_Composition"
  ) %>%
  add_residualized_score(
    outcome_col = x_col,
    covariates = qc_comp_covariates,
    new_col = "MyeloidTreg_resid_QC_Composition"
  ) %>%
  add_residualized_score(
    outcome_col = y_col,
    covariates = qc_comp_covariates,
    new_col = "DediffStromal_resid_QC_Composition"
  )

safe_write_csv(
  visium_comp_df,
  file.path(out_table_dir, "18E_Visium_spot_state_composition_residualized_scores.csv")
)

############################################################
## 14. Raw vs residualized dual-high enrichment
############################################################

safe_fisher_or <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
  
  if (is.null(ft)) {
    return(data.frame(
      OR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p_value = NA_real_
    ))
  }
  
  data.frame(
    OR = unname(ft$estimate),
    CI_low = unname(ft$conf.int[1]),
    CI_high = unname(ft$conf.int[2]),
    p_value = ft$p.value
  )
}

dual_or_for_scores <- function(df, x_score, y_score, label, cutoff_quantile = 0.75) {
  dat <- df %>%
    filter(is.finite(.data[[x_score]]), is.finite(.data[[y_score]]))
  
  x_cut <- quantile(dat[[x_score]], probs = cutoff_quantile, na.rm = TRUE, names = FALSE)
  y_cut <- quantile(dat[[y_score]], probs = cutoff_quantile, na.rm = TRUE, names = FALSE)
  
  x_high <- dat[[x_score]] >= x_cut
  y_high <- dat[[y_score]] >= y_cut
  
  both <- sum(x_high & y_high)
  x_only <- sum(x_high & !y_high)
  y_only <- sum(!x_high & y_high)
  neither <- sum(!x_high & !y_high)
  
  or_tbl <- safe_fisher_or(
    a = both,
    b = x_only,
    c = y_only,
    d = neither
  )
  
  cbind(
    data.frame(
      Analysis = label,
      x_score = x_score,
      y_score = y_score,
      n_spots = nrow(dat),
      cutoff_quantile = cutoff_quantile,
      x_cutoff = x_cut,
      y_cutoff = y_cut,
      BothHigh = both,
      MyeloidTregHighOnly = x_only,
      DediffStromalHighOnly = y_only,
      NeitherHigh = neither,
      stringsAsFactors = FALSE
    ),
    or_tbl
  )
}

dual_or_tbl <- bind_rows(
  dual_or_for_scores(
    visium_comp_df,
    x_score = x_col,
    y_score = y_col,
    label = "Raw state scores"
  ),
  dual_or_for_scores(
    visium_comp_df,
    x_score = "MyeloidTreg_resid_QC",
    y_score = "DediffStromal_resid_QC",
    label = "QC-residualized state scores"
  ),
  dual_or_for_scores(
    visium_comp_df,
    x_score = "MyeloidTreg_resid_Composition",
    y_score = "DediffStromal_resid_Composition",
    label = "Composition-residualized state scores"
  ),
  dual_or_for_scores(
    visium_comp_df,
    x_score = "MyeloidTreg_resid_QC_Composition",
    y_score = "DediffStromal_resid_QC_Composition",
    label = "QC/composition-residualized state scores"
  )
)

dual_or_tbl <- dual_or_tbl %>%
  mutate(
    OR_attenuation_vs_raw = OR[Analysis == "Raw state scores"][1] - OR,
    log2OR = log2(OR)
  )

safe_write_csv(
  dual_or_tbl,
  file.path(out_table_dir, "18E_raw_vs_residualized_dual_high_Fisher_OR.csv")
)

p_or <- dual_or_tbl %>%
  mutate(Analysis = factor(Analysis, levels = rev(Analysis))) %>%
  ggplot(aes(x = Analysis, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.6) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.15, linewidth = 0.5) +
  coord_flip() +
  scale_y_log10() +
  labs(
    title = "Raw and residualized dual-high co-enrichment in Visium",
    subtitle = "Composition-aware sensitivity analysis; not functional validation",
    x = NULL,
    y = "Fisher OR, log10 scale"
  ) +
  theme_revision(base_size = 12.0)

safe_ggsave(
  file.path(out_fig_dir, "18E_raw_vs_residualized_dual_high_OR.pdf"),
  p_or,
  width = 8,
  height = 4.8
)

safe_ggsave(
  file.path(out_fig_dir, "18E_raw_vs_residualized_dual_high_OR.png"),
  p_or,
  width = 8,
  height = 4.8
)

############################################################
## 15. Raw vs residualized continuous correlation
############################################################

cor_raw_resid <- bind_rows(
  cor_test_pair(visium_comp_df[[x_col]], visium_comp_df[[y_col]], method = "spearman") %>%
    mutate(Analysis = "Raw state scores"),
  cor_test_pair(visium_comp_df$MyeloidTreg_resid_QC, visium_comp_df$DediffStromal_resid_QC, method = "spearman") %>%
    mutate(Analysis = "QC-residualized state scores"),
  cor_test_pair(visium_comp_df$MyeloidTreg_resid_Composition, visium_comp_df$DediffStromal_resid_Composition, method = "spearman") %>%
    mutate(Analysis = "Composition-residualized state scores"),
  cor_test_pair(visium_comp_df$MyeloidTreg_resid_QC_Composition, visium_comp_df$DediffStromal_resid_QC_Composition, method = "spearman") %>%
    mutate(Analysis = "QC/composition-residualized state scores")
) %>%
  select(Analysis, n, rho = estimate, p_value)

safe_write_csv(
  cor_raw_resid,
  file.path(out_table_dir, "18E_raw_vs_residualized_state_score_Spearman_correlation.csv")
)

p_cor_resid <- cor_raw_resid %>%
  mutate(Analysis = factor(Analysis, levels = rev(Analysis))) %>%
  ggplot(aes(x = Analysis, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 2.6) +
  coord_flip() +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    title = "Raw and residualized continuous state-score correlation",
    subtitle = "Spearman correlation between myeloid–Treg and\ntumor-dedifferentiation/stromal-remodeling scores",
    x = NULL,
    y = "Spearman rho"
  ) +
  theme_revision(base_size = 12.0) +
  theme(plot.subtitle = element_text(size = 10.5, lineheight = 1.0))

safe_ggsave(
  file.path(out_fig_dir, "18E_raw_vs_residualized_state_score_correlation.pdf"),
  p_cor_resid,
  width = 8,
  height = 4.8
)

safe_ggsave(
  file.path(out_fig_dir, "18E_raw_vs_residualized_state_score_correlation.png"),
  p_cor_resid,
  width = 8,
  height = 4.8
)

############################################################
## 16. Spatial maps of raw and residualized categories
############################################################

visium_comp_df <- assign_dual_category(
  visium_comp_df,
  x_col = "MyeloidTreg_resid_QC_Composition",
  y_col = "DediffStromal_resid_QC_Composition",
  cutoff_quantile = 0.75,
  prefix = "Residualized_"
)

map_cat_df <- bind_rows(
  visium_comp_df %>%
    transmute(
      Spot,
      spatial_x,
      spatial_y,
      Category = Raw_DualHighCategory,
      Analysis = "Raw state scores"
    ),
  visium_comp_df %>%
    transmute(
      Spot,
      spatial_x,
      spatial_y,
      Category = Residualized_DualHighCategory,
      Analysis = "QC/composition-residualized state scores"
    )
)

p_cat_map <- ggplot(
  map_cat_df,
  aes(x = spatial_x, y = spatial_y, color = Category)
) +
  geom_point(size = 1.0, alpha = 0.9) +
  scale_y_reverse() +
  coord_equal() +
  facet_wrap(~ Analysis, ncol = 2) +
  labs(
    title = "Raw and composition-residualized high-state categories",
    subtitle = "Top-quartile operational categories before and after QC/composition residualization",
    x = "Spatial x",
    y = "Spatial y",
    color = "Category"
  ) +
  theme_revision(base_size = 11.5) +
  theme(plot.subtitle = element_text(size = 10.5))

safe_ggsave(
  file.path(out_fig_dir, "18E_raw_vs_QCcomposition_residualized_category_maps.pdf"),
  p_cat_map,
  width = 10,
  height = 5
)

safe_ggsave(
  file.path(out_fig_dir, "18E_raw_vs_QCcomposition_residualized_category_maps.png"),
  p_cat_map,
  width = 10,
  height = 5
)

############################################################
## 17. Interpretation helper table
############################################################

interpretation_summary <- data.frame(
  Item = c(
    "Purpose",
    "Composition module interpretation",
    "Both-high composition comparison",
    "Residualized dual-high OR",
    "If residualized OR remains > 1",
    "If residualized OR is strongly attenuated",
    "Interpretation boundary"
  ),
  RecommendedInterpretation = c(
    "Evaluate whether Visium dual-high co-enrichment is partly explained by local mixed-cell composition.",
    "Marker-module scores provide composition-aware context for Visium spots; they are not formal cell-type deconvolution.",
    "Both-high spots can be described as enriched or depleted for CAF/stromal, myeloid, T/NK, B/plasma, melanoma-lineage, or endothelial marker modules.",
    "Raw and residualized Fisher ORs quantify whether the dual-high relationship persists after removing QC and marker-module composition effects.",
    "State that dual-high co-enrichment was not fully explained by QC/composition covariates, while still interpreting the result as mixed-cell spatial support.",
    "State that the raw dual-high pattern was substantially composition-associated and should be interpreted primarily as a mixed-cell tissue-level niche.",
    "These analyses do not establish functional interaction, causal signaling, or population-level spatial generalizability."
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  interpretation_summary,
  file.path(out_table_dir, "18E_interpretation_summary_for_manuscript_revision.csv")
)

############################################################
## 18. Output inventory and session info
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
    "figure",
    "figure"
  ),
  File = c(
    file.path(out_table_dir, "18E_selected_visium_expression_source_summary.csv"),
    file.path(out_table_dir, "18E_Visium_composition_marker_gene_presence_audit.csv"),
    file.path(out_table_dir, "18E_Visium_spot_state_QC_composition_table.csv"),
    file.path(out_table_dir, "18E_composition_module_summary_by_raw_dual_high_category.csv"),
    file.path(out_table_dir, "18E_composition_BothHigh_vs_NeitherHigh_Wilcoxon.csv"),
    file.path(out_table_dir, "18E_Visium_state_composition_Spearman_correlation.csv"),
    file.path(out_table_dir, "18E_raw_vs_residualized_dual_high_Fisher_OR.csv"),
    file.path(out_table_dir, "18E_raw_vs_residualized_state_score_Spearman_correlation.csv"),
    file.path(out_fig_dir, "18E_composition_modules_by_raw_dual_high_category.pdf"),
    file.path(out_fig_dir, "18E_spatial_maps_of_composition_modules.pdf"),
    file.path(out_fig_dir, "18E_Visium_state_composition_correlation_heatmap.pdf"),
    file.path(out_fig_dir, "18E_raw_vs_residualized_dual_high_OR.pdf"),
    file.path(out_fig_dir, "18E_raw_vs_QCcomposition_residualized_category_maps.pdf")
  ),
  stringsAsFactors = FALSE
)

output_inventory$Exists <- file.exists(output_inventory$File)

safe_write_csv(
  output_inventory,
  file.path(out_table_dir, "18E_output_inventory.csv")
)

sink(file.path(log_dir, "sessionInfo_18E_Visium_composition_aware_spatial_analysis_revision.txt"))
print(sessionInfo())
sink()

message("============================================================")
message("18E Visium composition-aware spatial analysis completed.")
message("Selected expression source: ", visium_expr_source)
message("Main table 1: ", file.path(out_table_dir, "18E_composition_BothHigh_vs_NeitherHigh_Wilcoxon.csv"))
message("Main table 2: ", file.path(out_table_dir, "18E_Visium_state_composition_Spearman_correlation.csv"))
message("Main table 3: ", file.path(out_table_dir, "18E_raw_vs_residualized_dual_high_Fisher_OR.csv"))
message("Main table 4: ", file.path(out_table_dir, "18E_raw_vs_residualized_state_score_Spearman_correlation.csv"))
message("============================================================")
## Public sequential end gate
.step05_out <- file.path(out_table_dir, "18E_Visium_spot_state_composition_residualized_scores.csv")
if (!file.exists(.step05_out)) stop("STEP 05 output missing: ", .step05_out, call. = FALSE)
cat("\nSTEP 05 PASS — composition-aware Visium analysis completed\n")
