############################################################
## 08b_GSE244983_malignant_subclustering.R
##
## v1.4.2 FINAL hardening:
##   - preserves the frozen malignant subclustering and Figure-5 compact programs
##   - adds the frozen five-program malignant analysis used by Figure 6/S10
##   - exports canonical Supplementary Figure S10A-D source CSVs directly
##   - no downstream S10 plotting script needs RDS fallback or fuzzy columns
##
## Purpose
## -------
## Re-cluster the canonical malignant-lineage compartment in GSE244983
## and generate analysis/source tables for Main Figure 5.
##
## Canonical analysis lineage
## --------------------------
## Script 07 freezes the primary GSE244983 major-cell annotation.
## Script 10 therefore consumes ONLY:
##
##   results/intermediate/GSE244983/
##     GSE244983_seurat_major_annotated.rds
##
## and subsets:
##
##   MajorCellType %in% c("Malignant", "Cycling malignant")
##
## The current canonical lineage contains 19,596 malignant-lineage cells:
##   - Malignant:          17,708
##   - Cycling malignant:  1,888
##
## The older manuscript/reporting value of 19,599 is NOT used as an
## analysis gate. Reporting harmonization is handled separately from the
## reproducible canonical analysis.
##
## Historical mother-script logic retained
## ---------------------------------------
##   - set.seed(20260430) before malignant-only reprocessing
##   - RNA assay
##   - 2,000 variable features (vst)
##   - ScaleData on variable features
##   - PCA: 40 components
##   - Harmony by Sample, dimensions 1:30
##   - FindNeighbors on Harmony dimensions 1:30
##   - FindClusters resolution = 0.6
##   - malignant labels Mal_0 ... Mal_9
##   - compact four-program marker/module analysis for Figure 5
##
## Public-pipeline hardening
## -------------------------
##   - no recursive file discovery
##   - no historical RDS fallback
##   - no package auto-install
##   - Harmony is required; no silent PCA fallback
##   - explicit cluster and UMAP random seeds matching Seurat defaults
##   - analysis/source tables only; no manuscript figures
##   - canonical object is saved only after all structural gates PASS
##   - input/output MD5s, cell-set MD5s and sessionInfo are recorded
############################################################

options(stringsAsFactors = FALSE)

SCRIPT_ID <- "10_GSE244983_malignant_subclustering"
SCRIPT_PATCH <- "v1.3_LF_hash_escape_fix"

############################################################
## 0. Project root
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")

if (!nzchar(project_dir)) {
  cwd <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )
  
  probe <- file.path(
    cwd,
    "results",
    "intermediate",
    "GSE244983",
    "GSE244983_seurat_major_annotated.rds"
  )
  
  if (file.exists(probe)) {
    project_dir <- cwd
  } else {
    stop(
      "Set ICB_PROJECT_DIR to the repository root, or run from the root ",
      "containing results/intermediate/GSE244983/",
      "GSE244983_seurat_major_annotated.rds."
    )
  }
}

project_dir <- normalizePath(
  project_dir,
  winslash = "/",
  mustWork = TRUE
)

message(
  "Project root: ",
  project_dir
)

############################################################
## 1. Packages -- NO auto-install
############################################################

required_pkgs <- c(
  "Seurat",
  "Matrix",
  "dplyr",
  "tidyr",
  "tibble",
  "harmony"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(
      missing_pkgs,
      collapse = ", "
    ),
    ". Restore the repository environment and rerun."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

############################################################
## 2. Canonical paths
############################################################

input_seurat_file <- file.path(
  project_dir,
  "results",
  "intermediate",
  "GSE244983",
  "GSE244983_seurat_major_annotated.rds"
)

out_table_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "GSE244983",
  "malignant_subclustering"
)

out_intermediate_dir <- file.path(
  project_dir,
  "results",
  "intermediate",
  "GSE244983"
)

log_dir <- file.path(
  project_dir,
  "logs"
)

dir.create(
  out_table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  out_intermediate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_seurat_file)) {
  stop(
    "Missing canonical script-07 input: ",
    input_seurat_file
  )
}

############################################################
## 3. Helpers
############################################################

safe_write_csv <- function(
    x,
    path,
    row.names = FALSE
) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  utils::write.csv(
    x,
    file = path,
    row.names = row.names,
    na = ""
  )
  
  message(
    "Saved: ",
    path
  )
  
  invisible(path)
}

safe_save_rds <- function(
    x,
    path
) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    x,
    file = path
  )
  
  message(
    "Saved: ",
    path
  )
  
  invisible(path)
}

md5_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  
  unname(
    tools::md5sum(path)
  )
}

cellset_md5 <- function(ids) {
  ## Cross-platform deterministic hash:
  ## sort unique IDs, encode as UTF-8, join with LF ("\n"),
  ## retain a final LF, and write bytes in binary mode.
  ##
  ## Do NOT use writeLines() here: on Windows it writes CRLF and
  ## therefore produces a different MD5 from Linux/macOS.
  ids <- sort(
    unique(
      enc2utf8(
        as.character(ids)
      )
    )
  )
  
  payload <- paste0(
    paste(
      ids,
      collapse = "\n"
    ),
    "\n"
  )
  
  tf <- tempfile(
    fileext = ".txt"
  )
  
  con <- file(
    tf,
    open = "wb"
  )
  
  on.exit(
    {
      try(
        close(con),
        silent = TRUE
      )
      unlink(tf)
    },
    add = TRUE
  )
  
  writeBin(
    charToRaw(payload),
    con
  )
  
  close(con)
  
  unname(
    tools::md5sum(tf)
  )
}

relative_to_project <- function(path) {
  p <- normalizePath(
    path,
    winslash = "/",
    mustWork = FALSE
  )
  
  root <- paste0(
    normalizePath(
      project_dir,
      winslash = "/",
      mustWork = TRUE
    ),
    "/"
  )
  
  if (startsWith(p, root)) {
    substring(
      p,
      nchar(root) + 1L
    )
  } else {
    p
  }
}

zscore_safe <- function(x) {
  x <- as.numeric(x)
  
  if (all(is.na(x))) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }
  
  s <- stats::sd(
    x,
    na.rm = TRUE
  )
  
  m <- mean(
    x,
    na.rm = TRUE
  )
  
  if (
    is.na(s) ||
    s == 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }
  
  (x - m) / s
}

get_normalized_rna <- function(obj) {
  if (!"RNA" %in% names(obj@assays)) {
    stop(
      "Canonical script-07 object does not contain an RNA assay."
    )
  }
  
  out <- tryCatch(
    Seurat::GetAssayData(
      obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e_layer) {
      tryCatch(
        Seurat::GetAssayData(
          obj,
          assay = "RNA",
          slot = "data"
        ),
        error = function(e_slot) {
          stop(
            "Normalized RNA data layer is unavailable. ",
            "Script 10 will not fall back to raw counts or silently ",
            "change normalization."
          )
        }
      )
    }
  )
  
  if (
    nrow(out) == 0 ||
    ncol(out) == 0
  ) {
    stop(
      "Normalized RNA data layer is empty."
    )
  }
  
  out
}

match_features_case_insensitive <- function(
    feature_names,
    requested_genes
) {
  feature_upper <- toupper(
    feature_names
  )
  
  result <- data.frame(
    RequestedGene = requested_genes,
    MatchedFeature = NA_character_,
    Present = FALSE,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(requested_genes)) {
    hit <- which(
      feature_upper ==
        toupper(
          requested_genes[[i]]
        )
    )
    
    if (length(hit) > 0) {
      result$MatchedFeature[[i]] <-
        feature_names[
          hit[[1]]
        ]
      
      result$Present[[i]] <- TRUE
    }
  }
  
  result
}

score_module_simple <- function(
    expr_data,
    requested_genes,
    cell_order,
    score_name
) {
  feature_match <- match_features_case_insensitive(
    rownames(expr_data),
    requested_genes
  )
  
  present_features <- feature_match$MatchedFeature[
    feature_match$Present
  ]
  
  if (length(present_features) < 2) {
    stop(
      "Too few genes present for compact malignant program ",
      score_name,
      ": ",
      length(present_features),
      "/",
      length(requested_genes)
    )
  }
  
  score <- Matrix::colMeans(
    expr_data[
      present_features,
      ,
      drop = FALSE
    ]
  )
  
  score <- score[
    cell_order
  ]
  
  list(
    score = as.numeric(score),
    presence = feature_match
  )
}

## Reproduce the numeric data used by Seurat::DotPlot without plotting.
##
## AverageExpressionLinear:
##   mean(expm1(log-normalized expression)) within a subcluster.
##
## PercentExpressed:
##   percentage of cells with normalized expression > 0.
##
## ScaledAverageExpression_DotPlot:
##   log1p(AverageExpressionLinear), z-standardized across subclusters
##   within each gene, then clipped to Seurat DotPlot's default
##   color range [-2.5, 2.5].
build_dotplot_source <- function(
    marker_expr_dense,
    marker_features,
    marker_labels,
    marker_programs,
    groups,
    group_levels
) {
  rows <- vector(
    "list",
    length(group_levels)
  )
  
  for (i in seq_along(group_levels)) {
    g <- group_levels[[i]]
    
    idx <- which(
      as.character(groups) ==
        g
    )
    
    if (length(idx) == 0) {
      next
    }
    
    sub_mat <- marker_expr_dense[
      ,
      idx,
      drop = FALSE
    ]
    
    avg_linear <- rowMeans(
      expm1(
        sub_mat
      )
    )
    
    pct_exp <- rowMeans(
      sub_mat >
        0
    ) *
      100
    
    rows[[i]] <- data.frame(
      MalignantSubcluster = g,
      N_cells = length(idx),
      Gene = marker_labels,
      MatchedFeature = marker_features,
      MarkerProgram = marker_programs,
      AverageExpressionLinear = avg_linear,
      AverageExpressionLog1p = log1p(
        avg_linear
      ),
      PercentExpressed = pct_exp,
      stringsAsFactors = FALSE
    )
  }
  
  out <- dplyr::bind_rows(
    rows
  )
  
  out <- out |>
    dplyr::group_by(
      Gene
    ) |>
    dplyr::mutate(
      ScaledAverageExpressionRaw =
        zscore_safe(
          AverageExpressionLog1p
        ),
      ScaledAverageExpression_DotPlot =
        pmax(
          pmin(
            ScaledAverageExpressionRaw,
            2.5
          ),
          -2.5
        )
    ) |>
    dplyr::ungroup()
  
  out
}

############################################################
## 4. Frozen lineage and preprocessing parameters
############################################################

## Script-07 primary lineage.
expected_parent_n_cells <- 26053L
expected_parent_n_features <- 11616L

expected_parent_cellset_md5 <-
  "04a351e47fcc11bdd28cae81976d1133"

## Canonical malignant-lineage subset from script 07:
## Malignant (17,708) + Cycling malignant (1,888) = 19,596.
expected_malignant_n_cells <- 19596L

## Cross-platform LF hash of the 19,596-cell malignant subset.
expected_malignant_cellset_md5 <-
  "4094e7ece0a55fc8579564b14237c03d"

expected_malignant_major_counts <- c(
  "Malignant" = 17708L,
  "Cycling malignant" = 1888L
)

expected_malignant_sample_counts <- c(
  "Pat_ICBnaive1" = 6372L,
  "Pat_ICBnaive2" = 7369L,
  "Pat42" = 5833L,
  "Pat5" = 22L
)

## Parent script-07 Seurat clusters represented in the malignant lineage.
expected_parent_cluster_counts <- c(
  "0" = 4547L,
  "1" = 3896L,
  "2" = 2705L,
  "3" = 2343L,
  "4" = 2205L,
  "5" = 2012L,
  "6" = 1888L
)

malignant_labels <- c(
  "Malignant",
  "Cycling malignant"
)

## Historical mother-script parameters.
global_seed <- 20260430L
n_variable_features <- 2000L
n_pcs_compute <- 40L
dims_use <- 1:30
cluster_resolution <- 0.6

## Explicit Seurat defaults used by the historical code path.
## Making them explicit improves version-to-version reproducibility without
## intentionally changing the historical clustering definition.
pca_seed <- 42L
findclusters_random_seed <- 0L
umap_seed <- 42L

expected_malignant_subclusters <-
  paste0(
    "Mal_",
    0:9
  )

parameters <- data.frame(
  Parameter = c(
    "malignant_labels",
    "global_seed",
    "n_variable_features",
    "n_pcs_compute",
    "dims_use",
    "cluster_resolution",
    "pca_seed",
    "findclusters_random_seed",
    "umap_seed",
    "cellset_md5_convention"
  ),
  Value = c(
    paste(
      malignant_labels,
      collapse = ";"
    ),
    as.character(
      global_seed
    ),
    as.character(
      n_variable_features
    ),
    as.character(
      n_pcs_compute
    ),
    paste(
      dims_use,
      collapse = ","
    ),
    as.character(
      cluster_resolution
    ),
    as.character(
      pca_seed
    ),
    as.character(
      findclusters_random_seed
    ),
    as.character(
      umap_seed
    ),
    "sorted unique UTF-8 cell IDs joined by LF with final LF; binary-byte MD5"
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  parameters,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subclustering_parameters.csv"
  )
)

############################################################
## 5. Compact Figure-5 marker/program definitions
############################################################

malignant_marker_sets <- list(
  "Tumor plasticity/dedifferentiation" = c(
    "AXL",
    "NGFR",
    "WNT5A",
    "EGFR",
    "SOX9",
    "ZEB1"
  ),
  "Melanocytic differentiation" = c(
    "MLANA",
    "PMEL",
    "TYR",
    "DCT",
    "MITF",
    "SOX10"
  ),
  "Antigen presentation" = c(
    "B2M",
    "HLA-A",
    "HLA-B",
    "TAP1",
    "PSMB9",
    "STAT1",
    "IRF1"
  ),
  "IFN response" = c(
    "IFIT1",
    "ISG15",
    "IFIT3",
    "CXCL10",
    "GBP1",
    "OAS1"
  )
)

program_score_map <- c(
  "Tumor plasticity/dedifferentiation" =
    "Tumor_plasticity_dedifferentiation_score",
  "Melanocytic differentiation" =
    "Melanocytic_differentiation_score",
  "Antigen presentation" =
    "Antigen_presentation_score",
  "IFN response" =
    "IFN_response_score"
)

marker_definition <- dplyr::bind_rows(
  lapply(
    names(
      malignant_marker_sets
    ),
    function(program_name) {
      data.frame(
        MarkerProgram = program_name,
        Gene = malignant_marker_sets[[program_name]],
        ScoreColumn =
          unname(
            program_score_map[[program_name]]
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

marker_definition$Order <-
  seq_len(
    nrow(
      marker_definition
    )
  )

safe_write_csv(
  marker_definition,
  file.path(
    out_table_dir,
    "GSE244983_malignant_compact_program_definitions.csv"
  )
)

############################################################
## 6. Input audit and parent-lineage gates
############################################################

input_audit <- data.frame(
  Input = "script07_major_annotated_Seurat",
  RelativePath = relative_to_project(
    input_seurat_file
  ),
  Exists = file.exists(
    input_seurat_file
  ),
  MD5 = md5_file(
    input_seurat_file
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_audit,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subclustering_input_audit.csv"
  )
)

message(
  "Loading canonical script-07 Seurat object..."
)

obj <- readRDS(
  input_seurat_file
)

if (!inherits(obj, "Seurat")) {
  stop(
    "Canonical script-07 input is not a Seurat object."
  )
}

DefaultAssay(obj) <- "RNA"

required_parent_meta <- c(
  "Sample",
  "SeuratCluster",
  "MajorCellType"
)

missing_parent_meta <- setdiff(
  required_parent_meta,
  colnames(
    obj@meta.data
  )
)

if (length(missing_parent_meta) > 0) {
  stop(
    "Canonical script-07 object lacks required metadata column(s): ",
    paste(
      missing_parent_meta,
      collapse = ", "
    )
  )
}

parent_cellset_md5 <- cellset_md5(
  colnames(
    obj
  )
)

parent_gates <- data.frame(
  Gate = c(
    "Script-07 object has 26,053 cells",
    "Script-07 object has 11,616 features",
    "Script-07 cell-set MD5 matches frozen primary lineage",
    "Sample has no missing labels",
    "SeuratCluster has no missing labels",
    "MajorCellType has no missing labels"
  ),
  Observed = c(
    as.character(
      ncol(
        obj
      )
    ),
    as.character(
      nrow(
        obj
      )
    ),
    parent_cellset_md5,
    as.character(
      sum(
        is.na(
          obj$Sample
        )
      )
    ),
    as.character(
      sum(
        is.na(
          obj$SeuratCluster
        )
      )
    ),
    as.character(
      sum(
        is.na(
          obj$MajorCellType
        )
      )
    )
  ),
  Expected = c(
    as.character(
      expected_parent_n_cells
    ),
    as.character(
      expected_parent_n_features
    ),
    expected_parent_cellset_md5,
    "0",
    "0",
    "0"
  ),
  Pass = c(
    ncol(
      obj
    ) ==
      expected_parent_n_cells,
    nrow(
      obj
    ) ==
      expected_parent_n_features,
    identical(
      parent_cellset_md5,
      expected_parent_cellset_md5
    ),
    sum(
      is.na(
        obj$Sample
      )
    ) ==
      0,
    sum(
      is.na(
        obj$SeuratCluster
      )
    ) ==
      0,
    sum(
      is.na(
        obj$MajorCellType
      )
    ) ==
      0
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  parent_gates,
  file.path(
    out_table_dir,
    "GSE244983_script07_parent_lineage_gates.csv"
  )
)

if (
  any(
    !parent_gates$Pass
  )
) {
  session_file <- file.path(
    log_dir,
    paste0(
      "sessionInfo_",
      SCRIPT_ID,
      ".txt"
    )
  )
  
  writeLines(
    capture.output(
      sessionInfo()
    ),
    session_file
  )
  
  failed_parent <- parent_gates[
    !parent_gates$Pass,
    ,
    drop = FALSE
  ]
  
  stop(
    "Script-07 parent-lineage gate failure:\n",
    paste0(
      failed_parent$Gate,
      " | observed=",
      failed_parent$Observed,
      " | expected=",
      failed_parent$Expected,
      collapse = "\n"
    ),
    "\nMalignant subclustering was not run."
  )
}

############################################################
## 7. Select canonical malignant-lineage cells
############################################################

major_values <- as.character(
  obj$MajorCellType
)

malignant_cells <- colnames(
  obj
)[
  major_values %in%
    malignant_labels
]

message(
  "Canonical malignant-lineage cells selected: ",
  length(
    malignant_cells
  )
)

malignant_subset_meta <- obj@meta.data[
  malignant_cells,
  ,
  drop = FALSE
]

malignant_subset_md5 <- cellset_md5(
  malignant_cells
)

observed_malignant_major_raw <- table(
  as.character(
    malignant_subset_meta$MajorCellType
  )
)

observed_malignant_major <- setNames(
  integer(
    length(
      expected_malignant_major_counts
    )
  ),
  names(
    expected_malignant_major_counts
  )
)

major_common <- intersect(
  names(
    observed_malignant_major_raw
  ),
  names(
    expected_malignant_major_counts
  )
)

observed_malignant_major[
  major_common
] <- as.integer(
  observed_malignant_major_raw[
    major_common
  ]
)

observed_malignant_sample_raw <- table(
  as.character(
    malignant_subset_meta$Sample
  )
)

observed_malignant_sample <- setNames(
  integer(
    length(
      expected_malignant_sample_counts
    )
  ),
  names(
    expected_malignant_sample_counts
  )
)

sample_common <- intersect(
  names(
    observed_malignant_sample_raw
  ),
  names(
    expected_malignant_sample_counts
  )
)

observed_malignant_sample[
  sample_common
] <- as.integer(
  observed_malignant_sample_raw[
    sample_common
  ]
)

observed_parent_cluster_raw <- table(
  as.character(
    malignant_subset_meta$SeuratCluster
  )
)

observed_parent_cluster <- setNames(
  integer(
    length(
      expected_parent_cluster_counts
    )
  ),
  names(
    expected_parent_cluster_counts
  )
)

parent_cluster_common <- intersect(
  names(
    observed_parent_cluster_raw
  ),
  names(
    expected_parent_cluster_counts
  )
)

observed_parent_cluster[
  parent_cluster_common
] <- as.integer(
  observed_parent_cluster_raw[
    parent_cluster_common
  ]
)

subset_gates <- data.frame(
  Gate = c(
    "Canonical malignant-lineage subset has 19,596 cells",
    "Canonical malignant-lineage cell-set MD5 matches frozen script-07 subset",
    "Malignant/Cycling malignant counts match frozen script-07 lineage",
    "Malignant-lineage sample counts match frozen script-07 lineage",
    "Parent Seurat-cluster counts match frozen script-07 malignant lineage",
    "Only malignant-lineage major labels are retained"
  ),
  Observed = c(
    as.character(
      length(
        malignant_cells
      )
    ),
    malignant_subset_md5,
    paste0(
      names(
        expected_malignant_major_counts
      ),
      "=",
      observed_malignant_major[
        names(
          expected_malignant_major_counts
        )
      ],
      collapse = ";"
    ),
    paste0(
      names(
        expected_malignant_sample_counts
      ),
      "=",
      observed_malignant_sample[
        names(
          expected_malignant_sample_counts
        )
      ],
      collapse = ";"
    ),
    paste0(
      names(
        expected_parent_cluster_counts
      ),
      "=",
      observed_parent_cluster[
        names(
          expected_parent_cluster_counts
        )
      ],
      collapse = ";"
    ),
    paste(
      sort(
        unique(
          as.character(
            malignant_subset_meta$MajorCellType
          )
        )
      ),
      collapse = ";"
    )
  ),
  Expected = c(
    as.character(
      expected_malignant_n_cells
    ),
    expected_malignant_cellset_md5,
    paste0(
      names(
        expected_malignant_major_counts
      ),
      "=",
      expected_malignant_major_counts,
      collapse = ";"
    ),
    paste0(
      names(
        expected_malignant_sample_counts
      ),
      "=",
      expected_malignant_sample_counts,
      collapse = ";"
    ),
    paste0(
      names(
        expected_parent_cluster_counts
      ),
      "=",
      expected_parent_cluster_counts,
      collapse = ";"
    ),
    paste(
      sort(
        malignant_labels
      ),
      collapse = ";"
    )
  ),
  Pass = c(
    length(
      malignant_cells
    ) ==
      expected_malignant_n_cells,
    identical(
      malignant_subset_md5,
      expected_malignant_cellset_md5
    ),
    identical(
      as.integer(
        observed_malignant_major[
          names(
            expected_malignant_major_counts
          )
        ]
      ),
      as.integer(
        expected_malignant_major_counts
      )
    ),
    identical(
      as.integer(
        observed_malignant_sample[
          names(
            expected_malignant_sample_counts
          )
        ]
      ),
      as.integer(
        expected_malignant_sample_counts
      )
    ),
    identical(
      as.integer(
        observed_parent_cluster[
          names(
            expected_parent_cluster_counts
          )
        ]
      ),
      as.integer(
        expected_parent_cluster_counts
      )
    ),
    setequal(
      unique(
        as.character(
          malignant_subset_meta$MajorCellType
        )
      ),
      malignant_labels
    )
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  subset_gates,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subset_lineage_gates.csv"
  )
)

if (
  any(
    !subset_gates$Pass
  )
) {
  session_file <- file.path(
    log_dir,
    paste0(
      "sessionInfo_",
      SCRIPT_ID,
      ".txt"
    )
  )
  
  writeLines(
    capture.output(
      sessionInfo()
    ),
    session_file
  )
  
  failed_subset <- subset_gates[
    !subset_gates$Pass,
    ,
    drop = FALSE
  ]
  
  stop(
    "Canonical malignant-lineage gate failure:\n",
    paste0(
      failed_subset$Gate,
      " | observed=",
      failed_subset$Observed,
      " | expected=",
      failed_subset$Expected,
      collapse = "\n"
    ),
    "\nMalignant re-clustering was not run."
  )
}

mal_obj <- subset(
  obj,
  cells = malignant_cells
)

DefaultAssay(
  mal_obj
) <- "RNA"

## Preserve parent script-07 labels explicitly before malignant-only clustering.
mal_obj@meta.data$ParentSeuratCluster <-
  as.character(
    mal_obj@meta.data$SeuratCluster
  )

mal_obj@meta.data$ParentMajorCellType <-
  as.character(
    mal_obj@meta.data$MajorCellType
  )

rm(
  obj,
  malignant_subset_meta
)

gc(
  verbose = FALSE
)

############################################################
## 8. Confirm normalized RNA layer
############################################################

rna_data_before <- get_normalized_rna(
  mal_obj
)

if (
  ncol(
    rna_data_before
  ) !=
  expected_malignant_n_cells
) {
  stop(
    "Normalized RNA layer does not contain all canonical malignant-lineage cells."
  )
}

############################################################
## 9. Reprocess malignant-lineage cells
############################################################

set.seed(
  global_seed
)

message(
  "Finding 2,000 malignant-lineage variable features..."
)

mal_obj <- Seurat::FindVariableFeatures(
  mal_obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = n_variable_features,
  verbose = FALSE
)

if (
  length(
    Seurat::VariableFeatures(
      mal_obj,
      assay = "RNA"
    )
  ) <
  n_variable_features
) {
  stop(
    "Fewer than 2,000 variable features were returned."
  )
}

mal_obj <- Seurat::ScaleData(
  mal_obj,
  assay = "RNA",
  features = Seurat::VariableFeatures(
    mal_obj,
    assay = "RNA"
  ),
  verbose = FALSE
)

mal_obj <- Seurat::RunPCA(
  mal_obj,
  assay = "RNA",
  features = Seurat::VariableFeatures(
    mal_obj,
    assay = "RNA"
  ),
  npcs = n_pcs_compute,
  seed.use = pca_seed,
  verbose = FALSE
)

available_pcs <- ncol(
  Seurat::Embeddings(
    mal_obj,
    "pca"
  )
)

if (
  available_pcs <
  max(
    dims_use
  )
) {
  stop(
    "PCA returned only ",
    available_pcs,
    " components; expected at least ",
    max(
      dims_use
    ),
    "."
  )
}

message(
  "Running malignant-lineage Harmony using Sample..."
)

mal_obj <- harmony::RunHarmony(
  object = mal_obj,
  group.by.vars = "Sample",
  reduction.use = "pca",
  dims.use = dims_use,
  reduction.save = "harmony",
  verbose = FALSE
)

if (
  !"harmony" %in%
  Seurat::Reductions(
    mal_obj
  )
) {
  stop(
    "Harmony reduction was not created. ",
    "Script 10 does not silently fall back to PCA."
  )
}

available_harmony_dims <- ncol(
  Seurat::Embeddings(
    mal_obj,
    "harmony"
  )
)

if (
  available_harmony_dims <
  max(
    dims_use
  )
) {
  stop(
    "Harmony returned only ",
    available_harmony_dims,
    " dimensions; expected at least ",
    max(
      dims_use
    ),
    "."
  )
}

mal_obj <- Seurat::FindNeighbors(
  mal_obj,
  reduction = "harmony",
  dims = dims_use,
  verbose = FALSE
)

mal_obj <- Seurat::FindClusters(
  mal_obj,
  resolution = cluster_resolution,
  random.seed = findclusters_random_seed,
  verbose = FALSE
)

mal_obj <- Seurat::RunUMAP(
  mal_obj,
  reduction = "harmony",
  dims = dims_use,
  seed.use = umap_seed,
  verbose = FALSE
)

malignant_cluster_id <-
  as.character(
    Seurat::Idents(
      mal_obj
    )
  )

mal_obj@meta.data$MalignantClusterID <-
  malignant_cluster_id

cluster_numeric <- suppressWarnings(
  as.numeric(
    unique(
      malignant_cluster_id
    )
  )
)

if (
  any(
    is.na(
      cluster_numeric
    )
  )
) {
  stop(
    "Malignant cluster identities are not numeric and cannot be mapped deterministically to Mal_# labels."
  )
}

cluster_levels_numeric <- sort(
  unique(
    as.numeric(
      malignant_cluster_id
    )
  )
)

cluster_levels_chr <- as.character(
  cluster_levels_numeric
)

malignant_subcluster_map <- setNames(
  paste0(
    "Mal_",
    cluster_levels_chr
  ),
  cluster_levels_chr
)

mal_obj@meta.data$MalignantSubcluster <- factor(
  malignant_subcluster_map[
    malignant_cluster_id
  ],
  levels = paste0(
    "Mal_",
    cluster_levels_chr
  )
)

message(
  "Malignant subclusters detected: ",
  paste(
    levels(
      mal_obj$MalignantSubcluster
    ),
    collapse = ", "
  )
)

############################################################
## 10. Compact malignant program scores for Figure 5
############################################################

expr_data <- get_normalized_rna(
  mal_obj
)

program_presence_rows <- list()

for (
  program_name in names(
    malignant_marker_sets
  )
) {
  score_col <-
    unname(
      program_score_map[[program_name]]
    )
  
  scored <- score_module_simple(
    expr_data = expr_data,
    requested_genes =
      malignant_marker_sets[[program_name]],
    cell_order = rownames(
      mal_obj@meta.data
    ),
    score_name = score_col
  )
  
  mal_obj@meta.data[[score_col]] <- scored$score
  
  pp <- scored$presence
  pp$MarkerProgram <- program_name
  pp$ScoreColumn <- score_col
  
  program_presence_rows[[program_name]] <- pp
  
  message(
    "Program ",
    program_name,
    ": ",
    sum(
      pp$Present
    ),
    "/",
    nrow(
      pp
    ),
    " genes present"
  )
}

program_gene_presence <- dplyr::bind_rows(
  program_presence_rows
)

safe_write_csv(
  program_gene_presence,
  file.path(
    out_table_dir,
    "GSE244983_malignant_compact_program_gene_presence.csv"
  )
)

program_score_cols <- unname(
  program_score_map[
    names(
      malignant_marker_sets
    )
  ]
)

############################################################
## 11. Cell-level and subcluster composition tables
############################################################

cell_table <- data.frame(
  Cell = rownames(
    mal_obj@meta.data
  ),
  Sample = as.character(
    mal_obj$Sample
  ),
  ParentSeuratCluster = as.character(
    mal_obj$ParentSeuratCluster
  ),
  ParentMajorCellType = as.character(
    mal_obj$ParentMajorCellType
  ),
  MalignantClusterID = as.character(
    mal_obj$MalignantClusterID
  ),
  MalignantSubcluster = as.character(
    mal_obj$MalignantSubcluster
  ),
  stringsAsFactors = FALSE
)

for (
  cc in program_score_cols
) {
  cell_table[[cc]] <- as.numeric(
    mal_obj@meta.data[[cc]]
  )
}

safe_write_csv(
  cell_table,
  file.path(
    out_table_dir,
    "GSE244983_malignant_cells_with_subclusters_and_compact_program_scores.csv"
  )
)

subcluster_counts <- cell_table |>
  dplyr::count(
    MalignantSubcluster,
    name = "N_cells"
  ) |>
  dplyr::mutate(
    Proportion =
      N_cells /
      sum(
        N_cells
      )
  ) |>
  dplyr::arrange(
    suppressWarnings(
      as.numeric(
        sub(
          "^Mal_",
          "",
          MalignantSubcluster
        )
      )
    )
  )

safe_write_csv(
  subcluster_counts,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subcluster_counts.csv"
  )
)

sample_composition <- cell_table |>
  dplyr::count(
    Sample,
    MalignantSubcluster,
    name = "N_cells"
  ) |>
  dplyr::group_by(
    Sample
  ) |>
  dplyr::mutate(
    Proportion =
      N_cells /
      sum(
        N_cells
      )
  ) |>
  dplyr::ungroup()

safe_write_csv(
  sample_composition,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subcluster_composition_by_sample.csv"
  )
)

parent_origin <- cell_table |>
  dplyr::count(
    ParentSeuratCluster,
    ParentMajorCellType,
    MalignantSubcluster,
    name = "N_cells"
  ) |>
  dplyr::group_by(
    ParentSeuratCluster,
    ParentMajorCellType
  ) |>
  dplyr::mutate(
    ProportionWithinParent =
      N_cells /
      sum(
        N_cells
      )
  ) |>
  dplyr::ungroup()

safe_write_csv(
  parent_origin,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subcluster_parent_origin.csv"
  )
)

############################################################
## 12. Compact program summaries and Figure-5 heatmap source
############################################################

program_summary <- cell_table |>
  dplyr::group_by(
    MalignantSubcluster
  ) |>
  dplyr::summarise(
    N_cells = dplyr::n(),
    dplyr::across(
      dplyr::all_of(
        program_score_cols
      ),
      list(
        mean = ~ mean(
          .x,
          na.rm = TRUE
        ),
        median = ~ stats::median(
          .x,
          na.rm = TRUE
        )
      )
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    suppressWarnings(
      as.numeric(
        sub(
          "^Mal_",
          "",
          MalignantSubcluster
        )
      )
    )
  )

safe_write_csv(
  program_summary,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subcluster_compact_program_summary.csv"
  )
)

mean_cols <- paste0(
  program_score_cols,
  "_mean"
)

program_matrix <- as.matrix(
  program_summary[
    ,
    mean_cols,
    drop = FALSE
  ]
)

rownames(
  program_matrix
) <- program_summary$MalignantSubcluster

program_matrix_z <- scale(
  program_matrix
)

program_matrix_z[
  !is.finite(
    program_matrix_z
  )
] <- 0

program_heatmap <- as.data.frame(
  program_matrix_z
)

program_heatmap$MalignantSubcluster <-
  rownames(
    program_heatmap
  )

program_heatmap <- program_heatmap |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      mean_cols
    ),
    names_to = "ProgramScoreColumn",
    values_to = "Average_z_score"
  ) |>
  dplyr::mutate(
    ProgramScoreColumn =
      sub(
        "_mean$",
        "",
        ProgramScoreColumn
      ),
    MarkerProgram =
      names(
        program_score_map
      )[
        match(
          ProgramScoreColumn,
          unname(
            program_score_map
          )
        )
      ]
  ) |>
  dplyr::select(
    MalignantSubcluster,
    MarkerProgram,
    ProgramScoreColumn,
    Average_z_score
  )

safe_write_csv(
  program_heatmap,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subcluster_compact_program_heatmap_matrix.csv"
  )
)

top_program_by_subcluster <- program_heatmap |>
  dplyr::group_by(
    MalignantSubcluster
  ) |>
  dplyr::slice_max(
    order_by = Average_z_score,
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

safe_write_csv(
  top_program_by_subcluster,
  file.path(
    out_table_dir,
    "GSE244983_top_compact_program_by_malignant_subcluster.csv"
  )
)

############################################################
## 13. Malignant-only UMAP source for Figure 5A/D
############################################################

umap_mat <- as.data.frame(
  Seurat::Embeddings(
    mal_obj,
    "umap"
  )
)

if (
  ncol(
    umap_mat
  ) <
  2
) {
  stop(
    "Malignant-only UMAP has fewer than two dimensions."
  )
}

umap_source <- data.frame(
  Cell = rownames(
    umap_mat
  ),
  UMAP_1 = umap_mat[[1]],
  UMAP_2 = umap_mat[[2]],
  stringsAsFactors = FALSE
)

umap_source <- merge(
  umap_source,
  cell_table,
  by = "Cell",
  all.x = TRUE,
  sort = FALSE
)

umap_source <- umap_source[
  match(
    rownames(
      umap_mat
    ),
    umap_source$Cell
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  umap_source,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subcluster_UMAP_source.csv"
  )
)

############################################################
## 14. DotPlot-compatible source for Figure 5B
############################################################

marker_order <- marker_definition$Gene

marker_match <- match_features_case_insensitive(
  rownames(
    expr_data
  ),
  marker_order
)

marker_presence <- merge(
  marker_definition,
  marker_match,
  by.x = "Gene",
  by.y = "RequestedGene",
  all.x = TRUE,
  sort = FALSE
)

marker_presence <- marker_presence[
  match(
    marker_order,
    marker_presence$Gene
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  marker_presence,
  file.path(
    out_table_dir,
    "GSE244983_malignant_marker_gene_presence.csv"
  )
)

present_marker_rows <- marker_presence[
  marker_presence$Present,
  ,
  drop = FALSE
]

marker_expr_dense <- as.matrix(
  expr_data[
    present_marker_rows$MatchedFeature,
    ,
    drop = FALSE
  ]
)

dotplot_source <- build_dotplot_source(
  marker_expr_dense = marker_expr_dense,
  marker_features =
    present_marker_rows$MatchedFeature,
  marker_labels =
    present_marker_rows$Gene,
  marker_programs =
    present_marker_rows$MarkerProgram,
  groups = as.character(
    mal_obj$MalignantSubcluster
  ),
  group_levels =
    levels(
      mal_obj$MalignantSubcluster
    )
)

dotplot_source$MalignantSubcluster <- factor(
  dotplot_source$MalignantSubcluster,
  levels =
    levels(
      mal_obj$MalignantSubcluster
    )
)

dotplot_source$Gene <- factor(
  dotplot_source$Gene,
  levels =
    present_marker_rows$Gene
)

dotplot_source <- dotplot_source |>
  dplyr::arrange(
    MalignantSubcluster,
    Gene
  )

dotplot_source$MalignantSubcluster <-
  as.character(
    dotplot_source$MalignantSubcluster
  )

dotplot_source$Gene <-
  as.character(
    dotplot_source$Gene
  )

safe_write_csv(
  dotplot_source,
  file.path(
    out_table_dir,
    "GSE244983_malignant_marker_dotplot_source.csv"
  )
)


############################################################
## 14B. Frozen five-program malignant analysis and
##      Supplementary Figure S10 reporting sources
############################################################

## IMPORTANT:
## The original compact Figure-5 program analysis above is preserved exactly.
## This additional block implements the frozen five-program split used by
## Main Figure 6 / Supplementary Figure S10 and exports reporting-only sources.
##
## No downstream figure script should rediscover metadata columns or recompute
## these summaries from an RDS.

five_program_gene_sets <- list(
  Tumor_plasticity_dedifferentiation = c(
    "AXL", "NGFR", "WNT5A", "EGFR", "SOX9", "ZEB1", "FN1", "VIM"
  ),
  Stromal_ECM_remodeling = c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2",
    "COL6A3", "DCN", "LUM", "FAP", "ACTA2", "PDPN", "THY1", "VCAN",
    "TNC", "TGFBI", "TIMP1", "INHBA", "MRC2", "ADAM12"
  ),
  Melanocytic_differentiation = c(
    "MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"
  ),
  Antigen_presentation = c(
    "B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2",
    "PSMB8", "PSMB9", "STAT1", "IRF1"
  ),
  IFN_response = c(
    "IFIT1", "IFIT2", "IFIT3", "ISG15", "CXCL10", "GBP1",
    "OAS1", "IFI6", "MX1"
  )
)

five_program_order <- names(five_program_gene_sets)

five_program_display <- c(
  Tumor_plasticity_dedifferentiation =
    "Tumor plasticity/dedifferentiation",
  Stromal_ECM_remodeling =
    "Stromal/ECM remodeling",
  Melanocytic_differentiation =
    "Melanocytic differentiation",
  Antigen_presentation =
    "Antigen presentation",
  IFN_response =
    "IFN response"
)

five_raw_cols <- paste0(five_program_order, "_score")
five_z_cols <- paste0(five_raw_cols, "_z")

five_program_definition <- dplyr::bind_rows(
  lapply(
    five_program_order,
    function(pg) {
      data.frame(
        Program = pg,
        ProgramDisplay = unname(five_program_display[[pg]]),
        Gene = five_program_gene_sets[[pg]],
        RawScoreColumn = paste0(pg, "_score"),
        ZScoreColumn = paste0(pg, "_score_z"),
        stringsAsFactors = FALSE
      )
    }
  )
)

safe_write_csv(
  five_program_definition,
  file.path(
    out_table_dir,
    "GSE244983_malignant_five_program_definitions.csv"
  )
)

five_presence_rows <- vector("list", length(five_program_order))

for (i in seq_along(five_program_order)) {
  pg <- five_program_order[[i]]
  
  scored <- score_module_simple(
    expr_data = expr_data,
    requested_genes = five_program_gene_sets[[pg]],
    cell_order = rownames(mal_obj@meta.data),
    score_name = pg
  )
  
  raw_col <- five_raw_cols[[i]]
  z_col <- five_z_cols[[i]]
  
  mal_obj@meta.data[[raw_col]] <- scored$score
  mal_obj@meta.data[[z_col]] <- zscore_safe(scored$score)
  
  pp <- scored$presence
  pp$Program <- pg
  pp$ProgramDisplay <- unname(five_program_display[[pg]])
  pp$RawScoreColumn <- raw_col
  pp$ZScoreColumn <- z_col
  five_presence_rows[[i]] <- pp
}

five_program_presence <- dplyr::bind_rows(five_presence_rows)

safe_write_csv(
  five_program_presence,
  file.path(
    out_table_dir,
    "GSE244983_malignant_five_program_gene_presence.csv"
  )
)

five_cell_table <- data.frame(
  Cell = rownames(mal_obj@meta.data),
  MalignantSubcluster = as.character(mal_obj$MalignantSubcluster),
  stringsAsFactors = FALSE
)

for (cc in c(five_raw_cols, five_z_cols)) {
  five_cell_table[[cc]] <- as.numeric(mal_obj@meta.data[[cc]])
}

safe_write_csv(
  five_cell_table,
  file.path(
    out_table_dir,
    "GSE244983_malignant_cells_with_five_program_scores.csv"
  )
)

## S10B: exact cell-level z-score reporting source.
s10b_rows <- vector("list", length(five_program_order))

for (i in seq_along(five_program_order)) {
  pg <- five_program_order[[i]]
  
  s10b_rows[[i]] <- data.frame(
    Cell = rownames(mal_obj@meta.data),
    MalignantSubcluster = as.character(mal_obj$MalignantSubcluster),
    Program = pg,
    ProgramDisplay = unname(five_program_display[[pg]]),
    ProgramZ = as.numeric(mal_obj@meta.data[[five_z_cols[[i]]]]),
    stringsAsFactors = FALSE
  )
}

s10b <- dplyr::bind_rows(s10b_rows)

safe_write_csv(
  s10b,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS10B_program_z_distribution_source.csv"
  )
)

## S10C: mean cell-level program z-score by frozen malignant subcluster.
s10c <- s10b |>
  dplyr::group_by(
    MalignantSubcluster,
    Program,
    ProgramDisplay
  ) |>
  dplyr::summarise(
    MeanProgramZ = mean(ProgramZ, na.rm = TRUE),
    N_cells = dplyr::n(),
    .groups = "drop"
  )

safe_write_csv(
  s10c,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS10C_mean_program_z_source.csv"
  )
)

## S10D: cell-level Spearman matrix across the same five program-z columns.
s10_wide <- s10b |>
  dplyr::select(
    Cell,
    Program,
    ProgramZ
  ) |>
  tidyr::pivot_wider(
    names_from = Program,
    values_from = ProgramZ
  )

s10_cor <- suppressWarnings(
  stats::cor(
    as.matrix(
      s10_wide[
        ,
        five_program_order,
        drop = FALSE
      ]
    ),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
)

s10d_rows <- vector(
  "list",
  length(five_program_order) * length(five_program_order)
)

kk_s10 <- 0L
for (row_pg in five_program_order) {
  for (col_pg in five_program_order) {
    kk_s10 <- kk_s10 + 1L
    
    s10d_rows[[kk_s10]] <- data.frame(
      ProgramRow = row_pg,
      ProgramRowDisplay = unname(five_program_display[[row_pg]]),
      ProgramColumn = col_pg,
      ProgramColumnDisplay = unname(five_program_display[[col_pg]]),
      SpearmanRho = s10_cor[row_pg, col_pg],
      stringsAsFactors = FALSE
    )
  }
}

s10d <- dplyr::bind_rows(s10d_rows)

safe_write_csv(
  s10d,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS10D_program_correlation_source.csv"
  )
)

## S10A: frozen five-program representative-marker panel.
s10_marker_sets <- list(
  Tumor_plasticity_dedifferentiation = c(
    "AXL", "NGFR", "WNT5A", "EGFR", "SOX9", "ZEB1", "VIM", "FN1"
  ),
  Stromal_ECM_remodeling = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "THY1", "VCAN"
  ),
  Melanocytic_differentiation = c(
    "MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10"
  ),
  Antigen_presentation = c(
    "B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2"
  ),
  IFN_response = c(
    "IFIT1", "IFIT2", "IFIT3", "ISG15", "CXCL10", "GBP1", "OAS1"
  )
)

s10_marker_definition <- dplyr::bind_rows(
  lapply(
    names(s10_marker_sets),
    function(pg) {
      data.frame(
        MarkerGroup = pg,
        Gene = s10_marker_sets[[pg]],
        stringsAsFactors = FALSE
      )
    }
  )
)

s10_marker_match <- match_features_case_insensitive(
  rownames(expr_data),
  s10_marker_definition$Gene
)

s10_marker_presence <- merge(
  s10_marker_definition,
  s10_marker_match,
  by.x = "Gene",
  by.y = "RequestedGene",
  all.x = TRUE,
  sort = FALSE
)

s10_marker_presence <- s10_marker_presence[
  match(
    s10_marker_definition$Gene,
    s10_marker_presence$Gene
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  s10_marker_definition,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS10_marker_definitions.csv"
  )
)

safe_write_csv(
  s10_marker_presence,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS10_marker_gene_presence.csv"
  )
)

s10_present <- s10_marker_presence[
  s10_marker_presence$Present,
  ,
  drop = FALSE
]

s10_marker_expr_dense <- as.matrix(
  expr_data[
    s10_present$MatchedFeature,
    ,
    drop = FALSE
  ]
)

s10a0 <- build_dotplot_source(
  marker_expr_dense = s10_marker_expr_dense,
  marker_features = s10_present$MatchedFeature,
  marker_labels = s10_present$Gene,
  marker_programs = s10_present$MarkerGroup,
  groups = as.character(mal_obj$MalignantSubcluster),
  group_levels = levels(mal_obj$MalignantSubcluster)
)

s10a <- s10a0 |>
  dplyr::transmute(
    MarkerGroup = MarkerProgram,
    Gene = Gene,
    MalignantSubcluster = MalignantSubcluster,
    PercentExpressed = PercentExpressed,
    AverageExpressionScaledDisplay = ScaledAverageExpression_DotPlot
  )

safe_write_csv(
  s10a,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS10A_marker_dotplot_source.csv"
  )
)

## Explicit S10 figure-source manifest.
s10_reporting_manifest <- data.frame(
  Target = c(
    "Supplementary Figure S10A",
    "Supplementary Figure S10B",
    "Supplementary Figure S10C",
    "Supplementary Figure S10D"
  ),
  SourceCSV = c(
    "GSE244983_SuppFigureS10A_marker_dotplot_source.csv",
    "GSE244983_SuppFigureS10B_program_z_distribution_source.csv",
    "GSE244983_SuppFigureS10C_mean_program_z_source.csv",
    "GSE244983_SuppFigureS10D_program_correlation_source.csv"
  ),
  Use = c(
    "Five-program malignant marker dot plot",
    "Five-program cell-level z-score distributions",
    "Mean cell-level program z-score by malignant subcluster",
    "Cell-level five-program Spearman correlation matrix"
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  s10_reporting_manifest,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS10_figure_source_manifest.csv"
  )
)

five_presence_summary <- five_program_presence |>
  dplyr::group_by(Program) |>
  dplyr::summarise(
    N_genes = dplyr::n(),
    N_present = sum(Present),
    .groups = "drop"
  )

s10_reporting_gates <- data.frame(
  Gate = c(
    "All five frozen malignant programs contain at least two detected genes",
    "All five cell-level z-score columns are finite for 19,596 cells",
    "S10A covers Mal_0 through Mal_9",
    "S10B contains 97,980 cell-program rows",
    "S10B contains exactly 19,596 unique cells",
    "S10C contains exactly 50 subcluster-program rows",
    "S10D contains exactly 25 program-pair rows",
    "S10D is symmetric",
    "S10D diagonal is exactly one within floating-point tolerance"
  ),
  Observed = c(
    paste0(
      five_presence_summary$Program,
      "=",
      five_presence_summary$N_present,
      "/",
      five_presence_summary$N_genes,
      collapse = ";"
    ),
    paste0(
      five_z_cols,
      "=",
      vapply(
        five_z_cols,
        function(cc) sum(is.finite(mal_obj@meta.data[[cc]])),
        integer(1)
      ),
      collapse = ";"
    ),
    paste(sort(unique(s10a$MalignantSubcluster)), collapse = ";"),
    as.character(nrow(s10b)),
    as.character(length(unique(s10b$Cell))),
    as.character(nrow(s10c)),
    as.character(nrow(s10d)),
    as.character(max(abs(s10_cor - t(s10_cor)), na.rm = TRUE)),
    as.character(max(abs(diag(s10_cor) - 1), na.rm = TRUE))
  ),
  Expected = c(
    ">=2 detected genes in each of 5 programs",
    "19596 finite scores for each of 5 programs",
    paste(sort(expected_malignant_subclusters), collapse = ";"),
    as.character(expected_malignant_n_cells * length(five_program_order)),
    as.character(expected_malignant_n_cells),
    "50",
    "25",
    "<=1e-12",
    "<=1e-12"
  ),
  Pass = c(
    nrow(five_presence_summary) == 5L &&
      all(five_presence_summary$N_present >= 2L),
    all(
      vapply(
        five_z_cols,
        function(cc) sum(is.finite(mal_obj@meta.data[[cc]])),
        integer(1)
      ) == expected_malignant_n_cells
    ),
    setequal(unique(s10a$MalignantSubcluster), expected_malignant_subclusters),
    nrow(s10b) == expected_malignant_n_cells * length(five_program_order),
    length(unique(s10b$Cell)) == expected_malignant_n_cells,
    nrow(s10c) == 50L,
    nrow(s10d) == 25L,
    max(abs(s10_cor - t(s10_cor)), na.rm = TRUE) <= 1e-12,
    max(abs(diag(s10_cor) - 1), na.rm = TRUE) <= 1e-12
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  s10_reporting_gates,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS10_reporting_contract_gates.csv"
  )
)

if (any(!s10_reporting_gates$Pass)) {
  stop(
    "Supplementary Figure S10 reporting-contract gate failure: ",
    paste(
      s10_reporting_gates$Gate[!s10_reporting_gates$Pass],
      collapse = " | "
    )
  )
}

rm(
  s10_marker_expr_dense
)


rm(
  marker_expr_dense,
  expr_data,
  rna_data_before
)

gc(
  verbose = FALSE
)

############################################################
## 15. Figure/source-data manifest
############################################################

source_manifest <- data.frame(
  Artifact = c(
    "Main Figure 5A",
    "Main Figure 5B",
    "Main Figure 5C",
    "Main Figure 5D",
    "Supplementary Table S9 compact-program component"
  ),
  SourceTable = c(
    "GSE244983_malignant_subcluster_UMAP_source.csv",
    "GSE244983_malignant_marker_dotplot_source.csv",
    "GSE244983_malignant_subcluster_compact_program_heatmap_matrix.csv",
    "GSE244983_malignant_subcluster_UMAP_source.csv",
    "GSE244983_malignant_subcluster_compact_program_summary.csv"
  ),
  PlotRole = c(
    "Mal_0 to Mal_9 malignant-only UMAP",
    "Representative compact malignant marker programs across Mal_0 to Mal_9",
    "Subcluster-level standardized compact program heatmap",
    "Malignant-only UMAP with compact program feature scores",
    "Subcluster counts and compact program means/medians"
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  source_manifest,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subclustering_figure_source_manifest.csv"
  )
)

############################################################
## 16. Reproducibility gates
############################################################

observed_subclusters <- levels(
  mal_obj$MalignantSubcluster
)

observed_cluster_ids <- sort(
  unique(
    as.character(
      mal_obj$MalignantClusterID
    )
  )
)

program_presence_summary <- program_gene_presence |>
  dplyr::group_by(
    MarkerProgram
  ) |>
  dplyr::summarise(
    N_genes = dplyr::n(),
    N_present = sum(
      Present
    ),
    .groups = "drop"
  )

program_scores_finite <- vapply(
  program_score_cols,
  function(cc) {
    sum(
      is.finite(
        mal_obj@meta.data[[cc]]
      )
    )
  },
  integer(1)
)

malignant_object_cellset_md5 <- cellset_md5(
  colnames(
    mal_obj
  )
)

umap_cellset_md5 <- cellset_md5(
  umap_source$Cell
)

subcluster_total <- sum(
  subcluster_counts$N_cells
)

sample_composition_total <- sum(
  sample_composition$N_cells
)

expected_dotplot_rows <-
  length(
    observed_subclusters
  ) *
  nrow(
    present_marker_rows
  )

analysis_gates <- data.frame(
  Gate = c(
    "Malignant object retains exactly 19,596 cells",
    "Malignant object retains exact canonical malignant cell set",
    "Exactly 2,000 variable features are present",
    "PCA contains at least 30 dimensions",
    "Harmony contains at least 30 dimensions",
    "Exactly 10 malignant subclusters are detected",
    "Malignant numeric cluster IDs are exactly 0 through 9",
    "Malignant subcluster labels are exactly Mal_0 through Mal_9",
    "No malignant subcluster labels are missing",
    "Subcluster counts sum to 19,596",
    "Sample-by-subcluster counts sum to 19,596",
    "All four compact malignant programs contain at least two detected genes",
    "All four compact malignant program scores are finite for all 19,596 cells",
    "Malignant-only UMAP source retains exact canonical malignant cell set",
    "Figure-5 marker DotPlot source has expected rows",
    "Figure-5 marker DotPlot source has no missing scaled averages",
    "Figure-5 compact program heatmap has exactly 40 rows"
  ),
  Observed = c(
    as.character(
      ncol(
        mal_obj
      )
    ),
    malignant_object_cellset_md5,
    as.character(
      length(
        Seurat::VariableFeatures(
          mal_obj,
          assay = "RNA"
        )
      )
    ),
    as.character(
      ncol(
        Seurat::Embeddings(
          mal_obj,
          "pca"
        )
      )
    ),
    as.character(
      ncol(
        Seurat::Embeddings(
          mal_obj,
          "harmony"
        )
      )
    ),
    as.character(
      length(
        observed_subclusters
      )
    ),
    paste(
      observed_cluster_ids,
      collapse = ";"
    ),
    paste(
      observed_subclusters,
      collapse = ";"
    ),
    as.character(
      sum(
        is.na(
          mal_obj$MalignantSubcluster
        )
      )
    ),
    as.character(
      subcluster_total
    ),
    as.character(
      sample_composition_total
    ),
    paste0(
      program_presence_summary$MarkerProgram,
      "=",
      program_presence_summary$N_present,
      "/",
      program_presence_summary$N_genes,
      collapse = ";"
    ),
    paste0(
      program_score_cols,
      "=",
      program_scores_finite,
      collapse = ";"
    ),
    umap_cellset_md5,
    as.character(
      nrow(
        dotplot_source
      )
    ),
    as.character(
      sum(
        is.na(
          dotplot_source$ScaledAverageExpression_DotPlot
        )
      )
    ),
    as.character(
      nrow(
        program_heatmap
      )
    )
  ),
  Expected = c(
    as.character(
      expected_malignant_n_cells
    ),
    expected_malignant_cellset_md5,
    as.character(
      n_variable_features
    ),
    ">=30",
    ">=30",
    "10",
    paste(
      as.character(
        0:9
      ),
      collapse = ";"
    ),
    paste(
      expected_malignant_subclusters,
      collapse = ";"
    ),
    "0",
    as.character(
      expected_malignant_n_cells
    ),
    as.character(
      expected_malignant_n_cells
    ),
    ">=2 detected genes in each of 4 programs",
    "19596 finite scores for each of 4 programs",
    expected_malignant_cellset_md5,
    as.character(
      expected_dotplot_rows
    ),
    "0",
    "40"
  ),
  Pass = c(
    ncol(
      mal_obj
    ) ==
      expected_malignant_n_cells,
    identical(
      malignant_object_cellset_md5,
      expected_malignant_cellset_md5
    ),
    length(
      Seurat::VariableFeatures(
        mal_obj,
        assay = "RNA"
      )
    ) ==
      n_variable_features,
    ncol(
      Seurat::Embeddings(
        mal_obj,
        "pca"
      )
    ) >=
      30,
    ncol(
      Seurat::Embeddings(
        mal_obj,
        "harmony"
      )
    ) >=
      30,
    length(
      observed_subclusters
    ) ==
      10,
    identical(
      observed_cluster_ids,
      as.character(
        0:9
      )
    ),
    identical(
      observed_subclusters,
      expected_malignant_subclusters
    ),
    sum(
      is.na(
        mal_obj$MalignantSubcluster
      )
    ) ==
      0,
    subcluster_total ==
      expected_malignant_n_cells,
    sample_composition_total ==
      expected_malignant_n_cells,
    nrow(
      program_presence_summary
    ) ==
      4 &&
      all(
        program_presence_summary$N_present >=
          2
      ),
    all(
      program_scores_finite ==
        expected_malignant_n_cells
    ),
    identical(
      umap_cellset_md5,
      expected_malignant_cellset_md5
    ),
    nrow(
      dotplot_source
    ) ==
      expected_dotplot_rows,
    sum(
      is.na(
        dotplot_source$ScaledAverageExpression_DotPlot
      )
    ) ==
      0,
    nrow(
      program_heatmap
    ) ==
      40
  ),
  stringsAsFactors = FALSE
)

all_gates <- dplyr::bind_rows(
  data.frame(
    GateClass = "script07_parent_lineage",
    parent_gates,
    stringsAsFactors = FALSE
  ),
  data.frame(
    GateClass = "canonical_malignant_subset",
    subset_gates,
    stringsAsFactors = FALSE
  ),
  data.frame(
    GateClass = "malignant_subclustering",
    analysis_gates,
    stringsAsFactors = FALSE
  )
)

safe_write_csv(
  all_gates,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subclustering_reproducibility_gates.csv"
  )
)

############################################################
## 17. Session info
############################################################

session_file <- file.path(
  log_dir,
  paste0(
    "sessionInfo_",
    SCRIPT_ID,
    ".txt"
  )
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_file
)

message(
  "Saved: ",
  session_file
)

if (
  any(
    !all_gates$Pass
  )
) {
  stop(
    "10 reproducibility gate failure: ",
    paste(
      all_gates$Gate[
        !all_gates$Pass
      ],
      collapse = " | "
    ),
    ". Canonical malignant-subclustered Seurat object was NOT saved."
  )
}

############################################################
## 18. Save canonical malignant-subclustered object
############################################################

canonical_object_file <- file.path(
  out_intermediate_dir,
  "GSE244983_malignant_subclustered.rds"
)

safe_save_rds(
  mal_obj,
  canonical_object_file
)

canonical_object_audit <- data.frame(
  Object = "GSE244983_malignant_subclustered.rds",
  RelativePath = relative_to_project(
    canonical_object_file
  ),
  MD5 = md5_file(
    canonical_object_file
  ),
  N_cells = ncol(
    mal_obj
  ),
  N_features = nrow(
    mal_obj
  ),
  CellSetMD5 = cellset_md5(
    colnames(
      mal_obj
    )
  ),
  N_subclusters =
    length(
      levels(
        mal_obj$MalignantSubcluster
      )
    ),
  Subclusters =
    paste(
      levels(
        mal_obj$MalignantSubcluster
      ),
      collapse = ";"
    ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  canonical_object_audit,
  file.path(
    out_table_dir,
    "GSE244983_canonical_malignant_subclustered_object_audit.csv"
  )
)

############################################################
## 19. Output inventory
############################################################

table_files <- list.files(
  out_table_dir,
  full.names = TRUE,
  recursive = FALSE
)

inventory_paths <- c(
  table_files,
  canonical_object_file,
  session_file
)

output_inventory <- data.frame(
  Type = c(
    rep(
      "table",
      length(
        table_files
      )
    ),
    "RDS",
    "sessionInfo"
  ),
  RelativePath = vapply(
    inventory_paths,
    relative_to_project,
    character(1)
  ),
  File = basename(
    inventory_paths
  ),
  Size_bytes = file.info(
    inventory_paths
  )$size,
  MD5 = vapply(
    inventory_paths,
    md5_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  output_inventory,
  file.path(
    out_table_dir,
    "GSE244983_malignant_subclustering_output_inventory.csv"
  )
)

message(
  "08b_GSE244983_malignant_subclustering.R finished successfully."
)
