############################################################
## 08_GSE244983_CAF_stromal_analysis.R
##
## Purpose
## -------
## Recompute the CAF/stromal marker analysis in the canonical GSE244983
## single-cell lineage produced by script 07, and generate analysis/source
## tables for:
##   - Main Figure 4A-D;
##   - Supplementary Figure S9;
##   - revised Supplementary Table S8.
##
## Public-pipeline design
## ----------------------
## This script:
##   1) consumes ONLY the canonical script-07 major-annotated Seurat object;
##   2) uses the frozen 21-gene CAF/stromal marker panel;
##   3) implements the final Supplementary Methods definition:
##        mean normalized expression of PRESENT marker genes
##        -> z-standardization across cells;
##   4) does NOT use the tumor-dedifferentiation/stromal-remodeling state
##      score as an input, avoiding circular source attribution;
##   5) creates explicit DotPlot-compatible source tables rather than
##      generating manuscript figures inside the analysis script;
##   6) does not recursively discover files or auto-install packages;
##   7) records input MD5, cell-set MD5, hard reproducibility gates,
##      output inventory, and sessionInfo().
##
## Canonical input
## ---------------
## results/intermediate/GSE244983/
##   GSE244983_seurat_major_annotated.rds
##
## Important provenance note
## -------------------------
## Historical script 03 embedded this CAF analysis inside preprocessing,
## while the final Supplementary Figure S9 plotting script independently
## re-read the historical major-annotated object. Public script 09 separates
## the CAF analysis into a single canonical analysis layer. Figure scripts
## should only plot the tables written here.
############################################################

options(stringsAsFactors = FALSE)

SCRIPT_ID <- "09_GSE244983_CAF_stromal_analysis"

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
  "tibble"
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
  "CAF_stromal_analysis"
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

md5_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  unname(
    tools::md5sum(path)
  )
}

cellset_md5 <- function(ids) {
  ids <- sort(
    unique(
      as.character(ids)
    )
  )

  tf <- tempfile(
    fileext = ".txt"
  )

  on.exit(
    unlink(tf),
    add = TRUE
  )

  writeLines(
    ids,
    tf,
    useBytes = TRUE
  )

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
            "Script 09 will not fall back to raw counts because that would ",
            "change the frozen CAF/stromal core-score definition."
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

## Reproduce the numeric data behind Seurat::DotPlot without creating a plot.
##
## For each group and marker:
##   AverageExpressionLinear = mean(expm1(log-normalized expression))
##   PercentExpressed = fraction of cells with normalized expression > 0
##
## As in Seurat DotPlot(scale = TRUE), the average expression is transformed
## by log1p and z-standardized across groups within each gene, then clipped
## to the default DotPlot color range [-2.5, 2.5].
##
## The submitted S9 plotting script subsequently displayed this quantity
## with a visual scale limited to [-2, 2]. Both values are retained here.
build_dotplot_source <- function(
  marker_expr_dense,
  marker_features,
  marker_labels,
  groups,
  group_levels,
  group_col_name
) {
  if (
    nrow(marker_expr_dense) !=
      length(marker_features)
  ) {
    stop(
      "marker_expr_dense row count does not match marker_features."
    )
  }

  if (
    ncol(marker_expr_dense) !=
      length(groups)
  ) {
    stop(
      "marker_expr_dense column count does not match grouping vector."
    )
  }

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
      Group = g,
      N_cells = length(idx),
      Gene = marker_labels,
      MatchedFeature = marker_features,
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

  if (nrow(out) == 0) {
    stop(
      "DotPlot source table has zero rows."
    )
  }

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
        ),
      ScaledAverageExpression_S9Display =
        pmax(
          pmin(
            ScaledAverageExpression_DotPlot,
            2
          ),
          -2
        )
    ) |>
    dplyr::ungroup()

  colnames(out)[
    colnames(out) ==
      "Group"
  ] <- group_col_name

  out
}

############################################################
## 4. Frozen lineage and CAF/stromal definitions
############################################################

expected_n_cells <- 26053L
expected_n_features <- 11616L

expected_cellset_md5 <-
  "a574b6a23054be747b74e2e09f5ddea9"

major_celltype_levels <- c(
  "Malignant",
  "Cycling malignant",
  "CAF/stromal-like cells",
  "Endothelial",
  "Myeloid cells",
  "T/NK cells",
  "T/NK/Treg-like cells",
  "B/Plasma cells",
  "Cycling cells"
)

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

expected_cluster_counts <- c(
  "0" = 4547L,
  "1" = 3896L,
  "2" = 2705L,
  "3" = 2343L,
  "4" = 2205L,
  "5" = 2012L,
  "6" = 1888L,
  "7" = 1565L,
  "8" = 1541L,
  "9" = 1103L,
  "10" = 616L,
  "11" = 382L,
  "12" = 333L,
  "13" = 305L,
  "14" = 231L,
  "15" = 193L,
  "16" = 143L,
  "17" = 45L
)

caf_markers <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL5A1",
  "COL6A1",
  "COL6A2",
  "COL6A3",
  "DCN",
  "LUM",
  "FAP",
  "ACTA2",
  "PDPN",
  "THY1",
  "VCAN",
  "TNC",
  "TGFB1",
  "FBN1",
  "MRC2",
  "ADAM12",
  "TIMP1",
  "INHBA"
)

safe_write_csv(
  data.frame(
    MarkerPanel = "CAF/stromal core",
    Gene = caf_markers,
    Order = seq_along(
      caf_markers
    ),
    stringsAsFactors = FALSE
  ),
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_marker_definitions.csv"
  )
)

############################################################
## 5. Input audit and script-07 gates
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
    "GSE244983_CAF_stromal_input_audit.csv"
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

required_meta <- c(
  "Sample",
  "SeuratCluster",
  "MajorCellType"
)

missing_meta <- setdiff(
  required_meta,
  colnames(
    obj@meta.data
  )
)

if (length(missing_meta) > 0) {
  stop(
    "Canonical script-07 object lacks required metadata column(s): ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

if (
  !"umap" %in%
    names(
      obj@reductions
    )
) {
  stop(
    "Canonical script-07 object lacks UMAP reduction."
  )
}

input_cellset_md5 <- cellset_md5(
  colnames(
    obj
  )
)

observed_major_raw <- table(
  as.character(
    obj$MajorCellType
  )
)

observed_major_counts <- setNames(
  integer(
    length(
      expected_major_counts
    )
  ),
  names(
    expected_major_counts
  )
)

major_common <- intersect(
  names(
    observed_major_raw
  ),
  names(
    expected_major_counts
  )
)

observed_major_counts[
  major_common
] <- as.integer(
  observed_major_raw[
    major_common
  ]
)

observed_cluster_raw <- table(
  as.character(
    obj$SeuratCluster
  )
)

observed_cluster_counts <- setNames(
  integer(
    length(
      expected_cluster_counts
    )
  ),
  names(
    expected_cluster_counts
  )
)

cluster_common <- intersect(
  names(
    observed_cluster_raw
  ),
  names(
    expected_cluster_counts
  )
)

observed_cluster_counts[
  cluster_common
] <- as.integer(
  observed_cluster_raw[
    cluster_common
  ]
)

cluster12_cells <- colnames(
  obj
)[
  as.character(
    obj$SeuratCluster
  ) ==
    "12"
]

cluster12_major_types <- unique(
  as.character(
    obj$MajorCellType[
      as.character(
        obj$SeuratCluster
      ) ==
        "12"
    ]
  )
)

lineage_gates <- data.frame(
  Gate = c(
    "Script-07 object has 26,053 cells",
    "Script-07 object has 11,616 features",
    "Script-07 cell-set MD5 matches frozen primary lineage",
    "MajorCellType has no missing labels",
    "Major-cell-type counts match frozen script-07 lineage",
    "Seurat cluster counts match frozen script-07 lineage",
    "Exactly 18 Seurat clusters are present",
    "Cluster 12 contains 333 cells",
    "Cluster 12 is uniformly CAF/stromal-like",
    "UMAP reduction present"
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
    input_cellset_md5,
    as.character(
      sum(
        is.na(
          obj$MajorCellType
        )
      )
    ),
    paste0(
      names(
        expected_major_counts
      ),
      "=",
      observed_major_counts[
        names(
          expected_major_counts
        )
      ],
      collapse = ";"
    ),
    paste0(
      names(
        expected_cluster_counts
      ),
      "=",
      observed_cluster_counts[
        names(
          expected_cluster_counts
        )
      ],
      collapse = ";"
    ),
    as.character(
      length(
        unique(
          as.character(
            obj$SeuratCluster
          )
        )
      )
    ),
    as.character(
      length(
        cluster12_cells
      )
    ),
    paste(
      cluster12_major_types,
      collapse = ";"
    ),
    as.character(
      "umap" %in%
        names(
          obj@reductions
        )
    )
  ),
  Expected = c(
    as.character(
      expected_n_cells
    ),
    as.character(
      expected_n_features
    ),
    expected_cellset_md5,
    "0",
    paste0(
      names(
        expected_major_counts
      ),
      "=",
      expected_major_counts,
      collapse = ";"
    ),
    paste0(
      names(
        expected_cluster_counts
      ),
      "=",
      expected_cluster_counts,
      collapse = ";"
    ),
    "18",
    "333",
    "CAF/stromal-like cells",
    "TRUE"
  ),
  Pass = c(
    ncol(
      obj
    ) ==
      expected_n_cells,
    nrow(
      obj
    ) ==
      expected_n_features,
    identical(
      input_cellset_md5,
      expected_cellset_md5
    ),
    sum(
      is.na(
        obj$MajorCellType
      )
    ) ==
      0,
    identical(
      as.integer(
        observed_major_counts[
          names(
            expected_major_counts
          )
        ]
      ),
      as.integer(
        expected_major_counts
      )
    ),
    identical(
      as.integer(
        observed_cluster_counts[
          names(
            expected_cluster_counts
          )
        ]
      ),
      as.integer(
        expected_cluster_counts
      )
    ),
    length(
      unique(
        as.character(
          obj$SeuratCluster
        )
      )
    ) ==
      18,
    length(
      cluster12_cells
    ) ==
      333,
    identical(
      cluster12_major_types,
      "CAF/stromal-like cells"
    ),
    "umap" %in%
      names(
        obj@reductions
      )
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  lineage_gates,
  file.path(
    out_table_dir,
    "GSE244983_script07_lineage_gates.csv"
  )
)

if (
  any(
    !lineage_gates$Pass
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

  stop(
    "Script-07 lineage gate failure: ",
    paste(
      lineage_gates$Gate[
        !lineage_gates$Pass
      ],
      collapse = " | "
    ),
    ". CAF/stromal scoring was not run."
  )
}

############################################################
## 6. Marker presence and normalized expression
############################################################

expr_data <- get_normalized_rna(
  obj
)

marker_match <- match_features_case_insensitive(
  rownames(
    expr_data
  ),
  caf_markers
)

safe_write_csv(
  marker_match,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_marker_gene_presence.csv"
  )
)

caf_markers_present <- marker_match$MatchedFeature[
  marker_match$Present
]

caf_marker_labels_present <- marker_match$RequestedGene[
  marker_match$Present
]

message(
  "CAF/stromal markers detected: ",
  length(
    caf_markers_present
  ),
  "/",
  length(
    caf_markers
  )
)

if (
  length(
    caf_markers_present
  ) <
    2
) {
  stop(
    "Too few CAF/stromal markers are present to calculate the core score."
  )
}

############################################################
## 7. CAF/stromal core score
##
## Final Supplementary Methods definition:
## mean normalized expression of PRESENT markers -> z-standardize across cells.
############################################################

caf_core_raw <- Matrix::colMeans(
  expr_data[
    caf_markers_present,
    ,
    drop = FALSE
  ]
)

caf_core_z <- zscore_safe(
  caf_core_raw
)

if (
  length(
    caf_core_z
  ) !=
    ncol(
      obj
    )
) {
  stop(
    "CAF/stromal core score length does not match object cell count."
  )
}

caf_core_cell_scores <- data.frame(
  Cell = colnames(
    obj
  ),
  Sample = as.character(
    obj$Sample
  ),
  SeuratCluster = as.character(
    obj$SeuratCluster
  ),
  MajorCellType = as.character(
    obj$MajorCellType
  ),
  CAF_stromal_core_score_raw =
    as.numeric(
      caf_core_raw
    ),
  CAF_stromal_core_score_z =
    as.numeric(
      caf_core_z
    ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  caf_core_cell_scores,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_core_score_cell_level.csv"
  )
)

score_distribution_audit <- data.frame(
  N_cells = length(
    caf_core_z
  ),
  N_finite = sum(
    is.finite(
      caf_core_z
    )
  ),
  Mean = mean(
    caf_core_z,
    na.rm = TRUE
  ),
  SD = stats::sd(
    caf_core_z,
    na.rm = TRUE
  ),
  Min = min(
    caf_core_z,
    na.rm = TRUE
  ),
  Q25 = as.numeric(
    stats::quantile(
      caf_core_z,
      0.25,
      na.rm = TRUE
    )
  ),
  Median = stats::median(
    caf_core_z,
    na.rm = TRUE
  ),
  Q75 = as.numeric(
    stats::quantile(
      caf_core_z,
      0.75,
      na.rm = TRUE
    )
  ),
  Max = max(
    caf_core_z,
    na.rm = TRUE
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  score_distribution_audit,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_core_score_distribution_audit.csv"
  )
)

############################################################
## 8. Major-cell-type and cluster summaries
############################################################

major_summary <- caf_core_cell_scores |>
  dplyr::mutate(
    MajorCellType = factor(
      MajorCellType,
      levels = major_celltype_levels
    )
  ) |>
  dplyr::group_by(
    MajorCellType,
    .drop = TRUE
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    CAF_stromal_core_score_z_mean =
      mean(
        CAF_stromal_core_score_z,
        na.rm = TRUE
      ),
    CAF_stromal_core_score_z_median =
      stats::median(
        CAF_stromal_core_score_z,
        na.rm = TRUE
      ),
    CAF_stromal_core_score_z_q25 =
      as.numeric(
        stats::quantile(
          CAF_stromal_core_score_z,
          0.25,
          na.rm = TRUE
        )
      ),
    CAF_stromal_core_score_z_q75 =
      as.numeric(
        stats::quantile(
          CAF_stromal_core_score_z,
          0.75,
          na.rm = TRUE
        )
      ),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    match(
      as.character(
        MajorCellType
      ),
      major_celltype_levels
    )
  )

major_summary$MajorCellType <-
  as.character(
    major_summary$MajorCellType
  )

safe_write_csv(
  major_summary,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_core_score_by_major_celltype.csv"
  )
)

cluster_summary <- caf_core_cell_scores |>
  dplyr::group_by(
    SeuratCluster,
    MajorCellType
  ) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    CAF_stromal_core_score_z_mean =
      mean(
        CAF_stromal_core_score_z,
        na.rm = TRUE
      ),
    CAF_stromal_core_score_z_median =
      stats::median(
        CAF_stromal_core_score_z,
        na.rm = TRUE
      ),
    CAF_stromal_core_score_z_q25 =
      as.numeric(
        stats::quantile(
          CAF_stromal_core_score_z,
          0.25,
          na.rm = TRUE
        )
      ),
    CAF_stromal_core_score_z_q75 =
      as.numeric(
        stats::quantile(
          CAF_stromal_core_score_z,
          0.75,
          na.rm = TRUE
        )
      ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    ClusterNumeric =
      suppressWarnings(
        as.numeric(
          SeuratCluster
        )
      )
  ) |>
  dplyr::arrange(
    ClusterNumeric,
    SeuratCluster
  ) |>
  dplyr::select(
    -ClusterNumeric
  )

safe_write_csv(
  cluster_summary,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_core_score_by_Seurat_cluster.csv"
  )
)

cluster_identity_audit <- cluster_summary |>
  dplyr::arrange(
    suppressWarnings(
      as.numeric(
        SeuratCluster
      )
    )
  )

safe_write_csv(
  cluster_identity_audit,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_cluster_identity_audit.csv"
  )
)

############################################################
## 9. UMAP source for Main Figure 4B
############################################################

umap_mat <- as.data.frame(
  Seurat::Embeddings(
    obj,
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
    "UMAP reduction has fewer than two dimensions."
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
  caf_core_cell_scores,
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
    "GSE244983_CAF_stromal_core_score_UMAP_source.csv"
  )
)

############################################################
## 10. DotPlot-compatible source tables for Figure 4A / S9
############################################################

## Only 21 genes x 26,053 cells are converted to dense form.
## This is intentionally small and makes the DotPlot statistics explicit.
marker_expr_dense <- as.matrix(
  expr_data[
    caf_markers_present,
    ,
    drop = FALSE
  ]
)

rownames(
  marker_expr_dense
) <- caf_markers_present

major_levels_observed <- major_celltype_levels[
  major_celltype_levels %in%
    unique(
      as.character(
        obj$MajorCellType
      )
    )
]

cluster_levels_observed <- names(
  expected_cluster_counts
)

dotplot_major <- build_dotplot_source(
  marker_expr_dense = marker_expr_dense,
  marker_features = caf_markers_present,
  marker_labels = caf_marker_labels_present,
  groups = as.character(
    obj$MajorCellType
  ),
  group_levels = major_levels_observed,
  group_col_name = "MajorCellType"
)

dotplot_major$MajorCellType <- factor(
  dotplot_major$MajorCellType,
  levels = major_levels_observed
)

dotplot_major$Gene <- factor(
  dotplot_major$Gene,
  levels = caf_marker_labels_present
)

dotplot_major <- dotplot_major |>
  dplyr::arrange(
    MajorCellType,
    Gene
  )

dotplot_major$MajorCellType <-
  as.character(
    dotplot_major$MajorCellType
  )

dotplot_major$Gene <-
  as.character(
    dotplot_major$Gene
  )

safe_write_csv(
  dotplot_major,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_dotplot_by_major_celltype.csv"
  )
)

dotplot_cluster <- build_dotplot_source(
  marker_expr_dense = marker_expr_dense,
  marker_features = caf_markers_present,
  marker_labels = caf_marker_labels_present,
  groups = as.character(
    obj$SeuratCluster
  ),
  group_levels = cluster_levels_observed,
  group_col_name = "SeuratCluster"
)

dotplot_cluster$SeuratCluster <- factor(
  dotplot_cluster$SeuratCluster,
  levels = cluster_levels_observed
)

dotplot_cluster$Gene <- factor(
  dotplot_cluster$Gene,
  levels = caf_marker_labels_present
)

dotplot_cluster <- dotplot_cluster |>
  dplyr::arrange(
    SeuratCluster,
    Gene
  )

dotplot_cluster$SeuratCluster <-
  as.character(
    dotplot_cluster$SeuratCluster
  )

dotplot_cluster$Gene <-
  as.character(
    dotplot_cluster$Gene
  )

safe_write_csv(
  dotplot_cluster,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_dotplot_by_Seurat_cluster.csv"
  )
)

rm(
  marker_expr_dense
)

gc(
  verbose = FALSE
)

############################################################
## 11. Figure/source-data manifest
############################################################

source_manifest <- data.frame(
  Artifact = c(
    "Main Figure 4A",
    "Main Figure 4B",
    "Main Figure 4C",
    "Main Figure 4D",
    "Supplementary Figure S9A",
    "Supplementary Figure S9B",
    "Supplementary Table S8"
  ),
  SourceTable = c(
    "GSE244983_CAF_stromal_dotplot_by_Seurat_cluster.csv",
    "GSE244983_CAF_stromal_core_score_UMAP_source.csv",
    "GSE244983_CAF_stromal_core_score_cell_level.csv",
    "GSE244983_CAF_stromal_core_score_cell_level.csv",
    "GSE244983_CAF_stromal_dotplot_by_major_celltype.csv",
    "GSE244983_CAF_stromal_dotplot_by_Seurat_cluster.csv",
    "GSE244983_CAF_stromal_core_score_by_major_celltype.csv"
  ),
  PlotRole = c(
    "CAF/stromal marker dot plot across Seurat clusters",
    "CAF/stromal core-score UMAP",
    "CAF/stromal core-score distribution by major cell type",
    "CAF/stromal core-score distribution by Seurat cluster",
    "CAF/stromal marker dot plot across major cell types",
    "CAF/stromal marker dot plot across Seurat clusters",
    "CAF/stromal core-score summary by major cell type"
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  source_manifest,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_figure_source_manifest.csv"
  )
)

############################################################
## 12. Reproducibility gates
############################################################

major_mean_order <- major_summary[
  order(
    -major_summary$CAF_stromal_core_score_z_mean
  ),
  ,
  drop = FALSE
]

cluster_mean_order <- cluster_summary[
  order(
    -cluster_summary$CAF_stromal_core_score_z_mean
  ),
  ,
  drop = FALSE
]

top_major_type <- major_mean_order$MajorCellType[[1]]
top_cluster <- cluster_mean_order$SeuratCluster[[1]]

major_n_cells_total <- sum(
  major_summary$n_cells
)

cluster_n_cells_total <- sum(
  cluster_summary$n_cells
)

cell_score_md5 <- cellset_md5(
  caf_core_cell_scores$Cell
)

umap_score_md5 <- cellset_md5(
  umap_source$Cell
)

expected_dotplot_major_rows <-
  length(
    major_levels_observed
  ) *
  length(
    caf_markers_present
  )

expected_dotplot_cluster_rows <-
  length(
    cluster_levels_observed
  ) *
  length(
    caf_markers_present
  )

analysis_gates <- data.frame(
  Gate = c(
    "All 21 frozen CAF/stromal markers are present",
    "CAF/stromal cell-level score covers all 26,053 cells",
    "CAF/stromal cell-score cell set matches script-07 lineage",
    "All CAF/stromal z-scores are finite",
    "CAF/stromal z-score mean is zero",
    "CAF/stromal z-score SD is one",
    "Major-cell-type summary accounts for all cells",
    "Cluster summary accounts for all cells",
    "CAF/stromal-like cells have the highest mean CAF/stromal core score",
    "Cluster 12 has the highest mean CAF/stromal core score",
    "UMAP score source covers the exact script-07 cell set",
    "Major-cell-type DotPlot source has expected rows",
    "Cluster DotPlot source has expected rows",
    "Major-cell-type DotPlot source has no missing scaled averages",
    "Cluster DotPlot source has no missing scaled averages"
  ),
  Observed = c(
    as.character(
      sum(
        marker_match$Present
      )
    ),
    as.character(
      nrow(
        caf_core_cell_scores
      )
    ),
    cell_score_md5,
    as.character(
      sum(
        is.finite(
          caf_core_cell_scores$CAF_stromal_core_score_z
        )
      )
    ),
    as.character(
      mean(
        caf_core_cell_scores$CAF_stromal_core_score_z
      )
    ),
    as.character(
      stats::sd(
        caf_core_cell_scores$CAF_stromal_core_score_z
      )
    ),
    as.character(
      major_n_cells_total
    ),
    as.character(
      cluster_n_cells_total
    ),
    top_major_type,
    top_cluster,
    umap_score_md5,
    as.character(
      nrow(
        dotplot_major
      )
    ),
    as.character(
      nrow(
        dotplot_cluster
      )
    ),
    as.character(
      sum(
        is.na(
          dotplot_major$ScaledAverageExpression_DotPlot
        )
      )
    ),
    as.character(
      sum(
        is.na(
          dotplot_cluster$ScaledAverageExpression_DotPlot
        )
      )
    )
  ),
  Expected = c(
    "21",
    as.character(
      expected_n_cells
    ),
    expected_cellset_md5,
    as.character(
      expected_n_cells
    ),
    "0",
    "1",
    as.character(
      expected_n_cells
    ),
    as.character(
      expected_n_cells
    ),
    "CAF/stromal-like cells",
    "12",
    expected_cellset_md5,
    as.character(
      expected_dotplot_major_rows
    ),
    as.character(
      expected_dotplot_cluster_rows
    ),
    "0",
    "0"
  ),
  Pass = c(
    sum(
      marker_match$Present
    ) ==
      21,
    nrow(
      caf_core_cell_scores
    ) ==
      expected_n_cells,
    identical(
      cell_score_md5,
      expected_cellset_md5
    ),
    sum(
      is.finite(
        caf_core_cell_scores$CAF_stromal_core_score_z
      )
    ) ==
      expected_n_cells,
    abs(
      mean(
        caf_core_cell_scores$CAF_stromal_core_score_z
      )
    ) <
      1e-8,
    abs(
      stats::sd(
        caf_core_cell_scores$CAF_stromal_core_score_z
      ) -
        1
    ) <
      1e-8,
    major_n_cells_total ==
      expected_n_cells,
    cluster_n_cells_total ==
      expected_n_cells,
    identical(
      top_major_type,
      "CAF/stromal-like cells"
    ),
    identical(
      top_cluster,
      "12"
    ),
    identical(
      umap_score_md5,
      expected_cellset_md5
    ),
    nrow(
      dotplot_major
    ) ==
      expected_dotplot_major_rows,
    nrow(
      dotplot_cluster
    ) ==
      expected_dotplot_cluster_rows,
    sum(
      is.na(
        dotplot_major$ScaledAverageExpression_DotPlot
      )
    ) ==
      0,
    sum(
      is.na(
        dotplot_cluster$ScaledAverageExpression_DotPlot
      )
    ) ==
      0
  ),
  stringsAsFactors = FALSE
)

all_gates <- dplyr::bind_rows(
  data.frame(
    GateClass = "script07_lineage",
    lineage_gates,
    stringsAsFactors = FALSE
  ),
  data.frame(
    GateClass = "CAF_stromal_analysis",
    analysis_gates,
    stringsAsFactors = FALSE
  )
)

safe_write_csv(
  all_gates,
  file.path(
    out_table_dir,
    "GSE244983_CAF_stromal_reproducibility_gates.csv"
  )
)

############################################################
## 13. Session info
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
    "09 reproducibility gate failure: ",
    paste(
      all_gates$Gate[
        !all_gates$Pass
      ],
      collapse = " | "
    ),
    ". No downstream figure script should be frozen from these outputs."
  )
}

############################################################
## 14. Output inventory
############################################################

table_files <- list.files(
  out_table_dir,
  full.names = TRUE,
  recursive = FALSE
)

inventory_paths <- c(
  table_files,
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
    "GSE244983_CAF_stromal_output_inventory.csv"
  )
)

message(
  "08_GSE244983_CAF_stromal_analysis.R finished successfully."
)
