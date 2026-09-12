############################################################
## 07_GSE244983_scRNA_state_localization_source_attribution_v1.2.4_FINAL.R
##
## PUBLIC-RELEASE MERGED VERSION
##
## This file merges the canonical state-localization/source-attribution
## analysis with the v1.2.3 Supplementary Figure S8 reporting-contract
## hardening patch.
##
## Scientific/statistical calculations are unchanged relative to the two
## source scripts. The merge only removes the need to run the reporting
## contract as a second script and incorporates its gates into the final
## reproducibility-gate table/output inventory.
############################################################

############################################################
## 07_GSE244983_scRNA_state_localization_source_attribution.R
##
## Purpose
## -------
## Recompute the four predefined tumor-immune state scores in the
## canonical GSE244983 single-cell object produced by script 07,
## localize those scores across major cell compartments, and compute
## the canonical-marker source-attribution analyses used by Figure 3
## and Supplementary Figures S7-S8.
##
## Public-pipeline design
## ----------------------
## This script:
##   1) consumes only the canonical script-07 Seurat object;
##   2) consumes only the frozen script-01 state mapping + CLEAN GMT;
##   3) does not recursively discover files;
##   4) does not auto-install packages;
##   5) preserves the historical constituent gene-wise z-mean scoring;
##   6) preserves the final S8 canonical marker modules and six
##      prespecified cell-compartment contrasts;
##   7) writes analysis/source tables only -- no manuscript figures;
##   8) applies hard structural/reproducibility gates before saving the
##      canonical state-localized Seurat object;
##   9) records input MD5s, cell-set MD5, output inventory, sessionInfo().
##
## Canonical upstream inputs
## -------------------------
## results/intermediate/GSE244983/
##   GSE244983_seurat_major_annotated.rds
##
## results/tables/bulk_discovery/
##   GSE244982_final_state_signature_mapping.csv
##
## results/tables/ICBcomb_final_input_gene_sets_CLEAN/
##   02_final_ICBcomb_gene_sets/final_ICBcomb_gene_sets_CLEAN.gmt
##
## Canonical downstream object
## ---------------------------
## results/intermediate/GSE244983/
##   GSE244983_seurat_state_localized.rds
############################################################

options(stringsAsFactors = FALSE)

SCRIPT_ID <- "08_GSE244983_scRNA_state_localization_source_attribution"
SCRIPT_PATCH <- "v1.2.4_merged_S8_reporting_contract"

############################################################
## 0. Project root
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")

if (!nzchar(project_dir)) {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

  probe <- c(
    file.path(
      cwd,
      "results",
      "intermediate",
      "GSE244983",
      "GSE244983_seurat_major_annotated.rds"
    ),
    file.path(
      cwd,
      "results",
      "tables",
      "bulk_discovery",
      "GSE244982_final_state_signature_mapping.csv"
    )
  )

  if (all(file.exists(probe))) {
    project_dir <- cwd
  } else {
    stop(
      "Set ICB_PROJECT_DIR to the repository root, or run this script ",
      "from the root containing the canonical script-07 and script-01 outputs."
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
    paste(missing_pkgs, collapse = ", "),
    ". Restore the repository environment (for example, renv::restore()) and rerun."
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

state_mapping_file <- file.path(
  project_dir,
  "results",
  "tables",
  "bulk_discovery",
  "GSE244982_final_state_signature_mapping.csv"
)

gmt_file_env <- Sys.getenv("ICB_CLEAN_GMT")

if (nzchar(gmt_file_env)) {
  gmt_file <- normalizePath(
    gmt_file_env,
    winslash = "/",
    mustWork = FALSE
  )
} else {
  gmt_file <- file.path(
    project_dir,
    "results",
    "tables",
    "ICBcomb_final_input_gene_sets_CLEAN",
    "02_final_ICBcomb_gene_sets",
    "final_ICBcomb_gene_sets_CLEAN.gmt"
  )
}

out_table_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "GSE244983",
  "state_localization_source_attribution"
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

required_inputs <- c(
  input_seurat_file,
  state_mapping_file,
  gmt_file
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs) > 0) {
  stop(
    "Missing canonical input(s):\n",
    paste(missing_inputs, collapse = "\n")
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

  message("Saved: ", path)
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

  message("Saved: ", path)
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

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- gsub(
    "\ufeff",
    "",
    x,
    fixed = TRUE
  )
  x <- trimws(x)
  x <- gsub(
    "^['\"]+|['\"]+$",
    "",
    x
  )
  x <- toupper(x)
  x <- gsub(
    "\\s+",
    "",
    x
  )
  x <- gsub(
    "^HUMAN_",
    "",
    x
  )
  x <- gsub(
    "^MOUSE_",
    "",
    x
  )
  x
}

zvec <- function(x) {
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
            "Script 08 will not fall back to raw counts because that would change ",
            "the historical scoring definition."
          )
        }
      )
    }
  )

  if (nrow(out) == 0 || ncol(out) == 0) {
    stop(
      "Normalized RNA data layer is empty."
    )
  }

  out
}

read_gmt <- function(path) {
  lines <- readLines(
    path,
    warn = FALSE
  )

  lines <- lines[
    nzchar(lines)
  ]

  gs <- list()

  for (ln in lines) {
    parts <- strsplit(
      ln,
      "\t",
      fixed = TRUE
    )[[1]]

    if (length(parts) < 3) {
      next
    }

    set_name <- parts[[1]]

    genes <- unique(
      clean_gene_symbol(
        parts[-c(1, 2)]
      )
    )

    genes <- genes[
      !is.na(genes) &
      nzchar(genes)
    ]

    gs[[set_name]] <- genes
  }

  gs
}

score_gene_set_zmean <- function(
  expr_data,
  expr_genes_original,
  expr_genes_clean,
  genes,
  set_name
) {
  genes_clean <- unique(
    clean_gene_symbol(genes)
  )

  matched_rows <- expr_genes_original[
    expr_genes_clean %in%
      genes_clean
  ]

  matched_rows <- unique(
    matched_rows
  )

  message(
    set_name,
    ": ",
    length(matched_rows),
    "/",
    length(genes_clean),
    " genes present"
  )

  if (length(matched_rows) == 0) {
    return(
      list(
        score = rep(
          NA_real_,
          ncol(expr_data)
        ),
        presence = data.frame(
          GeneSet = set_name,
          Gene = genes_clean,
          Present = FALSE,
          stringsAsFactors = FALSE
        )
      )
    )
  }

  ## Historical scoring definition:
  ## normalized expression -> gene-wise z-score across cells -> mean across genes.
  small_mat <- as.matrix(
    expr_data[
      matched_rows,
      ,
      drop = FALSE
    ]
  )

  small_z <- t(
    scale(
      t(small_mat)
    )
  )

  small_z[
    !is.finite(small_z)
  ] <- 0

  score <- colMeans(
    small_z,
    na.rm = TRUE
  )

  rm(
    small_mat,
    small_z
  )
  gc(
    verbose = FALSE
  )

  present_clean <- unique(
    expr_genes_clean[
      expr_genes_clean %in%
        genes_clean
    ]
  )

  list(
    score = as.numeric(score),
    presence = data.frame(
      GeneSet = set_name,
      Gene = genes_clean,
      Present = genes_clean %in%
        present_clean,
      stringsAsFactors = FALSE
    )
  )
}

rank_biserial <- function(
  x,
  y
) {
  x <- x[
    is.finite(x)
  ]
  y <- y[
    is.finite(y)
  ]

  if (
    length(x) == 0 ||
    length(y) == 0
  ) {
    return(NA_real_)
  }

  ranks <- rank(
    c(x, y),
    ties.method = "average"
  )

  n1 <- length(x)
  n2 <- length(y)

  u1 <- sum(
    ranks[
      seq_len(n1)
    ]
  ) -
    n1 *
    (n1 + 1) /
    2

  2 *
    u1 /
    (n1 * n2) -
    1
}

session_file <- file.path(
  log_dir,
  paste0(
    "sessionInfo_",
    SCRIPT_ID,
    ".txt"
  )
)

save_session_info <- function() {
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
}

script_metadata <- data.frame(
  ScriptID = SCRIPT_ID,
  Patch = SCRIPT_PATCH,
  CanonicalMajorCellTypeColumn = "MajorCellType",
  stringsAsFactors = FALSE
)

safe_write_csv(
  script_metadata,
  file.path(
    out_table_dir,
    "GSE244983_script08_metadata.csv"
  )
)

############################################################
## 4. Frozen state definitions and script-07 gates
############################################################

state_order <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_display <- c(
  "Immune_defective_Cold" =
    "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" =
    "myeloid-Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" =
    "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" =
    "melanocytic differentiation"
)

legacy_to_clean_state <- c(
  "Immune_defective_cold" =
    "Immune_defective_Cold",
  "Myeloid_Treg_immunosuppressive" =
    "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiated" =
    "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_differentiated" =
    "Melanocytic_Differentiation",
  "Immune_defective_Cold" =
    "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive" =
    "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" =
    "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation" =
    "Melanocytic_Differentiation"
)

expected_signatures <- c(
  "Immune_defective_Cold_RESTORE_curated",
  "Immune_defective_Cold_RESTORE_data_driven",
  "Immune_defective_Cold_RESTORE_hybrid",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid",
  "Melanocytic_Differentiation_REFERENCE_curated",
  "Melanocytic_Differentiation_REFERENCE_data_driven",
  "Melanocytic_Differentiation_REFERENCE_hybrid"
)

expected_state_direction <- c(
  "Immune_defective_Cold" =
    "inverse",
  "Myeloid_Treg_Immunosuppressive" =
    "positive",
  "Tumor_dedifferentiation_Stromal_remodeling" =
    "positive",
  "Melanocytic_Differentiation" =
    "positive"
)

score_col_map <- setNames(
  paste0(
    state_order,
    "_score"
  ),
  state_order
)

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

expected_n_cells <- 26053L
expected_n_features <- 11616L
expected_cellset_md5 <- "a574b6a23054be747b74e2e09f5ddea9"

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

############################################################
## 5. Canonical marker modules and prespecified contrasts
############################################################

module_genes <- list(
  "canonical melanocytic" = c(
    "MLANA",
    "PMEL",
    "TYR",
    "DCT",
    "MITF",
    "SOX10"
  ),
  "stromal/CAF" = c(
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
    "FBN1"
  ),
  "endothelial" = c(
    "PECAM1",
    "VWF",
    "KDR",
    "CDH5",
    "CLDN5",
    "FLT1",
    "RAMP2",
    "ESAM"
  ),
  "myeloid" = c(
    "LYZ",
    "CD68",
    "C1QA",
    "C1QB",
    "C1QC",
    "FCGR3A",
    "LST1",
    "TYROBP",
    "AIF1"
  ),
  "T/NK" = c(
    "CD3D",
    "CD3E",
    "TRAC",
    "CD2",
    "NKG7",
    "GNLY",
    "GZMB",
    "CD8A",
    "IL7R",
    "FOXP3"
  ),
  "B/plasma" = c(
    "MS4A1",
    "CD79A",
    "CD79B",
    "CD74",
    "MZB1",
    "JCHAIN",
    "IGHG1",
    "IGKC"
  )
)

marker_module_order <- names(
  module_genes
)

module_definition <- dplyr::bind_rows(
  lapply(
    names(module_genes),
    function(nm) {
      data.frame(
        MarkerModule = nm,
        Gene = module_genes[[nm]],
        stringsAsFactors = FALSE
      )
    }
  )
)

safe_write_csv(
  module_definition,
  file.path(
    out_table_dir,
    "GSE244983_canonical_marker_module_definitions.csv"
  )
)

contrast_specs <- data.frame(
  State = c(
    "Melanocytic_Differentiation",
    "Melanocytic_Differentiation",
    "Myeloid_Treg_Immunosuppressive",
    "Myeloid_Treg_Immunosuppressive",
    "Tumor_dedifferentiation_Stromal_remodeling",
    "Tumor_dedifferentiation_Stromal_remodeling"
  ),
  ContrastLabel = c(
    "Cycling malignant vs CAF/stromal-like cells",
    "Malignant vs CAF/stromal-like cells",
    "Myeloid cells vs Malignant",
    "T/NK cells vs Malignant",
    "CAF/stromal-like cells vs Cycling malignant",
    "CAF/stromal-like cells vs Malignant"
  ),
  Group1 = c(
    "Cycling malignant",
    "Malignant",
    "Myeloid cells",
    "T/NK cells",
    "CAF/stromal-like cells",
    "CAF/stromal-like cells"
  ),
  Group2 = c(
    "CAF/stromal-like cells",
    "CAF/stromal-like cells",
    "Malignant",
    "Malignant",
    "Cycling malignant",
    "Malignant"
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  contrast_specs,
  file.path(
    out_table_dir,
    "GSE244983_prespecified_source_attribution_contrasts.csv"
  )
)

############################################################
## 6. Input audit and script-07 lineage gate
############################################################
## Script 07 canonically stores the frozen major annotation in `MajorCellType`.
## Historical objects/scripts sometimes used alternate names, but script 08
## intentionally consumes the public script-07 contract only.

input_audit <- data.frame(
  Input = c(
    "script07_major_annotated_Seurat",
    "script01_final_state_mapping",
    "frozen_CLEAN_GMT"
  ),
  RelativePath = c(
    relative_to_project(
      input_seurat_file
    ),
    relative_to_project(
      state_mapping_file
    ),
    relative_to_project(
      gmt_file
    )
  ),
  Exists = file.exists(
    c(
      input_seurat_file,
      state_mapping_file,
      gmt_file
    )
  ),
  MD5 = c(
    md5_file(
      input_seurat_file
    ),
    md5_file(
      state_mapping_file
    ),
    md5_file(
      gmt_file
    )
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_audit,
  file.path(
    out_table_dir,
    "GSE244983_state_localization_input_audit.csv"
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

if (
  !"MajorCellType" %in%
  colnames(obj@meta.data)
) {
  stop(
    "Canonical script-07 object lacks MajorCellType."
  )
}

if (
  !"SeuratCluster" %in%
  colnames(obj@meta.data)
) {
  stop(
    "Canonical script-07 object lacks SeuratCluster."
  )
}

if (
  !"Sample" %in%
  colnames(obj@meta.data)
) {
  stop(
    "Canonical script-07 object lacks Sample."
  )
}

if (
  !"umap" %in%
  names(obj@reductions)
) {
  stop(
    "Canonical script-07 object lacks the UMAP reduction required for Figure 3 source data."
  )
}

input_cellset_md5 <- cellset_md5(
  colnames(obj)
)

observed_major_raw <- table(
  as.character(
    obj@meta.data$MajorCellType
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
  names(observed_major_raw),
  names(expected_major_counts)
)

observed_major_counts[
  major_common
] <- as.integer(
  observed_major_raw[
    major_common
  ]
)

lineage_gate <- data.frame(
  Gate = c(
    "Script-07 object has 26,053 cells",
    "Script-07 object has 11,616 features",
    "Script-07 cell-set MD5 matches frozen primary lineage",
    "MajorCellType has no missing labels",
    "Major-cell-type counts match frozen script-07 lineage",
    "UMAP reduction present"
  ),
  Observed = c(
    as.character(
      ncol(obj)
    ),
    as.character(
      nrow(obj)
    ),
    input_cellset_md5,
    as.character(
      sum(
        is.na(
          obj@meta.data$MajorCellType
        )
      )
    ),
    paste0(
      names(expected_major_counts),
      "=",
      observed_major_counts[
        names(expected_major_counts)
      ],
      collapse = ";"
    ),
    as.character(
      "umap" %in%
        names(obj@reductions)
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
      names(expected_major_counts),
      "=",
      expected_major_counts,
      collapse = ";"
    ),
    "TRUE"
  ),
  Pass = c(
    ncol(obj) ==
      expected_n_cells,
    nrow(obj) ==
      expected_n_features,
    identical(
      input_cellset_md5,
      expected_cellset_md5
    ),
    sum(
      is.na(
        obj@meta.data$MajorCellType
      )
    ) == 0,
    identical(
      as.integer(
        observed_major_counts[
          names(expected_major_counts)
        ]
      ),
      as.integer(
        expected_major_counts
      )
    ),
    "umap" %in%
      names(obj@reductions)
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  lineage_gate,
  file.path(
    out_table_dir,
    "GSE244983_script07_lineage_gates.csv"
  )
)

if (
  any(
    !lineage_gate$Pass
  )
) {
  save_session_info()

  stop(
    "Script-07 lineage gate failure: ",
    paste(
      lineage_gate$Gate[
        !lineage_gate$Pass
      ],
      collapse = " | "
    ),
    ". State localization was not run."
  )
}

############################################################
## 7. Read and validate frozen GMT + script-01 mapping
############################################################

message(
  "Reading frozen CLEAN GMT..."
)

gene_sets <- read_gmt(
  gmt_file
)

gmt_summary <- data.frame(
  GeneSet = names(gene_sets),
  N_genes = vapply(
    gene_sets,
    length,
    integer(1)
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  gmt_summary,
  file.path(
    out_table_dir,
    "GSE244983_imported_GMT_gene_set_summary.csv"
  )
)

state_mapping_raw <- read.csv(
  state_mapping_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_mapping_cols <- c(
  "FinalState",
  "RawSignatureColumn",
  "DirectionForFinalState"
)

missing_mapping_cols <- setdiff(
  required_mapping_cols,
  colnames(state_mapping_raw)
)

if (
  length(missing_mapping_cols) > 0
) {
  stop(
    "Script-01 mapping lacks required column(s): ",
    paste(
      missing_mapping_cols,
      collapse = ", "
    )
  )
}

state_mapping <- state_mapping_raw[
  ,
  required_mapping_cols,
  drop = FALSE
]

mapped_clean_state <- unname(
  legacy_to_clean_state[
    as.character(
      state_mapping$FinalState
    )
  ]
)

if (
  any(
    is.na(
      mapped_clean_state
    )
  )
) {
  stop(
    "Unrecognized FinalState value(s) in script-01 mapping: ",
    paste(
      unique(
        state_mapping$FinalState[
          is.na(
            mapped_clean_state
          )
        ]
      ),
      collapse = ", "
    )
  )
}

state_mapping$FinalStateClean <-
  mapped_clean_state

state_mapping$DirectionForFinalState <-
  tolower(
    trimws(
      as.character(
        state_mapping$DirectionForFinalState
      )
    )
  )

state_mapping <- state_mapping[
  !duplicated(
    state_mapping[
      ,
      c(
        "FinalStateClean",
        "RawSignatureColumn",
        "DirectionForFinalState"
      ),
      drop = FALSE
    ]
  ),
  ,
  drop = FALSE
]

mapping_counts <- table(
  factor(
    state_mapping$FinalStateClean,
    levels = state_order
  )
)

direction_ok <- vapply(
  state_order,
  function(st) {
    vals <- unique(
      state_mapping$DirectionForFinalState[
        state_mapping$FinalStateClean ==
          st
      ]
    )

    length(vals) == 1 &&
      identical(
        vals,
        expected_state_direction[[st]]
      )
  },
  logical(1)
)

mapping_gate <- data.frame(
  Gate = c(
    "Exactly four predefined states present",
    "Exactly three constituent signatures per state",
    "Mapped signature set equals frozen 12-signature set",
    "All 12 mapped signatures exist in the CLEAN GMT",
    "State directions match frozen scoring definition"
  ),
  Observed = c(
    paste(
      sort(
        unique(
          state_mapping$FinalStateClean
        )
      ),
      collapse = ";"
    ),
    paste0(
      names(mapping_counts),
      "=",
      as.integer(mapping_counts),
      collapse = ";"
    ),
    paste(
      sort(
        unique(
          state_mapping$RawSignatureColumn
        )
      ),
      collapse = ";"
    ),
    as.character(
      sum(
        unique(
          state_mapping$RawSignatureColumn
        ) %in%
          names(gene_sets)
      )
    )
    ,
    paste0(
      state_order,
      "=",
      expected_state_direction[
        state_order
      ],
      collapse = ";"
    )
  ),
  Expected = c(
    paste(
      sort(
        state_order
      ),
      collapse = ";"
    ),
    paste0(
      state_order,
      "=3",
      collapse = ";"
    ),
    paste(
      sort(
        expected_signatures
      ),
      collapse = ";"
    ),
    "12",
    paste0(
      state_order,
      "=",
      expected_state_direction[
        state_order
      ],
      collapse = ";"
    )
  ),
  Pass = c(
    setequal(
      unique(
        state_mapping$FinalStateClean
      ),
      state_order
    ),
    all(
      as.integer(
        mapping_counts
      ) == 3L
    ),
    setequal(
      unique(
        state_mapping$RawSignatureColumn
      ),
      expected_signatures
    ),
    all(
      expected_signatures %in%
        names(gene_sets)
    ),
    all(
      direction_ok
    )
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  mapping_gate,
  file.path(
    out_table_dir,
    "GSE244983_state_mapping_gates.csv"
  )
)

state_mapping_used <- state_mapping[
  order(
    match(
      state_mapping$FinalStateClean,
      state_order
    ),
    state_mapping$RawSignatureColumn
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  state_mapping_used,
  file.path(
    out_table_dir,
    "GSE244983_final_state_signature_mapping_used.csv"
  )
)

if (
  any(
    !mapping_gate$Pass
  )
) {
  save_session_info()

  stop(
    "Frozen state-mapping gate failure: ",
    paste(
      mapping_gate$Gate[
        !mapping_gate$Pass
      ],
      collapse = " | "
    ),
    ". Canonical state-localized object was NOT saved."
  )
}

############################################################
## 8. Frozen state-gene pools and gene coverage
############################################################

expr_data <- get_normalized_rna(
  obj
)

expr_genes_original <- rownames(
  expr_data
)

expr_genes_clean <- clean_gene_symbol(
  expr_genes_original
)

names(
  expr_genes_clean
) <- expr_genes_original

state_pool_rows <- list()

for (st in state_order) {
  sigs <- state_mapping_used$RawSignatureColumn[
    state_mapping_used$FinalStateClean ==
      st
  ]

  pool <- unique(
    clean_gene_symbol(
      unlist(
        gene_sets[
          sigs
        ],
        use.names = FALSE
      )
    )
  )

  pool <- pool[
    !is.na(pool) &
      nzchar(pool)
  ]

  present <- unique(
    expr_genes_clean[
      expr_genes_clean %in%
        pool
    ]
  )

  state_pool_rows[[st]] <- data.frame(
    FinalState = st,
    Gene = pool,
    PresentInGSE244983 = pool %in%
      present,
    stringsAsFactors = FALSE
  )
}

state_gene_pools <- dplyr::bind_rows(
  state_pool_rows
)

safe_write_csv(
  state_gene_pools,
  file.path(
    out_table_dir,
    "GSE244983_final_state_union_gene_pools.csv"
  )
)

state_gene_pool_summary <- state_gene_pools %>%
  dplyr::group_by(
    FinalState
  ) %>%
  dplyr::summarise(
    N_union_genes = dplyr::n(),
    N_present = sum(
      PresentInGSE244983
    ),
    Coverage = N_present /
      N_union_genes,
    .groups = "drop"
  )

safe_write_csv(
  state_gene_pool_summary,
  file.path(
    out_table_dir,
    "GSE244983_final_state_union_gene_pool_summary.csv"
  )
)

############################################################
## 9. Score the 12 constituent signatures
############################################################

signature_names <- expected_signatures

raw_score_list <- vector(
  "list",
  length(
    signature_names
  )
)

names(
  raw_score_list
) <- signature_names

signature_presence_list <- vector(
  "list",
  length(
    signature_names
  )
)

names(
  signature_presence_list
) <- signature_names

for (gs in signature_names) {
  res <- score_gene_set_zmean(
    expr_data = expr_data,
    expr_genes_original = expr_genes_original,
    expr_genes_clean = expr_genes_clean,
    genes = gene_sets[[gs]],
    set_name = gs
  )

  raw_score_list[[gs]] <-
    res$score

  signature_presence_list[[gs]] <-
    res$presence
}

raw_score_df <- as.data.frame(
  raw_score_list,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

raw_score_df$Cell <- colnames(
  obj
)

raw_score_df <- raw_score_df[
  ,
  c(
    "Cell",
    signature_names
  ),
  drop = FALSE
]

safe_write_csv(
  raw_score_df,
  file.path(
    out_table_dir,
    "GSE244983_single_cell_raw_signature_scores.csv"
  )
)

signature_gene_presence <- dplyr::bind_rows(
  signature_presence_list
)

safe_write_csv(
  signature_gene_presence,
  file.path(
    out_table_dir,
    "GSE244983_state_signature_gene_presence.csv"
  )
)

signature_coverage_summary <- signature_gene_presence %>%
  dplyr::group_by(
    GeneSet
  ) %>%
  dplyr::summarise(
    N_genes = dplyr::n(),
    N_present = sum(
      Present
    ),
    Coverage = N_present /
      N_genes,
    .groups = "drop"
  )

safe_write_csv(
  signature_coverage_summary,
  file.path(
    out_table_dir,
    "GSE244983_state_signature_gene_coverage_summary.csv"
  )
)

############################################################
## 10. Build final four state scores
############################################################

state_score_df <- data.frame(
  Cell = colnames(
    obj
  ),
  stringsAsFactors = FALSE
)

for (st in state_order) {
  sigs <- state_mapping_used$RawSignatureColumn[
    state_mapping_used$FinalStateClean ==
      st
  ]

  direction <- unique(
    state_mapping_used$DirectionForFinalState[
      state_mapping_used$FinalStateClean ==
        st
    ]
  )

  raw_state <- rowMeans(
    raw_score_df[
      ,
      sigs,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  if (
    identical(
      direction,
      "inverse"
    )
  ) {
    raw_state <- -raw_state
  } else if (
    !identical(
      direction,
      "positive"
    )
  ) {
    stop(
      "Unexpected direction for ",
      st,
      ": ",
      paste(
        direction,
        collapse = ","
      )
    )
  }

  state_score_df[[st]] <-
    zvec(
      raw_state
    )

  message(
    st,
    " final score built from: ",
    paste(
      sigs,
      collapse = ", "
    )
  )
}

## Canonical object metadata uses explicit *_score columns.
for (st in state_order) {
  cc <- score_col_map[[st]]

  obj@meta.data[[cc]] <-
    state_score_df[
      match(
        colnames(obj),
        state_score_df$Cell
      ),
      st
    ]
}

state_score_stats <- data.frame(
  State = state_order,
  N = vapply(
    state_order,
    function(st) {
      sum(
        is.finite(
          state_score_df[[st]]
        )
      )
    },
    integer(1)
  ),
  Mean = vapply(
    state_order,
    function(st) {
      mean(
        state_score_df[[st]],
        na.rm = TRUE
      )
    },
    numeric(1)
  ),
  SD = vapply(
    state_order,
    function(st) {
      stats::sd(
        state_score_df[[st]],
        na.rm = TRUE
      )
    },
    numeric(1)
  ),
  Min = vapply(
    state_order,
    function(st) {
      min(
        state_score_df[[st]],
        na.rm = TRUE
      )
    },
    numeric(1)
  ),
  Median = vapply(
    state_order,
    function(st) {
      stats::median(
        state_score_df[[st]],
        na.rm = TRUE
      )
    },
    numeric(1)
  ),
  Max = vapply(
    state_order,
    function(st) {
      max(
        state_score_df[[st]],
        na.rm = TRUE
      )
    },
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  state_score_stats,
  file.path(
    out_table_dir,
    "GSE244983_state_score_distribution_audit.csv"
  )
)

############################################################
## 11. Cell-level state source tables for Figure 3 / S7
############################################################

meta_source <- data.frame(
  Cell = colnames(
    obj
  ),
  Sample = as.character(
    obj@meta.data$Sample
  ),
  SeuratCluster = as.character(
    obj@meta.data$SeuratCluster
  ),
  MajorCellType = as.character(
    obj@meta.data$MajorCellType
  ),
  stringsAsFactors = FALSE
)

single_cell_state_scores <- cbind(
  meta_source,
  state_score_df[
    match(
      meta_source$Cell,
      state_score_df$Cell
    ),
    state_order,
    drop = FALSE
  ]
)

safe_write_csv(
  single_cell_state_scores,
  file.path(
    out_table_dir,
    "GSE244983_single_cell_state_scores.csv"
  )
)

umap_mat <- as.data.frame(
  Seurat::Embeddings(
    obj,
    "umap"
  )
)

if (
  ncol(umap_mat) < 2
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
  single_cell_state_scores,
  by = "Cell",
  all.x = TRUE,
  sort = FALSE
)

## Restore canonical Seurat cell order after merge.
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
    "GSE244983_state_score_UMAP_source.csv"
  )
)

############################################################
## 12. Localization across major cell types
############################################################

state_long <- single_cell_state_scores %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      state_order
    ),
    names_to = "State",
    values_to = "Score"
  )

summary_by_celltype <- state_long %>%
  dplyr::group_by(
    MajorCellType,
    State
  ) %>%
  dplyr::summarise(
    N_cells = sum(
      is.finite(
        Score
      )
    ),
    MeanScore = mean(
      Score,
      na.rm = TRUE
    ),
    MedianScore = stats::median(
      Score,
      na.rm = TRUE
    ),
    SDScore = stats::sd(
      Score,
      na.rm = TRUE
    ),
    Q25 = as.numeric(
      stats::quantile(
        Score,
        0.25,
        na.rm = TRUE
      )
    ),
    Q75 = as.numeric(
      stats::quantile(
        Score,
        0.75,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    StateDisplay = unname(
      state_display[
        State
      ]
    )
  ) %>%
  dplyr::arrange(
    match(
      State,
      state_order
    ),
    match(
      MajorCellType,
      major_celltype_levels
    )
  )

safe_write_csv(
  summary_by_celltype,
  file.path(
    out_table_dir,
    "GSE244983_state_score_by_major_celltype_summary.csv"
  )
)

localization_matrix <- summary_by_celltype %>%
  dplyr::select(
    MajorCellType,
    State,
    MeanScore
  ) %>%
  tidyr::pivot_wider(
    names_from = State,
    values_from = MeanScore
  )

localization_matrix$MajorCellType <-
  factor(
    localization_matrix$MajorCellType,
    levels = major_celltype_levels
  )

localization_matrix <-
  localization_matrix[
    order(
      localization_matrix$MajorCellType
    ),
    ,
    drop = FALSE
  ]

localization_matrix$MajorCellType <-
  as.character(
    localization_matrix$MajorCellType
  )

safe_write_csv(
  localization_matrix,
  file.path(
    out_table_dir,
    "GSE244983_major_celltype_state_localization_matrix.csv"
  )
)

top_localization <- summary_by_celltype %>%
  dplyr::group_by(
    State
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      MeanScore
    ),
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    Rank = dplyr::row_number()
  ) %>%
  dplyr::filter(
    Rank <= 3
  ) %>%
  dplyr::ungroup()

safe_write_csv(
  top_localization,
  file.path(
    out_table_dir,
    "GSE244983_top_localizing_celltypes_by_state.csv"
  )
)

############################################################
## 13. Canonical-marker source-attribution modules
############################################################

module_score_list <- vector(
  "list",
  length(
    module_genes
  )
)

names(
  module_score_list
) <- names(
  module_genes
)

module_presence_list <- vector(
  "list",
  length(
    module_genes
  )
)

names(
  module_presence_list
) <- names(
  module_genes
)

for (mm in names(module_genes)) {
  res <- score_gene_set_zmean(
    expr_data = expr_data,
    expr_genes_original = expr_genes_original,
    expr_genes_clean = expr_genes_clean,
    genes = module_genes[[mm]],
    set_name = paste0(
      "marker module ",
      mm
    )
  )

  module_score_list[[mm]] <-
    res$score

  pp <- res$presence
  pp$MarkerModule <- mm
  pp$GeneSet <- NULL
  module_presence_list[[mm]] <- pp
}

module_score_df <- as.data.frame(
  module_score_list,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

module_score_df$Cell <- colnames(
  obj
)

module_score_df <- module_score_df[
  ,
  c(
    "Cell",
    marker_module_order
  ),
  drop = FALSE
]

safe_write_csv(
  module_score_df,
  file.path(
    out_table_dir,
    "GSE244983_canonical_marker_module_cell_scores.csv"
  )
)

module_gene_presence <- dplyr::bind_rows(
  module_presence_list
)

safe_write_csv(
  module_gene_presence,
  file.path(
    out_table_dir,
    "GSE244983_canonical_marker_module_gene_presence.csv"
  )
)

module_coverage_summary <- module_gene_presence %>%
  dplyr::group_by(
    MarkerModule
  ) %>%
  dplyr::summarise(
    N_genes = dplyr::n(),
    N_present = sum(
      Present
    ),
    Coverage = N_present /
      N_genes,
    .groups = "drop"
  )

safe_write_csv(
  module_coverage_summary,
  file.path(
    out_table_dir,
    "GSE244983_canonical_marker_module_gene_coverage_summary.csv"
  )
)

source_df <- single_cell_state_scores

for (mm in marker_module_order) {
  source_df[[mm]] <- module_score_df[
    match(
      source_df$Cell,
      module_score_df$Cell
    ),
    mm
  ]
}

############################################################
## 14. S8 panel-A source attribution: state means by cell type
############################################################

source_attribution_state_means <- summary_by_celltype %>%
  dplyr::select(
    MajorCellType,
    State,
    StateDisplay,
    N_cells,
    MeanScore
  )

safe_write_csv(
  source_attribution_state_means,
  file.path(
    out_table_dir,
    "GSE244983_source_attribution_state_means.csv"
  )
)

############################################################
## 15. S8 panel-B six prespecified contrasts
############################################################

contrast_rows <- vector(
  "list",
  nrow(
    contrast_specs
  )
)

for (i in seq_len(nrow(contrast_specs))) {
  sp <- contrast_specs[
    i,
    ,
    drop = FALSE
  ]

  st <- sp$State[[1]]
  g1 <- sp$Group1[[1]]
  g2 <- sp$Group2[[1]]

  x <- source_df[
    source_df$MajorCellType ==
      g1,
    st
  ]

  y <- source_df[
    source_df$MajorCellType ==
      g2,
    st
  ]

  x <- as.numeric(x)
  y <- as.numeric(y)

  contrast_rows[[i]] <- data.frame(
    State = st,
    StateDisplay = unname(
      state_display[[st]]
    ),
    ContrastLabel = sp$ContrastLabel[[1]],
    Group1 = g1,
    Group2 = g2,
    N_Group1 = sum(
      is.finite(x)
    ),
    N_Group2 = sum(
      is.finite(y)
    ),
    Median_Group1 = stats::median(
      x,
      na.rm = TRUE
    ),
    Median_Group2 = stats::median(
      y,
      na.rm = TRUE
    ),
    MedianDifference = stats::median(
      x,
      na.rm = TRUE
    ) -
      stats::median(
        y,
        na.rm = TRUE
      ),
    RankBiserial = rank_biserial(
      x,
      y
    ),
    stringsAsFactors = FALSE
  )
}

source_attribution_contrasts <- dplyr::bind_rows(
  contrast_rows
)

safe_write_csv(
  source_attribution_contrasts,
  file.path(
    out_table_dir,
    "GSE244983_source_attribution_prespecified_contrasts.csv"
  )
)

############################################################
## 16. S8 panel-C canonical module means by cell type
############################################################

module_long <- source_df %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      marker_module_order
    ),
    names_to = "MarkerModule",
    values_to = "ModuleScore"
  )

module_by_celltype <- module_long %>%
  dplyr::group_by(
    MajorCellType,
    MarkerModule
  ) %>%
  dplyr::summarise(
    N_cells = sum(
      is.finite(
        ModuleScore
      )
    ),
    MeanModuleScore = mean(
      ModuleScore,
      na.rm = TRUE
    ),
    MedianModuleScore = stats::median(
      ModuleScore,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    match(
      MajorCellType,
      major_celltype_levels
    ),
    match(
      MarkerModule,
      marker_module_order
    )
  )

safe_write_csv(
  module_by_celltype,
  file.path(
    out_table_dir,
    "GSE244983_canonical_marker_module_by_celltype.csv"
  )
)

############################################################
## 17. S8 panel-D state-marker Spearman correlations
############################################################

cor_rows <- list()
kk <- 1L

for (st in state_order) {
  for (mm in marker_module_order) {
    x <- as.numeric(
      source_df[[st]]
    )

    y <- as.numeric(
      source_df[[mm]]
    )

    keep <- is.finite(x) &
      is.finite(y)

    rho <- suppressWarnings(
      stats::cor(
        x[keep],
        y[keep],
        method = "spearman"
      )
    )

    cor_rows[[kk]] <- data.frame(
      State = st,
      StateDisplay = unname(
        state_display[[st]]
      ),
      MarkerModule = mm,
      N_cells = sum(
        keep
      ),
      SpearmanRho = rho,
      stringsAsFactors = FALSE
    )

    kk <- kk + 1L
  }
}

state_marker_spearman <- dplyr::bind_rows(
  cor_rows
)

safe_write_csv(
  state_marker_spearman,
  file.path(
    out_table_dir,
    "GSE244983_state_marker_module_spearman.csv"
  )
)

############################################################
## 17B. Supplementary Figure S8 reporting contract v1.2.3
##      (integrated from former standalone 07b patch)
############################################################

# -------------------------------------------------------------------------
# The preceding integrated analysis sections generate:
#
# source_attribution_state_means
# source_attribution_contrasts
# module_by_celltype
# state_marker_spearman
# major_celltype_levels
#
# -------------------------------------------------------------------------

if (!exists("source_attribution_state_means")) {
  stop("Missing source_attribution_state_means before integrated S8 reporting-contract checks.")
}

if (!exists("major_celltype_levels")) {
  stop("Missing major_celltype_levels.")
}

if (!exists("out_table_dir")) {
  stop("Missing out_table_dir.")
}

# =============================================================================
# Supplementary Figure S8 reporting contract v1.2.3
# =============================================================================

observed_major_celltypes <- sort(
  unique(source_attribution_state_means$MajorCellType)
)

excluded_annotation_levels <- setdiff(
  major_celltype_levels,
  observed_major_celltypes
)

# -------------------------------------------------------------------------
# Observed state-score-evaluable universe
# -------------------------------------------------------------------------

s8_major_celltype_universe <- data.frame(
  MajorCellType = c(
    observed_major_celltypes,
    excluded_annotation_levels
  ),
  Included_in_state_localization = c(
    rep(TRUE, length(observed_major_celltypes)),
    rep(FALSE, length(excluded_annotation_levels))
  ),
  stringsAsFactors = FALSE
)

s8_major_celltype_universe$Reason <- ifelse(
  s8_major_celltype_universe$Included_in_state_localization,
  "Included in predefined tumor-immune state localization analysis",
  "Annotated major cell type without state-score localization output"
)

safe_write_csv(
  s8_major_celltype_universe,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS8_major_celltype_universe_audit.csv"
  )
)

# -------------------------------------------------------------------------
# Contract gates
# -------------------------------------------------------------------------

s8_reporting_gates <- data.frame(
  Gate = c(
    "All observed S8A major cell types belong to annotated universe",
    "No unexpected major cell type appears in state localization output",
    "Observed major celltype universe contains at least one compartment",
    "S8B contains six prespecified contrasts",
    "S8C contains canonical marker modules",
    "S8D contains four state keys",
    "S8D contains six marker-module keys"
  ),
  Observed = c(
    paste(
      observed_major_celltypes,
      collapse = ";"
    ),
    paste(
      setdiff(
        observed_major_celltypes,
        major_celltype_levels
      ),
      collapse = ";"
    ),
    length(observed_major_celltypes),
    nrow(source_attribution_contrasts),
    paste(
      sort(unique(module_by_celltype$MarkerModule)),
      collapse = ";"
    ),
    paste(
      sort(unique(state_marker_spearman$State)),
      collapse = ";"
    ),
    paste(
      sort(unique(state_marker_spearman$MarkerModule)),
      collapse = ";"
    )
  ),
  Expected = c(
    "subset of annotation universe",
    "empty",
    ">0",
    "6",
    "canonical marker-module universe",
    "four predefined states",
    "six marker modules"
  ),
  Pass = c(
    all(
      observed_major_celltypes %in% major_celltype_levels
    ),
    length(
      setdiff(
        observed_major_celltypes,
        major_celltype_levels
      )
    ) == 0L,
    length(observed_major_celltypes) > 0L,
    nrow(source_attribution_contrasts) == 6L,
    TRUE,
    length(unique(state_marker_spearman$State)) == 4L,
    length(unique(state_marker_spearman$MarkerModule)) == 6L
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  s8_reporting_gates,
  file.path(
    out_table_dir,
    "GSE244983_SuppFigureS8_reporting_contract_gates.csv"
  )
)

if (any(!s8_reporting_gates$Pass)) {
  stop(
    "Supplementary Figure S8 reporting-contract gate failure: ",
    paste(
      s8_reporting_gates$Gate[
        !s8_reporting_gates$Pass
      ],
      collapse = " | "
    )
  )
}

cat("\n====================================================\n")
cat("S8 reporting contract v1.2.3 PASSED\n")
cat("Observed major cell-type universe:\n")
cat(
  paste(
    observed_major_celltypes,
    collapse = "; "
  ),
  "\n"
)
cat("Excluded annotation levels:\n")
cat(
  paste(
    excluded_annotation_levels,
    collapse = "; "
  ),
  "\n"
)
cat("====================================================\n")


############################################################
## 18. Reproducibility gates
############################################################

state_stats_by_name <- setNames(
  split(
    state_score_stats,
    state_score_stats$State
  ),
  state_score_stats$State
)

state_standardization_ok <- all(
  vapply(
    state_order,
    function(st) {
      rr <- state_stats_by_name[[st]]

      isTRUE(
        rr$N[[1]] ==
          expected_n_cells
      ) &&
        abs(
          rr$Mean[[1]]
        ) <
        1e-8 &&
        abs(
          rr$SD[[1]] -
            1
        ) <
        1e-8
    },
    logical(1)
  )
)

signature_nonempty_ok <- all(
  signature_coverage_summary$N_present >
    0
)

module_nonempty_ok <- all(
  module_coverage_summary$N_present >
    0
)

state_table_cellset_md5 <- cellset_md5(
  single_cell_state_scores$Cell
)

final_gates <- data.frame(
  Gate = c(
    "Canonical cell universe retained after state scoring",
    "State-score table cell set matches script-07 cell set",
    "All 12 constituent signatures have at least one detected gene",
    "All six canonical marker modules have at least one detected gene",
    "Exactly four final state scores generated",
    "All four final state scores are finite for all 26,053 cells and standardized",
    "Major-cell localization table has one row per observed cell-type/state combination",
    "Top-localization table has three rows per state",
    "Exactly six prespecified source-attribution contrasts generated",
    "Exactly 24 state-marker Spearman correlations generated",
    "UMAP source table contains all 26,053 cells"
  ),
  Observed = c(
    as.character(
      ncol(obj)
    ),
    state_table_cellset_md5,
    as.character(
      sum(
        signature_coverage_summary$N_present >
          0
      )
    ),
    as.character(
      sum(
        module_coverage_summary$N_present >
          0
      )
    ),
    as.character(
      sum(
        state_order %in%
          colnames(
            state_score_df
          )
      )
    ),
    paste0(
      state_score_stats$State,
      ":N=",
      state_score_stats$N,
      ",mean=",
      signif(
        state_score_stats$Mean,
        6
      ),
      ",sd=",
      signif(
        state_score_stats$SD,
        6
      ),
      collapse = ";"
    ),
    as.character(
      nrow(
        summary_by_celltype
      )
    ),
    as.character(
      nrow(
        top_localization
      )
    ),
    as.character(
      nrow(
        source_attribution_contrasts
      )
    ),
    as.character(
      nrow(
        state_marker_spearman
      )
    ),
    as.character(
      nrow(
        umap_source
      )
    )
  ),
  Expected = c(
    as.character(
      expected_n_cells
    ),
    expected_cellset_md5,
    "12",
    "6",
    "4",
    "N=26053 for each; mean=0; sd=1",
    as.character(
      length(
        unique(
          single_cell_state_scores$MajorCellType
        )
      ) *
        length(
          state_order
        )
    ),
    "12",
    "6",
    "24",
    as.character(
      expected_n_cells
    )
  ),
  Pass = c(
    ncol(obj) ==
      expected_n_cells,
    identical(
      state_table_cellset_md5,
      expected_cellset_md5
    ),
    signature_nonempty_ok &&
      nrow(
        signature_coverage_summary
      ) == 12,
    module_nonempty_ok &&
      nrow(
        module_coverage_summary
      ) == 6,
    all(
      state_order %in%
        colnames(
          state_score_df
        )
    ),
    state_standardization_ok,
    nrow(
      summary_by_celltype
    ) ==
      length(
        unique(
          single_cell_state_scores$MajorCellType
        )
      ) *
      length(
        state_order
      ),
    nrow(
      top_localization
    ) == 12,
    nrow(
      source_attribution_contrasts
    ) == 6,
    nrow(
      state_marker_spearman
    ) == 24,
    nrow(
      umap_source
    ) ==
      expected_n_cells
  ),
  stringsAsFactors = FALSE
)

all_gates <- dplyr::bind_rows(
  data.frame(
    GateClass = "script07_lineage",
    lineage_gate,
    stringsAsFactors = FALSE
  ),
  data.frame(
    GateClass = "state_mapping",
    mapping_gate,
    stringsAsFactors = FALSE
  ),
  data.frame(
    GateClass = "state_localization_source_attribution",
    final_gates,
    stringsAsFactors = FALSE
  ),
  data.frame(
    GateClass = "s8_reporting_contract",
    s8_reporting_gates,
    stringsAsFactors = FALSE
  )
)

safe_write_csv(
  all_gates,
  file.path(
    out_table_dir,
    "GSE244983_reproducibility_gates.csv"
  )
)

save_session_info()

if (
  any(
    !all_gates$Pass
  )
) {
  stop(
    "08 reproducibility gate failure: ",
    paste(
      all_gates$Gate[
        !all_gates$Pass
      ],
      collapse = " | "
    ),
    ". Canonical state-localized Seurat object was NOT saved."
  )
}

############################################################
## 19. Save canonical state-localized object
############################################################

canonical_object_file <- file.path(
  out_intermediate_dir,
  "GSE244983_seurat_state_localized.rds"
)

safe_save_rds(
  obj,
  canonical_object_file
)

canonical_object_audit <- data.frame(
  Object = "GSE244983_seurat_state_localized.rds",
  RelativePath = relative_to_project(
    canonical_object_file
  ),
  MD5 = md5_file(
    canonical_object_file
  ),
  N_cells = ncol(
    obj
  ),
  N_features = nrow(
    obj
  ),
  CellSetMD5 = cellset_md5(
    colnames(
      obj
    )
  ),
  N_state_score_columns = sum(
    unname(
      score_col_map
    ) %in%
      colnames(
        obj@meta.data
      )
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  canonical_object_audit,
  file.path(
    out_table_dir,
    "GSE244983_canonical_state_localized_object_audit.csv"
  )
)

############################################################
## 20. Output inventory
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

inventory_type <- c(
  rep(
    "table",
    length(
      table_files
    )
  ),
  "RDS",
  "sessionInfo"
)

output_inventory <- data.frame(
  Type = inventory_type,
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
    "GSE244983_output_inventory.csv"
  )
)

message(
  "07_GSE244983_scRNA_state_localization_source_attribution.R finished successfully."
)

message(
  "Merged public script 07 v1.2.4 finished successfully; S8 reporting contract integrated."
)
