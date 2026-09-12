############################################################
## 10_CellChat_LR_compare_with_frozen.R
##
## SELF-CONTAINED BLOCKING REPRODUCTION AUDIT
##
## Frozen numerical authorities are bundled under:
##   frozen_reference/
##
## They were extracted from the final historical Supplementary
## Tables S24-S26 and are SHA-tracked by:
##   frozen_reference/FROZEN_REFERENCE_MANIFEST.csv
##
## This script does NOT rerun CellChat or CosMx.
############################################################

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR")

if (!nzchar(PROJECT_DIR)) {
  PROJECT_DIR <- "D:/ICB_resistance_project"
}

PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)

SCRIPT_DIR <- file.path(
  PROJECT_DIR,
  "scripts",
  "analysis",
  "14_CellChat_LR"
)

REF_DIR <- file.path(
  SCRIPT_DIR,
  "frozen_reference"
)

LOCKED_INPUT_DIR <- file.path(
  SCRIPT_DIR,
  "locked_input"
)

CUR16 <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "CellChat_candidate_LR"
)

CUR17 <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "CosMx_spatial_LR_validation"
)

REPORT_TABLE <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "CellChat_LR_reporting"
)

REPORT_FIG <- file.path(
  PROJECT_DIR,
  "results",
  "figures",
  "CellChat_LR_reporting"
)

AUDIT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "audit",
  "14_CellChat_LR_sequential_reproduction_check"
)

dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
## 1. Helpers
############################################################

read_csv_required <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing required file: ",
      path,
      call. = FALSE
    )
  }

  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  if (requireNamespace("digest", quietly = TRUE)) {
    return(
      digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }

  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, "rb")
    on.exit(
      try(close(con), silent = TRUE),
      add = TRUE
    )

    return(
      paste0(
        as.character(
          openssl::sha256(con)
        ),
        collapse = ""
      )
    )
  }

  stop(
    "Package 'digest' or 'openssl' is required for SHA256 verification.",
    call. = FALSE
  )
}

safe_num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

logical_robust <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  y <- tolower(
    trimws(
      as.character(x)
    )
  )

  y %in% c(
    "true",
    "t",
    "1",
    "1.0",
    "yes",
    "y"
  )
}

normalize_chr <- function(x) {
  y <- as.character(x)

  y <- gsub(
    "\u2013|\u2014",
    "-",
    y
  )

  y <- gsub(
    "[[:space:]]+",
    " ",
    y
  )

  trimws(y)
}

normalize_chr_empty_na <- function(x) {
  y <- normalize_chr(x)

  y[
    is.na(y) |
      !nzchar(y) |
      tolower(y) %in% c(
        "na",
        "nan"
      )
  ] <- NA_character_

  y
}

character_exact_empty_na <- function(
  a,
  b
) {
  aa <- normalize_chr_empty_na(a)
  bb <- normalize_chr_empty_na(b)

  (
    is.na(aa) &
      is.na(bb)
  ) |
    (
      !is.na(aa) &
        !is.na(bb) &
        aa == bb
    )
}

numeric_exact <- function(
  a,
  b,
  tol = 1e-12
) {
  aa <- safe_num(a)
  bb <- safe_num(b)

  (
    is.na(aa) &
    is.na(bb)
  ) |
    (
      is.finite(aa) &
      is.finite(bb) &
      abs(
        aa -
        bb
      ) <=
        tol
    )
}

character_exact <- function(
  a,
  b
) {
  aa <- normalize_chr(a)
  bb <- normalize_chr(b)

  (
    is.na(aa) &
    is.na(bb)
  ) |
    (
      !is.na(aa) &
      !is.na(bb) &
      aa ==
        bb
    )
}

max_numeric_diff <- function(
  a,
  b
) {
  aa <- safe_num(a)
  bb <- safe_num(b)

  d <- abs(
    aa -
      bb
  )

  if (
    length(d) == 0L ||
    all(
      is.na(d)
    )
  ) {
    return(0)
  }

  max(
    d,
    na.rm = TRUE
  )
}

align_by_pair <- function(
  frozen,
  current
) {
  if (
    !setequal(
      as.character(
        frozen$pair_family
      ),
      as.character(
        current$pair_family
      )
    )
  ) {
    stop(
      "pair_family universe differs between current and frozen tables.",
      call. = FALSE
    )
  }

  current[
    match(
      as.character(
        frozen$pair_family
      ),
      as.character(
        current$pair_family
      )
    ),
    ,
    drop = FALSE
  ]
}

make_k_pair_contrast_key <- function(
  d
) {
  paste(
    as.integer(
      safe_num(
        d$k
      )
    ),
    normalize_chr(
      d$pair_family
    ),
    normalize_chr(
      d$contrast
    ),
    sep = "||"
  )
}

gate <- function(
  id,
  section,
  description,
  pass,
  observed,
  expected
) {
  data.frame(
    gate_id = id,
    section = section,
    description = description,
    observed = as.character(
      observed
    ),
    expected = as.character(
      expected
    ),
    pass = isTRUE(
      pass
    ),
    stringsAsFactors = FALSE
  )
}

############################################################
## 2. Read bundled frozen authorities
############################################################

FROZEN_S24 <- file.path(
  REF_DIR,
  "S24_CellChat_candidate_summary_FROZEN.csv"
)

FROZEN_S25 <- file.path(
  REF_DIR,
  "S25_CosMx_K20_summary_FROZEN.csv"
)

FROZEN_S26 <- file.path(
  REF_DIR,
  "S26_CosMx_all_kNN_180_tests_FROZEN.csv"
)

FROZEN_MANIFEST <- file.path(
  REF_DIR,
  "FROZEN_REFERENCE_MANIFEST.csv"
)

f24 <- read_csv_required(
  FROZEN_S24
)

f25 <- read_csv_required(
  FROZEN_S25
)

f26 <- read_csv_required(
  FROZEN_S26
)

manifest <- read_csv_required(
  FROZEN_MANIFEST
)

frozen_reference_paths <- c(
  S24_CellChat_candidate_summary_FROZEN.csv =
    FROZEN_S24,
  S25_CosMx_K20_summary_FROZEN.csv =
    FROZEN_S25,
  S26_CosMx_all_kNN_180_tests_FROZEN.csv =
    FROZEN_S26
)

frozen_reference_sha_detail <- data.frame()

for (nm in names(frozen_reference_paths)) {
  row_manifest <- manifest[
    as.character(
      manifest$derived_reference
    ) == nm,
    ,
    drop = FALSE
  ]

  expected_sha <- if (
    nrow(row_manifest) == 1L
  ) {
    as.character(
      row_manifest$derived_sha256[1]
    )
  } else {
    NA_character_
  }

  observed_sha <- sha256_file(
    frozen_reference_paths[[nm]]
  )

  frozen_reference_sha_detail <- rbind(
    frozen_reference_sha_detail,
    data.frame(
      reference = nm,
      expected_sha256 = expected_sha,
      observed_sha256 = observed_sha,
      pass =
        !is.na(expected_sha) &&
        identical(
          expected_sha,
          observed_sha
        ),
      stringsAsFactors = FALSE
    )
  )
}

utils::write.csv(
  frozen_reference_sha_detail,
  file.path(
    AUDIT_DIR,
    "00_FROZEN_REFERENCE_SHA256_CHECK.csv"
  ),
  row.names = FALSE
)

frozen_contract_pass <- (
  nrow(f24) == 12L &&
  nrow(f25) == 12L &&
  nrow(f26) == 180L &&
  nrow(manifest) == 3L &&
  all(
    frozen_reference_sha_detail$pass
  )
)

LOCKED_SELECTION_FILE <- file.path(
  LOCKED_INPUT_DIR,
  "CellChat_historical_sampling_group_manifest.csv"
)

LOCKED_SELECTION_MANIFEST_FILE <- file.path(
  LOCKED_INPUT_DIR,
  "LOCKED_INPUT_MANIFEST.csv"
)

locked_selection <- read_csv_required(
  LOCKED_SELECTION_FILE
)

locked_selection_manifest <- read_csv_required(
  LOCKED_SELECTION_MANIFEST_FILE
)

locked_selection_sha <- sha256_file(
  LOCKED_SELECTION_FILE
)

expected_locked_selection_sha <-
  "b4e1dd0cd5f3da4808d9284e23de7491dd4779186021530971c7dad43f579620"

locked_selection_groups <- sort(
  unique(
    as.character(
      locked_selection$cellchat_group
    )
  )
)

EXPECTED_CELLCHAT_GROUPS <- sort(
  c(
    "Malignant",
    "Cycling malignant",
    "CAF/stromal-like cells",
    "Endothelial",
    "Myeloid cells",
    "T/NK cells",
    "T/NK/Treg-like cells",
    "B/plasma cells"
  )
)

locked_input_contract_pass <- (
  nrow(locked_selection) == 3219L &&
  length(
    unique(
      as.character(
        locked_selection$cell_id
      )
    )
  ) == 3219L &&
  setequal(
    locked_selection_groups,
    EXPECTED_CELLCHAT_GROUPS
  ) &&
  identical(
    locked_selection_sha,
    expected_locked_selection_sha
  ) &&
  nrow(locked_selection_manifest) == 1L &&
  identical(
    as.character(
      locked_selection_manifest$sha256[1]
    ),
    expected_locked_selection_sha
  )
)

############################################################
## 3. Current Step01 outputs
############################################################

cur_summary <- read_csv_required(
  file.path(
    CUR16,
    "CellChat_candidate_LR_summary.csv"
  )
)

cur_pairs <- read_csv_required(
  file.path(
    CUR16,
    "CellChat_candidate_pairs_for_CosMx_validation.csv"
  )
)

cur_nom <- read_csv_required(
  file.path(
    CUR16,
    "CellChat_candidate_LR_nomination.csv"
  )
)

cur_input_summary_path <- file.path(
  CUR16,
  "CellChat_input_summary.csv"
)

cur_input_summary <- if (
  file.exists(
    cur_input_summary_path
  )
) {
  read_csv_required(
    cur_input_summary_path
  )
} else {
  NULL
}

############################################################
## 4. Step01 frozen pair-axis contract
############################################################

cur_pairs_aligned <- align_by_pair(
  f24,
  cur_pairs
)

pair_family_exact <- all(
  character_exact(
    f24$pair_family,
    cur_pairs_aligned$pair_family
  )
)

ligand_pattern_exact <- all(
  character_exact(
    f24$ligand_pattern,
    cur_pairs_aligned$ligand_pattern
  )
)

receptor_pattern_exact <- all(
  character_exact(
    f24$receptor_pattern,
    cur_pairs_aligned$receptor_pattern
  )
)

pair_nomination_exact <- all(
  logical_robust(
    f24$nominated_by_cellchat
  ) ==
  logical_robust(
    cur_pairs_aligned$nominated_by_cellchat
  )
)

axis_contract_pass <- (
  nrow(cur_pairs_aligned) == 12L &&
  pair_family_exact &&
  ligand_pattern_exact &&
  receptor_pattern_exact &&
  pair_nomination_exact
)

############################################################
## 5. Step01 frozen CellChat summary comparison
############################################################

cur_summary_aligned <- align_by_pair(
  f24,
  cur_summary
)

summary_fields_numeric <- c(
  "n_interactions",
  "n_source_target_pairs",
  "max_prob",
  "median_prob",
  "min_p_value"
)

summary_fields_character <- c(
  "top_source_target",
  "direction_classes"
)

s24_detail <- data.frame()

for (nm in summary_fields_numeric) {
  exact <- numeric_exact(
    f24[[nm]],
    cur_summary_aligned[[nm]]
  )

  row <- data.frame(
    field = nm,
    exact_rows = sum(exact),
    total_rows = length(exact),
    max_numeric_abs_diff =
      max_numeric_diff(
        f24[[nm]],
        cur_summary_aligned[[nm]]
      ),
    pass = all(exact),
    stringsAsFactors = FALSE
  )

  s24_detail <- rbind(
    s24_detail,
    row
  )
}

for (nm in summary_fields_character) {
  exact <- character_exact_empty_na(
    f24[[nm]],
    cur_summary_aligned[[nm]]
  )

  row <- data.frame(
    field = nm,
    exact_rows = sum(exact),
    total_rows = length(exact),
    max_numeric_abs_diff = NA_real_,
    pass = all(exact),
    stringsAsFactors = FALSE
  )

  s24_detail <- rbind(
    s24_detail,
    row
  )
}

nom_exact <- (
  logical_robust(
    f24$nominated_by_cellchat
  ) ==
  logical_robust(
    cur_summary_aligned$nominated_by_cellchat
  )
)

s24_detail <- rbind(
  s24_detail,
  data.frame(
    field = "nominated_by_cellchat",
    exact_rows = sum(nom_exact),
    total_rows = length(nom_exact),
    max_numeric_abs_diff = NA_real_,
    pass = all(nom_exact),
    stringsAsFactors = FALSE
  )
)

utils::write.csv(
  s24_detail,
  file.path(
    AUDIT_DIR,
    "01_STEP01_CURRENT_VS_FROZEN_S24_DETAIL.csv"
  ),
  row.names = FALSE
)

s24_numeric_exact_pass <- all(
  s24_detail$pass
)

############################################################
## 6. Step01 current nomination -> current summary rebuild
############################################################

EXPECTED_AXES <- as.character(
  f24$pair_family
)

rebuilt_rows <- lapply(
  EXPECTED_AXES,
  function(pf) {
    dd <- cur_nom[
      as.character(
        cur_nom$pair_family
      ) ==
        pf,
      ,
      drop = FALSE
    ]

    if (nrow(dd) == 0L) {
      return(
        data.frame(
          pair_family = pf,
          n_interactions = 0L,
          n_source_target_pairs = 0L,
          max_prob = NA_real_,
          median_prob = NA_real_,
          min_p_value = NA_real_,
          top_source_target = NA_character_,
          direction_classes = NA_character_,
          nominated_by_cellchat = FALSE,
          stringsAsFactors = FALSE
        )
      )
    }

    prob <- safe_num(
      dd$prob
    )

    pval <- safe_num(
      dd$pval
    )

    idx_max <- which.max(
      prob
    )

    data.frame(
      pair_family = pf,
      n_interactions = nrow(dd),
      n_source_target_pairs =
        length(
          unique(
            as.character(
              dd$source_target
            )
          )
        ),
      max_prob =
        max(
          prob,
          na.rm = TRUE
        ),
      median_prob =
        stats::median(
          prob,
          na.rm = TRUE
        ),
      min_p_value =
        min(
          pval,
          na.rm = TRUE
        ),
      top_source_target =
        as.character(
          dd$source_target[
            idx_max
          ]
        ),
      direction_classes =
        paste(
          sort(
            unique(
              as.character(
                dd$direction_class
              )
            )
          ),
          collapse = ";"
        ),
      nominated_by_cellchat =
        any(
          pval <
            0.05,
          na.rm = TRUE
        ),
      stringsAsFactors = FALSE
    )
  }
)

rebuilt <- do.call(
  rbind,
  rebuilt_rows
)

rebuilt_aligned <- align_by_pair(
  cur_summary_aligned,
  rebuilt
)

current_rebuild_pass <- TRUE

for (nm in summary_fields_numeric) {
  current_rebuild_pass <- (
    current_rebuild_pass &&
    all(
      numeric_exact(
        cur_summary_aligned[[nm]],
        rebuilt_aligned[[nm]]
      )
    )
  )
}

for (nm in summary_fields_character) {
  current_rebuild_pass <- (
    current_rebuild_pass &&
    all(
      character_exact_empty_na(
        cur_summary_aligned[[nm]],
        rebuilt_aligned[[nm]]
      )
    )
  )
}

current_rebuild_pass <- (
  current_rebuild_pass &&
  all(
    logical_robust(
      cur_summary_aligned$nominated_by_cellchat
    ) ==
    logical_robust(
      rebuilt_aligned$nominated_by_cellchat
    )
  )
)

############################################################
## 7. Step02-04 frozen S25 comparison
############################################################

cur_k20 <- read_csv_required(
  file.path(
    REPORT_TABLE,
    "CosMx_LR_K20_summary.csv"
  )
)

cur_k20_aligned <- align_by_pair(
  f25,
  cur_k20
)

s25_numeric <- c(
  "n_evaluable_contrasts",
  "n_enriched_contrasts",
  "max_log2_spatial_enrichment",
  "max_neighbor_fraction_percent",
  "min_FDR"
)

s25_character <- c(
  "ligand_gene",
  "receptor_genes",
  "validation_priority",
  "top_contrast",
  "overall_support"
)

s25_detail <- data.frame()

for (nm in s25_numeric) {
  exact <- numeric_exact(
    f25[[nm]],
    cur_k20_aligned[[nm]]
  )

  s25_detail <- rbind(
    s25_detail,
    data.frame(
      field = nm,
      exact_rows = sum(exact),
      total_rows = length(exact),
      max_numeric_abs_diff =
        max_numeric_diff(
          f25[[nm]],
          cur_k20_aligned[[nm]]
        ),
      pass = all(exact),
      stringsAsFactors = FALSE
    )
  )
}

for (nm in s25_character) {
  exact <- character_exact(
    f25[[nm]],
    cur_k20_aligned[[nm]]
  )

  s25_detail <- rbind(
    s25_detail,
    data.frame(
      field = nm,
      exact_rows = sum(exact),
      total_rows = length(exact),
      max_numeric_abs_diff = NA_real_,
      pass = all(exact),
      stringsAsFactors = FALSE
    )
  )
}

s25_nom_exact <- (
  logical_robust(
    f25$nominated_by_cellchat
  ) ==
  logical_robust(
    cur_k20_aligned$nominated_by_cellchat
  )
)

s25_detail <- rbind(
  s25_detail,
  data.frame(
    field = "nominated_by_cellchat",
    exact_rows = sum(s25_nom_exact),
    total_rows = length(s25_nom_exact),
    max_numeric_abs_diff = NA_real_,
    pass = all(s25_nom_exact),
    stringsAsFactors = FALSE
  )
)

utils::write.csv(
  s25_detail,
  file.path(
    AUDIT_DIR,
    "02_STEP02_04_CURRENT_VS_FROZEN_S25_DETAIL.csv"
  ),
  row.names = FALSE
)

s25_pass <- all(
  s25_detail$pass
)

############################################################
## 8. Step02-04 frozen S26 180-test comparison
############################################################

cur_all <- read_csv_required(
  file.path(
    REPORT_TABLE,
    "CosMx_LR_all_kNN.csv"
  )
)

if (
  nrow(cur_all) != 180L ||
  nrow(f26) != 180L
) {
  stop(
    "S26 test-universe row count is not 180/180.",
    call. = FALSE
  )
}

f26$key <- make_k_pair_contrast_key(
  f26
)

cur_all$key <- make_k_pair_contrast_key(
  cur_all
)

if (
  !setequal(
    f26$key,
    cur_all$key
  )
) {
  stop(
    "Current and frozen S26 key universes differ.",
    call. = FALSE
  )
}

cur_all_aligned <- cur_all[
  match(
    f26$key,
    cur_all$key
  ),
  ,
  drop = FALSE
]

s26_numeric <- c(
  "log2_spatial_enrichment",
  "FDR",
  "neighbor_fraction_percent",
  "expected_fraction_percent"
)

s26_character <- c(
  "validation_priority",
  "status",
  "interpretation"
)

s26_detail <- data.frame()

for (nm in s26_numeric) {
  exact <- numeric_exact(
    f26[[nm]],
    cur_all_aligned[[nm]]
  )

  s26_detail <- rbind(
    s26_detail,
    data.frame(
      field = nm,
      exact_rows = sum(exact),
      total_rows = length(exact),
      max_numeric_abs_diff =
        max_numeric_diff(
          f26[[nm]],
          cur_all_aligned[[nm]]
        ),
      pass = all(exact),
      stringsAsFactors = FALSE
    )
  )
}

for (nm in s26_character) {
  exact <- character_exact(
    f26[[nm]],
    cur_all_aligned[[nm]]
  )

  s26_detail <- rbind(
    s26_detail,
    data.frame(
      field = nm,
      exact_rows = sum(exact),
      total_rows = length(exact),
      max_numeric_abs_diff = NA_real_,
      pass = all(exact),
      stringsAsFactors = FALSE
    )
  )
}

utils::write.csv(
  s26_detail,
  file.path(
    AUDIT_DIR,
    "03_STEP02_04_CURRENT_VS_FROZEN_S26_DETAIL.csv"
  ),
  row.names = FALSE
)

s26_pass <- all(
  s26_detail$pass
)

############################################################
## 9. k-specific BH-FDR reconstruction
############################################################

fdr_detail <- data.frame()

for (k in c(10, 20, 30)) {
  d <- read_csv_required(
    file.path(
      CUR17,
      paste0(
        "KNN",
        k
      ),
      "CosMx_spatial_LR_validation_all.csv"
    )
  )

  if (
    nrow(d) != 60L
  ) {
    stop(
      "k=",
      k,
      " current table does not contain 60 tests.",
      call. = FALSE
    )
  }

  rec <- stats::p.adjust(
    safe_num(
      d$empirical_p
    ),
    method = "BH"
  )

  exact <- numeric_exact(
    rec,
    d$FDR
  )

  fdr_detail <- rbind(
    fdr_detail,
    data.frame(
      k = k,
      exact_rows = sum(exact),
      total_rows = length(exact),
      max_numeric_abs_diff =
        max_numeric_diff(
          rec,
          d$FDR
        ),
      pass = all(exact),
      stringsAsFactors = FALSE
    )
  )
}

utils::write.csv(
  fdr_detail,
  file.path(
    AUDIT_DIR,
    "04_K_SPECIFIC_BH_FDR_REPRODUCTION.csv"
  ),
  row.names = FALSE
)

fdr_pass <- all(
  fdr_detail$pass
)

############################################################
## 10. Reporting contracts / artifacts
############################################################

pair_by_k <- read_csv_required(
  file.path(
    REPORT_TABLE,
    "CosMx_LR_pair_summary_by_k.csv"
  )
)

robustness <- read_csv_required(
  file.path(
    REPORT_TABLE,
    "CosMx_LR_kNN_robustness.csv"
  )
)

consolidation_pass <- (
  nrow(cur_all) == 180L &&
  nrow(pair_by_k) == 36L &&
  nrow(robustness) == 12L
)

s24_files <- file.path(
  REPORT_FIG,
  c(
    "Supplementary_Figure_S24_CellChat_candidate_LR_axes.png",
    "Supplementary_Figure_S24_CellChat_candidate_LR_axes.jpg",
    "Supplementary_Figure_S24_CellChat_candidate_LR_axes.pdf",
    "Supplementary_Figure_S24_CellChat_candidate_LR_axes_audit.csv"
  )
)

s25_files <- file.path(
  REPORT_FIG,
  c(
    "Supplementary_Figure_S25_CosMx_spatial_proximity_candidate_LR_axes.png",
    "Supplementary_Figure_S25_CosMx_spatial_proximity_candidate_LR_axes.jpg",
    "Supplementary_Figure_S25_CosMx_spatial_proximity_candidate_LR_axes.pdf",
    "Supplementary_Figure_S25_CosMx_spatial_proximity_candidate_LR_axes_audit.csv",
    "Supplementary_Figure_S25_ranked_axis_contrast_source.csv"
  )
)

s26_files <- file.path(
  REPORT_FIG,
  c(
    "Supplementary_Figure_S26_CosMx_LR_kNN_sensitivity.png",
    "Supplementary_Figure_S26_CosMx_LR_kNN_sensitivity.jpg",
    "Supplementary_Figure_S26_CosMx_LR_kNN_sensitivity.pdf",
    "Supplementary_Figure_S26_CosMx_LR_kNN_sensitivity_audit.csv"
  )
)

table_files <- file.path(
  REPORT_TABLE,
  c(
    "Supplementary_Table_S24_CellChat_candidate_LR.xlsx",
    "Supplementary_Table_S25_CosMx_LR_K20.xlsx",
    "Supplementary_Table_S26_CosMx_LR_kNN_sensitivity.xlsx"
  )
)

s24_artifacts_pass <- all(
  file.exists(
    s24_files
  )
)

s25_artifacts_pass <- all(
  file.exists(
    s25_files
  )
)

s26_artifacts_pass <- all(
  file.exists(
    s26_files
  )
)

tables_pass <- all(
  file.exists(
    table_files
  )
)

############################################################
## 11. Input lineage record
############################################################

input_lineage <- data.frame(
  item = c(
    "current_input_summary",
    "selected_scRNA_object",
    "n_cells_after_filter_downsample",
    "groups_used",
    "selection_manifest",
    "selection_manifest_sha256",
    "selection_manifest_n",
    "manifest_cells_missing_from_public"
  ),
  value = c(
    cur_input_summary_path,
    if (
      !is.null(cur_input_summary) &&
      "selected_scRNA_object" %in%
        names(cur_input_summary)
    ) {
      as.character(
        cur_input_summary$selected_scRNA_object[1]
      )
    } else {
      NA_character_
    },
    if (
      !is.null(cur_input_summary) &&
      "n_cells_after_filter_downsample" %in%
        names(cur_input_summary)
    ) {
      as.character(
        cur_input_summary$n_cells_after_filter_downsample[1]
      )
    } else {
      NA_character_
    },
    if (
      !is.null(cur_input_summary) &&
      "groups_used" %in%
        names(cur_input_summary)
    ) {
      as.character(
        cur_input_summary$groups_used[1]
      )
    } else {
      NA_character_
    },
    if (
      !is.null(cur_input_summary) &&
      "selection_manifest" %in%
        names(cur_input_summary)
    ) {
      as.character(
        cur_input_summary$selection_manifest[1]
      )
    } else {
      NA_character_
    },
    if (
      !is.null(cur_input_summary) &&
      "selection_manifest_sha256" %in%
        names(cur_input_summary)
    ) {
      as.character(
        cur_input_summary$selection_manifest_sha256[1]
      )
    } else {
      NA_character_
    },
    if (
      !is.null(cur_input_summary) &&
      "selection_manifest_n" %in%
        names(cur_input_summary)
    ) {
      as.character(
        cur_input_summary$selection_manifest_n[1]
      )
    } else {
      NA_character_
    },
    if (
      !is.null(cur_input_summary) &&
      "manifest_cells_missing_from_public" %in%
        names(cur_input_summary)
    ) {
      as.character(
        cur_input_summary$manifest_cells_missing_from_public[1]
      )
    } else {
      NA_character_
    }
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  input_lineage,
  file.path(
    AUDIT_DIR,
    "05_CURRENT_STEP01_INPUT_LINEAGE.csv"
  ),
  row.names = FALSE
)

current_input_manifest_pass <- (
  !is.null(cur_input_summary) &&
  "selection_manifest_sha256" %in%
    names(cur_input_summary) &&
  "selection_manifest_n" %in%
    names(cur_input_summary) &&
  "manifest_cells_missing_from_public" %in%
    names(cur_input_summary) &&
  identical(
    as.character(
      cur_input_summary$selection_manifest_sha256[1]
    ),
    expected_locked_selection_sha
  ) &&
  safe_num(
    cur_input_summary$selection_manifest_n[1]
  ) == 3219 &&
  safe_num(
    cur_input_summary$manifest_cells_missing_from_public[1]
  ) == 0 &&
  safe_num(
    cur_input_summary$n_cells_after_filter_downsample[1]
  ) == 3219
)

############################################################
## 12. Blocking gates
############################################################

gates <- rbind(
  gate(
    "R01",
    "Frozen_reference",
    "Bundled frozen S24/S25/S26 authority contract is complete",
    frozen_contract_pass,
    paste0(
      nrow(f24),
      ";",
      nrow(f25),
      ";",
      nrow(f26)
    ),
    "12;12;180"
  ),
  gate(
    "R01B",
    "Frozen_reference",
    "Locked historical Step16 cell/group selection manifest is complete and SHA256-verified",
    locked_input_contract_pass,
    paste0(
      nrow(locked_selection),
      " cells; SHA=",
      locked_selection_sha
    ),
    "3219 cells; locked SHA256"
  ),
  gate(
    "R02",
    "Step01_CellChat",
    "Frozen 12-axis ligand/receptor pattern and nominated-status universe reproduces",
    axis_contract_pass,
    axis_contract_pass,
    TRUE
  ),
  gate(
    "R02B",
    "Step01_CellChat",
    "Current Step01 consumes the locked 3,219-cell historical selection manifest with zero missing public cells",
    current_input_manifest_pass,
    current_input_manifest_pass,
    TRUE
  ),
  gate(
    "R03",
    "Step01_CellChat",
    "Current CellChat 12-axis numerical summary exactly reproduces frozen S24",
    s24_numeric_exact_pass,
    paste0(
      sum(
        s24_detail$pass
      ),
      "/",
      nrow(
        s24_detail
      ),
      " fields exact"
    ),
    "8/8 fields exact"
  ),
  gate(
    "R04",
    "Step01_CellChat",
    "Current nomination rows exactly reconstruct current 12-axis summary",
    current_rebuild_pass,
    current_rebuild_pass,
    TRUE
  ),
  gate(
    "R05",
    "Step02_04_CosMx_LR",
    "Current K20 12-axis summary exactly reproduces frozen S25",
    s25_pass,
    paste0(
      sum(
        s25_detail$pass
      ),
      "/",
      nrow(
        s25_detail
      ),
      " fields exact"
    ),
    paste0(
      nrow(
        s25_detail
      ),
      "/",
      nrow(
        s25_detail
      )
    )
  ),
  gate(
    "R06",
    "Step02_04_CosMx_LR",
    "Current 180 contrast-level tests exactly reproduce frozen S26",
    s26_pass,
    paste0(
      sum(
        s26_detail$pass
      ),
      "/",
      nrow(
        s26_detail
      ),
      " fields exact"
    ),
    paste0(
      nrow(
        s26_detail
      ),
      "/",
      nrow(
        s26_detail
      )
    )
  ),
  gate(
    "R07",
    "Step02_04_CosMx_LR",
    "BH-FDR exactly recomputes within each k-specific 60-test universe",
    fdr_pass,
    paste0(
      sum(
        fdr_detail$pass
      ),
      "/3"
    ),
    "3/3"
  ),
  gate(
    "R08",
    "Step05_consolidation",
    "Reporting contract is 180 tests / 36 pair-by-k rows / 12 axes",
    consolidation_pass,
    ifelse(
      consolidation_pass,
      "180;36;12",
      "failed"
    ),
    "180;36;12"
  ),
  gate(
    "R09",
    "Step06_09_reporting",
    "S24 PNG/JPG/PDF/audit exist",
    s24_artifacts_pass,
    sum(
      file.exists(
        s24_files
      )
    ),
    4
  ),
  gate(
    "R10",
    "Step06_09_reporting",
    "S25 PNG/JPG/PDF/source/audit exist",
    s25_artifacts_pass,
    sum(
      file.exists(
        s25_files
      )
    ),
    5
  ),
  gate(
    "R11",
    "Step06_09_reporting",
    "S26 PNG/JPG/PDF/audit exist",
    s26_artifacts_pass,
    sum(
      file.exists(
        s26_files
      )
    ),
    4
  ),
  gate(
    "R12",
    "Step06_09_reporting",
    "Supplementary Tables S24-S26 exist",
    tables_pass,
    sum(
      file.exists(
        table_files
      )
    ),
    3
  )
)

utils::write.csv(
  gates,
  file.path(
    AUDIT_DIR,
    "06_REPRODUCTION_GATES.csv"
  ),
  row.names = FALSE
)

section_summary <- do.call(
  rbind,
  lapply(
    unique(
      gates$section
    ),
    function(sec) {
      dd <- gates[
        gates$section ==
          sec,
        ,
        drop = FALSE
      ]

      data.frame(
        section = sec,
        passed = sum(
          dd$pass
        ),
        total = nrow(
          dd
        ),
        status =
          ifelse(
            all(
              dd$pass
            ),
            "PASS",
            "FAIL"
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

utils::write.csv(
  section_summary,
  file.path(
    AUDIT_DIR,
    "07_REPRODUCTION_SUMMARY.csv"
  ),
  row.names = FALSE
)

blocking_failures <- sum(
  !gates$pass
)

only_cellchat_numeric_drift <- (
  blocking_failures == 1L &&
  !gates$pass[
    gates$gate_id ==
      "R03"
  ]
)

decision <- if (
  blocking_failures ==
    0L
) {
  "READY_FOR_PUBLIC_FREEZE"
} else if (
  isTRUE(
    only_cellchat_numeric_drift
  )
) {
  "BLOCKED_STEP01_CELLCHAT_NUMERICAL_DRIFT"
} else {
  "BLOCKED"
}

writeLines(
  c(
    "14_CellChat_LR self-contained reproduction audit",
    paste0(
      "Blocking failures: ",
      blocking_failures
    ),
    paste0(
      "Overall: ",
      ifelse(
        blocking_failures ==
          0L,
        "PASS",
        "FAIL"
      )
    ),
    paste0(
      "Decision: ",
      decision
    ),
    "",
    "Important interpretation:",
    "If R03 is the only failure, the frozen 12-axis candidate universe and all downstream CosMx K20/K10/K30 results reproduce, but the rerun CellChat communication counts/probabilities differ from the frozen S24 numerical summary.",
    "Do not change downstream CosMx logic to repair an isolated R03 failure."
  ),
  file.path(
    AUDIT_DIR,
    "DECISION.txt"
  )
)

cat("\n============================================================\n")
cat("14_CellChat_LR SELF-CONTAINED REPRODUCTION AUDIT\n")
cat("============================================================\n")

for (
  i in seq_len(
    nrow(
      section_summary
    )
  )
) {
  cat(
    section_summary$section[i],
    ": ",
    section_summary$passed[i],
    "/",
    section_summary$total[i],
    " ",
    section_summary$status[i],
    "\n",
    sep = ""
  )
}

cat(
  "Blocking failures: ",
  blocking_failures,
  "\n",
  sep = ""
)

cat(
  "Overall: ",
  ifelse(
    blocking_failures ==
      0L,
    "PASS",
    "FAIL"
  ),
  "\n",
  sep = ""
)

cat(
  "Decision: ",
  decision,
  "\n",
  sep = ""
)

cat(
  "Audit: ",
  AUDIT_DIR,
  "\n",
  sep = ""
)

cat("============================================================\n")

if (
  blocking_failures !=
    0L
) {
  stop(
    paste0(
      "14_CellChat_LR reproduction audit not fully closed: ",
      decision
    ),
    call. = FALSE
  )
}
