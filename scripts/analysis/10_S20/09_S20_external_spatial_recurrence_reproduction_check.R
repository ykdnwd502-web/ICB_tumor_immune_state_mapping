############################################################
## 10_S20_external_spatial_recurrence_reproduction_check.R
## CHECKER PATCH: v1.0.2 — robust Figure S20 audit comparison
##
## PURPOSE
## -------
## Reproduction audit ONLY for the complete 10_S20 external
## melanoma spatial recurrence module:
##
##   01 GSE250636 raw -> historical 18J
##   02 Thrane2018 raw -> historical 18J
##   03 historical 18J -> historical S28 reporting layer
##   04 GSE250636 -> source-locked canonical S20
##   05 Thrane2018 -> source-locked canonical S20
##   06 canonical S20 aggregate
##   07 Supplementary Figure S20
##
## This script:
##   - does NOT rerun any spatial analysis
##   - does NOT modify any result
##   - does NOT redefine state signatures
##   - verifies file lineage, state-gene locks, historical-to-
##     canonical numerical reproduction, aggregate consistency,
##     and final Figure-S20 reporting artifacts
##
## Interpretation boundary:
##   External spatial recurrence audit only.
##   NOT matched ICB-response validation.
##   NOT biomarker validation.
##   NOT functional signaling / niche validation.
############################################################

options(stringsAsFactors = FALSE)

############################################################
## 0. Project root
############################################################

PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR"
)

if (!nzchar(PROJECT_DIR)) {
  PROJECT_DIR <- "D:/ICB_resistance_project"
}

PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)

TOL <- 1e-12

############################################################
## 1. Canonical paths
############################################################

HIST_ROOT <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "revision_external_spatial_recurrence_18J"
)

HIST_GSE_DIR <- file.path(
  HIST_ROOT,
  "GSE250636"
)

HIST_THR_DIR <- file.path(
  HIST_ROOT,
  "Thrane2018_legacyST"
)

HIST_S28_DIR <- file.path(
  HIST_ROOT,
  "S28_external_spatial_recurrence_audits"
)

CANON_ROOT <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "S20_external_spatial_recurrence_CANONICAL_v1.0"
)

CANON_GSE_DIR <- file.path(
  CANON_ROOT,
  "GSE250636"
)

CANON_THR_DIR <- file.path(
  CANON_ROOT,
  "Thrane2018_legacyST"
)

COMBINED_DIR <- file.path(
  CANON_ROOT,
  "combined"
)

FIG_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "figures",
  "supplementary"
)

LOCK_FILE <- file.path(
  PROJECT_DIR,
  "public_release",
  "S20_external_spatial_recurrence",
  "config",
  "S20_locked_state_genes.csv"
)

AUDIT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "audit",
  "10_S20_external_spatial_recurrence_reproduction_check"
)

dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
## 2. Required package
############################################################

if (!requireNamespace(
  "data.table",
  quietly = TRUE
)) {
  stop(
    "Required package missing: data.table",
    call. = FALSE
  )
}

############################################################
## 3. Helpers
############################################################

read_csv0 <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Missing required file: ",
      path,
      call. = FALSE
    )
  }

  as.data.frame(
    data.table::fread(
      path,
      data.table = FALSE,
      check.names = FALSE
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
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

  y <- trimws(y)

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

safe_num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

standardize_state <- function(x) {
  xx <- tolower(
    as.character(x)
  )

  xx <- gsub(
    "[^a-z0-9]+",
    "_",
    xx
  )

  out <- rep(
    NA_character_,
    length(xx)
  )

  out[
    grepl(
      "immune",
      xx
    ) &
      (
        grepl(
          "cold",
          xx
        ) |
          grepl(
            "defective",
            xx
          )
      )
  ] <- "Immune_defective_Cold"

  out[
    grepl(
      "myeloid",
      xx
    ) &
      (
        grepl(
          "treg",
          xx
        ) |
          grepl(
            "immunosuppress",
            xx
          )
      )
  ] <- "Myeloid_Treg_Immunosuppressive"

  out[
    (
      grepl(
        "dediff",
        xx
      ) |
        grepl(
          "dedifferentiation",
          xx
        )
    ) &
      (
        grepl(
          "stromal",
          xx
        ) |
          grepl(
            "remodel",
            xx
          )
      )
  ] <- "Tumor_dedifferentiation_Stromal_remodeling"

  out[
    grepl(
      "melanocytic",
      xx
    ) |
      (
        grepl(
          "differentiation",
          xx
        ) &
          !grepl(
            "dediff|stromal|remodel",
            xx
          )
      )
  ] <- "Melanocytic_Differentiation"

  out
}

parse_gene_cell <- function(x) {
  x <- as.character(x)
  x <- x[
    !is.na(x)
  ]

  if (length(x) == 0L) {
    return(
      character(0)
    )
  }

  x <- paste(
    x,
    collapse = ";"
  )

  x <- gsub(
    "[\r\n\t]+",
    ";",
    x
  )

  x <- unlist(
    strsplit(
      x,
      "[,;|]+|\\s+"
    )
  )

  x <- toupper(
    trimws(x)
  )

  x <- sub(
    "\\..*$",
    "",
    x
  )

  x <- x[
    !is.na(x) &
      nzchar(x) &
      x != "NA"
  ]

  unique(x)
}

extract_state_gene <- function(path) {
  df <- read_csv0(
    path
  )

  nms <- colnames(
    df
  )

  state_candidates <- c(
    "state",
    "State",
    "FinalState",
    "final_state",
    "module",
    "Module"
  )

  gene_candidates <- c(
    "gene",
    "Gene",
    "GeneSymbol",
    "gene_symbol",
    "symbol",
    "Symbol",
    "retained_genes",
    "RetainedGenes",
    "genes",
    "Genes"
  )

  state_col <- state_candidates[
    state_candidates %in%
      nms
  ]

  gene_col <- gene_candidates[
    gene_candidates %in%
      nms
  ]

  rows <- list()
  kk <- 1L

  if (
    length(state_col) > 0L &&
      length(gene_col) > 0L
  ) {
    state_col <- state_col[[1]]
    gene_col <- gene_col[[1]]

    for (i in seq_len(nrow(df))) {
      st <- standardize_state(
        df[[state_col]][[i]]
      )

      genes <- parse_gene_cell(
        df[[gene_col]][[i]]
      )

      if (
        !is.na(st) &&
          length(genes) > 0L
      ) {
        rows[[kk]] <- data.frame(
          state = st,
          gene = genes,
          stringsAsFactors = FALSE
        )

        kk <- kk + 1L
      }
    }
  } else {
    for (cc in nms) {
      st <- standardize_state(
        cc
      )

      if (is.na(st)) {
        next
      }

      genes <- parse_gene_cell(
        df[[cc]]
      )

      if (length(genes) > 0L) {
        rows[[kk]] <- data.frame(
          state = st,
          gene = genes,
          stringsAsFactors = FALSE
        )

        kk <- kk + 1L
      }
    }
  }

  if (length(rows) == 0L) {
    stop(
      "Could not parse four-state gene lock from: ",
      path,
      call. = FALSE
    )
  }

  out <- unique(
    do.call(
      rbind,
      rows
    )
  )

  out <- out[
    order(
      out$state,
      out$gene
    ),
    ,
    drop = FALSE
  ]

  rownames(
    out
  ) <- NULL

  out
}

state_gene_equal <- function(
  path_a,
  path_b
) {
  a <- extract_state_gene(
    path_a
  )

  b <- extract_state_gene(
    path_b
  )

  identical(
    a,
    b
  )
}

file_inventory <- function(
  paths,
  labels
) {
  data.frame(
    Item = labels,
    Path = paths,
    Exists = file.exists(paths),
    stringsAsFactors = FALSE
  )
}

compare_scientific_csv <- function(
  historical_file,
  canonical_file,
  key_candidates = character(0),
  ignore_regex = "(?i)(^.*path$|^.*dir$|selected.*file|gene.*set.*file)"
) {
  result <- data.frame(
    Historical = historical_file,
    Canonical = canonical_file,
    Rows_historical = NA_integer_,
    Rows_canonical = NA_integer_,
    Common_columns = NA_integer_,
    Key_columns = NA_character_,
    Max_numeric_abs_diff = NA_real_,
    Character_exact = FALSE,
    Numeric_exact = FALSE,
    Pass = FALSE,
    stringsAsFactors = FALSE
  )

  if (
    !file.exists(historical_file) ||
      !file.exists(canonical_file)
  ) {
    return(result)
  }

  a <- read_csv0(
    historical_file
  )

  b <- read_csv0(
    canonical_file
  )

  result$Rows_historical <- nrow(a)
  result$Rows_canonical <- nrow(b)

  common_cols <- intersect(
    colnames(a),
    colnames(b)
  )

  if (
    length(common_cols) == 0L ||
      nrow(a) != nrow(b)
  ) {
    return(result)
  }

  ignore_cols <- common_cols[
    grepl(
      ignore_regex,
      common_cols,
      perl = TRUE
    )
  ]

  compare_cols <- setdiff(
    common_cols,
    ignore_cols
  )

  if (length(compare_cols) == 0L) {
    return(result)
  }

  result$Common_columns <- length(
    compare_cols
  )

  keys <- intersect(
    key_candidates,
    compare_cols
  )

  if (length(keys) == 0L) {
    keys <- compare_cols[
      vapply(
        a[
          compare_cols
        ],
        function(x) {
          is.character(x) ||
            is.factor(x) ||
            is.logical(x)
        },
        logical(1)
      )
    ]
  }

  make_key <- function(df, key_cols) {
    if (length(key_cols) == 0L) {
      return(
        as.character(
          seq_len(
            nrow(df)
          )
        )
      )
    }

    vv <- lapply(
      df[
        key_cols
      ],
      function(x) {
        y <- normalize_chr(x)
        y[
          is.na(y)
        ] <- "<NA>"
        y
      }
    )

    do.call(
      paste,
      c(
        vv,
        sep = "\u241f"
      )
    )
  }

  key_a <- make_key(
    a,
    keys
  )

  key_b <- make_key(
    b,
    keys
  )

  result$Key_columns <- paste(
    keys,
    collapse = ";"
  )

  ## If keys are unique, align directly.
  ## If they are not, retain file row order. Historical and canonical
  ## producers use the same scientific aggregation order.
  if (
    length(keys) > 0L &&
      !anyDuplicated(key_a) &&
      !anyDuplicated(key_b)
  ) {
    if (!setequal(key_a, key_b)) {
      return(result)
    }

    b <- b[
      match(
        key_a,
        key_b
      ),
      ,
      drop = FALSE
    ]
  }

  numeric_cols <- compare_cols[
    vapply(
      compare_cols,
      function(cc) {
        is.numeric(a[[cc]]) &&
          is.numeric(b[[cc]])
      },
      logical(1)
    )
  ]

  character_cols <- setdiff(
    compare_cols,
    numeric_cols
  )

  numeric_ok <- TRUE
  max_diff <- 0

  for (cc in numeric_cols) {
    aa <- safe_num(
      a[[cc]]
    )

    bb <- safe_num(
      b[[cc]]
    )

    both_na <- is.na(aa) &
      is.na(bb)

    both_finite <- is.finite(aa) &
      is.finite(bb)

    exact <- both_na |
      (
        both_finite &
          abs(
            aa - bb
          ) <= TOL
      )

    if (!all(exact)) {
      numeric_ok <- FALSE
    }

    dd <- abs(
      aa - bb
    )

    dd <- dd[
      is.finite(dd)
    ]

    if (length(dd) > 0L) {
      max_diff <- max(
        max_diff,
        max(dd)
      )
    }
  }

  character_ok <- TRUE

  for (cc in character_cols) {
    aa <- normalize_chr(
      a[[cc]]
    )

    bb <- normalize_chr(
      b[[cc]]
    )

    exact <- (
      is.na(aa) &
        is.na(bb)
    ) |
      (
        !is.na(aa) &
          !is.na(bb) &
          aa == bb
      )

    if (!all(exact)) {
      character_ok <- FALSE
    }
  }

  result$Max_numeric_abs_diff <- max_diff
  result$Character_exact <- character_ok
  result$Numeric_exact <- numeric_ok
  result$Pass <- numeric_ok &&
    character_ok

  result
}

compare_subset_to_parent <- function(
  parent_file,
  subset_file,
  readout_value
) {
  if (
    !file.exists(parent_file) ||
      !file.exists(subset_file)
  ) {
    return(FALSE)
  }

  parent <- read_csv0(
    parent_file
  )

  subset <- read_csv0(
    subset_file
  )

  if (!"readout" %in% colnames(parent)) {
    return(FALSE)
  }

  expected <- parent[
    as.character(
      parent$readout
    ) ==
      readout_value,
    ,
    drop = FALSE
  ]

  common_cols <- intersect(
    colnames(expected),
    colnames(subset)
  )

  if (length(common_cols) == 0L) {
    return(FALSE)
  }

  ## compare a stable scientific row serialization
  canon_rows <- function(df, cols) {
    df <- df[
      ,
      cols,
      drop = FALSE
    ]

    vals <- lapply(
      df,
      function(x) {
        if (is.numeric(x)) {
          out <- rep(
            "<NA>",
            length(x)
          )

          ok <- is.finite(x)

          out[ok] <- sprintf(
            "%.12g",
            x[ok]
          )

          out
        } else {
          out <- normalize_chr(x)

          out[
            is.na(out)
          ] <- "<NA>"

          out
        }
      }
    )

    sort(
      do.call(
        paste,
        c(
          vals,
          sep = "\u241f"
        )
      )
    )
  }

  identical(
    canon_rows(
      expected,
      common_cols
    ),
    canon_rows(
      subset,
      common_cols
    )
  )
}

############################################################
## 4. Historical 18J output inventory
############################################################

hist_gse_required <- c(
  "18J_GSE250636_state_genes_used_long.csv",
  "18J_GSE250636_gene_coverage.csv",
  "18J_GSE250636_gene_coverage_summary.csv",
  "18J_GSE250636_sample_input_log.csv",
  "18J_GSE250636_key_readout.csv",
  "18J_GSE250636_primary_pair_dual_high_OR.csv",
  "18J_GSE250636_primary_pair_Moran_descriptive.csv",
  "18J_GSE250636_primary_pair_raw_vs_residualized_correlation.csv",
  "18J_GSE250636_spot_scores_with_residualized_primary_pair.csv",
  "18J_GSE250636_spot_state_composition_scores.csv"
)

hist_thr_required <- c(
  "18J_Thrane2018_state_genes_used_long.csv",
  "18J_Thrane2018_gene_coverage.csv",
  "18J_Thrane2018_gene_coverage_summary.csv",
  "18J_Thrane2018_sample_input_log.csv",
  "18J_Thrane2018_key_readout.csv",
  "18J_Thrane2018_primary_pair_dual_high_OR.csv",
  "18J_Thrane2018_primary_pair_Moran_descriptive.csv",
  "18J_Thrane2018_primary_pair_raw_vs_residualized_correlation.csv",
  "18J_Thrane2018_spot_scores_with_residualized_primary_pair.csv",
  "18J_Thrane2018_spot_state_composition_scores.csv"
)

hist_inventory <- rbind(
  file_inventory(
    file.path(
      HIST_GSE_DIR,
      hist_gse_required
    ),
    paste0(
      "GSE250636/",
      hist_gse_required
    )
  ),
  file_inventory(
    file.path(
      HIST_THR_DIR,
      hist_thr_required
    ),
    paste0(
      "Thrane2018/",
      hist_thr_required
    )
  )
)

write.csv(
  hist_inventory,
  file.path(
    AUDIT_DIR,
    "10_S20_historical_18J_output_inventory.csv"
  ),
  row.names = FALSE
)

############################################################
## 5. Historical S28 reporting-layer inventory
############################################################

s28_required <- c(
  "S28_plot_input_overview.csv",
  "S28_dataset_overview.csv",
  "S28_combined_gene_coverage_summary.csv",
  "S28_primary_and_composition_gene_coverage.csv",
  "S28_plot_input_coverage.csv",
  "S28_plot_input_spearman.csv",
  "S28_primary_pair_focused_summary.csv",
  "S28_combined_key_readout.csv",
  "S28_plot_input_dual_high_OR.csv",
  "S28_plot_input_moran.csv"
)

s28_inventory <- file_inventory(
  file.path(
    HIST_S28_DIR,
    s28_required
  ),
  s28_required
)

write.csv(
  s28_inventory,
  file.path(
    AUDIT_DIR,
    "10_S20_historical_S28_output_inventory.csv"
  ),
  row.names = FALSE
)

############################################################
## 6. Canonical S20 dataset output inventory
############################################################

canon_gse_required <- c(
  "S20_GSE250636_state_genes_used_long.csv",
  "S20_GSE250636_gene_coverage.csv",
  "S20_GSE250636_gene_coverage_summary.csv",
  "S20_GSE250636_sample_input_log.csv",
  "S20_GSE250636_key_readout.csv",
  "S20_GSE250636_primary_pair_dual_high_OR.csv",
  "S20_GSE250636_primary_pair_Moran_descriptive.csv",
  "S20_GSE250636_primary_pair_raw_vs_residualized_correlation.csv",
  "S20_GSE250636_spot_scores_with_residualized_primary_pair.csv",
  "S20_GSE250636_spot_state_composition_scores.csv"
)

canon_thr_required <- c(
  "S20_Thrane2018_state_genes_used_long.csv",
  "S20_Thrane2018_gene_coverage.csv",
  "S20_Thrane2018_gene_coverage_summary.csv",
  "S20_Thrane2018_sample_input_log.csv",
  "S20_Thrane2018_key_readout.csv",
  "S20_Thrane2018_primary_pair_dual_high_OR.csv",
  "S20_Thrane2018_primary_pair_Moran_descriptive.csv",
  "S20_Thrane2018_primary_pair_raw_vs_residualized_correlation.csv",
  "S20_Thrane2018_spot_scores_with_residualized_primary_pair.csv",
  "S20_Thrane2018_spot_state_composition_scores.csv",
  "S20_Thrane2018_inclusion_decision_helper.csv"
)

canon_inventory <- rbind(
  file_inventory(
    file.path(
      CANON_GSE_DIR,
      canon_gse_required
    ),
    paste0(
      "GSE250636/",
      canon_gse_required
    )
  ),
  file_inventory(
    file.path(
      CANON_THR_DIR,
      canon_thr_required
    ),
    paste0(
      "Thrane2018/",
      canon_thr_required
    )
  )
)

write.csv(
  canon_inventory,
  file.path(
    AUDIT_DIR,
    "10_S20_canonical_dataset_output_inventory.csv"
  ),
  row.names = FALSE
)

############################################################
## 7. State-gene lock audit
############################################################

state_files <- c(
  public_lock = LOCK_FILE,
  historical_GSE250636 = file.path(
    HIST_GSE_DIR,
    "18J_GSE250636_state_genes_used_long.csv"
  ),
  historical_Thrane2018 = file.path(
    HIST_THR_DIR,
    "18J_Thrane2018_state_genes_used_long.csv"
  ),
  canonical_GSE250636 = file.path(
    CANON_GSE_DIR,
    "S20_GSE250636_state_genes_used_long.csv"
  ),
  canonical_Thrane2018 = file.path(
    CANON_THR_DIR,
    "S20_Thrane2018_state_genes_used_long.csv"
  )
)

state_exists <- file.exists(
  state_files
)

state_lock_pass <- all(
  state_exists
)

reference_state_file <- state_files[
  "historical_GSE250636"
]

state_lock_rows <- list()

for (nm in names(state_files)) {
  pp <- state_files[[nm]]

  identical_to_reference <- FALSE
  n_rows <- NA_integer_
  n_states <- NA_integer_

  if (
    file.exists(pp) &&
      file.exists(reference_state_file)
  ) {
    parsed <- tryCatch(
      extract_state_gene(
        pp
      ),
      error = function(e) {
        NULL
      }
    )

    if (!is.null(parsed)) {
      n_rows <- nrow(
        parsed
      )

      n_states <- length(
        unique(
          parsed$state
        )
      )

      identical_to_reference <- tryCatch(
        state_gene_equal(
          reference_state_file,
          pp
        ),
        error = function(e) {
          FALSE
        }
      )
    }
  }

  state_lock_rows[[nm]] <- data.frame(
    Source = nm,
    Path = pp,
    Exists = file.exists(pp),
    N_state_gene_rows = n_rows,
    N_states = n_states,
    Identical_to_historical_GSE250636 = identical_to_reference,
    stringsAsFactors = FALSE
  )
}

state_lock_audit <- do.call(
  rbind,
  state_lock_rows
)

rownames(
  state_lock_audit
) <- NULL

state_lock_pass <- all(
  state_lock_audit$Exists
) &&
  all(
    state_lock_audit$Identical_to_historical_GSE250636
  ) &&
  all(
    state_lock_audit$N_states == 4L
  )

write.csv(
  state_lock_audit,
  file.path(
    AUDIT_DIR,
    "10_S20_state_gene_lock_audit.csv"
  ),
  row.names = FALSE
)

############################################################
## 8. Historical 18J -> canonical S20 numerical comparison
############################################################

comparison_specs <- list(
  GSE250636_key_readout = list(
    hist = file.path(
      HIST_GSE_DIR,
      "18J_GSE250636_key_readout.csv"
    ),
    canon = file.path(
      CANON_GSE_DIR,
      "S20_GSE250636_key_readout.csv"
    ),
    keys = c(
      "readout",
      "dataset",
      "sample_id",
      "adjustment",
      "statistic_label"
    )
  ),
  GSE250636_gene_coverage_summary = list(
    hist = file.path(
      HIST_GSE_DIR,
      "18J_GSE250636_gene_coverage_summary.csv"
    ),
    canon = file.path(
      CANON_GSE_DIR,
      "S20_GSE250636_gene_coverage_summary.csv"
    ),
    keys = c(
      "module"
    )
  ),
  GSE250636_primary_pair_correlation = list(
    hist = file.path(
      HIST_GSE_DIR,
      "18J_GSE250636_primary_pair_raw_vs_residualized_correlation.csv"
    ),
    canon = file.path(
      CANON_GSE_DIR,
      "S20_GSE250636_primary_pair_raw_vs_residualized_correlation.csv"
    ),
    keys = c(
      "dataset",
      "sample_id",
      "adjustment"
    )
  ),
  GSE250636_dual_high_OR = list(
    hist = file.path(
      HIST_GSE_DIR,
      "18J_GSE250636_primary_pair_dual_high_OR.csv"
    ),
    canon = file.path(
      CANON_GSE_DIR,
      "S20_GSE250636_primary_pair_dual_high_OR.csv"
    ),
    keys = c(
      "dataset",
      "sample_id",
      "adjustment"
    )
  ),
  GSE250636_Moran = list(
    hist = file.path(
      HIST_GSE_DIR,
      "18J_GSE250636_primary_pair_Moran_descriptive.csv"
    ),
    canon = file.path(
      CANON_GSE_DIR,
      "S20_GSE250636_primary_pair_Moran_descriptive.csv"
    ),
    keys = c(
      "dataset",
      "sample_id",
      "adjustment"
    )
  ),
  Thrane2018_key_readout = list(
    hist = file.path(
      HIST_THR_DIR,
      "18J_Thrane2018_key_readout.csv"
    ),
    canon = file.path(
      CANON_THR_DIR,
      "S20_Thrane2018_key_readout.csv"
    ),
    keys = c(
      "readout",
      "dataset",
      "sample_id",
      "adjustment",
      "statistic_label"
    )
  ),
  Thrane2018_gene_coverage_summary = list(
    hist = file.path(
      HIST_THR_DIR,
      "18J_Thrane2018_gene_coverage_summary.csv"
    ),
    canon = file.path(
      CANON_THR_DIR,
      "S20_Thrane2018_gene_coverage_summary.csv"
    ),
    keys = c(
      "module"
    )
  ),
  Thrane2018_primary_pair_correlation = list(
    hist = file.path(
      HIST_THR_DIR,
      "18J_Thrane2018_primary_pair_raw_vs_residualized_correlation.csv"
    ),
    canon = file.path(
      CANON_THR_DIR,
      "S20_Thrane2018_primary_pair_raw_vs_residualized_correlation.csv"
    ),
    keys = c(
      "dataset",
      "sample_id",
      "adjustment"
    )
  ),
  Thrane2018_dual_high_OR = list(
    hist = file.path(
      HIST_THR_DIR,
      "18J_Thrane2018_primary_pair_dual_high_OR.csv"
    ),
    canon = file.path(
      CANON_THR_DIR,
      "S20_Thrane2018_primary_pair_dual_high_OR.csv"
    ),
    keys = c(
      "dataset",
      "sample_id",
      "adjustment"
    )
  ),
  Thrane2018_Moran = list(
    hist = file.path(
      HIST_THR_DIR,
      "18J_Thrane2018_primary_pair_Moran_descriptive.csv"
    ),
    canon = file.path(
      CANON_THR_DIR,
      "S20_Thrane2018_primary_pair_Moran_descriptive.csv"
    ),
    keys = c(
      "dataset",
      "sample_id",
      "adjustment"
    )
  )
)

comparison_rows <- list()

for (nm in names(comparison_specs)) {
  sp <- comparison_specs[[nm]]

  rr <- compare_scientific_csv(
    historical_file = sp$hist,
    canonical_file = sp$canon,
    key_candidates = sp$keys
  )

  rr$Comparison <- nm

  comparison_rows[[nm]] <- rr
}

numeric_comparison <- do.call(
  rbind,
  comparison_rows
)

numeric_comparison <- numeric_comparison[
  ,
  c(
    "Comparison",
    setdiff(
      colnames(
        numeric_comparison
      ),
      "Comparison"
    )
  ),
  drop = FALSE
]

rownames(
  numeric_comparison
) <- NULL

historical_to_canonical_pass <- all(
  numeric_comparison$Pass
)

write.csv(
  numeric_comparison,
  file.path(
    AUDIT_DIR,
    "10_S20_historical_to_canonical_numeric_comparison.csv"
  ),
  row.names = FALSE
)

############################################################
## 9. Canonical aggregate inventory
############################################################

combined_required <- c(
  "S20_dataset_overview.csv",
  "S20_plot_input_overview.csv",
  "S20_combined_gene_coverage_summary.csv",
  "S20_primary_and_composition_gene_coverage.csv",
  "S20_plot_input_coverage.csv",
  "S20_combined_key_readout.csv",
  "S20_primary_pair_focused_summary.csv",
  "S20_plot_input_spearman.csv",
  "S20_plot_input_dual_high_OR.csv",
  "S20_plot_input_moran.csv"
)

combined_inventory <- file_inventory(
  file.path(
    COMBINED_DIR,
    combined_required
  ),
  combined_required
)

write.csv(
  combined_inventory,
  file.path(
    AUDIT_DIR,
    "10_S20_combined_output_inventory.csv"
  ),
  row.names = FALSE
)

############################################################
## 10. Aggregate consistency checks
############################################################

aggregate_gates <- data.frame(
  Gate = c(
    "Combined key readout row count equals two canonical datasets",
    "Combined coverage row count equals two canonical datasets",
    "Dataset overview contains exactly two external datasets",
    "Spearman plot input is exact readout subset",
    "Dual-high OR plot input is exact readout subset",
    "Moran plot input is exact readout subset",
    "Overview plot input row count matches dataset overview"
  ),
  Pass = FALSE,
  stringsAsFactors = FALSE
)

if (all(combined_inventory$Exists)) {
  gse_key <- read_csv0(
    file.path(
      CANON_GSE_DIR,
      "S20_GSE250636_key_readout.csv"
    )
  )

  thr_key <- read_csv0(
    file.path(
      CANON_THR_DIR,
      "S20_Thrane2018_key_readout.csv"
    )
  )

  gse_cov <- read_csv0(
    file.path(
      CANON_GSE_DIR,
      "S20_GSE250636_gene_coverage_summary.csv"
    )
  )

  thr_cov <- read_csv0(
    file.path(
      CANON_THR_DIR,
      "S20_Thrane2018_gene_coverage_summary.csv"
    )
  )

  combined_key_file <- file.path(
    COMBINED_DIR,
    "S20_combined_key_readout.csv"
  )

  combined_cov_file <- file.path(
    COMBINED_DIR,
    "S20_combined_gene_coverage_summary.csv"
  )

  overview_file <- file.path(
    COMBINED_DIR,
    "S20_dataset_overview.csv"
  )

  overview_plot_file <- file.path(
    COMBINED_DIR,
    "S20_plot_input_overview.csv"
  )

  combined_key <- read_csv0(
    combined_key_file
  )

  combined_cov <- read_csv0(
    combined_cov_file
  )

  overview <- read_csv0(
    overview_file
  )

  overview_plot <- read_csv0(
    overview_plot_file
  )

  aggregate_gates$Pass[[1]] <-
    nrow(combined_key) ==
      nrow(gse_key) +
        nrow(thr_key)

  aggregate_gates$Pass[[2]] <-
    nrow(combined_cov) ==
      nrow(gse_cov) +
        nrow(thr_cov)

  aggregate_gates$Pass[[3]] <-
    nrow(overview) == 2L &&
      "dataset_id" %in%
        colnames(
          overview
        ) &&
      setequal(
        as.character(
          overview$dataset_id
        ),
        c(
          "GSE250636",
          "Thrane2018_legacyST"
        )
      )

  aggregate_gates$Pass[[4]] <-
    compare_subset_to_parent(
      combined_key_file,
      file.path(
        COMBINED_DIR,
        "S20_plot_input_spearman.csv"
      ),
      "primary_pair_spearman"
    )

  aggregate_gates$Pass[[5]] <-
    compare_subset_to_parent(
      combined_key_file,
      file.path(
        COMBINED_DIR,
        "S20_plot_input_dual_high_OR.csv"
      ),
      "primary_pair_dual_high_OR"
    )

  aggregate_gates$Pass[[6]] <-
    compare_subset_to_parent(
      combined_key_file,
      file.path(
        COMBINED_DIR,
        "S20_plot_input_moran.csv"
      ),
      "primary_pair_bivariate_Moran_descriptive"
    )

  aggregate_gates$Pass[[7]] <-
    nrow(overview_plot) ==
      nrow(overview)
}

write.csv(
  aggregate_gates,
  file.path(
    AUDIT_DIR,
    "10_S20_aggregate_consistency_gates.csv"
  ),
  row.names = FALSE
)

############################################################
## 11. Final Supplementary Figure S20 artifact inventory
############################################################

figure_required <- c(
  "Supplementary_Figure_S20_external_spatial_recurrence.png",
  "Supplementary_Figure_S20_external_spatial_recurrence.jpg",
  "Supplementary_Figure_S20_external_spatial_recurrence.pdf",
  "Supplementary_Figure_S20_external_spatial_recurrence_audit.csv",
  "Supplementary_Figure_S20_panelA_external_dataset_counts.csv",
  "Supplementary_Figure_S20_panelB_signature_gene_coverage.csv",
  "Supplementary_Figure_S20_panelC_spearman_correlation.csv",
  "Supplementary_Figure_S20_panelD_dual_high_enrichment.csv",
  "Supplementary_Figure_S20_panelE_bivariate_moran.csv"
)

figure_inventory <- file_inventory(
  file.path(
    FIG_DIR,
    figure_required
  ),
  figure_required
)

write.csv(
  figure_inventory,
  file.path(
    AUDIT_DIR,
    "10_S20_figure_output_inventory.csv"
  ),
  row.names = FALSE
)

############################################################
## 12. Figure audit consistency
############################################################

figure_audit_gate <- FALSE
panel_row_counts_ok <- FALSE
panel_labels_ok <- FALSE

figure_audit_file <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S20_external_spatial_recurrence_audit.csv"
)

figure_audit_detail <- data.frame(
  Item = character(),
  Observed = character(),
  Expected = character(),
  Pass = logical(),
  stringsAsFactors = FALSE
)

if (
  all(
    figure_inventory$Exists
  ) &&
    file.exists(
      figure_audit_file
    )
) {
  fig_audit <- read_csv0(
    figure_audit_file
  )

  panel_files <- c(
    A = file.path(
      FIG_DIR,
      "Supplementary_Figure_S20_panelA_external_dataset_counts.csv"
    ),
    B = file.path(
      FIG_DIR,
      "Supplementary_Figure_S20_panelB_signature_gene_coverage.csv"
    ),
    C = file.path(
      FIG_DIR,
      "Supplementary_Figure_S20_panelC_spearman_correlation.csv"
    ),
    D = file.path(
      FIG_DIR,
      "Supplementary_Figure_S20_panelD_dual_high_enrichment.csv"
    ),
    E = file.path(
      FIG_DIR,
      "Supplementary_Figure_S20_panelE_bivariate_moran.csv"
    )
  )

  actual_panel_rows <- vapply(
    unname(
      panel_files
    ),
    function(pp) {
      nrow(
        read_csv0(
          pp
        )
      )
    },
    integer(1)
  )

  audit_row_cols <- c(
    "n_panelA_rows",
    "n_panelB_rows",
    "n_panelC_rows",
    "n_panelD_rows",
    "n_panelE_rows"
  )

  required_audit_cols <- c(
    audit_row_cols,
    "panel_labels"
  )

  audit_schema_ok <-
    nrow(fig_audit) >= 1L &&
    all(
      required_audit_cols %in%
        colnames(
          fig_audit
        )
    )

  if (audit_schema_ok) {

    expected_panel_rows <- suppressWarnings(
      as.integer(
        unlist(
          fig_audit[
            1,
            audit_row_cols,
            drop = FALSE
          ],
          use.names = FALSE
        )
      )
    )

    panel_row_counts_ok <-
      length(actual_panel_rows) == 5L &&
      length(expected_panel_rows) == 5L &&
      all(
        is.finite(
          actual_panel_rows
        )
      ) &&
      all(
        is.finite(
          expected_panel_rows
        )
      ) &&
      all(
        as.integer(
          actual_panel_rows
        ) ==
          as.integer(
            expected_panel_rows
          )
      )

    observed_labels <- trimws(
      unlist(
        strsplit(
          as.character(
            fig_audit$panel_labels[[1]]
          ),
          ";",
          fixed = TRUE
        )
      )
    )

    panel_labels_ok <- identical(
      observed_labels,
      c(
        "A",
        "B",
        "C",
        "D",
        "E"
      )
    )

    figure_audit_detail <- rbind(
      data.frame(
        Item = paste0(
          "Panel ",
          names(panel_files),
          " row count"
        ),
        Observed = as.character(
          actual_panel_rows
        ),
        Expected = as.character(
          expected_panel_rows
        ),
        Pass = as.integer(
          actual_panel_rows
        ) ==
          as.integer(
            expected_panel_rows
          ),
        stringsAsFactors = FALSE
      ),
      data.frame(
        Item = "Panel labels",
        Observed = paste(
          observed_labels,
          collapse = ";"
        ),
        Expected = "A;B;C;D;E",
        Pass = panel_labels_ok,
        stringsAsFactors = FALSE
      )
    )

    figure_audit_gate <-
      panel_row_counts_ok &&
      panel_labels_ok
  } else {
    figure_audit_detail <- data.frame(
      Item = "Figure audit schema",
      Observed = paste(
        colnames(
          fig_audit
        ),
        collapse = ";"
      ),
      Expected = paste(
        required_audit_cols,
        collapse = ";"
      ),
      Pass = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

write.csv(
  figure_audit_detail,
  file.path(
    AUDIT_DIR,
    "10_S20_figure_audit_consistency_detail.csv"
  ),
  row.names = FALSE
)

cat(
  "Figure panel row-count match: ",
  panel_row_counts_ok,
  "\n",
  sep = ""
)

cat(
  "Figure panel-label match: ",
  panel_labels_ok,
  "\n",
  sep = ""
)

############################################################
## 13. Final blocking summary
############################################################

final_summary <- data.frame(
  Gate = c(
    "Historical 18J outputs",
    "Historical S28 reporting layer",
    "Canonical S20 dataset outputs",
    "Four-state gene lock",
    "Historical-to-canonical numerical reproduction",
    "Canonical aggregate outputs",
    "Canonical aggregate consistency",
    "Supplementary Figure S20 outputs",
    "Supplementary Figure S20 audit consistency"
  ),
  Observed = c(
    paste0(
      sum(
        hist_inventory$Exists
      ),
      "/",
      nrow(
        hist_inventory
      )
    ),
    paste0(
      sum(
        s28_inventory$Exists
      ),
      "/",
      nrow(
        s28_inventory
      )
    ),
    paste0(
      sum(
        canon_inventory$Exists
      ),
      "/",
      nrow(
        canon_inventory
      )
    ),
    paste0(
      sum(
        state_lock_audit$Identical_to_historical_GSE250636
      ),
      "/",
      nrow(
        state_lock_audit
      ),
      " sources identical"
    ),
    paste0(
      sum(
        numeric_comparison$Pass
      ),
      "/",
      nrow(
        numeric_comparison
      ),
      " comparisons PASS; max diff=",
      format(
        max(
          numeric_comparison$Max_numeric_abs_diff,
          na.rm = TRUE
        ),
        scientific = TRUE
      )
    ),
    paste0(
      sum(
        combined_inventory$Exists
      ),
      "/",
      nrow(
        combined_inventory
      )
    ),
    paste0(
      sum(
        aggregate_gates$Pass
      ),
      "/",
      nrow(
        aggregate_gates
      )
    ),
    paste0(
      sum(
        figure_inventory$Exists
      ),
      "/",
      nrow(
        figure_inventory
      )
    ),
    as.character(
      figure_audit_gate
    )
  ),
  Pass = c(
    all(
      hist_inventory$Exists
    ),
    all(
      s28_inventory$Exists
    ),
    all(
      canon_inventory$Exists
    ),
    state_lock_pass,
    historical_to_canonical_pass,
    all(
      combined_inventory$Exists
    ),
    all(
      aggregate_gates$Pass
    ),
    all(
      figure_inventory$Exists
    ),
    isTRUE(
      figure_audit_gate
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  final_summary,
  file.path(
    AUDIT_DIR,
    "10_S20_final_gate_summary.csv"
  ),
  row.names = FALSE
)

blocking <- sum(
  !final_summary$Pass
)

max_numeric_diff <- suppressWarnings(
  max(
    numeric_comparison$Max_numeric_abs_diff,
    na.rm = TRUE
  )
)

if (!is.finite(max_numeric_diff)) {
  max_numeric_diff <- NA_real_
}

cat(
  "\n============================================================\n"
)

cat(
  "10 S20 EXTERNAL SPATIAL RECURRENCE REPRODUCTION CHECK\n"
)

cat(
  "============================================================\n"
)

cat(
  "Historical 18J outputs: ",
  sum(
    hist_inventory$Exists
  ),
  "/",
  nrow(
    hist_inventory
  ),
  "\n",
  sep = ""
)

cat(
  "Historical S28 core outputs: ",
  sum(
    s28_inventory$Exists
  ),
  "/",
  nrow(
    s28_inventory
  ),
  "\n",
  sep = ""
)

cat(
  "Canonical dataset outputs: ",
  sum(
    canon_inventory$Exists
  ),
  "/",
  nrow(
    canon_inventory
  ),
  "\n",
  sep = ""
)

cat(
  "State-gene lock sources identical: ",
  sum(
    state_lock_audit$Identical_to_historical_GSE250636
  ),
  "/",
  nrow(
    state_lock_audit
  ),
  "\n",
  sep = ""
)

cat(
  "Historical -> canonical numeric tables: ",
  sum(
    numeric_comparison$Pass
  ),
  "/",
  nrow(
    numeric_comparison
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum numeric absolute difference: ",
  format(
    max_numeric_diff,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Aggregate outputs: ",
  sum(
    combined_inventory$Exists
  ),
  "/",
  nrow(
    combined_inventory
  ),
  "\n",
  sep = ""
)

cat(
  "Aggregate consistency gates: ",
  sum(
    aggregate_gates$Pass
  ),
  "/",
  nrow(
    aggregate_gates
  ),
  "\n",
  sep = ""
)

cat(
  "Figure/reporting outputs: ",
  sum(
    figure_inventory$Exists
  ),
  "/",
  nrow(
    figure_inventory
  ),
  "\n",
  sep = ""
)

cat(
  "Figure audit consistency: ",
  figure_audit_gate,
  "\n",
  sep = ""
)

cat(
  "Blocking failures: ",
  blocking,
  "\n",
  sep = ""
)

cat(
  "Overall: ",
  ifelse(
    blocking == 0L,
    "PASS",
    "FAIL"
  ),
  "\n",
  sep = ""
)

cat(
  "Decision: ",
  ifelse(
    blocking == 0L,
    "READY_FOR_S20_ANALYSIS_FREEZE",
    "BLOCKED"
  ),
  "\n",
  sep = ""
)

cat(
  "Audit: ",
  AUDIT_DIR,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)

if (blocking > 0L) {
  failed <- final_summary[
    !final_summary$Pass,
    ,
    drop = FALSE
  ]

  stop(
    paste0(
      "10 S20 reproduction audit failed.\n",
      paste(
        paste0(
          failed$Gate,
          " | observed=",
          failed$Observed
        ),
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}
