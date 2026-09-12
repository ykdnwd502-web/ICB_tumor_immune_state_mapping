## 07_Visium_primary_canonical_analysis.R
## Canonical numerical producer for Primary Visium and S15-S19 panel sources.
## No reporting style is changed here because this script draws no figures.

############################################################
## 21_Primary_Visium_spatial_analysis_CANONICAL_v1.0.R
##
## PURPOSE
##   Canonical analytical mother script for the primary Visium
##   melanoma section after 21A3 confirmed the historical Step11 edge
##   vectorization artifact.
##
## ANALYTICAL PRINCIPLES
##   1. Canonical spot/state/QC lineage remains fixed at n = 3,458.
##   2. All kNN graphs are reconstructed directly from spatial_x/spatial_y
##      with FNN::get.knn; no historical malformed edge RDS is used for
##      canonical graph-based statistics.
##   3. Historical 18D observed bivariate statistic is preserved exactly:
##        I_sym = (I_xy + I_yx) / 2
##      with equal k-neighbor weights.
##   4. Canonical permutation null for symmetric bivariate co-variation
##      uses ONE permutation of y per iteration and recomputes the full
##      symmetric statistic. This is more internally coherent than the
##      historical 18D independent-direction permutation implementation.
##   5. A NEW composition-stratified sensitivity null is implemented
##      explicitly as within-joint-quartile permutation of y within
##      CAF_Stromal x Endothelial strata.
##   6. Fisher ORs are descriptive effect-size summaries, not the primary
##      inferential evidence for spatial co-variation.
##
## THIS SCRIPT DOES NOT:
##   - modify the frozen state definitions;
##   - modify the 18E composition-module definitions;
##   - regenerate manuscript text or figure legends;
##   - treat the primary Visium section as independent cohort validation.
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################################################
## 0. Project configuration
############################################################

project_dir <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)

project_dir <- normalizePath(
  project_dir,
  winslash = "/",
  mustWork = TRUE
)

EXPECTED_N <- 3458L
K_VALUES <- c(4L, 6L, 8L, 12L)
MAIN_K <- 6L
N_PERM <- 999L

SEED_UNIVARIATE <- 20260503L
SEED_BIVARIATE <- 20260608L
SEED_COMPOSITION_STRATIFIED <- 20260629L
SEED_DUAL_HIGH <- 20260701L

STATE_COLS <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

PRIMARY_X <- "Myeloid_Treg_Immunosuppressive"
PRIMARY_Y <- "Tumor_dedifferentiation_Stromal_remodeling"

QC_COLS <- c(
  "nCount_Spatial",
  "nFeature_Spatial",
  "percent_mt"
)

COMPOSITION_COLS <- c(
  "CAF_Stromal",
  "Myeloid",
  "T_NK",
  "B_Plasma",
  "Melanoma_Lineage",
  "Endothelial"
)

############################################################
## 1. Inputs
############################################################

step11_meta_rds <- file.path(
  project_dir,
  "results",
  "intermediate",
  "spatial_melanoma_validation",
  "Step11_spatial_analysis_metadata_with_coordinates_and_scores.rds"
)

step11_qc_rds <- file.path(
  project_dir,
  "results",
  "intermediate",
  "spatial_melanoma_validation",
  "Step11_spatial_analysis_metadata_with_QCresidual_scores.rds"
)

legacy_edge_rds <- file.path(
  project_dir,
  "results",
  "intermediate",
  "spatial_melanoma_validation",
  "Step11_spatial_knn6_edge_list.rds"
)

composition_file <- file.path(
  project_dir,
  "results",
  "tables",
  "revision_visium_composition",
  "18E_Visium_spot_state_composition_residualized_scores.csv"
)

historical_18D_moran_file <- file.path(
  project_dir,
  "results",
  "tables",
  "revision_spatial_covariation",
  "Visium_bivariate_spatial_Moran_kNN_sensitivity.csv"
)

required_inputs <- c(
  step11_meta_rds,
  step11_qc_rds,
  legacy_edge_rds,
  composition_file
)

if (any(!file.exists(required_inputs))) {
  stop(
    "Missing required 21C input(s):\n",
    paste(
      required_inputs[
        !file.exists(required_inputs)
      ],
      collapse = "\n"
    ),
    call. = FALSE
  )
}

############################################################
## 2. Output
############################################################

canonical_table_root <- file.path(
  project_dir,
  "results",
  "tables",
  "primary_visium_spatial"
)

canonical_intermediate_root <- file.path(
  project_dir,
  "results",
  "intermediate",
  "primary_visium_spatial"
)

canonical_reporting_root <- file.path(
  project_dir,
  "results",
  "reporting",
  "primary_visium_spatial"
)

table_dir <- canonical_table_root
intermediate_dir <- canonical_intermediate_root
log_dir <- file.path(
  canonical_reporting_root,
  "logs"
)

panel_root <- file.path(
  table_dir,
  "panel_sources"
)

for (
  dd in c(
    table_dir,
    intermediate_dir,
    log_dir,
    panel_root,
    file.path(panel_root, "S15"),
    file.path(panel_root, "S16"),
    file.path(panel_root, "S17"),
    file.path(panel_root, "S18"),
    file.path(panel_root, "S19"),
    canonical_reporting_root
  )
) {
  dir.create(
    dd,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

## Keep this alias only so downstream provenance code from the audited 21C
## implementation remains mechanically simple. It is NOT an audit output path.
audit_root <- canonical_reporting_root

############################################################
## 3. Package gate
############################################################

if (!requireNamespace("FNN", quietly = TRUE)) {
  stop(
    "Package FNN is required for 21C.",
    call. = FALSE
  )
}

############################################################
## 4. Helpers
############################################################

safe_csv <- function(x, file) {
  utils::write.csv(
    x,
    file,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )

  message(
    "Saved: ",
    normalizePath(
      file,
      winslash = "/",
      mustWork = FALSE
    )
  )
}

md5_file <- function(file) {
  if (!file.exists(file)) {
    return(NA_character_)
  }

  unname(
    tools::md5sum(file)
  )
}

as_df <- function(x, label) {
  if (is.data.frame(x)) {
    return(x)
  }

  if (is.matrix(x)) {
    return(
      as.data.frame(
        x,
        stringsAsFactors = FALSE
      )
    )
  }

  stop(
    label,
    " cannot be converted safely to data.frame. Class: ",
    paste(
      class(x),
      collapse = "/"
    ),
    call. = FALSE
  )
}

resolve_ids <- function(df) {
  candidates <- c(
    "Spot",
    "barcode",
    "Barcode",
    "spot",
    "spot_id",
    "SpotID"
  )

  for (cc in candidates) {
    if (
      cc %in%
      colnames(df)
    ) {
      ids <- as.character(
        df[[cc]]
      )

      if (
        length(ids) ==
        nrow(df) &&
        all(
          !is.na(ids) &
          nzchar(ids)
        )
      ) {
        return(
          list(
            ids = ids,
            source = cc
          )
        )
      }
    }
  }

  rn <- rownames(df)

  if (
    !is.null(rn) &&
    length(rn) ==
    nrow(df) &&
    !identical(
      rn,
      as.character(
        seq_len(
          nrow(df)
        )
      )
    )
  ) {
    return(
      list(
        ids = as.character(rn),
        source = "rownames"
      )
    )
  }

  stop(
    "Could not resolve spot IDs.",
    call. = FALSE
  )
}

canonicalize_percent_mt <- function(
  df,
  label
) {
  source_col <- NA_character_
  action <- "UNRESOLVED"

  if (
    "percent_mt" %in%
    colnames(df)
  ) {
    source_col <- "percent_mt"
    action <- "ALREADY_CANONICAL"
  } else if (
    "percent.mt" %in%
    colnames(df)
  ) {
    df$percent_mt <-
      suppressWarnings(
        as.numeric(
          df[["percent.mt"]]
        )
      )

    source_col <- "percent.mt"
    action <- "ALIASED_percent.mt_TO_percent_mt"
  }

  list(
    data = df,
    audit = data.frame(
      Object = label,
      SourceColumn =
        source_col,
      CanonicalColumn =
        "percent_mt",
      Action =
        action,
      CanonicalPresent =
        "percent_mt" %in%
        colnames(df),
      N_finite =
        if (
          "percent_mt" %in%
          colnames(df)
        ) {
          sum(
            is.finite(
              suppressWarnings(
                as.numeric(
                  df[["percent_mt"]]
                )
              )
            )
          )
        } else {
          0L
        },
      stringsAsFactors = FALSE
    )
  )
}

zvec <- function(x) {
  x <- as.numeric(x)

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

make_knn_index <- function(
  coords,
  k
) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"

  FNN::get.knn(
    coords,
    k = k
  )$nn.index
}

knn_index_to_edges <- function(
  knn_idx
) {
  n <- nrow(knn_idx)
  k <- ncol(knn_idx)

  data.frame(
    from =
      rep(
        seq_len(n),
        each = k
      ),
    to =
      as.integer(
        as.vector(
          t(knn_idx)
        )
      ),
    stringsAsFactors = FALSE
  )
}

graph_structure <- function(
  edges,
  n_expected,
  k_expected
) {
  key <- paste(
    edges$from,
    edges$to,
    sep = "->"
  )

  per_source_n <-
    table(
      edges$from
    )

  per_source_unique <-
    vapply(
      split(
        edges$to,
        edges$from
      ),
      function(v) {
        length(
          unique(v)
        )
      },
      integer(1)
    )

  data.frame(
    n_edges =
      nrow(edges),
    n_sources =
      length(
        unique(
          edges$from
        )
      ),
    n_self_loops =
      sum(
        edges$from ==
        edges$to
      ),
    n_duplicate_directed_rows =
      sum(
        duplicated(key)
      ),
    min_rows_per_source =
      min(per_source_n),
    max_rows_per_source =
      max(per_source_n),
    min_unique_neighbors =
      min(per_source_unique),
    max_unique_neighbors =
      max(per_source_unique),
    exact_structure =
      nrow(edges) ==
      n_expected *
      k_expected &&
      length(
        unique(
          edges$from
        )
      ) ==
      n_expected &&
      sum(
        edges$from ==
        edges$to
      ) ==
      0L &&
      sum(
        duplicated(key)
      ) ==
      0L &&
      all(
        per_source_unique ==
        k_expected
      ),
    stringsAsFactors = FALSE
  )
}

neighbor_mean <- function(
  v,
  knn_idx
) {
  m <- matrix(
    v[
      knn_idx
    ],
    nrow =
      nrow(
        knn_idx
      ),
    ncol =
      ncol(
        knn_idx
      )
  )

  rowMeans(
    m,
    na.rm = TRUE
  )
}

univariate_moran_from_index <- function(
  x,
  knn_idx
) {
  x <- as.numeric(x)

  ok <- is.finite(x)

  if (!all(ok)) {
    stop(
      "21C univariate Moran expects a complete analysis vector.",
      call. = FALSE
    )
  }

  z <- x -
    mean(
      x
    )

  denom <- sum(
    z^2
  )

  if (
    !is.finite(denom) ||
    denom == 0
  ) {
    return(NA_real_)
  }

  lag_z <- neighbor_mean(
    z,
    knn_idx
  )

  ## Equivalent to historical binary-edge Moran with S0 = n*k,
  ## but computed directly from the valid row-aligned kNN index.
  sum(
    z *
    lag_z
  ) /
  denom
}

univariate_moran_permutation <- function(
  x,
  knn_idx,
  n_perm,
  seed
) {
  obs <-
    univariate_moran_from_index(
      x,
      knn_idx
    )

  set.seed(seed)

  perm <- vapply(
    seq_len(n_perm),
    function(i) {
      xp <- sample(
        x,
        length(x),
        replace = FALSE
      )

      univariate_moran_from_index(
        xp,
        knn_idx
      )
    },
    numeric(1)
  )

  p_two <-
    (
      sum(
        abs(perm) >=
        abs(obs),
        na.rm = TRUE
      ) + 1
    ) /
    (
      sum(
        is.finite(
          perm
        )
      ) + 1
    )

  list(
    observed = obs,
    permutation = perm,
    p_two_sided = p_two
  )
}

symmetric_bivariate_from_z <- function(
  zx,
  zy,
  knn_idx
) {
  lag_y <-
    neighbor_mean(
      zy,
      knn_idx
    )

  lag_x <-
    neighbor_mean(
      zx,
      knn_idx
    )

  I_xy <- mean(
    zx *
    lag_y
  )

  I_yx <- mean(
    zy *
    lag_x
  )

  c(
    I_xy = I_xy,
    I_yx = I_yx,
    I_symmetric =
      mean(
        c(
          I_xy,
          I_yx
        )
      )
  )
}

symmetric_bivariate <- function(
  x,
  y,
  knn_idx
) {
  symmetric_bivariate_from_z(
    zvec(x),
    zvec(y),
    knn_idx
  )
}

rank_quantile_bin <- function(
  x,
  n_bin = 4L
) {
  x <- as.numeric(x)

  if (
    any(
      !is.finite(x)
    )
  ) {
    stop(
      "Composition stratum variable contains non-finite values.",
      call. = FALSE
    )
  }

  r <- rank(
    x,
    ties.method = "average"
  )

  n <- length(r)

  out <-
    floor(
      (
        r - 1
      ) /
      n *
      n_bin
    ) + 1L

  out[
    out < 1L
  ] <- 1L

  out[
    out > n_bin
  ] <- n_bin

  as.integer(out)
}

build_composition_strata <- function(
  caf,
  endothelial,
  n_bin = 4L
) {
  caf_bin <-
    rank_quantile_bin(
      caf,
      n_bin
    )

  endo_bin <-
    rank_quantile_bin(
      endothelial,
      n_bin
    )

  interaction(
    paste0(
      "CAFQ",
      caf_bin
    ),
    paste0(
      "EndoQ",
      endo_bin
    ),
    drop = TRUE,
    lex.order = TRUE
  )
}

permute_within_strata <- function(
  v,
  strata
) {
  out <- v

  split_idx <-
    split(
      seq_along(v),
      strata
    )

  for (
    idx in split_idx
  ) {
    if (
      length(idx) > 1L
    ) {
      out[idx] <-
        sample(
          v[idx],
          length(idx),
          replace = FALSE
        )
    }
  }

  out
}

bivariate_permutation <- function(
  x,
  y,
  knn_idx,
  n_perm,
  seed,
  strata = NULL
) {
  zx <- zvec(x)
  zy <- zvec(y)

  obs_vec <-
    symmetric_bivariate_from_z(
      zx,
      zy,
      knn_idx
    )

  obs <-
    unname(
      obs_vec[
        "I_symmetric"
      ]
    )

  set.seed(seed)

  perm_vals <- vapply(
    seq_len(n_perm),
    function(i) {
      zy_perm <-
        if (
          is.null(strata)
        ) {
          sample(
            zy,
            length(zy),
            replace = FALSE
          )
        } else {
          permute_within_strata(
            zy,
            strata
          )
        }

      unname(
        symmetric_bivariate_from_z(
          zx,
          zy_perm,
          knn_idx
        )[
          "I_symmetric"
        ]
      )
    },
    numeric(1)
  )

  perm_mean <- mean(
    perm_vals,
    na.rm = TRUE
  )

  p_greater <-
    (
      sum(
        perm_vals >=
        obs,
        na.rm = TRUE
      ) + 1
    ) /
    (
      n_perm + 1
    )

  p_two_centered <-
    (
      sum(
        abs(
          perm_vals -
          perm_mean
        ) >=
        abs(
          obs -
          perm_mean
        ),
        na.rm = TRUE
      ) + 1
    ) /
    (
      n_perm + 1
    )

  list(
    observed = obs_vec,
    perm = perm_vals,
    perm_mean = perm_mean,
    perm_sd =
      stats::sd(
        perm_vals,
        na.rm = TRUE
      ),
    p_greater =
      p_greater,
    p_two_centered =
      p_two_centered
  )
}

partial_spearman_residual <- function(
  df,
  x_col,
  y_col,
  covariates
) {
  covariates <-
    covariates[
      covariates %in%
      colnames(df)
    ]

  vars <- c(
    x_col,
    y_col,
    covariates
  )

  dat <-
    df[
      ,
      vars,
      drop = FALSE
    ]

  dat <-
    dat[
      stats::complete.cases(
        dat
      ),
      ,
      drop = FALSE
    ]

  rank_df <-
    as.data.frame(
      lapply(
        dat,
        function(v) {
          rank(
            as.numeric(v),
            ties.method = "average",
            na.last = "keep"
          )
        }
      )
    )

  fml_x <-
    as.formula(
      paste(
        x_col,
        "~",
        paste(
          covariates,
          collapse = " + "
        )
      )
    )

  fml_y <-
    as.formula(
      paste(
        y_col,
        "~",
        paste(
          covariates,
          collapse = " + "
        )
      )
    )

  rx <-
    residuals(
      lm(
        fml_x,
        data = rank_df
      )
    )

  ry <-
    residuals(
      lm(
        fml_y,
        data = rank_df
      )
    )

  ct <-
    suppressWarnings(
      stats::cor.test(
        rx,
        ry,
        method = "pearson"
      )
    )

  data.frame(
    n = nrow(dat),
    estimate =
      unname(
        ct$estimate
      ),
    p_value =
      ct$p.value,
    method =
      "partial Spearman by rank residualization",
    covariates =
      paste(
        covariates,
        collapse = " + "
      ),
    stringsAsFactors = FALSE
  )
}

safe_fisher_or <- function(
  a,
  b,
  c,
  d
) {
  mat <- matrix(
    c(
      a,
      b,
      c,
      d
    ),
    nrow = 2,
    byrow = TRUE
  )

  ft <-
    suppressWarnings(
      fisher.test(
        mat
      )
    )

  data.frame(
    OR =
      unname(
        ft$estimate
      ),
    CI_low =
      ft$conf.int[[1]],
    CI_high =
      ft$conf.int[[2]],
    p_value =
      ft$p.value,
    stringsAsFactors = FALSE
  )
}

dual_high_counts <- function(
  x,
  y,
  q = 0.75
) {
  x_cut <-
    as.numeric(
      stats::quantile(
        x,
        probs = q,
        na.rm = TRUE,
        names = FALSE
      )
    )

  y_cut <-
    as.numeric(
      stats::quantile(
        y,
        probs = q,
        na.rm = TRUE,
        names = FALSE
      )
    )

  x_high <-
    x >=
    x_cut

  y_high <-
    y >=
    y_cut

  both <-
    sum(
      x_high &
      y_high,
      na.rm = TRUE
    )

  x_only <-
    sum(
      x_high &
      !y_high,
      na.rm = TRUE
    )

  y_only <-
    sum(
      !x_high &
      y_high,
      na.rm = TRUE
    )

  neither <-
    sum(
      !x_high &
      !y_high,
      na.rm = TRUE
    )

  list(
    x_cut = x_cut,
    y_cut = y_cut,
    x_high = x_high,
    y_high = y_high,
    both = both,
    x_only = x_only,
    y_only = y_only,
    neither = neither
  )
}

dual_high_permutation <- function(
  x,
  y,
  q,
  n_perm,
  seed,
  strata = NULL
) {
  obs <-
    dual_high_counts(
      x,
      y,
      q = q
    )

  set.seed(seed)

  perm_both <- vapply(
    seq_len(n_perm),
    function(i) {
      yp <-
        if (
          is.null(strata)
        ) {
          sample(
            y,
            length(y),
            replace = FALSE
          )
        } else {
          permute_within_strata(
            y,
            strata
          )
        }

      dual_high_counts(
        x,
        yp,
        q = q
      )$both
    },
    numeric(1)
  )

  p_greater <-
    (
      sum(
        perm_both >=
        obs$both
      ) + 1
    ) /
    (
      n_perm + 1
    )

  list(
    observed = obs,
    perm_both = perm_both,
    p_greater = p_greater
  )
}

threshold_grid_or <- function(
  x,
  y
) {
  defs <- data.frame(
    CutoffName = c(
      "top15",
      "top20",
      "top25",
      "top30",
      "top35",
      "median"
    ),
    Quantile = c(
      0.85,
      0.80,
      0.75,
      0.70,
      0.65,
      0.50
    ),
    stringsAsFactors = FALSE
  )

  rows <- lapply(
    seq_len(
      nrow(defs)
    ),
    function(i) {
      q <-
        defs$Quantile[[i]]

      cc <-
        dual_high_counts(
          x,
          y,
          q = q
        )

      or <-
        safe_fisher_or(
          cc$both,
          cc$x_only,
          cc$y_only,
          cc$neither
        )

      cbind(
        data.frame(
          CutoffName =
            defs$CutoffName[[i]],
          QuantileHighDefinition =
            q,
          x_cutoff =
            cc$x_cut,
          y_cutoff =
            cc$y_cut,
          BothHigh =
            cc$both,
          XHighOnly =
            cc$x_only,
          YHighOnly =
            cc$y_only,
          NeitherHigh =
            cc$neither,
          stringsAsFactors = FALSE
        ),
        or
      )
    }
  )

  x_high <-
    x > 0

  y_high <-
    y > 0

  both <-
    sum(
      x_high &
      y_high
    )

  x_only <-
    sum(
      x_high &
      !y_high
    )

  y_only <-
    sum(
      !x_high &
      y_high
    )

  neither <-
    sum(
      !x_high &
      !y_high
    )

  or <-
    safe_fisher_or(
      both,
      x_only,
      y_only,
      neither
    )

  rows[[length(rows) + 1L]] <- cbind(
    data.frame(
      CutoffName =
        "z_gt_0",
      QuantileHighDefinition =
        NA_real_,
      x_cutoff = 0,
      y_cutoff = 0,
      BothHigh = both,
      XHighOnly = x_only,
      YHighOnly = y_only,
      NeitherHigh = neither,
      stringsAsFactors = FALSE
    ),
    or
  )

  do.call(
    rbind,
    rows
  )
}

gate_row <- function(
  Gate,
  Observed,
  Expected,
  Pass,
  Classification = "BLOCKING"
) {
  data.frame(
    Gate =
      as.character(Gate),
    Observed =
      as.character(Observed),
    Expected =
      as.character(Expected),
    Pass =
      isTRUE(Pass),
    Classification =
      as.character(
        Classification
      ),
    stringsAsFactors = FALSE
  )
}

############################################################
## 5. Load and align canonical spot-level inputs
############################################################

meta <-
  as_df(
    readRDS(
      step11_meta_rds
    ),
    "Step11 metadata"
  )

qc <-
  as_df(
    readRDS(
      step11_qc_rds
    ),
    "Step11 QC-residual metadata"
  )

legacy_edges <-
  as_df(
    readRDS(
      legacy_edge_rds
    ),
    "Historical kNN6 edge list"
  )

comp <-
  utils::read.csv(
    composition_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

meta_qc_alias <-
  canonicalize_percent_mt(
    meta,
    "Step11_metadata"
  )

meta <-
  meta_qc_alias$data

qc_qc_alias <-
  canonicalize_percent_mt(
    qc,
    "Step11_QCresidual_metadata"
  )

qc <-
  qc_qc_alias$data

comp_qc_alias <-
  canonicalize_percent_mt(
    comp,
    "18E_composition_table"
  )

comp <-
  comp_qc_alias$data

qc_alias_audit <-
  do.call(
    rbind,
    list(
      meta_qc_alias$audit,
      qc_qc_alias$audit,
      comp_qc_alias$audit
    )
  )

safe_csv(
  qc_alias_audit,
  file.path(
    table_dir,
    "21_QC_percent_mt_alias_audit.csv"
  )
)

meta_id <-
  resolve_ids(meta)

qc_id <-
  resolve_ids(qc)

comp_id <-
  resolve_ids(comp)

meta_ids <-
  meta_id$ids

qc_ids <-
  qc_id$ids

comp_ids <-
  comp_id$ids

if (
  length(
    unique(
      meta_ids
    )
  ) !=
  EXPECTED_N
) {
  stop(
    "Canonical Step11 spot count is not 3458.",
    call. = FALSE
  )
}

qc_match <-
  match(
    meta_ids,
    qc_ids
  )

comp_match <-
  match(
    meta_ids,
    comp_ids
  )

if (
  anyNA(qc_match) ||
  anyNA(comp_match)
) {
  stop(
    "QC or composition table does not recover the exact canonical spot universe.",
    call. = FALSE
  )
}

qc <-
  qc[
    qc_match,
    ,
    drop = FALSE
  ]

comp <-
  comp[
    comp_match,
    ,
    drop = FALSE
  ]

############################################################
## 6. Lineage and column gates
############################################################

required_meta_cols <- c(
  "spatial_x",
  "spatial_y",
  QC_COLS,
  STATE_COLS
)

missing_meta <-
  setdiff(
    required_meta_cols,
    colnames(meta)
  )

if (
  length(
    missing_meta
  ) > 0L
) {
  stop(
    "Step11 metadata missing columns: ",
    paste(
      missing_meta,
      collapse = ", "
    ),
    call. = FALSE
  )
}

qc_state_cols <-
  paste0(
    STATE_COLS,
    "_QCresid"
  )

missing_qc <-
  setdiff(
    qc_state_cols,
    colnames(qc)
  )

if (
  length(
    missing_qc
  ) > 0L
) {
  stop(
    "QC-residual table missing columns: ",
    paste(
      missing_qc,
      collapse = ", "
    ),
    call. = FALSE
  )
}

score_pairs <- list(
  Raw = c(
    x =
      PRIMARY_X,
    y =
      PRIMARY_Y
  ),
  QC_residualized = c(
    x =
      "MyeloidTreg_resid_QC",
    y =
      "DediffStromal_resid_QC"
  ),
  Composition_residualized = c(
    x =
      "MyeloidTreg_resid_Composition",
    y =
      "DediffStromal_resid_Composition"
  ),
  QC_Composition_residualized = c(
    x =
      "MyeloidTreg_resid_QC_Composition",
    y =
      "DediffStromal_resid_QC_Composition"
  )
)

required_comp_cols <- unique(
  c(
    STATE_COLS,
    COMPOSITION_COLS,
    unname(
      unlist(
        score_pairs
      )
    )
  )
)

missing_comp <-
  setdiff(
    required_comp_cols,
    colnames(comp)
  )

if (
  length(
    missing_comp
  ) > 0L
) {
  stop(
    "18E composition table missing required columns: ",
    paste(
      missing_comp,
      collapse = ", "
    ),
    call. = FALSE
  )
}

state_equivalence <- do.call(
  rbind,
  lapply(
    STATE_COLS,
    function(st) {
      a <-
        as.numeric(
          meta[[st]]
        )

      b <-
        as.numeric(
          comp[[st]]
        )

      d <-
        abs(
          a - b
        )

      data.frame(
        State = st,
        N =
          length(a),
        N_finite_meta =
          sum(
            is.finite(a)
          ),
        N_finite_composition =
          sum(
            is.finite(b)
          ),
        MaxAbsDiff =
          max(
            d,
            na.rm = TRUE
          ),
        Exact_1e10 =
          max(
            d,
            na.rm = TRUE
          ) <=
          1e-10,
        stringsAsFactors = FALSE
      )
    }
  )
)

safe_csv(
  state_equivalence,
  file.path(
    table_dir,
    "21_state_score_equivalence_Step11_vs_18E.csv"
  )
)

coords <-
  as.matrix(
    meta[
      ,
      c(
        "spatial_x",
        "spatial_y"
      )
    ]
  )

storage.mode(coords) <-
  "numeric"

coord_complete <-
  apply(
    coords,
    1,
    function(v) {
      all(
        is.finite(v)
      )
    }
  )

coord_key <-
  paste(
    format(
      coords[, 1],
      digits = 17,
      scientific = FALSE,
      trim = TRUE
    ),
    format(
      coords[, 2],
      digits = 17,
      scientific = FALSE,
      trim = TRUE
    ),
    sep = "||"
  )

############################################################
## 7. Common complete universe for primary-pair sensitivity
############################################################

pair_complete <-
  coord_complete

for (
  nm in names(
    score_pairs
  )
) {
  pair <-
    score_pairs[[nm]]

  pair_complete <-
    pair_complete &
    is.finite(
      as.numeric(
        comp[[pair[["x"]]]]
      )
    ) &
    is.finite(
      as.numeric(
        comp[[pair[["y"]]]]
      )
    )
}

for (
  cc in COMPOSITION_COLS
) {
  pair_complete <-
    pair_complete &
    is.finite(
      as.numeric(
        comp[[cc]]
      )
    )
}

common_idx <-
  which(
    pair_complete
  )

common_n <-
  length(
    common_idx
  )

common_meta <-
  meta[
    common_idx,
    ,
    drop = FALSE
  ]

common_comp <-
  comp[
    common_idx,
    ,
    drop = FALSE
  ]

common_coords <-
  coords[
    common_idx,
    ,
    drop = FALSE
  ]

############################################################
## 8. Rebuild corrected kNN graphs
############################################################

graph_index <- list()
graph_edges <- list()
graph_audit_rows <- list()

for (
  k in K_VALUES
) {
  idx <-
    make_knn_index(
      common_coords,
      k = k
    )

  ed <-
    knn_index_to_edges(
      idx
    )

  ga <-
    graph_structure(
      ed,
      n_expected =
        common_n,
      k_expected =
        k
    )

  ga$k <- k

  graph_index[[as.character(k)]] <- idx

  graph_edges[[as.character(k)]] <- ed

  graph_audit_rows[[length(graph_audit_rows) + 1L]] <- ga
}

graph_audit <-
  do.call(
    rbind,
    graph_audit_rows
  )

graph_audit <-
  graph_audit[
    ,
    c(
      "k",
      setdiff(
        colnames(
          graph_audit
        ),
        "k"
      )
    )
  ]

safe_csv(
  graph_audit,
  file.path(
    table_dir,
    "21_corrected_knn_graph_audit.csv"
  )
)

saveRDS(
  graph_index,
  file.path(
    intermediate_dir,
    "21_corrected_knn_indices_k4_k6_k8_k12.rds"
  )
)

saveRDS(
  graph_edges[[as.character(MAIN_K)]],
  file.path(
    intermediate_dir,
    "21_corrected_knn6_edge_list.rds"
  )
)

############################################################
## 9. Historical malformed graph vs corrected k6 pair-set audit
############################################################

legacy_from <-
  as.integer(
    legacy_edges$from
  )

legacy_to <-
  as.integer(
    legacy_edges$to
  )

correct_k6_edges <-
  graph_edges[[as.character(MAIN_K)]]

legacy_key <-
  unique(
    paste(
      legacy_from,
      legacy_to,
      sep = "->"
    )
  )

correct_key <-
  unique(
    paste(
      correct_k6_edges$from,
      correct_k6_edges$to,
      sep = "->"
    )
  )

graph_pair_comparison <- data.frame(
  Historical_unique_pairs =
    length(
      legacy_key
    ),
  Corrected_unique_pairs =
    length(
      correct_key
    ),
  Intersection =
    length(
      intersect(
        legacy_key,
        correct_key
      )
    ),
  Historical_minus_corrected =
    length(
      setdiff(
        legacy_key,
        correct_key
      )
    ),
  Corrected_minus_historical =
    length(
      setdiff(
        correct_key,
        legacy_key
      )
    ),
  Jaccard =
    length(
      intersect(
        legacy_key,
        correct_key
      )
    ) /
    length(
      union(
        legacy_key,
        correct_key
      )
    ),
  Historical_self_loops =
    sum(
      legacy_from ==
      legacy_to
    ),
  Historical_duplicate_rows =
    sum(
      duplicated(
        paste(
          legacy_from,
          legacy_to,
          sep = "->"
        )
      )
    ),
  stringsAsFactors = FALSE
)

safe_csv(
  graph_pair_comparison,
  file.path(
    table_dir,
    "21_historical_vs_corrected_knn6_pair_set_audit.csv"
  )
)

############################################################
## 10. Corrected univariate Moran for four states, raw + QC
##     Main k = 6; historical permutation definition retained.
############################################################

k6_idx <-
  graph_index[[as.character(MAIN_K)]]

univariate_rows <- list()
univariate_null_rows <- list()

for (
  s_idx in seq_along(
    STATE_COLS
  )
) {
  st <-
    STATE_COLS[[s_idx]]

  raw_x <-
    as.numeric(
      common_meta[[st]]
    )

  qc_col <-
    paste0(
      st,
      "_QCresid"
    )

  qc_x <-
    as.numeric(
      qc[
        common_idx,
        qc_col
      ]
    )

  variants <- list(
    Raw = raw_x,
    QC_residualized = qc_x
  )

  for (
    v_idx in seq_along(
      variants
    )
  ) {
    label <-
      names(
        variants
      )[[v_idx]]

    x <-
      variants[[v_idx]]

    seed <-
      SEED_UNIVARIATE +
      s_idx *
      100L +
      v_idx

    rr <-
      univariate_moran_permutation(
        x,
        k6_idx,
        n_perm = N_PERM,
        seed = seed
      )

    univariate_rows[[length(univariate_rows) + 1L]] <- data.frame(
      State = st,
      ScoreType = label,
      k = MAIN_K,
      n = length(x),
      Moran_I =
        rr$observed,
      p_empirical_two_sided =
        rr$p_two_sided,
      n_perm = N_PERM,
      seed = seed,
      Graph =
        "Corrected row-aligned FNN kNN",
      stringsAsFactors = FALSE
    )

    univariate_null_rows[[length(univariate_null_rows) + 1L]] <- data.frame(
      State = st,
      ScoreType = label,
      k = MAIN_K,
      permutation =
        seq_len(
          N_PERM
        ),
      Moran_I_perm =
        rr$permutation,
      seed = seed,
      stringsAsFactors = FALSE
    )
  }
}

univariate_moran <-
  do.call(
    rbind,
    univariate_rows
  )

univariate_moran$FDR_BH <-
  p.adjust(
    univariate_moran$
      p_empirical_two_sided,
    method = "BH"
  )

univariate_null <-
  do.call(
    rbind,
    univariate_null_rows
  )

safe_csv(
  univariate_moran,
  file.path(
    table_dir,
    "21_corrected_univariate_Moran_raw_QC_k6.csv"
  )
)

safe_csv(
  univariate_null,
  file.path(
    table_dir,
    "21_corrected_univariate_Moran_permutation_null_k6.csv"
  )
)

############################################################
## 11. Corrected descriptive high-state neighbor enrichment
############################################################

neighbor_rows <- list()

for (
  st in STATE_COLS
) {
  x <-
    as.numeric(
      common_meta[[st]]
    )

  cutoff <-
    as.numeric(
      stats::quantile(
        x,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE
      )
    )

  high <-
    x >=
    cutoff

  a <-
    high[
      correct_k6_edges$from
    ]

  b <-
    high[
      correct_k6_edges$to
    ]

  tab <-
    table(
      factor(
        a,
        levels = c(
          TRUE,
          FALSE
        )
      ),
      factor(
        b,
        levels = c(
          TRUE,
          FALSE
        )
      )
    )

  m <-
    matrix(
      as.numeric(
        tab
      ),
      nrow = 2
    )

  or <-
    safe_fisher_or(
      m[1, 1],
      m[1, 2],
      m[2, 1],
      m[2, 2]
    )

  neighbor_rows[[length(neighbor_rows) + 1L]] <- cbind(
    data.frame(
      State = st,
      k = MAIN_K,
      cutoff_quantile =
        0.75,
      cutoff =
        cutoff,
      n_source_high_neighbor_high =
        m[1, 1],
      n_source_high_neighbor_not_high =
        m[1, 2],
      n_source_not_high_neighbor_high =
        m[2, 1],
      n_source_not_high_neighbor_not_high =
        m[2, 2],
      InferentialStatus =
        "DESCRIPTIVE_EDGE_LEVEL_EFFECT_SIZE",
      stringsAsFactors = FALSE
    ),
    or
  )
}

neighbor_enrichment <-
  do.call(
    rbind,
    neighbor_rows
  )

safe_csv(
  neighbor_enrichment,
  file.path(
    table_dir,
    "21_corrected_high_state_neighbor_enrichment_top25_k6_DESCRIPTIVE.csv"
  )
)

############################################################
## 12. Primary pair continuous association
############################################################

raw_x <-
  as.numeric(
    common_meta[[PRIMARY_X]]
  )

raw_y <-
  as.numeric(
    common_meta[[PRIMARY_Y]]
  )

sp <-
  suppressWarnings(
    stats::cor.test(
      raw_x,
      raw_y,
      method = "spearman",
      exact = FALSE
    )
  )

pe <-
  suppressWarnings(
    stats::cor.test(
      raw_x,
      raw_y,
      method = "pearson"
    )
  )

partial <-
  partial_spearman_residual(
    common_meta,
    PRIMARY_X,
    PRIMARY_Y,
    QC_COLS
  )

continuous_summary <- rbind(
  data.frame(
    Analysis =
      "Raw Spearman",
    n =
      length(
        raw_x
      ),
    estimate =
      unname(
        sp$estimate
      ),
    p_value =
      sp$p.value,
    method =
      "Spearman",
    covariates =
      "None",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis =
      "Raw Pearson",
    n =
      length(
        raw_x
      ),
    estimate =
      unname(
        pe$estimate
      ),
    p_value =
      pe$p.value,
    method =
      "Pearson",
    covariates =
      "None",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis =
      "QC-adjusted partial Spearman",
    n =
      partial$n,
    estimate =
      partial$estimate,
    p_value =
      partial$p_value,
    method =
      partial$method,
    covariates =
      partial$covariates,
    stringsAsFactors = FALSE
  )
)

safe_csv(
  continuous_summary,
  file.path(
    table_dir,
    "21_primary_pair_continuous_association.csv"
  )
)

############################################################
## 13. Raw threshold-grid Fisher OR sensitivity
############################################################

threshold_table <-
  threshold_grid_or(
    raw_x,
    raw_y
  )

threshold_table$InferentialStatus <-
  "DESCRIPTIVE_OPERATIONAL_SENSITIVITY"

safe_csv(
  threshold_table,
  file.path(
    table_dir,
    "21_primary_pair_threshold_grid_Fisher_OR_DESCRIPTIVE.csv"
  )
)

############################################################
## 14. Composition strata
############################################################

strata <-
  build_composition_strata(
    caf =
      as.numeric(
        common_comp$
          CAF_Stromal
      ),
    endothelial =
      as.numeric(
        common_comp$
          Endothelial
      ),
    n_bin = 4L
  )

strata_tab <-
  as.data.frame(
    table(
      strata
    ),
    stringsAsFactors = FALSE
  )

colnames(
  strata_tab
) <- c(
  "Stratum",
  "N"
)

strata_audit <- data.frame(
  n_spots =
    length(
      strata
    ),
  n_strata =
    length(
      unique(
        strata
      )
    ),
  min_stratum_n =
    min(
      strata_tab$N
    ),
  median_stratum_n =
    stats::median(
      strata_tab$N
    ),
  max_stratum_n =
    max(
      strata_tab$N
    ),
  Definition =
    "Joint rank-quartile strata: CAF_Stromal x Endothelial",
  Status =
    "NEW_CANONICAL_COMPOSITION_STRATIFIED_SENSITIVITY",
  stringsAsFactors = FALSE
)

safe_csv(
  strata_tab,
  file.path(
    table_dir,
    "21_composition_strata_counts.csv"
  )
)

safe_csv(
  strata_audit,
  file.path(
    table_dir,
    "21_composition_strata_definition.csv"
  )
)

############################################################
## 15. Bivariate observed statistics across all methods and k
############################################################

biv_observed_rows <- list()

for (
  method_i in seq_along(
    score_pairs
  )
) {
  method <-
    names(
      score_pairs
    )[[method_i]]

  pair <-
    score_pairs[[method_i]]

  x <-
    as.numeric(
      common_comp[[pair[["x"]]]]
    )

  y <-
    as.numeric(
      common_comp[[pair[["y"]]]]
    )

  for (
    k in K_VALUES
  ) {
    obs <-
      symmetric_bivariate(
        x,
        y,
        graph_index[[as.character(k)]]
      )

    biv_observed_rows[[length(biv_observed_rows) + 1L]] <- data.frame(
      Method = method,
      k = k,
      n = length(x),
      I_xy =
        unname(
          obs[
            "I_xy"
          ]
        ),
      I_yx =
        unname(
          obs[
            "I_yx"
          ]
        ),
      I_symmetric =
        unname(
          obs[
            "I_symmetric"
          ]
        ),
      Graph =
        "Direct row-aligned FNN::get.knn index",
      stringsAsFactors = FALSE
    )
  }
}

biv_observed <-
  do.call(
    rbind,
    biv_observed_rows
  )

safe_csv(
  biv_observed,
  file.path(
    table_dir,
    "21_bivariate_Moran_observed_all_methods_k4_k6_k8_k12.csv"
  )
)

############################################################
## 16. Raw bivariate k sensitivity with unrestricted permutation
############################################################

raw_k_perm_rows <- list()
raw_k_null_rows <- list()

for (
  k in K_VALUES
) {
  rr <-
    bivariate_permutation(
      x =
        as.numeric(
          common_comp[[score_pairs$Raw[["x"]]]]
        ),
      y =
        as.numeric(
          common_comp[[score_pairs$Raw[["y"]]]]
        ),
      knn_idx =
        graph_index[[as.character(k)]],
      n_perm =
        N_PERM,
      seed =
        SEED_BIVARIATE +
        k,
      strata =
        NULL
    )

  raw_k_perm_rows[[length(raw_k_perm_rows) + 1L]] <- data.frame(
    Method = "Raw",
    k = k,
    n = common_n,
    I_xy =
      unname(
        rr$observed[
          "I_xy"
        ]
      ),
    I_yx =
      unname(
        rr$observed[
          "I_yx"
        ]
      ),
    I_symmetric =
      unname(
        rr$observed[
          "I_symmetric"
        ]
      ),
    perm_mean =
      rr$perm_mean,
    perm_sd =
      rr$perm_sd,
    p_greater =
      rr$p_greater,
    p_two_centered =
      rr$p_two_centered,
    n_perm =
      N_PERM,
    seed =
      SEED_BIVARIATE +
      k,
    Null =
      "UNRESTRICTED_SINGLE_Y_PERMUTATION",
    stringsAsFactors = FALSE
  )

  raw_k_null_rows[[length(raw_k_null_rows) + 1L]] <- data.frame(
    Method = "Raw",
    k = k,
    permutation =
      seq_len(
        N_PERM
      ),
    I_symmetric_perm =
      rr$perm,
    stringsAsFactors = FALSE
  )
}

raw_k_perm <-
  do.call(
    rbind,
    raw_k_perm_rows
  )

raw_k_null <-
  do.call(
    rbind,
    raw_k_null_rows
  )

safe_csv(
  raw_k_perm,
  file.path(
    table_dir,
    "21_raw_bivariate_Moran_k_sensitivity_unrestricted_permutation.csv"
  )
)

safe_csv(
  raw_k_null,
  file.path(
    table_dir,
    "21_raw_bivariate_Moran_k_sensitivity_nulls.csv"
  )
)

############################################################
## 17. Main k=6 composition-aware bivariate permutation
##     unrestricted + NEW CAF/endothelial-stratified null
############################################################

composition_biv_rows <- list()
composition_biv_null_rows <- list()

for (
  method_i in seq_along(
    score_pairs
  )
) {
  method <-
    names(
      score_pairs
    )[[method_i]]

  pair <-
    score_pairs[[method_i]]

  x <-
    as.numeric(
      common_comp[[pair[["x"]]]]
    )

  y <-
    as.numeric(
      common_comp[[pair[["y"]]]]
    )

  null_defs <- list(
    Unrestricted = NULL,
    CAF_Endothelial_joint_quartile_stratified =
      strata
  )

  for (
    null_i in seq_along(
      null_defs
    )
  ) {
    null_name <-
      names(
        null_defs
      )[[null_i]]

    seed <-
      if (
        null_name ==
        "Unrestricted"
      ) {
        SEED_BIVARIATE +
        1000L +
        method_i
      } else {
        SEED_COMPOSITION_STRATIFIED +
        1000L +
        method_i
      }

    rr <-
      bivariate_permutation(
        x,
        y,
        k6_idx,
        n_perm =
          N_PERM,
        seed =
          seed,
        strata =
          null_defs[[null_i]]
      )

    composition_biv_rows[[length(composition_biv_rows) + 1L]] <- data.frame(
      Method = method,
      k = MAIN_K,
      n = length(x),
      Null = null_name,
      I_xy =
        unname(
          rr$observed[
            "I_xy"
          ]
        ),
      I_yx =
        unname(
          rr$observed[
            "I_yx"
          ]
        ),
      I_symmetric =
        unname(
          rr$observed[
            "I_symmetric"
          ]
        ),
      perm_mean =
        rr$perm_mean,
      perm_sd =
        rr$perm_sd,
      p_greater =
        rr$p_greater,
      p_two_centered =
        rr$p_two_centered,
      n_perm =
        N_PERM,
      seed =
        seed,
      stringsAsFactors = FALSE
    )

    composition_biv_null_rows[[length(composition_biv_null_rows) + 1L]] <- data.frame(
      Method = method,
      k = MAIN_K,
      Null = null_name,
      permutation =
        seq_len(
          N_PERM
        ),
      I_symmetric_perm =
        rr$perm,
      seed =
        seed,
      stringsAsFactors = FALSE
    )
  }
}

composition_biv <-
  do.call(
    rbind,
    composition_biv_rows
  )

composition_biv_null <-
  do.call(
    rbind,
    composition_biv_null_rows
  )

safe_csv(
  composition_biv,
  file.path(
    table_dir,
    "21_composition_aware_bivariate_Moran_k6_permutation_summary.csv"
  )
)

safe_csv(
  composition_biv_null,
  file.path(
    table_dir,
    "21_composition_aware_bivariate_Moran_k6_permutation_nulls.csv"
  )
)

############################################################
## 18. Raw/residualized top-quartile overlap permutation
##     unrestricted + composition-stratified
############################################################

dual_perm_rows <- list()
dual_perm_null_rows <- list()

for (
  method_i in seq_along(
    score_pairs
  )
) {
  method <-
    names(
      score_pairs
    )[[method_i]]

  pair <-
    score_pairs[[method_i]]

  x <-
    as.numeric(
      common_comp[[pair[["x"]]]]
    )

  y <-
    as.numeric(
      common_comp[[pair[["y"]]]]
    )

  obs <-
    dual_high_counts(
      x,
      y,
      q = 0.75
    )

  or <-
    safe_fisher_or(
      obs$both,
      obs$x_only,
      obs$y_only,
      obs$neither
    )

  null_defs <- list(
    Unrestricted = NULL,
    CAF_Endothelial_joint_quartile_stratified =
      strata
  )

  for (
    null_i in seq_along(
      null_defs
    )
  ) {
    null_name <-
      names(
        null_defs
      )[[null_i]]

    seed <-
      SEED_DUAL_HIGH +
      method_i *
      100L +
      null_i

    rr <-
      dual_high_permutation(
        x,
        y,
        q = 0.75,
        n_perm =
          N_PERM,
        seed =
          seed,
        strata =
          null_defs[[null_i]]
      )

    dual_perm_rows[[length(dual_perm_rows) + 1L]] <- cbind(
      data.frame(
        Method = method,
        cutoff_quantile =
          0.75,
        n = length(x),
        BothHigh =
          obs$both,
        XHighOnly =
          obs$x_only,
        YHighOnly =
          obs$y_only,
        NeitherHigh =
          obs$neither,
        Null =
          null_name,
        PermutationStatistic =
          "BothHigh_count",
        permutation_mean =
          mean(
            rr$perm_both
          ),
        permutation_sd =
          stats::sd(
            rr$perm_both
          ),
        p_overlap_greater =
          rr$p_greater,
        n_perm =
          N_PERM,
        seed =
          seed,
        FisherStatus =
          "DESCRIPTIVE_EFFECT_SIZE_ONLY",
        stringsAsFactors = FALSE
      ),
      or
    )

    dual_perm_null_rows[[length(dual_perm_null_rows) + 1L]] <- data.frame(
      Method = method,
      Null =
        null_name,
      permutation =
        seq_len(
          N_PERM
        ),
      BothHigh_perm =
        rr$perm_both,
      seed =
        seed,
      stringsAsFactors = FALSE
    )
  }
}

dual_perm <-
  do.call(
    rbind,
    dual_perm_rows
  )

dual_perm_null <-
  do.call(
    rbind,
    dual_perm_null_rows
  )

safe_csv(
  dual_perm,
  file.path(
    table_dir,
    "21_top25_dual_high_overlap_permutation_summary.csv"
  )
)

safe_csv(
  dual_perm_null,
  file.path(
    table_dir,
    "21_top25_dual_high_overlap_permutation_nulls.csv"
  )
)

############################################################
## 19. Historical 18D observed-statistic replication audit
##     Important: 18D used direct FNN indices and is NOT affected by
##     the Step11 edge-vectorization artifact for observed bivariate I.
############################################################

historical_replication <- data.frame(
  HistoricalFileExists =
    file.exists(
      historical_18D_moran_file
    ),
  NHistoricalRows =
    NA_integer_,
  MaxAbsDiff_Ixy =
    NA_real_,
  MaxAbsDiff_Iyx =
    NA_real_,
  MaxAbsDiff_Isymmetric =
    NA_real_,
  ExactObserved_1e10 =
    NA,
  stringsAsFactors = FALSE
)

if (
  file.exists(
    historical_18D_moran_file
  )
) {
  hist <-
    utils::read.csv(
      historical_18D_moran_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

  can <-
    biv_observed[
      biv_observed$Method ==
      "Raw",
      ,
      drop = FALSE
    ]

  mm <-
    merge(
      hist,
      can,
      by = "k",
      suffixes =
        c(
          "_historical",
          "_21C"
        )
    )

  if (
    nrow(mm) > 0L
  ) {
    historical_replication$NHistoricalRows <-
      nrow(hist)

    historical_replication$MaxAbsDiff_Ixy <-
      max(
        abs(
          mm$I_xy_historical -
          mm$I_xy_21C
        ),
        na.rm = TRUE
      )

    historical_replication$MaxAbsDiff_Iyx <-
      max(
        abs(
          mm$I_yx_historical -
          mm$I_yx_21C
        ),
        na.rm = TRUE
      )

    historical_replication$MaxAbsDiff_Isymmetric <-
      max(
        abs(
          mm$I_symmetric_historical -
          mm$I_symmetric_21C
        ),
        na.rm = TRUE
      )

    historical_replication$ExactObserved_1e10 <-
      all(
        c(
          historical_replication$
            MaxAbsDiff_Ixy,
          historical_replication$
            MaxAbsDiff_Iyx,
          historical_replication$
            MaxAbsDiff_Isymmetric
        ) <=
        1e-10
      )
  }
}

safe_csv(
  historical_replication,
  file.path(
    table_dir,
    "21_historical_18D_observed_bivariate_replication_audit.csv"
  )
)

############################################################
## 20. Method contract
############################################################

method_contract <- data.frame(
  Analysis = c(
    "Univariate Moran raw/QC",
    "Primary-pair symmetric bivariate Moran",
    "Primary-pair unrestricted permutation",
    "Composition-stratified bivariate permutation",
    "Raw threshold-grid dual-high",
    "Raw/residualized top25 dual-high overlap permutation",
    "Fisher OR",
    "QC-adjusted partial Spearman"
  ),
  CanonicalDefinition = c(
    "Corrected direct FNN k=6; historical Moran formula retained; 999 two-sided global score permutations",
    "I_sym=(I_xy+I_yx)/2; I_xy=mean(zx_i * mean_k(zy_neighbors)); direct row-aligned FNN graph",
    "One permutation of y per iteration; recompute full symmetric statistic; upper-tail p primary; centered two-sided p sensitivity",
    "Permute y within joint rank-quartile CAF_Stromal x Endothelial strata; 999 permutations; NEW canonical sensitivity",
    "top15/top20/top25/top30/top35/median/z>0; Fisher OR descriptive",
    "Top25; permute y unrestricted or within CAF/endothelial strata; BothHigh count as permutation statistic",
    "Descriptive effect-size summary only because spots/edges are not treated as independent inferential replicates",
    "Rank all variables, residualize ranked x and y on ranked QC covariates, Pearson correlation of residuals"
  ),
  Status = c(
    "CORRECTED_CANONICAL",
    "CANONICAL_OBSERVED_FORMULA_FROM_18D",
    "CANONICAL_HARDENED_NULL",
    "NEW_CANONICAL_SENSITIVITY",
    "CANONICAL_DESCRIPTIVE_SENSITIVITY",
    "NEW_CANONICAL_PERMUTATION_SENSITIVITY",
    "DESCRIPTIVE_ONLY",
    "CANONICAL_FROM_18D"
  ),
  stringsAsFactors = FALSE
)

safe_csv(
  method_contract,
  file.path(
    table_dir,
    "21_method_contract.csv"
  )
)

############################################################
## 21. Gates
############################################################

gates <- do.call(
  rbind,
  list(
    gate_row(
      "Step11_percent_mt_canonicalized",
      "percent_mt" %in%
        colnames(meta),
      TRUE,
      "percent_mt" %in%
        colnames(meta)
    ),
    gate_row(
      "Canonical_spot_count",
      length(
        unique(
          meta_ids
        )
      ),
      EXPECTED_N,
      length(
        unique(
          meta_ids
        )
      ) ==
      EXPECTED_N
    ),
    gate_row(
      "QC_spot_set_exact",
      setequal(
        meta_ids,
        qc_ids
      ),
      TRUE,
      setequal(
        meta_ids,
        qc_ids
      )
    ),
    gate_row(
      "Composition_spot_set_exact",
      setequal(
        meta_ids,
        comp_ids
      ),
      TRUE,
      setequal(
        meta_ids,
        comp_ids
      )
    ),
    gate_row(
      "Step11_vs_18E_four_state_scores_exact",
      paste0(
        sum(
          state_equivalence$
            Exact_1e10
        ),
        "/4"
      ),
      "4/4",
      all(
        state_equivalence$
          Exact_1e10
      )
    ),
    gate_row(
      "Coordinate_rows_complete",
      sum(
        coord_complete
      ),
      EXPECTED_N,
      sum(
        coord_complete
      ) ==
      EXPECTED_N
    ),
    gate_row(
      "Coordinate_pairs_unique",
      length(
        unique(
          coord_key
        )
      ),
      EXPECTED_N,
      length(
        unique(
          coord_key
        )
      ) ==
      EXPECTED_N
    ),
    gate_row(
      "Primary_pair_common_complete_universe",
      common_n,
      EXPECTED_N,
      common_n ==
      EXPECTED_N
    ),
    gate_row(
      "All_corrected_graphs_structurally_clean",
      paste(
        graph_audit$
          exact_structure,
        collapse = ";"
      ),
      "all TRUE",
      all(
        graph_audit$
          exact_structure
      )
    ),
    gate_row(
      "Historical_k6_graph_not_equal_corrected",
      graph_pair_comparison$
        Jaccard,
      "<1",
      is.finite(
        graph_pair_comparison$
          Jaccard
      ) &&
      graph_pair_comparison$
        Jaccard <
      1,
      "DIAGNOSTIC"
    ),
    gate_row(
      "Historical_18D_observed_bivariate_replication",
      historical_replication$
        ExactObserved_1e10,
      TRUE,
      isTRUE(
        historical_replication$
          ExactObserved_1e10
      ),
      "DIAGNOSTIC"
    ),
    gate_row(
      "Composition_strata_nontrivial",
      strata_audit$
        n_strata,
      ">1",
      strata_audit$
        n_strata >
      1L
    ),
    gate_row(
      "Composition_strata_min_size",
      strata_audit$
        min_stratum_n,
      ">=10",
      strata_audit$
        min_stratum_n >=
      10L
    )
  )
)

safe_csv(
  gates,
  file.path(
    table_dir,
    "21_reproducibility_gates.csv"
  )
)

blocking_fail <-
  sum(
    !gates$Pass &
    gates$Classification ==
      "BLOCKING"
  )

diagnostic_fail <-
  sum(
    !gates$Pass &
    gates$Classification ==
      "DIAGNOSTIC"
  )

############################################################
## 22. Decision summary
############################################################

main_raw_k6 <-
  raw_k_perm[
    raw_k_perm$k ==
      MAIN_K,
    ,
    drop = FALSE
  ]

main_comp_biv <-
  composition_biv[
    composition_biv$k ==
      MAIN_K,
    ,
    drop = FALSE
  ]

decision <- data.frame(
  Item = c(
    "Stage1_lineage_status",
    "Historical_Step11_edge_graph_status",
    "Historical_18D_observed_bivariate_status",
    "Corrected_common_spots",
    "Corrected_k6_graph_clean",
    "Raw_k6_I_symmetric",
    "Raw_k6_unrestricted_p_greater",
    "Composition_aware_methods_tested",
    "New_composition_stratified_null",
    "Blocking_failures",
    "Diagnostic_failures",
    "21_status"
  ),
  Value = c(
    "FINAL_FREEZE",
    "LEGACY_INVALID_FOR_CANONICAL_GRAPH_INFERENCE",
    "DIRECT_FNN_OBSERVED_STATISTIC_CAN_BE_RETAINED_IF_REPLICATION_GATE_PASSES",
    common_n,
    graph_audit$
      exact_structure[
        graph_audit$k ==
          MAIN_K
      ][[1]],
    if (
      nrow(
        main_raw_k6
      ) == 1L
    ) {
      main_raw_k6$
        I_symmetric[[1]]
    } else {
      NA
    },
    if (
      nrow(
        main_raw_k6
      ) == 1L
    ) {
      main_raw_k6$
        p_greater[[1]]
    } else {
      NA
    },
    paste(
      names(
        score_pairs
      ),
      collapse = "; "
    ),
    "CAF_Stromal x Endothelial joint rank quartiles; within-stratum y permutation",
    blocking_fail,
    diagnostic_fail,
    if (
      blocking_fail == 0L
    ) {
      "STAGE2_CORRECTED_REBUILD_PASS_PENDING_NUMERICAL_INTERPRETATION"
    } else {
      "STAGE2_CORRECTED_REBUILD_BLOCKED"
    }
  ),
  stringsAsFactors = FALSE
)

safe_csv(
  decision,
  file.path(
    table_dir,
    "21_decision_summary.csv"
  )
)


############################################################
## 23. Canonical S15-S19 panel-source tables
##
## IMPORTANT
##   These are reporting DATA CONTRACTS only. No figure is drawn here.
##   S15-S18 are fully determined by the canonical primary Visium analysis.
##   S19 additionally requires a unique spot-level humoral/TLS-like module
##   table with exact 3,458-spot identity.
############################################################

derive_dual_category <- function(
  x,
  y,
  q = 0.75
) {
  x_cut <- as.numeric(
    stats::quantile(
      x,
      probs = q,
      na.rm = TRUE,
      names = FALSE
    )
  )

  y_cut <- as.numeric(
    stats::quantile(
      y,
      probs = q,
      na.rm = TRUE,
      names = FALSE
    )
  )

  xh <- x >= x_cut
  yh <- y >= y_cut

  code <- ifelse(
    xh & yh,
    "BothHigh",
    ifelse(
      xh & !yh,
      "MyeloidTregOnly",
      ifelse(
        !xh & yh,
        "TumorStromalOnly",
        "NeitherHigh"
      )
    )
  )

  display <- c(
    NeitherHigh = "neither high",
    MyeloidTregOnly = "myeloid-Treg high only",
    TumorStromalOnly = "tumor-dedifferentiation/stromal-remodeling high only",
    BothHigh = "both high"
  )

  list(
    code = code,
    display = unname(
      display[code]
    ),
    x_cut = x_cut,
    y_cut = y_cut
  )
}

display_x_col <- if (
  "imagecol" %in%
  colnames(common_meta)
) {
  "imagecol"
} else {
  "spatial_x"
}

display_y_col <- if (
  "imagerow" %in%
  colnames(common_meta)
) {
  "imagerow"
} else {
  "spatial_y"
}

spot_id_common <- meta_ids[
  common_idx
]

raw_category <- derive_dual_category(
  as.numeric(
    common_comp[[score_pairs$Raw[["x"]]]]
  ),
  as.numeric(
    common_comp[[score_pairs$Raw[["y"]]]]
  ),
  q = 0.75
)

combined_category <- derive_dual_category(
  as.numeric(
    common_comp[[score_pairs$QC_Composition_residualized[["x"]]]]
  ),
  as.numeric(
    common_comp[[score_pairs$QC_Composition_residualized[["y"]]]]
  ),
  q = 0.75
)

state_display <- c(
  Immune_defective_Cold = "immune-defective/cold",
  Myeloid_Treg_Immunosuppressive = "myeloid-Treg immunosuppressive",
  Tumor_dedifferentiation_Stromal_remodeling =
    "tumor-dedifferentiation/stromal-remodeling",
  Melanocytic_Differentiation = "melanocytic differentiation"
)

## -------------------------
## S15 Panel A
## -------------------------

state_matrix <- as.matrix(
  common_meta[
    ,
    STATE_COLS,
    drop = FALSE
  ]
)

storage.mode(
  state_matrix
) <- "numeric"

dominant_idx <- max.col(
  state_matrix,
  ties.method = "first"
)

s15A <- data.frame(
  Spot = spot_id_common,
  DisplayX =
    as.numeric(
      common_meta[[display_x_col]]
    ),
  DisplayY =
    as.numeric(
      common_meta[[display_y_col]]
    ),
  DisplayYInvert =
    TRUE,
  DominantStateCode =
    STATE_COLS[
      dominant_idx
    ],
  DominantState =
    unname(
      state_display[
        STATE_COLS[
          dominant_idx
        ]
      ]
    ),
  stringsAsFactors = FALSE
)

for (
  st in STATE_COLS
) {
  s15A[[st]] <- as.numeric(
    common_meta[[st]]
  )
}

safe_csv(
  s15A,
  file.path(
    panel_root,
    "S15",
    "21_S15A_dominant_state_spatial_map.csv"
  )
)

## S15 Panel B = corrected raw/QC univariate Moran.
safe_csv(
  univariate_moran,
  file.path(
    panel_root,
    "S15",
    "21_S15B_corrected_univariate_Moran_raw_QC_k6.csv"
  )
)

## S15 Panel C = descriptive edge-level OR only.
safe_csv(
  neighbor_enrichment,
  file.path(
    panel_root,
    "S15",
    "21_S15C_corrected_high_state_neighbor_OR_DESCRIPTIVE.csv"
  )
)

## -------------------------
## S16
## -------------------------

s16A <- data.frame(
  Spot = spot_id_common,
  MyeloidTreg =
    as.numeric(
      common_comp[[score_pairs$Raw[["x"]]]]
    ),
  TumorStromal =
    as.numeric(
      common_comp[[score_pairs$Raw[["y"]]]]
    ),
  SpatialX =
    as.numeric(
      common_meta$spatial_x
    ),
  SpatialY =
    as.numeric(
      common_meta$spatial_y
    ),
  stringsAsFactors = FALSE
)

safe_csv(
  s16A,
  file.path(
    panel_root,
    "S16",
    "21_S16A_continuous_primary_pair_spot_scores.csv"
  )
)

safe_csv(
  continuous_summary,
  file.path(
    panel_root,
    "S16",
    "21_S16A_continuous_primary_pair_summary.csv"
  )
)

s16B <- merge(
  biv_observed[
    biv_observed$Method ==
      "Raw",
    ,
    drop = FALSE
  ],
  raw_k_perm[
    ,
    c(
      "k",
      "perm_mean",
      "perm_sd",
      "p_greater",
      "p_two_centered",
      "n_perm",
      "seed"
    )
  ],
  by = "k",
  all.x = TRUE,
  sort = TRUE
)

safe_csv(
  s16B,
  file.path(
    panel_root,
    "S16",
    "21_S16B_symmetric_bivariate_Moran_k4_k6_k8_k12.csv"
  )
)

safe_csv(
  threshold_table,
  file.path(
    panel_root,
    "S16",
    "21_S16C_threshold_grid_dual_high_OR_DESCRIPTIVE.csv"
  )
)

## -------------------------
## S17
## -------------------------

s17A_rows <- list()

for (
  cc in COMPOSITION_COLS
) {
  s17A_rows[[length(s17A_rows) + 1L]] <- data.frame(
    Spot = spot_id_common,
    CategoryCode =
      raw_category$code,
    Category =
      raw_category$display,
    Module = cc,
    Score =
      as.numeric(
        common_comp[[cc]]
      ),
    stringsAsFactors = FALSE
  )
}

s17A <- do.call(
  rbind,
  s17A_rows
)

safe_csv(
  s17A,
  file.path(
    panel_root,
    "S17",
    "21_S17A_composition_modules_by_raw_category_long.csv"
  )
)

s17B_rows <- list()

for (
  st in STATE_COLS
) {
  for (
    cc in COMPOSITION_COLS
  ) {
    ct <- suppressWarnings(
      stats::cor.test(
        as.numeric(
          common_comp[[st]]
        ),
        as.numeric(
          common_comp[[cc]]
        ),
        method = "spearman",
        exact = FALSE
      )
    )

    s17B_rows[[length(s17B_rows) + 1L]] <- data.frame(
      StateCode = st,
      State =
        unname(
          state_display[[st]]
        ),
      Module = cc,
      rho =
        unname(
          ct$estimate
        ),
      p_nominal =
        ct$p.value,
      n =
        common_n,
      InferentialStatus =
        "DESCRIPTIVE_SPOT_LEVEL_ASSOCIATION",
      stringsAsFactors = FALSE
    )
  }
}

s17B <- do.call(
  rbind,
  s17B_rows
)

s17B$FDR_BH_descriptive <-
  p.adjust(
    s17B$p_nominal,
    method = "BH"
  )

safe_csv(
  s17B,
  file.path(
    panel_root,
    "S17",
    "21_S17B_state_composition_Spearman.csv"
  )
)

s17C <- dual_perm[
  !duplicated(
    dual_perm$Method
  ),
  c(
    "Method",
    "cutoff_quantile",
    "n",
    "BothHigh",
    "XHighOnly",
    "YHighOnly",
    "NeitherHigh",
    "OR",
    "CI_low",
    "CI_high",
    "p_value",
    "FisherStatus"
  ),
  drop = FALSE
]

safe_csv(
  s17C,
  file.path(
    panel_root,
    "S17",
    "21_S17C_raw_residualized_top25_Fisher_OR_DESCRIPTIVE.csv"
  )
)

s17D <- rbind(
  data.frame(
    Spot = spot_id_common,
    DisplayX =
      as.numeric(
        common_meta[[display_x_col]]
      ),
    DisplayY =
      as.numeric(
        common_meta[[display_y_col]]
      ),
    Method =
      "Raw",
    CategoryCode =
      raw_category$code,
    Category =
      raw_category$display,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Spot = spot_id_common,
    DisplayX =
      as.numeric(
        common_meta[[display_x_col]]
      ),
    DisplayY =
      as.numeric(
        common_meta[[display_y_col]]
      ),
    Method =
      "QC_Composition_residualized",
    CategoryCode =
      combined_category$code,
    Category =
      combined_category$display,
    stringsAsFactors = FALSE
  )
)

safe_csv(
  s17D,
  file.path(
    panel_root,
    "S17",
    "21_S17D_raw_vs_QCcomposition_spatial_categories.csv"
  )
)

## -------------------------
## S18
## -------------------------

safe_csv(
  s17B,
  file.path(
    panel_root,
    "S18",
    "21_S18A_state_composition_Spearman.csv"
  )
)

s18B_rows <- list()

for (
  method_i in seq_along(
    score_pairs
  )
) {
  method <-
    names(
      score_pairs
    )[[method_i]]

  pair <-
    score_pairs[[method_i]]

  x <- as.numeric(
    common_comp[[pair[["x"]]]]
  )

  y <- as.numeric(
    common_comp[[pair[["y"]]]]
  )

  ct <- suppressWarnings(
    stats::cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )

  s18B_rows[[length(s18B_rows) + 1L]] <- data.frame(
    Method = method,
    n =
      length(x),
    rho =
      unname(
        ct$estimate
      ),
    p_nominal =
      ct$p.value,
    StatisticalRole =
      "CONTINUOUS_EFFECT_SIZE_SUMMARY",
    stringsAsFactors = FALSE
  )
}

s18B <- do.call(
  rbind,
  s18B_rows
)

safe_csv(
  s18B,
  file.path(
    panel_root,
    "S18",
    "21_S18B_raw_residualized_primary_pair_Spearman.csv"
  )
)

safe_csv(
  dual_perm,
  file.path(
    panel_root,
    "S18",
    "21_S18C_top25_dual_high_OR_and_permutation.csv"
  )
)

safe_csv(
  composition_biv,
  file.path(
    panel_root,
    "S18",
    "21_S18D_bivariate_Moran_and_permutation_k6.csv"
  )
)

## -------------------------
## S19 dependency resolver
## -------------------------

tls_module_aliases <- list(
  B_Cell =
    c(
      "B_Cell",
      "B-cell",
      "B_cell"
    ),
  Plasma_Cell =
    c(
      "Plasma_Cell",
      "plasma cell",
      "Plasma_cell"
    ),
  T_Cell_Zone =
    c(
      "T_Cell_Zone",
      "T-cell zone",
      "T_cell_zone"
    ),
  Tfh_Like =
    c(
      "Tfh_Like",
      "Tfh-like",
      "Tfh_like"
    ),
  FDC_HEV_Like =
    c(
      "FDC_HEV_Like",
      "FDC/HEV-like",
      "FDC_HEV_like"
    ),
  IG_Humoral =
    c(
      "IG_Humoral",
      "Ig-humoral",
      "IG_humoral"
    ),
  TLS_Chemokine_Core =
    c(
      "TLS_Chemokine_Core",
      "TLS chemokine core",
      "TLS_chemokine_core"
    )
)

resolve_tls_module_columns <- function(df) {
  mapping <- character(
    length(
      tls_module_aliases
    )
  )

  names(mapping) <-
    names(
      tls_module_aliases
    )

  for (
    nm in names(
      tls_module_aliases
    )
  ) {
    hits <- intersect(
      tls_module_aliases[[nm]],
      colnames(df)
    )

    if (
      length(hits) == 0L
    ) {
      return(NULL)
    }

    mapping[[nm]] <-
      hits[[1]]
  }

  mapping
}

read_tls_candidate <- function(file) {
  ext <- tolower(
    tools::file_ext(file)
  )

  out <- tryCatch(
    {
      if (
        ext == "csv"
      ) {
        utils::read.csv(
          file,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      } else if (
        ext %in%
        c(
          "tsv",
          "txt"
        )
      ) {
        utils::read.delim(
          file,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      } else if (
        ext == "rds"
      ) {
        as_df(
          readRDS(file),
          basename(file)
        )
      } else {
        NULL
      }
    },
    error = function(e) {
      NULL
    }
  )

  out
}

manual_tls_file <-
  Sys.getenv(
    "ICB_VISIUM_TLS_MODULE_FILE",
    unset = ""
  )

## Public sequential source-lineage lock for S19.
## Step 06 deliberately writes both a dedicated module-score table and a
## broader joined spot/state/composition/TLS table. Both contain the same
## seven TLS-like module columns and exact 3,458-spot lineage, so an
## unrestricted recursive scan can yield >1 valid path and leave the S19
## dependency formally unresolved. Prefer the dedicated Step-06 module-score
## table as the canonical S19 source. The environment-variable override is
## retained for an explicit user-supplied source; recursive discovery remains
## a fallback only when the preferred source is absent. No score, threshold,
## state definition, graph, permutation, or statistical calculation changes.
preferred_tls_file <- file.path(
  project_dir,
  "results",
  "tables",
  "revision_humoral_TLS_like",
  "18F_humoral_TLS_like_module_scores.csv"
)

tls_candidates <- character(0)

if (
  nzchar(
    manual_tls_file
  )
) {
  if (
    !file.exists(
      manual_tls_file
    )
  ) {
    stop(
      "ICB_VISIUM_TLS_MODULE_FILE was set but does not exist: ",
      manual_tls_file,
      call. = FALSE
    )
  }

  tls_candidates <-
    normalizePath(
      manual_tls_file,
      winslash = "/",
      mustWork = TRUE
    )
} else if (
  file.exists(
    preferred_tls_file
  )
) {
  tls_candidates <-
    normalizePath(
      preferred_tls_file,
      winslash = "/",
      mustWork = TRUE
    )
} else {
  all_table_candidates <- list.files(
    file.path(
      project_dir,
      "results",
      "tables"
    ),
    pattern =
      "\\.(csv|tsv|txt|rds)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  all_table_candidates <-
    all_table_candidates[
      !grepl(
        "primary_visium_spatial",
        all_table_candidates,
        fixed = TRUE
      )
    ]

  for (
    f in all_table_candidates
  ) {
    dd <- read_tls_candidate(f)

    if (
      is.null(dd)
    ) {
      next
    }

    mp <- resolve_tls_module_columns(dd)

    if (
      !is.null(mp)
    ) {
      tls_candidates <-
        c(
          tls_candidates,
          f
        )
    }
  }
}

tls_candidate_rows <- list()
tls_valid_paths <- character(0)
tls_loaded <- list()

for (
  f in unique(
    tls_candidates
  )
) {
  dd <- read_tls_candidate(f)
  mp <- if (
    is.null(dd)
  ) {
    NULL
  } else {
    resolve_tls_module_columns(dd)
  }

  id_exact <- FALSE
  id_source <- NA_character_
  n_rows <- if (
    is.null(dd)
  ) {
    NA_integer_
  } else {
    nrow(dd)
  }

  if (
    !is.null(dd) &&
    !is.null(mp)
  ) {
    id_info <- tryCatch(
      resolve_ids(dd),
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(id_info)
    ) {
      id_source <-
        id_info$source

      id_exact <-
        nrow(dd) ==
        EXPECTED_N &&
        setequal(
          id_info$ids,
          meta_ids
        )

      if (
        id_exact
      ) {
        tls_valid_paths <-
          c(
            tls_valid_paths,
            f
          )

        tls_loaded[[f]] <- list(
          data = dd,
          ids = id_info$ids,
          mapping = mp
        )
      }
    }
  }

  tls_candidate_rows[[length(tls_candidate_rows) + 1L]] <- data.frame(
    Path =
      normalizePath(
        f,
        winslash = "/",
        mustWork = FALSE
      ),
    MD5 =
      md5_file(f),
    N_rows =
      n_rows,
    ID_source =
      id_source,
    SpotSetExact =
      id_exact,
    ModuleMapping =
      if (
        is.null(mp)
      ) {
        NA_character_
      } else {
        paste(
          names(mp),
          mp,
          sep = "=",
          collapse = "; "
        )
      },
    stringsAsFactors = FALSE
  )
}

tls_candidate_inventory <- if (
  length(
    tls_candidate_rows
  ) > 0L
) {
  do.call(
    rbind,
    tls_candidate_rows
  )
} else {
  data.frame(
    Path = character(0),
    MD5 = character(0),
    N_rows = integer(0),
    ID_source = character(0),
    SpotSetExact = logical(0),
    ModuleMapping = character(0),
    stringsAsFactors = FALSE
  )
}

safe_csv(
  tls_candidate_inventory,
  file.path(
    panel_root,
    "S19",
    "21_S19_TLS_module_source_candidate_inventory.csv"
  )
)

tls_source_resolved <-
  length(
    unique(
      tls_valid_paths
    )
  ) ==
  1L

tls_source_file <-
  if (
    tls_source_resolved
  ) {
    unique(
      tls_valid_paths
    )[[1]]
  } else {
    NA_character_
  }

if (
  tls_source_resolved
) {
  tls_obj <-
    tls_loaded[[tls_source_file]]

  tls_df <-
    tls_obj$data

  tls_match <-
    match(
      meta_ids,
      tls_obj$ids
    )

  tls_df <-
    tls_df[
      tls_match[
        common_idx
      ],
      ,
      drop = FALSE
    ]

  tls_map <-
    tls_obj$mapping

  tls_display_map <- c(
    B_Cell = "B-cell",
    FDC_HEV_Like = "FDC/HEV-like",
    IG_Humoral = "Ig-humoral",
    Plasma_Cell = "plasma cell",
    T_Cell_Zone = "T-cell zone",
    Tfh_Like = "Tfh-like",
    TLS_Chemokine_Core = "TLS chemokine core"
  )

  tls_matrix <- data.frame(
    Spot = spot_id_common,
    stringsAsFactors = FALSE
  )

  for (
    nm in names(
      tls_display_map
    )
  ) {
    tls_matrix[[tls_display_map[[nm]]]] <- as.numeric(
      tls_df[[tls_map[[nm]]]]
    )
  }

  tls_display_order <- unname(
    tls_display_map
  )

  s19A_rows <- list()
  s19B_rows <- list()

  for (
    mm in tls_display_order
  ) {
    s19A_rows[[length(s19A_rows) + 1L]] <- data.frame(
      Spot = spot_id_common,
      CategoryCode =
        raw_category$code,
      Category =
        raw_category$display,
      Module = mm,
      Score =
        tls_matrix[[mm]],
      stringsAsFactors = FALSE
    )

    s19B_rows[[length(s19B_rows) + 1L]] <- data.frame(
      Spot = spot_id_common,
      SpatialX =
        as.numeric(
          common_meta$spatial_x
        ),
      SpatialY =
        as.numeric(
          common_meta$spatial_y
        ),
      Module = mm,
      Score =
        tls_matrix[[mm]],
      stringsAsFactors = FALSE
    )
  }

  s19A <- do.call(
    rbind,
    s19A_rows
  )

  s19B <- do.call(
    rbind,
    s19B_rows
  )

  safe_csv(
    s19A,
    file.path(
      panel_root,
      "S19",
      "21_S19A_humoral_TLS_modules_by_raw_category_long.csv"
    )
  )

  safe_csv(
    s19B,
    file.path(
      panel_root,
      "S19",
      "21_S19B_humoral_TLS_module_spatial_scores_long.csv"
    )
  )

  module_autocorr_from_index <- function(
    x,
    idx
  ) {
    xz <- zvec(x)

    mean(
      xz *
      neighbor_mean(
        xz,
        idx
      )
    )
  }

  s19C_rows <- list()

  for (
    k in K_VALUES
  ) {
    idx <-
      graph_index[[as.character(k)]]

    for (
      mm in tls_display_order
    ) {
      s19C_rows[[length(s19C_rows) + 1L]] <- data.frame(
        Module = mm,
        k = k,
        n =
          common_n,
        SpatialAutocorrelation =
          module_autocorr_from_index(
            tls_matrix[[mm]],
            idx
          ),
        StatisticalRole =
          "DESCRIPTIVE_DIRECT_FNN_AUTOCORRELATION",
        stringsAsFactors = FALSE
      )
    }
  }

  s19C <- do.call(
    rbind,
    s19C_rows
  )

  safe_csv(
    s19C,
    file.path(
      panel_root,
      "S19",
      "21_S19C_humoral_TLS_module_spatial_autocorrelation_k_sensitivity.csv"
    )
  )

  q_tls <- 0.75

  cutoff_tls <- vapply(
    tls_display_order,
    function(mm) {
      as.numeric(
        stats::quantile(
          tls_matrix[[mm]],
          probs = q_tls,
          na.rm = TRUE,
          names = FALSE
        )
      )
    },
    numeric(1)
  )

  b_high <-
    tls_matrix[["B-cell"]] >=
    cutoff_tls[["B-cell"]]

  t_high <-
    tls_matrix[["T-cell zone"]] >=
    cutoff_tls[["T-cell zone"]]

  ig_high <-
    tls_matrix[["Ig-humoral"]] >=
    cutoff_tls[["Ig-humoral"]]

  tls_high <-
    tls_matrix[["TLS chemokine core"]] >=
    cutoff_tls[["TLS chemokine core"]]

  feature <- ifelse(
    tls_high &
    b_high &
    t_high,
    "TLS-chemokine/B-cell/T-cell co-high",
    ifelse(
      ig_high &
      t_high,
      "Ig-humoral/T-cell co-high",
      ifelse(
        tls_high,
        "TLS chemokine high/other",
        ifelse(
          ig_high,
          "Ig-humoral high/other",
          "other"
        )
      )
    )
  )

  s19D <- data.frame(
    Spot = spot_id_common,
    SpatialX =
      as.numeric(
        common_meta$spatial_x
      ),
    SpatialY =
      as.numeric(
        common_meta$spatial_y
      ),
    Feature =
      feature,
    HighQuantile =
      q_tls,
    InterpretiveStatus =
      "OPERATIONAL_TRANSCRIPTIONAL_CONTEXT_NOT_MATURE_TLS_CALLING",
    stringsAsFactors = FALSE
  )

  safe_csv(
    s19D,
    file.path(
      panel_root,
      "S19",
      "21_S19D_operational_humoral_TLS_cohigh_features.csv"
    )
  )
}

############################################################
## 24. S15-S19 rebuild contract and reporting gates
############################################################

rebuild_contract <- data.frame(
  Figure = c(
    "S15","S15","S15",
    "S16","S16","S16",
    "S17","S17","S17","S17",
    "S18","S18","S18","S18",
    "S19","S19","S19","S19"
  ),
  Panel = c(
    "A","B","C",
    "A","B","C",
    "A","B","C","D",
    "A","B","C","D",
    "A","B","C","D"
  ),
  CanonicalSource = c(
    "21_S15A_dominant_state_spatial_map.csv",
    "21_S15B_corrected_univariate_Moran_raw_QC_k6.csv",
    "21_S15C_corrected_high_state_neighbor_OR_DESCRIPTIVE.csv",
    "21_S16A_continuous_primary_pair_spot_scores.csv + summary",
    "21_S16B_symmetric_bivariate_Moran_k4_k6_k8_k12.csv",
    "21_S16C_threshold_grid_dual_high_OR_DESCRIPTIVE.csv",
    "21_S17A_composition_modules_by_raw_category_long.csv",
    "21_S17B_state_composition_Spearman.csv",
    "21_S17C_raw_residualized_top25_Fisher_OR_DESCRIPTIVE.csv",
    "21_S17D_raw_vs_QCcomposition_spatial_categories.csv",
    "21_S18A_state_composition_Spearman.csv",
    "21_S18B_raw_residualized_primary_pair_Spearman.csv",
    "21_S18C_top25_dual_high_OR_and_permutation.csv",
    "21_S18D_bivariate_Moran_and_permutation_k6.csv",
    "21_S19A_humoral_TLS_modules_by_raw_category_long.csv",
    "21_S19B_humoral_TLS_module_spatial_scores_long.csv",
    "21_S19C_humoral_TLS_module_spatial_autocorrelation_k_sensitivity.csv",
    "21_S19D_operational_humoral_TLS_cohigh_features.csv"
  ),
  RebuildRule = c(
    "REBUILD_FROM_CANONICAL_RAW_STATE_SCORES",
    "REPLACE_LEGACY_EDGE_MORAN_WITH_CORRECTED_FNN",
    "REPLACE_LEGACY_EDGE_OR; DISPLAY_AS_DESCRIPTIVE",
    "REBUILD_FROM_FROZEN_3458_SPOT_SCORES",
    "DIRECT_FNN_OBSERVED + HARDENED_PERMUTATION",
    "REBUILD_THRESHOLD_GRID; FISHER_DESCRIPTIVE",
    "REBUILD_FROM_18E_FROZEN_COMPOSITION_MODULES",
    "REBUILD_DESCRIPTIVE_STATE_MODULE_CORRELATIONS",
    "REBUILD_OR_ACROSS_RAW/QC/COMP/QC+COMP; DESCRIPTIVE",
    "REBUILD_RAW_VS_QC+COMP_OPERATIONAL_MAP",
    "SAME_CANONICAL_CORRELATION_SOURCE_AS_S17B",
    "REBUILD_CONTINUOUS_EFFECT_SIZE_ACROSS_METHODS",
    "REBUILD_WITH_UNRESTRICTED_AND_COMPOSITION_STRATIFIED_NULLS",
    "REBUILD_WITH_CORRECTED_DIRECT_FNN + HARDENED_NULLS",
    "REBUILD_IF_UNIQUE_3458_SPOT_TLS_SOURCE_RESOLVED",
    "REBUILD_IF_UNIQUE_3458_SPOT_TLS_SOURCE_RESOLVED",
    "DIRECT_FNN_DESCRIPTIVE_AUTOCORRELATION; NO_UNSUPPORTED_PERMUTATION_CLAIM",
    "OPERATIONAL_TOP25_FEATURES; NOT_MATURE_TLS_CALLING"
  ),
  InferentialRole = c(
    "DESCRIPTIVE_MAP",
    "PERMUTATION_SUPPORTED_SPATIAL_AUTOCORRELATION",
    "DESCRIPTIVE_EDGE_EFFECT_SIZE",
    "CONTINUOUS_EFFECT_SIZE",
    "PERMUTATION_SUPPORTED_SPATIAL_COVARIATION",
    "DESCRIPTIVE_THRESHOLD_SENSITIVITY",
    "DESCRIPTIVE_COMPOSITION_CONTEXT",
    "DESCRIPTIVE_SPOT_LEVEL_ASSOCIATION",
    "DESCRIPTIVE_EFFECT_SIZE",
    "DESCRIPTIVE_MAP",
    "DESCRIPTIVE_SPOT_LEVEL_ASSOCIATION",
    "CONTINUOUS_EFFECT_SIZE",
    "PERMUTATION_SUPPORTED_OVERLAP_SENSITIVITY",
    "PERMUTATION_SUPPORTED_SPATIAL_COVARIATION",
    "DESCRIPTIVE_HUMORAL_TLS_LIKE_CONTEXT",
    "DESCRIPTIVE_SPATIAL_MAP",
    "DESCRIPTIVE_SPATIAL_AUTOCORRELATION",
    "OPERATIONAL_CONTEXT_MAP"
  ),
  HistoricalStatus = c(
    "RETAIN_CONCEPT_REBUILD_SOURCE",
    "LEGACY_GRAPH_INVALID",
    "LEGACY_GRAPH_INVALID",
    "RETAIN_CONCEPT_REBUILD_SOURCE",
    "HISTORICAL_OBSERVED_DIRECT_FNN_VALID; USE_CANONICAL_SOURCE",
    "RETAIN_CONCEPT_REBUILD_SOURCE",
    "18E_COMPOSITION_ANALYSIS_RETAINED",
    "18E_COMPOSITION_ANALYSIS_RETAINED",
    "18E_COMPOSITION_ANALYSIS_RETAINED",
    "18E_COMPOSITION_ANALYSIS_RETAINED",
    "RETAIN_CONCEPT_REBUILD_SOURCE",
    "RETAIN_CONCEPT_REBUILD_SOURCE",
    "HISTORICAL_PERMUTATION_WORDING_HARDENED",
    "HISTORICAL_NULL_HARDENED",
    "DEPENDENCY_REQUIRES_EXACT_SOURCE_LINEAGE",
    "DEPENDENCY_REQUIRES_EXACT_SOURCE_LINEAGE",
    "REPORTING_DIRECT_FNN_VALID; REBUILD_CANONICAL",
    "RETAIN_OPERATIONAL_ONLY"
  ),
  stringsAsFactors = FALSE
)

safe_csv(
  rebuild_contract,
  file.path(
    canonical_reporting_root,
    "21D_S15_S19_rebuild_contract.csv"
  )
)

reporting_gates <- data.frame(
  Gate = c(
    "S15_panel_sources_complete",
    "S16_panel_sources_complete",
    "S17_panel_sources_complete",
    "S18_panel_sources_complete",
    "S19_TLS_source_unique_exact_3458_spots",
    "S19_panel_sources_complete_if_resolved"
  ),
  Pass = c(
    all(
      file.exists(
        file.path(
          panel_root,
          "S15",
          c(
            "21_S15A_dominant_state_spatial_map.csv",
            "21_S15B_corrected_univariate_Moran_raw_QC_k6.csv",
            "21_S15C_corrected_high_state_neighbor_OR_DESCRIPTIVE.csv"
          )
        )
      )
    ),
    all(
      file.exists(
        file.path(
          panel_root,
          "S16",
          c(
            "21_S16A_continuous_primary_pair_spot_scores.csv",
            "21_S16A_continuous_primary_pair_summary.csv",
            "21_S16B_symmetric_bivariate_Moran_k4_k6_k8_k12.csv",
            "21_S16C_threshold_grid_dual_high_OR_DESCRIPTIVE.csv"
          )
        )
      )
    ),
    all(
      file.exists(
        file.path(
          panel_root,
          "S17",
          c(
            "21_S17A_composition_modules_by_raw_category_long.csv",
            "21_S17B_state_composition_Spearman.csv",
            "21_S17C_raw_residualized_top25_Fisher_OR_DESCRIPTIVE.csv",
            "21_S17D_raw_vs_QCcomposition_spatial_categories.csv"
          )
        )
      )
    ),
    all(
      file.exists(
        file.path(
          panel_root,
          "S18",
          c(
            "21_S18A_state_composition_Spearman.csv",
            "21_S18B_raw_residualized_primary_pair_Spearman.csv",
            "21_S18C_top25_dual_high_OR_and_permutation.csv",
            "21_S18D_bivariate_Moran_and_permutation_k6.csv"
          )
        )
      )
    ),
    tls_source_resolved,
    if (
      tls_source_resolved
    ) {
      all(
        file.exists(
          file.path(
            panel_root,
            "S19",
            c(
              "21_S19A_humoral_TLS_modules_by_raw_category_long.csv",
              "21_S19B_humoral_TLS_module_spatial_scores_long.csv",
              "21_S19C_humoral_TLS_module_spatial_autocorrelation_k_sensitivity.csv",
              "21_S19D_operational_humoral_TLS_cohigh_features.csv"
            )
          )
        )
      )
    } else {
      FALSE
    }
  ),
  Classification = c(
    "BLOCKING",
    "BLOCKING",
    "BLOCKING",
    "BLOCKING",
    "S19_DEPENDENCY",
    "S19_DEPENDENCY"
  ),
  stringsAsFactors = FALSE
)

safe_csv(
  reporting_gates,
  file.path(
    canonical_reporting_root,
    "21D_S15_S19_reporting_source_gates.csv"
  )
)

s15_s18_ready <-
  all(
    reporting_gates$Pass[
      reporting_gates$Gate %in%
      c(
        "S15_panel_sources_complete",
        "S16_panel_sources_complete",
        "S17_panel_sources_complete",
        "S18_panel_sources_complete"
      )
    ]
  )

s19_ready <-
  all(
    reporting_gates$Pass[
      reporting_gates$Classification ==
        "S19_DEPENDENCY"
    ]
  )

reporting_decision <- data.frame(
  Item = c(
    "CanonicalAnalysisBlockingFailures",
    "S15_S18_reporting_sources_ready",
    "S19_TLS_source_resolved",
    "S19_reporting_sources_ready",
    "ReportingRebuildStatus"
  ),
  Value = c(
    blocking_fail,
    s15_s18_ready,
    tls_source_resolved,
    s19_ready,
    if (
      blocking_fail == 0L &&
      s15_s18_ready &&
      s19_ready
    ) {
      "S15_S19_PANEL_SOURCES_READY_FOR_REPORTING_REBUILD"
    } else if (
      blocking_fail == 0L &&
      s15_s18_ready
    ) {
      "S15_S18_READY__S19_SOURCE_LINEAGE_REQUIRES_CLOSURE"
    } else {
      "REPORTING_SOURCE_BUILD_BLOCKED"
    }
  ),
  stringsAsFactors = FALSE
)

safe_csv(
  reporting_decision,
  file.path(
    canonical_reporting_root,
    "21D_reporting_rebuild_decision.csv"
  )
)


############################################################
## 25. Provenance
############################################################

input_manifest <- data.frame(
  Input = c(
    "Step11 metadata RDS",
    "Step11 QC-residual RDS",
    "Historical malformed kNN6 edge RDS",
    "18E composition/residualized score table",
    "Historical 18D Moran table"
  ),
  Path = c(
    step11_meta_rds,
    step11_qc_rds,
    legacy_edge_rds,
    composition_file,
    historical_18D_moran_file
  ),
  Exists = file.exists(
    c(
      step11_meta_rds,
      step11_qc_rds,
      legacy_edge_rds,
      composition_file,
      historical_18D_moran_file
    )
  ),
  MD5 = vapply(
    c(
      step11_meta_rds,
      step11_qc_rds,
      legacy_edge_rds,
      composition_file,
      historical_18D_moran_file
    ),
    md5_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

safe_csv(
  input_manifest,
  file.path(
    table_dir,
    "21_input_manifest.csv"
  )
)

session_file <- file.path(
  log_dir,
  "sessionInfo_21_Primary_Visium_spatial_analysis.txt"
)

sink(
  session_file
)
print(
  sessionInfo()
)
sink()

manifest_file <- file.path(
  table_dir,
  "21_output_MD5_manifest.csv"
)

output_files <- c(
  list.files(
    table_dir,
    full.names = TRUE
  ),
  list.files(
    intermediate_dir,
    full.names = TRUE
  ),
  session_file
)

output_files <-
  output_files[
    normalizePath(
      output_files,
      winslash = "/",
      mustWork = FALSE
    ) !=
    normalizePath(
      manifest_file,
      winslash = "/",
      mustWork = FALSE
    )
  ]

output_manifest <- data.frame(
  File =
    basename(
      output_files
    ),
  Path =
    output_files,
  SizeBytes =
    file.info(
      output_files
    )$size,
  MD5 =
    vapply(
      output_files,
      md5_file,
      character(1)
    ),
  stringsAsFactors = FALSE
)

safe_csv(
  output_manifest,
  manifest_file
)

############################################################
## 26. Console summary
############################################################

cat("\n============================================================\n")
cat("21D PRIMARY VISIUM CANONICAL MOTHER SCRIPT COMPLETED\n")
cat("============================================================\n")
cat("Canonical spot universe: ", length(unique(meta_ids)), "\n", sep = "")
cat("Common complete Stage-2 spots: ", common_n, "\n", sep = "")
cat("Corrected graphs clean: ", all(graph_audit$exact_structure), "\n", sep = "")
cat(
  "Historical k6 vs corrected Jaccard: ",
  sprintf(
    "%.8f",
    graph_pair_comparison$Jaccard
  ),
  "\n",
  sep = ""
)
cat(
  "Historical 18D observed bivariate replication: ",
  historical_replication$ExactObserved_1e10,
  "\n",
  sep = ""
)

if (
  nrow(
    main_raw_k6
  ) == 1L
) {
  cat(
    "Raw k6 symmetric bivariate I: ",
    sprintf(
      "%.8f",
      main_raw_k6$I_symmetric[[1]]
    ),
    "\n",
    sep = ""
  )

  cat(
    "Raw k6 unrestricted permutation P(greater): ",
    sprintf(
      "%.6f",
      main_raw_k6$p_greater[[1]]
    ),
    "\n",
    sep = ""
  )
}

cat(
  "Composition strata: ",
  strata_audit$n_strata,
  " (min n=",
  strata_audit$min_stratum_n,
  ")\n",
  sep = ""
)
cat("Blocking canonical-analysis gates failed: ", blocking_fail, "\n", sep = "")
cat("Diagnostic canonical-analysis gates failed: ", diagnostic_fail, "\n", sep = "")
cat(
  "Canonical analysis status: ",
  if (
    blocking_fail == 0L
  ) {
    "STAGE2_CORRECTED_REBUILD_PASS_PENDING_NUMERICAL_INTERPRETATION"
  } else {
    "STAGE2_CORRECTED_REBUILD_BLOCKED"
  },
  "\n",
  sep = ""
)
cat(
  "S15-S18 panel sources ready: ",
  s15_s18_ready,
  "\n",
  sep = ""
)
cat(
  "S19 TLS source resolved: ",
  tls_source_resolved,
  "\n",
  sep = ""
)
cat(
  "S19 panel sources ready: ",
  s19_ready,
  "\n",
  sep = ""
)
cat(
  "Reporting rebuild status: ",
  reporting_decision$Value[
    reporting_decision$Item ==
      "ReportingRebuildStatus"
  ],
  "\n",
  sep = ""
)
cat(
  "Canonical tables: ",
  normalizePath(
    table_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)
cat(
  "Canonical reporting: ",
  normalizePath(
    canonical_reporting_root,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")

## Public sequential end gate
.step07_out <- file.path(panel_root, "S18", "21_S18D_bivariate_Moran_and_permutation_k6.csv")
.step07_gate <- file.path(canonical_reporting_root, "21D_reporting_rebuild_decision.csv")
if (!file.exists(.step07_out)) stop("STEP 07 panel source missing: ", .step07_out, call. = FALSE)
if (!file.exists(.step07_gate)) stop("STEP 07 reporting decision missing: ", .step07_gate, call. = FALSE)
cat("\nSTEP 07 PASS — canonical Primary Visium panel sources rebuilt\n")
