############################################################
## 05_bulk_response_sensitivity_analysis.R
##
## Purpose:
##   Public reproducibility script for the bulk-response sensitivity analyses
##   underlying Supplementary Figures S11-S13.
##
## Consolidates three historical revision modules:
##   - 18A: composition-aware audit (S12)
##   - 18C: GSE78220 fragility audit (S11)
##   - 18G: formal composition-adjusted response analysis (S13)
##
## Public-clean principles:
##   - Analysis only: no manuscript figure rendering in this script.
##   - No recursive input discovery and no local-machine path fallback.
##   - Canonical upstream outputs from scripts 01, 02, and 04 are required.
##   - Duplicate raw GSE78220 AUC point estimates are computed once and audited.
##   - S12 and S13 composition marker definitions remain separate because the
##     historical 18A and 18G analyses used different frozen marker panels.
##   - GSE91061 formal analysis is explicitly pretreatment-only, matching the
##     manuscript definition and script 04.
##   - FINAL legacy-CI provenance patch:
##       * submitted S11/18C and S12/18A bootstrap CIs are frozen verbatim for
##         exact submitted-figure/result reproduction;
##       * those historical pROC bootstrap CI calls used 5,000 stratified
##         resamples but did not record an RNG seed;
##       * deterministic seeded reanalysis CIs are retained in separate audit
##         tables and never silently substituted for the submitted CIs;
##       * all non-CI S11/S12 statistics and the entire S13/18G analysis remain
##         unchanged from the previously validated public-clean candidate.
##
## Interpretation:
##   These analyses are descriptive robustness / sensitivity analyses.
##   They are not predictive model validation or biomarker validation.
############################################################

options(stringsAsFactors = FALSE)

############################################################
## 0. Project root and output directories
############################################################

resolve_project_dir <- function() {
  env_dir <- Sys.getenv("ICB_PROJECT_DIR")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, winslash = "/", mustWork = FALSE))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  looks_like_root <- (
    dir.exists(file.path(wd, "data_processed")) &&
    dir.exists(file.path(wd, "results"))
  ) || file.exists(file.path(wd, "renv.lock"))

  if (looks_like_root) {
    return(wd)
  }

  stop(
    "Cannot determine project root. Run from the repository root or set ",
    'Sys.setenv(ICB_PROJECT_DIR = "D:/ICB_resistance_project").'
  )
}

project_dir <- resolve_project_dir()
message("Project directory: ", project_dir)

data_processed_dir <- file.path(project_dir, "data_processed")
results_dir <- file.path(project_dir, "results")
table_root <- file.path(results_dir, "tables")
intermediate_root <- file.path(results_dir, "intermediate")
out_dir <- file.path(table_root, "bulk_response_sensitivity")
out_intermediate_dir <- file.path(intermediate_root, "bulk_response_sensitivity")
log_dir <- file.path(project_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "stringr",
  "pROC",
  "tibble",
  "matrixStats"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Restore the repository environment with renv::restore()."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(pROC)
  library(tibble)
  library(matrixStats)
})

HAS_MCP <- requireNamespace("MCPcounter", quietly = TRUE)
HAS_ESTIMATE <- requireNamespace("estimate", quietly = TRUE)

############################################################
## 2. Frozen state definitions and random-number settings
############################################################

state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

state_labels <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" =
    "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

## S11 / historical 18C.
S11_BOOT_N <- 5000L
S11_PERM_N <- 10000L
S11_BOOT_SEED_BASE <- 20260608L
S11_PERM_SEED_BASE <- 20260680L
S11_MAX_AUC_PERM_SEED <- 20260701L

## Seeded script-02 reanalysis used only as an audit comparator for S11 raw CI.
## The submitted S11/18C CI itself is frozen below because historical 18C used
## pROC bootstrap n=5000 without a recorded RNG seed.
GSE78220_PRIMARY_AUC_SEED <- 20260506L
GSE78220_PRIMARY_AUC_BOOT_N <- 10000L

## S12 / historical 18A.
## Historical 18A used pROC bootstrap n=5000 without an explicit RNG seed.
## The submitted legacy CIs are frozen below. Deterministic seeded reanalysis
## CIs are still calculated and saved separately as a reproducibility audit.
## Point AUCs and residualized scores are unaffected by this provenance patch.
S12_RESIDUAL_CI_BOOT_N <- 5000L
S12_RESIDUAL_CI_SEED_BASE <- 20260608L

## S13 / historical 18G: preserve exactly.
S13_BOOT_N <- 10000L
S13_PERM_N <- 10000L
S13_SEED_BASE <- 20260622L

############################################################
## 3. Canonical upstream inputs
############################################################

input_files <- c(
  GSE244982_expression =
    file.path(data_processed_dir, "GSE244982_bulk_matrix.rds"),

  GSE244982_state_scores =
    file.path(
      table_root,
      "bulk_discovery",
      "GSE244982_final_tumor_immune_state_scores.csv"
    ),

  GSE78220_expression =
    file.path(
      data_processed_dir,
      "GSE78220_external_expression_matrix_for_scoring.rds"
    ),

  GSE78220_state_scores =
    file.path(
      table_root,
      "external_validation",
      "GSE78220_external_state_scores.csv"
    ),

  GSE78220_primary_auc =
    file.path(
      table_root,
      "external_validation",
      "GSE78220_nonresponse_AUC_with_95CI.csv"
    ),

  GSE91061_expression =
    file.path(
      intermediate_root,
      "GSE91061",
      "GSE91061_expression_symbol_matrix_for_state_scoring.rds"
    ),

  GSE91061_pretreatment_state_scores =
    file.path(
      table_root,
      "GSE91061",
      "GSE91061_state_scores_pretreatment.csv"
    )
)

############################################################
## 4. General helper functions
############################################################

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.frame(x), file)
  message("Saved: ", file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
  message("Saved: ", file)
}

relative_to_project <- function(path) {
  path_n <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_n <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root_n, "/")

  if (startsWith(path_n, prefix)) {
    substring(path_n, nchar(prefix) + 1L)
  } else {
    path_n
  }
}

safe_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

zvec <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\|.*$", "", x)
  x <- gsub("\\..*$", "", x)
  x <- gsub("\\s+", "", x)
  toupper(x)
}

make_sample_key <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub(
    "\\.BASELINE$|\\.PRE$|\\.PRETREATMENT$|\\.ON$|\\.ON_TREATMENT$",
    "",
    x,
    ignore.case = TRUE
  )
  gsub("[^A-Z0-9]", "", x)
}

collapse_duplicate_genes <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- clean_gene_symbol(rownames(mat))

  keep <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep, , drop = FALSE]

  if (!anyDuplicated(rownames(mat))) {
    return(mat)
  }

  dt <- data.table::as.data.table(mat)
  dt[, GeneSymbol := rownames(mat)]

  dt2 <- dt[
    ,
    lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
    by = GeneSymbol
  ]

  genes <- dt2$GeneSymbol
  mat2 <- as.matrix(dt2[, -1, with = FALSE])
  rownames(mat2) <- genes
  storage.mode(mat2) <- "numeric"
  mat2
}

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  row_mean <- rowMeans(mat, na.rm = TRUE)
  row_sd <- apply(mat, 1, stats::sd, na.rm = TRUE)

  z <- sweep(mat, 1, row_mean, "-")
  z <- sweep(z, 1, row_sd, "/")
  z[!is.finite(z)] <- NA_real_
  z
}

standardize_state_table <- function(df) {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)

  aliases <- list(
    Immune_defective_Cold = c(
      "Immune_defective_Cold",
      "Immune_defective_cold",
      "Immune_defective_Cold_score",
      "Immune_defective_cold_score"
    ),
    Myeloid_Treg_Immunosuppressive = c(
      "Myeloid_Treg_Immunosuppressive",
      "Myeloid_Treg_immunosuppressive",
      "Myeloid_Treg_Immunosuppressive_score",
      "Myeloid_Treg_immunosuppressive_score"
    ),
    Tumor_dedifferentiation_Stromal_remodeling = c(
      "Tumor_dedifferentiation_Stromal_remodeling",
      "Tumor_dedifferentiation_stromal_remodeling",
      "Tumor_dedifferentiation_Stromal_remodeling_score",
      "Tumor_dedifferentiated",
      "Tumor_dedifferentiated_score",
      "Dediff_Stromal",
      "Dediff_Stromal_score"
    ),
    Melanocytic_Differentiation = c(
      "Melanocytic_Differentiation",
      "Melanocytic_differentiation",
      "Melanocytic_Differentiation_score",
      "Melanocytic_differentiated",
      "Melanocytic_differentiated_score"
    )
  )

  for (canonical in names(aliases)) {
    if (canonical %in% colnames(df)) next
    hit <- intersect(aliases[[canonical]], colnames(df))
    if (length(hit) > 0) {
      colnames(df)[match(hit[1], colnames(df))] <- canonical
    }
  }

  if (!"Sample" %in% colnames(df)) {
    sample_aliases <- c(
      "sample",
      "sample_id",
      "SampleID",
      "GSM",
      "geo_accession",
      "Patient"
    )
    hit <- intersect(sample_aliases, colnames(df))
    if (length(hit) == 0) {
      stop("Could not identify Sample column in state-score table.")
    }
    colnames(df)[match(hit[1], colnames(df))] <- "Sample"
  }

  missing_state <- setdiff(state_cols, colnames(df))
  if (length(missing_state) > 0) {
    stop("Missing state-score columns: ", paste(missing_state, collapse = ", "))
  }

  for (st in state_cols) {
    df[[st]] <- as.numeric(df[[st]])
  }

  df$Sample <- as.character(df$Sample)
  df
}

strict_response_group <- function(df) {
  if ("BinaryResponseGroup" %in% colnames(df)) {
    x <- as.character(df$BinaryResponseGroup)
  } else if ("ResponseGroup" %in% colnames(df)) {
    x <- as.character(df$ResponseGroup)
  } else if ("ResponseRaw" %in% colnames(df)) {
    x <- as.character(df$ResponseRaw)
  } else {
    stop("No response field found.")
  }

  y <- toupper(trimws(x))
  y <- gsub("-", "_", y)
  y <- gsub("\\s+", "_", y)

  out <- rep(NA_character_, length(y))

  out[
    y %in% c(
      "RESPONDER",
      "RESPONSE",
      "R",
      "CR",
      "PR",
      "CR/PR",
      "COMPLETE_RESPONSE",
      "PARTIAL_RESPONSE"
    )
  ] <- "Responder"

  out[
    y %in% c(
      "NONRESPONDER",
      "NON_RESPONDER",
      "NON_RESPONSE",
      "NR",
      "PD",
      "PROGRESSIVE_DISEASE",
      "PROGRESSION"
    )
  ] <- "NonResponder"

  out
}

align_expression_to_samples <- function(expr, sample_vec, dataset_name) {
  expr <- collapse_duplicate_genes(expr)
  expr_cols <- colnames(expr)
  sample_vec <- as.character(sample_vec)

  expr_key <- make_sample_key(expr_cols)
  sample_key <- make_sample_key(sample_vec)

  if (anyDuplicated(expr_key)) {
    stop(dataset_name, ": duplicated normalized expression sample keys.")
  }

  idx <- match(sample_key, expr_key)
  matched <- !is.na(idx)

  audit <- data.frame(
    Dataset = dataset_name,
    StateSample = sample_vec,
    StateSampleKey = sample_key,
    MatchedExpressionColumn = ifelse(matched, expr_cols[idx], NA_character_),
    Matched = matched,
    stringsAsFactors = FALSE
  )

  if (!all(matched)) {
    stop(
      dataset_name,
      ": not all state-score samples matched expression columns. ",
      "See matching audit."
    )
  }

  expr2 <- expr[, idx, drop = FALSE]
  colnames(expr2) <- sample_vec

  list(expr = expr2, audit = audit)
}

compute_auc_fixed_direction <- function(y, pred) {
  y <- as.numeric(y)
  pred <- as.numeric(pred)

  keep <- is.finite(y) & is.finite(pred)
  y <- y[keep]
  pred <- pred[keep]

  if (length(y) < 6 || length(unique(y)) != 2) {
    return(NA_real_)
  }

  roc_obj <- pROC::roc(
    response = y,
    predictor = pred,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  as.numeric(pROC::auc(roc_obj))
}

bootstrap_auc_distribution <- function(y, pred, n_boot, seed) {
  set.seed(seed)

  y <- as.numeric(y)
  pred <- as.numeric(pred)

  keep <- is.finite(y) & is.finite(pred)
  y <- y[keep]
  pred <- pred[keep]

  idx0 <- which(y == 0)
  idx1 <- which(y == 1)

  if (length(idx0) < 2 || length(idx1) < 2) {
    return(rep(NA_real_, n_boot))
  }

  out <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    idx <- c(
      sample(idx0, length(idx0), replace = TRUE),
      sample(idx1, length(idx1), replace = TRUE)
    )
    out[b] <- compute_auc_fixed_direction(y[idx], pred[idx])
  }

  out
}

permutation_auc <- function(y, pred, n_perm, seed) {
  set.seed(seed)

  y <- as.numeric(y)
  pred <- as.numeric(pred)

  keep <- is.finite(y) & is.finite(pred)
  y <- y[keep]
  pred <- pred[keep]

  observed_auc <- compute_auc_fixed_direction(y, pred)
  perm <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    perm[i] <- compute_auc_fixed_direction(sample(y), pred)
  }

  p_greater <- (
    sum(perm >= observed_auc, na.rm = TRUE) + 1
  ) / (
    sum(is.finite(perm)) + 1
  )

  list(
    observed_auc = observed_auc,
    p_greater = p_greater,
    perm = perm
  )
}

loso_auc <- function(y, pred, sample_ids) {
  y <- as.numeric(y)
  pred <- as.numeric(pred)
  sample_ids <- as.character(sample_ids)

  keep <- is.finite(y) & is.finite(pred)
  y <- y[keep]
  pred <- pred[keep]
  sample_ids <- sample_ids[keep]

  full_auc <- compute_auc_fixed_direction(y, pred)

  out <- vector("list", length(y))

  for (i in seq_along(y)) {
    out[[i]] <- data.frame(
      RemovedSample = sample_ids[i],
      RemovedClass = y[i],
      Full_AUC = full_auc,
      LeaveOneOut_AUC = compute_auc_fixed_direction(y[-i], pred[-i]),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      Delta_AUC_minus_full = LeaveOneOut_AUC - Full_AUC
    )
}

calc_pROC_bootstrap_ci_seeded <- function(
  response_group,
  predictor,
  boot_n,
  seed
) {
  response_group <- as.character(response_group)
  predictor <- as.numeric(predictor)

  keep <- response_group %in% c("Responder", "NonResponder") &
    is.finite(predictor)

  response_group <- factor(
    response_group[keep],
    levels = c("Responder", "NonResponder")
  )
  predictor <- predictor[keep]

  roc_obj <- pROC::roc(
    response = response_group,
    predictor = predictor,
    levels = c("Responder", "NonResponder"),
    direction = "<",
    quiet = TRUE
  )

  set.seed(seed)

  ci <- as.numeric(
    pROC::ci.auc(
      roc_obj,
      conf.level = 0.95,
      method = "bootstrap",
      boot.n = boot_n,
      boot.stratified = TRUE
    )
  )

  c(
    AUC = as.numeric(pROC::auc(roc_obj)),
    CI_low = ci[1],
    CI_mid = ci[2],
    CI_high = ci[3]
  )
}

save_session_info <- function() {
  f <- file.path(
    log_dir,
    "sessionInfo_05_bulk_response_sensitivity_analysis.txt"
  )
  sink(f)
  print(sessionInfo())
  sink()
  message("Saved: ", f)
}

############################################################
## 5. Input audit and load upstream objects
############################################################

input_audit <- data.frame(
  Input = names(input_files),
  RelativePath = vapply(input_files, relative_to_project, character(1)),
  Exists = file.exists(input_files),
  MD5 = vapply(input_files, safe_md5, character(1)),
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_audit,
  file.path(out_dir, "05_input_file_audit.csv")
)

if (any(!input_audit$Exists)) {
  stop(
    "Missing required upstream input(s): ",
    paste(
      input_audit$RelativePath[!input_audit$Exists],
      collapse = "; "
    )
  )
}

expr244 <- collapse_duplicate_genes(
  readRDS(input_files[["GSE244982_expression"]])
)

state244 <- standardize_state_table(
  data.table::fread(
    input_files[["GSE244982_state_scores"]],
    data.table = FALSE,
    check.names = FALSE
  )
)

expr782 <- collapse_duplicate_genes(
  readRDS(input_files[["GSE78220_expression"]])
)

state782 <- standardize_state_table(
  data.table::fread(
    input_files[["GSE78220_state_scores"]],
    data.table = FALSE,
    check.names = FALSE
  )
)

auc02_raw <- data.table::fread(
  input_files[["GSE78220_primary_auc"]],
  data.table = FALSE,
  check.names = FALSE
)

expr910 <- collapse_duplicate_genes(
  readRDS(input_files[["GSE91061_expression"]])
)

state910_pre <- standardize_state_table(
  data.table::fread(
    input_files[["GSE91061_pretreatment_state_scores"]],
    data.table = FALSE,
    check.names = FALSE
  )
)

state782$StrictResponseGroup <- strict_response_group(state782)
state910_pre$StrictResponseGroup <- strict_response_group(state910_pre)

binary782 <- state782 %>%
  dplyr::filter(
    StrictResponseGroup %in% c("Responder", "NonResponder")
  ) %>%
  dplyr::mutate(
    NonResponse = ifelse(StrictResponseGroup == "NonResponder", 1L, 0L)
  )

binary910 <- state910_pre %>%
  dplyr::filter(
    StrictResponseGroup %in% c("Responder", "NonResponder")
  ) %>%
  dplyr::mutate(
    NonResponse = ifelse(StrictResponseGroup == "NonResponder", 1L, 0L)
  )

############################################################
## 6. Consolidated GSE78220 raw AUC audit
############################################################

auc02_required <- c(
  "State",
  "AUC_raw_score_for_nonresponse",
  "AUC_95CI_low",
  "AUC_95CI_high"
)

if (!all(auc02_required %in% colnames(auc02_raw))) {
  stop(
    "Script-02 AUC table is missing required columns: ",
    paste(setdiff(auc02_required, colnames(auc02_raw)), collapse = ", ")
  )
}

auc02 <- auc02_raw %>%
  dplyr::mutate(
    State = as.character(State)
  ) %>%
  dplyr::filter(State %in% state_cols) %>%
  dplyr::transmute(
    Dataset = "GSE78220",
    State,
    StateLabel = unname(state_labels[State]),
    AUC = as.numeric(AUC_raw_score_for_nonresponse),
    CI_low = as.numeric(AUC_95CI_low),
    CI_high = as.numeric(AUC_95CI_high),
    n_samples = if ("n_samples" %in% colnames(auc02_raw)) {
      as.integer(n_samples)
    } else {
      nrow(binary782)
    },
    n_responders = if ("n_responders" %in% colnames(auc02_raw)) {
      as.integer(n_responders)
    } else {
      sum(binary782$NonResponse == 0)
    },
    n_nonresponders = if ("n_nonresponders" %in% colnames(auc02_raw)) {
      as.integer(n_nonresponders)
    } else {
      sum(binary782$NonResponse == 1)
    },
    BootstrapN = GSE78220_PRIMARY_AUC_BOOT_N,
    BootstrapSeed = GSE78220_PRIMARY_AUC_SEED,
    Source = "02_GSE78220_state_scoring_response_survival.R"
  )

if (nrow(auc02) != length(state_cols)) {
  stop("Expected four states in script-02 AUC table.")
}

direct_auc782 <- data.frame(
  State = state_cols,
  DirectAUC = vapply(
    state_cols,
    function(st) {
      compute_auc_fixed_direction(
        binary782$NonResponse,
        binary782[[st]]
      )
    },
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

auc_consistency <- auc02 %>%
  dplyr::left_join(direct_auc782, by = "State") %>%
  dplyr::mutate(
    AbsoluteDifference = abs(AUC - DirectAUC),
    Pass = AbsoluteDifference <= 1e-12
  )

safe_write_csv(
  auc_consistency,
  file.path(out_dir, "GSE78220_raw_auc_cross_script_consistency_audit.csv")
)

if (any(!auc_consistency$Pass)) {
  stop("GSE78220 raw AUC differs between script 02 and direct recalculation.")
}

## Preserve the deterministic script-02 CI as an audit-only reanalysis.
## Do NOT use it as the canonical S11 figure-source CI because the submitted
## S11 was generated by historical 18C using 5,000 stratified pROC bootstrap
## resamples without a recorded RNG seed.
s11_auc_seeded_reanalysis <- auc02 %>%
  dplyr::mutate(
    CIProvenance =
      "Seeded script-02 reanalysis; audit only; not the submitted S11 CI"
  )

safe_write_csv(
  s11_auc_seeded_reanalysis,
  file.path(
    out_dir,
    "S11_GSE78220_raw_auc_summary_seeded_reanalysis.csv"
  )
)

## Frozen submitted S11 / historical 18C bootstrap CIs.
## These values were recovered from the original 18C result table and are
## intentionally treated as immutable provenance because the historical
## pROC::ci.auc() call used boot.n=5000, stratified=TRUE, with no recorded seed.
s11_legacy_ci <- data.frame(
  State = state_cols,
  Legacy_CI_low = c(
    0.261111111111111,
    0.316666666666667,
    0.600000000000000,
    0.122222222222222
  ),
  Legacy_CI_high = c(
    0.722222222222222,
    0.777777777777778,
    0.944444444444444,
    0.538888888888889
  ),
  stringsAsFactors = FALSE
)

s11_auc <- auc02 %>%
  dplyr::rename(
    SeededReanalysis_CI_low = CI_low,
    SeededReanalysis_CI_high = CI_high,
    SeededReanalysis_BootstrapN = BootstrapN,
    SeededReanalysis_BootstrapSeed = BootstrapSeed,
    SeededReanalysis_Source = Source
  ) %>%
  dplyr::left_join(
    s11_legacy_ci,
    by = "State"
  ) %>%
  dplyr::mutate(
    CI_low = Legacy_CI_low,
    CI_high = Legacy_CI_high,
    BootstrapN = 5000L,
    BootstrapSeed = NA_integer_,
    CIProvenance =
      paste0(
        "Frozen submitted S11 / historical 18C pROC bootstrap CI; ",
        "5000 stratified resamples; RNG seed not recorded"
      ),
    Source =
      "18C_GSE78220_response_association_fragility_audit_revision.R"
  ) %>%
  dplyr::select(
    Dataset,
    State,
    StateLabel,
    AUC,
    CI_low,
    CI_high,
    n_samples,
    n_responders,
    n_nonresponders,
    BootstrapN,
    BootstrapSeed,
    CIProvenance,
    Source,
    SeededReanalysis_CI_low,
    SeededReanalysis_CI_high,
    SeededReanalysis_BootstrapN,
    SeededReanalysis_BootstrapSeed,
    SeededReanalysis_Source
  )

safe_write_csv(
  s11_auc,
  file.path(out_dir, "S11_GSE78220_raw_auc_summary.csv")
)

############################################################
## 7. S11: historical 18C fragility audit
############################################################

## 7.1 Stratified bootstrap distributions.
s11_boot_all <- list()
s11_boot_summary <- list()

for (i in seq_along(state_cols)) {
  st <- state_cols[i]
  seed_i <- S11_BOOT_SEED_BASE + i

  boot_vec <- bootstrap_auc_distribution(
    y = binary782$NonResponse,
    pred = binary782[[st]],
    n_boot = S11_BOOT_N,
    seed = seed_i
  )

  s11_boot_all[[st]] <- data.frame(
    State = st,
    StateLabel = state_labels[[st]],
    BootstrapIteration = seq_len(S11_BOOT_N),
    AUC = boot_vec,
    BootstrapSeed = seed_i,
    stringsAsFactors = FALSE
  )

  s11_boot_summary[[st]] <- data.frame(
    State = st,
    StateLabel = state_labels[[st]],
    n_boot = S11_BOOT_N,
    BootstrapSeed = seed_i,
    median_AUC = median(boot_vec, na.rm = TRUE),
    mean_AUC = mean(boot_vec, na.rm = TRUE),
    CI_low_percentile = as.numeric(
      stats::quantile(boot_vec, 0.025, na.rm = TRUE)
    ),
    CI_high_percentile = as.numeric(
      stats::quantile(boot_vec, 0.975, na.rm = TRUE)
    ),
    min_AUC = min(boot_vec, na.rm = TRUE),
    max_AUC = max(boot_vec, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

s11_boot_all <- dplyr::bind_rows(s11_boot_all)
s11_boot_summary <- dplyr::bind_rows(s11_boot_summary)

safe_write_csv(
  s11_boot_all,
  file.path(out_dir, "S11_GSE78220_bootstrap_auc_values.csv")
)

safe_write_csv(
  s11_boot_summary,
  file.path(out_dir, "S11_GSE78220_bootstrap_auc_summary.csv")
)

## 7.2 Per-state label permutation.
s11_perm_all <- list()
s11_perm_summary <- list()

for (i in seq_along(state_cols)) {
  st <- state_cols[i]
  seed_i <- S11_PERM_SEED_BASE + i

  perm_res <- permutation_auc(
    y = binary782$NonResponse,
    pred = binary782[[st]],
    n_perm = S11_PERM_N,
    seed = seed_i
  )

  s11_perm_all[[st]] <- data.frame(
    State = st,
    StateLabel = state_labels[[st]],
    Permutation = seq_len(S11_PERM_N),
    AUC = perm_res$perm,
    Observed_AUC = perm_res$observed_auc,
    PermutationSeed = seed_i,
    stringsAsFactors = FALSE
  )

  s11_perm_summary[[st]] <- data.frame(
    State = st,
    StateLabel = state_labels[[st]],
    Observed_AUC = perm_res$observed_auc,
    Permutation_mean_AUC = mean(perm_res$perm, na.rm = TRUE),
    Permutation_sd_AUC = stats::sd(perm_res$perm, na.rm = TRUE),
    p_greater = perm_res$p_greater,
    n_perm = S11_PERM_N,
    PermutationSeed = seed_i,
    stringsAsFactors = FALSE
  )
}

s11_perm_all <- dplyr::bind_rows(s11_perm_all)
s11_perm_summary <- dplyr::bind_rows(s11_perm_summary)

safe_write_csv(
  s11_perm_all,
  file.path(out_dir, "S11_GSE78220_permutation_auc_values.csv")
)

safe_write_csv(
  s11_perm_summary,
  file.path(out_dir, "S11_GSE78220_permutation_auc_summary.csv")
)

## 7.3 Leave-one-sample-out.
s11_loso <- dplyr::bind_rows(
  lapply(
    state_cols,
    function(st) {
      loso_auc(
        y = binary782$NonResponse,
        pred = binary782[[st]],
        sample_ids = binary782$Sample
      ) %>%
        dplyr::mutate(
          State = st,
          StateLabel = state_labels[[st]],
          RemovedResponseGroup = ifelse(
            RemovedClass == 1,
            "Non-responder",
            "Responder"
          )
        )
    }
  )
)

s11_loso_summary <- s11_loso %>%
  dplyr::group_by(State, StateLabel) %>%
  dplyr::summarise(
    Full_AUC = dplyr::first(Full_AUC),
    n_leave_one_out = dplyr::n(),
    median_LOSO_AUC = median(LeaveOneOut_AUC, na.rm = TRUE),
    min_LOSO_AUC = min(LeaveOneOut_AUC, na.rm = TRUE),
    max_LOSO_AUC = max(LeaveOneOut_AUC, na.rm = TRUE),
    max_abs_delta_AUC = max(
      abs(Delta_AUC_minus_full),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

safe_write_csv(
  s11_loso,
  file.path(out_dir, "S11_GSE78220_loso_auc_values.csv")
)

safe_write_csv(
  s11_loso_summary,
  file.path(out_dir, "S11_GSE78220_loso_auc_summary.csv")
)

## 7.4 Max-AUC state-selection permutation.
set.seed(S11_MAX_AUC_PERM_SEED)

observed_state_auc <- vapply(
  state_cols,
  function(st) {
    compute_auc_fixed_direction(
      binary782$NonResponse,
      binary782[[st]]
    )
  },
  numeric(1)
)

observed_max_auc <- max(observed_state_auc, na.rm = TRUE)
observed_best_state <- names(which.max(observed_state_auc))

max_perm_vec <- numeric(S11_PERM_N)
max_perm_best_state <- character(S11_PERM_N)

for (i in seq_len(S11_PERM_N)) {
  y_perm <- sample(binary782$NonResponse)

  auc_i <- vapply(
    state_cols,
    function(st) {
      compute_auc_fixed_direction(
        y_perm,
        binary782[[st]]
      )
    },
    numeric(1)
  )

  max_perm_vec[i] <- max(auc_i, na.rm = TRUE)
  max_perm_best_state[i] <- names(which.max(auc_i))
}

s11_max_perm_null <- data.frame(
  Permutation = seq_len(S11_PERM_N),
  Max_AUC_across_states = max_perm_vec,
  BestState_under_permutation = max_perm_best_state,
  Observed_max_AUC = observed_max_auc,
  PermutationSeed = S11_MAX_AUC_PERM_SEED,
  stringsAsFactors = FALSE
)

s11_max_perm_summary <- data.frame(
  Observed_best_state = observed_best_state,
  Observed_best_state_label = state_labels[[observed_best_state]],
  Observed_max_AUC = observed_max_auc,
  Permutation_mean_max_AUC = mean(max_perm_vec, na.rm = TRUE),
  Permutation_sd_max_AUC = stats::sd(max_perm_vec, na.rm = TRUE),
  p_max_AUC_selection_adjusted = (
    sum(max_perm_vec >= observed_max_auc, na.rm = TRUE) + 1
  ) / (
    S11_PERM_N + 1
  ),
  n_perm = S11_PERM_N,
  PermutationSeed = S11_MAX_AUC_PERM_SEED,
  stringsAsFactors = FALSE
)

safe_write_csv(
  s11_max_perm_null,
  file.path(out_dir, "S11_GSE78220_max_auc_permutation_null.csv")
)

safe_write_csv(
  s11_max_perm_summary,
  file.path(out_dir, "S11_GSE78220_max_auc_permutation_summary.csv")
)

############################################################
## 8. S12: historical 18A compact composition-aware audit
############################################################

## IMPORTANT:
## These compact marker sets are the frozen 18A/S12 definitions.
## Do not replace them with the expanded 18G/S13 marker sets.
s12_marker_modules <- list(
  CAF_Stromal = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN",
    "LUM", "FAP", "ACTA2", "PDGFRA"
  ),
  Myeloid = c(
    "LYZ", "CD68", "CSF1R", "LST1", "AIF1"
  ),
  T_NK = c(
    "CD3D", "CD3E", "CD8A", "NKG7", "GZMB"
  ),
  B_Plasma = c(
    "MS4A1", "CD79A", "MZB1", "JCHAIN", "IGKC"
  ),
  Melanoma_Lineage = c(
    "MLANA", "PMEL", "TYR", "MITF", "SOX10"
  ),
  Endothelial = c(
    "PECAM1", "VWF", "KDR", "ENG"
  )
)

s12_marker_table <- dplyr::bind_rows(
  lapply(
    names(s12_marker_modules),
    function(module_name) {
      data.frame(
        Module = module_name,
        Gene = s12_marker_modules[[module_name]],
        stringsAsFactors = FALSE
      )
    }
  )
)

safe_write_csv(
  s12_marker_table,
  file.path(out_dir, "S12_composition_marker_gene_sets.csv")
)

compute_s12_module_scores <- function(expr, modules, dataset_name) {
  expr <- collapse_duplicate_genes(expr)
  gene_z <- zscore_rows(expr)

  score_df <- data.frame(
    Sample = colnames(expr),
    stringsAsFactors = FALSE
  )

  coverage <- list()

  for (module_name in names(modules)) {
    genes <- clean_gene_symbol(modules[[module_name]])
    present <- intersect(genes, rownames(gene_z))

    score_df[[module_name]] <- if (length(present) == 0) {
      rep(NA_real_, ncol(gene_z))
    } else {
      colMeans(
        gene_z[present, , drop = FALSE],
        na.rm = TRUE
      )
    }

    coverage[[module_name]] <- data.frame(
      Dataset = dataset_name,
      Module = module_name,
      Gene = genes,
      Present = genes %in% rownames(gene_z),
      stringsAsFactors = FALSE
    )
  }

  coverage_df <- dplyr::bind_rows(coverage) %>%
    dplyr::group_by(Dataset, Module) %>%
    dplyr::mutate(
      n_requested = dplyr::n(),
      n_present = sum(Present),
      n_missing = sum(!Present)
    ) %>%
    dplyr::ungroup()

  list(scores = score_df, coverage = coverage_df)
}

align244 <- align_expression_to_samples(
  expr244,
  state244$Sample,
  "GSE244982"
)

align782 <- align_expression_to_samples(
  expr782,
  state782$Sample,
  "GSE78220"
)

safe_write_csv(
  align244$audit,
  file.path(out_dir, "S12_GSE244982_expression_state_matching_audit.csv")
)

safe_write_csv(
  align782$audit,
  file.path(out_dir, "S12_GSE78220_expression_state_matching_audit.csv")
)

s12_comp244 <- compute_s12_module_scores(
  align244$expr,
  s12_marker_modules,
  "GSE244982"
)

s12_comp782 <- compute_s12_module_scores(
  align782$expr,
  s12_marker_modules,
  "GSE78220"
)

s12_comp_scores <- dplyr::bind_rows(
  s12_comp244$scores %>% dplyr::mutate(Dataset = "GSE244982"),
  s12_comp782$scores %>% dplyr::mutate(Dataset = "GSE78220")
)

s12_coverage <- dplyr::bind_rows(
  s12_comp244$coverage,
  s12_comp782$coverage
)

safe_write_csv(
  s12_comp_scores,
  file.path(out_dir, "S12_composition_scores_by_sample.csv")
)

safe_write_csv(
  s12_coverage,
  file.path(out_dir, "S12_composition_marker_gene_coverage.csv")
)

s12_merged244 <- state244 %>%
  dplyr::left_join(s12_comp244$scores, by = "Sample")

s12_merged782 <- state782 %>%
  dplyr::left_join(s12_comp782$scores, by = "Sample")

compute_s12_correlations <- function(df, dataset_name) {
  out <- list()

  for (st in state_cols) {
    for (module_name in names(s12_marker_modules)) {
      x <- as.numeric(df[[st]])
      y <- as.numeric(df[[module_name]])

      ok <- is.finite(x) & is.finite(y)

      if (sum(ok) < 4) {
        rho <- NA_real_
        p_value <- NA_real_
      } else {
        ct <- suppressWarnings(
          stats::cor.test(
            x[ok],
            y[ok],
            method = "spearman",
            exact = FALSE
          )
        )
        rho <- unname(ct$estimate)
        p_value <- ct$p.value
      }

      out[[length(out) + 1L]] <- data.frame(
        Dataset = dataset_name,
        State = st,
        StateLabel = state_labels[[st]],
        CompositionModule = module_name,
        n = sum(ok),
        rho = rho,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    }
  }

  dplyr::bind_rows(out) %>%
    dplyr::group_by(Dataset) %>%
    dplyr::mutate(
      p_adj_BH = stats::p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    dplyr::ungroup()
}

s12_correlations <- dplyr::bind_rows(
  compute_s12_correlations(s12_merged244, "GSE244982"),
  compute_s12_correlations(s12_merged782, "GSE78220")
)

safe_write_csv(
  s12_correlations,
  file.path(out_dir, "S12_state_composition_correlations.csv")
)

add_residualized_score <- function(
  df,
  outcome_col,
  covariates,
  new_col
) {
  out <- df
  vars <- c(outcome_col, covariates)

  missing_vars <- setdiff(vars, colnames(out))
  if (length(missing_vars) > 0) {
    stop(
      "Residualization missing variables: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  complete <- stats::complete.cases(
    out[, vars, drop = FALSE]
  )

  resid_vec <- rep(NA_real_, nrow(out))

  if (sum(complete) >= length(covariates) + 3) {
    fml <- stats::as.formula(
      paste(
        outcome_col,
        "~",
        paste(covariates, collapse = " + ")
      )
    )

    fit <- stats::lm(
      fml,
      data = out[complete, , drop = FALSE]
    )

    resid_vec[complete] <- zvec(
      stats::residuals(fit)
    )
  }

  out[[new_col]] <- resid_vec
  out
}

s12_binary782 <- s12_merged782
s12_binary782$StrictResponseGroup <- strict_response_group(s12_binary782)

s12_binary782 <- s12_binary782 %>%
  dplyr::filter(
    StrictResponseGroup %in% c("Responder", "NonResponder")
  )

tumor_state <- "Tumor_dedifferentiation_Stromal_remodeling"

s12_binary782 <- s12_binary782 %>%
  add_residualized_score(
    tumor_state,
    c("CAF_Stromal"),
    "Tumor_state_resid_CAF_Stromal"
  ) %>%
  add_residualized_score(
    tumor_state,
    c("CAF_Stromal", "Melanoma_Lineage"),
    "Tumor_state_resid_CAF_MelanomaLineage"
  ) %>%
  add_residualized_score(
    tumor_state,
    c("CAF_Stromal", "Myeloid", "B_Plasma"),
    "Tumor_state_resid_CAF_Myeloid_BPlasma"
  ) %>%
  add_residualized_score(
    tumor_state,
    c(
      "CAF_Stromal",
      "Myeloid",
      "B_Plasma",
      "Melanoma_Lineage"
    ),
    "Tumor_state_resid_CAF_Myeloid_BPlasma_MelanomaLineage"
  ) %>%
  add_residualized_score(
    tumor_state,
    names(s12_marker_modules),
    "Tumor_state_resid_AllCompositionModules"
  )

safe_write_csv(
  s12_binary782,
  file.path(out_dir, "S12_GSE78220_residualized_scores.csv")
)

s12_models <- data.frame(
  Analysis = c(
    "Raw score",
    "CAF/stromal residualized",
    "CAF + melanoma-lineage residualized",
    "CAF + myeloid + B/plasma residualized",
    "CAF + myeloid + B/plasma + melanoma-lineage residualized",
    "All-composition residualized"
  ),
  Predictor = c(
    tumor_state,
    "Tumor_state_resid_CAF_Stromal",
    "Tumor_state_resid_CAF_MelanomaLineage",
    "Tumor_state_resid_CAF_Myeloid_BPlasma",
    "Tumor_state_resid_CAF_Myeloid_BPlasma_MelanomaLineage",
    "Tumor_state_resid_AllCompositionModules"
  ),
  stringsAsFactors = FALSE
)

s12_auc_rows <- vector("list", nrow(s12_models))

for (i in seq_len(nrow(s12_models))) {
  model <- s12_models[i, , drop = FALSE]

  if (i == 1L) {
    raw_row <- auc02 %>%
      dplyr::filter(State == tumor_state)

    s12_auc_rows[[i]] <- data.frame(
      Analysis = model$Analysis,
      Predictor = model$Predictor,
      n_total = raw_row$n_samples,
      n_responder = raw_row$n_responders,
      n_nonresponder = raw_row$n_nonresponders,
      AUC = raw_row$AUC,
      CI_low = raw_row$CI_low,
      CI_high = raw_row$CI_high,
      BootstrapN = raw_row$BootstrapN,
      BootstrapSeed = raw_row$BootstrapSeed,
      CI_source = "Consolidated from script 02",
      stringsAsFactors = FALSE
    )
  } else {
    seed_i <- S12_RESIDUAL_CI_SEED_BASE + i

    ci <- calc_pROC_bootstrap_ci_seeded(
      response_group = s12_binary782$StrictResponseGroup,
      predictor = s12_binary782[[model$Predictor]],
      boot_n = S12_RESIDUAL_CI_BOOT_N,
      seed = seed_i
    )

    s12_auc_rows[[i]] <- data.frame(
      Analysis = model$Analysis,
      Predictor = model$Predictor,
      n_total = nrow(s12_binary782),
      n_responder = sum(
        s12_binary782$StrictResponseGroup == "Responder"
      ),
      n_nonresponder = sum(
        s12_binary782$StrictResponseGroup == "NonResponder"
      ),
      AUC = unname(ci["AUC"]),
      CI_low = unname(ci["CI_low"]),
      CI_high = unname(ci["CI_high"]),
      BootstrapN = S12_RESIDUAL_CI_BOOT_N,
      BootstrapSeed = seed_i,
      CI_source =
        "18A method preserved; explicit public-release RNG seed added",
      stringsAsFactors = FALSE
    )
  }
}

## Deterministic seeded reanalysis table (audit-only for CI provenance).
s12_auc_seeded_reanalysis <- dplyr::bind_rows(s12_auc_rows)
raw_s12_auc <- s12_auc_seeded_reanalysis$AUC[
  s12_auc_seeded_reanalysis$Analysis == "Raw score"
][1]

s12_auc_seeded_reanalysis$AUC_attenuation_vs_raw <-
  raw_s12_auc - s12_auc_seeded_reanalysis$AUC

safe_write_csv(
  s12_auc_seeded_reanalysis,
  file.path(
    out_dir,
    "S12_GSE78220_raw_vs_composition_residualized_auc_seeded_reanalysis.csv"
  )
)

## Frozen submitted S12 / historical 18A bootstrap CIs.
## Point estimates and residualized scores are recalculated above and must match;
## only the stochastic historical CI layer is frozen for submitted-figure
## reproduction because the original 18A pROC bootstrap used no recorded seed.
s12_legacy_ci <- data.frame(
  Analysis = c(
    "Raw score",
    "CAF/stromal residualized",
    "CAF + melanoma-lineage residualized",
    "CAF + myeloid + B/plasma residualized",
    "CAF + myeloid + B/plasma + melanoma-lineage residualized",
    "All-composition residualized"
  ),
  Legacy_CI_low = c(
    0.600000000000000,
    0.383333333333333,
    0.394444444444444,
    0.388888888888889,
    0.383333333333333,
    0.361111111111111
  ),
  Legacy_CI_high = c(
    0.944444444444444,
    0.833333333333333,
    0.827777777777778,
    0.838888888888889,
    0.822222222222222,
    0.800000000000000
  ),
  stringsAsFactors = FALSE
)

s12_auc <- s12_auc_seeded_reanalysis %>%
  dplyr::rename(
    SeededReanalysis_CI_low = CI_low,
    SeededReanalysis_CI_high = CI_high,
    SeededReanalysis_BootstrapN = BootstrapN,
    SeededReanalysis_BootstrapSeed = BootstrapSeed,
    SeededReanalysis_CI_source = CI_source
  ) %>%
  dplyr::left_join(
    s12_legacy_ci,
    by = "Analysis"
  ) %>%
  dplyr::mutate(
    CI_low = Legacy_CI_low,
    CI_high = Legacy_CI_high,
    BootstrapN = 5000L,
    BootstrapSeed = NA_integer_,
    CI_source =
      paste0(
        "Frozen submitted S12 / historical 18A pROC bootstrap CI; ",
        "5000 stratified resamples; RNG seed not recorded"
      )
  ) %>%
  dplyr::select(
    Analysis,
    Predictor,
    n_total,
    n_responder,
    n_nonresponder,
    AUC,
    CI_low,
    CI_high,
    BootstrapN,
    BootstrapSeed,
    CI_source,
    AUC_attenuation_vs_raw,
    SeededReanalysis_CI_low,
    SeededReanalysis_CI_high,
    SeededReanalysis_BootstrapN,
    SeededReanalysis_BootstrapSeed,
    SeededReanalysis_CI_source
  )

safe_write_csv(
  s12_auc,
  file.path(
    out_dir,
    "S12_GSE78220_raw_vs_composition_residualized_auc.csv"
  )
)

############################################################
## 9. S13: historical 18G expanded formal composition analysis
############################################################

## IMPORTANT:
## These expanded marker sets are the frozen 18G/S13 definitions.
## They intentionally differ from the compact S12/18A marker sets.
s13_marker_modules <- list(
  CAF_Stromal = c(
    "COL1A1", "COL1A2", "COL3A1", "COL5A1",
    "COL6A1", "COL6A2", "COL6A3",
    "DCN", "LUM", "FAP", "ACTA2", "PDPN",
    "THY1", "VCAN", "TNC", "TGFB1",
    "FBN1", "MRC2", "ADAM12", "TIMP1",
    "INHBA", "SPARC", "FN1", "MMP2"
  ),
  Myeloid = c(
    "LYZ", "LST1", "TYROBP", "CD14", "FCGR3A",
    "CD68", "CD163", "MSR1",
    "C1QA", "C1QB", "C1QC",
    "S100A8", "S100A9", "IL1B", "FCER1G", "AIF1"
  ),
  T_NK = c(
    "PTPRC", "CD3D", "CD3E", "CD2", "TRAC",
    "CD8A", "CD8B", "NKG7", "GNLY",
    "GZMB", "PRF1", "FOXP3", "IL2RA",
    "CTLA4", "TIGIT", "LAG3", "PDCD1"
  ),
  B_Plasma = c(
    "MS4A1", "CD79A", "CD79B", "CD74",
    "MZB1", "JCHAIN", "IGHG1", "IGHG3",
    "IGKC", "IGLC1", "SDC1", "XBP1", "TNFRSF17"
  ),
  Melanoma_Lineage = c(
    "MLANA", "PMEL", "TYR", "DCT", "MITF",
    "SOX10", "S100B", "MIA",
    "TFAP2A", "TYRP1", "GPR143", "EDNRB"
  ),
  Endothelial = c(
    "PECAM1", "VWF", "KDR", "FLT1", "CLDN5",
    "RAMP2", "ESAM", "ENG", "EMCN", "PLVAP", "ACKR1"
  )
)

s13_marker_table <- dplyr::bind_rows(
  lapply(
    names(s13_marker_modules),
    function(module_name) {
      data.frame(
        Module = module_name,
        Gene = s13_marker_modules[[module_name]],
        stringsAsFactors = FALSE
      )
    }
  )
)

safe_write_csv(
  s13_marker_table,
  file.path(out_dir, "S13_composition_marker_gene_sets.csv")
)

score_s13_gene_module <- function(expr, genes) {
  genes <- clean_gene_symbol(genes)
  present <- intersect(genes, rownames(expr))

  if (length(present) < 2) {
    return(
      list(
        score = rep(NA_real_, ncol(expr)),
        n_present = length(present),
        n_total = length(unique(genes))
      )
    )
  }

  mat <- expr[present, , drop = FALSE]
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- NA_real_

  score <- colMeans(mat_z, na.rm = TRUE)
  score <- zvec(score)

  list(
    score = score,
    n_present = length(present),
    n_total = length(unique(genes))
  )
}

compute_s13_composition <- function(expr, cohort) {
  score_df <- data.frame(
    sample_id = colnames(expr),
    stringsAsFactors = FALSE
  )
  coverage <- list()

  for (module_name in names(s13_marker_modules)) {
    res <- score_s13_gene_module(
      expr,
      s13_marker_modules[[module_name]]
    )

    score_df[[module_name]] <- res$score

    coverage[[module_name]] <- data.frame(
      cohort = cohort,
      module = module_name,
      n_present = res$n_present,
      n_total = res$n_total,
      coverage = res$n_present / res$n_total,
      stringsAsFactors = FALSE
    )
  }

  score_df <- score_df %>%
    dplyr::mutate(
      Immune_Composite = zvec(
        rowMeans(
          dplyr::across(c(Myeloid, T_NK, B_Plasma)),
          na.rm = TRUE
        )
      ),
      Tumor_Purity_Proxy = zvec(
        Melanoma_Lineage -
          rowMeans(
            dplyr::across(
              c(
                CAF_Stromal,
                Myeloid,
                T_NK,
                B_Plasma,
                Endothelial
              )
            ),
            na.rm = TRUE
          )
      )
    )

  list(
    scores = score_df,
    coverage = dplyr::bind_rows(coverage)
  )
}

create_marker_pc_scores <- function(
  comp_df,
  module_cols,
  prefix = "MarkerComp"
) {
  out <- data.frame(
    sample_id = comp_df$sample_id,
    stringsAsFactors = FALSE
  )

  available <- intersect(
    module_cols,
    colnames(comp_df)
  )

  dat <- comp_df[, available, drop = FALSE]

  valid <- vapply(
    dat,
    function(x) {
      sum(is.finite(x)) >= 5 &&
        stats::sd(x, na.rm = TRUE) > 0
    },
    logical(1)
  )

  dat <- dat[, valid, drop = FALSE]

  if (ncol(dat) < 2) {
    return(out)
  }

  dat_imp <- as.data.frame(dat)

  for (cn in colnames(dat_imp)) {
    m <- mean(dat_imp[[cn]], na.rm = TRUE)
    dat_imp[[cn]][!is.finite(dat_imp[[cn]])] <- m
  }

  pc <- stats::prcomp(
    dat_imp,
    center = TRUE,
    scale. = TRUE
  )

  n_pc <- min(2L, ncol(pc$x))

  for (i in seq_len(n_pc)) {
    out[[paste0(prefix, "_PC", i)]] <- zvec(pc$x[, i])
  }

  out
}

safe_residualize <- function(score, covar_df) {
  score <- as.numeric(score)

  if (ncol(covar_df) == 0) {
    return(zvec(score))
  }

  dat <- data.frame(
    score = score,
    covar_df,
    check.names = FALSE
  )

  complete <- stats::complete.cases(dat)
  dat_complete <- dat[complete, , drop = FALSE]

  out <- rep(NA_real_, length(score))

  if (nrow(dat_complete) < 8) {
    return(out)
  }

  form <- stats::as.formula(
    paste(
      "score ~",
      paste(colnames(covar_df), collapse = " + ")
    )
  )

  mm <- stats::model.matrix(
    form,
    data = dat_complete
  )

  if (qr(mm)$rank < ncol(mm)) {
    return(out)
  }

  if ((nrow(dat_complete) - qr(mm)$rank) < 3) {
    return(out)
  }

  fit <- tryCatch(
    stats::lm(form, data = dat_complete),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(out)
  }

  out[which(complete)] <- stats::residuals(fit)
  zvec(out)
}

run_mcp_counter <- function(expr) {
  if (!HAS_MCP) {
    return(
      list(
        scores = NULL,
        log = "MCPcounter not installed; skipped."
      )
    )
  }

  res <- tryCatch(
    {
      mcp <- MCPcounter::MCPcounter.estimate(
        expr,
        featuresType = "HUGO_symbols"
      )

      mcp <- as.data.frame(
        t(mcp),
        check.names = FALSE
      )

      mcp$sample_id <- rownames(mcp)

      rename_map <- c(
        "Fibroblasts" = "MCP_Fibroblasts",
        "Monocytic lineage" = "MCP_Monocytic_lineage",
        "T cells" = "MCP_T_cells",
        "CD8 T cells" = "MCP_CD8_T_cells",
        "B lineage" = "MCP_B_lineage",
        "NK cells" = "MCP_NK_cells",
        "Endothelial cells" = "MCP_Endothelial_cells"
      )

      for (old in names(rename_map)) {
        if (old %in% colnames(mcp)) {
          colnames(mcp)[
            colnames(mcp) == old
          ] <- rename_map[[old]]
        }
      }

      numeric_cols <- setdiff(
        colnames(mcp),
        "sample_id"
      )

      mcp[numeric_cols] <- lapply(
        mcp[numeric_cols],
        zvec
      )

      list(
        scores = mcp,
        log = "MCPcounter completed."
      )
    },
    error = function(e) {
      list(
        scores = NULL,
        log = paste(
          "MCPcounter failed:",
          conditionMessage(e)
        )
      )
    }
  )

  res
}

write_gct <- function(expr, path) {
  df <- data.frame(
    Name = rownames(expr),
    Description = rownames(expr),
    as.data.frame(expr, check.names = FALSE),
    check.names = FALSE
  )

  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("#1.2", con)
  writeLines(
    paste(nrow(df), ncol(df) - 2L, sep = "\t"),
    con
  )

  utils::write.table(
    df,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

run_estimate_scores <- function(
  expr,
  cohort_name
) {
  if (!HAS_ESTIMATE) {
    return(
      list(
        scores = NULL,
        log = "estimate not installed; skipped."
      )
    )
  }

  res <- tryCatch(
    {
      est_dir <- file.path(
        out_intermediate_dir,
        paste0("estimate_", cohort_name)
      )

      dir.create(
        est_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      input_gct <- file.path(
        est_dir,
        paste0(cohort_name, "_input.gct")
      )

      common_gct <- file.path(
        est_dir,
        paste0(cohort_name, "_common.gct")
      )

      score_gct <- file.path(
        est_dir,
        paste0(cohort_name, "_estimate_scores.gct")
      )

      write_gct(expr, input_gct)

      estimate::filterCommonGenes(
        input.f = input_gct,
        output.f = common_gct,
        id = "GeneSymbol"
      )

      estimate::estimateScore(
        input.ds = common_gct,
        output.ds = score_gct,
        platform = "illumina"
      )

      score_raw <- data.table::fread(
        score_gct,
        skip = 2,
        data.table = FALSE,
        check.names = FALSE
      )

      rownames(score_raw) <- score_raw$NAME
      score_raw$NAME <- NULL

      if ("Description" %in% colnames(score_raw)) {
        score_raw$Description <- NULL
      }

      score_df <- as.data.frame(
        t(score_raw),
        check.names = FALSE
      )

      score_df$sample_id <- rownames(score_df)

      for (cn in setdiff(colnames(score_df), "sample_id")) {
        score_df[[cn]] <- zvec(
          as.numeric(score_df[[cn]])
        )
      }

      colnames(score_df) <- gsub(
        "StromalScore",
        "ESTIMATE_StromalScore",
        colnames(score_df)
      )

      colnames(score_df) <- gsub(
        "ImmuneScore",
        "ESTIMATE_ImmuneScore",
        colnames(score_df)
      )

      colnames(score_df) <- gsub(
        "ESTIMATEScore",
        "ESTIMATE_ESTIMATEScore",
        colnames(score_df)
      )

      colnames(score_df) <- gsub(
        "TumorPurity",
        "ESTIMATE_TumorPurity",
        colnames(score_df)
      )

      list(
        scores = score_df,
        log = "ESTIMATE completed."
      )
    },
    error = function(e) {
      list(
        scores = NULL,
        log = paste(
          "ESTIMATE failed:",
          conditionMessage(e)
        )
      )
    }
  )

  res
}

prepare_formal_cohort <- function(
  cohort,
  expr,
  state_df
) {
  state_df <- standardize_state_table(state_df)

  align <- align_expression_to_samples(
    expr,
    state_df$Sample,
    cohort
  )

  expr_aligned <- align$expr

  state_work <- state_df %>%
    dplyr::arrange(
      match(Sample, colnames(expr_aligned))
    )

  state_score_df <- data.frame(
    sample_id = state_work$Sample,
    stringsAsFactors = FALSE
  )

  for (st in state_cols) {
    state_score_df[[st]] <- zvec(
      state_work[[st]]
    )
  }

  response_group <- strict_response_group(state_work)

  clinical_df <- data.frame(
    sample_id = state_work$Sample,
    response_group = response_group,
    non_response = ifelse(
      response_group == "NonResponder",
      1L,
      ifelse(
        response_group == "Responder",
        0L,
        NA_integer_
      )
    ),
    stringsAsFactors = FALSE
  )

  comp_res <- compute_s13_composition(
    expr_aligned,
    cohort
  )

  comp_df <- comp_res$scores

  mcp_res <- run_mcp_counter(expr_aligned)
  if (!is.null(mcp_res$scores)) {
    comp_df <- comp_df %>%
      dplyr::left_join(
        mcp_res$scores,
        by = "sample_id"
      )
  }

  estimate_res <- run_estimate_scores(
    expr_aligned,
    cohort
  )
  if (!is.null(estimate_res$scores)) {
    comp_df <- comp_df %>%
      dplyr::left_join(
        estimate_res$scores,
        by = "sample_id"
      )
  }

  marker_cols <- c(
    "CAF_Stromal",
    "Myeloid",
    "T_NK",
    "B_Plasma",
    "Melanoma_Lineage",
    "Endothelial"
  )

  pc_df <- create_marker_pc_scores(
    comp_df,
    marker_cols,
    prefix = "MarkerComp"
  )

  comp_df <- comp_df %>%
    dplyr::left_join(
      pc_df,
      by = "sample_id"
    )

  dat <- clinical_df %>%
    dplyr::left_join(
      state_score_df,
      by = "sample_id"
    ) %>%
    dplyr::left_join(
      comp_df,
      by = "sample_id"
    )

  method_log <- data.frame(
    cohort = cohort,
    method = c(
      "MCPcounter",
      "ESTIMATE",
      "CompositionPC"
    ),
    status = c(
      mcp_res$log,
      estimate_res$log,
      if (
        all(
          c("MarkerComp_PC1", "MarkerComp_PC2") %in%
            colnames(comp_df)
        )
      ) {
        "Composition PC1/PC2 available."
      } else {
        "Composition PC1/PC2 missing."
      }
    ),
    stringsAsFactors = FALSE
  )

  list(
    data = dat,
    composition_scores = comp_df,
    coverage = comp_res$coverage,
    method_log = method_log,
    matching_audit = align$audit
  )
}

formal782 <- prepare_formal_cohort(
  "GSE78220",
  expr782,
  state782
)

formal910 <- prepare_formal_cohort(
  "GSE91061",
  expr910,
  state910_pre
)

s13_matching_audit <- dplyr::bind_rows(
  formal782$matching_audit,
  formal910$matching_audit
)

safe_write_csv(
  s13_matching_audit,
  file.path(out_dir, "S13_expression_state_matching_audit.csv")
)

s13_comp_scores <- dplyr::bind_rows(
  formal782$composition_scores %>%
    dplyr::mutate(cohort = "GSE78220"),
  formal910$composition_scores %>%
    dplyr::mutate(cohort = "GSE91061")
)

s13_coverage <- dplyr::bind_rows(
  formal782$coverage,
  formal910$coverage
)

s13_method_log <- dplyr::bind_rows(
  formal782$method_log,
  formal910$method_log
)

safe_write_csv(
  s13_comp_scores,
  file.path(out_dir, "S13_composition_scores_by_sample.csv")
)

safe_write_csv(
  s13_coverage,
  file.path(out_dir, "S13_composition_marker_gene_coverage.csv")
)

safe_write_csv(
  s13_method_log,
  file.path(out_dir, "S13_optional_method_availability_log.csv")
)

compute_formal_state_comp_cor <- function(
  dat,
  cohort
) {
  comp_cols <- setdiff(
    colnames(dat),
    c(
      "sample_id",
      "response_group",
      "non_response",
      state_cols
    )
  )

  out <- list()

  for (st in state_cols) {
    for (cc in comp_cols) {
      x <- as.numeric(dat[[st]])
      y <- suppressWarnings(
        as.numeric(dat[[cc]])
      )

      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 6) next

      ct <- suppressWarnings(
        stats::cor.test(
          x[ok],
          y[ok],
          method = "spearman",
          exact = FALSE
        )
      )

      out[[length(out) + 1L]] <- data.frame(
        cohort = cohort,
        state = state_labels[[st]],
        state_internal = st,
        composition_feature = cc,
        n = sum(ok),
        spearman_rho = unname(ct$estimate),
        p_value = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  dplyr::bind_rows(out)
}

s13_state_comp_cor <- dplyr::bind_rows(
  compute_formal_state_comp_cor(
    formal782$data,
    "GSE78220"
  ),
  compute_formal_state_comp_cor(
    formal910$data,
    "GSE91061"
  )
) %>%
  dplyr::group_by(cohort) %>%
  dplyr::mutate(
    q_value = stats::p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

safe_write_csv(
  s13_state_comp_cor,
  file.path(out_dir, "S13_state_composition_correlations.csv")
)

run_formal_auc_analysis <- function(
  dat,
  cohort
) {
  dat_bin <- dat %>%
    dplyr::filter(!is.na(non_response))

  adjustment_sets <- list(
    Raw_score = character(0),
    CAF_stromal_residualized = c("CAF_Stromal"),
    CAF_plus_melanoma_lineage_residualized = c(
      "CAF_Stromal",
      "Melanoma_Lineage"
    ),
    All_marker_composition_residualized = c(
      "CAF_Stromal",
      "Myeloid",
      "T_NK",
      "B_Plasma",
      "Melanoma_Lineage",
      "Endothelial"
    ),
    Purity_proxy_residualized = c(
      "Tumor_Purity_Proxy"
    ),
    Marker_composition_PC1_PC2_residualized = c(
      "MarkerComp_PC1",
      "MarkerComp_PC2"
    )
  )

  if (
    all(
      c(
        "ESTIMATE_StromalScore",
        "ESTIMATE_ImmuneScore"
      ) %in% colnames(dat_bin)
    )
  ) {
    adjustment_sets[["ESTIMATE_residualized"]] <- intersect(
      c(
        "ESTIMATE_StromalScore",
        "ESTIMATE_ImmuneScore",
        "ESTIMATE_TumorPurity",
        "ESTIMATE_ESTIMATEScore"
      ),
      colnames(dat_bin)
    )
  }

  mcp_covars <- intersect(
    c(
      "MCP_Fibroblasts",
      "MCP_Monocytic_lineage",
      "MCP_T_cells",
      "MCP_CD8_T_cells",
      "MCP_B_lineage",
      "MCP_Endothelial_cells"
    ),
    colnames(dat_bin)
  )

  if (length(mcp_covars) >= 2) {
    adjustment_sets[["MCPcounter_residualized"]] <- mcp_covars
  }

  summary_rows <- list()
  loso_rows <- list()
  boot_rows <- list()
  perm_rows <- list()

  for (state_index in seq_along(state_cols)) {
    st <- state_cols[state_index]

    for (
      method_index in seq_along(adjustment_sets)
    ) {
      adjustment <- names(adjustment_sets)[method_index]
      covars <- adjustment_sets[[method_index]]
      covars <- intersect(
        covars,
        colnames(dat_bin)
      )

      if (
        adjustment != "Raw_score" &&
        length(covars) == 0
      ) {
        next
      }

      pred_raw <- dat_bin[[st]]

      pred <- if (adjustment == "Raw_score") {
        zvec(pred_raw)
      } else {
        safe_residualize(
          pred_raw,
          dat_bin[, covars, drop = FALSE]
        )
      }

      y <- dat_bin$non_response
      sid <- dat_bin$sample_id

      auc_obs <- compute_auc_fixed_direction(
        y,
        pred
      )

      boot_seed <-
        S13_SEED_BASE +
        state_index +
        method_index * 100L

      boot_vec <- bootstrap_auc_distribution(
        y = y,
        pred = pred,
        n_boot = S13_BOOT_N,
        seed = boot_seed
      )

      ci_low <- as.numeric(
        stats::quantile(
          boot_vec,
          0.025,
          na.rm = TRUE
        )
      )

      ci_high <- as.numeric(
        stats::quantile(
          boot_vec,
          0.975,
          na.rm = TRUE
        )
      )

      loso_df <- loso_auc(
        y = y,
        pred = pred,
        sample_ids = sid
      ) %>%
        dplyr::transmute(
          cohort = cohort,
          state = state_labels[[st]],
          state_internal = st,
          adjustment = adjustment,
          removed_sample = RemovedSample,
          removed_class = RemovedClass,
          auc = LeaveOneOut_AUC,
          Full_AUC = Full_AUC
        )

      perm_seed <-
        S13_SEED_BASE +
        state_index * 10L +
        method_index * 1000L

      perm_res <- permutation_auc(
        y = y,
        pred = pred,
        n_perm = S13_PERM_N,
        seed = perm_seed
      )

      summary_rows[[length(summary_rows) + 1L]] <-
        data.frame(
          cohort = cohort,
          outcome_mode = "strict",
          n_total = sum(is.finite(y) & is.finite(pred)),
          n_responders = sum(
            y == 0 & is.finite(pred)
          ),
          n_nonresponders = sum(
            y == 1 & is.finite(pred)
          ),
          state = state_labels[[st]],
          state_internal = st,
          adjustment = adjustment,
          covariates = if (
            length(covars) == 0
          ) {
            "none"
          } else {
            paste(covars, collapse = "; ")
          },
          auc = auc_obs,
          bootstrap_ci_low = ci_low,
          bootstrap_ci_high = ci_high,
          bootstrap_n = S13_BOOT_N,
          bootstrap_seed = boot_seed,
          loso_auc_median = median(
            loso_df$auc,
            na.rm = TRUE
          ),
          loso_auc_min = min(
            loso_df$auc,
            na.rm = TRUE
          ),
          loso_auc_max = max(
            loso_df$auc,
            na.rm = TRUE
          ),
          permutation_p = perm_res$p_greater,
          permutation_n = S13_PERM_N,
          permutation_seed = perm_seed,
          stringsAsFactors = FALSE
        )

      loso_rows[[length(loso_rows) + 1L]] <-
        loso_df

      boot_rows[[length(boot_rows) + 1L]] <-
        data.frame(
          cohort = cohort,
          state = state_labels[[st]],
          state_internal = st,
          adjustment = adjustment,
          bootstrap_id = seq_len(S13_BOOT_N),
          auc = boot_vec,
          bootstrap_seed = boot_seed,
          stringsAsFactors = FALSE
        )

      perm_rows[[length(perm_rows) + 1L]] <-
        data.frame(
          cohort = cohort,
          state = state_labels[[st]],
          state_internal = st,
          adjustment = adjustment,
          permutation_id = seq_len(S13_PERM_N),
          auc = perm_res$perm,
          permutation_seed = perm_seed,
          stringsAsFactors = FALSE
        )
    }
  }

  list(
    summary = dplyr::bind_rows(summary_rows),
    loso = dplyr::bind_rows(loso_rows),
    bootstrap = dplyr::bind_rows(boot_rows),
    permutation = dplyr::bind_rows(perm_rows)
  )
}

s13_auc782 <- run_formal_auc_analysis(
  formal782$data,
  "GSE78220"
)

s13_auc910 <- run_formal_auc_analysis(
  formal910$data,
  "GSE91061"
)

s13_auc_summary <- dplyr::bind_rows(
  s13_auc782$summary,
  s13_auc910$summary
) %>%
  dplyr::mutate(
    auc_ci = sprintf(
      "%.3f (%.3f–%.3f)",
      auc,
      bootstrap_ci_low,
      bootstrap_ci_high
    ),
    loso_range = sprintf(
      "%.3f–%.3f",
      loso_auc_min,
      loso_auc_max
    )
  )

s13_loso <- dplyr::bind_rows(
  s13_auc782$loso,
  s13_auc910$loso
)

s13_boot <- dplyr::bind_rows(
  s13_auc782$bootstrap,
  s13_auc910$bootstrap
)

s13_perm <- dplyr::bind_rows(
  s13_auc782$permutation,
  s13_auc910$permutation
)

safe_write_csv(
  s13_auc_summary,
  file.path(out_dir, "S13_formal_composition_adjusted_auc_summary.csv")
)

safe_write_csv(
  s13_loso,
  file.path(out_dir, "S13_loso_auc_values.csv")
)

safe_write_csv(
  s13_boot,
  file.path(out_dir, "S13_bootstrap_auc_values.csv")
)

safe_write_csv(
  s13_perm,
  file.path(out_dir, "S13_permutation_auc_values.csv")
)

############################################################
## 9B. Frozen legacy-CI provenance record for submitted S11/S12
############################################################

legacy_ci_provenance <- dplyr::bind_rows(
  s11_auc %>%
    dplyr::transmute(
      Figure = "Supplementary Figure S11",
      HistoricalModule = "18C",
      HistoricalResult =
        "18C_GSE78220_raw_AUC_all_states_bootstrap_CI.csv",
      Target = StateLabel,
      PointEstimate = AUC,
      Frozen_CI_low = CI_low,
      Frozen_CI_high = CI_high,
      HistoricalBootstrapN = BootstrapN,
      HistoricalBootstrapSeed = BootstrapSeed,
      SeedStatus = "Not recorded in historical CI call",
      SeededReanalysis_CI_low = SeededReanalysis_CI_low,
      SeededReanalysis_CI_high = SeededReanalysis_CI_high,
      PublicUse =
        "Frozen legacy CI is canonical for exact submitted S11 reproduction"
    ),
  s12_auc %>%
    dplyr::transmute(
      Figure = "Supplementary Figure S12",
      HistoricalModule = "18A",
      HistoricalResult =
        "GSE78220_raw_vs_composition_residualized_AUC.csv",
      Target = Analysis,
      PointEstimate = AUC,
      Frozen_CI_low = CI_low,
      Frozen_CI_high = CI_high,
      HistoricalBootstrapN = BootstrapN,
      HistoricalBootstrapSeed = BootstrapSeed,
      SeedStatus = "Not recorded in historical CI call",
      SeededReanalysis_CI_low = SeededReanalysis_CI_low,
      SeededReanalysis_CI_high = SeededReanalysis_CI_high,
      PublicUse =
        "Frozen legacy CI is canonical for exact submitted S12 reproduction"
    )
)

safe_write_csv(
  legacy_ci_provenance,
  file.path(out_dir, "05_legacy_ci_provenance.csv")
)

############################################################
## 10. Reproducibility gates
############################################################

expected_auc782 <- c(
  "Immune_defective_Cold" = 0.488888888888889,
  "Myeloid_Treg_Immunosuppressive" = 0.555555555555556,
  "Tumor_dedifferentiation_Stromal_remodeling" = 0.794444444444444,
  "Melanocytic_Differentiation" = 0.316666666666667
)

expected_auc910 <- c(
  "Immune_defective_Cold" = 0.660869565217391,
  "Myeloid_Treg_Immunosuppressive" = 0.365217391304348,
  "Tumor_dedifferentiation_Stromal_remodeling" = 0.434782608695652,
  "Melanocytic_Differentiation" = 0.621739130434783
)

gate_rows <- list()

add_gate <- function(
  gate,
  observed,
  expected,
  tolerance = 0
) {
  data.frame(
    Gate = gate,
    Observed = as.character(observed),
    Expected = as.character(expected),
    Tolerance = as.character(tolerance),
    Pass = if (
      is.numeric(observed) &&
      is.numeric(expected)
    ) {
      is.finite(observed) &&
        abs(observed - expected) <= tolerance
    } else {
      identical(
        as.character(observed),
        as.character(expected)
      )
    },
    stringsAsFactors = FALSE
  )
}

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE244982_state_samples",
    nrow(state244),
    41
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE78220_state_samples",
    nrow(state782),
    28
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE78220_binary_samples",
    nrow(binary782),
    27
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE78220_responders",
    sum(binary782$NonResponse == 0),
    15
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE78220_nonresponders",
    sum(binary782$NonResponse == 1),
    12
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE91061_pretreatment_samples",
    nrow(state910_pre),
    51
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE91061_binary_samples",
    nrow(binary910),
    33
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE91061_responders",
    sum(binary910$NonResponse == 0),
    10
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "GSE91061_nonresponders",
    sum(binary910$NonResponse == 1),
    23
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "S11_LOSO_rows",
    nrow(s11_loso),
    27 * 4
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "S11_max_AUC_permutations",
    nrow(s11_max_perm_null),
    S11_PERM_N
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "S12_state_composition_correlation_rows",
    nrow(s12_correlations),
    2 * 4 * 6
  )

gate_rows[[length(gate_rows) + 1L]] <-
  add_gate(
    "S12_residualization_models",
    nrow(s12_auc),
    6
  )

## Legacy submitted-CI provenance gates.
expected_s11_legacy_low <- setNames(
  s11_legacy_ci$Legacy_CI_low,
  s11_legacy_ci$State
)
expected_s11_legacy_high <- setNames(
  s11_legacy_ci$Legacy_CI_high,
  s11_legacy_ci$State
)

for (st in state_cols) {
  row_st <- s11_auc %>%
    dplyr::filter(State == st)

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0("S11_legacy_CI_low_", st),
      row_st$CI_low[1],
      expected_s11_legacy_low[[st]],
      1e-12
    )

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0("S11_legacy_CI_high_", st),
      row_st$CI_high[1],
      expected_s11_legacy_high[[st]],
      1e-12
    )

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0("S11_point_AUC_unchanged_", st),
      row_st$AUC[1],
      auc02$AUC[match(st, auc02$State)],
      1e-12
    )
}

for (analysis_name in s12_legacy_ci$Analysis) {
  row_model <- s12_auc %>%
    dplyr::filter(Analysis == analysis_name)

  legacy_row <- s12_legacy_ci %>%
    dplyr::filter(Analysis == analysis_name)

  reanalysis_row <- s12_auc_seeded_reanalysis %>%
    dplyr::filter(Analysis == analysis_name)

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0(
        "S12_legacy_CI_low_",
        gsub("[^A-Za-z0-9]+", "_", analysis_name)
      ),
      row_model$CI_low[1],
      legacy_row$Legacy_CI_low[1],
      1e-12
    )

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0(
        "S12_legacy_CI_high_",
        gsub("[^A-Za-z0-9]+", "_", analysis_name)
      ),
      row_model$CI_high[1],
      legacy_row$Legacy_CI_high[1],
      1e-12
    )

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0(
        "S12_point_AUC_unchanged_",
        gsub("[^A-Za-z0-9]+", "_", analysis_name)
      ),
      row_model$AUC[1],
      reanalysis_row$AUC[1],
      1e-12
    )
}

for (st in state_cols) {
  obs <- direct_auc782$DirectAUC[
    match(st, direct_auc782$State)
  ]

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0("GSE78220_raw_AUC_", st),
      obs,
      expected_auc782[[st]],
      1e-12
    )
}

direct_auc910 <- data.frame(
  State = state_cols,
  DirectAUC = vapply(
    state_cols,
    function(st) {
      compute_auc_fixed_direction(
        binary910$NonResponse,
        binary910[[st]]
      )
    },
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

for (st in state_cols) {
  obs <- direct_auc910$DirectAUC[
    match(st, direct_auc910$State)
  ]

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0("GSE91061_raw_AUC_", st),
      obs,
      expected_auc910[[st]],
      1e-12
    )
}

target_label <- state_labels[[tumor_state]]

for (cohort in c("GSE78220", "GSE91061")) {
  core_methods <- c(
    "Raw_score",
    "CAF_stromal_residualized",
    "CAF_plus_melanoma_lineage_residualized",
    "All_marker_composition_residualized",
    "Purity_proxy_residualized",
    "Marker_composition_PC1_PC2_residualized"
  )

  target_methods <- s13_auc_summary %>%
    dplyr::filter(
      cohort == .env$cohort,
      state == target_label
    ) %>%
    dplyr::pull(adjustment)

  missing_core <- setdiff(
    core_methods,
    target_methods
  )

  gate_rows[[length(gate_rows) + 1L]] <-
    add_gate(
      paste0(
        cohort,
        "_S13_core_target_state_methods_present"
      ),
      length(missing_core),
      0
    )
}

reproducibility_gates <- dplyr::bind_rows(
  gate_rows
)

safe_write_csv(
  reproducibility_gates,
  file.path(out_dir, "05_reproducibility_gates.csv")
)

if (any(!reproducibility_gates$Pass)) {
  failed <- reproducibility_gates %>%
    dplyr::filter(!Pass)

  stop(
    "One or more public-release reproducibility gates failed:\n",
    paste(
      paste0(
        failed$Gate,
        " [observed=",
        failed$Observed,
        "; expected=",
        failed$Expected,
        "]"
      ),
      collapse = "\n"
    )
  )
}

############################################################
## 11. Source/output dependency audit
############################################################

dependency_audit <- data.frame(
  HistoricalModule = c(
    "18C / S11",
    "18C / S11",
    "18C / S11",
    "18A / S12",
    "18A / S12",
    "18G / S13",
    "18G / S13",
    "18G / S13"
  ),
  HistoricalOutput = c(
    "18C_GSE78220_raw_AUC_all_states_bootstrap_CI.csv",
    "18C_GSE78220_leave_one_sample_out_AUC_influence_all_states.csv",
    "18C_GSE78220_max_AUC_state_selection_permutation_null.csv",
    "ALL_state_composition_spearman_correlation.csv",
    "GSE78220_raw_vs_composition_residualized_AUC.csv",
    "18G_formal_composition_adjusted_AUC_summary.csv",
    "18G_LOSO_AUC_values.csv",
    "18G_state_composition_correlations.csv"
  ),
  PublicOutput = c(
    "S11_GSE78220_raw_auc_summary.csv",
    "S11_GSE78220_loso_auc_values.csv",
    "S11_GSE78220_max_auc_permutation_null.csv",
    "S12_state_composition_correlations.csv",
    "S12_GSE78220_raw_vs_composition_residualized_auc.csv",
    "S13_formal_composition_adjusted_auc_summary.csv",
    "S13_loso_auc_values.csv",
    "S13_state_composition_correlations.csv"
  ),
  PublicDecision = c(
    paste0(
      "AUC point estimate preserved; submitted 18C 5000-bootstrap CI frozen ",
      "verbatim because historical CI RNG seed was not recorded; seeded ",
      "script-02 CI retained separately as audit"
    ),
    "Preserved; same 27-sample binary subset",
    "Preserved; 10,000 permutations; seed 20260701",
    "Preserved compact 18A marker panels; 48 correlations",
    paste0(
      "Point/residual definitions preserved; submitted 18A 5000-bootstrap ",
      "CIs frozen verbatim because historical CI RNG seed was not recorded; ",
      "seeded reanalysis retained separately as audit"
    ),
    "Preserved expanded 18G marker panels and seed schedule",
    "Preserved",
    "Preserved with canonical state labels"
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  dependency_audit,
  file.path(out_dir, "05_source_output_dependency_audit.csv")
)

############################################################
## 12. Save compact analysis object, output inventory, session info
############################################################

safe_save_rds(
  list(
    S11_raw_auc = s11_auc,
    S11_raw_auc_seeded_reanalysis = s11_auc_seeded_reanalysis,
    S11_loso_summary = s11_loso_summary,
    S11_max_auc_permutation_summary = s11_max_perm_summary,
    S12_correlations = s12_correlations,
    S12_residualized_auc = s12_auc,
    S12_residualized_auc_seeded_reanalysis = s12_auc_seeded_reanalysis,
    S13_auc_summary = s13_auc_summary,
    S13_method_log = s13_method_log,
    legacy_ci_provenance = legacy_ci_provenance,
    reproducibility_gates = reproducibility_gates
  ),
  file.path(
    out_intermediate_dir,
    "bulk_response_sensitivity_summary.rds"
  )
)

save_session_info()

canonical_outputs <- c(
  file.path(out_dir, "05_input_file_audit.csv"),
  file.path(out_dir, "05_source_output_dependency_audit.csv"),
  file.path(out_dir, "05_reproducibility_gates.csv"),
  file.path(out_dir, "05_legacy_ci_provenance.csv"),
  file.path(out_dir, "GSE78220_raw_auc_cross_script_consistency_audit.csv"),
  file.path(out_dir, "S11_GSE78220_raw_auc_summary.csv"),
  file.path(out_dir, "S11_GSE78220_raw_auc_summary_seeded_reanalysis.csv"),
  file.path(out_dir, "S11_GSE78220_bootstrap_auc_values.csv"),
  file.path(out_dir, "S11_GSE78220_bootstrap_auc_summary.csv"),
  file.path(out_dir, "S11_GSE78220_permutation_auc_values.csv"),
  file.path(out_dir, "S11_GSE78220_permutation_auc_summary.csv"),
  file.path(out_dir, "S11_GSE78220_loso_auc_values.csv"),
  file.path(out_dir, "S11_GSE78220_loso_auc_summary.csv"),
  file.path(out_dir, "S11_GSE78220_max_auc_permutation_null.csv"),
  file.path(out_dir, "S11_GSE78220_max_auc_permutation_summary.csv"),
  file.path(out_dir, "S12_composition_marker_gene_sets.csv"),
  file.path(out_dir, "S12_composition_scores_by_sample.csv"),
  file.path(out_dir, "S12_composition_marker_gene_coverage.csv"),
  file.path(out_dir, "S12_state_composition_correlations.csv"),
  file.path(out_dir, "S12_GSE78220_residualized_scores.csv"),
  file.path(out_dir, "S12_GSE78220_raw_vs_composition_residualized_auc.csv"),
  file.path(out_dir, "S12_GSE78220_raw_vs_composition_residualized_auc_seeded_reanalysis.csv"),
  file.path(out_dir, "S13_composition_marker_gene_sets.csv"),
  file.path(out_dir, "S13_composition_scores_by_sample.csv"),
  file.path(out_dir, "S13_composition_marker_gene_coverage.csv"),
  file.path(out_dir, "S13_optional_method_availability_log.csv"),
  file.path(out_dir, "S13_state_composition_correlations.csv"),
  file.path(out_dir, "S13_formal_composition_adjusted_auc_summary.csv"),
  file.path(out_dir, "S13_loso_auc_values.csv"),
  file.path(out_dir, "S13_bootstrap_auc_values.csv"),
  file.path(out_dir, "S13_permutation_auc_values.csv"),
  file.path(
    out_intermediate_dir,
    "bulk_response_sensitivity_summary.rds"
  ),
  file.path(
    log_dir,
    "sessionInfo_05_bulk_response_sensitivity_analysis.txt"
  )
)

output_inventory <- data.frame(
  File = basename(canonical_outputs),
  RelativePath = vapply(
    canonical_outputs,
    relative_to_project,
    character(1)
  ),
  Exists = file.exists(canonical_outputs),
  SizeBytes = ifelse(
    file.exists(canonical_outputs),
    file.info(canonical_outputs)$size,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  output_inventory,
  file.path(out_dir, "05_output_inventory.csv")
)

message("============================================================")
message("05_bulk_response_sensitivity_analysis.R completed.")
message("S11 fragility module: completed.")
message("S12 composition-aware audit: completed.")
message("S13 formal composition-adjusted module: completed.")
message("Submitted S11/S12 legacy bootstrap CIs frozen with explicit provenance.")
message("Seeded S11/S12 CI reanalysis outputs retained separately for audit.")
message("All public-release reproducibility gates passed.")
message("============================================================")
