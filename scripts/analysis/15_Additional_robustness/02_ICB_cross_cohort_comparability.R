############################################################
## 02_ICB_cross_cohort_comparability.R
##
## Project: ICB_resistance_project
##
## ROLE
##   Canonical public-repository analysis script for the cross-cohort
##   response-comparability sensitivity analysis combining GSE78220
##   and pretreatment GSE91061 in a common expression space.
##
## CANONICAL CONTRACT
##   1) GSE78220: strict binary response subset, n = 27.
##   2) GSE91061: pretreatment strict CR/PR-vs-PD subset, n = 33.
##   3) Combined response-analysis universe: n = 60.
##   4) Common HGNC-symbol universe: exactly 19,799 genes.
##   5) Batch harmonization:
##        sva::ComBat(raw common-expression matrix,
##                    batch = cohort,
##                    mod = NULL,
##                    par.prior = TRUE)
##      Response is NOT included in the ComBat model.
##   6) Post-ComBat gene-wise z-standardization across the 60 samples.
##   7) Four-state scoring uses the same frozen GMT/signature mapping
##      as scripts 02 and 04.
##   8) AUC direction is frozen prospectively:
##        higher state score -> greater non-response tendency
##      with NonResponse = 1 and fixed ROC direction.
##   9) Deterministic DeLong 95% CI for the point AUC.
##  10) Seeded sensitivity layers:
##        - 2,000 bootstrap resamples, stratified by cohort x response.
##        - 5,000 within-cohort response-label permutations.
##        - permutation test is two-sided around AUC = 0.5.
##
## HISTORICAL 18H PROVENANCE
##   The historical 18H expression/ComBat/row-z/state-score lineage was
##   previously source-locked and numerically audited.
##
##   One historical endpoint inconsistency was identified:
##     Immune_defective_Cold historical pooled AUC = 0.437714...
##     canonical fixed-direction AUC               = 0.562286...
##     and the historical CI is the exact reflected interval.
##
##   Therefore the historical 18H AUC/bootstrap/permutation summary is
##   retained only as a LEGACY provenance snapshot. It is not used as
##   the canonical inferential output of this public script.
##
## OUTPUT
##   results/tables/cross_cohort_harmonization/
##   results/intermediate/cross_cohort_harmonization/
##   logs/
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################################################
## 0. Project root
############################################################

resolve_project_dir <- function() {
  env <- Sys.getenv("ICB_PROJECT_DIR", unset = "")
  if (nzchar(env) && dir.exists(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }
  
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (dir.exists(file.path(wd, "scripts")) &&
      dir.exists(file.path(wd, "results"))) {
    return(wd)
  }
  
  stop(
    "Cannot determine project root. Set:\n",
    'Sys.setenv(ICB_PROJECT_DIR = "D:/ICB_resistance_project")'
  )
}

project_dir <- resolve_project_dir()

table_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "cross_cohort_harmonization"
)
intermediate_dir <- file.path(
  project_dir,
  "results",
  "intermediate",
  "cross_cohort_harmonization"
)
log_dir <- file.path(project_dir, "logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", project_dir)

############################################################
## 1. Packages
############################################################

required_pkgs <- c("sva", "pROC")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Restore the frozen project environment before running this cross-cohort analysis."
  )
}

############################################################
## 2. Frozen parameters
############################################################

EXPECTED_G78220_BINARY_N <- 27L
EXPECTED_G91061_BINARY_N <- 33L
EXPECTED_COMBINED_N <- 60L
EXPECTED_COMMON_GENE_N <- 19799L

STATE_COLS <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

STATE_DISPLAY <- c(
  "Immune_defective_Cold" =
    "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" =
    "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" =
    "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" =
    "melanocytic differentiation"
)

## Frozen state-score signature definitions.
IMMUNE_SIGNATURES <- c(
  "Immune_defective_Cold_RESTORE_hybrid",
  "Immune_defective_Cold_RESTORE_data_driven",
  "Immune_defective_Cold_RESTORE_curated"
)

MYELOID_SIGNATURES <- c(
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated"
)

TUMOR_SIGNATURES <- c(
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated"
)

MEL_SIGNATURES <- c(
  "Melanocytic_Differentiation_REFERENCE_hybrid",
  "Melanocytic_Differentiation_REFERENCE_data_driven",
  "Melanocytic_Differentiation_REFERENCE_curated"
)

## Deterministic public-release RNG plan.
BOOT_N <- 2000L
BOOT_SEED_BASE <- 202608290L

PERM_N <- 5000L
PERM_SEED_BASE <- 202608390L

## Canonical deterministic point estimates established in 19B.
EXPECTED_CANONICAL_AUC <- c(
  "Immune_defective_Cold" = 0.5622857142857143,
  "Myeloid_Treg_Immunosuppressive" = 0.4720000000000000,
  "Tumor_dedifferentiation_Stromal_remodeling" = 0.6057142857142858,
  "Melanocytic_Differentiation" = 0.4548571428571429
)

############################################################
## 3. Helpers
############################################################

resolve_existing_file <- function(paths) {
  for (p in paths) {
    if (file.exists(p)) {
      return(p)
    }
  }
  # Fallback to the first path so upstream logic can report the missing primary target
  return(paths[1])
}

safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    x,
    file = file,
    row.names = FALSE,
    na = ""
  )
  message(
    "Saved: ",
    normalizePath(file, winslash = "/", mustWork = FALSE)
  )
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
  message(
    "Saved: ",
    normalizePath(file, winslash = "/", mustWork = FALSE)
  )
}

md5_file <- function(file) {
  if (!file.exists(file)) return(NA_character_)
  unname(tools::md5sum(file))
}

assert_true <- function(x, msg) {
  if (!isTRUE(x)) stop(msg, call. = FALSE)
}

zvec <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  
  (x - m) / s
}

row_z <- function(m) {
  m <- as.matrix(m)
  mu <- rowMeans(m, na.rm = TRUE)
  ss <- apply(m, 1, stats::sd, na.rm = TRUE)
  
  out <- sweep(m, 1, mu, "-")
  good <- is.finite(ss) & ss > 0
  
  out[good, ] <- sweep(
    out[good, , drop = FALSE],
    1,
    ss[good],
    "/"
  )
  
  out[!good, ] <- 0
  out
}

read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  
  gs <- list()
  
  for (ln in lines) {
    sp <- strsplit(ln, "\t")[[1]]
    if (length(sp) < 3L) next
    
    nm <- sp[[1]]
    genes <- unique(
      toupper(
        trimws(
          sp[-c(1, 2)]
        )
      )
    )
    genes <- genes[nzchar(genes)]
    
    gs[[nm]] <- genes
  }
  
  gs
}

auc_exact_components <- function(y, score) {
  keep <- !is.na(y) & is.finite(score)
  
  y <- as.integer(y[keep])
  score <- as.numeric(score[keep])
  
  assert_true(
    identical(sort(unique(y)), c(0L, 1L)),
    "AUC requires both binary outcome classes."
  )
  
  pos <- score[y == 1L]
  neg <- score[y == 0L]
  
  cmp <- outer(pos, neg, "-")
  
  n_pairs <- length(pos) * length(neg)
  
  ## Twice the Mann-Whitney concordance numerator:
  ##   2 * concordant pairs + tied pairs.
  ## This integer-valued representation avoids floating-point ambiguity
  ## when a permutation statistic lies exactly on the observed AUC
  ## distance from 0.5.
  concordance_numerator_x2 <-
    2L * sum(cmp > 0) +
    sum(cmp == 0)
  
  auc <-
    concordance_numerator_x2 /
    (2 * n_pairs)
  
  ## Exact two-sided distance from AUC = 0.5 on the same integer scale.
  centered_stat_x2 <-
    abs(
      concordance_numerator_x2 -
        n_pairs
    )
  
  list(
    auc = auc,
    n_pairs = n_pairs,
    concordance_numerator_x2 =
      concordance_numerator_x2,
    centered_stat_x2 =
      centered_stat_x2
  )
}

rank_auc <- function(y, score) {
  auc_exact_components(y, score)$auc
}

bootstrap_auc_cohort_response <- function(
    df,
    score_col,
    n_boot,
    seed
) {
  set.seed(seed)
  
  strata <- interaction(
    df$Cohort,
    df$ResponseGroup,
    drop = TRUE
  )
  
  strata_index <- split(
    seq_len(nrow(df)),
    strata
  )
  
  out <- rep(NA_real_, n_boot)
  
  for (b in seq_len(n_boot)) {
    idx <- unlist(
      lapply(
        strata_index,
        function(ii) {
          sample(ii, length(ii), replace = TRUE)
        }
      ),
      use.names = FALSE
    )
    
    out[b] <- rank_auc(
      y = df$NonResponse[idx],
      score = df[[score_col]][idx]
    )
  }
  
  out
}

permutation_auc_within_cohort <- function(
    df,
    score_col,
    n_perm,
    seed
) {
  set.seed(seed)
  
  obs_comp <- auc_exact_components(
    df$NonResponse,
    df[[score_col]]
  )
  
  cohort_index <- split(
    seq_len(nrow(df)),
    df$Cohort
  )
  
  perm_auc <- rep(NA_real_, n_perm)
  perm_centered_stat_x2 <- rep(NA_integer_, n_perm)
  
  for (i in seq_len(n_perm)) {
    y_perm <- df$NonResponse
    
    for (ii in cohort_index) {
      y_perm[ii] <- sample(
        y_perm[ii],
        replace = FALSE
      )
    }
    
    perm_comp <- auc_exact_components(
      y_perm,
      df[[score_col]]
    )
    
    ## The within-cohort permutation preserves the global class counts,
    ## so the Mann-Whitney pair denominator must remain identical.
    assert_true(
      perm_comp$n_pairs == obs_comp$n_pairs,
      "Permutation changed the AUC pair denominator unexpectedly."
    )
    
    perm_auc[i] <-
      perm_comp$auc
    
    perm_centered_stat_x2[i] <-
      perm_comp$centered_stat_x2
  }
  
  ## Exact discrete two-sided permutation test around AUC = 0.5.
  ## Compare the integer Mann-Whitney centered statistic rather than
  ## floating-point abs(AUC - 0.5), so ties at the observed boundary
  ## are counted reproducibly.
  n_extreme <- sum(
    perm_centered_stat_x2 >=
      obs_comp$centered_stat_x2
  )
  
  p_two_sided <- (
    1 + n_extreme
  ) /
    (n_perm + 1)
  
  list(
    observed = obs_comp$auc,
    observed_n_pairs =
      obs_comp$n_pairs,
    observed_concordance_numerator_x2 =
      obs_comp$concordance_numerator_x2,
    observed_centered_stat_x2 =
      obs_comp$centered_stat_x2,
    permutation = perm_auc,
    permutation_centered_stat_x2 =
      perm_centered_stat_x2,
    n_extreme = n_extreme,
    p_two_sided = p_two_sided
  )
}

############################################################
## 4. Canonical upstream inputs
############################################################

g782_expr_file <- resolve_existing_file(c(
  file.path(
    project_dir,
    "data_processed",
    "GSE78220_external_expression_matrix_for_scoring.rds"
  ),
  file.path(
    project_dir,
    "data_processed",
    "GSE78220_external_expression_matrix_for_scoring_CLEAN.rds"
  )
))

g782_state_file <- resolve_existing_file(c(
  file.path(
    project_dir,
    "results",
    "tables",
    "external_validation",
    "GSE78220_external_state_scores.csv"
  ),
  file.path(
    project_dir,
    "results",
    "tables",
    "external_validation",
    "GSE78220_external_state_scores_CLEAN_STANDARD.csv"
  )
))

g910_expr_file <- file.path(
  project_dir,
  "results",
  "intermediate",
  "GSE91061",
  "GSE91061_expression_symbol_matrix_for_state_scoring.rds"
)

g910_pre_file <- file.path(
  project_dir,
  "results",
  "tables",
  "GSE91061",
  "GSE91061_state_scores_pretreatment.csv"
)

gmt_file <- file.path(
  project_dir,
  "results",
  "tables",
  "ICBcomb_final_input_gene_sets_CLEAN",
  "02_final_ICBcomb_gene_sets",
  "final_ICBcomb_gene_sets_CLEAN.gmt"
)

## Historical 18H is audit-only.
legacy_root <- file.path(
  project_dir,
  "results",
  "revision_cross_cohort_comparability_18H"
)

legacy_common_file <- file.path(
  legacy_root,
  "intermediate",
  "18H_common_expression_space_ComBat_harmonized.rds"
)

legacy_rowz_file <- file.path(
  legacy_root,
  "intermediate",
  "18H_ComBat_row_z_expression_matrix.rds"
)

legacy_score_file <- file.path(
  legacy_root,
  "tables",
  "18H_harmonized_four_state_scores.csv"
)

legacy_auc_file <- file.path(
  legacy_root,
  "tables",
  "18H_harmonized_AUC_bootstrap_permutation_summary.csv"
)

required_inputs <- c(
  g782_expr_file,
  g782_state_file,
  g910_expr_file,
  g910_pre_file,
  gmt_file
)

if (any(!file.exists(required_inputs))) {
  stop(
    "Missing canonical input(s):\n",
    paste(
      required_inputs[!file.exists(required_inputs)],
      collapse = "\n"
    )
  )
}

input_manifest <- data.frame(
  Input = c(
    "GSE78220 canonical expression",
    "GSE78220 canonical state/response table",
    "GSE91061 canonical expression",
    "GSE91061 pretreatment state/response table",
    "Frozen CLEAN GMT"
  ),
  Path = required_inputs,
  MD5 = vapply(
    required_inputs,
    md5_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_manifest,
  file.path(
    table_dir,
    "cross_cohort_input_manifest.csv"
  )
)

############################################################
## 5. Read canonical upstream objects
############################################################

g782_expr <- readRDS(g782_expr_file)
g910_expr <- readRDS(g910_expr_file)

if (!is.matrix(g782_expr)) {
  g782_expr <- as.matrix(g782_expr)
}
if (!is.matrix(g910_expr)) {
  g910_expr <- as.matrix(g910_expr)
}

storage.mode(g782_expr) <- "numeric"
storage.mode(g910_expr) <- "numeric"

g782_state <- utils::read.csv(
  g782_state_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

g910_pre <- utils::read.csv(
  g910_pre_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

############################################################
## 6. Freeze the response-analysis sample universe
############################################################

assert_true(
  all(
    c(
      "Sample",
      "ResponseGroup"
    ) %in% colnames(g782_state)
  ),
  "GSE78220 state table is missing Sample/ResponseGroup."
)

assert_true(
  all(
    c(
      "Sample",
      "BinaryResponseGroup"
    ) %in% colnames(g910_pre)
  ),
  "GSE91061 pretreatment table is missing Sample/BinaryResponseGroup."
)

g782_bin <- g782_state[
  g782_state$ResponseGroup %in%
    c("Responder", "NonResponder"),
  ,
  drop = FALSE
]

g910_bin <- g910_pre[
  g910_pre$BinaryResponseGroup %in%
    c("Responder", "NonResponder"),
  ,
  drop = FALSE
]

assert_true(
  nrow(g782_bin) == EXPECTED_G78220_BINARY_N,
  paste0(
    "GSE78220 strict-binary universe drift: observed ",
    nrow(g782_bin),
    "; expected ",
    EXPECTED_G78220_BINARY_N
  )
)

assert_true(
  nrow(g910_bin) == EXPECTED_G91061_BINARY_N,
  paste0(
    "GSE91061 strict-binary pretreatment universe drift: observed ",
    nrow(g910_bin),
    "; expected ",
    EXPECTED_G91061_BINARY_N
  )
)

sample_manifest <- rbind(
  data.frame(
    Sample = as.character(g782_bin$Sample),
    Cohort = "GSE78220",
    ResponseGroup = as.character(g782_bin$ResponseGroup),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Sample = as.character(g910_bin$Sample),
    Cohort = "GSE91061",
    ResponseGroup = as.character(g910_bin$BinaryResponseGroup),
    stringsAsFactors = FALSE
  )
)

sample_manifest$NonResponse <- ifelse(
  sample_manifest$ResponseGroup == "NonResponder",
  1L,
  0L
)

assert_true(
  nrow(sample_manifest) == EXPECTED_COMBINED_N,
  "Combined response-analysis universe is not n=60."
)

assert_true(
  !anyDuplicated(sample_manifest$Sample),
  "Duplicate sample IDs detected across the combined n=60 universe."
)

sample_summary <- do.call(
  rbind,
  lapply(
    split(sample_manifest, sample_manifest$Cohort),
    function(x) {
      data.frame(
        Cohort = unique(x$Cohort),
        N = nrow(x),
        Responders =
          sum(x$ResponseGroup == "Responder"),
        NonResponders =
          sum(x$ResponseGroup == "NonResponder"),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(sample_summary) <- NULL

safe_write_csv(
  sample_manifest,
  file.path(
    table_dir,
    "cross_cohort_sample_manifest.csv"
  )
)

safe_write_csv(
  sample_summary,
  file.path(
    table_dir,
    "cross_cohort_sample_summary.csv"
  )
)

############################################################
## 7. Exact common-gene universe
############################################################

common_genes <- intersect(
  rownames(g782_expr),
  rownames(g910_expr)
)

assert_true(
  length(common_genes) == EXPECTED_COMMON_GENE_N,
  paste0(
    "Common gene universe drift: observed ",
    length(common_genes),
    "; expected ",
    EXPECTED_COMMON_GENE_N
  )
)

common_gene_manifest <- data.frame(
  Gene = common_genes,
  stringsAsFactors = FALSE
)

safe_write_csv(
  common_gene_manifest,
  file.path(
    table_dir,
    "cross_cohort_common_gene_manifest.csv"
  )
)

############################################################
## 8. Construct exact 19,799 x 60 raw matrix
############################################################

raw_common <- matrix(
  NA_real_,
  nrow = length(common_genes),
  ncol = nrow(sample_manifest),
  dimnames = list(
    common_genes,
    sample_manifest$Sample
  )
)

idx782 <- which(
  sample_manifest$Cohort == "GSE78220"
)

idx910 <- which(
  sample_manifest$Cohort == "GSE91061"
)

raw_common[, idx782] <-
  g782_expr[
    common_genes,
    sample_manifest$Sample[idx782],
    drop = FALSE
  ]

raw_common[, idx910] <-
  g910_expr[
    common_genes,
    sample_manifest$Sample[idx910],
    drop = FALSE
  ]

assert_true(
  identical(
    dim(raw_common),
    c(
      EXPECTED_COMMON_GENE_N,
      EXPECTED_COMBINED_N
    )
  ),
  "Raw common-expression matrix dimension gate failed."
)

############################################################
## 9. Canonical ComBat: cohort batch only
############################################################

batch <- factor(
  sample_manifest$Cohort,
  levels = c(
    "GSE78220",
    "GSE91061"
  )
)

combat_common <- sva::ComBat(
  dat = raw_common,
  batch = batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE
)

assert_true(
  identical(
    dim(combat_common),
    c(
      EXPECTED_COMMON_GENE_N,
      EXPECTED_COMBINED_N
    )
  ),
  "ComBat matrix dimension gate failed."
)

combat_rowz <- row_z(combat_common)

############################################################
## 10. State scoring
############################################################

gene_sets <- read_gmt(gmt_file)

get_present <- function(signature_name) {
  assert_true(
    signature_name %in% names(gene_sets),
    paste0(
      "Frozen signature missing from GMT: ",
      signature_name
    )
  )
  
  intersect(
    gene_sets[[signature_name]],
    rownames(combat_rowz)
  )
}

all_signatures <- c(
  IMMUNE_SIGNATURES,
  MYELOID_SIGNATURES,
  TUMOR_SIGNATURES,
  MEL_SIGNATURES
)

raw_signature_scores <- data.frame(
  Sample = colnames(combat_rowz),
  stringsAsFactors = FALSE
)

coverage_rows <- list()

for (sig in all_signatures) {
  total_genes <- gene_sets[[sig]]
  present_genes <- get_present(sig)
  
  assert_true(
    length(present_genes) >= 3L,
    paste0(
      "Too few common genes for signature ",
      sig,
      ": ",
      length(present_genes)
    )
  )
  
  raw_signature_scores[[sig]] <-
    colMeans(
      combat_rowz[
        present_genes,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  
  coverage_rows[[sig]] <- data.frame(
    Signature = sig,
    N_total = length(total_genes),
    N_present_common = length(present_genes),
    CoverageFraction =
      length(present_genes) /
      length(total_genes),
    stringsAsFactors = FALSE
  )
}

signature_coverage <- do.call(
  rbind,
  coverage_rows
)
rownames(signature_coverage) <- NULL

safe_write_csv(
  signature_coverage,
  file.path(
    table_dir,
    "cross_cohort_signature_gene_coverage.csv"
  )
)

state_scores <- data.frame(
  Sample =
    raw_signature_scores$Sample,
  
  Immune_defective_Cold =
    zvec(
      -rowMeans(
        raw_signature_scores[
          ,
          IMMUNE_SIGNATURES,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    ),
  
  Myeloid_Treg_Immunosuppressive =
    zvec(
      rowMeans(
        raw_signature_scores[
          ,
          MYELOID_SIGNATURES,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    ),
  
  Tumor_dedifferentiation_Stromal_remodeling =
    zvec(
      rowMeans(
        raw_signature_scores[
          ,
          TUMOR_SIGNATURES,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    ),
  
  Melanocytic_Differentiation =
    zvec(
      rowMeans(
        raw_signature_scores[
          ,
          MEL_SIGNATURES,
          drop = FALSE
        ],
        na.rm = TRUE
      )
    ),
  
  stringsAsFactors = FALSE
)

final_df <- merge(
  sample_manifest,
  state_scores,
  by = "Sample",
  all.x = TRUE,
  sort = FALSE
)

final_df <- final_df[
  match(
    sample_manifest$Sample,
    final_df$Sample
  ),
  ,
  drop = FALSE
]

assert_true(
  nrow(final_df) == EXPECTED_COMBINED_N,
  "Final state/response table must contain exactly 60 rows."
)

for (st in STATE_COLS) {
  assert_true(
    all(is.finite(final_df[[st]])),
    paste0(
      "Non-finite state score detected for ",
      st
    )
  )
}

safe_write_csv(
  raw_signature_scores,
  file.path(
    table_dir,
    "cross_cohort_raw_signature_scores.csv"
  )
)

safe_write_csv(
  state_scores,
  file.path(
    table_dir,
    "cross_cohort_harmonized_four_state_scores.csv"
  )
)

safe_write_csv(
  final_df,
  file.path(
    table_dir,
    "cross_cohort_harmonized_state_scores_with_response.csv"
  )
)

############################################################
## 11. Canonical fixed-direction AUC + DeLong CI
############################################################

auc_rows <- list()

for (i in seq_along(STATE_COLS)) {
  st <- STATE_COLS[[i]]
  
  y <- final_df$NonResponse
  x <- final_df[[st]]
  
  auc_rank <- rank_auc(y, x)
  
  roc_obj <- pROC::roc(
    response = y,
    predictor = x,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  auc_proc <- as.numeric(
    pROC::auc(roc_obj)
  )
  
  assert_true(
    abs(auc_rank - auc_proc) <= 1e-12,
    paste0(
      "Rank AUC and pROC fixed-direction AUC disagree for ",
      st
    )
  )
  
  assert_true(
    abs(
      auc_rank -
        unname(
          EXPECTED_CANONICAL_AUC[[st]]
        )
    ) <= 1e-12,
    paste0(
      "Canonical AUC drift for ",
      st,
      ": observed ",
      signif(auc_rank, 15),
      "; expected ",
      signif(
        EXPECTED_CANONICAL_AUC[[st]],
        15
      )
    )
  )
  
  ci <- as.numeric(
    pROC::ci.auc(
      roc_obj,
      conf.level = 0.95,
      method = "delong"
    )
  )
  
  auc_rows[[st]] <- data.frame(
    State = st,
    StateDisplay = unname(
      STATE_DISPLAY[[st]]
    ),
    N = nrow(final_df),
    Responders =
      sum(
        final_df$ResponseGroup ==
          "Responder"
      ),
    NonResponders =
      sum(
        final_df$ResponseGroup ==
          "NonResponder"
      ),
    AUC = auc_rank,
    CI_method = "DeLong",
    CI_low = ci[[1]],
    CI_high = ci[[3]],
    Direction =
      "higher score -> greater non-response tendency",
    OutcomeCoding =
      "NonResponse=1; Responder=0",
    stringsAsFactors = FALSE
  )
}

auc_summary <- do.call(
  rbind,
  auc_rows
)
rownames(auc_summary) <- NULL

safe_write_csv(
  auc_summary,
  file.path(
    table_dir,
    "cross_cohort_fixed_direction_AUC_DeLong_summary.csv"
  )
)

############################################################
## 12. Seeded bootstrap sensitivity
##
## Stratification preserves BOTH cohort and response composition.
############################################################

bootstrap_values <- list()
bootstrap_summary <- list()

for (i in seq_along(STATE_COLS)) {
  st <- STATE_COLS[[i]]
  seed_i <- BOOT_SEED_BASE + i
  
  boot <- bootstrap_auc_cohort_response(
    df = final_df,
    score_col = st,
    n_boot = BOOT_N,
    seed = seed_i
  )
  
  bootstrap_values[[st]] <- data.frame(
    State = st,
    Iteration = seq_len(BOOT_N),
    AUC = boot,
    BootstrapSeed = seed_i,
    stringsAsFactors = FALSE
  )
  
  bootstrap_summary[[st]] <- data.frame(
    State = st,
    StateDisplay =
      unname(
        STATE_DISPLAY[[st]]
      ),
    N_boot = BOOT_N,
    BootstrapSeed = seed_i,
    Mean_AUC =
      mean(
        boot,
        na.rm = TRUE
      ),
    Median_AUC =
      stats::median(
        boot,
        na.rm = TRUE
      ),
    Percentile_low =
      as.numeric(
        stats::quantile(
          boot,
          0.025,
          na.rm = TRUE
        )
      ),
    Percentile_high =
      as.numeric(
        stats::quantile(
          boot,
          0.975,
          na.rm = TRUE
        )
      ),
    BootstrapStrata =
      "Cohort x ResponseGroup",
    AUCDirection =
      "fixed: higher score -> greater non-response tendency",
    stringsAsFactors = FALSE
  )
}

bootstrap_values <- do.call(
  rbind,
  bootstrap_values
)
bootstrap_summary <- do.call(
  rbind,
  bootstrap_summary
)

rownames(bootstrap_values) <- NULL
rownames(bootstrap_summary) <- NULL

safe_write_csv(
  bootstrap_values,
  file.path(
    table_dir,
    "cross_cohort_seeded_bootstrap_AUC_values.csv"
  )
)

safe_write_csv(
  bootstrap_summary,
  file.path(
    table_dir,
    "cross_cohort_seeded_bootstrap_AUC_summary.csv"
  )
)

############################################################
## 13. Seeded within-cohort permutation sensitivity
##
## Permutation preserves cohort membership and cohort-specific
## response prevalence. The test is two-sided around AUC = 0.5.
############################################################

perm_values <- list()
perm_summary <- list()

for (i in seq_along(STATE_COLS)) {
  st <- STATE_COLS[[i]]
  seed_i <- PERM_SEED_BASE + i
  
  pp <- permutation_auc_within_cohort(
    df = final_df,
    score_col = st,
    n_perm = PERM_N,
    seed = seed_i
  )
  
  perm_values[[st]] <- data.frame(
    State = st,
    Iteration = seq_len(PERM_N),
    PermutationAUC = pp$permutation,
    PermutationCenteredStatX2 =
      pp$permutation_centered_stat_x2,
    ObservedAUC = pp$observed,
    ObservedCenteredStatX2 =
      pp$observed_centered_stat_x2,
    ObservedConcordanceNumeratorX2 =
      pp$observed_concordance_numerator_x2,
    AUC_pair_denominator =
      pp$observed_n_pairs,
    PermutationSeed = seed_i,
    stringsAsFactors = FALSE
  )
  
  perm_summary[[st]] <- data.frame(
    State = st,
    StateDisplay =
      unname(
        STATE_DISPLAY[[st]]
      ),
    Observed_AUC =
      pp$observed,
    Permutation_mean_AUC =
      mean(
        pp$permutation,
        na.rm = TRUE
      ),
    Permutation_sd_AUC =
      stats::sd(
        pp$permutation,
        na.rm = TRUE
      ),
    P_two_sided_from_0_5 =
      pp$p_two_sided,
    N_extreme_exact =
      pp$n_extreme,
    ObservedCenteredStatX2 =
      pp$observed_centered_stat_x2,
    AUC_pair_denominator =
      pp$observed_n_pairs,
    N_perm = PERM_N,
    PermutationSeed = seed_i,
    PermutationScheme =
      "response labels permuted within cohort; exact discrete two-sided Mann-Whitney distance from AUC=0.5",
    stringsAsFactors = FALSE
  )
}

perm_values <- do.call(
  rbind,
  perm_values
)
perm_summary <- do.call(
  rbind,
  perm_summary
)

rownames(perm_values) <- NULL
rownames(perm_summary) <- NULL

safe_write_csv(
  perm_values,
  file.path(
    table_dir,
    "cross_cohort_seeded_within_cohort_permutation_values.csv"
  )
)

safe_write_csv(
  perm_summary,
  file.path(
    table_dir,
    "cross_cohort_seeded_within_cohort_permutation_summary.csv"
  )
)

############################################################
## 14. Historical 18H equivalence / discrepancy audit
############################################################

legacy_available <- all(
  file.exists(
    c(
      legacy_common_file,
      legacy_rowz_file,
      legacy_score_file,
      legacy_auc_file
    )
  )
)

legacy_audit <- data.frame()

if (legacy_available) {
  legacy_common <- readRDS(
    legacy_common_file
  )
  legacy_rowz <- readRDS(
    legacy_rowz_file
  )
  
  legacy_scores <- utils::read.csv(
    legacy_score_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  legacy_auc <- utils::read.csv(
    legacy_auc_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  ## Align matrices by exact names.
  combat_aligned <- combat_common[
    rownames(legacy_common),
    colnames(legacy_common),
    drop = FALSE
  ]
  
  rowz_aligned <- combat_rowz[
    rownames(legacy_rowz),
    colnames(legacy_rowz),
    drop = FALSE
  ]
  
  combat_diff <- max(
    abs(
      combat_aligned -
        legacy_common
    ),
    na.rm = TRUE
  )
  
  rowz_diff <- max(
    abs(
      rowz_aligned -
        legacy_rowz
    ),
    na.rm = TRUE
  )
  
  score_diffs <- vapply(
    STATE_COLS,
    function(st) {
      ls <- legacy_scores[
        match(
          state_scores$Sample,
          legacy_scores$Sample
        ),
        st
      ]
      
      max(
        abs(
          state_scores[[st]] -
            ls
        ),
        na.rm = TRUE
      )
    },
    numeric(1)
  )
  
  legacy_auc_cmp <- merge(
    auc_summary[
      ,
      c(
        "State",
        "AUC",
        "CI_low",
        "CI_high"
      )
    ],
    legacy_auc[
      ,
      intersect(
        c(
          "State",
          "AUC",
          "CI_low",
          "CI_high",
          "Bootstrap_mean",
          "Bootstrap_low",
          "Bootstrap_high",
          "Permutation_P",
          "Permutation_N"
        ),
        colnames(legacy_auc)
      )
    ],
    by = "State",
    suffixes = c(
      "_Canonical",
      "_Historical"
    ),
    all.x = TRUE
  )
  
  legacy_auc_cmp$HistoricalPointRelationship <-
    ifelse(
      abs(
        legacy_auc_cmp$AUC_Canonical -
          legacy_auc_cmp$AUC_Historical
      ) <= 1e-12,
      "exact",
      ifelse(
        abs(
          legacy_auc_cmp$AUC_Historical -
            (1 - legacy_auc_cmp$AUC_Canonical)
        ) <= 1e-12,
        "exact complement (direction reversal)",
        "other discrepancy"
      )
    )
  
  safe_write_csv(
    legacy_auc_cmp,
    file.path(
      table_dir,
      "cross_cohort_historical_18H_AUC_discrepancy_audit.csv"
    )
  )
  
  legacy_audit <- data.frame(
    Check = c(
      "Historical ComBat matrix exact",
      "Historical row-z matrix exact",
      paste0(
        "Historical state scores exact: ",
        STATE_COLS
      )
    ),
    MaxAbsoluteDifference = c(
      combat_diff,
      rowz_diff,
      unname(score_diffs)
    ),
    Pass = c(
      combat_diff <= 1e-10,
      rowz_diff <= 1e-10,
      score_diffs <= 1e-10
    ),
    stringsAsFactors = FALSE
  )
  
  safe_write_csv(
    legacy_audit,
    file.path(
      table_dir,
      "cross_cohort_historical_18H_expression_score_equivalence.csv"
    )
  )
}

############################################################
## 15. Save canonical intermediate objects
############################################################

safe_save_rds(
  raw_common,
  file.path(
    intermediate_dir,
    "cross_cohort_common_expression_raw.rds"
  )
)

safe_save_rds(
  combat_common,
  file.path(
    intermediate_dir,
    "cross_cohort_common_expression_ComBat_batch_only.rds"
  )
)

safe_save_rds(
  combat_rowz,
  file.path(
    intermediate_dir,
    "cross_cohort_ComBat_gene_row_z.rds"
  )
)

safe_save_rds(
  final_df,
  file.path(
    intermediate_dir,
    "cross_cohort_final_analysis_table.rds"
  )
)

############################################################
## 16. Canonical method manifest
############################################################

method_manifest <- data.frame(
  Parameter = c(
    "Role",
    "GSE78220 response universe",
    "GSE91061 response universe",
    "Combined sample universe",
    "Common-gene universe",
    "Batch factor",
    "ComBat model covariates",
    "ComBat par.prior",
    "Post-ComBat scaling",
    "State-score mapping",
    "AUC endpoint",
    "AUC direction",
    "AUC CI",
    "Bootstrap",
    "Permutation",
    "Historical 18H role"
  ),
  Value = c(
    "Cross-cohort response-comparability sensitivity analysis",
    "27 strict Responder/NonResponder samples",
    "33 pretreatment strict CR/PR-vs-PD samples",
    "60 samples",
    "19,799 exact common HGNC symbols",
    "Cohort: GSE78220 vs GSE91061",
    "None; response status is NOT included",
    "TRUE",
    "Gene-wise z-standardization across all 60 samples",
    "Frozen scripts 02/04 CLEAN GMT mapping",
    "Non-response",
    "Higher state score -> greater non-response tendency",
    "DeLong 95% CI",
    paste0(
      BOOT_N,
      " seeded resamples, stratified by cohort x response"
    ),
    paste0(
      PERM_N,
      " seeded permutations within cohort; exact discrete two-sided Mann-Whitney distance around AUC=0.5"
    ),
    "Legacy provenance only; not canonical inferential authority"
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  method_manifest,
  file.path(
    table_dir,
    "cross_cohort_method_manifest.csv"
  )
)

############################################################
## 17. Reproducibility gates
############################################################

gate_rows <- list()

add_gate <- function(
    gate,
    observed,
    expected,
    pass
) {
  gate_rows[[length(gate_rows) + 1L]] <<-
    data.frame(
      Gate = gate,
      Observed = as.character(observed),
      Expected = as.character(expected),
      Pass = isTRUE(pass),
      stringsAsFactors = FALSE
    )
}

add_gate(
  "GSE78220_binary_n",
  nrow(g782_bin),
  EXPECTED_G78220_BINARY_N,
  nrow(g782_bin) ==
    EXPECTED_G78220_BINARY_N
)

add_gate(
  "GSE91061_binary_n",
  nrow(g910_bin),
  EXPECTED_G91061_BINARY_N,
  nrow(g910_bin) ==
    EXPECTED_G91061_BINARY_N
)

add_gate(
  "combined_n",
  nrow(final_df),
  EXPECTED_COMBINED_N,
  nrow(final_df) ==
    EXPECTED_COMBINED_N
)

add_gate(
  "common_gene_n",
  length(common_genes),
  EXPECTED_COMMON_GENE_N,
  length(common_genes) ==
    EXPECTED_COMMON_GENE_N
)

add_gate(
  "all_state_scores_finite",
  all(
    vapply(
      STATE_COLS,
      function(st) {
        all(
          is.finite(
            final_df[[st]]
          )
        )
      },
      logical(1)
    )
  ),
  TRUE,
  all(
    vapply(
      STATE_COLS,
      function(st) {
        all(
          is.finite(
            final_df[[st]]
          )
        )
      },
      logical(1)
    )
  )
)

auc_gate <- all(
  abs(
    auc_summary$AUC -
      unname(
        EXPECTED_CANONICAL_AUC[
          auc_summary$State
        ]
      )
  ) <= 1e-12
)

add_gate(
  "four_fixed_direction_AUCs",
  paste(
    signif(
      auc_summary$AUC,
      8
    ),
    collapse = ";"
  ),
  paste(
    signif(
      unname(
        EXPECTED_CANONICAL_AUC[
          auc_summary$State
        ]
      ),
      8
    ),
    collapse = ";"
  ),
  auc_gate
)

add_gate(
  "permutation_exact_integer_statistics",
  all(
    perm_summary$N_extreme_exact >= 0 &
      perm_summary$N_extreme_exact <= PERM_N &
      perm_summary$ObservedCenteredStatX2 >= 0
  ),
  TRUE,
  all(
    perm_summary$N_extreme_exact >= 0 &
      perm_summary$N_extreme_exact <= PERM_N &
      perm_summary$ObservedCenteredStatX2 >= 0
  )
)

if (legacy_available) {
  add_gate(
    "legacy_ComBat_exact",
    max(
      legacy_audit$MaxAbsoluteDifference[
        legacy_audit$Check ==
          "Historical ComBat matrix exact"
      ]
    ),
    "<=1e-10",
    legacy_audit$Pass[
      legacy_audit$Check ==
        "Historical ComBat matrix exact"
    ]
  )
  
  add_gate(
    "legacy_rowz_exact",
    max(
      legacy_audit$MaxAbsoluteDifference[
        legacy_audit$Check ==
          "Historical row-z matrix exact"
      ]
    ),
    "<=1e-10",
    legacy_audit$Pass[
      legacy_audit$Check ==
        "Historical row-z matrix exact"
    ]
  )
  
  add_gate(
    "legacy_state_scores_4of4_exact",
    sum(
      legacy_audit$Pass[
        grepl(
          "^Historical state scores exact:",
          legacy_audit$Check
        )
      ]
    ),
    "4",
    sum(
      legacy_audit$Pass[
        grepl(
          "^Historical state scores exact:",
          legacy_audit$Check
        )
      ]
    ) == 4L
  )
}

gates <- do.call(
  rbind,
  gate_rows
)

safe_write_csv(
  gates,
  file.path(
    table_dir,
    "cross_cohort_reproducibility_gates.csv"
  )
)

if (any(!gates$Pass)) {
  stop(
    "Cross-cohort reproducibility gate failure: ",
    paste(
      gates$Gate[!gates$Pass],
      collapse = " | "
    )
  )
}

############################################################
## 18. Session information
############################################################

session_file <- file.path(
  log_dir,
  "sessionInfo_15_02_ICB_cross_cohort_comparability.txt"
)

sink(session_file)
print(sessionInfo())
sink()

############################################################
## 19. Output MD5 manifest
############################################################

output_manifest_file <- file.path(
  table_dir,
  "cross_cohort_output_MD5_manifest.csv"
)

output_files <- c(
  list.files(
    table_dir,
    full.names = TRUE,
    recursive = FALSE
  ),
  list.files(
    intermediate_dir,
    full.names = TRUE,
    recursive = FALSE
  ),
  session_file
)

## A file cannot contain a stable MD5 hash of its own final bytes.
## On reruns, an older output manifest may already exist in table_dir;
## exclude it explicitly before computing the new manifest.
output_files <- output_files[
  normalizePath(
    output_files,
    winslash = "/",
    mustWork = FALSE
  ) !=
    normalizePath(
      output_manifest_file,
      winslash = "/",
      mustWork = FALSE
    )
]

assert_true(
  !any(
    basename(output_files) ==
      basename(output_manifest_file)
  ),
  "Output MD5 manifest self-exclusion gate failed."
)

output_manifest <- data.frame(
  File = basename(output_files),
  Path = output_files,
  SizeBytes = file.info(output_files)$size,
  MD5 = vapply(
    output_files,
    md5_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  output_manifest,
  output_manifest_file
)

############################################################
## 20. Console summary
############################################################

cat("\n============================================================\n")
cat("SCRIPT 19 CROSS-COHORT HARMONIZATION: PASS\n")
cat("============================================================\n")
cat("GSE78220 strict-binary n: ", nrow(g782_bin), "\n", sep = "")
cat("GSE91061 strict-binary pretreatment n: ", nrow(g910_bin), "\n", sep = "")
cat("Combined n: ", nrow(final_df), "\n", sep = "")
cat("Common genes: ", length(common_genes), "\n", sep = "")
cat("ComBat model: cohort batch only; response covariate = NO\n")
cat("Canonical fixed-direction AUCs:\n")
print(
  auc_summary[
    ,
    c(
      "State",
      "AUC",
      "CI_low",
      "CI_high"
    )
  ]
)
cat("Seeded bootstrap n/state: ", BOOT_N, "\n", sep = "")
cat("Seeded within-cohort permutation n/state: ", PERM_N, "\n", sep = "")
cat("All reproducibility gates: PASS\n")
cat(
  "Output table directory: ",
  normalizePath(
    table_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")