############################################################
## 11_CosMx_reproduction_check_v1.3.R
##
## CosMx reproduction audit v1.3
##
## Patch from v1.2:
## Only modifies S22 threshold sensitivity comparison.
##
## Reason:
## Step02 uses cutoff_mode = pooled/per_fov and numeric fractions
## (0.20/0.25/0.30), whereas S22 uses Mode = Pooled/Per-FOV
## and TopFraction = 20%/25%/30%.
##
## v1.3 normalizes both Mode and TopFraction, matches rows by
## Mode + TopFraction, and compares OR / lower CI / upper CI.
##
## Step02 sensitivity is calculated on 86600 scored cells;
## S22 is recalculated from the 86572 spatially resolved cells.
## Therefore reporting-level tolerance is 0.02, not bit identity.
##
## No analytical output is changed.
############################################################

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)

PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash="/",
  mustWork=TRUE
)

AUDIT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "audit",
  "11_CosMx_reproduction_check"
)
dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

EXPECTED_RAW_CELLS <- 86600L
EXPECTED_SPATIAL_CELLS <- 86572L

TOL_EXACT <- 1e-12
TOL_REPORT <- 0.02

read_csv0 <- function(path){
  if(!file.exists(path))
    stop("Missing file: ", path, call.=FALSE)

  if(requireNamespace("data.table", quietly=TRUE)){
    return(as.data.frame(
      data.table::fread(
        path,
        data.table=FALSE,
        check.names=FALSE,
        showProgress=FALSE
      )
    ))
  }

  read.csv(
    path,
    check.names=FALSE,
    stringsAsFactors=FALSE
  )
}

maxdiff <- function(a,b){
  a <- suppressWarnings(as.numeric(a))
  b <- suppressWarnings(as.numeric(b))

  if(length(a)!=length(b))
    return(Inf)

  max(
    abs(a-b),
    na.rm=TRUE
  )
}

normalize_mode <- function(x){

  xx <- tolower(
    trimws(
      as.character(x)
    )
  )

  xx <- gsub(
    "-",
    "_",
    xx,
    fixed = TRUE
  )

  xx <- gsub(
    "[[:space:]]+",
    "_",
    xx
  )

  xx <- gsub(
    "_+",
    "_",
    xx
  )

  out <- xx

  out[
    grepl(
      "^pooled$|^global$",
      xx
    )
  ] <- "pooled"

  out[
    grepl(
      "^per_fov$|^perfov$|^fov$",
      xx
    )
  ] <- "per_fov"

  out
}

normalize_fraction <- function(x){

  xx <- trimws(
    as.character(x)
  )

  has_percent <- grepl(
    "%",
    xx,
    fixed = TRUE
  )

  xx <- gsub(
    "%",
    "",
    xx,
    fixed = TRUE
  )

  vv <- suppressWarnings(
    as.numeric(xx)
  )

  vv[
    has_percent &
      is.finite(vv)
  ] <- vv[
    has_percent &
      is.finite(vv)
  ] / 100

  vv[
    !has_percent &
      is.finite(vv) &
      vv > 1
  ] <- vv[
    !has_percent &
      is.finite(vv) &
      vv > 1
  ] / 100

  vv
}

make_key <- function(mode, frac){

  mode2 <- normalize_mode(mode)
  frac2 <- normalize_fraction(frac)

  if (
    any(
      is.na(mode2)
    ) ||
      any(
        !is.finite(frac2)
      )
  ) {
    stop(
      "Could not normalize Mode/TopFraction for S22 comparison.",
      call. = FALSE
    )
  }

  paste0(
    mode2,
    "__",
    sprintf(
      "%.6f",
      frac2
    )
  )
}

############################################################
## Inputs
############################################################

STEP01 <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "IF_mIF_validation",
  "CosMx_cell_state_scores.csv"
)

STEP02 <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "CosMx_spatial_validation",
  "CosMx_state_dual_high_with_coordinates.csv"
)

STEP02_SENS <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "CosMx_spatial_validation",
  "CosMx_dual_high_OR_threshold_sensitivity.csv"
)

STEP03 <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "CosMx_metadata_QC",
  "CosMx_state_dual_high_with_coordinates_AND_metadata_QC.csv"
)

S22_SENS <- file.path(
  PROJECT_DIR,
  "results",
  "canonical",
  "S21_S22_CosMx_v1.0",
  "figures",
  "Supplementary_Figure_S22_threshold_sensitivity.csv"
)

state_cols <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

############################################################
## Load and lineage checks
############################################################

s1 <- read_csv0(STEP01)
s2 <- read_csv0(STEP02)
s3 <- read_csv0(STEP03)

s1$key <- paste0(
  s1$fov,
  "__",
  s1$cell_ID
)

s2$key <- paste0(
  s2$fov,
  "__",
  s2$cell_ID
)

s3$key <- paste0(
  s3$fov,
  "__",
  s3$cell_ID
)

lineage_diff <- max(
  vapply(
    state_cols,
    function(x)
      maxdiff(
        s1[
          match(s2$key,s1$key),
          x
        ],
        s2[[x]]
      ),
    numeric(1)
  )
)

step03_diff <- max(
  vapply(
    c(
      state_cols,
      "spatial_x",
      "spatial_y"
    ),
    function(x)
      maxdiff(
        s2[[x]],
        s3[
          match(s2$key,s3$key),
          x
        ]
      ),
    numeric(1)
  )
)

############################################################
## S22 threshold sensitivity
############################################################

s22_diff <- NA_real_
s22_match_n <- 0L

s22_detail <- data.frame(
  key = character(),
  Mode = character(),
  TopFraction = numeric(),
  Step02_OR = numeric(),
  S22_OR = numeric(),
  OR_abs_diff = numeric(),
  Step02_low = numeric(),
  S22_low = numeric(),
  low_abs_diff = numeric(),
  Step02_high = numeric(),
  S22_high = numeric(),
  high_abs_diff = numeric(),
  stringsAsFactors = FALSE
)

if (
  file.exists(STEP02_SENS) &&
    file.exists(S22_SENS)
) {

  a <- read_csv0(
    STEP02_SENS
  )

  b <- read_csv0(
    S22_SENS
  )

  required_a <- c(
    "cutoff_mode",
    "top_fraction",
    "odds_ratio",
    "conf_low",
    "conf_high"
  )

  required_b <- c(
    "Mode",
    "TopFraction",
    "OR",
    "low",
    "high"
  )

  if (
    !all(required_a %in% colnames(a))
  ) {
    stop(
      "Step02 sensitivity table lacks required columns: ",
      paste(
        setdiff(
          required_a,
          colnames(a)
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (
    !all(required_b %in% colnames(b))
  ) {
    stop(
      "S22 sensitivity table lacks required columns: ",
      paste(
        setdiff(
          required_b,
          colnames(b)
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  a2 <- data.frame(
    Mode = normalize_mode(
      a$cutoff_mode
    ),
    TopFraction = normalize_fraction(
      a$top_fraction
    ),
    OR = suppressWarnings(
      as.numeric(
        a$odds_ratio
      )
    ),
    low = suppressWarnings(
      as.numeric(
        a$conf_low
      )
    ),
    high = suppressWarnings(
      as.numeric(
        a$conf_high
      )
    ),
    stringsAsFactors = FALSE
  )

  b2 <- data.frame(
    Mode = normalize_mode(
      b$Mode
    ),
    TopFraction = normalize_fraction(
      b$TopFraction
    ),
    OR = suppressWarnings(
      as.numeric(
        b$OR
      )
    ),
    low = suppressWarnings(
      as.numeric(
        b$low
      )
    ),
    high = suppressWarnings(
      as.numeric(
        b$high
      )
    ),
    stringsAsFactors = FALSE
  )

  a2$key <- make_key(
    a2$Mode,
    a2$TopFraction
  )

  b2$key <- make_key(
    b2$Mode,
    b2$TopFraction
  )

  if (
    anyDuplicated(a2$key) ||
      anyDuplicated(b2$key)
  ) {
    stop(
      "Duplicated Mode + TopFraction key in S22 sensitivity comparison.",
      call. = FALSE
    )
  }

  expected_keys <- as.vector(
    outer(
      c(
        "pooled",
        "per_fov"
      ),
      c(
        0.20,
        0.25,
        0.30
      ),
      FUN = function(m, f) {
        paste0(
          m,
          "__",
          sprintf(
            "%.6f",
            f
          )
        )
      }
    )
  )

  if (
    !setequal(
      a2$key,
      expected_keys
    )
  ) {
    stop(
      "Step02 sensitivity table does not contain the expected 2 x 3 key grid.",
      call. = FALSE
    )
  }

  if (
    !setequal(
      b2$key,
      expected_keys
    )
  ) {
    stop(
      "S22 sensitivity table does not contain the expected 2 x 3 key grid.",
      call. = FALSE
    )
  }

  b2 <- b2[
    match(
      a2$key,
      b2$key
    ),
    ,
    drop = FALSE
  ]

  s22_match_n <- nrow(
    a2
  )

  s22_detail <- data.frame(
    key = a2$key,
    Mode = a2$Mode,
    TopFraction = a2$TopFraction,
    Step02_OR = a2$OR,
    S22_OR = b2$OR,
    OR_abs_diff = abs(
      a2$OR -
        b2$OR
    ),
    Step02_low = a2$low,
    S22_low = b2$low,
    low_abs_diff = abs(
      a2$low -
        b2$low
    ),
    Step02_high = a2$high,
    S22_high = b2$high,
    high_abs_diff = abs(
      a2$high -
        b2$high
    ),
    stringsAsFactors = FALSE
  )

  all_diffs <- c(
    s22_detail$OR_abs_diff,
    s22_detail$low_abs_diff,
    s22_detail$high_abs_diff
  )

  if (
    length(all_diffs) > 0L &&
      all(
        is.finite(
          all_diffs
        )
      )
  ) {
    s22_diff <- max(
      all_diffs
    )
  }
}

write.csv(
  s22_detail,
  file.path(
    AUDIT_DIR,
    "11_CosMx_v1.3_S22_threshold_comparison_detail.csv"
  ),
  row.names = FALSE
)

############################################################
## Summary
############################################################

summary <- data.frame(
  Gate=c(
    "Step01 raw universe",
    "Step01->Step02 score lineage",
    "Step02->Step03 lineage",
    "S22 threshold sensitivity"
  ),
  Pass=c(
    nrow(s1)==EXPECTED_RAW_CELLS,
    lineage_diff<=TOL_EXACT,
    step03_diff<=TOL_EXACT,
    is.finite(s22_diff) &&
      s22_diff<=TOL_REPORT
  ),
  Observed=c(
    nrow(s1),
    format(lineage_diff,scientific=TRUE),
    format(step03_diff,scientific=TRUE),
    paste0(
      "matched=",
      s22_match_n,
      "/6; max_diff=",
      format(s22_diff,scientific=TRUE)
    )
  )
)

write.csv(
  summary,
  file.path(
    AUDIT_DIR,
    "11_CosMx_v1.3_final_gate_detail.csv"
  ),
  row.names=FALSE
)

blocking <- sum(
  !summary$Pass
)

cat("\n============================================================\n")
cat("11 COSMX REPRODUCTION CHECK v1.3\n")
cat("============================================================\n")

cat(
  "Step01 -> Step02 max diff: ",
  format(lineage_diff,scientific=TRUE),
  "\n",
  sep=""
)

cat(
  "Step02 -> Step03 max diff: ",
  format(step03_diff,scientific=TRUE),
  "\n",
  sep=""
)

cat(
  "S22 threshold matched rows: ",
  s22_match_n,
  "/6\n",
  sep=""
)

cat(
  "S22 threshold max diff: ",
  format(s22_diff,scientific=TRUE),
  "\n",
  sep=""
)

cat(
  "Blocking failures: ",
  blocking,
  "\n",
  sep=""
)

cat(
  "Overall: ",
  ifelse(blocking==0,"PASS","FAIL"),
  "\n",
  sep=""
)

if(blocking>0){
  stop(
    "11 CosMx reproduction audit v1.3 failed.",
    call.=FALSE
  )
}
