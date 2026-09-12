############################################################
##
## 06_verify_S20_locked_state_genes.R
##
## Verify ONLY the frozen S20 four-state gene lock.
##
## Frozen input:
##   public_release/
##   S20_external_spatial_recurrence/
##   config/
##   S20_locked_state_genes.csv
##
## Verification targets:
##   1) GSE250636 canonical state-gene table
##   2) Thrane2018 canonical state-gene table
##
## No gene selection.
## No gene filtering.
## No state reconstruction.
## No file writing to frozen configuration.
##
############################################################


options(stringsAsFactors = FALSE)


PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)


PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)



############################################################
# 1. Paths
############################################################


LOCK_FILE <- file.path(
  PROJECT_DIR,
  "public_release",
  "S20_external_spatial_recurrence",
  "config",
  "S20_locked_state_genes.csv"
)


GSE_CANONICAL_FILE <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "S20_external_spatial_recurrence_CANONICAL_v1.0",
  "GSE250636",
  "S20_GSE250636_state_genes_used_long.csv"
)


THRANE_CANONICAL_FILE <- file.path(
  PROJECT_DIR,
  "results",
  "tables",
  "S20_external_spatial_recurrence_CANONICAL_v1.0",
  "Thrane2018_legacyST",
  "S20_Thrane2018_state_genes_used_long.csv"
)



############################################################
# 2. Reader
############################################################


read_lock <- function(path, label) {
  
  
  if (!file.exists(path)) {
    
    stop(
      "Missing ",
      label,
      ": ",
      path,
      call. = FALSE
    )
    
  }
  
  
  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  if (!all(c("state", "gene") %in% colnames(x))) {
    
    stop(
      "Required columns state/gene missing in: ",
      path,
      call. = FALSE
    )
    
  }
  
  
  x <- unique(
    x[
      ,
      c(
        "state",
        "gene"
      ),
      drop = FALSE
    ]
  )
  
  
  x$state <- as.character(
    x$state
  )
  
  
  x$gene <- toupper(
    trimws(
      as.character(
        x$gene
      )
    )
  )
  
  
  x <- x[
    !is.na(x$state) &
      nzchar(x$state) &
      !is.na(x$gene) &
      nzchar(x$gene),
    ,
    drop = FALSE
  ]
  
  
  x <- x[
    order(
      x$state,
      x$gene
    ),
    ,
    drop = FALSE
  ]
  
  
  rownames(x) <- NULL
  
  
  x
  
}



############################################################
# 3. Load frozen configuration and canonical outputs
############################################################


cat(
  "\n============================================================\n",
  "S20 LOCK VERIFICATION START\n",
  "============================================================\n"
)



locked <- read_lock(
  LOCK_FILE,
  "frozen S20 locked state-gene configuration"
)


gse <- read_lock(
  GSE_CANONICAL_FILE,
  "GSE250636 canonical state-gene table"
)


thrane <- read_lock(
  THRANE_CANONICAL_FILE,
  "Thrane2018 canonical state-gene table"
)



############################################################
# 4. Canonical consistency checks
############################################################


if (!identical(gse, thrane)) {
  
  stop(
    paste0(
      "BLOCKED: GSE250636 and Thrane2018 canonical ",
      "state-gene tables are not identical."
    ),
    call. = FALSE
  )
  
}



if (!identical(locked, gse)) {
  
  stop(
    paste0(
      "BLOCKED: frozen S20_locked_state_genes.csv ",
      "does not match canonical rebuilt state-gene table."
    ),
    call. = FALSE
  )
  
}



############################################################
# 5. Four-state check
############################################################


required_states <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)


if (
  length(
    unique(locked$state)
  ) != 4L
) {
  
  stop(
    "BLOCKED: frozen S20 lock does not contain exactly four states.",
    call. = FALSE
  )
  
}



if (
  !all(
    required_states %in% unique(locked$state)
  )
) {
  
  stop(
    "BLOCKED: frozen S20 lock is missing one or more required states.",
    call. = FALSE
  )
  
}



############################################################
# 6. PASS summary
############################################################


cat(
  "\n============================================================\n",
  "S20 LOCK VERIFICATION PASS\n",
  "============================================================\n",
  "Frozen lock:\n",
  LOCK_FILE,
  "\n\n",
  "Canonical GSE250636:\n",
  GSE_CANONICAL_FILE,
  "\n\n",
  "Canonical Thrane2018:\n",
  THRANE_CANONICAL_FILE,
  "\n\n",
  "States: ",
  length(unique(locked$state)),
  "\n",
  "Rows: ",
  nrow(locked),
  "\n",
  "============================================================\n",
  sep = ""
)