############################################################
## 00_run_data_acquisition.R
##
## One-click raw-data preparation runner.
##
## IMPORTANT:
## CosMx Dryad Slide_4.zip may require manual download unless
## COSMX_SLIDE4_URL is supplied. If step 03 stops for that reason:
##   1) download Slide_4.zip from the Dryad DOI;
##   2) place it in the exact path reported;
##   3) rerun this runner.
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

SCRIPT_DIR <- file.path(
  PROJECT_DIR,
  "scripts",
  "data_acquisition"
)

steps <- c(
  "01_download_GEO_inputs.R",
  "02_prepare_primary_Visium_10x.R",
  "03_prepare_CosMx_Dryad_Slide4.R",
  "04_prepare_external_spatial_inputs.R",
  "05_validate_raw_data_contract.R"
)

paths <- file.path(
  SCRIPT_DIR,
  steps
)

missing <- paths[
  !file.exists(paths)
]

if (length(missing) > 0L) {
  stop(
    "Missing acquisition script(s):\n",
    paste(missing, collapse = "\n"),
    call. = FALSE
  )
}

log_dir <- file.path(
  PROJECT_DIR,
  "release_metadata",
  "raw_data"
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

run_log <- data.frame(
  Step = steps,
  Status = NA_character_,
  Start = NA_character_,
  End = NA_character_,
  Elapsed_minutes = NA_real_,
  Message = NA_character_,
  stringsAsFactors = FALSE
)

old_wd <- getwd()

on.exit(
  setwd(old_wd),
  add = TRUE
)

setwd(PROJECT_DIR)

for (i in seq_along(paths)) {

  cat(
    "\n============================================================\n",
    "RAW-DATA STEP ", i, "/", length(paths), "\n",
    steps[[i]], "\n",
    "============================================================\n",
    sep = ""
  )

  t0 <- Sys.time()
  run_log$Start[[i]] <- format(
    t0,
    "%Y-%m-%d %H:%M:%S"
  )

  env <- new.env(
    parent = globalenv()
  )

  err <- tryCatch(
    {
      sys.source(
        paths[[i]],
        envir = env,
        chdir = FALSE
      )
      NULL
    },
    error = function(e) e
  )

  t1 <- Sys.time()

  run_log$End[[i]] <- format(
    t1,
    "%Y-%m-%d %H:%M:%S"
  )

  run_log$Elapsed_minutes[[i]] <- round(
    as.numeric(
      difftime(
        t1,
        t0,
        units = "mins"
      )
    ),
    3
  )

  if (is.null(err)) {
    run_log$Status[[i]] <- "PASS"
    run_log$Message[[i]] <- ""
  } else {
    run_log$Status[[i]] <- "FAIL"
    run_log$Message[[i]] <- conditionMessage(err)

    write.csv(
      run_log,
      file.path(
        log_dir,
        "RAW_DATA_ACQUISITION_RUN_LOG.csv"
      ),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    stop(
      "Data-acquisition step failed: ",
      steps[[i]],
      "\n",
      conditionMessage(err),
      call. = FALSE
    )
  }

  write.csv(
    run_log,
    file.path(
      log_dir,
      "RAW_DATA_ACQUISITION_RUN_LOG.csv"
    ),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(
    log_dir,
    "sessionInfo_data_acquisition.txt"
  ),
  useBytes = TRUE
)

cat("\n============================================================\n")
cat("RAW DATA ACQUISITION / PREPARATION COMPLETE\n")
cat("============================================================\n")
cat("All steps PASS: ", all(run_log$Status == "PASS"), "\n", sep = "")
cat("Next: run the frozen analysis runners from the project root.\n")
cat("============================================================\n")
