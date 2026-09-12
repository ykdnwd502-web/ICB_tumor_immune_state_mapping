############################################################
## run_11_CosMx_full_sequential_reproduction.R
##
## One-click public sequential rerun:
##   01 -> 02 -> 03 -> 04 -> 05 -> reproduction check
##
## Each child script is sourced into a separate environment.
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

Sys.setenv(
  ICB_PROJECT_DIR = PROJECT_DIR
)

SCRIPT_DIR <- file.path(
  PROJECT_DIR,
  "scripts",
  "analysis",
  "11_CosMx"
)

steps <- c(
  "01_CosMx_state_scoring.R",
  "02_CosMx_cell_level_spatial_validation.R",
  "03_CosMx_metadata_QC_and_technical_sensitivity.R",
  "04_CosMx_make_S21.R",
  "05_CosMx_make_S22.R",
  "11_CosMx_reproduction_check.R"
)

step_paths <- file.path(
  SCRIPT_DIR,
  steps
)

missing_scripts <- step_paths[
  !file.exists(step_paths)
]

if (length(missing_scripts) > 0L) {
  stop(
    "Missing 11_CosMx script(s):\n",
    paste(
      missing_scripts,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

log_dir <- file.path(
  PROJECT_DIR,
  "results",
  "audit",
  "11_CosMx_full_sequential_runner"
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

setwd(
  PROJECT_DIR
)

############################################################
# SessionInfo logging
############################################################

session_log <- file.path(
  PROJECT_DIR,
  "logs",
  "sessionInfo_11_CosMx_full_sequential_reproduction.txt"
)

dir.create(
  dirname(session_log),
  recursive = TRUE,
  showWarnings = FALSE
)

capture.output(
  sessionInfo(),
  file = session_log
)

cat(
  "\n============================================================\n",
  "11_CosMx FULL SEQUENTIAL REPRODUCTION START\n",
  "Project: ", PROJECT_DIR, "\n",
  "Script dir: ", SCRIPT_DIR, "\n",
  "============================================================\n",
  sep = ""
)

for (i in seq_along(step_paths)) {

  ff <- step_paths[[i]]
  nm <- steps[[i]]

  cat(
    "\n\n============================================================\n",
    "RUNNING STEP ", i, "/", length(step_paths), "\n",
    nm, "\n",
    "============================================================\n",
    sep = ""
  )

  t0 <- Sys.time()

  run_log$Start[[i]] <- format(
    t0,
    "%Y-%m-%d %H:%M:%S"
  )

  child_env <- new.env(
    parent = globalenv()
  )

  err <- NULL

  tryCatch(
    {
      sys.source(
        ff,
        envir = child_env,
        chdir = FALSE,
        keep.source = TRUE
      )
    },
    error = function(e) {
      err <<- conditionMessage(e)
    }
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

    cat(
      "\n>>> STEP PASS: ",
      nm,
      " | elapsed ",
      run_log$Elapsed_minutes[[i]],
      " min\n",
      sep = ""
    )

  } else {

    run_log$Status[[i]] <- "FAIL"
    run_log$Message[[i]] <- err

    write.csv(
      run_log,
      file.path(
        log_dir,
        "11_CosMx_full_sequential_runner_log.csv"
      ),
      row.names = FALSE
    )

    cat(
      "\n>>> STEP FAIL: ",
      nm,
      "\n",
      err,
      "\n",
      sep = ""
    )

    stop(
      "11_CosMx sequential reproduction stopped at: ",
      nm,
      call. = FALSE
    )
  }
}

write.csv(
  run_log,
  file.path(
    log_dir,
    "11_CosMx_full_sequential_runner_log.csv"
  ),
  row.names = FALSE
)

cat(
  "\n============================================================\n",
  "11_CosMx FULL SEQUENTIAL REPRODUCTION COMPLETED\n",
  "Steps PASS: ",
  sum(
    run_log$Status ==
      "PASS"
  ),
  "/",
  nrow(
    run_log
  ),
  "\n",
  "Runner log: ",
  file.path(
    log_dir,
    "11_CosMx_full_sequential_runner_log.csv"
  ),
  "\n",
  "============================================================\n",
  sep = ""
)

print(
  run_log[
    ,
    c(
      "Step",
      "Status",
      "Elapsed_minutes"
    ),
    drop = FALSE
  ]
)

capture.output(
  sessionInfo(),
  file = session_log
)
