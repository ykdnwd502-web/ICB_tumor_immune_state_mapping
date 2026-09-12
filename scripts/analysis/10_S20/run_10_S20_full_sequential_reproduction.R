############################################################
##
## run_10_S20_full_sequential_reproduction.R
##
## One-click sequential runner for the frozen S20 module.
##
## Sequence:
##
## 01 -> 02 -> 03 -> 04 -> 05
## -> 06 lock verification
## -> 07 aggregation
## -> 08 figure
## -> 09 final reproduction check
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


Sys.setenv(
  ICB_PROJECT_DIR = PROJECT_DIR
)


SCRIPT_DIR <- file.path(
  PROJECT_DIR,
  "scripts",
  "analysis",
  "10_S20"
)


steps <- c(
  
  "01_GSE250636_build_18J.R",
  
  "02_Thrane2018_build_18J.R",
  
  "03_build_historical_S28_from_18J.R",
  
  "04_GSE250636_build_S20_canonical.R",
  
  "05_Thrane2018_build_S20_canonical.R",
  
  "06_verify_S20_locked_state_genes.R",
  
  "07_S20_aggregate.R",
  
  "08_S20_make_figure.R",
  
  "09_S20_external_spatial_recurrence_reproduction_check.R"
  
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
    "Missing S20 script(s):\n",
    paste(
      missing_scripts,
      collapse = "\n"
    ),
    call. = FALSE
  )
  
}


old_wd <- getwd()

on.exit(
  setwd(old_wd),
  add = TRUE
)


setwd(PROJECT_DIR)

############################################################
# SessionInfo logging
############################################################

LOG_DIR <- file.path(
  PROJECT_DIR,
  "logs"
)

dir.create(
  LOG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


SESSION_LOG <- file.path(
  LOG_DIR,
  "sessionInfo_10_S20_full_sequential_reproduction.txt"
)


capture.output(
  sessionInfo(),
  file = SESSION_LOG
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



cat(
  "\n============================================================\n",
  "S20 FULL SEQUENTIAL REPRODUCTION START\n",
  "Project: ",
  PROJECT_DIR,
  "\n",
  "============================================================\n",
  sep = ""
)



for (i in seq_along(step_paths)) {
  
  
  ff <- step_paths[[i]]
  
  nm <- steps[[i]]
  
  
  cat(
    "\n============================================================\n",
    "RUNNING STEP ",
    i,
    "/",
    length(step_paths),
    "\n",
    nm,
    "\n",
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
    
    error = function(e){
      
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
  
  
  if (is.null(err)){
    
    
    run_log$Status[[i]] <- "PASS"
    
    run_log$Message[[i]] <- ""
    
    
    cat(
      "\n>>> STEP PASS: ",
      nm,
      "\n",
      sep=""
    )
    
    
  } else {
    
    
    run_log$Status[[i]] <- "FAIL"
    
    run_log$Message[[i]] <- err
    
    
    cat(
      "\n>>> STEP FAIL: ",
      nm,
      "\n",
      err,
      "\n",
      sep=""
    )
    
    
    stop(
      "S20 sequential reproduction stopped at: ",
      nm,
      call.=FALSE
    )
    
  }
  
}



log_dir <- file.path(
  PROJECT_DIR,
  "results",
  "audit",
  "10_S20_full_sequential_runner"
)


dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  run_log,
  file.path(
    log_dir,
    "10_S20_full_sequential_runner_log.csv"
  ),
  row.names = FALSE
)



cat(
  "\n============================================================\n",
  "S20 FULL SEQUENTIAL REPRODUCTION COMPLETED\n",
  "PASS: ",
  sum(run_log$Status=="PASS"),
  "/",
  nrow(run_log),
  "\n",
  "============================================================\n",
  sep=""
)


print(
  run_log[
    ,
    c(
      "Step",
      "Status",
      "Elapsed_minutes"
    ),
    drop=FALSE
  ]
)

############################################################
# Final sessionInfo
############################################################

capture.output(
  sessionInfo(),
  file = SESSION_LOG
)


cat(
  "\nSessionInfo saved:\n",
  SESSION_LOG,
  "\n"
)