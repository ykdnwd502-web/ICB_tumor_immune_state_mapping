############################################################
## run_12_Primary_Visium_full_sequential_reproduction.R
##
## One-click public sequential rerun:
##
## 01 -> 02 -> 03 -> 03B -> 04 -> 05
## -> 06 -> 07 -> 08 -> 09 -> 09A
##
## Each child script is sourced into a separate environment.
##
############################################################


options(stringsAsFactors = FALSE)


############################################################
# Project configuration
############################################################


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
  "12_Primary_Visium"
)



############################################################
# Sequential execution order
############################################################


steps <- c(
  
  "01_Visium_build_object_and_state_scoring.R",
  
  "02_Visium_spatial_QC_Moran_colocalization.R",
  
  "03_Visium_dual_high_niche_analysis.R",
  
  "03B_Visium_make_Figure8.R",
  
  "04_Visium_threshold_independent_covariation.R",
  
  "05_Visium_composition_aware_analysis.R",
  
  "06_Visium_humoral_TLS_context.R",
  
  "07_Visium_primary_canonical_analysis.R",
  
  "08_Visium_make_S15_S19_figures.R",
  
  "09_Visium_compare_with_frozen.R",
  
  "09A_Visium_diagnose_reproduction_mismatches.R"
  
)



step_paths <- file.path(
  SCRIPT_DIR,
  steps
)



############################################################
# Script existence check
############################################################


missing_scripts <- step_paths[
  !file.exists(step_paths)
]


if (length(missing_scripts) > 0L) {
  
  stop(
    "Missing 12_Primary_Visium script(s):\n",
    paste(
      missing_scripts,
      collapse = "\n"
    ),
    call. = FALSE
  )
  
}



############################################################
# Audit log
############################################################


log_dir <- file.path(
  PROJECT_DIR,
  "results",
  "audit",
  "12_Primary_Visium_full_sequential_runner"
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



############################################################
# Session information
############################################################


session_log <- file.path(
  PROJECT_DIR,
  "logs",
  "sessionInfo_12_Primary_Visium_full_sequential_reproduction.txt"
)


dir.create(
  dirname(session_log),
  recursive = TRUE,
  showWarnings = FALSE
)



############################################################
# Working directory
############################################################


old_wd <- getwd()

on.exit(
  setwd(old_wd),
  add = TRUE
)


setwd(
  PROJECT_DIR
)



############################################################
# Start message
############################################################


cat(
  "\n============================================================\n",
  "12_Primary_Visium FULL SEQUENTIAL REPRODUCTION START\n",
  "Project: ",
  PROJECT_DIR,
  "\n",
  "Script dir: ",
  SCRIPT_DIR,
  "\n",
  "============================================================\n",
  sep = ""
)



############################################################
# Execute pipeline
############################################################


for (i in seq_along(step_paths)) {
  
  
  ff <- step_paths[[i]]
  
  nm <- steps[[i]]
  
  
  
  cat(
    "\n\n============================================================\n",
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
        "12_Primary_Visium_full_sequential_runner_log.csv"
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
    
    
    
    capture.output(
      sessionInfo(),
      file = session_log
    )
    
    
    
    stop(
      "12_Primary_Visium sequential reproduction stopped at: ",
      nm,
      call. = FALSE
    )
    
  }
  
}



############################################################
# Save final logs
############################################################


write.csv(
  run_log,
  file.path(
    log_dir,
    "12_Primary_Visium_full_sequential_runner_log.csv"
  ),
  row.names = FALSE
)



capture.output(
  sessionInfo(),
  file = session_log
)



############################################################
# Final summary
############################################################


cat(
  "\n============================================================\n",
  "12_Primary_Visium FULL SEQUENTIAL REPRODUCTION COMPLETED\n",
  "Steps PASS: ",
  sum(
    run_log$Status == "PASS"
  ),
  "/",
  nrow(
    run_log
  ),
  "\n",
  "Runner log: ",
  file.path(
    log_dir,
    "12_Primary_Visium_full_sequential_runner_log.csv"
  ),
  "\n",
  "SessionInfo: ",
  session_log,
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