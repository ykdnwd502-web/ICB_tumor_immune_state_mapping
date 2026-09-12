############################################################
## 00_RUN_IF_mIF_SEQUENTIAL.R
##
## Complete public IF/mIF execution chain.
##
## 01  raw CosMx expression + gene-set discovery -> 4-state scores
## 02  IF metadata overlap selection -> matched-cell table -> Spearman/FDR
## 03  Supplementary Figure S23
## 04  Supplementary Table S23
## 05  compare all key rebuilt outputs with frozen Step14B v4 / Step14A v8
##
## Reproducibility:
##   - isolated child environments
##   - console log
##   - execution summary CSV
############################################################


options(
  stringsAsFactors = FALSE
)


############################################################
## Project paths
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


Sys.setenv(
  ICB_PROJECT_DIR = PROJECT_DIR
)


SCRIPT_DIR <- file.path(
  PROJECT_DIR,
  "scripts",
  "analysis",
  "13_IF_mIF"
)



############################################################
## Script list
############################################################

scripts <- c(
  "01_IF_mIF_build_CosMx_state_scores.R",
  "02_IF_mIF_match_and_correlation.R",
  "03_IF_mIF_make_S23.R",
  "04_IF_mIF_make_Supplementary_Table_S23.R",
  "05_IF_mIF_compare_with_frozen.R"
)


script_paths <- file.path(
  SCRIPT_DIR,
  scripts
)



############################################################
## Check scripts
############################################################

missing_scripts <- script_paths[
  !file.exists(script_paths)
]


if (length(missing_scripts) > 0L) {
  
  stop(
    "Missing script(s):\n",
    paste(
      missing_scripts,
      collapse = "\n"
    ),
    call. = FALSE
  )
  
}



############################################################
## Logging setup
############################################################

log_dir <- file.path(
  PROJECT_DIR,
  "results",
  "audit",
  "13_IF_mIF_runner"
)


dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


console_log <- file.path(
  log_dir,
  "13_IF_mIF_sequential_runner_console.log"
)


execution_csv <- file.path(
  log_dir,
  "13_IF_mIF_sequential_runner_execution_log.csv"
)



## Create single logging connection

log_con <- file(
  console_log,
  open = "wt"
)


sink(
  log_con,
  split = TRUE
)


sink(
  log_con,
  type = "message"
)



############################################################
## Execution table
############################################################

execution_log <- data.frame(
  Script = scripts,
  Status = NA_character_,
  Start = NA_character_,
  End = NA_character_,
  Elapsed_minutes = NA_real_,
  Message = NA_character_,
  stringsAsFactors = FALSE
)



############################################################
## Start message
############################################################

cat(
  "\n============================================================\n",
  "13_IF_mIF COMPLETE SEQUENTIAL CHAIN START\n",
  "Project: ",
  PROJECT_DIR,
  "\n",
  "Script directory: ",
  SCRIPT_DIR,
  "\n",
  "============================================================\n",
  sep = ""
)



############################################################
## Run sequential scripts
############################################################


for (i in seq_along(script_paths)) {
  
  
  ff <- script_paths[[i]]
  
  nm <- scripts[[i]]
  
  
  
  cat(
    "\n============================================================\n",
    "RUNNING STEP ",
    i,
    "/",
    length(script_paths),
    "\n",
    nm,
    "\n",
    "============================================================\n",
    sep = ""
  )
  
  
  
  t0 <- Sys.time()
  
  
  execution_log$Start[[i]] <-
    format(
      t0,
      "%Y-%m-%d %H:%M:%S"
    )
  
  
  
  err <- NULL
  
  
  
  child_env <- new.env(
    parent = globalenv()
  )
  
  
  
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
  
  
  
  execution_log$End[[i]] <-
    format(
      t1,
      "%Y-%m-%d %H:%M:%S"
    )
  
  
  
  execution_log$Elapsed_minutes[[i]] <-
    round(
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
    
    
    execution_log$Status[[i]] <- "PASS"
    
    execution_log$Message[[i]] <- ""
    
    
    
    cat(
      "\n>>> STEP PASS: ",
      nm,
      " | elapsed ",
      execution_log$Elapsed_minutes[[i]],
      " min\n",
      sep = ""
    )
    
    
  } else {
    
    
    execution_log$Status[[i]] <- "FAIL"
    
    execution_log$Message[[i]] <- err
    
    
    
    write.csv(
      execution_log,
      execution_csv,
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
    
    
    
    ## close logging safely
    
    sink(
      type = "message"
    )
    
    sink()
    
    close(
      log_con
    )
    
    
    
    stop(
      "13_IF_mIF sequential reproduction stopped at: ",
      nm,
      call. = FALSE
    )
    
  }
  
}



############################################################
## Save final execution log
############################################################

write.csv(
  execution_log,
  execution_csv,
  row.names = FALSE
)



############################################################
## Final summary
############################################################


cat(
  "\n============================================================\n",
  "13_IF_mIF COMPLETE SEQUENTIAL CHAIN FINISHED\n",
  "Steps PASS: ",
  sum(
    execution_log$Status == "PASS"
  ),
  "/",
  nrow(
    execution_log
  ),
  "\n",
  "Execution log:\n",
  execution_csv,
  "\n",
  "Console log:\n",
  console_log,
  "\n",
  "============================================================\n",
  sep = ""
)



print(
  execution_log[
    ,
    c(
      "Script",
      "Status",
      "Elapsed_minutes"
    ),
    drop = FALSE
  ]
)



############################################################
## Close logging
############################################################


sink(
  type = "message"
)


sink()


close(
  log_con
)

############################################################
## Save sessionInfo
############################################################

session_file <- file.path(
  PROJECT_DIR,
  "logs",
  "sessionInfo_13_IF_mIF_full_sequential_reproduction.txt"
)


sink(session_file)

cat("============================================================\n")
cat("13_IF_mIF FULL SEQUENTIAL REPRODUCTION SESSION INFO\n")
cat("============================================================\n\n")

cat("Project:\n")
cat(PROJECT_DIR)

cat("\n\nModule:\n")
cat("13_IF_mIF\n")

cat("\n\nExecution date:\n")
cat(as.character(Sys.Date()))

cat("\n\nR session information:\n\n")

print(sessionInfo())

sink()


cat(
  "\nSessionInfo saved:\n",
  session_file,
  "\n"
)