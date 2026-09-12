############################################################
## 00_RUN_CellChat_LR_SEQUENTIAL.R
##
## Complete public sequential execution chain:
##   CellChat nomination -> CosMx k=20/10/30 proximity ->
##   reporting sources -> S24/S25/S26 ->
##   Tables S24-S26 ->
##   current-vs-frozen reproduction audit.
############################################################


############################################################
## Project paths
############################################################

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR")

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
  "14_CellChat_LR"
)



############################################################
## Script runner
############################################################

run_script <- function(name) {
  
  p <- file.path(
    SCRIPT_DIR,
    name
  )
  
  if (!file.exists(p)) {
    stop(
      "Missing CellChat/LR script: ",
      p,
      call. = FALSE
    )
  }
  
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "RUNNING: ",
    p,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================================\n"
  )
  
  
  sys.source(
    p,
    envir = .GlobalEnv
  )
  
}



############################################################
## Sequential execution
############################################################


run_script(
  "01_CellChat_candidate_LR_nomination.R"
)

run_script(
  "02_CosMx_spatial_LR_primary_K20.R"
)

run_script(
  "03_CosMx_spatial_LR_sensitivity_K10.R"
)

run_script(
  "04_CosMx_spatial_LR_sensitivity_K30.R"
)

run_script(
  "05_CellChat_LR_build_reporting_sources.R"
)

run_script(
  "06_CellChat_LR_make_S24.R"
)

run_script(
  "07_CellChat_LR_make_S25.R"
)

run_script(
  "08_CellChat_LR_make_S26.R"
)

run_script(
  "09_CellChat_LR_make_Supplementary_Tables_S24_S26.R"
)

run_script(
  "10_CellChat_LR_compare_with_frozen.R"
)



############################################################
## Completion message
############################################################


cat(
  "\n============================================================\n"
)

cat(
  "14_CellChat_LR COMPLETE SEQUENTIAL CHAIN FINISHED\n"
)

cat(
  "============================================================\n"
)



############################################################
## Save sessionInfo
############################################################


audit_dir <- file.path(
  PROJECT_DIR,
  "logs",
  "14_CellChat_LR_runner"
)


dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


session_file <- file.path(
  audit_dir,
  "sessionInfo_14_CellChat_LR_full_sequential_reproduction.txt"
)


sink(
  session_file
)


cat(
  "============================================================\n"
)

cat(
  "14_CellChat_LR FULL SEQUENTIAL REPRODUCTION SESSION INFO\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Project:\n"
)

cat(
  PROJECT_DIR
)


cat(
  "\n\nModule:\n"
)

cat(
  "14_CellChat_LR\n"
)


cat(
  "\nExecution date:\n"
)

cat(
  as.character(Sys.Date())
)


cat(
  "\n\nR session information:\n\n"
)


print(
  sessionInfo()
)


sink()



cat(
  "\nSessionInfo saved:\n",
  session_file,
  "\n"
)



############################################################
## Final
############################################################


cat(
  "\n============================================================\n"
)

cat(
  "14_CellChat_LR SESSION INFO COMPLETED\n"
)

cat(
  "============================================================\n"
)
