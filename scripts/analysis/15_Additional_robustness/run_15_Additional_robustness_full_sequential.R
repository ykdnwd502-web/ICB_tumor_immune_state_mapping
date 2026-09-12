############################################################
## run_15_Additional_robustness_full_sequential.R
############################################################
options(stringsAsFactors = FALSE)

# 1. 设置工作路径
PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR", unset = "D:/ICB_resistance_project")
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
Sys.setenv(ICB_PROJECT_DIR = PROJECT_DIR)

# 2. 确保 logs 文件夹存在
log_dir <- file.path(PROJECT_DIR, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# 3. 定义执行步骤 (已将所有的 "figures" 路径修正为 "analysis")
steps <- c(
  file.path(PROJECT_DIR,"scripts","analysis","15_Additional_robustness","01_GSE244983_patient_level_summary.R"),
  file.path(PROJECT_DIR,"scripts","analysis","15_Additional_robustness","01_make_S29_patient_level_structure.R"),
  file.path(PROJECT_DIR,"scripts","analysis","15_Additional_robustness","02_ICB_cross_cohort_comparability.R"),
  file.path(PROJECT_DIR,"scripts","analysis","15_Additional_robustness","02_make_S30_cross_cohort_robustness.R"),
  file.path(PROJECT_DIR,"scripts","analysis","15_Additional_robustness","03_statistical_testing_inventory.R"),
  file.path(PROJECT_DIR,"scripts","analysis","15_Additional_robustness","03_make_Table_S29_statistical_analysis_summary.R"),
  file.path(PROJECT_DIR,"scripts","analysis","15_Additional_robustness","04_make_Table_S30_cross_cohort_comparability.R"),
  file.path(PROJECT_DIR,"scripts","analysis","15_Additional_robustness","15_Additional_robustness_reproduction_check.R")
)

missing <- steps[!file.exists(steps)]
if(length(missing) > 0) {
  stop("Missing script(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
}

cat("\n============================================================\n")
cat("15 ADDITIONAL ROBUSTNESS SEQUENTIAL RUN START\n")
cat("============================================================\n")

# 4. 循环执行脚本
for(i in seq_along(steps)){
  cat("\n--- STEP ", i, "/", length(steps), ": ", basename(steps[[i]]), " ---\n", sep = "")
  source(steps[[i]], local = new.env(parent = globalenv()), echo = FALSE, chdir = FALSE)
}

# 5. 统一保存 Runner 级别的 sessionInfo 到 logs 文件夹
session_file <- file.path(log_dir, "sessionInfo_15_Runner_Sequential.txt")
sink(session_file)
cat("=== RUNNER SESSION INFO ===\n")
cat("Time:", as.character(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("15 ADDITIONAL ROBUSTNESS SEQUENTIAL RUN COMPLETE\n")
cat("Runner sessionInfo saved to: logs/sessionInfo_15_Runner_Sequential.txt\n")
cat("============================================================\n")