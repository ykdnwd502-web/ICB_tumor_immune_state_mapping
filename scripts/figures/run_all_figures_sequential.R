############################################################
## run_all_figures_sequential.R
############################################################
options(stringsAsFactors = FALSE)

# 1. 设置工作路径
PROJECT_DIR <- "D:/ICB_resistance_project"
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = FALSE)
Sys.setenv(ICB_PROJECT_DIR = PROJECT_DIR)
setwd(PROJECT_DIR)

# 2. 确保 logs 文件夹存在
log_dir <- file.path(PROJECT_DIR, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 2.5 Load shared figure plotting helpers
##
## Required for clean reproducibility.
## No figure script should depend on previous R sessions.
############################################################

helper_file <- file.path(
  PROJECT_DIR,
  "scripts",
  "figures",
  "00_plot_helpers.R"
)


if(!file.exists(helper_file)){
  
  stop(
    "Missing shared plotting helper:\n",
    helper_file,
    call.=FALSE
  )
  
}


source(
  helper_file,
  local = FALSE
)


cat(
  "Loaded shared plotting helper:\n",
  helper_file,
  "\n"
)


# 3. 自定义执行顺序：先主 Figure (从小到大)，再 Supplementary Figure (从小到大)
steps <- c(
  # --- 1. Main Figures ---
  "scripts/figures/Figure_2.R",
  "scripts/figures/Figure_3.R",
  "scripts/figures/Figure_4.R",
  "scripts/figures/Figure_5.R",
  "scripts/figures/Figure_6.R",
  "scripts/figures/Figure 7.R",
  "scripts/figures/Figure 8.R",
  # Figure 9 结合了 S27 和 S28，作为最后一个主图脚本执行
  "scripts/figures/Figure_9_Supplementary Figure S27_Supplementary Figure S28.R",
  
  # --- 2. Supplementary Figures ---
  "scripts/figures/Supplementary Figure S2.R",
  "scripts/figures/Supplementary Figure S3.R",
  "scripts/figures/Supplementary Figure S4.R",
  "scripts/figures/Supplementary Figure S5.R",
  "scripts/figures/Supplementary Figure S6.R",
  "scripts/figures/Supplementary Figure S7.R",
  "scripts/figures/Supplementary Figure S8.R",
  "scripts/figures/Supplementary Figure S9.R",
  "scripts/figures/Supplementary Figure S10.R",
  "scripts/figures/Supplementary Figure S11.R",
  "scripts/figures/Supplementary Figure S12.R",
  "scripts/figures/Supplementary Figure S13.R",
  "scripts/figures/Supplementary Figure S14.R",
  "scripts/figures/Supplementary Figure S15.R",
  "scripts/figures/Supplementary Figure S16.R",
  "scripts/figures/Supplementary Figure S17.R",
  "scripts/figures/Supplementary Figure S18.R",
  "scripts/figures/Supplementary Figure S19.R",
  "scripts/figures/Supplementary Figure S20.R",
  "scripts/figures/Supplementary Figure S21.R",
  "scripts/figures/Supplementary Figure S22.R",
  "scripts/figures/Supplementary Figure S23.R",
  "scripts/figures/Supplementary Figure S24.R",
  "scripts/figures/Supplementary Figure S25.R",
  "scripts/figures/Supplementary Figure S26.R",
  "scripts/figures/Supplementary Figure S29.R",
  "scripts/figures/Supplementary Figure S30.R"
)

# 转换为绝对路径
steps_full <- file.path(PROJECT_DIR, steps)

# 检查文件是否存在（预检）
missing <- steps_full[!file.exists(steps_full)]
if(length(missing) > 0) {
  stop("Missing script(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
}

# 4. 核心执行逻辑：双重 tryCatch
tryCatch({
  cat("\n============================================================\n")
  cat("ALL FIGURES SEQUENTIAL RUN START\n")
  cat("============================================================\n")
  
  for(i in seq_along(steps_full)){
    script_path <- steps_full[[i]]
    script_name <- basename(script_path)
    
    cat("\n--- STEP ", i, "/", length(steps_full), ": ", script_name, " ---\n", sep = "")
    
    # 针对每一个单独的子脚本运行内置 tryCatch
    tryCatch({
      # 核心运行：使用 local = new.env() 隔离环境，防止变量相互污染
      source(script_path, local = new.env(parent = globalenv()), echo = FALSE, chdir = FALSE)
      cat("-> 成功执行: ", script_name, "\n")
      
    }, error = function(e) {
      cat("\n[!!!] 子脚本执行出错: ", script_name, "\n")
      cat("错误信息: ", e$message, "\n")
      # 如果您希望某个图报错后立刻停止后续所有脚本，请取消下方这行的注释：
      # stop(e) 
      
    }, finally = {
      # 无论该子脚本成功还是失败，都为它保存专属的 sessionInfo
      # 日志命名规则例如：sessionInfo_Figure_2.txt
      log_name <- paste0("sessionInfo_", tools::file_path_sans_ext(script_name), ".txt")
      log_path <- file.path(log_dir, log_name)
      
      sink(log_path)
      cat("=== SESSION INFO FOR:", script_name, "===\n")
      cat("Time:", as.character(Sys.time()), "\n\n")
      print(sessionInfo())
      sink()
      
      cat("-> 已保存本脚本独立环境日志至: logs/", log_name, "\n", sep = "")
    })
  }
  
  cat("\n============================================================\n")
  cat("ALL FIGURES SEQUENTIAL RUN COMPLETE\n")
  cat("============================================================\n")
  
}, error = function(e) {
  cat("\n[!!!] Runner 主进程捕获到致命错误:\n")
  cat(e$message, "\n\n")
  
}, finally = {
  # 5. 最后，为整个 Runner 脚本流水线保存一份全局的 sessionInfo
  global_session_file <- file.path(log_dir, "sessionInfo_Runner_Global.txt")
  
  sink(global_session_file)
  cat("=== GLOBAL RUNNER SESSION INFO ===\n")
  cat("Time:", as.character(Sys.time()), "\n\n")
  print(sessionInfo())
  sink()
  
  cat("\n-> 全局 Runner sessionInfo 已经强制保存至: logs/sessionInfo_Runner_Global.txt\n")
})