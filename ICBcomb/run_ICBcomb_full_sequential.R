############################################################
## run_ICBcomb_full_sequential.R
## ICBcomb path-locked true clean-run runner v1.0.3
##
## This runner is designed for the ACTUAL module tree:
##   <ICBcomb root>/run_ICBcomb_full_sequential.R
##   <ICBcomb root>/scripts/audit/01_ICBcomb_source_lock_audit.R
##   <ICBcomb root>/scripts/analysis/02_ICBcomb_build_source_tables.R
##   <ICBcomb root>/scripts/figures/03_make_Figure9_ICBcomb.R
##   <ICBcomb root>/scripts/figures/04_make_Supplementary_Figure_S27.R
##   <ICBcomb root>/scripts/figures/05_make_Supplementary_Figure_S28.R
##   <ICBcomb root>/scripts/tables/06_make_Supplementary_Tables_S27A_S27B.R
##   <ICBcomb root>/scripts/audit/07_ICBcomb_post_rebuild_audit.R
##
## IMPORTANT PATH FIX
##   The module root is resolved from THIS runner file when source() is used.
##   It does NOT assume getwd() is the ICBcomb folder and does NOT expect an
##   extra scripts/*/ICBcomb/ directory level.
############################################################

options(stringsAsFactors = FALSE)

.normalize_existing_dir <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(x)) return(NA_character_)
  if (!dir.exists(x)) return(NA_character_)
  normalizePath(x, winslash = "/", mustWork = TRUE)
}

.runner_file <- tryCatch({
  of <- sys.frame(1)$ofile
  if (is.null(of) || !nzchar(of)) NA_character_ else normalizePath(of, winslash = "/", mustWork = TRUE)
}, error = function(e) NA_character_)

## Resolution order:
## 1) explicit ICBCOMB_ROOT
## 2) directory containing this sourced runner
## 3) explicit ICB_PUBLIC_ROOT (legacy compatibility)
## 4) getwd(), but only if it itself is a valid ICBcomb module root
root_candidates <- c(
  Sys.getenv("ICBCOMB_ROOT"),
  if (!is.na(.runner_file)) dirname(.runner_file) else NA_character_,
  Sys.getenv("ICB_PUBLIC_ROOT"),
  getwd()
)
root_candidates <- unique(root_candidates[!is.na(root_candidates) & nzchar(root_candidates)])

.is_module_root <- function(x) {
  dir.exists(x) &&
    dir.exists(file.path(x, "inputs", "ICBcomb_INPUT_FREEZE_v1.0")) &&
    dir.exists(file.path(x, "scripts"))
}

valid_roots <- root_candidates[vapply(root_candidates, .is_module_root, logical(1))]
if (length(valid_roots) == 0L) {
  stop(
    "Could not resolve the ICBcomb module root.\n",
    "The runner must reside in the folder containing inputs/ and scripts/.\n",
    "Runner file seen as: ", ifelse(is.na(.runner_file), "<unavailable>", .runner_file), "\n",
    "getwd(): ", getwd(), "\n",
    "You may also set ICBCOMB_ROOT explicitly.",
    call. = FALSE
  )
}
ROOT <- normalizePath(valid_roots[[1]], winslash = "/", mustWork = TRUE)

## Force every downstream module script to use the same resolved root.
Sys.setenv(ICBCOMB_ROOT = ROOT)
Sys.setenv(ICB_PUBLIC_ROOT = ROOT)

steps <- c(
  "scripts/audit/01_ICBcomb_source_lock_audit.R",
  "scripts/analysis/02_ICBcomb_build_source_tables.R",
  "scripts/figures/03_make_Figure9_ICBcomb.R",
  "scripts/figures/04_make_Supplementary_Figure_S27.R",
  "scripts/figures/05_make_Supplementary_Figure_S28.R",
  "scripts/tables/06_make_Supplementary_Tables_S27A_S27B.R",
  "scripts/audit/07_ICBcomb_post_rebuild_audit.R"
)
step_paths <- file.path(ROOT, steps)

INPUT <- file.path(ROOT, "inputs", "ICBcomb_INPUT_FREEZE_v1.0")
preflight_inputs <- c(
  file.path(INPUT, "SHA256SUMS.csv"),
  file.path(INPUT, "data", "raw_platform_results", "ICBcomb_Immune_cold_RESTORE_hybrid_results.csv"),
  file.path(INPUT, "data", "raw_platform_results", "ICBcomb_Myeloid_Treg_SUPPRESS_hybrid_results.csv"),
  file.path(INPUT, "data", "raw_platform_results", "ICBcomb_Tumor_dedifferentiation_SUPPRESS_hybrid_results.csv"),
  file.path(INPUT, "data", "query_gene_set_bundle", "02_final_ICBcomb_gene_sets", "final_ICBcomb_gene_set_summary_CLEAN.csv"),
  file.path(INPUT, "data", "query_gene_set_bundle", "02_final_ICBcomb_gene_sets", "Immune_defective_Cold_RESTORE_hybrid.txt"),
  file.path(INPUT, "data", "query_gene_set_bundle", "02_final_ICBcomb_gene_sets", "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid.txt"),
  file.path(INPUT, "data", "query_gene_set_bundle", "02_final_ICBcomb_gene_sets", "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid.txt"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_combined_raw_results_CLEAN.csv"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_positive_reversal_rows_CLEAN.csv"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_query_gene_set_counts_CLEAN.csv"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_query_gene_set_counts_FINAL_LOCKED.csv"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_state_strategy_summary_CLEAN.csv"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_strategy_mapping_audit_CLEAN.csv"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_top_positive_reversal_strategies_FOR_FIGURE8B_FINAL_LOCKED.csv"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_strategy_state_heatmap_table_FOR_FIGURE8C_FINAL_LOCKED.csv"),
  file.path(INPUT, "data", "historical_reference", "ICBcomb_reversal_analysis", "ICBcomb_state_level_maximum_positive_reversal_signal_FOR_FIGURE8D_FINAL_LOCKED.csv")
)

missing_scripts <- step_paths[!file.exists(step_paths)]
missing_inputs <- preflight_inputs[!file.exists(preflight_inputs)]

cat("============================================================\n")
cat("ICBcomb CLEAN RUN v1.0.3 PATH PREFLIGHT\n")
cat("Resolved module root: ", ROOT, "\n", sep = "")
cat("Runner file: ", ifelse(is.na(.runner_file), "<unavailable>", .runner_file), "\n", sep = "")
cat("Working directory: ", normalizePath(getwd(), winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Expected script layout: scripts/audit|analysis|figures|tables directly (NO extra /ICBcomb/ level)\n")
cat("Named-source path mapping hardening: ENABLED for Figure 9 and S27A/S27B rebuilds\n")
cat("============================================================\n")

if (length(missing_scripts) > 0L || length(missing_inputs) > 0L) {
  if (length(missing_scripts) > 0L) {
    cat("\nMISSING REQUIRED SCRIPT(S):\n")
    cat(paste0("  - ", missing_scripts), sep = "\n")
    cat("\n")
  }
  if (length(missing_inputs) > 0L) {
    cat("\nMISSING REQUIRED FROZEN INPUT(S):\n")
    cat(paste0("  - ", missing_inputs), sep = "\n")
    cat("\n")
  }
  stop("ICBcomb path/input preflight failed. No analysis step was started.", call. = FALSE)
}

required_packages <- c(
  "digest", "readr", "dplyr", "stringr", "tidyr",
  "ggplot2", "patchwork", "scales", "ggrepel", "openxlsx"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install/restore these packages before rerunning.",
    call. = FALSE
  )
}

RUN_AUDIT <- file.path(ROOT, "results", "audit", "ICBcomb_cleanrun_runner")
dir.create(RUN_AUDIT, recursive = TRUE, showWarnings = FALSE)
run_log <- data.frame(
  step = steps,
  status = "NOT_RUN",
  elapsed_seconds = NA_real_,
  stringsAsFactors = FALSE
)

cat("\nICBcomb path/input/package preflight: PASS\n")
cat("Starting 7-step clean rebuild.\n")

# --- 新增: 在工作文件夹(getwd())下创建 logs 文件夹 ---
log_dir <- file.path(getwd(), "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# --- 新增: 将核心运行逻辑包裹在 tryCatch 中 ---
tryCatch({
  
  for (i in seq_along(steps)) {
    s <- steps[[i]]
    cat("\n============================================================\n")
    cat(sprintf("STEP %02d/%02d: %s\n", i, length(steps), s))
    cat("============================================================\n")
    t0 <- proc.time()[["elapsed"]]
    
    ok <- tryCatch(
      {
        source(file.path(ROOT, s), local = new.env(parent = globalenv()))
        TRUE
      },
      error = function(e) {
        run_log$status[[i]] <<- paste0("ERROR: ", conditionMessage(e))
        run_log$elapsed_seconds[[i]] <<- proc.time()[["elapsed"]] - t0
        utils::write.csv(
          run_log,
          file.path(RUN_AUDIT, "ICBcomb_cleanrun_step_status.csv"),
          row.names = FALSE
        )
        stop(
          "ICBcomb runner stopped at: ", s, "\n",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    
    if (isTRUE(ok)) {
      run_log$status[[i]] <- "PASS"
      run_log$elapsed_seconds[[i]] <- proc.time()[["elapsed"]] - t0
      utils::write.csv(
        run_log,
        file.path(RUN_AUDIT, "ICBcomb_cleanrun_step_status.csv"),
        row.names = FALSE
      )
    }
  }
  
  cat("\n============================================================\n")
  cat("ICBcomb TRUE CLEAN-RUN REBUILD COMPLETED\n")
  cat("Resolved module root: ", ROOT, "\n", sep = "")
  cat("Steps passed: ", sum(run_log$status == "PASS"), "/", nrow(run_log), "\n", sep = "")
  cat("Decision: READY_FOR_MODULE_FREEZE\n")
  cat("============================================================\n")
  
  # --- 新增: 无论前面的代码是顺利跑完还是遇到 stop() 中断，都会强制执行 finally 模块 ---
}, finally = {
  session_file <- file.path(log_dir, "sessionInfo_ICBcomb_Runner.txt")
  
  sink(session_file)
  cat("=== RUNNER SESSION INFO ===\n")
  cat("Time:", as.character(Sys.time()), "\n\n")
  print(sessionInfo())
  sink()
  
  cat("\n-> Runner sessionInfo 已经强制保存至:", normalizePath(session_file, winslash = "/", mustWork = FALSE), "\n")
})
