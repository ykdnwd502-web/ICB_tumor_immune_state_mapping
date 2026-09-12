############################################################
## 09A-04. Primary doublet-removal robustness analysis
##
## Compares all 26,053 cells with the primary 7.5%-rate pANN
## singlet subset, using the already frozen four state-score columns.
## No state scores are recalculated here.
############################################################

options(stringsAsFactors = FALSE)
SCRIPT_ID <- "09A_04_primary_doublet_robustness"

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

required_pkgs <- c("readr", "dplyr", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)

out_table_dir <- file.path(project_dir, "results", "tables", "GSE244983", "doublet_CNV_sensitivity")
out_obj_dir <- file.path(project_dir, "results", "objects", "GSE244983", "doublet_CNV_sensitivity")
all_rds <- file.path(out_obj_dir, "GSE244983_all_cells_with_pANN_doublet_calls.rds")
singlet_rds <- file.path(out_obj_dir, "GSE244983_singlets_after_pANN_doublet_filter.rds")

write_csv <- function(x, name) readr::write_csv(x, file.path(out_table_dir, name), na = "")

state_cols <- c(
  "Immune_defective_Cold_score",
  "Myeloid_Treg_Immunosuppressive_score",
  "Tumor_dedifferentiation_Stromal_remodeling_score",
  "Melanocytic_Differentiation_score"
)
state_display <- c(
  "Immune_defective_Cold_score" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive_score" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling_score" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation_score" = "melanocytic differentiation"
)

if (!file.exists(all_rds) || !file.exists(singlet_rds)) stop("09A-03 objects are missing.", call. = FALSE)
obj_all <- readRDS(all_rds)
obj_singlet <- readRDS(singlet_rds)
if (!all(state_cols %in% colnames(obj_all@meta.data))) stop("State-score columns missing from all-cell object.", call. = FALSE)
if (!all(state_cols %in% colnames(obj_singlet@meta.data))) stop("State-score columns missing from singlet object.", call. = FALSE)

summarize_localization <- function(obj, group_col, universe) {
  md <- obj@meta.data
  md$Cell <- rownames(md)
  long <- tidyr::pivot_longer(
    md[, c("Cell", group_col, state_cols), drop = FALSE],
    cols = dplyr::all_of(state_cols), names_to = "StateColumn", values_to = "Score"
  )
  out <- dplyr::summarise(
    dplyr::group_by(long, .data[[group_col]], StateColumn),
    N = sum(is.finite(Score)),
    MeanScore = mean(Score, na.rm = TRUE),
    MedianScore = stats::median(Score, na.rm = TRUE),
    .groups = "drop"
  )
  names(out)[names(out) == group_col] <- "Group"
  out$Grouping <- group_col
  out$Universe <- universe
  out$State <- unname(state_display[out$StateColumn])
  out
}

all_major <- summarize_localization(obj_all, "MajorCellType", "All cells")
sing_major <- summarize_localization(obj_singlet, "MajorCellType", "Retained singlets")
all_cluster <- summarize_localization(obj_all, "SeuratCluster", "All cells")
sing_cluster <- summarize_localization(obj_singlet, "SeuratCluster", "Retained singlets")

localization <- rbind(all_major, sing_major, all_cluster, sing_cluster)
write_csv(localization, "D04_state_localization_all_vs_singlets.csv")

compare_one <- function(all_df, sing_df, grouping) {
  states <- unique(all_df$StateColumn)
  do.call(rbind, lapply(states, function(st) {
    a <- all_df[all_df$StateColumn == st, c("Group", "MeanScore"), drop = FALSE]
    s <- sing_df[sing_df$StateColumn == st, c("Group", "MeanScore"), drop = FALSE]
    names(a)[2] <- "MeanAll"
    names(s)[2] <- "MeanSinglet"
    d <- merge(a, s, by = "Group", all = TRUE)
    d <- d[is.finite(d$MeanAll) & is.finite(d$MeanSinglet), , drop = FALSE]
    top_all <- d$Group[which.max(d$MeanAll)]
    top_singlet <- d$Group[which.max(d$MeanSinglet)]
    data.frame(
      Grouping = grouping,
      StateColumn = st,
      State = unname(state_display[[st]]),
      N_groups = nrow(d),
      Spearman = suppressWarnings(stats::cor(d$MeanAll, d$MeanSinglet, method = "spearman")),
      Pearson = suppressWarnings(stats::cor(d$MeanAll, d$MeanSinglet, method = "pearson")),
      MaxAbsDelta = max(abs(d$MeanSinglet - d$MeanAll), na.rm = TRUE),
      MeanAbsDelta = mean(abs(d$MeanSinglet - d$MeanAll), na.rm = TRUE),
      TopGroupAll = top_all,
      TopGroupSinglet = top_singlet,
      TopGroupUnchanged = identical(as.character(top_all), as.character(top_singlet)),
      stringsAsFactors = FALSE
    )
  }))
}

concordance <- rbind(
  compare_one(all_major, sing_major, "MajorCellType"),
  compare_one(all_cluster, sing_cluster, "SeuratCluster")
)
write_csv(concordance, "D04_state_localization_concordance.csv")

## High-state doublet enrichment: top quartile of each state among all cells.
md <- obj_all@meta.data
is_doublet <- as.character(md$pANN_doublet_primary) == "Doublet"
enrichment <- do.call(rbind, lapply(state_cols, function(st) {
  x <- md[[st]]
  threshold <- as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, type = 7))
  high <- is.finite(x) & x >= threshold
  tab <- table(
    HighState = factor(high, levels = c(FALSE, TRUE)),
    PredictedDoublet = factor(is_doublet, levels = c(FALSE, TRUE))
  )
  ft <- stats::fisher.test(tab)
  est <- if (length(ft$estimate)) unname(ft$estimate[[1]]) else NA_real_
  ci <- if (!is.null(ft$conf.int)) as.numeric(ft$conf.int) else c(NA_real_, NA_real_)
  data.frame(
    StateColumn = st,
    State = unname(state_display[[st]]),
    HighStateThresholdQ75 = threshold,
    N_high = sum(high),
    N_doublet_high = sum(high & is_doublet),
    N_doublet_not_high = sum(!high & is_doublet),
    OddsRatio = est,
    CI_low = ci[[1]],
    CI_high = ci[[2]],
    P_value = ft$p.value,
    stringsAsFactors = FALSE
  )
}))
enrichment$FDR_BH <- stats::p.adjust(enrichment$P_value, method = "BH")
write_csv(enrichment, "D04_high_state_doublet_enrichment.csv")

summary <- data.frame(
  item = c(
    "major_state_comparisons",
    "cluster_state_comparisons",
    "major_top_group_unchanged_n",
    "major_spearman_min",
    "major_spearman_max",
    "major_max_abs_delta_max",
    "interpretation"
  ),
  value = c(
    sum(concordance$Grouping == "MajorCellType"),
    sum(concordance$Grouping == "SeuratCluster"),
    sum(concordance$TopGroupUnchanged[concordance$Grouping == "MajorCellType"]),
    min(concordance$Spearman[concordance$Grouping == "MajorCellType"], na.rm = TRUE),
    max(concordance$Spearman[concordance$Grouping == "MajorCellType"], na.rm = TRUE),
    max(concordance$MaxAbsDelta[concordance$Grouping == "MajorCellType"], na.rm = TRUE),
    "technical_sensitivity_only_no_primary_lineage_replacement"
  ),
  stringsAsFactors = FALSE
)
write_csv(summary, "D04_primary_doublet_robustness_summary.csv")
writeLines(capture.output(sessionInfo()), file.path(out_table_dir, paste0(SCRIPT_ID, "_sessionInfo.txt")))
message("09A-04 complete. Scientific differences, if any, are reported rather than hard-gated.")
