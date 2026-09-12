############################################################
## 09A-05. Doublet-rate sensitivity (5%, 7.5%, 10%)
##
## Reuses the same per-cell pANN ranking and varies only the
## prespecified expected-doublet-rate call. This avoids repeated
## PCA/kNN computation and isolates rate-choice sensitivity.
############################################################

options(stringsAsFactors = FALSE)
SCRIPT_ID <- "09A_05_doublet_rate_sensitivity"

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- getwd()
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

required_pkgs <- c("readr", "dplyr", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)

out_table_dir <- file.path(project_dir, "results", "tables", "GSE244983", "doublet_CNV_sensitivity")
out_obj_dir <- file.path(project_dir, "results", "objects", "GSE244983", "doublet_CNV_sensitivity")
input_rds <- file.path(out_obj_dir, "GSE244983_all_cells_with_pANN_doublet_calls.rds")
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
rate_map <- data.frame(
  ExpectedRate = c(0.05, 0.075, 0.10),
  LabelColumn = c("pANN_doublet_5pct", "pANN_doublet_primary", "pANN_doublet_10pct"),
  stringsAsFactors = FALSE
)

if (!file.exists(input_rds)) stop("09A-02 labeled object missing.", call. = FALSE)
obj <- readRDS(input_rds)
md <- obj@meta.data

## Structural nestedness must hold because all calls use one pANN ranking.
d5 <- rownames(md)[as.character(md$pANN_doublet_5pct) == "Doublet"]
d75 <- rownames(md)[as.character(md$pANN_doublet_primary) == "Doublet"]
d10 <- rownames(md)[as.character(md$pANN_doublet_10pct) == "Doublet"]
nestedness <- data.frame(
  gate = c("5pct_subset_of_7.5pct", "7.5pct_subset_of_10pct"),
  pass = c(all(d5 %in% d75), all(d75 %in% d10)),
  stringsAsFactors = FALSE
)
write_csv(nestedness, "D05_doublet_rate_nestedness_gates.csv")
if (!all(nestedness$pass)) stop("Doublet-rate calls are not nested; inspect D02 classification logic.", call. = FALSE)

all_major <- do.call(rbind, lapply(state_cols, function(st) {
  tmp <- aggregate(md[[st]], by = list(MajorCellType = as.character(md$MajorCellType)), FUN = mean, na.rm = TRUE)
  names(tmp)[2] <- "MeanAll"
  tmp$StateColumn <- st
  tmp
}))

results <- list()
composition_results <- list()
for (i in seq_len(nrow(rate_map))) {
  rate <- rate_map$ExpectedRate[[i]]
  lab <- rate_map$LabelColumn[[i]]
  singlet <- as.character(md[[lab]]) == "Singlet"
  n_doublet <- sum(!singlet)

  for (st in state_cols) {
    tmp <- aggregate(md[[st]][singlet], by = list(MajorCellType = as.character(md$MajorCellType[singlet])), FUN = mean, na.rm = TRUE)
    names(tmp)[2] <- "MeanSinglet"
    base <- all_major[all_major$StateColumn == st, c("MajorCellType", "MeanAll"), drop = FALSE]
    d <- merge(base, tmp, by = "MajorCellType", all = TRUE)
    d <- d[is.finite(d$MeanAll) & is.finite(d$MeanSinglet), , drop = FALSE]
    top_all <- d$MajorCellType[which.max(d$MeanAll)]
    top_sing <- d$MajorCellType[which.max(d$MeanSinglet)]
    results[[length(results) + 1L]] <- data.frame(
      ExpectedRate = rate,
      LabelColumn = lab,
      N_predicted_doublets = n_doublet,
      N_retained_singlets = sum(singlet),
      StateColumn = st,
      State = unname(state_display[[st]]),
      Spearman = suppressWarnings(stats::cor(d$MeanAll, d$MeanSinglet, method = "spearman")),
      Pearson = suppressWarnings(stats::cor(d$MeanAll, d$MeanSinglet, method = "pearson")),
      MaxAbsDelta = max(abs(d$MeanSinglet - d$MeanAll), na.rm = TRUE),
      MeanAbsDelta = mean(abs(d$MeanSinglet - d$MeanAll), na.rm = TRUE),
      TopGroupAll = top_all,
      TopGroupSinglet = top_sing,
      TopGroupUnchanged = identical(as.character(top_all), as.character(top_sing)),
      stringsAsFactors = FALSE
    )
  }

  all_tt <- prop.table(table(as.character(md$MajorCellType)))
  sing_tt <- prop.table(table(as.character(md$MajorCellType[singlet])))
  lev <- union(names(all_tt), names(sing_tt))
  av <- as.numeric(all_tt[lev]); av[is.na(av)] <- 0
  sv <- as.numeric(sing_tt[lev]); sv[is.na(sv)] <- 0
  composition_results[[length(composition_results) + 1L]] <- data.frame(
    ExpectedRate = rate,
    N_predicted_doublets = n_doublet,
    N_retained_singlets = sum(singlet),
    MaxAbsMajorCompositionShift = max(abs(sv - av)),
    MeanAbsMajorCompositionShift = mean(abs(sv - av)),
    stringsAsFactors = FALSE
  )
}

result_df <- do.call(rbind, results)
comp_df <- do.call(rbind, composition_results)
write_csv(result_df, "D05_doublet_rate_state_localization_sensitivity.csv")
write_csv(comp_df, "D05_doublet_rate_composition_sensitivity.csv")

summary <- dplyr::summarise(
  dplyr::group_by(result_df, ExpectedRate),
  N_predicted_doublets = dplyr::first(N_predicted_doublets),
  N_retained_singlets = dplyr::first(N_retained_singlets),
  N_states_top_group_unchanged = sum(TopGroupUnchanged),
  MinSpearman = min(Spearman, na.rm = TRUE),
  MaxAbsDeltaAcrossStates = max(MaxAbsDelta, na.rm = TRUE),
  .groups = "drop"
)
summary <- dplyr::left_join(summary, comp_df, by = c("ExpectedRate", "N_predicted_doublets", "N_retained_singlets"))
write_csv(summary, "D05_doublet_rate_sensitivity_summary.csv")
writeLines(capture.output(sessionInfo()), file.path(out_table_dir, paste0(SCRIPT_ID, "_sessionInfo.txt")))
message("09A-05 complete: 5%, 7.5%, and 10% expected-rate sensitivity summarized.")
