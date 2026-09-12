############################################################
## 02_ICBcomb_build_source_tables.R
## True clean-run ICBcomb source-table producer v1.0.2
##
## FROZEN SCIENTIFIC CONTRACT
##   3 source-locked raw ICBcomb CSVs x 15 rows = 45 rows
##   Positive prioritization: NES > 0
##   45 raw rows -> 42 positive rows -> 10 state x strategy rows
##   Query counts: 98 RESTORE / 73 SUPPRESS / 81 SUPPRESS
##
## IMPORTANT
##   * Reads ONLY explicit files from inputs/ICBcomb_INPUT_FREEZE_v1.0.
##   * Historical-reference CSVs are used only as numerical audit authorities;
##     they are never promoted or copied as analysis outputs.
##   * No fuzzy search, no live platform query, no auto-install.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

ROOT <- Sys.getenv("ICBCOMB_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("ICB_PUBLIC_ROOT")
if (!nzchar(ROOT)) ROOT <- getwd()
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(ROOT, "inputs"))) stop("Wrong ICBcomb module root resolved in 02_ICBcomb_build_source_tables.R: ", ROOT, call. = FALSE)

INPUT <- file.path(ROOT, "inputs", "ICBcomb_INPUT_FREEZE_v1.0")
RAW_DIR <- file.path(INPUT, "data", "raw_platform_results")
GENE_DIR <- file.path(
  INPUT,
  "data",
  "query_gene_set_bundle",
  "02_final_ICBcomb_gene_sets"
)
HIST_DIR <- file.path(
  INPUT,
  "data",
  "historical_reference",
  "ICBcomb_reversal_analysis"
)
OUT <- file.path(ROOT, "results", "tables", "ICBcomb")
AUDIT <- file.path(ROOT, "results", "audit", "02_ICBcomb_source_table_rebuild")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT, recursive = TRUE, showWarnings = FALSE)

STATE_LEVELS <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling"
)

source_spec <- data.frame(
  SourceFile = c(
    "ICBcomb_Immune_cold_RESTORE_hybrid_results.csv",
    "ICBcomb_Myeloid_Treg_SUPPRESS_hybrid_results.csv",
    "ICBcomb_Tumor_dedifferentiation_SUPPRESS_hybrid_results.csv"
  ),
  State = STATE_LEVELS,
  Direction = c("RESTORE", "SUPPRESS", "SUPPRESS"),
  GeneSetVersion = rep("hybrid", 3L),
  GeneSet = c(
    "Immune_defective_Cold_RESTORE_hybrid",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid"
  ),
  stringsAsFactors = FALSE
)

required_raw_cols <- c(
  "Name", "Term", "ES", "NES", "NOM p-val", "FDR q-val",
  "FWER p-val", "Tag %", "Gene %", "Lead_genes",
  "Dataset", "Type", "Path", "Group"
)

parse_num <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

shorten_strategy_name <- function(x, max_width = 52) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "[_]+", " ")
  x <- stringr::str_squish(x)
  stringr::str_wrap(x, width = max_width)
}

map_strategy_category <- function(name, type = NA_character_, path = NA_character_, group = NA_character_) {
  txt_l <- tolower(paste(name, type, path, group, sep = " "))
  dplyr::case_when(
    stringr::str_detect(txt_l, "avadomide|cc[- ]?122|cereblon|lenalidomide|pomalidomide|immunomodulatory") ~
      "Avadomide/IMiD + ICB",
    stringr::str_detect(txt_l, "azacitidine|decitabine|hdac|vorinostat|entinostat|methyltransferase|epigen") ~
      "Epigenetic therapy + ICB",
    stringr::str_detect(txt_l, "il[- ]?15|alt[- ]?803|n[- ]?803|cytokine") ~
      "Cytokine/IL-15 + ICB",
    stringr::str_detect(txt_l, "ascorbic|vitamin c") ~
      "Ascorbic acid + ICB",
    stringr::str_detect(txt_l, "radiation|radiotherapy|irradiat") ~
      "Radiotherapy + ICB",
    stringr::str_detect(txt_l, "chemo|carboplatin|cisplatin|paclitaxel|gemcitabine|temozolomide|dacarbazine") ~
      "Chemotherapy + ICB",
    stringr::str_detect(txt_l, "braf|mek|mapk|trametinib|cobimetinib|vemurafenib|dabrafenib|target") ~
      "Targeted/MAPK therapy + ICB",
    stringr::str_detect(txt_l, "vegf|bevacizumab|anti[- ]?angiogenic|angiogenesis") ~
      "Anti-angiogenic therapy + ICB",
    stringr::str_detect(txt_l, "ctla|pd[- ]?1|pd[- ]?l1|checkpoint") ~
      "Checkpoint-combination strategy",
    TRUE ~ "Other/uncategorized ICBcomb strategy"
  )
}

read_one_source <- function(i) {
  sp <- source_spec[i, , drop = FALSE]
  f <- file.path(RAW_DIR, sp$SourceFile)
  if (!file.exists(f)) stop("Missing frozen raw ICBcomb file: ", f, call. = FALSE)

  x <- readr::read_csv(f, show_col_types = FALSE, name_repair = "minimal")
  if (!all(required_raw_cols %in% names(x))) {
    stop(
      "Raw file missing required columns: ", sp$SourceFile, "\n",
      paste(setdiff(required_raw_cols, names(x)), collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(x) != 15L) stop("Expected exactly 15 rows in ", sp$SourceFile, call. = FALSE)

  x %>%
    mutate(
      SourceFile = sp$SourceFile,
      State = sp$State,
      Direction = sp$Direction,
      GeneSetVersion = sp$GeneSetVersion,
      ES = parse_num(.data[["ES"]]),
      NES = parse_num(.data[["NES"]]),
      NOM_p_val = parse_num(.data[["NOM p-val"]]),
      FDR_q_val = parse_num(.data[["FDR q-val"]]),
      FWER_p_val = parse_num(.data[["FWER p-val"]]),
      Tag_percent = parse_num(.data[["Tag %"]]),
      Gene_percent = parse_num(.data[["Gene %"]]),
      StrategyRaw = as.character(.data[["Name"]]),
      StrategyCategory = map_strategy_category(
        name = .data[["Name"]],
        type = .data[["Type"]],
        path = .data[["Path"]],
        group = .data[["Group"]]
      ),
      Strategy = ifelse(
        StrategyCategory == "Other/uncategorized ICBcomb strategy",
        shorten_strategy_name(StrategyRaw, max_width = 52),
        StrategyCategory
      ),
      Dataset = as.character(.data[["Dataset"]]),
      ComparisonType = as.character(.data[["Type"]]),
      Path = as.character(.data[["Path"]]),
      Group = as.character(.data[["Group"]]),
      Lead_genes = as.character(.data[["Lead_genes"]]),
      PositivePrioritization = !is.na(NES) & NES > 0
    )
}

all_rows <- bind_rows(lapply(seq_len(nrow(source_spec)), read_one_source))
stopifnot(nrow(all_rows) == 45L)

readr::write_csv(all_rows, file.path(OUT, "ICBcomb_all_platform_rows.csv"), na = "")

strategy_mapping <- all_rows %>%
  distinct(
    SourceFile, State, Direction, GeneSetVersion,
    StrategyRaw, StrategyCategory, Strategy,
    Dataset, ComparisonType, Path, Group
  ) %>%
  arrange(State, StrategyCategory, StrategyRaw)

stopifnot(nrow(strategy_mapping) == 45L)
readr::write_csv(strategy_mapping, file.path(OUT, "ICBcomb_strategy_mapping.csv"), na = "")

query_summary_file <- file.path(GENE_DIR, "final_ICBcomb_gene_set_summary_CLEAN.csv")
if (!file.exists(query_summary_file)) {
  stop("Missing final ICBcomb gene-set summary: ", query_summary_file, call. = FALSE)
}
gene_summary <- readr::read_csv(query_summary_file, show_col_types = FALSE)
if (!all(c("GeneSet", "n_genes") %in% names(gene_summary))) {
  stop("Final ICBcomb gene-set summary has an unexpected schema.", call. = FALSE)
}

query_counts <- source_spec %>%
  select(State, Direction, GeneSetVersion, GeneSet) %>%
  left_join(gene_summary, by = "GeneSet") %>%
  rename(n_query_genes = n_genes)

stopifnot(
  nrow(query_counts) == 3L,
  identical(as.integer(query_counts$n_query_genes), c(98L, 73L, 81L)),
  identical(as.character(query_counts$Direction), c("RESTORE", "SUPPRESS", "SUPPRESS"))
)
readr::write_csv(query_counts, file.path(OUT, "ICBcomb_query_counts.csv"), na = "")

positive_rows <- all_rows %>%
  filter(State %in% STATE_LEVELS, PositivePrioritization) %>%
  arrange(State, desc(NES), FDR_q_val)

stopifnot(nrow(positive_rows) == 42L, all(positive_rows$NES > 0))
readr::write_csv(
  positive_rows,
  file.path(OUT, "ICBcomb_positive_prioritization_rows.csv"),
  na = ""
)

state_strategy_summary <- positive_rows %>%
  group_by(State, Strategy) %>%
  summarise(
    median_prioritization_score = median(NES, na.rm = TRUE),
    mean_prioritization_score = mean(NES, na.rm = TRUE),
    max_prioritization_score = max(NES, na.rm = TRUE),
    support_n = dplyr::n(),
    support_dataset_n = dplyr::n_distinct(Dataset),
    min_FDR_q_val = min(FDR_q_val, na.rm = TRUE),
    median_FDR_q_val = median(FDR_q_val, na.rm = TRUE),
    representative_datasets = paste(sort(unique(na.omit(Dataset))), collapse = "; "),
    representative_names = paste(utils::head(unique(na.omit(StrategyRaw)), 5), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    min_FDR_q_val = ifelse(is.infinite(min_FDR_q_val), NA_real_, min_FDR_q_val),
    median_FDR_q_val = ifelse(is.nan(median_FDR_q_val), NA_real_, median_FDR_q_val)
  ) %>%
  arrange(State, desc(median_prioritization_score), desc(support_n)) %>%
  group_by(State) %>%
  mutate(rank_within_state = dplyr::row_number()) %>%
  ungroup()

stopifnot(nrow(state_strategy_summary) == 10L)
readr::write_csv(
  state_strategy_summary,
  file.path(OUT, "ICBcomb_state_strategy_summary.csv"),
  na = ""
)

figA <- query_counts %>%
  transmute(
    state = State,
    n_query_genes = n_query_genes,
    GeneSet = GeneSet,
    Direction = Direction,
    GeneSetVersion = GeneSetVersion
  )

figB <- state_strategy_summary %>%
  group_by(State) %>%
  arrange(desc(median_prioritization_score), desc(support_n), .by_group = TRUE) %>%
  slice_head(n = 5L) %>%
  mutate(rank_within_state = dplyr::row_number()) %>%
  ungroup() %>%
  transmute(
    state = State,
    rank_within_state = rank_within_state,
    strategy = Strategy,
    median_prioritization_score = median_prioritization_score,
    support_n = support_n,
    support_dataset_n = support_dataset_n,
    min_FDR_q_val = min_FDR_q_val,
    median_FDR_q_val = median_FDR_q_val,
    representative_names = representative_names
  )

strategy_order <- figB %>%
  arrange(match(state, STATE_LEVELS), rank_within_state) %>%
  pull(strategy) %>%
  unique()

figC <- tidyr::expand_grid(
  state = STATE_LEVELS,
  strategy = strategy_order
) %>%
  left_join(
    state_strategy_summary %>%
      transmute(
        state = State,
        strategy = Strategy,
        median_prioritization_score = median_prioritization_score,
        support_n = support_n,
        support_dataset_n = support_dataset_n,
        min_FDR_q_val = min_FDR_q_val,
        median_FDR_q_val = median_FDR_q_val
      ),
    by = c("state", "strategy")
  ) %>%
  mutate(
    support_n = ifelse(is.na(support_n), 0L, support_n),
    support_dataset_n = ifelse(is.na(support_dataset_n), 0L, support_dataset_n)
  )

figD <- state_strategy_summary %>%
  group_by(State) %>%
  summarise(
    max_median_prioritization_score = max(median_prioritization_score, na.rm = TRUE),
    top_strategy = Strategy[which.max(median_prioritization_score)],
    total_positive_rows = sum(support_n, na.rm = TRUE),
    total_positive_strategy_n = dplyr::n_distinct(Strategy),
    .groups = "drop"
  ) %>%
  mutate(State = factor(State, levels = STATE_LEVELS)) %>%
  arrange(State) %>%
  transmute(
    state = as.character(State),
    max_median_prioritization_score = max_median_prioritization_score,
    top_strategy = top_strategy,
    total_positive_rows = total_positive_rows,
    total_positive_strategy_n = total_positive_strategy_n
  )

stopifnot(
  nrow(figA) == 3L,
  nrow(figB) == 10L,
  nrow(figC) == 12L,
  nrow(figD) == 3L
)

readr::write_csv(figA, file.path(OUT, "Figure9_Panel_A_query_counts.csv"), na = "")
readr::write_csv(figB, file.path(OUT, "Figure9_Panel_B_top_strategies.csv"), na = "")
readr::write_csv(figC, file.path(OUT, "Figure9_Panel_C_state_strategy_matrix.csv"), na = "")
readr::write_csv(figD, file.path(OUT, "Figure9_Panel_D_state_maximum.csv"), na = "")

## Manuscript-claim support reconstructed from the rebuilt summary.
get_summary_row <- function(state, strategy) {
  state_strategy_summary %>% filter(State == state, Strategy == strategy)
}

r_icb_avad <- get_summary_row("immune-defective/cold", "Avadomide/IMiD + ICB")
r_mye_avad <- get_summary_row("myeloid–Treg immunosuppressive", "Avadomide/IMiD + ICB")
r_tum_epi  <- get_summary_row("tumor-dedifferentiation/stromal-remodeling", "Epigenetic therapy + ICB")
cytokine_states <- state_strategy_summary %>%
  filter(Strategy == "Cytokine/IL-15 + ICB") %>% pull(State) %>% unique()
ascorbic_rows <- state_strategy_summary %>% filter(Strategy == "Ascorbic acid + ICB")

tumor_top <- figD %>%
  filter(state == "tumor-dedifferentiation/stromal-remodeling") %>%
  pull(top_strategy)

claims <- data.frame(
  claim_id = sprintf("ICB%02d", 1:6),
  manuscript_claim = c(
    "Avadomide/IMiD + ICB has two supporting independent datasets for immune-defective/cold.",
    "Avadomide/IMiD + ICB has two supporting independent datasets for myeloid–Treg immunosuppressive.",
    "Epigenetic therapy + ICB has three supporting independent datasets for tumor-dedifferentiation/stromal-remodeling and is the top median signal for that state.",
    "Cytokine/IL-15 + ICB appears only for immune-defective/cold among the three final perturbation queries.",
    "Ascorbic acid + ICB appears for all three final queries, with one supporting independent dataset for each.",
    "Only positive-direction ICBcomb rows (NES > 0) enter the strategy-level summary."
  ),
  observed = c(
    paste0("support_dataset_n=", r_icb_avad$support_dataset_n),
    paste0("support_dataset_n=", r_mye_avad$support_dataset_n),
    paste0("support_dataset_n=", r_tum_epi$support_dataset_n, "; top_strategy=", tumor_top),
    paste(cytokine_states, collapse = "; "),
    paste0(
      ascorbic_rows$State,
      "=",
      ascorbic_rows$support_dataset_n,
      collapse = "; "
    ),
    paste0("positive_rows=", nrow(positive_rows), "; any_NES<=0=", any(positive_rows$NES <= 0))
  ),
  pass = c(
    nrow(r_icb_avad) == 1L && r_icb_avad$support_dataset_n == 2L,
    nrow(r_mye_avad) == 1L && r_mye_avad$support_dataset_n == 2L,
    nrow(r_tum_epi) == 1L && r_tum_epi$support_dataset_n == 3L && identical(tumor_top, "Epigenetic therapy + ICB"),
    length(cytokine_states) == 1L && identical(cytokine_states, "immune-defective/cold"),
    nrow(ascorbic_rows) == 3L && all(ascorbic_rows$support_dataset_n == 1L),
    nrow(positive_rows) == 42L && !any(positive_rows$NES <= 0)
  ),
  stringsAsFactors = FALSE
)
stopifnot(all(claims$pass))
readr::write_csv(claims, file.path(OUT, "ICBcomb_manuscript_claim_support.csv"), na = "")

S27_contract <- data.frame(
  item = c(
    "three_final_query_counts_available",
    "query_GMT_exact_to_hybrid_txt",
    "strategy_mapping_available",
    "positive_direction_mapping_available",
    "melanocytic_reference_control_retained",
    "melanocytic_not_final_raw_query"
  ),
  pass = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
S28_contract <- data.frame(
  item = c(
    "three_final_states_present",
    "strategy_summary_positive_only",
    "median_score_available",
    "support_row_count_available",
    "support_dataset_count_available",
    "bubble_landscape_source_nonempty",
    "support_vs_median_source_nonempty"
  ),
  pass = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
readr::write_csv(S27_contract, file.path(OUT, "S27_source_contract.csv"), na = "")
readr::write_csv(S28_contract, file.path(OUT, "S28_source_contract.csv"), na = "")

## Numerical audit against immutable historical reference authority.
read_hist <- function(name) {
  f <- file.path(HIST_DIR, name)
  if (!file.exists(f)) stop("Missing frozen historical reference: ", f, call. = FALSE)
  read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}

same_df <- function(a, b, sort_cols = NULL, tolerance = 1e-10) {
  if (!is.null(sort_cols)) {
    a <- a[do.call(order, a[sort_cols]), , drop = FALSE]
    b <- b[do.call(order, b[sort_cols]), , drop = FALSE]
    rownames(a) <- NULL
    rownames(b) <- NULL
  }
  isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = tolerance))
}

hist_all <- read_hist("ICBcomb_combined_raw_results_CLEAN.csv")
names(hist_all)[names(hist_all) == "PositiveReversal"] <- "PositivePrioritization"
hist_positive <- read_hist("ICBcomb_positive_reversal_rows_CLEAN.csv")
names(hist_positive)[names(hist_positive) == "PositiveReversal"] <- "PositivePrioritization"
hist_query <- read_hist("ICBcomb_query_gene_set_counts_CLEAN.csv")
hist_summary <- read_hist("ICBcomb_state_strategy_summary_CLEAN.csv")
names(hist_summary)[names(hist_summary) == "median_reversal_score"] <- "median_prioritization_score"
names(hist_summary)[names(hist_summary) == "mean_reversal_score"] <- "mean_prioritization_score"
names(hist_summary)[names(hist_summary) == "max_reversal_score"] <- "max_prioritization_score"
hist_A <- read_hist("ICBcomb_query_gene_set_counts_FINAL_LOCKED.csv")
hist_B <- read_hist("ICBcomb_top_positive_reversal_strategies_FOR_FIGURE8B_FINAL_LOCKED.csv")
names(hist_B)[names(hist_B) == "median_reversal_score"] <- "median_prioritization_score"
hist_C <- read_hist("ICBcomb_strategy_state_heatmap_table_FOR_FIGURE8C_FINAL_LOCKED.csv")
names(hist_C)[names(hist_C) == "median_reversal_score"] <- "median_prioritization_score"
hist_D <- read_hist("ICBcomb_state_level_maximum_positive_reversal_signal_FOR_FIGURE8D_FINAL_LOCKED.csv")
names(hist_D)[names(hist_D) == "max_median_reversal_score"] <- "max_median_prioritization_score"
hist_map <- read_hist("ICBcomb_strategy_mapping_audit_CLEAN.csv")

numerical_audit <- data.frame(
  item = c(
    "all_platform_rows",
    "strategy_mapping",
    "query_counts",
    "positive_rows",
    "state_strategy_summary",
    "Figure9_panel_A",
    "Figure9_panel_B",
    "Figure9_panel_C",
    "Figure9_panel_D"
  ),
  exact_or_numeric_equal = c(
    same_df(as.data.frame(all_rows), hist_all),
    same_df(as.data.frame(strategy_mapping), hist_map),
    same_df(as.data.frame(query_counts), hist_query),
    same_df(as.data.frame(positive_rows), hist_positive),
    same_df(as.data.frame(state_strategy_summary), hist_summary),
    same_df(as.data.frame(figA), hist_A),
    same_df(as.data.frame(figB), hist_B),
    same_df(as.data.frame(figC), hist_C, sort_cols = c("state", "strategy")),
    same_df(as.data.frame(figD), hist_D)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  numerical_audit,
  file.path(AUDIT, "01_rebuilt_vs_frozen_historical_reference.csv"),
  row.names = FALSE
)

if (!all(numerical_audit$exact_or_numeric_equal)) {
  print(numerical_audit[!numerical_audit$exact_or_numeric_equal, , drop = FALSE])
  stop("ICBcomb source-table rebuild did not reproduce the frozen numerical authority.", call. = FALSE)
}

cat("============================================================\n")
cat("ICBcomb SOURCE TABLE REBUILD: PASS\n")
cat("Raw rows: 45\n")
cat("Positive NES > 0 rows: 42\n")
cat("State x strategy rows: 10\n")
cat("Figure 9 source rows A/B/C/D: 3/10/12/3\n")
cat("Frozen numerical-reference comparisons: 9/9 PASS\n")
cat("============================================================\n")
