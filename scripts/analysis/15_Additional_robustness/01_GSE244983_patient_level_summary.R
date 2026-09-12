############################################################
## 01_GSE244983_patient_level_summary.R
##
## Public revision analysis for Supplementary Figure S29.
## Produces patient-level source tables from the frozen primary
## GSE244983 26,053-cell state-localization lineage.
##
## IMPORTANT:
##   - reads only canonical 06/07 public outputs;
##   - does NOT use the historical 25,972-cell doublet/inferCNV branch;
##   - performs no new hypothesis test and changes no state definition.
############################################################

options(stringsAsFactors = FALSE)

required_pkgs <- c("readr", "dplyr", "tidyr", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) stop("Missing package(s): ", paste(missing_pkgs, collapse=", "), call.=FALSE)

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR", unset = "D:/ICB_resistance_project")
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)

PREPROCESS_DIR <- file.path(PROJECT_DIR, "results", "tables", "GSE244983", "preprocess_major_annotation")
STATE_DIR <- file.path(PROJECT_DIR, "results", "tables", "GSE244983", "state_localization_source_attribution")
OUT_DIR <- file.path(PROJECT_DIR, "results", "tables", "additional_robustness", "S29_patient_level")
LOG_DIR <- file.path(PROJECT_DIR, "logs")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(LOG_DIR, recursive=TRUE, showWarnings=FALSE)

CELL_FILE <- file.path(STATE_DIR, "GSE244983_single_cell_state_scores.csv")
COUNT_FILE <- file.path(PREPROCESS_DIR, "GSE244983_sample_cell_counts.csv")

EXPECTED_N <- 26053L
EXPECTED_PATIENTS <- 4L
EXPECTED_COUNTS <- c(
  "Pat_ICBnaive1"=7414L,
  "Pat_ICBnaive2"=9081L,
  "Pat42"=7188L,
  "Pat5"=2370L
)
PATIENT_ORDER <- c("Pat_ICBnaive2", "Pat_ICBnaive1", "Pat42", "Pat5")
STATE_COLS <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)
STATE_LABELS <- c(
  "Immune_defective_Cold"="immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive"="myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling"="tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation"="melanocytic differentiation"
)
TARGET_STATE <- "Tumor_dedifferentiation_Stromal_remodeling"

assert_true <- function(x, msg) if (!isTRUE(x)) stop(msg, call.=FALSE)
md5_file <- function(path) if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_
safe_csv <- function(x, name) readr::write_csv(x, file.path(OUT_DIR, name), na="")

input_manifest <- data.frame(
  Input=c("Canonical state-localization cell table", "Canonical preprocessing sample-count table"),
  Path=c(CELL_FILE, COUNT_FILE),
  Exists=file.exists(c(CELL_FILE, COUNT_FILE)),
  MD5=vapply(c(CELL_FILE, COUNT_FILE), md5_file, character(1)),
  stringsAsFactors=FALSE
)
safe_csv(input_manifest, "S29_input_manifest.csv")
assert_true(all(input_manifest$Exists), "Missing canonical S29 input(s).")

cell_df <- readr::read_csv(CELL_FILE, show_col_types=FALSE, progress=FALSE) |> as.data.frame(check.names=FALSE)
required <- c("Cell", "Sample", "SeuratCluster", "MajorCellType", STATE_COLS)
assert_true(all(required %in% colnames(cell_df)), paste("Missing columns:", paste(setdiff(required, colnames(cell_df)), collapse=", ")))
for (st in STATE_COLS) cell_df[[st]] <- suppressWarnings(as.numeric(cell_df[[st]]))
cell_df$Cell <- as.character(cell_df$Cell)
cell_df$Sample <- as.character(cell_df$Sample)

count_df <- readr::read_csv(COUNT_FILE, show_col_types=FALSE, progress=FALSE) |> as.data.frame(check.names=FALSE)
assert_true(all(c("Sample","N_cells") %in% colnames(count_df)), "Sample-count file requires Sample and N_cells.")
count_df$Sample <- as.character(count_df$Sample)
count_df$N_cells <- as.integer(count_df$N_cells)

obs <- table(cell_df$Sample)
obs <- setNames(as.integer(obs), names(obs))
obs07 <- setNames(count_df$N_cells, count_df$Sample)
finite_state <- vapply(STATE_COLS, function(st) all(is.finite(cell_df[[st]])), logical(1))

lineage_gates <- data.frame(
  Gate=c(
    "Primary cell count",
    "Unique cell identifiers",
    "Exactly four patients",
    "Patient labels",
    "Patient counts from cell table",
    "Patient counts equal preprocessing source",
    "All four state scores finite",
    "Historical 25,972-cell support object not used"
  ),
  Observed=c(
    nrow(cell_df),
    length(unique(cell_df$Cell)),
    length(unique(cell_df$Sample)),
    paste(sort(unique(cell_df$Sample)), collapse=";"),
    paste(paste0(names(EXPECTED_COUNTS), "=", obs[names(EXPECTED_COUNTS)]), collapse=";"),
    paste(paste0(names(EXPECTED_COUNTS), "=", obs07[names(EXPECTED_COUNTS)]), collapse=";"),
    paste(paste0(STATE_COLS, "=", finite_state), collapse=";"),
    "TRUE"
  ),
  Expected=c(
    EXPECTED_N,
    EXPECTED_N,
    EXPECTED_PATIENTS,
    paste(sort(names(EXPECTED_COUNTS)), collapse=";"),
    paste(paste0(names(EXPECTED_COUNTS), "=", EXPECTED_COUNTS), collapse=";"),
    paste(paste0(names(EXPECTED_COUNTS), "=", EXPECTED_COUNTS), collapse=";"),
    "TRUE for all states",
    "TRUE"
  ),
  Pass=c(
    nrow(cell_df)==EXPECTED_N,
    length(unique(cell_df$Cell))==EXPECTED_N,
    length(unique(cell_df$Sample))==EXPECTED_PATIENTS,
    setequal(unique(cell_df$Sample), names(EXPECTED_COUNTS)),
    setequal(names(obs), names(EXPECTED_COUNTS)) && identical(as.integer(obs[names(EXPECTED_COUNTS)]), as.integer(EXPECTED_COUNTS)),
    setequal(names(obs07), names(EXPECTED_COUNTS)) && identical(as.integer(obs07[names(EXPECTED_COUNTS)]), as.integer(EXPECTED_COUNTS)),
    all(finite_state),
    TRUE
  ),
  stringsAsFactors=FALSE
)
safe_csv(lineage_gates, "S29_lineage_gates.csv")
if (any(!lineage_gates$Pass)) stop("S29 primary-lineage gate failure: ", paste(lineage_gates$Gate[!lineage_gates$Pass], collapse=" | "), call.=FALSE)

patient_counts <- cell_df |>
  dplyr::count(Sample, name="n_cells") |>
  dplyr::mutate(
    proportion=n_cells/sum(n_cells),
    proportion_pct=scales::percent(proportion, accuracy=0.1),
    Sample=factor(Sample, levels=PATIENT_ORDER)
  ) |>
  dplyr::arrange(Sample) |>
  dplyr::mutate(Sample=as.character(Sample))

patient_means <- cell_df |>
  dplyr::group_by(Sample) |>
  dplyr::summarise(
    n_cells=dplyr::n(),
    Immune_defective_Cold_mean=mean(Immune_defective_Cold, na.rm=TRUE),
    Myeloid_Treg_Immunosuppressive_mean=mean(Myeloid_Treg_Immunosuppressive, na.rm=TRUE),
    Tumor_dedifferentiation_Stromal_remodeling_mean=mean(Tumor_dedifferentiation_Stromal_remodeling, na.rm=TRUE),
    Melanocytic_Differentiation_mean=mean(Melanocytic_Differentiation, na.rm=TRUE),
    .groups="drop"
  ) |>
  dplyr::mutate(Sample=factor(Sample, levels=PATIENT_ORDER)) |>
  dplyr::arrange(Sample) |>
  dplyr::mutate(Sample=as.character(Sample))

patient_means_long <- patient_means |>
  tidyr::pivot_longer(
    cols=dplyr::ends_with("_mean"),
    names_to="State",
    values_to="MeanScore"
  ) |>
  dplyr::mutate(
    State=sub("_mean$", "", State),
    StateLabel=unname(STATE_LABELS[State])
  )

panelC <- cell_df |>
  dplyr::transmute(Cell=Cell, Sample=Sample, Tumor_dedifferentiation_Stromal_remodeling_score=.data[[TARGET_STATE]])

safe_csv(patient_counts, "S29_patient_counts.csv")
safe_csv(patient_means, "S29_patient_state_means.csv")
safe_csv(patient_means_long, "S29_patient_state_means_long.csv")
safe_csv(panelC, "S29_target_state_cell_scores.csv")

method_record <- data.frame(
  Item=c("Role","Primary lineage","Patients","Patient-level summaries","Hypothesis testing","Historical support branch"),
  Value=c(
    "Reviewer-requested patient-level structure audit",
    "26,053-cell primary GSE244983 annotation/state-localization universe",
    "4",
    "cell contribution; four-state patient means; target-state cell distribution",
    "None; descriptive structure audit",
    "25,972-cell doublet/inferCNV branch is intentionally excluded"
  ),
  stringsAsFactors=FALSE
)
safe_csv(method_record, "S29_method_record.csv")

writeLines(capture.output(sessionInfo()), file.path(LOG_DIR, "sessionInfo_15_01_GSE244983_patient_level_summary.txt"))

cat("\n============================================================\n")
cat("15/01 GSE244983 PATIENT-LEVEL SUMMARY\n")
cat("Primary cells: ", nrow(cell_df), "\n", sep="")
cat("Patients: ", length(unique(cell_df$Sample)), "\n", sep="")
cat("Lineage gates: ", sum(lineage_gates$Pass), "/", nrow(lineage_gates), "\n", sep="")
cat("Historical 25,972-cell support branch used: NO\n")
cat("Overall: PASS\n")
cat("============================================================\n")
