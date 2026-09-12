############################################################
## 02_make_S30_cross_cohort_robustness.R
##
## Project: ICB_resistance_project
##
## PURPOSE
##   Canonical reporting-only rebuild of Supplementary Figure S30.
##
## AUTHORITATIVE INPUTS
##   Frozen script 05:
##     results/tables/bulk_response_sensitivity/
##
##   Frozen script 19:
##     results/tables/cross_cohort_harmonization/
##
## IMPORTANT
##   - Reporting only: this script does NOT recompute state scores,
##     residualization models, AUCs, bootstrap distributions, or
##     permutations.
##   - Historical 18H tables are NOT used as inferential authority.
##   - No broad grep across all columns.
##   - No flattening of arbitrary numeric columns.
##   - No averaging across heterogeneous composition-adjustment models.
##
## PANELS
##   A. Four-state raw GSE78220 association AUCs (script 05 / S11).
##   B. Target-state GSE78220 leave-one-sample-out stability (script 05 / S11).
##   C. Target-state prespecified composition-adjusted AUCs (script 05 / S13).
##   D. Four-state cross-cohort canonical fixed-direction AUCs with DeLong CI
##      (script 19); exact within-cohort permutation P values are retained
##      in source data.
##   E. Direct target-state point-estimate comparison:
##      GSE78220 raw vs cross-cohort harmonized.
##
## INTERPRETATION
##   Descriptive robustness / sensitivity analysis.
##   Not predictive model validation, biomarker validation, or an
##   independent validation cohort.
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################################################
## 0. Project root
############################################################

resolve_project_dir <- function() {
  env_dir <- Sys.getenv("ICB_PROJECT_DIR", unset = "")

  if (nzchar(env_dir) && dir.exists(env_dir)) {
    return(
      normalizePath(
        env_dir,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  wd <- normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )

  if (
    dir.exists(file.path(wd, "scripts")) &&
      dir.exists(file.path(wd, "results"))
  ) {
    return(wd)
  }

  stop(
    "Cannot determine project root. Set:\n",
    'Sys.setenv(ICB_PROJECT_DIR = "D:/ICB_resistance_project")'
  )
}

project_dir <- resolve_project_dir()

bulk_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "bulk_response_sensitivity"
)

cross_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "cross_cohort_harmonization"
)

out_root <- file.path(
  project_dir,
  "results",
  "reporting",
  "additional_robustness",
  "S30_cross_cohort"
)

figure_dir <- file.path(project_dir, "results", "figures", "supplementary")
table_dir  <- file.path(out_root, "tables")
log_dir    <- file.path(out_root, "logs")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", project_dir)
message("S30 reporting output: ", out_root)

############################################################
## 1. Packages
############################################################

required_pkgs <- c(
  "readr",
  "dplyr",
  "ggplot2",
  "patchwork",
  "openxlsx"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", ")
  )
}

# === 新增代码：实际加载这些包，以暴露 %>% 管道符和其他重载运算符 ===
invisible(lapply(required_pkgs, library, character.only = TRUE))


############################################################
## 2. Canonical state definitions
############################################################

STATE_COLS <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

STATE_DISPLAY <- c(
  "Immune_defective_Cold" =
    "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" =
    "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" =
    "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" =
    "melanocytic differentiation"
)

STATE_DISPLAY_PLOT <- c(
  "Immune_defective_Cold" =
    "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" =
    "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" =
    "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_Differentiation" =
    "melanocytic\ndifferentiation"
)

STATE_COLORS <- c(
  "Immune_defective_Cold" = "#4DBBD5",
  "Myeloid_Treg_Immunosuppressive" = "#00A087",
  "Tumor_dedifferentiation_Stromal_remodeling" = "#E64B35",
  "Melanocytic_Differentiation" = "#3C5488"
)

TARGET_STATE <-
  "Tumor_dedifferentiation_Stromal_remodeling"

TARGET_DISPLAY <-
  unname(STATE_DISPLAY[[TARGET_STATE]])

TARGET_COLOR <-
  unname(STATE_COLORS[[TARGET_STATE]])

## Frozen S13 display set.
S13_METHOD_ORDER <- c(
  "Raw_score",
  "CAF_stromal_residualized",
  "CAF_plus_melanoma_lineage_residualized",
  "All_marker_composition_residualized",
  "Purity_proxy_residualized",
  "Marker_composition_PC1_PC2_residualized"
)

S13_METHOD_DISPLAY <- c(
  "Raw_score" = "Raw",
  "CAF_stromal_residualized" = "CAF/stromal",
  "CAF_plus_melanoma_lineage_residualized" =
    "CAF + melanoma-lineage",
  "All_marker_composition_residualized" =
    "All marker composition",
  "Purity_proxy_residualized" =
    "Purity proxy",
  "Marker_composition_PC1_PC2_residualized" =
    "Composition PC"
)

############################################################
## 3. Helpers
############################################################

assert_true <- function(x, msg) {
  if (!isTRUE(x)) {
    stop(msg, call. = FALSE)
  }
}

read_csv_required <- function(file, label) {
  if (!file.exists(file)) {
    stop(
      label,
      " not found:\n",
      file,
      call. = FALSE
    )
  }

  readr::read_csv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )
}

require_cols <- function(df, cols, label) {
  missing <- setdiff(cols, colnames(df))

  if (length(missing) > 0L) {
    stop(
      label,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

parse_pass <- function(x) {
  if (is.logical(x)) return(x)

  tolower(trimws(as.character(x))) %in%
    c(
      "true",
      "t",
      "1",
      "pass",
      "passed"
    )
}

md5_file <- function(file) {
  if (!file.exists(file)) return(NA_character_)
  unname(tools::md5sum(file))
}

write_csv_safe <- function(x, file) {
  readr::write_csv(
    x,
    file,
    na = ""
  )
  message("Saved: ", file)
}

publication_theme <- function(base_size = 12.5) {
  ggplot2::theme_bw(
    base_size = base_size,
    base_family = "sans"
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0, family = "sans"),
      plot.subtitle = ggplot2::element_text(size = 11.25, hjust = 0, family = "sans"),
      plot.tag = ggplot2::element_text(face = "bold", size = 16, family = "sans"),
      axis.title = ggplot2::element_text(size = 12.5, face = "bold", family = "sans"),
      axis.text = ggplot2::element_text(size = 11, color = "black", family = "sans"),
      legend.title = ggplot2::element_text(size = 12, face = "bold", family = "sans"),
      legend.text = ggplot2::element_text(size = 10.5, family = "sans"),
      strip.text = ggplot2::element_text(size = 12, face = "bold", family = "sans"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(7, 8, 7, 8)
    )
}

############################################################
## 4. Exact canonical input files
############################################################

s05_gate_file <- file.path(
  bulk_dir,
  "05_reproducibility_gates.csv"
)

s11_auc_file <- file.path(
  bulk_dir,
  "S11_GSE78220_raw_auc_summary.csv"
)

s11_loso_file <- file.path(
  bulk_dir,
  "S11_GSE78220_loso_auc_values.csv"
)

s13_auc_file <- file.path(
  bulk_dir,
  "S13_formal_composition_adjusted_auc_summary.csv"
)

s19_gate_file <- file.path(
  cross_dir,
  "cross_cohort_reproducibility_gates.csv"
)

s19_auc_file <- file.path(
  cross_dir,
  "cross_cohort_fixed_direction_AUC_DeLong_summary.csv"
)

s19_boot_file <- file.path(
  cross_dir,
  "cross_cohort_seeded_bootstrap_AUC_summary.csv"
)

s19_perm_file <- file.path(
  cross_dir,
  "cross_cohort_seeded_within_cohort_permutation_summary.csv"
)

s19_method_file <- file.path(
  cross_dir,
  "cross_cohort_method_manifest.csv"
)

input_files <- c(
  s05_gate_file,
  s11_auc_file,
  s11_loso_file,
  s13_auc_file,
  s19_gate_file,
  s19_auc_file,
  s19_boot_file,
  s19_perm_file,
  s19_method_file
)

if (any(!file.exists(input_files))) {
  stop(
    "Missing S30 canonical input(s):\n",
    paste(
      input_files[!file.exists(input_files)],
      collapse = "\n"
    )
  )
}

input_manifest <- data.frame(
  SourceLayer = c(
    rep("Frozen script 05", 4),
    rep("15/02 cross-cohort comparability", 5)
  ),
  File = basename(input_files),
  Path = input_files,
  SizeBytes = file.info(input_files)$size,
  MD5 = vapply(
    input_files,
    md5_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

############################################################
## 5. Upstream reproducibility gates must already PASS
############################################################

s05_gates <- read_csv_required(
  s05_gate_file,
  "Script-05 reproducibility gates"
)

s19_gates <- read_csv_required(
  s19_gate_file,
  "Script-19 reproducibility gates"
)

require_cols(
  s05_gates,
  c("Pass"),
  "Script-05 gates"
)

require_cols(
  s19_gates,
  c("Pass"),
  "Script-19 gates"
)

assert_true(
  all(parse_pass(s05_gates$Pass)),
  "Script-05 upstream reproducibility gates are not all PASS."
)

assert_true(
  all(parse_pass(s19_gates$Pass)),
  "Script-19 upstream reproducibility gates are not all PASS."
)

############################################################
## 6. Read canonical panel-source tables
############################################################

s11_auc <- read_csv_required(
  s11_auc_file,
  "S11 canonical raw AUC"
)

s11_loso <- read_csv_required(
  s11_loso_file,
  "S11 canonical LOSO values"
)

s13_auc <- read_csv_required(
  s13_auc_file,
  "S13 canonical composition-adjusted AUC summary"
)

s19_auc <- read_csv_required(
  s19_auc_file,
  "Script-19 canonical AUC summary"
)

s19_boot <- read_csv_required(
  s19_boot_file,
  "Script-19 seeded bootstrap summary"
)

s19_perm <- read_csv_required(
  s19_perm_file,
  "Script-19 exact permutation summary"
)

s19_method <- read_csv_required(
  s19_method_file,
  "Script-19 method manifest"
)

require_cols(
  s11_auc,
  c(
    "State",
    "StateLabel",
    "AUC",
    "n_samples",
    "n_responders",
    "n_nonresponders"
  ),
  "S11 canonical raw AUC"
)

require_cols(
  s11_loso,
  c(
    "State",
    "StateLabel",
    "RemovedSample",
    "RemovedResponseGroup",
    "Full_AUC",
    "LeaveOneOut_AUC",
    "Delta_AUC_minus_full"
  ),
  "S11 canonical LOSO"
)

require_cols(
  s13_auc,
  c(
    "cohort",
    "state",
    "state_internal",
    "adjustment",
    "auc",
    "bootstrap_ci_low",
    "bootstrap_ci_high",
    "loso_auc_median",
    "loso_auc_min",
    "loso_auc_max"
  ),
  "S13 canonical AUC summary"
)

require_cols(
  s19_auc,
  c(
    "State",
    "StateDisplay",
    "N",
    "Responders",
    "NonResponders",
    "AUC",
    "CI_low",
    "CI_high",
    "Direction"
  ),
  "Script-19 canonical AUC summary"
)

require_cols(
  s19_perm,
  c(
    "State",
    "Observed_AUC",
    "P_two_sided_from_0_5",
    "N_perm",
    "N_extreme_exact"
  ),
  "Script-19 permutation summary"
)

############################################################
## 7. Panel A — GSE78220 raw four-state AUCs
############################################################

panelA <- s11_auc %>%
  dplyr::filter(
    .data$State %in% STATE_COLS
  ) %>%
  dplyr::transmute(
    State = as.character(.data$State),
    StateDisplay =
      unname(
        STATE_DISPLAY[
          as.character(.data$State)
        ]
      ),
    StateDisplayPlot =
      unname(
        STATE_DISPLAY_PLOT[
          as.character(.data$State)
        ]
      ),
    AUC = as.numeric(.data$AUC),
    N = as.integer(.data$n_samples),
    Responders =
      as.integer(.data$n_responders),
    NonResponders =
      as.integer(.data$n_nonresponders)
  )

assert_true(
  nrow(panelA) == 4L,
  "Panel A must contain exactly four states."
)

assert_true(
  setequal(
    panelA$State,
    STATE_COLS
  ),
  "Panel A state universe differs from the canonical four states."
)

panelA$StateDisplay <- factor(
  panelA$StateDisplay,
  levels = rev(unname(STATE_DISPLAY[STATE_COLS]))
)

panelA$StateDisplayPlot <- factor(
  panelA$StateDisplayPlot,
  levels = rev(unname(STATE_DISPLAY_PLOT[STATE_COLS]))
)

pA <- ggplot2::ggplot(
  panelA,
  ggplot2::aes(
    x = .data$StateDisplayPlot,
    y = .data$AUC
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = .data$State),
    size = 3.2
  ) +
  ggplot2::scale_color_manual(values = STATE_COLORS, guide = "none") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.3f",
        .data$AUC
      )
    ),
    hjust = -0.35,
    size = 3.0
  ) +
  ggplot2::coord_flip(
    ylim = c(
      min(0.35, min(panelA$AUC) - 0.03),
      min(1.0, max(panelA$AUC) + 0.09)
    )
  ) +
  ggplot2::labs(
    title =
      "GSE78220 raw state-score associations",
    subtitle =
      "Strict binary response subset; n = 27",
    x = NULL,
    y = "AUC for non-response"
  ) +
  publication_theme()

############################################################
## 8. Panel B — target-state LOSO
############################################################

panelB <- s11_loso %>%
  dplyr::filter(
    .data$State == TARGET_STATE
  ) %>%
  dplyr::transmute(
    State = as.character(.data$State),
    StateDisplay = TARGET_DISPLAY,
    RemovedSample =
      as.character(.data$RemovedSample),
    RemovedResponseGroup =
      as.character(
        .data$RemovedResponseGroup
      ),
    Full_AUC =
      as.numeric(.data$Full_AUC),
    LOSO_AUC =
      as.numeric(
        .data$LeaveOneOut_AUC
      ),
    Delta_AUC =
      as.numeric(
        .data$Delta_AUC_minus_full
      )
  )

assert_true(
  nrow(panelB) == 27L,
  paste0(
    "Panel B target-state LOSO must contain 27 leave-one-out values; observed ",
    nrow(panelB),
    "."
  )
)

assert_true(
  length(
    unique(panelB$RemovedSample)
  ) == 27L,
  "Panel B RemovedSample IDs are not unique."
)

assert_true(
  length(
    unique(panelB$Full_AUC)
  ) == 1L,
  "Panel B Full_AUC is not constant across LOSO rows."
)

panelB$DisplayGroup <-
  "tumor-dedifferentiation/\nstromal-remodeling"

pB <- ggplot2::ggplot(
  panelB,
  ggplot2::aes(
    x = .data$DisplayGroup,
    y = .data$LOSO_AUC
  )
) +
  ggplot2::geom_hline(
    yintercept =
      unique(panelB$Full_AUC),
    linewidth = 0.45,
    linetype = "solid"
  ) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  ggplot2::geom_boxplot(
    width = 0.28,
    outlier.shape = NA,
    color = TARGET_COLOR
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(
      shape =
        .data$RemovedResponseGroup
    ),
    width = 0.06,
    height = 0,
    size = 1.9,
    alpha = 0.85,
    color = TARGET_COLOR
  ) +
  ggplot2::scale_shape_discrete(
    name = "Removed sample"
  ) +
  ggplot2::labs(
    title =
      "Leave-one-sample-out stability",
    subtitle =
      paste0(
        "Solid line = full-sample AUC ",
        sprintf(
          "%.3f",
          unique(panelB$Full_AUC)
        )
      ),
    x = NULL,
    y = "LOSO AUC"
  ) +
  publication_theme()

############################################################
## 9. Panel C — GSE78220 composition-adjusted target state
############################################################

panelC <- s13_auc %>%
  dplyr::filter(
    .data$cohort == "GSE78220",
    .data$state_internal ==
      TARGET_STATE,
    .data$adjustment %in%
      S13_METHOD_ORDER
  ) %>%
  dplyr::transmute(
    Cohort =
      as.character(.data$cohort),
    State =
      as.character(
        .data$state_internal
      ),
    StateDisplay = TARGET_DISPLAY,
    Adjustment =
      as.character(
        .data$adjustment
      ),
    MethodDisplay =
      unname(
        S13_METHOD_DISPLAY[
          as.character(
            .data$adjustment
          )
        ]
      ),
    AUC =
      as.numeric(.data$auc),
    CI_low =
      as.numeric(
        .data$bootstrap_ci_low
      ),
    CI_high =
      as.numeric(
        .data$bootstrap_ci_high
      ),
    LOSO_median =
      as.numeric(
        .data$loso_auc_median
      ),
    LOSO_min =
      as.numeric(
        .data$loso_auc_min
      ),
    LOSO_max =
      as.numeric(
        .data$loso_auc_max
      )
  )

assert_true(
  nrow(panelC) ==
    length(S13_METHOD_ORDER),
  paste0(
    "Panel C must contain exactly ",
    length(S13_METHOD_ORDER),
    " prespecified S13 methods; observed ",
    nrow(panelC),
    "."
  )
)

assert_true(
  setequal(
    panelC$Adjustment,
    S13_METHOD_ORDER
  ),
  "Panel C is missing one or more prespecified S13 methods."
)

panelC$MethodDisplay <- factor(
  panelC$MethodDisplay,
  levels = rev(
    unname(
      S13_METHOD_DISPLAY[
        S13_METHOD_ORDER
      ]
    )
  )
)

pC <- ggplot2::ggplot(
  panelC,
  ggplot2::aes(
    x = .data$MethodDisplay,
    y = .data$AUC
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = .data$CI_low,
      ymax = .data$CI_high
    ),
    width = 0.18,
    linewidth = 0.45,
    color = TARGET_COLOR
  ) +
  ggplot2::geom_point(
    size = 2.8,
    color = TARGET_COLOR
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title =
      "Composition-aware sensitivity",
    subtitle =
      "GSE78220 tumor-dedifferentiation/stromal-remodeling state",
    x = NULL,
    y = "AUC for non-response"
  ) +
  publication_theme()

############################################################
## 10. Panel D — canonical cross-cohort harmonization
############################################################

panelD <- s19_auc %>%
  dplyr::filter(
    .data$State %in% STATE_COLS
  ) %>%
  dplyr::left_join(
    s19_perm %>%
      dplyr::transmute(
        State =
          as.character(.data$State),
        Permutation_P =
          as.numeric(
            .data$P_two_sided_from_0_5
          ),
        Permutation_N =
          as.integer(.data$N_perm),
        N_extreme_exact =
          as.integer(
            .data$N_extreme_exact
          )
      ),
    by = "State"
  ) %>%
  dplyr::transmute(
    State =
      as.character(.data$State),
    StateDisplay =
      unname(
        STATE_DISPLAY[
          as.character(.data$State)
        ]
      ),
    StateDisplayPlot =
      unname(
        STATE_DISPLAY_PLOT[
          as.character(.data$State)
        ]
      ),
    N =
      as.integer(.data$N),
    Responders =
      as.integer(.data$Responders),
    NonResponders =
      as.integer(.data$NonResponders),
    AUC =
      as.numeric(.data$AUC),
    CI_low =
      as.numeric(.data$CI_low),
    CI_high =
      as.numeric(.data$CI_high),
    Direction =
      as.character(.data$Direction),
    Permutation_P =
      .data$Permutation_P,
    Permutation_N =
      .data$Permutation_N,
    N_extreme_exact =
      .data$N_extreme_exact
  )

assert_true(
  nrow(panelD) == 4L,
  "Panel D must contain exactly four canonical states."
)

assert_true(
  setequal(
    panelD$State,
    STATE_COLS
  ),
  "Panel D state universe differs from the canonical four states."
)

assert_true(
  all(panelD$N == 60L),
  "Panel D cross-cohort sample size must be n=60 for every state."
)

assert_true(
  all(
    panelD$Permutation_N ==
      5000L
  ),
  "Panel D exact permutation N must be 5,000 for every state."
)

panelD$StateDisplay <- factor(
  panelD$StateDisplay,
  levels = rev(unname(STATE_DISPLAY[STATE_COLS]))
)

panelD$StateDisplayPlot <- factor(
  panelD$StateDisplayPlot,
  levels = rev(unname(STATE_DISPLAY_PLOT[STATE_COLS]))
)

pD <- ggplot2::ggplot(
  panelD,
  ggplot2::aes(
    x = .data$StateDisplayPlot,
    y = .data$AUC
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = .data$CI_low,
      ymax = .data$CI_high,
      color = .data$State
    ),
    width = 0.18,
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = .data$State),
    size = 3.2
  ) +
  ggplot2::scale_color_manual(values = STATE_COLORS, guide = "none") +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.3f",
        .data$AUC
      )
    ),
    hjust = -0.32,
    size = 2.8
  ) +
  ggplot2::coord_flip(
    ylim = c(
      min(0.25, min(panelD$CI_low) - 0.03),
      min(1.0, max(panelD$CI_high) + 0.07)
    )
  ) +
  ggplot2::labs(
    title =
      "Cross-cohort harmonized associations",
    subtitle =
      "GSE78220 n=27 + GSE91061 n=33; cohort-only ComBat; DeLong 95% CI",
    x = NULL,
    y = "Fixed-direction AUC for non-response"
  ) +
  publication_theme()

############################################################
## 11. Panel E — direct target-state attenuation
##
## IMPORTANT:
##   No averaging across heterogeneous composition-adjustment
##   models and no averaging of cross-cohort state AUCs.
############################################################

raw_target <- panelA %>%
  dplyr::filter(
    .data$State ==
      TARGET_STATE
  )

harm_target <- panelD %>%
  dplyr::filter(
    .data$State ==
      TARGET_STATE
  )

assert_true(
  nrow(raw_target) == 1L,
  "Panel E raw target-state row is not unique."
)

assert_true(
  nrow(harm_target) == 1L,
  "Panel E harmonized target-state row is not unique."
)

panelE <- data.frame(
  Stage = c(
    "GSE78220 raw",
    "Cross-cohort harmonized"
  ),
  AUC = c(
    raw_target$AUC[[1]],
    harm_target$AUC[[1]]
  ),
  N = c(
    raw_target$N[[1]],
    harm_target$N[[1]]
  ),
  stringsAsFactors = FALSE
)

panelE$Stage <- factor(
  panelE$Stage,
  levels = panelE$Stage
)

pE <- ggplot2::ggplot(
  panelE,
  ggplot2::aes(
    x = .data$Stage,
    y = .data$AUC,
    group = 1
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_line(
    linewidth = 0.65,
    color = TARGET_COLOR
  ) +
  ggplot2::geom_point(
    size = 3.4,
    color = TARGET_COLOR
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%.3f",
        .data$AUC
      )
    ),
    vjust = -0.9,
    size = 3.1
  ) +
  ggplot2::labs(
    title =
      "Attenuation of the target-state association after harmonization",
    subtitle =
      "Direct point-estimate comparison; no averaging across heterogeneous models",
    x = NULL,
    y = "AUC for non-response"
  ) +
  publication_theme()

############################################################
## 12. Cross-panel reporting gates
############################################################

target_raw_auc <-
  raw_target$AUC[[1]]

target_harmonized_auc <-
  harm_target$AUC[[1]]

reporting_gates <- dplyr::bind_rows(
  data.frame(
    Gate = "Script05_all_upstream_gates_PASS",
    Observed = as.character(sum(parse_pass(s05_gates$Pass))),
    Expected = as.character(nrow(s05_gates)),
    Pass = all(parse_pass(s05_gates$Pass)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "Cross_cohort_all_upstream_gates_PASS",
    Observed = as.character(sum(parse_pass(s19_gates$Pass))),
    Expected = as.character(nrow(s19_gates)),
    Pass = all(parse_pass(s19_gates$Pass)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "PanelA_four_states",
    Observed = as.character(nrow(panelA)),
    Expected = "4",
    Pass = nrow(panelA) == 4L,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "PanelB_target_LOSO_n",
    Observed = as.character(nrow(panelB)),
    Expected = "27",
    Pass = nrow(panelB) == 27L,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "PanelC_six_prespecified_methods",
    Observed = as.character(nrow(panelC)),
    Expected = as.character(length(S13_METHOD_ORDER)),
    Pass = nrow(panelC) == length(S13_METHOD_ORDER),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "PanelD_four_states_n60",
    Observed = paste0(nrow(panelD), " states; N=", paste(sort(unique(panelD$N)), collapse = ";")),
    Expected = "4 states; N=60",
    Pass = nrow(panelD) == 4L && all(panelD$N == 60L),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "PanelD_target_canonical_AUC",
    Observed = as.character(target_harmonized_auc),
    Expected = "0.6057142857142858",
    Pass = abs(target_harmonized_auc - 0.6057142857142858) <= 1e-12,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "PanelD_target_exact_permutation_P",
    Observed = as.character(harm_target$Permutation_P[[1]]),
    Expected = "0.1569686062787443",
    Pass = abs(harm_target$Permutation_P[[1]] - 0.1569686062787443) <= 1e-12,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "PanelE_no_heterogeneous_AUC_average",
    Observed = "raw target point + harmonized target point",
    Expected = "raw target point + harmonized target point",
    Pass = TRUE,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gate = "Historical_18H_not_used",
    Observed = "No legacy 18H input path",
    Expected = "No legacy 18H input path",
    Pass = TRUE,
    stringsAsFactors = FALSE
  )
)

if (any(!reporting_gates$Pass)) {
  stop(
    "S30 reporting gate failure:\n",
    paste(
      reporting_gates$Gate[
        !reporting_gates$Pass
      ],
      collapse = "\n"
    )
  )
}

############################################################
## 13. Assemble final figure
############################################################

fig_s30 <-
  (
    pA |
      pB
  ) /
  (
    pC |
      pD
  ) /
  pE +
  patchwork::plot_layout(
    heights = c(
      1,
      1.15,
      0.85
    )
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title =
      paste0(
        "Supplementary Figure S30. Integrated fragility, composition-aware sensitivity, ",
        "and cross-cohort robustness of tumor–immune state associations"
      ),
    subtitle =
      paste0(
        "Integrated reporting of prespecified sensitivity and cross-cohort comparability analyses; ",
        "descriptive robustness analysis, not biomarker validation"
      ),
    theme =
      ggplot2::theme(
        plot.title =
          ggplot2::element_text(
            face = "bold",
            size = 15,
            hjust = 0,
            family = "sans"
          ),
        plot.subtitle =
          ggplot2::element_text(
            size = 11.25,
            hjust = 0,
            family = "sans"
          ),
        plot.tag =
          ggplot2::element_text(
            face = "bold",
            size = 16,
            family = "sans"
          )
      )
  )

print(fig_s30)

############################################################
## 14. Save figure
############################################################

figure_base <- file.path(
  figure_dir,
  "Supplementary_Figure_S30_integrated_fragility_crosscohort"
)

ggplot2::ggsave(
  filename =
    paste0(
      figure_base,
      ".pdf"
    ),
  plot = fig_s30,
  width = 12.5,
  height = 14.5,
  units = "in"
)

ggplot2::ggsave(
  filename =
    paste0(
      figure_base,
      ".png"
    ),
  plot = fig_s30,
  width = 12.5,
  height = 14.5,
  units = "in",
  dpi = 300
)

ggplot2::ggsave(
  filename =
    paste0(
      figure_base,
      ".jpg"
    ),
  plot = fig_s30,
  width = 12.5,
  height = 14.5,
  units = "in",
  dpi = 300
)

############################################################
## 15. Export source-data workbook
############################################################

readme_df <- data.frame(
  Item = c(
    "Figure",
    "Role",
    "Panel A source",
    "Panel B source",
    "Panel C source",
    "Panel D source",
    "Panel D permutation source",
    "Panel E definition",
    "Historical 18H used",
    "Interpretation"
  ),
  Value = c(
    "Supplementary Figure S30",
    "Integrated robustness / sensitivity reporting",
    basename(s11_auc_file),
    basename(s11_loso_file),
    basename(s13_auc_file),
    basename(s19_auc_file),
    basename(s19_perm_file),
    paste0(
      "Direct ",
      TARGET_DISPLAY,
      " point-estimate comparison: GSE78220 raw vs cross-cohort harmonized"
    ),
    "No",
    "Descriptive robustness/sensitivity; not predictive model or biomarker validation"
  ),
  stringsAsFactors = FALSE
)

wb <- openxlsx::createWorkbook()

sheet_list <- list(
  README = readme_df,
  PanelA_raw_AUC = panelA,
  PanelB_LOSO = panelB,
  PanelC_composition = panelC,
  PanelD_cross_cohort = panelD,
  PanelE_attenuation = panelE,
  Script19_bootstrap = s19_boot,
  Script19_permutation = s19_perm,
  Input_manifest = input_manifest,
  Reporting_gates = reporting_gates
)

header_style <- openxlsx::createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE
)

for (sheet_name in names(sheet_list)) {
  dat <- sheet_list[[sheet_name]]

  openxlsx::addWorksheet(
    wb,
    sheetName = sheet_name,
    gridLines = FALSE
  )

  openxlsx::writeData(
    wb,
    sheet = sheet_name,
    x = dat,
    headerStyle = header_style
  )

  openxlsx::freezePane(
    wb,
    sheet = sheet_name,
    firstActiveRow = 2
  )

  openxlsx::setColWidths(
    wb,
    sheet = sheet_name,
    cols = seq_len(ncol(dat)),
    widths = "auto"
  )
}

source_xlsx <- file.path(
  table_dir,
  "Supplementary_Figure_S30_source_data.xlsx"
)

openxlsx::saveWorkbook(
  wb,
  source_xlsx,
  overwrite = TRUE
)

############################################################
## 16. Export compact CSV sources / audit
############################################################

write_csv_safe(
  panelA,
  file.path(
    table_dir,
    "S30_PanelA_raw_AUC.csv"
  )
)

write_csv_safe(
  panelB,
  file.path(
    table_dir,
    "S30_PanelB_LOSO.csv"
  )
)

write_csv_safe(
  panelC,
  file.path(
    table_dir,
    "S30_PanelC_composition_adjusted.csv"
  )
)

write_csv_safe(
  panelD,
  file.path(
    table_dir,
    "S30_PanelD_cross_cohort.csv"
  )
)

write_csv_safe(
  panelE,
  file.path(
    table_dir,
    "S30_PanelE_attenuation.csv"
  )
)

write_csv_safe(
  input_manifest,
  file.path(
    table_dir,
    "S30_input_manifest.csv"
  )
)

write_csv_safe(
  reporting_gates,
  file.path(
    table_dir,
    "S30_reporting_gates.csv"
  )
)

############################################################
## 17. Session information
############################################################

session_file <- file.path(
  log_dir,
  "sessionInfo_15_S30_cross_cohort_robustness.txt"
)

sink(session_file)
print(sessionInfo())
sink()

############################################################
## 18. Output MD5 manifest
############################################################

manifest_file <- file.path(
  table_dir,
  "S30_output_MD5_manifest.csv"
)

output_files <- c(
  list.files(
    figure_dir,
    full.names = TRUE
  ),
  list.files(
    table_dir,
    full.names = TRUE
  ),
  session_file
)

## Exclude manifest from its own hash inventory.
output_files <- output_files[
  normalizePath(
    output_files,
    winslash = "/",
    mustWork = FALSE
  ) !=
    normalizePath(
      manifest_file,
      winslash = "/",
      mustWork = FALSE
    )
]

output_manifest <- data.frame(
  File = basename(output_files),
  Path = output_files,
  SizeBytes = file.info(output_files)$size,
  MD5 = vapply(
    output_files,
    md5_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_safe(
  output_manifest,
  manifest_file
)

############################################################
## 19. Console summary
############################################################

cat("\n============================================================\n")
cat("SUPPLEMENTARY FIGURE S30 REPORTING REBUILD: PASS\n")
cat("============================================================\n")
cat("Panel A states: ", nrow(panelA), "\n", sep = "")
cat("Panel B target-state LOSO values: ", nrow(panelB), "\n", sep = "")
cat("Panel C composition methods: ", nrow(panelC), "\n", sep = "")
cat("Panel D cross-cohort states: ", nrow(panelD), "\n", sep = "")
cat(
  "Target raw GSE78220 AUC: ",
  sprintf("%.6f", target_raw_auc),
  "\n",
  sep = ""
)
cat(
  "Target harmonized AUC: ",
  sprintf("%.6f", target_harmonized_auc),
  "\n",
  sep = ""
)
cat(
  "Target exact within-cohort permutation P: ",
  sprintf(
    "%.6f",
    harm_target$Permutation_P[[1]]
  ),
  "\n",
  sep = ""
)
cat("Historical 18H inferential table used: NO\n")
cat("All S30 reporting gates: PASS\n")
cat(
  "Output root: ",
  normalizePath(
    out_root,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")
