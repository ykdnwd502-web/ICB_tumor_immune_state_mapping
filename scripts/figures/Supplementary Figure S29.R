################--------------------------------------------
## S29_patient_level_structure.R
##
## Reporting-only rebuild of Supplementary Figure S29 from
## 15/01 patient-level source tables.
##
## STYLE CONTRACT
##   - sans-serif publication theme
##   - base text 12.5 pt
##   - panel titles 15 pt; panel tags 16 pt
##   - subtitles ~11.25 pt; axis titles 12.5 pt; ticks 11 pt
##   - dense angled labels 10.5 pt
##   - semantic wrapping for long state labels
################--------------------------------------------

options(stringsAsFactors = FALSE)
required_pkgs <- c("readr", "dplyr", "ggplot2", "patchwork", "openxlsx")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR", unset = "D:/ICB_resistance_project")
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
SRC_DIR <- file.path(PROJECT_DIR, "results", "tables", "additional_robustness", "S29_patient_level")
FIG_DIR <- file.path(PROJECT_DIR, "results", "figures", "supplementary")
TAB_DIR <- file.path(PROJECT_DIR, "results", "reporting", "additional_robustness", "S29_patient_level")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

readc <- function(name) {
  readr::read_csv(file.path(SRC_DIR, name), show_col_types = FALSE, progress = FALSE) |>
    as.data.frame(check.names = FALSE)
}

counts <- readc("S29_patient_counts.csv")
means <- readc("S29_patient_state_means.csv")
means_long <- readc("S29_patient_state_means_long.csv")
panelC <- readc("S29_target_state_cell_scores.csv")
gates <- readc("S29_lineage_gates.csv")
manifest <- readc("S29_input_manifest.csv")

PATIENT_ORDER <- c("Pat_ICBnaive2", "Pat_ICBnaive1", "Pat42", "Pat5")
EXPECTED_COUNTS <- c(Pat_ICBnaive2 = 9081L, Pat_ICBnaive1 = 7414L, Pat42 = 7188L, Pat5 = 2370L)
EXPECTED_MEANS <- rbind(
  Pat_ICBnaive2 = c(0.3325608419486311, -0.006194506141643475, 0.1759250596840239, 0.29067310298207455),
  Pat_ICBnaive1 = c(0.39356295782044515, -0.2730601224886434, -0.22499937648125515, 0.16393561070755508),
  Pat42 = c(-0.2653437647613925, -0.19702662352376465, 0.026172049120634644, 0.037924389913327694),
  Pat5 = c(-1.7006623603001734, 1.4755052440050156, -0.04960117250536579, -1.7416124812924185)
)
MEAN_COLS <- c(
  "Immune_defective_Cold_mean",
  "Myeloid_Treg_Immunosuppressive_mean",
  "Tumor_dedifferentiation_Stromal_remodeling_mean",
  "Melanocytic_Differentiation_mean"
)
TOL <- 1e-12
TARGET_COLOR <- "#E64B35"

STATE_DISPLAY_PLOT <- c(
  "immune-defective/cold" = "immune-defective/cold",
  "myeloid–Treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "melanocytic differentiation" = "melanocytic\ndifferentiation"
)

publication_theme <- function() {
  ggplot2::theme_bw(base_size = 12.5, base_family = "sans") +
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
      plot.margin = ggplot2::margin(7, 8, 7, 8)
    )
}

stopifnot(all(as.logical(gates$Pass)))
counts <- counts[match(PATIENT_ORDER, counts$Sample), , drop = FALSE]
stopifnot(identical(as.integer(counts$n_cells), as.integer(EXPECTED_COUNTS[PATIENT_ORDER])))
means <- means[match(PATIENT_ORDER, means$Sample), , drop = FALSE]
max_mean_diff <- max(abs(as.matrix(means[, MEAN_COLS]) - EXPECTED_MEANS[PATIENT_ORDER, ]), na.rm = TRUE)
if (!is.finite(max_mean_diff) || max_mean_diff > TOL) stop("S29 patient-state mean drift: ", max_mean_diff, call. = FALSE)

counts$Sample <- factor(counts$Sample, levels = PATIENT_ORDER)
means_long$Sample <- factor(means_long$Sample, levels = PATIENT_ORDER)
panelC$Sample <- factor(panelC$Sample, levels = PATIENT_ORDER)
state_levels <- names(STATE_DISPLAY_PLOT)
means_long$StateLabel <- factor(means_long$StateLabel, levels = state_levels)
means_long$StateLabelPlot <- factor(
  unname(STATE_DISPLAY_PLOT[as.character(means_long$StateLabel)]),
  levels = unname(STATE_DISPLAY_PLOT[state_levels])
)

pA <- ggplot2::ggplot(counts, ggplot2::aes(x = Sample, y = n_cells)) +
  ggplot2::geom_col(width = 0.70) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(n_cells, "\n(", proportion_pct, ")")),
    vjust = -0.30, size = 3.6, family = "sans"
  ) +
  ggplot2::labs(title = "Per-patient cell contribution", x = NULL, y = "Number of cells") +
  ggplot2::coord_cartesian(ylim = c(0, max(counts$n_cells) * 1.18), clip = "off") +
  publication_theme()

pB <- ggplot2::ggplot(means_long, ggplot2::aes(x = StateLabelPlot, y = Sample, fill = MeanScore)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", MeanScore)), size = 3.1, family = "sans") +
  ggplot2::scale_y_discrete(limits = rev(PATIENT_ORDER), drop = FALSE) +
  ggplot2::scale_fill_gradient2(low = "#3B7DDD", mid = "white", high = "#D84B4B", midpoint = 0) +
  ggplot2::labs(
    title = "Patient-level mean tumor–immune state scores",
    x = NULL, y = NULL,
    fill = "Mean score across cells\nfrom each patient"
  ) +
  publication_theme() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 10.5))

pC <- ggplot2::ggplot(
  panelC,
  ggplot2::aes(x = Sample, y = Tumor_dedifferentiation_Stromal_remodeling_score)
) +
  ggplot2::geom_violin(trim = FALSE, fill = TARGET_COLOR, color = TARGET_COLOR, alpha = 0.22) +
  ggplot2::geom_boxplot(width = 0.15, fill = "white", color = TARGET_COLOR, outlier.size = 0.3) +
  ggplot2::labs(
    title = "Distribution of tumor-dedifferentiation/\nstromal-remodeling state scores across patients",
    x = NULL,
    y = "Tumor-dedifferentiation/stromal-remodeling score"
  ) +
  publication_theme()

fig <- (pA | pB) / pC +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Patient-level structure of predefined tumor–immune states in GSE244983",
    subtitle = "GSE244983 single-cell cohort: n = 26,053 cells from 4 patients",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0, family = "sans"),
      plot.subtitle = ggplot2::element_text(size = 11.25, hjust = 0, family = "sans"),
      plot.tag = ggplot2::element_text(face = "bold", size = 16, family = "sans")
    )
  )

base <- file.path(FIG_DIR, "Supplementary Figure S29. Patient-level structure of predefined tumor–immune states in GSE244983")
for (ext in c("pdf", "png", "jpg")) {
  ggplot2::ggsave(
    paste0(base, ".", ext), fig,
    width = 13, height = 11, dpi = 300,
    limitsize = FALSE, bg = "white"
  )
}

readme <- data.frame(
  Item = c("Figure", "Primary lineage", "Patients", "Role", "Historical support branch", "Max patient-state mean diff"),
  Value = c(
    "Supplementary Figure S29", "26,053 cells", "4",
    "Patient-level structure audit requested during revision",
    "Not used", format(max_mean_diff, scientific = TRUE)
  ),
  stringsAsFactors = FALSE
)
openxlsx::write.xlsx(
  list(
    README = readme,
    Patient_counts = counts,
    State_means = means,
    State_means_long = means_long,
    Tumor_stromal_scores = panelC,
    Lineage_gates = gates,
    Input_manifest = manifest
  ),
  file.path(TAB_DIR, "Supplementary_Figure_S29_source_data.xlsx"),
  overwrite = TRUE,
  asTable = FALSE
)

figure_gate <- data.frame(
  Gate = c(
    "Primary lineage gates PASS", "Four patient counts exact",
    "Patient-state means exact within 1e-12",
    "Figure PDF exists", "Figure PNG exists", "Figure JPG exists"
  ),
  Observed = c(
    all(as.logical(gates$Pass)), sum(counts$n_cells), max_mean_diff,
    file.exists(paste0(base, ".pdf")), file.exists(paste0(base, ".png")), file.exists(paste0(base, ".jpg"))
  ),
  Expected = c(TRUE, 26053, "<=1e-12", TRUE, TRUE, TRUE),
  Pass = c(
    all(as.logical(gates$Pass)), sum(counts$n_cells) == 26053, max_mean_diff <= TOL,
    file.exists(paste0(base, ".pdf")), file.exists(paste0(base, ".png")), file.exists(paste0(base, ".jpg"))
  ),
  stringsAsFactors = FALSE
)
readr::write_csv(figure_gate, file.path(TAB_DIR, "S29_reporting_gates.csv"))
if (any(!figure_gate$Pass)) stop("S29 reporting gate failure.", call. = FALSE)
cat("S29 reporting PASS; max patient-state mean diff: ", format(max_mean_diff, scientific = TRUE), "\n", sep = "")