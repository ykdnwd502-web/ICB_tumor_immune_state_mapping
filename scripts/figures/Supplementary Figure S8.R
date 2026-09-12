# =============================================================================
# Supplementary_Figure_S8_build_final_v2.1.4.2_FINAL_visual_polish.R
#
# Reporting-only builder for Supplementary Figure S8
#
# Frozen architecture:
#   canonical CSV sources -> plotting layer -> final figure
#
# No Seurat object.
# No RDS.
# No score recalculation.
#
# v2.1.4.2 visual polish:
#   - Panel B compartment strip labels wrapped to two lines where needed
#   - Panel D long Y-axis label wrapped to two lines
#   - Panel D title kept on one line
#   - Layout widened on right column to avoid clipping
#   - Panel D whitespace reduced by shorter labels and tighter margins
#   - JPG and PNG export at 300 dpi
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) {
  project_dir <- "D:/ICB_resistance_project"
}

source_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "GSE244983",
  "state_localization_source_attribution"
)

figure_dir <- file.path(
  project_dir,
  "results",
  "figures",
  "all_supplementary_figures"
)

audit_dir <- file.path(
  project_dir,
  "results",
  "audit",
  "Supplementary_Figure_S8_v2.1.4.2_visual_polish"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

read_source <- function(filename) {
  read.csv(
    file.path(source_dir, filename),
    check.names = FALSE
  )
}

pick_col <- function(df, candidates, what) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop(
      "Cannot identify column for ", what,
      ". Available columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  hit[1]
}

# =============================================================================
# Read frozen canonical sources
# =============================================================================

s8a_raw <- read_source("GSE244983_source_attribution_state_means.csv")
s8b_raw <- read_source("GSE244983_source_attribution_prespecified_contrasts.csv")
s8c_raw <- read_source("GSE244983_canonical_marker_module_by_celltype.csv")
s8d_raw <- read_source("GSE244983_state_marker_module_spearman.csv")

# =============================================================================
# Canonical column mapping
# =============================================================================

s8a <- s8a_raw %>%
  rename(
    MajorCellType = all_of(pick_col(
      s8a_raw,
      c("MajorCellType", "MajorCellType...MajorCellType", "CellType"),
      "S8A MajorCellType"
    )),
    State = all_of(pick_col(
      s8a_raw,
      c("State", "State...State"),
      "S8A State"
    )),
    MeanScore = all_of(pick_col(
      s8a_raw,
      c("MeanScore", "MeanScore...MeanScore", "StateScore"),
      "S8A MeanScore"
    ))
  )

if ("StateDisplay" %in% names(s8a_raw)) {
  s8a$StateDisplay <- s8a_raw[[pick_col(
    s8a_raw,
    c("StateDisplay"),
    "S8A StateDisplay"
  )]]
} else {
  s8a$StateDisplay <- s8a$State
}

s8b <- s8b_raw %>%
  rename(
    ContrastLabel = all_of(pick_col(
      s8b_raw,
      c("ContrastLabel", "Contrast", "Comparison"),
      "S8B ContrastLabel"
    )),
    MedianDifference = all_of(pick_col(
      s8b_raw,
      c("MedianDifference", "Difference", "MeanDifference", "Delta"),
      "S8B MedianDifference"
    )),
    RankBiserial = all_of(pick_col(
      s8b_raw,
      c("RankBiserial", "r_rb", "Rank_Biserial"),
      "S8B RankBiserial"
    ))
  )

if ("StateDisplay" %in% names(s8b_raw)) {
  s8b$StateDisplay <- s8b_raw[["StateDisplay"]]
} else if ("State" %in% names(s8b_raw)) {
  s8b$StateDisplay <- s8b_raw[["State"]]
} else {
  stop("S8B requires StateDisplay or State")
}

s8c <- s8c_raw %>%
  rename(
    MajorCellType = all_of(pick_col(
      s8c_raw,
      c("MajorCellType", "MajorCellType...MajorCellType", "CellType"),
      "S8C MajorCellType"
    )),
    MarkerModule = all_of(pick_col(
      s8c_raw,
      c("MarkerModule", "Module"),
      "S8C MarkerModule"
    )),
    MeanModuleScore = all_of(pick_col(
      s8c_raw,
      c("MeanModuleScore", "MeanScore", "ModuleScore"),
      "S8C MeanModuleScore"
    ))
  )

s8d <- s8d_raw %>%
  rename(
    State = all_of(pick_col(
      s8d_raw,
      c("State", "State...State"),
      "S8D State"
    )),
    MarkerModule = all_of(pick_col(
      s8d_raw,
      c("MarkerModule", "Module"),
      "S8D MarkerModule"
    )),
    SpearmanRho = all_of(pick_col(
      s8d_raw,
      c("SpearmanRho", "rho", "Spearman_rho"),
      "S8D SpearmanRho"
    ))
  )

# =============================================================================
# Controlled display labels
# =============================================================================

state_map <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid-Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation",
  "immune-defective/cold" = "immune-defective/cold",
  "myeloid-Treg immunosuppressive" = "myeloid-Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation" = "melanocytic differentiation"
)

module_map <- c(
  "canonical melanocytic" = "canonical melanocytic",
  "stromal/CAF" = "stromal/CAF",
  "endothelial" = "endothelial",
  "myeloid" = "myeloid",
  "T/NK" = "T/NK",
  "B/plasma" = "B/plasma"
)

celltype_order <- c(
  "B/Plasma cells",
  "CAF/stromal-like cells",
  "Cycling malignant",
  "Endothelial",
  "Malignant",
  "Myeloid cells",
  "T/NK cells",
  "T/NK/Treg-like cells"
)

state_order_A <- c(
  "immune-defective/cold",
  "melanocytic differentiation",
  "myeloid-Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling"
)

state_order_B <- c(
  "myeloid-Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation"
)

state_order_D <- c(
  "tumor-dedifferentiation/stromal-remodeling",
  "myeloid-Treg immunosuppressive",
  "melanocytic differentiation",
  "immune-defective/cold"
)

state_facet_labels <- c(
  "myeloid-Treg immunosuppressive" = "myeloid-Treg\nimmunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "melanocytic differentiation" = "melanocytic differentiation"
)

state_axis_labels_D <- c(
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "myeloid-Treg immunosuppressive" = "myeloid-Treg immunosuppressive",
  "melanocytic differentiation" = "melanocytic differentiation",
  "immune-defective/cold" = "immune-defective/cold"
)

# =============================================================================
# Harmonize data
# =============================================================================

s8a$StateDisplay <- ifelse(
  s8a$StateDisplay %in% names(state_map),
  unname(state_map[s8a$StateDisplay]),
  s8a$StateDisplay
)
s8a$MajorCellType <- factor(s8a$MajorCellType, levels = celltype_order)
s8a$StateDisplay <- factor(s8a$StateDisplay, levels = state_order_A)

s8b$StateDisplay <- ifelse(
  s8b$StateDisplay %in% names(state_map),
  unname(state_map[s8b$StateDisplay]),
  s8b$StateDisplay
)
s8b$StateDisplay <- factor(s8b$StateDisplay, levels = state_order_B)
s8b$ContrastDisplay <- gsub(" vs ", "\nvs ", s8b$ContrastLabel, fixed = TRUE)

s8c$MajorCellType <- factor(s8c$MajorCellType, levels = rev(celltype_order))
s8c$MarkerModule <- factor(
  s8c$MarkerModule,
  levels = c(
    "B/plasma",
    "canonical melanocytic",
    "endothelial",
    "myeloid",
    "stromal/CAF",
    "T/NK"
  )
)

s8d$StateDisplay <- ifelse(
  s8d$State %in% names(state_map),
  unname(state_map[s8d$State]),
  s8d$State
)
s8d$StateDisplay <- factor(s8d$StateDisplay, levels = state_order_D)
s8d$MarkerModule <- factor(
  s8d$MarkerModule,
  levels = c(
    "B/plasma",
    "canonical melanocytic",
    "endothelial",
    "myeloid",
    "stromal/CAF",
    "T/NK"
  )
)

# =============================================================================
# Plots
# =============================================================================

pA <- ggplot(
  s8a,
  aes(x = StateDisplay, y = MajorCellType, fill = MeanScore)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", MeanScore)), size = 3) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0
  ) +
  labs(
    title = "Single-cell source attribution",
    subtitle = "Mean z-standardized predefined state score by cell type",
    x = NULL,
    y = NULL,
    fill = "Mean\nz-score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(face = "bold")
  )

pB <- ggplot(
  s8b,
  aes(x = ContrastDisplay, y = MedianDifference)
) +
  geom_col(fill = "grey45", width = 0.82) +
  geom_text(
    aes(
      label = paste0(
        "\u0394=", sprintf("%.2f", MedianDifference),
        "\nr_rb=", sprintf("%.2f", RankBiserial)
      )
    ),
    vjust = -0.15,
    size = 3.2,
    lineheight = 0.95
  ) +
  facet_wrap(
    ~ StateDisplay,
    nrow = 1,
    scales = "free_x",
    labeller = as_labeller(state_facet_labels)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.30))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Pre-specified cell-type contrasts",
    subtitle = "Positive values indicate higher score in the first listed compartment",
    x = NULL,
    y = "Median difference"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 9.5, lineheight = 0.95),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.5),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(t = 8, r = 18, b = 12, l = 8)
  )

pC <- ggplot(
  s8c,
  aes(x = MarkerModule, y = MajorCellType, fill = MeanModuleScore)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", MeanModuleScore)), size = 3) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0
  ) +
  labs(
    title = "Canonical marker-module assessment",
    subtitle = "Mean z-standardized canonical module score by cell type",
    x = NULL,
    y = NULL,
    fill = "Mean module\nz-score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(face = "bold")
  )

pD <- ggplot(
  s8d,
  aes(x = MarkerModule, y = StateDisplay, fill = SpearmanRho)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", SpearmanRho)), size = 3) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0
  ) +
  scale_y_discrete(labels = state_axis_labels_D) +
  labs(
    title = "Composite state versus canonical marker modules",
    subtitle = "Spearman correlation across single cells",
    x = NULL,
    y = NULL,
    fill = "Spearman\nrho"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9.5, lineheight = 0.95),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t = 8, r = 14, b = 8, l = 2)
  )

fig <- (pA | pB) / (pC | pD) +
  plot_layout(
    widths = c(1.00, 1.15),
    heights = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 18, face = "bold")
  )

# =============================================================================
# Export
# =============================================================================

base_name <- "Supplementary Figure S8. Single-cell source attribution and canonical marker-module assessment of predefined tumor–immune state scores in GSE244983"

jpg_file <- file.path(figure_dir, paste0(base_name, ".jpg"))
png_file <- file.path(figure_dir, paste0(base_name, ".png"))

ggsave(
  filename = jpg_file,
  plot = fig,
  width = 15,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = fig,
  width = 15,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

write.csv(
  data.frame(
    Figure = "Supplementary Figure S8",
    Version = "v2.1.4.2_FINAL_visual_polish",
    JPG = jpg_file,
    PNG = png_file,
    stringsAsFactors = FALSE
  ),
  file.path(audit_dir, "S8_v2.1.4.2_export_manifest.csv"),
  row.names = FALSE
)

cat("Supplementary Figure S8 v2.1.4.2 export completed\n")
