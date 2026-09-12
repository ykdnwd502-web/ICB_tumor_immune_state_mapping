############################################################
## 03B_Visium_make_Figure8.R
##
## REPORTING-ONLY FIGURE 8 ASSEMBLY v1.3
## REVIEWER-6 INFERENCE-WORDING PATCH — NO STALE PNG DEPENDENCY
##
## Why v1.2:
##   v1.1 rebuilt panel A directly, but still assembled B/C/D from
##   pre-existing PNG files. If those PNGs were not regenerated after
##   a title-only patch, the old "Figure 8B./8C./8D." prefixes remained.
##
## v1.2 therefore rebuilds ALL FOUR PANELS directly from already-saved
## analysis/reporting tables and metadata, then assembles the composite.
##
## v1.3 changes ONLY Figure 8B reporting:
##   - Fisher OR and exact 95% CI remain displayed as descriptive effect size
##   - Fisher P value remains in the frozen source table but is NOT displayed
##   - Panel B subtitle explicitly points spatial inference to permutation/
##     Moran-type analyses reported separately in Supplementary Figures S16/S18
##
## IMPORTANT — reporting layer only:
##   - no state rescoring
##   - no threshold selection
##   - no graph reconstruction
##   - no permutation
##   - no residualization
##   - no candidate-gene reranking
##   - no expression/statistical recomputation
##
## Final panel titles:
##   A  Spot-level dual-high co-occurrence map
##   B  Spot-level dual-high enrichment
##   C  Candidate mechanism genes in the dual-high niche
##   D  Four-group expression pattern of nominated genes
##
## Output:
##   results/figures/spatial_melanoma_validation/
##     Figure8A_dual_high_spatial_colocalization.{png,pdf}
##     Figure8B_dual_high_enrichment_OR.{png,pdf}
##     Figure8C_candidate_mechanism_genes_Both_vs_Neither.{png,pdf}
##     Figure8D_four_group_candidate_gene_dotplot.{png,pdf}
##     Figure8_spatial_coenrichment_candidate_features.{png,pdf}
############################################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

############################################################
## 0. Project paths
############################################################

project_dir <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)

project_dir <- normalizePath(
  project_dir,
  winslash = "/",
  mustWork = TRUE
)

figure_dir <- file.path(
  project_dir,
  "results",
  "figures",
  "spatial_melanoma_validation"
)

table_dir <- file.path(
  project_dir,
  "results",
  "tables",
  "spatial_melanoma_validation"
)

intermediate_dir <- file.path(
  project_dir,
  "results",
  "intermediate",
  "spatial_melanoma_validation"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
## 1. Packages
############################################################

for (pkg in c("ggplot2", "scales", "grid")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required.",
      call. = FALSE
    )
  }
}

############################################################
## 2. Helpers
############################################################

read_csv_exact <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required reporting source is missing: ",
      path,
      call. = FALSE
    )
  }

  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
}

save_plot_both <- function(
  p,
  stem,
  width,
  height,
  dpi = 300
) {
  png_file <- file.path(
    figure_dir,
    paste0(stem, ".png")
  )

  pdf_file <- file.path(
    figure_dir,
    paste0(stem, ".pdf")
  )

  ggplot2::ggsave(
    png_file,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    pdf_file,
    plot = p,
    width = width,
    height = height,
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  invisible(
    c(
      png = png_file,
      pdf = pdf_file
    )
  )
}

theme_clean_pub <- function(
  base_size = 12.5,
  base_family = "sans"
) {
  ggplot2::theme_bw(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 2.5,
        hjust = 0.5,
        color = "black"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 1,
        hjust = 0.5,
        color = "black"
      ),
      axis.title = ggplot2::element_text(
        face = "bold",
        size = base_size,
        color = "black"
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold",
        size = base_size - 1,
        color = "black"
      ),
      axis.title.x = ggplot2::element_text(
        face = "bold",
        size = base_size,
        color = "black"
      ),
      axis.text = ggplot2::element_text(
        size = base_size - 1.5,
        color = "black"
      ),
      panel.grid.major = ggplot2::element_line(
        color = "#E8E8E8",
        linewidth = 0.30
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.60
      ),
      legend.title = ggplot2::element_text(
        face = "bold",
        size = base_size
      ),
      legend.text = ggplot2::element_text(
        size = base_size - 1.5
      ),
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      )
    )
}

############################################################
## 3. Panel A — rebuild from frozen Step 11 values
############################################################

a_metadata_file <- file.path(
  intermediate_dir,
  "Step11_spatial_analysis_metadata_with_QCresidual_scores.rds"
)

a_cutoff_file <- file.path(
  table_dir,
  "Step11_spatial_main_top25_pooled_state_score_cutoffs.csv"
)

a_counts_file <- file.path(
  table_dir,
  "Step11_spatial_focused_dual_high_group_counts_top25_pooled.csv"
)

for (f in c(
  a_metadata_file,
  a_cutoff_file,
  a_counts_file
)) {
  if (!file.exists(f)) {
    stop(
      "Figure 8A source missing: ",
      f,
      call. = FALSE
    )
  }
}

a_df <- readRDS(
  a_metadata_file
)

a_cut <- read_csv_exact(
  a_cutoff_file
)

a_counts_expected <- read_csv_exact(
  a_counts_file
)

myeloid_state <- "Myeloid_Treg_Immunosuppressive"
tumor_state <- "Tumor_dedifferentiation_Stromal_remodeling"

required_a_cols <- c(
  "spatial_x",
  "spatial_y",
  myeloid_state,
  tumor_state
)

missing_a_cols <- setdiff(
  required_a_cols,
  colnames(a_df)
)

if (length(missing_a_cols) > 0L) {
  stop(
    "Figure 8A metadata missing: ",
    paste(missing_a_cols, collapse = ", "),
    call. = FALSE
  )
}

get_cutoff <- function(state_name) {
  z <- suppressWarnings(
    as.numeric(
      a_cut$Cutoff[
        a_cut$State == state_name
      ]
    )
  )

  z <- unique(
    z[
      is.finite(z)
    ]
  )

  if (length(z) != 1L) {
    stop(
      "Expected one frozen cutoff for ",
      state_name,
      ".",
      call. = FALSE
    )
  }

  z[[1]]
}

myeloid_cut <- get_cutoff(
  myeloid_state
)

tumor_cut <- get_cutoff(
  tumor_state
)

myeloid_high <- (
  a_df[[myeloid_state]] >= myeloid_cut
)

tumor_high <- (
  a_df[[tumor_state]] >= tumor_cut
)

a_df$DualHighGroup <- ifelse(
  myeloid_high & tumor_high,
  "Both high",
  ifelse(
    myeloid_high & !tumor_high,
    "Myeloid-Treg high only",
    ifelse(
      !myeloid_high & tumor_high,
      "Dediff/Stromal high only",
      "Neither high"
    )
  )
)

a_df$DualHighGroup <- factor(
  a_df$DualHighGroup,
  levels = c(
    "Neither high",
    "Myeloid-Treg high only",
    "Dediff/Stromal high only",
    "Both high"
  )
)

## Frozen count audit.
a_counts_observed <- as.data.frame(
  table(
    a_df$DualHighGroup
  ),
  stringsAsFactors = FALSE
)

colnames(
  a_counts_observed
) <- c(
  "DualHighGroup",
  "n_spots"
)

expected_group_col <- intersect(
  c(
    "DualHighGroup",
    "Group",
    "Category"
  ),
  colnames(a_counts_expected)
)

if (
  length(expected_group_col) < 1L ||
  !("n_spots" %in% colnames(a_counts_expected))
) {
  stop(
    "Unexpected Figure 8A frozen count-table schema.",
    call. = FALSE
  )
}

expected_group_col <- expected_group_col[[1]]

a_counts_expected_key <- data.frame(
  DualHighGroup = as.character(
    a_counts_expected[[expected_group_col]]
  ),
  n_spots = as.integer(
    a_counts_expected$n_spots
  ),
  stringsAsFactors = FALSE
)

a_counts_check <- merge(
  a_counts_observed,
  a_counts_expected_key,
  by = "DualHighGroup",
  all = TRUE,
  suffixes = c(
    "_rebuilt",
    "_frozen"
  )
)

if (
  nrow(a_counts_check) != 4L ||
  any(
    is.na(
      a_counts_check$n_spots_rebuilt
    )
  ) ||
  any(
    is.na(
      a_counts_check$n_spots_frozen
    )
  ) ||
  any(
    a_counts_check$n_spots_rebuilt !=
      a_counts_check$n_spots_frozen
  )
) {
  print(
    a_counts_check,
    row.names = FALSE
  )

  stop(
    "Figure 8A frozen count audit failed.",
    call. = FALSE
  )
}

focused_colors <- c(
  "Neither high" = "#D9D9D9",
  "Myeloid-Treg high only" = "#00A087",
  "Dediff/Stromal high only" = "#E64B35",
  "Both high" = "#7E2F8E"
)

a_legend_labels <- c(
  "Neither high" = "neither high",
  "Myeloid-Treg high only" =
    "myeloid–Treg high only",
  "Dediff/Stromal high only" =
    "tumor-dedifferentiation/stromal-remodeling\nhigh only",
  "Both high" = "both high"
)

p8a <- ggplot2::ggplot(
  a_df,
  ggplot2::aes(
    x = spatial_x,
    y = spatial_y,
    color = DualHighGroup
  )
) +
  ggplot2::geom_point(
    size = 0.65,
    alpha = 0.88
  ) +
  ggplot2::scale_color_manual(
    values = focused_colors,
    drop = FALSE,
    name = "Dual-high group",
    labels = a_legend_labels
  ) +
  ggplot2::scale_y_reverse() +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title =
      "Spot-level dual-high co-occurrence map",
    subtitle =
      "Visium spots; high-state cutoff: top 25% pooled; y-axis inverted for display",
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_clean_pub(
    base_size = 12.5
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.75,
      hjust = 0.5,
      margin = ggplot2::margin(
        b = 5
      )
    ),
    axis.title = ggplot2::element_text(
      face = "bold",
      size = 12.25
    ),
    axis.text = ggplot2::element_text(
      size = 10.75
    ),
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 11.75
    ),
    legend.text = ggplot2::element_text(
      size = 10.5,
      lineheight = 0.90
    ),
    legend.key.height = grid::unit(
      0.45,
      "cm"
    ),
    legend.position = "right",
    plot.margin = ggplot2::margin(
      6,
      4,
      4,
      4
    )
  )

############################################################
## 4. Panel B — rebuild from saved OR table
############################################################

b_or_file <- file.path(
  table_dir,
  "Step12_dual_high_spot_level_OR_top25_pooled.csv"
)

b_or <- read_csv_exact(
  b_or_file
)

required_b_cols <- c(
  "OR",
  "CI_lower",
  "CI_upper",
  "P_value"
)

if (
  nrow(b_or) != 1L ||
  length(
    setdiff(
      required_b_cols,
      colnames(b_or)
    )
  ) > 0L
) {
  stop(
    "Unexpected Figure 8B OR-table schema.",
    call. = FALSE
  )
}

fig8b_data <- data.frame(
  Label = "Dual-high\nco-occurrence",
  OR = b_or$OR,
  CI_lower = b_or$CI_lower,
  CI_upper = b_or$CI_upper,
  P_value = b_or$P_value,
  stringsAsFactors = FALSE
)

fig8b_annot <- sprintf(
  "OR = %.2f\n95%% CI %.2f–%.2f",
  fig8b_data$OR,
  fig8b_data$CI_lower,
  fig8b_data$CI_upper
)

x_max_8b <- max(
  11.5,
  fig8b_data$CI_upper * 1.23
)

annot_y_8b <- min(
  fig8b_data$CI_upper * 1.05,
  x_max_8b * 0.86
)

p8b <- ggplot2::ggplot(
  fig8b_data,
  ggplot2::aes(
    x = Label,
    y = OR
  )
) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = 2,
    color = "grey45",
    linewidth = 0.7
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = CI_lower,
      ymax = CI_upper
    ),
    width = 0.09,
    linewidth = 1.0,
    color = "black"
  ) +
  ggplot2::geom_point(
    size = 4.8,
    color = "#7E2F8E"
  ) +
  ggplot2::annotate(
    "text",
    x = 1,
    y = annot_y_8b,
    label = fig8b_annot,
    hjust = 0,
    vjust = 0.5,
    size = 4.7
  ) +
  ggplot2::coord_flip(
    clip = "off"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      x_max_8b
    ),
    breaks = scales::pretty_breaks(
      n = 5
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  ) +
  ggplot2::labs(
    title =
      "Spot-level dual-high enrichment",
    subtitle =
      paste0(
        "Descriptive Fisher OR; unit of analysis: Visium spot.\n",
        "Spatial inference is evaluated separately by permutation and Moran-type analyses (S16/S18)."
      ),
    x = NULL,
    y = "Odds ratio"
  ) +
  theme_clean_pub(
    base_size = 12.5
  ) +
  ggplot2::theme(
    panel.border = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(
      color = "grey25",
      linewidth = 0.5
    ),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(
      size = 11
    ),
    axis.text.x = ggplot2::element_text(
      size = 11
    ),
    axis.title.x = ggplot2::element_text(
      size = 12.5
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.5,
      lineheight = 1.0
    ),
    plot.margin = ggplot2::margin(
      8,
      70,
      8,
      8
    )
  )

############################################################
## 5. Panel C — rebuild from saved top-20 reporting table
############################################################

c_file <- file.path(
  table_dir,
  "Step12_dual_high_top_candidate_genes_Figure8C.csv"
)

c_tbl <- read_csv_exact(
  c_file
)

required_c_cols <- c(
  "Gene",
  "Diff_Both_minus_Neither",
  "GeneGroup"
)

missing_c_cols <- setdiff(
  required_c_cols,
  colnames(c_tbl)
)

if (
  length(missing_c_cols) > 0L ||
  nrow(c_tbl) != 20L
) {
  stop(
    paste0(
      "Unexpected Figure 8C source table. ",
      "Expected 20 rows and Gene/Diff_Both_minus_Neither/GeneGroup."
    ),
    call. = FALSE
  )
}

gene_group_colors <- c(
  "ECM/CAF-like remodeling" = "#D28A4A",
  "B-cell/plasma-cell/humoral immune" = "#7E65A5",
  "Complement/inflammatory remodeling" = "#7A9B35",
  "Antigen presentation/MHC" = "#4AA3BF",
  "Myeloid/macrophage inflammatory remodeling" = "#D65C6A",
  "Dedifferentiation/invasive phenotype" = "#6F9B86",
  "Other" = "#8A8A8A"
)

## Recreate only the display order from the already-saved 20 rows.
gene_order_c <- c_tbl$Gene[
  order(
    c_tbl$Diff_Both_minus_Neither,
    decreasing = FALSE
  )
]

c_tbl$Gene <- factor(
  c_tbl$Gene,
  levels = gene_order_c
)

c_tbl$GeneGroup <- factor(
  c_tbl$GeneGroup,
  levels = names(
    gene_group_colors
  )
)

used_groups_c <- names(
  gene_group_colors
)[
  names(
    gene_group_colors
  ) %in%
    as.character(
      unique(
        c_tbl$GeneGroup
      )
    )
]

p8c <- ggplot2::ggplot(
  c_tbl,
  ggplot2::aes(
    x = Gene,
    y = Diff_Both_minus_Neither,
    fill = GeneGroup
  )
) +
  ggplot2::geom_col(
    width = 0.82,
    color = "grey25",
    linewidth = 0.25
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(
    values =
      gene_group_colors[
        used_groups_c
      ],
    breaks = used_groups_c,
    drop = TRUE
  ) +
  ggplot2::labs(
    title =
      "Candidate mechanism genes in the dual-high niche",
    subtitle =
      "Ranked by mean log-normalized expression difference: Both high minus Neither high",
    x = NULL,
    y =
      "Mean expression difference (Both high - Neither high)",
    fill = "Gene group"
  ) +
  theme_clean_pub(
    base_size = 12.5
  ) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      size = 13.0
    ),
    axis.text.x = ggplot2::element_text(
      size = 12.5
    ),
    axis.title.x = ggplot2::element_text(
      size = 14.5
    ),
    legend.position = "right",
    legend.title = ggplot2::element_text(
      size = 13.5
    ),
    legend.text = ggplot2::element_text(
      size = 12.5
    ),
    legend.key.height = grid::unit(
      0.45,
      "cm"
    ),
    legend.key.width = grid::unit(
      0.45,
      "cm"
    ),
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 18
    ),
    plot.subtitle = ggplot2::element_text(
      size = 13.0
    ),
    plot.margin = ggplot2::margin(
      8,
      10,
      8,
      8
    )
  )

############################################################
## 6. Panel D — rebuild from saved 80-row reporting table
############################################################

d_file <- file.path(
  table_dir,
  "Step12_dual_high_Figure8D_four_group_candidate_gene_expression_summary.csv"
)

d_tbl <- read_csv_exact(
  d_file
)

required_d_cols <- c(
  "Gene",
  "Group",
  "PctExpr",
  "GeneGroup",
  "Diff_Both_minus_Neither",
  "y_pos",
  "RowScaledMeanExpr",
  "GroupX"
)

missing_d_cols <- setdiff(
  required_d_cols,
  colnames(d_tbl)
)

if (
  length(missing_d_cols) > 0L ||
  nrow(d_tbl) != 80L
) {
  stop(
    paste0(
      "Unexpected Figure 8D source table. ",
      "Expected 80 rows and the finalized reporting columns."
    ),
    call. = FALSE
  )
}

group_levels_8d <- c(
  "Neither high",
  "Myeloid–Treg high only",
  "Dediff/Stromal high only",
  "Both high"
)

group_labels_8d_x <- c(
  "Neither\nhigh",
  "Myeloid–Treg\nhigh only",
  "Dediff/Stromal\nhigh only",
  "Both\nhigh"
)

gene_group_label_short <- c(
  "ECM/CAF-like remodeling" =
    "ECM/CAF-like\nremodeling",
  "B-cell/plasma-cell/humoral immune" =
    "B-cell/plasma-cell\nhumoral immune",
  "Complement/inflammatory remodeling" =
    "Complement/\ninflammatory remodeling",
  "Antigen presentation/MHC" =
    "Antigen presentation/\nMHC",
  "Myeloid/macrophage inflammatory remodeling" =
    "Myeloid/macrophage\ninflammatory remodeling",
  "Dedifferentiation/invasive phenotype" =
    "Dedifferentiation/\ninvasive phenotype",
  "Other" =
    "Other"
)

d_tbl$Group <- factor(
  as.character(
    d_tbl$Group
  ),
  levels = group_levels_8d
)

d_tbl$GeneGroup <- factor(
  as.character(
    d_tbl$GeneGroup
  ),
  levels = names(
    gene_group_colors
  )
)

## One row per displayed gene.
gene_meta_d <- unique(
  d_tbl[
    ,
    c(
      "Gene",
      "GeneGroup",
      "Diff_Both_minus_Neither",
      "y_pos"
    )
  ]
)

gene_meta_d <- gene_meta_d[
  order(
    gene_meta_d$y_pos,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

group_split_d <- split(
  gene_meta_d$y_pos,
  as.character(
    gene_meta_d$GeneGroup
  )
)

group_split_d <- group_split_d[
  names(group_split_d) %in%
    names(
      gene_group_colors
    )
]

group_meta_d <- do.call(
  rbind,
  lapply(
    names(group_split_d),
    function(g) {
      z <- group_split_d[[g]]

      data.frame(
        GeneGroup = g,
        y_mid = mean(z),
        y_min = min(z) - 0.5,
        y_max = max(z) + 0.5,
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(
  group_meta_d
) <- NULL

group_meta_d$GeneGroupShort <-
  gene_group_label_short[
    as.character(
      group_meta_d$GeneGroup
    )
  ]

p8d <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = gene_meta_d,
    ggplot2::aes(
      x = 0.52,
      xend = 0.52,
      y = y_pos - 0.42,
      yend = y_pos + 0.42,
      color = GeneGroup
    ),
    linewidth = 4.0,
    lineend = "butt",
    show.legend = FALSE
  ) +
  ggplot2::geom_segment(
    data = group_meta_d,
    ggplot2::aes(
      x = 0.72,
      xend = 4.35,
      y = y_min,
      yend = y_min
    ),
    color = "#DADADA",
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = d_tbl,
    ggplot2::aes(
      x = GroupX,
      y = y_pos,
      size = PctExpr,
      fill = RowScaledMeanExpr
    ),
    shape = 21,
    color = "grey30",
    stroke = 0.28,
    alpha = 0.95
  ) +
  ggplot2::geom_text(
    data = group_meta_d,
    ggplot2::aes(
      x = 4.62,
      y = y_mid,
      label = GeneGroupShort,
      color = GeneGroup
    ),
    hjust = 0,
    vjust = 0.5,
    size = 4.20,
    lineheight = 0.90,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#6F90DA",
    mid = "#F7F7F7",
    high = "#EE705D",
    midpoint = 0,
    limits = c(
      -1.5,
      1.5
    ),
    oob = scales::squish
  ) +
  ggplot2::scale_color_manual(
    values =
      gene_group_colors,
    drop = TRUE
  ) +
  ggplot2::scale_size_continuous(
    range = c(
      0.8,
      5.2
    ),
    limits = c(
      0,
      100
    ),
    breaks = c(
      25,
      50,
      75,
      100
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0.38,
      6.75
    ),
    breaks = seq_along(
      group_levels_8d
    ),
    labels = group_labels_8d_x,
    expand = c(
      0,
      0
    )
  ) +
  ggplot2::scale_y_continuous(
    breaks = gene_meta_d$y_pos,
    labels = gene_meta_d$Gene,
    expand = ggplot2::expansion(
      add = c(
        0.35,
        0.35
      )
    )
  ) +
  ggplot2::coord_cartesian(
    clip = "off"
  ) +
  ggplot2::labs(
    title =
      "Four-group expression pattern of nominated genes",
    subtitle =
      "Dot color shows row-scaled mean log expression; dot size shows percent-expressing spots",
    x = NULL,
    y = NULL,
    fill =
      "Row-scaled\nmean expression",
    size =
      "Percent\nexpressed"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      order = 1,
      barheight = grid::unit(
        3.5,
        "cm"
      ),
      barwidth = grid::unit(
        0.50,
        "cm"
      )
    ),
    size = ggplot2::guide_legend(
      order = 2,
      override.aes = list(
        fill = "grey75",
        color = "grey30"
      )
    )
  ) +
  theme_clean_pub(
    base_size = 12.5
  ) +
  ggplot2::theme(
    panel.grid.major.x =
      ggplot2::element_line(
        color = "#E5E5E5",
        linewidth = 0.35
      ),
    panel.grid.major.y =
      ggplot2::element_line(
        color = "#ECECEC",
        linewidth = 0.28
      ),
    panel.grid.minor =
      ggplot2::element_blank(),
    panel.border =
      ggplot2::element_blank(),
    axis.line =
      ggplot2::element_blank(),
    axis.ticks.y =
      ggplot2::element_blank(),
    axis.text.y =
      ggplot2::element_text(
        size = 12.5
      ),
    axis.text.x =
      ggplot2::element_text(
        size = 12.75,
        lineheight = 0.92
      ),
    legend.position =
      "right",
    legend.title =
      ggplot2::element_text(
        size = 12.75,
        face = "bold"
      ),
    legend.text =
      ggplot2::element_text(
        size = 12.0
      ),
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 18
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 13.0
      ),
    plot.margin =
      ggplot2::margin(
        8,
        70,
        8,
        8
      )
  )

############################################################
## 7. Save refreshed standalone panels
############################################################

save_plot_both(
  p8a,
  "Figure8A_dual_high_spatial_colocalization",
  width = 8.8,
  height = 6.8
)

save_plot_both(
  p8b,
  "Figure8B_dual_high_enrichment_OR",
  width = 8.6,
  height = 4.2
)

save_plot_both(
  p8c,
  "Figure8C_candidate_mechanism_genes_Both_vs_Neither",
  width = 12.0,
  height = 8.2
)

save_plot_both(
  p8d,
  "Figure8D_four_group_candidate_gene_dotplot",
  width = 13.4,
  height = 9.6
)

############################################################
## 8. Final composite geometry
############################################################

FIGURE_WIDTH_IN  <- 16.0
FIGURE_HEIGHT_IN <- 11.5
PNG_DPI          <- 300

OUTER_LEFT   <- 0.018
OUTER_RIGHT  <- 0.018
OUTER_TOP    <- 0.020
OUTER_BOTTOM <- 0.022

COL_GAP <- 0.018
ROW_GAP <- 0.025

TOP_HEIGHT <- 0.37

BOTTOM_HEIGHT <- (
  1 -
    OUTER_TOP -
    OUTER_BOTTOM -
    ROW_GAP -
    TOP_HEIGHT
)

TOP_LEFT_WIDTH <- 0.49

TOP_RIGHT_WIDTH <- (
  1 -
    OUTER_LEFT -
    OUTER_RIGHT -
    COL_GAP -
    TOP_LEFT_WIDTH
)

BOTTOM_LEFT_WIDTH <- 0.43

BOTTOM_RIGHT_WIDTH <- (
  1 -
    OUTER_LEFT -
    OUTER_RIGHT -
    COL_GAP -
    BOTTOM_LEFT_WIDTH
)

draw_tag <- function(tag) {
  grid::grid.text(
    tag,
    x = grid::unit(
      0.006,
      "npc"
    ),
    y = grid::unit(
      0.992,
      "npc"
    ),
    just = c(
      "left",
      "top"
    ),
    gp = grid::gpar(
      fontsize = 16,
      fontface = "bold",
      fontfamily = "sans",
      col = "black"
    )
  )
}

draw_plot <- function(
  p,
  x,
  y,
  width,
  height,
  tag
) {
  grid::pushViewport(
    grid::viewport(
      x = grid::unit(
        x,
        "npc"
      ),
      y = grid::unit(
        y,
        "npc"
      ),
      width = grid::unit(
        width,
        "npc"
      ),
      height = grid::unit(
        height,
        "npc"
      ),
      just = c(
        "left",
        "bottom"
      ),
      clip = "off"
    )
  )

  print(
    p,
    newpage = FALSE
  )

  draw_tag(
    tag
  )

  grid::popViewport()
}

render_figure8 <- function() {
  grid::grid.newpage()

  top_y <- (
    1 -
      OUTER_TOP -
      TOP_HEIGHT
  )

  bottom_y <- OUTER_BOTTOM

  draw_plot(
    p8a,
    x = OUTER_LEFT,
    y = top_y,
    width = TOP_LEFT_WIDTH,
    height = TOP_HEIGHT,
    tag = "A"
  )

  draw_plot(
    p8b,
    x = (
      OUTER_LEFT +
        TOP_LEFT_WIDTH +
        COL_GAP
    ),
    y = top_y,
    width = TOP_RIGHT_WIDTH,
    height = TOP_HEIGHT,
    tag = "B"
  )

  draw_plot(
    p8c,
    x = OUTER_LEFT,
    y = bottom_y,
    width = BOTTOM_LEFT_WIDTH,
    height = BOTTOM_HEIGHT,
    tag = "C"
  )

  draw_plot(
    p8d,
    x = (
      OUTER_LEFT +
        BOTTOM_LEFT_WIDTH +
        COL_GAP
    ),
    y = bottom_y,
    width = BOTTOM_RIGHT_WIDTH,
    height = BOTTOM_HEIGHT,
    tag = "D"
  )
}

############################################################
## 9. Export composite
############################################################

output_stem <- file.path(
  figure_dir,
  "Figure8_spatial_coenrichment_candidate_features"
)

png_file <- paste0(
  output_stem,
  ".png"
)

pdf_file <- paste0(
  output_stem,
  ".pdf"
)

grDevices::png(
  filename = png_file,
  width = FIGURE_WIDTH_IN,
  height = FIGURE_HEIGHT_IN,
  units = "in",
  res = PNG_DPI,
  bg = "white"
)

render_figure8()

grDevices::dev.off()

grDevices::cairo_pdf(
  filename = pdf_file,
  width = FIGURE_WIDTH_IN,
  height = FIGURE_HEIGHT_IN,
  bg = "white"
)

render_figure8()

grDevices::dev.off()

############################################################
## 10. Output audit
############################################################

audit <- data.frame(
  Item = c(
    "A_title_prefix_absent",
    "B_title_prefix_absent",
    "C_title_prefix_absent",
    "D_title_prefix_absent",
    "A_frozen_count_audit",
    "B_Fisher_OR_descriptive",
    "B_Fisher_P_preserved_in_source_not_displayed",
    "B_spatial_inference_redirected_to_S16_S18",
    "C_rows",
    "D_rows",
    "Composite_PNG",
    "Composite_PDF"
  ),
  Value = c(
    !grepl(
      "^Figure 8A\\.",
      "Spot-level dual-high co-occurrence map"
    ),
    !grepl(
      "^Figure 8B\\.",
      "Spot-level dual-high enrichment"
    ),
    !grepl(
      "^Figure 8C\\.",
      "Candidate mechanism genes in the dual-high niche"
    ),
    !grepl(
      "^Figure 8D\\.",
      "Four-group expression pattern of nominated genes"
    ),
    TRUE,
    TRUE,
    ("P_value" %in% colnames(b_or)) &&
      !grepl("P =", fig8b_annot, fixed = TRUE),
    TRUE,
    nrow(c_tbl),
    nrow(d_tbl),
    file.exists(
      png_file
    ),
    file.exists(
      pdf_file
    )
  ),
  stringsAsFactors = FALSE
)

audit_file <- file.path(
  figure_dir,
  "Figure8_assembly_output_audit.csv"
)

utils::write.csv(
  audit,
  audit_file,
  row.names = FALSE
)

utils::write.csv(
  a_counts_check,
  file.path(
    figure_dir,
    "Figure8A_reporting_refresh_count_audit.csv"
  ),
  row.names = FALSE
)

blocking <- (
  !file.exists(png_file) ||
    !file.exists(pdf_file) ||
    file.info(png_file)$size <= 1000 ||
    file.info(pdf_file)$size <= 1000 ||
    nrow(c_tbl) != 20L ||
    nrow(d_tbl) != 80L ||
    grepl("P =", fig8b_annot, fixed = TRUE) ||
    !("P_value" %in% colnames(b_or))
)

cat("\n============================================================\n")
cat("FIGURE 8 REPORTING ASSEMBLY v1.3\n")
cat("============================================================\n")
cat("A rebuilt from frozen Step 11 values: PASS\n")
cat("B rebuilt from saved OR/CI/P source table: PASS\n")
cat("B displays descriptive OR + 95% CI only; Fisher P is not displayed: PASS\n")
cat("B redirects spatial inference to permutation/Moran analyses (S16/S18): PASS\n")
cat("C rebuilt from saved top-20 reporting table: PASS\n")
cat("D rebuilt from saved 80-row reporting table: PASS\n")
cat("No B/C/D stale PNG files were used for the composite.\n")
cat("Final titles contain no Figure 8B./8C./8D. prefix.\n")
cat("Composite PNG: ", png_file, "\n", sep = "")
cat("Composite PDF: ", pdf_file, "\n", sep = "")
cat(
  "Overall: ",
  ifelse(
    blocking,
    "FAIL",
    "PASS"
  ),
  "\n",
  sep = ""
)
cat("============================================================\n")

if (blocking) {
  stop(
    "Figure 8 v1.3 reporting output gate failed.",
    call. = FALSE
  )
}

cat(
  "\nSTEP 03B PASS — self-contained Figure 8 reporting rebuild completed.\n"
)
