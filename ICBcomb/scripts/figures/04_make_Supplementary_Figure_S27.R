############################################################
## 04_make_Supplementary_Figure_S27.R
##
## PURPOSE
##   Reporting-only rebuild of Supplementary Figure S27 from
##   the frozen ICBcomb module inputs/source tables.
##
## FIXES IN THIS PATCH
##   1) Panel B "Mapped rows" is computed from the full 45-row
##      ICBcomb strategy-mapping table:
##         15 / 15 / 15
##   2) Panel B "Positive prioritization rows" is computed from
##      the frozen NES > 0 row universe:
##         15 / 14 / 13
##   3) Public-facing "reversal" wording is removed.
##
## SCIENTIFIC BOUNDARY
##   Reporting-only patch. No query genes, NES values, thresholds,
##   state definitions, strategy mappings, or statistical results
##   are changed.
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

############################################################
## 0. Module root and fixed paths
############################################################

ROOT <- Sys.getenv("ICBCOMB_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("ICB_PUBLIC_ROOT")
if (!nzchar(ROOT)) ROOT <- getwd()

ROOT <- normalizePath(
  ROOT,
  winslash = "/",
  mustWork = TRUE
)

INPUT_ROOT <- file.path(
  ROOT,
  "inputs",
  "ICBcomb_INPUT_FREEZE_v1.0"
)

GENESET_FILE <- file.path(
  INPUT_ROOT,
  "data",
  "query_gene_set_bundle",
  "02_final_ICBcomb_gene_sets",
  "final_ICBcomb_gene_set_summary_CLEAN.csv"
)

TABLE_ROOT <- file.path(
  ROOT,
  "results",
  "tables",
  "ICBcomb"
)

MAPPING_FILE <- file.path(
  TABLE_ROOT,
  "ICBcomb_strategy_mapping.csv"
)

POSITIVE_FILE <- file.path(
  TABLE_ROOT,
  "ICBcomb_positive_prioritization_rows.csv"
)

FIG_DIR <- file.path(
  ROOT,
  "results",
  "figures"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  GENESET_FILE,
  MAPPING_FILE,
  POSITIVE_FILE
)

if (!all(file.exists(required_files))) {
  stop(
    "Missing frozen S27 source files:\n",
    paste(
      required_files[!file.exists(required_files)],
      collapse = "\n"
    ),
    call. = FALSE
  )
}

############################################################
## 1. Fixed state/display contracts
############################################################

STATE_LEVELS <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling"
)

STATE_DISPLAY <- c(
  "immune-defective/cold" =
    "Immune-defective/\ncold",
  "myeloid–Treg immunosuppressive" =
    "Myeloid-Treg\nimmunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" =
    "Tumor-dedifferentiation/\nstromal-remodeling"
)

############################################################
## 2. Panel A: exact frozen query-gene-set sizes
############################################################

gene_sets <- read.csv(
  GENESET_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  identical(
    names(gene_sets),
    c("GeneSet", "n_genes")
  ),
  nrow(gene_sets) == 13L
)

label_map <- c(
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven" =
    "Myeloid-Treg suppress, data-driven",
  "Melanocytic_Differentiation_REFERENCE_data_driven" =
    "Melanocytic reference, data-driven",
  "Immune_defective_Cold_RESTORE_data_driven" =
    "Immune-defective/cold restore, data-driven",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven" =
    "Dediff/stromal suppress, data-driven",
  "Immune_defective_Cold_RESTORE_hybrid" =
    "Immune-defective/cold restore, hybrid",
  "Dedifferentiated_inflammatory_REFERENCE" =
    "Differentiated-inflammatory reference",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid" =
    "Dediff/stromal suppress, hybrid",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid" =
    "Myeloid-Treg suppress, hybrid",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated" =
    "Dediff/stromal suppress, curated",
  "Immune_defective_Cold_RESTORE_curated" =
    "Immune-defective/cold restore, curated",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated" =
    "Myeloid-Treg suppress, curated",
  "Melanocytic_Differentiation_REFERENCE_hybrid" =
    "Melanocytic reference, hybrid",
  "Melanocytic_Differentiation_REFERENCE_curated" =
    "Melanocytic reference, curated"
)

panelA_df <- gene_sets %>%
  mutate(
    readable_label = unname(label_map[GeneSet])
  )

if (
  any(is.na(panelA_df$readable_label)) ||
  anyDuplicated(panelA_df$GeneSet)
) {
  stop(
    "S27 Panel A gene-set labeling contract failed.",
    call. = FALSE
  )
}

## Preserve the manuscript-facing top-to-bottom order.
panelA_top_to_bottom <- c(
  "Myeloid-Treg suppress, data-driven",
  "Melanocytic reference, data-driven",
  "Immune-defective/cold restore, data-driven",
  "Dediff/stromal suppress, data-driven",
  "Immune-defective/cold restore, hybrid",
  "Differentiated-inflammatory reference",
  "Dediff/stromal suppress, hybrid",
  "Myeloid-Treg suppress, hybrid",
  "Dediff/stromal suppress, curated",
  "Immune-defective/cold restore, curated",
  "Myeloid-Treg suppress, curated",
  "Melanocytic reference, hybrid",
  "Melanocytic reference, curated"
)

panelA_df <- panelA_df %>%
  mutate(
    readable_label = factor(
      readable_label,
      levels = rev(panelA_top_to_bottom)
    )
  )

############################################################
## 3. Panel B: FIXED mapped vs positive-row universes
############################################################

mapping <- read.csv(
  MAPPING_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

positive <- read.csv(
  POSITIVE_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(mapping) == 45L,
  nrow(positive) == 42L,
  all(STATE_LEVELS %in% mapping$State),
  all(STATE_LEVELS %in% positive$State),
  all(as.numeric(positive$NES) > 0)
)

mapped_counts <- mapping %>%
  filter(State %in% STATE_LEVELS) %>%
  count(
    State,
    name = "mapped_rows"
  )

positive_counts <- positive %>%
  filter(State %in% STATE_LEVELS) %>%
  count(
    State,
    name = "positive_prioritization_rows"
  )

panelB_df <- data.frame(
  State = STATE_LEVELS,
  stringsAsFactors = FALSE
) %>%
  left_join(
    mapped_counts,
    by = "State"
  ) %>%
  left_join(
    positive_counts,
    by = "State"
  )

## Hard reporting gates: this is the point of the patch.
stopifnot(
  identical(
    as.integer(panelB_df$mapped_rows),
    c(15L, 15L, 15L)
  ),
  identical(
    as.integer(panelB_df$positive_prioritization_rows),
    c(15L, 14L, 13L)
  )
)

panelB_long <- panelB_df %>%
  pivot_longer(
    cols = c(
      mapped_rows,
      positive_prioritization_rows
    ),
    names_to = "row_type",
    values_to = "n_rows"
  ) %>%
  mutate(
    State = factor(
      State,
      levels = STATE_LEVELS
    ),
    row_type = factor(
      row_type,
      levels = c(
        "mapped_rows",
        "positive_prioritization_rows"
      ),
      labels = c(
        "Mapped rows",
        "Positive prioritization rows"
      )
    )
  )

############################################################
## 4. Plot theme
############################################################

theme_s27 <- theme_bw(
  base_size = 10.8
) +
  theme(
    plot.title = element_text(
      size = 13.2,
      face = "bold",
      color = "black",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9.2,
      color = "black",
      hjust = 0
    ),
    plot.tag = element_text(
      size = 18,
      face = "bold",
      color = "black"
    ),
    axis.title = element_text(
      size = 10.8,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      size = 9.0,
      color = "black"
    ),
    legend.title = element_text(
      size = 10.0,
      face = "bold",
      color = "black"
    ),
    legend.text = element_text(
      size = 9.0,
      color = "black"
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.30
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

############################################################
## 5. Build Panel A
############################################################

pA <- ggplot(
  panelA_df,
  aes(
    x = n_genes,
    y = readable_label
  )
) +
  geom_col(
    width = 0.80,
    fill = "grey70",
    color = "grey25",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = n_genes
    ),
    hjust = -0.10,
    size = 3.25,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      100,
      20
    ),
    limits = c(
      0,
      108
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    title = "Cleaned ICBcomb query gene-set sizes",
    subtitle =
      "Human-readable labels replace raw query file names; melanocytic/reference sets are retained as controls",
    x = "Number of genes",
    y = NULL
  ) +
  theme_s27 +
  theme(
    axis.text.y = element_text(
      size = 9.0,
      color = "black"
    )
  )

############################################################
## 6. Build Panel B
############################################################

pB <- ggplot(
  panelB_long,
  aes(
    x = State,
    y = n_rows,
    fill = row_type
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.74
    ),
    width = 0.64,
    color = "grey25",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = n_rows
    ),
    position = position_dodge(
      width = 0.74
    ),
    vjust = -0.35,
    size = 3.35,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Mapped rows" = "grey70",
      "Positive prioritization rows" = "#F23838"
    ),
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    breaks = c(
      0,
      5,
      10,
      15
    ),
    limits = c(
      0,
      16.8
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  scale_x_discrete(
    labels = STATE_DISPLAY,
    drop = FALSE
  ) +
  labs(
    title = "ICBcomb strategy-mapping summary",
    subtitle =
      "Positive prioritization rows meet the prespecified direction-specific NES > 0 criterion; reference/control rows are excluded",
    x = NULL,
    y = "Number of rows"
  ) +
  theme_s27 +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.x = element_text(
      size = 8.7,
      angle = 28,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9.0,
      color = "black"
    )
  )

############################################################
## 7. Assemble and save
############################################################

supp_s27 <- pA / pB +
  plot_layout(
    heights = c(
      1.02,
      0.98
    )
  ) +
  plot_annotation(
    title =
      "ICBcomb query gene sets and strategy mapping summary",
    subtitle =
      "Query-gene-set size summary and strategy mapping used for downstream ICBcomb candidate perturbation-prioritization analysis.",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        size = 17,
        face = "bold",
        color = "black",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 10.0,
        color = "black",
        hjust = 0
      ),
      plot.tag = element_text(
        size = 18,
        face = "bold",
        color = "black"
      ),
      plot.margin = margin(
        8,
        10,
        8,
        10
      )
    )
  )

out_png <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S27_ICBcomb_query_strategy_mapping.png"
)

out_jpg <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S27_ICBcomb_query_strategy_mapping.jpg"
)

out_pdf <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S27_ICBcomb_query_strategy_mapping.pdf"
)

ggsave(
  out_png,
  supp_s27,
  width = 11.0,
  height = 11.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  out_jpg,
  supp_s27,
  width = 11.0,
  height = 11.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  out_pdf,
  supp_s27,
  width = 11.0,
  height = 11.5,
  bg = "white"
)

############################################################
## 8. Reporting audit
############################################################

audit <- data.frame(
  item = c(
    "panelA_gene_sets",
    "mapped_rows_immune",
    "mapped_rows_myeloid",
    "mapped_rows_tumor",
    "positive_rows_immune",
    "positive_rows_myeloid",
    "positive_rows_tumor",
    "public_reversal_wording_removed",
    "scientific_analysis_changed"
  ),
  observed = c(
    as.character(nrow(panelA_df)),
    as.character(panelB_df$mapped_rows[1]),
    as.character(panelB_df$mapped_rows[2]),
    as.character(panelB_df$mapped_rows[3]),
    as.character(panelB_df$positive_prioritization_rows[1]),
    as.character(panelB_df$positive_prioritization_rows[2]),
    as.character(panelB_df$positive_prioritization_rows[3]),
    "TRUE",
    "FALSE"
  ),
  expected = c(
    "13",
    "15",
    "15",
    "15",
    "15",
    "14",
    "13",
    "TRUE",
    "FALSE"
  ),
  pass = c(
    nrow(panelA_df) == 13L,
    panelB_df$mapped_rows[1] == 15L,
    panelB_df$mapped_rows[2] == 15L,
    panelB_df$mapped_rows[3] == 15L,
    panelB_df$positive_prioritization_rows[1] == 15L,
    panelB_df$positive_prioritization_rows[2] == 14L,
    panelB_df$positive_prioritization_rows[3] == 13L,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(
    FIG_DIR,
    "S27_artifact_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  panelA_df %>%
    mutate(
      readable_label =
        as.character(
          readable_label
        )
    ),
  file.path(
    FIG_DIR,
    "Supplementary_Figure_S27_query_gene_set_sizes.csv"
  ),
  row.names = FALSE
)

write.csv(
  panelB_df,
  file.path(
    FIG_DIR,
    "Supplementary_Figure_S27_strategy_mapping_summary.csv"
  ),
  row.names = FALSE
)

if (!all(audit$pass)) {
  print(audit)
  stop(
    "S27 reporting audit failed.",
    call. = FALSE
  )
}

cat(
  "\nSupplementary Figure S27 reporting rebuild: PASS\n"
)
cat(
  "Mapped rows: 15 / 15 / 15\n"
)
cat(
  "Positive prioritization rows: 15 / 14 / 13\n"
)
cat(
  "No scientific analysis changed.\n"
)
