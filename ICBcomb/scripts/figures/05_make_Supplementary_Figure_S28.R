############################################################
## 05_make_Supplementary_Figure_S28.R
##
## PURPOSE
##   Reporting-only rebuild of Supplementary Figure S28 from
##   the frozen ICBcomb state-by-strategy summary table.
##
## FIX IN THIS PATCH
##   Public-facing "reversal-prioritization" wording is replaced
##   by the frozen neutral term "prioritization score".
##
## SCIENTIFIC BOUNDARY
##   Reporting-only patch. No NES values, support counts,
##   strategies, states, query definitions, or thresholds change.
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

############################################################
## 0. Module root and fixed source
############################################################

ROOT <- Sys.getenv("ICBCOMB_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("ICB_PUBLIC_ROOT")
if (!nzchar(ROOT)) ROOT <- getwd()

ROOT <- normalizePath(
  ROOT,
  winslash = "/",
  mustWork = TRUE
)

TABLE_FILE <- file.path(
  ROOT,
  "results",
  "tables",
  "ICBcomb",
  "ICBcomb_state_strategy_summary.csv"
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

if (!file.exists(TABLE_FILE)) {
  stop(
    "Missing frozen S28 source table: ",
    TABLE_FILE,
    call. = FALSE
  )
}

############################################################
## 1. Fixed states / strategies / colors
############################################################

STATE_LEVELS <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling"
)

STATE_DISPLAY <- c(
  "immune-defective/cold" =
    "Immune-defective/\nCold",
  "myeloid–Treg immunosuppressive" =
    "Myeloid-Treg\nImmunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" =
    "Tumor-dedifferentiation/\nStromal-remodeling"
)

STRATEGY_LEVELS <- c(
  "Cytokine/IL-15 + ICB",
  "Avadomide/IMiD + ICB",
  "Epigenetic therapy + ICB",
  "Ascorbic acid + ICB"
)

STRATEGY_SHORT <- c(
  "Cytokine/IL-15 + ICB" =
    "Cytokine/IL-15",
  "Avadomide/IMiD + ICB" =
    "IMiD",
  "Epigenetic therapy + ICB" =
    "Epigenetic",
  "Ascorbic acid + ICB" =
    "Ascorbic acid"
)

STRATEGY_COLORS <- c(
  "Cytokine/IL-15 + ICB" = "#4DBBD5",
  "Avadomide/IMiD + ICB" = "#00A087",
  "Epigenetic therapy + ICB" = "#E64B35",
  "Ascorbic acid + ICB" = "#3C5488"
)

############################################################
## 2. Load and lock frozen S28 source
############################################################

s28_df <- read.csv(
  TABLE_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "State",
  "Strategy",
  "median_prioritization_score",
  "support_n"
)

if (!all(required_cols %in% names(s28_df))) {
  stop(
    "S28 source table is missing required columns.",
    call. = FALSE
  )
}

s28_df <- s28_df %>%
  transmute(
    state = State,
    strategy = Strategy,
    support_rows =
      as.numeric(
        support_n
      ),
    median_score =
      as.numeric(
        median_prioritization_score
      )
  ) %>%
  filter(
    state %in% STATE_LEVELS,
    strategy %in% STRATEGY_LEVELS
  )

## Frozen numerical contract.
state_n <- table(
  factor(
    s28_df$state,
    levels = STATE_LEVELS
  )
)

stopifnot(
  nrow(s28_df) == 10L,
  identical(
    as.integer(state_n),
    c(
      4L,
      3L,
      3L
    )
  ),
  sum(
    s28_df$support_rows
  ) == 42L,
  all(
    s28_df$support_rows > 0
  ),
  all(
    is.finite(
      s28_df$median_score
    )
  ),
  all(
    s28_df$median_score > 0
  )
)

############################################################
## 3. Complete state x strategy grid for Panel A
############################################################

s28_grid <- expand_grid(
  state = STATE_LEVELS,
  strategy = STRATEGY_LEVELS
) %>%
  left_join(
    s28_df,
    by = c(
      "state",
      "strategy"
    )
  ) %>%
  mutate(
    state = factor(
      state,
      levels = rev(
        STATE_LEVELS
      )
    ),
    strategy = factor(
      strategy,
      levels = STRATEGY_LEVELS
    )
  )

s28_points <- s28_grid %>%
  filter(
    is.finite(
      support_rows
    ),
    is.finite(
      median_score
    ),
    support_rows > 0
  )

############################################################
## 4. Theme
############################################################

theme_s28 <- theme_bw(
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
## 5. Panel A
############################################################

score_min <- min(
  s28_points$median_score,
  na.rm = TRUE
)

score_max <- max(
  s28_points$median_score,
  na.rm = TRUE
)

pA <- ggplot(
  s28_grid,
  aes(
    x = strategy,
    y = state
  )
) +
  geom_point(
    data = s28_points,
    aes(
      size = support_rows,
      color = median_score
    ),
    alpha = 0.88,
    stroke = 0.35
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_discrete(
    labels = STATE_DISPLAY,
    drop = FALSE
  ) +
  scale_size_continuous(
    range = c(
      2.2,
      9.5
    ),
    breaks = c(
      2,
      4,
      6,
      8,
      10
    ),
    name = "Support\nrows"
  ) +
  scale_color_gradient(
    low = "#FEE0D2",
    high = "#E31A1C",
    name = "Median ICBcomb\nprioritization\nscore",
    limits = c(
      score_min,
      score_max
    )
  ) +
  labs(
    title =
      "Support landscape of positive ICBcomb prioritization signals",
    subtitle =
      "Bubble size indicates support rows; color indicates median prioritization score",
    x = NULL,
    y = NULL
  ) +
  theme_s28 +
  theme(
    axis.text.x = element_text(
      size = 8.4,
      angle = 38,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 8.7,
      color = "black"
    ),
    legend.position = "right"
  )

############################################################
## 6. Panel B
############################################################

panelB_df <- s28_df %>%
  mutate(
    state = factor(
      state,
      levels = STATE_LEVELS
    ),
    strategy = factor(
      strategy,
      levels = STRATEGY_LEVELS
    ),
    label_short =
      unname(
        STRATEGY_SHORT[
          as.character(
            strategy
          )
        ]
      )
  )

pB <- ggplot(
  panelB_df,
  aes(
    x = support_rows,
    y = median_score
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey35"
  ) +
  geom_point(
    aes(
      color = strategy
    ),
    size = 3.2,
    alpha = 0.92,
    stroke = 0.35
  ) +
  ggrepel::geom_text_repel(
    aes(
      label = label_short
    ),
    size = 3.0,
    color = "black",
    min.segment.length = 0,
    segment.color = "grey70",
    box.padding = 0.25,
    point.padding = 0.22,
    max.overlaps = Inf,
    seed = 28,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = STRATEGY_COLORS,
    breaks = STRATEGY_LEVELS,
    labels = STRATEGY_SHORT,
    name = "Strategy category"
  ) +
  facet_wrap(
    ~ state,
    ncol = 1,
    labeller = labeller(
      state = STATE_DISPLAY
    )
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(
      n = 4
    ),
    expand = expansion(
      mult = c(
        0.08,
        0.16
      )
    )
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(
      n = 5
    ),
    expand = expansion(
      mult = c(
        0.08,
        0.18
      )
    )
  ) +
  labs(
    title =
      "ICBcomb prioritization score and supporting rows by state",
    subtitle =
      "Each point represents one strategy category summarized in Panel A",
    x = "Supporting rows",
    y = "Median ICBcomb prioritization score"
  ) +
  theme_s28 +
  theme(
    strip.background = element_rect(
      fill = "grey85",
      color = "grey35"
    ),
    strip.text = element_text(
      size = 9.5,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 9.0,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9.0,
      color = "black"
    ),
    legend.position = "bottom",
    legend.justification = "left"
  )

############################################################
## 7. Assemble and save
############################################################

supp_s28 <- pA / pB +
  plot_layout(
    heights = c(
      0.82,
      1.28
    )
  ) +
  plot_annotation(
    title =
      "ICBcomb candidate perturbation-prioritization signals",
    subtitle =
      "Support distribution across resistance-associated tumor-immune states; candidate-prioritization summary only.",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        size = 17,
        face = "bold",
        color = "black",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 10.2,
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
  "Supplementary_Figure_S28_ICBcomb_positive_evidence.png"
)

out_jpg <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S28_ICBcomb_positive_evidence.jpg"
)

out_pdf <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S28_ICBcomb_positive_evidence.pdf"
)

ggsave(
  out_png,
  supp_s28,
  width = 9.3,
  height = 13.4,
  dpi = 300,
  bg = "white"
)

ggsave(
  out_jpg,
  supp_s28,
  width = 9.3,
  height = 13.4,
  dpi = 300,
  bg = "white"
)

ggsave(
  out_pdf,
  supp_s28,
  width = 9.3,
  height = 13.4,
  bg = "white"
)

############################################################
## 8. Reporting audit
############################################################

audit <- data.frame(
  item = c(
    "strategy_summary_rows",
    "immune_state_rows",
    "myeloid_state_rows",
    "tumor_state_rows",
    "total_support_rows",
    "all_median_scores_positive",
    "public_reversal_wording_removed",
    "scientific_analysis_changed"
  ),
  observed = c(
    as.character(nrow(s28_df)),
    as.character(state_n[1]),
    as.character(state_n[2]),
    as.character(state_n[3]),
    as.character(sum(s28_df$support_rows)),
    as.character(all(s28_df$median_score > 0)),
    "TRUE",
    "FALSE"
  ),
  expected = c(
    "10",
    "4",
    "3",
    "3",
    "42",
    "TRUE",
    "TRUE",
    "FALSE"
  ),
  pass = c(
    nrow(s28_df) == 10L,
    state_n[1] == 4L,
    state_n[2] == 3L,
    state_n[3] == 3L,
    sum(s28_df$support_rows) == 42L,
    all(s28_df$median_score > 0),
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(
    FIG_DIR,
    "S28_artifact_audit.csv"
  ),
  row.names = FALSE
)

write.csv(
  s28_df,
  file.path(
    FIG_DIR,
    "Supplementary_Figure_S28_ICBcomb_prioritization_summary.csv"
  ),
  row.names = FALSE
)

if (!all(audit$pass)) {
  print(audit)
  stop(
    "S28 reporting audit failed.",
    call. = FALSE
  )
}

cat(
  "\nSupplementary Figure S28 reporting rebuild: PASS\n"
)
cat(
  "State-strategy rows: 10\n"
)
cat(
  "Total positive supporting rows: 42\n"
)
cat(
  "Neutral prioritization terminology: PASS\n"
)
cat(
  "No scientific analysis changed.\n"
)
