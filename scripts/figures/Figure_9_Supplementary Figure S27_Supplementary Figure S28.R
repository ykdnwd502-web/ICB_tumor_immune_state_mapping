################--------------------------------------------
## ICBcomb Integrated Analysis and Visualization Pipeline
################--------------------------------------------

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR", unset = "")
if (!nzchar(PROJECT_DIR)) PROJECT_DIR <- "D:/ICB_resistance_project"
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)

ICB_ROOT_ENV <- Sys.getenv("ICB_ICBCOMB_ROOT", unset = "")
ROOT_CANDIDATES <- c(
  if (nzchar(ICB_ROOT_ENV)) ICB_ROOT_ENV else character(0),
  file.path(PROJECT_DIR, "modules", "ICBcomb"),
  file.path(PROJECT_DIR, "ICBcomb")
)
ROOT_CANDIDATES <- unique(ROOT_CANDIDATES[dir.exists(ROOT_CANDIDATES)])
if (length(ROOT_CANDIDATES) == 0L) {
  stop(
    "ICBcomb module directory not found. Please verify environment paths.",
    call. = FALSE
  )
}
ROOT <- normalizePath(ROOT_CANDIDATES[[1]], winslash = "/", mustWork = TRUE)
Sys.setenv(ICB_PUBLIC_ROOT = ROOT)
setwd(ROOT)

TABLE_ROOT <- file.path(ROOT, "results", "tables", "ICBcomb")
FIG_DIR    <- file.path(ROOT, "results", "figures")
RPT_DIR    <- file.path(ROOT, "results", "reports")
INPUT_ROOT <- file.path(ROOT, "inputs", "ICBcomb_INPUT_FREEZE_v1.0")

dir.create(TABLE_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(ggrepel)
  library(stringr)
})

################--------------------------------------------
## 1. Load and Standardize Input Data
################--------------------------------------------

files <- c(
  all_rows = "ICBcomb_all_platform_rows.csv",
  mapping  = "ICBcomb_strategy_mapping.csv",
  queries  = "ICBcomb_query_counts.csv",
  positive = "ICBcomb_positive_prioritization_rows.csv",
  summary  = "ICBcomb_state_strategy_summary.csv",
  figA     = "Figure9_Panel_A_query_counts.csv",
  figB     = "Figure9_Panel_B_top_strategies.csv",
  figC     = "Figure9_Panel_C_state_strategy_matrix.csv",
  figD     = "Figure9_Panel_D_state_maximum.csv"
)

R_tables <- lapply(file.path(TABLE_ROOT, files), function(x) {
  if (file.exists(x)) {
    df <- read.csv(x, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- tolower(colnames(df))
    if ("state" %in% colnames(df)) {
      df$state <- trimws(tolower(df$state))
      df$state <- gsub("–", "-", df$state)
      df$state <- gsub("strmal", "stromal", df$state)
    }
    return(df)
  } else {
    return(NULL)
  }
})
names(R_tables) <- names(files)

state_levels_main <- c(
  "immune-defective/cold",
  "myeloid-treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling"
)

state_cols <- c(
  "immune-defective/cold" = "#4DBBD5",
  "myeloid-treg immunosuppressive" = "#00A087",
  "tumor-dedifferentiation/stromal-remodeling" = "#E64B35"
)

state_label <- c(
  "immune-defective/cold" = "immune-defective/\ncold",
  "myeloid-treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling"
)

state_label_one_line <- c(
  "immune-defective/cold" = "immune-defective/cold",
  "myeloid-treg immunosuppressive" = "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/stromal-remodeling"
)

theme_icb <- function(base_size = 10.5) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", color = "black", hjust = 0),
      strip.background = element_rect(fill = "#D9D9D9", color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      plot.margin = margin(5, 6, 5, 6)
    )
}

FONT_SCALE_MAIN <- 1.40
scale_font <- function(x) x * FONT_SCALE_MAIN
PANEL_TITLE_SIZE <- scale_font(10.2)
PANEL_A_X_TITLE_SIZE <- scale_font(9.2)

################--------------------------------------------
## 2. Generate Figure 9
################--------------------------------------------

dfA <- R_tables$figA
if (!is.null(dfA)) {
  dfA <- dfA %>% 
    mutate(
      state = trimws(tolower(state)),
      state = gsub("–", "-", state),
      state = gsub("strmal", "stromal", state),
      state = factor(state, levels = rev(state_levels_main))
    )
  
  pA <- ggplot(dfA, aes(x = n_query_genes, y = state, fill = state)) +
    geom_col(width = 0.58) +
    geom_text(aes(label = n_query_genes), hjust = -0.15, size = scale_font(3.8), color = "black") +
    scale_fill_manual(values = state_cols, drop = FALSE) +
    scale_y_discrete(labels = state_label) +
    scale_x_continuous(limits = c(0, max(dfA$n_query_genes, na.rm = TRUE) * 1.15), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "A  ICBcomb query gene sets", x = "Number of query genes", y = NULL) +
    theme_icb(base_size = scale_font(10.5)) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(size = scale_font(10.0)),
      axis.text.y = element_text(size = scale_font(10.2)),
      axis.title.x = element_text(size = PANEL_A_X_TITLE_SIZE, margin = margin(t = 4)),
      plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0)
    )
}

dfB <- R_tables$figB
if (!is.null(dfB)) {
  score_col <- intersect(c("median_reversal_score", "median_prioritization_score", "nes"), colnames(dfB))[1]
  dfB$score_val <- as.numeric(dfB[[score_col]])
  
  if (!"support_n" %in% colnames(dfB)) dfB$support_n <- 1
  if (!"rank_within_state" %in% colnames(dfB)) {
    dfB <- dfB %>% group_by(state) %>% mutate(rank_within_state = row_number()) %>% ungroup()
  }
  
  FONT_SCALE_PANEL_B <- 1.15
  scale_font_b <- function(x) x * FONT_SCALE_PANEL_B
  
  B_plot <- dfB %>%
    mutate(
      state = trimws(tolower(state)),
      state = gsub("–", "-", state),
      state = gsub("strmal", "stromal", state),
      state = factor(state, levels = state_levels_main),
      strategy_label = paste0(rank_within_state, ". ", strategy),
      value_label = paste0(sprintf("%.2f", score_val), " (n=", support_n, ")")
    ) %>%
    arrange(state, desc(rank_within_state)) %>%
    mutate(
      y_id = paste(state, strategy_label, sep = "___"),
      y_id = factor(y_id, levels = unique(y_id))
    )
  
  B_xmax <- max(B_plot$score_val, na.rm = TRUE) * 1.33
  
  pB <- ggplot(B_plot, aes(x = score_val, y = y_id, fill = state)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey45") +
    geom_col(width = 0.62, color = "black", linewidth = 0.25) +
    geom_text(aes(x = score_val + B_xmax * 0.025, label = value_label), hjust = 0, size = scale_font_b(3.1), color = "black") +
    facet_wrap(~ state, ncol = 1, scales = "free_y", labeller = labeller(state = state_label_one_line)) +
    scale_y_discrete(labels = function(x) str_wrap(str_replace(x, "^.*___", ""), width = 34)) +
    scale_fill_manual(values = state_cols, drop = FALSE) +
    coord_cartesian(xlim = c(0, B_xmax), clip = "off") +
    labs(title = "B  Top positive candidate reversal strategies", x = "Median prioritization score (NES)", y = NULL) +
    theme_icb(base_size = 10.2) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.35),
      axis.text.y = element_text(size = scale_font_b(9.2), lineheight = 0.92),
      axis.text.x = element_text(size = scale_font_b(10)),
      axis.title.x = element_text(size = scale_font_b(10.2), face = "bold", margin = margin(t = 4)),
      strip.text = element_text(size = 10.4, face = "bold", margin = margin(t = 1.5, r = 2, b = 1.5, l = 2), lineheight = 0.85),
      plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0),
      plot.margin = margin(5, 42, 5, 5)
    )
}

dfC <- R_tables$figC
if (!is.null(dfC)) {
  strategy_order <- unique(dfC$strategy)
  score_col_c <- intersect(c("median_reversal_score", "median_prioritization_score", "nes"), colnames(dfC))[1]
  dfC$score_val <- as.numeric(dfC[[score_col_c]])
  
  C_plot <- dfC %>%
    mutate(
      state = trimws(tolower(state)),
      state = gsub("–", "-", state),
      state = gsub("strmal", "stromal", state),
      state = factor(state, levels = rev(state_levels_main)),
      strategy = factor(strategy, levels = strategy_order),
      label = ifelse(is.na(score_val), "—", paste0(sprintf("%.2f", score_val), "\n", "n=", support_n))
    )
  
  C_abs <- max(abs(C_plot$score_val), na.rm = TRUE)
  if (!is.finite(C_abs) || C_abs <= 0) C_abs <- 1
  
  pC <- ggplot(C_plot, aes(x = strategy, y = state, fill = score_val)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = label), size = scale_font(3.35), lineheight = 0.88, color = "black") +
    scale_y_discrete(labels = state_label) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 18)) +
    scale_fill_gradient2(
      low = "#3B82F6", mid = "white", high = "#EF4444", 
      midpoint = 0, limits = c(-C_abs, C_abs), na.value = "#F3F4F6", 
      name = "Median\nprioritization\nscore"
    ) +
    labs(title = "C  State-specific matrix", x = NULL, y = NULL) +
    theme_icb(base_size = scale_font(10.2)) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = scale_font(8.9)),
      axis.text.y = element_text(size = scale_font(9.8)),
      legend.key.height = unit(0.43, "cm"),
      legend.key.width = unit(0.32, "cm"),
      legend.title = element_text(size = scale_font(9.7), face = "bold"),
      legend.text = element_text(size = scale_font(8.6)),
      plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0)
    )
}

dfD <- R_tables$figD
if (!is.null(dfD)) {
  score_col_d <- intersect(c("max_median_reversal_score", "max_median_prioritization_score", "max_nes"), colnames(dfD))[1]
  dfD$score_val <- as.numeric(dfD[[score_col_d]])
  dfD <- dfD %>% 
    mutate(
      state = trimws(tolower(state)),
      state = gsub("–", "-", state),
      state = gsub("strmal", "stromal", state),
      state = factor(state, levels = state_levels_main)
    )
  
  D_ymax <- max(dfD$score_val, na.rm = TRUE) * 1.20
  
  pD <- ggplot(dfD, aes(x = state, y = score_val, fill = state)) +
    geom_col(width = 0.54, color = NA) +
    geom_text(aes(label = sprintf("%.2f", score_val)), vjust = -0.35, size = scale_font(4.0) * 0.85) +
    scale_fill_manual(values = state_cols, drop = FALSE) +
    scale_x_discrete(labels = state_label) +
    scale_y_continuous(limits = c(0, D_ymax), expand = expansion(mult = c(0, 0.02))) +
    labs(title = "D  Maximum positive signal", x = NULL, y = "Maximum median score") +
    theme_icb(base_size = scale_font(10.6)) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = scale_font(9.3), lineheight = 0.9),
      axis.text.y = element_text(size = scale_font(10.0)),
      axis.title.y = element_text(size = scale_font(10.5), margin = margin(r = 5)),
      plot.title = element_text(size = PANEL_TITLE_SIZE, hjust = 0)
    )
}

fig9 <- ((pA | pB) / (pC | pD)) +
  plot_layout(heights = c(0.98, 1.28), widths = c(1.05, 1.20)) +
  plot_annotation(
    title = "ICBcomb nomination of candidate perturbation-prioritization signals for selected candidate resistance-relevant tumor–immune state queries",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13.5, margin = margin(b = 5), family = "sans"))
  )

ggsave(file.path(FIG_DIR, "Figure 9. ICBcomb nomination of candidate perturbation-prioritization signals for selected candidate resistance-relevant tumor–immune state queries.png"), fig9, width = 16.8, height = 10.0, dpi = 300, bg = "white")
ggsave(file.path(FIG_DIR, "Figure 9. ICBcomb nomination of candidate perturbation-prioritization signals for selected candidate resistance-relevant tumor–immune state queries.jpg"), fig9, width = 16.8, height = 10.0, dpi = 300, bg = "white")
ggsave(file.path(FIG_DIR, "Figure 9. ICBcomb nomination of candidate perturbation-prioritization signals for selected candidate resistance-relevant tumor–immune state queries.pdf"), fig9, width = 16.8, height = 10.0, device = cairo_pdf, bg = "white")

################--------------------------------------------
## 3. Generate Supplementary Figure S27
################--------------------------------------------

GENESET_FILE <- file.path(
  INPUT_ROOT, "data", "query_gene_set_bundle", 
  "02_final_ICBcomb_gene_sets", "final_ICBcomb_gene_set_summary_CLEAN.csv"
)

STATE_DISPLAY_SUPP <- c(
  "immune-defective/cold" = "immune-defective/\ncold",
  "myeloid-treg immunosuppressive" = "myeloid-Treg\nimmunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling"
)

if (file.exists(GENESET_FILE)) {
  gene_sets <- read.csv(GENESET_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  
  label_map <- c(
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven" = "Myeloid-Treg suppress, data-driven",
    "Melanocytic_Differentiation_REFERENCE_data_driven" = "Melanocytic reference, data-driven",
    "Immune_defective_Cold_RESTORE_data_driven" = "immune-defective/cold restore, data-driven",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven" = "tumor-dedifferentiation/stromal suppress, data-driven",
    "Immune_defective_Cold_RESTORE_hybrid" = "immune-defective/cold restore, hybrid",
    "Dedifferentiated_inflammatory_REFERENCE" = "Differentiated-inflammatory reference",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid" = "tumor-dedifferentiation/stromal suppress, hybrid",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid" = "Myeloid-Treg suppress, hybrid",
    "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated" = "tumor-dedifferentiation/stromal suppress, curated",
    "Immune_defective_Cold_RESTORE_curated" = "immune-defective/cold restore, curated",
    "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated" = "Myeloid-Treg suppress, curated",
    "Melanocytic_Differentiation_REFERENCE_hybrid" = "Melanocytic reference, hybrid",
    "Melanocytic_Differentiation_REFERENCE_curated" = "Melanocytic reference, curated"
  )
  
  panelA_df <- gene_sets %>% mutate(readable_label = unname(label_map[GeneSet]))
  panelA_top_to_bottom <- c(
    "Myeloid-Treg suppress, data-driven", "Melanocytic reference, data-driven", "immune-defective/cold restore, data-driven",
    "tumor-dedifferentiation/stromal suppress, data-driven", "immune-defective/cold restore, hybrid", "Differentiated-inflammatory reference",
    "tumor-dedifferentiation/stromal suppress, hybrid", "Myeloid-Treg suppress, hybrid", "tumor-dedifferentiation/stromal suppress, curated",
    "immune-defective/cold restore, curated", "Myeloid-Treg suppress, curated", "Melanocytic reference, hybrid", "Melanocytic reference, curated"
  )
  panelA_df <- panelA_df %>% mutate(readable_label = factor(readable_label, levels = rev(panelA_top_to_bottom)))
  
  mapped_counts <- R_tables$mapping %>% 
    mutate(state = trimws(tolower(state)), state = gsub("–", "-", state), state = gsub("strmal", "stromal", state)) %>%
    filter(state %in% state_levels_main) %>% 
    count(state, name = "mapped_rows")
  
  positive_counts <- R_tables$positive %>% 
    mutate(state = trimws(tolower(state)), state = gsub("–", "-", state), state = gsub("strmal", "stromal", state)) %>%
    filter(state %in% state_levels_main) %>% 
    count(state, name = "positive_prioritization_rows")
  
  panelB_df <- data.frame(state = state_levels_main, stringsAsFactors = FALSE) %>%
    left_join(mapped_counts, by = "state") %>%
    left_join(positive_counts, by = "state")
  
  panelB_long <- panelB_df %>%
    pivot_longer(cols = c(mapped_rows, positive_prioritization_rows), names_to = "row_type", values_to = "n_rows") %>%
    mutate(
      state = factor(state, levels = state_levels_main),
      row_type = factor(row_type, levels = c("mapped_rows", "positive_prioritization_rows"), labels = c("Mapped rows", "Positive prioritization rows"))
    )
  
  theme_supp <- theme_bw(base_size = 10.8) + theme(
    plot.title = element_text(size = 13.2, face = "bold", color = "black", hjust = 0),
    plot.subtitle = element_text(size = 9.2, color = "black", hjust = 0),
    plot.tag = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 10.8, face = "bold", color = "black"),
    axis.text = element_text(size = 9.0, color = "black"),
    legend.title = element_text(size = 10.0, face = "bold", color = "black"),
    legend.text = element_text(size = 9.0, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )
  
  pA_s27 <- ggplot(panelA_df, aes(x = n_genes, y = readable_label)) +
    geom_col(width = 0.80, fill = "grey70", color = "grey25", linewidth = 0.35) +
    geom_text(aes(label = n_genes), hjust = -0.10, size = 3.25, color = "black") +
    scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 108), expand = expansion(mult = c(0, 0))) +
    labs(title = "Cleaned ICBcomb query gene-set sizes", subtitle = "Human-readable labels replace raw query file names", x = "Number of genes", y = NULL) +
    theme_supp
  
  pB_s27 <- ggplot(panelB_long, aes(x = state, y = n_rows, fill = row_type)) +
    geom_col(position = position_dodge(width = 0.74), width = 0.64, color = "grey25", linewidth = 0.35) +
    geom_text(aes(label = n_rows), position = position_dodge(width = 0.74), vjust = -0.35, size = 3.35, color = "black") +
    scale_fill_manual(values = c("Mapped rows" = "grey70", "Positive prioritization rows" = "#F23838"), drop = FALSE, name = NULL) +
    scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 16.8), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(labels = STATE_DISPLAY_SUPP, drop = FALSE) +
    labs(title = "ICBcomb strategy-mapping summary", subtitle = "Positive prioritization rows meet NES > 0 criterion", x = NULL, y = "Number of rows") +
    theme_supp + theme(legend.position = "top", legend.justification = "left", axis.text.x = element_text(size = 8.7, angle = 28, hjust = 1, vjust = 1, color = "black"))
  
  supp_s27 <- pA_s27 / pB_s27 + plot_layout(heights = c(1.02, 0.98)) + plot_annotation(
    title = "ICBcomb query gene sets and strategy mapping summary",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 15.5, face = "bold", color = "black", hjust = 0))
  )
  
  ggsave(file.path(FIG_DIR, "Supplementary Figure S27. ICBcomb query gene sets and strategy mapping summary.png"), supp_s27, width = 11.0, height = 11.5, dpi = 300, bg = "white")
  ggsave(file.path(FIG_DIR, "Supplementary Figure S27. ICBcomb query gene sets and strategy mapping summary.jpg"), supp_s27, width = 11.0, height = 11.5, dpi = 300, bg = "white")
  ggsave(file.path(FIG_DIR, "Supplementary Figure S27. ICBcomb query gene sets and strategy mapping summary.pdf"), supp_s27, width = 11.0, height = 11.5, device = cairo_pdf, bg = "white")
}

################--------------------------------------------
## 4. Generate Supplementary Figure S28
################--------------------------------------------

TABLE_FILE_S28 <- file.path(TABLE_ROOT, "ICBcomb_state_strategy_summary.csv")

if (file.exists(TABLE_FILE_S28)) {
  STRATEGY_LEVELS <- c("Cytokine/IL-15 + ICB", "Avadomide/IMiD + ICB", "Epigenetic therapy + ICB", "Ascorbic acid + ICB")
  STRATEGY_SHORT <- c("Cytokine/IL-15 + ICB" = "Cytokine/IL-15", "Avadomide/IMiD + ICB" = "IMiD", "Epigenetic therapy + ICB" = "Epigenetic", "Ascorbic acid + ICB" = "Ascorbic acid")
  STRATEGY_COLORS <- c("Cytokine/IL-15 + ICB" = "#4DBBD5", "Avadomide/IMiD + ICB" = "#00A087", "Epigenetic therapy + ICB" = "#E64B35", "Ascorbic acid + ICB" = "#3C5488")
  
  s28_raw <- read.csv(TABLE_FILE_S28, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(s28_raw) <- tolower(colnames(s28_raw))
  
  s28_df <- s28_raw %>%
    mutate(
      state = trimws(tolower(state)), 
      state = gsub("–", "-", state),
      state = gsub("strmal", "stromal", state)
    ) %>%
    transmute(state = state, strategy = strategy, support_rows = as.numeric(support_n), median_score = as.numeric(median_prioritization_score)) %>%
    filter(state %in% state_levels_main, strategy %in% STRATEGY_LEVELS)
  
  s28_grid <- expand_grid(state = state_levels_main, strategy = STRATEGY_LEVELS) %>%
    left_join(s28_df, by = c("state", "strategy")) %>%
    mutate(
      state = factor(state, levels = rev(state_levels_main)), 
      strategy = factor(strategy, levels = STRATEGY_LEVELS)
    )
  
  s28_points <- s28_grid %>% filter(is.finite(support_rows), is.finite(median_score), support_rows > 0)
  
  theme_supp_s28 <- theme_bw(base_size = 10.8) + theme(
    plot.title = element_text(size = 13.2, face = "bold", color = "black", hjust = 0),
    plot.subtitle = element_text(size = 9.2, color = "black", hjust = 0),
    plot.tag = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 10.8, face = "bold", color = "black"),
    axis.text = element_text(size = 9.0, color = "black"),
    legend.title = element_text(size = 10.0, face = "bold", color = "black"),
    legend.text = element_text(size = 9.0, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )
  
  score_min <- min(s28_points$median_score, na.rm = TRUE)
  score_max <- max(s28_points$median_score, na.rm = TRUE)
  
  pA_s28 <- ggplot(s28_grid, aes(x = strategy, y = state)) +
    geom_point(data = s28_points, aes(size = support_rows, color = median_score), alpha = 0.88, stroke = 0.35) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(labels = STATE_DISPLAY_SUPP, drop = FALSE) +
    scale_size_continuous(range = c(2.2, 9.5), breaks = c(2, 4, 6, 8, 10), name = "Support\nrows") +
    scale_color_gradient(low = "#FEE0D2", high = "#E31A1C", name = "Median ICBcomb\nprioritization\nscore", limits = c(score_min, score_max)) +
    labs(title = "Support landscape of positive ICBcomb prioritization signals", subtitle = "Bubble size indicates support rows; color indicates median score", x = NULL, y = NULL) +
    theme_supp_s28 + theme(axis.text.x = element_text(size = 8.4, angle = 38, hjust = 1, vjust = 1, color = "black"), legend.position = "right")
  
  panelB_df <- s28_df %>% 
    mutate(
      state = factor(state, levels = state_levels_main), 
      strategy = factor(strategy, levels = STRATEGY_LEVELS), 
      label_short = unname(STRATEGY_SHORT[as.character(strategy)])
    )
  
  pB_s28 <- ggplot(panelB_df, aes(x = support_rows, y = median_score)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey35") +
    geom_point(aes(color = strategy), size = 3.2, alpha = 0.92, stroke = 0.35) +
    ggrepel::geom_text_repel(aes(label = label_short), size = 3.0, color = "black", min.segment.length = 0, segment.color = "grey70", seed = 28, show.legend = FALSE) +
    scale_color_manual(values = STRATEGY_COLORS, breaks = STRATEGY_LEVELS, labels = STRATEGY_SHORT, name = "Strategy category") +
    facet_wrap(~ state, ncol = 1, labeller = labeller(state = STATE_DISPLAY_SUPP)) +
    scale_x_continuous(breaks = pretty_breaks(n = 4), expand = expansion(mult = c(0.08, 0.16))) +
    scale_y_continuous(breaks = pretty_breaks(n = 5), expand = expansion(mult = c(0.08, 0.18))) +
    labs(title = "ICBcomb prioritization score and supporting rows by state", x = "Supporting rows", y = "Median ICBcomb prioritization score") +
    theme_supp_s28 + theme(strip.background = element_rect(fill = "grey85", color = "grey35"), strip.text = element_text(size = 9.5, face = "bold", color = "black"), legend.position = "bottom", legend.justification = "left")
  
  supp_s28 <- pA_s28 / pB_s28 + plot_layout(heights = c(0.82, 1.28)) + plot_annotation(
    title = "ICBcomb candidate perturbation-prioritization signals across resistance-associated tumor-immune states",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 12.0, face = "bold", color = "black", hjust = 0, margin = margin(b = 6))
    )
  )
  
  ggsave(file.path(FIG_DIR, "Supplementary Figure S28. ICBcomb candidate perturbation-prioritization signals across resistance-associated tumor-immune states.png"), supp_s28, width = 9.3, height = 13.4, dpi = 300, bg = "white")
  ggsave(file.path(FIG_DIR, "Supplementary Figure S28. ICBcomb candidate perturbation-prioritization signals across resistance-associated tumor-immune states.jpg"), supp_s28, width = 9.3, height = 13.4, dpi = 300, bg = "white")
  ggsave(file.path(FIG_DIR, "Supplementary Figure S28. ICBcomb candidate perturbation-prioritization signals across resistance-associated tumor-immune states.pdf"), supp_s28, width = 9.3, height = 13.4, device = cairo_pdf, bg = "white")
}
