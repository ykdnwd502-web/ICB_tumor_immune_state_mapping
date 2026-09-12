############################################################
## Supplementary Figure S4. Robustness of predefined tumor–immune state scoring
##
## Purpose:
##   Render Supplementary Figure S4 with text size uniformly scaled up by 15%:
##   A) Gene-pool Jaccard overlap among predefined states
##   B) Alternative scoring-method concordance
##   C) Gene-subsampling robustness across fractions
##
## Output:
##   D:/ICB_resistance_project/results/figures/state_scoring_robustness/
##     - Supplementary Figure S4. Robustness of predefined tumor–immune state scoring.png
##     - Supplementary Figure S4. Robustness of predefined tumor–immune state scoring.jpg
##     - Supplementary Figure S4. Robustness of predefined tumor–immune state scoring.pdf
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(patchwork)
})

############################################################
## 0. User settings
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

table_dir <- file.path(project_dir, "results", "tables", "state_scoring_robustness")
fig_dir   <- file.path(project_dir, "results", "figures", "state_scoring_robustness")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

############################################################
## 1. Helpers
############################################################

safe_ggsave <- function(file, plot, width, height, dpi = 300, device = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE,
    bg = "white",
    device = device
  )
  message("Saved: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

## 基础字体大小从 10 放大 15% 调整为 11.5
theme_s4 <- function(base_size = 11.5) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#E9E9E9", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(hjust = 0.5, color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      strip.background = element_rect(fill = "#D9D9D9", color = "grey50", linewidth = 0.4),
      strip.text = element_text(face = "bold", color = "black"),
      plot.tag = element_text(face = "bold", size = base_size * 1.7, color = "black"),
      plot.margin = margin(5, 8, 5, 8)
    )
}

state_display <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Immune_defective_cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Myeloid_Treg_immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Tumor_dedifferentiated" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation",
  "Melanocytic_differentiated" = "melanocytic differentiation"
)

state_order_candidates <- c(
  "Immune_defective_Cold",
  "Immune_defective_cold",
  "Myeloid_Treg_Immunosuppressive",
  "Myeloid_Treg_immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Tumor_dedifferentiated",
  "Melanocytic_Differentiation",
  "Melanocytic_differentiated"
)

display_state <- function(x, wrap_width = NULL) {
  out <- unname(state_display[as.character(x)])
  out[is.na(out)] <- str_replace_all(as.character(x[is.na(out)]), "_", " ")
  if (!is.null(wrap_width)) out <- str_wrap(out, width = wrap_width)
  out
}

method_display <- c(
  "GSVA_sSGSEA_Optional" = "GSVA/ssGSEA",
  "GSVA_ssGSEA_Optional" = "GSVA/ssGSEA",
  "GSVA_ssGSEA" = "GSVA/ssGSEA",
  "rank_mean_singscore_like" = "rank-mean score",
  "union_gene_zmean" = "union-gene z-mean",
  "zmean_union_gene_pool" = "union-gene z-mean"
)

display_method <- function(x) {
  out <- unname(method_display[as.character(x)])
  out[is.na(out)] <- str_replace_all(as.character(x[is.na(out)]), "_", " ")
  out
}

############################################################
## 2. Input files
############################################################

overlap_long_file   <- file.path(table_dir, "final_state_gene_pool_overlap_long.csv")
jaccard_matrix_file <- file.path(table_dir, "final_state_gene_pool_Jaccard_matrix.csv")
method_file         <- file.path(table_dir, "scoring_method_concordance_vs_primary.csv")
sub_iter_file       <- file.path(table_dir, "gene_subsampling_all_iterations.csv")
sub_summary_file    <- file.path(table_dir, "gene_subsampling_stability_summary.csv")

if (!file.exists(overlap_long_file) && !file.exists(jaccard_matrix_file)) {
  stop("找不到Jaccard重叠文件: ", table_dir)
}
if (!file.exists(method_file)) {
  stop("找不到打分方法一致性文件: ", method_file)
}
if (!file.exists(sub_iter_file) && !file.exists(sub_summary_file)) {
  stop("找不到基因亚抽样稳健性文件: ", table_dir)
}

message("Using method concordance table: ", method_file)

############################################################
## 3. Panel A: Jaccard overlap heatmap
############################################################

if (file.exists(overlap_long_file)) {
  overlap_long <- readr::read_csv(overlap_long_file, show_col_types = FALSE)
} else {
  jac_mat <- read.csv(jaccard_matrix_file, check.names = FALSE)
  colnames(jac_mat)[1] <- "State1"
  overlap_long <- jac_mat %>%
    pivot_longer(
      cols = -State1,
      names_to = "State2",
      values_to = "Jaccard"
    )
}

state_order <- unique(c(as.character(overlap_long$State1), as.character(overlap_long$State2)))
state_order <- state_order_candidates[state_order_candidates %in% state_order]
if (length(state_order) == 0) {
  state_order <- unique(as.character(overlap_long$State1))
}

overlap_plot <- overlap_long %>%
  mutate(
    State1 = factor(State1, levels = state_order),
    State2 = factor(State2, levels = state_order),
    State1Label = factor(display_state(State1, wrap_width = 34), levels = rev(display_state(state_order, wrap_width = 34))),
    State2Label = factor(display_state(State2, wrap_width = 34), levels = display_state(state_order, wrap_width = 34)),
    label = ifelse(is.na(Jaccard), "", sprintf("%.2f", Jaccard))
  )

pA <- ggplot(overlap_plot, aes(x = State2Label, y = State1Label, fill = Jaccard)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 3.8, family = "sans", color = "black") + # 对应放大
  scale_fill_gradient(
    low = "white",
    high = "#2166AC",
    limits = c(0, 1),
    na.value = "grey90",
    name = "Jaccard\nindex"
  ) +
  labs(
    title = "Gene-pool overlap among predefined tumor–immune states",
    subtitle = "Jaccard index based on final-state union gene pools",
    x = NULL,
    y = NULL
  ) +
  theme_s4(base_size = 11.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 10.0),
    axis.text.y = element_text(size = 10.2),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.6)
  )

############################################################
## 4. Panel B: alternative scoring-method concordance
############################################################

method_conc <- readr::read_csv(method_file, show_col_types = FALSE)

method_plot <- method_conc %>%
  mutate(
    Dataset = factor(Dataset, levels = c("GSE244982", "GSE78220")),
    FinalState = factor(FinalState, levels = state_order),
    StateLabel = factor(display_state(FinalState, wrap_width = 38), levels = rev(display_state(state_order, wrap_width = 38))),
    MethodLabel = factor(
      display_method(AlternativeScoringMethod),
      levels = c("GSVA/ssGSEA", "rank-mean score", "union-gene z-mean")
    ),
    label = ifelse(is.na(Spearman_rho_vs_primary), "", sprintf("%.2f", Spearman_rho_vs_primary)),
    Spearman_rho_display = ifelse(
      is.na(Spearman_rho_vs_primary),
      NA_real_,
      pmin(pmax(Spearman_rho_vs_primary, 0.90), 1.00)
    )
  ) %>%
  filter(!is.na(MethodLabel))

pB <- ggplot(method_plot, aes(x = MethodLabel, y = StateLabel, fill = Spearman_rho_display)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 3.7, family = "sans", color = "black") + # 对应放大
  facet_wrap(~ Dataset, ncol = 1) +
  scale_fill_gradient(
    low = "#FEE0D2",
    high = "#EF4444",
    limits = c(0.90, 1.00),
    breaks = c(0.90, 0.95, 1.00),
    labels = c("0.90", "0.95", "1.00"),
    na.value = "grey90",
    name = "Spearman\nrho"
  ) +
  labs(
    title = "Alternative scoring-method concordance",
    subtitle = "Near-perfect concordance with primary constituent z-mean score",
    x = "Alternative scoring method",
    y = NULL
  ) +
  theme_s4(base_size = 11.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9.5),
    axis.text.y = element_text(size = 9.0, lineheight = 0.92),
    axis.title.x = element_text(size = 10.5, face = "bold"),
    legend.title = element_text(size = 10.3, face = "bold"),
    legend.text = element_text(size = 9.5),
    strip.text = element_text(size = 9.8, face = "bold")
  )

############################################################
## 5. Panel C: gene-subsampling robustness
############################################################

if (file.exists(sub_iter_file)) {
  sub_all <- readr::read_csv(sub_iter_file, show_col_types = FALSE)
  if (!"Status" %in% colnames(sub_all)) sub_all$Status <- "Completed"
  
  sub_plot <- sub_all %>%
    filter(Status == "Completed") %>%
    mutate(
      Dataset = factor(Dataset, levels = c("GSE244982", "GSE78220")),
      FinalState = factor(FinalState, levels = state_order),
      StateLabel = factor(display_state(FinalState, wrap_width = 34), levels = display_state(state_order, wrap_width = 34)),
      FractionRetained = factor(
        paste0(round(FractionRetained * 100), "% retained"),
        levels = c("70% retained", "80% retained", "90% retained")
      )
    )
  
  pC <- ggplot(sub_plot, aes(x = FractionRetained, y = Spearman_rho_vs_full)) +
    geom_hline(yintercept = 0.80, linetype = "dashed", linewidth = 0.38, color = "black") +
    geom_boxplot(outlier.size = 0.35, width = 0.58, linewidth = 0.35) +
    facet_grid(Dataset ~ StateLabel) +
    coord_cartesian(ylim = c(0.75, 1.02)) +
    scale_y_continuous(breaks = c(0.75, 0.80, 0.90, 1.00)) +
    labs(
      title = "Gene-subsampling robustness of predefined state scores",
      subtitle = "One hundred random subsampling iterations per retained-gene fraction",
      x = "Gene-retention fraction",
      y = "Spearman rho vs full score"
    ) +
    theme_s4(base_size = 10.8) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9.0),
      axis.text.y = element_text(size = 9.2),
      axis.title.x = element_text(size = 10.5, face = "bold"),
      axis.title.y = element_text(size = 10.5, face = "bold"),
      strip.text.x = element_text(size = 8.4, face = "bold", lineheight = 0.90),
      strip.text.y = element_text(size = 9.2, face = "bold"),
      panel.grid.minor = element_blank()
    )
} else {
  sub_summary <- readr::read_csv(sub_summary_file, show_col_types = FALSE)
  sub_plot <- sub_summary %>%
    mutate(
      Dataset = factor(Dataset, levels = c("GSE244982", "GSE78220")),
      FinalState = factor(FinalState, levels = state_order),
      StateLabel = factor(display_state(FinalState, wrap_width = 34), levels = display_state(state_order, wrap_width = 34)),
      FractionRetained = factor(
        paste0(round(FractionRetained * 100), "% retained"),
        levels = c("70% retained", "80% retained", "90% retained")
      )
    )
  
  pC <- ggplot(sub_plot, aes(x = FractionRetained, y = median_rho)) +
    geom_hline(yintercept = 0.80, linetype = "dashed", linewidth = 0.38, color = "black") +
    geom_pointrange(aes(ymin = IQR_low, ymax = IQR_high), linewidth = 0.35, size = 0.9) +
    facet_grid(Dataset ~ StateLabel) +
    coord_cartesian(ylim = c(0.75, 1.02)) +
    scale_y_continuous(breaks = c(0.75, 0.80, 0.90, 1.00)) +
    labs(
      title = "Gene-subsampling robustness of predefined state scores",
      subtitle = "Median and interquartile range across one hundred random subsampling iterations",
      x = "Gene-retention fraction",
      y = "Spearman rho vs full score"
    ) +
    theme_s4(base_size = 10.8) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9.0),
      axis.text.y = element_text(size = 9.2),
      axis.title.x = element_text(size = 10.5, face = "bold"),
      axis.title.y = element_text(size = 10.5, face = "bold"),
      strip.text.x = element_text(size = 8.4, face = "bold", lineheight = 0.90),
      strip.text.y = element_text(size = 9.2, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

############################################################
## 6. Combined Supplementary Figure S4
############################################################

supp_s4 <- (pA | pB) / pC +
  plot_layout(heights = c(0.95, 1.30), widths = c(1.02, 1.18)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 28, family = "sans"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

base_filename <- "Supplementary Figure S4. Robustness of predefined tumor–immune state scoring"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

safe_ggsave(out_png, supp_s4, width = 15.2, height = 11.6, dpi = dpi_out) # 略微调大画布以适应放大后的文字
safe_ggsave(out_jpg, supp_s4, width = 15.2, height = 11.6, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s4, width = 15.2, height = 11.6, dpi = dpi_out, device = cairo_pdf)

cat("\nSupplementary Figure S4 generated successfully with text scaled up by 15% in PNG, JPG, and PDF formats (300 DPI)!\n")