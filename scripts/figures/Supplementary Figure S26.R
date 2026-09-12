################--------------------------------------------
## Supplementary Figure S26. kNN sensitivity summary for CosMx spatial ligand–receptor proximity analysis
##
## Purpose:
##   1) Load K10, K20, and K30 CosMx spatial validation summary tables from KNN10, KNN20, and KNN30 directories
##   2) Generate and export Supplementary Figure S26 (Panel A kNN sensitivity heatmap, Panel B robustness classification summary)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S26. kNN sensitivity summary for CosMx spatial ligand–receptor proximity analysis.png
##     - Supplementary Figure S26. kNN sensitivity summary for CosMx spatial ligand–receptor proximity analysis.jpg
##     - Supplementary Figure S26. kNN sensitivity summary for CosMx spatial ligand–receptor proximity analysis.pdf
################--------------------------------------------

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
})

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
transmute <- dplyr::transmute
bind_rows <- dplyr::bind_rows

################--------------------------------------------
## 0. Paths and user controls
################--------------------------------------------

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

validation_root <- file.path(project_dir, "results", "tables", "CosMx_spatial_LR_validation")
fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300
base_filename <- "Supplementary Figure S26. kNN sensitivity summary for CosMx spatial ligand–receptor proximity analysis"
base_family <- "sans"

## Using the newly provided folder and file structure for KNN10, KNN20, and KNN30
file_k10 <- file.path(validation_root, "KNN10", "CosMx_spatial_LR_validation_summary.csv")
file_k20 <- file.path(validation_root, "KNN20", "CosMx_spatial_LR_validation_summary.csv")
file_k30 <- file.path(validation_root, "KNN30", "CosMx_spatial_LR_validation_summary.csv")

required_files <- c(file_k10, file_k20, file_k30)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required summary input file(s):\n",
    paste(missing_files, collapse = "\n"),
    "\n\nPlease ensure the CSV files are placed in their respective KNN10, KNN20, and KNN30 directories."
  )
}

################--------------------------------------------
## 1. Read and standardize input data
################--------------------------------------------

read_one_summary <- function(path, k_label) {
  dat <- readr::read_csv(path, show_col_types = FALSE)
  
  required_cols <- c(
    "pair_family",
    "validation_priority",
    "n_evaluable_contrasts",
    "n_enriched_contrasts",
    "max_log2_spatial_enrichment",
    "max_neighbor_fraction_percent",
    "min_FDR",
    "top_contrast",
    "overall_support"
  )
  
  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop("File lacks required column(s): ", basename(path), "\n", paste(missing_cols, collapse = ", "))
  }
  
  dat %>%
    mutate(
      kNN = k_label,
      kNN_num = as.integer(str_remove(k_label, "kNN=")),
      pair_family = as.character(pair_family),
      validation_priority = as.character(validation_priority),
      top_contrast = as.character(top_contrast),
      significant = min_FDR < 0.05,
      positive_enrichment = max_log2_spatial_enrichment > 0,
      heatmap_label = sprintf("%.2f%s", max_log2_spatial_enrichment, ifelse(significant, "*", ""))
    )
}

summary_long <- bind_rows(
  read_one_summary(file_k10, "kNN=10"),
  read_one_summary(file_k20, "kNN=20"),
  read_one_summary(file_k30, "kNN=30")
)

pair_order <- summary_long %>%
  group_by(pair_family) %>%
  summarise(
    priority = dplyr::first(validation_priority),
    mean_log2 = mean(max_log2_spatial_enrichment, na.rm = TRUE),
    k20_log2 = max_log2_spatial_enrichment[kNN == "kNN=20"][1],
    .groups = "drop"
  ) %>%
  mutate(
    k20_log2 = ifelse(is.na(k20_log2), mean_log2, k20_log2),
    priority_rank = case_when(
      priority == "High" ~ 1L,
      priority == "Medium" ~ 2L,
      TRUE ~ 3L
    )
  ) %>%
  arrange(priority_rank, desc(k20_log2)) %>%
  pull(pair_family)

summary_long <- summary_long %>%
  mutate(
    pair_family = factor(pair_family, levels = rev(pair_order)),
    kNN = factor(kNN, levels = c("kNN=10", "kNN=20", "kNN=30"))
  )

robust_tbl <- summary_long %>%
  group_by(pair_family) %>%
  summarise(
    validation_priority = dplyr::first(validation_priority),
    positive_kNN_n = sum(positive_enrichment, na.rm = TRUE),
    significant_kNN_n = sum(significant, na.rm = TRUE),
    min_enriched_contrasts = min(n_enriched_contrasts, na.rm = TRUE),
    mean_neighbor_fraction_percent = mean(max_neighbor_fraction_percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    robustness_class = case_when(
      positive_kNN_n == 3 & significant_kNN_n >= 2 & min_enriched_contrasts >= 3 ~ "Robust positive",
      positive_kNN_n >= 2 & significant_kNN_n >= 1 ~ "Moderate positive",
      TRUE ~ "Limited/negative"
    ),
    robustness_class = factor(robustness_class, levels = c("Robust positive", "Moderate positive", "Limited/negative")),
    significance_label = paste0(significant_kNN_n, "/3"),
    pair_family = factor(as.character(pair_family), levels = levels(summary_long$pair_family))
  ) %>%
  arrange(pair_family)

################--------------------------------------------
## 2. Build panels
################--------------------------------------------

max_abs <- max(abs(summary_long$max_log2_spatial_enrichment), na.rm = TRUE)
max_abs <- max(2.5, ceiling(max_abs * 10) / 10)

p_heat <- ggplot(summary_long, aes(x = kNN, y = pair_family, fill = max_log2_spatial_enrichment)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = heatmap_label), size = 3.2, family = base_family) +
  scale_fill_gradient2(
    low = "#3B5BA9",
    mid = "white",
    high = "#E64B35",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    name = "Max log2\nspatial\nenrichment"
  ) +
  labs(
    title = "kNN sensitivity heatmap of candidate L–R axes",
    subtitle = "Values show maximum log2 spatial enrichment across tested source–target contrasts; * indicates FDR < 0.05",
    x = "Spatial kNN setting",
    y = "Candidate L–R axis"
  ) +
  theme_bw(base_family = base_family) +
  theme(
    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 9.2),
    axis.title = element_text(face = "bold", size = 10.5),
    axis.text.x = element_text(size = 9.0, face = "bold", color = "black"),
    axis.text.y = element_text(size = 9.0, color = "black"),
    legend.title = element_text(face = "bold", size = 9.8),
    legend.text = element_text(size = 8.7),
    panel.grid = element_blank()
  )

class_colors <- c(
  "Robust positive" = "#E64B35",
  "Moderate positive" = "#12B8B8",
  "Limited/negative" = "#BDBDBD"
)

p_robust <- ggplot(robust_tbl, aes(x = significant_kNN_n, y = pair_family, color = robustness_class, size = mean_neighbor_fraction_percent)) +
  geom_point(alpha = 0.9) +
  geom_text(aes(label = significance_label), nudge_x = 0.18, size = 3, color = "black", family = base_family) +
  scale_x_continuous(breaks = 0:3, limits = c(-0.2, 3.45), expand = expansion(mult = c(0.02, 0.08))) +
  scale_color_manual(values = class_colors, name = "Robustness\nclass") +
  scale_size_continuous(range = c(2.5, 7), name = "Mean max\nneighbor\nfraction (%)") +
  labs(
    title = "Robustness classification across kNN settings",
    subtitle = "x-axis shows the number of kNN settings with FDR < 0.05; labels show significant kNN settings out of 3",
    x = "Number of significant kNN settings",
    y = "Candidate L–R axis"
  ) +
  theme_bw(base_family = base_family) +
  theme(
    plot.title = element_text(face = "bold", size = 13.5),
    plot.subtitle = element_text(size = 9.2),
    axis.title = element_text(face = "bold", size = 10.5),
    axis.text = element_text(size = 9.0, color = "black"),
    legend.title = element_text(face = "bold", size = 9.8),
    legend.text = element_text(size = 8.7),
    panel.grid.minor = element_blank()
  )

################--------------------------------------------
## 3. Assemble and save
################--------------------------------------------

supp_s26 <- p_heat / p_robust +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 18, face = "bold", color = "black"),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 8.2
fig_height <- 11.5

ggplot2::ggsave(out_png, plot = supp_s26, width = fig_width, height = fig_height, units = "in", dpi = dpi_out, bg = "white", limitsize = FALSE)
ggplot2::ggsave(out_jpg, plot = supp_s26, width = fig_width, height = fig_height, units = "in", dpi = dpi_out, bg = "white", limitsize = FALSE)
ggplot2::ggsave(out_pdf, plot = supp_s26, width = fig_width, height = fig_height, units = "in", dpi = dpi_out, bg = "white", limitsize = FALSE, device = cairo_pdf)

################--------------------------------------------
## 4. Audit outputs
################--------------------------------------------

audit <- data.frame(
  file_k10 = file_k10,
  file_k20 = file_k20,
  file_k30 = file_k30,
  n_pairs = nrow(robust_tbl),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S26_kNN_sensitivity_summary_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)
write.csv(robust_tbl, file.path(fig_dir, "Supplementary_Figure_S26_robustness_summary_table.csv"), row.names = FALSE)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S26 script finished successfully.")