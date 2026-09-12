################--------------------------------------------
## Supplementary Figure S25. CosMx spatial proximity-constrained analysis of CellChat-nominated ligand–receptor axes
##
## Purpose:
##   1) Load CosMx spatial proximity-constrained validation table from KNN20 directory
##   2) Generate and export Supplementary Figure S25 (Panel A bubble plot across spatial contrasts, Panel B ranked spatially enriched L-R axis/contrast pairs)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S25. CosMx spatial proximity-constrained analysis of CellChat-nominated ligand–receptor axes.png
##     - Supplementary Figure S25. CosMx spatial proximity-constrained analysis of CellChat-nominated ligand–receptor axes.jpg
##     - Supplementary Figure S25. CosMx spatial proximity-constrained analysis of CellChat-nominated ligand–receptor axes.pdf
################--------------------------------------------

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(scales)
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

################--------------------------------------------
## 0. Paths and user controls
################--------------------------------------------

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300
base_filename <- "Supplementary Figure S25. CosMx spatial proximity-constrained analysis of CellChat-nominated ligand–receptor axes"

## Manual override utilizing the provided file path.
manual_s25_file <- file.path(
  project_dir,
  "results", "tables",
  "CosMx_spatial_LR_validation", "KNN20",
  "CosMx_spatial_LR_validation_all.csv"
)

n_top_panelB <- 18
fdr_cutoff <- 0.05

contrast_order <- c(
  "Dediff/Stromal high -> Myeloid-Treg high",
  "Myeloid-Treg high -> Dediff/Stromal high",
  "Dediff-only -> Myeloid-only",
  "Myeloid-only -> Dediff-only",
  "Within Both-high niche"
)

tracked_axes_order <- c(
  "MIF-CD74",
  "MIF-CD44",
  "MIF-CXCR4",
  "SPP1-CD44",
  "SPP1-ITGAV/ITGB1",
  "SPP1-ITGAV/ITGB5",
  "FN1-CD44",
  "APP-CD74",
  "GDF15-TGFBR2",
  "COL1A1-CD44",
  "COL1A2-CD44",
  "COL3A1-CD44"
)

################--------------------------------------------
## 1. Helper functions
################--------------------------------------------

read_any_table <- function(file) {
  if (!file.exists(file)) stop("File does not exist: ", file)
  
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("csv", "txt")) {
    out <- suppressMessages(readr::read_csv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "tsv") {
    out <- suppressMessages(readr::read_tsv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "rds") {
    out <- readRDS(file)
    if (!is.data.frame(out)) stop("RDS object is not a data.frame: ", file)
  } else {
    stop("Unsupported file extension: ", file)
  }
  
  as.data.frame(out, stringsAsFactors = FALSE)
}

safe_ggsave <- function(file, plot, width, height, dpi = 300, device = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE,
    device = device
  )
  message("Saved: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

find_col <- function(df, patterns, label, required = TRUE) {
  cols <- colnames(df)
  for (pat in patterns) {
    hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE, perl = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  
  if (required) {
    stop(
      "Could not find column for ", label, ". Tried patterns:\n",
      paste(patterns, collapse = "\n"),
      "\nAvailable columns:\n",
      paste(cols, collapse = ", ")
    )
  }
  
  NA_character_
}

pick_numeric_col <- function(df, preferred, label, required = TRUE) {
  for (cc in preferred) {
    if (cc %in% colnames(df)) {
      vv <- suppressWarnings(as.numeric(df[[cc]]))
      if (sum(is.finite(vv)) > 0) return(cc)
    }
  }
  
  if (required) {
    stop(
      "Could not find numeric column for ", label, ". Tried:\n",
      paste(preferred, collapse = ", "),
      "\nAvailable columns:\n",
      paste(colnames(df), collapse = ", ")
    )
  }
  
  NA_character_
}

standardize_lr_axis <- function(x) {
  x0 <- as.character(x)
  x1 <- toupper(x0)
  x1 <- gsub("\\s+", "", x1)
  x1 <- gsub("–|—", "-", x1)
  x1 <- gsub("_", "-", x1)
  x1 <- gsub("\\+", "-", x1)
  x1 <- gsub("\\(|\\)|\\[|\\]", "-", x1)
  
  dplyr::case_when(
    grepl("MIF.*CD74", x1) ~ "MIF-CD74",
    grepl("MIF.*CD44", x1) ~ "MIF-CD44",
    grepl("MIF.*CXCR4", x1) ~ "MIF-CXCR4",
    grepl("SPP1.*CD44", x1) ~ "SPP1-CD44",
    grepl("SPP1.*ITGAV.*ITGB1", x1) ~ "SPP1-ITGAV/ITGB1",
    grepl("SPP1.*ITGAV.*ITGB5", x1) ~ "SPP1-ITGAV/ITGB5",
    grepl("FN1.*CD44", x1) ~ "FN1-CD44",
    grepl("APP.*CD74", x1) ~ "APP-CD74",
    grepl("GDF15.*TGFBR2", x1) ~ "GDF15-TGFBR2",
    grepl("COL1A1.*CD44", x1) ~ "COL1A1-CD44",
    grepl("COL1A2.*CD44", x1) ~ "COL1A2-CD44",
    grepl("COL3A1.*CD44", x1) ~ "COL3A1-CD44",
    TRUE ~ gsub("–|—", "-", x0)
  )
}

clean_side <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("_", " ", x)
  x <- gsub("-", " ", x)
  x <- gsub("/", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

classify_side <- function(side) {
  s <- clean_side(side)
  lineage <- dplyr::case_when(
    grepl("dediff|stromal|dedifferentiation", s) ~ "dediff",
    grepl("myeloid|treg", s) ~ "myeloid",
    TRUE ~ NA_character_
  )
  specificity <- ifelse(grepl("strict|only", s), "only", "high")
  list(lineage = lineage, specificity = specificity)
}

standardize_contrast <- function(x) {
  x0 <- as.character(x)
  vapply(x0, function(z) {
    zl <- tolower(z)
    if (grepl("within", zl) && grepl("both", zl)) {
      return("Within Both-high niche")
    }
    
    z1 <- z
    z1 <- gsub("→|➜|=>|-->|—>|–>", "->", z1)
    z1 <- gsub("_to_|\\bto\\b", "->", z1, ignore.case = TRUE)
    z1 <- gsub("_", " ", z1)
    z1 <- gsub("\\s+", " ", z1)
    z1 <- trimws(z1)
    
    parts <- strsplit(z1, "->", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      src <- paste(parts[1], collapse = " ")
      tgt <- paste(parts[2:length(parts)], collapse = " ")
      src_info <- classify_side(src)
      tgt_info <- classify_side(tgt)
      
      if (!is.na(src_info$lineage) && !is.na(tgt_info$lineage)) {
        if (src_info$lineage == "dediff" && tgt_info$lineage == "myeloid") {
          if (src_info$specificity == "only" || tgt_info$specificity == "only") {
            return("Dediff-only -> Myeloid-only")
          } else {
            return("Dediff/Stromal high -> Myeloid-Treg high")
          }
        }
        if (src_info$lineage == "myeloid" && tgt_info$lineage == "dediff") {
          if (src_info$specificity == "only" || tgt_info$specificity == "only") {
            return("Myeloid-only -> Dediff-only")
          } else {
            return("Myeloid-Treg high -> Dediff/Stromal high")
          }
        }
      }
    }
    
    zc <- clean_side(z)
    pos_d <- regexpr("dediff|stromal|dedifferentiation", zc, perl = TRUE)[1]
    pos_m <- regexpr("myeloid|treg", zc, perl = TRUE)[1]
    has_only <- grepl("strict|only", zc)
    
    if (pos_d > 0 && pos_m > 0) {
      if (pos_d < pos_m) {
        if (has_only) return("Dediff-only -> Myeloid-only")
        return("Dediff/Stromal high -> Myeloid-Treg high")
      } else {
        if (has_only) return("Myeloid-only -> Dediff-only")
        return("Myeloid-Treg high -> Dediff/Stromal high")
      }
    }
    z
  }, character(1))
}

short_contrast <- function(x) {
  dplyr::recode(
    as.character(x),
    "Dediff/Stromal high -> Myeloid-Treg high" = "D/S high -> M/T high",
    "Myeloid-Treg high -> Dediff/Stromal high" = "M/T high -> D/S high",
    "Dediff-only -> Myeloid-only" = "D-only -> M-only",
    "Myeloid-only -> Dediff-only" = "M-only -> D-only",
    "Within Both-high niche" = "Both-high niche",
    .default = as.character(x)
  )
}

standardize_priority <- function(x, fdr, log2e) {
  if (!missing(x) && !all(is.na(x))) {
    xl <- tolower(as.character(x))
    out <- dplyr::case_when(
      grepl("high", xl) ~ "High",
      grepl("medium|moderate|mid", xl) ~ "Medium",
      grepl("low", xl) ~ "Low",
      TRUE ~ NA_character_
    )
    if (any(!is.na(out))) return(out)
  }
  
  dplyr::case_when(
    is.finite(fdr) & fdr < 0.05 & is.finite(log2e) & log2e >= 0.75 ~ "High",
    is.finite(fdr) & fdr < 0.10 & is.finite(log2e) & log2e >= 0.50 ~ "Medium",
    TRUE ~ "Low"
  )
}

theme_s25 <- theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(size = 13.2, face = "bold", color = "black", hjust = 0),
    plot.subtitle = element_text(size = 9.2, color = "black", hjust = 0),
    plot.tag = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 10.5, face = "bold", color = "black"),
    axis.text = element_text(size = 8.7, color = "black"),
    legend.title = element_text(size = 9.8, face = "bold", color = "black"),
    legend.text = element_text(size = 8.7, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

################--------------------------------------------
## 2. Load and harmonize data
################--------------------------------------------

s25_file <- manual_s25_file
if (!nzchar(s25_file) || !file.exists(s25_file)) {
  stop("Spatial validation file does not exist: ", manual_s25_file)
}
s25_file <- normalizePath(s25_file, winslash = "/", mustWork = TRUE)
s25_raw <- read_any_table(s25_file)

axis_col <- find_col(s25_raw, c("^pair_family$", "lr_axis", "L_R_axis", "candidate.*axis", "axis", "ligand.*receptor", "interaction", "pair"), "candidate L-R axis")
contrast_col <- find_col(s25_raw, c("^contrast$", "source.*target.*contrast", "spatial.*contrast", "comparison", "direction", "neighborhood", "niche"), "spatial source-target contrast")
log2_col <- pick_numeric_col(s25_raw, c("log2_spatial_enrichment", "log2_enrichment", "log2FC", "log2_fc", "log2OR", "log2_odds_ratio", "enrichment_log2", "effect_log2", "log2"), "log2 spatial enrichment")
fdr_col <- pick_numeric_col(s25_raw, c("FDR", "fdr", "q_value", "qvalue", "qval", "padj", "adj_p", "adj_pval", "p_adj"), "FDR / adjusted P value")
fraction_col <- pick_numeric_col(s25_raw, c("neighbor_fraction_percent", "neighbor_fraction", "lr_neighbor_fraction", "L_high_R_high_neighbor_fraction", "neighbor_frac", "fraction", "frac", "percent", "pct", "neighbor_percent"), "neighbor fraction", required = FALSE)
priority_col <- find_col(s25_raw, c("validation_priority", "support_priority", "priority", "support", "tier"), "support priority", required = FALSE)

raw_contrast_map <- data.frame(
  raw_contrast = unique(as.character(s25_raw[[contrast_col]])),
  standardized_contrast = standardize_contrast(unique(as.character(s25_raw[[contrast_col]]))),
  stringsAsFactors = FALSE
)

s25_df <- s25_raw %>%
  transmute(
    lr_axis = standardize_lr_axis(.data[[axis_col]]),
    raw_contrast = as.character(.data[[contrast_col]]),
    contrast = standardize_contrast(.data[[contrast_col]]),
    log2_enrichment = as.numeric(.data[[log2_col]]),
    fdr = as.numeric(.data[[fdr_col]]),
    neighbor_fraction_raw = if (!is.na(fraction_col)) as.numeric(.data[[fraction_col]]) else NA_real_,
    support_priority_raw = if (!is.na(priority_col)) as.character(.data[[priority_col]]) else NA_character_
  ) %>%
  mutate(
    neighbor_fraction = dplyr::case_when(
      is.finite(neighbor_fraction_raw) & neighbor_fraction_raw <= 1 ~ neighbor_fraction_raw * 100,
      is.finite(neighbor_fraction_raw) ~ neighbor_fraction_raw,
      TRUE ~ NA_real_
    ),
    support_priority = standardize_priority(support_priority_raw, fdr, log2_enrichment),
    significance = ifelse(is.finite(fdr) & fdr < fdr_cutoff, "FDR < 0.05", "NS")
  ) %>%
  filter(
    lr_axis %in% tracked_axes_order,
    contrast %in% contrast_order,
    is.finite(log2_enrichment),
    is.finite(fdr)
  ) %>%
  group_by(lr_axis, contrast) %>%
  arrange(fdr, desc(log2_enrichment), .by_group = TRUE) %>%
  summarise(
    raw_contrast_example = dplyr::first(raw_contrast),
    log2_enrichment = dplyr::first(log2_enrichment),
    fdr = dplyr::first(fdr),
    neighbor_fraction = dplyr::first(neighbor_fraction),
    support_priority = dplyr::first(support_priority),
    significance = dplyr::first(significance),
    .groups = "drop"
  )

if (nrow(s25_df) == 0) stop("S25 harmonized table has zero rows after filtering.")

if (all(!is.finite(s25_df$neighbor_fraction))) {
  s25_df$neighbor_fraction <- 5
}

s25_df <- s25_df %>%
  mutate(
    neighbor_fraction = pmax(neighbor_fraction, 0),
    neighbor_fraction = pmin(neighbor_fraction, max(20, max(neighbor_fraction, na.rm = TRUE)))
  )

mapped_contrasts <- sort(unique(s25_df$contrast))
missing_contrasts <- setdiff(contrast_order, mapped_contrasts)

################--------------------------------------------
## 3. Panel A and Panel B data
################--------------------------------------------

panelA_grid <- tidyr::expand_grid(
  lr_axis = tracked_axes_order,
  contrast = contrast_order
) %>%
  left_join(s25_df, by = c("lr_axis", "contrast")) %>%
  mutate(
    lr_axis = factor(lr_axis, levels = rev(tracked_axes_order)),
    contrast = factor(contrast, levels = contrast_order),
    significance = ifelse(is.na(significance), "Not tested", significance)
  )

panelA_points <- panelA_grid %>%
  filter(is.finite(log2_enrichment))

panelB_df <- s25_df %>%
  filter(is.finite(log2_enrichment), log2_enrichment > 0, is.finite(fdr), fdr < fdr_cutoff) %>%
  mutate(
    contrast_short = short_contrast(contrast),
    axis_contrast = paste0(lr_axis, " | ", contrast_short),
    support_priority = factor(
      ifelse(support_priority %in% c("High", "Medium", "Low"), support_priority, "Low"),
      levels = c("High", "Medium", "Low")
    )
  ) %>%
  arrange(fdr, desc(log2_enrichment), support_priority)

if (nrow(panelB_df) == 0) {
  panelB_df <- s25_df %>%
    filter(is.finite(log2_enrichment), log2_enrichment > 0) %>%
    mutate(
      contrast_short = short_contrast(contrast),
      axis_contrast = paste0(lr_axis, " | ", contrast_short),
      support_priority = factor(
        ifelse(support_priority %in% c("High", "Medium", "Low"), support_priority, "Low"),
        levels = c("High", "Medium", "Low")
      )
    ) %>%
    arrange(desc(log2_enrichment), fdr)
}

if (is.finite(n_top_panelB)) {
  panelB_df <- utils::head(panelB_df, n_top_panelB)
}

panelB_df <- panelB_df %>%
  arrange(log2_enrichment) %>%
  mutate(axis_contrast = factor(axis_contrast, levels = axis_contrast))

################--------------------------------------------
## 4. Build panels
################--------------------------------------------

abs_lim <- max(abs(panelA_points$log2_enrichment), na.rm = TRUE)
abs_lim <- max(2, ceiling(abs_lim * 10) / 10)

max_neighbor_fraction <- max(panelA_points$neighbor_fraction, na.rm = TRUE)
if (!is.finite(max_neighbor_fraction)) max_neighbor_fraction <- 20
size_limits <- c(0, max(20, max_neighbor_fraction))

pA <- ggplot(panelA_grid, aes(x = contrast, y = lr_axis)) +
  geom_point(
    data = panelA_points,
    aes(color = log2_enrichment, size = neighbor_fraction, shape = significance),
    alpha = 0.86
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_color_gradient2(
    low = "#3B5BA9",
    mid = "white",
    high = "#E64B35",
    midpoint = 0,
    limits = c(-abs_lim, abs_lim),
    name = "log2 spatial\nenrichment"
  ) +
  scale_size_continuous(
    range = c(1.4, 6.0),
    breaks = c(5, 10, 15, 20),
    limits = size_limits,
    name = "L-high/R-high\nneighbor fraction (%)"
  ) +
  scale_shape_manual(
    values = c("FDR < 0.05" = 16, "NS" = 1, "Not tested" = 4),
    breaks = c("FDR < 0.05", "NS"),
    name = NULL
  ) +
  labs(
    title = "Spatial proximity-constrained candidate L-R analysis",
    subtitle = "CosMx neighborhoods; kNN = 20; permutations = 999; candidate-prioritization audit only",
    x = "Spatial source-target contrast",
    y = "Candidate L-R axis"
  ) +
  theme_s25 +
  theme(
    axis.text.x = element_text(size = 8.0, angle = 30, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(size = 8.5, color = "black"),
    legend.position = "right"
  )

priority_colors <- c(
  "High" = "#F8766D",
  "Medium" = "#00BFC4",
  "Low" = "grey70"
)

pB <- ggplot(panelB_df, aes(x = log2_enrichment, y = axis_contrast, fill = support_priority)) +
  geom_col(width = 0.76, color = NA) +
  scale_fill_manual(
    values = priority_colors,
    breaks = c("High", "Medium", "Low"),
    drop = TRUE,
    name = "Support priority"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  labs(
    title = "Spatially enriched candidate L-R axes in CosMx",
    subtitle = "Ranked by FDR and log2 spatial enrichment; proximity support only",
    x = "log2 spatial enrichment",
    y = NULL
  ) +
  theme_s25 +
  theme(
    axis.text.x = element_text(angle = 0, color = "black"),
    axis.text.y = element_text(size = 8.2, color = "black"),
    legend.position = "right"
  )

################--------------------------------------------
## 5. Assemble and save
################--------------------------------------------

supp_s25 <- pA / pB +
  plot_layout(heights = c(0.95, 1.05)) +
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

fig_width <- 10.8
fig_height <- 11.0

safe_ggsave(out_png, supp_s25, width = fig_width, height = fig_height, dpi = dpi_out, device = NULL)
safe_ggsave(out_jpg, supp_s25, width = fig_width, height = fig_height, dpi = dpi_out, device = NULL)
safe_ggsave(out_pdf, supp_s25, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 6. Audit outputs
################--------------------------------------------

audit <- data.frame(
  s25_file = s25_file,
  axis_col = axis_col,
  contrast_col = contrast_col,
  log2_col = log2_col,
  fdr_col = fdr_col,
  fraction_col = fraction_col,
  priority_col = priority_col,
  preferred_knn20_input = grepl("KNN20|KNN20-CosMx|K20", s25_file, ignore.case = TRUE),
  unicode_arrow_removed = TRUE,
  ascii_lr_axis_used = TRUE,
  direction_aware_contrast_parser_used = TRUE,
  all_five_contrasts_kept_in_panelA = TRUE,
  mapped_contrasts = paste(mapped_contrasts, collapse = "; "),
  missing_contrasts = paste(missing_contrasts, collapse = "; "),
  n_harmonized_rows = nrow(s25_df),
  n_panelA_points = nrow(panelA_points),
  n_panelB_rows = nrow(panelB_df),
  n_top_panelB = ifelse(is.finite(n_top_panelB), n_top_panelB, NA),
  fdr_cutoff = fdr_cutoff,
  log2_enrichment_min = min(s25_df$log2_enrichment, na.rm = TRUE),
  log2_enrichment_max = max(s25_df$log2_enrichment, na.rm = TRUE),
  neighbor_fraction_min = min(s25_df$neighbor_fraction, na.rm = TRUE),
  neighbor_fraction_max = max(s25_df$neighbor_fraction, na.rm = TRUE),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S25_CosMx_spatial_proximity_candidate_LR_axes_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

write.csv(raw_contrast_map, file.path(fig_dir, "Supplementary_Figure_S25_raw_to_standardized_contrast_mapping.csv"), row.names = FALSE)
write.csv(s25_df, file.path(fig_dir, "Supplementary_Figure_S25_spatial_proximity_candidate_LR_axes_harmonized.csv"), row.names = FALSE)
write.csv(panelB_df %>% mutate(axis_contrast = as.character(axis_contrast)), file.path(fig_dir, "Supplementary_Figure_S25_spatially_enriched_candidate_LR_axes_ranked.csv"), row.names = FALSE)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S25 script finished successfully.")