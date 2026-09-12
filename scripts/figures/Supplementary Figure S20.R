################--------------------------------------------
## Supplementary Figure S20. External melanoma spatial recurrence audits of composition-linked co-variation
##
## Purpose:
##   1) Load hard-locked S28 external spatial recurrence audit tables
##   2) Generate and export Supplementary Figure S20 (Panel A dataset spots/samples, Panel B gene coverage, Panel C Spearman correlation, Panel D dual-high OR, Panel E bivariate Moran-type statistic)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S20. External melanoma spatial recurrence audits of composition-linked co-variation.png
##     - Supplementary Figure S20. External melanoma spatial recurrence audits of composition-linked co-variation.jpg
##     - Supplementary Figure S20. External melanoma spatial recurrence audits of composition-linked co-variation.pdf
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
## 0. Paths and manual overrides
################--------------------------------------------

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

s28_dir <- file.path(
  project_dir,
  "results", "tables",
  "revision_external_spatial_recurrence_18J",
  "S28_external_spatial_recurrence_audits"
)

fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300
base_filename <- "Supplementary Figure S20. External melanoma spatial recurrence audits of composition-linked co-variation"

manual_panelA_file <- ""
manual_panelB_file <- ""
manual_panelC_file <- ""
manual_panelD_file <- ""
manual_panelE_file <- ""

################--------------------------------------------
## 1. Helper functions
################--------------------------------------------

read_any_table <- function(file) {
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  
  ext <- tolower(tools::file_ext(file))
  
  if (ext %in% c("csv", "txt")) {
    out <- suppressMessages(readr::read_csv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "tsv") {
    out <- suppressMessages(readr::read_tsv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "rds") {
    out <- readRDS(file)
    if (!is.data.frame(out)) {
      stop("RDS object is not a data.frame: ", file)
    }
  } else {
    stop("Unsupported file extension: ", file)
  }
  
  as.data.frame(out, stringsAsFactors = FALSE)
}

first_existing <- function(paths, label) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop(
      "No existing file found for ", label, ". Tried:\n",
      paste(paths, collapse = "\n")
    )
  }
  normalizePath(hit[1], winslash = "/", mustWork = FALSE)
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

find_col <- function(df, patterns, label, required = TRUE, exclude_patterns = character(0)) {
  cols <- colnames(df)
  
  for (pat in patterns) {
    hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE, perl = TRUE)
    
    if (length(hit) > 0 && length(exclude_patterns) > 0) {
      for (ep in exclude_patterns) {
        hit <- hit[!grepl(ep, hit, ignore.case = TRUE, perl = TRUE)]
      }
    }
    
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

pick_numeric_col <- function(df, preferred, label) {
  for (cc in preferred) {
    if (cc %in% colnames(df)) {
      vv <- suppressWarnings(as.numeric(df[[cc]]))
      if (sum(is.finite(vv)) > 0) return(cc)
    }
  }
  stop(
    "Could not find a numeric column for ", label, ". Tried:\n",
    paste(preferred, collapse = ", "),
    "\nAvailable columns:\n",
    paste(colnames(df), collapse = ", ")
  )
}

filter_readout_if_present <- function(df, pattern) {
  if (!"readout" %in% colnames(df)) return(df)
  
  hit <- df %>%
    filter(grepl(pattern, readout, ignore.case = TRUE, perl = TRUE))
  
  if (nrow(hit) > 0) hit else df
}

standardize_dataset <- function(x) {
  x0 <- as.character(x)
  xl <- tolower(x0)
  
  dplyr::case_when(
    grepl("gse250636|250636|visium", xl) ~ "GSE250636 Visium",
    grepl("thrane|legacy", xl) ~ "Thrane 2018 legacy-ST",
    TRUE ~ x0
  )
}

standardize_state <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("/", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("myeloid.*treg|treg.*immunosuppress", xl) ~ "Myeloid–Treg immunosuppressive",
    grepl("tumor.*dediff|dediff.*stromal|stromal.*remodel", xl) ~ "Tumor-dedifferentiation/stromal-remodeling",
    TRUE ~ x0
  )
}

standardize_method <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("^raw$|\\braw\\b", xl) ~ "Raw",
    grepl("stromal|vascular", xl) ~ "Stromal/vascular",
    grepl("broad.*composition|broad.*comp|broad", xl) ~ "Broad composition",
    grepl("tumor.*purity|purity", xl) ~ "Tumor-purity proxy",
    grepl("pc1|pc2|composition pc|comp.*pc|pc", xl) ~ "Composition PC",
    TRUE ~ x0
  )
}

dataset_order <- c("GSE250636 Visium", "Thrane 2018 legacy-ST")
state_order <- c("Myeloid–Treg immunosuppressive", "Tumor-dedifferentiation/stromal-remodeling")
method_order <- c("Raw", "Stromal/vascular", "Broad composition", "Tumor-purity proxy", "Composition PC")

wrap_state <- function(x) {
  dplyr::recode(
    as.character(x),
    "Myeloid–Treg immunosuppressive" = "Myeloid–Treg\nimmunosuppressive",
    "Tumor-dedifferentiation/stromal-remodeling" = "Tumor-dedifferentiation/\nstromal-remodeling",
    .default = as.character(x)
  )
}

wrap_method <- function(x) {
  dplyr::recode(
    as.character(x),
    "Raw" = "Raw",
    "Stromal/vascular" = "Stromal/\nvascular",
    "Broad composition" = "Broad\ncomposition",
    "Tumor-purity proxy" = "Tumor-purity\nproxy",
    "Composition PC" = "Composition PC",
    .default = as.character(x)
  )
}

## 核心修改：在上一次基础上再放大 10%（即总计相比最初版本放大约 21%，系数约为 1.21）
theme_s20 <- theme_bw(base_size = 13.31) +
  theme(
    plot.title = element_text(size = 14.76, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 10.89, hjust = 0, color = "black"),
    plot.tag = element_text(size = 21.78, face = "bold", color = "black"),
    plot.tag.position = c(0, 1.03),
    axis.title = element_text(size = 12.1, face = "bold", color = "black"),
    axis.text = element_text(size = 10.04, color = "black"),
    strip.text = element_text(size = 10.4, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey85", color = "grey45", linewidth = 0.30),
    legend.title = element_text(size = 10.89, face = "bold", color = "black"),
    legend.text = element_text(size = 9.68, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 5, 5, 5)
  )

################--------------------------------------------
## 2. Locate hard-locked S28 input files[cite: 1]
################--------------------------------------------

panelA_file <- if (nzchar(manual_panelA_file)) manual_panelA_file else first_existing(
  c(
    file.path(s28_dir, "S28_plot_input_overview.csv"),
    file.path(s28_dir, "S28_dataset_overview.csv")
  ),
  "Panel A overview"
)

panelB_file <- if (nzchar(manual_panelB_file)) manual_panelB_file else first_existing(
  c(
    file.path(s28_dir, "S28_combined_gene_coverage_summary.csv"),
    file.path(s28_dir, "S28_primary_and_composition_gene_coverage.csv"),
    file.path(s28_dir, "S28_plot_input_coverage.csv")
  ),
  "Panel B coverage"
)

panelC_file <- if (nzchar(manual_panelC_file)) manual_panelC_file else first_existing(
  c(
    file.path(s28_dir, "S28_plot_input_spearman.csv"),
    file.path(s28_dir, "S28_primary_pair_focused_summary.csv"),
    file.path(s28_dir, "S28_combined_key_readout.csv")
  ),
  "Panel C Spearman"
)

panelD_file <- if (nzchar(manual_panelD_file)) manual_panelD_file else first_existing(
  c(
    file.path(s28_dir, "S28_plot_input_dual_high_OR.csv"),
    file.path(s28_dir, "S28_primary_pair_focused_summary.csv"),
    file.path(s28_dir, "S28_combined_key_readout.csv")
  ),
  "Panel D dual-high OR"
)

panelE_file <- if (nzchar(manual_panelE_file)) manual_panelE_file else first_existing(
  c(
    file.path(s28_dir, "S28_plot_input_moran.csv"),
    file.path(s28_dir, "S28_primary_pair_focused_summary.csv"),
    file.path(s28_dir, "S28_combined_key_readout.csv")
  ),
  "Panel E Moran"
)

panelA_raw <- read_any_table(panelA_file)
panelB_raw <- read_any_table(panelB_file)
panelC_raw <- read_any_table(panelC_file)
panelD_raw <- read_any_table(panelD_file)
panelE_raw <- read_any_table(panelE_file)

################--------------------------------------------
## 3. Harmonize Panel A
################--------------------------------------------

ds_col_A <- find_col(panelA_raw, c("dataset_label", "dataset_id", "dataset", "cohort", "source", "platform"), "Panel A dataset")
spots_col_A <- find_col(panelA_raw, c("total_spots", "n_spots", "spot_count", "spots", "nspot", "total"), "Panel A total spots")
samples_col_A <- find_col(panelA_raw, c("n_samples", "sample_count", "samples", "nsample"), "Panel A sample count", required = FALSE)

panelA_df <- panelA_raw %>%
  transmute(
    Dataset = standardize_dataset(.data[[ds_col_A]]),
    n_spots = as.numeric(.data[[spots_col_A]]),
    n_samples = if (!is.na(samples_col_A)) as.numeric(.data[[samples_col_A]]) else NA_real_
  ) %>%
  filter(Dataset %in% dataset_order, is.finite(n_spots)) %>%
  group_by(Dataset) %>%
  summarise(
    n_spots = max(n_spots, na.rm = TRUE),
    n_samples = suppressWarnings(max(n_samples, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_order),
    label = ifelse(
      is.finite(n_samples),
      paste0("n=", comma(n_spots), "\n", n_samples, " samples"),
      paste0("n=", comma(n_spots))
    ),
    label_y = n_spots + 0.035 * max(n_spots, na.rm = TRUE)
  )

if (nrow(panelA_df) == 0) stop("Panel A harmonized table has zero rows.")

################--------------------------------------------
## 4. Harmonize Panel B
################--------------------------------------------

ds_col_B <- find_col(panelB_raw, c("dataset_label", "dataset_id", "dataset", "cohort", "source", "platform"), "Panel B dataset")
state_col_B <- find_col(panelB_raw, c("module_label", "module", "state", "signature", "gene_set"), "Panel B state/signature")

if ("median_coverage" %in% colnames(panelB_raw)) {
  coverage_col_B <- "median_coverage"
  low_col_B <- if ("min_coverage" %in% colnames(panelB_raw)) "min_coverage" else NA_character_
  high_col_B <- if ("max_coverage" %in% colnames(panelB_raw)) "max_coverage" else NA_character_
} else {
  coverage_col_B <- find_col(panelB_raw, c("coverage_fraction", "coverage", "detected.*coverage", "gene.*coverage", "fraction", "proportion", "pct", "percent"), "Panel B coverage")
  low_col_B <- find_col(panelB_raw, c("ci.*low", "low", "lower", "min_coverage", "coverage_low", "q025"), "Panel B lower interval", required = FALSE)
  high_col_B <- find_col(panelB_raw, c("ci.*high", "high", "upper", "max_coverage", "coverage_high", "q975"), "Panel B upper interval", required = FALSE)
}

panelB_df <- panelB_raw %>%
  transmute(
    Dataset = standardize_dataset(.data[[ds_col_B]]),
    State = standardize_state(.data[[state_col_B]]),
    coverage = as.numeric(.data[[coverage_col_B]]),
    low = if (!is.na(low_col_B)) as.numeric(.data[[low_col_B]]) else NA_real_,
    high = if (!is.na(high_col_B)) as.numeric(.data[[high_col_B]]) else NA_real_
  ) %>%
  filter(Dataset %in% dataset_order, State %in% state_order, is.finite(coverage)) %>%
  mutate(
    coverage = ifelse(coverage > 1.5, coverage / 100, coverage),
    low = ifelse(is.finite(low) & low > 1.5, low / 100, low),
    high = ifelse(is.finite(high) & high > 1.5, high / 100, high),
    low = ifelse(is.finite(low), low, coverage),
    high = ifelse(is.finite(high), high, coverage),
    Dataset = factor(Dataset, levels = dataset_order),
    State = factor(State, levels = state_order)
  ) %>%
  group_by(Dataset, State) %>%
  summarise(
    coverage = median(coverage, na.rm = TRUE),
    low = min(low, na.rm = TRUE),
    high = max(high, na.rm = TRUE),
    .groups = "drop"
  )

if (nrow(panelB_df) == 0) stop("Panel B harmonized table has zero rows.")

################--------------------------------------------
## 5. Harmonize Panel C
################--------------------------------------------

panelC_raw2 <- filter_readout_if_present(panelC_raw, "spearman|correlation|rho")

ds_col_C <- find_col(panelC_raw2, c("dataset_label", "dataset_id", "dataset", "cohort", "source", "platform"), "Panel C dataset")
method_col_C <- find_col(panelC_raw2, c("adjustment_label", "adjustment", "method", "analysis", "resid", "score_type", "model"), "Panel C method")
rho_col_C <- pick_numeric_col(panelC_raw2, c("statistic", "spearman_rho", "rho", "correlation", "cor"), "Panel C Spearman rho")

panelC_df <- panelC_raw2 %>%
  transmute(
    Dataset = standardize_dataset(.data[[ds_col_C]]),
    Method = standardize_method(.data[[method_col_C]]),
    rho = as.numeric(.data[[rho_col_C]])
  ) %>%
  filter(Dataset %in% dataset_order, Method %in% method_order, is.finite(rho)) %>%
  group_by(Dataset, Method) %>%
  summarise(rho = median(rho, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_order),
    Method = factor(Method, levels = method_order)
  ) %>%
  arrange(Dataset, Method)

if (nrow(panelC_df) == 0) stop("Panel C harmonized table has zero rows.")

################--------------------------------------------
## 6. Harmonize Panel D
################--------------------------------------------

panelD_raw2 <- filter_readout_if_present(panelD_raw, "dual|high|fisher|odds|or|enrichment")

ds_col_D <- find_col(panelD_raw2, c("dataset_label", "dataset_id", "dataset", "cohort", "source", "platform"), "Panel D dataset")
method_col_D <- find_col(panelD_raw2, c("adjustment_label", "adjustment", "method", "analysis", "resid", "score_type", "model"), "Panel D method")
or_col_D <- pick_numeric_col(panelD_raw2, c("OR", "or", "statistic", "odds_ratio", "fisher_or", "ratio"), "Panel D Fisher odds ratio")
low_col_D <- find_col(panelD_raw2, c("^lower$", "^low$", "ci.*low", "or_low", "conf.low", "lcl", "q025", "CI_low"), "Panel D lower CI", required = FALSE)
high_col_D <- find_col(panelD_raw2, c("^upper$", "^high$", "ci.*high", "or_high", "conf.high", "ucl", "q975", "CI_high"), "Panel D upper CI", required = FALSE)

panelD_df <- panelD_raw2 %>%
  transmute(
    Dataset = standardize_dataset(.data[[ds_col_D]]),
    Method = standardize_method(.data[[method_col_D]]),
    OR = as.numeric(.data[[or_col_D]]),
    low = if (!is.na(low_col_D)) as.numeric(.data[[low_col_D]]) else NA_real_,
    high = if (!is.na(high_col_D)) as.numeric(.data[[high_col_D]]) else NA_real_
  ) %>%
  filter(Dataset %in% dataset_order, Method %in% method_order, is.finite(OR), OR > 0) %>%
  mutate(
    low = ifelse(is.finite(low) & low > 0, low, OR),
    high = ifelse(is.finite(high) & high > 0, high, OR)
  ) %>%
  group_by(Dataset, Method) %>%
  summarise(
    OR = median(OR, na.rm = TRUE),
    low = min(low, na.rm = TRUE),
    high = max(high, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_order),
    Method = factor(Method, levels = method_order)
  ) %>%
  arrange(Dataset, Method)

if (nrow(panelD_df) == 0) stop("Panel D harmonized table has zero rows.")

################--------------------------------------------
## 7. Harmonize Panel E
################--------------------------------------------

panelE_raw2 <- filter_readout_if_present(panelE_raw, "moran|bivariate")

ds_col_E <- find_col(panelE_raw2, c("dataset_label", "dataset_id", "dataset", "cohort", "source", "platform"), "Panel E dataset")
method_col_E <- find_col(panelE_raw2, c("adjustment_label", "adjustment", "method", "analysis", "resid", "score_type", "model"), "Panel E method")
moran_col_E <- pick_numeric_col(
  panelE_raw2,
  c("statistic", "moran", "Moran", "bivariate_moran", "I", "I_symmetric", "Moran_type_I"),
  "Panel E bivariate Moran-type statistic"
)

panelE_df <- panelE_raw2 %>%
  transmute(
    Dataset = standardize_dataset(.data[[ds_col_E]]),
    Method = standardize_method(.data[[method_col_E]]),
    Moran = as.numeric(.data[[moran_col_E]])
  ) %>%
  filter(Dataset %in% dataset_order, Method %in% method_order, is.finite(Moran)) %>%
  group_by(Dataset, Method) %>%
  summarise(Moran = median(Moran, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Dataset = factor(Dataset, levels = dataset_order),
    Method = factor(Method, levels = method_order)
  ) %>%
  arrange(Dataset, Method)

if (nrow(panelE_df) == 0) stop("Panel E harmonized table has zero rows.")

################--------------------------------------------
## 8. Plots with explicit A–E labels[cite: 1]
################--------------------------------------------

pA <- ggplot(panelA_df, aes(x = Dataset, y = n_spots)) +
  geom_col(width = 0.62, fill = "grey35") +
  geom_text(
    aes(y = label_y, label = label),
    vjust = 0,
    size = 3.75,  # 3.41 * 1.10
    lineheight = 0.92
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    tag = "A",
    title = "External spatial datasets",
    x = NULL,
    y = "Total spatial spots"
  ) +
  theme_s20 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

pB <- ggplot(panelB_df, aes(x = State, y = coverage)) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.12, linewidth = 0.45, color = "black") +
  geom_point(size = 2.0, color = "black") +
  facet_wrap(~ Dataset, nrow = 1) +
  scale_x_discrete(labels = wrap_state) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.05),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    tag = "B",
    title = "Primary signature gene coverage",
    x = NULL,
    y = "Detected signature-gene coverage"
  ) +
  theme_s20 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8.83)) # 8.03 * 1.10

pC <- ggplot(panelC_df, aes(x = Method, y = rho, group = 1)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
  geom_line(linewidth = 0.45, color = "black") +
  geom_point(size = 2.0, color = "black") +
  facet_wrap(~ Dataset, nrow = 1) +
  scale_x_discrete(labels = wrap_method) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.10))) +
  labs(
    tag = "C",
    title = "Primary-pair Spearman correlation",
    x = NULL,
    y = "Spearman rho"
  ) +
  theme_s20 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8.71)) # 7.92 * 1.10

pD <- ggplot(panelD_df, aes(x = Method, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.12, linewidth = 0.45, color = "black") +
  geom_point(size = 2.0, color = "black") +
  facet_wrap(~ Dataset, nrow = 1) +
  scale_x_discrete(labels = wrap_method) +
  scale_y_log10(
    breaks = c(1, 3, 5, 10),
    labels = c("1", "3", "5", "10"),
    minor_breaks = NULL,
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  labs(
    tag = "D",
    title = "Primary-pair dual-high enrichment",
    subtitle = "Dashed line indicates odds ratio = 1",
    x = NULL,
    y = "Fisher odds ratio (log10-scaled axis)"
  ) +
  theme_s20 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8.71))

pE <- ggplot(panelE_df, aes(x = Method, y = Moran, group = 1)) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
  geom_line(linewidth = 0.45, color = "black") +
  geom_point(size = 2.0, color = "black") +
  facet_wrap(~ Dataset, nrow = 1) +
  scale_x_discrete(labels = wrap_method) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
  labs(
    tag = "E",
    title = "Descriptive bivariate Moran-type statistic",
    subtitle = "External spatial recurrence audit; positive values indicate spatial co-variation",
    x = NULL,
    y = "Bivariate Moran-type statistic"
  ) +
  theme_s20 +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 8.71))

################--------------------------------------------
## 9. Assemble and save
################--------------------------------------------

header <- ggplot() +
  annotate(
    "text",
    x = 0,
    y = 1,
    hjust = 0,
    vjust = 1,
    label = "External melanoma spatial recurrence audits of composition-linked co-variation",
    fontface = "bold",
    size = 5.32  # 4.84 * 1.10
  ) +
  annotate(
    "text",
    x = 0,
    y = 0.42,
    hjust = 0,
    vjust = 1,
    label = "Primary state pair: myeloid–Treg immunosuppressive vs tumor-dedifferentiation/stromal-remodeling. Analyses are spatial recurrence audits and do not use matched ICB-response annotation.",
    size = 3.75  # 3.41 * 1.10
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(0, 4, 0, 4))

main_grid <- ((pA | pB) / (pC | pD) / pE) +
  plot_layout(heights = c(0.95, 1.05, 1.0))

supp_s20 <- header / main_grid +
  plot_layout(heights = c(0.10, 1.0)) &
  theme(plot.margin = margin(4, 6, 4, 6))

out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 15.8
fig_height <- 12.2

safe_ggsave(out_png, supp_s20, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s20, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s20, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 10. Audit outputs
################--------------------------------------------

audit <- data.frame(
  panelA_file = panelA_file,
  panelB_file = panelB_file,
  panelC_file = panelC_file,
  panelD_file = panelD_file,
  panelE_file = panelE_file,
  panelA_dataset_col = ds_col_A,
  panelA_spots_col = spots_col_A,
  panelA_samples_col = samples_col_A,
  panelB_dataset_col = ds_col_B,
  panelB_state_col = state_col_B,
  panelB_coverage_col = coverage_col_B,
  panelB_low_col = low_col_B,
  panelB_high_col = high_col_B,
  panelC_dataset_col = ds_col_C,
  panelC_method_col = method_col_C,
  panelC_rho_col = rho_col_C,
  panelD_dataset_col = ds_col_D,
  panelD_method_col = method_col_D,
  panelD_or_col = or_col_D,
  panelD_low_col = low_col_D,
  panelD_high_col = high_col_D,
  panelE_dataset_col = ds_col_E,
  panelE_method_col = method_col_E,
  panelE_moran_col = moran_col_E,
  n_panelA_rows = nrow(panelA_df),
  n_panelB_rows = nrow(panelB_df),
  n_panelC_rows = nrow(panelC_df),
  n_panelD_rows = nrow(panelD_df),
  n_panelE_rows = nrow(panelE_df),
  dataset_order = paste(dataset_order, collapse = "; "),
  method_order = paste(method_order, collapse = "; "),
  panel_labels = "A; B; C; D; E",
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S20_external_spatial_recurrence_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

write.csv(panelA_df, file.path(fig_dir, "Supplementary_Figure_S20_panelA_external_dataset_counts_harmonized.csv"), row.names = FALSE)
write.csv(panelB_df, file.path(fig_dir, "Supplementary_Figure_S20_panelB_signature_gene_coverage_harmonized.csv"), row.names = FALSE)
write.csv(panelC_df, file.path(fig_dir, "Supplementary_Figure_S20_panelC_spearman_correlation_harmonized.csv"), row.names = FALSE)
write.csv(panelD_df, file.path(fig_dir, "Supplementary_Figure_S20_panelD_dual_high_enrichment_harmonized.csv"), row.names = FALSE)
write.csv(panelE_df, file.path(fig_dir, "Supplementary_Figure_S20_panelE_bivariate_moran_harmonized.csv"), row.names = FALSE)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S20 script finished successfully.")