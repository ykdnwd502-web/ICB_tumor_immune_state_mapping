############################################################
## Supplementary Figure S12. Composition-aware audit of bulk state scores and response association
##
## Purpose:
##   1) Load harmonized composition-aware correlation and AUC tables with precise file filtering
##   2) Generate and export Supplementary Figure S12 (Panel A composition correlation heatmap and Panel B residualized AUC)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S12. Composition-aware audit of bulk state scores and response association.png
##     - Supplementary Figure S12. Composition-aware audit of bulk state scores and response association.jpg
##     - Supplementary Figure S12. Composition-aware audit of bulk state scores and response association.pdf
############################################################

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

############################################################
## 0. Paths and settings
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

table_root <- file.path(project_dir, "results", "tables")
fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300

manual_corr_file <- ""
manual_auc_file  <- ""

HEATMAP_LOW  <- "#3B82F6"
HEATMAP_MID  <- "#FFFFFF"
HEATMAP_HIGH <- "#EF4444"

############################################################
## 1. Helper functions
############################################################

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

theme_s12 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 12.2, face = "bold", color = "black"),
    axis.text = element_text(size = 10.8, color = "black"),
    legend.title = element_text(size = 10.8, face = "bold", color = "black"),
    legend.text = element_text(size = 9.5, color = "black"),
    strip.background = element_rect(fill = "#D9D9D9", color = "grey50", linewidth = 0.40),
    strip.text = element_text(size = 10.5, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.34),
    panel.grid.minor = element_blank(),
    plot.margin = margin(6, 6, 6, 6)
  )

read_any_table <- function(file) {
  if (!file.exists(file)) stop("File does not exist: ", file)
  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("csv", "txt")) {
    out <- suppressMessages(readr::read_csv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext %in% c("tsv")) {
    out <- suppressMessages(readr::read_tsv(file, show_col_types = FALSE, guess_max = 100000))
  } else if (ext == "rds") {
    out <- readRDS(file)
    if (!is.data.frame(out)) stop("RDS is not a data.frame: ", file)
  } else {
    stop("Unsupported file extension: ", file)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

score_file <- function(files, include_patterns, exclude_patterns = character(0)) {
  f_low <- tolower(files)
  score <- rep(0, length(files))
  for (pat in include_patterns) {
    score <- score + ifelse(grepl(pat, f_low, ignore.case = TRUE), 25, 0)
  }
  for (pat in exclude_patterns) {
    score <- score - ifelse(grepl(pat, f_low, ignore.case = TRUE), 50, 0)
  }
  data.frame(file = files, score = score, stringsAsFactors = FALSE) %>%
    arrange(desc(score), file)
}

auto_find_file <- function(label, include_patterns, exclude_patterns = character(0), required_cols_any = NULL) {
  all_files <- list.files(
    table_root,
    pattern = "\\.(csv|tsv|txt|rds)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(all_files) == 0) stop("No table files found under: ", table_root)
  
  scored <- score_file(all_files, include_patterns, exclude_patterns) %>%
    filter(score > 0)
  
  if (nrow(scored) == 0) {
    stop("Could not find candidate table for ", label, " under: ", table_root)
  }
  
  for (f in scored$file) {
    dat <- tryCatch(read_any_table(f), error = function(e) NULL)
    if (is.null(dat)) next
    
    col_low <- tolower(colnames(dat))
    if (is.null(required_cols_any)) return(f)
    
    ok <- all(vapply(required_cols_any, function(pat) any(grepl(pat, col_low, ignore.case = TRUE)), logical(1)))
    if (ok) return(f)
  }
  
  stop("Candidate files found for ", label, ", but none matched required columns.")
}

find_col <- function(df, patterns, label, required = TRUE) {
  cols <- colnames(df)
  for (pat in patterns) {
    hit <- grep(pat, cols, value = TRUE, ignore.case = TRUE)
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

standardize_state_label <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  x1_low <- tolower(x1)
  
  dplyr::case_when(
    grepl("immune.*defective|defective.*cold|immune.*cold", x1_low) ~ "immune-defective/cold",
    grepl("myeloid.*treg|treg.*immunosuppress", x1_low) ~ "myeloid–Treg immunosuppressive",
    grepl("tumor.*dediff|dediff.*stromal|stromal.*remodel", x1_low) ~ "tumor-dedifferentiation/stromal-remodeling",
    grepl("melanocytic|melanocyte|melanoma.*lineage.*state", x1_low) ~ "melanocytic differentiation",
    TRUE ~ x0
  )
}

standardize_module_label <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("\\s+", " ", x1)
  x1_low <- tolower(x1)
  
  dplyr::case_when(
    grepl("caf|stromal", x1_low) ~ "CAF/stromal",
    grepl("^myeloid$|myeloid", x1_low) ~ "myeloid",
    grepl("t.?nk|t/nk|tnk", x1_low) ~ "T/NK",
    grepl("b.?plasma|b/plasma|plasma", x1_low) ~ "B/plasma",
    grepl("melanoma.*lineage|lineage|melanocytic|melanoma", x1_low) ~ "melanoma lineage",
    grepl("endothelial|endo", x1_low) ~ "endothelial",
    TRUE ~ x0
  )
}

standardize_auc_label <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("\\s+", " ", x1)
  x1_low <- tolower(x1)
  
  dplyr::case_when(
    grepl("^raw|raw score", x1_low) ~ "Raw score",
    grepl("all.*composition|composition.*all", x1_low) ~ "All-composition residualized",
    grepl("caf.*myeloid.*b.?plasma.*melanoma|caf.*myeloid.*b.?plasma.*lineage", x1_low) ~
      "CAF + myeloid + B/plasma + melanoma-lineage residualized",
    grepl("caf.*myeloid.*b.?plasma", x1_low) ~ "CAF + myeloid + B/plasma residualized",
    grepl("caf.*melanoma|caf.*lineage", x1_low) ~ "CAF + melanoma-lineage residualized",
    grepl("caf.*stromal|caf", x1_low) ~ "CAF/stromal residualized",
    TRUE ~ x0
  )
}

state_order <- c(
  "melanocytic differentiation",
  "tumor-dedifferentiation/stromal-remodeling",
  "myeloid–Treg immunosuppressive",
  "immune-defective/cold"
)

module_order <- c(
  "CAF/stromal",
  "myeloid",
  "T/NK",
  "B/plasma",
  "melanoma lineage",
  "endothelial"
)

auc_order <- c(
  "Raw score",
  "CAF + melanoma-lineage residualized",
  "CAF/stromal residualized",
  "CAF + myeloid + B/plasma residualized",
  "CAF + myeloid + B/plasma + melanoma-lineage residualized",
  "All-composition residualized"
)

wrap_state <- function(x) {
  dplyr::recode(
    as.character(x),
    "immune-defective/cold" = "immune-defective/cold",
    "myeloid–Treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
    "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
    "melanocytic differentiation" = "melanocytic\ndifferentiation",
    .default = as.character(x)
  )
}

wrap_module <- function(x) {
  dplyr::recode(
    as.character(x),
    "CAF/stromal" = "CAF/stromal",
    "myeloid" = "myeloid",
    "T/NK" = "T/NK",
    "B/plasma" = "B/plasma",
    "melanoma lineage" = "melanoma\nlineage",
    "endothelial" = "endothelial",
    .default = as.character(x)
  )
}

############################################################
## 2. Locate and load tables (精准匹配相关性与AUC表)
############################################################

corr_file <- manual_corr_file
if (!nzchar(corr_file)) {
  corr_file <- auto_find_file(
    label = "Panel A correlation",
    include_patterns = c("composition", "corr", "spearman", "audit", "state"),
    exclude_patterns = c("auc", "panelc", "paneld", "marker_module", "s8", "s9"),
    required_cols_any = c("state", "module", "rho")
  )
}

auc_file <- manual_auc_file
if (!nzchar(auc_file)) {
  auc_file <- auto_find_file(
    label = "Panel B AUC",
    include_patterns = c("gse78220", "residual", "auc"),
    exclude_patterns = c("loso", "leave", "perm", "cox"),
    required_cols_any = c("auc", "raw")
  )
}

corr_raw <- read_any_table(corr_file)
auc_raw  <- read_any_table(auc_file)

############################################################
## 3. Harmonize data tables
############################################################

dataset_col <- find_col(corr_raw, c("^Dataset$", "^dataset$", "cohort", "source"), "Dataset", required = FALSE)
state_col <- find_col(corr_raw, c("^State$", "state_name", "final.*state", "tumor.*immune.*state", "signature", "score"), "State")
module_col <- find_col(corr_raw, c("composition.*module", "marker.*module", "module", "marker", "cell.*type", "compartment"), "Module")
rho_col <- find_col(corr_raw, c("spearman.*rho", "rho", "spearman", "correlation", "^cor$", "estimate"), "Rho")

panelA_df <- corr_raw %>%
  transmute(
    Dataset = if (!is.na(dataset_col)) as.character(.data[[dataset_col]]) else "GSE244982",
    State = standardize_state_label(.data[[state_col]]),
    Module = standardize_module_label(.data[[module_col]]),
    SpearmanRho = as.numeric(.data[[rho_col]])
  ) %>%
  filter(
    Dataset %in% c("GSE244982", "GSE78220"),
    State %in% state_order,
    Module %in% module_order,
    is.finite(SpearmanRho)
  ) %>%
  group_by(Dataset, State, Module) %>%
  summarise(SpearmanRho = SpearmanRho[1], .groups = "drop") %>%
  mutate(
    Dataset = factor(Dataset, levels = c("GSE244982", "GSE78220")),
    State = factor(State, levels = rev(state_order)),
    Module = factor(Module, levels = module_order),
    label = sprintf("%.2f", SpearmanRho)
  )

auc_label_col <- find_col(auc_raw, c("analysis", "model", "score.*type", "predictor", "residual", "adjust", "method", "label"), "AUC Label")
auc_col <- find_col(auc_raw, c("^auc$", "auc_non", "non.*response.*auc", "auc.*nonresponse", "auc.*non.*response", "AUC_raw", "estimate"), "AUC")
auc_low_col <- find_col(auc_raw, c("ci.*low", "lower", "auc_low", "lo95", "lcl", "conf.low", "q025", "p025", "percentile_2.5", "ci2.5"), "Lower CI", required = FALSE)
auc_high_col <- find_col(auc_raw, c("ci.*high", "upper", "auc_high", "hi95", "ucl", "conf.high", "q975", "p975", "percentile_97.5", "ci97.5"), "Upper CI", required = FALSE)

panelB_df <- auc_raw %>%
  transmute(
    Model = standardize_auc_label(.data[[auc_label_col]]),
    AUC = as.numeric(.data[[auc_col]]),
    CI_low = if (!is.na(auc_low_col)) as.numeric(.data[[auc_low_col]]) else NA_real_,
    CI_high = if (!is.na(auc_high_col)) as.numeric(.data[[auc_high_col]]) else NA_real_
  ) %>%
  filter(Model %in% auc_order, is.finite(AUC)) %>%
  group_by(Model) %>%
  summarise(
    AUC = AUC[1],
    CI_low = CI_low[1],
    CI_high = CI_high[1],
    .groups = "drop"
  ) %>%
  mutate(
    CI_low = ifelse(is.finite(CI_low), CI_low, 0.5),
    CI_high = ifelse(is.finite(CI_high), CI_high, AUC),
    Model = factor(Model, levels = rev(auc_order)),
    label = sprintf("%.3f", AUC)
  )

############################################################
## 4. Plotting
############################################################

pA <- ggplot(panelA_df, aes(x = Module, y = State, fill = SpearmanRho)) +
  geom_tile(color = "white", linewidth = 0.40) +
  geom_text(aes(label = label), size = 3.5, color = "black") +
  facet_grid(Dataset ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = HEATMAP_LOW,
    mid = HEATMAP_MID,
    high = HEATMAP_HIGH,
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman\nrho"
  ) +
  scale_x_discrete(labels = wrap_module) +
  scale_y_discrete(labels = wrap_state) +
  labs(
    title = "State-score associations with marker-based composition modules",
    subtitle = "Spearman correlation; marker modules are used as a composition-aware audit, not formal deconvolution",
    x = "Composition marker module",
    y = "Tumor–immune state"
  ) +
  theme_s12 +
  theme(
    plot.title = element_text(size = 15.5, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10.8, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10.0),
    axis.text.y = element_text(size = 10.2),
    axis.title.x = element_text(size = 12.0, face = "bold"),
    axis.title.y = element_text(size = 12.0, face = "bold"),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.2),
    panel.grid = element_blank()
  )

pB <- ggplot(panelB_df, aes(x = AUC, y = Model)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.45, color = "black") +
  geom_segment(aes(x = CI_low, xend = CI_high, y = Model, yend = Model), linewidth = 0.65, color = "black") +
  geom_point(size = 2.7, color = "black") +
  geom_text(aes(label = label), nudge_x = 0.010, hjust = 0, size = 3.8, color = "black") +
  scale_x_continuous(limits = c(0.43, 0.90), breaks = seq(0.5, 0.8, 0.1), labels = number_format(accuracy = 0.01)) +
  labs(
    title = "GSE78220 raw versus composition-residualized AUC",
    subtitle = "Tumor-dedifferentiation/stromal-remodeling score; descriptive sensitivity analysis",
    x = "AUC for non-response",
    y = NULL
  ) +
  theme_s12 +
  theme(
    plot.title = element_text(size = 15.5, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11.2, hjust = 0.5),
    axis.title.x = element_text(size = 12.5, face = "bold"),
    axis.text.x = element_text(size = 11.0),
    axis.text.y = element_text(size = 10.8)
  )

supp_s12 <- pA / pB +
  plot_layout(heights = c(1.12, 0.88)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 21, face = "bold", family = "sans"),
      plot.margin = margin(4, 5, 4, 5)
    )
  )

############################################################
## 5. Export outputs
############################################################

base_filename <- "Supplementary Figure S12. Composition-aware audit of bulk state scores and response association"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 11.8
fig_height <- 12.2

safe_ggsave(out_png, supp_s12, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s12, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s12, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

audit <- data.frame(
  corr_file = corr_file,
  auc_file = auc_file,
  n_panelA_rows = nrow(panelA_df),
  n_panelB_rows = nrow(panelB_df),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S12_composition_aware_audit_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

cat("\nSupplementary Figure S12 generated successfully in PNG, JPG, and PDF formats (300 DPI)!\n")