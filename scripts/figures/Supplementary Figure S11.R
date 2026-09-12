############################################################
## Supplementary Figure S11. Fragility audit of GSE78220 response association
##
## Purpose:
##   1) Load GSE78220 raw state scores and binary response data via auto-matching
##   2) Dynamically compute leave-one-sample-out (LOSO) AUC influence for all 4 states
##   3) Generate and export Supplementary Figure S11 with accurate per-sample jittered scatter
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S11. Fragility audit of GSE78220 response association.png
##     - Supplementary Figure S11. Fragility audit of GSE78220 response association.jpg
##     - Supplementary Figure S11. Fragility audit of GSE78220 response association.pdf
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
  library(pROC)
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

manual_auc_file  <- ""
manual_perm_file <- ""

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

theme_s11 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 12.5, face = "bold", color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 10.5, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
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
    stop("Unsupported input extension: ", file)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

score_file <- function(files, include_patterns, exclude_patterns = character(0)) {
  f_low <- tolower(files)
  score <- rep(0, length(files))
  for (pat in include_patterns) {
    score <- score + ifelse(grepl(pat, f_low, ignore.case = TRUE), 20, 0)
  }
  for (pat in exclude_patterns) {
    score <- score - ifelse(grepl(pat, f_low, ignore.case = TRUE), 30, 0)
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
  x1 <- gsub("\\s+", " ", x1)
  x1_low <- tolower(x1)
  
  dplyr::case_when(
    grepl("immune.*defective|defective.*cold|immune.*cold", x1_low) ~ "immune-defective/cold",
    grepl("myeloid.*treg|treg.*immunosuppress", x1_low) ~ "myeloid–Treg immunosuppressive",
    grepl("tumor.*dediff|dediff.*stromal|stromal.*remodel", x1_low) ~ "tumor-dedifferentiation/stromal-remodeling",
    grepl("melanocytic|melanocyte", x1_low) ~ "melanocytic differentiation",
    TRUE ~ x0
  )
}

state_order <- c(
  "melanocytic differentiation",
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling"
)

state_order_B <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation"
)

wrap_state <- function(x) {
  dplyr::recode(
    as.character(x),
    "melanocytic differentiation" = "melanocytic\ndifferentiation",
    "immune-defective/cold" = "immune-defective/cold",
    "myeloid–Treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
    "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
    .default = as.character(x)
  )
}

standardize_response_group <- function(x) {
  x0 <- as.character(x)
  x1 <- tolower(x0)
  dplyr::case_when(
    grepl("non.?responder|non.?response|progress|pd|non responder", x1) ~ "Non-responder",
    grepl("responder|response|cr|pr", x1) ~ "Responder",
    TRUE ~ "Non-responder"
  )
}

calc_auc_fixed_direction <- function(response, predictor) {
  dat <- data.frame(
    ResponseGroup = response,
    Predictor = as.numeric(predictor),
    stringsAsFactors = FALSE
  ) %>%
    filter(ResponseGroup %in% c("Responder", "NonResponder"), is.finite(Predictor)) %>%
    mutate(ResponseGroup = factor(ResponseGroup, levels = c("Responder", "NonResponder")))
  
  if (nrow(dat) < 6 || sum(dat$ResponseGroup == "Responder") < 2 || sum(dat$ResponseGroup == "NonResponder") < 2) {
    return(NA_real_)
  }
  
  roc_obj <- pROC::roc(
    response = dat$ResponseGroup,
    predictor = dat$Predictor,
    levels = c("Responder", "NonResponder"),
    direction = "<",
    quiet = TRUE
  )
  as.numeric(pROC::auc(roc_obj))
}

############################################################
## 2. Locate and load tables & Compute LOSO dynamically
############################################################

auc_file <- manual_auc_file
if (!nzchar(auc_file)) {
  auc_file <- auto_find_file(
    label = "Panel A raw AUC summary",
    include_patterns = c("gse78220", "auc", "raw|binary|response|nonresponse|non.response|summary"),
    exclude_patterns = c("loso|leave|perm|permutation|subsampl|bootstrap"),
    required_cols_any = c("auc", "state|signature")
  )
}

perm_file <- manual_perm_file
if (!nzchar(perm_file)) {
  perm_file <- auto_find_file(
    label = "Panel C permutation maximum AUC",
    include_patterns = c("gse78220", "perm|permutation|null", "auc|max"),
    exclude_patterns = c("loso|leave"),
    required_cols_any = c("auc|max")
  )
}

auc_df_raw  <- read_any_table(auc_file)
perm_df_raw <- read_any_table(perm_file)

## 自动通过 auto_find_file 智能匹配任何带有 external_state_scores 的评分表
state_score_file <- auto_find_file(
  label = "GSE78220 external state scores",
  include_patterns = c("gse78220", "external_state_scores"),
  exclude_patterns = c("audit", "inventory"),
  required_cols_any = c("sample", "responsegroup")
)

gse78220_raw <- read_any_table(state_score_file)

raw_cols <- colnames(gse78220_raw)
state_cols_raw <- intersect(c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
), raw_cols)

gse78220_binary <- gse78220_raw %>%
  mutate(ResponseGroup = standardize_response_group(ResponseGroup)) %>%
  filter(ResponseGroup %in% c("Responder", "Non-responder")) %>%
  mutate(ResponseGroupInternal = ifelse(ResponseGroup == "Non-responder", "NonResponder", "Responder"))

loso_list <- list()
for (st in state_cols_raw) {
  dat <- gse78220_binary %>% filter(is.finite(.data[[st]]))
  for (i in seq_len(nrow(dat))) {
    dat_sub <- dat[-i, ]
    auc_i <- calc_auc_fixed_direction(dat_sub$ResponseGroupInternal, dat_sub[[st]])
    
    std_state_name <- standardize_state_label(st)
    loso_list[[length(loso_list) + 1]] <- data.frame(
      State = std_state_name,
      LOSO_AUC = auc_i,
      RemovedGroup = dat$ResponseGroup[i],
      RemovedSample = dat$Sample[i],
      stringsAsFactors = FALSE
    )
  }
}
panelB_df <- bind_rows(loso_list)

############################################################
## 3. Harmonize data tables
############################################################

auc_state_col <- find_col(auc_df_raw, c("^state$", "state_name", "final.*state", "signature", "program"), "State A")
auc_col <- find_col(auc_df_raw, c("^auc$", "^AUC$", "auc.*raw", "raw.*auc", "auc_raw", "AUC_raw_score_for_nonresponse", "auc.*nonresponse", "nonresponse.*auc", "non.*response.*auc", "estimate"), "AUC A")
auc_low_col <- find_col(auc_df_raw, c("AUC_95CI_low", "auc.*95.*ci.*low", "ci.*low", "lower", "auc_low", "lo95", "lcl", "conf.low", "q025", "p025", "percentile_2.5", "ci2.5"), "Lower CI A", required = FALSE)
auc_high_col <- find_col(auc_df_raw, c("AUC_95CI_high", "auc.*95.*ci.*high", "ci.*high", "upper", "auc_high", "hi95", "ucl", "conf.high", "q975", "p975", "percentile_97.5", "ci97.5"), "Upper CI A", required = FALSE)

panelA_df <- auc_df_raw %>%
  transmute(
    State = standardize_state_label(.data[[auc_state_col]]),
    AUC = as.numeric(.data[[auc_col]]),
    CI_low = if (!is.na(auc_low_col)) as.numeric(.data[[auc_low_col]]) else NA_real_,
    CI_high = if (!is.na(auc_high_col)) as.numeric(.data[[auc_high_col]]) else NA_real_
  ) %>%
  filter(State %in% state_order, is.finite(AUC)) %>%
  group_by(State) %>%
  summarise(
    AUC = AUC[1],
    CI_low = CI_low[1],
    CI_high = CI_high[1],
    .groups = "drop"
  ) %>%
  mutate(
    CI_low = ifelse(is.finite(CI_low), CI_low, AUC),
    CI_high = ifelse(is.finite(CI_high), CI_high, AUC),
    State = factor(State, levels = rev(state_order))
  )

panelB_df <- panelB_df %>%
  filter(State %in% state_order_B, is.finite(LOSO_AUC)) %>%
  mutate(
    State = factor(State, levels = state_order_B),
    RemovedGroup = factor(RemovedGroup, levels = c("Non-responder", "Responder"))
  )

perm_auc_col <- find_col(perm_df_raw, c("Max_AUC_across_states", "Observed_max_AUC", "max.*auc", "maximum.*auc", "auc_max", "max_auc", "^auc$"), "Perm AUC")
panelC_df <- perm_df_raw %>%
  transmute(MaxAUC = as.numeric(.data[[perm_auc_col]])) %>%
  filter(is.finite(MaxAUC))

observed_max_auc <- max(panelA_df$AUC, na.rm = TRUE)

############################################################
## 4. Plotting
############################################################

pA <- ggplot(panelA_df, aes(x = AUC, y = State)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.45, color = "black") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.18, linewidth = 0.55, color = "black") +
  geom_point(size = 2.6, color = "black") +
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), labels = number_format(accuracy = 0.01)) +
  scale_y_discrete(labels = wrap_state) +
  labs(
    title = "Raw non-response AUC across predefined states",
    subtitle = "GSE78220 binary response subset; descriptive analysis",
    x = "AUC for non-response",
    y = NULL
  ) +
  theme_s11 +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title.x = element_text(size = 13.5, face = "bold"),
    axis.text.x = element_text(size = 12.2),
    axis.text.y = element_text(size = 12.2)
  )

pB <- ggplot(panelB_df, aes(x = State, y = LOSO_AUC, color = RemovedGroup)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.35, color = "black") +
  geom_boxplot(
    aes(group = State),
    width = 0.20,
    outlier.shape = NA,
    color = "black",
    fill = "white",
    linewidth = 0.35
  ) +
  geom_jitter(width = 0.12, height = 0, size = 1.9, alpha = 0.9) +
  scale_color_manual(
    values = c(
      "Non-responder" = "#F8766D",
      "Responder" = "#00BFC4"
    ),
    name = "Removed group",
    drop = FALSE
  ) +
  scale_x_discrete(labels = wrap_state) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), labels = number_format(accuracy = 0.01)) +
  labs(
    title = "Leave-one-sample-out AUC influence",
    subtitle = "Each point recalculates AUC after removing one GSE78220 sample",
    x = "Final tumor–immune state",
    y = "Leave-one-sample-out AUC"
  ) +
  theme_s11 +
  theme(
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11.0, hjust = 0.5),
    axis.title.x = element_text(size = 11.8, face = "bold"),
    axis.title.y = element_text(size = 11.8, face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 10.2),
    axis.text.y = element_text(size = 10.5),
    legend.position = "right",
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 10)
  )

pC <- ggplot(panelC_df, aes(x = MaxAUC)) +
  geom_histogram(bins = 60, fill = "grey55", color = "grey55", linewidth = 0.15) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.45, color = "black") +
  geom_vline(xintercept = observed_max_auc, linetype = "solid", linewidth = 0.55, color = "black") +
  scale_x_continuous(limits = c(0.35, 1.0), breaks = seq(0.4, 1.0, 0.2), labels = number_format(accuracy = 0.01)) +
  labs(
    title = "Permutation null for maximum AUC across four states",
    subtitle = "Assesses small-cohort state-selection optimism",
    x = "Maximum AUC across four states under permuted labels",
    y = "Frequency"
  ) +
  theme_s11 +
  theme(
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11.0, hjust = 0.5),
    axis.title.x = element_text(size = 11.8, face = "bold"),
    axis.title.y = element_text(size = 11.8, face = "bold"),
    axis.text = element_text(size = 10.5)
  )

supp_s11 <- pA / (pB | pC) +
  plot_layout(heights = c(1.05, 0.95), widths = c(1.0, 1.15)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 21, face = "bold", family = "sans"),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

############################################################
## 5. Export outputs
############################################################

base_filename <- "Supplementary Figure S11. Fragility audit of GSE78220 response association"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 15.0
fig_height <- 10.2

safe_ggsave(out_png, supp_s11, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s11, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s11, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

audit <- data.frame(
  auc_file = auc_file,
  perm_file = perm_file,
  observed_max_auc = observed_max_auc,
  n_panelA_rows = nrow(panelA_df),
  n_panelB_rows = nrow(panelB_df),
  n_panelC_rows = nrow(panelC_df),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S11_Fragility_audit_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

cat("\nSupplementary Figure S11 generated successfully with smart auto-matching in PNG, JPG, and PDF formats (300 DPI)!\n")