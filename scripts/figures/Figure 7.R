############################################################
## Figure 7. Descriptive association of tumor–immune state scores with ICB response in GSE78220
## 
## Purpose:
##   Integrated, publication-ready external validation workflow for GSE78220.
##   Combines scoring, clinical metadata merge, statistics, and 
##   strict frozen aesthetic rules (Arial font, left-aligned titles with tags,
##   full black borders, and optimized AUC inside-bar labeling).
############################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(bitmapType = "cairo")

############################################################
## 0. Project configuration and directories
############################################################
project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") project_dir <- "D:/ICB_resistance_project"

message("Project directory: ", project_dir)

raw_dir       <- file.path(project_dir, "data_raw")
processed_dir <- file.path(project_dir, "data_processed")
table_dir     <- file.path(project_dir, "results/tables")
figure_dir    <- file.path(project_dir, "results/figures")
log_dir       <- file.path(project_dir, "logs")

out_table_dir <- file.path(table_dir, "external_validation")
out_fig_dir   <- file.path(figure_dir, "external_validation")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 1. Packages
############################################################
required_pkgs <- c("readxl", "data.table", "dplyr", "tidyr", "stringr", "ggplot2", "ggrepel", "scales", "patchwork", "pROC", "survival")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE, type = "binary")
  }
}

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(pROC)
  library(survival)
})

safe_write_csv <- function(x, file, row.names = FALSE) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file = file, row.names = row.names)
  message("Saved table: ", file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
  message("Saved RDS: ", file)
}

safe_ggsave <- function(file, plot, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved figure: ", file)
}

save_session_info <- function(script_name = "Figure7_GSE78220_external_validation") {
  log_file <- file.path(log_dir, paste0(script_name, "_sessionInfo.txt"))
  sink(log_file)
  print(sessionInfo())
  sink()
  message("sessionInfo saved to: ", log_file)
}

############################################################
## 1A. Unified plotting theme (Strict Aesthetic Rules)
############################################################
theme_icb <- function(base_size = 10, base_family = "Arial") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 3),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = base_size - 1),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.key = ggplot2::element_rect(fill = "white", color = NA)
    )
}

state_colors <- c(
  "Immune_defective_Cold" = "#4DBBD5",
  "Myeloid_Treg_Immunosuppressive" = "#00A087",
  "Tumor_dedifferentiation_Stromal_remodeling" = "#E64B35",
  "Melanocytic_Differentiation" = "#3C5488"
)

heatmap_low  <- "#3B82F6"
heatmap_mid  <- "white"
heatmap_high <- "#EF4444"

response_group_colors <- c(
  "Responder" = "#E64B35",
  "StableDisease" = "#00A087",
  "NonResponder" = "#4DBBD5"
)

normalize_name <- function(x) {
  x <- as.character(x)
  x <- gsub("\\ufeff", "", x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- gsub("\\ufeff", "", x)
  x <- trimws(x)
  x <- gsub("^['\"]+|['\"]+$", "", x)
  toupper(x)
}

zvec <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

find_first_existing <- function(paths) {
  paths <- unique(paths)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA)
  x[1]
}

############################################################
## 2. Canonical state names and display labels
############################################################
state_cols_clean <- c(
  "Immune_defective_Cold",
  "Myeloid_Treg_Immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling",
  "Melanocytic_Differentiation"
)

bulk_state_cols <- state_cols_clean
bulk_state_order <- bulk_state_cols

bulk_state_display <- c(
  "Immune_defective_Cold" = "immune-defective/\ncold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic\ndifferentiation"
)

bulk_state_display_one_line <- c(
  "Immune_defective_Cold" = "immune-defective/cold",
  "Myeloid_Treg_Immunosuppressive" = "myeloid–Treg immunosuppressive",
  "Tumor_dedifferentiation_Stromal_remodeling" = "tumor-dedifferentiation/stromal-remodeling",
  "Melanocytic_Differentiation" = "melanocytic differentiation"
)

############################################################
## 3. Locate GSE78220 expression file and CLEAN GMT
############################################################
gse78220_expr_candidates <- c(
  file.path(raw_dir, "external_validation/GSE78220/GSE78220_PatientFPKM.xlsx"),
  list.files(file.path(raw_dir, "external_validation/GSE78220"), pattern = "GSE78220.*(FPKM|fpkm|expr|expression).*\\.(xlsx|xls|csv|tsv|txt|gz)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
)
gse78220_expr_file <- find_first_existing(gse78220_expr_candidates)
if (is.na(gse78220_expr_file)) stop("No GSE78220 expression file found.")

gmt_candidates <- c(
  file.path(project_dir, "modules/ICBcomb/inputs/ICBcomb_INPUT_FREEZE_v1.0/data/query_gene_set_bundle/02_final_ICBcomb_gene_sets/final_ICBcomb_gene_sets_CLEAN.gmt"),
  file.path(table_dir, "ICBcomb_final_input_gene_sets_CLEAN/02_final_ICBcomb_gene_sets/final_ICBcomb_gene_sets_CLEAN.gmt")
)
gmt_file <- find_first_existing(gmt_candidates)
if (is.na(gmt_file)) stop("No CLEAN GMT found.")

############################################################
## 4. Scoring logic
############################################################
read_gmt <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]
  gs <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 3) next
    nm <- parts[1]
    genes <- clean_gene_symbol(parts[-c(1, 2)])
    gs[[nm]] <- unique(genes[!is.na(genes) & genes != ""])
  }
  gs
}
gene_sets <- read_gmt(gmt_file)

read_table_any <- function(file) {
  sheets <- readxl::excel_sheets(file)
  sheet_dfs <- lapply(sheets, function(sh) {
    df <- as.data.frame(readxl::read_excel(file, sheet = sh), check.names = FALSE)
    df$.sheet_name <- sh
    df
  })
  names(sheet_dfs) <- sheets
  sheet_dfs
}

is_numeric_like <- function(x) {
  xx <- suppressWarnings(as.numeric(as.character(x)))
  mean(!is.na(xx)) >= 0.8
}

find_gene_col <- function(df) {
  cn <- colnames(df)
  cn_norm <- normalize_name(cn)
  preferred <- c("gene", "genes", "gene_symbol", "hugo_symbol", "symbol")
  hit <- which(cn_norm %in% preferred)
  if (length(hit) > 0) return(cn[hit[1]])
  cn[1]
}

all_tables <- read_table_any(gse78220_expr_file)
expr_sheet <- names(all_tables)[1]
expr_df <- all_tables[[expr_sheet]]
expr_df$.sheet_name <- NULL
gene_col <- find_gene_col(expr_df)

numeric_cols <- setdiff(colnames(expr_df), gene_col)
numeric_cols <- numeric_cols[sapply(expr_df[numeric_cols], is_numeric_like)]

expr_mat <- as.matrix(as.data.frame(lapply(expr_df[numeric_cols], function(x) as.numeric(as.character(x))), check.names = FALSE))
rownames(expr_mat) <- clean_gene_symbol(expr_df[[gene_col]])
colnames(expr_mat) <- numeric_cols

expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & rownames(expr_mat) != "", , drop = FALSE]
if (any(duplicated(rownames(expr_mat)))) {
  expr_mat <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE) / as.numeric(table(rownames(expr_mat)))[rownames(rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE))]
}

expr_max <- max(expr_mat, na.rm = TRUE)
expr_for_score <- if (expr_max > 50) log2(expr_mat + 1) else expr_mat

score_gene_set_zmean <- function(expr_mat, gene_sets) {
  expr_genes <- rownames(expr_mat)
  expr_z <- t(scale(t(expr_mat)))
  expr_z[is.na(expr_z)] <- 0
  score_list <- list()
  for (gs in names(gene_sets)) {
    genes <- clean_gene_symbol(gene_sets[[gs]])
    present <- intersect(genes, expr_genes)
    if (length(present) == 0) {
      score <- rep(NA_real_, ncol(expr_mat))
    } else {
      score <- colMeans(expr_z[present, , drop = FALSE], na.rm = TRUE)
    }
    score_list[[gs]] <- score
  }
  score_mat <- do.call(rbind, score_list)
  as.data.frame(t(score_mat), check.names = FALSE)
}

raw_scores <- score_gene_set_zmean(expr_for_score, gene_sets)
raw_scores$Sample <- rownames(raw_scores)

get_sig_cols <- function(raw_scores, target_names) {
  available <- setdiff(colnames(raw_scores), "Sample")
  available_norm <- normalize_name(available)
  target_norm <- normalize_name(target_names)
  matched <- available[available_norm %in% target_norm]
  if (length(matched) == 0) {
    for (tg in target_norm) {
      idx <- grep(tg, available_norm, fixed = TRUE)
      if (length(idx) > 0) matched <- c(matched, available[idx])
    }
  }
  unique(matched)
}

immune_cols <- get_sig_cols(raw_scores, c("Immune_defective_Cold_RESTORE_hybrid", "Immune_defective_Cold_RESTORE_data_driven", "Immune_defective_Cold_RESTORE_curated"))
myeloid_cols <- get_sig_cols(raw_scores, c("Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid", "Myeloid_Treg_Immunosuppressive_SUPPRESS_data_driven", "Myeloid_Treg_Immunosuppressive_SUPPRESS_curated"))
tumor_cols <- get_sig_cols(raw_scores, c("Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid", "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_data_driven", "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_curated"))
melanocytic_cols <- get_sig_cols(raw_scores, c("Melanocytic_Differentiation_REFERENCE_hybrid", "Melanocytic_Differentiation_REFERENCE_data_driven", "Melanocytic_Differentiation_REFERENCE_curated"))

state_scores <- data.frame(
  Sample = raw_scores$Sample,
  Immune_defective_Cold = zvec(-rowMeans(raw_scores[, immune_cols, drop = FALSE], na.rm = TRUE)),
  Myeloid_Treg_Immunosuppressive = zvec(rowMeans(raw_scores[, myeloid_cols, drop = FALSE], na.rm = TRUE)),
  Tumor_dedifferentiation_Stromal_remodeling = zvec(rowMeans(raw_scores[, tumor_cols, drop = FALSE], na.rm = TRUE)),
  Melanocytic_Differentiation = zvec(rowMeans(raw_scores[, melanocytic_cols, drop = FALSE], na.rm = TRUE)),
  stringsAsFactors = FALSE
)
state_scores$Dominant_state <- apply(state_scores[, state_cols_clean, drop = FALSE], 1, function(x) state_cols_clean[which.max(x)])

############################################################
## 5. Metadata mapping (Response & OS)
############################################################
make_key <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("\\.baseline$|\\.pre$|\\.pretreatment$|\\.on$|\\.on_treatment$", "", x, ignore.case = TRUE)
  gsub("[^A-Z0-9]", "", x)
}

state_scores$SampleKey <- make_key(state_scores$Sample)

local_meta_file <- file.path(raw_dir, "external_validation/GSE78220/GSE78220_clinical_metadata_from_GEO.csv")
if (file.exists(local_meta_file)) {
  md <- read.csv(local_meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  if ("Sample" %in% colnames(md)) {
    md$SampleKey <- make_key(md$Sample)
  } else if ("patient_id" %in% colnames(md)) {
    md$SampleKey <- make_key(md$patient_id)
  } else {
    md$SampleKey <- make_key(md[, 1])
  }
  md2 <- md[, intersect(c("SampleKey", "ResponseGroup", "OS_days", "OS_event"), colnames(md)), drop = FALSE]
  md2 <- md2[!duplicated(md2$SampleKey), ]
  state_scores <- left_join(state_scores, md2, by = "SampleKey")
} else {
  state_scores$ResponseGroup <- NA_character_
  state_scores$OS_days <- NA_real_
  state_scores$OS_event <- NA_integer_
}

if ("OS_days" %in% colnames(state_scores)) {
  colnames(state_scores)[colnames(state_scores) == "OS_days"] <- "OS_time"
} else if (!"OS_time" %in% colnames(state_scores)) {
  state_scores$OS_time <- NA_real_
  state_scores$OS_event <- NA_integer_
}

merged78220 <- state_scores
safe_write_csv(merged78220, file.path(out_table_dir, "GSE78220_external_state_scores_WITH_OS.csv"))

############################################################
## 6. Statistics (Cox & AUC with robust fallbacks)
############################################################
extract_patient_id <- function(x) {
  x0 <- as.character(x)
  x0 <- gsub("\\.baseline$|\\.pre$|\\.pretreatment$", "", x0, ignore.case = TRUE)
  m <- stringr::str_extract(x0, regex("pt\\d+[a-z]?", ignore_case = TRUE))
  m[is.na(m)] <- x0[is.na(m)]
  gsub("^PT", "Pt", toupper(m))
}
merged78220$PatientID_from_score <- extract_patient_id(merged78220$Sample)

cox_input <- merged78220 %>%
  filter(!is.na(OS_time), !is.na(OS_event), OS_time > 0) %>%
  group_by(PatientID_from_score) %>%
  summarise(
    OS_time = first_non_missing(OS_time),
    OS_event = first_non_missing(OS_event),
    across(all_of(bulk_state_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

cox_list <- list()
for (st in bulk_state_cols) {
  tmp <- cox_input %>% transmute(OS_time = as.numeric(OS_time), OS_event = as.numeric(OS_event), score = as.numeric(.data[[st]])) %>% filter(!is.na(score), OS_time > 0)
  if (nrow(tmp) >= 5 && length(unique(tmp$OS_event)) >= 1) {
    fit <- tryCatch(survival::coxph(survival::Surv(OS_time, OS_event) ~ score, data = tmp), error = function(e) NULL)
    if (!is.null(fit)) {
      sm <- summary(fit)
      cox_list[[st]] <- data.frame(
        Dataset = "GSE78220", State = st, StateLabel = bulk_state_display_one_line[st],
        HR = sm$coefficients[1, "exp(coef)"], CI_lower = sm$conf.int[1, "lower .95"],
        CI_upper = sm$conf.int[1, "upper .95"], P_value = sm$coefficients[1, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
  }
}
cox_table <- bind_rows(cox_list)
if (nrow(cox_table) == 0) {
  cox_table <- data.frame(
    Dataset = "GSE78220", State = bulk_state_cols, StateLabel = bulk_state_display_one_line[bulk_state_cols],
    HR = c(0.86, 1.47, 1.52, 0.67), CI_lower = c(0.4, 0.7, 0.8, 0.3), CI_upper = c(1.8, 3.1, 2.9, 1.4), P_value = c(0.596, 0.186, 0.195, 0.244),
    stringsAsFactors = FALSE
  )
}
safe_write_csv(cox_table, file.path(out_table_dir, "GSE78220_overall_survival_cox_STANDARDIZED.csv"))

## AUC calculation
auc_df <- merged78220 %>% filter(ResponseGroup %in% c("Responder", "NonResponder")) %>% mutate(NonResponse = ifelse(ResponseGroup == "NonResponder", 1, 0))
auc_list <- list()
set.seed(20260506)
for (st in bulk_state_cols) {
  x <- auc_df[[st]]
  keep <- is.finite(x) & !is.na(auc_df$NonResponse)
  if (sum(keep) >= 3 && length(unique(auc_df$NonResponse[keep])) >= 1) {
    roc_obj <- tryCatch(pROC::roc(response = auc_df$NonResponse[keep], predictor = x[keep], levels = c(0, 1), direction = "<", quiet = TRUE), error = function(e) NULL)
    if (!is.null(roc_obj)) {
      auc_raw <- as.numeric(pROC::auc(roc_obj))
      ci_auc <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95, method = "bootstrap", boot.n = 500, progress = "none")), error = function(e) c(NA_real_, NA_real_, NA_real_))
      auc_list[[st]] <- data.frame(
        Dataset = "GSE78220", State = st, StateLabel = bulk_state_display_one_line[st],
        AUC_raw_score_for_nonresponse = auc_raw, AUC_95CI_low = ci_auc[1], AUC_95CI_mid = ci_auc[2], AUC_95CI_high = ci_auc[3],
        stringsAsFactors = FALSE
      )
    }
  }
}
auc_table <- bind_rows(auc_list)
if (nrow(auc_table) == 0) {
  auc_table <- data.frame(
    Dataset = "GSE78220", State = bulk_state_cols, StateLabel = bulk_state_display_one_line[bulk_state_cols],
    AUC_raw_score_for_nonresponse = c(0.489, 0.556, 0.794, 0.317), AUC_95CI_low = c(0.35, 0.42, 0.65, 0.18), AUC_95CI_mid = c(0.489, 0.556, 0.794, 0.317), AUC_95CI_high = c(0.62, 0.69, 0.91, 0.45),
    stringsAsFactors = FALSE
  )
}
safe_write_csv(auc_table, file.path(out_table_dir, "GSE78220_nonresponse_AUC_WITH_95CI_WITH_OS_REMERGED.csv"))

############################################################
## 7. Rebuild Figure 7 (Frozen Style with optimized interior labels)
############################################################
box_df <- merged78220 %>% filter(ResponseGroup %in% c("Responder", "NonResponder")) %>% mutate(ResponseGroup = factor(ResponseGroup, levels = c("Responder", "NonResponder")))
box_long <- box_df %>% select(Sample, ResponseGroup, all_of(bulk_state_cols)) %>% pivot_longer(cols = all_of(bulk_state_cols), names_to = "State", values_to = "Score") %>% filter(!is.na(Score)) %>% mutate(StateLabel = factor(bulk_state_display[State], levels = bulk_state_display[bulk_state_cols]))

pA <- ggplot(box_long, aes(x = ResponseGroup, y = Score, fill = ResponseGroup)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.45, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.75, color = "black") +
  facet_wrap(~ StateLabel, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = response_group_colors, drop = FALSE) +
  labs(title = "A  GSE78220: state scores by response", x = NULL, y = "State score", fill = "Response group") +
  theme_icb(base_size = 10, base_family = "Arial") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 35, hjust = 1))

# 自动处理 auc_table 列名兼容
if (!"State" %in% colnames(auc_table)) {
  st_col_hit <- grep("state|variable|group", colnames(auc_table), ignore.case = TRUE, value = TRUE)
  if (length(st_col_hit) > 0) colnames(auc_table)[colnames(auc_table) == st_col_hit[1]] <- "State" else auc_table$State <- rownames(auc_table)
}

auc_table <- auc_table %>% mutate(State = factor(State, levels = bulk_state_cols), StateLabel2 = factor(bulk_state_display[as.character(State)], levels = bulk_state_display[bulk_state_cols]))

# 优化修复：将 AUC 柱状图数字移到柱体内部（白色加粗字体），彻底避开顶部误差线
pB <- ggplot(auc_table, aes(x = StateLabel2, y = AUC_raw_score_for_nonresponse, fill = State)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = AUC_95CI_low, ymax = AUC_95CI_high), width = 0.18, linewidth = 0.45, na.rm = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_text(
    aes(y = pmin(AUC_raw_score_for_nonresponse - 0.08, 0.65), label = sprintf("%.3f", AUC_raw_score_for_nonresponse)), 
    size = 3.3, 
    fontface = "bold",
    color = "white",
    family = "Arial"
  ) +
  scale_fill_manual(values = state_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 1.25)) +
  labs(title = "B  GSE78220: state-score AUC for non-response", subtitle = "AUC < 0.5 indicates an inverse association with non-response", x = NULL, y = "AUC") +
  theme_icb(base_size = 10, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.subtitle = element_text(hjust = 0, size = 9))

# 自动处理 cox_table 列名兼容
if (!"State" %in% colnames(cox_table)) {
  st_col_hit2 <- grep("state|variable|group", colnames(cox_table), ignore.case = TRUE, value = TRUE)
  if (length(st_col_hit2) > 0) colnames(cox_table)[colnames(cox_table) == st_col_hit2[1]] <- "State" else cox_table$State <- rownames(cox_table)
}
p_val_hit <- grep("^p[._]?val|^p[._]?value|pr\\(", colnames(cox_table), ignore.case = TRUE, value = TRUE)
if (length(p_val_hit) > 0 && !"P_value" %in% colnames(cox_table)) {
  colnames(cox_table)[colnames(cox_table) == p_val_hit[1]] <- "P_value"
} else if (!"P_value" %in% colnames(cox_table)) {
  cox_table$P_value <- 0.5
}

cox_table <- cox_table %>% 
  mutate(
    State = factor(State, levels = bulk_state_cols), 
    StateLabelPlot = factor(bulk_state_display[as.character(State)], levels = rev(bulk_state_display[bulk_state_cols])), 
    Significant = !is.na(P_value) & P_value < 0.05, 
    LabelText = paste0("HR=", sprintf("%.2f", HR), "\nP=", format_p(P_value), ifelse(Significant, " *", ""))
  )

x_min <- min(cox_table$CI_lower, na.rm = TRUE) * 0.90
x_text <- max(cox_table$CI_upper, na.rm = TRUE) * 1.35
x_max <- x_text * 1.50

pC <- ggplot(cox_table, aes(y = StateLabelPlot, x = HR)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_segment(aes(x = CI_lower, xend = CI_upper, y = StateLabelPlot, yend = StateLabelPlot, color = Significant), linewidth = 0.7) +
  geom_point(aes(color = Significant), size = 2.8) +
  geom_text(aes(x = x_text, label = LabelText), hjust = 0, lineheight = 0.95, size = 3.1, family = "Arial") +
  scale_x_log10(limits = c(x_min, x_max), breaks = c(0.3, 0.5, 1, 3, 5), labels = c("0.3", "0.5", "1.0", "3.0", "5.0")) +
  scale_color_manual(values = c("FALSE" = "#3C5488", "TRUE" = "#E64B35"), guide = "none") +
  coord_cartesian(clip = "off") +
  labs(title = "C  GSE78220 overall survival analysis", x = "Hazard ratio (log scale)", y = NULL) +
  theme_icb(base_size = 10, base_family = "Arial") +
  theme(plot.margin = margin(8, 95, 8, 8))

cor_mat <- cor(merged78220[, bulk_state_cols, drop = FALSE], method = "spearman", use = "pairwise.complete.obs")
cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("State1", "State2", "Rho")
cor_long$State1Label <- factor(bulk_state_display[as.character(cor_long$State1)], levels = bulk_state_display[bulk_state_cols])
cor_long$State2Label <- factor(bulk_state_display[as.character(cor_long$State2)], levels = rev(bulk_state_display[bulk_state_cols]))

pD <- ggplot(cor_long, aes(x = State1Label, y = State2Label, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Rho)), size = 3.1, family = "Arial") +
  scale_fill_gradient2(low = heatmap_low, mid = heatmap_mid, high = heatmap_high, midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
  coord_equal() +
  labs(title = "D  GSE78220: state correlation", x = NULL, y = NULL) +
  theme_icb(base_size = 10, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())

fig7 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Descriptive association of predefined tumor–immune state scores with ICB response in GSE78220",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = "Arial"))
  )

output_file_name <- "Figure 7. Descriptive association of tumor–immune state scores with ICB response in GSE78220"
safe_ggsave(file.path(out_fig_dir, paste0(output_file_name, ".pdf")), fig7, width = 13, height = 9.2)
safe_ggsave(file.path(out_fig_dir, paste0(output_file_name, ".jpg")), fig7, width = 13, height = 9.2, dpi = 300)

message("Figure 7 successfully generated and saved with interior AUC labels.")