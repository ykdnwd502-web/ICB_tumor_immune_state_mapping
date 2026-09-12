################--------------------------------------------
## Supplementary Figure S24. CellChat-based nomination of candidate ligand–receptor axes
##
## Purpose:
##   1) Load CellChat all-communications table from CellChat_candidate_LR directory
##   2) Generate and export Supplementary Figure S24 (Panel A focused L-R axes, Panel B tracked axes support, Panel C aggregated interaction counts)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S24. CellChat-based nomination of candidate ligand–receptor axes.png
##     - Supplementary Figure S24. CellChat-based nomination of candidate ligand–receptor axes.jpg
##     - Supplementary Figure S24. CellChat-based nomination of candidate ligand–receptor axes.pdf
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
  library(grid)
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
## 0. Paths and manual overrides
################--------------------------------------------

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

table_root <- file.path(project_dir, "results", "tables")
fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300
base_filename <- "Supplementary Figure S24. CellChat-based nomination of candidate ligand–receptor axes"

## Manual override utilizing the provided file paths.
manual_lr_file <- file.path(
  project_dir,
  "results", "tables",
  "CellChat_candidate_LR",
  "CellChat_all_communications.csv"
)

max_panelA_pairs <- 18
compute_panelC_from_lr_file <- TRUE

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

standardize_group <- function(x) {
  x0 <- as.character(x)
  x1 <- x0
  x1 <- gsub("_", " ", x1)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("b.?plasma|plasma|b cell|b-cell", xl) ~ "B/plasma cells",
    grepl("caf|stromal|fibro", xl) ~ "CAF/stromal-like cells",
    grepl("cycling.*malig|malig.*cycling|cycling", xl) ~ "Cycling malignant",
    grepl("endothel", xl) ~ "Endothelial",
    grepl("malig|tumou?r|cancer", xl) ~ "Malignant",
    grepl("myeloid|macro|mono|dc", xl) ~ "Myeloid cells",
    grepl("t.?nk.?treg|treg", xl) ~ "T/NK/Treg-like cells",
    grepl("t.?nk|nk|t cell|t-cell", xl) ~ "T/NK cells",
    TRUE ~ x0
  )
}

short_group <- function(x) {
  dplyr::recode(
    as.character(x),
    "B/plasma cells" = "B/plasma",
    "CAF/stromal-like cells" = "CAF/stromal",
    "Cycling malignant" = "Cycling malig.",
    "Endothelial" = "Endothelial",
    "Malignant" = "Malignant",
    "Myeloid cells" = "Myeloid",
    "T/NK cells" = "T/NK",
    "T/NK/Treg-like cells" = "T/NK/Treg-like",
    .default = as.character(x)
  )
}

tracked_axes_order <- c(
  "MIF–CD74",
  "MIF–CD44",
  "MIF–CXCR4",
  "SPP1–CD44",
  "SPP1–ITGAV/ITGB1",
  "SPP1–ITGAV/ITGB5",
  "FN1–CD44",
  "APP–CD74",
  "GDF15–TGFBR2",
  "COL1A1–CD44",
  "COL1A2–CD44",
  "COL3A1–CD44"
)

group_order <- c(
  "B/plasma cells",
  "CAF/stromal-like cells",
  "Cycling malignant",
  "Endothelial",
  "Malignant",
  "Myeloid cells",
  "T/NK cells",
  "T/NK/Treg-like cells"
)

detect_tracked_axes <- function(key) {
  k <- toupper(as.character(key))
  k <- gsub("\\s+", "", k)
  k <- gsub("–|—", "-", k)
  k <- gsub("_", "-", k)
  k <- gsub("\\+", "-", k)
  k <- gsub("\\(|\\)|\\[|\\]", "-", k)
  
  out <- character(0)
  has <- function(pattern) grepl(pattern, k, perl = TRUE)
  
  if (has("MIF") && has("CD74")) out <- c(out, "MIF–CD74")
  if (has("MIF") && has("CD44")) out <- c(out, "MIF–CD44")
  if (has("MIF") && has("CXCR4")) out <- c(out, "MIF–CXCR4")
  
  if (has("SPP1") && has("CD44")) out <- c(out, "SPP1–CD44")
  if (has("SPP1") && has("ITGAV") && has("ITGB1")) out <- c(out, "SPP1–ITGAV/ITGB1")
  if (has("SPP1") && has("ITGAV") && has("ITGB5")) out <- c(out, "SPP1–ITGAV/ITGB5")
  
  if (has("FN1") && has("CD44")) out <- c(out, "FN1–CD44")
  if (has("APP") && has("CD74")) out <- c(out, "APP–CD74")
  if (has("GDF15") && has("TGFBR2")) out <- c(out, "GDF15–TGFBR2")
  
  if (has("COL1A1") && has("CD44")) out <- c(out, "COL1A1–CD44")
  if (has("COL1A2") && has("CD44")) out <- c(out, "COL1A2–CD44")
  if (has("COL3A1") && has("CD44")) out <- c(out, "COL3A1–CD44")
  
  unique(out)
}

theme_s24 <- theme_bw(base_size = 10.5) +
  theme(
    plot.title = element_text(size = 13.5, face = "bold", color = "black", hjust = 0),
    plot.subtitle = element_text(size = 9.2, color = "black", hjust = 0),
    plot.tag = element_text(size = 18, face = "bold", color = "black"),
    axis.title = element_text(size = 10.5, face = "bold", color = "black"),
    axis.text = element_text(size = 8.8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.title = element_text(size = 9.8, face = "bold", color = "black"),
    legend.text = element_text(size = 8.7, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

################--------------------------------------------
## 2. Load CellChat data
################--------------------------------------------

lr_file <- manual_lr_file
if (!nzchar(lr_file) || !file.exists(lr_file)) {
  stop("CellChat all communications file does not exist: ", manual_lr_file)
}
lr_file <- normalizePath(lr_file, winslash = "/", mustWork = TRUE)
lr_raw <- read_any_table(lr_file)

source_col <- find_col(lr_raw, patterns = c("^source$", "source.*group", "source.*cell", "sender", "^from$"), label = "source group")
target_col <- find_col(lr_raw, patterns = c("^target$", "target.*group", "target.*cell", "receiver", "^to$"), label = "target group")
prob_col <- pick_numeric_col(lr_raw, preferred = c("prob", "probability", "CellChat_probability", "cellchat_probability", "communication_probability", "score", "weight", "mean"), label = "CellChat probability/score")

ligand_col <- find_col(lr_raw, patterns = c("^ligand$", "ligand"), label = "ligand", required = FALSE)
receptor_col <- find_col(lr_raw, patterns = c("^receptor$", "receptor"), label = "receptor", required = FALSE)
interaction_col <- find_col(lr_raw, patterns = c("interaction_name_2", "interaction_name", "^interaction$", "lr_pair", "pair", "axis", "pathway_name"), label = "interaction / pathway", required = FALSE)
pathway_col <- find_col(lr_raw, patterns = c("pathway_name", "pathway"), label = "pathway", required = FALSE)

lr_base <- lr_raw %>%
  mutate(
    .source_group = standardize_group(.data[[source_col]]),
    .target_group = standardize_group(.data[[target_col]]),
    .prob = as.numeric(.data[[prob_col]]),
    .ligand = if (!is.na(ligand_col)) as.character(.data[[ligand_col]]) else "",
    .receptor = if (!is.na(receptor_col)) as.character(.data[[receptor_col]]) else "",
    .interaction = if (!is.na(interaction_col)) as.character(.data[[interaction_col]]) else "",
    .pathway = if (!is.na(pathway_col)) as.character(.data[[pathway_col]]) else "",
    .axis_key = paste(.ligand, .receptor, .interaction, .pathway, sep = "|"),
    .source_short = short_group(.source_group),
    .target_short = short_group(.target_group),
    .source_target = paste0(.source_short, "\u2192", .target_short)
  ) %>%
  filter(
    is.finite(.prob),
    .source_group %in% group_order,
    .target_group %in% group_order
  )

if (nrow(lr_base) == 0) stop("S24 CellChat base table has zero rows after harmonization.")

axis_list <- lapply(lr_base$.axis_key, detect_tracked_axes)
lr_expanded <- lr_base %>%
  mutate(.row_id = dplyr::row_number(), .lr_axis_list = axis_list) %>%
  tidyr::unnest_longer(.lr_axis_list, values_to = "lr_axis", keep_empty = FALSE) %>%
  filter(lr_axis %in% tracked_axes_order) %>%
  transmute(
    row_id = .row_id,
    source_group = .source_group,
    target_group = .target_group,
    source_short = .source_short,
    target_short = .target_short,
    source_target = .source_target,
    lr_axis = lr_axis,
    prob = .prob
  )

if (nrow(lr_expanded) == 0) stop("No tracked candidate L-R axes were detected after expanding CellChat rows.")

tracked_present <- intersect(unique(lr_expanded$lr_axis), tracked_axes_order)

################--------------------------------------------
## 3. Panel A and Panel B tables
################--------------------------------------------

pair_summary <- lr_expanded %>%
  group_by(source_target) %>%
  summarise(
    max_prob = max(prob, na.rm = TRUE),
    sum_prob = sum(prob, na.rm = TRUE),
    n_axes_positive = dplyr::n_distinct(lr_axis[prob > 0]),
    .groups = "drop"
  ) %>%
  filter(is.finite(max_prob), max_prob > 0) %>%
  arrange(desc(max_prob), desc(n_axes_positive), desc(sum_prob), source_target)

selected_pairs <- utils::head(pair_summary$source_target, max_panelA_pairs)

panelA_df <- lr_expanded %>%
  filter(source_target %in% selected_pairs, prob > 0) %>%
  group_by(source_target, lr_axis) %>%
  summarise(prob = max(prob, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    source_target = factor(source_target, levels = selected_pairs),
    lr_axis = factor(lr_axis, levels = rev(tracked_axes_order))
  ) %>%
  filter(!is.na(lr_axis))

panelB_df <- lr_expanded %>%
  filter(lr_axis %in% tracked_axes_order) %>%
  group_by(lr_axis) %>%
  summarise(max_prob = max(prob, na.rm = TRUE), .groups = "drop")

panelB_df <- data.frame(lr_axis = tracked_axes_order, stringsAsFactors = FALSE) %>%
  left_join(panelB_df, by = "lr_axis") %>%
  mutate(
    max_prob = ifelse(is.finite(max_prob), max_prob, 0),
    support_status = ifelse(max_prob > 0, "CellChat-supported", "Not detected/limited")
  ) %>%
  arrange(max_prob) %>%
  mutate(lr_axis = factor(lr_axis, levels = lr_axis))

################--------------------------------------------
## 4. Panel C counts
################--------------------------------------------

count_df <- lr_base %>%
  filter(.prob > 0) %>%
  group_by(source_group = .source_group, target_group = .target_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::complete(
    source_group = group_order,
    target_group = group_order,
    fill = list(count = 0)
  ) %>%
  mutate(
    source_group = factor(source_group, levels = rev(group_order)),
    target_group = factor(target_group, levels = group_order)
  )

################--------------------------------------------
## 5. Build panels
################--------------------------------------------

prob_max <- max(panelA_df$prob, na.rm = TRUE)
prob_mid <- stats::median(panelA_df$prob[panelA_df$prob > 0], na.rm = TRUE)
if (!is.finite(prob_mid)) prob_mid <- prob_max / 2

pA <- ggplot(panelA_df, aes(x = source_target, y = lr_axis)) +
  geom_point(aes(size = prob, color = prob), alpha = 0.86) +
  scale_size_continuous(
    range = c(1.0, 4.2),
    breaks = pretty(c(0, prob_max), n = 3),
    name = "CellChat-inferred\nscore"
  ) +
  scale_color_gradient2(
    low = "#3B7DDD",
    mid = "white",
    high = "#EF3B2C",
    midpoint = prob_mid,
    name = "CellChat-inferred\nscore"
  ) +
  labs(
    title = "Candidate L–R axes across focused source–target pairs",
    subtitle = "CellChat-inferred scores in focused tumor/stromal–immune directions; candidate nomination only",
    x = "Source \u2192 target",
    y = "Candidate L–R axis"
  ) +
  theme_s24 +
  theme(
    axis.text.x = element_text(size = 7.0, angle = 48, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8.2),
    legend.position = "right"
  )

pB <- ggplot(panelB_df, aes(x = max_prob, y = lr_axis, fill = support_status)) +
  geom_col(width = 0.72, color = NA) +
  geom_point(
    data = panelB_df %>% filter(max_prob == 0),
    aes(x = 0, y = lr_axis),
    inherit.aes = FALSE,
    color = "grey45",
    size = 1.8
  ) +
  scale_fill_manual(
    values = c("CellChat-supported" = "#E64B35", "Not detected/limited" = "grey78"),
    breaks = c("CellChat-supported", "Not detected/limited"),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.00, 0.05))) +
  labs(
    title = "CellChat support for tracked candidate L–R axes",
    subtitle = "Maximum CellChat-inferred score across source–target pairs",
    x = "Maximum CellChat-inferred score",
    y = "Candidate L–R axis"
  ) +
  theme_s24 +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "top",
    legend.justification = "left"
  )

pC <- ggplot(count_df, aes(x = target_group, y = source_group, fill = count)) +
  geom_tile(color = "white", linewidth = 0.50) +
  geom_text(aes(label = ifelse(count > 0, as.character(round(count)), "")), size = 3.0, color = "black") +
  scale_x_discrete(labels = short_group, drop = FALSE) +
  scale_y_discrete(labels = short_group, drop = FALSE) +
  scale_fill_gradient(
    low = "white",
    high = "#EF3B2C",
    name = "Count",
    limits = c(0, max(count_df$count, na.rm = TRUE)),
    oob = scales::squish
  ) +
  labs(
    title = "Aggregated CellChat-inferred interaction counts",
    subtitle = "Number of nonzero CellChat-inferred interactions; candidate-nomination context only",
    x = "Target group",
    y = "Source group"
  ) +
  coord_fixed(ratio = 0.78) +
  theme_s24 +
  theme(
    axis.text.x = element_text(size = 8.0, angle = 42, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8.3),
    legend.position = "right"
  )

################--------------------------------------------
## 6. Assemble and save
################--------------------------------------------

supp_s24 <- pA / pB / pC +
  plot_layout(heights = c(0.96, 0.82, 1.12)) +
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

fig_width <- 9.4
fig_height <- 14.0

safe_ggsave(out_png, supp_s24, width = fig_width, height = fig_height, dpi = dpi_out, device = NULL)
safe_ggsave(out_jpg, supp_s24, width = fig_width, height = fig_height, dpi = dpi_out, device = NULL)
safe_ggsave(out_pdf, supp_s24, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 7. Audit outputs
################--------------------------------------------

audit <- data.frame(
  lr_file = lr_file,
  count_source = "computed_directly_from_selected_CellChat_all_communications_table",
  source_col = source_col,
  target_col = target_col,
  prob_col = prob_col,
  ligand_col = ligand_col,
  receptor_col = receptor_col,
  interaction_col = interaction_col,
  pathway_col = pathway_col,
  n_lr_rows_raw = nrow(lr_raw),
  n_lr_rows_base_after_group_filter = nrow(lr_base),
  n_lr_rows_expanded_to_tracked_axes = nrow(lr_expanded),
  n_panelA_rows = nrow(panelA_df),
  n_panelA_pairs = length(selected_pairs),
  max_panelA_pairs = max_panelA_pairs,
  selected_panelA_pairs = paste(selected_pairs, collapse = "; "),
  tracked_axes_present = paste(sort(tracked_present), collapse = "; "),
  n_panelB_axes = nrow(panelB_df),
  panelB_supported_axes = paste(as.character(panelB_df$lr_axis[panelB_df$max_prob > 0]), collapse = "; "),
  n_panelC_rows = nrow(count_df),
  panelC_count_min = min(count_df$count, na.rm = TRUE),
  panelC_count_max = max(count_df$count, na.rm = TRUE),
  prevented_wrong_count_table_autoselection = TRUE,
  conservative_wording_used = TRUE,
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S24_CellChat_candidate_LR_axes_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

write.csv(
  panelA_df %>%
    mutate(
      source_target = as.character(source_target),
      lr_axis = as.character(lr_axis)
    ),
  file.path(fig_dir, "Supplementary_Figure_S24_panelA_candidate_LR_axes_harmonized.csv"),
  row.names = FALSE
)

write.csv(
  panelB_df %>%
    mutate(lr_axis = as.character(lr_axis)),
  file.path(fig_dir, "Supplementary_Figure_S24_panelB_tracked_candidate_axes_harmonized.csv"),
  row.names = FALSE
)

write.csv(
  count_df %>%
    mutate(
      source_group = as.character(source_group),
      target_group = as.character(target_group)
    ),
  file.path(fig_dir, "Supplementary_Figure_S24_panelC_aggregated_interaction_counts_harmonized.csv"),
  row.names = FALSE
)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S24 script finished successfully.")