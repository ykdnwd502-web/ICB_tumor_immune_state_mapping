################--------------------------------------------
## Supplementary Figure S22. CosMx cell-level dual-high co-localization and threshold sensitivity analysis
##
## Purpose:
##   1) Load CosMx cell-level coordinate, category, and score tables
##   2) Generate and export Supplementary Figure S22 (Panel A dual-high co-localization map, Panel B Both-high distribution map, Panel C threshold sensitivity analysis with visible error bars)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S22. CosMx cell-level dual-high co-localization and threshold sensitivity analysis.png
##     - Supplementary Figure S22. CosMx cell-level dual-high co-localization and threshold sensitivity analysis.jpg
##     - Supplementary Figure S22. CosMx cell-level dual-high co-localization and threshold sensitivity analysis.pdf
################--------------------------------------------

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
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
base_filename <- "Supplementary Figure S22. CosMx cell-level dual-high co-localization and threshold sensitivity analysis"

## Manual overrides. Leave blank to use preferred input or auto-search fallback.
manual_cosmx_file <- ""
manual_sensitivity_file <- ""

preferred_cosmx_file <- file.path(
  project_dir,
  "results", "tables",
  "CosMx_dual_high_QC_sensitivity_with_metadata_v2",
  "CosMx_state_dual_high_with_coordinates_AND_metadata_QC_ENHANCED_v2.csv"
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
      "Could not find a numeric column for ", label, ". Tried:\n",
      paste(preferred, collapse = ", "),
      "\nAvailable columns:\n",
      paste(colnames(df), collapse = ", ")
    )
  }
  
  NA_character_
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  xx <- tolower(trimws(as.character(x)))
  out <- xx %in% c("true", "t", "yes", "y", "1", "high", "positive", "pos")
  out[xx %in% c("false", "f", "no", "n", "0", "low", "negative", "neg", "", "na", "nan")] <- FALSE
  out
}

standardize_dual_category <- function(x) {
  x0 <- as.character(x)
  x1 <- x0
  x1 <- gsub("_", " ", x1)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("/", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("both", xl) ~ "both high",
    grepl("neither|none|other|not high|no high|low", xl) ~ "neither high",
    grepl("myeloid.*treg", xl) & !grepl("dediff|stromal|both", xl) ~ "myeloid–Treg high only",
    grepl("dediff|stromal", xl) & !grepl("myeloid|treg|both", xl) ~ "tumor-dedifferentiation/stromal-remodeling high only",
    grepl("myeloid.*only|treg.*only", xl) ~ "myeloid–Treg high only",
    grepl("dediff.*only|stromal.*only", xl) ~ "tumor-dedifferentiation/stromal-remodeling high only",
    TRUE ~ x0
  )
}

wrap_dual_category <- function(x) {
  dplyr::recode(
    as.character(x),
    "myeloid–Treg high only" = "Myeloid–Treg high only",
    "tumor-dedifferentiation/stromal-remodeling high only" = "Dediff/Stromal high only",
    "both high" = "Both high",
    "neither high" = "Neither high",
    .default = as.character(x)
  )
}

dual_order <- c(
  "neither high",
  "myeloid–Treg high only",
  "tumor-dedifferentiation/stromal-remodeling high only",
  "both high"
)

dual_colors <- c(
  "neither high" = "grey78",
  "myeloid–Treg high only" = "#00A087",
  "tumor-dedifferentiation/stromal-remodeling high only" = "#E64B35",
  "both high" = "#8E44AD"
)

standardize_mode <- function(x) {
  x0 <- as.character(x)
  xl <- tolower(x0)
  
  dplyr::case_when(
    grepl("pooled|global", xl) ~ "Pooled",
    grepl("per.*fov|perfov|per_fov|fov|per.*sample|per_sample|sample", xl) ~ "Per-FOV",
    TRUE ~ x0
  )
}

mode_order <- c("Pooled", "Per-FOV")
mode_colors <- c("Pooled" = "#F8766D", "Per-FOV" = "#00BFC4")

auto_find_cosmx_file <- function(root) {
  files <- list.files(
    file.path(root, "results", "tables"),
    pattern = "\\.(csv|tsv|txt|rds)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(files) == 0) stop("No table files found under results/tables.")
  
  score <- rep(0, length(files))
  f_low <- tolower(files)
  
  score <- score + ifelse(grepl("cosmx", f_low), 80, 0)
  score <- score + ifelse(grepl("dual_high|dual-high|both_high|both-high|state|dominant", f_low), 40, 0)
  score <- score + ifelse(grepl("coordinate|coordinates|metadata|polygon|global|spatial", f_low), 30, 0)
  score <- score + ifelse(grepl("enhanced|qc", f_low), 10, 0)
  score <- score - ifelse(grepl("audit|summary|moran|auc|permutation|null|coverage|sensitivity_or", f_low), 50, 0)
  
  cand <- data.frame(file = files, score = score, stringsAsFactors = FALSE) %>%
    arrange(desc(score), file) %>%
    filter(score > 0)
  
  message("Top candidates for S22 CosMx coordinate/category table:")
  print(utils::head(cand, 20))
  
  for (f in cand$file) {
    x <- tryCatch(read_any_table(f), error = function(e) NULL)
    if (is.null(x)) next
    
    cn <- colnames(x)
    has_x <- any(grepl("^spatial_x$|global.*x|coordinate.*x|coord.*x|^x$", cn, ignore.case = TRUE, perl = TRUE))
    has_y <- any(grepl("^spatial_y$|global.*y|coordinate.*y|coord.*y|^y$", cn, ignore.case = TRUE, perl = TRUE))
    has_cat <- any(grepl("dual.*high.*category|raw.*category|category|both.*high|myeloid.*high|dediff.*high", cn, ignore.case = TRUE, perl = TRUE))
    
    if (has_x && has_y && has_cat) {
      message("Selected S22 CosMx file: ", f)
      return(f)
    }
  }
  
  stop("No suitable CosMx coordinate/category table found. Set manual_cosmx_file.")
}

auto_find_sensitivity_file <- function(root) {
  files <- list.files(
    file.path(root, "results", "tables"),
    pattern = "\\.(csv|tsv|txt|rds)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(files) == 0) return(NA_character_)
  
  score <- rep(0, length(files))
  f_low <- tolower(files)
  
  score <- score + ifelse(grepl("cosmx", f_low), 80, 0)
  score <- score + ifelse(grepl("dual_high|dual-high|both_high|both-high", f_low), 50, 0)
  score <- score + ifelse(grepl("sensitivity|threshold|topfraction|top_fraction|cutoff", f_low), 60, 0)
  score <- score + ifelse(grepl("or|odds|fisher", f_low), 50, 0)
  score <- score - ifelse(grepl("plot_input|coordinate_state|counts|audit|moran|auc|coverage|permutation|null", f_low), 40, 0)
  
  cand <- data.frame(file = files, score = score, stringsAsFactors = FALSE) %>%
    arrange(desc(score), file) %>%
    filter(score > 0)
  
  message("Top candidates for S22 threshold sensitivity OR table:")
  print(utils::head(cand, 20))
  
  for (f in cand$file) {
    x <- tryCatch(read_any_table(f), error = function(e) NULL)
    if (is.null(x)) next
    
    cn <- colnames(x)
    has_frac <- any(grepl("top.*fraction|fraction|quantile|cutoff", cn, ignore.case = TRUE, perl = TRUE))
    has_mode <- any(grepl("mode|cutoff|pooled|fov|sample|stratum", cn, ignore.case = TRUE, perl = TRUE))
    has_or <- any(grepl("^OR$|odds|ratio|fisher|estimate", cn, ignore.case = TRUE, perl = TRUE))
    
    if (has_frac && has_mode && has_or) {
      message("Selected S22 threshold sensitivity table: ", f)
      return(f)
    }
  }
  
  NA_character_
}

theme_s22 <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 14.2, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 9.6, hjust = 0, color = "black"),
    axis.title = element_text(size = 11.0, face = "bold", color = "black"),
    axis.text = element_text(size = 9.0, color = "black"),
    legend.title = element_text(size = 10.0, face = "bold", color = "black"),
    legend.text = element_text(size = 8.7, color = "black"),
    legend.key.height = unit(0.42, "cm"),
    legend.key.width = unit(0.42, "cm"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.32),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 5, 5, 5)
  )

################--------------------------------------------
## 2. Load and harmonize CosMx coordinate/category data
################--------------------------------------------

cosmx_file <- manual_cosmx_file
if (!nzchar(cosmx_file)) {
  if (file.exists(preferred_cosmx_file)) {
    cosmx_file <- preferred_cosmx_file
  } else {
    cosmx_file <- auto_find_cosmx_file(project_dir)
  }
}
cosmx_file <- normalizePath(cosmx_file, winslash = "/", mustWork = TRUE)

cosmx_raw <- read_any_table(cosmx_file)

x_col <- find_col(
  cosmx_raw,
  patterns = c("^spatial_x$", "^global_x$", "^x_global$", "^global\\.x$", "polygon.*x", "coordinate.*x", "coord.*x", "^x$"),
  label = "CosMx global/spatial x"
)

y_col <- find_col(
  cosmx_raw,
  patterns = c("^spatial_y$", "^global_y$", "^y_global$", "^global\\.y$", "polygon.*y", "coordinate.*y", "coord.*y", "^y$"),
  label = "CosMx global/spatial y"
)

category_col <- find_col(
  cosmx_raw,
  patterns = c(
    "^Raw_DualHighCategory$",
    "^Raw_DualHigh_Category$",
    "^DualHighCategory$",
    "^Dual_High_Category$",
    "dual.*high.*category",
    "raw.*high.*category",
    "raw.*category",
    "^Category$",
    "category"
  ),
  label = "CosMx dual-high category",
  required = FALSE
)

myeloid_high_col <- find_col(
  cosmx_raw,
  patterns = c("Raw_MyeloidTreg_high", "MyeloidTreg.*high", "myeloid.*treg.*high", "myeloid.*high"),
  label = "myeloid–Treg high flag",
  required = FALSE
)

dediff_high_col <- find_col(
  cosmx_raw,
  patterns = c("Raw_DediffStromal_high", "DediffStromal.*high", "dediff.*stromal.*high", "dediff.*high", "stromal.*high"),
  label = "dediff/stromal high flag",
  required = FALSE
)

sample_col <- find_col(
  cosmx_raw,
  patterns = c("^fov$", "FOV", "sample", "slide", "roi", "section", "patient"),
  label = "FOV/sample identifier",
  required = FALSE
)

myeloid_score_col <- find_col(
  cosmx_raw,
  patterns = c(
    "^Myeloid_Treg_Immunosuppressive$",
    "^MyeloidTreg_score$",
    "Myeloid.*Treg.*score",
    "myeloid.*treg",
    "Mye.*Treg"
  ),
  label = "myeloid–Treg score",
  required = FALSE,
  exclude_patterns = c("high", "category", "flag")
)

dediff_score_col <- find_col(
  cosmx_raw,
  patterns = c(
    "^Tumor_dedifferentiation_Stromal_remodeling$",
    "^DediffStromal_score$",
    "Dediff.*Stromal.*score",
    "dediff.*stromal",
    "Tumor.*dediff"
  ),
  label = "dediff/stromal score",
  required = FALSE,
  exclude_patterns = c("high", "category", "flag")
)

if (!is.na(category_col)) {
  category_raw <- cosmx_raw[[category_col]]
  category_std <- standardize_dual_category(category_raw)
} else if (!is.na(myeloid_high_col) && !is.na(dediff_high_col)) {
  m_high <- as_bool(cosmx_raw[[myeloid_high_col]])
  d_high <- as_bool(cosmx_raw[[dediff_high_col]])
  
  category_std <- dplyr::case_when(
    m_high & d_high ~ "both high",
    m_high & !d_high ~ "myeloid–Treg high only",
    !m_high & d_high ~ "tumor-dedifferentiation/stromal-remodeling high only",
    TRUE ~ "neither high"
  )
} else {
  stop(
    "Could not derive dual-high category. Need a category column or both high-flag columns.\n",
    "Available columns:\n", paste(colnames(cosmx_raw), collapse = ", ")
  )
}

plot_df <- cosmx_raw %>%
  transmute(
    x = as.numeric(.data[[x_col]]),
    y = as.numeric(.data[[y_col]]),
    dual_category = factor(category_std, levels = dual_order),
    sample_id = if (!is.na(sample_col)) as.character(.data[[sample_col]]) else NA_character_,
    myeloid_score = if (!is.na(myeloid_score_col)) as.numeric(.data[[myeloid_score_col]]) else NA_real_,
    dediff_score = if (!is.na(dediff_score_col)) as.numeric(.data[[dediff_score_col]]) else NA_real_
  ) %>%
  filter(is.finite(x), is.finite(y), !is.na(dual_category))

if (nrow(plot_df) == 0) {
  stop("S22 harmonized CosMx table has zero rows after coordinate/category filtering.")
}

################--------------------------------------------
## 3. Panel C: compute threshold sensitivity with Fisher 95% CI
################--------------------------------------------

compute_or <- function(my_high, de_high) {
  my_high <- as.logical(my_high)
  de_high <- as.logical(de_high)
  
  tab <- table(
    factor(my_high, levels = c(FALSE, TRUE)),
    factor(de_high, levels = c(FALSE, TRUE))
  )
  
  ft <- suppressWarnings(stats::fisher.test(tab))
  
  data.frame(
    OR = as.numeric(ft$estimate),
    low = as.numeric(ft$conf.int[1]),
    high = as.numeric(ft$conf.int[2]),
    stringsAsFactors = FALSE
  )
}

compute_threshold_sensitivity_from_scores <- function(df, fractions = c(0.20, 0.25, 0.30)) {
  if (!all(is.finite(df$myeloid_score)) || !all(is.finite(df$dediff_score))) {
    stop("Score columns are not usable for Panel C direct computation.")
  }
  
  if (all(is.na(df$sample_id))) {
    df$sample_id <- "all"
  }
  
  out <- list()
  
  for (fr in fractions) {
    q_m <- stats::quantile(df$myeloid_score, probs = 1 - fr, na.rm = TRUE, names = FALSE)
    q_d <- stats::quantile(df$dediff_score, probs = 1 - fr, na.rm = TRUE, names = FALSE)
    
    or_pooled <- compute_or(df$myeloid_score >= q_m, df$dediff_score >= q_d) %>%
      mutate(Mode = "Pooled", TopFraction = fr)
    
    tmp <- df %>%
      group_by(sample_id) %>%
      mutate(
        q_m = stats::quantile(myeloid_score, probs = 1 - fr, na.rm = TRUE, names = FALSE),
        q_d = stats::quantile(dediff_score, probs = 1 - fr, na.rm = TRUE, names = FALSE),
        my_high = myeloid_score >= q_m,
        de_high = dediff_score >= q_d
      ) %>%
      ungroup()
    
    or_perfov <- compute_or(tmp$my_high, tmp$de_high) %>%
      mutate(Mode = "Per-FOV", TopFraction = fr)
    
    out[[paste0(fr, "_pooled")]] <- or_pooled
    out[[paste0(fr, "_perfov")]] <- or_perfov
  }
  
  bind_rows(out)
}

use_direct_fisher_ci <- !is.na(myeloid_score_col) && !is.na(dediff_score_col) &&
  all(is.finite(plot_df$myeloid_score)) && all(is.finite(plot_df$dediff_score))

if (use_direct_fisher_ci) {
  sensitivity_file <- "computed_directly_from_CosMx_score_columns"
  sens_df <- compute_threshold_sensitivity_from_scores(plot_df, fractions = c(0.20, 0.25, 0.30))
} else {
  sensitivity_file <- manual_sensitivity_file
  if (!nzchar(sensitivity_file)) {
    sensitivity_file <- auto_find_sensitivity_file(project_dir)
  }
  
  if (is.na(sensitivity_file) || !nzchar(sensitivity_file) || !file.exists(sensitivity_file)) {
    stop("No usable external threshold sensitivity table found.")
  }
  
  sensitivity_file <- normalizePath(sensitivity_file, winslash = "/", mustWork = TRUE)
  sens_raw <- read_any_table(sensitivity_file)
  
  frac_col <- find_col(sens_raw, c("^TopFraction$", "top.*fraction", "fraction", "quantile", "cutoff"), "threshold top fraction")
  mode_col <- find_col(sens_raw, c("^Mode$", "^mode$", "cutoff.*mode", "stratum", "pooled", "fov", "sample"), "cutoff mode")
  or_col <- pick_numeric_col(sens_raw, c("OR", "or", "odds_ratio", "Fisher_OR", "fisher_or", "estimate", "ratio"), "cell-level odds ratio")
  low_col <- pick_numeric_col(sens_raw, c("CI_low", "ci_low", "conf.low", "lower", "low", "lcl"), "OR lower CI")
  high_col <- pick_numeric_col(sens_raw, c("CI_high", "ci_high", "conf.high", "upper", "high", "ucl"), "OR upper CI")
  
  sens_df <- sens_raw %>%
    transmute(
      Mode = standardize_mode(.data[[mode_col]]),
      TopFraction = as.numeric(.data[[frac_col]]),
      OR = as.numeric(.data[[or_col]]),
      low = as.numeric(.data[[low_col]]),
      high = as.numeric(.data[[high_col]])
    ) %>%
    mutate(TopFraction = ifelse(TopFraction > 1.5, TopFraction / 100, TopFraction))
}

sens_df <- sens_df %>%
  filter(
    Mode %in% mode_order,
    is.finite(TopFraction),
    TopFraction %in% c(0.20, 0.25, 0.30),
    is.finite(OR), OR > 0,
    is.finite(low), low > 0,
    is.finite(high), high > 0
  ) %>%
  mutate(
    Mode = factor(Mode, levels = mode_order),
    TopFraction = factor(
      TopFraction,
      levels = c(0.20, 0.25, 0.30),
      labels = c("20%", "25%", "30%")
    )
  ) %>%
  arrange(Mode, TopFraction)

if (nrow(sens_df) == 0) {
  stop("Panel C threshold sensitivity table has zero rows after harmonization.")
}

################--------------------------------------------
## 4. Build panels
################--------------------------------------------

panel_high_df <- plot_df %>% filter(dual_category != "neither high")
panel_bg_df <- plot_df

pA <- ggplot() +
  geom_point(data = panel_bg_df, aes(x = x, y = y), size = 0.055, alpha = 0.08, color = "grey65") +
  geom_point(data = panel_high_df, aes(x = x, y = y, color = dual_category), size = 0.115, alpha = 0.83) +
  coord_equal() +
  scale_y_reverse() +
  scale_color_manual(
    values = dual_colors,
    breaks = c("myeloid–Treg high only", "tumor-dedifferentiation/stromal-remodeling high only", "both high"),
    labels = wrap_dual_category,
    name = NULL,
    drop = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(size = 2.0, alpha = 1), ncol = 1)) +
  labs(
    title = "CosMx cell-level dual-high co-localization map",
    subtitle = "Top 25% pooled cutoff; unit: CosMx cell; high-state categories shown; y-axis inverted for display",
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_s22 +
  theme(
    legend.position = c(0.79, 0.13),
    legend.background = element_rect(fill = scales::alpha("white", 0.65), color = NA),
    legend.text = element_text(size = 8.2)
  )

pB <- ggplot() +
  geom_point(data = panel_bg_df, aes(x = x, y = y), size = 0.055, alpha = 0.10, color = "grey70") +
  geom_point(data = plot_df %>% filter(dual_category == "both high"), aes(x = x, y = y), size = 0.115, alpha = 0.85, color = dual_colors[["both high"]]) +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "CosMx Both-high cell distribution map",
    subtitle = "Only Both-high cells highlighted; faint full-cell background retained for spatial context; y-axis inverted for display",
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_s22 +
  theme(legend.position = "none")

pC <- ggplot(sens_df, aes(x = TopFraction, y = OR, group = Mode, color = Mode)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_line(linewidth = 0.60) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.10, linewidth = 0.75, alpha = 0.95) +
  geom_point(size = 2.45) +
  scale_color_manual(values = mode_colors, breaks = mode_order, name = "Cutoff mode") +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5), expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    title = "CosMx dual-high threshold sensitivity analysis",
    subtitle = "Myeloid–Treg high and Dediff/Stromal high; Fisher odds ratios with 95% CIs",
    x = "Top fraction used to define high-state",
    y = "Cell-level odds ratio"
  ) +
  theme_s22 +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))

################--------------------------------------------
## 5. Assemble and save
################--------------------------------------------

supp_s22 <- pA / pB / pC +
  plot_layout(heights = c(1, 1, 0.78)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold", color = "black"),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 8.8
fig_height <- 13.4

safe_ggsave(out_png, supp_s22, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s22, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s22, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 6. Audit outputs
################--------------------------------------------

category_counts <- as.data.frame(table(plot_df$dual_category), stringsAsFactors = FALSE)
colnames(category_counts) <- c("dual_category", "n_cells")

audit <- data.frame(
  cosmx_file = cosmx_file,
  sensitivity_source = sensitivity_file,
  panelC_direct_fisher_from_scores = use_direct_fisher_ci,
  x_col = x_col,
  y_col = y_col,
  y_axis_inverted_for_display = TRUE,
  category_col = ifelse(is.na(category_col), "derived_from_high_flags", category_col),
  myeloid_high_col = myeloid_high_col,
  dediff_high_col = dediff_high_col,
  myeloid_score_col = myeloid_score_col,
  dediff_score_col = dediff_score_col,
  sample_col = sample_col,
  n_cells = nrow(plot_df),
  x_min = min(plot_df$x, na.rm = TRUE),
  x_max = max(plot_df$x, na.rm = TRUE),
  y_min = min(plot_df$y, na.rm = TRUE),
  y_max = max(plot_df$y, na.rm = TRUE),
  category_counts = paste(category_counts$dual_category, category_counts$n_cells, sep = "=", collapse = "; "),
  n_panelC_rows = nrow(sens_df),
  panelC_error_bars_restored = TRUE,
  panelC_min_ci_width = min(sens_df$high - sens_df$low, na.rm = TRUE),
  panelC_max_ci_width = max(sens_df$high - sens_df$low, na.rm = TRUE),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S22_CosMx_dual_high_colocalization_threshold_sensitivity_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

write.csv(category_counts, file.path(fig_dir, "Supplementary_Figure_S22_CosMx_dual_high_category_counts.csv"), row.names = FALSE)
write.csv(plot_df %>% select(x, y, dual_category, sample_id, myeloid_score, dediff_score), file.path(fig_dir, "Supplementary_Figure_S22_CosMx_coordinate_dual_high_plot_input.csv"), row.names = FALSE)
write.csv(sens_df, file.path(fig_dir, "Supplementary_Figure_S22_threshold_sensitivity_harmonized.csv"), row.names = FALSE)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S22 script finished successfully.")