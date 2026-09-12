################--------------------------------------------
## Supplementary Figure S19. Humoral and TLS-like spatial context in Visium dual-high regions
##
## Purpose:
##   1) Load spatial spot-level humoral/TLS module scores and coordinate/category context tables
##   2) Generate and export Supplementary Figure S19 (Panel A module score boxplots, Panel B spatial score distribution maps, Panel C spatial autocorrelation curves, Panel D operational co-high spatial features)
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S19. Humoral and TLS-like spatial context in Visium dual-high regions.png
##     - Supplementary Figure S19. Humoral and TLS-like spatial context in Visium dual-high regions.jpg
##     - Supplementary Figure S19. Humoral and TLS-like spatial context in Visium dual-high regions.pdf
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

table_root <- file.path(project_dir, "results", "tables")
fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300
high_quantile <- 0.75
knn_values <- c(4, 6, 8, 12)

manual_tls_module_file <- ""
manual_spot_context_file <- file.path(
  project_dir,
  "results", "tables", "revision_visium_composition",
  "18E_Visium_spot_state_composition_residualized_scores.csv"
)

base_filename <- "Supplementary Figure S19. Humoral and TLS-like spatial context in Visium dual-high regions"

################--------------------------------------------
## 1. Helpers
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

score_file <- function(files, include_patterns, prefer_patterns = character(0), exclude_patterns = character(0)) {
  f_low <- tolower(files)
  score <- rep(0, length(files))
  
  for (pat in include_patterns) {
    score <- score + ifelse(grepl(pat, f_low, ignore.case = TRUE, perl = TRUE), 20, 0)
  }
  for (pat in prefer_patterns) {
    score <- score + ifelse(grepl(pat, f_low, ignore.case = TRUE, perl = TRUE), 60, 0)
  }
  for (pat in exclude_patterns) {
    score <- score - ifelse(grepl(pat, f_low, ignore.case = TRUE, perl = TRUE), 100, 0)
  }
  
  data.frame(file = files, score = score, stringsAsFactors = FALSE) %>%
    arrange(desc(score), file)
}

auto_find_tls_module_file <- function(root) {
  files <- list.files(
    root,
    pattern = "\\.(csv|tsv|txt|rds)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(files) == 0) stop("No table files found under: ", root)
  
  scored <- score_file(
    files,
    include_patterns = c("tls|humoral|plasma|b_cell|t_cell|tfh|chemokine|fdc|hev|visium|spatial"),
    prefer_patterns = c("tls", "humoral", "module", "visium", "spatial", "18f", "sr", "figure_s19"),
    exclude_patterns = c("auc|loso|permutation|null|external|gse78220|gse91061|cosmx|inventory|audit")
  ) %>%
    filter(score > 0)
  
  if (nrow(scored) == 0) stop("No candidate TLS/humoral module tables found.")
  
  required_any <- c(
    "b.*cell|b_cell|bcell",
    "plasma",
    "t.*cell.*zone|t_cell_zone",
    "tfh",
    "fdc|hev",
    "ig.*humoral|humoral",
    "tls.*chemokine|chemokine"
  )
  
  for (f in scored$file) {
    dat <- tryCatch(read_any_table(f), error = function(e) NULL)
    if (is.null(dat)) next
    cn <- colnames(dat)
    
    ok <- all(vapply(required_any, function(pat) {
      any(grepl(pat, cn, ignore.case = TRUE, perl = TRUE))
    }, logical(1)))
    
    if (ok) {
      message("Selected S19 TLS/humoral module table: ", f)
      return(f)
    }
  }
  
  stop(
    "Candidate files were found, but none contained all seven TLS/humoral module columns.\n",
    "Set manual_tls_module_file to the correct table path."
  )
}

standardize_category <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("_", " ", x0)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("/", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)
  
  dplyr::case_when(
    grepl("both", xl) ~ "both high",
    grepl("neither|none|no high|low", xl) ~ "neither high",
    grepl("myeloid.*treg", xl) & !grepl("dediff|stromal|both", xl) ~ "myeloid–Treg high only",
    grepl("dediff|stromal", xl) & !grepl("myeloid|treg|both", xl) ~ "tumor-dedifferentiation/stromal-remodeling high only",
    grepl("myeloid.*only|treg.*only", xl) ~ "myeloid–Treg high only",
    grepl("dediff.*only|stromal.*only", xl) ~ "tumor-dedifferentiation/stromal-remodeling high only",
    TRUE ~ x0
  )
}

## 核心修改：将过长的分类标签采用精炼缩写与折行显示，防止 Panel A 发生文本互相挤压重叠
wrap_category <- function(x) {
  dplyr::recode(
    as.character(x),
    "neither high" = "neither high",
    "myeloid–Treg high only" = "myeloid–Treg\nhigh only",
    "tumor-dedifferentiation/stromal-remodeling high only" = "dediff/stromal\nhigh only",
    "both high" = "both high",
    .default = as.character(x)
  )
}

category_order <- c(
  "neither high",
  "myeloid–Treg high only",
  "tumor-dedifferentiation/stromal-remodeling high only",
  "both high"
)

module_display_order <- c(
  "B-cell",
  "FDC/HEV-like",
  "Ig-humoral",
  "plasma cell",
  "T-cell zone",
  "Tfh-like",
  "TLS chemokine core"
)

module_col_patterns <- list(
  "B-cell" = c("^B_Cell$", "^Bcell$", "B.*Cell"),
  "FDC/HEV-like" = c("^FDC_HEV_Like$", "FDC.*HEV", "FDC", "HEV"),
  "Ig-humoral" = c("^IG_Humoral$", "^Ig_Humoral$", "IG.*Humoral", "Ig.*Humoral", "Humoral"),
  "plasma cell" = c("^Plasma_Cell$", "Plasma.*Cell", "Plasma"),
  "T-cell zone" = c("^T_Cell_Zone$", "T.*Cell.*Zone", "T_cell_zone"),
  "Tfh-like" = c("^Tfh_Like$", "Tfh.*Like", "Tfh"),
  "TLS chemokine core" = c("^TLS_Chemokine_Core$", "TLS.*Chemokine", "Chemokine.*Core")
)

feature_order <- c(
  "other",
  "Ig-humoral high/other",
  "TLS chemokine high/other",
  "Ig-humoral/T-cell co-high",
  "TLS-chemokine/B-cell/T-cell co-high"
)

feature_colors <- c(
  "other" = "grey82",
  "Ig-humoral high/other" = "#E64B35",
  "TLS chemokine high/other" = "#00A087",
  "Ig-humoral/T-cell co-high" = "#4DBBD5",
  "TLS-chemokine/B-cell/T-cell co-high" = "#B26BFF"
)

wrap_feature <- function(x) {
  dplyr::recode(
    as.character(x),
    "other" = "other",
    "Ig-humoral high/other" = "Ig-humoral\nhigh/other",
    "TLS chemokine high/other" = "TLS chemokine\nhigh/other",
    "Ig-humoral/T-cell co-high" = "Ig-humoral/\nT-cell co-high",
    "TLS-chemokine/B-cell/T-cell co-high" = "TLS-chemokine/\nB-cell/T-cell co-high",
    .default = as.character(x)
  )
}

theme_s19 <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 12.8, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 9.1, hjust = 0, color = "black"),
    axis.title = element_text(size = 10.4, face = "bold", color = "black"),
    axis.text = element_text(size = 8.6, color = "black"),
    strip.text = element_text(size = 8.2, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey85", color = "grey45", linewidth = 0.30),
    legend.title = element_text(size = 9.0, face = "bold", color = "black"),
    legend.text = element_text(size = 7.9, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

################--------------------------------------------
## 2. Load and harmonize input data
################--------------------------------------------

tls_file <- manual_tls_module_file
if (!nzchar(tls_file)) {
  tls_file <- auto_find_tls_module_file(table_root)
}

tls_raw <- read_any_table(tls_file)
context_raw <- read_any_table(manual_spot_context_file)

module_cols <- character(0)
for (nm in names(module_col_patterns)) {
  module_cols[nm] <- find_col(
    tls_raw,
    patterns = module_col_patterns[[nm]],
    label = paste0("module score: ", nm),
    required = TRUE
  )
}

x_col <- find_col(
  tls_raw,
  patterns = c("^spatial_x$", "^imagecol$", "^pxl_col_in_fullres$", "^x$", "pixel.*x", "coord.*x", "array_col", "col"),
  label = "spatial x",
  required = FALSE
)

y_col <- find_col(
  tls_raw,
  patterns = c("^spatial_y$", "^imagerow$", "^pxl_row_in_fullres$", "^y$", "pixel.*y", "coord.*y", "array_row", "row"),
  label = "spatial y",
  required = FALSE
)

cat_col <- find_col(
  tls_raw,
  patterns = c("Raw.*Dual.*High.*Category", "Raw.*High.*Category", "Raw.*Category", "DualHighCategory", "Dual_High_Category", "High.*Category", "^Category$"),
  label = "raw high-state category",
  required = FALSE
)

need_context <- is.na(x_col) || is.na(y_col) || is.na(cat_col)
tls_joined <- tls_raw

join_key <- NA_character_
if (need_context) {
  possible_keys <- c("Spot", "barcode", "Barcode", "spot", "SpatialBarcode", "SpatialSampleID")
  key_hits <- possible_keys[possible_keys %in% colnames(tls_raw) & possible_keys %in% colnames(context_raw)]
  
  if (length(key_hits) > 0) {
    join_key <- key_hits[1]
    context_keep_cols <- unique(c(
      join_key,
      "spatial_x",
      "spatial_y",
      "imagecol",
      "imagerow",
      "Raw_DualHighCategory",
      "Raw_DualHigh_Category",
      "DualHighCategory",
      "Category"
    ))
    context_keep_cols <- context_keep_cols[context_keep_cols %in% colnames(context_raw)]
    
    tls_joined <- tls_raw %>%
      left_join(
        context_raw[, context_keep_cols, drop = FALSE],
        by = join_key,
        suffix = c("", ".context")
      )
  } else if (nrow(tls_raw) == nrow(context_raw)) {
    join_key <- "row_order"
    tls_joined <- bind_cols(
      tls_raw,
      context_raw %>%
        select(any_of(c(
          "spatial_x",
          "spatial_y",
          "imagecol",
          "imagerow",
          "Raw_DualHighCategory",
          "Raw_DualHigh_Category",
          "DualHighCategory",
          "Category"
        ))) %>%
        rename_with(~ paste0(.x, ".context"))
    )
  } else {
    stop(
      "TLS module table lacks coordinates/category and cannot be joined to context table.\n",
      "Set manual_tls_module_file to a spot-level table with module scores, coordinates, and Raw_DualHighCategory."
    )
  }
  
  x_col <- find_col(
    tls_joined,
    patterns = c("^spatial_x$", "^spatial_x\\.context$", "^imagecol$", "^imagecol\\.context$", "^pxl_col_in_fullres$", "^x$", "pixel.*x", "coord.*x", "array_col", "col"),
    label = "spatial x after context join"
  )
  
  y_col <- find_col(
    tls_joined,
    patterns = c("^spatial_y$", "^spatial_y\\.context$", "^imagerow$", "^imagerow\\.context$", "^pxl_row_in_fullres$", "^y$", "pixel.*y", "coord.*y", "array_row", "row"),
    label = "spatial y after context join"
  )
  
  cat_col <- find_col(
    tls_joined,
    patterns = c("Raw.*Dual.*High.*Category", "Raw.*High.*Category", "Raw.*Category", "DualHighCategory", "Dual_High_Category", "High.*Category", "^Category$", "Category\\.context"),
    label = "raw high-state category after context join"
  )
}

plot_df <- tls_joined %>%
  transmute(
    x = as.numeric(.data[[x_col]]),
    y = as.numeric(.data[[y_col]]),
    Category = standardize_category(.data[[cat_col]]),
    !!!setNames(
      lapply(unname(module_cols), function(cc) as.numeric(tls_joined[[cc]])),
      names(module_cols)
    )
  ) %>%
  filter(is.finite(x), is.finite(y)) %>%
  mutate(
    Category = factor(Category, levels = category_order)
  )

if (nrow(plot_df) == 0) stop("Harmonized S19 spot-level table has zero rows.")

################--------------------------------------------
## 3. Panel A: module score boxplots by raw high-state category
################--------------------------------------------

panelA_df <- plot_df %>%
  select(Category, all_of(module_display_order)) %>%
  filter(!is.na(Category)) %>%
  pivot_longer(
    cols = all_of(module_display_order),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(
    Module = factor(Module, levels = module_display_order),
    Category = factor(Category, levels = category_order)
  ) %>%
  filter(is.finite(Score))

if (nrow(panelA_df) == 0) stop("Panel A table has zero rows.")

################--------------------------------------------
## 4. Panel B: spatial distribution of TLS/humoral modules
################--------------------------------------------

panelB_df <- plot_df %>%
  select(x, y, all_of(module_display_order)) %>%
  pivot_longer(
    cols = all_of(module_display_order),
    names_to = "Module",
    values_to = "Score"
  ) %>%
  mutate(Module = factor(Module, levels = module_display_order)) %>%
  filter(is.finite(Score), is.finite(x), is.finite(y))

if (nrow(panelB_df) == 0) stop("Panel B table has zero rows.")

score_limits <- quantile(panelB_df$Score, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
score_abs <- max(abs(score_limits), na.rm = TRUE)
if (!is.finite(score_abs) || score_abs == 0) score_abs <- 1

################--------------------------------------------
## 5. Panel C: spatial autocorrelation statistics by module and k
################--------------------------------------------

compute_knn_indices <- function(coords, k) {
  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knn(coords, k = k)
    return(nn$nn.index)
  }
  
  d <- as.matrix(stats::dist(coords))
  diag(d) <- Inf
  idx <- t(apply(d, 1, function(v) order(v)[seq_len(k)]))
  idx
}

univariate_moran_from_index <- function(x, idx) {
  x <- as.numeric(x)
  keep <- is.finite(x)
  xz <- as.numeric(scale(x[keep]))
  
  if (!all(keep)) {
    idx <- idx[keep, , drop = FALSE]
  }
  
  mean(vapply(
    seq_along(xz),
    function(i) mean(xz[i] * xz[idx[i, ]]),
    numeric(1)
  ))
}

complete_rows <- which(is.finite(plot_df$x) & is.finite(plot_df$y))
for (m in module_display_order) {
  complete_rows <- intersect(complete_rows, which(is.finite(plot_df[[m]])))
}

if (length(complete_rows) < 100) stop("Too few complete rows for Moran analysis.")

coords <- as.matrix(plot_df[complete_rows, c("x", "y")])
coords <- apply(coords, 2, as.numeric)

panelC_df <- list()
for (k in knn_values) {
  idx <- compute_knn_indices(coords, k = k)
  
  tmp <- lapply(module_display_order, function(m) {
    data.frame(
      Module = m,
      k = k,
      Moran = univariate_moran_from_index(plot_df[[m]][complete_rows], idx),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()
  
  panelC_df[[as.character(k)]] <- tmp
}

panelC_df <- bind_rows(panelC_df) %>%
  mutate(
    Module = factor(Module, levels = module_display_order),
    k = as.numeric(k)
  )

if (nrow(panelC_df) == 0) stop("Panel C table has zero rows.")

################--------------------------------------------
## 6. Panel D: operational co-high spatial features
################--------------------------------------------

cutoff <- vapply(
  module_display_order,
  function(m) as.numeric(stats::quantile(plot_df[[m]], probs = high_quantile, na.rm = TRUE, names = FALSE)),
  numeric(1)
)

b_high <- plot_df[["B-cell"]] >= cutoff[["B-cell"]]
t_high <- plot_df[["T-cell zone"]] >= cutoff[["T-cell zone"]]
ig_high <- plot_df[["Ig-humoral"]] >= cutoff[["Ig-humoral"]]
tls_high <- plot_df[["TLS chemokine core"]] >= cutoff[["TLS chemokine core"]]

feature <- dplyr::case_when(
  tls_high & b_high & t_high ~ "TLS-chemokine/B-cell/T-cell co-high",
  ig_high & t_high ~ "Ig-humoral/T-cell co-high",
  tls_high ~ "TLS chemokine high/other",
  ig_high ~ "Ig-humoral high/other",
  TRUE ~ "other"
)

panelD_df <- plot_df %>%
  transmute(
    x = x,
    y = y,
    Feature = factor(feature, levels = feature_order)
  ) %>%
  filter(is.finite(x), is.finite(y), !is.na(Feature))

if (nrow(panelD_df) == 0) stop("Panel D table has zero rows.")

################--------------------------------------------
## 7. Plots
################--------------------------------------------

pA <- ggplot(panelA_df, aes(x = Category, y = Score)) +
  geom_boxplot(width = 0.62, outlier.size = 0.28, outlier.alpha = 0.45, linewidth = 0.35) +
  facet_wrap(~ Module, ncol = 3, scales = "free_y") +
  scale_x_discrete(labels = wrap_category) +
  labs(
    title = "Humoral and TLS-like transcriptional modules across Visium high-state categories",
    subtitle = "Used for cautious TLS-like spatial context, not mature TLS calling",
    x = "Raw high-state category",
    y = "Module score"
  ) +
  theme_s19 +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.8),
    strip.text = element_text(size = 8.1, face = "bold"),
    plot.margin = margin(5, 5, 12, 5)
  )

pB <- ggplot(panelB_df, aes(x = x, y = y, color = Score)) +
  geom_point(size = 0.28, alpha = 0.88) +
  coord_equal() +
  scale_y_reverse() +
  facet_wrap(~ Module, ncol = 3) +
  scale_color_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    limits = c(-score_abs, score_abs),
    oob = scales::squish,
    name = "Module\nscore"
  ) +
  labs(
    title = "Spatial distribution of humoral/TLS-like modules",
    subtitle = "Transcriptional spatial context; not mature TLS identification",
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_s19 +
  theme(
    axis.title = element_text(size = 9.0, face = "bold"),
    axis.text = element_text(size = 6.8),
    strip.text = element_text(size = 7.7, face = "bold"),
    legend.title = element_text(size = 8.0, face = "bold"),
    legend.text = element_text(size = 7.0)
  )

pC <- ggplot(panelC_df, aes(x = k, y = Moran)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_line(linewidth = 0.45, color = "black") +
  geom_point(size = 1.75, color = "black") +
  facet_wrap(~ Module, ncol = 3) +
  scale_x_continuous(
    breaks = knn_values,
    labels = as.character(knn_values)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  labs(
    title = "Spatial autocorrelation of humoral/TLS-like modules",
    subtitle = "Positive values indicate spatial autocorrelation within the Visium section",
    x = "k nearest neighbors",
    y = "Spatial autocorrelation statistic"
  ) +
  theme_s19 +
  theme(
    strip.text = element_text(size = 7.7, face = "bold"),
    axis.text = element_text(size = 7.4)
  )

pD <- ggplot(panelD_df, aes(x = x, y = y, color = Feature)) +
  geom_point(size = 0.42, alpha = 0.90) +
  coord_equal() +
  scale_y_reverse() +
  scale_color_manual(
    values = feature_colors,
    breaks = feature_order,
    labels = wrap_feature,
    name = "Operational feature",
    drop = FALSE
  ) +
  labs(
    title = "Operational humoral/TLS-like co-high spatial features",
    subtitle = "Top-quartile co-high categories; not mature TLS calling",
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_s19 +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7.0),
    legend.title = element_text(size = 8.0, face = "bold"),
    legend.key.height = unit(0.38, "cm"),
    legend.key.width = unit(0.38, "cm")
  )

################--------------------------------------------
## 8. Assemble and save
################--------------------------------------------

supp_s19 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.0, 1.0)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 19, face = "bold", family = "sans"),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

base_filename <- "Supplementary Figure S19. Humoral and TLS-like spatial context in Visium dual-high regions"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

fig_width <- 17.2
fig_height <- 11.2

safe_ggsave(out_png, supp_s19, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s19, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s19, width = fig_width, height = fig_height, dpi = dpi_out, device = cairo_pdf)

################--------------------------------------------
## 9. Audit outputs
################--------------------------------------------

audit <- data.frame(
  tls_module_file = tls_file,
  spot_context_file = manual_spot_context_file,
  join_key = join_key,
  x_col = x_col,
  y_col = y_col,
  category_col = cat_col,
  n_spots = nrow(plot_df),
  n_panelA_rows = nrow(panelA_df),
  n_panelB_rows = nrow(panelB_df),
  n_panelC_rows = nrow(panelC_df),
  n_panelD_rows = nrow(panelD_df),
  high_quantile = high_quantile,
  knn_values = paste(knn_values, collapse = "; "),
  module_columns = paste(names(module_cols), module_cols, sep = "=", collapse = "; "),
  feature_counts = paste(names(table(panelD_df$Feature)), as.integer(table(panelD_df$Feature)), sep = "=", collapse = "; "),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, "Supplementary_Figure_S19_humoral_TLS_like_spatial_context_audit.csv")
write.csv(audit, audit_file, row.names = FALSE)

write.csv(panelA_df, file.path(fig_dir, "Supplementary_Figure_S19_panelA_module_scores_by_high_state_category.csv"), row.names = FALSE)
write.csv(panelB_df, file.path(fig_dir, "Supplementary_Figure_S19_panelB_spatial_module_scores.csv"), row.names = FALSE)
write.csv(panelC_df, file.path(fig_dir, "Supplementary_Figure_S19_panelC_module_spatial_autocorrelation.csv"), row.names = FALSE)
write.csv(panelD_df, file.path(fig_dir, "Supplementary_Figure_S19_panelD_operational_cohigh_features.csv"), row.names = FALSE)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("Supplementary Figure S19 script finished successfully.")