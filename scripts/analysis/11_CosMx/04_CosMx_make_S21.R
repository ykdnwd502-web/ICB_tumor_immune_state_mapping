############################################################
## STYLE-ONLY v1.1.1 — LONG-SUBTITLE ADAPTIVE FIX
## S20-reference visual ratio retained.
## Panel B subtitle is locally reduced and wrapped at a semantic boundary
## to prevent horizontal clipping. All other figure proportions, canvas,
## analytical inputs, state definitions, and numerical outputs are unchanged.
############################################################

############################################################
## 04_CosMx_make_S21.R
##
## PUBLIC SEQUENTIAL ROLE
##   Step03 primary metadata-QC cell table -> final Supplementary Figure S21
##
## This is the frozen canonical S21 rebuild authority, renamed for the
## simple public sequence. Internal source-lock to metadata-v2 is retained.
############################################################

############################################################
## 21H2A S21 CosMx CANONICAL rebuild v1.0
## Based on Supplementary Figure S21 reviewer-fixed script v3 FINAL
## CosMx coordinate matching and RNA-defined state remapping overview
##
## Why v3 FINAL:
##   v2 fixed the y-axis orientation, but the A/B panel labels were still
##   clipped or not visible in the exported image. This final version uses
##   a safer patchwork tag system with visible A/B labels and no negative
##   tag coordinates.
##
## Final fixes:
##   1) Keeps scale_y_reverse() in both panels to match the original
##      image-coordinate display orientation.
##
##   2) Uses white background with light grey grid, aligned with S11-S20.
##
##   3) Uses final four-state names:
##        immune-defective/cold
##        myeloid–Treg immunosuppressive
##        tumor-dedifferentiation/stromal-remodeling
##        melanocytic differentiation
##
##   4) Uses patchwork plot_annotation(tag_levels = "A") so A/B labels are
##      visible and not clipped.
##
##   5) Exports PNG, JPG, and PDF; raster outputs are 300 dpi.
##
## Primary input:
##   D:/ICB_resistance_project/results/tables/
##   CosMx_metadata_QC/
##   CosMx_state_dual_high_with_coordinates_AND_metadata_QC.csv
##
## If needed, set manual_cosmx_file below.
############################################################

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
transmute <- dplyr::transmute

############################################################
## 0. Paths and manual override
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

canonical_root <- file.path(project_dir, "results", "canonical", "S21_S22_CosMx_v1.0")
fig_dir <- file.path(canonical_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dpi_out <- 300
out_prefix <- "Supplementary_Figure_S21_CosMx_spatial_state_distribution"

## Manual override. Leave blank to use preferred input or auto-search fallback.
manual_cosmx_file <- ""

preferred_cosmx_file <- file.path(
  project_dir,
  "results", "tables",
  "CosMx_metadata_QC",
  "CosMx_state_dual_high_with_coordinates_AND_metadata_QC.csv"
)

############################################################
## 1. Helper functions
############################################################

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

safe_ggsave <- function(file, plot, width, height, dpi = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
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
  score <- score + ifelse(grepl("coordinate|coordinates|metadata|polygon|global|spatial", f_low), 30, 0)
  score <- score + ifelse(grepl("dual_high|state|dominant", f_low), 30, 0)
  score <- score + ifelse(grepl("enhanced|qc", f_low), 10, 0)
  score <- score - ifelse(grepl("audit|summary|moran|auc|permutation|null|coverage", f_low), 50, 0)

  cand <- data.frame(file = files, score = score, stringsAsFactors = FALSE) %>%
    arrange(desc(score), file) %>%
    filter(score > 0)

  message("Top candidates for S21 CosMx coordinate/state table:")
  print(utils::head(cand, 20))

  for (f in cand$file) {
    x <- tryCatch(read_any_table(f), error = function(e) NULL)
    if (is.null(x)) next

    cn <- colnames(x)
    has_x <- any(grepl("^spatial_x$|global.*x|coordinate.*x|coord.*x|^x$", cn, ignore.case = TRUE, perl = TRUE))
    has_y <- any(grepl("^spatial_y$|global.*y|coordinate.*y|coord.*y|^y$", cn, ignore.case = TRUE, perl = TRUE))
    has_state <- any(grepl("dominant.*state|state|final.*state|category", cn, ignore.case = TRUE, perl = TRUE))

    if (has_x && has_y && has_state) {
      message("Selected S21 CosMx file: ", f)
      return(f)
    }
  }

  stop("No suitable CosMx coordinate/state table found. Set manual_cosmx_file.")
}

standardize_state <- function(x) {
  x0 <- as.character(x)
  x1 <- x0
  x1 <- gsub("_", " ", x1)
  x1 <- gsub("-", " ", x1)
  x1 <- gsub("\\s+", " ", x1)
  xl <- tolower(x1)

  dplyr::case_when(
    grepl("immune.*defective|cold", xl) ~ "immune-defective/cold",
    grepl("myeloid.*treg|treg.*immunosuppress", xl) ~ "myeloid–Treg immunosuppressive",
    grepl("dediff|stromal.*remodel|tumor.*stromal|stromal remodeling", xl) ~ "tumor-dedifferentiation/stromal-remodeling",
    grepl("melanocytic|melanoma.*differentiation|differentiation", xl) ~ "melanocytic differentiation",
    TRUE ~ x0
  )
}

wrap_state_legend <- function(x) {
  dplyr::recode(
    as.character(x),
    "immune-defective/cold" = "immune-defective/\ncold",
    "myeloid–Treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
    "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling",
    "melanocytic differentiation" = "melanocytic\ndifferentiation",
    .default = as.character(x)
  )
}

state_order <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling",
  "melanocytic differentiation"
)

state_colors <- c(
  "immune-defective/cold" = "#4DBBD5",
  "myeloid–Treg immunosuppressive" = "#00A087",
  "tumor-dedifferentiation/stromal-remodeling" = "#E64B35",
  "melanocytic differentiation" = "#3C5488"
)

theme_s21 <- theme_bw(base_size = 12.5, base_family = "sans") +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 12.5, face = "bold", color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12.5, face = "bold", color = "black"),
    legend.text = element_text(size = 10.75, color = "black", lineheight = 0.95),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width = unit(0.55, "cm"),
    panel.grid.major = element_line(color = "#E8E8E8", linewidth = 0.30),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.60),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(6, 6, 6, 6)
  )

############################################################
## 2. Load and harmonize CosMx data
############################################################

cosmx_file <- manual_cosmx_file
if (!nzchar(cosmx_file)) {
  if (file.exists(preferred_cosmx_file)) {
    cosmx_file <- preferred_cosmx_file
  } else {
    cosmx_file <- auto_find_cosmx_file(project_dir)
  }
}
cosmx_file <- normalizePath(cosmx_file, winslash = "/", mustWork = TRUE)

locked_primary_source <- normalizePath(
  file.path(
    project_dir, "results", "tables",
    "CosMx_metadata_QC",
    "CosMx_state_dual_high_with_coordinates_AND_metadata_QC.csv"
  ),
  winslash = "/",
  mustWork = TRUE
)

if (!identical(cosmx_file, locked_primary_source)) {
  stop(
    paste0(
      "Canonical source-lock violation. Expected primary metadata-QC authority: ",
      locked_primary_source,
      " ; observed: ",
      cosmx_file
    ),
    call. = FALSE
  )
}

cosmx_raw <- read_any_table(cosmx_file)

message("S21 CosMx input file: ", cosmx_file)
message("S21 CosMx input columns:")
message(paste(colnames(cosmx_raw), collapse = ", "))

x_col <- find_col(
  cosmx_raw,
  patterns = c(
    "^spatial_x$",
    "^global_x$",
    "^x_global$",
    "^global\\.x$",
    "polygon.*x",
    "coordinate.*x",
    "coord.*x",
    "^x$"
  ),
  label = "CosMx global/spatial x"
)

y_col <- find_col(
  cosmx_raw,
  patterns = c(
    "^spatial_y$",
    "^global_y$",
    "^y_global$",
    "^global\\.y$",
    "polygon.*y",
    "coordinate.*y",
    "coord.*y",
    "^y$"
  ),
  label = "CosMx global/spatial y"
)

state_col <- find_col(
  cosmx_raw,
  patterns = c(
    "^Dominant_state$",
    "^DominantState$",
    "dominant.*state",
    "^FinalState$",
    "final.*state",
    "rna.*state",
    "tumor.*immune.*state",
    "^state$",
    "category"
  ),
  label = "RNA-defined dominant tumor–immune state"
)

sample_col <- find_col(
  cosmx_raw,
  patterns = c("sample", "fov", "slide", "roi", "section", "patient"),
  label = "sample/FOV identifier",
  required = FALSE
)

plot_df <- cosmx_raw %>%
  transmute(
    x = as.numeric(.data[[x_col]]),
    y = as.numeric(.data[[y_col]]),
    state_raw = as.character(.data[[state_col]]),
    state = standardize_state(.data[[state_col]]),
    sample_id = if (!is.na(sample_col)) as.character(.data[[sample_col]]) else NA_character_
  ) %>%
  filter(is.finite(x), is.finite(y)) %>%
  mutate(
    state = factor(state, levels = state_order)
  ) %>%
  filter(!is.na(state))

if (nrow(plot_df) == 0) {
  stop("S21 harmonized CosMx table has zero rows after coordinate/state filtering.")
}

message("S21 harmonized rows: ", nrow(plot_df))
message("State counts:")
print(table(plot_df$state, useNA = "ifany"))

n_cells <- nrow(plot_df)
n_cells_label <- scales::comma(n_cells)

############################################################
## 3. Build panels
############################################################

pA <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_point(size = 0.075, alpha = 0.42, color = "grey40") +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "CosMx coordinate matching overview",
    subtitle = paste0("Matched CosMx cells with cell-level spatial/global coordinates; y-axis inverted for display; n = ", n_cells_label),
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_s21 +
  theme(legend.position = "none")

pB <- ggplot(plot_df, aes(x = x, y = y, color = state)) +
  geom_point(size = 0.105, alpha = 0.78) +
  coord_equal() +
  scale_y_reverse() +
  scale_color_manual(
    values = state_colors,
    breaks = state_order,
    labels = wrap_state_legend,
    name = "Dominant state",
    drop = FALSE
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 2.4, alpha = 1),
      title.position = "top",
      title.hjust = 0,
      byrow = TRUE
    )
  ) +
  labs(
    title = "CosMx remapping of RNA-defined tumor–immune states",
    subtitle = "Unit: CosMx cell; cell-level spatial/global coordinates; y-axis inverted for display;\nlocalization overview, not functional validation",
    x = "Spatial x",
    y = "Spatial y"
  ) +
  theme_s21 +
  theme(
    plot.subtitle = element_text(size = 10.75, lineheight = 1.0, margin = margin(b = 6)),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 8)
  )

############################################################
## 4. Assemble and save
############################################################

## This is the key fix: use patchwork tags at the composite level.
## No negative tag coordinate is used, so A/B labels remain visible.
supp_s21 <- pA / pB +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold", color = "black", family = "sans"),
      plot.tag.position = c(0.005, 0.995),
      plot.margin = margin(4, 6, 4, 6)
    )
  )

out_png <- file.path(fig_dir, paste0(out_prefix, ".png"))
out_jpg <- file.path(fig_dir, paste0(out_prefix, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(out_prefix, ".pdf"))

fig_width <- 9.3
fig_height <- 11.2

safe_ggsave(out_png, supp_s21, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_jpg, supp_s21, width = fig_width, height = fig_height, dpi = dpi_out)
safe_ggsave(out_pdf, supp_s21, width = fig_width, height = fig_height, dpi = dpi_out)

############################################################
## 5. Audit outputs
############################################################

state_counts <- as.data.frame(table(plot_df$state), stringsAsFactors = FALSE)
colnames(state_counts) <- c("state", "n_cells")

audit <- data.frame(
  cosmx_file = cosmx_file,
  x_col = x_col,
  y_col = y_col,
  y_axis_inverted_for_display = TRUE,
  panel_labels_generated_by = "patchwork plot_annotation(tag_levels = 'A')",
  state_col = state_col,
  sample_col = sample_col,
  n_cells = n_cells,
  x_min = min(plot_df$x, na.rm = TRUE),
  x_max = max(plot_df$x, na.rm = TRUE),
  y_min = min(plot_df$y, na.rm = TRUE),
  y_max = max(plot_df$y, na.rm = TRUE),
  state_counts = paste(state_counts$state, state_counts$n_cells, sep = "=", collapse = "; "),
  output_png = out_png,
  output_jpg = out_jpg,
  output_pdf = out_pdf,
  stringsAsFactors = FALSE
)

audit_file <- file.path(fig_dir, paste0(out_prefix, "_audit.csv"))
write.csv(audit, audit_file, row.names = FALSE)

write.csv(
  state_counts,
  file.path(fig_dir, "Supplementary_Figure_S21_CosMx_state_counts.csv"),
  row.names = FALSE
)

write.csv(
  plot_df %>%
    select(x, y, state_raw, state, sample_id),
  file.path(fig_dir, "Supplementary_Figure_S21_CosMx_coordinate_state_plot_input.csv"),
  row.names = FALSE
)

message("Saved audit: ", normalizePath(audit_file, winslash = "/", mustWork = FALSE))
message("21H2A S21 CosMx CANONICAL rebuild v1.0 finished successfully.")


## Public sequential end gate
.step04_src <- file.path(
  fig_dir,
  "Supplementary_Figure_S21_CosMx_coordinate_state_plot_input.csv"
)
if (!file.exists(.step04_src)) stop("STEP 06 S21 source table missing: ", .step04_src, call. = FALSE)
.step04_n <- nrow(utils::read.csv(.step04_src, check.names = FALSE))
cat("\n============================================================\n")
cat("STEP 04 PASS — final S21 rebuilt\n")
cat("S21 plot-source cells: ", .step04_n, "\n", sep = "")
cat("============================================================\n")
