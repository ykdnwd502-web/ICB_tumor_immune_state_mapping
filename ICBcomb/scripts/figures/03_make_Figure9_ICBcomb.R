############################################################
## 03_make_Figure9_ICBcomb.R
## Genuine Figure 9 reporting rebuild v1.0.3
##
## Reads ONLY the four fixed Figure 9 source tables produced by step 02.
## No analytical recalculation, fuzzy search, or result-table mutation.
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(patchwork)
  library(grid)
})

ROOT <- Sys.getenv("ICBCOMB_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("ICB_PUBLIC_ROOT")
if (!nzchar(ROOT)) ROOT <- getwd()
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(ROOT, "scripts"))) stop("Wrong ICBcomb module root resolved in 03_make_Figure9_ICBcomb.R: ", ROOT, call. = FALSE)

T <- file.path(ROOT, "results", "tables", "ICBcomb")
F <- file.path(ROOT, "results", "figures")
dir.create(F, recursive = TRUE, showWarnings = FALSE)

files <- c(
  A = "Figure9_Panel_A_query_counts.csv",
  B = "Figure9_Panel_B_top_strategies.csv",
  C = "Figure9_Panel_C_state_strategy_matrix.csv",
  D = "Figure9_Panel_D_state_maximum.csv"
)
paths <- stats::setNames(
  file.path(T, unname(files)),
  names(files)
)
if (!all(file.exists(paths))) {
  stop(
    "Missing Figure 9 source table(s):\n",
    paste(paths[!file.exists(paths)], collapse = "\n"),
    call. = FALSE
  )
}

A <- read.csv(paths[["A"]], check.names = FALSE, stringsAsFactors = FALSE)
B <- read.csv(paths[["B"]], check.names = FALSE, stringsAsFactors = FALSE)
C <- read.csv(paths[["C"]], check.names = FALSE, stringsAsFactors = FALSE)
D <- read.csv(paths[["D"]], check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- list(
  A = c("state", "n_query_genes"),
  B = c("state", "rank_within_state", "strategy", "median_prioritization_score", "support_n"),
  C = c("state", "strategy", "median_prioritization_score", "support_n"),
  D = c("state", "max_median_prioritization_score")
)
source_tables <- list(A = A, B = B, C = C, D = D)
for (nm in names(required_cols)) {
  miss <- setdiff(required_cols[[nm]], colnames(source_tables[[nm]]))
  if (length(miss) > 0L) {
    stop(
      "Figure 9 source schema mismatch for panel ", nm, ": missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
}

stopifnot(nrow(A) == 3L, nrow(B) == 10L, nrow(C) == 12L, nrow(D) == 3L)

STATE_LEVELS <- c(
  "immune-defective/cold",
  "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling"
)
STATE_COLS <- c(
  "immune-defective/cold" = "#4DBBD5",
  "myeloid–Treg immunosuppressive" = "#00A087",
  "tumor-dedifferentiation/stromal-remodeling" = "#E64B35"
)
STATE_LABEL <- c(
  "immune-defective/cold" = "immune-defective/cold",
  "myeloid–Treg immunosuppressive" = "myeloid–Treg\nimmunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/\nstromal-remodeling"
)
STATE_STRIP <- c(
  "immune-defective/cold" = "immune-defective/cold",
  "myeloid–Treg immunosuppressive" = "myeloid–Treg immunosuppressive",
  "tumor-dedifferentiation/stromal-remodeling" = "tumor-dedifferentiation/stromal-remodeling"
)

heatmap_low <- "#3B82F6"
heatmap_mid <- "white"
heatmap_high <- "#EF4444"

safe_ggsave <- function(path, plot, width, height, dpi = 300) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "pdf") {
    ggplot2::ggsave(
      filename = path,
      plot = plot,
      width = width,
      height = height,
      device = grDevices::cairo_pdf,
      limitsize = FALSE
    )
  } else {
    ggplot2::ggsave(
      filename = path,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
}

theme_icb <- function(base_size = 10.5) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", color = "black", hjust = 0),
      strip.background = element_rect(fill = "#D9D9D9", color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      plot.margin = margin(5, 6, 5, 6)
    )
}

## Panel A
A_plot <- A %>%
  mutate(state = factor(state, levels = rev(STATE_LEVELS)))
A_xmax <- max(A_plot$n_query_genes, na.rm = TRUE) * 1.15

pA <- ggplot(A_plot, aes(x = n_query_genes, y = state, fill = state)) +
  geom_col(width = 0.58) +
  geom_text(aes(label = n_query_genes), hjust = -0.15, size = 4.0, color = "black") +
  scale_fill_manual(values = STATE_COLS, drop = FALSE) +
  scale_y_discrete(labels = STATE_LABEL) +
  scale_x_continuous(limits = c(0, A_xmax), breaks = c(0, 30, 60, 90), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A  ICBcomb query gene sets", x = "Number of query genes", y = NULL) +
  theme_icb(base_size = 11.2) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10.2),
    plot.title = element_text(size = 13.2)
  )

## Panel B
B_plot <- B %>%
  mutate(
    state = factor(state, levels = STATE_LEVELS),
    strategy_label = paste0(rank_within_state, ". ", strategy),
    value_label = paste0(sprintf("%.2f", median_prioritization_score), " (n=", support_n, ")")
  ) %>%
  arrange(state, desc(rank_within_state)) %>%
  mutate(
    y_id = paste(state, strategy_label, sep = "___"),
    y_id = factor(y_id, levels = unique(y_id))
  )
B_xmax <- max(B_plot$median_prioritization_score, na.rm = TRUE) * 1.33

pB <- ggplot(B_plot, aes(x = median_prioritization_score, y = y_id, fill = state)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey45") +
  geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  geom_text(
    aes(x = median_prioritization_score + B_xmax * 0.025, label = value_label),
    hjust = 0,
    size = 3.1,
    color = "black"
  ) +
  facet_wrap(~ state, ncol = 1, scales = "free_y", labeller = labeller(state = STATE_STRIP)) +
  scale_y_discrete(labels = function(x) stringr::str_replace(x, "^.*___", "")) +
  scale_fill_manual(values = STATE_COLS, drop = FALSE) +
  coord_cartesian(xlim = c(0, B_xmax), clip = "off") +
  labs(
    title = "B  Top positive candidate reversal strategies",
    x = "Median prioritization score (NES)",
    y = NULL
  ) +
  theme_icb(base_size = 10.2) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.35),
    axis.text.y = element_text(size = 9.2, lineheight = 0.92),
    strip.text = element_text(size = 10.4, face = "bold"),
    plot.title = element_text(size = 13.2),
    plot.margin = margin(5, 42, 5, 5)
  )

## Panel C
strategy_order <- unique(B[order(match(B$state, STATE_LEVELS), B$rank_within_state), "strategy"])
C_plot <- C %>%
  mutate(
    state = factor(state, levels = rev(STATE_LEVELS)),
    strategy = factor(strategy, levels = strategy_order),
    label = ifelse(
      is.na(median_prioritization_score),
      "—",
      paste0(sprintf("%.2f", median_prioritization_score), "\n", "n=", support_n)
    )
  )
C_abs <- max(abs(C_plot$median_prioritization_score), na.rm = TRUE)
if (!is.finite(C_abs) || C_abs <= 0) C_abs <- 1

pC <- ggplot(C_plot, aes(x = strategy, y = state, fill = median_prioritization_score)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 3.4, lineheight = 0.88, color = "black") +
  scale_y_discrete(labels = STATE_LABEL) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    limits = c(-C_abs, C_abs),
    na.value = "#F3F4F6",
    name = "Median\nprioritization\nscore"
  ) +
  labs(title = "C  State-specific matrix", x = NULL, y = NULL) +
  theme_icb(base_size = 11.0) +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9.3),
    axis.text.y = element_text(size = 10.0),
    legend.key.height = grid::unit(0.43, "cm"),
    legend.key.width = grid::unit(0.32, "cm"),
    plot.title = element_text(size = 13.2)
  )

## Panel D
D_plot <- D %>% mutate(state = factor(state, levels = STATE_LEVELS))
D_ymax <- max(D_plot$max_median_prioritization_score, na.rm = TRUE) * 1.20

pD <- ggplot(D_plot, aes(x = state, y = max_median_prioritization_score, fill = state)) +
  geom_col(width = 0.54, color = NA) +
  geom_text(aes(label = sprintf("%.2f", max_median_prioritization_score)), vjust = -0.35, size = 4.0) +
  scale_fill_manual(values = STATE_COLS, drop = FALSE) +
  scale_x_discrete(labels = STATE_LABEL) +
  scale_y_continuous(limits = c(0, D_ymax), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "D  Maximum positive signal", x = NULL, y = "Maximum median score") +
  theme_icb(base_size = 11.0) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 9.5, lineheight = 0.9),
    plot.title = element_text(size = 13.2)
  )

fig9 <- ((pA | pB) / (pC | pD)) +
  plot_layout(heights = c(0.80, 1.35), widths = c(1.05, 1.20)) +
  plot_annotation(
    title = "ICBcomb-based candidate perturbation-prioritization analysis of selected resistance-associated state queries",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 15.8,
        margin = margin(b = 5),
        family = "sans"
      )
    )
  )

out_base <- file.path(F, "Figure9_ICBcomb_perturbation_prioritization")
safe_ggsave(paste0(out_base, ".png"), fig9, width = 16.8, height = 10.0, dpi = 300)
safe_ggsave(paste0(out_base, ".jpg"), fig9, width = 16.8, height = 10.0, dpi = 300)
safe_ggsave(paste0(out_base, ".pdf"), fig9, width = 16.8, height = 10.0)

manuscript_base <- paste0(
  "Figure 9. ICBcomb nomination of candidate perturbation-prioritization signals ",
  "for selected candidate resistance-relevant tumor–immune state queries"
)
for (ext in c("png", "jpg", "pdf")) {
  src <- paste0(out_base, ".", ext)
  dst <- file.path(F, paste0(manuscript_base, ".", ext))
  ok <- file.copy(src, dst, overwrite = TRUE)
  if (!ok) stop("Failed to create manuscript-facing Figure 9 alias: ", dst, call. = FALSE)
}

artifact_audit <- data.frame(
  figure = "Figure9",
  panel_A_rows = nrow(A),
  panel_B_rows = nrow(B),
  panel_C_rows = nrow(C),
  panel_D_rows = nrow(D),
  analytical_recalculation = FALSE,
  fuzzy_autosearch = FALSE,
  png_exists = file.exists(paste0(out_base, ".png")),
  jpg_exists = file.exists(paste0(out_base, ".jpg")),
  pdf_exists = file.exists(paste0(out_base, ".pdf")),
  stringsAsFactors = FALSE
)
write.csv(artifact_audit, file.path(F, "Figure9_artifact_audit.csv"), row.names = FALSE)

if (!all(artifact_audit[c("png_exists", "jpg_exists", "pdf_exists")])) {
  stop("Figure 9 artifact generation failed.", call. = FALSE)
}

cat("Figure 9 rebuild: PASS (3/10/12/3 source rows; PNG/JPG/PDF written)\n")
