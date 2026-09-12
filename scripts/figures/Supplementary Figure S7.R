############################################################
## Supplementary Figure S7. Distribution of tumor–immune state scores across major cell types in GSE244983
##
## Purpose:
##   1) Load Seurat object from the specified path
##   2) Detect major cell-type and state-score columns
##   3) Generate and export Supplementary Figure S7
##
## Output:
##   D:/ICB_resistance_project/results/figures/supplementary/
##     - Supplementary Figure S7. Distribution of tumor–immune state scores across major cell types in GSE244983.png
##     - Supplementary Figure S7. Distribution of tumor–immune state scores across major cell types in GSE244983.jpg
##     - Supplementary Figure S7. Distribution of tumor–immune state scores across major cell types in GSE244983.pdf
############################################################

options(stringsAsFactors = FALSE)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
})

############################################################
## 0. Paths
############################################################

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (project_dir == "") {
  project_dir <- "D:/ICB_resistance_project"
}
project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)

fig_dir <- file.path(project_dir, "results", "figures", "supplementary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## 优先指定用户要求的 Seurat 对象路径
manual_rds_path <- file.path(project_dir, "results", "intermediate", "GSE244983", "GSE244983_seurat_state_localized.rds")

############################################################
## 1. Helper functions
############################################################

find_first_col <- function(meta_cols, patterns, required = TRUE, label = "column") {
  for (pat in patterns) {
    hit <- grep(pat, meta_cols, value = TRUE, ignore.case = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  if (required) {
    stop(
      "Could not find ", label, ". Tried patterns:\n",
      paste(patterns, collapse = "\n"),
      "\nAvailable metadata columns:\n",
      paste(meta_cols, collapse = ", ")
    )
  }
  NA_character_
}

theme_s7 <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    strip.background = element_rect(fill = "#D9D9D9", color = "grey50", linewidth = 0.45),
    strip.text = element_text(size = 11.2, face = "bold", color = "black"),
    legend.title = element_text(size = 11, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

safe_ggsave <- function(file, plot, width, height, dpi = 300, device = NULL) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE,
    device = device
  )
  message("Saved: ", normalizePath(file, winslash = "/", mustWork = FALSE))
}

############################################################
## 2. Load Seurat object
############################################################

if (!file.exists(manual_rds_path)) {
  stop("manual_rds_path does not exist: ", manual_rds_path)
}

message("Using Seurat object: ", manual_rds_path)
obj <- readRDS(manual_rds_path)

if (!inherits(obj, "Seurat")) {
  stop("The selected RDS is not a Seurat object: ", manual_rds_path)
}

meta <- obj@meta.data
meta_cols <- colnames(meta)

############################################################
## 3. Detect major cell-type and state-score columns
############################################################

celltype_col <- find_first_col(
  meta_cols,
  patterns = c(
    "^MajorCellType_LOCKED$",
    "^MajorCellType$",
    "^major_cell_type$",
    "^major_celltype$",
    "^CellType$",
    "^cell_type$",
    "^Cell_type$"
  ),
  label = "major cell-type column"
)

immune_col <- find_first_col(
  meta_cols,
  patterns = c(
    "Immune.*defective.*Cold.*score",
    "Immune.*defective.*cold.*score",
    "immune.*defective.*cold",
    "Immune_defective_Cold",
    "Immune_defective_cold"
  ),
  label = "immune-defective/cold score column"
)

myeloid_col <- find_first_col(
  meta_cols,
  patterns = c(
    "Myeloid.*Treg.*Immunosuppressive.*score",
    "Myeloid.*Treg.*immunosuppressive.*score",
    "myeloid.*treg.*immunosuppressive",
    "Myeloid_Treg_Immunosuppressive",
    "Myeloid_Treg_immunosuppressive"
  ),
  label = "myeloid–Treg immunosuppressive score column"
)

tdsr_col <- find_first_col(
  meta_cols,
  patterns = c(
    "Tumor.*dedifferentiation.*Stromal.*remodeling.*score",
    "Tumor.*dedifferentiation.*stromal.*remodeling.*score",
    "tumor.*dediff.*stromal.*remodel",
    "Tumor_dedifferentiation_Stromal_remodeling",
    "Tumor_dedifferentiated",
    "Dediff.*Stromal",
    "dediff.*stromal"
  ),
  label = "tumor-dedifferentiation/stromal-remodeling score column"
)

mel_col <- find_first_col(
  meta_cols,
  patterns = c(
    "Melanocytic.*Differentiation.*score",
    "Melanocytic.*differentiation.*score",
    "melanocytic.*differentiation",
    "Melanocytic_Differentiation",
    "Melanocytic_differentiated"
  ),
  label = "melanocytic differentiation score column"
)

score_cols <- c(
  "immune-defective/cold" = immune_col,
  "myeloid–Treg immunosuppressive" = myeloid_col,
  "tumor-dedifferentiation/stromal-remodeling" = tdsr_col,
  "melanocytic differentiation" = mel_col
)

############################################################
## 4. Build plotting table
############################################################

plot_df <- meta %>%
  dplyr::select(
    MajorCellType = all_of(celltype_col),
    all_of(unname(score_cols))
  )

colnames(plot_df) <- c("MajorCellType", names(score_cols))

celltype_order <- c(
  "Malignant",
  "Cycling malignant",
  "CAF/stromal-like cells",
  "Endothelial",
  "Myeloid cells",
  "T/NK cells",
  "T/NK/Treg-like cells",
  "B/Plasma cells"
)

observed_celltypes <- unique(as.character(plot_df$MajorCellType))
final_celltype_order <- c(
  celltype_order[celltype_order %in% observed_celltypes],
  setdiff(sort(observed_celltypes), celltype_order)
)

plot_df <- plot_df %>%
  mutate(
    MajorCellType = factor(as.character(MajorCellType), levels = final_celltype_order)
  ) %>%
  pivot_longer(
    cols = all_of(names(score_cols)),
    names_to = "State",
    values_to = "Score"
  ) %>%
  mutate(
    State = factor(
      State,
      levels = c(
        "immune-defective/cold",
        "myeloid–Treg immunosuppressive",
        "tumor-dedifferentiation/stromal-remodeling",
        "melanocytic differentiation"
      )
    )
  )

############################################################
## 5. Colors
############################################################

celltype_colors <- c(
  "Malignant" = "#3F5A91",
  "Cycling malignant" = "#4DBBD5",
  "CAF/stromal-like cells" = "#E64B35",
  "Endothelial" = "#7E6148",
  "Myeloid cells" = "#00A087",
  "T/NK cells" = "#8491B4",
  "T/NK/Treg-like cells" = "#91D1C2",
  "B/Plasma cells" = "#F39B7F"
)

missing_colors <- setdiff(final_celltype_order, names(celltype_colors))
if (length(missing_colors) > 0) {
  extra_cols <- scales::hue_pal()(length(missing_colors))
  names(extra_cols) <- missing_colors
  celltype_colors <- c(celltype_colors, extra_cols)
}

############################################################
## 6. Plot
############################################################

p_s7 <- ggplot(plot_df, aes(x = MajorCellType, y = Score, fill = MajorCellType)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "black",
    linewidth = 0.35,
    alpha = 0.95
  ) +
  geom_boxplot(
    width = 0.10,
    outlier.shape = NA,
    color = "grey25",
    fill = "white",
    linewidth = 0.35
  ) +
  facet_wrap(~ State, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  labs(
    title = "Tumor–immune state scores across major cell types in GSE244983",
    x = NULL,
    y = "State score"
  ) +
  theme_s7 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 11.2, face = "bold")
  )

############################################################
## 7. Save outputs at 300 dpi (PNG, JPG, PDF)
############################################################

base_filename <- "Supplementary Figure S7. Distribution of tumor–immune state scores across major cell types in GSE244983"
out_png <- file.path(fig_dir, paste0(base_filename, ".png"))
out_jpg <- file.path(fig_dir, paste0(base_filename, ".jpg"))
out_pdf <- file.path(fig_dir, paste0(base_filename, ".pdf"))

safe_ggsave(out_png, p_s7, width = 13.5, height = 9.2, dpi = 300)
safe_ggsave(out_jpg, p_s7, width = 13.5, height = 9.2, dpi = 300)
safe_ggsave(out_pdf, p_s7, width = 13.5, height = 9.2, dpi = 300, device = cairo_pdf)

cat("\nSupplementary Figure S7 generated successfully in PNG, JPG, and PDF formats (300 DPI)!\n")