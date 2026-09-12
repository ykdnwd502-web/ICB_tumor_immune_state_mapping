############################################################
## 01_CellChat_candidate_LR_nomination.R
##
## Public sequential producer derived from the frozen Step16 v4
## source-locked/path-safe CellChat analysis.
##
## SCIENTIFIC CONTRACT — UNCHANGED
##   seed = 20260504
##   max 800 cells/group
##   min 30 cells/group
##   CellChatDB.human full database
##   computeCommunProb(type="truncatedMean", trim=0.1)
##   filterCommunication(min.cells=30)
##   12 tracked candidate ligand-receptor axes
##
## Source lock:
##   public expression object:
##     results/intermediate/GSE244983/GSE244983_seurat_major_annotated.rds
##   immutable historical Step16 cell/group selection manifest:
##     locked_input/CellChat_historical_sampling_group_manifest.csv
##
## The small manifest replaces the historical 2.5-GB Seurat object as the
## public reproducibility boundary. 10D verified that all manifest cells are
## present in the public expression object and reproduce frozen S24 exactly.
##
## Interpretation:
##   candidate nomination only; not causal or functional proof.
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260504)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ------------------------------
# 0. Config
# ------------------------------
project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir,winslash="/",mustWork=TRUE)

out_table_dir <- file.path(project_dir, "results", "tables", "CellChat_candidate_LR")
out_fig_dir   <- file.path(project_dir, "results", "figures", "CellChat_candidate_LR")
out_rds_dir   <- file.path(project_dir, "results", "rds", "CellChat_candidate_LR")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_rds_dir, recursive = TRUE, showWarnings = FALSE)

script_dir <- file.path(
  project_dir,
  "scripts",
  "analysis",
  "14_CellChat_LR"
)

public_seurat_file <- file.path(
  project_dir,
  "results",
  "intermediate",
  "GSE244983",
  "GSE244983_seurat_major_annotated.rds"
)

selection_manifest_file <- file.path(
  script_dir,
  "locked_input",
  "CellChat_historical_sampling_group_manifest.csv"
)

expected_selection_manifest_sha256 <-
  "b4e1dd0cd5f3da4808d9284e23de7491dd4779186021530971c7dad43f579620"

# Historical Step16 selection contract represented by the locked manifest.
max_cells_per_group <- 800
min_cells_per_group <- 30

# CellChat settings.
species <- "human"
# If the presto package is not installed, CellChat may fail unless
# identifyOverExpressedGenes() is run with do.fast = FALSE.
# Options: "auto", TRUE, FALSE. Keep "auto" for clean-room reproducibility.
use_presto_fast_wilcox <- "auto"
use_secreted_only <- FALSE   # FALSE keeps ECM-receptor and cell-cell contact axes such as FN1/COL-CD44.
cellchat_min_cells <- 30

# Candidate L-R axes to explicitly track.
candidate_lr_def <- data.frame(
  pair_family = c(
    "MIF-CD74", "MIF-CD44", "MIF-CXCR4",
    "SPP1-CD44", "SPP1-ITGAV/ITGB1", "SPP1-ITGAV/ITGB5",
    "FN1-CD44", "APP-CD74", "GDF15-TGFBR2",
    "COL1A1-CD44", "COL1A2-CD44", "COL3A1-CD44"
  ),
  ligand_pattern = c(
    "^MIF$", "^MIF$", "^MIF$",
    "^SPP1$", "^SPP1$", "^SPP1$",
    "^FN1$", "^APP$", "^GDF15$",
    "^COL1A1$", "^COL1A2$", "^COL3A1$"
  ),
  receptor_pattern = c(
    "CD74", "CD44", "CXCR4",
    "CD44", "ITGAV.*ITGB1|ITGB1.*ITGAV", "ITGAV.*ITGB5|ITGB5.*ITGAV",
    "CD44", "CD74", "TGFBR2",
    "CD44", "CD44", "CD44"
  ),
  stringsAsFactors = FALSE
)

# Compartments relevant to the dual-high niche.
tumor_stromal_groups <- c(
  "Malignant",
  "Cycling malignant",
  "CAF/stromal-like cells",
  "Endothelial"
)

immune_groups <- c(
  "Myeloid cells",
  "T/NK cells",
  "T/NK/Treg-like cells",
  "B/plasma cells"
)

focused_groups <- unique(c(tumor_stromal_groups, immune_groups))

# ------------------------------
# 1. Utilities
# ------------------------------
safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE)
  message("Saved table: ", path)
}

safe_num <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

sanitize_names <- function(df) {
  nm <- colnames(df)
  nm[is.na(nm) | nm == ""] <- paste0("blank_col_", which(is.na(nm) | nm == ""))
  colnames(df) <- make.unique(nm, sep = "_")
  df
}

clean_text <- function(x) {
  x <- as.character(x)
  trimws(x)
}

safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  if (requireNamespace("digest", quietly = TRUE)) {
    return(
      digest::digest(
        file = path,
        algo = "sha256",
        serialize = FALSE
      )
    )
  }

  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, "rb")
    on.exit(
      try(close(con), silent = TRUE),
      add = TRUE
    )

    return(
      paste0(
        as.character(
          openssl::sha256(con)
        ),
        collapse = ""
      )
    )
  }

  stop(
    "Package 'digest' or 'openssl' is required for locked-input SHA256 verification.",
    call. = FALSE
  )
}

get_meta <- function(obj) {
  if (!is.null(obj@meta.data)) return(obj@meta.data)
  stop("Object does not look like a Seurat object with @meta.data.")
}

guess_celltype_col <- function(meta) {
  candidates <- c(
    "major_celltype", "major_cell_type", "celltype_major", "major_annotation",
    "cell_type", "celltype", "CellType", "cell_type_final", "cell_type_clean",
    "annotation", "Annotation", "seurat_annotations", "predicted.celltype",
    "manual_annotation", "final_celltype"
  )

  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) > 0) return(hit[1])

  # Soft search fallback.
  soft <- colnames(meta)[grepl("cell.*type|celltype|annotation|major", colnames(meta), ignore.case = TRUE)]
  if (length(soft) > 0) return(soft[1])

  NA_character_
}

standardize_celltype <- function(x) {
  xx <- clean_text(x)
  xl <- tolower(xx)
  xl <- gsub("_", " ", xl)
  xl <- gsub("-", " ", xl)
  xl <- gsub("\\s+", " ", xl)

  out <- xx

  out[grepl("cycling.*malignant|malignant.*cycling|prolif.*malignant", xl)] <- "Cycling malignant"
  out[grepl("malignant|tumou?r|melanoma", xl) & !grepl("cycling|prolif", xl)] <- "Malignant"

  out[grepl("fibroblast/caf|caf|stromal|fibro|fibroblast", xl)] <- "CAF/stromal-like cells"
  out[grepl("endo", xl)] <- "Endothelial"
  out[grepl("myeloid|macrophage|mono|dendritic|dc", xl)] <- "Myeloid cells"

  # Treg-like before generic T/NK.
  out[grepl("treg|t/nk/treg|t nk treg|regulatory", xl)] <- "T/NK/Treg-like cells"
  out[grepl("t/nk|t nk|\\bt cell\\b|\\bt cells\\b|\\bnk\\b|natural killer|lymphocyte", xl) &
        !grepl("treg|regulatory", xl)] <- "T/NK cells"

  out[grepl("b/plasma|b plasma|plasma|\\bb cell\\b|\\bb cells\\b", xl)] <- "B/plasma cells"

  out
}

score_rds_candidate <- function(path, obj) {
  if (is.null(obj)) {
    return(data.frame(
      file = path, basename = basename(path), readable = FALSE,
      n_cells = NA_integer_, n_features = NA_integer_,
      has_meta = FALSE, guessed_celltype_col = NA_character_,
      n_celltypes = NA_integer_, score = -Inf,
      stringsAsFactors = FALSE
    ))
  }

  n_cells <- tryCatch(ncol(obj), error = function(e) NA_integer_)
  n_features <- tryCatch(nrow(obj), error = function(e) NA_integer_)
  meta <- tryCatch(get_meta(obj), error = function(e) NULL)

  if (is.null(meta)) {
    return(data.frame(
      file = path, basename = basename(path), readable = TRUE,
      n_cells = n_cells, n_features = n_features,
      has_meta = FALSE, guessed_celltype_col = NA_character_,
      n_celltypes = NA_integer_, score = 0,
      stringsAsFactors = FALSE
    ))
  }

  ct_col <- guess_celltype_col(meta)
  n_ct <- if (!is.na(ct_col)) length(unique(meta[[ct_col]])) else NA_integer_

  p_low <- tolower(path)
  score <- 0
  if (grepl("gse244983", p_low)) score <- score + 50
  if (grepl("scrna|single|seurat|annotation|major", p_low)) score <- score + 30
  if (!is.na(ct_col)) score <- score + 50
  if (!is.na(n_cells) && n_cells > 1000) score <- score + 20
  if (!is.na(n_ct) && n_ct >= 3) score <- score + 20
  if (grepl("malignant|cellchat|spatial|cosmx", p_low)) score <- score - 5

  data.frame(
    file = path,
    basename = basename(path),
    readable = TRUE,
    n_cells = n_cells,
    n_features = n_features,
    has_meta = TRUE,
    guessed_celltype_col = ct_col,
    n_celltypes = n_ct,
    score = score,
    stringsAsFactors = FALSE
  )
}

get_assay_for_cellchat <- function(obj) {
  assays <- tryCatch(Seurat::Assays(obj), error = function(e) character(0))
  if ("RNA" %in% assays) return("RNA")
  if ("SCT" %in% assays) return("SCT")
  tryCatch(Seurat::DefaultAssay(obj), error = function(e) assays[1])
}

safe_diet_seurat <- function(obj, assay_use) {
  # Keep only the assay needed by CellChat and drop reductions/graphs when possible.
  # Different Seurat versions expose slightly different DietSeurat arguments, so use tryCatch.
  tryCatch({
    Seurat::DefaultAssay(obj) <- assay_use
    Seurat::DietSeurat(
      obj,
      assays = assay_use,
      dimreducs = NULL,
      graphs = NULL,
      misc = FALSE
    )
  }, error = function(e) {
    message("DietSeurat skipped: ", conditionMessage(e))
    obj
  })
}

match_candidate_pair <- function(comm, lr_def) {
  if (!all(c("ligand", "receptor") %in% colnames(comm))) {
    comm$ligand <- NA_character_
    comm$receptor <- NA_character_
  }

  out <- list()
  k <- 1

  for (i in seq_len(nrow(lr_def))) {
    tmp <- comm %>%
      filter(
        grepl(lr_def$ligand_pattern[i], ligand, ignore.case = TRUE),
        grepl(lr_def$receptor_pattern[i], receptor, ignore.case = TRUE)
      )

    if (nrow(tmp) > 0) {
      tmp$pair_family <- lr_def$pair_family[i]
      tmp$tracked_ligand_pattern <- lr_def$ligand_pattern[i]
      tmp$tracked_receptor_pattern <- lr_def$receptor_pattern[i]
      out[[k]] <- tmp
      k <- k + 1
    }
  }

  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}


run_identify_overexpressed_genes <- function(cellchat, use_presto_fast_wilcox = "auto") {
  has_presto <- requireNamespace("presto", quietly = TRUE)
  fn_formals <- names(formals(CellChat::identifyOverExpressedGenes))
  supports_do_fast <- "do.fast" %in% fn_formals

  if (identical(use_presto_fast_wilcox, "auto")) {
    do_fast <- has_presto
  } else {
    do_fast <- isTRUE(use_presto_fast_wilcox)
  }

  if (supports_do_fast) {
    if (do_fast && has_presto) {
      message("Running CellChat::identifyOverExpressedGenes() with do.fast = TRUE using presto.")
      return(CellChat::identifyOverExpressedGenes(cellchat, do.fast = TRUE))
    } else {
      message("presto is not available or fast Wilcoxon disabled; running CellChat::identifyOverExpressedGenes() with do.fast = FALSE.")
      return(CellChat::identifyOverExpressedGenes(cellchat, do.fast = FALSE))
    }
  }

  message("This CellChat version does not expose do.fast; running identifyOverExpressedGenes() with default arguments.")
  CellChat::identifyOverExpressedGenes(cellchat)
}

# ------------------------------
# 2. Check required packages
# ------------------------------
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Package 'Seurat' is required but not installed.")
}
if (!requireNamespace("CellChat", quietly = TRUE)) {
  stop(
    "Package 'CellChat' is required but not installed. ",
    "Install CellChat before running Step16."
  )
}
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

# ------------------------------
# 3. Load public scRNA object + locked historical selection manifest
# ------------------------------
if (!file.exists(public_seurat_file)) {
  stop(
    "Public GSE244983 Seurat object is missing: ",
    public_seurat_file,
    call. = FALSE
  )
}

if (!file.exists(selection_manifest_file)) {
  stop(
    "Locked CellChat selection manifest is missing: ",
    selection_manifest_file,
    call. = FALSE
  )
}

observed_selection_manifest_sha256 <- sha256_file(
  selection_manifest_file
)

if (
  !identical(
    observed_selection_manifest_sha256,
    expected_selection_manifest_sha256
  )
) {
  stop(
    "Locked CellChat selection manifest SHA256 mismatch. Expected ",
    expected_selection_manifest_sha256,
    " but observed ",
    observed_selection_manifest_sha256,
    ".",
    call. = FALSE
  )
}

selection_manifest <- utils::read.csv(
  selection_manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_manifest_cols <- c(
  "cell_id",
  "cellchat_group",
  "historical_original_celltype"
)

if (
  !all(
    required_manifest_cols %in%
      names(selection_manifest)
  )
) {
  stop(
    "Locked CellChat selection manifest lacks required columns.",
    call. = FALSE
  )
}

if (
  nrow(selection_manifest) != 3219L ||
  length(
    unique(
      as.character(
        selection_manifest$cell_id
      )
    )
  ) != 3219L
) {
  stop(
    "Locked CellChat selection manifest must contain exactly 3,219 unique cell IDs.",
    call. = FALSE
  )
}

manifest_groups <- sort(
  unique(
    as.character(
      selection_manifest$cellchat_group
    )
  )
)

if (
  !setequal(
    manifest_groups,
    focused_groups
  )
) {
  stop(
    "Locked CellChat selection manifest does not contain the expected eight CellChat groups.",
    call. = FALSE
  )
}

obj <- readRDS(
  public_seurat_file
)

if (!inherits(obj, "Seurat")) {
  stop(
    "Public GSE244983 input is not a Seurat object.",
    call. = FALSE
  )
}

selected_file <- public_seurat_file
celltype_col <- "locked historical manifest: cellchat_group"

missing_manifest_cells <- setdiff(
  as.character(
    selection_manifest$cell_id
  ),
  colnames(obj)
)

if (length(missing_manifest_cells) > 0L) {
  stop(
    "Locked historical CellChat cells are missing from the public GSE244983 object. Missing n=",
    length(missing_manifest_cells),
    call. = FALSE
  )
}

object_audit <- data.frame(
  file = selected_file,
  basename = basename(selected_file),
  readable = TRUE,
  n_cells = ncol(obj),
  n_features = nrow(obj),
  has_meta = TRUE,
  guessed_celltype_col = "not used; locked historical manifest",
  n_celltypes = length(manifest_groups),
  score = NA_real_,
  locked_selection_manifest =
    selection_manifest_file,
  locked_selection_manifest_sha256 =
    observed_selection_manifest_sha256,
  manifest_cells = nrow(selection_manifest),
  manifest_cells_missing_from_public =
    length(missing_manifest_cells),
  stringsAsFactors = FALSE
)

safe_write_csv(
  object_audit,
  file.path(
    out_table_dir,
    "CellChat_input_object_audit.csv"
  )
)

message(
  "Selected public scRNA object: ",
  selected_file
)

message(
  "Locked historical CellChat selection manifest: ",
  selection_manifest_file
)

group_counts_all <- selection_manifest %>%
  count(
    cellchat_group,
    name = "n_cells_all"
  ) %>%
  arrange(
    desc(n_cells_all)
  )

safe_write_csv(
  group_counts_all,
  file.path(
    out_table_dir,
    "CellChat_group_counts_before_filter.csv"
  )
)

sampled_cells <- as.character(
  selection_manifest$cell_id
)

sample_audit <- selection_manifest %>%
  count(
    cellchat_group,
    name = "n_cells_available"
  ) %>%
  mutate(
    n_cells_sampled = n_cells_available,
    selection_mode =
      "locked historical Step16 cell/group manifest",
    historical_seed = 20260504L,
    historical_max_cells_per_group =
      max_cells_per_group,
    historical_min_cells_per_group =
      min_cells_per_group
  ) %>%
  arrange(
    desc(n_cells_available)
  )

safe_write_csv(
  sample_audit,
  file.path(
    out_table_dir,
    "CellChat_memory_safe_sampling_audit.csv"
  )
)

message(
  "Subsetting public Seurat object to locked historical Step16 cells; cells = ",
  length(sampled_cells)
)

obj_focus <- subset(
  x = obj,
  cells = sampled_cells
)

group_map <- setNames(
  as.character(
    selection_manifest$cellchat_group
  ),
  as.character(
    selection_manifest$cell_id
  )
)

original_celltype_map <- setNames(
  as.character(
    selection_manifest$historical_original_celltype
  ),
  as.character(
    selection_manifest$cell_id
  )
)

obj_focus$cellchat_group <- unname(
  group_map[
    colnames(obj_focus)
  ]
)

obj_focus$Step16_original_celltype <- unname(
  original_celltype_map[
    colnames(obj_focus)
  ]
)

if (
  anyNA(
    obj_focus$cellchat_group
  ) ||
  any(
    !nzchar(
      trimws(
        as.character(
          obj_focus$cellchat_group
        )
      )
    )
  )
) {
  stop(
    "Locked CellChat group assignment failed after public-object subsetting.",
    call. = FALSE
  )
}

# 10D validated the exact public-object + manifest path without DietSeurat.
# Keep that tested path unchanged here.
assay_use <- get_assay_for_cellchat(
  obj_focus
)

Seurat::DefaultAssay(
  obj_focus
) <- assay_use

gc()

group_counts_used <- data.frame(
  cell = colnames(obj_focus),
  cellchat_group =
    as.character(
      obj_focus$cellchat_group
    ),
  stringsAsFactors = FALSE
) %>%
  count(
    cellchat_group,
    name = "n_cells_used"
  ) %>%
  arrange(
    desc(n_cells_used)
  )

safe_write_csv(
  group_counts_used,
  file.path(
    out_table_dir,
    "CellChat_group_counts.csv"
  )
)

input_summary <- data.frame(
  step16_version =
    "historical_selection_manifest_bridge",
  interpretation_note =
    "Historical Step16 cell/group selection is source-locked as a small immutable manifest; CellChat remains candidate nomination only.",
  selected_scRNA_object =
    selected_file,
  celltype_col_used =
    celltype_col,
  assay_used =
    assay_use,
  max_cells_per_group =
    max_cells_per_group,
  min_cells_per_group =
    min_cells_per_group,
  n_cells_before_filter =
    ncol(obj),
  n_cells_after_filter_downsample =
    ncol(obj_focus),
  groups_used =
    paste(
      group_counts_used$cellchat_group,
      collapse = ";"
    ),
  selection_manifest =
    selection_manifest_file,
  selection_manifest_sha256 =
    observed_selection_manifest_sha256,
  selection_manifest_n =
    nrow(selection_manifest),
  manifest_cells_missing_from_public =
    length(missing_manifest_cells),
  selection_provenance =
    "Frozen historical Step16 broad-celltype standardization and seed=20260504/max800/min30 sampling boundary",
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_summary,
  file.path(
    out_table_dir,
    "CellChat_input_summary.csv"
  )
)

# ------------------------------
# 4. Run CellChat
# ------------------------------
message("Running CellChat with assay: ", assay_use)
message("Groups: ", paste(group_counts_used$cellchat_group, collapse = ", "))

cellchat <- CellChat::createCellChat(
  object = obj_focus,
  group.by = "cellchat_group",
  assay = assay_use
)

if (tolower(species) == "human") {
  data(CellChatDB.human, package = "CellChat")
  CellChatDB <- CellChat::CellChatDB.human
} else {
  data(CellChatDB.mouse, package = "CellChat")
  CellChatDB <- CellChat::CellChatDB.mouse
}

if (use_secreted_only) {
  CellChatDB.use <- CellChat::subsetDB(CellChatDB, search = "Secreted Signaling")
} else {
  CellChatDB.use <- CellChatDB
}

cellchat@DB <- CellChatDB.use

cellchat <- CellChat::subsetData(cellchat)
cellchat <- run_identify_overexpressed_genes(cellchat, use_presto_fast_wilcox = use_presto_fast_wilcox)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = cellchat_min_cells)
cellchat <- CellChat::computeCommunProbPathway(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)

saveRDS(cellchat, file.path(out_rds_dir, "CellChat_object.rds"))

comm <- CellChat::subsetCommunication(cellchat)
comm <- as.data.frame(comm)

safe_write_csv(comm, file.path(out_table_dir, "CellChat_all_communications.csv"))

# ------------------------------
# 5. Candidate L-R nomination
# ------------------------------
candidate_hits <- match_candidate_pair(comm, candidate_lr_def)

if (nrow(candidate_hits) > 0) {
  candidate_hits <- candidate_hits %>%
    mutate(
      source = as.character(source),
      target = as.character(target),
      direction_class = case_when(
        source %in% tumor_stromal_groups & target %in% immune_groups ~ "Tumor/stromal -> immune",
        source %in% immune_groups & target %in% tumor_stromal_groups ~ "Immune -> tumor/stromal",
        source %in% tumor_stromal_groups & target %in% tumor_stromal_groups ~ "Within tumor/stromal",
        source %in% immune_groups & target %in% immune_groups ~ "Within immune",
        TRUE ~ "Other"
      ),
      source_target = paste(source, "→", target)
    ) %>%
    arrange(pair_family, pval, desc(prob))
} else {
  candidate_hits <- data.frame()
}

safe_write_csv(candidate_hits, file.path(out_table_dir, "CellChat_candidate_LR_nomination.csv"))

# Summary for tracked candidate pairs.
if (nrow(candidate_hits) > 0) {
  candidate_summary <- candidate_hits %>%
    group_by(pair_family) %>%
    summarise(
      n_interactions = n(),
      n_source_target_pairs = n_distinct(source_target),
      max_prob = max(prob, na.rm = TRUE),
      median_prob = median(prob, na.rm = TRUE),
      min_p_value = min(pval, na.rm = TRUE),
      top_source_target = source_target[which.max(prob)],
      direction_classes = paste(sort(unique(direction_class)), collapse = ";"),
      nominated_by_cellchat = any(pval < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  candidate_summary <- candidate_lr_def %>%
    transmute(
      pair_family,
      n_interactions = 0,
      n_source_target_pairs = 0,
      max_prob = NA_real_,
      median_prob = NA_real_,
      min_p_value = NA_real_,
      top_source_target = NA_character_,
      direction_classes = NA_character_,
      nominated_by_cellchat = FALSE
    )
}

# Ensure all tracked pairs are represented.
candidate_summary <- candidate_lr_def %>%
  select(pair_family, ligand_pattern, receptor_pattern) %>%
  left_join(candidate_summary, by = "pair_family") %>%
  mutate(
    n_interactions = ifelse(is.na(n_interactions), 0, n_interactions),
    n_source_target_pairs = ifelse(is.na(n_source_target_pairs), 0, n_source_target_pairs),
    nominated_by_cellchat = ifelse(is.na(nominated_by_cellchat), FALSE, nominated_by_cellchat)
  ) %>%
  arrange(desc(nominated_by_cellchat), desc(max_prob), pair_family)

safe_write_csv(candidate_summary, file.path(out_table_dir, "CellChat_candidate_LR_summary.csv"))

# Candidate pairs to pass into CosMx spatial L-R validation.
# Keep all tracked pairs, but explicitly flag CellChat support.
candidate_pairs_for_cosmx <- candidate_lr_def %>%
  left_join(candidate_summary %>% select(pair_family, nominated_by_cellchat, max_prob, min_p_value, top_source_target), by = "pair_family") %>%
  mutate(
    validation_priority = case_when(
      nominated_by_cellchat & pair_family %in% c("MIF-CD74", "MIF-CD44", "MIF-CXCR4", "SPP1-CD44",
                                                 "SPP1-ITGAV/ITGB1", "SPP1-ITGAV/ITGB5",
                                                 "FN1-CD44", "APP-CD74") ~ "High",
      nominated_by_cellchat ~ "Medium",
      pair_family %in% c("MIF-CD74", "MIF-CD44", "SPP1-CD44", "FN1-CD44", "APP-CD74") ~ "Track despite limited CellChat support",
      TRUE ~ "Exploratory"
    ),
    intended_next_step = "CosMx spatial proximity-constrained LR validation"
  ) %>%
  arrange(factor(validation_priority, levels = c("High", "Medium", "Track despite limited CellChat support", "Exploratory")), pair_family)

safe_write_csv(candidate_pairs_for_cosmx, file.path(out_table_dir, "CellChat_candidate_pairs_for_CosMx_validation.csv"))

# ------------------------------
# 6. Figures
# ------------------------------
# 6A. Aggregated CellChat communication count heatmap.
if (!is.null(cellchat@net$count)) {
  net_count <- as.data.frame(as.table(cellchat@net$count))
  colnames(net_count) <- c("source", "target", "count")
  net_count$source <- as.character(net_count$source)
  net_count$target <- as.character(net_count$target)

  p_net <- ggplot(net_count, aes(x = target, y = source, fill = count)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(count > 0, count, "")), size = 3) +
    scale_fill_gradient(low = "white", high = "#E64B35", name = "Count") +
    labs(
      title = "CellChat aggregated communication counts",
      subtitle = "Candidate nomination only; not interpreted as causal communication",
      x = "Target group",
      y = "Source group"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      legend.title = element_text(face = "bold")
    )

  ggsave(file.path(out_fig_dir, "CellChat_communication_count_heatmap.png"),
         p_net, width = 8.5, height = 7.2, dpi = 400, bg = "white")
  ggsave(file.path(out_fig_dir, "CellChat_communication_count_heatmap.pdf"),
         p_net, width = 8.5, height = 7.2, device = cairo_pdf, bg = "white")
}

# 6B. Candidate L-R bubble plot for focused source/target directions.
if (nrow(candidate_hits) > 0) {
  plot_hits <- candidate_hits %>%
    filter(direction_class %in% c("Tumor/stromal -> immune", "Immune -> tumor/stromal")) %>%
    mutate(
      neg_log10_p = -log10(pmax(pval, .Machine$double.xmin)),
      pair_family = factor(pair_family, levels = rev(unique(candidate_lr_def$pair_family))),
      source_target = factor(source_target, levels = unique(source_target[order(direction_class, source, target)]))
    )

  # Keep the plot readable if too many source-target combinations.
  if (nrow(plot_hits) > 0) {
    top_dirs <- plot_hits %>%
      group_by(source_target) %>%
      summarise(max_prob = max(prob, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(max_prob)) %>%
      slice_head(n = 16) %>%
      pull(source_target)

    plot_hits <- plot_hits %>% filter(source_target %in% top_dirs)

    p_bubble <- ggplot(plot_hits, aes(x = source_target, y = pair_family)) +
      geom_point(aes(size = neg_log10_p, color = prob), alpha = 0.85) +
      scale_color_gradient2(low = "#3B82F6", mid = "white", high = "#EF4444",
                            midpoint = median(plot_hits$prob, na.rm = TRUE),
                            name = "CellChat\nprobability") +
      scale_size_continuous(range = c(1.5, 6), name = "-log10(P)") +
      labs(
        title = "CellChat nomination of candidate L-R axes",
        subtitle = "Focused directions between tumor/stromal-remodeling and immune compartments",
        x = "Source → target",
        y = "Candidate L-R axis"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )

    ggsave(file.path(out_fig_dir, "CellChat_candidate_LR_bubble.png"),
           p_bubble, width = 11, height = 7.5, dpi = 450, bg = "white")
    ggsave(file.path(out_fig_dir, "CellChat_candidate_LR_bubble.pdf"),
           p_bubble, width = 11, height = 7.5, device = cairo_pdf, bg = "white")
  }
}

# 6C. Candidate-pair summary bar plot.
plot_summary <- candidate_summary %>%
  mutate(
    max_prob_plot = ifelse(is.finite(max_prob), max_prob, 0),
    pair_family = factor(pair_family, levels = rev(pair_family)),
    support_label = ifelse(nominated_by_cellchat, "CellChat-supported", "Not detected/limited")
  )

p_summary <- ggplot(plot_summary, aes(x = max_prob_plot, y = pair_family, fill = support_label)) +
  geom_col(width = 0.75, color = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = c("CellChat-supported" = "#E64B35", "Not detected/limited" = "#BDBDBD")) +
  labs(
    title = "CellChat support for tracked candidate L-R axes",
    subtitle = "Maximum communication probability across source-target pairs",
    x = "Maximum CellChat probability",
    y = "Tracked candidate L-R axis",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_fig_dir, "CellChat_candidate_LR_summary.png"),
       p_summary, width = 8.5, height = 6.8, dpi = 400, bg = "white")
ggsave(file.path(out_fig_dir, "CellChat_candidate_LR_summary.pdf"),
       p_summary, width = 8.5, height = 6.8, device = cairo_pdf, bg = "white")

# ------------------------------
# 7. Console summary
# ------------------------------
message("\nCellChat candidate L-R nomination finished.")
message("Selected public scRNA object: ", selected_file)
message("CellChat grouping source: ", celltype_col)
message("Locked historical cells used: ", ncol(obj_focus))
message("Groups used: ", paste(group_counts_used$cellchat_group, collapse = ", "))
message("Selection manifest SHA256: ", observed_selection_manifest_sha256)
message("All communications: ", nrow(comm))
message("Candidate L-R hits: ", nrow(candidate_hits))
message("Tables: ", out_table_dir)
message("Figures: ", out_fig_dir)
message("CellChat object: ", file.path(out_rds_dir, "CellChat_object.rds"))
