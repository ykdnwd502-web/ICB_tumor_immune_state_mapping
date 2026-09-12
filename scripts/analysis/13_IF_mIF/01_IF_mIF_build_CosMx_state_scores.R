############################################################
## 01_IF_mIF_build_CosMx_state_scores.R
##
## PUBLIC SEQUENTIAL IF/mIF — UPSTREAM STATE-SCORE PRODUCER
##
## Scientific logic is inherited unchanged from the frozen Step14B v4
## composite-ID workflow:
##   - CosMx cell key = fov + cell_ID
##   - log1p transform = TRUE
##   - gene-wise z scaling before state-score averaging
##   - direction applied when available
##   - minimum matched genes/state thresholds unchanged
##
## This public version changes only:
##   - project-root resolution
##   - formal output filenames
##   - comments/reporting text
##
## It intentionally retains the original expression/gene-set discovery logic
## so the public rerun can be compared against the frozen historical outputs.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

project_dir <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(project_dir)) project_dir <- "D:/ICB_resistance_project"
project_dir <- normalizePath(project_dir,winslash="/",mustWork=TRUE)

out_table_dir <- file.path(project_dir, "results", "tables", "IF_mIF_validation")
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)

manual_expr_file <- NA_character_
manual_gene_set_file <- NA_character_

min_total_matched_genes <- 8
min_states_with_matches <- 2
min_genes_per_state_to_score <- 2

log1p_transform <- TRUE
use_direction_if_available <- TRUE

state_aliases <- list(
  Immune_defective_Cold = c(
    "Immune_defective_Cold", "Immune_defective_cold",
    "Immune-defective/Cold", "Immune defective Cold", "Cold",
    "Immune-defective_Cold", "Immune.defective.Cold"
  ),
  Myeloid_Treg_Immunosuppressive = c(
    "Myeloid_Treg_Immunosuppressive", "Myeloid_Treg_immunosuppressive",
    "Myeloid–Treg Immunosuppressive", "Myeloid-Treg Immunosuppressive",
    "Myeloid Treg Immunosuppressive", "Myeloid_Treg",
    "Myeloid.Treg.Immunosuppressive"
  ),
  Tumor_dedifferentiation_Stromal_remodeling = c(
    "Tumor_dedifferentiation_Stromal_remodeling",
    "Tumor_dedifferentiated",
    "Tumor-dedifferentiation/Stromal-remodeling",
    "Tumor dedifferentiation Stromal remodeling",
    "Dediff_Stromal", "Stromal_remodeling",
    "Tumor.dedifferentiation.Stromal.remodeling"
  ),
  Melanocytic_Differentiation = c(
    "Melanocytic_Differentiation", "Melanocytic_differentiated",
    "Melanocytic Differentiation", "Melanocytic",
    "Melanocytic.Differentiation"
  )
)

safe_write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE)
  message("Saved table: ", path)
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

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

make_merge_id <- function(fov, cell_id) {
  paste0(clean_text(fov), "__", clean_text(cell_id))
}

clean_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  unique(x)
}

standard_state_name <- function(x) {
  x0 <- as.character(x)
  x0_low <- tolower(gsub("[^A-Za-z0-9]+", "_", x0))
  for (nm in names(state_aliases)) {
    aliases <- unique(c(state_aliases[[nm]], nm))
    aliases_low <- tolower(gsub("[^A-Za-z0-9]+", "_", aliases))
    if (x0_low %in% aliases_low) return(nm)
    if (any(grepl(paste(aliases_low, collapse = "|"), x0_low))) return(nm)
  }
  NA_character_
}

guess_col <- function(cols, patterns) {
  for (p in patterns) {
    hit <- cols[grepl(p, cols, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

read_gmt_file <- function(f) {
  lines <- readLines(f, warn = FALSE)
  lines <- lines[nchar(lines) > 0]
  out <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) next
    raw_state <- parts[1]
    st <- standard_state_name(raw_state)
    if (is.na(st)) next
    genes <- clean_gene(parts[-c(1, 2)])
    if (length(genes) == 0) next
    out[[length(out) + 1]] <- data.frame(
      state = st,
      gene = genes,
      direction = 1,
      source_set_name = raw_state,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(NULL)
  bind_rows(out)
}

read_gs_file <- function(f) {
  ext <- tolower(tools::file_ext(f))
  if (ext == "gmt") return(read_gmt_file(f))
  if (ext == "rds") {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj)) return(NULL)
    if (is.data.frame(obj)) return(sanitize_names(obj))
    if (is.list(obj)) {
      if (!is.null(names(obj)) && all(vapply(obj, function(z) is.atomic(z) || is.factor(z), logical(1)))) {
        lens <- lengths(obj)
        if (sum(lens) > 0) {
          return(data.frame(
            state = rep(names(obj), lens),
            gene = unlist(obj, use.names = FALSE),
            direction = 1,
            stringsAsFactors = FALSE
          ))
        }
      }
      df_items <- obj[vapply(obj, is.data.frame, logical(1))]
      if (length(df_items) > 0) return(sanitize_names(df_items[[1]]))
    }
    return(NULL)
  }
  if (ext == "csv") {
    df <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
    if (!is.null(df)) df <- sanitize_names(df)
    return(df)
  }
  if (ext %in% c("tsv", "txt")) {
    df <- tryCatch(read.delim(f, check.names = FALSE), error = function(e) NULL)
    if (!is.null(df)) df <- sanitize_names(df)
    return(df)
  }
  NULL
}

extract_gene_sets <- function(df_or_gmt) {
  if (is.null(df_or_gmt) || nrow(df_or_gmt) == 0) return(NULL)
  if (all(c("state", "gene") %in% colnames(df_or_gmt))) {
    gs <- df_or_gmt
    if (!"direction" %in% colnames(gs)) gs$direction <- 1
    gs$state <- vapply(gs$state, standard_state_name, character(1))
    gs <- gs[!is.na(gs$state) & !is.na(gs$gene) & gs$gene != "", ]
    return(gs[, c("state", "gene", "direction")])
  }
  df <- df_or_gmt
  cols <- colnames(df)
  gene_col <- guess_col(cols, c("^gene$", "gene_symbol", "symbol", "^Gene$", "genes", "GeneSymbol", "HGNC"))
  state_col <- guess_col(cols, c("FinalState", "^state$", "State", "signature", "program", "group", "Final_State"))
  dir_col <- guess_col(cols, c("DirectionForFinalState", "Direction", "direction", "sign", "weight"))
  if (!is.na(gene_col) && !is.na(state_col)) {
    tmp <- data.frame(
      state_raw = as.character(df[[state_col]]),
      gene = as.character(df[[gene_col]]),
      direction_raw = if (!is.na(dir_col)) as.character(df[[dir_col]]) else "up",
      stringsAsFactors = FALSE
    )
    tmp$state <- vapply(tmp$state_raw, standard_state_name, character(1))
    tmp <- tmp[!is.na(tmp$state) & !is.na(tmp$gene) & tmp$gene != "", ]
    if (nrow(tmp) == 0) return(NULL)
    tmp$direction <- 1
    if (use_direction_if_available && !is.na(dir_col)) {
      tmp$direction[grepl("down|negative|neg|reverse|-1", tmp$direction_raw, ignore.case = TRUE)] <- -1
    }
    return(tmp[, c("state", "gene", "direction")])
  }
  out <- list()
  for (cc in cols) {
    st <- standard_state_name(cc)
    if (!is.na(st)) {
      genes <- clean_gene(df[[cc]])
      if (length(genes) > 0) {
        out[[length(out) + 1]] <- data.frame(state = st, gene = genes, direction = 1, stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out) == 0) return(NULL)
  bind_rows(out)
}

match_genes_to_panel <- function(genes, gene_cols) {
  genes <- clean_gene(genes)
  exact <- intersect(genes, gene_cols)
  gene_map <- setNames(gene_cols, toupper(gene_cols))
  upper_hit <- intersect(toupper(genes), names(gene_map))
  upper <- unname(gene_map[upper_hit])
  unique(c(exact, upper))
}

# 1. Read expression matrix
if (!is.na(manual_expr_file) && file.exists(manual_expr_file)) {
  expr_file <- manual_expr_file
} else {
  expr_candidates <- list.files(
    file.path(project_dir, "data_raw"),
    pattern = "exprMat.*\\.csv$|expr.*matrix.*\\.csv$|matrix\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  expr_candidates <- expr_candidates[grepl("CosMx|Run5611|Slide_4|Dryad", expr_candidates, ignore.case = TRUE)]
  if (length(expr_candidates) == 0) stop("No CosMx expression matrix found. Set manual_expr_file.")
  expr_file <- expr_candidates[order(file.info(expr_candidates)$size, decreasing = TRUE)][1]
}

message("Selected CosMx expression matrix: ", expr_file)
expr <- read.csv(expr_file, check.names = FALSE)
expr <- sanitize_names(expr)

fov_col <- guess_col(colnames(expr), c("^fov$", "FOV", "field"))
cell_col <- guess_col(colnames(expr), c("^cell_ID$", "cellid", "cell_id", "CellID", "Cell_ID", "^cell$"))
if (is.na(fov_col) || is.na(cell_col)) {
  stop("Expression matrix must contain both fov and cell_ID columns for composite merge key.")
}

expr$.merge_id <- make_merge_id(expr[[fov_col]], expr[[cell_col]])
gene_cols <- setdiff(colnames(expr), c(fov_col, cell_col, ".merge_id"))
num_like <- vapply(gene_cols, function(cc) {
  vals <- expr[[cc]]
  mean(is.finite(safe_num(vals[seq_len(min(length(vals), 2000))]))) > 0.80
}, logical(1))
gene_cols <- gene_cols[num_like]

message("Expression matrix rows: ", nrow(expr))
message("Unique composite cell IDs: ", length(unique(expr$.merge_id)))
message("Expression matrix numeric-like genes: ", length(gene_cols))

expr_gene_inventory <- data.frame(
  expr_file = expr_file,
  n_rows = nrow(expr),
  n_unique_composite_ids = length(unique(expr$.merge_id)),
  fov_col = fov_col,
  cell_col = cell_col,
  n_numeric_like_genes = length(gene_cols),
  first_80_genes = paste(head(gene_cols, 80), collapse = ";"),
  stringsAsFactors = FALSE
)
safe_write_csv(expr_gene_inventory, file.path(out_table_dir, "IF_mIF_CosMx_expression_gene_inventory.csv"))

# 2. Find gene sets
if (!is.na(manual_gene_set_file) && file.exists(manual_gene_set_file)) {
  gs_candidates <- manual_gene_set_file
} else {
  gs_candidates <- list.files(
    project_dir,
    pattern = "\\.(gmt|csv|tsv|txt|rds)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  final_like <- grepl(
    "final.*gene.*set|gene.*set.*final|ICBcomb.*gene.*set|final_ICBcomb|query_gene|GMT|gmt",
    gs_candidates,
    ignore.case = TRUE
  )
  other_like <- grepl(
    "gene.set|geneset|signature|ICBcomb|query|FinalState|state_gene|gene_sets|mapping|signature_mapping|final_input",
    gs_candidates,
    ignore.case = TRUE
  )
  gs_candidates <- unique(c(gs_candidates[final_like], gs_candidates[other_like]))
  gs_candidates <- gs_candidates[
    !grepl("exprMat|tx_file|metadata|polygons|processed_metadata|state_scores|raw_signature_scores|bulk_raw_signature_scores",
           gs_candidates, ignore.case = TRUE)
  ]
}

if (length(gs_candidates) == 0) stop("No candidate gene set/signature/GMT file found. Set manual_gene_set_file.")

gs_audit <- list()
parsed_gene_sets <- list()
for (f in unique(gs_candidates)) {
  raw <- read_gs_file(f)
  gs <- extract_gene_sets(raw)
  if (is.null(gs)) {
    gs_audit[[length(gs_audit) + 1]] <- data.frame(
      file = f, n_states = 0, n_signature_genes = 0,
      total_matched_genes = 0, states_with_matches = 0,
      matched_by_state = "", selected_candidate = FALSE,
      stringsAsFactors = FALSE
    )
    next
  }
  gs <- gs %>%
    mutate(
      gene = trimws(as.character(gene)),
      direction = suppressWarnings(as.numeric(direction))
    ) %>%
    mutate(direction = ifelse(is.na(direction), 1, direction)) %>%
    filter(!is.na(gene), gene != "", state %in% names(state_aliases)) %>%
    distinct(state, gene, direction)
  overlap_tbl <- lapply(names(state_aliases), function(st) {
    genes <- gs$gene[gs$state == st]
    matched <- match_genes_to_panel(genes, gene_cols)
    data.frame(
      state = st,
      signature_genes = length(unique(genes)),
      matched_genes_n = length(matched),
      matched_genes = paste(head(matched, 100), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  total_matched <- sum(overlap_tbl$matched_genes_n)
  states_with_matches <- sum(overlap_tbl$matched_genes_n >= min_genes_per_state_to_score)
  ok <- total_matched >= min_total_matched_genes && states_with_matches >= min_states_with_matches
  gs_audit[[length(gs_audit) + 1]] <- data.frame(
    file = f,
    n_states = length(unique(gs$state)),
    n_signature_genes = length(unique(gs$gene)),
    total_matched_genes = total_matched,
    states_with_matches = states_with_matches,
    matched_by_state = paste0(overlap_tbl$state, "=", overlap_tbl$matched_genes_n, collapse = "; "),
    selected_candidate = ok,
    stringsAsFactors = FALSE
  )
  attr(gs, "overlap_tbl") <- overlap_tbl
  parsed_gene_sets[[f]] <- gs
}

gs_audit_df <- bind_rows(gs_audit) %>% arrange(desc(total_matched_genes), desc(states_with_matches))
safe_write_csv(gs_audit_df, file.path(out_table_dir, "IF_mIF_gene_set_candidate_overlap_audit.csv"))

usable <- gs_audit_df %>% filter(selected_candidate)
if (nrow(usable) == 0) {
  stop("No usable gene-set file had sufficient overlap with the CosMx gene panel. Inspect IF_mIF_gene_set_candidate_overlap_audit.csv.")
}

selected_gs_file <- usable$file[1]
gene_sets <- parsed_gene_sets[[selected_gs_file]]
overlap_audit <- attr(gene_sets, "overlap_tbl")

message("Selected gene-set file: ", selected_gs_file)
safe_write_csv(gene_sets, file.path(out_table_dir, "IF_mIF_state_gene_sets_used.csv"))
safe_write_csv(overlap_audit, file.path(out_table_dir, "IF_mIF_state_score_gene_overlap_audit.csv"))

# 3. Compute scores
score_df <- data.frame(
  fov = clean_text(expr[[fov_col]]),
  cell_ID = clean_text(expr[[cell_col]]),
  .merge_id = expr$.merge_id,
  stringsAsFactors = FALSE
)

for (st in names(state_aliases)) {
  gs <- gene_sets %>% filter(state == st)
  overlap <- match_genes_to_panel(gs$gene, gene_cols)
  if (length(overlap) < min_genes_per_state_to_score) {
    warning("Too few matched genes for state ", st, ": ", length(overlap))
    score_df[[st]] <- NA_real_
    next
  }
  mat <- as.matrix(as.data.frame(lapply(expr[, overlap, drop = FALSE], safe_num)))
  if (log1p_transform) mat <- log1p(mat)
  mat_scaled <- scale(mat)
  mat_scaled[!is.finite(mat_scaled)] <- NA_real_

  dir_vec <- rep(1, length(overlap))
  if (use_direction_if_available) {
    dir_tbl <- gs %>% distinct(gene, direction)
    dir_vec <- vapply(overlap, function(g) {
      idx <- which(toupper(dir_tbl$gene) == toupper(g))
      if (length(idx) > 0) return(as.numeric(dir_tbl$direction[idx[1]]))
      1
    }, numeric(1))
    dir_vec[is.na(dir_vec)] <- 1
  }
  mat_scaled <- sweep(mat_scaled, 2, dir_vec, `*`)
  score_df[[st]] <- rowMeans(mat_scaled, na.rm = TRUE)
}

out_score <- file.path(out_table_dir, "IF_mIF_CosMx_cell_state_scores.csv")
safe_write_csv(score_df, out_score)

state_summary <- score_df %>%
  summarise(across(all_of(names(state_aliases)), list(
    n_finite = ~sum(is.finite(.x)),
    mean = ~mean(.x, na.rm = TRUE),
    sd = ~sd(.x, na.rm = TRUE)
  )))
safe_write_csv(state_summary, file.path(out_table_dir, "IF_mIF_CosMx_cell_state_score_summary.csv"))

input_selection_audit <- data.frame(
  expression_matrix = normalizePath(
    expr_file,
    winslash = "/",
    mustWork = TRUE
  ),
  gene_set_file = normalizePath(
    selected_gs_file,
    winslash = "/",
    mustWork = TRUE
  ),
  expression_rows = nrow(expr),
  unique_composite_cell_ids = length(unique(expr$.merge_id)),
  numeric_like_genes = length(gene_cols),
  min_total_matched_genes = min_total_matched_genes,
  min_states_with_matches = min_states_with_matches,
  min_genes_per_state_to_score = min_genes_per_state_to_score,
  log1p_transform = log1p_transform,
  use_direction_if_available = use_direction_if_available,
  stringsAsFactors = FALSE
)

safe_write_csv(
  input_selection_audit,
  file.path(
    out_table_dir,
    "IF_mIF_input_selection_audit.csv"
  )
)

message("IF/mIF Step 01 state-score build finished.")
message("State-score table saved to: ", out_score)
message("Next: run 02_IF_mIF_match_and_correlation.R.")
