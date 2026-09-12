############################################################
## 09A_Visium_diagnose_reproduction_mismatches.R
## DIAGNOSTIC ONLY — does not modify analytical results.
## Purpose:
##   1) resolve exact nonnumeric mismatches in Step02 / Step03;
##   2) resolve row/key differences in Figure 8C / 8D source tables;
##   3) distinguish label/encoding drift from analytical drift.
############################################################

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)

FROZEN_ROOT <- Sys.getenv(
  "ICB_FROZEN_ROOT",
  unset = "D:/ICB_resistance_project_FROZEN_BACKUP"
)

PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)
FROZEN_ROOT <- normalizePath(FROZEN_ROOT, winslash = "/", mustWork = TRUE)

AUDIT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "audit",
  "12_Primary_Visium_sequential_reproduction_check"
)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

pairs <- data.frame(
  group = c(
    rep("Step02_spatial_QC_Moran", 4),
    rep("Step03_dual_high_niche", 3)
  ),
  current = c(
    "results/tables/spatial_melanoma_validation/Step11_spatial_QC_state_Spearman_correlation.csv",
    "results/tables/spatial_melanoma_validation/Step11_spatial_raw_vs_QCresidual_state_correlation.csv",
    "results/tables/spatial_melanoma_validation/Step11_spatial_Morans_I_raw_and_QCresidual.csv",
    "results/tables/spatial_melanoma_validation/Step11_spatial_high_state_neighbor_enrichment_top25_pooled.csv",
    "results/tables/spatial_melanoma_validation/Step12_dual_high_group_counts.csv",
    "results/tables/spatial_melanoma_validation/Step12_dual_high_top_candidate_genes_Figure8C.csv",
    "results/tables/spatial_melanoma_validation/Step12_dual_high_Figure8D_four_group_candidate_gene_expression_summary.csv"
  ),
  frozen = c(
    "results/tables/spatial_melanoma_validation/Step11_spatial_QC_state_Spearman_correlation_CLEAN.csv",
    "results/tables/spatial_melanoma_validation/Step11_spatial_raw_vs_QCresidual_state_correlation_CLEAN.csv",
    "results/tables/spatial_melanoma_validation/Step11_spatial_Morans_I_raw_and_QCresidual_CLEAN.csv",
    "results/tables/spatial_melanoma_validation/Step11_spatial_high_state_neighbor_enrichment_top25_pooled_CLEAN.csv",
    "results/tables/spatial_melanoma_validation/Step12_dual_high_group_counts_CLEAN.csv",
    "results/tables/spatial_melanoma_validation/Step12_dual_high_top_candidate_genes_Figure8C_CLEAN.csv",
    "results/tables/spatial_melanoma_validation/Step12_dual_high_Figure8D_four_group_candidate_gene_expression_summary_CLEAN.csv"
  ),
  stringsAsFactors = FALSE
)

read_exact <- function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
}

norm_dash_ws <- function(x) {
  x <- as.character(x)
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}

norm_case_dash_ws <- function(x) {
  tolower(norm_dash_ws(x))
}

detail_rows <- list()
key_rows <- list()

for (ii in seq_len(nrow(pairs))) {

  cur_path <- file.path(PROJECT_DIR, pairs$current[ii])
  fro_path <- file.path(FROZEN_ROOT, pairs$frozen[ii])

  cat("\n============================================================\n")
  cat(pairs$group[ii], "\n")
  cat("Current: ", cur_path, "\n", sep = "")
  cat("Frozen : ", fro_path, "\n", sep = "")

  if (!file.exists(cur_path) || !file.exists(fro_path)) {
    cat("MISSING FILE — current=", file.exists(cur_path),
        ", frozen=", file.exists(fro_path), "\n", sep = "")
    next
  }

  a <- read_exact(cur_path)
  b <- read_exact(fro_path)

  cat("Rows: current=", nrow(a), " frozen=", nrow(b), "\n", sep = "")
  cat("Columns identical: ", identical(names(a), names(b)), "\n", sep = "")

  common_cols <- intersect(names(a), names(b))

  ## A. Column-level exact character diagnostics when row counts align
  if (nrow(a) == nrow(b)) {
    for (nm in common_cols) {
      aa <- a[[nm]]
      bb <- b[[nm]]

      if (!(is.numeric(aa) && is.numeric(bb))) {
        aa_chr <- as.character(aa)
        bb_chr <- as.character(bb)

        mismatch <- which(
          xor(is.na(aa_chr), is.na(bb_chr)) |
          (!is.na(aa_chr) & !is.na(bb_chr) & aa_chr != bb_chr)
        )

        if (length(mismatch) > 0) {
          dash_equal <- identical(
            norm_dash_ws(aa_chr),
            norm_dash_ws(bb_chr)
          )
          case_dash_equal <- identical(
            norm_case_dash_ws(aa_chr),
            norm_case_dash_ws(bb_chr)
          )

          detail_rows[[length(detail_rows) + 1L]] <- data.frame(
            group = pairs$group[ii],
            file = basename(cur_path),
            column = nm,
            n_mismatch = length(mismatch),
            dash_whitespace_normalized_equal = dash_equal,
            case_dash_whitespace_normalized_equal = case_dash_equal,
            first_row = mismatch[1],
            current_value = aa_chr[mismatch[1]],
            frozen_value = bb_chr[mismatch[1]],
            stringsAsFactors = FALSE
          )

          cat("\nCharacter mismatch column: ", nm,
              " | n=", length(mismatch), "\n", sep = "")
          show_n <- head(mismatch, 8)
          print(
            data.frame(
              row = show_n,
              current = aa_chr[show_n],
              frozen = bb_chr[show_n],
              stringsAsFactors = FALSE
            ),
            row.names = FALSE
          )
          cat("dash/ws normalized equal: ", dash_equal, "\n", sep = "")
          cat("case+dash/ws normalized equal: ", case_dash_equal, "\n", sep = "")
        }
      }
    }
  }

  ## B. Key-set diagnostics for row-count mismatches / derived display subsets
  key_candidates <- list(
    c("Gene", "Group"),
    c("Gene"),
    c("State", "QC_metric"),
    c("State", "ScoreType"),
    c("State"),
    c("Group")
  )

  usable_key <- NULL
  for (kk in key_candidates) {
    if (all(kk %in% common_cols)) {
      usable_key <- kk
      break
    }
  }

  if (!is.null(usable_key)) {
    make_key <- function(df, kk) {
      do.call(
        paste,
        c(lapply(df[kk], as.character), sep = " || ")
      )
    }

    ka <- make_key(a, usable_key)
    kb <- make_key(b, usable_key)

    only_cur <- setdiff(unique(ka), unique(kb))
    only_fro <- setdiff(unique(kb), unique(ka))

    key_rows[[length(key_rows) + 1L]] <- data.frame(
      group = pairs$group[ii],
      file = basename(cur_path),
      key_columns = paste(usable_key, collapse = " + "),
      current_rows = nrow(a),
      frozen_rows = nrow(b),
      current_unique_keys = length(unique(ka)),
      frozen_unique_keys = length(unique(kb)),
      keys_only_current = length(only_cur),
      keys_only_frozen = length(only_fro),
      example_only_current = if (length(only_cur)) only_cur[1] else NA_character_,
      example_only_frozen = if (length(only_fro)) only_fro[1] else NA_character_,
      stringsAsFactors = FALSE
    )

    cat("\nKey diagnostic: ", paste(usable_key, collapse = " + "), "\n", sep = "")
    cat("Unique keys current/frozen: ",
        length(unique(ka)), "/", length(unique(kb)), "\n", sep = "")
    cat("Keys only current: ", length(only_cur), "\n", sep = "")
    if (length(only_cur)) print(head(only_cur, 20))
    cat("Keys only frozen : ", length(only_fro), "\n", sep = "")
    if (length(only_fro)) print(head(only_fro, 20))
  }
}

detail_df <- if (length(detail_rows)) {
  do.call(rbind, detail_rows)
} else {
  data.frame(
    group = character(0),
    file = character(0),
    column = character(0),
    n_mismatch = integer(0),
    dash_whitespace_normalized_equal = logical(0),
    case_dash_whitespace_normalized_equal = logical(0),
    first_row = integer(0),
    current_value = character(0),
    frozen_value = character(0),
    stringsAsFactors = FALSE
  )
}

key_df <- if (length(key_rows)) {
  do.call(rbind, key_rows)
} else {
  data.frame()
}

out1 <- file.path(AUDIT_DIR, "05_DETAILED_CHARACTER_MISMATCH_DIAGNOSTIC.csv")
out2 <- file.path(AUDIT_DIR, "06_ROW_KEY_MISMATCH_DIAGNOSTIC.csv")

utils::write.csv(detail_df, out1, row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(key_df, out2, row.names = FALSE, fileEncoding = "UTF-8")

cat("\n============================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("Character mismatch audit: ", out1, "\n", sep = "")
cat("Row/key mismatch audit   : ", out2, "\n", sep = "")
cat("No analytical result was modified.\n")
cat("============================================================\n")
