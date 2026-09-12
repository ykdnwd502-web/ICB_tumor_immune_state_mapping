############################################################
## 05_CellChat_LR_build_reporting_sources.R
##
## Deterministic reporting-source consolidation from the newly
## rerun Step01-04 outputs.
##
## No CellChat or permutation analysis is rerun here.
## BH FDR is recomputed from empirical P values within each
## k-specific 60-test universe and must agree within 1e-12.
############################################################

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("ICB_PROJECT_DIR")
if (!nzchar(PROJECT_DIR)) {
  PROJECT_DIR <- "D:/ICB_resistance_project"
}
PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)

STEP16_DIR <- file.path(
  PROJECT_DIR,
  "results", "tables",
  "CellChat_candidate_LR"
)

STEP17_DIR <- file.path(
  PROJECT_DIR,
  "results", "tables",
  "CosMx_spatial_LR_validation"
)

TABLE_DIR <- file.path(
  PROJECT_DIR,
  "results", "tables",
  "CellChat_LR_reporting"
)

AUDIT_DIR <- file.path(
  PROJECT_DIR,
  "results", "audit",
  "14_CellChat_LR_sequential_reproduction_check"
)

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing sequential source: ", path, call. = FALSE)
  }

  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

as_logical_robust <- function(x) {
  if (is.logical(x)) return(x)

  tolower(
    trimws(
      as.character(x)
    )
  ) %in% c("true", "t", "1", "yes", "y")
}

EXPECTED_AXES <- c(
  "MIF-CD74",
  "MIF-CD44",
  "MIF-CXCR4",
  "SPP1-CD44",
  "SPP1-ITGAV/ITGB1",
  "SPP1-ITGAV/ITGB5",
  "FN1-CD44",
  "APP-CD74",
  "GDF15-TGFBR2",
  "COL1A1-CD44",
  "COL1A2-CD44",
  "COL3A1-CD44"
)

SRC_ALL_COMM <- file.path(
  STEP16_DIR,
  "CellChat_all_communications.csv"
)
SRC_NOM <- file.path(
  STEP16_DIR,
  "CellChat_candidate_LR_nomination.csv"
)
SRC_SUMMARY <- file.path(
  STEP16_DIR,
  "CellChat_candidate_LR_summary.csv"
)
SRC_PAIRS <- file.path(
  STEP16_DIR,
  "CellChat_candidate_pairs_for_CosMx_validation.csv"
)

src_k <- function(k) {
  file.path(
    STEP17_DIR,
    paste0("KNN", k),
    "CosMx_spatial_LR_validation_all.csv"
  )
}

all_comm <- read_csv(SRC_ALL_COMM)
nom <- read_csv(SRC_NOM)
input_summary <- read_csv(SRC_SUMMARY)
pairs <- read_csv(SRC_PAIRS)

if (
  !("pair_family" %in% names(pairs)) ||
  nrow(pairs) != 12L ||
  !setequal(as.character(pairs$pair_family), EXPECTED_AXES)
) {
  stop("Step01 candidate-pair 12-axis contract failed.", call. = FALSE)
}

req_nom <- c(
  "pair_family",
  "source_target",
  "prob",
  "pval",
  "direction_class"
)

if (!all(req_nom %in% names(nom))) {
  stop("CellChat nomination table lacks required columns.", call. = FALSE)
}

rebuild_step16_summary <- function(nom_df) {
  out <- lapply(
    EXPECTED_AXES,
    function(pf) {
      dd <- nom_df[
        as.character(nom_df$pair_family) == pf,
        ,
        drop = FALSE
      ]

      if (nrow(dd) == 0L) {
        return(
          data.frame(
            pair_family = pf,
            n_interactions = 0L,
            n_source_target_pairs = 0L,
            max_prob = NA_real_,
            median_prob = NA_real_,
            min_p_value = NA_real_,
            top_source_target = NA_character_,
            direction_classes = NA_character_,
            nominated_by_cellchat = FALSE,
            stringsAsFactors = FALSE
          )
        )
      }

      pr <- safe_num(dd$prob)
      pv <- safe_num(dd$pval)
      idx <- which.max(pr)

      data.frame(
        pair_family = pf,
        n_interactions = nrow(dd),
        n_source_target_pairs =
          length(unique(as.character(dd$source_target))),
        max_prob = max(pr, na.rm = TRUE),
        median_prob = stats::median(pr, na.rm = TRUE),
        min_p_value = min(pv, na.rm = TRUE),
        top_source_target =
          as.character(dd$source_target[idx]),
        direction_classes =
          paste(
            sort(unique(as.character(dd$direction_class))),
            collapse = ";"
          ),
        nominated_by_cellchat =
          any(pv < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, out)
}

step16_summary <- rebuild_step16_summary(nom)

if (!all(c("source", "target", "prob") %in% names(all_comm))) {
  stop("CellChat all-communications table lacks source/target/prob.", call. = FALSE)
}

nonzero <- (
  is.finite(safe_num(all_comm$prob)) &
  safe_num(all_comm$prob) > 0
)

all_comm_nz <- all_comm[
  nonzero,
  ,
  drop = FALSE
]

source_target_counts <- stats::aggregate(
  rep(1L, nrow(all_comm_nz)),
  by = list(
    source = as.character(all_comm_nz$source),
    target = as.character(all_comm_nz$target)
  ),
  FUN = sum
)

names(source_target_counts)[3] <- "nonzero_interaction_count"

rebuild_step17_summary <- function(res) {
  needed <- c(
    "pair_family",
    "ligand_gene",
    "receptor_genes",
    "nominated_by_cellchat",
    "validation_priority",
    "contrast",
    "status",
    "log2_spatial_enrichment",
    "neighbor_fraction_percent",
    "FDR"
  )

  if (!all(needed %in% names(res))) {
    stop("CosMx L-R all-table lacks required columns.", call. = FALSE)
  }

  pairs0 <- unique(as.character(res$pair_family))

  out <- lapply(
    pairs0,
    function(pf) {
      dd <- res[
        as.character(res$pair_family) == pf,
        ,
        drop = FALSE
      ]

      log2e <- safe_num(dd$log2_spatial_enrichment)
      neigh <- safe_num(dd$neighbor_fraction_percent)
      fdr <- safe_num(dd$FDR)

      interpretation <- if ("interpretation" %in% names(dd)) {
        as.character(dd$interpretation)
      } else {
        ifelse(
          as.character(dd$status) != "ok",
          "not evaluable",
          ifelse(
            is.finite(fdr) & fdr < 0.05 & log2e > 0,
            "spatially enriched",
            ifelse(
              is.finite(fdr) & fdr < 0.05 & log2e < 0,
              "spatially depleted",
              "not significant"
            )
          )
        )
      }

      max_log2 <- if (any(is.finite(log2e))) {
        max(log2e, na.rm = TRUE)
      } else {
        NA_real_
      }

      max_neigh <- if (any(is.finite(neigh))) {
        max(neigh, na.rm = TRUE)
      } else {
        NA_real_
      }

      min_fdr <- if (any(is.finite(fdr))) {
        min(fdr, na.rm = TRUE)
      } else {
        NA_real_
      }

      top_contrast <- if (any(is.finite(log2e))) {
        as.character(
          dd$contrast[
            which.max(
              ifelse(
                is.finite(log2e),
                log2e,
                -Inf
              )
            )
          ][1]
        )
      } else {
        NA_character_
      }

      priority <- as.character(dd$validation_priority[1])
      n_enriched <- sum(
        interpretation == "spatially enriched",
        na.rm = TRUE
      )

      overall_support <- if (
        n_enriched > 0 &&
        priority %in% c("High", "Medium")
      ) {
        "CellChat-nominated and spatially enriched"
      } else if (n_enriched > 0) {
        "Spatially enriched exploratory pair"
      } else {
        "No positive spatial enrichment"
      }

      data.frame(
        pair_family = pf,
        ligand_gene = as.character(dd$ligand_gene[1]),
        receptor_genes = as.character(dd$receptor_genes[1]),
        nominated_by_cellchat =
          as_logical_robust(dd$nominated_by_cellchat[1]),
        validation_priority = priority,
        n_evaluable_contrasts =
          sum(as.character(dd$status) == "ok"),
        n_enriched_contrasts = n_enriched,
        max_log2_spatial_enrichment = max_log2,
        max_neighbor_fraction_percent = max_neigh,
        min_FDR = min_fdr,
        top_contrast = top_contrast,
        overall_support = overall_support,
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, out)
}

all_by_k <- list()
summary_by_k <- list()
fdr_audit <- list()

for (k in c(10, 20, 30)) {
  dd <- read_csv(src_k(k))

  if (
    nrow(dd) != 60L ||
    !all(c("empirical_p", "FDR") %in% names(dd))
  ) {
    stop("k=", k, " sequential 60-test contract failed.", call. = FALSE)
  }

  rec_fdr <- stats::p.adjust(
    safe_num(dd$empirical_p),
    method = "BH"
  )

  old_fdr <- safe_num(dd$FDR)

  exact <- (
    (is.na(old_fdr) & is.na(rec_fdr)) |
    (
      is.finite(old_fdr) &
      is.finite(rec_fdr) &
      abs(old_fdr - rec_fdr) <= 1e-12
    )
  )

  fdr_audit[[as.character(k)]] <- data.frame(
    k = k,
    rows = nrow(dd),
    FDR_exact_1e12 = all(exact),
    stringsAsFactors = FALSE
  )

  if (!all(exact)) {
    stop("k=", k, " BH-FDR does not reproduce from empirical_p.", call. = FALSE)
  }

  dd$FDR <- rec_fdr
  dd$k <- k

  all_by_k[[as.character(k)]] <- dd
  summary_by_k[[as.character(k)]] <-
    rebuild_step17_summary(dd)
}

fdr_audit <- do.call(rbind, fdr_audit)

utils::write.csv(
  fdr_audit,
  file.path(AUDIT_DIR, "CellChat_LR_k_specific_FDR_audit.csv"),
  row.names = FALSE
)

combined_all <- do.call(rbind, all_by_k)

combined_summary <- do.call(
  rbind,
  lapply(
    c(10, 20, 30),
    function(k) {
      x <- summary_by_k[[as.character(k)]]
      x$k <- k
      x
    }
  )
)

robust_rows <- lapply(
  EXPECTED_AXES,
  function(pf) {
    rows <- lapply(
      c(10, 20, 30),
      function(k) {
        x <- summary_by_k[[as.character(k)]]
        x[
          as.character(x$pair_family) == pf,
          ,
          drop = FALSE
        ]
      }
    )

    max_log2 <- vapply(
      rows,
      function(x) {
        if (nrow(x) == 1L) {
          safe_num(x$max_log2_spatial_enrichment[1])
        } else {
          NA_real_
        }
      },
      numeric(1)
    )

    min_fdr <- vapply(
      rows,
      function(x) {
        if (nrow(x) == 1L) {
          safe_num(x$min_FDR[1])
        } else {
          NA_real_
        }
      },
      numeric(1)
    )

    max_neigh <- vapply(
      rows,
      function(x) {
        if (nrow(x) == 1L) {
          safe_num(x$max_neighbor_fraction_percent[1])
        } else {
          NA_real_
        }
      },
      numeric(1)
    )

    enriched <- vapply(
      rows,
      function(x) {
        if (nrow(x) == 1L) {
          safe_num(x$n_enriched_contrasts[1])
        } else {
          NA_real_
        }
      },
      numeric(1)
    )

    sig_pos <- (
      is.finite(max_log2) &
      max_log2 > 0 &
      is.finite(min_fdr) &
      min_fdr < 0.05
    )

    sig_any <- (
      is.finite(min_fdr) &
      min_fdr < 0.05
    )

    data.frame(
      pair_family = pf,
      k_settings_present =
        sum(vapply(rows, nrow, integer(1)) == 1L),
      positive_kNN_n =
        sum(is.finite(max_log2) & max_log2 > 0),
      significant_positive_kNN_n = sum(sig_pos),
      FDR_significant_kNN_n = sum(sig_any),
      total_enriched_contrasts =
        sum(enriched, na.rm = TRUE),
      mean_max_neighbor_fraction_percent =
        if (any(is.finite(max_neigh))) {
          mean(max_neigh, na.rm = TRUE)
        } else {
          NA_real_
        },
      min_FDR_across_k =
        if (any(is.finite(min_fdr))) {
          min(min_fdr, na.rm = TRUE)
        } else {
          NA_real_
        },
      k10_max_log2 = max_log2[1],
      k20_max_log2 = max_log2[2],
      k30_max_log2 = max_log2[3],
      stringsAsFactors = FALSE
    )
  }
)

robust_tbl <- do.call(rbind, robust_rows)

utils::write.csv(
  nom,
  file.path(TABLE_DIR, "CellChat_candidate_nomination.csv"),
  row.names = FALSE
)

utils::write.csv(
  step16_summary,
  file.path(TABLE_DIR, "CellChat_candidate_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  pairs,
  file.path(TABLE_DIR, "CellChat_candidate_pairs.csv"),
  row.names = FALSE
)

utils::write.csv(
  source_target_counts,
  file.path(TABLE_DIR, "CellChat_source_target_nonzero_counts.csv"),
  row.names = FALSE
)

utils::write.csv(
  combined_all,
  file.path(TABLE_DIR, "CosMx_LR_all_kNN.csv"),
  row.names = FALSE
)

utils::write.csv(
  combined_summary,
  file.path(TABLE_DIR, "CosMx_LR_pair_summary_by_k.csv"),
  row.names = FALSE
)

utils::write.csv(
  robust_tbl,
  file.path(TABLE_DIR, "CosMx_LR_kNN_robustness.csv"),
  row.names = FALSE
)

for (k in c(10, 20, 30)) {
  utils::write.csv(
    all_by_k[[as.character(k)]],
    file.path(
      TABLE_DIR,
      paste0("CosMx_LR_K", k, "_all.csv")
    ),
    row.names = FALSE
  )

  utils::write.csv(
    summary_by_k[[as.character(k)]],
    file.path(
      TABLE_DIR,
      paste0("CosMx_LR_K", k, "_summary.csv")
    ),
    row.names = FALSE
  )
}

source_manifest <- data.frame(
  role = c(
    "Step01 all communications",
    "Step01 candidate nomination",
    "Step01 candidate summary",
    "Step01 candidate pairs",
    "k10 all",
    "k20 all",
    "k30 all"
  ),
  path = c(
    SRC_ALL_COMM,
    SRC_NOM,
    SRC_SUMMARY,
    SRC_PAIRS,
    src_k(10),
    src_k(20),
    src_k(30)
  ),
  policy = "NEW_SEQUENTIAL_OUTPUT",
  stringsAsFactors = FALSE
)

utils::write.csv(
  source_manifest,
  file.path(AUDIT_DIR, "CellChat_LR_reporting_source_manifest.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("CELLCHAT/LR REPORTING SOURCE BUILD COMPLETED\n")
cat("============================================================\n")
cat("Tracked axes: ", nrow(step16_summary), "\n", sep = "")
cat("Combined rows: ", nrow(combined_all), "\n", sep = "")
cat("Pair-by-k rows: ", nrow(combined_summary), "\n", sep = "")
cat("Robustness axes: ", nrow(robust_tbl), "\n", sep = "")
cat("============================================================\n")
