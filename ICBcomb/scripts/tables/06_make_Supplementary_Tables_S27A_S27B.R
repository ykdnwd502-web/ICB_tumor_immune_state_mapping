############################################################
## 06_make_Supplementary_Tables_S27A_S27B.R
## Genuine Supplementary Tables S27A/S27B rebuild v1.0.3
##
## Rebuilds both XLSX workbooks from canonical CSV outputs of step 02.
## No analytical recalculation and no hidden pre-existing XLSX dependency.
############################################################

options(stringsAsFactors = FALSE)

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required to rebuild Supplementary Tables S27A/S27B.", call. = FALSE)
}

ROOT <- Sys.getenv("ICBCOMB_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("ICB_PUBLIC_ROOT")
if (!nzchar(ROOT)) ROOT <- getwd()
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(ROOT, "scripts"))) stop("Wrong ICBcomb module root resolved in 06_make_Supplementary_Tables_S27A_S27B.R: ", ROOT, call. = FALSE)

T <- file.path(ROOT, "results", "tables", "ICBcomb")
RPT <- file.path(ROOT, "results", "reports")
TABLE_OUT <- file.path(ROOT, "results", "tables")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_OUT, recursive = TRUE, showWarnings = FALSE)

required <- c(
  query = "ICBcomb_query_counts.csv",
  positive = "ICBcomb_positive_prioritization_rows.csv",
  mapping = "ICBcomb_strategy_mapping.csv",
  claims = "ICBcomb_manuscript_claim_support.csv",
  summary = "ICBcomb_state_strategy_summary.csv",
  A = "Figure9_Panel_A_query_counts.csv",
  B = "Figure9_Panel_B_top_strategies.csv",
  C = "Figure9_Panel_C_state_strategy_matrix.csv",
  D = "Figure9_Panel_D_state_maximum.csv"
)
paths <- stats::setNames(
  file.path(T, unname(required)),
  names(required)
)
if (!all(file.exists(paths))) {
  stop(
    "Missing canonical source table(s) for Supplementary Tables S27A/S27B:\n",
    paste(paths[!file.exists(paths)], collapse = "\n"),
    call. = FALSE
  )
}

X <- lapply(paths, function(x) read.csv(x, check.names = FALSE, stringsAsFactors = FALSE))

stopifnot(
  nrow(X$query) == 3L,
  identical(as.integer(X$query$n_query_genes), c(98L, 73L, 81L)),
  nrow(X$positive) == 42L,
  all(as.numeric(X$positive$NES) > 0),
  nrow(X$mapping) == 45L,
  nrow(X$claims) == 6L && all(as.logical(X$claims$pass)),
  nrow(X$summary) == 10L,
  c(nrow(X$A), nrow(X$B), nrow(X$C), nrow(X$D)) == c(3L, 10L, 12L, 3L)
)

header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#4F81BD",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)
subtle_style <- openxlsx::createStyle(
  fgFill = "#D9EAF7",
  textDecoration = "bold",
  valign = "top"
)
wrap_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")

write_sheet <- function(wb, sheet, df, freeze = TRUE) {
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df, withFilter = FALSE)
  openxlsx::addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
  if (nrow(df) > 0L) {
    openxlsx::addStyle(
      wb, sheet, wrap_style,
      rows = 2:(nrow(df) + 1L),
      cols = seq_len(ncol(df)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  if (freeze) openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = "auto")
  ## Cap extremely wide columns for readable public workbooks.
  if (ncol(df) > 0L) {
    widths <- pmin(35, pmax(10, vapply(df, function(z) {
      min(35, max(nchar(as.character(c(names(z), utils::head(z, 200)))), na.rm = TRUE) + 2)
    }, numeric(1))))
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = widths)
  }
}

readme_A <- data.frame(
  Field = c(
    "Table",
    "Scope",
    "Interpretation boundary",
    "Final query contract",
    "Positive-direction rule"
  ),
  Value = c(
    "Supplementary Table S27A",
    "Frozen query definitions and row-level positive ICBcomb prioritization evidence.",
    "Exploratory database-derived prioritization only; not therapeutic efficacy or causal validation.",
    "98 RESTORE / 73 SUPPRESS / 81 SUPPRESS; melanocytic differentiation is reference/control only.",
    "Only NES > 0 rows after direction-specific RESTORE/SUPPRESS query construction enter Positive_rows."
  ),
  stringsAsFactors = FALSE
)

wbA <- openxlsx::createWorkbook()
write_sheet(wbA, "README", readme_A, freeze = TRUE)
write_sheet(wbA, "Final_query_counts", X$query)
write_sheet(wbA, "Positive_rows", X$positive)
write_sheet(wbA, "Strategy_mapping", X$mapping)
write_sheet(wbA, "Claim_support", X$claims)
openxlsx::addStyle(wbA, "README", subtle_style, rows = 2:(nrow(readme_A) + 1L), cols = 1, gridExpand = TRUE, stack = TRUE)

readme_B <- data.frame(
  Field = c(
    "Table",
    "Scope",
    "Score",
    "Support metrics",
    "Interpretation boundary"
  ),
  Value = c(
    "Supplementary Table S27B",
    "State-by-strategy prioritization summaries and exact Figure 9 reporting-source tables.",
    "Median ICBcomb NES reported as prioritization score.",
    "support_n = positive platform rows; support_dataset_n = distinct archived datasets.",
    "Exploratory database-derived prioritization only; not therapeutic efficacy or clinical actionability."
  ),
  stringsAsFactors = FALSE
)

wbB <- openxlsx::createWorkbook()
write_sheet(wbB, "README", readme_B, freeze = TRUE)
write_sheet(wbB, "Strategy_summary", X$summary)
write_sheet(wbB, "Figure9_A", X$A)
write_sheet(wbB, "Figure9_B", X$B)
write_sheet(wbB, "Figure9_C", X$C)
write_sheet(wbB, "Figure9_D", X$D)
openxlsx::addStyle(wbB, "README", subtle_style, rows = 2:(nrow(readme_B) + 1L), cols = 1, gridExpand = TRUE, stack = TRUE)

outA <- file.path(RPT, "Supplementary_Table_S27A_ICBcomb_positive_rows.xlsx")
outB <- file.path(RPT, "Supplementary_Table_S27B_ICBcomb_strategy_summary.xlsx")
openxlsx::saveWorkbook(wbA, outA, overwrite = TRUE)
openxlsx::saveWorkbook(wbB, outB, overwrite = TRUE)

## Preserve the historical duplicate public locations used by downstream packaging.
copyA <- file.path(TABLE_OUT, basename(outA))
copyB <- file.path(TABLE_OUT, basename(outB))
if (!file.copy(outA, copyA, overwrite = TRUE)) stop("Failed to copy S27A workbook to results/tables.", call. = FALSE)
if (!file.copy(outB, copyB, overwrite = TRUE)) stop("Failed to copy S27B workbook to results/tables.", call. = FALSE)

if (!all(file.exists(c(outA, outB, copyA, copyB)))) {
  stop("Supplementary Tables S27A/S27B rebuild failed.", call. = FALSE)
}

cat("Supplementary Tables S27A/S27B rebuild: PASS\n")
cat("S27A sheets: README / Final_query_counts / Positive_rows / Strategy_mapping / Claim_support\n")
cat("S27B sheets: README / Strategy_summary / Figure9_A / Figure9_B / Figure9_C / Figure9_D\n")
