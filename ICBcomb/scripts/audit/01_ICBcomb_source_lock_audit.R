############################################################
## 01_ICBcomb_source_lock_audit.R
## ICBcomb input-freeze + source-lock audit (clean-run v1.0.2)
##
## PURPOSE
##   * Validate the ENTIRE immutable ICBcomb input freeze against SHA256SUMS.csv.
##   * Validate the three frozen raw platform downloads and the three final
##     perturbation-query gene sets before any analysis is started.
##   * No live ICBcomb query, no fuzzy file discovery, no auto-install.
############################################################

options(stringsAsFactors = FALSE)

ROOT <- Sys.getenv("ICBCOMB_ROOT")
if (!nzchar(ROOT)) ROOT <- Sys.getenv("ICB_PUBLIC_ROOT")
if (!nzchar(ROOT)) ROOT <- getwd()
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = TRUE)

if (!dir.exists(file.path(ROOT, "inputs", "ICBcomb_INPUT_FREEZE_v1.0"))) {
  stop("Step 01 resolved the wrong ICBcomb root: ", ROOT, "\nExpected inputs/ICBcomb_INPUT_FREEZE_v1.0 directly under this root.", call. = FALSE)
}

INPUT <- file.path(ROOT, "inputs", "ICBcomb_INPUT_FREEZE_v1.0")
AUDIT <- file.path(ROOT, "results", "audit", "01_ICBcomb_source_lock")
dir.create(AUDIT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Package 'digest' is required for ICBcomb SHA256 validation.", call. = FALSE)
}

manifest_file <- file.path(INPUT, "SHA256SUMS.csv")
if (!file.exists(manifest_file)) {
  stop("Missing frozen input manifest: ", manifest_file, call. = FALSE)
}

manifest <- read.csv(
  manifest_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_manifest_cols <- c("relative_path", "sha256")
if (!all(required_manifest_cols %in% names(manifest))) {
  stop(
    "SHA256SUMS.csv is missing required columns: ",
    paste(setdiff(required_manifest_cols, names(manifest)), collapse = ", "),
    call. = FALSE
  )
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

manifest_check <- manifest
manifest_check$absolute_path <- file.path(INPUT, manifest_check$relative_path)
manifest_check$exists <- file.exists(manifest_check$absolute_path)
manifest_check$sha256_observed <- vapply(
  manifest_check$absolute_path,
  function(x) if (file.exists(x)) sha256_file(x) else NA_character_,
  character(1)
)
manifest_check$pass <-
  manifest_check$exists &
  tolower(manifest_check$sha256_observed) == tolower(manifest_check$sha256)

write.csv(
  manifest_check,
  file.path(AUDIT, "01_input_freeze_sha256_validation.csv"),
  row.names = FALSE
)

if (!all(manifest_check$pass)) {
  print(manifest_check[!manifest_check$pass, , drop = FALSE])
  stop(
    "ICBcomb input-freeze SHA256 validation failed: ",
    sum(!manifest_check$pass),
    " blocking mismatch(es).",
    call. = FALSE
  )
}

RAW_DIR <- file.path(INPUT, "data", "raw_platform_results")
GENE_DIR <- file.path(
  INPUT,
  "data",
  "query_gene_set_bundle",
  "02_final_ICBcomb_gene_sets"
)

raw_expected <- c(
  "ICBcomb_Immune_cold_RESTORE_hybrid_results.csv",
  "ICBcomb_Myeloid_Treg_SUPPRESS_hybrid_results.csv",
  "ICBcomb_Tumor_dedifferentiation_SUPPRESS_hybrid_results.csv"
)
query_expected <- c(
  "Immune_defective_Cold_RESTORE_hybrid.txt",
  "Myeloid_Treg_Immunosuppressive_SUPPRESS_hybrid.txt",
  "Tumor_dedifferentiation_Stromal_remodeling_SUPPRESS_hybrid.txt"
)
query_expected_n <- c(98L, 73L, 81L)

raw_paths <- file.path(RAW_DIR, raw_expected)
query_paths <- file.path(GENE_DIR, query_expected)

if (!all(file.exists(raw_paths))) {
  stop(
    "Missing frozen raw ICBcomb platform file(s):\n",
    paste(raw_paths[!file.exists(raw_paths)], collapse = "\n"),
    call. = FALSE
  )
}
if (!all(file.exists(query_paths))) {
  stop(
    "Missing frozen final ICBcomb query file(s):\n",
    paste(query_paths[!file.exists(query_paths)], collapse = "\n"),
    call. = FALSE
  )
}

raw_rows <- vapply(
  raw_paths,
  function(x) nrow(read.csv(x, check.names = FALSE, stringsAsFactors = FALSE)),
  integer(1)
)
query_n <- vapply(
  query_paths,
  function(x) length(readLines(x, warn = FALSE)),
  integer(1)
)

source_contract <- data.frame(
  item = c(
    "manifest_entries",
    "manifest_pass",
    "raw_platform_files",
    "raw_rows_by_file",
    "raw_rows_total",
    "final_query_gene_counts",
    "final_query_directions"
  ),
  observed = c(
    as.character(nrow(manifest_check)),
    paste0(sum(manifest_check$pass), "/", nrow(manifest_check)),
    as.character(length(raw_paths)),
    paste(raw_rows, collapse = "/"),
    as.character(sum(raw_rows)),
    paste(query_n, collapse = "/"),
    "RESTORE/SUPPRESS/SUPPRESS"
  ),
  expected = c(
    "86",
    "86/86",
    "3",
    "15/15/15",
    "45",
    "98/73/81",
    "RESTORE/SUPPRESS/SUPPRESS"
  ),
  pass = c(
    nrow(manifest_check) == 86L,
    all(manifest_check$pass),
    length(raw_paths) == 3L,
    identical(as.integer(raw_rows), c(15L, 15L, 15L)),
    sum(raw_rows) == 45L,
    identical(as.integer(query_n), query_expected_n),
    TRUE
  ),
  stringsAsFactors = FALSE
)

write.csv(
  source_contract,
  file.path(AUDIT, "02_source_contract.csv"),
  row.names = FALSE
)

if (!all(source_contract$pass)) {
  print(source_contract[!source_contract$pass, , drop = FALSE])
  stop("ICBcomb frozen source contract failed.", call. = FALSE)
}

cat("============================================================\n")
cat("ICBcomb SOURCE-LOCK AUDIT: PASS\n")
cat("Input freeze SHA256: 86/86 PASS\n")
cat("Raw platform files: 3 x 15 = 45 rows\n")
cat("Final query genes: 98 / 73 / 81\n")
cat("Directions: RESTORE / SUPPRESS / SUPPRESS\n")
cat("============================================================\n")
