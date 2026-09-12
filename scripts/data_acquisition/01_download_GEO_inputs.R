############################################################
## ICB_resistance_project — public data acquisition
##
## Project root is controlled by:
##   Sys.getenv("ICB_PROJECT_DIR", unset = "D:/ICB_resistance_project")
##
## Safety:
##   - never deletes source data
##   - existing non-empty downloads are reused
##   - extraction is skipped when the required final raw contract exists
############################################################

options(stringsAsFactors = FALSE)
options(timeout = max(3600, getOption("timeout")))

PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)

PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)

norm_path <- function(x) gsub("\\\\", "/", x)

dir_create <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
  invisible(x)
}

file_ok <- function(x, min_bytes = 1L) {
  file.exists(x) &&
    isTRUE(file.info(x)$size >= min_bytes)
}

download_one <- function(url, dest, min_bytes = 1L) {

  dir_create(dirname(dest))

  if (file_ok(dest, min_bytes)) {
    message("REUSE: ", dest)
    return(invisible(dest))
  }

  if (file.exists(dest)) {
    stop(
      "Existing file is too small/invalid; remove it manually before rerun:\n",
      dest,
      call. = FALSE
    )
  }

  message("DOWNLOAD:\n  ", url, "\n-> ", dest)

  curl_bin <- Sys.which("curl")

  if (nzchar(curl_bin)) {
    status <- system2(
      curl_bin,
      args = c(
        "-L",
        "--fail",
        "--retry", "3",
        "--retry-delay", "5",
        "--connect-timeout", "60",
        "-o", shQuote(dest),
        shQuote(url)
      )
    )

    if (!identical(status, 0L)) {
      stop("curl download failed: ", url, call. = FALSE)
    }
  } else {
    utils::download.file(
      url = url,
      destfile = dest,
      mode = "wb",
      method = "libcurl",
      quiet = FALSE
    )
  }

  if (!file_ok(dest, min_bytes)) {
    stop("Downloaded file failed size validation: ", dest, call. = FALSE)
  }

  invisible(dest)
}

write_utf8_csv <- function(x, path) {
  dir_create(dirname(path))
  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

############################################################
## 01_download_GEO_inputs.R
##
## Downloads the exact GEO supplementary files represented in
## the validated clean-run raw directory.
##
## No decompression is performed for .gz supplementary files:
## the frozen analysis scripts read these compressed files directly.
############################################################

RAW_ROOT <- file.path(PROJECT_DIR, "data_raw")

geo_files <- data.frame(
  Dataset = c(
    "GSE244982",
    "GSE244983",
    "GSE244983",
    "GSE244983",
    "GSE78220",
    "GSE91061",
    "GSE91061",
    "GSE91061",
    "GSE91061"
  ),

  URL = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244982/suppl/GSE244982_ProcessedData_bulkRNAseq.txt.gz",

    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244983/suppl/GSE244983_NormalizedData_scRNAseq.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244983/suppl/GSE244983_RawCounts_scRNAseq.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244983/suppl/GSE244983_SingleCellAnnotations.txt.gz",

    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx",

    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt.gz"
  ),

  Target = c(
    file.path(RAW_ROOT, "GSE244982",
              "GSE244982_ProcessedData_bulkRNAseq.txt.gz"),

    file.path(RAW_ROOT, "GSE244983",
              "GSE244983_NormalizedData_scRNAseq.txt.gz"),
    file.path(RAW_ROOT, "GSE244983",
              "GSE244983_RawCounts_scRNAseq.txt.gz"),
    file.path(RAW_ROOT, "GSE244983",
              "GSE244983_SingleCellAnnotations.txt.gz"),

    file.path(RAW_ROOT, "external_validation", "GSE78220",
              "GSE78220_PatientFPKM.xlsx"),

    file.path(RAW_ROOT, "external_validation", "GSE91061",
              "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz"),
    file.path(RAW_ROOT, "external_validation", "GSE91061",
              "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz"),
    file.path(RAW_ROOT, "external_validation", "GSE91061",
              "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz"),
    file.path(RAW_ROOT, "external_validation", "GSE91061",
              "GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt.gz")
  ),

  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(geo_files))) {
  download_one(
    geo_files$URL[[i]],
    geo_files$Target[[i]],
    min_bytes = 1024L
  )
}

############################################################
## GEO phenotype metadata
##
## The validated raw tree contains pData exports for GSE78220
## and GSE91061. Rebuild them from GEO rather than freezing a
## hand-edited metadata table.
############################################################

pdata_targets <- c(
  GSE78220 = file.path(
    RAW_ROOT, "external_validation", "GSE78220",
    "GSE78220_pData_raw_from_GEO.csv"
  ),
  GSE91061 = file.path(
    RAW_ROOT, "external_validation", "GSE91061",
    "GSE91061_pData_raw_from_GEO.csv"
  )
)

need_pdata <- names(pdata_targets)[
  !vapply(pdata_targets, file_ok, logical(1), min_bytes = 100L)
]

if (length(need_pdata) > 0L) {

  if (!requireNamespace("GEOquery", quietly = TRUE) ||
      !requireNamespace("Biobase", quietly = TRUE)) {
    stop(
      "GEO phenotype metadata must be rebuilt, but GEOquery/Biobase is missing.\n",
      "Install with BiocManager::install(c('GEOquery','Biobase')) and rerun.",
      call. = FALSE
    )
  }

  for (gse in need_pdata) {

    message("FETCH GEO pData: ", gse)

    g <- GEOquery::getGEO(
      gse,
      GSEMatrix = TRUE,
      AnnotGPL = FALSE,
      getGPL = FALSE
    )

    if (!is.list(g) || length(g) < 1L) {
      stop("GEOquery returned no ExpressionSet for ", gse, call. = FALSE)
    }

    pd <- Biobase::pData(g[[1]])
    pd <- data.frame(
      geo_accession = rownames(pd),
      pd,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    write_utf8_csv(
      pd,
      pdata_targets[[gse]]
    )

    message("SAVED: ", pdata_targets[[gse]])
  }
}

geo_files$Exists <- file.exists(geo_files$Target)
geo_files$Size_bytes <- file.info(geo_files$Target)$size

manifest_path <- file.path(
  PROJECT_DIR,
  "release_metadata",
  "raw_data",
  "GEO_DOWNLOAD_MANIFEST.csv"
)

write_utf8_csv(
  geo_files,
  manifest_path
)

if (!all(geo_files$Exists)) {
  stop("One or more required GEO supplementary files are missing.", call. = FALSE)
}

cat("\n============================================================\n")
cat("01 GEO INPUT DOWNLOAD/PREPARATION PASS\n")
cat("============================================================\n")
cat("Supplementary files: ", sum(geo_files$Exists), "/", nrow(geo_files), "\n", sep = "")
cat("pData GSE78220: ", file_ok(pdata_targets[["GSE78220"]], 100L), "\n", sep = "")
cat("pData GSE91061: ", file_ok(pdata_targets[["GSE91061"]], 100L), "\n", sep = "")
cat("============================================================\n")
