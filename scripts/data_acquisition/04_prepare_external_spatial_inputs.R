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
## 04_prepare_external_spatial_inputs.R
##
## Prepares the two external-spatial recurrence sources:
##   A. GSE250636 Visium
##   B. Thrane et al. 2018 legacy ST
############################################################

############################################################
## A. GSE250636
############################################################

GSE_ROOT <- file.path(
  PROJECT_DIR,
  "data_external_spatial_18J",
  "GSE250636"
)

GSE_TAR <- file.path(
  GSE_ROOT,
  "raw_geo",
  "GSE250636",
  "GSE250636_RAW.tar"
)

GSE_TARGET <- file.path(
  GSE_ROOT,
  "unpacked_recursive",
  "GSE250636_RAW"
)

GSE_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE250nnn/GSE250636/suppl/GSE250636_RAW.tar"
)

download_one(
  GSE_URL,
  GSE_TAR,
  min_bytes = 1024L
)

samples <- data.frame(
  GSM = c(
    "GSM7983358",
    "GSM7983359",
    "GSM7983360",
    "GSM7983361",
    "GSM7983362",
    "GSM7983363",
    "GSM7983364",
    "GSM7983365",
    "GSM7983366"
  ),
  Sample = c(
    "sample2",
    "sample4",
    "sample7",
    "sample8",
    "sample12",
    "sample13",
    "sample14",
    "sample15",
    "sample16"
  ),
  stringsAsFactors = FALSE
)

sample_prefix <- paste0(
  samples$GSM,
  "_",
  samples$Sample
)

expected_h5 <- file.path(
  GSE_TARGET,
  paste0(
    sample_prefix,
    "_filtered_feature_bc_matrix.h5"
  )
)

expected_pos <- file.path(
  GSE_TARGET,
  paste0(
    sample_prefix,
    "_tissue_positions_list.csv.gz"
  )
)

gse_complete <- (
  all(file.exists(expected_h5)) &&
  all(file.exists(expected_pos))
)

if (!gse_complete) {

  if (dir.exists(GSE_TARGET) &&
      length(list.files(GSE_TARGET, all.files = FALSE)) > 0L) {
    stop(
      "Incomplete GSE250636 extraction directory already exists:\n",
      GSE_TARGET,
      "\nFor safety, inspect/remove it manually and rerun.",
      call. = FALSE
    )
  }

  dir_create(GSE_TARGET)

  message("UNTAR GSE250636 -> ", GSE_TARGET)

  utils::untar(
    GSE_TAR,
    exdir = GSE_TARGET
  )
}

gse_contract <- data.frame(
  Sample = sample_prefix,
  Filtered_H5 = norm_path(expected_h5),
  Filtered_H5_exists = file.exists(expected_h5),
  Position_file = norm_path(expected_pos),
  Position_exists = file.exists(expected_pos),
  stringsAsFactors = FALSE
)

if (!all(gse_contract$Filtered_H5_exists) ||
    !all(gse_contract$Position_exists)) {
  stop(
    "GSE250636 raw contract failed: expected 9 filtered H5 + 9 position files.",
    call. = FALSE
  )
}

write_utf8_csv(
  gse_contract,
  file.path(
    PROJECT_DIR,
    "release_metadata",
    "raw_data",
    "GSE250636_RAW_CONTRACT.csv"
  )
)

############################################################
## B. Thrane 2018 legacy ST
############################################################

THRANE_ROOT <- file.path(
  PROJECT_DIR,
  "data_external_spatial_18J",
  "Thrane2018_legacyST"
)

THRANE_RAW <- file.path(
  THRANE_ROOT,
  "raw_zenodo"
)

THRANE_UNPACK <- file.path(
  THRANE_ROOT,
  "unpacked"
)

THRANE_ZIP <- file.path(
  THRANE_RAW,
  "Thrane_et_al_2018_CAN_RES.zip"
)

THRANE_URL <- paste0(
  "https://zenodo.org/records/8215682/files/",
  "Thrane_et_al_2018_CAN_RES.zip?download=1"
)

dir_create(THRANE_RAW)
dir_create(THRANE_UNPACK)

## The frozen canonical scripts use the same conservative size floor.
download_one(
  THRANE_URL,
  THRANE_ZIP,
  min_bytes = 850 * 1024^2
)

COUNT_DIR <- file.path(
  THRANE_UNPACK,
  "Thrane_et_al_2018_CAN_RES",
  "Count_data"
)

expected_thrane <- file.path(
  COUNT_DIR,
  c(
    "ST_mel1_rep1_counts.tsv",
    "ST_mel1_rep2_counts.tsv",
    "ST_mel2_rep1_counts.tsv",
    "ST_mel2_rep2_counts.tsv",
    "ST_mel3_rep1_counts.tsv",
    "ST_mel3_rep2_counts.tsv",
    "ST_mel4_rep1_counts.tsv",
    "ST_mel4_rep2_counts.tsv"
  )
)

if (!all(file.exists(expected_thrane))) {

  if (length(list.files(
    THRANE_UNPACK,
    recursive = TRUE,
    all.files = FALSE
  )) > 0L) {
    stop(
      "Incomplete Thrane unpack directory exists:\n",
      THRANE_UNPACK,
      "\nFor safety, inspect/remove it manually and rerun.",
      call. = FALSE
    )
  }

  message("UNZIP Thrane 2018 -> ", THRANE_UNPACK)

  utils::unzip(
    THRANE_ZIP,
    exdir = THRANE_UNPACK
  )
}

thrane_contract <- data.frame(
  File = norm_path(expected_thrane),
  Exists = file.exists(expected_thrane),
  Size_bytes = ifelse(
    file.exists(expected_thrane),
    file.info(expected_thrane)$size,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

if (!all(thrane_contract$Exists)) {
  stop(
    "Thrane 2018 raw contract failed: expected 8 count matrices.",
    call. = FALSE
  )
}

write_utf8_csv(
  thrane_contract,
  file.path(
    PROJECT_DIR,
    "release_metadata",
    "raw_data",
    "THRANE2018_RAW_CONTRACT.csv"
  )
)

cat("\n============================================================\n")
cat("04 EXTERNAL SPATIAL INPUT PREPARATION PASS\n")
cat("============================================================\n")
cat("GSE250636 complete samples: 9/9\n")
cat("Thrane count matrices: 8/8\n")
cat("NOTE: standardized/*.rds is NOT downloaded or frozen; it is regenerated.\n")
cat("============================================================\n")
