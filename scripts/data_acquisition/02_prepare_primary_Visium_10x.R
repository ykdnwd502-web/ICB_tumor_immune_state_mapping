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
## 02_prepare_primary_Visium_10x.R
##
## Downloads and extracts the 10x Genomics public Human
## Melanoma, IF Stained (FFPE) CytAssist output into the exact
## directory contract consumed by 12_Primary_Visium.
############################################################

RAW_DIR <- file.path(
  PROJECT_DIR,
  "data_raw",
  "spatial_10x_human_melanoma_IF_FFPE"
)

dir_create(RAW_DIR)

BASE_URL <- paste0(
  "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/",
  "CytAssist_FFPE_Human_Skin_Melanoma/"
)

matrix_tar <- file.path(
  RAW_DIR,
  "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.tar.gz"
)

spatial_tar <- file.path(
  RAW_DIR,
  "CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz"
)

metrics_csv <- file.path(
  RAW_DIR,
  "CytAssist_FFPE_Human_Skin_Melanoma_metrics_summary.csv"
)

download_one(
  paste0(
    BASE_URL,
    "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.tar.gz"
  ),
  matrix_tar,
  min_bytes = 1024L
)

download_one(
  paste0(
    BASE_URL,
    "CytAssist_FFPE_Human_Skin_Melanoma_spatial.tar.gz"
  ),
  spatial_tar,
  min_bytes = 1024L
)

download_one(
  paste0(
    BASE_URL,
    "CytAssist_FFPE_Human_Skin_Melanoma_metrics_summary.csv"
  ),
  metrics_csv,
  min_bytes = 100L
)

matrix_dir <- file.path(
  RAW_DIR,
  "filtered_feature_bc_matrix"
)

matrix_extract_audit_dir <- file.path(
  RAW_DIR,
  "filtered_feature_bc_matrix_extracted"
)

spatial_dir <- file.path(
  RAW_DIR,
  "spatial"
)

matrix_required <- c(
  file.path(matrix_dir, "barcodes.tsv.gz"),
  file.path(matrix_dir, "features.tsv.gz"),
  file.path(matrix_dir, "matrix.mtx.gz")
)

matrix_audit_required <- c(
  file.path(matrix_extract_audit_dir,
            "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
  file.path(matrix_extract_audit_dir,
            "filtered_feature_bc_matrix", "features.tsv.gz"),
  file.path(matrix_extract_audit_dir,
            "filtered_feature_bc_matrix", "matrix.mtx.gz")
)

spatial_required <- c(
  file.path(spatial_dir, "tissue_positions.csv"),
  file.path(spatial_dir, "scalefactors_json.json")
)

if (!all(file.exists(matrix_required))) {
  message("EXTRACT matrix archive -> raw root")
  utils::untar(
    matrix_tar,
    exdir = RAW_DIR,
    compressed = "gzip"
  )
}

## Keep the second extraction because the validated frozen Visium
## builder explicitly audits this historical extraction directory.
if (!all(file.exists(matrix_audit_required))) {
  dir_create(matrix_extract_audit_dir)
  message("EXTRACT matrix archive -> filtered_feature_bc_matrix_extracted")
  utils::untar(
    matrix_tar,
    exdir = matrix_extract_audit_dir,
    compressed = "gzip"
  )
}

if (!all(file.exists(spatial_required))) {
  message("EXTRACT spatial archive -> raw root")
  utils::untar(
    spatial_tar,
    exdir = RAW_DIR,
    compressed = "gzip"
  )
}

required <- c(
  matrix_required,
  matrix_audit_required,
  spatial_required
)

contract <- data.frame(
  File = norm_path(required),
  Exists = file.exists(required),
  Size_bytes = ifelse(
    file.exists(required),
    file.info(required)$size,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

write_utf8_csv(
  contract,
  file.path(
    PROJECT_DIR,
    "release_metadata",
    "raw_data",
    "VISIUM_10X_RAW_CONTRACT.csv"
  )
)

if (!all(contract$Exists)) {
  stop(
    "10x Visium extraction incomplete. Inspect VISIUM_10X_RAW_CONTRACT.csv",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("02 PRIMARY VISIUM RAW PREPARATION PASS\n")
cat("============================================================\n")
cat("Matrix contract: 3/3\n")
cat("Historical extraction audit copy: 3/3\n")
cat("Spatial contract: 2/2\n")
cat("Expected downstream raw root:\n", RAW_DIR, "\n", sep = "")
cat("============================================================\n")
