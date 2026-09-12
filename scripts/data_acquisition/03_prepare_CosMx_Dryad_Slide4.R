############################################################
## ICB_resistance_project — public data acquisition
##
## 03_prepare_CosMx_Dryad_Slide4_v1.1.R
##
## Purpose
## -------
## Download (when programmatically permitted), stage, extract, and validate
## the exact CosMx Slide 4 raw-data contract used by the frozen public
## analysis chain.
##
## Canonical raw-data target:
##   data_raw/CosMx_melanocytic_tumors_Dryad/
##     Slide_4/
##       Run5611_MK3/
##
## Dryad dataset:
##   DOI: 10.5061/dryad.ksn02v7b1
##   File: Slide_4.zip
##
## v1.1 changes
## ------------
## 1. Hard-codes the current Dryad Slide_4.zip file-stream URL as the
##    default automatic download target.
## 2. Keeps COSMX_SLIDE4_URL as an optional override.
## 3. Supports COSMX_SLIDE4_ZIP as an optional path to an already
##    downloaded local Slide_4.zip.
## 4. Uses browser-like curl headers + Dryad landing-page cookie attempt.
## 5. Validates archive size and ZIP integrity before extraction.
## 6. Extracts into a temporary directory, discovers Run5611_MK3
##    robustly, then installs it at the exact frozen-script path.
## 7. If Dryad rejects programmatic download (e.g. HTTP 403), the script
##    reports the exact browser link and does not modify partial raw data.
##
## Safety
## ------
## - never deletes existing source data;
## - reuses an already complete canonical Slide_4 directory;
## - never overwrites a pre-existing canonical ZIP;
## - partial canonical Slide_4 directories cause a safe STOP;
## - temporary files created by this script may be removed on failure.
############################################################

options(stringsAsFactors = FALSE)
options(timeout = max(7200, getOption("timeout")))

############################################################
## 0. Project configuration
############################################################

PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR",
  unset = "D:/ICB_resistance_project"
)

PROJECT_DIR <- normalizePath(
  PROJECT_DIR,
  winslash = "/",
  mustWork = TRUE
)

norm_path <- function(x) {
  gsub("\\\\", "/", x)
}

dir_create <- function(x) {
  dir.create(
    x,
    recursive = TRUE,
    showWarnings = FALSE
  )
  invisible(x)
}

file_ok <- function(
  x,
  min_bytes = 1
) {
  file.exists(x) &&
    isTRUE(
      !is.na(file.info(x)$size) &&
        file.info(x)$size >= min_bytes
    )
}

write_utf8_csv <- function(
  x,
  path
) {
  dir_create(dirname(path))

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

############################################################
## 1. Paths and authoritative Dryad locations
############################################################

BASE_DIR <- file.path(
  PROJECT_DIR,
  "data_raw",
  "CosMx_melanocytic_tumors_Dryad"
)

ZIP_FILE <- file.path(
  BASE_DIR,
  "Slide_4.zip"
)

SLIDE_DIR <- file.path(
  BASE_DIR,
  "Slide_4"
)

RUN_DIR <- file.path(
  SLIDE_DIR,
  "Run5611_MK3"
)

dir_create(BASE_DIR)

DRYAD_DOI <- "10.5061/dryad.ksn02v7b1"

DRYAD_LANDING_URL <- paste0(
  "https://datadryad.org/dataset/doi:",
  DRYAD_DOI
)

## Current file-specific link exposed by the Dryad dataset page.
## Can be overridden without editing the script:
##   Sys.setenv(COSMX_SLIDE4_URL = "...")
DRYAD_FILE_URL_DEFAULT <-
  "https://datadryad.org/downloads/file_stream/3121158"

DRYAD_FILE_URL <- Sys.getenv(
  "COSMX_SLIDE4_URL",
  unset = DRYAD_FILE_URL_DEFAULT
)

## Optional local archive override:
##   Sys.setenv(COSMX_SLIDE4_ZIP = "D:/Downloads/Slide_4.zip")
LOCAL_ZIP_OVERRIDE <- Sys.getenv(
  "COSMX_SLIDE4_ZIP",
  unset = ""
)

## Dryad reports Slide_4.zip at ~3.68 GB.
## Use a conservative lower bound to reject HTML error pages and truncated
## downloads while avoiding dependence on an exact byte count.
MIN_ZIP_BYTES <- 3.0e9

############################################################
## 2. Frozen raw-data contract
############################################################

required_names <- c(
  "Run5611_MK3_exprMat_file.csv",
  "Run5611_MK3_metadata_file.csv",
  "Run5611_MK3_fov_positions_file.csv",
  "Run5611_MK3_tx_file.csv"
)

required <- file.path(
  RUN_DIR,
  required_names
)

canonical_complete <- function() {
  all(
    file.exists(required)
  )
}

############################################################
## 3. ZIP validation
############################################################

validate_zip <- function(
  path,
  min_bytes = MIN_ZIP_BYTES
) {

  if (!file.exists(path)) {
    return(FALSE)
  }

  size_ok <- isTRUE(
    !is.na(file.info(path)$size) &&
      file.info(path)$size >= min_bytes
  )

  if (!size_ok) {
    return(FALSE)
  }

  zip_list <- tryCatch(
    utils::unzip(
      path,
      list = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(zip_list) ||
      nrow(zip_list) == 0L) {
    return(FALSE)
  }

  ## Require evidence that the archive actually contains the frozen
  ## Run5611_MK3 contract, not merely any valid ZIP.
  zip_names <- gsub(
    "\\\\",
    "/",
    zip_list$Name
  )

  required_in_zip <- vapply(
    required_names,
    function(x) {
      any(
        endsWith(
          zip_names,
          paste0(
            "Run5611_MK3/",
            x
          )
        ) |
          basename(zip_names) == x
      )
    },
    logical(1)
  )

  all(required_in_zip)
}

############################################################
## 4. Optional local ZIP staging
############################################################

stage_local_zip_if_available <- function() {

  if (file.exists(ZIP_FILE)) {
    return(invisible(FALSE))
  }

  candidates <- character()

  if (nzchar(LOCAL_ZIP_OVERRIDE)) {
    candidates <- c(
      candidates,
      LOCAL_ZIP_OVERRIDE
    )
  }

  ## Convenient fallback after a browser download on Windows/macOS/Linux.
  candidates <- c(
    candidates,
    file.path(
      path.expand("~"),
      "Downloads",
      "Slide_4.zip"
    )
  )

  candidates <- unique(
    norm_path(candidates)
  )

  candidates <- candidates[
    file.exists(candidates)
  ]

  if (length(candidates) == 0L) {
    return(invisible(FALSE))
  }

  for (candidate in candidates) {

    message(
      "Checking local Slide_4.zip candidate: ",
      candidate
    )

    if (!validate_zip(candidate)) {
      message(
        "Local candidate rejected by size/ZIP/content validation: ",
        candidate
      )
      next
    }

    message(
      "Staging validated local archive:\n  ",
      candidate,
      "\n-> ",
      ZIP_FILE
    )

    ok <- file.copy(
      candidate,
      ZIP_FILE,
      overwrite = FALSE
    )

    if (!isTRUE(ok)) {
      stop(
        "Could not copy validated local Slide_4.zip to canonical path:\n",
        ZIP_FILE,
        call. = FALSE
      )
    }

    if (!validate_zip(ZIP_FILE)) {
      stop(
        "Canonical copied Slide_4.zip failed post-copy validation.",
        call. = FALSE
      )
    }

    return(invisible(TRUE))
  }

  invisible(FALSE)
}

############################################################
## 5. Programmatic Dryad download
############################################################

try_dryad_download <- function() {

  if (file.exists(ZIP_FILE)) {
    return(
      validate_zip(ZIP_FILE)
    )
  }

  curl_bin <- Sys.which("curl")

  if (!nzchar(curl_bin)) {
    message(
      "System curl not available; automatic Dryad download skipped."
    )
    return(FALSE)
  }

  part_file <- paste0(
    ZIP_FILE,
    ".part"
  )

  cookie_file <- tempfile(
    pattern = "dryad_cookie_"
  )

  landing_tmp <- tempfile(
    pattern = "dryad_landing_",
    fileext = ".html"
  )

  on.exit(
    {
      if (file.exists(cookie_file)) {
        unlink(cookie_file)
      }

      if (file.exists(landing_tmp)) {
        unlink(landing_tmp)
      }

      if (file.exists(part_file) &&
          !validate_zip(part_file)) {
        unlink(part_file)
      }
    },
    add = TRUE
  )

  if (file.exists(part_file)) {
    unlink(part_file)
  }

  user_agent <- paste(
    "Mozilla/5.0",
    "(Windows NT 10.0; Win64; x64)",
    "AppleWebKit/537.36",
    "(KHTML, like Gecko)",
    "Chrome/152.0.0.0 Safari/537.36"
  )

  ##########################################################
  ## 5A. Visit landing page first to obtain any public cookie
  ##########################################################

  message(
    "Dryad landing page:\n  ",
    DRYAD_LANDING_URL
  )

  invisible(
    system2(
      curl_bin,
      args = c(
        "-L",
        "--silent",
        "--show-error",
        "--compressed",
        "--connect-timeout", "60",
        "--user-agent", shQuote(user_agent),
        "--cookie-jar", shQuote(cookie_file),
        "-o", shQuote(landing_tmp),
        shQuote(DRYAD_LANDING_URL)
      )
    )
  )

  ##########################################################
  ## 5B. Download Slide_4.zip
  ##########################################################

  message(
    "Attempting automatic Dryad Slide_4.zip download:\n  ",
    DRYAD_FILE_URL,
    "\n-> ",
    ZIP_FILE
  )

  status <- system2(
    curl_bin,
    args = c(
      "-L",
      "--fail",
      "--show-error",
      "--compressed",
      "--retry", "3",
      "--retry-delay", "5",
      "--connect-timeout", "60",
      "--user-agent", shQuote(user_agent),
      "--referer", shQuote(DRYAD_LANDING_URL),
      "--cookie", shQuote(cookie_file),
      "--cookie-jar", shQuote(cookie_file),
      "-o", shQuote(part_file),
      shQuote(DRYAD_FILE_URL)
    )
  )

  if (!identical(status, 0L)) {

    message(
      "Automatic Dryad download was rejected or failed (curl status ",
      status,
      ")."
    )

    return(FALSE)
  }

  if (!validate_zip(part_file)) {

    message(
      "Downloaded object failed size/ZIP/content validation.\n",
      "It may be an HTML access-denied page or a truncated download."
    )

    return(FALSE)
  }

  ok <- file.rename(
    part_file,
    ZIP_FILE
  )

  if (!isTRUE(ok)) {

    ok <- file.copy(
      part_file,
      ZIP_FILE,
      overwrite = FALSE
    )

    if (!isTRUE(ok)) {
      stop(
        "Validated download could not be moved into canonical location.",
        call. = FALSE
      )
    }

    unlink(part_file)
  }

  validate_zip(ZIP_FILE)
}

############################################################
## 6. Acquire ZIP only when canonical extracted data are absent
############################################################

if (!canonical_complete()) {

  ##########################################################
  ## 6A. Existing partial canonical extraction = safe STOP
  ##########################################################

  if (dir.exists(SLIDE_DIR) &&
      length(
        list.files(
          SLIDE_DIR,
          recursive = TRUE,
          all.files = FALSE
        )
      ) > 0L) {

    stop(
      paste0(
        "Partial/incomplete canonical Slide_4 directory exists:\n",
        SLIDE_DIR,
        "\n\n",
        "For safety, this script will not merge new files into it.\n",
        "Inspect/remove or archive the incomplete directory manually, then rerun."
      ),
      call. = FALSE
    )
  }

  ##########################################################
  ## 6B. Check existing canonical ZIP
  ##########################################################

  if (file.exists(ZIP_FILE)) {

    if (!validate_zip(ZIP_FILE)) {
      stop(
        paste0(
          "Existing canonical Slide_4.zip failed validation:\n",
          ZIP_FILE,
          "\n\n",
          "The file may be incomplete or may not be the Dryad Slide_4 archive.\n",
          "For safety it was NOT deleted automatically."
        ),
        call. = FALSE
      )
    }

    message(
      "REUSE validated canonical archive: ",
      ZIP_FILE
    )
  }

  ##########################################################
  ## 6C. Try local Downloads / explicit local override
  ##########################################################

  if (!file.exists(ZIP_FILE)) {
    stage_local_zip_if_available()
  }

  ##########################################################
  ## 6D. Try direct Dryad download automatically
  ##########################################################

  if (!file.exists(ZIP_FILE)) {
    try_dryad_download()
  }

  ##########################################################
  ## 6E. Manual browser fallback
  ##########################################################

  if (!file.exists(ZIP_FILE) ||
      !validate_zip(ZIP_FILE)) {

    if (interactive()) {
      try(
        utils::browseURL(
          DRYAD_FILE_URL
        ),
        silent = TRUE
      )
    }

    stop(
      paste0(
        "Automatic Dryad download could not obtain a validated Slide_4.zip.\n\n",
        "This is usually caused by Dryad rejecting programmatic file-stream requests.\n\n",
        "Official Dryad dataset:\n",
        DRYAD_LANDING_URL,
        "\n\n",
        "Slide_4.zip direct browser link:\n",
        DRYAD_FILE_URL,
        "\n\n",
        "Dryad currently lists Slide_4.zip at approximately 3.68 GB.\n\n",
        "Option 1 — easiest:\n",
        "Download Slide_4.zip in your browser. On the next rerun, this script will\n",
        "also look automatically in ~/Downloads/Slide_4.zip.\n\n",
        "Option 2 — place it directly at:\n",
        ZIP_FILE,
        "\n\n",
        "Option 3 — point to an existing local file before rerunning:\n",
        "Sys.setenv(COSMX_SLIDE4_ZIP = 'D:/your/path/Slide_4.zip')\n\n",
        "No partial canonical raw data were installed."
      ),
      call. = FALSE
    )
  }

  ##########################################################
  ## 6F. Final pre-extraction ZIP gate
  ##########################################################

  if (!validate_zip(ZIP_FILE)) {
    stop(
      "Slide_4.zip did not pass the final pre-extraction validation.",
      call. = FALSE
    )
  }
}

############################################################
## 7. Robust extraction to exact frozen-script path
############################################################

if (!canonical_complete()) {

  EXTRACT_TMP <- file.path(
    BASE_DIR,
    ".Slide_4_extract_tmp"
  )

  if (dir.exists(EXTRACT_TMP)) {
    stop(
      paste0(
        "Temporary extraction directory already exists:\n",
        EXTRACT_TMP,
        "\nRemove this script-generated temporary directory and rerun."
      ),
      call. = FALSE
    )
  }

  dir_create(EXTRACT_TMP)

  extraction_success <- FALSE

  on.exit(
    {
      if (!extraction_success &&
          dir.exists(EXTRACT_TMP)) {
        unlink(
          EXTRACT_TMP,
          recursive = TRUE,
          force = TRUE
        )
      }
    },
    add = TRUE
  )

  message(
    "UNZIP validated Dryad archive:\n  ",
    ZIP_FILE,
    "\n-> ",
    EXTRACT_TMP
  )

  utils::unzip(
    ZIP_FILE,
    exdir = EXTRACT_TMP
  )

  ##########################################################
  ## Discover the extracted Run5611_MK3 directory robustly.
  ##########################################################

  candidate_dirs <- unique(
    dirname(
      list.files(
        EXTRACT_TMP,
        pattern = "^Run5611_MK3_exprMat_file\\.csv$",
        recursive = TRUE,
        full.names = TRUE
      )
    )
  )

  candidate_dirs <- candidate_dirs[
    basename(candidate_dirs) ==
      "Run5611_MK3"
  ]

  if (length(candidate_dirs) != 1L) {
    stop(
      paste0(
        "Could not uniquely locate Run5611_MK3 after extraction.\n",
        "Candidate count: ",
        length(candidate_dirs)
      ),
      call. = FALSE
    )
  }

  extracted_run_dir <- candidate_dirs[[1]]

  extracted_required <- file.path(
    extracted_run_dir,
    required_names
  )

  if (!all(file.exists(extracted_required))) {
    stop(
      "Extracted Run5611_MK3 directory is missing one or more required files.",
      call. = FALSE
    )
  }

  ##########################################################
  ## Install only after full extracted contract validation.
  ##########################################################

  dir_create(SLIDE_DIR)

  move_ok <- file.rename(
    extracted_run_dir,
    RUN_DIR
  )

  if (!isTRUE(move_ok)) {

    dir_create(RUN_DIR)

    items <- list.files(
      extracted_run_dir,
      full.names = TRUE,
      all.files = TRUE,
      no.. = TRUE
    )

    copy_ok <- file.copy(
      items,
      RUN_DIR,
      recursive = TRUE,
      overwrite = FALSE,
      copy.mode = TRUE,
      copy.date = TRUE
    )

    if (!all(copy_ok)) {
      stop(
        "Fallback copy of extracted Run5611_MK3 directory failed.",
        call. = FALSE
      )
    }
  }

  if (!canonical_complete()) {
    stop(
      "Canonical Run5611_MK3 installation failed post-copy validation.",
      call. = FALSE
    )
  }

  extraction_success <- TRUE

  ## Remove only this script's temporary extraction directory.
  if (dir.exists(EXTRACT_TMP)) {
    unlink(
      EXTRACT_TMP,
      recursive = TRUE,
      force = TRUE
    )
  }
}

############################################################
## 8. Final raw-data contract
############################################################

contract <- data.frame(
  Required_file = norm_path(required),
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
    "COSMX_SLIDE4_RAW_CONTRACT.csv"
  )
)

if (!all(contract$Exists)) {
  stop(
    "CosMx Slide 4 extraction does not match the frozen script contract.",
    call. = FALSE
  )
}

############################################################
## 9. Acquisition provenance
############################################################

provenance <- data.frame(
  Dataset = "CosMx melanoma RNA-SMI Slide 4",
  DOI = DRYAD_DOI,
  Landing_URL = DRYAD_LANDING_URL,
  File = "Slide_4.zip",
  File_URL = DRYAD_FILE_URL,
  Canonical_ZIP = norm_path(ZIP_FILE),
  Canonical_Run_Directory = norm_path(RUN_DIR),
  ZIP_exists = file.exists(ZIP_FILE),
  ZIP_size_bytes = ifelse(
    file.exists(ZIP_FILE),
    file.info(ZIP_FILE)$size,
    NA_real_
  ),
  Required_files_PASS = all(contract$Exists),
  stringsAsFactors = FALSE
)

write_utf8_csv(
  provenance,
  file.path(
    PROJECT_DIR,
    "release_metadata",
    "raw_data",
    "COSMX_SLIDE4_ACQUISITION_PROVENANCE.csv"
  )
)

############################################################
## 10. Final message
############################################################

cat("\n============================================================\n")
cat("03 COSMX DRYAD SLIDE-4 PREPARATION PASS\n")
cat("============================================================\n")
cat(
  "Required files: ",
  sum(contract$Exists),
  "/",
  nrow(contract),
  "\n",
  sep = ""
)
cat(
  "Expression matrix used by frozen scripts:\n",
  norm_path(required[[1]]),
  "\n",
  sep = ""
)
cat(
  "Metadata used by frozen scripts:\n",
  norm_path(required[[2]]),
  "\n",
  sep = ""
)
cat(
  "Canonical raw directory:\n",
  norm_path(RUN_DIR),
  "\n",
  sep = ""
)
cat(
  "NOTE: Slide_4_unzipped is NOT created; it was a redundant legacy copy.\n"
)
cat("============================================================\n")
