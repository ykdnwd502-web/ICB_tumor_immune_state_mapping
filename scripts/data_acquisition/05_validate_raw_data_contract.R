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
## 05_validate_raw_data_contract.R
##
## Final pre-analysis validation only.
## Does not download or modify source files.
############################################################

contract <- data.frame(
  Dataset = character(),
  Required_path = character(),
  Required = logical(),
  Purpose = character(),
  stringsAsFactors = FALSE
)

add_req <- function(dataset, path, purpose, required = TRUE) {
  contract <<- rbind(
    contract,
    data.frame(
      Dataset = dataset,
      Required_path = norm_path(path),
      Required = required,
      Purpose = purpose,
      stringsAsFactors = FALSE
    )
  )
}

## GEO discovery / external ICB
add_req(
  "GSE244982",
  file.path(PROJECT_DIR, "data_raw", "GSE244982",
            "GSE244982_ProcessedData_bulkRNAseq.txt.gz"),
  "Bulk discovery expression input"
)

add_req(
  "GSE244983",
  file.path(PROJECT_DIR, "data_raw", "GSE244983",
            "GSE244983_RawCounts_scRNAseq.txt.gz"),
  "scRNA raw counts"
)

add_req(
  "GSE244983",
  file.path(PROJECT_DIR, "data_raw", "GSE244983",
            "GSE244983_SingleCellAnnotations.txt.gz"),
  "Author single-cell annotation"
)

add_req(
  "GSE244983",
  file.path(PROJECT_DIR, "data_raw", "GSE244983",
            "GSE244983_NormalizedData_scRNAseq.txt.gz"),
  "Source-completeness copy; not the primary reconstruction input",
  required = FALSE
)

add_req(
  "GSE78220",
  file.path(PROJECT_DIR, "data_raw", "external_validation",
            "GSE78220", "GSE78220_PatientFPKM.xlsx"),
  "External ICB response/survival expression input"
)

add_req(
  "GSE78220",
  file.path(PROJECT_DIR, "data_raw", "external_validation",
            "GSE78220", "GSE78220_pData_raw_from_GEO.csv"),
  "GEO phenotype provenance",
  required = FALSE
)

add_req(
  "GSE91061",
  file.path(PROJECT_DIR, "data_raw", "external_validation",
            "GSE91061",
            "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz"),
  "Pretreatment RLD expression used by frozen analysis"
)

add_req(
  "GSE91061",
  file.path(PROJECT_DIR, "data_raw", "external_validation",
            "GSE91061", "GSE91061_pData_raw_from_GEO.csv"),
  "GEO phenotype metadata"
)

## Primary 10x Visium
for (f in c(
  file.path(PROJECT_DIR, "data_raw",
            "spatial_10x_human_melanoma_IF_FFPE",
            "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
  file.path(PROJECT_DIR, "data_raw",
            "spatial_10x_human_melanoma_IF_FFPE",
            "filtered_feature_bc_matrix", "features.tsv.gz"),
  file.path(PROJECT_DIR, "data_raw",
            "spatial_10x_human_melanoma_IF_FFPE",
            "filtered_feature_bc_matrix", "matrix.mtx.gz"),
  file.path(PROJECT_DIR, "data_raw",
            "spatial_10x_human_melanoma_IF_FFPE",
            "spatial", "tissue_positions.csv"),
  file.path(PROJECT_DIR, "data_raw",
            "spatial_10x_human_melanoma_IF_FFPE",
            "spatial", "scalefactors_json.json")
)) {
  add_req(
    "10x_Visium_Human_Melanoma_FFPE",
    f,
    "Primary Visium raw contract"
  )
}

## CosMx
for (f in c(
  "Run5611_MK3_exprMat_file.csv",
  "Run5611_MK3_metadata_file.csv",
  "Run5611_MK3_fov_positions_file.csv"
)) {
  add_req(
    "CosMx_Dryad_Slide4",
    file.path(
      PROJECT_DIR,
      "data_raw",
      "CosMx_melanocytic_tumors_Dryad",
      "Slide_4",
      "Run5611_MK3",
      f
    ),
    "CosMx raw expression/spatial metadata contract"
  )
}

## GSE250636
samples <- c(
  "GSM7983358_sample2",
  "GSM7983359_sample4",
  "GSM7983360_sample7",
  "GSM7983361_sample8",
  "GSM7983362_sample12",
  "GSM7983363_sample13",
  "GSM7983364_sample14",
  "GSM7983365_sample15",
  "GSM7983366_sample16"
)

gse250_dir <- file.path(
  PROJECT_DIR,
  "data_external_spatial_18J",
  "GSE250636",
  "unpacked_recursive",
  "GSE250636_RAW"
)

for (s in samples) {
  add_req(
    "GSE250636",
    file.path(
      gse250_dir,
      paste0(s, "_filtered_feature_bc_matrix.h5")
    ),
    "External Visium recurrence filtered H5"
  )

  add_req(
    "GSE250636",
    file.path(
      gse250_dir,
      paste0(s, "_tissue_positions_list.csv.gz")
    ),
    "External Visium recurrence coordinates"
  )
}

## Thrane
thrane_count_dir <- file.path(
  PROJECT_DIR,
  "data_external_spatial_18J",
  "Thrane2018_legacyST",
  "unpacked",
  "Thrane_et_al_2018_CAN_RES",
  "Count_data"
)

for (f in c(
  "ST_mel1_rep1_counts.tsv",
  "ST_mel1_rep2_counts.tsv",
  "ST_mel2_rep1_counts.tsv",
  "ST_mel2_rep2_counts.tsv",
  "ST_mel3_rep1_counts.tsv",
  "ST_mel3_rep2_counts.tsv",
  "ST_mel4_rep1_counts.tsv",
  "ST_mel4_rep2_counts.tsv"
)) {
  add_req(
    "Thrane2018_legacyST",
    file.path(thrane_count_dir, f),
    "External legacy-ST recurrence count matrix"
  )
}

contract$Exists <- file.exists(contract$Required_path)
contract$Size_bytes <- ifelse(
  contract$Exists,
  file.info(contract$Required_path)$size,
  NA_real_
)

contract$PASS <- ifelse(
  contract$Required,
  contract$Exists & contract$Size_bytes > 0,
  TRUE
)

out_csv <- file.path(
  PROJECT_DIR,
  "release_metadata",
  "raw_data",
  "RAW_DATA_CONTRACT.csv"
)

write_utf8_csv(
  contract,
  out_csv
)

required_n <- sum(contract$Required)
required_pass <- sum(contract$PASS[contract$Required])

summary_lines <- c(
  "============================================================",
  "RAW DATA ACQUISITION CONTRACT",
  "============================================================",
  paste0("Project: ", PROJECT_DIR),
  paste0("Required entries PASS: ", required_pass, "/", required_n),
  paste0(
    "Optional entries present: ",
    sum(contract$Exists[!contract$Required]),
    "/",
    sum(!contract$Required)
  ),
  "",
  paste0(
    "Decision: ",
    if (required_pass == required_n)
      "READY_FOR_ANALYSIS_CHAIN"
    else
      "BLOCKED"
  ),
  "============================================================"
)

writeLines(
  summary_lines,
  file.path(
    PROJECT_DIR,
    "release_metadata",
    "raw_data",
    "RAW_DATA_CONTRACT_SUMMARY.txt"
  ),
  useBytes = TRUE
)

cat(
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n"
)

if (required_pass != required_n) {
  bad <- contract[
    contract$Required &
      !contract$PASS,
    ,
    drop = FALSE
  ]

  print(
    bad[
      ,
      c(
        "Dataset",
        "Required_path",
        "Exists",
        "Size_bytes"
      )
    ]
  )

  stop(
    "Raw-data contract is incomplete. Inspect RAW_DATA_CONTRACT.csv.",
    call. = FALSE
  )
}
