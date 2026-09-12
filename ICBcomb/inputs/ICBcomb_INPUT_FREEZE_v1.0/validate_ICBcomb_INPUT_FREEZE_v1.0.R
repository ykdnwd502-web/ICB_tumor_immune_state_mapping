## ICBcomb input freeze v1.0 validator
## Run from the root of the extracted freeze package.

options(stringsAsFactors = FALSE)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
manifest_file <- file.path(root, "SHA256SUMS.csv")

if (!file.exists(manifest_file)) {
  stop("SHA256SUMS.csv not found. Set working directory to the freeze-package root.")
}

manifest <- read.csv(manifest_file, check.names = FALSE, stringsAsFactors = FALSE)

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for SHA256 validation.")
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

check <- manifest
check$exists <- file.exists(file.path(root, check$relative_path))
check$sha256_observed <- vapply(
  file.path(root, check$relative_path),
  function(x) if (file.exists(x)) sha256_file(x) else NA_character_,
  character(1)
)
check$pass <- check$exists & check$sha256_observed == check$sha256

cat("============================================================\n")
cat("ICBcomb INPUT FREEZE v1.0 VALIDATION\n")
cat("Files:", sum(check$pass), "/", nrow(check), "PASS\n")
cat("Blocking failures:", sum(!check$pass), "\n")
cat("Overall:", ifelse(all(check$pass), "PASS", "FAIL"), "\n")
cat("============================================================\n")

if (!all(check$pass)) {
  print(check[!check$pass, , drop = FALSE])
  stop("ICBcomb input freeze validation failed.")
}
