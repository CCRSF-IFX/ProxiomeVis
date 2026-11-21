#!/usr/bin/env Rscript
# Use trusted application bundles. The same engine runs in the app and offline.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
for (file in c("precomputed_contract.R", "preparation_cli.R", "proximity_store.R",
               "hybrid_analysis.R", "hybrid_module.R", "bundle_builder.R")) {
  source(file.path(root, "R", file))
}
usage <- paste(
  "Usage: Rscript scripts/run_comparison.R BUNDLE_RDS REQUEST_JSON NEW_COMPARISON_RDS [options]",
  "  --id TEXT             Unique analysis ID (default: an ID from the request and time)",
  "  --label TEXT          Display label (default: target vs reference)",
  "  --output-bundle RDS   Also append the result to a new application RDS beside the input bundle",
  sep = "\n")
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0L) }
parsed <- parse_preparation_cli(args, list(id = NULL, label = NULL, "output-bundle" = NULL), 3L, usage)
bundle <- normalizePath(parsed$positional[[1L]], mustWork = TRUE)
request_path <- normalizePath(parsed$positional[[2L]], mustWork = TRUE)
output <- parsed$positional[[3L]]
if (tolower(tools::file_ext(output)) != "rds" || file.exists(output)) {
  stop("Choose a new comparison .rds file; existing files are not overwritten.", call. = FALSE)
}
bundle_output <- parsed$options[["output-bundle"]]
if (!is.null(bundle_output) && (file.exists(bundle_output) ||
    tolower(tools::file_ext(bundle_output)) != "rds")) stop("Choose a new application .rds output.", call. = FALSE)
request <- jsonlite::read_json(request_path, simplifyVector = TRUE)
result <- run_hybrid_analysis(bundle, request)
request <- result$provenance$request
invisible(gc())
data <- read_precomputed_rds(bundle)
id <- parsed$options$id
if (is.null(id)) id <- paste0("prepared-", substr(digest::digest(list(request, result$provenance$created_utc),
  algo = "sha256"), 1L, 16L))
label <- parsed$options$label
if (is.null(label)) label <- paste(request$target, "vs", request$reference)
entry <- list(id = id, label = label, request = request, result = result)
adapted <- data
adapted$contract <- hybrid_result_contract(data, entry)
adapted$contract$defaults[c("sample_column", "donor_column", "condition_column", "celltype_column", "cluster_column")] <-
  list("sample_id", "donor", "condition", "celltype", "cluster")
validate_precomputed_data(adapted)
parent <- normalizePath(dirname(output), mustWork = TRUE)
output <- file.path(parent, basename(output))
temporary <- tempfile(".proxiomevis-comparison-", tmpdir = parent)
tryCatch({
  saveRDS(list(result_version = 1L, entry = entry), temporary, compress = "xz")
  if (file.exists(output) || !file.rename(temporary, output)) stop("Could not publish the comparison RDS.", call. = FALSE)
}, finally = if (file.exists(temporary)) unlink(temporary))
message("Prepared comparison: ", output)
if (!is.null(bundle_output)) {
  rm(data, adapted, result, entry, request)
  invisible(gc())
  import_user_comparison(bundle, output, bundle_output)
  message("Application bundle with appended results: ", bundle_output)
}
