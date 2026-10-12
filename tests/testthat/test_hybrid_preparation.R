hybrid_preparation_fixture <- function() {
  root <- tempfile("hybrid-preparation-")
  dir.create(root)
  withr::defer(unlink(root, recursive = TRUE), envir = parent.frame())
  input <- file.path(root, "source.rds")
  output <- file.path(root, "prepared")
  source_dir <- file.path(root, "pxl")
  proximity_dir <- file.path(output, "proximity")
  dir.create(source_dir)
  dir.create(proximity_dir, recursive = TRUE)
  object <- precomputed_fixture()
  data <- precomputed_app_data(object)
  samples <- sort(unique(data$metadata$sample_id))
  sources <- setNames(file.path(source_dir, paste0(samples, ".pxl")), samples)
  # These tiny payloads exercise checksums only; no PXL/Parquet reader is invoked.
  for (sample in samples) saveRDS(list(synthetic_sample = sample), sources[[sample]])
  data$contract$proximity_source$pxl_by_sample <- sources
  object@misc$proxiomevis <- data$contract
  saveRDS(object, input)
  partitions <- sprintf("sample-%03d.parquet", seq_along(samples))
  for (partition in partitions) saveRDS(list(synthetic_partition = partition), file.path(proximity_dir, partition))
  manifest <- list(dataset_fingerprint = hybrid_dataset_fingerprint(data), partitions = data.frame(
    sample_id = samples, file = partitions,
    md5 = unname(tools::md5sum(file.path(proximity_dir, partitions))),
    source_pxl_md5 = unname(tools::md5sum(sources))))
  manifest_path <- file.path(proximity_dir, "manifest.rds")
  saveRDS(manifest, manifest_path)
  bundle_path <- file.path(output, "analysis.rds")
  saveRDS(list(existing_bundle = TRUE), bundle_path)

  environment <- new.env(parent = environment())
  repository <- normalizePath(file.path("..", ".."), mustWork = TRUE)
  withr::with_dir(repository, sys.source("scripts/prepare_hybrid_bundle.R", envir = environment))
  calls <- new.env(parent = emptyenv())
  calls$export <- calls$graph <- NULL
  calls$prepare <- 0L
  environment$read_precomputed_rds <- function(path) {
    stopifnot(identical(path, input))
    data
  }
  environment$write_precomputed_bundle <- function(data, path, assets) {
    calls$export <- list(path = path, assets = assets, dataset_id = data$contract$dataset_id)
    saveRDS(list(exported_bundle = TRUE, assets = assets), path)
  }
  environment$prepare_proximity_store <- function(...) {
    calls$prepare <- calls$prepare + 1L
    stop("Matching sidecars must be reused.")
  }
  environment$prepare_cell_graph_assets <- function(object, pxl_dir, output_dir, cells, dataset_id) {
    calls$graph <- list(is_seurat = inherits(object, "Seurat"), pxl_dir = pxl_dir,
                        output_dir = output_dir, cells = cells, dataset_id = dataset_id)
  }
  environment$requireNamespace <- function(...) TRUE
  list(input = input, output = output, source_dir = source_dir, sources = sources,
       partitions = file.path(proximity_dir, partitions), manifest = manifest,
       manifest_path = manifest_path, bundle_path = bundle_path, calls = calls,
       build = environment$build_hybrid_bundle)
}

test_that("sourcing preparation respects the CLI guard and reuses matching sidecars", {
  fixture <- hybrid_preparation_fixture()
  expect_true(is.function(fixture$build))
  original <- tools::md5sum(c(fixture$input, fixture$sources, fixture$partitions, fixture$manifest_path))
  expect_message(path <- fixture$build(fixture$input, fixture$output), "Prepared application bundle")
  expect_identical(path, fixture$bundle_path)
  expect_identical(fixture$calls$prepare, 0L)
  expect_identical(fixture$calls$export$assets, list(
    proximity = list(path = "proximity", manifest = "manifest.rds"),
    graphs = list(path = "graphs", index = "index.rds")))
  expect_true(fixture$calls$graph$is_seurat)
  expect_identical(fixture$calls$graph$pxl_dir, fixture$source_dir)
  expect_identical(fixture$calls$graph$output_dir, file.path(fixture$output, "graphs"))
  expect_identical(fixture$calls$graph$cells, colnames(readRDS(fixture$input)))
  expect_identical(fixture$calls$graph$dataset_id, "fixture-v1")
  expect_identical(tools::md5sum(names(original)), original)
})

test_that("selected graph cells are forwarded without restricting other bundle assets", {
  fixture <- hybrid_preparation_fixture()
  selected <- c("cell2", "cell1")
  original <- tools::md5sum(c(fixture$input, fixture$sources, fixture$partitions, fixture$manifest_path))
  expect_message(fixture$build(fixture$input, fixture$output, graph_cells = selected), "Prepared application bundle")
  expect_identical(fixture$calls$graph$cells, selected)
  expect_identical(fixture$calls$prepare, 0L)
  expect_identical(fixture$calls$export$assets$proximity, list(path = "proximity", manifest = "manifest.rds"))
  expect_identical(tools::md5sum(names(original)), original)
})

test_that("invalid graph selections reject before writes and cannot conflict with graphs FALSE", {
  fixture <- hybrid_preparation_fixture()
  original <- tools::md5sum(c(fixture$input, fixture$bundle_path, fixture$sources,
                              fixture$partitions, fixture$manifest_path))
  for (selection in list(character(), c("cell1", "cell1"), "unknown", NA_character_, "", 42)) {
    expect_error(fixture$build(fixture$input, fixture$output, graph_cells = selection), "unique exact cell IDs")
  }
  expect_error(fixture$build(fixture$input, fixture$output, graphs = FALSE, graph_cells = "cell1"),
               "requires graphs = TRUE")
  expect_null(fixture$calls$export)
  expect_null(fixture$calls$graph)
  expect_identical(fixture$calls$prepare, 0L)
  expect_identical(tools::md5sum(names(original)), original)
})

test_that("the optional fourth CLI argument is read and validated before preparation", {
  fixture <- hybrid_preparation_fixture()
  selected_file <- file.path(dirname(fixture$input), "graph_cells.txt")
  writeLines("unknown-cell", selected_file)
  original <- tools::md5sum(c(fixture$input, fixture$bundle_path, fixture$sources, fixture$manifest_path))
  repository <- normalizePath(file.path("..", ".."), mustWork = TRUE)
  result <- callr::rscript(file.path(repository, "scripts", "prepare_hybrid_bundle.R"),
    cmdargs = c(fixture$input, fixture$output, fixture$source_dir, selected_file),
    wd = repository, env = c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", ARROW_NUM_THREADS = "1"),
    fail_on_status = FALSE, show = FALSE)
  expect_identical(result$status, 1L)
  expect_match(result$stderr, "unique exact cell IDs")
  expect_identical(tools::md5sum(names(original)), original)
})

test_that("preparation supports a proximity-only bundle without graph calls", {
  fixture <- hybrid_preparation_fixture()
  expect_message(fixture$build(fixture$input, fixture$output, fixture$source_dir, graphs = FALSE),
                   "Prepared application bundle")
  expect_null(fixture$calls$graph)
  expect_identical(fixture$calls$prepare, 0L)
  expect_identical(fixture$calls$export$assets, list(proximity = list(path = "proximity", manifest = "manifest.rds")))
  expect_true(readRDS(fixture$bundle_path)$exported_bundle)
})

test_that("legacy graph sidecars reject before overwriting the application bundle", {
  fixture <- hybrid_preparation_fixture()
  graph_dir <- file.path(fixture$output, "graphs")
  dir.create(graph_dir)
  graph_index <- file.path(graph_dir, "index.rds")
  saveRDS(list(schema_version = 1L), graph_index)
  original <- tools::md5sum(c(fixture$bundle_path, fixture$manifest_path, graph_index))
  expect_error(fixture$build(fixture$input, fixture$output), "new directory")
  expect_null(fixture$calls$export)
  expect_null(fixture$calls$graph)
  expect_identical(fixture$calls$prepare, 0L)
  expect_identical(tools::md5sum(names(original)), original)
})

test_that("graph preparation failure does not publish a new application bundle", {
  fixture <- hybrid_preparation_fixture()
  environment(fixture$build)$prepare_cell_graph_assets <- function(...) stop("Graph export failed.")
  original <- tools::md5sum(fixture$bundle_path)
  expect_error(fixture$build(fixture$input, fixture$output), "Graph export failed")
  expect_null(fixture$calls$export)
  expect_identical(tools::md5sum(fixture$bundle_path), original)
})

test_that("sidecar and source mismatches fail before export without rewriting files", {
  for (failure in c("dataset", "partition", "source")) {
    fixture <- hybrid_preparation_fixture()
    manifest <- fixture$manifest
    if (failure == "dataset") manifest$dataset_fingerprint <- "different-dataset"
    if (failure == "partition") manifest$partitions$md5[[1L]] <- "different-partition"
    if (failure == "source") manifest$partitions$source_pxl_md5[[1L]] <- "different-source"
    saveRDS(manifest, fixture$manifest_path)
    original <- tools::md5sum(c(fixture$input, fixture$bundle_path, fixture$sources,
                                fixture$partitions, fixture$manifest_path))
    expected <- switch(failure, dataset = "different dataset", partition = "checksum validation", source = "source PXL files changed")
    expect_error(fixture$build(fixture$input, fixture$output), expected)
    expect_null(fixture$calls$export)
    expect_null(fixture$calls$graph)
    expect_identical(fixture$calls$prepare, 0L)
    expect_identical(tools::md5sum(names(original)), original)
  }
})
