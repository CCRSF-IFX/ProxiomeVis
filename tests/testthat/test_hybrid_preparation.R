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

test_that("version 2 graph upgrade rejects before proximity preparation or bundle writes", {
  fixture <- hybrid_preparation_fixture()
  graph_dir <- file.path(fixture$output, "graphs")
  dir.create(graph_dir)
  graph_index <- file.path(graph_dir, "index.rds")
  saveRDS(list(schema_version = 2L), graph_index)
  environment(fixture$build)$read_cell_graph_index <- function(...) list(schema_version = 2L)
  original <- tools::md5sum(c(fixture$bundle_path, fixture$manifest_path, graph_index))
  expect_error(fixture$build(fixture$input, fixture$output), "new directory")
  expect_null(fixture$calls$export)
  expect_null(fixture$calls$graph)
  expect_identical(fixture$calls$prepare, 0L)
  expect_identical(tools::md5sum(names(original)), original)
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

test_that("ranked export writes only selected graph cells and accurate per-sample logs", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("arrow"); skip_if_not_installed("pixelatorR")
  fixture <- hybrid_preparation_fixture()
  env <- environment(fixture$build)
  data <- env$read_precomputed_rds(fixture$input)
  meta <- data$metadata
  manifest <- fixture$manifest
  manifest$schema_version <- 1L
  manifest$dataset_id <- data$contract$dataset_id
  for (i in seq_along(fixture$partitions)) {
    cells <- meta$component[meta$sample_id == manifest$partitions$sample_id[[i]]]
    arrow::write_parquet(rbind(
      data.frame(component = cells, marker_1 = "A", marker_2 = "B", log2_ratio = c(1, 2, NA, NA)),
      data.frame(component = cells, marker_1 = "A", marker_2 = "C", log2_ratio = seq_along(cells)),
      data.frame(component = cells, marker_1 = "A", marker_2 = "A", log2_ratio = seq_along(cells))),
      fixture$partitions[[i]])
  }
  manifest$partitions$md5 <- unname(tools::md5sum(fixture$partitions))
  saveRDS(manifest, fixture$manifest_path)
  mapping <- data.frame(pxl_file = unname(fixture$sources))
  mapping$id_map <- lapply(names(fixture$sources), function(sample) {
    cells <- meta$component[meta$sample_id == sample]
    data.frame(current_id = cells, original_id = cells)
  })
  db <- list(close = function() invisible(NULL), query = function(sql) {
    if (startsWith(sql, "SET ")) return(invisible(NULL))
    if (startsWith(sql, "SELECT DISTINCT")) return(data.frame(layout = "coarsened_pmds_3d", graph_projection = "full"))
    cell <- sub(".*SELECT '(cell[0-9]+)' AS cell_id.*", "\\1", sql)
    stopifnot(cell %in% meta$component)
    nodes <- data.frame(cell_id = cell, node_id = c("1-umi1", "2-umi2"), x = c(0, 1), y = c(1, 0), z = c(0, 1))
    nodes$markers <- list("A", "C")
    nodes
  })
  testthat::local_mocked_bindings(FSMap = function(object) mapping,
    PixelDB = list(new = function(path) db), .package = "pixelatorR")
  env$prepare_cell_graph_assets <- prepare_cell_graph_assets
  env$write_precomputed_bundle <- write_precomputed_bundle
  original <- tools::md5sum(c(fixture$input, fixture$sources, fixture$partitions, fixture$manifest_path))
  expect_message(fixture$build(fixture$input, fixture$output, graph_markers = c("A", "B")), "3D export")
  selection <- readRDS(file.path(fixture$output, "graph_selection.rds"))
  expect_length(selection$cells, 8L)
  expect_true(all(selection$selected$colocalization_top & selection$selected$colocalization_bottom))
  expect_identical(selection$parameters$markers, c("A", "B"))
  summary <- read.delim(file.path(fixture$output, "graph_export_summary.tsv"))
  expect_equal(summary$total_cells, rep(4, 4))
  expect_equal(summary$selected_cells, rep(2, 4))
  expect_identical(summary$exported_cells, summary$selected_cells)
  expect_setequal(read.delim(file.path(fixture$output, "graph_selected_cells.tsv"))$component, selection$cells)
  index <- read_cell_graph_index(file.path(fixture$output, "graphs"), data$contract$dataset_id, meta$component)
  expect_setequal(index$cells$cell_id, selection$cells)
  for (file in index$samples$file) {
    nodes <- arrow::read_parquet(file.path(fixture$output, "graphs", file))
    expect_true(all(nodes$cell_id %in% selection$cells))
    expect_false("node_id" %in% names(nodes))
    expect_true("C" %in% unlist(nodes$markers))
    expect_equal(nrow(nodes), 2L * length(unique(nodes$cell_id)))
  }
  expect_equal(nrow(read_precomputed_rds(fixture$bundle_path)$metadata), nrow(meta))
  expect_identical(tools::md5sum(names(original)), original)
  expect_message(fixture$build(fixture$input, fixture$output, graph_markers = c("A", "B")), "Reusing sample")
  expect_error(fixture$build(fixture$input, fixture$output), "selected graph export")
  expect_message(fixture$build(fixture$input, fixture$output, graph_markers = c("A", "B", "C")), "Prepared application bundle")
  bundle_hash <- tools::md5sum(fixture$bundle_path)
  expect_error(fixture$build(fixture$input, fixture$output, graph_markers = c("A", "B")), "outside the ranked selection")
  expect_identical(tools::md5sum(fixture$bundle_path), bundle_hash)
})

test_that("marker options reject conflicting or invalid requests before writing outputs", {
  fixture <- hybrid_preparation_fixture()
  expect_error(fixture$build(fixture$input, fixture$output, graph_markers = "A"), "at least two distinct")
  expect_error(fixture$build(fixture$input, fixture$output, graph_markers = c("A", "unknown")), "Unknown graph markers")
  expect_error(fixture$build("missing", fixture$output, graph_markers = c("A", "B"), graphs = FALSE), "requires graphs")
  expect_error(fixture$build("missing", fixture$output, graph_markers = c("A", "B"), graph_cells = "cell1"), "cannot be combined")
  expect_null(fixture$calls$graph)
  expect_null(fixture$calls$export)
})

test_that("the marker CLI parses lists and validates them before changing outputs", {
  fixture <- hybrid_preparation_fixture()
  repository <- normalizePath(file.path("..", ".."), mustWork = TRUE)
  original <- tools::md5sum(c(fixture$input, fixture$bundle_path, fixture$manifest_path))
  cases <- c("A, unknown" = "Unknown graph markers: unknown", "A,A" = "at least two distinct",
    "A,B," = "missing after the final comma", "A,,B" = "at least two distinct")
  for (markers in names(cases)) {
    result <- callr::rscript(file.path(repository, "scripts", "prepare_hybrid_bundle.R"),
      cmdargs = c(fixture$input, fixture$output, "--graph-markers", markers),
      wd = repository, env = c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", ARROW_NUM_THREADS = "1"),
      fail_on_status = FALSE, show = FALSE)
    expect_identical(result$status, 1L)
    expect_match(result$stderr, cases[[markers]])
  }
  expect_identical(tools::md5sum(names(original)), original)
})
