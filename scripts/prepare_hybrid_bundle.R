#!/usr/bin/env Rscript
# Run from the repository after prepare_precomputed_rds.R; sources are read-only.
source(file.path("R", "precomputed_contract.R"))
source(file.path("R", "proximity_store.R"))
source(file.path("R", "cell_graph_assets.R"))

build_hybrid_bundle <- function(input_rds, output_dir, pxl_dir = NULL, graphs = TRUE, graph_cells = NULL) {
  data <- read_precomputed_rds(input_rds)
  if (!is.null(graph_cells)) {
    if (!isTRUE(graphs)) stop("graph_cells requires graphs = TRUE.", call. = FALSE)
    if (!is.character(graph_cells) || !length(graph_cells) || anyNA(graph_cells) ||
        any(!nzchar(graph_cells)) || anyDuplicated(graph_cells) ||
        !all(graph_cells %in% data$metadata$component)) {
      stop("graph_cells must contain nonempty, unique exact cell IDs from the source dataset.", call. = FALSE)
    }
  }
  sources <- data$contract$proximity_source$pxl_by_sample
  if (is.null(sources)) stop("The source bundle has no PXL preparation manifest.", call. = FALSE)
  if (is.null(pxl_dir)) {
    directories <- unique(dirname(sources))
    if (length(directories) != 1L) stop("Specify a directory containing the source PXL files.", call. = FALSE)
    pxl_dir <- directories
  }
  sources <- stats::setNames(file.path(pxl_dir, basename(sources)), names(sources))
  if (!all(file.exists(sources))) stop("A source PXL file is missing from the preparation directory.", call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  path <- file.path(output_dir, "analysis.rds")
  if (identical(normalizePath(input_rds, mustWork = TRUE), path)) {
    stop("Write the application bundle separately from the source RDS.", call. = FALSE)
  }
  graph_dir <- file.path(output_dir, "graphs")
  if (graphs && file.exists(file.path(graph_dir, "index.rds"))) {
    read_cell_graph_index(graph_dir, data$contract$dataset_id, data$metadata$component)
  }
  proximity_dir <- file.path(output_dir, "proximity")
  manifest_path <- file.path(proximity_dir, "manifest.rds")
  if (file.exists(manifest_path)) {
    manifest <- readRDS(manifest_path)
    if (!identical(manifest$dataset_fingerprint, hybrid_dataset_fingerprint(data))) {
      stop("Existing proximity assets belong to a different dataset. Choose a new output directory.", call. = FALSE)
    }
    paths <- vapply(manifest$partitions$file, function(file) proximity_store_file(proximity_dir, file), character(1))
    if (!identical(unname(tools::md5sum(paths)), manifest$partitions$md5)) {
      stop("Existing proximity assets failed checksum validation.", call. = FALSE)
    }
    source_hashes <- unname(tools::md5sum(sources[manifest$partitions$sample_id]))
    if (!identical(source_hashes, manifest$partitions$source_pxl_md5)) {
      stop("The source PXL files changed since proximity preparation. Prepare a new, matching analysis bundle.", call. = FALSE)
    }
  } else {
    prepare_proximity_store(data, sources, proximity_dir)
  }
  assets <- list(proximity = list(path = "proximity", manifest = "manifest.rds"))
  if (graphs) assets$graphs <- list(path = "graphs", index = "index.rds")
  if (graphs) {
    if (!requireNamespace("pixelatorR", quietly = TRUE)) stop("Install pixelatorR for graph preparation.")
    object <- readRDS(input_rds)
    if (!inherits(object, "Seurat")) stop("Graph preparation needs the source Seurat RDS with its exact FSMap.")
    prepare_cell_graph_assets(object, pxl_dir, graph_dir,
      cells = if (is.null(graph_cells)) colnames(object) else graph_cells,
      dataset_id = data$contract$dataset_id)
  }
  write_precomputed_bundle(data, path, assets = assets)
  message("Prepared application bundle: ", path)
  invisible(path)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) %in% c(2L, 3L, 4L)) {
    stop("Usage: Rscript scripts/prepare_hybrid_bundle.R PRECOMPUTED_RDS OUTPUT_DIR [PXL_DIR [GRAPH_CELLS_TXT]]", call. = FALSE)
  }
  build_hybrid_bundle(args[[1L]], args[[2L]], if (length(args) >= 3L) args[[3L]] else NULL,
    graph_cells = if (length(args) == 4L) readLines(args[[4L]], warn = FALSE) else NULL)
}
