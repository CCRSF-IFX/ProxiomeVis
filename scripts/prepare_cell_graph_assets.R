#!/usr/bin/env Rscript
# Writes one node-only Parquet per sample; PXL files are needed only during preparation.
# An optional cells.txt extends existing sample files atomically (use CPU-safe batches).
source(file.path("R", "cell_graph_assets.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L)) {
  stop("Usage: Rscript scripts/prepare_cell_graph_assets.R source.rds pxl_directory output_directory dataset_id [cells.txt]", call. = FALSE)
}
if (!requireNamespace("pixelatorR", quietly = TRUE) || !requireNamespace("arrow", quietly = TRUE)) {
  stop("Install pixelatorR and arrow for preparation.")
}
invisible(loadNamespace("pixelatorR"))
object <- readRDS(args[[1]])
cells <- if (length(args) == 5L) readLines(args[[5]], warn = FALSE) else colnames(object)
prepare_cell_graph_assets(object, args[[2]], args[[3]], cells = cells, dataset_id = args[[4]])
