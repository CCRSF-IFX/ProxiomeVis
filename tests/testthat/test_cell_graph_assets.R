source(file.path("..", "..", "R", "cell_graph_assets.R"))

graph_asset_fixture <- function(cell = "cell1") {
  nodes <- data.frame(cell_id = cell, node_id = paste0("node", 1:4), x = as.double(1:4),
                      y = c(0, 1, 0, 1), z = c(1, 1, 2, 2))
  nodes$markers <- list(character(), "CD81", "CD82", c("CD81", "CD82"))
  list(schema_version = 2L, dataset_id = "fixture-v1", cell_id = cell,
       source_cell_id = sub("cell", "original", cell), sample = "sample1", source_pxl = "sample.pxl",
       source_pxl_md5 = paste(rep("a", 32), collapse = ""), layout = "coarsened_pmds_3d",
       graph_projection = "full", nodes = nodes,
       interpretation = "Abstract graph coordinates, not physical cell positions; all source layout nodes retained.")
}

write_graph_fixture <- function(root, assets = list(graph_asset_fixture())) {
  temporary <- tempfile(tmpdir = root)
  schema <- cell_graph_schema()
  sink <- arrow::FileOutputStream$create(temporary)
  writer <- arrow::ParquetFileWriter$create(schema, sink, arrow::ParquetWriterProperties$create(schema$names))
  for (asset in assets) writer$WriteTable(arrow::Table$create(asset$nodes, schema = schema), nrow(asset$nodes))
  writer$Close()
  sink$close()
  md5 <- unname(tools::md5sum(temporary))
  file <- paste0(md5, ".parquet")
  file.rename(temporary, file.path(root, file))
  first <- assets[[1]]
  cells <- do.call(rbind, lapply(assets, function(asset) data.frame(
    cell_id = asset$cell_id, source_cell_id = asset$source_cell_id, sample = asset$sample, n_nodes = nrow(asset$nodes))))
  index <- list(schema_version = 2L, dataset_id = "fixture-v1", n_total_cells = 2L,
                coverage = if (length(assets) == 2L) "all_cells" else "selected_cells", cells = cells,
                samples = data.frame(sample = first$sample, source_pxl = first$source_pxl,
                  source_pxl_md5 = first$source_pxl_md5, layout = first$layout, graph_projection = first$graph_projection,
                  file = file, md5 = md5, size_bytes = file.info(file.path(root, file))$size,
                  n_nodes = sum(cells$n_nodes), n_cells = nrow(cells)))
  saveRDS(index, file.path(root, "index.rds"))
  read_cell_graph_index(root, "fixture-v1")
}

testthat::test_that("Arrow filters the selected cell from one sample file and enforces limits", {
  root <- withr::local_tempdir()
  index <- write_graph_fixture(root, list(graph_asset_fixture(), graph_asset_fixture("cell2")))
  testthat::expect_equal(load_cell_graph_asset(index, root, "cell1"), graph_asset_fixture())
  testthat::expect_equal(load_cell_graph_asset(index, root, "cell2"), graph_asset_fixture("cell2"))
  testthat::expect_length(list.files(root, pattern = "parquet$"), 1L)
  reader <- arrow::ParquetFileReader$create(file.path(root, index$samples$file))
  testthat::expect_equal(reader$num_row_groups, 2L)
  testthat::expect_true(reader$GetSchema()$Equals(cell_graph_schema(), check_metadata = FALSE))
  testthat::expect_error(load_cell_graph_asset(index, root, "missing"), "No prepared")
  testthat::expect_error(load_cell_graph_asset(index, root, "cell1", max_nodes = 2), "no nodes were silently removed")
  index$cells$n_nodes[1] <- 1L
  testthat::expect_error(load_cell_graph_asset(index, root, "cell1", max_nodes = 2), "no nodes were silently removed")
  testthat::expect_error(load_cell_graph_asset(index, root, "cell1"), "row count")
  testthat::expect_error(load_cell_graph_asset(index, root, "cell1", max_nodes = NA_real_), "positive integer")
})

testthat::test_that("trusted indexes reject legacy schemas, bad totals, unknown cells and escaped paths", {
  root <- withr::local_tempdir()
  outside <- withr::local_tempdir()
  index <- write_graph_fixture(root)
  testthat::expect_error(read_cell_graph_index(root, "another-dataset"), "does not match")
  testthat::expect_error(read_cell_graph_index(root, "fixture-v1", "cell2"), "unknown cells")
  testthat::expect_error(cell_graph_path(root, "../anything.parquet"), "filename")
  testthat::expect_error(cell_graph_path(root, file.path(outside, index$samples$file)), "filename")
  if (.Platform$OS.type != "windows") {
    file.symlink(file.path(root, index$samples$file), file.path(outside, index$samples$file))
    testthat::expect_error(cell_graph_path(outside, index$samples$file), "outside")
  }
  index$samples$n_nodes <- 9
  saveRDS(index, file.path(root, "index.rds"))
  testthat::expect_error(read_cell_graph_index(root, "fixture-v1"), "totals")
  saveRDS(list(schema_version = 1L), file.path(root, "index.rds"))
  testthat::expect_error(read_cell_graph_index(root, "fixture-v1"), "new directory")
})

testthat::test_that("sample checksums are cached per index and invalidated after a file change", {
  root <- withr::local_tempdir()
  index <- write_graph_fixture(root)
  checks <- 0L
  hash <- tools::md5sum
  testthat::local_mocked_bindings(md5sum = function(files) { checks <<- checks + 1L; hash(files) }, .package = "tools")
  load_cell_graph_asset(index, root, "cell1")
  load_cell_graph_asset(index, root, "cell1")
  testthat::expect_equal(checks, 1L)
  file <- file.path(root, index$samples$file)
  Sys.setFileTime(file, Sys.time() + 10)
  load_cell_graph_asset(index, root, "cell1")
  testthat::expect_equal(checks, 2L)
  saveRDS(list(), file)
  testthat::expect_error(load_cell_graph_asset(index, root, "cell1"), "checksum mismatch")
})

testthat::test_that("coordinates and lossless marker identity lists are validated", {
  root <- withr::local_tempdir()
  asset <- graph_asset_fixture()
  asset$nodes$markers[[4]] <- c("comma,label", "semicolon;label", "<marker>")
  index <- write_graph_fixture(root, list(asset))
  testthat::expect_equal(load_cell_graph_asset(index, root, "cell1")$nodes$markers, asset$nodes$markers)
  asset$nodes$x[1] <- Inf
  testthat::expect_error(validate_cell_graph_asset(asset), "invalid coordinates")
  asset <- graph_asset_fixture()
  asset$nodes$markers[[1]] <- c("CD81", "CD81")
  testthat::expect_error(validate_cell_graph_asset(asset), "marker identities")
  asset$nodes$markers[[1]] <- NA_character_
  testthat::expect_error(validate_cell_graph_asset(asset), "marker identities")
  asset <- graph_asset_fixture()
  asset$nodes$node_id[1] <- asset$nodes$node_id[2]
  testthat::expect_error(validate_cell_graph_asset(asset), "node identities")
})

testthat::test_that("node-only Plotly preserves every node, combinations, identities and escaped hover", {
  asset <- graph_asset_fixture()
  original <- serialize(asset, NULL)
  for (markers in list(character(), "CD81", c("CD81", "CD82"))) {
    plot <- plotly::plotly_build(plot_cell_graph_asset(asset, markers))
    testthat::expect_true(isTRUE(plot$x$config$scrollZoom))
    testthat::expect_equal(sum(vapply(plot$x$data, function(trace) length(trace$x), integer(1))), 4L)
    text <- unlist(lapply(plot$x$data, `[[`, "text"))
    testthat::expect_true(any(grepl("Markers: CD81<br>CD82", text, fixed = TRUE)))
    testthat::expect_false(any(grepl("CD81: 1", text, fixed = TRUE)))
    testthat::expect_true(all(vapply(plot$x$data, function(trace) trace$mode == "markers", logical(1))))
  }
  both <- plotly::plotly_build(plot_cell_graph_asset(asset, c("CD81", "CD82")))
  testthat::expect_setequal(vapply(both$x$data, `[[`, character(1), "name"),
                           c("Other nodes", "CD81", "CD82", "CD81 + CD82"))
  testthat::expect_identical(serialize(asset, NULL), original)
  testthat::expect_error(plot_cell_graph_asset(asset, "unknown"), "up to three")
  testthat::expect_error(plot_cell_graph_asset(asset, rep("CD81", 4)), "up to three")
  asset$nodes$node_id[1] <- "<node>"
  asset$nodes$markers[[1]] <- "<marker>"
  text <- unlist(lapply(plotly::plotly_build(plot_cell_graph_asset(asset, "<marker>"))$x$data, `[[`, "text"))
  testthat::expect_true(any(grepl("&lt;node&gt;", text, fixed = TRUE)))
  testthat::expect_true(any(grepl("&lt;marker&gt;", text, fixed = TRUE)))
})

testthat::test_that("preparation extends samples atomically, resumes and respects frozen source fingerprints", {
  testthat::skip_if_not_installed("pixelatorR")
  root <- withr::local_tempdir()
  pxl_dir <- file.path(root, "source")
  dir.create(pxl_dir)
  saveRDS("synthetic source fingerprint", file.path(pxl_dir, "sample.pxl"))
  asset_dir <- file.path(root, "graphs")
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(matrix(c(1, 2), nrow = 1,
    dimnames = list("CD81", c("cell1", "cell2"))), sparse = TRUE), meta.data = data.frame(
      sample_alias = rep("sample1", 2), row.names = c("cell1", "cell2")))
  object@misc$proxiomevis <- list(proximity_source = list(pxl_manifest = data.frame(
    sample_alias = "sample1", pxl_md5 = unname(tools::md5sum(file.path(pxl_dir, "sample.pxl"))))))
  mapping <- data.frame(sample = 1L, pxl_file = "/original/sample.pxl")
  mapping$id_map <- list(data.frame(current_id = c("cell1", "cell2"), original_id = c("original1", "original2")))
  calls <- character()
  fail <- FALSE
  db <- list(close = function() invisible(NULL), query = function(sql) {
    if (startsWith(sql, "SET ")) return(invisible(NULL))
    if (startsWith(sql, "SELECT DISTINCT")) return(data.frame(layout = "coarsened_pmds_3d", graph_projection = "full"))
    if (fail) stop("synthetic interruption")
    cell <- if (grepl("original1", sql, fixed = TRUE)) "cell1" else "cell2"
    calls <<- c(calls, cell)
    graph_asset_fixture(cell)$nodes
  })
  testthat::local_mocked_bindings(FSMap = function(object) mapping,
    PixelDB = list(new = function(path) db), .package = "pixelatorR")
  first <- suppressMessages(prepare_cell_graph_assets(object, pxl_dir, asset_dir, "cell1", "fixture-v1"))
  first_hash <- unname(tools::md5sum(file.path(asset_dir, "index.rds")))
  testthat::expect_identical(first$coverage, "selected_cells")
  fail <- TRUE
  testthat::expect_error(suppressMessages(prepare_cell_graph_assets(object, pxl_dir, asset_dir, dataset_id = "fixture-v1")), "synthetic interruption")
  testthat::expect_identical(unname(tools::md5sum(file.path(asset_dir, "index.rds"))), first_hash)
  testthat::expect_length(list.files(asset_dir, pattern = "parquet$"), 1L)
  testthat::expect_equal(load_cell_graph_asset(first, asset_dir, "cell1")$nodes, graph_asset_fixture()$nodes)
  fail <- FALSE
  second <- suppressMessages(prepare_cell_graph_assets(object, pxl_dir, asset_dir, dataset_id = "fixture-v1"))
  testthat::expect_identical(second$coverage, "all_cells")
  testthat::expect_identical(calls, c("cell1", "cell2"))
  testthat::expect_false(file.exists(file.path(asset_dir, first$samples$file)))
  testthat::expect_equal(load_cell_graph_asset(second, asset_dir, "cell2")$source_cell_id, "original2")
  testthat::expect_equal(arrow::ParquetFileReader$create(file.path(asset_dir, second$samples$file))$num_row_groups, 2L)
  second_hash <- unname(tools::md5sum(file.path(asset_dir, "index.rds")))
  third <- suppressMessages(prepare_cell_graph_assets(object, pxl_dir, asset_dir, dataset_id = "fixture-v1"))
  testthat::expect_identical(third$cells, second$cells)
  testthat::expect_identical(unname(tools::md5sum(file.path(asset_dir, "index.rds"))), second_hash)
  testthat::expect_identical(calls, c("cell1", "cell2"))
  testthat::expect_error(prepare_cell_graph_assets(object, pxl_dir, asset_dir, dataset_id = "other"), "does not match")
  testthat::expect_error(prepare_cell_graph_assets(object, pxl_dir, asset_dir, "unknown", "fixture-v1"), "known unique cells")
  saveRDS("changed source must not mix with frozen reports", file.path(pxl_dir, "sample.pxl"))
  testthat::expect_error(prepare_cell_graph_assets(object, pxl_dir, asset_dir, dataset_id = "fixture-v1"), "frozen analysis")
  testthat::expect_identical(unname(tools::md5sum(file.path(asset_dir, "index.rds"))), second_hash)
})

testthat::test_that("real DuckDB extraction preserves large IDs, multi-marker memberships and unmarked nodes", {
  testthat::skip_if_not_installed("pixelatorR")
  testthat::skip_if_not_installed("duckdb")
  root <- withr::local_tempdir()
  pxl_dir <- file.path(root, "source")
  dir.create(pxl_dir)
  saveRDS("source fingerprint", file.path(pxl_dir, "sample.pxl"))
  object <- SeuratObject::CreateSeuratObject(Matrix::sparseMatrix(i = 1, j = 1, x = 1,
    dimnames = list("CD81", "cell'1")), meta.data = data.frame(
      sample_alias = "sample1", row.names = "cell'1"))
  mapping <- data.frame(pxl_file = "sample.pxl")
  mapping$id_map <- list(data.frame(current_id = "cell'1", original_id = "original'one"))
  connection <- DBI::dbConnect(duckdb::duckdb())
  withr::defer(DBI::dbDisconnect(connection, shutdown = TRUE))
  DBI::dbExecute(connection, paste0("CREATE TABLE layouts AS SELECT * FROM (VALUES ",
    "('original''one',9007199254740993::BIGINT,'A',1.0,2.0,3.0,'coarsened_pmds_3d','full'),",
    "('original''one',2::BIGINT,'B',4.0,5.0,6.0,'coarsened_pmds_3d','full'),",
    "('original''one',3::BIGINT,'A',7.0,8.0,9.0,'coarsened_pmds_3d','full')) ",
    "AS t(component,index,pixel_type,x,y,z,layout,graph_projection)"))
  DBI::dbExecute(connection, paste0("CREATE TABLE edgelist AS SELECT * FROM (VALUES ",
    "('original''one',9007199254740993::BIGINT,2::BIGINT,'CD,81','CD45'),",
    "('original''one',9007199254740993::BIGINT,2::BIGINT,'CD,81','CD45'),",
    "('original''one',9007199254740993::BIGINT,2::BIGINT,'CD;82','CD45')) ",
    "AS t(component,umi1,umi2,marker_1,marker_2)"))
  mutate_source <- FALSE
  db <- list(query = function(sql) {
    result <- DBI::dbGetQuery(connection, sql)
    if (mutate_source && startsWith(sql, "WITH ")) saveRDS("changed mid-extraction", file.path(pxl_dir, "sample.pxl"))
    result
  }, close = function() invisible(NULL))
  testthat::local_mocked_bindings(FSMap = function(object) mapping,
    PixelDB = list(new = function(path) db), .package = "pixelatorR")
  index <- suppressMessages(prepare_cell_graph_assets(object, pxl_dir, file.path(root, "graphs"), dataset_id = "fixture-v1"))
  asset <- load_cell_graph_asset(index, file.path(root, "graphs"), "cell'1")
  testthat::expect_identical(asset$nodes$node_id, c("2-umi2", "3-umi1", "9007199254740993-umi1"))
  testthat::expect_identical(asset$nodes$markers, list("CD45", character(), c("CD,81", "CD;82")))
  testthat::expect_identical(names(asset$nodes), c("cell_id", "node_id", "x", "y", "z", "markers"))
  testthat::expect_equal(asset$nodes$x, c(4, 7, 1))
  testthat::expect_identical(asset$layout, "coarsened_pmds_3d")
  testthat::expect_identical(asset$graph_projection, "full")
  mutate_source <- TRUE
  testthat::expect_error(suppressMessages(prepare_cell_graph_assets(object, pxl_dir,
    file.path(root, "changed-source-graphs"), dataset_id = "fixture-v1")), "changed during preparation")
  testthat::expect_length(list.files(file.path(root, "changed-source-graphs")), 0L)
})
