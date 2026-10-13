test_that("abundance distributions use exact cell values and comparison scope", {
  data <- precomputed_app_data(precomputed_fixture())
  data$abundance["A", "cell9"] <- data$abundance["A", "cell9"] + 0.2
  rows <- precomputed_feature_distribution_data(data, "abundance", "D1", "B", c("NC", "TRT"), c("A", "A"))
  expected <- with(data$metadata, component[donor == "D1" & celltype == "B"])
  expect_identical(rows$component, expected)
  expect_equal(rows$value, as.numeric(data$abundance["A", expected]))
  expect_identical(levels(rows$condition), c("NC", "TRT"))
  expect_true(all(precomputed_feature_distribution_data(data, "abundance", "D1", PRECOMPUTED_ALL,
    "TRT", c("A", "A"))$condition == "TRT"))
  expect_error(precomputed_feature_distribution_data(data, "abundance", "D1", "B", "NC", c("A", "B")), "matches")
  plot <- plot_precomputed_feature_distribution(rows, "abundance")
  expect_true(any(vapply(plot$layers, function(x) inherits(x$geom, "GeomViolin"), logical(1))))
  expect_true(any(vapply(plot$layers, function(x) inherits(x$geom, "GeomBoxplot"), logical(1))))
  expect_s3_class(plotly::plotly_build(plotly::ggplotly(plot)), "plotly")
  file <- tempfile(fileext = ".svg")
  on.exit(unlink(file), add = TRUE)
  precomputed_write_svg(plot, file, 10, 7)
  expect_gt(file.info(file)$size, 0)
  for (small in list(rows[1, ], transform(rows, value = 1), rows[FALSE, ])) {
    expect_s3_class(plotly::plotly_build(plotly::ggplotly(
      plot_precomputed_feature_distribution(small, "abundance"))), "plotly")
  }
})

test_that("spatial distributions use observed unordered pair scores without imputation", {
  root <- withr::local_tempdir()
  fixture <- ranking_fixture(root = file.path(root, "proximity"))
  data <- fixture$data
  data$bundle_path <- file.path(root, "analysis.rds")
  data$assets <- list(proximity = list(path = "proximity"))
  read <- function(pair = c("A", "B"), readout = "colocalization")
    precomputed_feature_distribution_data(data, readout, "D1", "B", c("NC", "TRT"), pair)
  rows <- read()
  expect_equal(rows, read(c("B", "A")))
  expect_equal(nrow(rows), 4L)
  expect_equal(rows$value, c(1, 1, 3, 3))
  expect_equal(read(c("A", "A"), "self_proximity")$value, rows$value)
  expect_equal(nrow(read(c("A", "C"))), 0L)
  partition <- file.path(fixture$root, fixture$manifest$partitions$file[1L])
  scores <- arrow::read_parquet(partition)
  scores$log2_ratio[scores$component == "cell1"] <- NA_real_
  arrow::write_parquet(scores, partition)
  fixture$manifest$partitions$md5[1L] <- unname(tools::md5sum(partition))
  saveRDS(fixture$manifest, file.path(fixture$root, "manifest.rds"))
  expect_equal(nrow(read()), 3L)
  expect_false("cell1" %in% read()$component)
  data$assets <- NULL
  expect_error(read(), "require the prepared proximity sidecar")
})

test_that("Tables shows the distribution, follows row selection, and exports the same plot", {
  current <- shiny::reactiveVal(precomputed_app_data(precomputed_fixture()))
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(readout = "abundance", donor = "D1", celltype = "B", comparison = "TRT vs NC",
      report_page = "Precomputed results", results_view = "Tables")
    select_feature(report_feature_key("A", "A"))
    expect_match(detail_ggplot()$labels$title, "cell-level distribution")
    expect_equal(nrow(detail_ggplot()$data), 4L)
    rows <- display_results()
    session$setInputs(results_table_rows_selected = which(rows$marker_1 == "B"))
    expect_equal(detail_ggplot()$data$value, as.numeric(current()$abundance["B", detail_ggplot()$data$component]))
    expect_gt(file.info(output$download_detail_svg)$size, 0)
    session$setInputs(celltype = "T")
    expect_true(all(detail_ggplot()$data$celltype == "T"))
    session$setInputs(readout = "self_proximity")
    expect_match(detail_ggplot()$labels$title, "unavailable")
    session$setInputs(results_view = "Summary")
    expect_match(detail_ggplot()$labels$title, "stored effects")
  })
})
