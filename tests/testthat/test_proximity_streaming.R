test_that("bounded proximity scans preserve filtered batches, empty schemas, and row limits", {
  skip_if_not_installed("arrow")
  data <- hybrid_fixture()
  root <- withr::local_tempdir()
  raw <- data.frame(component = rep(data$metadata$component, 3L),
    marker_1 = rep(c("A", "A", "B"), each = 84L),
    marker_2 = rep(c("A", "B", "B"), each = 84L),
    log2_ratio = rep(seq_len(84L), 3L), min_marker_fraction = 0.1)
  raw$log2_ratio[85L] <- Inf
  raw$min_marker_fraction[86L] <- 0
  path <- file.path(root, "scores.parquet")
  arrow::write_parquet(raw, path, chunk_size = 8L)
  manifest_path <- file.path(root, "manifest.rds")
  manifest <- list(donor_cutoffs = c(D1 = 0.001))
  saveRDS(manifest, manifest_path)
  reader <- read_proximity_subset
  environment(reader) <- list2env(list(validated_proximity_partitions = function(...) {
    list(paths = path, manifest_path = manifest_path, manifest = manifest)
  }), parent = environment(read_proximity_subset))
  request <- hybrid_request("colocalization")

  scores <- reader(data, root, request)
  expected <- raw[raw$marker_1 != raw$marker_2 & is.finite(raw$log2_ratio) &
    raw$min_marker_fraction >= 0.001, , drop = FALSE]
  expect_equal(nrow(scores), 82L)
  expect_equal(lapply(dplyr::arrange(scores, component), identity),
    lapply(dplyr::arrange(expected, component), identity))
  expect_equal(attr(scores, "provenance")$observed_rows, 82L)

  empty <- reader(data, root, request, marker_1 = "B")
  expect_identical(names(empty), names(raw))
  expect_equal(nrow(empty), 0L)
  expect_equal(attr(empty, "provenance")$observed_rows, 0L)
  expect_equal(nrow(reader(data, root, hybrid_request("self_proximity"))), 168L)

  environment(reader)$HYBRID_LIMITS <- modifyList(HYBRID_LIMITS, list(rows = 81L))
  expect_error(reader(data, root, request), "Too many proximity rows")
  environment(reader)$HYBRID_LIMITS$rows <- 82L
  expect_equal(nrow(reader(data, root, request)), 82L)
  expect_equal(nrow(reader(data, root, request, marker_1 = "B")), 0L)
})
