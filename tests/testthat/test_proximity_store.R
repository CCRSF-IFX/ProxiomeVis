test_that("asset resolvers reject traversal and escaped partitions", {
  root <- tempfile(); dir.create(root)
  dir.create(file.path(root, "proximity"))
  path <- file.path(root, "analysis.rds")
  expect_equal(resolve_hybrid_asset(path, list(path = "proximity")), normalizePath(file.path(root, "proximity")))
  expect_error(resolve_hybrid_asset(path, list(path = "../outside")), "relative path")
  expect_error(resolve_hybrid_asset(path, list(path = "/tmp")), "relative path")
  expect_error(proximity_store_file(root, "../file.parquet"), "Invalid")
  outside <- tempfile(); dir.create(outside)
  file.symlink(outside, file.path(root, "escape"))
  expect_error(resolve_hybrid_asset(path, list(path = "escape")), "inside")
})

test_that("sidecar background filtering matches FilterProximityScores at the cutoff", {
  skip_if_not_installed("pixelatorR")
  data <- hybrid_fixture()
  root <- withr::local_tempdir()
  connection <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1"))
  withr::defer(DBI::dbDisconnect(connection, shutdown = TRUE))
  raw <- data.frame(component = data$metadata$component, marker_1 = "A", marker_2 = "B",
    log2_ratio = seq_len(84), p1 = rep(c(0.0005, 0.001, 0.1), 28),
    p2 = rep(c(0.1, 0.01, 0.0005), 28))
  DBI::dbWriteTable(connection, "scores", raw)
  expected <- dplyr::collect(pixelatorR::FilterProximityScores(dplyr::tbl(connection, "scores"),
    background_threshold_pct = 0.001))
  stored <- raw[c("component", "marker_1", "marker_2", "log2_ratio")]
  stored$min_marker_fraction <- pmin(raw$p1, raw$p2)
  path <- file.path(root, "scores.parquet")
  arrow::write_parquet(stored, path)
  manifest_path <- file.path(root, "manifest.rds")
  manifest <- list(donor_cutoffs = c(D1 = 0.001))
  saveRDS(manifest, manifest_path)
  reader <- read_proximity_subset
  environment(reader) <- list2env(list(validated_proximity_partitions = function(...) {
    list(paths = path, manifest_path = manifest_path, manifest = manifest)
  }), parent = environment(read_proximity_subset))
  actual <- reader(data, root, hybrid_request("colocalization"))
  expect_setequal(actual$component, expected$component)
  expect_equal(nrow(actual), 28L)
  expect_true(attr(actual, "provenance")$background_filter_enabled)
  expect_equal(attr(actual, "provenance")$background_threshold_pct, 0.001)
})

test_that("PXL preparation maps cells and matches finite observed score filtering", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("arrow")
  data <- hybrid_fixture()
  # Overlapping ID sets must not trigger a second mapping after preparation.
  data$cell_sources$original_id <- rev(data$cell_sources$component)
  expected_scores <- stats::setNames(log2(seq_len(nrow(data$metadata)) + 1), data$metadata$component)
  root <- tempfile(); dir.create(root)
  paths <- setNames(file.path(root, c("NC.pxl", "TRT.pxl")), c("D1_NC", "D1_TRT"))
  for (sample in names(paths)) {
    ids <- data$cell_sources[data$cell_sources$sample_id == sample, ]
    raw <- data.frame(component = rep(ids$original_id, 3L),
                       marker_1 = rep(c("A", "A", "B"), each = nrow(ids)),
                       marker_2 = rep(c("A", "B", "B"), each = nrow(ids)),
                       join_count = rep(c(0, 8, 3), each = nrow(ids)), join_count_expected_mean = 2)
    raw$join_count[1L] <- NA_real_
    raw$join_count[raw$marker_1 != raw$marker_2] <- 2 * (match(ids$component, data$metadata$component) + 1)
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = paths[[sample]], config = list(threads = "1"))
    DBI::dbWriteTable(con, "proximity", raw)
    DBI::dbDisconnect(con, shutdown = TRUE)
  }
  out <- file.path(root, "proximity")
  manifest <- prepare_proximity_store(data, paths, out)
  expect_equal(sum(manifest$partitions$rows), 84 * 3)
  expect_identical(manifest$partitions$source_pxl_md5, unname(tools::md5sum(paths[manifest$partitions$sample_id])))
  expect_equal(unname(manifest$donor_cutoffs), .001)
  request <- validate_hybrid_request(data, hybrid_request("colocalization"))
  scores <- read_proximity_subset(data, out, request)
  expect_equal(nrow(scores), 84L)
  expect_true(all(scores$component %in% data$metadata$component))
  expect_true(all(scores$marker_1 == "A" & scores$marker_2 == "B"))
  expect_equal(scores$log2_ratio, unname(expected_scores[scores$component]))
  narrowed <- request
  narrowed$group_a_cells <- head(request$group_a_cells, 40L)
  narrowed$group_b_cells <- tail(request$group_b_cells, 40L)
  scoped <- read_proximity_subset(data, out, narrowed)
  expect_setequal(scoped$component, c(narrowed$group_a_cells, narrowed$group_b_cells))
  expect_equal(scoped$log2_ratio, unname(expected_scores[scoped$component]))
  expect_equal(read_proximity_subset(data, out, request, marker_1 = "A"), scores)
  expect_equal(nrow(read_proximity_subset(data, out, request, marker_1 = "B")), 0L)
  pinned <- data
  pinned$contract$proximity_source$pxl_manifest <- data.frame(
    sample_alias = manifest$partitions$sample_id, pxl_md5 = manifest$partitions$source_pxl_md5)
  expect_setequal(read_proximity_subset(pinned, out, request)$component, scores$component)
  changed_source <- pinned
  changed_source$contract$proximity_source$pxl_manifest$pxl_md5[[1L]] <- "new-scoring-source"
  expect_identical(hybrid_dataset_fingerprint(changed_source), hybrid_dataset_fingerprint(pinned))
  expect_error(read_proximity_subset(changed_source, out, request), "source PXL checksums do not match")
  for (failure in c("missing", "duplicate", "empty_hash", "missing_column")) {
    invalid <- pinned
    inventory <- invalid$contract$proximity_source$pxl_manifest
    if (failure == "missing") inventory <- inventory[-1L, , drop = FALSE]
    if (failure == "duplicate") inventory <- rbind(inventory, inventory[1L, ])
    if (failure == "empty_hash") inventory$pxl_md5[[1L]] <- NA_character_
    if (failure == "missing_column") inventory$pxl_md5 <- NULL
    invalid$contract$proximity_source$pxl_manifest <- inventory
    expect_error(read_proximity_subset(invalid, out, request), "source PXL")
  }
  incomplete <- manifest
  incomplete$partitions$source_pxl_md5 <- NULL
  saveRDS(incomplete, file.path(out, "manifest.rds"))
  expect_error(read_proximity_subset(pinned, out, request), "Missing source PXL checksums")
  saveRDS(manifest, file.path(out, "manifest.rds"))
  reader <- read_proximity_subset
  environment(reader) <- list2env(list(HYBRID_LIMITS = modifyList(HYBRID_LIMITS, list(rows = 80L))),
                                   parent = environment(read_proximity_subset))
  expect_error(reader(data, out, request), "Too many proximity rows")
  self <- read_proximity_subset(data, out, validate_hybrid_request(data, hybrid_request("self_proximity")))
  expect_equal(nrow(self), 166L)
  expect_true(all(is.finite(self$log2_ratio)))
  expect_equal(unique(self$log2_ratio[self$marker_1 == "A"]), -1)
  expect_error(prepare_proximity_store(data, paths, out), "not overwritten")
  mismatched <- data
  mismatched$contract$proximity_source$pxl_manifest <- data.frame(sample_alias = names(paths), pxl_md5 = "wrong-source")
  expect_error(prepare_proximity_store(mismatched, paths, file.path(root, "wrong-generation")), "frozen report provenance")
  altered <- data; altered$counts[1, 1] <- 999
  expect_error(read_proximity_subset(altered, out, request), "does not match")
  altered <- data; altered$metadata$isotype_fraction[1L] <- .02
  expect_error(read_proximity_subset(altered, out, request), "does not match")
  partition <- file.path(out, manifest$partitions$file[[1L]])
  changed <- arrow::read_parquet(partition); changed$log2_ratio[1L] <- 0
  arrow::write_parquet(changed, partition)
  expect_error(read_proximity_subset(data, out, request), "partition changed")
})
