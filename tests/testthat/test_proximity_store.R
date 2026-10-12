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

test_that("PXL preparation maps cells and matches finite observed score filtering", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("arrow")
  data <- hybrid_fixture()
  root <- tempfile(); dir.create(root)
  paths <- setNames(file.path(root, c("NC.pxl", "TRT.pxl")), c("D1_NC", "D1_TRT"))
  for (sample in names(paths)) {
    ids <- data$cell_sources[data$cell_sources$sample_id == sample, ]
    raw <- data.frame(component = rep(ids$original_id, 3L),
                       marker_1 = rep(c("A", "A", "B"), each = nrow(ids)),
                       marker_2 = rep(c("A", "B", "B"), each = nrow(ids)),
                       join_count = rep(c(0, 8, 3), each = nrow(ids)), join_count_expected_mean = 2)
    raw$join_count[1L] <- NA_real_
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
  expect_equal(scores$log2_ratio, rep(2, 84L))
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
