test_that("ranking uses the exact cohort, supports self pairs and either pair order without test cutoffs", {
  fixture <- ranking_fixture()
  rank <- function(a, b, donor = "D1", celltype = "B", condition = "NC") {
    rank_cell_proximity(fixture$data, fixture$root, donor, celltype, condition, a, b)
  }
  ab <- rank("A", "B")
  expect_identical(ab, rank("B", "A"))
  expect_identical(ab$component, c("cell9", "cell1"))
  expect_equal(ab$log2_ratio, c(3, 1))
  expect_identical(ab$rank, 1:2)
  expect_true(all(ab$marker_1 == "A" & ab$marker_2 == "B"))
  expect_identical(rank("A", "A")$component, ab$component)
  expect_identical(rank("A", "A", "D2", "T", "TRT")$component, c("cell16", "cell8"))
  expect_equal(attr(ab, "provenance")$cohort_cells, 2L)
  expect_match(attr(ab, "provenance")$filtering, "no imputation")
  expect_equal(nrow(rank("C", "C")), 0L)
  expect_error(rank("unknown", "A"), "known markers")
  expect_error(rank(c("A", "B"), "A"), "one marker")
  expect_error(rank("A", "A", condition = "unknown"), "No cells match")
})

test_that("missing scores stay missing, ties are deterministic, and duplicate reversed rows fail closed", {
  fixture <- ranking_fixture()
  path <- file.path(fixture$root, fixture$manifest$partitions$file[1L])
  original <- arrow::read_parquet(path)
  update <- function(rows) {
    arrow::write_parquet(rows, path)
    fixture$manifest$partitions$md5[1L] <- unname(tools::md5sum(path))
    saveRDS(fixture$manifest, file.path(fixture$root, "manifest.rds"))
  }
  rank <- function() rank_cell_proximity(fixture$data, fixture$root, "D1", "B", "NC", "A", "B")
  tied <- original
  tied$log2_ratio <- -2
  update(tied)
  expect_identical(rank()$component, c("cell1", "cell9"))
  missing <- tied
  missing$log2_ratio[missing$component == "cell1"] <- NA_real_
  update(missing)
  expect_identical(rank()$component, "cell9")
  missing$log2_ratio <- Inf
  update(missing)
  expect_equal(nrow(rank()), 0L)
  reversed <- original[original$component == "cell1" & original$marker_1 == "B", ]
  reversed$marker_1 <- "A"; reversed$marker_2 <- "B"
  update(rbind(original, reversed))
  expect_error(rank(), "Duplicate proximity observations")
  update(original)
  arrow::write_parquet(tied, path)
  expect_error(rank(), "partition changed")
  fixture$data$counts[1, 1] <- 999
  expect_error(rank(), "does not match")
})
