ranked_graph_fixture <- function() {
  root <- tempfile("ranked-graphs-")
  dir.create(root)
  data <- precomputed_app_data(precomputed_fixture())
  # Five cells in every donor/type/condition stratum; a separate sample has no scores.
  meta <- expand.grid(donor = c("D1", "D2"), celltype = c("B", "T"),
    condition = c("NC", "TRT"), replicate = 1:5, stringsAsFactors = FALSE)
  meta$component <- sprintf("cell%03d", seq_len(nrow(meta)))
  meta$sample_id <- paste(meta$donor, meta$condition, sep = "_")
  empty <- meta[1, ]; empty$component <- "empty"; empty$sample_id <- "empty_sample"
  data$metadata <- rbind(meta, empty)
  data$contract$proximity_source <- NULL
  scores <- do.call(rbind, lapply(1:3, function(pair) {
    data.frame(component = meta$component, marker_1 = c("A", "B", "C")[[pair]],
      marker_2 = c("A", "A", "A")[[pair]],
      log2_ratio = if (pair == 1L) c(7, 6, 9, 8, 5)[meta$replicate] else
        if (pair == 2L) c(1, 5, 3, 4, 2)[meta$replicate] else c(3, 4, 5, 1, 2)[meta$replicate])
  }))
  # Include both orientations in different cells, not duplicate observations.
  flip <- scores$marker_1 == "B" & match(scores$component, meta$component) %% 2L == 0L
  scores$marker_1[flip] <- "A"; scores$marker_2[flip] <- "B"
  samples <- unique(data$metadata$sample_id)
  files <- paste0(samples, ".parquet")
  for (i in seq_along(samples)) {
    ids <- data$metadata$component[data$metadata$sample_id == samples[[i]]]
    arrow::write_parquet(scores[scores$component %in% ids, ], file.path(root, files[[i]]))
  }
  manifest <- list(schema_version = 1L, dataset_id = data$contract$dataset_id,
    dataset_fingerprint = hybrid_dataset_fingerprint(data),
    partitions = data.frame(sample_id = samples,
      file = files, md5 = unname(tools::md5sum(file.path(root, files)))))
  saveRDS(manifest, file.path(root, "manifest.rds"))
  list(data = data, root = root, scores = scores, manifest = manifest)
}

test_that("marker selection unions top/bottom two different-marker cells in every donor/type/condition", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("arrow")
  fixture <- ranked_graph_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  result <- select_marker_graph_cells(fixture$data, fixture$root, c("A", "B"))
  meta <- fixture$data$metadata
  expected <- meta$component[meta$replicate %in% c(1, 2, 4, 5) & meta$component != "empty"]
  expect_setequal(result$cells, expected)
  expect_identical(anyDuplicated(result$cells), 0L)
  expect_equal(nrow(result$selected), 32L)
  expect_equal(result$summary$selected_cells, c(8, 8, 8, 8, 0))
  expect_equal(result$summary$colocalization_top_cells, c(4, 4, 4, 4, 0))
  expect_equal(result$summary$colocalization_bottom_cells, c(4, 4, 4, 4, 0))
  expect_setequal(select_marker_graph_cells(fixture$data, fixture$root, c("A", "B", "C"))$cells,
    setdiff(meta$component, "empty"))
  expect_identical(select_marker_graph_cells(fixture$data, fixture$root, c("B", "A")), result)
  expect_identical(result$parameters$markers, c("A", "B"))
  expect_identical(result$parameters$strata, c("donor", "celltype", "condition"))
})

test_that("graph ranking rejects reversed duplicates, ignores nonfinite scores, and breaks ties by ID", {
  skip_if_not_installed("duckdb"); skip_if_not_installed("arrow")
  fixture <- ranked_graph_fixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  update_scores <- function(scores) {
    for (i in seq_len(nrow(fixture$manifest$partitions))) {
      partition <- fixture$manifest$partitions[i, ]
      ids <- fixture$data$metadata$component[fixture$data$metadata$sample_id == partition$sample_id]
      path <- file.path(fixture$root, partition$file)
      arrow::write_parquet(scores[scores$component %in% ids, ], path)
      fixture$manifest$partitions$md5[[i]] <- unname(tools::md5sum(path))
    }
    saveRDS(fixture$manifest, file.path(fixture$root, "manifest.rds"))
  }
  scores <- fixture$scores
  duplicate <- scores[scores$marker_1 == "B", ][1, ]
  duplicate$marker_1 <- "A"; duplicate$marker_2 <- "B"
  update_scores(rbind(scores, duplicate))
  expect_error(select_marker_graph_cells(fixture$data, fixture$root, c("A", "B")), "Duplicate proximity")
  scores$log2_ratio <- 1
  scores$log2_ratio[scores$marker_1 == "C"] <- Inf
  update_scores(scores)
  result <- select_marker_graph_cells(fixture$data, fixture$root, c("A", "B", "C"))
  meta <- fixture$data$metadata
  expect_setequal(result$cells, meta$component[meta$replicate %in% 1:2 & meta$component != "empty"])
  expect_equal(sum(result$summary$selected_cells), 16)
  scores$log2_ratio <- NA_real_
  update_scores(scores)
  expect_error(select_marker_graph_cells(fixture$data, fixture$root, c("A", "B")), "No finite")
  scores$log2_ratio[scores$marker_1 == scores$marker_2] <- 100
  update_scores(scores)
  expect_error(select_marker_graph_cells(fixture$data, fixture$root, c("A", "B")), "No finite")
  scores$log2_ratio[scores$component == "cell001" & scores$marker_1 == "B"] <- 1
  update_scores(scores)
  expect_identical(select_marker_graph_cells(fixture$data, fixture$root, c("A", "B"))$cells, "cell001")
  for (markers in list(character(), "A", c("A", "A"), c("A", NA_character_), c("A", ""), 42)) {
    expect_error(select_marker_graph_cells(fixture$data, fixture$root, markers), "at least two distinct")
  }
  expect_error(select_marker_graph_cells(fixture$data, fixture$root, c("A", "unknown")), "Unknown graph markers")
})
