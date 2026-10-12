test_that("requests bind every cell to one donor, cell type and comparison", {
  data <- hybrid_fixture()
  request <- validate_hybrid_request(data, hybrid_request())
  expect_length(request$group_a_cells, 42L)
  expect_length(request$group_b_cells, 42L)
  expect_identical(request$dataset_id, data$contract$dataset_id)
  bad <- hybrid_request(); bad$group_a_cells <- request$group_b_cells
  expect_error(validate_hybrid_request(data, bad), "belong")
  bad <- hybrid_request(); bad$dataset_id <- "other"
  expect_error(validate_hybrid_request(data, bad), "dataset changed")
  data$bundle_md5 <- "current-snapshot"
  bad <- hybrid_request(); bad$bundle_md5 <- "old-snapshot"
  expect_error(validate_hybrid_request(data, bad), "bundle changed")
  expect_identical(validate_hybrid_request(data, hybrid_request())$bundle_md5, "current-snapshot")
  bad <- hybrid_request(); bad$reference <- "TRT"
  expect_error(validate_hybrid_request(data, bad), "different")
  bad <- hybrid_request(); bad$group_a_cells <- request$group_a_cells[1:2]
  expect_error(validate_hybrid_request(data, bad), "at least 3")
  bad <- hybrid_request(); bad$group_a_cells <- rep(request$group_a_cells[1], 3L)
  expect_error(validate_hybrid_request(data, bad), "unique")
  bad <- hybrid_request(); bad$donor <- "D2"
  expect_error(validate_hybrid_request(data, bad), "No cells")
  bad <- hybrid_request(); bad$markers <- "absent"
  expect_error(validate_hybrid_request(data, bad), "known markers")
  expect_error(validate_hybrid_request(data, hybrid_request(), modifyList(HYBRID_LIMITS, list(cells = 80L))), "at most 80")
  expect_error(validate_hybrid_request(data, hybrid_request("colocalization"), modifyList(HYBRID_LIMITS, list(pairs = 0L))), "marker pairs")
  selected <- hybrid_request(); selected$group_a_cells <- request$group_a_cells[1:3]
  selected$group_b_cells <- request$group_b_cells[1:3]
  expect_length(validate_hybrid_request(data, selected)$group_a_cells, 3L)
  selected$readout <- "self_proximity"
  expect_error(validate_hybrid_request(data, selected), "at least 20")
})

test_that("replacing a bundle with the same dataset ID rejects the pinned request", {
  path <- tempfile(fileext = ".rds")
  object <- precomputed_fixture()
  saveRDS(object, path)
  loaded <- read_precomputed_rds(path)
  request <- hybrid_request()
  request$dataset_id <- loaded$contract$dataset_id
  request$bundle_md5 <- loaded$bundle_md5
  values <- SeuratObject::LayerData(object, assay = "PNA", layer = "data")
  values[1L, 1L] <- values[1L, 1L] + 1
  SeuratObject::LayerData(object, assay = "PNA", layer = "data") <- values
  saveRDS(object, path)
  expect_error(run_hybrid_analysis(path, request), "bundle changed")
})

test_that("selected-cell abundance exactly reuses the report method", {
  skip_if_not_installed("pixelatorR")
  data <- hybrid_fixture()
  path <- tempfile(fileext = ".rds")
  saveRDS(data, path)
  runner <- run_hybrid_analysis
  environment(runner) <- list2env(list(read_precomputed_rds = function(path) readRDS(path)), parent = environment(run_hybrid_analysis))
  result <- runner(path, hybrid_request())
  request <- validate_hybrid_request(data, hybrid_request())
  cells <- c(request$group_a_cells, request$group_b_cells)
  object <- SeuratObject::CreateSeuratObject(data$counts[request$markers, cells], assay = "PNA",
    meta.data = data.frame(condition = rep(c("TRT", "NC"), each = 42L), row.names = cells))
  SeuratObject::LayerData(object, assay = "PNA", layer = "data") <- data$abundance[request$markers, cells]
  expected <- pixelatorR::RunDAA(object, contrast_column = "condition", targets = "TRT", reference = "NC",
    assay = "PNA", mean_fxn = Matrix::rowMeans, fc_name = "difference", p_adjust_method = "BH",
    slot = "data", fc.slot = "data", min.pct = 0, verbose = FALSE)
  expected <- expected[match(result$results$marker_1, expected$marker), ]
  expect_equal(result$results$effect_size, expected$difference)
  expect_equal(result$results$p_value, expected$p)
  expect_equal(result$results$p_adj, expected$p_adj)
  expect_true(all(result$results$effect_size < 0))
  expect_true(all(result$results$n_cells_a == 42L))
  expect_identical(result$provenance$request$group_a_cells, request$group_a_cells)
  expect_match(result$provenance$inference, "not biological replicates")
})

test_that("spatial tests retain missingness and report support separately", {
  skip_if_not_installed("pixelatorR")
  data <- hybrid_fixture()
  request <- validate_hybrid_request(data, hybrid_request("self_proximity"))
  cells <- c(request$group_a_cells, request$group_b_cells)
  scores <- data.frame(component = rep(cells, 2L), marker_1 = rep(c("A", "B"), each = length(cells)),
                       marker_2 = rep(c("A", "B"), each = length(cells)),
                       log2_ratio = rep(seq_along(cells) / 10, 2L))
  scores <- scores[!(scores$marker_1 == "B" & scores$component %in% request$group_a_cells[1:30]), ]
  path <- tempfile(fileext = ".rds"); saveRDS(data, path)
  runner <- run_hybrid_analysis
  environment(runner) <- list2env(list(read_precomputed_rds = function(path) readRDS(path),
    read_proximity_subset = function(data, asset_root, request) scores), parent = environment(run_hybrid_analysis))
  result <- runner(path, request, asset_root = tempdir())
  expect_identical(result$results$n_scored_a, c(42L, 12L))
  expect_identical(result$results$n_cells_a, c(42L, 42L))
  expect_identical(result$results$status, c("tested", "insufficient_observed_scores"))
  expect_true(is.na(result$results$p_value[2L]))
  scores$condition <- ifelse(scores$component %in% request$group_a_cells, "TRT", "NC")
  expected <- pixelatorR::DifferentialProximityAnalysis(scores, contrast_column = "condition", targets = "TRT",
    reference = "NC", proximity_metric = "log2_ratio", metric_type = "self", p_adjust_method = "BH",
    min_cells_per_group = 20, verbose = FALSE)
  expect_equal(result$results$p_value[1L], expected$p)
  expect_equal(result$results$effect_size[1L], expected$diff_median)
  expect_equal(result$results$p_adj[1L], expected$p_adj)
  expect_equal(result$results$p_value[1L], stats::wilcox.test(scores$log2_ratio[scores$marker_1 == "A" & scores$condition == "TRT"],
    scores$log2_ratio[scores$marker_1 == "A" & scores$condition == "NC"], exact = FALSE)$p.value)
})

test_that("selected colocalization uses observed pair scores and selected-family BH", {
  skip_if_not_installed("pixelatorR")
  data <- hybrid_fixture()
  input <- hybrid_request("colocalization"); input$markers <- c("A", "B", "C")
  request <- validate_hybrid_request(data, input)
  cells <- c(request$group_a_cells, request$group_b_cells)
  scores <- data.frame(component = rep(cells, 2L), marker_1 = "A",
    marker_2 = rep(c("B", "C"), each = length(cells)),
    log2_ratio = c(seq_along(cells) / 10, rep(c(0, 0, 1, 2), 21L)))
  path <- tempfile(fileext = ".rds"); saveRDS(data, path)
  runner <- run_hybrid_analysis
  environment(runner) <- list2env(list(read_precomputed_rds = function(path) readRDS(path),
    read_proximity_subset = function(data, asset_root, request) scores), parent = environment(run_hybrid_analysis))
  result <- runner(path, request, asset_root = tempdir())
  expected_p <- vapply(c("B", "C"), function(marker) {
    a <- scores$log2_ratio[scores$marker_2 == marker & scores$component %in% request$group_a_cells]
    b <- scores$log2_ratio[scores$marker_2 == marker & scores$component %in% request$group_b_cells]
    stats::wilcox.test(a, b, exact = FALSE)$p.value
  }, numeric(1))
  expect_equal(result$results$p_value[1:2], unname(expected_p), tolerance = 1e-12)
  expect_equal(result$results$p_adj[1:2], unname(p.adjust(expected_p, "BH")), tolerance = 1e-12)
  expect_identical(result$results$status, c("tested", "tested", "insufficient_observed_scores"))
  expect_identical(result$provenance$tested_features, 2L)
  expect_equal(result$results$n_scored_a, c(42L, 42L, 0L))
})
