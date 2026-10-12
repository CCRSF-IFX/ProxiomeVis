test_that("versioned input resolves defaults and works with unavailable PXLs", {
  object <- precomputed_fixture()
  expect_true(validate_precomputed_object(object))
  data <- precomputed_app_data(object)
  expect_identical(data$metadata$component, colnames(data$abundance))
  expect_identical(rownames(data$embeddings$umap), data$metadata$component)
  expect_s4_class(data$abundance, "sparseMatrix")
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(object, file)
  expect_identical(read_precomputed_rds(file)$contract$dataset_id, "fixture-v1")
})

test_that("invalid schema, annotation defaults and duplicated results are rejected", {
  object <- precomputed_fixture()
  object@misc$proxiomevis$schema_version <- 2L
  expect_error(validate_precomputed_object(object), "schema_version")
  object@misc$proxiomevis$schema_version <- 1.5
  expect_error(validate_precomputed_object(object), "schema_version")
  object@misc$proxiomevis$schema_version <- 1L
  object@misc$proxiomevis$defaults$celltype_column <- "absent"
  expect_error(validate_precomputed_object(object), "metadata column")
  object <- precomputed_fixture()
  object@misc$proxiomevis$differential_results$abundance <- rbind(
    object@misc$proxiomevis$differential_results$abundance,
    object@misc$proxiomevis$differential_results$abundance[1, ]
  )
  expect_error(validate_precomputed_object(object), "duplicate keys")
})

test_that("support, direction, pair types and descriptive significance are validated", {
  object <- precomputed_fixture()
  object@misc$proxiomevis$sample_summaries$self_proximity$n_scored_cells[1] <- 3L
  expect_error(validate_precomputed_object(object), "support counts")
  object <- precomputed_fixture()
  object@misc$proxiomevis$donor_contrasts$abundance$effect_size[1] <- -1
  expect_error(validate_precomputed_object(object), "contrast direction")
  object <- precomputed_fixture()
  object@misc$proxiomevis$differential_results$colocalization$marker_2[1] <- "A"
  expect_error(validate_precomputed_object(object), "wrong pair type")
  object <- precomputed_fixture()
  object@misc$proxiomevis$differential_results$abundance$p_value[1] <- 0.001
  expect_error(validate_precomputed_object(object), "descriptive effects")
  object <- precomputed_fixture()
  object@misc$proxiomevis$donor_contrasts$abundance$p_adj <- 2
  expect_error(validate_precomputed_object(object), "invalid probabilities")
})

test_that("missing scores stay missing and cell order is aligned explicitly", {
  object <- precomputed_fixture()
  object@misc$proxiomevis$sample_summaries$colocalization$value[1] <- NA_real_
  object@misc$proxiomevis$sample_summaries$colocalization$n_evaluable_cells <- NA_integer_
  object@reductions$umap@cell.embeddings <- object@reductions$umap@cell.embeddings[16:1, ]
  data <- precomputed_app_data(object)
  expect_true(is.na(data$contract$sample_summaries$colocalization$value[1]))
  expect_equal(data$embeddings$umap[1, 1], 1)
})
