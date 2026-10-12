test_that("compact bundles preserve the viewer data without analysis objects", {
  object <- precomputed_fixture()
  object@misc$proxiomevis$report <- list(source_results = list(original = data.frame(p_adj = c(0, NA_real_))))
  data <- precomputed_app_data(object)
  data$cell_sources <- data.frame(component = data$metadata$component,
                                  original_id = paste0("source-", seq_len(nrow(data$metadata))),
                                  sample_id = data$metadata$sample_id)
  data$discarded_analysis_object <- object
  assets <- list(proximity = list(path = "proximity"), graphs = list(path = "graphs"))
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  write_precomputed_bundle(data, file, assets)
  bundle <- readRDS(file)
  expect_identical(bundle$bundle_version, 2L)
  expect_identical(bundle$contract$schema_version, 1L)
  expect_null(bundle$discarded_analysis_object)
  expect_s4_class(bundle$abundance, "sparseMatrix")
  expect_s4_class(bundle$counts, "sparseMatrix")
  restored <- read_precomputed_rds(file)
  expect_identical(restored$bundle_path, normalizePath(file))
  expect_identical(restored$bundle_md5, unname(tools::md5sum(file)))
  expect_identical(restored$assets, assets)
  for (field in c("contract", "metadata", "embeddings", "abundance", "counts", "cell_sources")) {
    expect_identical(restored[[field]], data[[field]])
  }
  write_precomputed_bundle(restored, file)
  expect_false(any(c("bundle_path", "bundle_md5") %in% names(readRDS(file))))
})

test_that("bundle loads reject concurrent replacement instead of mislabeling the snapshot", {
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  write_precomputed_bundle(precomputed_app_data(precomputed_fixture()), file)
  loader <- read_precomputed_rds
  environment(loader) <- list2env(list(readRDS = function(path) {
    object <- base::readRDS(path)
    changed <- object
    changed$contract$title <- "Replaced while reading"
    saveRDS(changed, path)
    object
  }), parent = environment(read_precomputed_rds))
  expect_error(loader(file), "changed while loading")
})

test_that("compact and legacy inputs share strict data validation", {
  valid <- c(list(bundle_version = 2L), precomputed_app_data(precomputed_fixture()))
  invalid <- valid
  invalid$bundle_version <- 3L
  expect_error(precomputed_app_data(invalid), "bundle_version")
  invalid <- valid
  invalid$metadata$celltype[1] <- NA_character_
  expect_error(precomputed_app_data(invalid), "missing values")
  invalid <- valid
  invalid$metadata <- NULL
  expect_error(precomputed_app_data(invalid), "metadata must")
  invalid <- valid
  invalid$embeddings$umap[1, 1] <- Inf
  expect_error(precomputed_app_data(invalid), "invalid coordinates")
  invalid <- valid
  rownames(invalid$embeddings$pca)[1] <- rownames(invalid$embeddings$pca)[2]
  expect_error(precomputed_app_data(invalid), "invalid coordinates")
  invalid <- valid
  invalid$abundance@x[1] <- NA_real_
  expect_error(precomputed_app_data(invalid), "abundance has invalid")
  invalid <- valid
  colnames(invalid$abundance)[1] <- "unknown-cell"
  expect_error(precomputed_app_data(invalid), "abundance has invalid")
  invalid <- valid
  invalid$counts@x[1] <- -1
  expect_error(precomputed_app_data(invalid), "nonnegative integers")
  invalid <- valid
  invalid$contract$sample_summaries$abundance$n_evaluable_cells[1] <- 100L
  expect_error(precomputed_app_data(invalid), "support counts")
  invalid <- valid
  invalid$contract$differential_results$abundance$p_adj[1] <- 0.01
  expect_error(precomputed_app_data(invalid), "descriptive effects")
  invalid <- valid
  invalid$contract$differential_results$colocalization$marker_2[1] <- "unknown-marker"
  expect_error(precomputed_app_data(invalid), "invalid marker pairs")
  invalid <- valid
  invalid$contract$donor_contrasts$self_proximity <- rbind(invalid$contract$donor_contrasts$self_proximity,
                                                        invalid$contract$donor_contrasts$self_proximity[1, ])
  expect_error(precomputed_app_data(invalid), "duplicate keys")
})

test_that("source cell mappings are complete, unique, and sample-aligned", {
  data <- precomputed_app_data(precomputed_fixture())
  data$cell_sources <- data.frame(component = data$metadata$component,
                                  original_id = data$metadata$component,
                                  sample_id = data$metadata$sample_id)
  expect_true(validate_precomputed_data(data))
  invalid <- data
  invalid$cell_sources <- invalid$cell_sources[-1, ]
  expect_error(validate_precomputed_data(invalid), "every metadata cell")
  invalid <- data
  invalid$cell_sources$original_id[2] <- invalid$cell_sources$original_id[1]
  invalid$cell_sources$sample_id[2] <- invalid$cell_sources$sample_id[1]
  expect_error(validate_precomputed_data(invalid), "duplicate keys")
  invalid <- data
  invalid$cell_sources$sample_id[1] <- "wrong-sample"
  expect_error(validate_precomputed_data(invalid), "original sample")
})

test_that("compact loading requires no Seurat or pixelator namespace in a fresh R process", {
  data <- precomputed_app_data(precomputed_fixture())
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  write_precomputed_bundle(data, file)
  source_file <- normalizePath(file.path("..", "..", "R", "precomputed_contract.R"))
  expression <- paste0("source(", deparse(source_file), "); ",
                       "x <- read_precomputed_rds(", deparse(file), "); ",
                       "stopifnot(nrow(x$metadata) == 16L, inherits(x$counts, 'sparseMatrix'), ",
                       "!any(c('SeuratObject', 'Seurat', 'pixelatorR') %in% loadedNamespaces()))")
  output <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", shQuote(expression)),
                    stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
})

test_that("export rejects hidden analysis classes and preserves sparse storage", {
  data <- precomputed_app_data(precomputed_fixture())
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  invalid <- data
  invalid$contract$analysis_object <- precomputed_fixture()
  expect_error(write_precomputed_bundle(invalid, file), "analysis S4")
  expect_false(file.exists(file))
  data$abundance <- as.matrix(data$abundance)
  data$counts <- as.matrix(data$counts)
  write_precomputed_bundle(data, file)
  restored <- read_precomputed_rds(file)
  expect_s4_class(restored$abundance, "sparseMatrix")
  expect_equal(as.matrix(restored$abundance), data$abundance)
  expect_equal(as.matrix(restored$counts), data$counts)
})
