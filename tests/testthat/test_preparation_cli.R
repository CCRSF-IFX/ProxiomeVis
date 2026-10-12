source(testthat::test_path("..", "..", "R", "preparation_cli.R"))

test_that("preparation CLI preserves positional arguments and rejects malformed options", {
  options <- list(umap = NULL, "graph-markers" = NULL)
  parsed <- parse_preparation_cli(c("in", "--umap", "custom", "out", "--graph-markers", "A,B"),
    options, 2L, "usage")
  expect_identical(parsed$positional, c("in", "out"))
  expect_identical(parsed$options, list(umap = "custom", "graph-markers" = "A,B"))
  expect_error(parse_preparation_cli(c("in", "out", "--typo", "x"), options, 2L, "usage"), "Unknown option")
  expect_error(parse_preparation_cli(c("--umap", "a", "--umap", "b"), options, 2L, "usage"), "Repeated")
  expect_error(parse_preparation_cli(c("in", "out", "--umap"), options, 2L, "usage"), "Supply a value")
  expect_error(parse_preparation_cli(c("--umap", "--graph-markers", "A,B"), options, 2L, "usage"), "Supply a value")
  expect_error(parse_preparation_cli("in", options, 2L, "usage"), "usage")
})

test_that("all nine command-line settings reach the prepared RDS and application reader", {
  skip_if_not_installed("pixelatorR")
  repository <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  demo <- file.path(repository, "..", "pixelgen_demo_data", "demo_analysis")
  input <- file.path(demo, "cca_analysis", "results", "demo_cca_clustered.rds")
  skip_if_not(file.exists(input), "demo bundle is not mounted")
  env <- new.env(parent = globalenv())
  sys.source(file.path(repository, "scripts", "prepare_precomputed_rds.R"), envir = env)
  env$.preparation_script_file <- file.path(repository, "scripts", "prepare_precomputed_rds.R")
  root <- tempfile("custom-preparation-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  analysis <- file.path(root, "sources", "tertiary_analysis")
  original_paths <- env$demo_paths(file.path(demo, "tertiary_analysis"))
  for (role in names(original_paths)) {
    relative <- substring(original_paths[[role]], nchar(demo) + 2L)
    destination <- file.path(root, "sources", relative)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (role != "input_provenance") stopifnot(file.symlink(original_paths[[role]], destination))
  }
  loadNamespace("pixelatorR")
  object <- readRDS(input)
  columns <- c(celltype_broad = "custom_type", cluster = "custom_cluster", donor = "custom_donor",
    condition = "custom_condition", sample_alias = "custom_sample")
  for (column in names(columns)) object@meta.data[[columns[[column]]]] <- object@meta.data[[column]]
  object@reductions$custom_umap <- object@reductions$umap.cca
  object@reductions$custom_pca <- object@reductions$pca
  object@assays$custom_assay <- object@assays$PNA
  SeuratObject::Key(object@assays$custom_assay) <- "customassay_"
  SeuratObject::LayerData(object, assay = "custom_assay", layer = "custom_layer") <-
    SeuratObject::LayerData(object, assay = "PNA", layer = "data")
  custom_input <- file.path(root, "input.rds")
  saveRDS(object, custom_input)
  provenance <- env$read_source_table(original_paths[["input_provenance"]])
  provenance$source_md5 <- unname(tools::md5sum(custom_input))
  utils::write.table(provenance, file.path(analysis, "input", "input_provenance.tsv"),
    sep = "\t", quote = TRUE, row.names = FALSE)
  output <- file.path(root, "prepared.rds")
  flags <- c("--umap", "custom_umap", "--pca", "custom_pca", "--celltype-column", "custom_type",
    "--cluster-column", "custom_cluster", "--donor-column", "custom_donor",
    "--condition-column", "custom_condition", "--sample-column", "custom_sample",
    "--assay", "custom_assay", "--abundance-layer", "custom_layer")
  result <- callr::rscript(file.path(repository, "scripts", "prepare_precomputed_rds.R"),
    cmdargs = c(custom_input, analysis, output, flags), wd = repository,
    env = c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", ARROW_NUM_THREADS = "1"),
    fail_on_status = FALSE, show = FALSE)
  expect_equal(result$status, 0L, info = result$stderr)
  if (result$status != 0L) return(invisible(NULL))
  data <- read_precomputed_rds(output)
  expect_identical(data$contract$defaults, list(reduction = "custom_umap", pca_reduction = "custom_pca",
    celltype_column = "custom_type", cluster_column = "custom_cluster", donor_column = "custom_donor",
    condition_column = "custom_condition", sample_column = "custom_sample", assay = "custom_assay",
    abundance_layer = "custom_layer"))
  expect_equal(data$metadata$celltype, as.character(object@meta.data$celltype_broad))
  expect_equal(data$metadata$sample_id, as.character(object@meta.data$sample_alias))
  expect_equal(data$abundance, SeuratObject::LayerData(object, assay = "custom_assay", layer = "custom_layer"))
  first_id <- data$contract$dataset_id
  env$build_precomputed_rds(custom_input, analysis, output)
  expect_false(identical(read_precomputed_rds(output)$contract$dataset_id, first_id))
  expect_error(env$build_precomputed_rds(custom_input, analysis, output, assay = "absent"), "selected assay")
  expect_error(env$build_precomputed_rds(custom_input, analysis, output, abundance_layer = "absent"), "selected abundance layer")
  expect_error(env$build_precomputed_rds(custom_input, analysis, output, reduction = "absent"), "selected UMAP")
  expect_error(env$build_precomputed_rds(custom_input, analysis, output, celltype_column = "absent"), "missing metadata")
  expect_error(env$build_precomputed_rds(custom_input, analysis, output, assay = NA_character_), "nonempty strings")
})
