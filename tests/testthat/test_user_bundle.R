source(file.path("..", "..", "R", "bundle_builder.R"))

user_seurat_fixture <- function() {
  metadata <- data.frame(participant = "participant-7", treatment = rep(c("vehicle", "treated"), each = 4L),
    annotation = rep(c("B cell", "T cell"), 4L), cluster_id = rep(c("0", "1"), 4L),
    isotype_fraction = 0.001, row.names = paste0("app-cell-", 1:8))
  metadata$library <- paste0("library-", metadata$treatment)
  counts <- matrix(rep(c(2, 3, 5), 8L), nrow = 3L, dimnames = list(c("CD19", "CD81", "CD82"), rownames(metadata)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE), assay = "protein", meta.data = metadata)
  SeuratObject::LayerData(object, assay = "protein", layer = "normalized") <- log1p(counts)
  for (name in c("umap_custom", "pca_custom")) {
    key <- paste0(gsub("_", "", name), "_")
    coords <- matrix(seq_len(16L), ncol = 2L, dimnames = list(rownames(metadata), paste0(key, 1:2)))
    object[[name]] <- SeuratObject::CreateDimReducObject(coords, key = key, assay = "protein")
  }
  object
}

user_bundle_options <- function() list(title = "Own data", assay = "protein", abundance_layer = "normalized",
  reduction = "umap_custom", pca_reduction = "pca_custom", celltype_column = "annotation", cluster_column = "cluster_id",
  donor_column = "participant", condition_column = "treatment", sample_column = "library")

test_that("an annotated Seurat object builds a valid bundle without saved comparisons", {
  skip_if_not_installed("SeuratObject")
  object <- user_seurat_fixture()
  expect_null(object@misc$proxiomevis)
  data <- do.call(user_bundle_data, c(list(object = object), user_bundle_options()))
  expect_silent(validate_precomputed_data(data))
  for (readout in PROXIOMEVIS_READOUTS) {
    expect_equal(nrow(data$contract$donor_contrasts[[readout]]), 0L)
    expect_equal(nrow(data$contract$differential_results[[readout]]), 0L)
  }
  expect_equal(nrow(data$contract$sample_summaries$abundance), 12L)
  expect_equal(data$contract$sample_summaries$abundance$value[1:3], log1p(c(2, 3, 5)))
  expect_equal(data$contract$qc$cell_counts$n_cells, c(4, 4))
  expect_identical(data$contract$qc$cell_counts$stage, rep("00_loaded", 2L))
  normalized <- precomputed_app_data(c(list(bundle_version = 2L), data))
  expect_identical(unique(normalized$metadata$donor), "participant-7")
  expect_equal(as.matrix(normalized$counts), as.matrix(SeuratObject::LayerData(object, assay = "protein", layer = "counts")))
})

test_that("saved sections are optional, validated, and retain supplied values", {
  skip_if_not_installed("SeuratObject")
  object <- precomputed_fixture()
  saved <- object@misc$proxiomevis[c("sample_summaries", "donor_contrasts", "differential_results")]
  object@misc$proxiomevis <- NULL
  data <- user_bundle_data(object, cluster_column = "cluster", saved_results = saved)
  expect_identical(data$contract$donor_contrasts, saved$donor_contrasts)
  expect_identical(data$contract$differential_results, saved$differential_results)
  partial <- user_bundle_data(object, cluster_column = "cluster",
    saved_results = list(sample_summaries = list(colocalization = saved$sample_summaries$colocalization)))
  expect_identical(partial$contract$sample_summaries$colocalization, saved$sample_summaries$colocalization)
  expect_true(nrow(partial$contract$sample_summaries$abundance) > 0L)
  saved$donor_contrasts$abundance$effect_size[1] <- 100
  expect_error(user_bundle_data(object, cluster_column = "cluster", saved_results = saved), "invalid donor or contrast")
  expect_error(user_bundle_data(object, cluster_column = "cluster", saved_results = list(results = list())), "named result sections")
})

test_that("metadata, layers, cell mappings, and QC are checked before publication", {
  skip_if_not_installed("SeuratObject")
  options <- user_bundle_options()
  object <- user_seurat_fixture()
  expect_error(do.call(user_bundle_data, c(list(object = object), modifyList(options, list(donor_column = "missing")))), "metadata column")
  expect_error(do.call(user_bundle_data, c(list(object = object), modifyList(options, list(abundance_layer = "missing")))), "normalized layer")
  mapping <- data.frame(component = colnames(object), original_id = paste0("raw-", 1:8), sample_id = object$library)
  data <- do.call(user_bundle_data, c(list(object = object, cell_sources = mapping), options))
  expect_identical(data$cell_sources, mapping)
  mapping$sample_id[1] <- "wrong-sample"
  expect_error(do.call(user_bundle_data, c(list(object = object, cell_sources = mapping), options)), "original sample")
  qc <- data$contract$qc$cell_counts
  qc$n_cells[1] <- -1
  expect_error(do.call(user_bundle_data, c(list(object = object, qc = qc), options)), "nonnegative integer")
})

test_that("preparation publishes a new directory atomically and never overwrites input", {
  skip_if_not_installed("SeuratObject")
  root <- tempfile("user-bundle-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  input <- file.path(root, "source.rds")
  saveRDS(user_seurat_fixture(), input)
  hash <- unname(tools::md5sum(input))
  destination <- file.path(root, "prepared")
  output <- do.call(prepare_user_bundle, c(list(input_rds = input, output_dir = destination), user_bundle_options()))
  data <- read_precomputed_rds(output)
  expect_identical(readRDS(output)$bundle_version, 2L)
  expect_identical(data$contract$title, "Own data")
  expect_identical(unname(tools::md5sum(input)), hash)
  expect_error(prepare_user_bundle(input, destination), "new output directory")
  expect_error(prepare_user_bundle(input, input), "new output directory")
  invalid <- file.path(root, "invalid")
  expect_error(do.call(prepare_user_bundle, c(list(input_rds = input, output_dir = invalid),
    modifyList(user_bundle_options(), list(reduction = "absent")))), "missing reduction")
  expect_false(file.exists(invalid))
  expect_false(any(startsWith(list.files(root, all.files = TRUE), ".invalid-")))
})

test_that("the own-data CLI runs outside the repository with custom mappings", {
  skip_if_not_installed("SeuratObject")
  root <- tempfile("user-cli-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  input <- file.path(root, "source.rds")
  saveRDS(user_seurat_fixture(), input)
  script <- normalizePath(file.path("..", "..", "scripts", "prepare_user_bundle.R"))
  previous <- setwd(root)
  on.exit(setwd(previous), add = TRUE)
  options <- user_bundle_options()
  flags <- unlist(Map(function(name, value) c(paste0("--", gsub("_", "-", name)), value), names(options), options), use.names = FALSE)
  output <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(script, input, "prepared", flags)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  data <- read_precomputed_rds(file.path(root, "prepared", "analysis.rds"))
  expect_identical(data$contract$defaults$celltype_column, "annotation")
  expect_identical(data$contract$defaults$abundance_layer, "normalized")
  expect_identical(unique(data$metadata$donor), "participant-7")
})

test_that("own-data PXL preparation makes a matching proximity sidecar", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("duckdb")
  skip_if_not_installed("arrow")
  root <- tempfile("user-pxl-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  object <- user_seurat_fixture()
  input <- file.path(root, "source.rds")
  saveRDS(object, input)
  mapping <- data.frame(component = colnames(object), original_id = paste0("raw-", 1:8), sample_id = object$library)
  samples <- unique(mapping$sample_id)
  paths <- stats::setNames(file.path(root, paste0(samples, ".pxl")), samples)
  for (sample in samples) {
    connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = paths[[sample]])
    rows <- data.frame(component = mapping$original_id[mapping$sample_id == sample],
      marker_1 = "CD82", marker_2 = "CD19", join_count = 4, join_count_expected_mean = 2)
    rows <- rows[-1L, , drop = FALSE]
    DBI::dbWriteTable(connection, "proximity", rows)
    DBI::dbDisconnect(connection, shutdown = TRUE)
  }
  output <- do.call(prepare_user_bundle, c(list(input_rds = input, output_dir = file.path(root, "prepared"),
    cell_sources = mapping, pxl_by_sample = paths), user_bundle_options()))
  data <- read_precomputed_rds(output)
  sidecar <- resolve_hybrid_asset(output, data$assets$proximity)
  prepared <- validated_proximity_partitions(data, sidecar, data$metadata$component)
  expect_length(prepared$paths, 2L)
  rows <- as.data.frame(arrow::read_parquet(prepared$paths[[1L]]))
  expect_true(all(rows$component %in% colnames(object)))
  expect_true(all(rows$marker_1 == "CD19" & rows$marker_2 == "CD82" & rows$log2_ratio == 1))
  expect_identical(names(prepared$manifest$donor_cutoffs), "participant-7")
  profile <- data$contract$sample_summaries$colocalization
  expect_equal(nrow(profile), 4L)
  expect_true(all(profile$value[profile$celltype == "B cell"] == 0.5))
  expect_true(all(profile$value[profile$celltype == "T cell"] == 1))
  expect_true(all(profile$n_total_cells == 2))
  expect_true(all(profile$n_scored_cells[profile$celltype == "B cell"] == 1))
  expect_true(all(profile$metric == "population_mean_log2_ratio_missing_zero"))

  connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = paths[[1L]])
  DBI::dbExecute(connection, "DROP TABLE proximity")
  DBI::dbDisconnect(connection, shutdown = TRUE)
  failed <- file.path(root, "failed")
  expect_error(do.call(prepare_user_bundle, c(list(input_rds = input, output_dir = failed,
    cell_sources = mapping, pxl_by_sample = paths), user_bundle_options())), "proximity")
  expect_false(file.exists(failed))
  expect_false(any(startsWith(list.files(root, all.files = TRUE), ".failed-")))
})

test_that("a trusted exported comparison can be saved as a separate validated bundle", {
  skip_if_not_installed("SeuratObject")
  root <- tempfile("user-import-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  input <- file.path(root, "source.rds")
  saveRDS(user_seurat_fixture(), input)
  bundle <- do.call(prepare_user_bundle, c(list(input_rds = input, output_dir = file.path(root, "prepared")), user_bundle_options()))
  data <- read_precomputed_rds(bundle)
  request <- validate_hybrid_request(data, list(readout = "self_proximity", donor = "participant-7",
    celltype = "B cell", target = "treated", reference = "vehicle", markers = "CD82", min_cells_per_group = 1L))
  entry <- list(id = "session-1", label = "My comparison", request = request, result = list(results = data.frame(
    feature_id = "CD82", marker_1 = "CD82", marker_2 = "CD82", effect_size = 0.5,
    p_value = 0.01, p_adj = 0.03, n_cells_a = 2, n_cells_b = 2), provenance = list(method = "fixture",
      dataset_id = request$dataset_id, bundle_md5 = request$bundle_md5, request = request)))
  comparison <- file.path(root, "comparison.rds")
  saveRDS(list(result_version = 1L, entry = entry), comparison)
  output <- import_user_comparison(bundle, comparison, file.path(dirname(bundle), "with_comparison.rds"))
  imported <- read_precomputed_rds(output)
  expect_identical(imported$contract$donor_contrasts$self_proximity$p_adj, 0.03)
  expect_identical(imported$contract$sample_summaries, data$contract$sample_summaries)
  expect_error(import_user_comparison(bundle, comparison, output), "new output RDS")
  entry$request$bundle_md5 <- "wrong"
  saveRDS(list(result_version = 1L, entry = entry), comparison)
  expect_error(import_user_comparison(bundle, comparison, file.path(root, "bad.rds")), "bundle changed")
  expect_false(file.exists(file.path(root, "bad.rds")))
})

test_that("spatial population means retain zero-observation sample denominators", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("arrow")
  root <- tempfile("user-profile-zero-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  metadata <- data.frame(component = paste0("cell-", 1:69),
    sample_id = rep(c("observed", "zero", "other-condition", "other-celltype", "other-donor"), c(10, 30, 20, 5, 4)),
    donor = c(rep("donor-7", 65), rep("donor-8", 4)),
    condition = c(rep("vehicle", 40), rep("treated", 20), rep("vehicle", 9)),
    celltype = c(rep("B", 60), rep("T", 5), rep("B", 4)))
  observed <- data.frame(component = metadata$component[1:10], marker_1 = "CD19", marker_2 = "CD82",
    log2_ratio = 1, min_marker_fraction = 0.2)
  samples <- unique(metadata$sample_id)
  manifest <- list(partitions = data.frame(sample_id = samples, file = paste0(seq_along(samples), ".parquet")),
    donor_cutoffs = c("donor-7" = 0.001, "donor-8" = 0.001))
  for (i in seq_along(samples)) {
    rows <- if (i == 1L) observed else observed[FALSE, , drop = FALSE]
    arrow::write_parquet(rows, file.path(root, manifest$partitions$file[[i]]))
  }
  data <- list(metadata = metadata, contract = empty_user_results())
  profiles <- user_spatial_profiles(data, root, manifest)
  rows <- profiles$colocalization
  expect_setequal(rows$sample_id, c("observed", "zero", "other-condition"))
  expect_equal(rows$value[match(c("observed", "zero", "other-condition"), rows$sample_id)], c(1, 0, 0))
  expect_equal(rows$n_total_cells[match(c("observed", "zero", "other-condition"), rows$sample_id)], c(10, 30, 20))
  expect_equal(rows$n_scored_cells[match(c("observed", "zero", "other-condition"), rows$sample_id)], c(10, 0, 0))
  expect_equal(nrow(profiles$self_proximity), 0L)
  expect_true(all(rows$marker_1 == "CD19" & rows$marker_2 == "CD82"))
  data$contract$sample_summaries[names(profiles)] <- profiles
  alignment <- report_spatial_alignment(data$contract)
  pooled <- alignment$profiles[alignment$profiles$panel == "vehicle" & alignment$profiles$celltype == "B", ]
  expect_equal(pooled$value, 0.25)
  expect_equal(pooled$pct_detected, 0.25)
  treated <- alignment$profiles[alignment$profiles$panel == "treated", ]
  expect_equal(treated$value, 0)
  expect_equal(treated$pct_detected, 0)

  arrow::write_parquet(observed[FALSE, , drop = FALSE], file.path(root, manifest$partitions$file[[1L]]))
  unavailable <- user_spatial_profiles(data, root, manifest)
  expect_equal(nrow(unavailable$colocalization), 0L)
  expect_equal(nrow(unavailable$self_proximity), 0L)
})

test_that("comparison import round-trips overlapping metadata mappings", {
  skip_if_not_installed("SeuratObject")
  root <- tempfile("user-import-aliases-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  object <- user_seurat_fixture()
  object$sample_id <- object$participant
  object$donor <- object$library
  object$celltype <- object$cluster_id
  object$cluster <- object$annotation
  input <- file.path(root, "source.rds")
  saveRDS(object, input)
  mapping <- data.frame(component = colnames(object), original_id = paste0("raw-", 1:8), sample_id = object$library)
  options <- modifyList(user_bundle_options(), list(donor_column = "sample_id", sample_column = "donor",
    celltype_column = "cluster", cluster_column = "celltype"))
  bundle <- do.call(prepare_user_bundle, c(list(input_rds = input, output_dir = file.path(root, "prepared"),
    cell_sources = mapping), options))
  original_hash <- unname(tools::md5sum(c(input, bundle)))
  data <- read_precomputed_rds(bundle)
  request <- validate_hybrid_request(data, list(readout = "self_proximity", donor = "participant-7",
    celltype = "B cell", target = "treated", reference = "vehicle", markers = "CD82", min_cells_per_group = 1L))
  entry <- list(id = "session-1", label = "Mapped comparison", request = request, result = list(results = data.frame(
    feature_id = "CD82", marker_1 = "CD82", marker_2 = "CD82", effect_size = 0.5,
    p_value = 0.01, p_adj = 0.03, n_cells_a = 2, n_cells_b = 2), provenance = list(method = "fixture",
      dataset_id = request$dataset_id, bundle_md5 = request$bundle_md5, request = request)))
  comparison <- file.path(root, "comparison.rds")
  saveRDS(list(result_version = 1L, entry = entry), comparison)
  output <- import_user_comparison(bundle, comparison, file.path(dirname(bundle), "with_comparison.rds"))
  imported <- read_precomputed_rds(output)
  expect_identical(imported$metadata, data$metadata)
  expect_identical(imported$cell_sources, data$cell_sources)
  expect_identical(imported$contract$sample_summaries, data$contract$sample_summaries)
  expect_identical(imported$contract$donor_contrasts$self_proximity$donor, "participant-7")
  expect_identical(imported$contract$defaults$donor_column, "donor")
  expect_identical(imported$contract$defaults$sample_column, "sample_id")
  expect_identical(hybrid_dataset_fingerprint(imported), hybrid_dataset_fingerprint(data))
  expect_identical(unname(tools::md5sum(c(input, bundle))), original_hash)
})
