source("../../R/data_adapter.R")

test_that("default demo RDS path points to the Pixelator v4.1.1 Seurat object", {
  path <- default_demo_rds_path()

  expect_match(path, "pg_data_combined_fil\\.pixelator_v4\\.1\\.1\\.rds$")
  expect_true(file.exists(path))
})

test_that("default demo RDS path can be provided by environment", {
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(list(ok = TRUE), rds_path)
  withr::local_envvar(c(PROXIOME_DEMO_RDS = rds_path))

  expect_equal(default_demo_rds_path(), normalizePath(rds_path, mustWork = TRUE))
})

test_that("demo data can load from bundled cache without source RDS", {
  cache_path <- tempfile(fileext = ".rds")
  cached_data <- list(
    source = list(
      rds_path = "/missing/demo.rds",
      resolved_rds_path = "/missing/demo.rds",
      cache_schema_version = DEMO_CACHE_SCHEMA_VERSION
    ),
    qc = list(origin_metadata = data.frame(component = character(0)))
  )
  saveRDS(cached_data, cache_path)

  data <- load_demo_proxiome_data(
    rds_path = "/missing/demo.rds",
    cache_path = cache_path
  )

  expect_identical(data, cached_data)
})

test_that("demo data default can use cache when repository RDS is absent", {
  cache_path <- tempfile(fileext = ".rds")
  work_dir <- tempfile("cache-only-app-")
  dir.create(work_dir, recursive = TRUE)
  cached_data <- list(
    source = list(
      rds_path = "/missing/demo.rds",
      resolved_rds_path = "/missing/demo.rds",
      cache_schema_version = DEMO_CACHE_SCHEMA_VERSION
    ),
    qc = list(origin_metadata = data.frame(component = character(0)))
  )
  saveRDS(cached_data, cache_path)
  withr::local_envvar(c(PROXIOME_DEMO_RDS = NA))
  withr::local_dir(work_dir)

  data <- load_demo_proxiome_data(cache_path = cache_path)

  expect_identical(data, cached_data)
})

test_that("older cached data is upgraded with spatial sample aliases", {
  cache_path <- tempfile(fileext = ".rds")
  cached_data <- list(
    source = list(
      rds_path = "/missing/demo.rds",
      resolved_rds_path = "/missing/demo.rds",
      cache_schema_version = DEMO_CACHE_SCHEMA_VERSION - 1L
    ),
    metadata = data.frame(
      component = c("cell-a", "cell-b"),
      sample = c("sample-1", "sample-2"),
      condition = c("UNT", "PHA"),
      celltype_manual = c("CD8 T", "CD4 T"),
      stringsAsFactors = FALSE
    ),
    colocalization = data.frame(
      component = c("cell-a", "cell-b"),
      marker_1 = c("CD3e", "CD4"),
      marker_2 = c("CD8", "CD8"),
      log2_ratio = c(0.3, -0.2),
      stringsAsFactors = FALSE
    ),
    clustering = data.frame(
      component = c("cell-a", "cell-b"),
      marker_1 = c("CD3e", "CD4"),
      marker_2 = c("CD3e", "CD4"),
      log2_ratio = c(0.8, 0.5),
      stringsAsFactors = FALSE
    ),
    qc = list(origin_metadata = data.frame(component = character(0)))
  )
  saveRDS(cached_data, cache_path)

  data <- load_demo_proxiome_data(
    rds_path = "/missing/demo.rds",
    cache_path = cache_path
  )

  expect_equal(data$source$cache_schema_version, DEMO_CACHE_SCHEMA_VERSION)
  expect_equal(data$metadata$sample_alias, c("sample-1", "sample-2"))
  expect_equal(data$colocalization$sample_alias, c("sample-1", "sample-2"))
  expect_equal(data$clustering$sample_alias, c("sample-1", "sample-2"))
})

test_that("proxiome readout context names the three required readouts", {
  context <- proxiome_readout_context()

  expect_s3_class(context, "data.frame")
  expect_equal(context$readout, c("Abundance", "Clustering", "Colocalization"))
  expect_true(all(nzchar(context$description)))
})

test_that("demo marker selection preserves requested markers that are present", {
  markers <- select_demo_markers(
    available_markers = c("CD45", "CD8", "CD81", "CD82", "CD3e"),
    requested_markers = c("CD8", "CD4", "CD81", "CD82")
  )

  expect_equal(markers, c("CD8", "CD81", "CD82"))
})

test_that("marker selection can keep demo defaults or all user markers", {
  available_markers <- c("CD3e", "CD4", "CD8", "CD45", "M1", "M2")

  expect_equal(
    select_proxiome_markers(available_markers, marker_selection = "demo"),
    c("CD3e", "CD4", "CD8", "CD45")
  )
  expect_equal(
    select_proxiome_markers(available_markers, marker_selection = "all"),
    available_markers
  )
  expect_equal(
    select_proxiome_markers(available_markers, requested_markers = c("M2", "missing"), marker_selection = "all"),
    "M2"
  )
})

test_that("sparse abundance matrix conversion avoids dense assay coercion", {
  skip_if_not_installed("Matrix")
  if (!methods::isClass("ProxiomeNoDenseSparseMatrix")) {
    methods::setClass("ProxiomeNoDenseSparseMatrix", contains = "dgCMatrix")
  }
  methods::setMethod("as.matrix", "ProxiomeNoDenseSparseMatrix", function(x, ...) {
    stop("dense coercion blocked", call. = FALSE)
  })

  sparse <- Matrix::sparseMatrix(
    i = c(1, 1, 2),
    j = c(1, 2, 3),
    x = c(1, 2, 3),
    dims = c(2, 3),
    dimnames = list(c("CD3e", "CD8"), c("cell-a", "cell-b", "cell-c"))
  )
  sparse <- methods::as(sparse, "ProxiomeNoDenseSparseMatrix")

  abundance <- matrix_to_long(sparse, value_name = "abundance")

  expect_equal(
    abundance,
    data.frame(
      marker = rep(c("CD3e", "CD8"), times = 3),
      component = rep(c("cell-a", "cell-b", "cell-c"), each = 2),
      abundance = c(1, 0, 2, 0, 0, 3),
      stringsAsFactors = FALSE
    )
  )
})

test_that("abundance matrix conversion handles empty marker selections", {
  empty_layer <- matrix(
    numeric(),
    nrow = 0,
    ncol = 2,
    dimnames = list(character(), c("cell-a", "cell-b"))
  )

  expect_equal(
    matrix_to_long(empty_layer, value_name = "abundance"),
    data.frame(
      marker = character(),
      component = character(),
      abundance = numeric(),
      stringsAsFactors = FALSE
    )
  )
  expect_equal(
    matrix_layers_to_long(empty_layer, empty_layer),
    data.frame(
      marker = character(),
      component = character(),
      abundance = numeric(),
      count = numeric(),
      stringsAsFactors = FALSE
    )
  )
})

test_that("abundance and count layers are combined in component-major order", {
  data_layer <- matrix(
    c(1, 0, 2, 3),
    nrow = 2,
    dimnames = list(c("CD3e", "CD8"), c("cell-a", "cell-b"))
  )
  counts_layer <- matrix(
    c(10, 0, 20, 30),
    nrow = 2,
    dimnames = list(c("CD3e", "CD8"), c("cell-a", "cell-b"))
  )

  abundance <- matrix_layers_to_long(data_layer, counts_layer, max_block_entries = 2)

  expect_equal(
    abundance,
    data.frame(
      marker = rep(c("CD3e", "CD8"), times = 2),
      component = rep(c("cell-a", "cell-b"), each = 2),
      abundance = c(1, 0, 2, 3),
      count = c(10, 0, 20, 30),
      stringsAsFactors = FALSE
    )
  )
})

test_that("metadata gains a sample alias for spatial summaries", {
  metadata <- data.frame(
    component = c("cell-a", "cell-b"),
    sample = c("long_sample_1", "long_sample_2"),
    condition = c("UNT", "PHA"),
    stringsAsFactors = FALSE
  )

  normalized <- normalize_spatial_metadata(metadata)

  expect_equal(normalized$sample_alias, c("long_sample_1", "long_sample_2"))
})

test_that("proximity summarization separates self clustering from pair colocalization", {
  metadata <- data.frame(
    component = c("cell-a", "cell-b"),
    condition = c("UNT", "CD3CD28"),
    celltype_manual = c("CD8 T", "CD8 T"),
    sample_alias = c("sample-1", "sample-2"),
    stringsAsFactors = FALSE
  )

  proximity <- data.frame(
    component = c("cell-a", "cell-b", "cell-a", "cell-b"),
    marker_1 = c("CD81", "CD81", "CD53", "CD53"),
    marker_2 = c("CD81", "CD81", "CD58", "CD58"),
    log2_ratio = c(1.2, 0.4, 2.1, -0.3),
    join_count = c(10, 4, 12, 2),
    join_count_expected_mean = c(5, 3, 4, 2),
    stringsAsFactors = FALSE
  )

  readouts <- summarize_proximity_readouts(proximity, metadata)

  expect_named(readouts, c("clustering", "clustering_summary", "colocalization", "colocalization_summary"))
  expect_equal(unique(readouts$clustering$marker), "CD81")
  expect_equal(unique(readouts$colocalization$marker_pair), "CD53 / CD58")
  expect_equal(readouts$colocalization$sample_alias, c("sample-1", "sample-2"))
  expect_equal(readouts$clustering_summary$n_cells, c(1L, 1L))
  expect_equal(readouts$colocalization_summary$n_cells, c(1L, 1L))
})

test_that("proximity and abundance summaries avoid base R merge and aggregate hot paths", {
  summarize_abundance_source <- paste(deparse(body(summarize_abundance)), collapse = "\n")
  summarize_proximity_source <- paste(deparse(body(summarize_proximity_readouts)), collapse = "\n")
  aggregate_source <- paste(deparse(body(aggregate_numeric_readout)), collapse = "\n")

  expect_false(grepl("\\bmerge\\s*\\(", summarize_abundance_source))
  expect_false(grepl("\\bmerge\\s*\\(", summarize_proximity_source))
  expect_false(grepl("\\baggregate\\s*\\(", aggregate_source))
})

test_that("stored assay proximity is used without calling pixelatorR ProximityScores", {
  if (!methods::isClass("ProxiomeTestAssay")) {
    methods::setClass("ProxiomeTestAssay", slots = c(proximity = "data.frame"))
  }

  assay <- methods::new(
    "ProxiomeTestAssay",
    proximity = data.frame(
      component = c("cell-a", "cell-a", "cell-b"),
      marker_1 = c("CD3e", "CD3e", "CD4"),
      marker_2 = c("CD3e", "CD8", "CD4"),
      log2_ratio = c(1.2, -0.4, 0.8),
      join_count = c(4, 2, 3),
      stringsAsFactors = FALSE
    )
  )

  proximity <- stored_assay_proximity(assay, markers = c("CD3e", "CD8"))

  expect_equal(nrow(proximity), 2L)
  expect_equal(proximity$marker_1, c("CD3e", "CD3e"))
  expect_equal(proximity$marker_2, c("CD3e", "CD8"))

  adapter_source <- paste(readLines("../../R/data_adapter.R", warn = FALSE), collapse = "\n")
  expect_false(grepl('require_namespace("pixelatorR")', adapter_source, fixed = TRUE))
  expect_false(grepl("pixelatorR::ProximityScores", adapter_source, fixed = TRUE))
})

test_that("user RDS schema inspection reports assay, metadata, embeddings, proximity, dimensions, and cache estimate", {
  if (!methods::isClass("ProxiomeSchemaAssay")) {
    methods::setClass("ProxiomeSchemaAssay", slots = c(proximity = "data.frame"))
  }
  if (!methods::isClass("ProxiomeSchemaObject")) {
    methods::setClass(
      "ProxiomeSchemaObject",
      slots = c(
        assays = "list",
        meta.data = "data.frame",
        reductions = "list",
        active.assay = "character"
      )
    )
  }

  assay <- methods::new(
    "ProxiomeSchemaAssay",
    proximity = data.frame(
      component = c("cell-a", "cell-b"),
      marker_1 = c("CD3e", "CD4"),
      marker_2 = c("CD3e", "CD8"),
      log2_ratio = c(0.5, -0.2),
      stringsAsFactors = FALSE
    )
  )
  rownames(assay@proximity) <- NULL
  object <- methods::new(
    "ProxiomeSchemaObject",
    assays = list(PNA = assay),
    meta.data = data.frame(
      condition = c("UNT", "CD3CD28"),
      celltype_manual = c("CD8 T", "CD4 T"),
      row.names = c("cell-a", "cell-b"),
      stringsAsFactors = FALSE
    ),
    reductions = list(umap = matrix(1:4, ncol = 2, dimnames = list(c("cell-a", "cell-b"), c("UMAP_1", "UMAP_2")))),
    active.assay = "PNA"
  )
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(object, rds_path)

  schema <- inspect_user_rds_schema(rds_path)

  expect_equal(schema$assay$name, "PNA")
  expect_true(schema$assay$available)
  expect_equal(schema$marker_count, 3L)
  expect_equal(schema$cell_count, 2L)
  expect_true(all(c("condition", "celltype_manual") %in% schema$metadata$present))
  expect_true("sample_alias" %in% schema$metadata$missing)
  expect_equal(schema$embeddings$names, "umap")
  expect_true(schema$embeddings$has_two_dimensional)
  expect_true(schema$proximity$available)
  expect_equal(schema$proximity$row_count, 2L)
  expect_true(schema$estimated_cache_size_bytes > 0)

  report <- format_user_rds_schema_report(schema)
  expect_match(report, "Expected assay: PNA", fixed = TRUE)
  expect_match(report, "Metadata columns", fixed = TRUE)
  expect_match(report, "celltype_manual", fixed = TRUE)
  expect_match(report, "Embeddings: umap", fixed = TRUE)
  expect_match(report, "Stored proximity: available", fixed = TRUE)
  expect_match(report, "3 markers", fixed = TRUE)
  expect_match(report, "2 cells", fixed = TRUE)
  expect_match(report, "Estimated cache", fixed = TRUE)
})

test_that("QC payload preserves original metadata and filter counts", {
  filtered_metadata <- data.frame(
    sample = c("s1", "s2"),
    n_umi = c(12000, 18000),
    isotype_fraction = c(0.0002, 0.0004),
    row.names = c("cell-b", "cell-c"),
    stringsAsFactors = FALSE
  )
  origin_metadata <- data.frame(
    sample = c("s1", "s1", "s2"),
    n_umi = c(9000, 12000, 18000),
    row.names = c("cell-a", "cell-b", "cell-c"),
    stringsAsFactors = FALSE
  )
  filter_counts <- data.frame(
    step = c("00_loaded", "01_after_n_umi_filter", "02_after_tau_filter"),
    n_cells = c(3, 2, 2),
    stringsAsFactors = FALSE
  )

  qc <- build_qc_payload(filtered_metadata, origin_metadata, filter_counts)

  expect_named(qc, c("origin_metadata", "filtered_metadata", "filter_counts"))
  expect_equal(qc$origin_metadata$component, c("cell-a", "cell-b", "cell-c"))
  expect_equal(qc$filtered_metadata$component, c("cell-b", "cell-c"))
  expect_equal(qc$filter_counts$fraction_loaded, c(1, 2 / 3, 2 / 3))
  expect_equal(qc$filter_counts$step_label, c("Loaded", "After n umi filter", "After tau filter"))
})

test_that("QC payload falls back when notebook misc tables are absent", {
  filtered_metadata <- data.frame(
    sample = c("s1", "s2"),
    n_umi = c(12000, 18000),
    row.names = c("cell-b", "cell-c"),
    stringsAsFactors = FALSE
  )

  qc <- build_qc_payload(filtered_metadata)

  expect_equal(nrow(qc$origin_metadata), 2L)
  expect_equal(qc$filter_counts$step, c("00_loaded", "99_filtered"))
  expect_equal(qc$filter_counts$n_cells, c(2L, 2L))
})

test_that("QC distribution choices use available numeric metadata columns", {
  choices <- available_qc_distribution_choices(
    data.frame(
      sample = c("s1", "s2"),
      n_umi = c(1, 2),
      n_edges = c(3, 4),
      reads_in_component = c(5, 6),
      isotype_fraction = c(0.1, 0.2),
      tau = c(0.3, 0.4),
      condition = c("UNT", "PHA")
    )
  )

  expect_named(choices, c("UMIs", "Edges", "Reads in component", "Isotype fraction", "Tau"))
  expect_equal(unname(choices), c("n_umi", "n_edges", "reads_in_component", "isotype_fraction", "tau"))
})

test_that("differential numeric readout compares two conditions within cell types", {
  readout <- data.frame(
    marker = rep(c("CD3e", "CD8"), each = 8),
    condition = rep(c("UNT", "UNT", "PHA", "PHA"), times = 4),
    celltype_manual = rep(c("B", "B", "B", "B", "T", "T", "T", "T"), times = 2),
    abundance = c(
      0.2, 0.4, 1.8, 2.0,
      0.7, 0.8, 0.8, 0.9,
      1.2, 1.4, 0.6, 0.8,
      2.1, 2.3, 2.0, 2.2
    ),
    stringsAsFactors = FALSE
  )

  differential <- calculate_differential_readout(
    readout,
    feature_cols = "marker",
    value_col = "abundance",
    group_a = "UNT",
    group_b = "PHA",
    celltype_filter = "B",
    min_cells = 2
  )

  expect_s3_class(differential, "data.frame")
  expect_named(
    differential,
    c(
      "marker", "celltype_manual", "group_a", "group_b",
      "mean_a", "mean_b", "median_a", "median_b",
      "n_a", "n_b", "effect_size", "p_value", "p_adj",
      "direction", "is_significant"
    )
  )
  expect_equal(differential$celltype_manual, c("B", "B"))
  expect_equal(
    differential$effect_size[match("CD3e", differential$marker)],
    -1.6
  )
  expect_equal(
    differential$direction[match("CD8", differential$marker)],
    "Higher in UNT"
  )
})

test_that("differential readout pools selected cell types by default", {
  readout <- data.frame(
    marker = "CD3e",
    condition = rep(c("UNT", "UNT", "PHA", "PHA"), times = 2),
    celltype_manual = rep(c("B", "B", "B", "B", "T", "T", "T", "T"), each = 1),
    abundance = c(0.2, 0.4, 1.8, 2.0, 0.7, 0.8, 0.8, 0.9),
    stringsAsFactors = FALSE
  )

  differential <- calculate_differential_readout(
    readout,
    feature_cols = "marker",
    value_col = "abundance",
    group_a = "UNT",
    group_b = "PHA",
    celltype_filter = c("B", "T"),
    min_cells = 2
  )

  expect_equal(nrow(differential), 1L)
  expect_equal(differential$celltype_manual, "Pooled cell types")
  expect_equal(differential$n_a, 4)
  expect_equal(differential$n_b, 4)
  expect_equal(differential$effect_size, -0.8)
})

test_that("differential readout can stratify selected cell types when requested", {
  readout <- data.frame(
    marker = "CD3e",
    condition = rep(c("UNT", "UNT", "PHA", "PHA"), times = 2),
    celltype_manual = rep(c("B", "B", "B", "B", "T", "T", "T", "T"), each = 1),
    abundance = c(0.2, 0.4, 1.8, 2.0, 0.7, 0.8, 0.8, 0.9),
    stringsAsFactors = FALSE
  )

  differential <- calculate_differential_readout(
    readout,
    feature_cols = "marker",
    value_col = "abundance",
    group_a = "UNT",
    group_b = "PHA",
    celltype_filter = c("B", "T"),
    stratify_by_celltype = TRUE,
    min_cells = 2
  )

  expect_equal(nrow(differential), 2L)
  expect_equal(sort(differential$celltype_manual), c("B", "T"))
})

test_that("differential readout preserves marker-pair feature columns", {
  readout <- data.frame(
    marker_pair = rep(c("CD3e / CD8", "CD4 / CD8"), each = 4),
    marker_1 = rep(c("CD3e", "CD4"), each = 4),
    marker_2 = "CD8",
    condition = rep(c("UNT", "UNT", "PHA", "PHA"), times = 2),
    celltype_manual = "T",
    log2_ratio = c(0.1, 0.2, 1.1, 1.2, 2.2, 2.0, 1.4, 1.5),
    stringsAsFactors = FALSE
  )

  differential <- calculate_differential_readout(
    readout,
    feature_cols = c("marker_pair", "marker_1", "marker_2"),
    value_col = "log2_ratio",
    group_a = "UNT",
    group_b = "PHA",
    min_cells = 2
  )

  expect_equal(nrow(differential), 2L)
  expect_true(all(c("marker_pair", "marker_1", "marker_2") %in% names(differential)))
  expect_equal(
    differential$effect_size[match("CD3e / CD8", differential$marker_pair)],
    -1.0
  )
})

test_that("differential readout uses median difference for the plotted effect", {
  readout <- data.frame(
    marker = rep("CD81", 6),
    condition = c("UNT", "UNT", "UNT", "PHA", "PHA", "PHA"),
    celltype_manual = "CD8 T",
    log2_ratio = c(0, 0, 100, 2, 2, 2),
    stringsAsFactors = FALSE
  )

  differential <- calculate_differential_readout(
    readout,
    feature_cols = "marker",
    value_col = "log2_ratio",
    group_a = "UNT",
    group_b = "PHA",
    min_cells = 2
  )

  expect_equal(differential$mean_a - differential$mean_b, 31.333, tolerance = 0.001)
  expect_equal(differential$effect_size, -2)
})
