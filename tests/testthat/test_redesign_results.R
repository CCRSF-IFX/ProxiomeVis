redesign_results_fixture <- function() {
  data <- precomputed_app_data(precomputed_fixture())
  pairs <- expand.grid(marker_1 = c("A", "B"), marker_2 = c("A", "B"), stringsAsFactors = FALSE)
  pairs$donor <- "D1"
  pairs$celltype <- "B"
  pairs$panel <- "TRT vs NC"
  pairs$value <- c(.1, .2, .2, -.3)
  pairs$p_adj <- c(.01, .2, .2, .03)
  pairs$pct_detected <- 1
  pairs$status <- "supported"
  profiles <- rbind(transform(pairs, panel = "NC"), transform(pairs, panel = "TRT"))
  data$contract$report <- list(alignment = list(
    cfg = list(significance_size_cap = 4, significance_threshold = .05),
    markers = data.table::data.table(celltype = "B", marker = c("B", "A"),
      tested_clustering = TRUE, tested_colocalization = TRUE),
    conditions = c("NC", "TRT"), contrasts = "TRT vs NC", donors = "D1",
    profiles = data.table::as.data.table(profiles), differences = data.table::as.data.table(pairs),
    limits = list(profiles = 1, differences = .5), metric = "Mean log2 ratio"))
  data
}

test_that("shared spatial selection validates lists without changing prepared values", {
  expect_identical(precomputed_spatial_selection(c("A", "B"), "prepared", "", c("A", "B", "C"), c("B", "A")), c("B", "A"))
  expect_identical(precomputed_spatial_selection(c("B", "A"), "alphabetical", "", c("A", "B"), c("B", "A")), c("A", "B"))
  expect_identical(precomputed_spatial_selection(NULL, "custom", "C, B\nA", c("A", "B", "C"), c("B", "A")), c("C", "B", "A"))
  expect_error(precomputed_spatial_selection(NULL, "custom", "A,A", c("A", "B"), "A"), "only once")
  expect_error(precomputed_spatial_selection(NULL, "custom", "D", c("A", "B"), "A"), "Unknown proteins")
  expect_error(precomputed_spatial_selection(character(), "prepared", "", c("A", "B"), "A"), "at least one")
})

test_that("Apply aligns all spatial panels and exports; switching views preserves order", {
  fixture <- redesign_results_fixture()
  original <- fixture$contract
  current <- shiny::reactiveVal(fixture)
  manager <- new_hybrid_job_manager()
  shiny::testServer(precomputed_server, args = list(data = current, job_manager = manager), {
    session$setInputs(readout = "colocalization", donor = "D1", celltype = "B", comparison = "TRT vs NC",
      report_page = "Results", results_view = "Spatial patterns", spatial_kind = "differences")
    expect_identical(spatial_markers(), c("B", "A"))
    session$setInputs(spatial_order = "custom", spatial_custom_order = "C,B,A")
    expect_identical(spatial_markers(), c("B", "A"))
    session$setInputs(spatial_apply = 1L, spatial_layout = "compare")
    expect_identical(spatial_markers(), c("C", "B", "A"))
    plots <- spatial_comparison_plots()
    for (plot in plots) {
      expect_identical(plot$scales$get_scales("x")$limits, c("C", "B", "A"))
      expect_identical(plot$scales$get_scales("y")$limits, c("A", "B", "C"))
      expect_true(all(plot$data$status[plot$data$marker_1 == "C"] == "unavailable"))
    }
    expect_identical(plots$reference$scales$get_scales("fill")$limits, plots$target$scales$get_scales("fill")$limits)
    expect_equal(plots$difference$scales$get_scales("fill")$limits, c(-.5, .5))
    combined <- read.csv(output$download_spatial_compare_csv)
    expect_equal(nrow(combined), 27L)
    expect_identical(unique(combined$marker_1), c("C", "B", "A"))
    expect_setequal(combined$view, c("reference", "target", "difference"))
    expect_true(all(is.na(combined$value[combined$marker_1 == "C"])))
    expect_match(paste(readLines(output$download_spatial_compare_svg, warn = FALSE), collapse = ""), "<svg", fixed = TRUE)
    expect_identical(read.csv(output$download_spatial_order)$protein, c("C", "B", "A"))
    session$setInputs(spatial_kind = "profiles", spatial_panel = "NC", spatial_layout = "single")
    expect_identical(unique(spatial_rows()$marker_1), c("C", "B", "A"))
    expect_equal(nrow(read.csv(output$download_spatial_csv)), 9L)
    session$setInputs(celltype = "T", donor = "D2", readout = "self_proximity")
    expect_identical(spatial_markers(), c("C", "B", "A"))
    session$setInputs(celltype = PRECOMPUTED_ALL)
    expect_error(spatial_rows(), "Select one cell type", class = "shiny.silent.error")
    expect_error(spatial_comparison_rows(), "Select one cell type", class = "shiny.silent.error")
    session$setInputs(spatial_custom_order = "unknown", spatial_apply = 2L)
    expect_identical(spatial_markers(), c("C", "B", "A"))
    expect_match(output$spatial_selection_status, "Not applied")
    expect_identical(current()$contract, original)
    expect_identical(manager$active, 0L)
    current(fixture)
    session$flushReact()
    session$setInputs(donor = "D1", celltype = "B", spatial_reset = 1L)
    expect_identical(spatial_markers(), c("B", "A"))
  })
})

test_that("profiles are available without any stored differential result", {
  fixture <- redesign_results_fixture()
  for (section in c("donor_contrasts", "differential_results"))
    fixture$contract[[section]] <- lapply(fixture$contract[[section]], function(x) x[FALSE, ])
  fixture$contract$report$alignment$differences <- fixture$contract$report$alignment$differences[FALSE, ]
  fixture$contract$report$alignment$contrasts <- character()
  current <- shiny::reactiveVal(fixture)
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(readout = "colocalization", donor = "D1", celltype = "B", spatial_kind = "profiles", spatial_panel = "NC")
    expect_equal(nrow(selected_results()), 0L)
    expect_setequal(spatial_scope()$panels, c("NC", "TRT"))
    expect_equal(nrow(spatial_rows()), 4L)
    expect_true(any(spatial_rows()$status == "supported"))
    expect_identical(spatial_report_ggplot()$scales$get_scales("x")$limits, c("B", "A"))
    expect_true(nrow(selected_samples()) > 0L)
    expect_identical(spatial_scope()$reference, "NC")
    expect_identical(spatial_scope()$target, "TRT")
    expect_equal(nrow(spatial_comparison_plots()$reference$data), 4L)
    expect_equal(nrow(spatial_comparison_plots()$target$data), 4L)
    session$setInputs(spatial_reference_condition = "TRT", spatial_target_condition = "NC")
    expect_identical(spatial_scope()$reference, "TRT")
    expect_identical(spatial_scope()$target, "NC")
    select_feature('["B","A"]', from_spatial = TRUE)
    session$setInputs(explore_3d = 1L)
    expect_setequal(graph_selection()$conditions, c("NC", "TRT"))
  })
})

test_that("a profile condition need not belong to the selected differential comparison", {
  fixture <- redesign_results_fixture()
  fixture$metadata$condition[fixture$metadata$donor == "D1" & fixture$metadata$celltype == "B"][1L] <- "Third"
  alignment <- fixture$contract$report$alignment
  extra <- data.table::copy(alignment$profiles[panel == "NC"])
  extra$panel <- "Third"
  fixture$contract$report$alignment$profiles <- rbind(alignment$profiles, extra)
  current <- shiny::reactiveVal(fixture)
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(readout = "colocalization", donor = "D1", celltype = "B", comparison = "TRT vs NC",
      spatial_kind = "profiles", spatial_panel = "Third")
    expect_setequal(spatial_scope()$panels, c("NC", "TRT", "Third"))
    expect_identical(spatial_scope()$panel, "Third")
    expect_true(all(spatial_rows()$panel == "Third"))
    expect_identical(spatial_scope()$reference, "NC")
    expect_identical(spatial_scope()$target, "TRT")
  })
})

test_that("report eligibility cannot hide cell types with profiles but no saved test", {
  fixture <- redesign_results_fixture()
  fixture$contract$report$eligibility <- list(D1 = list(colocalization = data.frame(celltype = "B")))
  extra <- data.table::copy(fixture$contract$report$alignment$profiles)
  extra$celltype <- "T"
  fixture$contract$report$alignment$profiles <- rbind(fixture$contract$report$alignment$profiles, extra)
  current <- shiny::reactiveVal(fixture)
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(readout = "colocalization", donor = "D1", celltype = "T",
      spatial_kind = "profiles", spatial_panel = "NC")
    expect_identical(scope()$celltype, "T")
    expect_true("T" %in% scope()$celltypes)
    expect_true(any(spatial_rows()$status == "supported"))
  })
})

test_that("new comparison opens in the shared viewer without replacing the loaded dataset", {
  fixture <- hybrid_fixture()
  fixture$bundle_path <- "/trusted/analysis.rds"
  fixture$bundle_md5 <- paste(rep("a", 32), collapse = "")
  fixture$assets$proximity <- list(path = "proximity")
  fixture$embeddings <- list(umap = matrix(seq_len(168L), ncol = 2,
    dimnames = list(fixture$metadata$component, c("UMAP_1", "UMAP_2"))))
  fixture$contract$defaults$reduction <- "umap"
  pending <- new.env(parent = emptyenv())
  hybrid <- hybrid_analysis_server
  environment(hybrid) <- list2env(list(start_hybrid_job = function(path, request, code_dir, manager) {
    pending$request <- request
    promise <- promises::promise(function(resolve, reject) pending$resolve <- resolve)
    list(promise = promise, cancel = function() NULL)
  }), parent = environment(hybrid_analysis_server))
  server <- precomputed_server
  environment(server) <- list2env(list(hybrid_analysis_server = hybrid), parent = environment(precomputed_server))
  current <- shiny::reactiveVal(fixture)
  shiny::testServer(server, args = list(data = current), {
    session$setInputs(readout = "colocalization", donor = "D1", celltype = "B")
    select_feature('["A","B"]')
    session$setInputs(`analysis-donor` = "D1", `analysis-celltype` = "B", `analysis-readout` = "abundance",
      `analysis-target` = "TRT", `analysis-reference` = "NC", `analysis-acknowledge` = TRUE)
    session$setInputs(`analysis-run` = 1L)
    request <- pending$request
    result <- list(results = data.frame(feature_id = c("A", "B", "C"), marker_1 = c("A", "B", "C"),
      marker_2 = c("A", "B", "C"), effect_size = c(.7, -.2, .1), p_value = c(.01, .3, .5),
      p_adj = c(.03, .45, .5), n_cells_a = 42, n_cells_b = 42),
      provenance = list(dataset_id = request$dataset_id, bundle_md5 = request$bundle_md5, request = request))
    pending$resolve(result)
    for (i in seq_len(30L)) { later::run_now(.01); session$flushReact() }
    expect_length(analysis$history(), 1L)
    expect_identical(active_result_id(), "prepared")
    session$setInputs(`analysis-open_result` = 1L)
    expect_false(identical(active_result_id(), "prepared"))
    expect_identical(selected_results()$effect_size, result$results$effect_size)
    expect_identical(scope()$donor, "D1")
    expect_identical(scope()$celltype, "B")
    expect_identical(current()$contract, fixture$contract)
    expect_equal(nrow(report_volcano_ggplot()$data), 3L)
    expect_identical(selected_key(), '["A","A"]')
    expect_equal(nrow(detail_distribution()), 84L)
    expect_equal(nrow(read.csv(output$download_results_csv)), 3L)
    session$setInputs(`analysis-readout` = "colocalization", `analysis-run` = 2L)
    request <- pending$request
    pairs <- hybrid_feature_table(request)
    result$results[names(pairs)] <- pairs
    result$provenance <- list(dataset_id = request$dataset_id, bundle_md5 = request$bundle_md5,
      request = request, filters = list(donor_isotype_fraction_cutoff = .001))
    pending$resolve(result)
    for (i in seq_len(30L)) { later::run_now(.01); session$flushReact() }
    session$setInputs(`analysis-open_result` = 2L)
    expect_length(analysis$history(), 2L)
    expect_identical(scope()$readout, "colocalization")
    expect_identical(selected_key(), '["A","B"]')
    select_feature('["A","A"]', from_spatial = TRUE)
    session$setInputs(readout = "colocalization")
    expect_identical(selected_key(), '["A","A"]')
    expect_identical(scope()$readout, "colocalization")
    expect_equal(nrow(detail_contrasts()), 0L)
    expect_identical(current()$contract, fixture$contract)
    session$setInputs(result_source = "prepared")
    expect_identical(view_data()$contract, fixture$contract)
    replacement <- fixture
    replacement$bundle_md5 <- paste(rep("b", 32), collapse = "")
    current(replacement)
    session$flushReact()
    expect_length(analysis$history(), 0L)
    expect_identical(active_result_id(), "prepared")
  })
})
