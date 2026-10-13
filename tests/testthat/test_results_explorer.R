results_explorer_fixture <- function() {
  value <- precomputed_app_data(precomputed_fixture())
  value$metadata <- value$metadata[value$metadata$donor == "D1", , drop = FALSE]
  cells <- value$metadata$component
  value$abundance <- value$abundance[, cells, drop = FALSE]
  value$counts <- value$counts[, cells, drop = FALSE]
  value$embeddings <- lapply(value$embeddings, function(x) x[cells, , drop = FALSE])
  for (readout in PROXIOMEVIS_READOUTS) {
    rows <- value$contract$donor_contrasts[[readout]]
    rows <- rows[rows$donor == "D1", , drop = FALSE]
    rows$p_value <- ifelse(rows$marker_2 == "B", 0.1, 0.005)
    rows$p_adj <- rows$p_value * 2
    rows$method <- "stored component test"
    rows$replicate_unit <- "component"
    rows$inference <- "exploratory_component_level"
    rows$adjustment_family <- paste(readout, rows$celltype, sep = ":")
    value$contract$donor_contrasts[[readout]] <- rows
    value$contract$differential_results[[readout]] <- value$contract$differential_results[[readout]][FALSE, ]
    samples <- value$contract$sample_summaries[[readout]]
    value$contract$sample_summaries[[readout]] <- samples[samples$donor == "D1", ]
  }
  sources <- lapply(value$contract$donor_contrasts, function(rows) {
    rows$contrast <- "TRT_vs_NC"
    rows$target <- rows$group_a
    rows$reference <- rows$group_b
    rows
  })
  aligned <- expand.grid(marker_1 = c("A", "B"), marker_2 = c("A", "B"), stringsAsFactors = FALSE)
  aligned$donor <- "D1"
  aligned$celltype <- "B"
  aligned$panel <- "TRT_vs_NC"
  aligned$value <- c(0, 0.5, 0.5, NA_real_)
  aligned$p_adj <- c(0.01, 0.02, 0.02, NA_real_)
  aligned$status <- c("supported", "supported", "supported", "unavailable")
  profiles <- aligned
  profiles$panel <- "TRT"
  value$contract$report <- list(
    source_results = list(D1 = sources),
    eligibility = list(D1 = stats::setNames(lapply(PROXIOMEVIS_READOUTS,
      function(x) data.frame(celltype = c("B", "T", "Not tested"))), PROXIOMEVIS_READOUTS)),
    alignment = list(conditions = c("NC", "TRT", "Unrelated"), contrasts = "TRT_vs_NC",
      differences = aligned, profiles = profiles),
    settings = list(D1 = list(thresholds = list(fdr_threshold = 0.05, min_cells_per_group = 20L),
      spatial = list(profiling = list(summary_stat = "mean", missing_score_policy = "missing_not_zero"))))
  )
  value
}

test_that("single-donor explorer uses stored contrasts without cross-donor results", {
  fixture <- results_explorer_fixture()
  original <- fixture$contract
  current <- shiny::reactiveVal(fixture)
  manager <- new_hybrid_job_manager()
  shiny::testServer(precomputed_server, args = list(data = current, job_manager = manager), {
    session$setInputs(readout = "abundance", celltype = "B", comparison = "TRT vs NC",
      report_page = "Precomputed results", results_view = "Tables", include_controls = TRUE)
    expect_identical(scope()$donor, "D1")
    expect_identical(scope()$donors, "D1")
    expect_equal(nrow(readout_tables()$results), 0L)
    expect_equal(nrow(selected_results()), 2L)
    expect_equal(selected_results()$effect_size, c(1, 1))
    expect_identical(report_scope(), scope())
    expect_identical(spatial_scope()$donor, "D1")
    expect_identical(spatial_scope()$celltype, "B")
    expect_identical(spatial_scope()$panel, "TRT_vs_NC")
    donor_html <- as.character(output$donor_control$html)
    expect_match(donor_html, "Donor:", fixed = TRUE)
    expect_false(grepl("<select", donor_html, fixed = TRUE))
    expect_match(as.character(output$result_summary$html), "TRT − NC", fixed = TRUE)
    expect_match(as.character(output$result_summary$html), "not donor-population inference", fixed = TRUE)
    expect_match(as.character(output$methods_ui$html), "differences of observed medians", fixed = TRUE)
    expect_match(as.character(output$methods_ui$html), "fdr_threshold", fixed = TRUE)
    expect_match(as.character(output$methods_ui$html), "0.05", fixed = TRUE)
    expect_match(as.character(output$methods_ui$html), "min_cells_per_group", fixed = TRUE)
    expect_match(as.character(output$methods_ui$html), "missing_not_zero", fixed = TRUE)
    session$setInputs(results_view = "Volcano", celltype = "T")
    expect_identical(report_scope()$celltype, "T")
    expect_identical(spatial_scope()$celltype, "T")
    expect_true(all(selected_results()$celltype == "T"))
    session$setInputs(readout = "colocalization", results_view = "Spatial matrix",
      spatial_kind = "profiles", spatial_panel = "Unrelated")
    expect_setequal(spatial_scope()$panels, c("NC", "TRT"))
    expect_false(spatial_scope()$panel == "Unrelated")
    session$setInputs(celltype = "Not tested")
    expect_identical(scope()$celltype, "Not tested")
    expect_equal(nrow(selected_results()), 0L)
    expect_true("Not tested" %in% scope()$celltypes)
    expect_match(as.character(output$methods_ui$html), "No stored differential test", fixed = TRUE)
    expect_match(as.character(output$methods_ui$html), "fdr_threshold", fixed = TRUE)
    session$setInputs(celltype = PRECOMPUTED_ALL)
    expect_equal(nrow(selected_results()), 4L)
    expect_equal(manager$active, 0L)
    expect_identical(current()$contract, original)
    expect_false(session$isClosed())
  })
})

test_that("plot and table selections share exact pairs and transfer their scope to Cell 3D", {
  current <- shiny::reactiveVal(results_explorer_fixture())
  messages <- list()
  mock <- shiny::MockShinySession$new()
  mock$sendInputMessage <- function(id, message) messages[[sub("^.*-", "", id)]] <<- message
  shiny::testServer(precomputed_server, args = list(data = current), session = mock, {
    session$setInputs(readout = "colocalization", celltype = "B", comparison = "TRT vs NC",
      report_page = "Precomputed results", results_view = "Tables", include_controls = TRUE)
    select_feature('["C","A"]')
    expect_identical(selected_key(), report_feature_key("A", "C"))
    expect_identical(unique(detail_contrasts()$marker_2), "C")
    session$setInputs(results_view = "Volcano")
    expect_identical(selected_key(), report_feature_key("A", "C"))
    expect_identical(report_scope()$celltype, "B")
    rows <- display_results()
    session$setInputs(results_table_rows_selected = which(rows$marker_2 == "B"))
    expect_identical(selected_key(), report_feature_key("A", "B"))
    session$setInputs(results_table_rows_selected = which(rows$marker_2 == "C"))
    expect_identical(selected_key(), report_feature_key("A", "C"))
    select_feature('["A","A"]', from_spatial = TRUE)
    expect_identical(messages$readout$value, "self_proximity")
    expect_identical(selected_key(), report_feature_key("A", "A"))
    select_feature('["B","A"]', from_spatial = TRUE)
    expect_identical(messages$readout$value, "colocalization")
    expect_identical(selected_key(), report_feature_key("A", "B"))
    for (bad in c("not JSON", '["A"]', '["A","unknown"]', '[1,2]')) select_feature(bad)
    expect_identical(selected_key(), report_feature_key("A", "B"))
    session$setInputs(celltype = "T")
    expect_true(all(detail_contrasts()$celltype == "T"))
    expect_identical(selected_key(), report_feature_key("A", "B"))
    session$setInputs(explore_3d = 1L)
    request <- graph_selection()
    expect_identical(request$donor, "D1")
    expect_identical(request$celltype, "T")
    expect_identical(c(request$marker_1, request$marker_2), c("A", "B"))
    expect_setequal(request$conditions, c("TRT", "NC"))
    expect_null(input[["graph-display"]])
    expect_false(session$isClosed())
  })
})

test_that("filtered CSVs match FDR and table searches while all-row exports stay frozen", {
  fixture <- results_explorer_fixture()
  fixture$contract$donor_contrasts$abundance$is_control <-
    fixture$contract$donor_contrasts$abundance$marker_1 == "B"
  original <- fixture$contract
  current <- shiny::reactiveVal(fixture)
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(readout = "abundance", celltype = "B", comparison = "TRT vs NC",
      report_page = "Precomputed results", results_view = "Tables",
      include_controls = TRUE, fdr_filter = TRUE, fdr_cutoff = 0.05)
    expect_equal(nrow(display_results()), 1L)
    expect_equal(nrow(read.csv(output$download_filtered_results_csv)), 1L)
    expect_equal(nrow(read.csv(output$download_results_csv)), 2L)
    expect_equal(nrow(read.csv(output$download_report_volcano_csv)), 2L)
    session$setInputs(fdr_filter = FALSE, results_table_rows_all = 2L)
    filtered <- read.csv(output$download_filtered_results_csv)
    expect_identical(filtered$marker_1, display_results()$marker_1[2L])
    expect_equal(filtered$p_adj, display_results()$p_adj[2L])
    session$setInputs(results_table_rows_all = integer())
    empty <- read.csv(output$download_filtered_results_csv)
    expect_equal(nrow(empty), 0L)
    expect_identical(names(empty), names(display_results()))
    expect_equal(nrow(read.csv(output$download_results_csv)), 2L)
    expect_equal(nrow(read.csv(output$download_report_volcano_csv)), 2L)
    session$setInputs(results_table_rows_all = NULL, results_table_search = "no matching marker")
    expect_equal(nrow(read.csv(output$download_filtered_results_csv)), 0L)
    session$setInputs(results_table_search = "",
      results_table_search_columns = list("", "no matching marker"))
    expect_equal(nrow(read.csv(output$download_filtered_results_csv)), 0L)
    expect_equal(nrow(read.csv(output$download_results_csv)), 2L)
    expect_equal(nrow(read.csv(output$download_report_volcano_csv)), 2L)
    session$setInputs(results_table_search_columns = list("", ""))
    expect_equal(nrow(read.csv(output$download_filtered_results_csv)), 2L)
    session$setInputs(include_controls = FALSE)
    expect_equal(nrow(display_results()), 2L)
    expect_equal(nrow(read.csv(output$download_results_csv)), 2L)
    expect_equal(nrow(read.csv(output$download_report_volcano_csv)), 2L)
    session$setInputs(readout = "colocalization")
    session$setInputs(spatial_kind = "differences")
    spatial <- read.csv(output$download_spatial_csv)
    expect_equal(nrow(spatial), 4L)
    expect_identical(spatial$status, original$report$alignment$differences$status)
    expect_true(is.na(spatial$value[4L]))
    expect_identical(current()$contract, original)
  })
})

test_that("searched rows distinguish no search from zero matches and reject invalid indices", {
  rows <- data.frame(feature = c("A", "B", "C"))
  expect_identical(precomputed_searched_rows(rows, NULL), rows)
  expect_identical(precomputed_searched_rows(rows, NULL, search = "missing"), rows[FALSE, , drop = FALSE])
  expect_identical(precomputed_searched_rows(rows, NULL, column_search = list("", "missing")), rows[FALSE, , drop = FALSE])
  expect_identical(precomputed_searched_rows(rows, NULL, search = "", column_search = list("", "")), rows)
  expect_identical(precomputed_searched_rows(rows, integer()), rows[FALSE, , drop = FALSE])
  expect_identical(precomputed_searched_rows(rows, c(3L, 3L, NA, 0L, -1L, 8L)), rows[3L, , drop = FALSE])
})

test_that("Spatial matrix is disabled for abundance and requests a safe tab fallback", {
  current <- shiny::reactiveVal(results_explorer_fixture())
  messages <- list()
  availability <- NULL
  mock <- shiny::MockShinySession$new()
  mock$sendInputMessage <- function(id, message) messages[[id]] <<- message
  mock$sendCustomMessage <- function(type, message) {
    if (type == "precomputed-spatial-availability") availability <<- message
  }
  shiny::testServer(precomputed_server, args = list(data = current), session = mock, {
    session$setInputs(readout = "abundance", celltype = "B", comparison = "TRT vs NC",
      report_page = "Precomputed results", results_view = "Tables")
    expect_true(availability$disabled)
    expect_identical(availability$id, session$ns("results_view"))
    for (readout in c("self_proximity", "colocalization")) {
      session$setInputs(readout = readout)
      expect_false(availability$disabled)
    }
    session$setInputs(results_view = "Spatial matrix")
    session$setInputs(readout = "abundance")
    expect_true(availability$disabled)
    expect_identical(messages$results_view$value, "Tables")
    expect_match(spatial_report_ggplot()$labels$title, "Spatial")
    expect_match(ggplot2::ggplot_build(spatial_report_ggplot())$data[[1L]]$label, "[Ss]elf|[Cc]olo")
    expect_false(session$isClosed())
  })
})
