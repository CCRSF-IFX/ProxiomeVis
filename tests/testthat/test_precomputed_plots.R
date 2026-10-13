test_that("donor views preserve stored exploratory statistics", {
  contract <- precomputed_fixture()@misc$proxiomevis
  rows <- contract$donor_contrasts$abundance
  rows$p_value <- seq_len(nrow(rows)) / 100
  rows$p_adj <- seq_len(nrow(rows)) / 50
  rows$adjustment_family <- "stored-family"
  rows$method <- "stored component test"
  rows$replicate_unit <- "component"
  rows$inference <- "exploratory_component_level"
  contract$donor_contrasts$abundance <- rows

  selected <- precomputed_result_rows_for_donor(contract, "abundance", "D1")
  expect_true(nrow(selected) > 0L)
  expect_true(all(selected$donor == "D1"))
  expect_equal(selected$p_value, rows$p_value[rows$donor == "D1"])
  expect_equal(selected$p_adj, rows$p_adj[rows$donor == "D1"])
  expect_identical(unique(selected$replicate_unit), "component")
  expect_identical(unique(selected$inference), "exploratory_component_level")
})

test_that("heatmap completion keeps unavailable values missing and uses symmetric limits", {
  results <- precomputed_fixture()@misc$proxiomevis$differential_results$abundance
  second <- results
  second$analysis_id <- "abundance:ALT_vs_NC"
  second$comparison <- "ALT vs NC"
  second$group_a <- "TRT"
  second <- second[-1L, , drop = FALSE]
  combined <- rbind(results, second)

  prepared <- precomputed_prepare_heatmap(combined, top_n = 20L)
  expect_equal(nrow(prepared), length(unique(combined$feature_id)) * 4L)
  expect_true(anyNA(prepared$effect_size))

  plot <- plot_precomputed_heatmap(combined)
  expect_s3_class(plot, "ggplot")
  limit <- max(abs(combined$effect_size), na.rm = TRUE)
  expect_equal(plot$scales$scales[[1L]]$limits, c(-limit, limit))
})

test_that("heatmap clicks use row-by-column pair keys, including missing tiles", {
  contrasts <- precomputed_fixture()@misc$proxiomevis$donor_contrasts$colocalization
  contrasts <- contrasts[contrasts$donor == "D1", , drop = FALSE]
  for (columns in c(1L, 2L)) {
    rows <- if (columns == 1L) contrasts[contrasts$celltype == "B", , drop = FALSE] else {
      contrasts[!(contrasts$celltype == "T" & contrasts$marker_2 == "C"), , drop = FALSE]
    }
    rows <- rows[rev(seq_len(nrow(rows))), , drop = FALSE]
    plot <- plot_precomputed_heatmap(rows)
    original_data <- plot$data
    original_layers <- ggplot2::ggplot_build(plot)$data
    interactive <- plot
    i <- match(interactive$data$feature_id, rows$feature_id)
    interactive$data$key <- report_feature_key(rows$marker_1[i], rows$marker_2[i])
    widget <- precomputed_heatmap_plotly(interactive, source = "linked-summary")
    built <- plotly::plotly_build(plotly::plotly_build(widget))
    trace <- Filter(function(trace) identical(trace$type, "heatmap"), built$x$data)[[1L]]
    key_matrix <- do.call(rbind, lapply(trace$key, as.character))
    expect_equal(dim(key_matrix), c(2L, columns))
    positions <- cbind(match(as.integer(interactive$data$feature_label), trace$y),
                       match(as.integer(interactive$data$analysis_label), trace$x))
    expect_identical(unname(key_matrix[positions]), interactive$data$key)
    expect_true(all(vapply(as.vector(key_matrix), function(key) {
      pair <- jsonlite::fromJSON(key)
      is.character(pair) && length(pair) == 2L && pair[1L] == "A" && pair[2L] %in% c("B", "C")
    }, logical(1))))
    serialized <- jsonlite::fromJSON(plotly::plotly_json(built, jsonedit = FALSE), simplifyVector = FALSE)
    json_trace <- Filter(function(trace) identical(trace$type, "heatmap"), serialized$data)[[1L]]
    expect_length(json_trace$key, 2L)
    expect_true(all(lengths(json_trace$key) == columns))
    expect_identical(json_trace$key[[1L]][[1L]], key_matrix[1L, 1L])
    if (columns == 2L) {
      missing <- which(is.na(interactive$data$effect_size))
      expect_length(missing, 1L)
      expect_identical(key_matrix[positions[missing, , drop = FALSE]], report_feature_key("A", "C"))
      expect_true(is.na(plot$data$effect_size[missing]))
    }
    expect_identical(built$x$source, "linked-summary")
    expect_true("plotly_click" %in% built$x$shinyEvents)
    expect_identical(plot$data, original_data)
    expect_equal(ggplot2::ggplot_build(plot)$data, original_layers)
  }
  single <- contrasts[contrasts$marker_2 == "B", , drop = FALSE]
  for (columns in c(1L, 2L)) {
    plot <- plot_precomputed_heatmap(single[seq_len(columns), , drop = FALSE])
    plot$data$key <- rep(report_feature_key("A", "B"), nrow(plot$data))
    built <- plotly::plotly_build(plotly::plotly_build(precomputed_heatmap_plotly(plot)))
    trace <- Filter(function(trace) identical(trace$type, "heatmap"), built$x$data)[[1L]]
    expect_true(trace$`_isSimpleKey`)
    expect_identical(as.character(trace$key), report_feature_key("A", "B"))
    serialized <- jsonlite::fromJSON(plotly::plotly_json(built, jsonedit = FALSE), simplifyVector = FALSE)
    json_trace <- Filter(function(trace) identical(trace$type, "heatmap"), serialized$data)[[1L]]
    expect_identical(unlist(json_trace$key), report_feature_key("A", "B"))
  }
  expect_s3_class(precomputed_heatmap_plotly(precomputed_empty_plot("No stored results", "Summary")), "plotly")
})

test_that("embedding helpers align cells and support sparse marker abundance", {
  app_data <- precomputed_app_data(precomputed_fixture())
  prepared <- precomputed_embedding_data(
    app_data$metadata,
    app_data$embeddings$umap,
    color_by = "abundance",
    abundance = app_data$abundance,
    marker = "A"
  )
  expect_identical(prepared$component, app_data$metadata$component)
  expect_true(is.numeric(prepared$color_value))
  expect_s3_class(
    plot_precomputed_embedding(app_data$metadata, app_data$embeddings$umap, "celltype"),
    "ggplot"
  )
})

test_that("QC table places filtering stages in columns without filling missing counts", {
  counts <- data.frame(sample_id = c("s2", "s1", "s1", "s2", "s1"), donor = "D1",
    condition = c("TRT", "NC", "NC", "TRT", "NC"),
    stage = c("00_loaded", "02_after_tau_filter", "00_loaded", "03_after_isotype_filter", "01_after_n_umi_filter"),
    n_cells_display = c(100, 60, 90, 0, 70))
  wide <- precomputed_qc_table(counts)
  expect_identical(names(wide), c("sample_id", "donor", "condition", "00_loaded",
    "01_after_n_umi_filter", "02_after_tau_filter", "03_after_isotype_filter"))
  expect_equal(nrow(wide), 2L)
  expect_equal(wide$`00_loaded`, c(100, 90))
  expect_equal(wide$`01_after_n_umi_filter`, c(NA, 70))
  expect_equal(wide$`02_after_tau_filter`, c(NA, 60))
  expect_equal(wide$`03_after_isotype_filter`, c(0, NA))
  expect_error(precomputed_qc_table(rbind(counts, counts[1L, ])), "uniquely")
})

test_that("QC retention and unavailable paired source values have explicit plots", {
  retention <- data.frame(
    stage = rep(c("00 loaded", "01 retained"), 2),
    sample_id = rep(c("D1_unstim", "D1_PHA"), each = 2),
    donor = "D1",
    condition = rep(c("unstim", "PHA"), each = 2),
    n_cells_display = c(100, 80, 120, 90),
    celltype = "All",
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_precomputed_qc_counts(retention), "ggplot")

  contrast <- precomputed_fixture()@misc$proxiomevis$donor_contrasts$colocalization[1, ]
  contrast$value_a <- NA_real_
  contrast$value_b <- NA_real_
  plot <- plot_precomputed_condition_detail(contrast)
  expect_s3_class(plot, "ggplot")
  expect_match(ggplot2::ggplot_build(plot)$data[[1L]]$label, "not stored")
})

test_that("the report UI groups stored views under shared navigation and selection", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  ui <- as.character(precomputed_ui("report"))
  expect_match(ui, 'id="report-results_view"', fixed = TRUE)
  expect_match(ui, 'data-value="Tables"', fixed = TRUE)
  expect_false(grepl('data-value="Summary"', ui, fixed = TRUE))
  expect_false(grepl('id="report-effect_heatmap"', ui, fixed = TRUE))
  expect_match(ui, "Spatial patterns", fixed = TRUE)
  expect_match(ui, 'id="report-donor_control"', fixed = TRUE)
  expect_match(ui, 'id="report-explore_3d"', fixed = TRUE)
  expect_match(ui, 'id="report-download_filtered_results_csv"', fixed = TRUE)
  expect_match(ui, 'id="report-download_spatial_csv"', fixed = TRUE)
  expect_false(grepl('id="report-report_donor"', ui, fixed = TRUE))
  expect_false(grepl('id="report-spatial_donor"', ui, fixed = TRUE))
  for (id in c("effect_detail", "report_volcano_plot", "volcano_distribution", "spatial_report_plot",
               "results_table", "fdr_controls")) {
    expect_equal(sum(gregexpr(paste0('id="report-', id, '"'), ui, fixed = TRUE)[[1L]] > 0L), 1L)
  }
  expect_match(ui, "Exploratory component-level test", fixed = TRUE)
  expect_false(grepl('data-value="Volcano rows"', ui, fixed = TRUE))
  expect_false(grepl('id="report-report_volcano_table"', ui, fixed = TRUE))
  expect_false(grepl('id="report-download_filtered_report_csv"', ui, fixed = TRUE))
})

test_that("report cell-type choices retain source scopes that were not tested", {
  object <- precomputed_fixture()
  contract <- object@misc$proxiomevis
  contract$report <- list(
    eligibility = list(
      D1 = list(abundance = data.frame(celltype = c("B", "Not tested")))
    )
  )
  donor_rows <- contract$donor_contrasts$abundance[
    contract$donor_contrasts$abundance$donor == "D1",
    ,
    drop = FALSE
  ]
  choices <- precomputed_report_celltypes(contract, object@meta.data, "D1", "abundance", donor_rows)
  expect_identical(choices, c("B", "Not tested"))
})
