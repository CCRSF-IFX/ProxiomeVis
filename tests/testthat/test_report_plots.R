source(file.path("..", "..", "R", "report_plot_sources.R"))
source(file.path("..", "..", "R", "report_plots.R"))

report_plot_fixture <- function() {
  readouts <- c("abundance", "self_proximity", "colocalization")
  canonical <- function(readout) data.frame(
    donor = "D1", celltype = "B", comparison = "TRT vs NC",
    analysis_id = paste0(readout, ":TRT_vs_NC"),
    group_a = "TRT", group_b = "NC", stringsAsFactors = FALSE
  )
  abundance <- data.frame(
    marker = c("A", "B"), p = c(0.01, 0.2), p_adj = c(0.02, 0.3),
    difference = c(1.5, -0.5), target = "TRT", reference = "NC",
    contrast = "TRT_vs_NC", celltype = "B", stringsAsFactors = FALSE
  )
  self <- data.frame(
    marker_1 = c("A", "B"), marker_2 = c("A", "B"),
    p = c(0.005, 0.1), p_adj = c(0.01, 0.2), diff_median = c(1.2, -0.4),
    target = "TRT", reference = "NC", contrast = "TRT_vs_NC", celltype = "B",
    stringsAsFactors = FALSE
  )
  coloc <- data.frame(
    marker_1 = "A", marker_2 = "B", pair = "A/B", p = 0.01, p_adj = 0.025,
    diff_median = 0.8, target = "TRT", reference = "NC",
    contrast = "TRT_vs_NC", celltype = "B", stringsAsFactors = FALSE
  )
  abundance_support <- data.frame(
    contrast = "TRT_vs_NC", target = "TRT", reference = "NC",
    celltype = c("B", "T"), n_cells_target = c(30L, 2L),
    n_cells_reference = c(30L, 1L), status = c("tested", "insufficient_cells"),
    meets_configured_minimum = c(TRUE, FALSE), reason = c(NA, "too few cells"),
    stringsAsFactors = FALSE
  )
  spatial_support <- data.frame(
    contrast = "TRT_vs_NC", target = "TRT", reference = "NC",
    celltype = c("B", "T"), n_cells_target = c(30L, 2L),
    n_cells_reference = c(30L, 1L), min_cells_per_group = 20L,
    status = c("tested", "insufficient_cells"), reason = c(NA, "too few cells"),
    stringsAsFactors = FALSE
  )
  profile <- data.frame(
    celltype = rep("B", 4L), condition = rep(c("NC", "TRT"), each = 2L),
    marker_1 = rep(c("A", "B"), 2L), marker_2 = rep(c("A", "B"), 2L),
    mean_log2_ratio = c(0.1, 0.2, 0.4, 0.6), stringsAsFactors = FALSE
  )
  markers <- data.table::data.table(
    celltype = rep("B", 2L), marker = c("A", "B"),
    tested_clustering = c(TRUE, TRUE), tested_colocalization = c(TRUE, TRUE)
  )
  matrix_rows <- expand.grid(
    marker_1 = c("A", "B"), marker_2 = c("A", "B"), stringsAsFactors = FALSE
  )
  aligned <- data.table::as.data.table(matrix_rows)
  aligned[, `:=`(
    celltype = "B", donor = "D1", panel = "TRT_vs_NC", status = "supported",
    value = c(1.2, 0.8, 0.8, -0.4), p_adj = c(0.01, 0.025, 0.025, 0.2),
    pct_detected = 1
  )]
  aligned_profiles <- data.table::copy(aligned)
  aligned_profiles[, `:=`(panel = "TRT", value = c(0.4, 0.2, 0.2, 0.6))]
  alignment <- list(
    cfg = list(significance_size_cap = 6, significance_threshold = 0.05),
    conditions = "TRT", contrasts = "TRT_vs_NC", donors = "D1",
    markers = markers, profiles = aligned_profiles, differences = aligned,
    limits = list(profiles = 1, differences = 2), metric = "mean log2 ratio"
  )
  heatmap_matrix <- matrix(
    c(2, 1, 3, 1.5), nrow = 2L, byrow = TRUE,
    dimnames = list(c("A", "B"), c("D1_TRT", "D1_NC"))
  )
  list(
    donor_contrasts = stats::setNames(lapply(readouts, canonical), readouts),
    report = list(
      source_results = list(D1 = list(
        abundance = abundance, self_proximity = self, colocalization = coloc
      )),
      eligibility = list(D1 = list(
        abundance = abundance_support,
        self_proximity = spatial_support,
        colocalization = spatial_support
      )),
      profiles = list(D1 = list(
        abundance = data.frame(), self_proximity = profile, colocalization = profile
      )),
      settings = list(D1 = list(
        spatial = list(profiling = list(summary_stat = "mean")),
        clustering = list(
          proximity_metric = "log2_ratio", min_cells_per_group = 20L,
          fdr_threshold = 0.05, volcano_difference_threshold = 1,
          volcano_x_limit = NULL, volcano_label_max_overlaps = 5L,
          show_reference_lines = FALSE, reference_line_alpha = 0.5
        ),
        colocalization = list(
          min_cells_per_group = 20L, fdr_threshold = 0.05,
          volcano_difference_threshold = 1, volcano_x_limit = NULL,
          volcano_label_extremes = 2L, show_reference_lines = FALSE,
          reference_line_alpha = 0.5
        )
      )),
      alignment = alignment,
      abundance_heatmaps = list(D1 = list(
        manifest = data.frame(
          contrast = "TRT_vs_NC", celltype = "B", plot_status = "plotted",
          selection_samples = "D1_TRT;D1_NC", fdr_threshold = 0.05,
          difference_threshold = 0.5, color_min = 0, color_max = 3,
          stringsAsFactors = FALSE
        ),
        sample_counts = data.frame(
          celltype = "B", sample = c("D1_TRT", "D1_NC"), n_cells = c(30L, 30L)
        ),
        samples = data.frame(
          sample = c("D1_TRT", "D1_NC"), condition = c("TRT", "NC"), donor = "D1"
        ),
        matrices = list(TRT_vs_NC = list(B = heatmap_matrix))
      ))
    )
  )
}

test_that("report volcanoes use frozen rows and preserve stored probabilities", {
  contract <- report_plot_fixture()
  original <- contract$report$source_results
  for (readout in c("abundance", "self_proximity", "colocalization")) {
    plot <- plot_report_volcano(contract, readout, "D1", "B", "TRT vs NC")
    expect_s3_class(plot, "ggplot")
    built <- ggplot2::ggplot_build(plot)
    expect_gt(length(built$data), 0L)
  }
  expect_identical(contract$report$source_results, original)
  expect_equal(
    plot_report_volcano(contract, "self_proximity", "D1", "B", "TRT vs NC")$data$p_adj,
    contract$report$source_results$D1$self_proximity$p_adj
  )
})

test_that("interactive volcanoes expose every point without changing report exports", {
  contract <- report_plot_fixture()
  for (readout in c("abundance", "self_proximity", "colocalization")) {
    plot <- plot_report_volcano(contract, readout, "D1", "B", "TRT vs NC")
    original_data <- ggplot2::ggplot_build(plot)$data
    original_layers <- plot$layers
    expect_warning(widget <- report_volcano_plotly(plot, source = "results-volcano"), NA)
    built <- plotly::plotly_build(widget)
    points <- Filter(function(trace) identical(trace$mode, "markers"), built$x$data)
    hover <- unlist(lapply(points, function(trace) trace$text))
    expect_true(any(grepl("Marker.*A", hover)))
    expect_true(any(grepl("Marker.*B", hover)))
    expect_true(all(grepl("Adjusted P value:", hover, fixed = TRUE)))
    expect_true(all(grepl("Target: TRT", hover, fixed = TRUE)))
    keys <- unlist(lapply(points, function(trace) trace$key), use.names = FALSE)
    expected_keys <- if (identical(readout, "colocalization")) {
      report_feature_key("B", "A")
    } else report_feature_key(c("A", "B"), c("A", "B"))
    expect_setequal(keys, expected_keys)
    expect_setequal(unlist(lapply(points, function(trace) trace$customdata)), expected_keys)
    expect_identical(built$x$source, "results-volcano")
    expect_true("plotly_click" %in% built$x$shinyEvents)
    expect_true(built$x$config$scrollZoom)
    expect_true(built$x$config$displayModeBar)
    expect_identical(built$x$layout$dragmode, "zoom")
    labels <- vapply(built$x$layout$annotations, function(annotation) annotation$text, character(1))
    expect_true("-log10(adjusted P)" %in% labels)
    expect_identical(built$x$config$doubleClick, "reset")
    expect_identical(plot$layers, original_layers)
    expect_equal(ggplot2::ggplot_build(plot)$data, original_data)
    expect_false(any(vapply(plot$layers, function(layer) "text" %in% names(layer$mapping), logical(1))))
    empty <- plot_report_volcano(contract, readout, "D1", "T", "TRT vs NC")
    expect_warning(expect_s3_class(report_volcano_plotly(empty), "plotly"), NA)
  }
})

test_that("feature keys preserve exact unordered marker identities", {
  expect_identical(report_feature_key(c("B", "A"), c("A", "A")), c('["A","B"]', '["A","A"]'))
  expect_identical(report_feature_key(factor("B"), factor("A")), report_feature_key("A", "B"))
  expect_identical(report_feature_key(character(), character()), character())
  expect_identical(report_feature_key(NA_character_, "A"), NA_character_)
  expect_error(report_feature_key(c("A", "B"), "A"))
  expect_false(identical(report_feature_key("A/B", "C"), report_feature_key("A", "B/C")))
  marker <- 'CD3"\\<β>'
  expect_setequal(jsonlite::fromJSON(report_feature_key(marker, "CD45")), c(marker, "CD45"))
})

test_that("volcano hover escapes labels and shows stored zero P values", {
  contract <- report_plot_fixture()
  contract$report$source_results$D1$abundance$marker[1L] <- "A<B"
  contract$report$source_results$D1$abundance$p_adj[1L] <- 0
  plot <- plot_report_volcano(contract, "abundance", "D1", "B", "TRT vs NC")
  widget <- plotly::plotly_build(report_volcano_plotly(plot))
  points <- Filter(function(trace) identical(trace$mode, "markers"), widget$x$data)
  hover <- unlist(lapply(points, function(trace) trace$text))
  expect_true(any(grepl("Marker: A&lt;B", hover, fixed = TRUE)))
  expect_true(any(grepl("Adjusted P value: 0$", hover)))
})

test_that("unsupported source scopes keep their not-tested facet", {
  contract <- report_plot_fixture()
  plot <- plot_report_volcano(contract, "abundance", "D1", "T", "TRT vs NC")
  expect_s3_class(plot, "ggplot")
  labels <- unlist(lapply(ggplot2::ggplot_build(plot)$data, function(layer) layer$label), use.names = FALSE)
  expect_true(any(grepl("Test unavailable", labels, fixed = TRUE)))
})

test_that("aligned spatial plots retain shared marker order and scales", {
  contract <- report_plot_fixture()
  plot <- plot_report_spatial(contract, "D1", "B", "differences", "TRT_vs_NC")
  expect_s3_class(plot, "ggplot")
  expect_identical(levels(plot$data$marker_1), c("A", "B"))
  expect_equal(plot$scales$get_scales("fill")$limits, c(-2, 2))
  expect_equal(plot$data$p_adj, contract$report$alignment$differences$p_adj)
})

test_that("interactive spatial matrices retain stored values, scales, order and report exports", {
  contract <- report_plot_fixture()
  for (kind in c("profiles", "differences")) {
    panel <- if (identical(kind, "profiles")) "TRT" else "TRT_vs_NC"
    plot <- plot_report_spatial(contract, "D1", "B", kind, panel)
    original_data <- ggplot2::ggplot_build(plot)$data
    original_layers <- plot$layers
    expect_warning(widget <- report_spatial_plotly(plot, source = "results-spatial"), NA)
    built <- plotly::plotly_build(widget)
    points <- Filter(function(trace) identical(trace$mode, "markers") && !is.null(trace$key), built$x$data)
    hover <- unlist(lapply(points, function(trace) trace$text), use.names = FALSE)
    expect_length(hover, 4L)
    expect_true(all(grepl("Availability: supported", hover, fixed = TRUE)))
    expect_true(all(grepl("Fraction detected: 1", hover, fixed = TRUE)))
    expect_true(all(grepl("Adjusted P value:", hover, fixed = TRUE)))
    keys <- unlist(lapply(points, function(trace) trace$key), use.names = FALSE)
    expect_setequal(keys, report_feature_key(c("A", "A", "B"), c("A", "B", "B")))
    expect_setequal(unlist(lapply(points, function(trace) trace$customdata)), keys)
    expect_identical(built$x$source, "results-spatial")
    expect_true("plotly_click" %in% built$x$shinyEvents)
    expect_identical(built$x$layout$xaxis$ticktext, c("A", "B"))
    expect_identical(built$x$layout$yaxis$ticktext, c("B", "A"))
    expect_identical(built$x$layout$xaxis$constrain, "domain")
    expect_identical(built$x$layout$yaxis$constrain, "domain")
    expect_identical(built$x$layout$xaxis$scaleanchor, FALSE)
    expect_identical(built$x$layout$yaxis$scaleanchor, "x")
    expect_equal(built$x$layout$yaxis$scaleratio, 1)
    expect_equal(built$x$layout$xaxis$range, c(0.4, 2.6))
    expect_equal(built$x$layout$yaxis$range, c(0.4, 2.6))
    colorbars <- Filter(function(trace) !is.null(trace$marker$colorbar), built$x$data)
    limits <- contract$report$alignment$limits[[kind]]
    expect_equal(range(as.numeric(colorbars[[1L]]$marker$colorbar$ticktext)), c(-limits, limits))
    expect_true(built$x$config$scrollZoom)
    expect_identical(built$x$layout$dragmode, "zoom")
    expect_identical(plot$layers, original_layers)
    expect_equal(ggplot2::ggplot_build(plot)$data, original_data)
  }
})

test_that("spatial hover distinguishes unavailable pairs, zero effects and color saturation", {
  contract <- report_plot_fixture()
  rows <- contract$report$alignment$differences
  rows[marker_1 == "A" & marker_2 == "B", `:=`(status = "unavailable", value = NA_real_, p_adj = NA_real_)]
  rows[marker_1 == "B" & marker_2 == "B", value := 0]
  rows[marker_1 == "A" & marker_2 == "A", value := 5]
  plot <- plot_report_spatial(contract, "D1", "B", "differences", "TRT_vs_NC")
  built <- plotly::plotly_build(report_spatial_plotly(plot))
  missing <- Filter(function(trace) identical(trace$type, "heatmap"), built$x$data)[[1L]]
  expect_match(missing$text[is.finite(missing$z)], "Availability: unavailable", fixed = TRUE)
  expect_match(missing$text[is.finite(missing$z)], "Difference in medians: Not available", fixed = TRUE)
  expect_equal(missing$z[1, 1], 1)
  expect_equal(sum(is.finite(missing$z)), 1L)
  expect_equal(missing$customdata[[1]][1], report_feature_key("A", "B"))
  expect_equal(as.character(missing$key[[1]][1]), report_feature_key("A", "B"))
  serialized <- jsonlite::fromJSON(as.character(plotly::plotly_json(built, jsonedit = FALSE)), simplifyVector = FALSE)
  tile <- Filter(function(trace) identical(trace$type, "heatmap"), serialized$data)[[1L]]
  expect_equal(tile$key[[1]][[1]], report_feature_key("A", "B"))
  expect_equal(tile$text[[1]][[1]], missing$text[1, 1])
  expect_true(all(vapply(missing$colorscale, function(x) x[[2L]] == "#D1D1D1", logical(1))))
  points <- Filter(function(trace) identical(trace$mode, "markers") && !is.null(trace$key), built$x$data)
  zero <- Filter(function(trace) any(grepl("Difference in medians: 0<br>", trace$text, fixed = TRUE)), points)[[1L]]
  expect_identical(zero$marker$color, "rgba(255,255,255,1)")
  expect_gt(zero$marker$line$width, 0)
  expect_true(any(grepl("Difference in medians: 5<br>", unlist(lapply(points, function(trace) trace$text)), fixed = TRUE)))
  expect_equal(plot$scales$get_scales("fill")$limits, c(-2, 2))
  expect_equal(plot$data$p_adj, rows$p_adj)
})

test_that("interactive spatial matrices handle unavailable and empty scopes", {
  contract <- report_plot_fixture()
  contract$report$alignment$differences[, `:=`(status = "unavailable", value = NA_real_, p_adj = NA_real_)]
  plot <- plot_report_spatial(contract, "D1", "B", "differences", "TRT_vs_NC")
  expect_warning(widget <- report_spatial_plotly(plot), NA)
  expect_s3_class(widget, "plotly")
  empty <- plot_report_spatial(contract, "D1", "T", "differences", "TRT_vs_NC")
  expect_warning(expect_s3_class(report_spatial_plotly(empty), "plotly"), NA)
})

test_that("large unavailable matrices use one heatmap instead of per-cell polygons", {
  contract <- report_plot_fixture()
  markers <- paste0("M", seq_len(40L))
  contract$report$alignment$markers <- data.table::data.table(celltype = "B", marker = markers,
    tested_clustering = TRUE, tested_colocalization = TRUE)
  rows <- data.table::CJ(marker_1 = markers, marker_2 = markers)
  rows[, `:=`(celltype = "B", donor = "D1", panel = "TRT_vs_NC", status = "unavailable",
    value = NA_real_, p_adj = NA_real_, pct_detected = NA_real_)]
  rows[1L, `:=`(status = "supported", value = 0, p_adj = 0.02, pct_detected = 1)]
  contract$report$alignment$differences <- rows
  widget <- plotly::plotly_build(report_spatial_plotly(
    plot_report_spatial(contract, "D1", "B", "differences", "TRT_vs_NC")))
  expect_lte(length(widget$x$data), 5L)
  tiles <- Filter(function(x) identical(x$type, "heatmap"), widget$x$data)
  expect_length(tiles, 1L)
  expect_equal(sum(is.finite(tiles[[1L]]$z)), 1599L)
  expect_false(tiles[[1L]]$hoverongaps)
  expect_identical(widget$x$layout$xaxis$ticktext, markers)
  expect_identical(widget$x$layout$yaxis$ticktext, rev(markers))
})

test_that("Spatial patterns renders its default comparison on first entry", {
  value <- precomputed_app_data(precomputed_fixture())
  value$contract$report <- report_plot_fixture()$report
  current <- shiny::reactiveVal(value)
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(readout = "self_proximity", donor = "D1", celltype = "B",
      comparison = "TRT vs NC", results_view = "Spatial patterns")
    expect_identical(spatial_scope()$kind, "differences")
    expect_identical(spatial_scope()$panel, "TRT_vs_NC")
    expect_equal(nrow(spatial_report_ggplot()$data), 4L)
    expect_true(all(spatial_report_ggplot()$data$status == "supported"))
    session$setInputs(spatial_kind = "profiles")
    expect_identical(spatial_scope()$panel, "TRT")
    expect_equal(nrow(spatial_report_ggplot()$data), 4L)
  })
  html <- as.character(precomputed_ui("report"))
  expect_match(html, 'id="report-spatial_loading"', fixed = TRUE)
  expect_match(html, "Loading spatial matrix", fixed = TRUE)
  expect_match(html, 'role="status"', fixed = TRUE)
  expect_match(htmltools::renderTags(precomputed_ui("report"))$head, ".recalculating", fixed = TRUE)
})

test_that("abundance heatmaps use the frozen displayed matrix", {
  contract <- report_plot_fixture()
  grob <- plot_report_abundance_heatmap(contract, "D1", "B", "TRT vs NC")
  expect_s3_class(grob, "gtable")
  expect_identical(
    rownames(contract$report$abundance_heatmaps$D1$matrices$TRT_vs_NC$B),
    c("A", "B")
  )
})
