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
    expect_warning(widget <- report_volcano_plotly(plot), NA)
    built <- plotly::plotly_build(widget)
    points <- Filter(function(trace) identical(trace$mode, "markers"), built$x$data)
    hover <- unlist(lapply(points, function(trace) trace$text))
    expect_true(any(grepl("Marker.*A", hover)))
    expect_true(any(grepl("Marker.*B", hover)))
    expect_true(all(grepl("Adjusted P value:", hover, fixed = TRUE)))
    expect_true(all(grepl("Target: TRT", hover, fixed = TRUE)))
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

test_that("abundance heatmaps use the frozen displayed matrix", {
  contract <- report_plot_fixture()
  grob <- plot_report_abundance_heatmap(contract, "D1", "B", "TRT vs NC")
  expect_s3_class(grob, "gtable")
  expect_identical(
    rownames(contract$report$abundance_heatmaps$D1$matrices$TRT_vs_NC$B),
    c("A", "B")
  )
})
