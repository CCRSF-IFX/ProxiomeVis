spatial_display_fixture <- function() {
  profiles <- expand.grid(condition = c("NC", "TRT"), marker_1 = c("A", "B"),
    stringsAsFactors = FALSE)
  profiles$marker_2 <- ifelse(profiles$marker_1 == "A", "A", "C")
  profiles$donor <- "D1"
  profiles$celltype <- "B"
  profiles$sample_id <- paste0("D1_", profiles$condition)
  profiles$value <- c(0, 0.2, 0.3, 0.4)
  profiles$n_total_cells <- 10
  profiles$n_scored_cells <- 5
  profiles$metric <- "population_mean_log2_ratio_missing_zero"
  differences <- data.frame(marker_1 = "A", marker_2 = "B", donor = "D1", celltype = "B",
    comparison = "TRT vs NC", effect_size = 0.7, p_adj = 0.025, group_a = "TRT", group_b = "NC",
    metric = "delta_median_log2_ratio")
  list(sample_summaries = list(self_proximity = profiles[profiles$marker_1 == profiles$marker_2, ],
    colocalization = profiles[profiles$marker_1 != profiles$marker_2, ]),
    donor_contrasts = list(self_proximity = differences[FALSE, ], colocalization = differences))
}

test_that("vectorized pair keys preserve the original JSON format for all marker text", {
  markers <- c("CD82", "CD19", "é", "Ω", "A\"B", "C\\D", "E\rF", "", NA_character_)
  pairs <- expand.grid(marker_1 = markers, marker_2 = markers, stringsAsFactors = FALSE)
  expected <- vapply(seq_len(nrow(pairs)), function(i) {
    pair <- unlist(pairs[i, ], use.names = FALSE)
    if (anyNA(pair)) return(NA_character_)
    as.character(jsonlite::toJSON(sort(pair, method = "radix")))
  }, character(1))
  expect_identical(report_feature_key(pairs$marker_1, pairs$marker_2), expected)
})

test_that("all spatial panels keep one explicit protein selection and order", {
  contract <- spatial_display_fixture()
  order <- c("C", "A", "B")
  alignment <- report_spatial_alignment(contract)
  expect_identical(report_spatial_markers(alignment, "B"), c("A", "B", "C"))
  for (kind in c("profiles", "differences")) {
    panels <- if (kind == "profiles") c("NC", "TRT") else "TRT vs NC"
    for (panel in panels) {
      plot <- plot_report_spatial(contract, "D1", "B", kind, panel, markers = order)
      expect_identical(plot$scales$get_scales("x")$limits, order)
      expect_identical(plot$scales$get_scales("y")$limits, rev(order))
      expect_equal(nrow(plot$data), 9L)
      expect_equal(plot$scales$get_scales("fill")$limits, c(-1, 1) * alignment$limits[[kind]])
      missing <- plot$data[status == "unavailable"]
      expect_true(all(is.na(missing$value)))
    }
  }
  plot <- plot_report_spatial(contract, "D1", "B", "profiles", "NC", markers = order)
  expect_identical(plot$data[marker_1 == "A" & marker_2 == "A", value], 0)
  diff <- plot_report_spatial(contract, "D1", "B", "differences", "TRT vs NC", markers = order)
  expect_equal(diff$data[marker_1 == "B" & marker_2 == "A", value], 0.7)
  expect_equal(diff$data[marker_1 == "B" & marker_2 == "A", p_adj], 0.025)
  expect_identical(contract, spatial_display_fixture())
})

test_that("comparison panels reserve identical matrix geometry despite different legends", {
  contract <- spatial_display_fixture()
  plots <- list(
    plot_report_spatial(contract, "D1", "B", "profiles", "NC"),
    plot_report_spatial(contract, "D1", "B", "profiles", "TRT"),
    plot_report_spatial(contract, "D1", "B", "differences", "TRT vs NC")
  )
  widgets <- lapply(plots, function(plot)
    plotly::plotly_build(report_spatial_plotly(plot, comparison = TRUE)))
  for (widget in widgets) {
    expect_identical(widget$x$layout$margin, widgets[[1L]]$x$layout$margin)
    expect_false(widget$x$layout$margin$autoexpand)
    for (axis in c("xaxis", "yaxis")) {
      expect_equal(widget$x$layout[[axis]]$domain, c(0, 1))
      expect_equal(widget$x$layout[[axis]]$range, widgets[[1L]]$x$layout[[axis]]$range)
      expect_false(widget$x$layout[[axis]]$automargin)
    }
    expect_identical(widget$x$layout$yaxis$scaleanchor, "x")
    expect_equal(widget$x$layout$yaxis$scaleratio, 1)
    expect_identical(widget$x$layout$dragmode, "zoom")
    expect_true(widget$x$config$scrollZoom)
  }
  single <- plotly::plotly_build(report_spatial_plotly(plots[[3L]]))
  extract_values <- function(widget) lapply(widget$x$data, function(trace)
    trace[intersect(c("x", "y", "z", "text", "key", "customdata"), names(trace))])
  expect_equal(extract_values(widgets[[3L]]), extract_values(single))
  expect_false(identical(single$x$layout$margin, widgets[[3L]]$x$layout$margin))
})

test_that("significance size legend has drawable traces matching plotted point diameters", {
  contract <- spatial_display_fixture()
  rows <- contract$donor_contrasts$colocalization[rep(1L, 4L), ]
  rows$marker_1 <- c("A", "A", "A", "B")
  rows$marker_2 <- c("A", "B", "C", "C")
  rows$p_adj <- 10 ^ -c(0, 1, 2, 4)
  contract$donor_contrasts$colocalization <- rows
  plot <- plot_report_spatial(contract, "D1", "B", "differences", "TRT vs NC")
  widget <- plotly::plotly_build(report_spatial_plotly(plot, comparison = TRUE))
  legend <- Filter(function(trace) identical(trace$legendgroup, "significance-size"), widget$x$data)
  points <- Filter(function(trace) identical(trace$mode, "markers") && !is.null(trace$key), widget$x$data)
  keys <- unlist(lapply(points, `[[`, "key"))
  sizes <- unlist(lapply(points, function(trace) rep_len(trace$marker$size, length(trace$key))))
  expect_identical(vapply(legend, `[[`, character(1), "name"), c("0", "1", "2", "4 (capped)"))
  expect_equal(vapply(legend, function(trace) trace$marker$size, numeric(1)),
    unname(sizes[match(report_feature_key(rows$marker_1, rows$marker_2), keys)]))
  for (trace in legend) {
    expect_length(trace$x, 1L)
    expect_length(trace$y, 1L)
    expect_true(is.na(trace$x[[1L]]) && is.na(trace$y[[1L]]))
    expect_true(trace$showlegend)
  }
  serialized <- jsonlite::fromJSON(as.character(plotly::plotly_json(widget, jsonedit = FALSE)),
    simplifyVector = FALSE)
  serialized_legend <- Filter(function(trace) identical(trace$legendgroup, "significance-size"), serialized$data)
  expect_identical(serialized_legend[[1L]]$x, list(NULL))
  expect_identical(serialized_legend[[1L]]$y, list(NULL))
})

test_that("profile-only bundles need no differential test", {
  contract <- spatial_display_fixture()
  contract$donor_contrasts <- NULL
  alignment <- report_spatial_alignment(contract)
  expect_identical(alignment$contrasts, character())
  expect_setequal(alignment$conditions, c("NC", "TRT"))
  plot <- plot_report_spatial(contract, "D1", "B", "profiles", "NC")
  expect_match(plot$labels$subtitle, "Descriptive profiles")
  expect_setequal(plot$scales$get_scales("x")$limits, c("A", "B", "C"))
  expect_warning(plotly::plotly_build(report_spatial_plotly(plot)), NA)
  contract$sample_summaries$self_proximity$metric <- "median_log2_ratio"
  alignment <- report_spatial_alignment(contract)
  expect_false("A" %in% alignment$profiles$marker_1)
})

test_that("new comparisons retain prepared profiles without keeping old effects", {
  contract <- spatial_display_fixture()
  alignment <- report_spatial_alignment(contract)
  alignment$profiles$value <- alignment$profiles$value + 10
  alignment$differences <- alignment$contrasts <- NULL
  alignment$limits$profiles <- 12
  contract$report$alignment <- alignment
  rebuilt <- report_spatial_alignment(contract)
  expect_identical(rebuilt$profiles, alignment$profiles)
  expect_equal(rebuilt$limits$profiles, 12)
  expect_equal(rebuilt$differences$value, 0.7)
  expect_identical(rebuilt$contrasts, "TRT vs NC")
})

test_that("prepared protein order includes profile-only proteins in differential views", {
  contract <- spatial_display_fixture()
  alignment <- report_spatial_alignment(contract)
  alignment$markers <- alignment$markers[c(3L, 1L, 2L)]
  alignment$markers$tested_clustering <- FALSE
  alignment$markers$tested_colocalization <- c(FALSE, TRUE, TRUE)
  contract$report$alignment <- alignment
  plot <- plot_report_spatial(contract, "D1", "B", "differences", "TRT vs NC")
  expect_identical(plot$scales$get_scales("x")$limits, c("C", "A", "B"))
  expect_true(all(plot$data[marker_1 == "C", status] == "unavailable"))
  expect_warning(plotly::plotly_build(report_spatial_plotly(plot)), NA)
})

test_that("canonical profiles supplement older frozen matrices without replacing stored rows", {
  contract <- spatial_display_fixture()
  alignment <- report_spatial_alignment(contract)
  alignment$profiles <- alignment$profiles[marker_1 == "A"]
  alignment$profiles[panel == "NC", `:=`(status = "unavailable", value = NA_real_)]
  alignment$markers <- alignment$markers[marker == "A"]
  alignment$limits$profiles <- 2
  contract$report$alignment <- alignment
  restored <- report_spatial_alignment(contract)
  expect_setequal(restored$profiles$marker_2, c("A", "C"))
  expect_true(is.na(restored$profiles[marker_1 == "A" & panel == "NC", value]))
  expect_identical(restored$differences, alignment$differences)
  expect_equal(restored$limits$profiles, 2)
  expect_identical(report_spatial_markers(restored, "B"), c("A", "B", "C"))
  expect_identical(as.data.frame(report_spatial_alignment(contract)$profiles), as.data.frame(restored$profiles))
  plot <- plot_report_spatial(contract, "D1", "B", "profiles", "NC", alignment = restored)
  expect_equal(plot$data[marker_1 == "C" & marker_2 == "B", value], 0.3)
})

test_that("profile aggregation weights stored population means by sample cell count", {
  contract <- spatial_display_fixture()
  rows <- contract$sample_summaries$self_proximity
  extra <- rows[rows$condition == "NC", ]
  extra$sample_id <- "D1_NC_repeat"
  extra$value <- 1
  extra$n_total_cells <- 30
  extra$n_scored_cells <- 15
  contract$sample_summaries$self_proximity <- rbind(rows, extra)
  alignment <- report_spatial_alignment(contract)
  expect_equal(alignment$profiles[panel == "NC" & marker_1 == "A", value], 0.75)
  expect_equal(alignment$profiles[panel == "NC" & marker_1 == "A", pct_detected], 0.5)
  contract$sample_summaries$self_proximity$value[1L] <- NA_real_
  expect_true(is.na(report_spatial_alignment(contract)$profiles[panel == "NC" & marker_1 == "A", value]))
})

test_that("an omitted sample cannot silently inflate imported population means", {
  contract <- spatial_display_fixture()
  abundance <- unique(contract$sample_summaries$self_proximity[c("sample_id", "donor", "celltype", "condition")])
  extra <- abundance[abundance$condition == "TRT", ]
  extra$sample_id <- "D1_TRT_repeat"
  contract$sample_summaries$abundance <- rbind(abundance, extra)
  alignment <- report_spatial_alignment(contract)
  expect_true(all(is.na(alignment$profiles[panel == "TRT", value])))
  expect_true(all(alignment$profiles[panel == "TRT", status] == "unavailable"))
  expect_true(all(is.finite(alignment$profiles[panel == "NC", value])))
  for (readout in c("self_proximity", "colocalization")) {
    rows <- contract$sample_summaries[[readout]]
    extra <- rows[rows$condition == "TRT", ]
    extra$sample_id <- "D1_TRT_repeat"
    extra$value <- extra$n_scored_cells <- 0
    extra$n_total_cells <- 30
    contract$sample_summaries[[readout]] <- rbind(rows, extra)
  }
  alignment <- report_spatial_alignment(contract)
  expect_equal(alignment$profiles[panel == "TRT" & marker_1 == "A", value], 0.05)
  expect_equal(alignment$profiles[panel == "TRT" & marker_1 == "A", pct_detected], 0.125)
})

test_that("a stored effect without a test stays visible without claiming nonsignificance", {
  contract <- spatial_display_fixture()
  contract$donor_contrasts$colocalization$p_adj <- NA_real_
  plot <- plot_report_spatial(contract, "D1", "B", "differences", "TRT vs NC")
  points <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), plot$layers)[[1L]]
  expect_identical(unique(points$data$significance), "Not tested")
  expect_equal(unique(points$data$value), 0.7)
  expect_warning(plotly::plotly_build(report_spatial_plotly(plot)), NA)
})

test_that("spatial display validates proteins and keeps absent slots for export", {
  alignment <- report_spatial_alignment(spatial_display_fixture())
  expect_error(report_spatial_rows(alignment, "D1", "B", "profiles", "NC", c("A", "A")), "unique")
  expect_error(report_spatial_rows(alignment, "D1", "B", "profiles", "NC", c("A", "")), "nonempty")
  rows <- report_spatial_rows(alignment, "D1", "B", "differences", "TRT vs NC", c("B", "C"))
  expect_equal(nrow(rows), 4L)
  expect_true(all(rows$status == "unavailable"))
  expect_true(all(is.na(rows$value)))
  expect_identical(unique(rows$marker_1), c("B", "C"))
  expect_identical(unique(rows$marker_2), c("B", "C"))
  missing <- report_spatial_rows(alignment, "D1", "B", "profiles", "NC", c("B", "D"))
  expect_true(all(missing[marker_1 == "D" | marker_2 == "D", status] == "unavailable"))
  plot <- plot_report_spatial(spatial_display_fixture(), "D2", "B", "profiles", "unstored", markers = c("D", "A"))
  expect_identical(plot$scales$get_scales("x")$limits, c("D", "A"))
  expect_true(all(plot$data$status == "unavailable"))
  expect_warning(plotly::plotly_build(report_spatial_plotly(plot)), NA)
})

test_that("canonical volcanoes use stored effects and probabilities without a frozen report", {
  contract <- spatial_display_fixture()
  plot <- plot_report_volcano(contract, "colocalization", "D1", "B", "TRT vs NC")
  expect_equal(plot$data$diff_median, 0.7)
  expect_equal(plot$data$p_adj, 0.025)
  widget <- plotly::plotly_build(report_volcano_plotly(plot))
  points <- Filter(function(trace) identical(trace$mode, "markers"), widget$x$data)
  expect_true(any(grepl("Adjusted P value: 0.025", unlist(lapply(points, `[[`, "text")), fixed = TRUE)))
  expect_warning(plotly::plotly_build(report_volcano_plotly(
    plot_report_volcano(contract, "colocalization", "D1", "T", "TRT vs NC"))), NA)
})
