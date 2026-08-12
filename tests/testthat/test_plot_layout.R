source("../../R/plot_layout.R")

test_that("default Plotly margins reserve space for labels and card tabs", {
  margins <- proxiome_plot_margins()

  expect_equal(margins$l, 86)
  expect_equal(margins$t, 56)
  expect_equal(margins$b, 112)
  expect_equal(margins$r, 86)
})

test_that("default Plotly height leaves room for the table below the plot", {
  expect_equal(proxiome_plot_height(), "430px")
})

test_that("plot option dimensions use defaults and valid overrides", {
  defaults <- plot_options_dimensions(width_px = NULL, height_px = NULL)
  overrides <- plot_options_dimensions(width_px = 960, height_px = 720)
  invalid <- plot_options_dimensions(width_px = -1, height_px = NA)

  expect_equal(defaults$width, 832)
  expect_equal(defaults$height, 520)
  expect_equal(overrides$width, 960)
  expect_equal(overrides$height, 720)
  expect_equal(invalid$width, defaults$width)
  expect_equal(invalid$height, defaults$height)
})

test_that("plot option IDs keep view and export controls separate", {
  ids <- plot_options_control_ids("example_width", "example_height")

  expect_equal(ids$display, "example_display")
  expect_equal(ids$view_width, "example_view_width")
  expect_equal(ids$view_height, "example_view_height")
  expect_equal(ids$export_width, "example_width")
  expect_equal(ids$export_height, "example_height")
  expect_equal(ids$options, "example_width_options")
})

test_that("plot options separate view display from export dimensions", {
  input <- list(
    example_display = "fit",
    example_view_width = 1400,
    example_view_height = 650,
    example_width = 1200,
    example_height = 800
  )

  view <- plot_options_view_dimensions(input, "example")
  export <- plot_options_export_dimensions(input, "example")

  expect_equal(view$display, "fit")
  expect_equal(view$width, 1400)
  expect_equal(view$height, 650)
  expect_equal(export$width, 1200)
  expect_equal(export$height, 800)

  default_view <- plot_options_view_dimensions(list(), "example")
  expect_equal(default_view$display, "fit")
})

test_that("Plotly display dimensions fit by default and scroll when requested", {
  base_dimensions <- list(
    width = 960,
    height = 540,
    margin = list(l = 90, r = 90, t = 56, b = 112)
  )

  fit_dimensions <- base_dimensions
  fit_dimensions$display <- "fit"
  fit_display <- plotly_display_dimensions(fit_dimensions)

  expect_null(fit_display$width)
  expect_equal(fit_display$height, base_dimensions$height)
  expect_equal(fit_display$margin, base_dimensions$margin)

  scroll_dimensions <- base_dimensions
  scroll_dimensions$display <- "scroll"
  scroll_display <- plotly_display_dimensions(scroll_dimensions)

  expect_equal(scroll_display$width, base_dimensions$width)
  expect_equal(scroll_display$height, base_dimensions$height)
  expect_equal(base_dimensions$width, 960)
})

test_that("static downloads remove plotly-only aesthetics", {
  plot <- ggplot2::ggplot(
    data.frame(x = 1:3, y = 1:3, hover = letters[1:3], id = letters[1:3]),
    ggplot2::aes(x, y, text = hover, key = id)
  ) +
    suppressWarnings(ggplot2::geom_point(ggplot2::aes(customdata = id)))

  static_plot <- prepare_ggplot_download(plot)

  expect_false("text" %in% names(static_plot$mapping))
  expect_false("key" %in% names(static_plot$mapping))
  expect_false("customdata" %in% names(static_plot$layers[[1]]$mapping))
})

test_that("Plotly frame helper applies margins and colorbar title", {
  widget <- plotly::plot_ly(
    x = c(1, 2),
    y = c(1, 2),
    type = "scatter",
    mode = "markers",
    marker = list(
      color = c(0.1, 0.8),
      colorbar = list(title = list(text = "old"))
    )
  )

  framed <- apply_proxiome_plot_frame(widget, colorbar_title = "CD3e abundance")

  expect_equal(framed$x$layout$margin, proxiome_plot_margins())
  expect_equal(framed$x$data[[1]]$marker$colorbar$title$text, "CD3e abundance")
})

test_that("differential Plotly frame helper pins domains with content-aware margins", {
  widget <- plotly::plot_ly(
    x = c(-1, 0, 1),
    y = c(2, 4, 1),
    type = "scatter",
    mode = "markers"
  )

  framed <- apply_differential_plot_frame(widget)

  expect_equal(framed$x$layout$margin, proxiome_plot_margins())
  expect_true(isTRUE(framed$x$layout$xaxis$automargin))
  expect_true(isTRUE(framed$x$layout$yaxis$automargin))
  expect_equal(framed$x$layout$xaxis$domain, c(0, 1))
  expect_equal(framed$x$layout$yaxis$domain, c(0, 1))

  dimensions <- list(
    width = 720,
    height = 540,
    margin = list(l = 90, r = 150, t = 56, b = 140)
  )
  framed_with_dimensions <- apply_differential_plot_frame(widget, dimensions = dimensions)

  expect_equal(framed_with_dimensions$x$layout$width, 720)
  expect_equal(framed_with_dimensions$x$layout$height, 540)
  expect_equal(framed_with_dimensions$x$layout$margin, dimensions$margin)
})
