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
