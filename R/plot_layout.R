proxiome_plot_margins <- function() {
  list(l = 86, r = 86, t = 56, b = 112)
}

proxiome_plot_height <- function() {
  "430px"
}

coloc_heatmap_plot_margins <- function() {
  list(l = 104, r = 184, t = 72, b = 126)
}

coloc_heatmap_panel_px <- function() {
  340
}

coloc_heatmap_widget_height_px <- function() {
  margins <- coloc_heatmap_plot_margins()
  coloc_heatmap_panel_px() + margins$t + margins$b
}

coloc_heatmap_output_height <- function() {
  paste0(coloc_heatmap_widget_height_px(), "px")
}

apply_plotly_dimensions <- function(widget, dimensions = NULL) {
  if (is.null(dimensions)) {
    return(widget)
  }

  if (!is.null(dimensions$width)) {
    widget$x$layout$width <- dimensions$width
  }
  if (!is.null(dimensions$height)) {
    widget$x$layout$height <- dimensions$height
  }

  widget
}

plot_dimensions_margin <- function(dimensions, fallback) {
  if (!is.null(dimensions) && !is.null(dimensions$margin)) {
    return(dimensions$margin)
  }

  fallback
}

apply_proxiome_plot_frame <- function(widget, colorbar_title = NULL, dimensions = NULL) {
  widget <- plotly::plotly_build(widget)
  widget <- apply_plotly_dimensions(widget, dimensions)
  widget$x$layout$margin <- plot_dimensions_margin(dimensions, proxiome_plot_margins())
  widget$x$layout$xaxis$automargin <- TRUE
  widget$x$layout$yaxis$automargin <- TRUE

  if (!is.null(colorbar_title)) {
    widget <- set_plotly_colorbar_title(widget, colorbar_title)
  }

  plotly::config(widget, displayModeBar = FALSE)
}

apply_differential_plot_frame <- function(widget, dimensions = NULL) {
  widget <- plotly::plotly_build(widget)
  widget <- apply_plotly_dimensions(widget, dimensions)
  widget$x$layout$margin <- plot_dimensions_margin(dimensions, proxiome_plot_margins())
  if (is.null(widget$x$layout$xaxis)) {
    widget$x$layout$xaxis <- list()
  }
  if (is.null(widget$x$layout$yaxis)) {
    widget$x$layout$yaxis <- list()
  }
  widget$x$layout$xaxis$automargin <- TRUE
  widget$x$layout$yaxis$automargin <- TRUE
  widget$x$layout$xaxis$domain <- c(0, 1)
  widget$x$layout$yaxis$domain <- c(0, 1)

  plotly::config(widget, displayModeBar = FALSE)
}

set_plotly_colorbar_title <- function(widget, title) {
  if (is.null(widget$x$data)) {
    widget$x$data <- list()
  }

  for (i in seq_along(widget$x$data)) {
    trace <- widget$x$data[[i]]
    if (!is.null(trace$marker$colorbar)) {
      widget$x$data[[i]]$marker$colorbar$title <- list(text = title)
    }
    if (!is.null(trace$colorbar)) {
      widget$x$data[[i]]$colorbar$title <- list(text = title)
    }
  }

  widget
}
