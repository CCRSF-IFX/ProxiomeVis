proxiome_plot_margins <- function() {
  list(l = 86, r = 86, t = 56, b = 112)
}

proxiome_plot_height <- function() {
  "430px"
}

plot_download_controls <- function(ns, output_id) {
  if (is.null(ns)) {
    ns <- identity
  }

  div(
    class = "plot-download-controls",
    downloadButton(
      ns(paste0(output_id, "_download_png")),
      "PNG",
      class = "btn btn-outline-secondary btn-sm plot-download-button"
    ),
    downloadButton(
      ns(paste0(output_id, "_download_svg")),
      "SVG",
      class = "btn btn-outline-secondary btn-sm plot-download-button"
    )
  )
}

plot_options_controls <- function(
  ns,
  width_id,
  height_id,
  width_value,
  height_value,
  min_width = 420,
  min_height = 320,
  max_value = 2600,
  step = 50,
  point_controls = NULL
) {
  if (is.null(ns)) {
    ns <- identity
  }

  sections <- list(
    div(
      class = "plot-options-section",
      div("Canvas", class = "plot-options-section-title"),
      div(
        class = "plot-options-field",
        numericInput(ns(width_id), "Width", value = width_value, min = min_width, max = max_value, step = step, width = "170px"),
        span("px", class = "plot-options-unit")
      ),
      div(
        class = "plot-options-field",
        numericInput(ns(height_id), "Height", value = height_value, min = min_height, max = max_value, step = step, width = "170px"),
        span("px", class = "plot-options-unit")
      )
    )
  )
  if (!is.null(point_controls)) {
    sections <- c(
      sections,
      list(
        div(
          class = "plot-options-section",
          div("Points", class = "plot-options-section-title"),
          point_controls
        )
      )
    )
  }

  bslib::popover(
    actionButton(
      ns(paste0(width_id, "_options")),
      "Options",
      class = "btn btn-outline-secondary btn-sm plot-options-button"
    ),
    do.call(div, c(list(class = "plot-options-popover"), sections)),
    title = "Options",
    placement = "bottom"
  )
}

plot_download_filename <- function(prefix, format) {
  prefix <- as.character(prefix)[1]
  if (is.na(prefix) || !nzchar(prefix)) {
    prefix <- "proxiomevis-plot"
  }
  prefix <- gsub("[^A-Za-z0-9._-]+", "-", prefix)
  prefix <- gsub("^-+|-+$", "", prefix)
  if (!nzchar(prefix)) {
    prefix <- "proxiomevis-plot"
  }

  paste0(prefix, ".", format)
}

plot_download_value <- function(value, default) {
  if (is.function(value)) {
    value <- value()
  }
  value <- suppressWarnings(as.numeric(value[1]))
  if (!is.finite(value) || value <= 0) {
    return(default)
  }

  value
}

plot_download_inches_from_px <- function(value, default, pixels_per_inch = 96) {
  value <- plot_download_value(value, default * pixels_per_inch)
  value / pixels_per_inch
}

plot_download_size_from_dimensions <- function(
  dimensions,
  default_width = 8,
  default_height = 5,
  pixels_per_inch = 96
) {
  if (is.function(dimensions)) {
    dimensions <- dimensions()
  }
  if (is.null(dimensions)) {
    return(list(width = default_width, height = default_height))
  }

  list(
    width = plot_download_inches_from_px(dimensions$width, default_width, pixels_per_inch = pixels_per_inch),
    height = plot_download_inches_from_px(dimensions$height, default_height, pixels_per_inch = pixels_per_inch)
  )
}

save_ggplot_download <- function(
  plot,
  file,
  format = c("png", "svg"),
  width = 8,
  height = 5,
  dpi = 300
) {
  format <- match.arg(format)
  if (!inherits(plot, "ggplot")) {
    stop("Plot download requires a ggplot object.", call. = FALSE)
  }

  if (identical(format, "svg")) {
    svg_device <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else grDevices::svg
    ggplot2::ggsave(
      filename = file,
      plot = plot,
      device = svg_device,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      limitsize = FALSE
    )
    return(invisible(file))
  }

  png_device <- if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png
  } else {
    "png"
  }
  ggplot2::ggsave(
    filename = file,
    plot = plot,
    device = png_device,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE
  )
  invisible(file)
}

register_ggplot_downloads <- function(
  output,
  output_id,
  plot,
  filename_prefix = output_id,
  width = 8,
  height = 5,
  dpi = 300
) {
  if (!is.function(plot)) {
    stop("plot must be a function or reactive expression.", call. = FALSE)
  }

  resolved_filename <- function(format) {
    prefix <- if (is.function(filename_prefix)) filename_prefix() else filename_prefix
    plot_download_filename(prefix, format)
  }
  resolved_width <- function() plot_download_value(width, 8)
  resolved_height <- function() plot_download_value(height, 5)
  resolved_dpi <- function() plot_download_value(dpi, 300)

  output[[paste0(output_id, "_download_png")]] <- downloadHandler(
    filename = function() resolved_filename("png"),
    content = function(file) {
      save_ggplot_download(
        plot(),
        file,
        format = "png",
        width = resolved_width(),
        height = resolved_height(),
        dpi = resolved_dpi()
      )
    },
    contentType = "image/png"
  )

  output[[paste0(output_id, "_download_svg")]] <- downloadHandler(
    filename = function() resolved_filename("svg"),
    content = function(file) {
      save_ggplot_download(
        plot(),
        file,
        format = "svg",
        width = resolved_width(),
        height = resolved_height(),
        dpi = resolved_dpi()
      )
    },
    contentType = "image/svg+xml"
  )

  invisible(output)
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
