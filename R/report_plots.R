report_feature_key <- function(marker_1, marker_2) {
  stopifnot(length(marker_1) == length(marker_2))
  marker_1 <- as.character(marker_1)
  marker_2 <- as.character(marker_2)
  if (!length(marker_1)) return(character())
  # Encode each protein once, not once per cell/pair; retain the exact JSON click keys.
  markers <- sort(unique(c(marker_1, marker_2)), method = "radix")
  encoded <- vapply(markers, function(marker)
    as.character(jsonlite::toJSON(marker, auto_unbox = TRUE)), character(1), USE.NAMES = FALSE)
  a <- match(marker_1, markers)
  b <- match(marker_2, markers)
  keys <- paste0("[", encoded[pmin(a, b)], ",", encoded[pmax(a, b)], "]")
  keys[is.na(a) | is.na(b)] <- NA_character_
  keys
}

report_interactive_layer <- function(layer, plot_data, fields) {
  rows <- if (is.data.frame(layer$data)) layer$data else plot_data
  if (!is.data.frame(rows) || !nrow(rows)) return(layer)
  rows <- as.data.frame(rows)
  if (all(c("marker_1", "marker_2") %in% names(rows))) {
    rows$.feature_key <- report_feature_key(rows$marker_1, rows$marker_2)
  } else if ("marker" %in% names(rows)) {
    rows$.feature_key <- report_feature_key(rows$marker, rows$marker)
  } else return(layer)
  fields <- fields[names(fields) %in% names(rows)]
  rows$.report_hover <- Reduce(function(a, b) paste(a, b, sep = "<br>"),
    lapply(names(fields), function(name) {
      value <- rows[[name]]
      if (is.numeric(value)) value <- signif(value, 5)
      value <- ifelse(is.na(value), "Not available", as.character(value))
      paste0(fields[[name]], ": ", htmltools::htmlEscape(value))
    }))
  mapping <- layer$mapping
  mapping$text <- ggplot2::aes(text = .report_hover)$text
  mapping$key <- ggplot2::aes(key = .feature_key)$key
  mapping$customdata <- ggplot2::aes(customdata = .feature_key)$customdata
  # Clone the layer so the original ggplot and its SVG export stay untouched.
  ggplot2::ggproto(NULL, layer, data = rows, mapping = mapping)
}

report_volcano_plotly <- function(plot, source = "A") {
  # R plotmath labels do not reliably translate to Plotly axis titles.
  if (is.expression(plot$labels$y)) plot <- plot + ggplot2::labs(y = "-log10(adjusted P)")
  # Plotly cannot translate repel labels; point hover replaces them on screen.
  plot$layers <- Filter(function(layer) !inherits(layer$geom, "GeomTextRepel"), plot$layers)
  plot$layers <- lapply(plot$layers, function(layer) {
    if (!inherits(layer$geom, "GeomPoint")) return(layer)
    fields <- c(
      marker = "Marker", marker_1 = "Marker 1", marker_2 = "Marker 2",
      target = "Target", reference = "Reference", celltype = "Cell type",
      difference = "Abundance difference", diff_median = "Median difference",
      p = "P value", p_adj = "Adjusted P value"
    )
    report_interactive_layer(layer, plot$data, fields)
  })
  plotly::ggplotly(plot, tooltip = "text", source = source) |>
    plotly::event_register("plotly_click") |>
    plotly::layout(dragmode = "zoom", hovermode = "closest") |>
    plotly::config(scrollZoom = TRUE, displayModeBar = TRUE,
                   displaylogo = FALSE, doubleClick = "reset")
}

report_spatial_plotly <- function(plot, source = "A", comparison = FALSE) {
  fields <- c(
    marker_1 = "Marker 1", marker_2 = "Marker 2", donor = "Donor",
    celltype = "Cell type", panel = "Condition / comparison",
    value = if (is.null(plot$labels$fill)) "Stored value" else plot$labels$fill,
    p_adj = "Adjusted P value", status = "Availability", pct_detected = "Fraction detected"
  )
  # One heatmap replaces thousands of separate polygons for unavailable pairs.
  tiles <- Filter(function(layer) inherits(layer$geom, "GeomTile") &&
    is.data.frame(layer$data) && nrow(layer$data), plot$layers)
  plot$layers <- Filter(function(layer) {
    !inherits(layer$geom, "GeomTile")
  }, plot$layers)
  plot$layers <- lapply(plot$layers, function(layer) {
    if (inherits(layer$geom, "GeomLabel")) {
      return(ggplot2::ggproto(NULL, layer, geom = ggplot2::GeomText))
    }
    if (!inherits(layer$geom, c("GeomPoint", "GeomTile"))) return(layer)
    report_interactive_layer(layer, plot$data, fields)
  })
  # A converted ggplot border stays paper-sized when Plotly shrinks the axis domain.
  widget <- plotly::ggplotly(plot + ggplot2::theme(panel.border = ggplot2::element_blank()),
                            tooltip = "text", source = source)
  if (length(tiles)) {
    rows <- report_interactive_layer(tiles[[1L]], plot$data, fields)$data
    x <- plot$scales$get_scales("x")$limits
    y <- plot$scales$get_scales("y")$limits
    z <- matrix(NA_real_, nrow = length(y), ncol = length(x))
    text <- keys <- matrix(NA_character_, nrow = length(y), ncol = length(x))
    index <- cbind(match(rows$marker_2, y), match(rows$marker_1, x))
    keep <- stats::complete.cases(index)
    z[index[keep, , drop = FALSE]] <- 1
    text[index[keep, , drop = FALSE]] <- rows$.report_hover[keep]
    keys[index[keep, , drop = FALSE]] <- rows$.feature_key[keep]
    widget <- plotly::add_trace(widget, type = "heatmap", x = seq_along(x), y = seq_along(y),
      z = z, text = text, hoverinfo = "text", hoverongaps = FALSE,
      colorscale = list(c(0, "#E9EDF0"), c(1, "#E9EDF0")), zmin = 0, zmax = 1,
      showscale = FALSE, showlegend = FALSE, inherit = FALSE)
    widget <- plotly::plotly_build(widget)
    i <- which(vapply(widget$x$data, function(trace) identical(trace$type, "heatmap"), logical(1)))
    # Row lists preserve heatmap click keys through Plotly's key validation.
    widget$x$data[[i]]$key <- lapply(seq_len(nrow(keys)), function(j) unname(keys[j, ]))
    widget$x$data[[i]]$customdata <- widget$x$data[[i]]$key
    widget$x$data[[i]]$`_isSimpleKey` <- FALSE
    widget$x$data[[i]]$`_isNestedKey` <- TRUE
  }
  # Keep the marker bounds tight; shrink the panel rather than expanding empty ranges.
  if (length(plot$scales$get_scales("x")$limits)) {
    widget <- plotly::layout(widget,
      xaxis = list(constrain = "domain", constraintoward = "left", scaleanchor = FALSE),
      yaxis = list(constrain = "domain", scaleanchor = "x", scaleratio = 1))
  }
  fill_scale <- plot$scales$get_scales("fill")
  is_difference <- identical(plot$labels$fill, "Difference in medians") ||
    identical(fill_scale$name, "Difference in medians") ||
    grepl("log10(adj. P)", as.character(plot$labels$size), fixed = TRUE) ||
    any(vapply(plot$layers, function(layer) {
      is.data.frame(layer$data) && "point_size" %in% names(layer$data)
    }, logical(1)))
  size_legend <- list()
  if (isTRUE(is_difference)) {
    size_scale <- plot$scales$get_scales("size")
    size_limits <- size_scale$limits
    cap <- suppressWarnings(as.numeric(if (length(size_limits) >= 2L) size_limits[[2L]] else 4))
    if (!is.finite(cap)) cap <- 4
    breaks <- unique(pmin(c(0, 1, 2, cap), cap))
    labels <- ifelse(breaks == cap & cap > 2, paste0(cap, " (capped)"), as.character(breaks))
    # ggplotly converts the size scale's millimetres to CSS pixels.
    sizes <- size_scale$map(breaks) * 96 / 25.4
    for (i in seq_along(breaks)) {
      size_legend[[i]] <- list(type = "scatter", mode = "markers", name = labels[[i]],
        legendgroup = "significance-size", x = list(NA_real_), y = list(NA_real_), marker = list(
          size = sizes[[i]], color = "#D9D9D9",
          line = list(color = "#596873", width = 1)),
        showlegend = TRUE, hoverinfo = "skip")
    }
    widget <- plotly::layout(widget, legend = list(title = list(text =
      "Significance<br>−log10(adj. P)")))
  } else if (!is.null(plot$scales$get_scales("size"))) {
    breaks <- c(.25, .5, .75, 1)
    sizes <- plot$scales$get_scales("size")$map(breaks) * 96 / 25.4
    size_legend <- lapply(seq_along(breaks), function(i) list(type = "scatter", mode = "markers",
      name = paste0(breaks[[i]] * 100, "%"), legendgroup = "detection-size",
      x = list(NA_real_), y = list(NA_real_), marker = list(size = sizes[[i]], color = "#D9D9D9",
        line = list(color = "#596873", width = 1)), showlegend = TRUE, hoverinfo = "skip"))
    widget <- plotly::layout(widget, legend = list(title = list(text = "Cells with<br>a score")))
  }
  if (length(size_legend)) {
    widget <- plotly::plotly_build(widget)
    widget$x$data <- c(widget$x$data, size_legend)
    widget$x$layout$showlegend <- TRUE
  }
  if (length(plot$scales$get_scales("x")$limits)) {
    widget <- plotly::plotly_build(widget)
    # Reserve identical space for every panel; extra effect legends must not shrink its matrix.
    widget$x$layout$margin <- list(l = 110, r = 150, t = 65, b = 125, autoexpand = FALSE)
    for (axis in c("xaxis", "yaxis")) {
      widget$x$layout[[axis]]$domain <- c(0, 1)
      widget$x$layout[[axis]]$automargin <- FALSE
      widget$x$layout[[axis]]$tickfont <- list(size = 11, color = "#263238")
    }
    widget$x$layout$xaxis$tickangle <- -45
    widget$x$layout$yaxis$constraintoward <- "top"
    widget$x$layout$legend <- utils::modifyList(widget$x$layout$legend,
      list(x = 1.03, y = 0.47, xanchor = "left", yanchor = "top", font = list(size = 11)))
    widget$x$data <- lapply(widget$x$data, function(trace) {
      if (!any(trace$legendgroup %in% c("significance-size", "detection-size"))) trace$showlegend <- FALSE
      if (!is.null(trace$marker$colorbar)) {
        trace$marker$colorbar <- utils::modifyList(trace$marker$colorbar,
          list(x = 1.03, y = 1, yanchor = "top", len = 0.46, thickness = 12,
            xpad = 0, ypad = 0, tickfont = list(size = 11)))
        trace$marker$colorbar$title <- list(text = if (isTRUE(is_difference)) "Difference<br>in medians" else "Mean<br>log2 ratio",
          font = list(size = 11))
      }
      trace
    })
  }
  widget |>
    plotly::event_register("plotly_click") |>
    plotly::layout(dragmode = "zoom", hovermode = "closest") |>
    plotly::config(scrollZoom = TRUE, displayModeBar = "hover",
                   displaylogo = FALSE, doubleClick = "reset")
}

report_network_profiles <- function(contract, donor, celltype, comparison = NULL) {
  if (identical(celltype, PRECOMPUTED_ALL)) stop("Select one cell type for its network.", call. = FALSE)
  rows <- precomputed_filter_table(
    precomputed_contract_table(contract, "sample_summaries", "colocalization"),
    donor = donor, celltype = celltype, include_controls = TRUE)
  conditions <- sort(unique(as.character(rows$condition)), method = "radix")
  if (length(comparison) && !all(is.na(comparison))) {
    comparison <- report_plot_scalar(comparison, "comparison")
    contrasts <- precomputed_filter_table(
      precomputed_contract_table(contract, "donor_contrasts", "colocalization"),
      donor = donor, celltype = celltype, comparison = comparison, include_controls = TRUE)
    groups <- unique(contrasts[c("group_b", "group_a")])
    if (nrow(groups) != 1L) stop("No unique stored condition comparison is available.", call. = FALSE)
    conditions <- unname(unlist(groups[1L, ], use.names = FALSE))
  }
  rows <- rows[rows$condition %in% conditions & rows$marker_1 != rows$marker_2, , drop = FALSE]
  if (!nrow(rows)) stop("No stored colocalization profiles are available for this scope.", call. = FALSE)
  if (anyNA(rows$metric) || any(rows$metric != "population_mean_log2_ratio_missing_zero")) {
    stop("This network requires stored population mean log2 ratios with an explicit missing-zero policy.", call. = FALSE)
  }
  rows$key <- report_feature_key(rows$marker_1, rows$marker_2)
  if (anyDuplicated(rows[c("sample_id", "condition", "key")])) {
    stop("Multiple summaries exist for a sample and pair; prepare one profile per sample and pair.", call. = FALSE)
  }
  rows <- as.data.frame(report_population_profiles(rows, contract$sample_summaries$abundance))
  rows$key <- report_feature_key(rows$marker_1, rows$marker_2)
  rows$condition <- factor(rows$panel, levels = conditions)
  rows$panel <- NULL
  rows$metric <- "population_mean_log2_ratio_missing_zero"
  rows
}

report_network_default_markers <- function(rows) {
  rows <- rows[is.finite(rows$value) & is.finite(rows$pct_detected) & rows$pct_detected > 0.5, , drop = FALSE]
  rows <- rows[order(-abs(rows$value), rows$key), , drop = FALSE]
  head(unique(as.vector(t(as.matrix(rows[c("marker_1", "marker_2")])))), 10L)
}

plot_report_network <- function(rows, markers, min_detection = 0.5, min_mean = 0.1) {
  if (length(markers) < 2L || length(markers) > 20L || anyNA(markers) || anyDuplicated(markers) ||
      !all(markers %in% c(rows$marker_1, rows$marker_2))) {
    stop("Select 2 to 20 available network markers.", call. = FALSE)
  }
  for (threshold in list(min_detection, min_mean)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
      stop("Network thresholds must be finite, nonnegative numbers.", call. = FALSE)
    }
  }
  if (min_detection > 1) stop("Detection fraction must be between 0 and 1.", call. = FALSE)
  markers <- sort(markers, method = "radix")
  candidates <- rows[rows$marker_1 %in% markers & rows$marker_2 %in% markers &
    is.finite(rows$value) & is.finite(rows$pct_detected), , drop = FALSE]
  edges <- candidates[candidates$pct_detected > min_detection & abs(candidates$value) > min_mean, , drop = FALSE]
  graph <- igraph::graph_from_data_frame(unique(candidates[c("marker_1", "marker_2")]),
    directed = FALSE, vertices = data.frame(name = markers))
  graph <- igraph::simplify(graph)
  # A fixed initial layout and unfiltered union keep positions stable when thresholds change.
  initial <- igraph::layout_in_circle(graph)
  coords <- if (igraph::ecount(graph)) igraph::layout_with_kk(graph, coords = initial) else initial
  coords <- igraph::norm_coords(coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)
  nodes <- data.frame(marker = markers, x = coords[, 1L], y = coords[, 2L])
  nodes <- merge(nodes, data.frame(condition = factor(levels(rows$condition), levels = levels(rows$condition))), by = NULL)
  edges$x <- coords[match(edges$marker_1, markers), 1L]
  edges$y <- coords[match(edges$marker_1, markers), 2L]
  edges$xend <- coords[match(edges$marker_2, markers), 1L]
  edges$yend <- coords[match(edges$marker_2, markers), 2L]
  limit <- max(c(0.1, abs(candidates$value)))
  empty_conditions <- setdiff(levels(rows$condition), as.character(edges$condition))
  notes <- data.frame(condition = factor(empty_conditions, levels = levels(rows$condition)),
    label = rep("No pairs pass the display filters", length(empty_conditions)))
  ggplot2::ggplot(edges) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
      color = value, linewidth = abs(value))) +
    ggplot2::geom_point(ggplot2::aes(x = (x + xend) / 2, y = (y + yend) / 2,
      color = value), size = 3, alpha = 0.15, show.legend = FALSE) +
    ggplot2::geom_label(data = nodes, ggplot2::aes(x = x, y = y, label = marker), size = 4) +
    ggplot2::geom_text(data = notes, ggplot2::aes(x = 0, y = -1.25, label = label), size = 3) +
    ggplot2::facet_wrap(~condition, nrow = 1, drop = FALSE) +
    ggplot2::scale_color_gradient2(low = "#1F395F", mid = "#FFFFFF", high = "#B54868",
      midpoint = 0, limits = c(-limit, limit), name = "Mean log2 ratio") +
    ggplot2::scale_linewidth_continuous(range = c(0.4, 5), limits = c(0, limit), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-1.4, 1.4)) +
    ggplot2::scale_y_continuous(limits = c(-1.4, 1.4)) +
    ggplot2::labs(title = "Colocalization network", subtitle = paste(unique(rows$donor), unique(rows$celltype),
      paste0("detection > ", min_detection, "; |mean| > ", min_mean), sep = " · "),
      caption = "Stored population means include missing scores as zero. Display thresholds are not significance tests.") +
    ggplot2::theme_void() + ggplot2::theme(aspect.ratio = 1,
      strip.text = ggplot2::element_text(size = 12), panel.border = ggplot2::element_rect(fill = NA, color = "grey70"))
}

report_network_plotly <- function(plot, source = "A") {
  plot$theme$aspect.ratio <- NULL
  if (is.data.frame(plot$data) && !nrow(plot$data)) {
    plot$layers <- Filter(function(layer) !inherits(layer$geom, c("GeomSegment", "GeomPoint")), plot$layers)
  }
  fields <- c(marker_1 = "Marker 1", marker_2 = "Marker 2", condition = "Condition",
    value = "Mean log2 ratio (missing scores zero)", pct_detected = "Fraction detected",
    n_scored_cells = "Detected cells", n_total_cells = "Total cells")
  plot$layers <- lapply(plot$layers, function(layer) {
    if (inherits(layer$geom, "GeomLabel")) return(ggplot2::ggproto(NULL, layer, geom = ggplot2::GeomText))
    if (inherits(layer$geom, c("GeomSegment", "GeomPoint"))) return(report_interactive_layer(layer, plot$data, fields))
    layer
  })
  plotly::ggplotly(plot, tooltip = "text", source = source) |>
    plotly::event_register("plotly_click") |>
    plotly::layout(dragmode = "zoom", hovermode = "closest") |>
    plotly::config(scrollZoom = TRUE, displaylogo = FALSE, doubleClick = "reset")
}

report_plot_error <- function(...) {
  stop(paste0(...), call. = FALSE)
}

report_plot_scalar <- function(value, label, choices = NULL) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value)) {
    report_plot_error(label, " must be one nonempty value.")
  }
  if (!is.null(choices) && !value %in% choices) {
    report_plot_error(label, " must be one of: ", paste(choices, collapse = ", "), ".")
  }
  value
}

report_plot_bundle <- function(contract) {
  if (!is.list(contract) || !is.list(contract$report)) {
    report_plot_error("The precomputed contract does not contain report source data.")
  }
  contract$report
}

report_plot_nested <- function(report, section, donor, readout = NULL) {
  value <- report[[section]]
  if (!is.list(value) || is.null(value[[donor]])) {
    report_plot_error("Missing report$", section, " data for donor ", donor, ".")
  }
  value <- value[[donor]]
  if (!is.null(readout)) {
    if (!is.list(value) || is.null(value[[readout]])) {
      report_plot_error("Missing report$", section, " data for ", donor, " / ", readout, ".")
    }
    value <- value[[readout]]
  }
  value
}

report_plot_table <- function(value, label) {
  if (!is.data.frame(value)) report_plot_error(label, " must be a data frame.")
  data.table::as.data.table(value)
}

report_plot_require_columns <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing)) {
    report_plot_error(label, " is missing columns: ", paste(missing, collapse = ", "), ".")
  }
  invisible(table)
}

report_plot_contrast <- function(contract, readout, donor, celltype, comparison, source) {
  canonical <- contract$donor_contrasts[[readout]]
  if (!is.data.frame(canonical)) report_plot_error("Missing canonical donor contrasts for ", readout, ".")
  canonical <- data.table::as.data.table(canonical)
  report_plot_require_columns(
    canonical,
    c("donor", "celltype", "comparison", "analysis_id", "group_a", "group_b"),
    paste("donor_contrasts", readout)
  )
  comparison_slug <- gsub("[^[:alnum:]]+", "_", comparison)
  comparison_slug <- gsub("^_+|_+$", "", comparison_slug)
  # Resolve the comparison at donor level. Unsupported cell types can have no
  # canonical feature rows but still need the source eligibility scaffold.
  keep <- as.character(canonical$donor) == donor &
    (as.character(canonical$comparison) == comparison |
     as.character(canonical$analysis_id) == comparison |
     endsWith(as.character(canonical$analysis_id), paste0(":", comparison_slug)))
  rows <- canonical[keep]
  if (!nrow(rows)) {
    report_plot_error("No stored contrast matches ", donor, " / ", celltype, " / ", comparison, ".")
  }
  groups <- unique(rows[, .(numerator = as.character(group_a), denominator = as.character(group_b))])
  if (nrow(groups) != 1L) report_plot_error("The requested comparison has ambiguous stored groups.")
  report_plot_require_columns(source, c("contrast", "target", "reference"), "report source results")
  keep <- as.character(source$target) == groups$numerator &
    as.character(source$reference) == groups$denominator
  names <- unique(as.character(source$contrast[keep]))
  if (length(names) != 1L) {
    report_plot_error("Could not resolve one source contrast for ", groups$numerator, " vs ", groups$denominator, ".")
  }
  list(name = names[[1L]], numerator = groups$numerator, denominator = groups$denominator)
}

report_plot_settings <- function(report, donor, section) {
  settings <- report_plot_nested(report, "settings", donor)
  value <- settings[[section]]
  if (!is.list(value)) report_plot_error("Missing ", section, " settings for donor ", donor, ".")
  value
}

report_plot_auto_limit <- function(values, threshold, configured = NULL) {
  finite <- abs(as.numeric(values)[is.finite(as.numeric(values))])
  maximum <- max(c(0, finite))
  configured <- suppressWarnings(as.numeric(configured))
  if (!length(configured) || is.na(configured[[1L]])) {
    return(1.05 * max(maximum, as.numeric(threshold)))
  }
  if (length(configured) != 1L || !is.finite(configured) || configured <= 0 || maximum > configured) {
    report_plot_error("The stored volcano x limit is invalid or would clip a stored result.")
  }
  configured
}

# Recreate a donor/cell-type volcano from frozen source rows.  This function
# transforms stored values for display only; it never reruns a statistical test.
plot_report_volcano <- function(contract, readout, donor, celltype, comparison) {
  readout <- report_plot_scalar(
    readout, "readout", c("abundance", "self_proximity", "colocalization")
  )
  donor <- report_plot_scalar(donor, "donor")
  celltype <- report_plot_scalar(celltype, "celltype")
  comparison <- report_plot_scalar(comparison, "comparison")
  if (is.null(contract$report$source_results[[donor]][[readout]])) {
    rows <- as.data.frame(contract$donor_contrasts[[readout]])
    report_plot_require_columns(rows, c("donor", "celltype", "comparison", "marker_1",
      "marker_2", "effect_size", "group_a", "group_b"), "stored donor contrasts")
    rows <- rows[rows$donor == donor & rows$celltype == celltype & rows$comparison == comparison, , drop = FALSE]
    if (!"p_adj" %in% names(rows)) rows$p_adj <- rep(NA_real_, nrow(rows))
    if (!"p_value" %in% names(rows)) rows$p_value <- rep(NA_real_, nrow(rows))
    rows$p <- rows$p_value
    rows$target <- rows$group_a
    rows$reference <- rows$group_b
    rows$marker <- rows$marker_1
    effect <- if (readout == "abundance") "difference" else "diff_median"
    rows[[effect]] <- rows$effect_size
    rows <- rows[is.finite(rows$effect_size) & is.finite(rows$p_adj), , drop = FALSE]
    rows$.significance <- -log10(pmax(rows$p_adj, .Machine$double.xmin))
    plot <- ggplot2::ggplot(rows, ggplot2::aes(x = .data[[effect]], y = .significance)) +
      ggplot2::geom_point(color = "#3B6EA8", alpha = 0.7) +
      ggplot2::labs(title = paste(donor, celltype, comparison, sep = " | "),
        x = if (readout == "abundance") "Abundance difference" else "Difference in medians",
        y = "-log10(adjusted P)", subtitle = "Stored results: target minus reference") +
      ggplot2::theme_minimal(base_size = 12)
    if (!nrow(rows)) plot <- plot + ggplot2::annotate("text", x = 0, y = 0,
      label = "No finite effects with adjusted P values for this scope")
    return(plot)
  }
  report <- report_plot_bundle(contract)
  source <- report_plot_table(
    report_plot_nested(report, "source_results", donor, readout),
    paste("source results", donor, readout)
  )
  contrast <- report_plot_contrast(contract, readout, donor, celltype, comparison, source)
  contrast_name <- contrast$name
  source_keep <- as.character(source$contrast) == contrast_name &
    as.character(source$celltype) == celltype
  rows <- source[source_keep]

  plot <- if (identical(readout, "abundance")) {
    support <- report_plot_table(
      report_plot_nested(report, "eligibility", donor, readout),
      paste("eligibility", donor, readout)
    )
    report_plot_require_columns(
      source, c("marker", "p_adj", "difference", "celltype"), "abundance source results"
    )
    support_keep <- as.character(support$contrast) == contrast_name &
      as.character(support$celltype) == celltype
    support <- support[support_keep]
    if (nrow(support) != 1L) report_plot_error("Expected one abundance eligibility row for the requested scope.")
    value <- .report_abundance_source$make_volcano_plot(
      rows, support, contrast$name, contrast$numerator, contrast$denominator,
      celltype, 1L
    )
    extent <- 1.1 * max(1, abs(source$difference[is.finite(source$difference)]))
    suppressMessages(value + ggplot2::coord_cartesian(xlim = c(-extent, extent), clip = "off"))
  } else if (identical(readout, "self_proximity")) {
    cfg <- report_plot_settings(report, donor, "clustering")
    profiles <- report_plot_table(
      report_plot_nested(report, "profiles", donor, readout),
      paste("profiles", donor, readout)
    )
    profile_keep <- as.character(profiles$celltype) == celltype
    profiles <- profiles[profile_keep]
    report_plot_require_columns(
      source, c("marker_1", "marker_2", "p_adj", "diff_median", "celltype"),
      "self-proximity source results"
    )
    summary_stat <- report_plot_settings(report, donor, "spatial")$profiling$summary_stat
    metric_column <- paste0(summary_stat, "_", cfg$proximity_metric)
    x_limit <- report_plot_auto_limit(
      source$diff_median, cfg$volcano_difference_threshold, cfg$volcano_x_limit
    )
    .report_self_proximity_source$make_volcano_plot(
      rows, profiles, contrast$name, contrast$numerator, contrast$denominator,
      celltype, metric_column, cfg, x_limit, cfg$min_cells_per_group
    )
  } else {
    cfg <- report_plot_settings(report, donor, "colocalization")
    report_plot_require_columns(
      source, c("marker_1", "marker_2", "pair", "p_adj", "diff_median", "celltype"),
      "colocalization source results"
    )
    x_limit <- report_plot_auto_limit(
      source$diff_median, cfg$volcano_difference_threshold, cfg$volcano_x_limit
    )
    .report_colocalization_source$make_volcano_plot(
      rows, contrast$name, contrast$numerator, contrast$denominator,
      cfg$volcano_label_extremes, cfg$volcano_difference_threshold,
      cfg$fdr_threshold, c(-x_limit, x_limit), cfg, celltype,
      cfg$min_cells_per_group
    )
  }

  plot + ggplot2::labs(
    title = paste(donor, celltype, contrast$name, sep = " | "),
    subtitle = paste(strwrap(plot$labels$subtitle, width = 85), collapse = "\n")
  )
}

# Matrices and networks use the same cell-weighted population means and coverage checks.
report_population_profiles <- function(samples, abundance = NULL) {
  profiles <- data.table::data.table(marker_1 = character(), marker_2 = character(),
    celltype = character(), donor = character(), panel = character(), status = character(),
    value = numeric(), p_adj = numeric(), pct_detected = numeric(),
    n_total_cells = numeric(), n_scored_cells = numeric())
  if (nrow(samples)) {
    report_plot_require_columns(samples, c("sample_id", "marker_1", "marker_2", "celltype", "donor",
      "condition", "value", "metric", "n_total_cells", "n_scored_cells"), "spatial sample summaries")
    # Do not label medians or an unknown missing-score policy as population means.
    samples <- data.table::as.data.table(samples)[metric == "population_mean_log2_ratio_missing_zero"]
    if (nrow(samples)) {
      design_columns <- c("sample_id", "donor", "celltype", "condition")
      design <- unique(samples[, ..design_columns])
      if (is.data.frame(abundance) && all(design_columns %in% names(abundance))) {
        design <- unique(data.table::rbindlist(list(design, data.table::as.data.table(abundance)[, ..design_columns])))
      }
      expected <- design[, .(expected_samples = .N), by = .(celltype, donor, panel = condition)]
      samples[, key := report_feature_key(marker_1, marker_2)]
      samples <- unique(samples, by = c("sample_id", "donor", "celltype", "condition", "key"))
      profiles <- samples[, .(
        marker_1 = marker_1[[1L]], marker_2 = marker_2[[1L]],
        value = if (all(is.finite(value) & is.finite(n_total_cells) & n_total_cells > 0))
          stats::weighted.mean(value, n_total_cells) else NA_real_,
        pct_detected = if (all(is.finite(n_total_cells) & n_total_cells > 0 &
          is.finite(n_scored_cells) & n_scored_cells >= 0 & n_scored_cells <= n_total_cells))
          sum(n_scored_cells) / sum(n_total_cells) else NA_real_,
        n_total_cells = if (all(is.finite(n_total_cells) & n_total_cells > 0)) sum(n_total_cells) else NA_real_,
        n_scored_cells = if (all(is.finite(n_total_cells) & n_total_cells > 0 &
          is.finite(n_scored_cells) & n_scored_cells >= 0 & n_scored_cells <= n_total_cells))
          sum(n_scored_cells) else NA_real_,
        p_adj = NA_real_, observed_samples = .N
      ), by = .(celltype, donor, panel = condition, key)]
      profiles[expected, on = .(celltype, donor, panel), expected_samples := i.expected_samples]
      # Absent imported sample rows are not proof of zero observations in that sample.
      profiles[observed_samples != expected_samples, `:=`(value = NA_real_, pct_detected = NA_real_,
        n_total_cells = NA_real_, n_scored_cells = NA_real_)]
      profiles[, `:=`(key = NULL, observed_samples = NULL, expected_samples = NULL,
        status = ifelse(is.finite(value), "supported", "unavailable"))]
    }
  }
  profiles
}

# Canonical bundles do not need the original report scaffolding to show stored profiles.
report_spatial_alignment <- function(contract) {
  frozen <- contract$report$alignment
  empty <- data.table::data.table(marker_1 = character(), marker_2 = character(),
    celltype = character(), donor = character(), panel = character(), status = character(),
    value = numeric(), p_adj = numeric(), pct_detected = numeric())
  collect <- function(section) {
    data.table::rbindlist(lapply(c("self_proximity", "colocalization"), function(readout) {
      rows <- contract[[section]][[readout]]
      if (is.data.frame(rows)) rows else NULL
    }), fill = TRUE)
  }
  profiles <- report_population_profiles(collect("sample_summaries"), contract$sample_summaries$abundance)
  profiles[, c("n_total_cells", "n_scored_cells") := NULL]
  contrasts <- if (is.data.frame(frozen$differences)) data.table::data.table() else collect("donor_contrasts")
  differences <- data.table::copy(empty)
  if (nrow(contrasts)) {
    report_plot_require_columns(contrasts, c("marker_1", "marker_2", "celltype", "donor",
      "comparison", "effect_size", "metric"), "spatial donor contrasts")
    contrasts <- contrasts[metric == "delta_median_log2_ratio"]
    if (!"p_adj" %in% names(contrasts)) contrasts[, p_adj := NA_real_]
    differences <- contrasts[, .(marker_1, marker_2, celltype, donor, panel = comparison,
      status = ifelse(is.finite(effect_size), "supported", "unavailable"),
      value = effect_size, p_adj, pct_detected = NA_real_)]
  }
  rows <- data.table::rbindlist(list(profiles, differences), fill = TRUE)
  markers <- unique(data.table::rbindlist(list(rows[, .(celltype, marker = marker_1)],
    rows[, .(celltype, marker = marker_2)])))
  data.table::setorderv(markers, c("celltype", "marker"))
  markers[, `:=`(tested_clustering = FALSE, tested_colocalization = FALSE)]
  if (nrow(differences)) {
    tested <- differences[is.finite(p_adj)]
    self <- unique(tested[marker_1 == marker_2, .(celltype, marker = marker_1)])
    co <- tested[marker_1 != marker_2]
    co <- unique(data.table::rbindlist(list(co[, .(celltype, marker = marker_1)], co[, .(celltype, marker = marker_2)])))
    markers[self, on = .(celltype, marker), tested_clustering := TRUE]
    markers[co, on = .(celltype, marker), tested_colocalization := TRUE]
  }
  alignment <- list(cfg = list(significance_size_cap = 4, significance_threshold = 0.05),
    conditions = unique(profiles$panel), contrasts = unique(differences$panel),
    donors = unique(rows$donor), markers = markers, profiles = profiles, differences = differences,
    limits = list(profiles = max(c(0.1, abs(profiles$value[is.finite(profiles$value)]))),
      differences = max(c(0.1, abs(differences$value[is.finite(differences$value)])))),
    metric = "Mean log2 ratio (missing scores zero)")
  if (is.list(frozen)) {
    for (kind in c("profiles", "differences")) {
      if (is.data.frame(frozen[[kind]])) {
        retained <- data.table::as.data.table(frozen[[kind]])
        if (kind == "profiles" && nrow(alignment$profiles)) {
          key <- function(rows) data.table::data.table(donor = as.character(rows$donor),
            celltype = as.character(rows$celltype), panel = as.character(rows$panel),
            pair = report_feature_key(rows$marker_1, rows$marker_2))
          new_keys <- key(alignment$profiles)
          matched <- new_keys[key(retained), on = .(donor, celltype, panel, pair), which = TRUE, nomatch = 0L]
          missing <- !seq_len(nrow(new_keys)) %in% matched
          retained <- data.table::rbindlist(list(retained, alignment$profiles[missing]), fill = TRUE)
        }
        alignment[[kind]] <- retained
        alignment$limits[[kind]] <- frozen$limits[[kind]]
      }
    }
    alignment$conditions <- unique(c(frozen$conditions, as.character(alignment$profiles$panel)))
    alignment$contrasts <- unique(c(frozen$contrasts, as.character(alignment$differences$panel)))
    alignment$donors <- unique(c(as.character(alignment$profiles$donor), as.character(alignment$differences$donor)))
    alignment$markers <- unique(data.table::rbindlist(list(frozen$markers, alignment$markers), fill = TRUE),
      by = c("celltype", "marker"))
    if (is.list(frozen$cfg)) alignment$cfg <- frozen$cfg
    if (is.character(frozen$metric)) alignment$metric <- frozen$metric
  }
  alignment
}

report_spatial_markers <- function(alignment, celltype = NULL) {
  select <- function(rows, fields) {
    if (!is.data.frame(rows) || !all(fields %in% names(rows))) return(character())
    rows <- as.data.frame(rows)
    if (!is.null(celltype)) rows <- rows[rows$celltype %in% celltype, , drop = FALSE]
    as.character(unlist(rows[fields], use.names = FALSE))
  }
  markers <- unique(c(select(alignment$markers, "marker"),
    select(alignment$profiles, c("marker_1", "marker_2")),
    select(alignment$differences, c("marker_1", "marker_2"))))
  markers[!is.na(markers) & nzchar(markers)]
}

# Union result proteins across donors/contrasts, not just the visible or significant rows.
report_spatial_protein_set <- function(contract, celltype, alignment, history = list()) {
  tested <- character()
  add <- function(rows, effect = "effect_size") {
    if (!is.data.frame(rows) || !all(c("marker_1", "marker_2", effect) %in% names(rows))) return(character())
    if ("celltype" %in% names(rows)) rows <- rows[rows$celltype %in% celltype, , drop = FALSE]
    rows <- rows[is.finite(rows[[effect]]), , drop = FALSE]
    if ("is_control" %in% names(rows)) rows <- rows[!rows$is_control %in% TRUE, , drop = FALSE]
    c(as.character(rows$marker_1), as.character(rows$marker_2))
  }
  for (readout in c("self_proximity", "colocalization")) {
    for (section in c("donor_contrasts", "differential_results")) {
      tested <- c(tested, add(as.data.frame(contract[[section]][[readout]])))
    }
  }
  tested <- c(tested, add(as.data.frame(alignment$differences), "value"))
  for (entry in history) {
    if (identical(entry$request$celltype, celltype) &&
        entry$request$readout %in% c("self_proximity", "colocalization")) {
      tested <- c(tested, add(as.data.frame(entry$result$results)))
    }
  }
  tested <- sort(unique(tested[!is.na(tested) & nzchar(tested)]), method = "radix")
  if (length(tested)) return(list(markers = tested, source = "Union of spatial comparison results"))
  profiles <- as.data.frame(alignment$profiles)
  keep <- profiles$celltype %in% celltype & is.finite(profiles$value) & profiles$status %in% "supported"
  if ("pct_detected" %in% names(profiles)) keep <- keep & is.finite(profiles$pct_detected) & profiles$pct_detected > 0
  profiles <- profiles[keep, , drop = FALSE]
  markers <- sort(unique(c(as.character(profiles$marker_1), as.character(profiles$marker_2))), method = "radix")
  list(markers = markers, source = "Available condition profiles (no spatial comparison results)")
}

report_spatial_reference <- function(alignment, celltype) {
  profiles <- as.data.frame(alignment$profiles)
  profiles <- profiles[profiles$celltype %in% celltype & is.finite(profiles$value) &
    profiles$status %in% "supported", , drop = FALSE]
  choices <- unique(profiles[c("donor", "panel")])
  choices <- choices[order(choices$donor, choices$panel, method = "radix"), , drop = FALSE]
  preferred <- which(choices$donor == alignment$cfg$ordering_reference_donor &
    choices$panel == alignment$cfg$ordering_reference_condition)
  index <- if (length(preferred)) preferred[[1L]] else if (nrow(choices)) 1L else integer()
  list(choices = choices, donor = choices$donor[index], condition = choices$panel[index])
}

# Order from summaries only; never load cell scores or change plotted missing values.
report_spatial_order <- function(alignment, celltype, markers, donor, condition) {
  markers <- sort(unique(as.character(markers)), method = "radix")
  prepared <- as.data.frame(alignment$markers)
  prepared <- prepared[prepared$celltype %in% celltype, , drop = FALSE]
  reference <- report_spatial_reference(alignment, celltype)
  default_reference <- identical(donor, reference$donor) && identical(condition, reference$condition)
  if (length(alignment$cfg$ordering_reference_donor) == 1L && length(alignment$cfg$ordering_reference_condition) == 1L) {
    default_reference <- identical(donor, alignment$cfg$ordering_reference_donor) &&
      identical(condition, alignment$cfg$ordering_reference_condition)
  }
  if (nrow(prepared) && default_reference) {
    if ("marker_rank" %in% names(prepared)) prepared <- prepared[order(prepared$marker_rank), , drop = FALSE]
    stored <- unique(as.character(prepared$marker))
    # Reuse an explicit source order; the alphabetical canonical scaffold is not a prepared order.
    if ("marker_rank" %in% names(prepared) || !identical(stored, sort(stored, method = "radix"))) {
      ordered <- stored[stored %in% markers]
      return(list(markers = c(ordered, setdiff(markers, ordered)), method = "Prepared reference order"))
    }
  }
  if (length(donor) != 1L || length(condition) != 1L || is.na(donor) || is.na(condition) || length(markers) < 2L) {
    return(list(markers = markers, method = "Alphabetical (insufficient reference data)"))
  }
  grid <- as.data.frame(report_spatial_rows(alignment, donor, celltype, "profiles", condition, markers))
  matrix <- matrix(grid$value, nrow = length(markers), dimnames = list(markers, markers))
  present <- rowSums(is.finite(matrix)) > 0L
  matrix <- matrix[present, present, drop = FALSE]
  # Ward distances require a common observed feature space; do not impute unavailable scores.
  complete_columns <- colSums(!is.finite(matrix)) == 0L
  signature <- matrix[, complete_columns, drop = FALSE]
  if (nrow(signature) < 2L || ncol(signature) < 2L) {
    return(list(markers = markers, method = "Alphabetical (insufficient shared reference values)"))
  }
  tree <- stats::hclust(stats::dist(signature), method = "ward.D2")
  ordered <- rownames(signature)[tree$order]
  list(markers = c(ordered, setdiff(markers, ordered)), method = "Reference profile · Ward.D2")
}

# Preserve explicit missing entries; use the reverse pair only when that entry is absent.
report_spatial_rows <- function(alignment, donor, celltype, kind, panel, markers = NULL) {
  report_plot_scalar(kind, "kind", c("profiles", "differences"))
  if (is.null(markers)) markers <- report_spatial_markers(alignment, celltype)
  if (!is.character(markers) || anyNA(markers) || anyDuplicated(markers) ||
      any(!nzchar(markers))) {
    report_plot_error("Choose unique, nonempty spatial proteins in the required order.")
  }
  rows <- as.data.frame(alignment[[kind]])
  rows <- rows[rows$donor == donor & rows$celltype == celltype & rows$panel == panel, , drop = FALSE]
  grid <- expand.grid(marker_1 = markers, marker_2 = markers, stringsAsFactors = FALSE)
  keys <- paste(rows$marker_1, rows$marker_2, sep = "\r")
  if (anyDuplicated(keys)) report_plot_error("Multiple spatial results exist for this scope and protein pair.")
  index <- match(paste(grid$marker_1, grid$marker_2, sep = "\r"), keys)
  reverse <- match(paste(grid$marker_2, grid$marker_1, sep = "\r"), keys)
  index[is.na(index)] <- reverse[is.na(index)]
  fields <- setdiff(names(rows), c("marker_1", "marker_2", "celltype", "donor", "panel"))
  for (field in fields) grid[[field]] <- rows[[field]][index]
  grid$donor <- rep(donor, nrow(grid))
  grid$celltype <- rep(celltype, nrow(grid))
  grid$panel <- rep(panel, nrow(grid))
  for (field in c("value", "p_adj", "pct_detected")) {
    if (!field %in% names(grid)) grid[[field]] <- rep(NA_real_, nrow(grid))
  }
  if (!"status" %in% names(grid)) grid$status <- rep("supported", nrow(grid))
  grid$status[is.na(grid$status) | !is.finite(grid$value)] <- "unavailable"
  grid$value[grid$status != "supported"] <- NA_real_
  grid$p_adj[grid$status != "supported"] <- NA_real_
  data.table::as.data.table(grid)
}

# Display filters never recompute effects, adjusted P values, or color limits.
plot_report_spatial <- function(
  contract,
  donor,
  celltype,
  kind = "differences",
  panel,
  markers = NULL,
  alignment = NULL
) {
  donor <- report_plot_scalar(donor, "donor")
  celltype <- report_plot_scalar(celltype, "celltype")
  kind <- report_plot_scalar(kind, "kind", c("differences", "profiles"))
  panel <- report_plot_scalar(panel, "panel")
  if (is.null(alignment)) alignment <- report_spatial_alignment(contract)
  allowed_panels <- if (identical(kind, "profiles")) alignment$conditions else alignment$contrasts
  if (is.null(markers)) report_plot_scalar(panel, "panel", as.character(allowed_panels))
  # Explicit display selections are validated against dataset metadata by the caller.
  # Keep an unavailable grid when that valid selection has no stored results.
  panel_field <- if (kind == "profiles") "conditions" else "contrasts"
  alignment[[panel_field]] <- unique(c(alignment[[panel_field]], panel))
  alignment$donors <- unique(c(alignment$donors, donor))
  if (is.null(markers)) markers <- report_spatial_markers(alignment, celltype)
  alignment[[kind]] <- report_spatial_rows(alignment, donor, celltype, kind, panel, markers)
  .report_comparable_source$comparison_plot(
    alignment, celltype, kind, donor, panel, markers = markers
  )
}

report_abundance_matrix <- function(bundle, contrast, celltype) {
  matrices <- bundle$matrices
  if (!is.list(matrices)) report_plot_error("The frozen abundance heatmap matrices are missing.")
  matrix <- NULL
  if (!is.null(matrices[[contrast]]) && is.list(matrices[[contrast]])) {
    matrix <- matrices[[contrast]][[celltype]]
  }
  if (is.null(matrix)) matrix <- matrices[[paste(contrast, celltype, sep = "\r")]]
  if (is.null(matrix)) report_plot_error("No frozen abundance heatmap matrix matches the requested scope.")
  matrix <- as.matrix(matrix)
  storage.mode(matrix) <- "double"
  matrix
}

# Return a grid grob made with the source pheatmap palette and annotations.
# The matrix row order is already the displayed source order, so it is not
# reclustered at runtime.
plot_report_abundance_heatmap <- function(contract, donor, celltype, comparison) {
  donor <- report_plot_scalar(donor, "donor")
  celltype <- report_plot_scalar(celltype, "celltype")
  comparison <- report_plot_scalar(comparison, "comparison")
  report <- report_plot_bundle(contract)
  source <- report_plot_table(
    report_plot_nested(report, "source_results", donor, "abundance"),
    paste("source results", donor, "abundance")
  )
  contrast <- report_plot_contrast(contract, "abundance", donor, celltype, comparison, source)
  contrast_name <- contrast$name
  bundle <- report_plot_nested(report, "abundance_heatmaps", donor)
  manifest <- report_plot_table(bundle$manifest, paste("abundance heatmap manifest", donor))
  manifest_keep <- as.character(manifest$contrast) == contrast_name &
    as.character(manifest$celltype) == celltype
  manifest <- manifest[manifest_keep]
  if (nrow(manifest) != 1L) report_plot_error("Expected one frozen abundance heatmap manifest row.")
  matrix <- report_abundance_matrix(bundle, contrast$name, celltype)
  title <- paste(celltype, contrast$name, sep = " | ")
  if (!nrow(matrix)) {
    message <- if ("plot_status" %in% names(manifest)) {
      paste("No abundance heatmap:", gsub("_", " ", manifest$plot_status[[1L]]))
    } else {
      "No significant abundance hit detected"
    }
    return(ggplot2::ggplotGrob(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = message) +
        ggplot2::labs(title = title) +
        ggplot2::theme_void()
    ))
  }
  samples <- as.data.frame(bundle$samples, stringsAsFactors = FALSE)
  if (!"sample" %in% names(samples)) report_plot_error("Abundance heatmap samples need a sample column.")
  index <- match(colnames(matrix), as.character(samples$sample))
  if (anyNA(index)) report_plot_error("Abundance heatmap columns do not match the frozen sample annotation.")
  samples <- samples[index, , drop = FALSE]
  annotation <- samples[setdiff(names(samples), "sample")]
  rownames(annotation) <- samples$sample
  counts <- report_plot_table(bundle$sample_counts, paste("abundance heatmap counts", donor))
  report_plot_require_columns(counts, c("celltype", "sample", "n_cells"), "abundance heatmap counts")
  counts_keep <- as.character(counts$celltype) == celltype
  counts <- counts[counts_keep]
  n_cells <- counts[match(colnames(matrix), as.character(sample)), as.numeric(n_cells)]
  if (anyNA(n_cells)) report_plot_error("Abundance heatmap sample counts are incomplete.")
  limits <- as.numeric(manifest[1L, .(color_min, color_max)])
  if (length(limits) != 2L || any(!is.finite(limits)) || limits[[1L]] >= limits[[2L]]) {
    report_plot_error("The frozen abundance heatmap color limits are invalid.")
  }
  fdr <- manifest$fdr_threshold[[1L]]
  difference <- manifest$difference_threshold[[1L]]
  selection <- manifest$selection_samples[[1L]]
  subtitle <- paste0(
    "FDR <= ", fdr, "; |difference| > ", difference,
    "\nHits selected from: ", gsub(";", ", ", selection, fixed = TRUE)
  )
  heatmap <- pheatmap::pheatmap(
    matrix,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(c("#FFFFFF", "#B2182B"))(100),
    breaks = seq(limits[[1L]], limits[[2L]], length.out = 101),
    annotation_col = annotation,
    na_col = "grey85",
    border_color = NA,
    labels_col = paste0(samples$sample, " (n=", n_cells, ")"),
    main = paste(title, subtitle, "Mean normalized abundance (data; no scaling)", sep = "\n"),
    fontsize = 9,
    fontsize_row = 8,
    angle_col = 90,
    silent = TRUE
  )
  heatmap$gtable
}
