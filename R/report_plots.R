report_feature_key <- function(marker_1, marker_2) {
  stopifnot(length(marker_1) == length(marker_2))
  marker_1 <- as.character(marker_1)
  marker_2 <- as.character(marker_2)
  vapply(seq_along(marker_1), function(i) {
    pair <- c(marker_1[[i]], marker_2[[i]])
    if (anyNA(pair)) return(NA_character_)
    as.character(jsonlite::toJSON(sort(pair, method = "radix")))
  }, character(1))
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

report_spatial_plotly <- function(plot, source = "A") {
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
  widget <- plotly::ggplotly(plot, tooltip = "text", source = source)
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
      colorscale = list(c(0, "#D1D1D1"), c(1, "#D1D1D1")), zmin = 0, zmax = 1,
      showscale = FALSE, showlegend = FALSE, inherit = FALSE)
    widget <- plotly::plotly_build(widget)
    i <- which(vapply(widget$x$data, function(trace) identical(trace$type, "heatmap"), logical(1)))
    # Row lists preserve heatmap click keys through Plotly's key validation.
    widget$x$data[[i]]$key <- lapply(seq_len(nrow(keys)), function(j) unname(keys[j, ]))
    widget$x$data[[i]]$customdata <- widget$x$data[[i]]$key
    widget$x$data[[i]]$`_isSimpleKey` <- FALSE
    widget$x$data[[i]]$`_isNestedKey` <- TRUE
  }
  widget |>
    plotly::event_register("plotly_click") |>
    plotly::layout(dragmode = "zoom", hovermode = "closest") |>
    plotly::config(scrollZoom = TRUE, displayModeBar = TRUE,
                   displaylogo = FALSE, doubleClick = "reset")
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

# Recreate the existing cross-donor aligned spatial matrix without changing its
# marker order, shared scales, stored effects, or stored adjusted P values.
plot_report_spatial <- function(
  contract,
  donor,
  celltype,
  kind = "differences",
  panel
) {
  donor <- report_plot_scalar(donor, "donor")
  celltype <- report_plot_scalar(celltype, "celltype")
  kind <- report_plot_scalar(kind, "kind", c("differences", "profiles"))
  panel <- report_plot_scalar(panel, "panel")
  alignment <- report_plot_bundle(contract)$alignment
  if (!is.list(alignment)) report_plot_error("The report alignment bundle is missing.")
  allowed_panels <- if (identical(kind, "profiles")) alignment$conditions else alignment$contrasts
  report_plot_scalar(panel, "panel", as.character(allowed_panels))
  .report_comparable_source$comparison_plot(
    alignment, celltype, kind, donor, panel
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
