PRECOMPUTED_ALL <- "__all__"

PRECOMPUTED_READOUTS <- c(
  "Abundance" = "abundance",
  "Self-proximity" = "self_proximity",
  "Colocalization" = "colocalization"
)

precomputed_empty_table <- function() {
  data.frame(stringsAsFactors = FALSE)
}

precomputed_as_table <- function(x) {
  if (is.null(x)) {
    return(precomputed_empty_table())
  }
  if (!is.data.frame(x)) {
    stop("Precomputed contract tables must be data frames.", call. = FALSE)
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

precomputed_contract_table <- function(contract, section, readout) {
  section <- as.character(section)[1L]
  readout <- as.character(readout)[1L]
  if (!length(section) || is.na(section) || !nzchar(section) ||
      !length(readout) || is.na(readout) || !nzchar(readout)) {
    return(precomputed_empty_table())
  }
  tables <- contract[[section]]
  if (is.null(tables) || !is.list(tables) || is.null(tables[[readout]])) {
    return(precomputed_empty_table())
  }
  precomputed_as_table(tables[[readout]])
}

precomputed_nonempty_text <- function(x, fallback = "") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- fallback
  x
}

precomputed_feature_labels <- function(data) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    return(character())
  }

  feature_id <- if ("feature_id" %in% names(data)) {
    precomputed_nonempty_text(data$feature_id)
  } else {
    rep("", nrow(data))
  }
  marker_1 <- if ("marker_1" %in% names(data)) {
    precomputed_nonempty_text(data$marker_1)
  } else {
    rep("", nrow(data))
  }
  marker_2 <- if ("marker_2" %in% names(data)) {
    precomputed_nonempty_text(data$marker_2)
  } else {
    rep("", nrow(data))
  }

  labels <- feature_id
  self <- nzchar(marker_1) & (!nzchar(marker_2) | marker_1 == marker_2)
  pair <- nzchar(marker_1) & nzchar(marker_2) & marker_1 != marker_2
  labels[self] <- marker_1[self]
  labels[pair] <- paste(marker_1[pair], marker_2[pair], sep = " ↔ ")
  labels[!nzchar(labels)] <- "Unlabeled feature"
  labels
}

precomputed_analysis_labels <- function(data, include_donor = TRUE) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    return(character())
  }

  donor <- if ("donor" %in% names(data)) precomputed_nonempty_text(data$donor, "All donors") else "All donors"
  comparison <- if ("comparison" %in% names(data)) precomputed_nonempty_text(data$comparison, "Stored comparison") else "Stored comparison"
  celltype <- if ("celltype" %in% names(data)) precomputed_nonempty_text(data$celltype, "All cells") else "All cells"
  if (isTRUE(include_donor)) {
    paste(donor, comparison, celltype, sep = " · ")
  } else {
    paste(comparison, celltype, sep = " · ")
  }
}

precomputed_scope_value <- function(value, fallback = PRECOMPUTED_ALL) {
  value <- as.character(value)[1L]
  if (length(value) == 0L || is.na(value) || !nzchar(value)) fallback else value
}

precomputed_filter_table <- function(
  data,
  donor = PRECOMPUTED_ALL,
  comparison = PRECOMPUTED_ALL,
  celltype = PRECOMPUTED_ALL,
  include_controls = FALSE
) {
  data <- precomputed_as_table(data)
  if (nrow(data) == 0L) {
    return(data)
  }

  donor <- precomputed_scope_value(donor)
  comparison <- precomputed_scope_value(comparison)
  celltype <- precomputed_scope_value(celltype)
  keep <- rep(TRUE, nrow(data))
  if (!identical(donor, PRECOMPUTED_ALL) && "donor" %in% names(data)) {
    keep <- keep & !is.na(data$donor) & as.character(data$donor) == donor
  }
  if (!identical(comparison, PRECOMPUTED_ALL) && "comparison" %in% names(data)) {
    keep <- keep & !is.na(data$comparison) & as.character(data$comparison) == comparison
  }
  if (!identical(celltype, PRECOMPUTED_ALL) && "celltype" %in% names(data)) {
    keep <- keep & !is.na(data$celltype) & as.character(data$celltype) == celltype
  }
  if (!isTRUE(include_controls) && "is_control" %in% names(data)) {
    is_control <- as.logical(data$is_control)
    is_control[is.na(is_control)] <- FALSE
    keep <- keep & !is_control
  }
  data[keep, , drop = FALSE]
}

precomputed_available_values <- function(data, column) {
  if (!is.data.frame(data) || !column %in% names(data)) {
    return(character())
  }
  values <- sort(unique(precomputed_nonempty_text(data[[column]])))
  values[nzchar(values)]
}

precomputed_keep_selection <- function(selected, choices, fallback = NULL) {
  selected <- as.character(selected)[1]
  if (length(selected) && !is.na(selected) && selected %in% choices) {
    return(selected)
  }
  if (!is.null(fallback) && fallback %in% choices) {
    return(fallback)
  }
  choices[1]
}

precomputed_colors <- function(contract, key, values) {
  values <- unique(precomputed_nonempty_text(values, "Unknown"))
  values <- values[nzchar(values)]
  if (length(values) == 0L) {
    return(character())
  }

  key <- as.character(key)[1L]
  configured <- contract$colors
  if (is.list(configured)) {
    configured <- if (length(key) == 1L && !is.na(key) && !is.null(configured[[key]])) configured[[key]] else NULL
  }
  configured <- unlist(configured, use.names = TRUE)
  if (length(configured) > 0L && !is.null(names(configured))) {
    configured <- as.character(configured)
    matched <- configured[values]
    missing <- is.na(matched) | !nzchar(matched)
    if (any(missing)) {
      fallback <- grDevices::hcl.colors(sum(missing), palette = "Dark 3")
      matched[missing] <- fallback
    }
    names(matched) <- values
    return(matched)
  }

  palette <- grDevices::hcl.colors(max(3L, length(values)), palette = "Dark 3")[seq_along(values)]
  stats::setNames(palette, values)
}

precomputed_empty_plot <- function(message, title = NULL) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, color = "#5f6b6a", size = 4.2) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

precomputed_embedding_data <- function(metadata, embedding, color_by, abundance = NULL, marker = NULL) {
  metadata <- precomputed_as_table(metadata)
  color_by <- as.character(color_by)[1L]
  if (length(color_by) == 0L || is.na(color_by) || !nzchar(color_by)) {
    color_by <- "celltype"
  }
  if (is.null(embedding) || length(dim(embedding)) != 2L || ncol(embedding) < 2L || nrow(embedding) == 0L) {
    return(precomputed_empty_table())
  }
  if (!"component" %in% names(metadata)) {
    stop("Embedding metadata must contain component.", call. = FALSE)
  }

  component <- rownames(embedding)
  if (is.null(component) || any(!nzchar(component))) {
    if (nrow(embedding) != nrow(metadata)) {
      stop("Embedding rows require cell row names when their count differs from metadata.", call. = FALSE)
    }
    component <- as.character(metadata$component)
  }
  metadata_index <- match(component, as.character(metadata$component))
  dimensions <- colnames(embedding)
  if (is.null(dimensions) || length(dimensions) < 2L) {
    dimensions <- c("Dimension 1", "Dimension 2")
  }

  result <- data.frame(
    component = component,
    embedding_x = as.numeric(embedding[, 1L]),
    embedding_y = as.numeric(embedding[, 2L]),
    axis_x = dimensions[1L],
    axis_y = dimensions[2L],
    stringsAsFactors = FALSE
  )

  if (identical(color_by, "abundance")) {
    marker <- as.character(marker)[1]
    values <- rep(NA_real_, length(component))
    if (!is.null(abundance) && length(marker) == 1L && !is.na(marker) && nzchar(marker)) {
      if (!is.null(rownames(abundance)) && marker %in% rownames(abundance)) {
        cell_index <- match(component, colnames(abundance))
        present <- !is.na(cell_index)
        if (any(present)) {
          values[present] <- as.numeric(abundance[marker, cell_index[present], drop = TRUE])
        }
      } else if (!is.null(colnames(abundance)) && marker %in% colnames(abundance)) {
        cell_index <- match(component, rownames(abundance))
        present <- !is.na(cell_index)
        if (any(present)) {
          values[present] <- as.numeric(abundance[cell_index[present], marker, drop = TRUE])
        }
      }
    }
    result$color_value <- values
    result$color_label <- if (length(marker) && !is.na(marker) && nzchar(marker)) marker else "Abundance"
  } else {
    if (!color_by %in% names(metadata)) {
      result$color_value <- "Unknown"
    } else {
      result$color_value <- precomputed_nonempty_text(metadata[[color_by]][metadata_index], "Unknown")
    }
    result$color_label <- switch(
      color_by,
      celltype = "Cell type",
      cluster = "Cluster",
      condition = "Condition",
      donor = "Donor",
      color_by
    )
  }
  result
}

plot_precomputed_embedding <- function(
  metadata,
  embedding,
  color_by = "celltype",
  abundance = NULL,
  marker = NULL,
  colors = NULL,
  title = NULL
) {
  plot_data <- precomputed_embedding_data(metadata, embedding, color_by, abundance, marker)
  if (nrow(plot_data) == 0L) {
    return(precomputed_empty_plot("No embedding coordinates are available.", title))
  }

  numeric_color <- is.numeric(plot_data$color_value)
  plot_data$hover <- paste0(
    "Cell: ", plot_data$component,
    "<br>", plot_data$color_label, ": ",
    if (numeric_color) format(round(plot_data$color_value, 3), trim = TRUE) else plot_data$color_value
  )
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = embedding_x, y = embedding_y, color = color_value, text = hover)
  ) +
    ggplot2::geom_point(size = 0.8, alpha = 0.78, na.rm = TRUE) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      x = plot_data$axis_x[1L],
      y = plot_data$axis_y[1L],
      color = plot_data$color_label[1L]
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))

  if (numeric_color) {
    p + ggplot2::scale_color_gradientn(
      colors = c("#edf7f4", "#78aeb2", "#f0b45b", "#c7503e"),
      na.value = "#d5dadd"
    )
  } else {
    plot_data$color_value <- factor(plot_data$color_value)
    p$data <- plot_data
    requested <- colors
    if (is.null(requested) || length(requested) == 0L) {
      requested <- stats::setNames(
        grDevices::hcl.colors(max(3L, nlevels(plot_data$color_value)), "Dark 3")[seq_len(nlevels(plot_data$color_value))],
        levels(plot_data$color_value)
      )
    }
    p + ggplot2::scale_color_manual(values = requested, na.value = "#d5dadd", drop = FALSE)
  }
}

precomputed_qc_counts <- function(data) {
  contract <- data$contract
  counts <- precomputed_as_table(contract$qc$cell_counts)
  if (nrow(counts) == 0L) {
    metadata <- precomputed_as_table(data$metadata)
    needed <- c("donor", "condition", "celltype")
    if (nrow(metadata) == 0L || !all(needed %in% names(metadata))) {
      return(counts)
    }
    counts <- stats::aggregate(
      rep(1L, nrow(metadata)),
      metadata[needed],
      sum
    )
    names(counts)[ncol(counts)] <- "n_cells"
  }

  count_column <- intersect(c("n_cells", "cell_count", "cells", "n"), names(counts))[1]
  if (is.na(count_column)) {
    return(precomputed_empty_table())
  }
  counts$n_cells_display <- as.numeric(counts[[count_column]])
  for (column in c("donor", "condition", "celltype")) {
    if (!column %in% names(counts)) {
      counts[[column]] <- "All"
    }
    counts[[column]] <- precomputed_nonempty_text(counts[[column]], "All")
  }
  counts
}

plot_precomputed_qc_counts <- function(counts, colors = NULL) {
  counts <- precomputed_as_table(counts)
  if (nrow(counts) == 0L || !"n_cells_display" %in% names(counts)) {
    return(precomputed_empty_plot("No precomputed cell counts are available.", "Cell counts"))
  }
  x_column <- if ("stage" %in% names(counts)) "stage" else "condition"
  counts$plot_x <- precomputed_nonempty_text(counts[[x_column]], "All")
  counts$hover <- paste0(
    "Donor: ", counts$donor,
    "<br>", if (identical(x_column, "stage")) "Stage" else "Condition", ": ", counts$plot_x,
    "<br>Cell type: ", counts$celltype,
    "<br>Cells: ", format(counts$n_cells_display, big.mark = ",", trim = TRUE)
  )
  if (identical(x_column, "stage") && "sample_id" %in% names(counts)) {
    groups <- unique(counts$condition)
    if (is.null(colors) || length(colors) == 0L || !all(groups %in% names(colors))) {
      colors <- stats::setNames(grDevices::hcl.colors(max(3L, length(groups)), "Dark 3")[seq_along(groups)], groups)
    }
    return(
      ggplot2::ggplot(
        counts,
        ggplot2::aes(x = plot_x, y = n_cells_display, group = sample_id, color = condition, text = hover)
      ) +
        ggplot2::geom_line(linewidth = 0.8, alpha = 0.8) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::facet_wrap(~donor, scales = "free_x") +
        ggplot2::scale_color_manual(values = colors, na.value = "#8a9699", drop = FALSE) +
        ggplot2::labs(title = "Cell retention across processing stages", x = NULL, y = "Cells", color = "Condition") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          plot.title = ggplot2::element_text(face = "bold")
        )
    )
  }
  if (is.null(colors) || length(colors) == 0L) {
    levels <- unique(counts$celltype)
    colors <- stats::setNames(grDevices::hcl.colors(max(3L, length(levels)), "Dark 3")[seq_along(levels)], levels)
  }

  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = plot_x, y = n_cells_display, fill = celltype, text = hover)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~donor, scales = "free_x") +
    ggplot2::scale_fill_manual(values = colors, na.value = "#d5dadd", drop = FALSE) +
    ggplot2::labs(title = "Cells represented in the report", x = NULL, y = "Cells", fill = "Cell type") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

precomputed_result_rows_for_donor <- function(contract, readout, donor = PRECOMPUTED_ALL) {
  donor <- precomputed_scope_value(donor)
  if (identical(donor, PRECOMPUTED_ALL)) {
    return(precomputed_contract_table(contract, "differential_results", readout))
  }

  rows <- precomputed_contract_table(contract, "donor_contrasts", readout)
  if (nrow(rows) == 0L || !"donor" %in% names(rows)) {
    return(rows)
  }
  rows <- rows[as.character(rows$donor) == donor, , drop = FALSE]
  if (nrow(rows) == 0L) {
    return(rows)
  }

  # Keep any stored exploratory component-level statistics intact. These fields
  # are filled only for older descriptive contracts that do not carry them.
  defaults <- list(
    min_effect = rows$effect_size,
    max_effect = rows$effect_size,
    n_pairs = rep(1L, nrow(rows)),
    direction_consistent = rep(NA, nrow(rows)),
    p_value = rep(NA_real_, nrow(rows)),
    p_adj = rep(NA_real_, nrow(rows)),
    adjustment_family = rep(NA_character_, nrow(rows)),
    method = rep("stored donor contrast", nrow(rows)),
    replicate_unit = rep("donor", nrow(rows)),
    inference = rep("descriptive", nrow(rows))
  )
  for (field in names(defaults)) {
    if (!field %in% names(rows)) rows[[field]] <- defaults[[field]]
  }
  rows
}

precomputed_prepare_heatmap <- function(results, top_n = 20L) {
  results <- precomputed_as_table(results)
  required <- c("feature_id", "effect_size")
  if (nrow(results) == 0L || !all(required %in% names(results))) {
    return(precomputed_empty_table())
  }

  results$feature_id <- as.character(results$feature_id)
  results$feature_label <- precomputed_feature_labels(results)
  results$analysis_label <- precomputed_analysis_labels(results, include_donor = "donor" %in% names(results))
  if ("analysis_id" %in% names(results)) {
    results$analysis_key <- paste(
      precomputed_nonempty_text(results$analysis_id, "analysis"),
      if ("donor" %in% names(results)) precomputed_nonempty_text(results$donor, "All donors") else "All donors",
      precomputed_nonempty_text(results$celltype, "All cells"),
      sep = "\r"
    )
  } else {
    results$analysis_key <- results$analysis_label
  }

  # A contract key identifies one stored value. Never average duplicate rows in the viewer.
  duplicate <- duplicated(results[c("analysis_key", "feature_id")])
  if (any(duplicate)) {
    stop("Precomputed heatmap received duplicate analysis/feature rows.", call. = FALSE)
  }

  finite <- is.finite(as.numeric(results$effect_size))
  scores <- if (any(finite)) {
    tapply(
      abs(as.numeric(results$effect_size[finite])),
      results$feature_id[finite],
      max,
      na.rm = TRUE
    )
  } else {
    numeric()
  }
  top_n <- suppressWarnings(as.integer(top_n)[1])
  if (!is.finite(top_n) || top_n < 1L) {
    top_n <- 20L
  }
  if (length(scores) > 0L) {
    selected <- names(sort(scores, decreasing = TRUE))[seq_len(min(top_n, length(scores)))]
  } else {
    selected <- head(unique(results$feature_id), top_n)
  }
  selected <- unique(selected)
  if (length(selected) == 0L) {
    return(precomputed_empty_table())
  }

  analyses <- unique(results[c("analysis_key", "analysis_label")])
  feature_rows <- results[match(selected, results$feature_id), c("feature_id", "feature_label"), drop = FALSE]
  analyses$analysis_label <- make.unique(as.character(analyses$analysis_label), sep = " · ")
  feature_rows$feature_label <- make.unique(as.character(feature_rows$feature_label), sep = " · ")
  grid <- merge(
    feature_rows,
    analyses,
    by = NULL,
    sort = FALSE
  )
  values <- results[c("analysis_key", "feature_id", "effect_size")]
  for (column in intersect(c("p_adj", "p_value", "metric", "inference"), names(results))) {
    values[[column]] <- results[[column]]
  }
  plot_data <- merge(
    grid,
    values,
    by = c("analysis_key", "feature_id"),
    all.x = TRUE,
    sort = FALSE
  )
  plot_data$feature_label <- factor(plot_data$feature_label, levels = rev(feature_rows$feature_label))
  plot_data$analysis_label <- factor(plot_data$analysis_label, levels = analyses$analysis_label)
  plot_data$hover <- paste0(
    "Feature: ", plot_data$feature_label,
    "<br>Analysis: ", plot_data$analysis_label,
    "<br>Stored effect: ", ifelse(
      is.finite(plot_data$effect_size),
      format(round(plot_data$effect_size, 3), trim = TRUE),
      "Not available"
    ),
    if ("p_adj" %in% names(plot_data)) paste0(
      "<br>Stored FDR: ",
      ifelse(is.finite(plot_data$p_adj), format(plot_data$p_adj, scientific = TRUE, digits = 3), "Not available")
    ) else ""
  )
  plot_data
}

plot_precomputed_heatmap <- function(results, top_n = 20L, title = "Stored differential effects") {
  plot_data <- precomputed_prepare_heatmap(results, top_n = top_n)
  if (nrow(plot_data) == 0L) {
    return(precomputed_empty_plot("No stored effects match the selected scope.", title))
  }
  limit <- max(abs(as.numeric(plot_data$effect_size)), na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) {
    limit <- 1
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = analysis_label, y = feature_label, fill = effect_size, text = hover)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::scale_fill_gradient2(
      low = "#176d73",
      mid = "#f7f8f7",
      high = "#c7503e",
      midpoint = 0,
      limits = c(-limit, limit),
      na.value = "#d5dadd",
      name = "Stored effect"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Grey = unavailable in the stored analysis; missing values are not replaced by zero.",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

precomputed_heatmap_plotly <- function(plot, source = "A") {
  widget <- plotly::ggplotly(plot, tooltip = "text", source = source)
  rows <- plot$data
  if (is.data.frame(rows) && nrow(rows) &&
      all(c("feature_label", "analysis_label", "key") %in% names(rows))) {
    widget$x$data <- lapply(widget$x$data, function(trace) {
      if (!identical(trace$type, "heatmap")) return(trace)
      # Plotly's heatmap click handler indexes keys by [row][column], like z.
      keys <- matrix(NA_character_, nrow = length(trace$y), ncol = length(trace$x))
      positions <- cbind(match(as.integer(rows$feature_label), trace$y),
                         match(as.integer(rows$analysis_label), trace$x))
      mapped <- stats::complete.cases(positions)
      keys[positions[mapped, , drop = FALSE]] <- rows$key[mapped]
      unique_keys <- unique(as.vector(keys))
      trace$`_isSimpleKey` <- length(unique_keys) == 1L
      trace$`_isNestedKey` <- !trace$`_isSimpleKey`
      # Row lists survive Plotly's key validation; an R matrix gets flattened.
      trace$key <- if (trace$`_isSimpleKey`) unique_keys else {
        lapply(seq_len(nrow(keys)), function(i) unname(keys[i, ]))
      }
      trace
    })
  }
  widget |>
    plotly::layout(autosize = TRUE, margin = list(l = 140, r = 50, t = 100, b = 150)) |>
    plotly::event_register("plotly_click")
}

precomputed_feature_distribution_data <- function(data, readout, donor, celltype, conditions, pair) {
  if (!is.character(pair) || length(pair) != 2L || anyNA(pair) ||
      !all(pair %in% rownames(data$abundance))) stop("Select a known marker or pair.", call. = FALSE)
  if (!readout %in% c("abundance", "self_proximity", "colocalization") ||
      (readout != "colocalization" && pair[1L] != pair[2L]) ||
      (readout == "colocalization" && pair[1L] == pair[2L])) {
    stop("Select a marker or pair that matches the current readout.", call. = FALSE)
  }
  meta <- data$metadata
  rows <- meta[meta$donor == donor & meta$condition %in% conditions &
    (celltype == PRECOMPUTED_ALL | meta$celltype == celltype),
    c("component", "sample_id", "donor", "celltype", "condition"), drop = FALSE]
  rows$value <- numeric(nrow(rows))
  if (nrow(rows) && readout == "abundance") {
    rows$value <- as.numeric(data$abundance[pair[1L], rows$component])
  } else if (nrow(rows)) {
    if (is.null(data$bundle_path) || is.null(data$assets$proximity)) {
      stop("Cell-level spatial distributions require the prepared proximity sidecar.", call. = FALSE)
    }
    root <- resolve_hybrid_asset(data$bundle_path, data$assets$proximity)
    scores <- read_proximity_feature(data, root, rows$component, pair)
    rows$value <- scores$log2_ratio[match(rows$component, scores$component)]
  }
  rows <- rows[is.finite(rows$value), , drop = FALSE]
  rows$condition <- factor(rows$condition, levels = unique(conditions))
  rows
}

plot_precomputed_feature_distribution <- function(rows, readout, colors = NULL, title = "Cell-level distribution") {
  if (!nrow(rows)) return(precomputed_empty_plot("No finite cell-level values are available for this selection.", title))
  groups <- split(seq_len(nrow(rows)), interaction(rows$celltype, rows$condition, drop = TRUE))
  violin_rows <- unlist(groups[vapply(groups, function(i) length(i) >= 2L &&
    length(unique(rows$value[i])) >= 2L, logical(1))], use.names = FALSE)
  counts <- do.call(rbind, lapply(groups, function(i) {
    row <- rows[i[1L], c("celltype", "condition"), drop = FALSE]
    row$n <- length(i)
    row
  }))
  p <- ggplot2::ggplot(rows, ggplot2::aes(x = condition, y = value, fill = condition))
  if (length(violin_rows)) {
    p <- p + ggplot2::geom_violin(data = rows[violin_rows, , drop = FALSE], trim = TRUE,
      scale = "width", width = 0.85, alpha = 0.5, color = "#59696B")
    if (length(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  p <- p +
    ggplot2::geom_boxplot(width = 0.16, fill = "white", alpha = 0.9, outlier.size = 1) +
    ggplot2::geom_text(data = counts, ggplot2::aes(x = condition, y = Inf, label = paste0("n=", n)),
      inherit.aes = FALSE, vjust = 1.4, size = 3) +
    ggplot2::facet_wrap(~celltype) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(title = title, x = "Condition", fill = "Condition",
      y = if (readout == "abundance") "Stored normalized abundance" else "Observed proximity (log2 ratio)",
      subtitle = if (readout == "abundance") "Cell-level values; boxes show median and interquartile range. No new test is run."
        else "Finite observed cell scores; no zero imputation or abundance cutoff. No new test is run.",
      caption = "Violin density is shown only for groups with at least two distinct values. Cells are not independent donors.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold"))
  p
}

precomputed_effect_detail_data <- function(results, donor_contrasts) {
  results <- precomputed_as_table(results)
  donor_contrasts <- precomputed_as_table(donor_contrasts)

  aggregate <- results
  if (nrow(aggregate) > 0L) {
    aggregate$analysis_label <- precomputed_analysis_labels(aggregate, include_donor = FALSE)
    aggregate$effect_size <- as.numeric(aggregate$effect_size)
    aggregate$min_effect <- if ("min_effect" %in% names(aggregate)) as.numeric(aggregate$min_effect) else aggregate$effect_size
    aggregate$max_effect <- if ("max_effect" %in% names(aggregate)) as.numeric(aggregate$max_effect) else aggregate$effect_size
    aggregate$min_effect[!is.finite(aggregate$min_effect)] <- aggregate$effect_size[!is.finite(aggregate$min_effect)]
    aggregate$max_effect[!is.finite(aggregate$max_effect)] <- aggregate$effect_size[!is.finite(aggregate$max_effect)]
    aggregate$hover <- paste0(
      "Stored aggregate: ", format(round(aggregate$effect_size, 3), trim = TRUE),
      "<br>Range: ", format(round(aggregate$min_effect, 3), trim = TRUE),
      " to ", format(round(aggregate$max_effect, 3), trim = TRUE),
      if ("n_pairs" %in% names(aggregate)) paste0("<br>Donor pairs: ", aggregate$n_pairs) else ""
    )
  }

  donors <- donor_contrasts
  if (nrow(donors) > 0L) {
    donors$analysis_label <- precomputed_analysis_labels(donors, include_donor = FALSE)
    donors$donor <- precomputed_nonempty_text(donors$donor, "Unknown donor")
    donors$effect_size <- as.numeric(donors$effect_size)
    donors$hover <- paste0(
      "Donor: ", donors$donor,
      "<br>Stored effect: ", format(round(donors$effect_size, 3), trim = TRUE),
      if (all(c("value_a", "value_b", "group_a", "group_b") %in% names(donors))) paste0(
        "<br>", donors$group_a, ": ", format(round(donors$value_a, 3), trim = TRUE),
        "<br>", donors$group_b, ": ", format(round(donors$value_b, 3), trim = TRUE)
      ) else ""
    )
  }
  list(aggregate = aggregate, donors = donors)
}

plot_precomputed_effect_detail <- function(
  results,
  donor_contrasts,
  colors = NULL,
  title = "Stored donor effects"
) {
  detail <- precomputed_effect_detail_data(results, donor_contrasts)
  aggregate <- detail$aggregate
  donors <- detail$donors
  if (nrow(aggregate) == 0L && nrow(donors) == 0L) {
    return(precomputed_empty_plot("No stored detail is available for this feature.", title))
  }

  p <- ggplot2::ggplot(mapping = ggplot2::aes(text = hover)) +
    ggplot2::geom_hline(yintercept = 0, color = "#8a9699", linewidth = 0.5) +
    ggplot2::labs(
      title = title,
      subtitle = if (nrow(aggregate)) "Dots are stored donor contrasts; diamonds and ranges are stored cross-donor summaries."
        else "Dots are stored effects for the selected donor; no cross-donor inference is shown.",
      x = NULL,
      y = "Stored effect",
      color = "Donor"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (nrow(aggregate) > 0L) {
    p <- p +
      ggplot2::geom_linerange(
        data = aggregate,
        ggplot2::aes(x = analysis_label, ymin = min_effect, ymax = max_effect),
        linewidth = 1.1,
        color = "#263238"
      ) +
      ggplot2::geom_point(
        data = aggregate,
        ggplot2::aes(x = analysis_label, y = effect_size),
        shape = 23,
        size = 4,
        fill = "white",
        color = "#263238",
        stroke = 1
      )
  }
  if (nrow(donors) > 0L) {
    if (is.null(colors) || length(colors) == 0L) {
      donor_values <- unique(donors$donor)
      colors <- stats::setNames(grDevices::hcl.colors(max(3L, length(donor_values)), "Dark 3")[seq_along(donor_values)], donor_values)
    }
    p <- p +
      ggplot2::geom_point(
        data = donors,
        ggplot2::aes(x = analysis_label, y = effect_size, color = donor),
        position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1),
        size = 3,
        alpha = 0.9
      ) +
      ggplot2::scale_color_manual(values = colors, na.value = "#8a9699", drop = FALSE)
  }
  p
}

precomputed_condition_detail_data <- function(donor_contrasts) {
  donor_contrasts <- precomputed_as_table(donor_contrasts)
  required <- c("donor", "group_a", "group_b", "value_a", "value_b")
  if (nrow(donor_contrasts) == 0L || !all(required %in% names(donor_contrasts))) {
    return(precomputed_empty_table())
  }

  id <- if ("analysis_id" %in% names(donor_contrasts)) as.character(donor_contrasts$analysis_id) else seq_len(nrow(donor_contrasts))
  common <- data.frame(
    row_id = seq_len(nrow(donor_contrasts)),
    analysis_id = id,
    donor = precomputed_nonempty_text(donor_contrasts$donor, "Unknown donor"),
    analysis_label = precomputed_analysis_labels(donor_contrasts, include_donor = FALSE),
    stringsAsFactors = FALSE
  )
  group_b <- common
  group_b$group <- as.character(donor_contrasts$group_b)
  group_b$value <- as.numeric(donor_contrasts$value_b)
  group_b$order <- 1L
  group_a <- common
  group_a$group <- as.character(donor_contrasts$group_a)
  group_a$value <- as.numeric(donor_contrasts$value_a)
  group_a$order <- 2L
  detail <- rbind(group_b, group_a)
  detail$pair_id <- paste(detail$analysis_id, detail$donor, detail$row_id, sep = "\r")
  detail$hover <- paste0(
    "Donor: ", detail$donor,
    "<br>Condition: ", detail$group,
    "<br>Stored value: ", ifelse(
      is.finite(detail$value),
      format(round(detail$value, 3), trim = TRUE),
      "Not available"
    )
  )
  detail
}

plot_precomputed_condition_detail <- function(
  donor_contrasts,
  colors = NULL,
  title = "Stored paired condition values"
) {
  detail <- precomputed_condition_detail_data(donor_contrasts)
  if (nrow(detail) == 0L) {
    return(precomputed_empty_plot("No paired condition values are available for this feature.", title))
  }
  if (!any(is.finite(detail$value))) {
    return(precomputed_empty_plot(
      "Condition-level source values were not stored for this feature; use the stored-effect view.",
      title
    ))
  }
  donor_values <- unique(detail$donor)
  if (is.null(colors) || length(colors) == 0L) {
    colors <- stats::setNames(grDevices::hcl.colors(max(3L, length(donor_values)), "Dark 3")[seq_along(donor_values)], donor_values)
  }

  ggplot2::ggplot(
    detail,
    ggplot2::aes(x = group, y = value, group = pair_id, color = donor, text = hover)
  ) +
    ggplot2::geom_line(alpha = 0.55, linewidth = 0.7, na.rm = FALSE) +
    ggplot2::geom_point(size = 3, alpha = 0.9, na.rm = FALSE) +
    ggplot2::facet_wrap(~analysis_label, scales = "free_x") +
    ggplot2::scale_color_manual(values = colors, na.value = "#8a9699", drop = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = "Values and donor pairings are read directly from the stored contrast table.",
      x = NULL,
      y = "Stored value",
      color = "Donor"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

precomputed_write_svg <- function(plot, path, width = 11, height = 7) {
  grDevices::svg(path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}

precomputed_draw_plot_object <- function(plot) {
  if (inherits(plot, "ggplot")) {
    print(plot)
  } else if (inherits(plot, "grob")) {
    grid::grid.newpage()
    grid::grid.draw(plot)
  } else {
    print(plot)
  }
  invisible(plot)
}

precomputed_write_plot_object_svg <- function(plot, path, width = 11, height = 7) {
  grDevices::svg(path, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  precomputed_draw_plot_object(plot)
  invisible(path)
}

precomputed_safe_filename <- function(...) {
  value <- paste(..., sep = "-")
  value <- gsub("[^A-Za-z0-9._-]+", "-", value)
  value <- gsub("-+", "-", value)
  value <- sub("^-", "", sub("-$", "", value))
  if (!nzchar(value)) "proxiomevis" else value
}

precomputed_report_values <- function(x, column = NULL) {
  if (is.null(x)) return(character())
  if (is.data.frame(x)) {
    if (!is.null(column) && column %in% names(x)) x <- x[[column]]
    else {
      candidate <- names(x)[vapply(x, function(value) is.character(value) || is.factor(value), logical(1))][1L]
      if (is.na(candidate)) return(character())
      x <- x[[candidate]]
    }
  } else if (is.list(x)) {
    if (!is.null(column) && !is.null(x[[column]])) {
      return(precomputed_report_values(x[[column]]))
    }
    simple <- vapply(x, function(value) length(value) == 1L && is.atomic(value), logical(1))
    if (length(simple) && all(simple)) x <- unlist(x, use.names = FALSE)
    else if (length(names(x)) && all(nzchar(names(x)))) x <- names(x)
    else return(character())
  }
  values <- unique(precomputed_nonempty_text(x))
  values[nzchar(values)]
}

precomputed_report_celltypes <- function(contract, metadata, donor, readout, donor_rows) {
  donor <- as.character(donor)[1L]
  readout <- as.character(readout)[1L]
  eligibility <- NULL
  valid_scope <- length(donor) == 1L && !is.na(donor) && nzchar(donor) &&
    length(readout) == 1L && !is.na(readout) && nzchar(readout)
  if (valid_scope && is.list(contract$report$eligibility) && !is.null(contract$report$eligibility[[donor]]) &&
      is.list(contract$report$eligibility[[donor]])) {
    eligibility <- contract$report$eligibility[[donor]][[readout]]
  }
  celltypes <- precomputed_report_values(eligibility, "celltype")
  if (length(celltypes) == 0L && valid_scope && is.data.frame(metadata) && all(c("donor", "celltype") %in% names(metadata))) {
    celltypes <- sort(unique(precomputed_nonempty_text(
      metadata$celltype[as.character(metadata$donor) == donor]
    )))
    celltypes <- celltypes[nzchar(celltypes)]
  }
  if (length(celltypes) == 0L) {
    celltypes <- precomputed_available_values(donor_rows, "celltype")
  }
  celltypes
}

precomputed_report_plot <- function(function_name, ..., title) {
  if (!exists(function_name, mode = "function", inherits = TRUE)) {
    return(precomputed_empty_plot("The stored report adapter is not available in this build.", title))
  }
  plot_function <- get(function_name, mode = "function", inherits = TRUE)
  tryCatch(
    plot_function(...),
    error = function(error) precomputed_empty_plot(conditionMessage(error), title)
  )
}

precomputed_download_group <- function(...) {
  shiny::div(class = "d-flex flex-wrap gap-2 mb-2", ...)
}

precomputed_metric_card <- function(label, output_id) {
  bslib::card(
    bslib::card_body(
      shiny::tags$div(class = "text-muted small text-uppercase", label),
      shiny::tags$div(class = "fs-3 fw-semibold", shiny::textOutput(output_id, inline = TRUE))
    )
  )
}


precomputed_explorer_scope <- function(data, readout = NULL, donor = NULL, celltype = NULL, comparison = NULL) {
  readouts <- unname(PRECOMPUTED_READOUTS)
  readout <- precomputed_keep_selection(readout, readouts)
  rows <- precomputed_contract_table(data$contract, "donor_contrasts", readout)
  donors <- sort(unique(as.character(data$metadata$donor)))
  donor <- precomputed_keep_selection(donor, donors)
  donor_rows <- precomputed_filter_table(rows, donor = donor, include_controls = TRUE)
  celltypes <- precomputed_report_celltypes(data$contract, data$metadata, donor, readout, donor_rows)
  if (!length(celltypes)) celltypes <- sort(unique(as.character(data$metadata$celltype)))
  celltype <- precomputed_keep_selection(celltype, c(celltypes, PRECOMPUTED_ALL))
  comparisons <- precomputed_available_values(donor_rows, "comparison")
  comparison <- precomputed_keep_selection(comparison, comparisons)
  selected <- if (length(comparison) && !is.na(comparison)) {
    precomputed_filter_table(donor_rows, celltype = celltype, comparison = comparison, include_controls = TRUE)
  } else donor_rows[FALSE, , drop = FALSE]
  list(readout = readout, donor = donor, donors = donors, celltype = celltype, celltypes = celltypes,
       comparison = comparison, comparisons = comparisons, rows = selected)
}

precomputed_display_rows <- function(rows, filter = FALSE, cutoff = 0.05) {
  if (!nrow(rows)) return(rows)
  rows$feature <- precomputed_feature_labels(rows)
  if ("p_adj" %in% names(rows) && any(is.finite(rows$p_adj))) {
    cutoff <- suppressWarnings(as.numeric(cutoff)[1L])
    if (!length(cutoff) || !is.finite(cutoff) || cutoff < 0 || cutoff > 1) cutoff <- 0.05
    rows$stored_fdr_status <- ifelse(!is.finite(rows$p_adj), "Not available",
      ifelse(rows$p_adj <= cutoff, paste0("FDR ≤ ", cutoff), paste0("FDR > ", cutoff)))
    if (isTRUE(filter)) rows <- rows[is.finite(rows$p_adj) & rows$p_adj <= cutoff, , drop = FALSE]
  }
  rows
}

precomputed_searched_rows <- function(rows, indices, search = "", column_search = character()) {
  # Shiny decodes DataTables' zero-match [] as NULL, just like an unbound table.
  if (is.null(indices)) {
    searching <- any(nzchar(c(search, unlist(column_search))), na.rm = TRUE)
    return(if (searching) rows[FALSE, , drop = FALSE] else rows)
  }
  indices <- as.integer(indices)
  indices <- indices[!is.na(indices) & indices > 0L & indices <= nrow(rows)]
  rows[unique(indices), , drop = FALSE]
}

precomputed_server <- function(id, data, job_manager = new_hybrid_job_manager(), code_dir = "R") {
  shiny::moduleServer(id, function(input, output, session) {
    current_data <- shiny::reactive({
      value <- data()
      shiny::req(value, value$contract, value$metadata, value$embeddings)
      value
    })
    hybrid_analysis_server("analysis", current_data, job_manager, code_dir,
      active = shiny::reactive(identical(input$report_page, "Custom analysis")))
    graph_selection <- shiny::reactiveVal(NULL)
    cell_graph_server("graph", current_data,
      active = shiny::reactive(identical(input$report_page, "Cell 3D")), selection = graph_selection)

    scope <- shiny::reactive(precomputed_explorer_scope(current_data(),
      input$readout, input$donor, input$celltype, input$comparison))
    output$donor_control <- shiny::renderUI({
      donors <- sort(unique(as.character(current_data()$metadata$donor)))
      if (length(donors) == 1L) return(shiny::tags$p(class = "pt-2", shiny::tags$strong("Donor: "), donors))
      shiny::selectInput(session$ns("donor"), "Donor", choices = donors,
        selected = precomputed_keep_selection(shiny::isolate(input$donor), donors))
    })

    readout_tables <- shiny::reactive({
      value <- current_data()
      readout <- input$readout
      if (is.null(readout) || !readout %in% unname(PRECOMPUTED_READOUTS)) {
        readout <- unname(PRECOMPUTED_READOUTS)[1L]
      }
      list(
        readout = readout,
        results = precomputed_contract_table(value$contract, "differential_results", readout),
        contrasts = precomputed_contract_table(value$contract, "donor_contrasts", readout),
        samples = precomputed_contract_table(value$contract, "sample_summaries", readout)
      )
    })

    shiny::observeEvent(current_data(), {
      value <- current_data()
      embedding_names <- names(value$embeddings)
      defaults <- value$contract$defaults
      umap_selected <- precomputed_keep_selection(defaults$reduction, embedding_names)
      pca_selected <- precomputed_keep_selection(defaults$pca_reduction, embedding_names, umap_selected)
      shiny::updateSelectInput(session, "umap_reduction", choices = embedding_names, selected = umap_selected)
      shiny::updateSelectInput(session, "pca_reduction", choices = embedding_names, selected = pca_selected)

      color_choices <- c(
        "Cell type" = "celltype",
        "Cluster" = "cluster",
        "Condition" = "condition",
        "Donor" = "donor"
      )
      color_choices <- color_choices[unname(color_choices) %in% names(value$metadata)]
      if (!is.null(value$abundance)) {
        color_choices <- c(color_choices, "Marker abundance" = "abundance")
      }
      shiny::updateSelectInput(
        session,
        "embedding_color",
        choices = color_choices,
        selected = precomputed_keep_selection(shiny::isolate(input$embedding_color), unname(color_choices), "celltype")
      )

      markers <- rownames(value$abundance)
      if (is.null(markers) || length(markers) == 0L) {
        markers <- colnames(value$abundance)
      }
      markers <- as.character(markers)
      shiny::updateSelectizeInput(
        session,
        "embedding_marker",
        choices = markers,
        selected = precomputed_keep_selection(shiny::isolate(input$embedding_marker), markers),
        server = TRUE
      )

      available_readouts <- unname(PRECOMPUTED_READOUTS)[vapply(
        unname(PRECOMPUTED_READOUTS),
        function(readout) nrow(precomputed_contract_table(value$contract, "differential_results", readout)) > 0L ||
          nrow(precomputed_contract_table(value$contract, "donor_contrasts", readout)) > 0L,
        logical(1)
      )]
      choices <- PRECOMPUTED_READOUTS[PRECOMPUTED_READOUTS %in% available_readouts]
      shiny::updateSelectInput(
        session,
        "readout",
        choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$readout), unname(choices))
      )
    }, ignoreInit = FALSE)

    shiny::observeEvent(list(current_data(), input$readout, input$donor), {
      selection <- scope()
      shiny::updateSelectInput(session, "celltype",
        choices = c(stats::setNames(selection$celltypes, selection$celltypes), "All cell types (summary)" = PRECOMPUTED_ALL),
        selected = selection$celltype)
      shiny::updateSelectInput(session, "comparison", choices = selection$comparisons, selected = selection$comparison)
    })

    selected_results <- shiny::reactive(scope()$rows)
    selected_contrasts <- shiny::reactive(selected_results())

    selected_samples <- shiny::reactive({
      selected <- scope()
      samples <- precomputed_filter_table(readout_tables()$samples, donor = selected$donor,
        celltype = selected$celltype, include_controls = TRUE)
      rows <- precomputed_filter_table(readout_tables()$contrasts, donor = selected$donor,
        comparison = selected$comparison, include_controls = TRUE)
      conditions <- unique(c(as.character(rows$group_a), as.character(rows$group_b)))
      samples[samples$condition %in% conditions, , drop = FALSE]
    })

    selected_key <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$readout, {
      pair <- tryCatch(jsonlite::fromJSON(selected_key()), error = function(e) character())
      if (length(pair) == 2L &&
          xor(scope()$readout == "colocalization", pair[1L] != pair[2L])) {
        keys <- feature_rows()$key
        selected_key(if (length(keys)) keys[1L] else NULL)
      }
    })
    feature_rows <- shiny::reactive({
      rows <- scope()$rows
      if (!nrow(rows)) {
        rows$key <- character()
        return(rows)
      }
      rows$key <- report_feature_key(rows$marker_1, rows$marker_2)
      rows[!duplicated(rows$key), , drop = FALSE]
    })
    shiny::observe({
      rows <- feature_rows()
      keys <- rows$key
      chosen <- selected_key()
      if (!length(chosen) && length(keys)) {
        selected_key(keys[1L])
        chosen <- keys[1L]
      }
      choices <- stats::setNames(keys, precomputed_feature_labels(rows))
      if (length(chosen) && !chosen %in% keys) {
        pair <- tryCatch(jsonlite::fromJSON(chosen), error = function(e) character())
        if (length(pair) == 2L) choices <- c(choices, stats::setNames(chosen,
          paste0(paste(unique(pair), collapse = " ↔ "), " (no stored test in this scope)")))
      }
      shiny::updateSelectizeInput(session, "detail_feature", choices = choices,
        selected = chosen, server = TRUE)
    })
    shiny::observeEvent(input$detail_feature, {
      if (input$detail_feature %in% c(feature_rows()$key, selected_key())) selected_key(input$detail_feature)
    })
    shiny::observeEvent(current_data(), { selected_key(NULL) }, ignoreInit = TRUE)
    detail_results <- shiny::reactive(precomputed_empty_table())
    detail_contrasts <- shiny::reactive({
      rows <- scope()$rows
      if (!nrow(rows) || !length(selected_key())) return(rows[FALSE, , drop = FALSE])
      rows[report_feature_key(rows$marker_1, rows$marker_2) == selected_key(), , drop = FALSE]
    })

    select_feature <- function(key, from_spatial = FALSE) {
      pair <- tryCatch(jsonlite::fromJSON(as.character(key)[1L]), error = function(e) character())
      if (!is.character(pair) || length(pair) != 2L || anyNA(pair) ||
          !all(pair %in% rownames(current_data()$abundance))) return(invisible(NULL))
      selected_key(report_feature_key(pair[1L], pair[2L]))
      if (from_spatial) shiny::updateSelectInput(session, "readout",
        selected = if (pair[1L] == pair[2L]) "self_proximity" else "colocalization")
    }
    for (view in c("summary", "volcano", "spatial")) local({
      plot_view <- view
      shiny::observeEvent({
        shiny::req(session$rootScope()$input[[paste0("plotly_click-", session$ns(plot_view))]])
        plotly::event_data("plotly_click", source = session$ns(plot_view), priority = "event")
      }, {
        event <- plotly::event_data("plotly_click", source = session$ns(plot_view), priority = "event")
        key <- if (!is.null(event$key)) event$key else event$customdata
        if (length(key)) select_feature(key[1L], from_spatial = plot_view == "spatial")
      })
    })
    shiny::observeEvent(input$results_table_rows_selected, {
      rows <- display_results()
      i <- input$results_table_rows_selected
      if (length(i) == 1L && i > 0 && i <= nrow(rows))
        select_feature(report_feature_key(rows$marker_1[i], rows$marker_2[i]))
    })
    shiny::observeEvent(input$report_volcano_table_rows_selected, {
      rows <- report_volcano_display()
      i <- input$report_volcano_table_rows_selected
      if (length(i) == 1L && i > 0 && i <= nrow(rows))
        select_feature(report_feature_key(rows$marker_1[i], rows$marker_2[i]))
    })
    output$feature_context <- shiny::renderUI({
      pair <- tryCatch(jsonlite::fromJSON(selected_key()), error = function(e) character())
      if (length(pair) != 2L) return(shiny::tags$p("Click a plot or table row to inspect a feature."))
      rows <- detail_contrasts()
      fields <- intersect(c("celltype", "effect_size", "metric", "p_adj", "n_cells_a", "n_cells_b"), names(rows))
      shiny::tagList(
        shiny::tags$strong(paste(unique(pair), collapse = " ↔ ")),
        if (!nrow(rows)) shiny::tags$p("No stored differential test for this feature in the selected scope."),
        if (nrow(rows)) shiny::tags$ul(lapply(seq_len(nrow(rows)), function(i)
          shiny::tags$li(paste(paste(fields, unlist(rows[i, fields, drop = FALSE]), sep = ": "), collapse = " · ")))),
        shiny::tags$p(class = "small", "Cell-level exploration is separate from these frozen statistics. ",
          "Choose a condition and one highlight marker in Cell 3D, then click Display."))
    })
    shiny::observeEvent(input$explore_3d, {
      shared <- scope()
      pair <- tryCatch(jsonlite::fromJSON(selected_key()), error = function(e) character())
      if (shared$celltype == PRECOMPUTED_ALL || length(pair) != 2L) {
        shiny::showNotification("Choose one cell type and a feature before exploring cells.", type = "warning")
        return()
      }
      rows <- precomputed_filter_table(readout_tables()$contrasts, donor = shared$donor,
        comparison = shared$comparison, include_controls = TRUE)
      graph_selection(list(donor = shared$donor, celltype = shared$celltype,
        marker_1 = pair[1L], marker_2 = pair[2L],
        conditions = unique(c(as.character(rows$group_a), as.character(rows$group_b))),
        token = input$explore_3d))
      bslib::nav_select("report_page", "Cell 3D", session = session)
    })

    output$dataset_heading <- shiny::renderUI({
      contract <- current_data()$contract
      shiny::tagList(
        shiny::tags$h3(contract$title),
        shiny::tags$p(
          class = "text-muted",
          paste0("Dataset ", contract$dataset_id, " · schema version ", contract$schema_version,
                 " · report-only views; PXL files are not required.")
        )
      )
    })

    output$metric_cells <- shiny::renderText(format(nrow(current_data()$metadata), big.mark = ","))
    output$metric_samples <- shiny::renderText(format(length(unique(current_data()$metadata$sample_id)), big.mark = ","))
    output$metric_donors <- shiny::renderText(format(length(unique(current_data()$metadata$donor)), big.mark = ","))
    output$metric_celltypes <- shiny::renderText(format(length(unique(current_data()$metadata$celltype)), big.mark = ","))

    qc_counts <- shiny::reactive(precomputed_qc_counts(current_data()))
    qc_ggplot <- shiny::reactive({
      value <- current_data()
      counts <- qc_counts()
      color_key <- if ("stage" %in% names(counts) && "sample_id" %in% names(counts)) "condition" else "celltype"
      color_values <- if (color_key %in% names(counts)) counts[[color_key]] else character()
      plot_precomputed_qc_counts(counts, precomputed_colors(value$contract, color_key, color_values))
    })
    output$qc_plot <- plotly::renderPlotly({
      plotly::ggplotly(qc_ggplot(), tooltip = "text") |>
        plotly::layout(autosize = TRUE, margin = list(l = 70, r = 30, t = 80, b = 90))
    })
    output$qc_table <- DT::renderDT({
      shiny::validate(shiny::need(nrow(qc_counts()) > 0L, "No stored QC counts are available."))
      DT::datatable(qc_counts(), rownames = FALSE, filter = "top", options = list(pageLength = 12, scrollX = TRUE))
    })

    provenance_json <- shiny::reactive({
      jsonlite::toJSON(
        current_data()$contract$provenance,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null"
      )
    })
    output$provenance_text <- shiny::renderText(provenance_json())
    output$download_provenance <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(current_data()$contract$dataset_id, "provenance"), ".json"),
      content = function(file) writeLines(provenance_json(), file, useBytes = TRUE)
    )

    embedding_ggplot <- function(reduction_input, title) {
      shiny::reactive({
        value <- current_data()
        reduction <- precomputed_keep_selection(input[[reduction_input]], names(value$embeddings))
        embedding <- value$embeddings[[reduction]]
        color_by <- precomputed_scope_value(input$embedding_color, "celltype")
        color_values <- if (!is.null(color_by) && color_by %in% names(value$metadata)) value$metadata[[color_by]] else character()
        colors <- if (!identical(color_by, "abundance")) {
          precomputed_colors(value$contract, color_by, color_values)
        } else {
          NULL
        }
        plot_precomputed_embedding(
          metadata = value$metadata,
          embedding = embedding,
          color_by = color_by,
          abundance = value$abundance,
          marker = input$embedding_marker,
          colors = colors,
          title = paste(title, reduction)
        )
      })
    }
    umap_ggplot <- embedding_ggplot("umap_reduction", "Embedding:")
    pca_ggplot <- embedding_ggplot("pca_reduction", "PCA:")
    output$umap_plot <- plotly::renderPlotly({
      plotly::ggplotly(umap_ggplot(), tooltip = "text") |>
        plotly::layout(autosize = TRUE, margin = list(l = 60, r = 40, t = 70, b = 60))
    })
    output$pca_plot <- plotly::renderPlotly({
      plotly::ggplotly(pca_ggplot(), tooltip = "text") |>
        plotly::layout(autosize = TRUE, margin = list(l = 60, r = 40, t = 70, b = 60))
    })
    output$download_umap_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename("embedding", input$umap_reduction, input$embedding_color), ".svg"),
      content = function(file) precomputed_write_svg(umap_ggplot(), file, 9, 7)
    )
    output$download_pca_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename("embedding", input$pca_reduction, input$embedding_color), ".svg"),
      content = function(file) precomputed_write_svg(pca_ggplot(), file, 9, 7)
    )

    report_scope <- scope

    report_volcano_ggplot <- shiny::reactive({
      scope <- report_scope()
      if (identical(scope$celltype, PRECOMPUTED_ALL)) return(precomputed_empty_plot(
        "Select one cell type for its stored volcano.", "Original report volcano"))
      if (!length(scope$readout) || is.na(scope$readout) || !length(scope$donor) || is.na(scope$donor) ||
          !length(scope$celltype) || is.na(scope$celltype) || !length(scope$comparison) || is.na(scope$comparison)) {
        return(precomputed_empty_plot("No stored donor-level volcano is available.", "Original report volcano"))
      }
      precomputed_report_plot(
        "plot_report_volcano",
        contract = current_data()$contract,
        readout = scope$readout,
        donor = scope$donor,
        celltype = scope$celltype,
        comparison = scope$comparison,
        title = "Original report volcano"
      )
    })
    output$report_volcano_plot <- plotly::renderPlotly({
      report_volcano_plotly(report_volcano_ggplot(), source = session$ns("volcano"))
    })

    report_volcano_display <- shiny::reactive(
      precomputed_display_rows(report_scope()$rows, input$fdr_filter, input$fdr_cutoff))
    output$report_volcano_table <- DT::renderDT({
      shiny::validate(shiny::need(nrow(report_volcano_display()) > 0L, "No stored rows match this display scope."))
      DT::datatable(report_volcano_display(), rownames = FALSE, filter = "top",
        selection = "single", options = list(pageLength = 15, scrollX = TRUE))
    })
    output$download_report_volcano_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(
        report_scope()$readout, report_scope()$donor, report_scope()$celltype,
        report_scope()$comparison, "report-volcano"
      ), ".svg"),
      content = function(file) precomputed_write_svg(report_volcano_ggplot(), file, 10, 8)
    )
    output$download_report_volcano_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(
        report_scope()$readout, report_scope()$donor, report_scope()$celltype,
        report_scope()$comparison, "stored-results"
      ), ".csv"),
      content = function(file) utils::write.csv(report_scope()$rows, file, row.names = FALSE, na = "")
    )

    report_abundance_heatmap <- shiny::reactive({
      scope <- report_scope()
      contract <- current_data()$contract
      if (!identical(scope$readout, "abundance")) {
        return(precomputed_empty_plot("Select Abundance to view its frozen report heatmap.", "Original abundance heatmap"))
      }
      if (is.null(contract$report$abundance_heatmaps) ||
          !exists("plot_report_abundance_heatmap", mode = "function", inherits = TRUE)) {
        return(precomputed_empty_plot("No frozen abundance heatmap is stored for this dataset.", "Original abundance heatmap"))
      }
      tryCatch(
        get("plot_report_abundance_heatmap", mode = "function", inherits = TRUE)(
          contract = contract,
          donor = scope$donor,
          celltype = scope$celltype,
          comparison = scope$comparison
        ),
        error = function(error) precomputed_empty_plot(conditionMessage(error), "Original abundance heatmap")
      )
    })
    output$report_abundance_heatmap_plot <- shiny::renderPlot({
      precomputed_draw_plot_object(report_abundance_heatmap())
    }, res = 110)
    output$download_report_abundance_heatmap_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(
        report_scope()$donor, report_scope()$celltype, report_scope()$comparison,
        "abundance-report-heatmap"
      ), ".svg"),
      content = function(file) precomputed_write_plot_object_svg(report_abundance_heatmap(), file, 12, 9)
    )

    spatial_scope <- shiny::reactive({
      shared <- scope()
      contract <- current_data()$contract
      alignment <- contract$report$alignment
      kind <- precomputed_scope_value(input$spatial_kind, "differences")
      if (!kind %in% c("differences", "profiles")) kind <- "differences"
      donor_rows <- precomputed_filter_table(
        precomputed_contract_table(contract, "donor_contrasts", shared$readout),
        donor = shared$donor, comparison = shared$comparison, include_controls = TRUE)
      groups <- unique(c(as.character(donor_rows$group_a), as.character(donor_rows$group_b)))
      if (identical(kind, "profiles")) {
        panels <- intersect(precomputed_report_values(alignment$conditions, "condition"), groups)
        panel <- precomputed_keep_selection(input$spatial_panel, panels)
      } else {
        candidates <- unlist(lapply(contract$report$source_results[[shared$donor]], function(rows) {
          if (!all(c("contrast", "target", "reference") %in% names(rows))) return(character())
          matched <- as.character(rows$target) %in% donor_rows$group_a &
                     as.character(rows$reference) %in% donor_rows$group_b
          unique(as.character(rows$contrast[matched]))
        }), use.names = FALSE)
        panels <- intersect(precomputed_report_values(alignment$contrasts, "comparison"), candidates)
        panel <- if (length(panels) == 1L) panels else NA_character_
      }
      list(kind = kind, donor = shared$donor, celltype = shared$celltype, panels = panels, panel = panel)
    })
    shiny::observeEvent(list(scope(), input$spatial_kind), {
      selected <- spatial_scope()
      shiny::updateSelectInput(session, "spatial_panel", choices = selected$panels, selected = selected$panel)
    })
    spatial_report_ggplot <- shiny::reactive({
      scope <- spatial_scope()
      if (identical(scope$celltype, PRECOMPUTED_ALL)) return(precomputed_empty_plot(
        "Select one cell type for its spatial matrix.", "Spatial matrix"))
      if (!length(scope$donor) || is.na(scope$donor) || !length(scope$celltype) || is.na(scope$celltype) ||
          !length(scope$panel) || is.na(scope$panel)) {
        return(precomputed_empty_plot("No aligned spatial report panel is stored.", "Original spatial report view"))
      }
      precomputed_report_plot(
        "plot_report_spatial",
        contract = current_data()$contract,
        donor = scope$donor,
        celltype = scope$celltype,
        kind = scope$kind,
        panel = scope$panel,
        title = "Original spatial report view"
      )
    })
    output$spatial_report_plot <- plotly::renderPlotly(
      report_spatial_plotly(spatial_report_ggplot(), source = session$ns("spatial")))
    output$download_spatial_report_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(
        spatial_scope()$donor, spatial_scope()$celltype, spatial_scope()$kind,
        spatial_scope()$panel, "spatial-report"
      ), ".svg"),
      content = function(file) precomputed_write_svg(spatial_report_ggplot(), file, 12, 9)
    )

    output$fdr_controls <- shiny::renderUI({
      results <- selected_results()
      has_fdr <- "p_adj" %in% names(results) && any(is.finite(results$p_adj))
      if (!has_fdr) {
        return(shiny::tags$div(
          class = "alert alert-info small mt-2",
          shiny::tags$strong("Stored probabilities unavailable. "),
          "No adjusted P values were stored, so FDR controls are disabled."
        ))
      }
      shiny::tagList(
        shiny::numericInput(session$ns("fdr_cutoff"), "Stored FDR display cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        shiny::checkboxInput(session$ns("fdr_filter"), "Show only rows at or below this stored FDR", value = FALSE),
        shiny::tags$p(class = "text-muted small", "Table filters affect the tables and Filtered CSV downloads only. Plots and Original report SVG retain frozen thresholds; All stored rows ignores table searches and FDR filters.")
      )
    })

    output$result_summary <- shiny::renderUI({
      selected <- scope()
      rows <- selected$rows
      groups <- if (nrow(rows)) unique(paste(rows$group_a, "−", rows$group_b)) else character()
      direction <- if (length(groups)) paste(groups, collapse = "; ") else
        paste(selected$comparison, "— no stored effects in this scope")
      shiny::div(class = "alert alert-light border",
        shiny::tags$strong(paste("Donor:", selected$donor, "·",
          if (selected$celltype == PRECOMPUTED_ALL) "All cell types" else selected$celltype, "·",
          direction, "·", names(PRECOMPUTED_READOUTS)[match(selected$readout, PRECOMPUTED_READOUTS)])),
        shiny::tags$p(class = "small mb-1", "Positive effects mean target exceeds reference. ",
          "Stored single-donor results; no analysis was recalculated."),
        shiny::tags$p(class = "small mb-0",
          "Exploratory component-level tests, where stored, are not donor-population inference. ",
          nrow(rows), " stored effects. Missing or untested values are not zero.")
      )
    })

    heatmap_ggplot <- shiny::reactive({
      readout_label <- names(PRECOMPUTED_READOUTS)[match(readout_tables()$readout, PRECOMPUTED_READOUTS)]
      plot_precomputed_heatmap(
        precomputed_filter_table(selected_results(), include_controls = input$include_controls),
        top_n = input$top_features,
        title = paste(readout_label, "stored differential effects")
      )
    })
    output$effect_heatmap <- plotly::renderPlotly({
      plot <- heatmap_ggplot()
      rows <- scope()$rows
      if (is.data.frame(plot$data) && nrow(plot$data) && nrow(rows)) {
        i <- match(plot$data$feature_id, rows$feature_id)
        plot$data$key <- report_feature_key(rows$marker_1[i], rows$marker_2[i])
      }
      precomputed_heatmap_plotly(plot, source = session$ns("summary"))
    })

    detail_distribution <- shiny::reactive({
      selected <- scope()
      rows <- precomputed_filter_table(readout_tables()$contrasts, donor = selected$donor,
        comparison = selected$comparison, include_controls = TRUE)
      conditions <- unique(c(as.character(rows$group_b), as.character(rows$group_a)))
      pair <- tryCatch(jsonlite::fromJSON(selected_key()), error = function(e) character())
      precomputed_feature_distribution_data(current_data(), selected$readout, selected$donor,
        selected$celltype, conditions, pair)
    })
    detail_ggplot <- shiny::reactive({
      value <- current_data()
      contrasts <- detail_contrasts()
      donor_colors <- precomputed_colors(value$contract, "donor", contrasts$donor)
      pair <- tryCatch(jsonlite::fromJSON(selected_key()), error = function(e) character())
      feature_label <- paste(unique(pair), collapse = " ↔ ")
      if (nrow(contrasts) > 0L) {
        feature_label <- precomputed_feature_labels(contrasts)[1L]
      } else if (nrow(detail_results()) > 0L) {
        feature_label <- precomputed_feature_labels(detail_results())[1L]
      }
      if (is.null(input$results_view) || identical(input$results_view, "Tables")) {
        tryCatch({
          rows <- detail_distribution()
          plot_precomputed_feature_distribution(rows, scope()$readout,
            colors = precomputed_colors(value$contract, "condition", as.character(rows$condition)),
            title = paste(feature_label, "cell-level distribution"))
        }, error = function(e) precomputed_empty_plot(conditionMessage(e), "Cell-level distribution unavailable"))
      } else if (identical(input$detail_display, "conditions")) {
        plot_precomputed_condition_detail(
          contrasts,
          colors = donor_colors,
          title = paste(feature_label, "paired condition values")
        )
      } else {
        plot_precomputed_effect_detail(
          detail_results(),
          contrasts,
          colors = donor_colors,
          title = paste(feature_label, "stored effects")
        )
      }
    })
    output$effect_detail <- plotly::renderPlotly({
      plotly::ggplotly(detail_ggplot(), tooltip = if (is.null(input$results_view) ||
        identical(input$results_view, "Tables")) "all" else "text") |>
        plotly::layout(autosize = TRUE, margin = list(l = 80, r = 50, t = 100, b = 100))
    })

    display_results <- shiny::reactive(
      precomputed_display_rows(selected_results(), input$fdr_filter, input$fdr_cutoff))
    output$results_table <- DT::renderDT({
      shiny::validate(shiny::need(nrow(display_results()) > 0L, "No stored result rows match this display scope."))
      DT::datatable(display_results(), rownames = FALSE, filter = "top",
        selection = "single", options = list(pageLength = 15, scrollX = TRUE))
    })
    output$samples_table <- DT::renderDT({
      shiny::validate(shiny::need(nrow(selected_samples()) > 0L, "No stored sample support rows match this scope."))
      DT::datatable(selected_samples(), rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
    })

    output$download_results_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(scope()$readout, scope()$donor, scope()$comparison, scope()$celltype, "results"), ".csv"),
      content = function(file) utils::write.csv(selected_results(), file, row.names = FALSE, na = "")
    )
    output$download_filtered_results_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(scope()$donor, scope()$celltype,
        scope()$comparison, scope()$readout, "filtered-results"), ".csv"),
      content = function(file) utils::write.csv(
        precomputed_searched_rows(display_results(), input$results_table_rows_all,
          input$results_table_search, input$results_table_search_columns),
        file, row.names = FALSE, na = "")
    )
    output$download_filtered_report_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(scope()$donor, scope()$celltype,
        scope()$comparison, scope()$readout, "filtered-report-rows"), ".csv"),
      content = function(file) utils::write.csv(
        precomputed_searched_rows(report_volcano_display(), input$report_volcano_table_rows_all,
          input$report_volcano_table_search, input$report_volcano_table_search_columns),
        file, row.names = FALSE, na = "")
    )
    spatial_rows <- shiny::reactive({
      selected <- spatial_scope()
      rows <- current_data()$contract$report$alignment[[selected$kind]]
      if (!is.data.frame(rows)) return(precomputed_empty_table())
      rows <- as.data.frame(rows)
      rows[!is.na(rows$panel) & rows$donor == selected$donor & rows$celltype == selected$celltype &
        rows$panel %in% selected$panel, , drop = FALSE]
    })
    output$download_spatial_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(spatial_scope()$donor,
        spatial_scope()$celltype, spatial_scope()$kind, spatial_scope()$panel, "spatial"), ".csv"),
      content = function(file) utils::write.csv(spatial_rows(), file, row.names = FALSE, na = "")
    )
    output$methods_ui <- shiny::renderUI({
      selected <- scope()
      rows <- selected$rows
      fields <- if (nrow(rows)) intersect(c("metric", "method", "replicate_unit", "inference", "adjustment_family"), names(rows)) else character()
      settings <- current_data()$contract$report$settings[[selected$donor]]
      section <- switch(selected$readout, abundance = "abundance", self_proximity = "clustering",
        colocalization = "colocalization")
      shiny::tagList(
        if (!nrow(rows)) shiny::tags$p("No stored differential test is available in this scope."),
        shiny::tags$p("These are frozen single-donor contrasts. Cells are not independent biological donors; ",
          "stored component-level P values are exploratory. Display filters never rerun or readjust tests."),
        shiny::tags$dl(lapply(fields, function(field) shiny::tagList(shiny::tags$dt(field),
          shiny::tags$dd(paste(unique(as.character(rows[[field]])), collapse = "; "))))),
        shiny::tags$p("Spatial condition profiles use the stored profile summary (mean in this demo); ",
          "spatial differential effects use differences of observed medians. Profile differences must not be substituted for stored effects. ",
          "Condition values are shown only if persisted. Missing stored entries are not replaced with zero; ",
          "profile means retain the upstream missing-observation policy shown in the frozen settings."),
        shiny::tags$details(shiny::tags$summary("Frozen plot and filtering settings"),
          shiny::tags$pre(as.character(jsonlite::toJSON(list(
            analysis = settings[[section]], thresholds = settings$thresholds, spatial = settings$spatial),
            auto_unbox = TRUE, pretty = TRUE, null = "null"))))
      )
    })
    output$download_samples_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(scope()$readout, scope()$donor, scope()$comparison, scope()$celltype, "sample-support"), ".csv"),
      content = function(file) utils::write.csv(selected_samples(), file, row.names = FALSE, na = "")
    )
    output$download_heatmap_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(scope()$readout, scope()$donor, scope()$comparison, scope()$celltype, "heatmap"), ".svg"),
      content = function(file) precomputed_write_svg(heatmap_ggplot(), file, 12, 8)
    )
    output$download_detail_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(readout_tables()$readout, input$detail_feature,
        if (is.null(input$results_view) || identical(input$results_view, "Tables")) "distribution" else input$detail_display), ".svg"),
      content = function(file) precomputed_write_svg(detail_ggplot(), file, 10, 7)
    )

    invisible(list(
      results = selected_results,
      donor_contrasts = selected_contrasts,
      sample_support = selected_samples,
      heatmap_plot = heatmap_ggplot,
      detail_plot = detail_ggplot
    ))
  })
}
