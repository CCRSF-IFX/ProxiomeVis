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

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, color = "#8a9699", linewidth = 0.5) +
    ggplot2::labs(
      title = title,
      subtitle = "Dots are stored donor contrasts; diamonds and ranges are stored cross-donor summaries.",
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
        ggplot2::aes(x = analysis_label, ymin = min_effect, ymax = max_effect, text = hover),
        linewidth = 1.1,
        color = "#263238"
      ) +
      ggplot2::geom_point(
        data = aggregate,
        ggplot2::aes(x = analysis_label, y = effect_size, text = hover),
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
        ggplot2::aes(x = analysis_label, y = effect_size, color = donor, text = hover),
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
  shiny::div(
    class = "col-6 col-lg-3 mb-3",
    bslib::card(
      class = "h-100",
      bslib::card_body(
        shiny::tags$div(class = "text-muted small text-uppercase", label),
        shiny::tags$div(class = "fs-3 fw-semibold", shiny::textOutput(output_id, inline = TRUE))
      )
    )
  )
}

precomputed_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::navset_card_underline(
    id = ns("report_page"),
    full_screen = TRUE,
    bslib::nav_panel(
      "Overview / QC",
      shiny::uiOutput(ns("dataset_heading")),
      shiny::div(
        class = "row",
        precomputed_metric_card("Cells", ns("metric_cells")),
        precomputed_metric_card("Samples", ns("metric_samples")),
        precomputed_metric_card("Donors", ns("metric_donors")),
        precomputed_metric_card("Cell types", ns("metric_celltypes"))
      ),
      shiny::fluidRow(
        shiny::column(
          width = 7,
          bslib::card(
            bslib::card_header("QC cell counts"),
            bslib::card_body(plotly::plotlyOutput(ns("qc_plot"), height = "480px"))
          )
        ),
        shiny::column(
          width = 5,
          bslib::card(
            bslib::card_header("QC table"),
            bslib::card_body(DT::DTOutput(ns("qc_table")))
          )
        )
      ),
      bslib::card(
        class = "mt-3",
        bslib::card_header("Analysis provenance"),
        bslib::card_body(
          precomputed_download_group(
            shiny::downloadButton(ns("download_provenance"), "Provenance JSON", class = "btn-sm btn-outline-primary")
          ),
          shiny::tags$pre(class = "small", shiny::textOutput(ns("provenance_text")))
        )
      )
    ),
    bslib::nav_panel(
      "Embeddings",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Embedding controls",
          width = 300,
          shiny::selectInput(ns("umap_reduction"), "UMAP / annotation reduction", choices = character()),
          shiny::selectInput(ns("pca_reduction"), "PCA reduction", choices = character()),
          shiny::selectInput(
            ns("embedding_color"),
            "Color cells by",
            choices = c(
              "Cell type" = "celltype",
              "Cluster" = "cluster",
              "Condition" = "condition",
              "Donor" = "donor",
              "Marker abundance" = "abundance"
            ),
            selected = "celltype"
          ),
          shiny::conditionalPanel(
            condition = "input.embedding_color == 'abundance'",
            ns = ns,
            shiny::selectizeInput(ns("embedding_marker"), "Abundance marker", choices = character())
          ),
          shiny::tags$p(
            class = "text-muted small",
            "Embedding colors are display annotations or normalized abundance from the loaded RDS; no spatial file is read."
          )
        ),
        shiny::fluidRow(
          shiny::column(
            width = 6,
            bslib::card(
              bslib::card_header("UMAP / annotation embedding"),
              bslib::card_body(
                precomputed_download_group(
                  shiny::downloadButton(ns("download_umap_svg"), "SVG", class = "btn-sm btn-outline-primary")
                ),
                plotly::plotlyOutput(ns("umap_plot"), height = "620px")
              )
            )
          ),
          shiny::column(
            width = 6,
            bslib::card(
              bslib::card_header("PCA embedding"),
              bslib::card_body(
                precomputed_download_group(
                  shiny::downloadButton(ns("download_pca_svg"), "SVG", class = "btn-sm btn-outline-primary")
                ),
                plotly::plotlyOutput(ns("pca_plot"), height = "620px")
              )
            )
          )
        )
      )
    ),
    bslib::nav_panel(
      "Precomputed results",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Stored analysis scope",
          width = 320,
          shiny::selectInput(ns("readout"), "Readout", choices = PRECOMPUTED_READOUTS),
          shiny::selectInput(ns("donor"), "Donor view", choices = c("All stored donors" = PRECOMPUTED_ALL)),
          shiny::selectInput(ns("comparison"), "Stored comparison", choices = c("All comparisons" = PRECOMPUTED_ALL)),
          shiny::selectInput(ns("celltype"), "Cell type", choices = c("All cell types" = PRECOMPUTED_ALL)),
          shiny::checkboxInput(ns("include_controls"), "Include control markers", value = FALSE),
          shiny::numericInput(ns("top_features"), "Heatmap top features", value = 20, min = 2, max = 80, step = 1),
          shiny::tags$hr(),
          shiny::selectizeInput(ns("detail_feature"), "Detail feature", choices = character()),
          shiny::radioButtons(
            ns("detail_display"),
            "Detail display",
            choices = c("Donor effects and stored range" = "effects", "Paired condition values" = "conditions"),
            selected = "effects"
          ),
          shiny::uiOutput(ns("fdr_controls")),
          shiny::tags$div(
            class = "alert alert-secondary small mt-3",
            "All donors shows the stored cross-donor summary. Selecting one donor switches to that donor's stored contrast; it does not recompute or pool an analysis."
          )
        ),
        shiny::uiOutput(ns("result_summary")),
        bslib::card(
          bslib::card_header("Differential-effect heatmap"),
          bslib::card_body(
            precomputed_download_group(
              shiny::downloadButton(ns("download_results_csv"), "Results CSV", class = "btn-sm btn-outline-primary"),
              shiny::downloadButton(ns("download_heatmap_svg"), "Heatmap SVG", class = "btn-sm btn-outline-primary")
            ),
            plotly::plotlyOutput(ns("effect_heatmap"), height = "720px")
          )
        ),
        bslib::card(
          class = "mt-3",
          bslib::card_header("Selected-feature detail"),
          bslib::card_body(
            precomputed_download_group(
              shiny::downloadButton(ns("download_detail_svg"), "Detail SVG", class = "btn-sm btn-outline-primary")
            ),
            plotly::plotlyOutput(ns("effect_detail"), height = "540px")
          )
        ),
        bslib::navset_card_tab(
          id = ns("result_tables"),
          bslib::nav_panel("Results", DT::DTOutput(ns("results_table"))),
          bslib::nav_panel(
            "Donor contrasts",
            precomputed_download_group(
              shiny::downloadButton(ns("download_contrasts_csv"), "Donor contrasts CSV", class = "btn-sm btn-outline-primary")
            ),
            DT::DTOutput(ns("contrasts_table"))
          ),
          bslib::nav_panel(
            "Sample support",
            precomputed_download_group(
              shiny::downloadButton(ns("download_samples_csv"), "Sample support CSV", class = "btn-sm btn-outline-primary")
            ),
            DT::DTOutput(ns("samples_table"))
          )
        )
      )
    ),
    bslib::nav_panel(
      "Report volcano",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Stored report scope",
          width = 320,
          shiny::selectInput(ns("report_readout"), "Readout", choices = PRECOMPUTED_READOUTS),
          shiny::selectInput(ns("report_donor"), "Donor", choices = character()),
          shiny::selectInput(ns("report_comparison"), "Comparison", choices = character()),
          shiny::selectInput(ns("report_celltype"), "Cell type", choices = character()),
          shiny::uiOutput(ns("report_fdr_controls")),
          shiny::tags$div(
            class = "alert alert-warning small mt-3",
            shiny::tags$strong("Exploratory component-level test. "),
            "The volcano and P values are the original stored result for one donor; they are not donor-population inference. The cutoff affects the table only."
          )
        ),
        bslib::card(
          bslib::card_header("Original report volcano"),
          bslib::card_body(
            precomputed_download_group(
              shiny::downloadButton(ns("download_report_volcano_svg"), "Volcano SVG", class = "btn-sm btn-outline-primary"),
              shiny::downloadButton(ns("download_report_volcano_csv"), "Stored rows CSV", class = "btn-sm btn-outline-primary")
            ),
            shiny::plotOutput(ns("report_volcano_plot"), height = "720px")
          )
        ),
        bslib::card(
          class = "mt-3",
          bslib::card_header("Stored component-level results"),
          bslib::card_body(DT::DTOutput(ns("report_volcano_table")))
        ),
        bslib::card(
          class = "mt-3",
          bslib::card_header("Original abundance heatmap"),
          bslib::card_body(
            precomputed_download_group(
              shiny::downloadButton(ns("download_report_abundance_heatmap_svg"), "Abundance heatmap SVG", class = "btn-sm btn-outline-primary")
            ),
            shiny::plotOutput(ns("report_abundance_heatmap_plot"), height = "820px")
          )
        )
      )
    ),
    bslib::nav_panel(
      "Spatial report",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Stored spatial panel",
          width = 320,
          shiny::selectInput(ns("spatial_donor"), "Donor", choices = character()),
          shiny::selectInput(ns("spatial_celltype"), "Cell type", choices = character()),
          shiny::radioButtons(
            ns("spatial_kind"),
            "Display",
            choices = c("Condition differences" = "differences", "Condition profiles" = "profiles"),
            selected = "differences"
          ),
          shiny::selectInput(ns("spatial_panel"), "Comparison / condition", choices = character()),
          shiny::tags$p(
            class = "text-muted small",
            "Marker order and color/size scales come from the frozen report alignment bundle; the app does not recalculate proximity."
          )
        ),
        bslib::card(
          bslib::card_header("Original spatial report view"),
          bslib::card_body(
            precomputed_download_group(
              shiny::downloadButton(ns("download_spatial_report_svg"), "Spatial report SVG", class = "btn-sm btn-outline-primary")
            ),
            shiny::plotOutput(ns("spatial_report_plot"), height = "820px")
          )
        )
      )
    )
  )
}

precomputed_server <- function(id, data) {
  shiny::moduleServer(id, function(input, output, session) {
    current_data <- shiny::reactive({
      value <- data()
      shiny::req(value, value$contract, value$metadata, value$embeddings)
      value
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

    shiny::observeEvent(list(current_data(), input$readout), {
      tables <- readout_tables()
      donors <- precomputed_available_values(tables$contrasts, "donor")
      comparisons <- sort(unique(c(
        precomputed_available_values(tables$results, "comparison"),
        precomputed_available_values(tables$contrasts, "comparison")
      )))
      celltypes <- sort(unique(c(
        precomputed_available_values(tables$results, "celltype"),
        precomputed_available_values(tables$contrasts, "celltype")
      )))
      donor_choices <- c("All stored donors" = PRECOMPUTED_ALL, stats::setNames(donors, donors))
      comparison_choices <- c("All comparisons" = PRECOMPUTED_ALL, stats::setNames(comparisons, comparisons))
      celltype_choices <- c("All cell types" = PRECOMPUTED_ALL, stats::setNames(celltypes, celltypes))
      shiny::updateSelectInput(
        session,
        "donor",
        choices = donor_choices,
        selected = precomputed_keep_selection(shiny::isolate(input$donor), unname(donor_choices), PRECOMPUTED_ALL)
      )
      shiny::updateSelectInput(
        session,
        "comparison",
        choices = comparison_choices,
        selected = precomputed_keep_selection(shiny::isolate(input$comparison), unname(comparison_choices), PRECOMPUTED_ALL)
      )
      shiny::updateSelectInput(
        session,
        "celltype",
        choices = celltype_choices,
        selected = precomputed_keep_selection(shiny::isolate(input$celltype), unname(celltype_choices), PRECOMPUTED_ALL)
      )
    }, ignoreInit = FALSE)

    selected_results <- shiny::reactive({
      value <- current_data()
      tables <- readout_tables()
      rows <- precomputed_result_rows_for_donor(value$contract, tables$readout, input$donor)
      precomputed_filter_table(
        rows,
        donor = input$donor,
        comparison = input$comparison,
        celltype = input$celltype,
        include_controls = input$include_controls
      )
    })

    selected_contrasts <- shiny::reactive({
      tables <- readout_tables()
      precomputed_filter_table(
        tables$contrasts,
        donor = input$donor,
        comparison = input$comparison,
        celltype = input$celltype,
        include_controls = input$include_controls
      )
    })

    selected_samples <- shiny::reactive({
      tables <- readout_tables()
      samples <- precomputed_filter_table(
        tables$samples,
        donor = input$donor,
        celltype = input$celltype,
        include_controls = input$include_controls
      )
      if (!identical(precomputed_scope_value(input$comparison), PRECOMPUTED_ALL) && nrow(samples) > 0L && "condition" %in% names(samples)) {
        contrast_scope <- precomputed_filter_table(
          tables$contrasts,
          donor = input$donor,
          comparison = input$comparison,
          celltype = input$celltype,
          include_controls = input$include_controls
        )
        conditions <- unique(c(as.character(contrast_scope$group_a), as.character(contrast_scope$group_b)))
        conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
        if (length(conditions) > 0L) {
          samples <- samples[as.character(samples$condition) %in% conditions, , drop = FALSE]
        }
      }
      samples
    })

    shiny::observe({
      results <- selected_results()
      if (nrow(results) == 0L || !"feature_id" %in% names(results)) {
        shiny::updateSelectizeInput(session, "detail_feature", choices = character(), selected = character(), server = TRUE)
        return()
      }
      feature_rows <- results[!duplicated(results$feature_id), , drop = FALSE]
      choices <- stats::setNames(as.character(feature_rows$feature_id), precomputed_feature_labels(feature_rows))
      selected <- precomputed_keep_selection(shiny::isolate(input$detail_feature), unname(choices))
      shiny::freezeReactiveValue(input, "detail_feature")
      shiny::updateSelectizeInput(session, "detail_feature", choices = choices, selected = selected, server = TRUE)
    })

    detail_results <- shiny::reactive({
      if (!identical(precomputed_scope_value(input$donor), PRECOMPUTED_ALL)) {
        return(precomputed_empty_table())
      }
      results <- selected_results()
      feature <- input$detail_feature
      if (is.null(feature) || !"feature_id" %in% names(results)) {
        return(precomputed_empty_table())
      }
      results[as.character(results$feature_id) == feature, , drop = FALSE]
    })

    detail_contrasts <- shiny::reactive({
      contrasts <- selected_contrasts()
      feature <- input$detail_feature
      if (is.null(feature) || !"feature_id" %in% names(contrasts)) {
        return(precomputed_empty_table())
      }
      contrasts[as.character(contrasts$feature_id) == feature, , drop = FALSE]
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

    report_scope <- shiny::reactive({
      value <- current_data()
      readouts <- unname(PRECOMPUTED_READOUTS)[vapply(
        unname(PRECOMPUTED_READOUTS),
        function(readout) nrow(precomputed_contract_table(value$contract, "donor_contrasts", readout)) > 0L,
        logical(1)
      )]
      readout <- precomputed_keep_selection(input$report_readout, readouts)
      rows <- precomputed_contract_table(value$contract, "donor_contrasts", readout)
      donors <- precomputed_available_values(rows, "donor")
      donor <- precomputed_keep_selection(input$report_donor, donors)
      donor_rows <- precomputed_filter_table(rows, donor = donor, include_controls = TRUE)
      celltypes <- precomputed_report_celltypes(
        value$contract, value$metadata, donor, readout, donor_rows
      )
      celltype <- precomputed_keep_selection(input$report_celltype, celltypes)
      comparisons <- precomputed_available_values(donor_rows, "comparison")
      comparison <- precomputed_keep_selection(input$report_comparison, comparisons)
      list(
        readouts = readouts,
        readout = readout,
        donors = donors,
        donor = donor,
        celltypes = celltypes,
        celltype = celltype,
        comparisons = comparisons,
        comparison = comparison,
        rows = precomputed_filter_table(
          donor_rows,
          comparison = comparison,
          celltype = celltype,
          include_controls = TRUE
        )
      )
    })

    shiny::observe({
      scope <- report_scope()
      readout_choices <- PRECOMPUTED_READOUTS[PRECOMPUTED_READOUTS %in% scope$readouts]
      shiny::updateSelectInput(session, "report_readout", choices = readout_choices, selected = scope$readout)
      shiny::updateSelectInput(session, "report_donor", choices = scope$donors, selected = scope$donor)
      shiny::updateSelectInput(session, "report_celltype", choices = scope$celltypes, selected = scope$celltype)
      shiny::updateSelectInput(session, "report_comparison", choices = scope$comparisons, selected = scope$comparison)
    })

    output$report_fdr_controls <- shiny::renderUI({
      rows <- report_scope()$rows
      has_fdr <- "p_adj" %in% names(rows) && any(is.finite(rows$p_adj))
      if (!has_fdr) {
        return(shiny::tags$div(
          class = "alert alert-info small mt-2",
          "No adjusted P values were stored for this report scope; FDR filtering is disabled."
        ))
      }
      shiny::tagList(
        shiny::numericInput(session$ns("report_fdr_cutoff"), "Stored FDR table cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        shiny::checkboxInput(session$ns("report_fdr_filter"), "Show only table rows at or below cutoff", value = FALSE)
      )
    })

    report_volcano_ggplot <- shiny::reactive({
      scope <- report_scope()
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
    output$report_volcano_plot <- shiny::renderPlot(report_volcano_ggplot(), res = 110)

    report_volcano_display <- shiny::reactive({
      rows <- report_scope()$rows
      if (nrow(rows) == 0L) return(rows)
      rows$feature <- precomputed_feature_labels(rows)
      if ("p_adj" %in% names(rows) && any(is.finite(rows$p_adj))) {
        cutoff <- suppressWarnings(as.numeric(input$report_fdr_cutoff)[1L])
        if (!is.finite(cutoff)) cutoff <- 0.05
        rows$stored_fdr_status <- ifelse(
          is.finite(rows$p_adj),
          ifelse(rows$p_adj <= cutoff, paste0("FDR ≤ ", cutoff), paste0("FDR > ", cutoff)),
          "Not available"
        )
        if (isTRUE(input$report_fdr_filter)) {
          rows <- rows[is.finite(rows$p_adj) & rows$p_adj <= cutoff, , drop = FALSE]
        }
      }
      rows
    })
    output$report_volcano_table <- DT::renderDT({
      shiny::validate(shiny::need(nrow(report_volcano_display()) > 0L, "No stored rows match this display scope."))
      DT::datatable(report_volcano_display(), rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
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
      contract <- current_data()$contract
      alignment <- contract$report$alignment
      kind <- precomputed_scope_value(input$spatial_kind, "differences")
      if (!kind %in% c("differences", "profiles")) kind <- "differences"
      donors <- precomputed_report_values(alignment$donors, "donor")
      celltypes <- precomputed_report_values(alignment$markers, "celltype")
      if (length(donors) == 0L || length(celltypes) == 0L) {
        contrast_tables <- lapply(
          unname(PRECOMPUTED_READOUTS),
          function(readout) precomputed_contract_table(contract, "donor_contrasts", readout)
        )
        if (length(donors) == 0L) {
          donors <- sort(unique(unlist(lapply(contrast_tables, precomputed_available_values, column = "donor"), use.names = FALSE)))
        }
        if (length(celltypes) == 0L) {
          celltypes <- sort(unique(unlist(lapply(contrast_tables, precomputed_available_values, column = "celltype"), use.names = FALSE)))
        }
      }
      panels <- if (identical(kind, "profiles")) {
        precomputed_report_values(alignment$conditions, "condition")
      } else {
        precomputed_report_values(alignment$contrasts, "comparison")
      }
      list(
        kind = kind,
        donors = donors,
        donor = precomputed_keep_selection(input$spatial_donor, donors),
        celltypes = celltypes,
        celltype = precomputed_keep_selection(input$spatial_celltype, celltypes),
        panels = panels,
        panel = precomputed_keep_selection(input$spatial_panel, panels)
      )
    })
    shiny::observe({
      scope <- spatial_scope()
      shiny::updateSelectInput(session, "spatial_donor", choices = scope$donors, selected = scope$donor)
      shiny::updateSelectInput(session, "spatial_celltype", choices = scope$celltypes, selected = scope$celltype)
      shiny::updateSelectInput(session, "spatial_panel", choices = scope$panels, selected = scope$panel)
    })
    spatial_report_ggplot <- shiny::reactive({
      scope <- spatial_scope()
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
    output$spatial_report_plot <- shiny::renderPlot(spatial_report_ggplot(), res = 110)
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
          shiny::tags$strong("Descriptive result. "),
          "No adjusted P values were stored, so FDR controls are disabled."
        ))
      }
      shiny::tagList(
        shiny::numericInput(session$ns("fdr_cutoff"), "Stored FDR display cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        shiny::checkboxInput(session$ns("fdr_filter"), "Show only rows at or below this stored FDR", value = FALSE),
        shiny::tags$p(class = "text-muted small", "This filters the on-screen table only; it does not adjust or rerun tests, alter plots, or change downloads.")
      )
    })

    output$result_summary <- shiny::renderUI({
      results <- selected_results()
      contrasts <- selected_contrasts()
      donor <- precomputed_scope_value(input$donor)
      donor_message <- if (identical(donor, PRECOMPUTED_ALL)) {
        "Showing stored cross-donor summaries and their contributing donor contrasts."
      } else {
        paste("Showing stored contrasts for", donor,
              "only; component-level tests are exploratory and are not donor-population inference. No analysis was recalculated.")
      }
      features <- if ("feature_id" %in% names(results)) length(unique(results$feature_id)) else 0L
      shiny::div(
        class = "alert alert-light border",
        shiny::tags$strong(donor_message),
        shiny::tags$div(
          class = "small mt-1",
          paste(format(nrow(results), big.mark = ","), "stored effect rows ·",
                format(features, big.mark = ","), "features ·",
                format(nrow(contrasts), big.mark = ","), "donor contrast rows")
        )
      )
    })

    heatmap_ggplot <- shiny::reactive({
      readout_label <- names(PRECOMPUTED_READOUTS)[match(readout_tables()$readout, PRECOMPUTED_READOUTS)]
      plot_precomputed_heatmap(
        selected_results(),
        top_n = input$top_features,
        title = paste(readout_label, "stored differential effects")
      )
    })
    output$effect_heatmap <- plotly::renderPlotly({
      plotly::ggplotly(heatmap_ggplot(), tooltip = "text") |>
        plotly::layout(autosize = TRUE, margin = list(l = 140, r = 50, t = 100, b = 150))
    })

    detail_ggplot <- shiny::reactive({
      value <- current_data()
      contrasts <- detail_contrasts()
      donor_colors <- precomputed_colors(value$contract, "donor", contrasts$donor)
      feature_label <- input$detail_feature
      if (nrow(contrasts) > 0L) {
        feature_label <- precomputed_feature_labels(contrasts)[1L]
      } else if (nrow(detail_results()) > 0L) {
        feature_label <- precomputed_feature_labels(detail_results())[1L]
      }
      if (identical(input$detail_display, "conditions")) {
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
      plotly::ggplotly(detail_ggplot(), tooltip = "text") |>
        plotly::layout(autosize = TRUE, margin = list(l = 80, r = 50, t = 100, b = 100))
    })

    display_results <- shiny::reactive({
      rows <- selected_results()
      if (nrow(rows) == 0L) {
        return(rows)
      }
      rows$feature <- precomputed_feature_labels(rows)
      if ("p_adj" %in% names(rows) && any(is.finite(rows$p_adj))) {
        cutoff <- suppressWarnings(as.numeric(input$fdr_cutoff)[1])
        if (!is.finite(cutoff)) cutoff <- 0.05
        rows$stored_fdr_status <- ifelse(
          is.finite(rows$p_adj),
          ifelse(rows$p_adj <= cutoff, paste0("FDR ≤ ", cutoff), paste0("FDR > ", cutoff)),
          "Not available"
        )
        if (isTRUE(input$fdr_filter)) {
          rows <- rows[is.finite(rows$p_adj) & rows$p_adj <= cutoff, , drop = FALSE]
        }
      }
      rows
    })
    output$results_table <- DT::renderDT({
      shiny::validate(shiny::need(nrow(display_results()) > 0L, "No stored result rows match this display scope."))
      DT::datatable(display_results(), rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
    })
    output$contrasts_table <- DT::renderDT({
      shiny::validate(shiny::need(nrow(selected_contrasts()) > 0L, "No stored donor contrasts match this scope."))
      DT::datatable(selected_contrasts(), rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
    })
    output$samples_table <- DT::renderDT({
      shiny::validate(shiny::need(nrow(selected_samples()) > 0L, "No stored sample support rows match this scope."))
      DT::datatable(selected_samples(), rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
    })

    output$download_results_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(readout_tables()$readout, input$donor, input$comparison, input$celltype, "results"), ".csv"),
      content = function(file) utils::write.csv(selected_results(), file, row.names = FALSE, na = "")
    )
    output$download_contrasts_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(readout_tables()$readout, input$donor, input$comparison, input$celltype, "donor-contrasts"), ".csv"),
      content = function(file) utils::write.csv(selected_contrasts(), file, row.names = FALSE, na = "")
    )
    output$download_samples_csv <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(readout_tables()$readout, input$donor, input$comparison, input$celltype, "sample-support"), ".csv"),
      content = function(file) utils::write.csv(selected_samples(), file, row.names = FALSE, na = "")
    )
    output$download_heatmap_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(readout_tables()$readout, input$donor, input$comparison, input$celltype, "heatmap"), ".svg"),
      content = function(file) precomputed_write_svg(heatmap_ggplot(), file, 12, 8)
    )
    output$download_detail_svg <- shiny::downloadHandler(
      filename = function() paste0(precomputed_safe_filename(readout_tables()$readout, input$detail_feature, input$detail_display), ".svg"),
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
