differential_plot_row <- function(volcano_output_id, detail_output_id) {
  div(
    class = "differential-plot-row",
    plot_pane(
      size = "compact",
      extra_class = "differential-plot-pane",
      plotlyOutput(volcano_output_id, width = "auto", height = "auto")
    ),
    plot_pane(
      size = "standard",
      extra_class = "differential-plot-pane",
      plotlyOutput(detail_output_id, width = "auto", height = "auto")
    )
  )
}

make_differential_config <- function(
  group_a,
  group_b,
  celltype_filter,
  stratify_by_celltype = FALSE,
  min_cells = 3,
  fdr_cutoff = 0.05,
  effect_cutoff = 0.25,
  anchor_marker = NULL
) {
  list(
    group_a = group_a,
    group_b = group_b,
    celltype_filter = celltype_filter,
    stratify_by_celltype = isTRUE(stratify_by_celltype),
    min_cells = min_cells,
    fdr_cutoff = fdr_cutoff,
    effect_cutoff = effect_cutoff,
    anchor_marker = anchor_marker
  )
}

default_differential_config <- function(
  conditions,
  cell_types,
  anchor_marker = NULL,
  min_cells = 3,
  fdr_cutoff = 0.05,
  effect_cutoff = 0.25
) {
  make_differential_config(
    group_a = conditions[1],
    group_b = conditions[min(2, length(conditions))],
    celltype_filter = cell_types,
    stratify_by_celltype = FALSE,
    min_cells = min_cells,
    fdr_cutoff = fdr_cutoff,
    effect_cutoff = effect_cutoff,
    anchor_marker = anchor_marker
  )
}

differential_volcano_dimensions <- function(x_label, legend = TRUE) {
  label_width <- nchar(as.character(x_label %||% ""), type = "width")
  width <- bounded_integer(360 + label_width * 7, lower = 620, upper = 1100)
  bottom_margin <- bounded_integer(86 + label_width * 1.2, lower = 130, upper = 180)
  right_margin <- if (isTRUE(legend)) 160 else 96

  list(
    width = width,
    height = 500,
    margin = list(l = 92, r = right_margin, t = 56, b = bottom_margin)
  )
}

differential_detail_dimensions <- function(plot_data, stratify_by_celltype = FALSE, y_label = NULL) {
  group_count <- if ("condition" %in% names(plot_data)) count_observed_values(plot_data$condition) else 2L
  facet_count <- if (isTRUE(stratify_by_celltype) && "celltype_manual" %in% names(plot_data)) {
    count_observed_values(plot_data$celltype_manual)
  } else {
    1L
  }
  facet_cols <- if (facet_count <= 1) 1L else ceiling(sqrt(facet_count))
  facet_rows <- ceiling(facet_count / facet_cols)
  y_label_width <- nchar(as.character(y_label %||% ""), type = "width")
  left_margin <- bounded_integer(70 + y_label_width * 1.2, lower = 90, upper = 180)
  panel_width <- max(260, group_count * 150)
  panel_height <- 250

  list(
    width = max(560, panel_width * facet_cols + left_margin + 72),
    height = max(500, panel_height * facet_rows + 170),
    margin = list(l = left_margin, r = 72, t = 56, b = 110),
    facet_count = facet_count,
    facet_cols = facet_cols,
    facet_rows = facet_rows,
    group_count = group_count
  )
}

differential_summary_row <- function(result, fdr_cutoff, effect_cutoff) {
  hits <- filter_differential_hits(result, fdr_cutoff = fdr_cutoff, effect_cutoff = effect_cutoff)
  up <- sum(hits$effect_size > 0, na.rm = TRUE)
  down <- sum(hits$effect_size < 0, na.rm = TRUE)
  eligible <- sum(!is.na(result$p_adj))

  metric_row(
    metric_tile("Tested", format(nrow(result), big.mark = ",")),
    metric_tile("With FDR", format(eligible, big.mark = ",")),
    metric_tile("Higher in A", format(up, big.mark = ",")),
    metric_tile("Higher in B", format(down, big.mark = ","))
  )
}

filter_differential_hits <- function(result, fdr_cutoff, effect_cutoff) {
  if (nrow(result) == 0) {
    return(result)
  }

  keep <- !is.na(result$p_adj) &
    result$p_adj <= fdr_cutoff &
    !is.na(result$effect_size) &
    abs(result$effect_size) >= effect_cutoff

  result <- result[keep, , drop = FALSE]
  order_differential_results(result)
}

order_differential_results <- function(result) {
  if (nrow(result) == 0) {
    return(result)
  }

  result[order(is.na(result$p_adj), result$p_adj, -abs(result$effect_size)), , drop = FALSE]
}

differential_volcano_plot <- function(
  result,
  label_col,
  x_label,
  fdr_cutoff,
  effect_cutoff,
  source = NULL,
  dimensions = differential_volcano_dimensions(x_label)
) {
  plot_data <- prepare_differential_plot_data(result, label_col, fdr_cutoff, effect_cutoff)
  fdr_line <- -log10(max(fdr_cutoff, .Machine$double.xmin))

  p <- ggplot(plot_data, aes(effect_size, neg_log10_fdr, color = threshold_status, text = hover, key = plot_key)) +
    geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), color = "#b3bdbf", linetype = "dashed", linewidth = 0.45) +
    geom_hline(yintercept = fdr_line, color = "#b3bdbf", linetype = "dashed", linewidth = 0.45) +
    geom_point(size = 2.3, alpha = 0.78) +
    scale_color_manual(
      values = c(
        "Higher in A" = "#c7503e",
        "Higher in B" = "#176d73",
        "Not significant" = "#9aa5a8"
      ),
      name = NULL
    ) +
    labs(x = x_label, y = "-log10(FDR)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ggplotly(p, tooltip = "text", source = source, width = dimensions$width, height = dimensions$height) |>
    apply_differential_plot_frame(dimensions = dimensions)
}

prepare_differential_plot_data <- function(result, label_col, fdr_cutoff, effect_cutoff) {
  plot_data <- result
  plot_data$p_adj_for_plot <- ifelse(is.na(plot_data$p_adj), 1, pmax(plot_data$p_adj, .Machine$double.xmin))
  plot_data$neg_log10_fdr <- -log10(plot_data$p_adj_for_plot)
  plot_data$threshold_status <- differential_threshold_status(plot_data, fdr_cutoff, effect_cutoff)
  plot_data$plot_label <- paste(plot_data[[label_col]], plot_data$celltype_manual, sep = " | ")
  plot_data$plot_key <- plot_data[[label_col]]
  plot_data$hover <- paste0(
    plot_data[[label_col]],
    "<br>Cell type: ", plot_data$celltype_manual,
    "<br>Effect: ", round(plot_data$effect_size, 3),
    "<br>FDR: ", format_p_value(plot_data$p_adj),
    "<br>", plot_data$group_a, " mean: ", round(plot_data$mean_a, 3),
    "<br>", plot_data$group_b, " reference mean: ", round(plot_data$mean_b, 3),
    "<br>", plot_data$direction
  )
  plot_data
}

differential_threshold_status <- function(result, fdr_cutoff, effect_cutoff) {
  is_hit <- !is.na(result$p_adj) &
    result$p_adj <= fdr_cutoff &
    !is.na(result$effect_size) &
    abs(result$effect_size) >= effect_cutoff

  ifelse(
    is_hit,
    ifelse(result$effect_size > 0, "Higher in A", "Higher in B"),
    "Not significant"
  )
}

format_differential_table <- function(result, effect_label, max_rows = 30) {
  result <- order_differential_results(result)
  result <- head(result, max_rows)

  numeric_cols <- intersect(c("mean_a", "mean_b", "median_a", "median_b", "effect_size"), names(result))
  for (col in numeric_cols) {
    result[[col]] <- round(result[[col]], 3)
  }
  result$p_value <- format_p_value(result$p_value)
  result$p_adj <- format_p_value(result$p_adj)
  names(result)[names(result) == "effect_size"] <- effect_label
  result
}

format_p_value <- function(values) {
  ifelse(
    is.na(values),
    NA_character_,
    formatC(values, format = "e", digits = 2)
  )
}
