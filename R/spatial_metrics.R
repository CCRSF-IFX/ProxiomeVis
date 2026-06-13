summarize_spatial_heatmap_by_sample <- function(
  proximity,
  selected_markers,
  min_cells_detected = 1L
) {
  summarize_spatial_heatmap(
    proximity = proximity,
    selected_markers = selected_markers,
    group_cols = c("sample_alias", "condition"),
    min_cells_detected = min_cells_detected
  )
}

summarize_spatial_heatmap_by_celltype <- function(
  proximity,
  selected_markers,
  min_cells_detected = 1L
) {
  summarize_spatial_heatmap(
    proximity = proximity,
    selected_markers = selected_markers,
    group_cols = c("sample_alias", "condition", "celltype_manual"),
    min_cells_detected = min_cells_detected
  )
}

summarize_spatial_heatmap <- function(
  proximity,
  selected_markers,
  group_cols,
  value_col = "log2_ratio",
  component_col = "component",
  marker1_col = "marker_1",
  marker2_col = "marker_2",
  min_cells_detected = 1L
) {
  required_cols <- unique(c(group_cols, component_col, marker1_col, marker2_col, value_col))
  missing_cols <- setdiff(required_cols, names(proximity))
  if (length(missing_cols) > 0) {
    stop("Missing columns for spatial heatmap summary: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  selected_markers <- as.character(selected_markers)
  require_spatial_namespace("data.table")
  proximity_dt <- data.table::as.data.table(proximity)
  rows <- proximity_dt[
    get(marker1_col) %in% selected_markers &
      get(marker2_col) %in% selected_markers &
      get(marker1_col) != get(marker2_col) &
      is.finite(as.numeric(get(value_col)))
  ]

  if (nrow(rows) == 0) {
    return(empty_spatial_heatmap_summary(group_cols, marker1_col, marker2_col))
  }

  split_cols <- unique(c(group_cols, marker1_col, marker2_col))
  result <- rows[, .(
    mean_log2_ratio = mean(as.numeric(get(value_col)), na.rm = TRUE),
    n_detected = data.table::uniqueN(get(component_col))
  ), keyby = split_cols]
  totals <- rows[, .(
    n_total = data.table::uniqueN(get(component_col))
  ), keyby = group_cols]

  result <- totals[result, on = group_cols]
  result <- result[n_detected >= min_cells_detected]
  result[, pct_detected := ifelse(n_total > 0, n_detected / n_total, NA_real_)]
  output_cols <- c(split_cols, "mean_log2_ratio", "n_detected", "n_total", "pct_detected")
  as.data.frame(result[, ..output_cols])
}

complete_spatial_marker_pairs <- function(
  summary,
  selected_markers,
  group_cols
) {
  if (nrow(summary) == 0) {
    return(summary)
  }

  selected_markers <- as.character(selected_markers)
  require_spatial_namespace("data.table")
  summary_dt <- symmetrize_spatial_marker_pair_summary(summary, group_cols = group_cols)
  pair_grid <- data.table::CJ(
    marker_1 = selected_markers,
    marker_2 = selected_markers,
    unique = FALSE
  )
  pair_grid <- pair_grid[marker_1 != marker_2]

  group_grid <- unique(summary_dt[, ..group_cols])
  group_grid[, proxiome_cross_key := 1L]
  pair_grid[, proxiome_cross_key := 1L]
  grid <- pair_grid[group_grid, on = "proxiome_cross_key", allow.cartesian = TRUE]
  grid[, proxiome_cross_key := NULL]
  grid <- grid[, c(group_cols, "marker_1", "marker_2"), with = FALSE]

  completed <- summary_dt[grid, on = c(group_cols, "marker_1", "marker_2")]
  completed$mean_log2_ratio[is.na(completed$mean_log2_ratio)] <- 0
  completed$pct_detected[is.na(completed$pct_detected)] <- 0
  completed$n_detected[is.na(completed$n_detected)] <- 0L
  completed$n_total[is.na(completed$n_total)] <- 0L
  as.data.frame(completed)
}

symmetrize_spatial_marker_pair_summary <- function(summary, group_cols) {
  require_spatial_namespace("data.table")
  summary_dt <- data.table::as.data.table(summary)
  if (nrow(summary_dt) == 0) {
    return(summary_dt)
  }

  marker_cols <- c("marker_1", "marker_2")
  key_cols <- unique(c(group_cols, marker_cols))
  missing_cols <- setdiff(key_cols, names(summary_dt))
  if (length(missing_cols) > 0) {
    stop("Missing columns for spatial marker-pair symmetry: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  reversed <- data.table::copy(summary_dt)
  reversed_marker_1 <- reversed$marker_2
  reversed_marker_2 <- reversed$marker_1
  reversed[, marker_1 := reversed_marker_1]
  reversed[, marker_2 := reversed_marker_2]

  combined <- data.table::rbindlist(list(summary_dt, reversed), fill = TRUE)
  value_cols <- setdiff(names(combined), key_cols)
  if (length(value_cols) == 0) {
    return(unique(combined[, ..key_cols]))
  }

  combined[, lapply(.SD, average_spatial_pair_values), by = key_cols, .SDcols = value_cols]
}

average_spatial_pair_values <- function(values) {
  if (is.numeric(values)) {
    if (all(is.na(values))) {
      return(NA_real_)
    }
    return(mean(values, na.rm = TRUE))
  }

  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(NA)
  }
  values[1]
}

select_spatial_heatmap_markers <- function(
  summary,
  available_markers,
  n_markers = 20L,
  min_pct_detected = 0.25,
  min_range = 0.2
) {
  available_markers <- unique(as.character(available_markers))
  if (length(available_markers) == 0) {
    return(character(0))
  }

  n_markers <- max(2L, min(as.integer(n_markers), length(available_markers)))
  if (nrow(summary) == 0) {
    return(head(available_markers, n_markers))
  }

  markers <- intersect(unique(as.character(c(summary$marker_1, summary$marker_2))), available_markers)
  scores <- do.call(rbind, lapply(markers, function(marker) {
    values <- summary$mean_log2_ratio[
      (summary$marker_1 == marker | summary$marker_2 == marker) &
        summary$pct_detected >= min_pct_detected
    ]
    values <- as.numeric(values)
    values <- values[is.finite(values)]
    if (length(values) == 0) {
      return(data.frame(marker = marker, sd = 0, range = 0, stringsAsFactors = FALSE))
    }
    data.frame(
      marker = marker,
      sd = stats::sd(values),
      range = max(values) - min(values),
      stringsAsFactors = FALSE
    )
  }))

  scores$sd[is.na(scores$sd)] <- 0
  scores <- scores[scores$range >= min_range, , drop = FALSE]
  scores <- scores[order(scores$sd, decreasing = TRUE), , drop = FALSE]
  selected <- intersect(scores$marker, available_markers)

  if (length(selected) < 2) {
    selected <- head(available_markers, n_markers)
  }

  head(unique(selected), n_markers)
}

empty_spatial_heatmap_summary <- function(group_cols, marker1_col = "marker_1", marker2_col = "marker_2") {
  result <- data.frame(
    mean_log2_ratio = numeric(),
    n_detected = integer(),
    n_total = integer(),
    pct_detected = numeric(),
    stringsAsFactors = FALSE
  )

  for (group_col in rev(group_cols)) {
    result[[group_col]] <- character()
    result <- result[, c(group_col, setdiff(names(result), group_col)), drop = FALSE]
  }

  result[[marker1_col]] <- character()
  result[[marker2_col]] <- character()
  result
}

summarize_clustering_heatmap <- function(
  clustering,
  selected_markers,
  group_cols = c("condition", "celltype_manual"),
  value_col = "log2_ratio",
  component_col = "component",
  marker_col = "marker"
) {
  required_cols <- unique(c(group_cols, marker_col, value_col))
  missing_cols <- setdiff(required_cols, names(clustering))
  if (length(missing_cols) > 0) {
    stop("Missing columns for clustering heatmap summary: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  selected_markers <- unique(as.character(selected_markers))
  require_spatial_namespace("data.table")
  clustering_dt <- data.table::as.data.table(clustering)
  rows <- clustering_dt[
    as.character(get(marker_col)) %in% selected_markers &
      is.finite(as.numeric(get(value_col)))
  ]

  if (nrow(rows) == 0) {
    return(empty_clustering_heatmap_summary(group_cols, marker_col))
  }

  split_cols <- unique(c(marker_col, group_cols))
  if (component_col %in% names(rows)) {
    result <- rows[, .(
      mean_log2_ratio = mean(as.numeric(get(value_col)), na.rm = TRUE),
      median_log2_ratio = stats::median(as.numeric(get(value_col)), na.rm = TRUE),
      n_cells = data.table::uniqueN(get(component_col))
    ), keyby = split_cols]
  } else {
    result <- rows[, .(
      mean_log2_ratio = mean(as.numeric(get(value_col)), na.rm = TRUE),
      median_log2_ratio = stats::median(as.numeric(get(value_col)), na.rm = TRUE),
      n_cells = .N
    ), keyby = split_cols]
  }

  as.data.frame(result)
}

empty_clustering_heatmap_summary <- function(group_cols, marker_col = "marker") {
  result <- data.frame(
    marker = character(),
    mean_log2_ratio = numeric(),
    median_log2_ratio = numeric(),
    n_cells = integer(),
    stringsAsFactors = FALSE
  )
  names(result)[1] <- marker_col

  for (group_col in rev(group_cols)) {
    result[[group_col]] <- character()
    result <- result[, c(group_col, setdiff(names(result), group_col)), drop = FALSE]
  }

  result
}

select_clustering_heatmap_markers <- function(summary, n_markers = 20L, marker_col = "marker", value_col = "mean_log2_ratio") {
  if (nrow(summary) == 0 || !marker_col %in% names(summary) || !value_col %in% names(summary)) {
    return(character(0))
  }

  n_markers <- suppressWarnings(as.integer(n_markers[1]))
  if (!is.finite(n_markers) || n_markers < 1) {
    n_markers <- 20L
  }

  markers <- unique(as.character(summary[[marker_col]]))
  scores <- do.call(rbind, lapply(markers, function(marker) {
    values <- as.numeric(summary[[value_col]][as.character(summary[[marker_col]]) == marker])
    values <- values[is.finite(values)]
    data.frame(
      marker = marker,
      range = if (length(values) == 0) 0 else max(values) - min(values),
      sd = if (length(values) <= 1) 0 else stats::sd(values),
      max_abs = if (length(values) == 0) 0 else max(abs(values)),
      stringsAsFactors = FALSE
    )
  }))

  scores <- scores[order(scores$range, scores$sd, scores$max_abs, decreasing = TRUE), , drop = FALSE]
  head(scores$marker, min(n_markers, nrow(scores)))
}

require_spatial_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package is not installed: ", package, call. = FALSE)
  }
}
