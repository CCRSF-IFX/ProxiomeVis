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
  rows <- proximity[
    proximity[[marker1_col]] %in% selected_markers &
      proximity[[marker2_col]] %in% selected_markers &
      proximity[[marker1_col]] != proximity[[marker2_col]] &
      is.finite(as.numeric(proximity[[value_col]])),
    ,
    drop = FALSE
  ]

  if (nrow(rows) == 0) {
    return(empty_spatial_heatmap_summary(group_cols, marker1_col, marker2_col))
  }

  split_cols <- unique(c(group_cols, marker1_col, marker2_col))
  row_groups <- split(
    seq_len(nrow(rows)),
    interaction(rows[split_cols], drop = TRUE, lex.order = TRUE)
  )

  result <- do.call(rbind, lapply(row_groups, function(row_index) {
    chunk <- rows[row_index, , drop = FALSE]
    out <- chunk[1, split_cols, drop = FALSE]
    out$mean_log2_ratio <- mean(as.numeric(chunk[[value_col]]), na.rm = TRUE)
    out$n_detected <- length(unique(chunk[[component_col]]))
    out
  }))
  rownames(result) <- NULL

  total_groups <- split(
    seq_len(nrow(rows)),
    interaction(rows[group_cols], drop = TRUE, lex.order = TRUE)
  )
  totals <- do.call(rbind, lapply(total_groups, function(row_index) {
    chunk <- rows[row_index, , drop = FALSE]
    out <- chunk[1, group_cols, drop = FALSE]
    out$n_total <- length(unique(chunk[[component_col]]))
    out
  }))
  rownames(totals) <- NULL

  result <- merge(result, totals, by = group_cols, all.x = TRUE, sort = FALSE)
  result <- result[result$n_detected >= min_cells_detected, , drop = FALSE]
  result$pct_detected <- ifelse(result$n_total > 0, result$n_detected / result$n_total, NA_real_)
  result
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
  pair_grid <- expand.grid(
    marker_1 = selected_markers,
    marker_2 = selected_markers,
    stringsAsFactors = FALSE
  )
  pair_grid <- pair_grid[pair_grid$marker_1 != pair_grid$marker_2, , drop = FALSE]

  group_grid <- unique(summary[group_cols])
  grid <- merge(group_grid, pair_grid, all = TRUE)

  completed <- merge(
    grid,
    summary,
    by = c(group_cols, "marker_1", "marker_2"),
    all.x = TRUE,
    sort = FALSE
  )
  completed$mean_log2_ratio[is.na(completed$mean_log2_ratio)] <- 0
  completed$pct_detected[is.na(completed$pct_detected)] <- 0
  completed$n_detected[is.na(completed$n_detected)] <- 0L
  completed$n_total[is.na(completed$n_total)] <- 0L
  completed
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
  rows <- clustering[
    as.character(clustering[[marker_col]]) %in% selected_markers &
      is.finite(as.numeric(clustering[[value_col]])),
    ,
    drop = FALSE
  ]

  if (nrow(rows) == 0) {
    return(empty_clustering_heatmap_summary(group_cols, marker_col))
  }

  split_cols <- unique(c(marker_col, group_cols))
  row_groups <- split(
    seq_len(nrow(rows)),
    interaction(rows[split_cols], drop = TRUE, lex.order = TRUE)
  )

  result <- do.call(rbind, lapply(row_groups, function(row_index) {
    chunk <- rows[row_index, , drop = FALSE]
    out <- chunk[1, split_cols, drop = FALSE]
    out$mean_log2_ratio <- mean(as.numeric(chunk[[value_col]]), na.rm = TRUE)
    out$median_log2_ratio <- stats::median(as.numeric(chunk[[value_col]]), na.rm = TRUE)
    out$n_cells <- if (component_col %in% names(chunk)) {
      length(unique(chunk[[component_col]]))
    } else {
      nrow(chunk)
    }
    out
  }))
  rownames(result) <- NULL
  result
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
