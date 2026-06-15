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

pixelator_raw_component_id <- function(component, sample) {
  component <- as.character(component)[1]
  sample <- as.character(sample)[1]
  if (is.na(component)) {
    return("")
  }
  if (is.na(sample) || !nzchar(sample)) {
    return(component)
  }

  prefix <- paste0(sample, "_")
  if (startsWith(component, prefix)) {
    return(substr(component, nchar(prefix) + 1L, nchar(component)))
  }

  component
}

pixelator_layout_pxl_path <- function(
  sample,
  source = NULL,
  layout_dir = Sys.getenv("PROXIOME_LAYOUT_DIR", unset = "")
) {
  sample <- as.character(sample)[1]
  if (is.na(sample) || !nzchar(sample)) {
    return("")
  }

  filename <- paste0(sample, ".layout.pxl")
  layout_dir <- as.character(layout_dir)[1]
  if (!is.na(layout_dir) && nzchar(layout_dir)) {
    candidate <- file.path(path.expand(layout_dir), filename)
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  for (root in pixelator_source_roots(source)) {
    for (relative_dir in pixelator_layout_relative_dirs()) {
      candidate <- file.path(root, relative_dir, filename)
      if (file.exists(candidate)) {
        return(normalizePath(candidate, mustWork = TRUE))
      }
    }
  }

  ""
}

pixelator_source_roots <- function(source) {
  paths <- character()
  if (!is.null(source)) {
    paths <- unlist(source[c("resolved_rds_path", "rds_path")], use.names = FALSE)
  }
  paths <- unique(as.character(paths))
  paths <- paths[!is.na(paths) & nzchar(paths)]

  roots <- character()
  for (path in paths) {
    parts <- strsplit(normalizePath(path, mustWork = FALSE), .Platform$file.sep, fixed = TRUE)[[1]]
    notebook_pos <- which(parts == "notebooks")
    if (length(notebook_pos) > 0 && notebook_pos[1] > 1) {
      root <- paste(parts[seq_len(notebook_pos[1] - 1L)], collapse = .Platform$file.sep)
      if (startsWith(path, .Platform$file.sep)) {
        root <- paste0(.Platform$file.sep, sub(paste0("^", .Platform$file.sep, "+"), "", root))
      }
      roots <- c(roots, root)
    }
    roots <- c(roots, dirname(dirname(dirname(path))))
  }

  unique(roots[nzchar(roots)])
}

pixelator_layout_relative_dirs <- function() {
  c(
    file.path("results", "run_pixelator-4.1.1_merged_pixelator_v0.27.2", "pixelator"),
    file.path("results", "nf-core_pixelator", "pixelator"),
    file.path("results", "run_pixelator-3.0.1_merged", "nf-core_pixelator_merged", "pixelator"),
    file.path("results", "run_pixelator-3.0.1", "nf-core_pixelator", "pixelator")
  )
}

read_pixelator_3d_layout <- function(pxl_path, component, layout_method = "wpmds_3d") {
  require_spatial_namespace("DBI")
  require_spatial_namespace("duckdb")

  pxl_path <- as.character(pxl_path)[1]
  component <- as.character(component)[1]
  if (is.na(pxl_path) || !file.exists(pxl_path)) {
    stop("Pixelator layout file does not exist: ", pxl_path, call. = FALSE)
  }
  if (is.na(component) || !nzchar(component)) {
    stop("Pixelator component id is required.", call. = FALSE)
  }

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = pxl_path, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  sql <- '
    with node_labels as (
      select cast(umi1 as varchar) as node_id, marker_1 as marker
      from edgelist
      where component = ?
      union all
      select cast(umi2 as varchar) as node_id, marker_2 as marker
      from edgelist
      where component = ?
    ),
    collapsed_labels as (
      select node_id, min(marker) as marker
      from node_labels
      where marker is not null
      group by node_id
    )
    select
      cast(l."index" as varchar) as node_id,
      cast(l.x as double) as x,
      cast(l.y as double) as y,
      cast(l.z as double) as z,
      coalesce(c.marker, \'unlabeled\') as marker
    from layouts l
    left join collapsed_labels c on cast(l."index" as varchar) = c.node_id
    where l.component = ? and l.layout = ?
  '

  DBI::dbGetQuery(con, sql, params = list(component, component, component, layout_method))
}

prepare_pixelator_3d_layout <- function(
  layout,
  highlighted_markers,
  max_background_nodes = 7000L,
  seed = 1L
) {
  required_cols <- c("node_id", "x", "y", "z", "marker")
  missing_cols <- setdiff(required_cols, names(layout))
  if (length(missing_cols) > 0) {
    stop("Missing columns for 3D layout: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  nodes <- layout[is.finite(layout$x) & is.finite(layout$y) & is.finite(layout$z), required_cols, drop = FALSE]
  if (nrow(nodes) == 0) {
    return(data.frame(
      node_id = character(),
      x = numeric(),
      y = numeric(),
      z = numeric(),
      marker = character(),
      x_scaled = numeric(),
      y_scaled = numeric(),
      z_scaled = numeric(),
      marker_group = character(),
      stringsAsFactors = FALSE
    ))
  }

  highlighted_markers <- unique(as.character(highlighted_markers))
  highlighted_markers <- highlighted_markers[!is.na(highlighted_markers) & nzchar(highlighted_markers)]
  nodes$marker <- as.character(nodes$marker)
  nodes$marker[is.na(nodes$marker) | !nzchar(nodes$marker)] <- "unlabeled"
  nodes$x_scaled <- scale_pixelator_layout_axis(nodes$x)
  nodes$y_scaled <- scale_pixelator_layout_axis(nodes$y)
  nodes$z_scaled <- scale_pixelator_layout_axis(nodes$z)
  nodes$marker_group <- ifelse(nodes$marker %in% highlighted_markers, nodes$marker, "Other")

  background <- nodes$marker_group == "Other"
  max_background_nodes <- suppressWarnings(as.integer(max_background_nodes[1]))
  if (!is.finite(max_background_nodes) || max_background_nodes < 0) {
    max_background_nodes <- 7000L
  }
  background_index <- which(background)
  if (length(background_index) > max_background_nodes) {
    keep_background <- background_index[round(seq(1, length(background_index), length.out = max_background_nodes))]
    nodes <- nodes[sort(c(which(!background), keep_background)), , drop = FALSE]
  }

  nodes
}

scale_pixelator_layout_axis <- function(values) {
  values <- as.numeric(values)
  values <- values - mean(values, na.rm = TRUE)
  denom <- max(abs(values), na.rm = TRUE)
  if (!is.finite(denom) || denom == 0) {
    return(rep(0, length(values)))
  }
  values / denom
}

pixelator_3d_layout_plot <- function(
  nodes,
  highlighted_markers,
  title = NULL,
  dimensions = list(width = 832, height = 620)
) {
  require_spatial_namespace("plotly")

  highlighted_markers <- unique(as.character(highlighted_markers))
  highlighted_markers <- highlighted_markers[!is.na(highlighted_markers) & nzchar(highlighted_markers)]
  groups <- c("Other", intersect(highlighted_markers, unique(nodes$marker_group)))
  colors <- pixelator_3d_marker_colors()
  plot <- plotly::plot_ly(width = dimensions$width, height = dimensions$height)

  for (group in groups) {
    rows <- nodes[nodes$marker_group == group, , drop = FALSE]
    if (nrow(rows) == 0) {
      next
    }
    plot <- plotly::add_trace(
      plot,
      data = rows,
      x = ~x_scaled,
      y = ~y_scaled,
      z = ~z_scaled,
      type = "scatter3d",
      mode = "markers",
      name = group,
      text = ~paste0("Node: ", node_id, "<br>Marker: ", marker),
      hoverinfo = "text",
      marker = list(
        size = if (identical(group, "Other")) 2 else 4,
        opacity = if (identical(group, "Other")) 0.12 else 0.9,
        color = if (group %in% names(colors)) colors[[group]] else "#6c757d"
      ),
      inherit = FALSE
    )
  }

  plotly::layout(
    plot,
    title = title,
    margin = list(l = 0, r = 0, t = 40, b = 0),
    scene = list(
      xaxis = pixelator_3d_axis(),
      yaxis = pixelator_3d_axis(),
      zaxis = pixelator_3d_axis()
    ),
    legend = list(orientation = "h", x = 0, y = 1.05)
  )
}

pixelator_3d_axis <- function() {
  list(title = "", showgrid = FALSE, nticks = 0, showticklabels = FALSE, zeroline = FALSE)
}

pixelator_3d_marker_colors <- function() {
  c(
    "CD54" = "#D62728",
    "ICAM-1" = "#D62728",
    "CD40" = "#54A24B",
    "CD8" = "#F58518",
    "CD3e" = "#4C78A8",
    "CD81" = "#9C755F",
    "CD82" = "#B279A2",
    "Other" = "#9aa3a6"
  )
}

require_spatial_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package is not installed: ", package, call. = FALSE)
  }
}
