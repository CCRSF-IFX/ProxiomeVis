DEFAULT_DEMO_RDS_RELATIVE <- file.path(
  "RnD_CS041188_BaoTran_XiaolinWu_3_Pixelgen_042126",
  "notebooks",
  "r",
  "pg_data_combined_fil.pixelator_v4.1.1.rds"
)

DEFAULT_DEMO_MARKERS <- c(
  "CD3e", "CD4", "CD8", "CD45", "CD81", "CD82", "CD53", "CD58",
  "CD69", "CD25", "CD279", "CD28", "CD40", "CD44", "HLA-DR", "B2M"
)

DEMO_CACHE_SCHEMA_VERSION <- 3L
ABUNDANCE_LONG_MAX_BLOCK_ENTRIES <- 500000L

default_demo_rds_path <- function(repo_root = NULL, must_work = TRUE) {
  configured_rds <- Sys.getenv("PROXIOME_DEMO_RDS", unset = "")
  if (nzchar(configured_rds)) {
    return(normalizePath(path.expand(configured_rds), mustWork = must_work))
  }

  repo_root <- repo_root %||% find_repo_root(required = must_work)
  if (is.null(repo_root)) {
    return(file.path(normalizePath(getwd(), mustWork = FALSE), DEFAULT_DEMO_RDS_RELATIVE))
  }

  file.path(repo_root, DEFAULT_DEMO_RDS_RELATIVE)
}

find_repo_root <- function(start = getwd(), required = TRUE) {
  current <- normalizePath(start, mustWork = required)

  while (TRUE) {
    candidate <- file.path(current, DEFAULT_DEMO_RDS_RELATIVE)
    if (file.exists(candidate)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      if (isTRUE(required)) {
        stop(
          "Could not find repository root containing ",
          DEFAULT_DEMO_RDS_RELATIVE,
          call. = FALSE
        )
      }
      return(NULL)
    }
    current <- parent
  }
}

proxiome_readout_context <- function() {
  data.frame(
    readout = c("Abundance", "Clustering", "Colocalization"),
    description = c(
      "Per-cell marker abundance from the PNA assay, shown on the abundance embedding and summarized by biological condition.",
      "Self-proximity for a marker within each cell graph, using marker_1 == marker_2 proximity scores to show whether molecules of the same protein cluster together.",
      "Marker-pair proximity between two different proteins, using marker_1 != marker_2 proximity scores to summarize co-organization across cell groups."
    ),
    stringsAsFactors = FALSE
  )
}

select_demo_markers <- function(
  available_markers,
  requested_markers = DEFAULT_DEMO_MARKERS,
  fallback_n = 12
) {
  selected <- requested_markers[requested_markers %in% available_markers]

  if (length(selected) == 0) {
    selected <- head(available_markers, fallback_n)
  }

  unique(selected)
}

select_proxiome_markers <- function(
  available_markers,
  requested_markers = NULL,
  marker_selection = c("demo", "all")
) {
  marker_selection <- match.arg(marker_selection)
  available_markers <- unique(as.character(available_markers))

  if (is.null(requested_markers)) {
    if (identical(marker_selection, "all")) {
      return(available_markers)
    }

    return(select_demo_markers(available_markers))
  }

  select_demo_markers(
    available_markers = available_markers,
    requested_markers = requested_markers,
    fallback_n = length(available_markers)
  )
}

load_demo_proxiome_data <- function(
  rds_path = default_demo_rds_path(must_work = FALSE),
  markers = NULL,
  marker_selection = c("demo", "all"),
  cache_path = NULL,
  force = FALSE,
  progress_callback = NULL
) {
  marker_selection <- match.arg(marker_selection)
  progress_callback <- progress_callback %||% function(...) invisible(NULL)

  if (!is.null(cache_path) && file.exists(cache_path) && !force) {
    cache <- readRDS(cache_path)
    cache_resolved_path <- cache$source$resolved_rds_path %||% cache$source$rds_path
    source_matches_cache <- identical(
      normalizePath(cache_resolved_path, mustWork = FALSE),
      normalizePath(rds_path, mustWork = FALSE)
    )
    source_is_missing <- !file.exists(rds_path)
    if (source_matches_cache || source_is_missing) {
      cache_needs_save <- !identical(cache$source$cache_schema_version, DEMO_CACHE_SCHEMA_VERSION)
      cache <- upgrade_cached_demo_data(cache)
      if (is.null(cache$qc)) {
        if (source_is_missing) {
          saveRDS(cache, cache_path)
          return(cache)
        }
        object <- readRDS(rds_path)
        cache$qc <- build_qc_data(object)
        cache$source$cache_schema_version <- DEMO_CACHE_SCHEMA_VERSION
        saveRDS(cache, cache_path)
      } else if (isTRUE(cache_needs_save)) {
        saveRDS(cache, cache_path)
      }
      return(cache)
    }
  }

  require_namespace("Seurat")
  require_namespace("SeuratObject")

  if (!file.exists(rds_path)) {
    stop(
      paste(
        "Could not find the demo RDS.",
        sprintf("Resolved path: %s", rds_path),
        "Set PROXIOME_DEMO_RDS to a readable RDS path or deploy cache/demo_proxiome_data.rds with the app.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  progress_callback("read_rds", "Reading RDS file...", 0.18)
  object <- readRDS(rds_path)
  progress_callback("build", "Building app data tables...", 0.28)
  data <- build_demo_proxiome_data(
    object,
    rds_path = rds_path,
    markers = markers,
    marker_selection = marker_selection,
    progress_callback = progress_callback
  )

  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(data, cache_path)
  }

  data
}

build_demo_proxiome_data <- function(
  object,
  rds_path,
  markers = NULL,
  marker_selection = c("demo", "all"),
  progress_callback = NULL
) {
  require_namespace("Seurat")

  progress_callback <- progress_callback %||% function(...) invisible(NULL)
  marker_selection <- match.arg(marker_selection)
  assay <- Seurat::DefaultAssay(object)
  available_markers <- rownames(object)
  markers <- select_proxiome_markers(
    available_markers = available_markers,
    requested_markers = markers,
    marker_selection = marker_selection
  )

  progress_callback("metadata", "Building metadata and QC tables...", 0.36)
  metadata <- build_embedding_table(object)
  qc <- build_qc_data(object)

  progress_callback("abundance", "Building abundance tables...", 0.48)
  abundance <- build_abundance_long(object, markers = markers, assay = assay)
  abundance_summary <- summarize_abundance(abundance, metadata)

  progress_callback("proximity", "Reading stored proximity scores...", 0.62)
  proximity <- stored_object_proximity(object, assay = assay, markers = markers)
  progress_callback("summarize", "Summarizing clustering and colocalization...", 0.88)
  readouts <- summarize_proximity_readouts(proximity, metadata)

  progress_callback("finalize", "Finalizing app data...", 0.96)
  list(
    source = list(
      rds_path = rds_path,
      resolved_rds_path = normalizePath(rds_path, mustWork = FALSE),
      assay = assay,
      n_cells = ncol(object),
      n_markers = nrow(object),
      cache_schema_version = DEMO_CACHE_SCHEMA_VERSION
    ),
    readout_context = proxiome_readout_context(),
    marker_options = markers,
    qc = qc,
    metadata = metadata,
    abundance = abundance,
    abundance_summary = abundance_summary,
    clustering = readouts$clustering,
    clustering_summary = readouts$clustering_summary,
    colocalization = readouts$colocalization,
    colocalization_summary = readouts$colocalization_summary
  )
}

stored_object_proximity <- function(object, assay = Seurat::DefaultAssay(object), markers = NULL) {
  if (!"assays" %in% slotNames(object) || !assay %in% names(object@assays)) {
    stop("Assay is not available for stored proximity extraction: ", assay, call. = FALSE)
  }

  stored_assay_proximity(object@assays[[assay]], markers = markers)
}

stored_assay_proximity <- function(assay_object, markers = NULL) {
  if (!"proximity" %in% slotNames(assay_object)) {
    stop("The selected assay does not contain a stored proximity slot.", call. = FALSE)
  }

  proximity <- methods::slot(assay_object, "proximity")
  if (is.null(proximity) || !is.data.frame(proximity) || nrow(proximity) == 0) {
    stop("The selected assay does not contain stored proximity rows.", call. = FALSE)
  }

  required_cols <- c("component", "marker_1", "marker_2", "log2_ratio")
  missing_cols <- setdiff(required_cols, names(proximity))
  if (length(missing_cols) > 0) {
    stop(
      "Stored proximity data is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(markers) && length(markers) > 0) {
    proximity <- proximity[
      proximity$marker_1 %in% markers & proximity$marker_2 %in% markers,
      ,
      drop = FALSE
    ]
  }

  proximity
}

build_qc_data <- function(object) {
  require_namespace("Seurat")

  filtered_metadata <- object[[]]
  origin_metadata <- object@misc$origin.meta.data %||% NULL
  filter_counts <- object@misc$qc_cell_counts_by_step %||% NULL

  build_qc_payload(
    filtered_metadata = filtered_metadata,
    origin_metadata = origin_metadata,
    filter_counts = filter_counts
  )
}

build_qc_payload <- function(filtered_metadata, origin_metadata = NULL, filter_counts = NULL) {
  filtered_metadata <- metadata_with_component(filtered_metadata)
  if (is.null(origin_metadata) || !is.data.frame(origin_metadata)) {
    origin_metadata <- filtered_metadata
  } else {
    origin_metadata <- metadata_with_component(origin_metadata)
  }

  filter_counts <- standardize_qc_filter_counts(
    filter_counts,
    loaded_n = nrow(origin_metadata),
    filtered_n = nrow(filtered_metadata)
  )

  list(
    origin_metadata = origin_metadata,
    filtered_metadata = filtered_metadata,
    filter_counts = filter_counts
  )
}

metadata_with_component <- function(metadata) {
  metadata <- as.data.frame(metadata)
  if (!"component" %in% names(metadata)) {
    component <- rownames(metadata)
    if (is.null(component) || any(!nzchar(component))) {
      component <- as.character(seq_len(nrow(metadata)))
    }
    metadata$component <- component
  }

  metadata <- metadata[, c("component", setdiff(names(metadata), "component")), drop = FALSE]
  rownames(metadata) <- NULL
  metadata
}

normalize_spatial_metadata <- function(metadata) {
  metadata <- as.data.frame(metadata)

  if (!"sample_alias" %in% names(metadata)) {
    sample_col <- first_matching_column(
      metadata,
      c("sample", "source_sample", "orig.ident", "orig_ident", "donor", "patient")
    )
    if (is.null(sample_col)) {
      metadata$sample_alias <- "sample"
    } else {
      metadata$sample_alias <- as.character(metadata[[sample_col]])
    }
  } else {
    metadata$sample_alias <- as.character(metadata$sample_alias)
  }
  metadata$sample_alias[is.na(metadata$sample_alias) | !nzchar(metadata$sample_alias)] <- "sample"

  if (!"celltype_manual" %in% names(metadata)) {
    metadata$celltype_manual <- "unannotated"
  }

  metadata
}

upgrade_cached_demo_data <- function(cache) {
  if (is.null(cache) || !is.list(cache)) {
    return(cache)
  }

  if (!is.null(cache$metadata) && is.data.frame(cache$metadata)) {
    cache$metadata <- normalize_spatial_metadata(cache$metadata)
  }

  cache$colocalization <- attach_spatial_metadata_to_readout(cache$colocalization, cache$metadata)
  cache$clustering <- attach_spatial_metadata_to_readout(cache$clustering, cache$metadata)

  if (is.null(cache$source) || !is.list(cache$source)) {
    cache$source <- list()
  }
  cache$source$cache_schema_version <- DEMO_CACHE_SCHEMA_VERSION

  cache
}

attach_spatial_metadata_to_readout <- function(readout, metadata) {
  if (
    is.null(readout) ||
      !is.data.frame(readout) ||
      is.null(metadata) ||
      !is.data.frame(metadata) ||
      !"component" %in% names(readout) ||
      !"component" %in% names(metadata)
  ) {
    return(readout)
  }

  metadata_cols <- intersect(c("component", "condition", "celltype_manual", "sample_alias"), names(metadata))
  metadata <- metadata[!duplicated(metadata$component), metadata_cols, drop = FALSE]
  merged <- merge(
    readout,
    metadata,
    by = "component",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("", "_metadata")
  )

  for (metadata_col in c("condition", "celltype_manual", "sample_alias")) {
    suffixed_col <- paste0(metadata_col, "_metadata")
    if (suffixed_col %in% names(merged)) {
      merged[[metadata_col]] <- merged[[metadata_col]] %||% merged[[suffixed_col]]
      merged[[suffixed_col]] <- NULL
    }
  }

  merged
}

standardize_qc_filter_counts <- function(filter_counts, loaded_n, filtered_n) {
  if (!is.null(filter_counts) && is.data.frame(filter_counts) && nrow(filter_counts) > 0) {
    filter_counts <- as.data.frame(filter_counts)
    step_col <- first_matching_column(filter_counts, c("step", "filter_step", "qc_step"))
    count_col <- first_matching_column(filter_counts, c("n_cells", "cells", "cell_count", "n", "count"))

    if (!is.null(step_col) && !is.null(count_col)) {
      names(filter_counts)[names(filter_counts) == step_col] <- "step"
      names(filter_counts)[names(filter_counts) == count_col] <- "n_cells"
      filter_counts$n_cells <- as.integer(filter_counts$n_cells)
      if (!"fraction_loaded" %in% names(filter_counts)) {
        filter_counts$fraction_loaded <- qc_fraction_loaded(filter_counts, loaded_n)
      }
      filter_counts$step_label <- format_qc_step_label(filter_counts$step)
      return(filter_counts)
    }
  }

  loaded_n <- as.integer(loaded_n)
  filtered_n <- as.integer(filtered_n)
  data.frame(
    step = c("00_loaded", "99_filtered"),
    n_cells = c(loaded_n, filtered_n),
    fraction_loaded = c(1, ifelse(loaded_n > 0, filtered_n / loaded_n, NA_real_)),
    step_label = c("Loaded", "Filtered"),
    stringsAsFactors = FALSE
  )
}

first_matching_column <- function(data, candidates) {
  matches <- candidates[candidates %in% names(data)]
  if (length(matches) == 0) {
    return(NULL)
  }
  matches[1]
}

qc_fraction_loaded <- function(filter_counts, loaded_n) {
  if ("sample" %in% names(filter_counts)) {
    loaded_by_sample <- tapply(
      filter_counts$n_cells,
      filter_counts$sample,
      function(values) values[1]
    )
    loaded <- loaded_by_sample[as.character(filter_counts$sample)]
  } else {
    loaded <- rep(filter_counts$n_cells[1] %||% loaded_n, nrow(filter_counts))
  }

  loaded <- as.numeric(loaded)
  ifelse(loaded > 0, filter_counts$n_cells / loaded, NA_real_)
}

format_qc_step_label <- function(step) {
  label <- sub("^\\d+_", "", as.character(step))
  label <- gsub("_", " ", label)
  paste0(toupper(substr(label, 1, 1)), substring(label, 2))
}

available_qc_distribution_choices <- function(metadata) {
  candidates <- c(
    "UMIs" = "n_umi",
    "Edges" = "n_edges",
    "Reads in component" = "reads_in_component",
    "Isotype fraction" = "isotype_fraction",
    "Tau" = "tau"
  )

  available <- candidates[candidates %in% names(metadata)]
  available[vapply(available, function(column) is.numeric(metadata[[column]]), logical(1))]
}

build_embedding_table <- function(object) {
  require_namespace("Seurat")

  metadata <- object[[]]
  metadata$component <- rownames(metadata)

  for (reduction in Seurat::Reductions(object)) {
    embedding <- as.data.frame(Seurat::Embeddings(object, reduction))
    embedding$component <- rownames(embedding)
    names(embedding) <- c(
      paste0(reduction, "_", seq_len(ncol(embedding) - 1)),
      "component"
    )
    metadata <- merge(metadata, embedding, by = "component", all.x = TRUE, sort = FALSE)
  }

  if (!"celltype_manual" %in% names(metadata)) {
    metadata$celltype_manual <- "unannotated"
  }
  if (!"condition" %in% names(metadata)) {
    metadata$condition <- "unknown"
  }
  if (!"seurat_cluster" %in% names(metadata) && "seurat_clusters" %in% names(metadata)) {
    metadata$seurat_cluster <- metadata$seurat_clusters
  }
  metadata <- normalize_spatial_metadata(metadata)

  metadata
}

build_abundance_long <- function(object, markers, assay = Seurat::DefaultAssay(object)) {
  require_namespace("Seurat")
  require_namespace("SeuratObject")

  data_layer <- SeuratObject::LayerData(object, assay = assay, layer = "data")
  counts_layer <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  markers <- markers[markers %in% rownames(data_layer)]

  matrix_layers_to_long(
    data_layer[markers, , drop = FALSE],
    counts_layer[markers, , drop = FALSE],
  )
}

matrix_layers_to_long <- function(
  data_layer,
  counts_layer,
  max_block_entries = ABUNDANCE_LONG_MAX_BLOCK_ENTRIES
) {
  layer_dim <- dim(data_layer)
  if (!identical(layer_dim, dim(counts_layer))) {
    stop("Abundance data and count layers must have the same dimensions.", call. = FALSE)
  }

  n_markers <- layer_dim[1]
  n_components <- layer_dim[2]
  markers <- rownames(data_layer) %||% as.character(seq_len(n_markers))
  components <- colnames(data_layer) %||% as.character(seq_len(n_components))
  total_values <- n_markers * n_components
  abundance <- numeric(total_values)
  count <- numeric(total_values)

  for (columns in matrix_column_chunks(n_components, n_markers, max_block_entries)) {
    value_index <- ((columns[1] - 1L) * n_markers + 1L):(tail(columns, 1) * n_markers)
    abundance[value_index] <- matrix_block_values(
      data_layer[, columns, drop = FALSE],
      n_rows = n_markers,
      n_columns = length(columns)
    )
    count[value_index] <- matrix_block_values(
      counts_layer[, columns, drop = FALSE],
      n_rows = n_markers,
      n_columns = length(columns)
    )
  }

  data.frame(
    marker = rep(markers, times = n_components),
    component = rep(components, each = n_markers),
    abundance = abundance,
    count = count,
    stringsAsFactors = FALSE
  )
}

matrix_to_long <- function(
  x,
  value_name,
  max_block_entries = ABUNDANCE_LONG_MAX_BLOCK_ENTRIES
) {
  x_dim <- dim(x)
  n_markers <- x_dim[1]
  n_components <- x_dim[2]
  markers <- rownames(x) %||% as.character(seq_len(n_markers))
  components <- colnames(x) %||% as.character(seq_len(n_components))
  values <- numeric(n_markers * n_components)

  for (columns in matrix_column_chunks(n_components, n_markers, max_block_entries)) {
    value_index <- ((columns[1] - 1L) * n_markers + 1L):(tail(columns, 1) * n_markers)
    values[value_index] <- matrix_block_values(
      x[, columns, drop = FALSE],
      n_rows = n_markers,
      n_columns = length(columns)
    )
  }

  result <- data.frame(
    marker = rep(markers, times = n_components),
    component = rep(components, each = n_markers),
    value = values,
    stringsAsFactors = FALSE
  )
  names(result) <- c("marker", "component", value_name)
  result
}

matrix_column_chunks <- function(
  n_columns,
  n_rows,
  max_block_entries = ABUNDANCE_LONG_MAX_BLOCK_ENTRIES
) {
  if (n_columns == 0 || n_rows == 0) {
    return(list())
  }

  max_block_entries <- suppressWarnings(as.integer(max_block_entries[1]))
  if (!is.finite(max_block_entries) || max_block_entries < 1) {
    max_block_entries <- ABUNDANCE_LONG_MAX_BLOCK_ENTRIES
  }

  columns_per_chunk <- max(1L, floor(max_block_entries / max(1L, n_rows)))
  split(seq_len(n_columns), ceiling(seq_len(n_columns) / columns_per_chunk))
}

matrix_block_values <- function(block, n_rows = nrow(block), n_columns = ncol(block)) {
  if (inherits(block, "sparseMatrix")) {
    require_namespace("Matrix")
    values <- numeric(n_rows * n_columns)
    triplet <- methods::as(block, "TsparseMatrix")
    if (length(triplet@x) > 0) {
      values[(triplet@i + 1L) + triplet@j * n_rows] <- triplet@x
    }
    return(values)
  }

  as.numeric(block)
}

summarize_abundance <- function(abundance, metadata) {
  joined <- join_metadata_by_component(
    abundance,
    metadata,
    metadata_cols = c("component", "condition", "celltype_manual")
  )

  aggregate_numeric_readout(
    joined,
    group_cols = c("marker", "condition", "celltype_manual"),
    value_col = "abundance"
  )
}

summarize_proximity_readouts <- function(proximity, metadata) {
  proximity <- join_metadata_by_component(
    proximity,
    metadata,
    metadata_cols = c("component", "condition", "celltype_manual", "sample_alias")
  )

  for (metadata_col in c("condition_metadata", "celltype_manual_metadata", "sample_alias_metadata")) {
    if (metadata_col %in% names(proximity)) {
      proximity[[metadata_col]] <- NULL
    }
  }

  clustering <- proximity[
    proximity$marker_1 == proximity$marker_2,
    ,
    drop = FALSE
  ]
  clustering$marker <- clustering$marker_1

  colocalization <- proximity[
    proximity$marker_1 != proximity$marker_2,
    ,
    drop = FALSE
  ]
  colocalization$marker_pair <- paste(colocalization$marker_1, "/", colocalization$marker_2)

  list(
    clustering = as.data.frame(clustering),
    clustering_summary = aggregate_numeric_readout(
      clustering,
      group_cols = c("marker", "condition", "celltype_manual"),
      value_col = "log2_ratio"
    ),
    colocalization = as.data.frame(colocalization),
    colocalization_summary = aggregate_numeric_readout(
      colocalization,
      group_cols = c("marker_pair", "marker_1", "marker_2", "condition", "celltype_manual"),
      value_col = "log2_ratio"
    )
  )
}

join_metadata_by_component <- function(data, metadata, metadata_cols) {
  require_namespace("data.table")
  data_dt <- data.table::as.data.table(data)
  metadata_cols <- intersect(metadata_cols, names(metadata))
  if (!"component" %in% names(data_dt) || !"component" %in% metadata_cols) {
    return(data_dt)
  }

  metadata_dt <- data.table::as.data.table(metadata[, metadata_cols, drop = FALSE])
  conflicting_cols <- setdiff(intersect(names(metadata_dt), names(data_dt)), "component")
  if (length(conflicting_cols) > 0) {
    data.table::setnames(metadata_dt, conflicting_cols, paste0(conflicting_cols, "_metadata"))
  }

  metadata_dt[data_dt, on = "component"]
}

aggregate_numeric_readout <- function(data, group_cols, value_col) {
  require_namespace("data.table")
  group_cols <- group_cols[group_cols %in% names(data)]

  if (nrow(data) == 0) {
    return(data.frame())
  }

  data_dt <- data.table::as.data.table(data)
  if ("component" %in% names(data_dt)) {
    result <- data_dt[, .(
      mean_value = mean(as.numeric(get(value_col)), na.rm = TRUE),
      median_value = stats::median(as.numeric(get(value_col)), na.rm = TRUE),
      n_cells = data.table::uniqueN(get("component"))
    ), keyby = group_cols]
  } else {
    result <- data_dt[, .(
      mean_value = mean(as.numeric(get(value_col)), na.rm = TRUE),
      median_value = stats::median(as.numeric(get(value_col)), na.rm = TRUE),
      n_cells = .N
    ), keyby = group_cols]
  }

  as.data.frame(result)
}

calculate_differential_readout <- function(
  data,
  feature_cols,
  value_col,
  group_a,
  group_b,
  condition_col = "condition",
  celltype_col = "celltype_manual",
  celltype_filter = NULL,
  stratify_by_celltype = FALSE,
  pooled_celltype_label = "Pooled cell types",
  min_cells = 3,
  fdr_cutoff = 0.05
) {
  required_cols <- unique(c(feature_cols, value_col, condition_col))
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns for differential readout: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (!celltype_col %in% names(data)) {
    data[[celltype_col]] <- "All cells"
  }

  output_cols <- c(
    feature_cols,
    celltype_col,
    "group_a", "group_b",
    "mean_a", "mean_b", "median_a", "median_b",
    "n_a", "n_b",
    "effect_size", "p_value", "p_adj",
    "direction", "is_significant"
  )

  if (nrow(data) == 0 || is.null(group_a) || is.null(group_b) || identical(group_a, group_b)) {
    return(empty_differential_readout(output_cols))
  }

  data <- data[data[[condition_col]] %in% c(group_a, group_b), , drop = FALSE]
  if (!is.null(celltype_filter) && length(celltype_filter) > 0) {
    data <- data[data[[celltype_col]] %in% celltype_filter, , drop = FALSE]
  }

  if (nrow(data) == 0) {
    return(empty_differential_readout(output_cols))
  }

  if (!isTRUE(stratify_by_celltype)) {
    selected_celltypes <- unique(as.character(data[[celltype_col]]))
    data[[celltype_col]] <- if (length(selected_celltypes) == 1) {
      selected_celltypes
    } else {
      pooled_celltype_label
    }
  }

  split_cols <- unique(c(feature_cols, celltype_col))
  row_groups <- split(
    seq_len(nrow(data)),
    interaction(data[split_cols], drop = TRUE, lex.order = TRUE)
  )

  results <- lapply(row_groups, function(row_index) {
    chunk <- data[row_index, , drop = FALSE]
    values_a <- finite_values(chunk[[value_col]][chunk[[condition_col]] == group_a])
    values_b <- finite_values(chunk[[value_col]][chunk[[condition_col]] == group_b])

    feature_row <- chunk[1, split_cols, drop = FALSE]
    feature_row$group_a <- group_a
    feature_row$group_b <- group_b
    feature_row$mean_a <- mean_or_na(values_a)
    feature_row$mean_b <- mean_or_na(values_b)
    feature_row$median_a <- median_or_na(values_a)
    feature_row$median_b <- median_or_na(values_b)
    feature_row$n_a <- length(values_a)
    feature_row$n_b <- length(values_b)
    feature_row$effect_size <- feature_row$median_a - feature_row$median_b
    feature_row$p_value <- differential_p_value(values_a, values_b, min_cells = min_cells)
    feature_row
  })

  result <- do.call(rbind, results)
  rownames(result) <- NULL

  result$p_adj <- NA_real_
  valid_p <- !is.na(result$p_value)
  if (any(valid_p)) {
    result$p_adj[valid_p] <- stats::p.adjust(result$p_value[valid_p], method = "BH")
  }

  result$direction <- ifelse(
    is.na(result$effect_size) | result$effect_size == 0,
    "No change",
    ifelse(result$effect_size > 0, paste("Higher in", group_a), paste("Higher in", group_b))
  )
  result$is_significant <- !is.na(result$p_adj) & result$p_adj <= fdr_cutoff
  result[output_cols]
}

empty_differential_readout <- function(columns) {
  result <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  names(result) <- columns
  result
}

finite_values <- function(values) {
  values <- as.numeric(values)
  values[is.finite(values)]
}

mean_or_na <- function(values) {
  if (length(values) == 0) {
    return(NA_real_)
  }
  mean(values)
}

median_or_na <- function(values) {
  if (length(values) == 0) {
    return(NA_real_)
  }
  stats::median(values)
}

differential_p_value <- function(values_a, values_b, min_cells) {
  if (length(values_a) < min_cells || length(values_b) < min_cells) {
    return(NA_real_)
  }
  if (length(unique(c(values_a, values_b))) < 2) {
    return(NA_real_)
  }

  tryCatch(
    stats::wilcox.test(values_b, values_a, exact = FALSE)$p.value,
    error = function(error) NA_real_
  )
}

require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package is not installed: ", package, call. = FALSE)
  }
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
