cell_graph_error <- function(message) stop(message, call. = FALSE)

cell_graph_schema <- function(version = 3L) {
  if (identical(version, 2L)) return(arrow::schema(
    cell_id = arrow::utf8(), node_id = arrow::utf8(), x = arrow::float64(),
    y = arrow::float64(), z = arrow::float64(), markers = arrow::list_of(arrow::utf8())))
  arrow::schema(cell_id = arrow::utf8(), x = arrow::float64(),
                y = arrow::float64(), z = arrow::float64(), markers = arrow::list_of(arrow::utf8()))
}

cell_graph_path <- function(asset_dir, file) {
  root <- normalizePath(asset_dir, mustWork = TRUE)
  if (length(file) != 1L || is.na(file) ||
      !grepl("^[a-f0-9]{32}\\.parquet$", file)) cell_graph_error("Invalid cell graph asset filename.")
  path <- normalizePath(file.path(root, file), mustWork = TRUE)
  if (!startsWith(path, paste0(root, .Platform$file.sep))) {
    cell_graph_error("Cell graph asset is outside the configured directory.")
  }
  path
}

read_cell_graph_index <- function(asset_dir, dataset_id, cell_ids = NULL) {
  root <- normalizePath(asset_dir, mustWork = TRUE)
  path <- normalizePath(file.path(root, "index.rds"), mustWork = TRUE)
  if (!startsWith(path, paste0(root, .Platform$file.sep))) cell_graph_error("Cell graph index is outside its directory.")
  index <- readRDS(path)
  if (!is.list(index) || !(identical(index$schema_version, 2L) || identical(index$schema_version, 3L))) {
    cell_graph_error("Unsupported cell graph schema. Rebuild sample Parquet assets into a new directory; legacy RDS assets are not modified.")
  }
  cell_fields <- c("cell_id", "source_cell_id", "sample", "n_nodes")
  sample_fields <- c("sample", "file", "md5", "size_bytes", "source_pxl", "source_pxl_md5", "layout", "n_nodes", "n_cells")
  if (!identical(index$dataset_id, dataset_id) || !is.data.frame(index$cells) ||
      !is.data.frame(index$samples) || !all(cell_fields %in% names(index$cells)) ||
      !all(sample_fields %in% names(index$samples))) {
    cell_graph_error("Cell graph index does not match this dataset or schema.")
  }
  cells <- index$cells
  samples <- index$samples
  positive_integers <- function(value) is.numeric(value) && all(is.finite(value) & value > 0 & value == floor(value))
  strings <- function(table, fields) all(vapply(table[fields], function(value) {
    is.character(value) && !anyNA(value) && all(nzchar(value))
  }, logical(1)))
  if (!nrow(cells) || !nrow(samples) || !strings(cells, cell_fields[1:3]) ||
      !strings(samples, setdiff(sample_fields, c("n_nodes", "n_cells", "size_bytes"))) ||
      anyDuplicated(cells$cell_id) || anyDuplicated(samples$sample) || anyDuplicated(samples$file) ||
      !positive_integers(cells$n_nodes) || !positive_integers(samples$n_nodes) ||
      !positive_integers(samples$n_cells) || !positive_integers(samples$size_bytes) ||
      any(!grepl("^[a-f0-9]{32}$", samples$md5)) ||
      any(samples$file != paste0(samples$md5, ".parquet")) ||
      any(!grepl("^[a-f0-9]{32}$", samples$source_pxl_md5)) ||
      !setequal(cells$sample, samples$sample) ||
      (!is.null(cell_ids) && !all(cells$cell_id %in% cell_ids)) ||
      !positive_integers(index$n_total_cells) || length(index$n_total_cells) != 1L ||
      nrow(cells) > index$n_total_cells ||
      !identical(index$coverage, if (nrow(cells) == index$n_total_cells) "all_cells" else "selected_cells")) {
    cell_graph_error("Cell graph index contains invalid or unknown cells.")
  }
  for (i in seq_len(nrow(samples))) {
    selected <- cells[cells$sample == samples$sample[i], , drop = FALSE]
    if (nrow(selected) != samples$n_cells[i] || sum(selected$n_nodes) != samples$n_nodes[i]) {
      cell_graph_error("Cell graph index sample totals do not match its cells.")
    }
  }
  # Cache is scoped to this loaded index, never global; changed files are verified again.
  attr(index, "verified_samples") <- new.env(parent = emptyenv())
  index
}

validate_cell_graph_asset <- function(asset) {
  fields <- c("cell_id", "x", "y", "z", "markers")
  if (!is.list(asset) || !identical(asset$schema_version, 3L) ||
      !is.character(asset$cell_id) || length(asset$cell_id) != 1L || is.na(asset$cell_id) || !nzchar(asset$cell_id) ||
      !is.data.frame(asset$nodes) || !identical(names(asset$nodes), fields) ||
      !nrow(asset$nodes) ||
      !is.character(asset$nodes$cell_id) || anyNA(asset$nodes$cell_id) ||
      !all(asset$nodes$cell_id == asset$cell_id) ||
      !all(vapply(asset$nodes[c("x", "y", "z")], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(asset$nodes[c("x", "y", "z")]))) ||
      !is.list(asset$nodes$markers) ||
      !all(vapply(asset$nodes$markers, function(markers) {
        is.character(markers) && !anyNA(markers) && all(nzchar(markers)) && !anyDuplicated(markers)
      }, logical(1)))) {
    cell_graph_error("Cell graph asset has invalid coordinates, cell identities, or marker identities.")
  }
  invisible(asset)
}

cell_graph_verified_file <- function(index, asset_dir, sample) {
  file <- cell_graph_path(asset_dir, sample$file)
  observed <- file.info(file)
  stamp <- list(path = file, md5 = sample$md5, size = observed$size, mtime = observed$mtime, ctime = observed$ctime)
  cache <- attr(index, "verified_samples")
  cached <- if (is.environment(cache)) cache[[sample$file]] else NULL
  if (!identical(cached, stamp)) {
    if (observed$size != sample$size_bytes || !identical(unname(tools::md5sum(file)), sample$md5)) {
      cell_graph_error("Cell graph sample checksum mismatch.")
    }
    if (is.environment(cache)) cache[[sample$file]] <- stamp
  }
  file
}

load_cell_graph_asset <- function(index, asset_dir, cell_id, max_nodes = 200000L) {
  if (!requireNamespace("arrow", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE)) {
    cell_graph_error("Install arrow and dplyr to read prepared sample Parquet files.")
  }
  if (length(cell_id) != 1L || is.na(cell_id)) cell_graph_error("Select one prepared cell.")
  if (length(max_nodes) != 1L || !is.numeric(max_nodes) || !is.finite(max_nodes) ||
      max_nodes < 1 || max_nodes != floor(max_nodes)) cell_graph_error("The viewer node limit must be a positive integer.")
  entry <- index$cells[index$cells$cell_id == cell_id, , drop = FALSE]
  if (nrow(entry) != 1L) cell_graph_error("No prepared 3D asset is available for this cell.")
  limit_error <- function() cell_graph_error("This cell exceeds the configured 3D viewer size limit; no nodes were silently removed.")
  if (entry$n_nodes > max_nodes) limit_error()
  sample <- index$samples[index$samples$sample == entry$sample, , drop = FALSE]
  if (nrow(sample) != 1L) cell_graph_error("Cell graph sample provenance is missing.")
  file <- cell_graph_verified_file(index, asset_dir, sample)
  dataset <- arrow::open_dataset(file, format = "parquet")
  if (!dataset$schema$Equals(cell_graph_schema(index$schema_version), check_metadata = FALSE)) {
    cell_graph_error("Sample Parquet does not match its graph schema version; version 3 requires only cell_id, x, y, z, and list<string> markers.")
  }
  selected_cell <- cell_id
  query <- dplyr::filter(dataset, .data$cell_id == .env$selected_cell)
  query <- dplyr::select(query, dplyr::all_of(cell_graph_schema()$names))
  nodes <- as.data.frame(dplyr::collect(utils::head(query, max_nodes + 1)))
  nodes$markers <- lapply(nodes$markers, as.character)
  if (nrow(nodes) > max_nodes) limit_error()
  asset <- list(schema_version = 3L, dataset_id = index$dataset_id, cell_id = cell_id,
                source_cell_id = entry$source_cell_id, sample = entry$sample,
                source_pxl = sample$source_pxl, source_pxl_md5 = sample$source_pxl_md5,
                layout = sample$layout, graph_projection = sample$graph_projection, nodes = nodes,
                interpretation = "Abstract graph coordinates, not physical cell positions; all source layout nodes retained.")
  validate_cell_graph_asset(asset)
  if (nrow(nodes) != entry$n_nodes) cell_graph_error("Cell graph row count does not match its index provenance.")
  asset
}

plot_cell_graph_asset <- function(asset, markers = character()) {
  validate_cell_graph_asset(asset)
  if (length(markers) > 3L || anyNA(markers) || anyDuplicated(markers) ||
      !all(markers %in% unlist(asset$nodes$markers, use.names = FALSE))) {
    cell_graph_error("Choose up to three markers present in the prepared cell.")
  }
  nodes <- asset$nodes
  code <- rep(0L, nrow(nodes))
  if (length(markers)) {
    for (i in seq_along(markers)) {
      code <- code + 2^(i - 1L) * vapply(nodes$markers, function(value) markers[i] %in% value, logical(1))
    }
  }
  colors <- c("#DEDEDE", "#D73027", "#2166AC", "#984EA3", "#1B9E77", "#E69F00", "#00BFC4", "#111111")
  node_markers <- vapply(nodes$markers, function(value) {
    if (!length(value)) "None" else paste(htmltools::htmlEscape(value), collapse = "<br>")
  }, character(1))
  hover <- paste0("Markers: ", node_markers)
  coordinates <- as.matrix(nodes[c("x", "y", "z")])
  ranges <- apply(coordinates, 2, range)
  span <- max(ranges[2, ] - ranges[1, ])
  # One common scale preserves geometry and also handles flat/single-node layouts.
  coordinates <- sweep(coordinates, 2, colMeans(ranges)) / if (span > 0) span else 1
  plot <- threejs::scatterplot3js(coordinates, color = paste0(colors[code + 1L], "D9"),
    labels = hover, axis = FALSE, grid = FALSE, size = 0.1, center = TRUE,
    xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5), zlim = c(-0.5, 0.5),
    main = htmltools::htmlEscape(asset$cell_id), font.main = "14px sans-serif")
  plot$elementId <- NULL
  plot
}

prepare_cell_graph_assets <- function(object, pxl_dir, output_dir, cells = colnames(object),
                                      dataset_id = object@misc$proxiomevis$dataset_id) {
  if (!requireNamespace("pixelatorR", quietly = TRUE) || !requireNamespace("arrow", quietly = TRUE)) {
    cell_graph_error("Install pixelatorR and arrow to prepare sample Parquet assets.")
  }
  if (!inherits(object, "Seurat") || length(dataset_id) != 1L || is.na(dataset_id) || !nzchar(dataset_id) ||
      !length(cells) || anyNA(cells) || anyDuplicated(cells) || !all(cells %in% colnames(object))) {
    cell_graph_error("Preparation requires a source Seurat object, dataset ID, and known unique cells.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  pxl_dir <- normalizePath(pxl_dir, mustWork = TRUE)
  if (identical(output_dir, pxl_dir) || startsWith(output_dir, paste0(pxl_dir, .Platform$file.sep))) {
    cell_graph_error("Write graph assets outside the source PXL directory.")
  }
  index <- if (file.exists(file.path(output_dir, "index.rds"))) {
    read_cell_graph_index(output_dir, dataset_id, colnames(object))
  } else list(schema_version = 3L, dataset_id = dataset_id, cells = NULL, samples = NULL)
  if (!identical(index$schema_version, 3L)) {
    cell_graph_error("Existing version 2 graph assets are read-only. Prepare version 3 assets into a new directory; existing files are not modified.")
  }
  mapping <- pixelatorR::FSMap(object)
  defaults <- object@misc$proxiomevis$defaults
  sample_column <- if (is.null(defaults$sample_column)) "sample_alias" else defaults$sample_column
  if (!sample_column %in% names(object@meta.data)) cell_graph_error("Source object is missing its sample metadata column.")
  manifest <- object@misc$proxiomevis$proximity_source$pxl_manifest
  if (!is.null(manifest) && (!is.data.frame(manifest) ||
      !all(c("sample_alias", "pxl_md5") %in% names(manifest)))) {
    cell_graph_error("Source PXL manifest must contain sample_alias and pxl_md5.")
  }
  sql_string <- function(value) paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
  seen_samples <- character()
  for (source in seq_len(nrow(mapping))) {
    ids <- as.data.frame(mapping$id_map[[source]])
    ids <- ids[ids$current_id %in% cells, , drop = FALSE]
    if (!nrow(ids)) next
    if (anyNA(ids) || anyDuplicated(ids$current_id) || anyDuplicated(ids$original_id)) {
      cell_graph_error("Source cell mappings must be unique and complete.")
    }
    sample <- unique(as.character(object@meta.data[ids$current_id, sample_column]))
    if (length(sample) != 1L || is.na(sample) || !nzchar(sample) || sample %in% seen_samples) {
      cell_graph_error("Each sample must map to exactly one PXL source.")
    }
    seen_samples <- c(seen_samples, sample)
    pxl_file <- normalizePath(file.path(pxl_dir, basename(mapping$pxl_file[[source]])), mustWork = TRUE)
    if (!startsWith(pxl_file, paste0(pxl_dir, .Platform$file.sep))) cell_graph_error("PXL file is outside the configured directory.")
    fingerprint <- unname(tools::md5sum(pxl_file))
    if (!is.null(manifest)) {
      expected <- as.character(manifest$pxl_md5[manifest$sample_alias == sample])
      if (length(expected) != 1L || !identical(expected, fingerprint)) {
        cell_graph_error("PXL checksum does not match the frozen analysis source manifest.")
      }
    }
    previous <- if (is.null(index$samples)) NULL else index$samples[index$samples$sample == sample, , drop = FALSE]
    previous_cells <- if (is.null(index$cells)) NULL else index$cells[index$cells$sample == sample, , drop = FALSE]
    old_file <- NULL
    if (NROW(previous)) {
      all_ids <- as.data.frame(mapping$id_map[[source]])
      if (!identical(previous$source_pxl_md5, fingerprint) ||
          !identical(previous$source_pxl, basename(pxl_file)) ||
          !identical(as.character(previous_cells$source_cell_id),
                     as.character(all_ids$original_id[match(previous_cells$cell_id, all_ids$current_id)]))) {
        cell_graph_error("Existing sample provenance changed; rebuild into a new directory.")
      }
      old_file <- cell_graph_verified_file(index, output_dir, previous)
      ids <- ids[!ids$current_id %in% previous_cells$cell_id, , drop = FALSE]
      if (!nrow(ids)) {
        message("Reusing sample Parquet: ", sample)
        next
      }
    }
    db <- pixelatorR::PixelDB$new(pxl_file)
    temporary <- tempfile("sample-", tmpdir = output_dir, fileext = ".parquet")
    writer <- sink <- NULL
    finished <- FALSE
    new_cells <- NULL
    tryCatch({
      suppressWarnings(db$query("SET threads = 1"))
      suppressWarnings(db$query("SET memory_limit = '256MB'"))
      suppressWarnings(db$query(paste0("SET temp_directory = ", sql_string(tempdir()))))
      layouts <- db$query("SELECT DISTINCT layout, graph_projection FROM layouts")
      if (nrow(layouts) != 1L || anyNA(layouts) || any(!nzchar(as.matrix(layouts)))) {
        cell_graph_error("Source must contain one unambiguous layout and graph projection.")
      }
      layout <- as.character(layouts$layout)
      projection <- as.character(layouts$graph_projection)
      if (NROW(previous) && (!identical(previous$layout, layout) ||
                            !identical(previous$graph_projection, projection))) {
        cell_graph_error("Existing sample layout changed; rebuild into a new directory.")
      }
      schema <- cell_graph_schema()
      sink <- arrow::FileOutputStream$create(temporary)
      writer <- arrow::ParquetFileWriter$create(schema, sink,
        arrow::ParquetWriterProperties$create(schema$names, compression = "zstd", write_statistics = TRUE))
      if (!is.null(old_file)) {
        reader <- arrow::ParquetFileReader$create(old_file)
        if (!reader$GetSchema()$Equals(schema, check_metadata = FALSE)) cell_graph_error("Existing sample Parquet schema is invalid.")
        for (group in seq_len(reader$num_row_groups)) {
          table <- reader$ReadRowGroup(group - 1L)
          writer$WriteTable(table, chunk_size = table$num_rows)
          if (group %% 25L == 0L) invisible(gc())
        }
        rm(reader, table)
      }
      for (i in seq_len(nrow(ids))) {
        cell <- as.character(ids$current_id[i])
        id <- sql_string(as.character(ids$original_id[i]))
        message("Preparing sample ", sample, " cell ", i, "/", nrow(ids), ": ", cell)
        query <- paste0(
          "WITH e AS MATERIALIZED (SELECT umi1, umi2, marker_1, marker_2 FROM edgelist WHERE component = ", id,
          "), endpoints AS (SELECT CAST(umi1 AS VARCHAR) || '-umi1' AS node_id, marker_1 AS marker FROM e ",
          "UNION ALL SELECT CAST(umi2 AS VARCHAR) || '-umi2', marker_2 FROM e",
          "), markers AS (SELECT node_id, list(DISTINCT marker ORDER BY marker) AS markers FROM endpoints GROUP BY node_id) ",
          "SELECT ", sql_string(cell), " AS cell_id, CAST(l.index AS VARCHAR) || CASE WHEN l.pixel_type IN ('A','umi1') THEN '-umi1' ",
          "WHEN l.pixel_type IN ('B','umi2') THEN '-umi2' ELSE NULL END AS node_id, ",
          "CAST(l.x AS DOUBLE) AS x, CAST(l.y AS DOUBLE) AS y, CAST(l.z AS DOUBLE) AS z, ",
          "coalesce(m.markers, []::VARCHAR[]) AS markers FROM layouts l LEFT JOIN markers m ON m.node_id = ",
          "CAST(l.index AS VARCHAR) || CASE WHEN l.pixel_type IN ('A','umi1') THEN '-umi1' ELSE '-umi2' END ",
          "WHERE l.component = ", id, " AND l.layout = ", sql_string(layout),
          " AND l.graph_projection = ", sql_string(projection), " ORDER BY node_id")
        nodes <- as.data.frame(db$query(query))
        # Source IDs are needed for marker joins and uniqueness checks, not display.
        if (!is.character(nodes$node_id) || anyNA(nodes$node_id) ||
            anyDuplicated(nodes$node_id) || any(!nzchar(nodes$node_id))) {
          cell_graph_error("Source layout contains invalid or duplicate node identities.")
        }
        nodes$node_id <- NULL
        validate_cell_graph_asset(list(schema_version = 3L, cell_id = cell, nodes = nodes))
        table <- arrow::Table$create(nodes, schema = schema)
        writer$WriteTable(table, chunk_size = nrow(nodes))
        new_cells <- rbind(new_cells, data.frame(cell_id = cell, source_cell_id = as.character(ids$original_id[i]),
                                                sample = sample, n_nodes = nrow(nodes)))
        rm(nodes, table)
        if (i %% 25L == 0L) invisible(gc())
      }
      writer$Close()
      writer <- NULL
      sink$close()
      sink <- NULL
      finished <- TRUE
    }, finally = {
      if (!is.null(writer)) try(writer$Close(), silent = TRUE)
      if (!is.null(sink)) try(sink$close(), silent = TRUE)
      db$close()
      if (!finished && file.exists(temporary)) unlink(temporary)
    })
    if (!identical(unname(tools::md5sum(pxl_file)), fingerprint)) {
      unlink(temporary)
      cell_graph_error("Source PXL changed during preparation; existing published assets were not modified.")
    }
    checksum <- unname(tools::md5sum(temporary))
    filename <- paste0(checksum, ".parquet")
    destination <- file.path(output_dir, filename)
    if (file.exists(destination)) {
      if (!identical(unname(tools::md5sum(destination)), checksum)) cell_graph_error("Existing sample output checksum mismatch.")
      unlink(temporary)
    } else if (!file.rename(temporary, destination)) cell_graph_error("Could not finalize sample Parquet.")
    retained_cells <- if (is.null(index$cells)) NULL else index$cells[index$cells$sample != sample, , drop = FALSE]
    retained_samples <- if (is.null(index$samples)) NULL else index$samples[index$samples$sample != sample, , drop = FALSE]
    sample_cells <- rbind(previous_cells, new_cells)
    index$cells <- rbind(retained_cells, sample_cells)
    index$samples <- rbind(retained_samples, data.frame(sample = sample, file = filename, md5 = checksum,
      size_bytes = file.info(destination)$size, source_pxl = basename(pxl_file), source_pxl_md5 = fingerprint,
      layout = layout, graph_projection = projection, n_nodes = sum(sample_cells$n_nodes), n_cells = nrow(sample_cells)))
    index$n_total_cells <- ncol(object)
    index$coverage <- if (nrow(index$cells) == ncol(object)) "all_cells" else "selected_cells"
    index$created_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    attr(index, "verified_samples") <- NULL
    temporary_index <- tempfile("index-", tmpdir = output_dir, fileext = ".rds")
    saveRDS(index, temporary_index)
    if (!file.rename(temporary_index, file.path(output_dir, "index.rds"))) cell_graph_error("Could not finalize graph index.")
    # Only superseded output is removed, after its replacement index is durable.
    if (!is.null(old_file)) unlink(old_file)
    index <- read_cell_graph_index(output_dir, dataset_id, colnames(object))
  }
  if (!all(cells %in% index$cells$cell_id)) cell_graph_error("Not all selected cells have graph mappings.")
  invisible(read_cell_graph_index(output_dir, dataset_id, colnames(object)))
}
