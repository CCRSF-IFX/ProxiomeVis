hybrid_dataset_fingerprint <- function(data) {
  digest::digest(list(dataset_id = data$contract$dataset_id,
                       metadata = data$metadata, counts = data$counts, cell_sources = data$cell_sources,
                       filter_settings = lapply(data$contract$report$settings, function(settings) {
                         settings$spatial$low_abundance_filter
                       })), algo = "sha256")
}

resolve_hybrid_asset <- function(bundle_path, descriptor) {
  if (!is.list(descriptor) || !is.character(descriptor$path) || length(descriptor$path) != 1L ||
      is.na(descriptor$path) || !nzchar(descriptor$path) ||
      grepl("(^[/\\\\])|(^|[/\\\\])\\.\\.([/\\\\]|$)|^[A-Za-z]:", descriptor$path)) {
    stop("This bundle has no valid relative path for the requested prepared asset.", call. = FALSE)
  }
  root <- normalizePath(dirname(bundle_path), mustWork = TRUE)
  target <- normalizePath(file.path(root, descriptor$path), mustWork = TRUE)
  if (!dir.exists(target) || !startsWith(target, paste0(root, "/"))) {
    stop("Prepared asset paths must stay inside the bundle directory.", call. = FALSE)
  }
  target
}

proximity_store_file <- function(root, relative) {
  if (length(relative) != 1L || is.na(relative) || !grepl("^[A-Za-z0-9_.-]+[.]parquet$", relative)) {
    stop("Invalid proximity partition path.", call. = FALSE)
  }
  root <- normalizePath(root, mustWork = TRUE)
  path <- normalizePath(file.path(root, relative), mustWork = TRUE)
  if (!startsWith(path, paste0(root, "/")) || dir.exists(path)) stop("Invalid proximity partition path.", call. = FALSE)
  path
}

prepare_proximity_store <- function(data, pxl_by_sample, output_dir) {
  metadata <- data$metadata
  mapping <- data$cell_sources
  if (is.null(data$counts) || !is.data.frame(mapping) ||
      !all(c("component", "original_id", "sample_id") %in% names(mapping)) ||
      anyNA(mapping) || anyDuplicated(mapping$component) || anyDuplicated(mapping[c("sample_id", "original_id")]) ||
      !setequal(mapping$component, metadata$component) ||
      !all(mapping$sample_id == metadata$sample_id[match(mapping$component, metadata$component)])) {
    stop("Preparation requires aligned counts and exact per-sample original-to-current cell IDs.", call. = FALSE)
  }
  samples <- sort(unique(metadata$sample_id))
  if (!is.character(pxl_by_sample) || is.null(names(pxl_by_sample)) || anyDuplicated(names(pxl_by_sample)) ||
      !all(samples %in% names(pxl_by_sample))) stop("Supply a PXL path for every sample.", call. = FALSE)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("The output directory must be empty; existing assets are not overwritten.", call. = FALSE)
  }
  donor_cutoffs <- vapply(sort(unique(metadata$donor)), function(donor) {
    settings <- data$contract$report$settings[[donor]]$spatial$low_abundance_filter
    if (!is.null(settings) && (!isTRUE(settings$enabled) || !identical(settings$fraction_mode, "isotype") ||
        !is.null(settings$minimum_marker_count) || !is.null(settings$minimum_cells_per_pair))) {
      stop("This preparation supports the report's isotype-fraction filter without additional count filters.", call. = FALSE)
    }
    column <- if (is.null(settings$isotype_column)) "isotype_fraction" else settings$isotype_column
    if (!column %in% names(metadata)) stop("The isotype-fraction column is required for spatial preparation.", call. = FALSE)
    values <- metadata[[column]][metadata$donor == donor]
    values <- values[is.finite(values) & values >= 0]
    if (!length(values)) stop("Missing finite donor isotype fractions.", call. = FALSE)
    multiplier <- if (is.null(settings$mad_multiplier)) 2 else settings$mad_multiplier
    floor <- if (is.null(settings$minimum_fraction)) 0.001 else settings$minimum_fraction
    cutoff <- max(stats::median(values) + multiplier * stats::mad(values), floor)
    if (!is.finite(cutoff) || cutoff < 0 || cutoff > 1) stop("Invalid donor isotype cutoff.", call. = FALSE)
    cutoff
  }, numeric(1))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  partitions <- lapply(seq_along(samples), function(i) {
    sample <- samples[[i]]
    source <- normalizePath(pxl_by_sample[[sample]], mustWork = TRUE)
    source_md5 <- unname(tools::md5sum(source))
    inventory <- data$contract$proximity_source$pxl_manifest
    if (!is.null(inventory)) {
      expected <- inventory$pxl_md5[match(sample, inventory$sample_alias)]
      if (length(expected) != 1L || is.na(expected) || !identical(source_md5, as.character(expected))) {
        stop("A source PXL does not match the frozen report provenance.", call. = FALSE)
      }
    }
    ids <- mapping[mapping$sample_id == sample, , drop = FALSE]
    counts <- data$counts[, ids$component, drop = FALSE]
    totals <- Matrix::colSums(counts)
    if (any(!is.finite(totals) | totals <= 0)) stop("Counts must have positive finite cell totals.", call. = FALSE)
    long_counts <- data.frame(component = rep(ids$original_id, each = nrow(counts)),
                              marker = rep(rownames(counts), times = ncol(counts)),
                              fraction = as.numeric(counts) / rep(totals, each = nrow(counts)))
    spill_dir <- tempfile("proxiomevis-proximity-")
    dir.create(spill_dir)
    on.exit(unlink(spill_dir, recursive = TRUE), add = TRUE)
    connection <- DBI::dbConnect(duckdb::duckdb(), dbdir = source, read_only = TRUE,
                                  config = list(threads = "1", memory_limit = "512MB", temp_directory = spill_dir))
    on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
    DBI::dbWriteTable(connection, "cell_map", ids, temporary = TRUE)
    DBI::dbWriteTable(connection, "marker_fractions", long_counts, temporary = TRUE)
    query <- paste(
      "SELECT m.component, least(p.marker_1,p.marker_2) AS marker_1, greatest(p.marker_1,p.marker_2) AS marker_2,",
      "CASE WHEN p.join_count IS NULL OR p.join_count_expected_mean IS NULL THEN NULL ELSE",
      "log2(greatest(p.join_count,1)/greatest(p.join_count_expected_mean,1)) END AS log2_ratio,",
      "least(c1.fraction,c2.fraction) AS min_marker_fraction FROM proximity p",
      "JOIN cell_map m ON p.component=m.original_id",
      "JOIN marker_fractions c1 ON p.component=c1.component AND p.marker_1=c1.marker",
      "JOIN marker_fractions c2 ON p.component=c2.component AND p.marker_2=c2.marker"
    )
    filename <- sprintf("sample-%03d.parquet", i)
    path <- file.path(output_dir, filename)
    DBI::dbExecute(connection, paste0("COPY (", query, ") TO ", DBI::dbQuoteString(connection, path),
                                      " (FORMAT PARQUET, COMPRESSION ZSTD)"))
    path_sql <- DBI::dbQuoteString(connection, path)
    duplicate <- DBI::dbGetQuery(connection, paste0("SELECT 1 FROM read_parquet(", path_sql,
      ") GROUP BY component, marker_1, marker_2 HAVING count(*) > 1 LIMIT 1"))
    if (nrow(duplicate)) stop("Duplicate cell/marker-pair rows in PXL input.", call. = FALSE)
    nrows <- DBI::dbGetQuery(connection, paste0("SELECT count(*) AS n FROM read_parquet(", path_sql, ")"))$n
    if (!identical(source_md5, unname(tools::md5sum(source)))) {
      stop("A source PXL changed during proximity preparation.", call. = FALSE)
    }
    data.frame(sample_id = sample, file = filename, rows = as.numeric(nrows),
               md5 = unname(tools::md5sum(path)), source_pxl = source, source_pxl_md5 = source_md5)
  })
  manifest <- list(schema_version = 1L, dataset_id = data$contract$dataset_id,
                   dataset_fingerprint = hybrid_dataset_fingerprint(data),
                   created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                   partitions = do.call(rbind, partitions), donor_cutoffs = donor_cutoffs,
                   metric = "log2(pmax(join_count,1)/pmax(join_count_expected_mean,1))",
                   missing_score_policy = "observed_only; no imputation; finite scores required at testing")
  saveRDS(manifest, file.path(output_dir, "manifest.rds"))
  invisible(manifest)
}

validated_proximity_partitions <- function(data, asset_root, cells) {
  if (!is.character(cells) || !length(cells) || anyNA(cells) || anyDuplicated(cells) ||
      !all(cells %in% data$metadata$component)) {
    stop("Choose known, unique cells for proximity lookup.", call. = FALSE)
  }
  asset_root <- normalizePath(asset_root, mustWork = TRUE)
  manifest_path <- file.path(asset_root, "manifest.rds")
  if (!file.exists(manifest_path)) stop("Prepare the proximity sidecar before running spatial tests.", call. = FALSE)
  manifest_path <- normalizePath(manifest_path, mustWork = TRUE)
  if (!startsWith(manifest_path, paste0(asset_root, "/"))) stop("Invalid proximity manifest path.", call. = FALSE)
  manifest <- readRDS(manifest_path)
  if (!identical(manifest$schema_version, 1L) || !identical(manifest$dataset_id, data$contract$dataset_id) ||
      !identical(manifest$dataset_fingerprint, hybrid_dataset_fingerprint(data))) {
    stop("The proximity sidecar does not match this dataset.", call. = FALSE)
  }
  samples <- unique(data$metadata$sample_id[match(cells, data$metadata$component)])
  partitions <- manifest$partitions
  if (!is.data.frame(partitions) || !all(c("sample_id", "file", "md5") %in% names(partitions)) ||
      anyDuplicated(partitions$sample_id) || !all(samples %in% partitions$sample_id)) {
    stop("Incomplete proximity partition manifest.", call. = FALSE)
  }
  partitions <- partitions[match(samples, partitions$sample_id), , drop = FALSE]
  inventory <- data$contract$proximity_source$pxl_manifest
  if (!is.null(inventory)) {
    if (!is.data.frame(inventory) || !all(c("sample_alias", "pxl_md5") %in% names(inventory)) ||
        !"source_pxl_md5" %in% names(partitions)) {
      stop("Missing source PXL checksums for the frozen analysis or proximity sidecar.", call. = FALSE)
    }
    expected <- inventory[inventory$sample_alias %in% samples, , drop = FALSE]
    if (anyDuplicated(expected$sample_alias) || !setequal(expected$sample_alias, samples)) {
      stop("Frozen source PXL provenance must identify each requested sample exactly once.", call. = FALSE)
    }
    expected_hashes <- as.character(expected$pxl_md5[match(partitions$sample_id, expected$sample_alias)])
    if (anyNA(expected_hashes) || any(!nzchar(expected_hashes)) ||
        !identical(as.character(partitions$source_pxl_md5), expected_hashes)) {
      stop("The proximity sidecar source PXL checksums do not match the frozen analysis.", call. = FALSE)
    }
  }
  paths <- vapply(partitions$file, function(file) proximity_store_file(asset_root, file), character(1))
  if (!identical(unname(tools::md5sum(paths)), as.character(partitions$md5))) {
    stop("A proximity partition changed; prepare a matching sidecar again.", call. = FALSE)
  }
  list(paths = unname(paths), manifest = manifest, manifest_path = manifest_path)
}

read_proximity_subset <- function(data, asset_root, request) {
  request <- validate_hybrid_request(data, request)
  cells <- c(request$group_a_cells, request$group_b_cells)
  prepared <- validated_proximity_partitions(data, asset_root, cells)
  cutoff <- prepared$manifest$donor_cutoffs[[request$donor]]
  if (is.null(cutoff) || length(cutoff) != 1L || !is.finite(cutoff) || cutoff < 0 || cutoff > 1) {
    stop("Missing donor filtering cutoff in proximity manifest.", call. = FALSE)
  }
  markers <- request$markers
  query <- arrow::open_dataset(prepared$paths)
  query <- dplyr::filter(query, component %in% cells, marker_1 %in% markers, marker_2 %in% markers,
                         min_marker_fraction >= cutoff, is.finite(log2_ratio))
  query <- if (request$readout == "self_proximity") {
    dplyr::filter(query, marker_1 == marker_2)
  } else {
    dplyr::filter(query, marker_1 != marker_2)
  }
  scores <- as.data.frame(dplyr::collect(utils::head(query, HYBRID_LIMITS$rows + 1L)))
  if (nrow(scores) > HYBRID_LIMITS$rows) stop("Too many proximity rows; select fewer cells or markers.", call. = FALSE)
  if (anyDuplicated(scores[c("component", "marker_1", "marker_2")])) stop("Duplicate proximity observations.", call. = FALSE)
  attr(scores, "provenance") <- list(manifest_md5 = unname(tools::md5sum(prepared$manifest_path)),
    observed_rows = nrow(scores), donor_isotype_fraction_cutoff = cutoff,
    filtering = "Both markers pass a donor-wide pinned isotype threshold; only finite observed log2 ratios are tested.")
  scores
}

validate_graph_markers <- function(markers, available) {
  if (!is.character(markers) || length(markers) < 2L || anyNA(markers) ||
      any(!nzchar(trimws(markers))) || anyDuplicated(markers)) {
    stop("graph_markers must contain at least two distinct, nonempty marker names.", call. = FALSE)
  }
  unknown <- setdiff(markers, available)
  if (length(unknown)) stop("Unknown graph markers: ", paste(unknown, collapse = ", "), call. = FALSE)
  sort(markers, method = "radix")
}

select_marker_graph_cells <- function(data, asset_root, markers) {
  markers <- validate_graph_markers(markers, rownames(data$abundance))
  prepared <- validated_proximity_partitions(data, asset_root, data$metadata$component)
  spill <- tempfile("proxiomevis-graph-ranking-")
  dir.create(spill)
  on.exit(unlink(spill, recursive = TRUE), add = TRUE)
  connection <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1",
    memory_limit = "512MB", temp_directory = spill))
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE, after = FALSE)
  DBI::dbWriteTable(connection, "metadata", data$metadata[c("component", "donor", "celltype", "condition")])
  paths <- paste(DBI::dbQuoteString(connection, prepared$paths), collapse = ",")
  marker_sql <- paste(DBI::dbQuoteString(connection, markers), collapse = ",")
  # Rank on disk; only the selected cell union enters R memory.
  DBI::dbExecute(connection, paste0("CREATE VIEW scores AS SELECT p.component, m.donor, m.celltype, m.condition, ",
    "least(p.marker_1,p.marker_2) AS marker_1, greatest(p.marker_1,p.marker_2) AS marker_2, p.log2_ratio ",
    "FROM read_parquet([", paths, "]) p JOIN metadata m USING (component) ",
    "WHERE p.marker_1 IN (", marker_sql, ") AND p.marker_2 IN (", marker_sql, ") ",
    "AND p.marker_1 <> p.marker_2"))
  duplicates <- DBI::dbGetQuery(connection, paste(
    "SELECT 1 FROM scores GROUP BY component, marker_1, marker_2 HAVING count(*) > 1 LIMIT 1"))
  if (nrow(duplicates)) stop("Duplicate proximity observations for an unordered marker pair.", call. = FALSE)
  message("Ranking graph cells: top 2 and bottom 2 per different-marker pair among ",
    paste(markers, collapse = ", "), ", within donor × cell type × condition.")
  selected <- DBI::dbGetQuery(connection, paste0(
    "WITH ranked AS (SELECT *, ",
    "row_number() OVER (PARTITION BY donor, celltype, condition, marker_1, marker_2 ",
    "ORDER BY log2_ratio DESC, component ASC) AS top_rank, ",
    "row_number() OVER (PARTITION BY donor, celltype, condition, marker_1, marker_2 ",
    "ORDER BY log2_ratio ASC, component ASC) AS bottom_rank ",
    "FROM scores WHERE isfinite(log2_ratio)) ",
    "SELECT component, bool_or(top_rank <= 2) AS colocalization_top, ",
    "bool_or(bottom_rank <= 2) AS colocalization_bottom ",
    "FROM ranked WHERE top_rank <= 2 OR bottom_rank <= 2 ",
    "GROUP BY component ORDER BY component"))
  if (!nrow(selected)) stop("No finite different-marker proximity scores are available for the selected graph markers.", call. = FALSE)
  selected <- merge(selected, data$metadata[c("component", "sample_id", "donor", "celltype", "condition")],
    by = "component", sort = TRUE)
  summary <- do.call(rbind, lapply(sort(unique(data$metadata$sample_id)), function(sample) {
    rows <- selected[selected$sample_id == sample, , drop = FALSE]
    data.frame(sample_id = sample, total_cells = sum(data$metadata$sample_id == sample),
      colocalization_top_cells = sum(rows$colocalization_top),
      colocalization_bottom_cells = sum(rows$colocalization_bottom), selected_cells = nrow(rows))
  }))
  list(cells = selected$component, selected = selected, summary = summary,
    parameters = list(markers = markers, strata = c("donor", "celltype", "condition"),
      colocalization = "Top 2 and bottom 2 per unordered different-marker pair; self-pairs excluded",
      metric = "log2_ratio", ties = "Exact cell ID ascending in each tail; at most 2 per tail",
      filtering = "Finite observed scores only; no imputation or abundance cutoff",
      proximity_manifest_md5 = unname(tools::md5sum(prepared$manifest_path))))
}

read_proximity_feature <- function(data, asset_root, cells, pair) {
  if (!requireNamespace("arrow", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install arrow and dplyr to read cell-level proximity scores.", call. = FALSE)
  }
  if (!is.character(pair) || length(pair) != 2L || anyNA(pair) ||
      !all(pair %in% rownames(data$abundance))) stop("Choose known markers.", call. = FALSE)
  prepared <- validated_proximity_partitions(data, asset_root, cells)
  pair_a <- pair[1L]
  pair_b <- pair[2L]
  query <- arrow::open_dataset(prepared$paths)
  # Match both stored orientations; never synthesize a missing score as zero.
  query <- dplyr::filter(query, .data$component %in% .env$cells, is.finite(.data$log2_ratio),
    (.data$marker_1 == .env$pair_a & .data$marker_2 == .env$pair_b) |
      (.data$marker_1 == .env$pair_b & .data$marker_2 == .env$pair_a))
  scores <- as.data.frame(dplyr::collect(utils::head(query, length(cells) + 1L)))
  if (nrow(scores) > length(cells) || anyDuplicated(scores$component)) {
    stop("Duplicate proximity observations for the selected unordered marker pair.", call. = FALSE)
  }
  attr(scores, "manifest_md5") <- unname(tools::md5sum(prepared$manifest_path))
  scores
}

rank_cell_proximity <- function(data, asset_root, donor, celltype, condition, marker_1, marker_2) {
  for (value in list(donor, celltype, condition, marker_1, marker_2)) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("Choose one donor, cell type, condition, and one marker in each pair selector.", call. = FALSE)
    }
  }
  pair <- sort(c(marker_1, marker_2), method = "radix")
  meta <- data$metadata
  cells <- meta$component[meta$donor == donor & meta$celltype == celltype & meta$condition == condition]
  if (!length(cells)) stop("No cells match this donor, cell type, and condition.", call. = FALSE)
  scores <- read_proximity_feature(data, asset_root, cells, pair)
  manifest_md5 <- attr(scores, "manifest_md5")
  scores <- scores[order(-scores$log2_ratio, scores$component, method = "radix"), , drop = FALSE]
  result <- data.frame(rank = seq_len(nrow(scores)), component = scores$component,
    sample_id = meta$sample_id[match(scores$component, meta$component)],
    donor = rep(donor, nrow(scores)), celltype = rep(celltype, nrow(scores)),
    condition = rep(condition, nrow(scores)), marker_1 = rep(pair[1L], nrow(scores)),
    marker_2 = rep(pair[2L], nrow(scores)), log2_ratio = scores$log2_ratio)
  attr(result, "provenance") <- list(manifest_md5 = manifest_md5,
    cohort_cells = length(cells), scored_cells = nrow(scores),
    ranking = "Descending observed log2_ratio; ties ordered by exact cell ID.",
    filtering = "Finite observed scores only; no imputation, abundance cutoff, or differential-test group-size threshold.")
  result
}
