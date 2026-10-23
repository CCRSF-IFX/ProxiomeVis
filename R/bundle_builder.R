# Source precomputed_contract.R before this file. Preparation runs offline.
empty_user_results <- function() {
  feature <- data.frame(celltype = character(), feature_id = character(), marker_1 = character(),
    marker_2 = character(), metric = character(), is_control = logical())
  comparison <- data.frame(analysis_id = character(), group_a = character(), group_b = character(),
    comparison = character())
  tables <- list(
    sample_summaries = cbind(feature, data.frame(sample_id = character(), donor = character(),
      condition = character(), value = numeric(), n_total_cells = numeric(),
      n_evaluable_cells = numeric(), n_scored_cells = numeric())),
    donor_contrasts = cbind(feature, comparison, data.frame(donor = character(), value_a = numeric(),
      value_b = numeric(), effect_size = numeric(), n_cells_a = numeric(), n_cells_b = numeric())),
    differential_results = cbind(feature, comparison, data.frame(effect_size = numeric(),
      min_effect = numeric(), max_effect = numeric(), n_pairs = numeric(), direction_consistent = logical(),
      p_value = numeric(), p_adj = numeric(), adjustment_family = character(), method = character(),
      replicate_unit = character(), inference = character()))
  )
  lapply(tables, function(table) stats::setNames(rep(list(table), length(PROXIOMEVIS_READOUTS)), PROXIOMEVIS_READOUTS))
}

user_bundle_data <- function(object, title = "My proximity dataset", dataset_id = NULL,
                             assay = "PNA", abundance_layer = "data", reduction = "umap", pca_reduction = "pca",
                             celltype_column = "celltype", cluster_column = "seurat_clusters",
                             donor_column = "donor", condition_column = "condition", sample_column = "sample_id",
                             saved_results = NULL, qc = NULL, cell_sources = NULL, control_markers = character()) {
  if (!inherits(object, "Seurat")) stop("The input must be an annotated Seurat object.", call. = FALSE)
  defaults <- list(reduction = reduction, pca_reduction = pca_reduction, celltype_column = celltype_column,
    cluster_column = cluster_column, donor_column = donor_column, condition_column = condition_column,
    sample_column = sample_column, assay = assay, abundance_layer = abundance_layer)
  validate_precomputed_header(list(schema_version = PROXIOMEVIS_SCHEMA_VERSION,
    dataset_id = if (is.null(dataset_id)) "pending" else dataset_id, title = title, defaults = defaults))
  metadata <- object@meta.data
  if (!identical(rownames(metadata), colnames(object))) stop("Seurat cell identifiers are not aligned.", call. = FALSE)
  if (!assay %in% names(object@assays)) stop("The configured assay does not exist: ", assay, call. = FALSE)
  layers <- SeuratObject::Layers(object[[assay]])
  if (!abundance_layer %in% layers) {
    stop("The configured normalized layer does not exist. Normalize the assay before preparation.", call. = FALSE)
  }
  abundance <- SeuratObject::LayerData(object[[assay]], layer = abundance_layer)
  counts <- if ("counts" %in% layers) SeuratObject::LayerData(object[[assay]], layer = "counts") else NULL
  contract <- c(list(schema_version = PROXIOMEVIS_SCHEMA_VERSION,
    dataset_id = if (is.null(dataset_id)) paste0("user-", substr(digest::digest(list(defaults, metadata, abundance),
      algo = "sha256"), 1L, 16L)) else dataset_id, title = title, defaults = defaults), empty_user_results())
  contract$qc <- list(cell_counts = data.frame())
  contract$provenance <- list(builder = "prepare_user_bundle", created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    abundance_layer = abundance_layer, normalization = "Stored values; no normalization or differential test is run.")
  if (is.null(cell_sources)) {
    assay_object <- object[[assay]]
    if ("fs_map" %in% methods::slotNames(assay_object) && nrow(assay_object@fs_map)) {
      cell_sources <- do.call(rbind, lapply(assay_object@fs_map$id_map, function(mapping) {
        data.frame(component = as.character(mapping$current_id), original_id = as.character(mapping$original_id))
      }))
      cell_sources <- cell_sources[cell_sources$component %in% rownames(metadata), , drop = FALSE]
      cell_sources$sample_id <- as.character(metadata[[sample_column]][match(cell_sources$component, rownames(metadata))])
      rownames(cell_sources) <- NULL
    }
  }
  data <- list(contract = contract, metadata = metadata,
    embeddings = lapply(object@reductions, function(x) x@cell.embeddings), abundance = abundance,
    counts = counts, cell_sources = cell_sources)
  validate_precomputed_data(data)
  # Normalize aliases through the same reader used by the app, including overlapping column names.
  normalized <- precomputed_app_data(c(list(bundle_version = PROXIOMEVIS_BUNDLE_VERSION), data))
  meta <- normalized$metadata
  if (!is.character(control_markers) || anyNA(control_markers) || anyDuplicated(control_markers) ||
      !all(control_markers %in% rownames(abundance))) stop("Choose known, unique control markers.", call. = FALSE)
  groups <- unique(meta[c("sample_id", "donor", "condition", "celltype")])
  profiles <- lapply(seq_len(nrow(groups)), function(i) {
    group <- groups[i, , drop = FALSE]
    cells <- meta$component[meta$sample_id == group$sample_id & meta$celltype == group$celltype]
    markers <- rownames(abundance)
    data.frame(sample_id = group$sample_id, donor = group$donor, condition = group$condition,
      celltype = group$celltype, feature_id = markers, marker_1 = markers, marker_2 = markers,
      metric = "mean_normalized_abundance", is_control = markers %in% control_markers,
      value = as.numeric(Matrix::rowMeans(abundance[, cells, drop = FALSE])),
      n_total_cells = as.numeric(length(cells)), n_evaluable_cells = as.numeric(length(cells)),
      n_scored_cells = as.numeric(length(cells)))
  })
  data$contract$sample_summaries$abundance <- do.call(rbind, profiles)
  design <- unique(meta[c("sample_id", "donor", "condition")])
  if (is.null(qc)) {
    qc <- design
    qc$stage <- "00_loaded"
    qc$n_cells <- as.numeric(table(meta$sample_id)[qc$sample_id])
    data$contract$provenance$qc <- "00_loaded counts the cells in the input object. Earlier filter stages are not available."
  }
  require_contract_columns(qc, c("sample_id", "donor", "condition", "stage", "n_cells"), "QC counts", c("sample_id", "stage"))
  at <- match(qc$sample_id, design$sample_id)
  if (anyNA(at) || anyNA(qc[c("donor", "condition", "stage")]) ||
      !all(qc$donor == design$donor[at] & qc$condition == design$condition[at]) ||
      !is.numeric(qc$n_cells) || any(!is.finite(qc$n_cells) | qc$n_cells < 0 | qc$n_cells != floor(qc$n_cells)) ||
      !all(qc$stage %in% c("00_loaded", "01_after_n_umi_filter", "02_after_tau_filter", "03_after_isotype_filter"))) {
    stop("QC counts must use the sample design, supported stages, and nonnegative integer counts.", call. = FALSE)
  }
  data$contract$qc$cell_counts <- qc
  if (!is.null(saved_results)) {
    sections <- c("sample_summaries", "donor_contrasts", "differential_results")
    if (!is.list(saved_results) || is.null(names(saved_results)) || anyDuplicated(names(saved_results)) ||
        any(!names(saved_results) %in% c(sections, "report", "provenance"))) {
      stop("Saved results must contain named result sections, with optional report and provenance entries.", call. = FALSE)
    }
    for (section in intersect(sections, names(saved_results))) {
      values <- saved_results[[section]]
      if (!is.list(values) || is.null(names(values)) || anyDuplicated(names(values)) ||
          any(!names(values) %in% PROXIOMEVIS_READOUTS)) stop("Invalid saved-result readouts in ", section, ".", call. = FALSE)
      for (readout in names(values)) data$contract[[section]][[readout]] <- values[[readout]]
    }
    if (!is.null(saved_results$report)) data$contract$report <- saved_results$report
    if (!is.null(saved_results$provenance)) data$contract$provenance$imported_results <- saved_results$provenance
  }
  validate_precomputed_data(data)
  data
}

user_spatial_profiles <- function(data, proximity_dir, manifest) {
  metadata <- data$metadata
  spill_dir <- tempfile("proxiomevis-profile-")
  dir.create(spill_dir)
  on.exit(unlink(spill_dir, recursive = TRUE), add = TRUE)
  connection <- DBI::dbConnect(duckdb::duckdb(), config = list(threads = "1", memory_limit = "512MB", temp_directory = spill_dir))
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  summaries <- lapply(seq_len(nrow(manifest$partitions)), function(i) {
    partition <- manifest$partitions[i, , drop = FALSE]
    cells <- metadata[metadata$sample_id == partition$sample_id, c("component", "celltype"), drop = FALSE]
    donor <- unique(metadata$donor[metadata$sample_id == partition$sample_id])
    condition <- unique(metadata$condition[metadata$sample_id == partition$sample_id])
    DBI::dbWriteTable(connection, "sample_cells", cells, temporary = TRUE, overwrite = TRUE)
    path <- proximity_store_file(proximity_dir, partition$file)
    query <- paste0("SELECT s.celltype, p.marker_1, p.marker_2, sum(p.log2_ratio) AS score_sum, ",
      "count(*) AS n_scored_cells FROM read_parquet(", DBI::dbQuoteString(connection, path), ") p ",
      "JOIN sample_cells s ON s.component=p.component WHERE isfinite(p.log2_ratio) AND p.min_marker_fraction >= ?",
      " GROUP BY s.celltype, p.marker_1, p.marker_2")
    rows <- DBI::dbGetQuery(connection, query, params = list(manifest$donor_cutoffs[[donor]]))
    if (!nrow(rows)) return(NULL)
    totals <- table(cells$celltype)
    rows$n_total_cells <- rows$n_evaluable_cells <- as.numeric(totals[rows$celltype])
    rows$n_scored_cells <- as.numeric(rows$n_scored_cells)
    rows$value <- rows$score_sum / rows$n_total_cells
    rows$score_sum <- NULL
    rows$sample_id <- partition$sample_id
    rows$donor <- donor
    rows$condition <- condition
    rows$feature_id <- ifelse(rows$marker_1 == rows$marker_2, rows$marker_1, paste(rows$marker_1, rows$marker_2, sep = " × "))
    rows$metric <- "population_mean_log2_ratio_missing_zero"
    abundance <- data$contract$sample_summaries$abundance
    controls <- abundance$marker_1[abundance$is_control %in% TRUE]
    rows$is_control <- rows$marker_1 %in% controls | rows$marker_2 %in% controls
    rows
  })
  rows <- do.call(rbind, summaries)
  empty <- empty_user_results()$sample_summaries
  if (is.null(rows)) return(empty[c("self_proximity", "colocalization")])
  # Missing-zero means need every sample's denominator, even when that sample has no passing scores.
  populations <- stats::aggregate(rep(1L, nrow(metadata)),
    metadata[c("sample_id", "donor", "condition", "celltype")], sum)
  names(populations)[ncol(populations)] <- "n_total_cells"
  pairs <- unique(rows[c("donor", "celltype", "marker_1", "marker_2", "feature_id", "metric", "is_control")])
  grid <- merge(populations, pairs, by = c("donor", "celltype"), sort = FALSE)
  keys <- c("sample_id", "donor", "condition", "celltype", "marker_1", "marker_2")
  rows <- merge(grid, rows[c(keys, "value", "n_scored_cells")], by = keys, all.x = TRUE, sort = FALSE)
  rows$value[is.na(rows$value)] <- 0
  rows$n_scored_cells[is.na(rows$n_scored_cells)] <- 0
  rows$n_evaluable_cells <- rows$n_total_cells
  list(self_proximity = rows[rows$marker_1 == rows$marker_2, , drop = FALSE],
    colocalization = rows[rows$marker_1 != rows$marker_2, , drop = FALSE])
}

prepare_user_bundle <- function(input_rds, output_dir, ..., pxl_by_sample = NULL,
                                prepare_graphs = FALSE, isotype_column = "isotype_fraction") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Install SeuratObject for preparation.", call. = FALSE)
  input_rds <- normalizePath(input_rds, mustWork = TRUE)
  if (dir.exists(input_rds)) stop("Choose a trusted Seurat RDS file.", call. = FALSE)
  if (file.exists(output_dir)) stop("Choose a new output directory; existing data are not overwritten.", call. = FALSE)
  if (length(prepare_graphs) != 1L || is.na(prepare_graphs) || !is.logical(prepare_graphs)) {
    stop("prepare_graphs must be TRUE or FALSE.", call. = FALSE)
  }
  if (prepare_graphs && is.null(pxl_by_sample)) stop("Graph preparation requires PXL files.", call. = FALSE)
  if (requireNamespace("pixelatorR", quietly = TRUE)) loadNamespace("pixelatorR")
  source_md5 <- unname(tools::md5sum(input_rds))
  object <- readRDS(input_rds)
  data <- user_bundle_data(object, ...)
  data$contract$provenance$input_rds <- list(file = basename(input_rds), md5 = source_md5)
  if (!is.null(pxl_by_sample)) {
    samples <- unique(as.character(data$metadata[[data$contract$defaults$sample_column]]))
    if (!is.character(pxl_by_sample) || is.null(names(pxl_by_sample)) || anyNA(pxl_by_sample) ||
        anyDuplicated(names(pxl_by_sample)) || !setequal(names(pxl_by_sample), samples) ||
        !all(file.exists(pxl_by_sample)) || any(dir.exists(pxl_by_sample))) {
      stop("Supply one existing PXL file for every sample, named by sample ID.", call. = FALSE)
    }
    pxl_by_sample <- stats::setNames(normalizePath(pxl_by_sample, mustWork = TRUE), names(pxl_by_sample))
    data$contract$proximity_source <- list(pxl_by_sample = pxl_by_sample,
      pxl_manifest = data.frame(sample_alias = names(pxl_by_sample), pxl_md5 = unname(tools::md5sum(pxl_by_sample))))
    donors <- unique(as.character(data$metadata[[data$contract$defaults$donor_column]]))
    for (donor in donors) {
      if (is.null(data$contract$report$settings[[donor]]$spatial$low_abundance_filter)) {
        data$contract$report$settings[[donor]]$spatial$low_abundance_filter <- list(enabled = TRUE,
          fraction_mode = "isotype", isotype_column = isotype_column, mad_multiplier = 2, minimum_fraction = 0.001)
      }
    }
  }
  dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
  output_dir <- file.path(normalizePath(dirname(output_dir), mustWork = TRUE), basename(output_dir))
  staging <- tempfile(paste0(".", basename(output_dir), "-"), tmpdir = dirname(output_dir))
  dir.create(staging)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
  if (!is.null(pxl_by_sample)) {
    normalized <- precomputed_app_data(c(list(bundle_version = PROXIOMEVIS_BUNDLE_VERSION), data))
    manifest <- prepare_proximity_store(normalized, pxl_by_sample, file.path(staging, "proximity"))
    profiles <- user_spatial_profiles(normalized, file.path(staging, "proximity"), manifest)
    for (readout in names(profiles)) {
      if (!nrow(data$contract$sample_summaries[[readout]])) data$contract$sample_summaries[[readout]] <- profiles[[readout]]
    }
    data$contract$provenance$spatial_profiles <- paste("Prepared per-sample population means use all cells in each cell type as the denominator.",
      "Only finite observed scores passing the donor isotype threshold contribute to the sum; missing scores contribute zero.",
      "Each pair observed in a donor and cell type has a summary for every sample in that stratum, across conditions; samples with no passing scores have zero mean and zero observed support.",
      "Pairs with no passing observations anywhere in that donor and cell type remain absent. Supplied nonempty spatial summaries are retained.")
    data$assets$proximity <- list(path = "proximity", manifest = "manifest.rds")
    if (prepare_graphs) {
      directories <- unique(dirname(pxl_by_sample))
      if (length(directories) != 1L) stop("Graph preparation requires the PXL files in one directory.", call. = FALSE)
      object@misc$proxiomevis <- data$contract
      prepare_cell_graph_assets(object, directories, file.path(staging, "graphs"),
        cells = colnames(object), dataset_id = data$contract$dataset_id)
      data$assets$graphs <- list(path = "graphs", index = "index.rds")
    }
  }
  write_precomputed_bundle(data, file.path(staging, "analysis.rds"))
  if (!identical(source_md5, unname(tools::md5sum(input_rds)))) stop("The input RDS changed during preparation.", call. = FALSE)
  if (file.exists(output_dir) || !file.rename(staging, output_dir)) stop("Could not publish the new output directory.", call. = FALSE)
  invisible(file.path(output_dir, "analysis.rds"))
}

import_user_comparison <- function(bundle_path, comparison_path, output_rds) {
  if (file.exists(output_rds)) stop("Choose a new output RDS; existing files are not overwritten.", call. = FALSE)
  if (tolower(tools::file_ext(output_rds)) != "rds") stop("The output file must use the .rds extension.", call. = FALSE)
  data <- read_precomputed_rds(bundle_path)
  exported <- readRDS(comparison_path)
  if (!is.list(exported) || !identical(exported$result_version, 1L) || !is.list(exported$entry)) {
    stop("Choose a version 1 comparison exported by ProxiomeVis.", call. = FALSE)
  }
  entry <- exported$entry
  validate_hybrid_request(data, entry$request)
  data$contract <- hybrid_result_contract(data, entry)
  # The reader has resolved aliases; original overlapping mappings no longer describe this metadata.
  mappings <- list(sample_column = "sample_id", donor_column = "donor", condition_column = "condition",
    celltype_column = "celltype", cluster_column = "cluster")
  data$contract$defaults[names(mappings)] <- mappings
  validate_precomputed_data(data)
  parent <- normalizePath(dirname(output_rds), mustWork = TRUE)
  if (length(data$assets) && !identical(parent, dirname(normalizePath(bundle_path, mustWork = TRUE)))) {
    stop("Save the imported RDS beside the original bundle so its sidecar paths remain valid.", call. = FALSE)
  }
  path <- file.path(parent, basename(output_rds))
  temporary <- tempfile(paste0(".", basename(path), "-"), tmpdir = parent)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  write_precomputed_bundle(data, temporary)
  if (file.exists(path) || !file.rename(temporary, path)) stop("Could not publish the imported comparison RDS.", call. = FALSE)
  invisible(path)
}
