PROXIOMEVIS_SCHEMA_VERSION <- 1L
PROXIOMEVIS_BUNDLE_VERSION <- 2L
PROXIOMEVIS_READOUTS <- c("abundance", "self_proximity", "colocalization")

contract_error <- function(message) {
  stop("Invalid ProxiomeVis RDS: ", message, call. = FALSE)
}

require_contract_columns <- function(table, columns, label, keys) {
  if (!is.data.frame(table)) contract_error(paste(label, "must be a data frame."))
  missing <- setdiff(columns, names(table))
  if (length(missing)) {
    contract_error(paste(label, "is missing", paste(missing, collapse = ", ")))
  }
  if (anyNA(table[keys]) || any(vapply(table[keys], function(x) any(!nzchar(as.character(x))), logical(1)))) {
    contract_error(paste(label, "contains empty keys."))
  }
  if (anyDuplicated(table[keys])) contract_error(paste(label, "contains duplicate keys."))
}

validate_precomputed_header <- function(contract) {
  if (!is.list(contract) || !is.numeric(contract$schema_version) ||
      length(contract$schema_version) != 1L || is.na(contract$schema_version) ||
      contract$schema_version != PROXIOMEVIS_SCHEMA_VERSION) {
    contract_error("expected an analysis contract with schema_version = 1.")
  }
  for (field in c("dataset_id", "title")) {
    if (!is.character(contract[[field]]) || length(contract[[field]]) != 1L ||
        is.na(contract[[field]]) || !nzchar(contract[[field]])) {
      contract_error(paste(field, "must be a nonempty string."))
    }
  }
  defaults <- contract$defaults
  required_defaults <- c("reduction", "pca_reduction", "celltype_column", "cluster_column",
                         "donor_column", "condition_column", "sample_column", "assay", "abundance_layer")
  if (!is.list(defaults) || !all(required_defaults %in% names(defaults)) ||
      any(!vapply(defaults[required_defaults], function(x) is.character(x) && length(x) == 1L &&
                    !is.na(x) && nzchar(x), logical(1)))) {
    contract_error("defaults must name the reductions, metadata columns, assay, and abundance layer.")
  }
  invisible(TRUE)
}

validate_precomputed_data <- function(data) {
  if (!is.list(data)) contract_error("expected application data as a list.")
  contract <- data$contract
  validate_precomputed_header(contract)
  defaults <- contract$defaults
  metadata <- data$metadata
  if (!is.data.frame(metadata)) contract_error("metadata must be a data frame.")
  mappings <- unlist(defaults[c("celltype_column", "cluster_column", "donor_column",
                               "condition_column", "sample_column")], use.names = FALSE)
  if (!all(mappings %in% names(metadata))) contract_error("a default metadata column does not exist.")
  if (nrow(metadata) == 0L || anyDuplicated(rownames(metadata)) ||
      anyNA(rownames(metadata)) || any(!nzchar(rownames(metadata)))) contract_error("cell identifiers are invalid.")
  if (anyNA(metadata[mappings]) || any(vapply(metadata[mappings], function(x) any(!nzchar(as.character(x))), logical(1)))) {
    contract_error("mapped cell annotations contain missing values.")
  }
  if (!is.list(data$embeddings) || is.null(names(data$embeddings)) ||
      anyNA(names(data$embeddings)) || any(!nzchar(names(data$embeddings))) ||
      anyDuplicated(names(data$embeddings))) contract_error("embeddings must be a named list.")
  for (reduction in unique(c(defaults$reduction, defaults$pca_reduction))) {
    if (!reduction %in% names(data$embeddings)) contract_error(paste("missing reduction", reduction))
  }
  for (reduction in names(data$embeddings)) {
    coordinates <- data$embeddings[[reduction]]
    if (!is.matrix(coordinates) || !is.numeric(coordinates) || ncol(coordinates) < 2L ||
        anyDuplicated(rownames(coordinates)) || !setequal(rownames(coordinates), rownames(metadata)) ||
        any(!is.finite(coordinates))) contract_error(paste(reduction, "has invalid coordinates."))
  }
  abundance <- data$abundance
  if (!(is.matrix(abundance) && is.numeric(abundance)) && !inherits(abundance, "dMatrix")) {
    contract_error("abundance must be a numeric matrix.")
  }
  if (inherits(abundance, "Matrix") && !isTRUE(methods::validObject(abundance, test = TRUE))) {
    contract_error("abundance is not a valid Matrix object.")
  }
  abundance_values <- if (inherits(abundance, "Matrix")) abundance@x else abundance
  if (nrow(abundance) == 0L || is.null(rownames(abundance)) || anyNA(rownames(abundance)) ||
      any(!nzchar(rownames(abundance))) || anyDuplicated(rownames(abundance)) ||
      anyDuplicated(colnames(abundance)) || !setequal(colnames(abundance), rownames(metadata)) ||
      any(!is.finite(abundance_values))) {
    contract_error("abundance has invalid values, markers, or cell identifiers.")
  }
  if (!is.null(data$counts)) {
    counts <- data$counts
    if (!(is.matrix(counts) && is.numeric(counts)) && !inherits(counts, "dMatrix")) {
      contract_error("counts must be a numeric matrix.")
    }
    if (inherits(counts, "Matrix") && !isTRUE(methods::validObject(counts, test = TRUE))) {
      contract_error("counts is not a valid Matrix object.")
    }
    values <- if (inherits(counts, "Matrix")) counts@x else counts
    if (!identical(rownames(counts), rownames(abundance)) || anyDuplicated(colnames(counts)) ||
        !setequal(colnames(counts), rownames(metadata)) || any(!is.finite(values)) ||
        any(values < 0 | values != floor(values))) {
      contract_error("counts must contain nonnegative integers aligned to abundance markers and metadata cells.")
    }
  }
  if (!is.null(data$cell_sources)) {
    source <- data$cell_sources
    require_contract_columns(source, c("component", "original_id", "sample_id"), "cell_sources", "component")
    require_contract_columns(source, c("component", "original_id", "sample_id"), "cell_sources", c("sample_id", "original_id"))
    index <- match(source$component, rownames(metadata))
    if (!setequal(source$component, rownames(metadata)) || anyNA(index) ||
        !identical(as.character(source$sample_id), as.character(metadata[[defaults$sample_column]][index]))) {
      contract_error("cell_sources must map every metadata cell to its original sample.")
    }
  }
  sample_design <- unique(data.frame(
    sample_id = as.character(metadata[[defaults$sample_column]]),
    donor = as.character(metadata[[defaults$donor_column]]),
    condition = as.character(metadata[[defaults$condition_column]])
  ))
  if (anyDuplicated(sample_design$sample_id)) contract_error("a sample maps to multiple donors or conditions.")
  allowed_cells <- unique(as.character(metadata[[defaults$celltype_column]]))
  feature_columns <- c("celltype", "feature_id", "marker_1", "marker_2", "metric", "is_control")
  group_columns <- c("analysis_id", "group_a", "group_b", "comparison")
  sample_columns <- c(feature_columns, "sample_id", "donor", "condition", "value",
                      "n_total_cells", "n_evaluable_cells", "n_scored_cells")
  contrast_columns <- c(feature_columns, group_columns, "donor", "value_a", "value_b",
                        "effect_size", "n_cells_a", "n_cells_b")
  result_columns <- c(feature_columns, group_columns, "effect_size", "min_effect", "max_effect",
                      "n_pairs", "direction_consistent", "p_value", "p_adj", "adjustment_family",
                      "method", "replicate_unit", "inference")
  for (section in c("sample_summaries", "donor_contrasts", "differential_results")) {
    if (!is.list(contract[[section]]) || !all(PROXIOMEVIS_READOUTS %in% names(contract[[section]]))) {
      contract_error(paste(section, "must contain all three readouts."))
    }
    for (readout in PROXIOMEVIS_READOUTS) {
      table <- contract[[section]][[readout]]
      label <- paste(section, readout)
      columns <- switch(section, sample_summaries = sample_columns,
                        donor_contrasts = contrast_columns, differential_results = result_columns)
      keys <- switch(section, sample_summaries = c("sample_id", "celltype", "feature_id"),
                     donor_contrasts = c("analysis_id", "celltype", "feature_id", "donor"),
                     differential_results = c("analysis_id", "celltype", "feature_id"))
      require_contract_columns(table, columns, label, keys)
      if (!all(table$celltype %in% allowed_cells)) contract_error(paste(label, "has an unknown cell type."))
      if (!is.logical(table$is_control) || anyNA(table$is_control)) contract_error(paste(label, "has invalid control flags."))
      if (anyNA(table[c("marker_1", "marker_2", "metric")]) ||
          any(table$marker_1 > table$marker_2) ||
          !all(c(table$marker_1, table$marker_2) %in% rownames(abundance))) {
        contract_error(paste(label, "has invalid marker pairs."))
      }
      if ((readout == "colocalization" && any(table$marker_1 == table$marker_2)) ||
          (readout != "colocalization" && any(table$marker_1 != table$marker_2))) {
        contract_error(paste(label, "has the wrong pair type."))
      }
      number_columns <- intersect(c("value", "value_a", "value_b", "effect_size", "min_effect", "max_effect",
                                    "n_total_cells", "n_evaluable_cells", "n_scored_cells", "n_cells_a", "n_cells_b",
                                    "n_pairs", "p_value", "p_adj"), names(table))
      for (field in number_columns) {
        if (!is.numeric(table[[field]]) || any(is.infinite(table[[field]]))) contract_error(paste(label, field, "must be numeric and finite or NA."))
        if (startsWith(field, "n_") && any(table[[field]] < 0 | table[[field]] != floor(table[[field]]), na.rm = TRUE)) {
          contract_error(paste(label, field, "must contain nonnegative counts."))
        }
      }
      for (field in intersect(c("p_value", "p_adj"), names(table))) {
        if (any(table[[field]] < 0 | table[[field]] > 1, na.rm = TRUE)) {
          contract_error(paste(label, "has invalid probabilities."))
        }
      }
      if (section == "sample_summaries") {
        index <- match(table$sample_id, sample_design$sample_id)
        if (anyNA(index) || !all(table$donor == sample_design$donor[index]) ||
            !all(table$condition == sample_design$condition[index])) contract_error(paste(label, "does not match the sample design."))
        if (any(table$n_scored_cells > table$n_total_cells |
                table$n_evaluable_cells > table$n_total_cells |
                table$n_scored_cells > table$n_evaluable_cells, na.rm = TRUE)) {
          contract_error(paste(label, "has inconsistent support counts."))
        }
      } else {
        if (anyNA(table[c("group_a", "group_b")]) || any(table$group_a == table$group_b) ||
            !all(c(table$group_a, table$group_b) %in% sample_design$condition)) {
          contract_error(paste(label, "has an invalid contrast."))
        }
        if (section == "donor_contrasts") {
          if (!all(table$donor %in% sample_design$donor) ||
              any(abs(table$effect_size - (table$value_a - table$value_b)) > 1e-8, na.rm = TRUE)) {
            contract_error(paste(label, "has an invalid donor or contrast direction."))
          }
        } else {
          descriptive <- table$inference == "descriptive"
          if (anyNA(descriptive) || any(descriptive & (!is.na(table$p_value) | !is.na(table$p_adj)))) {
            contract_error(paste(label, "descriptive effects cannot contain inferential p-values."))
          }
          if (any(table$min_effect > table$effect_size | table$effect_size > table$max_effect, na.rm = TRUE)) {
            contract_error(paste(label, "has an inconsistent donor-effect range."))
          }
        }
      }
    }
  }
  if (!is.list(contract$provenance) || !length(contract$provenance)) contract_error("analysis provenance is missing.")
  if (!is.list(contract$qc) || !is.data.frame(contract$qc$cell_counts)) contract_error("QC cell counts are missing.")
  invisible(TRUE)
}

extract_precomputed_object <- function(object) {
  if (!inherits(object, "Seurat")) contract_error("expected a Seurat object or version 2 application bundle.")
  contract <- object@misc$proxiomevis
  validate_precomputed_header(contract)
  defaults <- contract$defaults
  metadata <- object@meta.data
  if (!identical(rownames(metadata), colnames(object))) contract_error("cell identifiers are not aligned.")
  if (!defaults$assay %in% names(object@assays)) contract_error("the default assay is missing.")
  for (assay in object@assays) {
    if ("proximity" %in% methods::slotNames(assay) && NROW(methods::slot(assay, "proximity")) > 0L) {
      contract_error("remove cell-level proximity rows before loading precomputed results.")
    }
  }
  layers <- SeuratObject::Layers(object[[defaults$assay]])
  if (!defaults$abundance_layer %in% layers || !"counts" %in% layers) {
    contract_error("the assay must contain counts and the configured abundance layer.")
  }
  data <- list(contract = contract, metadata = metadata,
               embeddings = lapply(object@reductions, function(x) x@cell.embeddings),
               abundance = SeuratObject::LayerData(object[[defaults$assay]], layer = defaults$abundance_layer),
               counts = SeuratObject::LayerData(object[[defaults$assay]], layer = "counts"))
  assay <- object[[defaults$assay]]
  if ("fs_map" %in% methods::slotNames(assay) && nrow(assay@fs_map) > 0L) {
    source <- do.call(rbind, lapply(assay@fs_map$id_map, function(mapping) {
      data.frame(component = unname(as.character(mapping$current_id)),
                 original_id = as.character(mapping$original_id))
    }))
    source <- source[source$component %in% rownames(metadata), , drop = FALSE]
    source$sample_id <- as.character(metadata[[defaults$sample_column]][match(source$component, rownames(metadata))])
    data$cell_sources <- source
  }
  data
}

validate_precomputed_object <- function(object) {
  validate_precomputed_data(extract_precomputed_object(object))
}

precomputed_app_data <- function(object) {
  if (is.list(object) && !is.null(object$bundle_version)) {
    if (!identical(object$bundle_version, PROXIOMEVIS_BUNDLE_VERSION)) {
      contract_error("expected bundle_version = 2L.")
    }
    data <- object[intersect(c("contract", "metadata", "embeddings", "abundance", "counts", "cell_sources", "assets"), names(object))]
  } else {
    data <- extract_precomputed_object(object)
  }
  validate_precomputed_data(data)
  defaults <- data$contract$defaults
  metadata <- data$metadata
  # Resolve all mappings before assigning aliases that may also be source column names.
  mapped <- lapply(c(sample_id = "sample_column", donor = "donor_column", condition = "condition_column",
                     celltype = "celltype_column", cluster = "cluster_column"),
                   function(key) as.character(metadata[[defaults[[key]]]]))
  metadata$component <- rownames(metadata)
  metadata[names(mapped)] <- mapped
  data$metadata <- metadata
  data$embeddings <- lapply(data$embeddings, function(x) x[metadata$component, , drop = FALSE])
  data$abundance <- data$abundance[, metadata$component, drop = FALSE]
  if (!is.null(data$counts)) data$counts <- data$counts[, metadata$component, drop = FALSE]
  if (!is.null(data$cell_sources)) {
    data$cell_sources <- data$cell_sources[match(metadata$component, data$cell_sources$component), , drop = FALSE]
    rownames(data$cell_sources) <- NULL
  }
  data
}

write_precomputed_bundle <- function(data, path, assets = data$assets) {
  validate_precomputed_data(data)
  bundle <- c(list(bundle_version = PROXIOMEVIS_BUNDLE_VERSION),
              data[intersect(c("contract", "metadata", "embeddings", "abundance", "counts", "cell_sources"), names(data))])
  bundle$abundance <- Matrix::Matrix(bundle$abundance, sparse = TRUE)
  if (!is.null(bundle$counts)) bundle$counts <- Matrix::Matrix(bundle$counts, sparse = TRUE)
  if (!is.null(assets)) bundle$assets <- assets
  # A viewer bundle must not carry analysis class definitions inside its contract.
  check_classes <- function(x) {
    if (isS4(x) && !inherits(x, "Matrix")) contract_error("application bundles cannot contain analysis S4 objects.")
    if (is.list(x)) invisible(lapply(x, check_classes))
  }
  check_classes(bundle)
  saveRDS(bundle, path, compress = "xz")
  invisible(path)
}

read_precomputed_rds <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path) || dir.exists(path)) {
    stop("Choose an existing precomputed .rds file on the server.", call. = FALSE)
  }
  if (tolower(tools::file_ext(path)) != "rds") stop("Expected an .rds file.", call. = FALSE)
  path <- normalizePath(path, mustWork = TRUE)
  checksum <- unname(tools::md5sum(path))
  if (is.na(checksum)) stop("Could not fingerprint the application bundle.", call. = FALSE)
  # Register only sparse-matrix classes for plain viewer bundles.
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install Matrix to load this RDS.", call. = FALSE)
  object <- readRDS(path)
  if (!identical(checksum, unname(tools::md5sum(path)))) {
    stop("The application bundle changed while loading; reload the dataset.", call. = FALSE)
  }
  if (inherits(object, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Install SeuratObject to load a legacy analysis RDS, or export a version 2 viewer bundle.", call. = FALSE)
    }
    if (requireNamespace("pixelatorR", quietly = TRUE)) loadNamespace("pixelatorR")
  }
  data <- precomputed_app_data(object)
  data$bundle_path <- path
  data$bundle_md5 <- checksum
  data
}
