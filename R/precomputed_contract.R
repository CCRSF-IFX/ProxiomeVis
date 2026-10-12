PROXIOMEVIS_SCHEMA_VERSION <- 1L
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

validate_precomputed_object <- function(object) {
  if (!inherits(object, "Seurat")) contract_error("expected a Seurat object.")
  contract <- object@misc$proxiomevis
  if (!is.list(contract) || !is.numeric(contract$schema_version) ||
      length(contract$schema_version) != 1L || is.na(contract$schema_version) ||
      contract$schema_version != PROXIOMEVIS_SCHEMA_VERSION) {
    contract_error("expected misc$proxiomevis with schema_version = 1.")
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
  metadata <- object@meta.data
  mappings <- unlist(defaults[c("celltype_column", "cluster_column", "donor_column",
                               "condition_column", "sample_column")], use.names = FALSE)
  if (!all(mappings %in% names(metadata))) contract_error("a default metadata column does not exist.")
  if (nrow(metadata) == 0L || anyDuplicated(rownames(metadata)) ||
      !identical(rownames(metadata), colnames(object))) contract_error("cell identifiers are not aligned.")
  if (anyNA(metadata[mappings]) || any(vapply(metadata[mappings], function(x) any(!nzchar(as.character(x))), logical(1)))) {
    contract_error("mapped cell annotations contain missing values.")
  }
  for (reduction in unique(c(defaults$reduction, defaults$pca_reduction))) {
    if (!reduction %in% names(object@reductions)) contract_error(paste("missing reduction", reduction))
    coordinates <- object@reductions[[reduction]]@cell.embeddings
    if (ncol(coordinates) < 2L || !setequal(rownames(coordinates), rownames(metadata)) ||
        any(!is.finite(coordinates))) contract_error(paste(reduction, "has invalid coordinates."))
  }
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
          !all(c(table$marker_1, table$marker_2) %in% rownames(object[[defaults$assay]]))) {
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

precomputed_app_data <- function(object) {
  validate_precomputed_object(object)
  contract <- object@misc$proxiomevis
  defaults <- contract$defaults
  metadata <- object@meta.data
  # Resolve all mappings before assigning aliases that may also be source column names.
  mapped <- lapply(c(sample_id = "sample_column", donor = "donor_column", condition = "condition_column",
                     celltype = "celltype_column", cluster = "cluster_column"),
                   function(key) as.character(metadata[[defaults[[key]]]]))
  metadata$component <- rownames(metadata)
  metadata[names(mapped)] <- mapped
  embeddings <- lapply(object@reductions, function(x) x@cell.embeddings[rownames(metadata), , drop = FALSE])
  abundance <- SeuratObject::LayerData(object[[defaults$assay]], layer = defaults$abundance_layer)
  if (!setequal(colnames(abundance), metadata$component)) contract_error("abundance cells do not match metadata.")
  list(contract = contract, metadata = metadata, embeddings = embeddings,
       abundance = abundance[, metadata$component, drop = FALSE])
}

read_precomputed_rds <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path) || dir.exists(path)) {
    stop("Choose an existing precomputed .rds file on the server.", call. = FALSE)
  }
  if (tolower(tools::file_ext(path)) != "rds") stop("Expected an .rds file.", call. = FALSE)
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Install SeuratObject to load this RDS.", call. = FALSE)
  # PNA assays need their class definition; ordinary Seurat fixtures do not.
  if (requireNamespace("pixelatorR", quietly = TRUE)) loadNamespace("pixelatorR")
  object <- readRDS(path)
  precomputed_app_data(object)
}
