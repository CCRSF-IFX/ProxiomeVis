#!/usr/bin/env Rscript

PRECOMPUTED_BUILDER_VERSION <- 1L
DEMO_DONORS <- c("d1", "d2")
DEMO_READOUTS <- c("abundance", "self_proximity", "colocalization")

.preparation_script_file <- local({
  command_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  frame_files <- unlist(lapply(sys.frames(), function(frame) frame$ofile), use.names = FALSE)
  candidates <- c(frame_files, command_file, file.path(getwd(), "scripts", "prepare_precomputed_rds.R"))
  candidates <- candidates[basename(candidates) == "prepare_precomputed_rds.R" & file.exists(candidates)]
  if (length(candidates)) normalizePath(tail(candidates, 1L), mustWork = TRUE) else NA_character_
})

preparation_error <- function(...) stop(..., call. = FALSE)

require_source_columns <- function(table, columns, label) {
  missing <- setdiff(columns, names(table))
  if (length(missing)) preparation_error(label, " is missing: ", paste(missing, collapse = ", "))
}

read_source_table <- function(path) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection), add = TRUE)
  utils::read.delim(
    connection, stringsAsFactors = FALSE, check.names = FALSE, quote = "\"",
    comment.char = "", na.strings = "NA"
  )
}

contract_file <- function() {
  candidates <- c(
    if (!is.na(.preparation_script_file)) {
      file.path(dirname(.preparation_script_file), "..", "R", "precomputed_contract.R")
    },
    file.path(getwd(), "R", "precomputed_contract.R")
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) preparation_error("Cannot find R/precomputed_contract.R.")
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

load_contract <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(contract_file(), envir = environment)
  if (!exists("validate_precomputed_object", environment, inherits = FALSE)) {
    preparation_error("R/precomputed_contract.R does not define validate_precomputed_object().")
  }
  environment
}

source_fingerprint <- function(path, role, root = NULL) {
  info <- file.info(path)
  normalized <- normalizePath(path, mustWork = TRUE)
  relative <- if (!is.null(root) && startsWith(normalized, paste0(root, .Platform$file.sep))) {
    substring(normalized, nchar(root) + 2L)
  } else {
    NA_character_
  }
  data.frame(
    role = role,
    relative_path = relative,
    source_path = normalized,
    size_bytes = as.numeric(info$size),
    modified_utc = format(info$mtime, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    md5 = unname(tools::md5sum(path)),
    stringsAsFactors = FALSE
  )
}

md5_text <- function(lines) {
  path <- tempfile("proxiomevis-fingerprint-")
  on.exit(unlink(path), add = TRUE)
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
  unname(tools::md5sum(path))
}

canonicalize_pairs <- function(table) {
  marker_1 <- as.character(table$marker_1)
  marker_2 <- as.character(table$marker_2)
  swap <- marker_1 > marker_2
  table$marker_1 <- ifelse(swap, marker_2, marker_1)
  table$marker_2 <- ifelse(swap, marker_1, marker_2)
  table
}

feature_id_for <- function(marker_1, marker_2) {
  ifelse(marker_1 == marker_2, marker_1, paste(marker_1, marker_2, sep = " × "))
}

readout_for_layer <- function(layer) {
  unname(c(abundance = "abundance", clustering = "self_proximity",
           colocalization = "colocalization")[as.character(layer)])
}

metric_for_readout <- function(readout) {
  ifelse(readout == "abundance", "delta_mean_clr", "delta_median_log2_ratio")
}

is_control_pair <- function(marker_1, marker_2, features) {
  controls <- intersect(c("mIgG1", "mIgG2a", "mIgG2b"), features)
  marker_1 %in% controls | marker_2 %in% controls
}

read_color_map <- function(path) {
  map <- jsonlite::fromJSON(path)$mapping
  stats::setNames(as.character(map$color), as.character(map$label))
}

demo_paths <- function(analysis_dir) {
  cca_dir <- file.path(dirname(analysis_dir), "cca_analysis")
  donor_root <- file.path(analysis_dir, "donor_specific")
  result_rel <- c(
    abundance = "04_differential_abundance/run_daa_results.tsv.gz",
    self_proximity = "05_differential_clustering/differential_proximity_clustering.tsv.gz",
    colocalization = "06_differential_colocalization/differential_proximity_colocalization.tsv.gz"
  )
  eligibility_rel <- c(
    abundance = "04_differential_abundance/run_daa_celltype_eligibility.tsv",
    self_proximity = "05_differential_clustering/clustering_celltype_eligibility.tsv",
    colocalization = "06_differential_colocalization/colocalization_celltype_eligibility.tsv"
  )
  profile_rel <- c(
    self_proximity = "05_differential_clustering/clustering_proximity_profiles.tsv.gz",
    colocalization = "06_differential_colocalization/colocalization_proximity_profiles.tsv.gz"
  )
  paths <- c(
    cell_counts = file.path(cca_dir, "tables", "celltype_composition.csv"),
    input_provenance = file.path(analysis_dir, "input", "input_provenance.tsv"),
    pxl_inventory = file.path(analysis_dir, "input", "pxl_provenance.tsv"),
    donor_results_long = file.path(analysis_dir, "cross_donor_review", "donor_specific_results_long.tsv"),
    cross_donor_alignment = file.path(analysis_dir, "cross_donor_review", "cross_donor_feature_alignment.tsv"),
    report_alignment = file.path(analysis_dir, "comparable_colocalization", "alignment.rds"),
    analysis_status = file.path(analysis_dir, "tertiary_analysis.ok"),
    validation_status = file.path(analysis_dir, "outputs_validated.ok"),
    result_summary = file.path(analysis_dir, "result_summary.tsv"),
    color_celltype = file.path(cca_dir, "colors", "celltype_color_map.json"),
    color_cluster = file.path(cca_dir, "colors", "cluster_color_map.json"),
    color_condition = file.path(cca_dir, "colors", "condition_color_map.json"),
    color_donor = file.path(cca_dir, "colors", "donor_color_map.json"),
    color_sample = file.path(cca_dir, "colors", "sample_color_map.json")
  )
  for (donor in DEMO_DONORS) {
    paths[paste(donor, "settings", sep = "_")] <- file.path(analysis_dir, "config", paste0(donor, ".yaml"))
    abundance_dir <- file.path(donor_root, donor, "04_differential_abundance")
    paths[paste(donor, "abundance_heatmap_manifest", sep = "_")] <-
      file.path(abundance_dir, "run_daa_heatmap_manifest.tsv")
    paths[paste(donor, "abundance_heatmap_sample_counts", sep = "_")] <-
      file.path(abundance_dir, "run_daa_heatmap_sample_cell_counts.tsv")
    matrix_files <- list.files(file.path(abundance_dir, "run_daa_heatmaps"), pattern = "[.]tsv$", full.names = TRUE)
    for (i in seq_along(matrix_files)) {
      paths[paste(donor, "abundance_heatmap_matrix", i, sep = "_")] <- matrix_files[[i]]
    }
    for (readout in DEMO_READOUTS) {
      paths[paste(donor, readout, "results", sep = "_")] <-
        file.path(donor_root, donor, result_rel[[readout]])
      paths[paste(donor, readout, "eligibility", sep = "_")] <-
        file.path(donor_root, donor, eligibility_rel[[readout]])
      if (readout != "abundance") {
        paths[paste(donor, readout, "profiles", sep = "_")] <-
          file.path(donor_root, donor, profile_rel[[readout]])
      }
    }
  }
  paths
}

make_cell_counts <- function(path) {
  source <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  require_source_columns(
    source,
    c("sample_alias", "celltype_broad", "n_components", "fraction", "donor", "condition"),
    "celltype_composition.csv"
  )
  total <- ave(source$n_components, source$sample_alias, FUN = sum)
  data.frame(
    sample_id = as.character(source$sample_alias), donor = as.character(source$donor),
    condition = as.character(source$condition), celltype = as.character(source$celltype_broad),
    n_cells = as.numeric(source$n_components), total_cells = as.numeric(total),
    fraction = as.numeric(source$fraction), stringsAsFactors = FALSE
  )
}

make_qc_stages <- function(object) {
  source <- as.data.frame(object@misc$qc_cell_counts_by_step)
  require_source_columns(
    source,
    c("step", "sample_alias", "donor", "condition", "n_components", "loaded_components", "fraction_loaded"),
    "object@misc$qc_cell_counts_by_step"
  )
  data.frame(
    stage = as.character(source$step), sample_id = as.character(source$sample_alias),
    donor = as.character(source$donor), condition = as.character(source$condition),
    n_cells = as.numeric(source$n_components), loaded_cells = as.numeric(source$loaded_components),
    fraction_loaded = as.numeric(source$fraction_loaded), stringsAsFactors = FALSE
  )
}

make_abundance_samples <- function(object, defaults) {
  metadata <- object@meta.data
  values <- SeuratObject::LayerData(object[[defaults$assay]], layer = defaults$abundance_layer)
  samples <- unique(data.frame(
    sample_id = as.character(metadata[[defaults$sample_column]]),
    donor = as.character(metadata[[defaults$donor_column]]),
    condition = as.character(metadata[[defaults$condition_column]]), stringsAsFactors = FALSE
  ))
  if (anyDuplicated(samples$sample_id)) preparation_error("Each demo sample must map to one donor and condition.")
  celltypes <- sort(unique(as.character(metadata[[defaults$celltype_column]])))
  markers <- rownames(values)
  rows <- vector("list", nrow(samples) * length(celltypes))
  index <- 0L
  for (sample_index in seq_len(nrow(samples))) {
    for (celltype in celltypes) {
      index <- index + 1L
      cells <- rownames(metadata)[
        as.character(metadata[[defaults$sample_column]]) == samples$sample_id[[sample_index]] &
          as.character(metadata[[defaults$celltype_column]]) == celltype
      ]
      value <- if (length(cells)) Matrix::rowMeans(values[, cells, drop = FALSE]) else rep(NA_real_, length(markers))
      rows[[index]] <- data.frame(
        sample_id = samples$sample_id[[sample_index]], donor = samples$donor[[sample_index]],
        condition = samples$condition[[sample_index]], celltype = celltype,
        feature_id = markers, marker_1 = markers, marker_2 = markers, value = as.numeric(value),
        n_total_cells = as.numeric(length(cells)), n_evaluable_cells = as.numeric(length(cells)),
        n_scored_cells = as.numeric(length(cells)), metric = "mean_clr",
        is_control = is_control_pair(markers, markers, markers), stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

make_spatial_samples <- function(source, donor, readout, features) {
  require_source_columns(
    source,
    c("celltype", "condition", "marker_1", "marker_2", "n_cells_detected", "n_cells", "mean_log2_ratio"),
    paste(donor, readout, "profiles")
  )
  source <- canonicalize_pairs(source)
  source <- if (readout == "self_proximity") {
    source[source$marker_1 == source$marker_2, , drop = FALSE]
  } else {
    source[source$marker_1 != source$marker_2, , drop = FALSE]
  }
  data.frame(
    sample_id = paste(donor, source$condition, sep = "_"), donor = donor,
    condition = as.character(source$condition), celltype = as.character(source$celltype),
    feature_id = feature_id_for(source$marker_1, source$marker_2),
    marker_1 = source$marker_1, marker_2 = source$marker_2,
    value = as.numeric(source$mean_log2_ratio), n_total_cells = as.numeric(source$n_cells),
    n_evaluable_cells = as.numeric(source$n_cells), n_scored_cells = as.numeric(source$n_cells_detected),
    metric = "population_mean_log2_ratio_missing_zero",
    is_control = is_control_pair(source$marker_1, source$marker_2, features), stringsAsFactors = FALSE
  )
}

make_donor_contrasts <- function(source, features, abundance_samples) {
  require_source_columns(
    source,
    c("donor", "layer", "contrast", "target", "reference", "celltype", "feature", "marker_1",
      "marker_2", "effect", "p", "p_adj", "n_target", "n_reference"),
    "donor_specific_results_long.tsv"
  )
  readout <- readout_for_layer(source$layer)
  if (anyNA(readout)) preparation_error("Unknown readout in donor-specific results.")
  marker_1 <- as.character(source$marker_1)
  marker_2 <- as.character(source$marker_2)
  abundance <- readout == "abundance"
  colocalization <- readout == "colocalization"
  marker_1[abundance] <- as.character(source$feature[abundance])
  marker_2[abundance] <- marker_1[abundance]
  # The frozen cross-donor table stores co-marker 1 in feature, co-marker 2 in marker_1,
  # and the display pair in marker_2.
  marker_2[colocalization] <- marker_1[colocalization]
  marker_1[colocalization] <- as.character(source$feature[colocalization])
  pairs <- canonicalize_pairs(data.frame(marker_1 = marker_1, marker_2 = marker_2))
  value_a <- value_b <- rep(NA_real_, nrow(source))
  abundance_key <- paste(
    abundance_samples$donor, abundance_samples$condition, abundance_samples$celltype,
    abundance_samples$feature_id, sep = "\r"
  )
  if (any(abundance)) {
    key_a <- paste(source$donor, source$target, source$celltype, source$feature, sep = "\r")
    key_b <- paste(source$donor, source$reference, source$celltype, source$feature, sep = "\r")
    value_a[abundance] <- abundance_samples$value[match(key_a[abundance], abundance_key)]
    value_b[abundance] <- abundance_samples$value[match(key_b[abundance], abundance_key)]
    if (anyNA(value_a[abundance]) || anyNA(value_b[abundance]) ||
        any(abs(value_a[abundance] - value_b[abundance] - source$effect[abundance]) > 1e-8)) {
      preparation_error("Frozen abundance means do not reproduce the source donor effect.")
    }
  }
  data.frame(
    analysis_id = paste(readout, source$contrast, sep = ":"), donor = as.character(source$donor),
    celltype = as.character(source$celltype), feature_id = feature_id_for(pairs$marker_1, pairs$marker_2),
    marker_1 = pairs$marker_1, marker_2 = pairs$marker_2,
    group_a = as.character(source$target), group_b = as.character(source$reference),
    comparison = paste(source$target, "vs", source$reference),
    value_a = value_a, value_b = value_b,
    effect_size = as.numeric(source$effect), n_cells_a = as.numeric(source$n_target),
    n_cells_b = as.numeric(source$n_reference), metric = metric_for_readout(readout),
    is_control = is_control_pair(pairs$marker_1, pairs$marker_2, features),
    p_value = as.numeric(source$p), p_adj = as.numeric(source$p_adj),
    adjustment_family = paste(source$donor, readout, source$contrast, source$celltype, "BH", sep = ":"),
    method = ifelse(abundance, "RunDAA (Wilcoxon/FindMarkers)", "DifferentialProximityAnalysis"),
    replicate_unit = "component", inference = "exploratory_component_level",
    stringsAsFactors = FALSE
  )
}

make_combined_results <- function(source, features) {
  require_source_columns(
    source,
    c("layer", "contrast", "target", "reference", "celltype", "feature", "marker_1", "marker_2",
      "effect_d1", "effect_d2", "tested_d1", "tested_d2", "same_direction", "support_status"),
    "cross_donor_feature_alignment.tsv"
  )
  readout <- readout_for_layer(source$layer)
  marker_1 <- as.character(source$marker_1)
  marker_2 <- as.character(source$marker_2)
  abundance <- readout == "abundance"
  colocalization <- readout == "colocalization"
  marker_1[abundance] <- as.character(source$feature[abundance])
  marker_2[abundance] <- marker_1[abundance]
  marker_2[colocalization] <- marker_1[colocalization]
  marker_1[colocalization] <- as.character(source$feature[colocalization])
  pairs <- canonicalize_pairs(data.frame(marker_1 = marker_1, marker_2 = marker_2))
  effects <- cbind(as.numeric(source$effect_d1), as.numeric(source$effect_d2))
  n_pairs <- rowSums(is.finite(effects))
  effect_size <- rowMeans(effects, na.rm = TRUE)
  effect_size[n_pairs == 0L] <- NA_real_
  min_effect <- apply(effects, 1L, function(x) if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_)
  max_effect <- apply(effects, 1L, function(x) if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_)
  data.frame(
    analysis_id = paste(readout, source$contrast, sep = ":"), celltype = as.character(source$celltype),
    feature_id = feature_id_for(pairs$marker_1, pairs$marker_2),
    marker_1 = pairs$marker_1, marker_2 = pairs$marker_2,
    group_a = as.character(source$target), group_b = as.character(source$reference),
    comparison = paste(source$target, "vs", source$reference), effect_size = effect_size,
    min_effect = min_effect, max_effect = max_effect, n_pairs = as.numeric(n_pairs),
    direction_consistent = as.logical(source$same_direction), metric = metric_for_readout(readout),
    is_control = is_control_pair(pairs$marker_1, pairs$marker_2, features),
    p_value = rep(NA_real_, nrow(source)), p_adj = rep(NA_real_, nrow(source)),
    adjustment_family = "not_applicable_descriptive",
    method = "mean of available donor-specific component-level effects",
    replicate_unit = "donor", inference = "descriptive",
    support_status = as.character(source$support_status), stringsAsFactors = FALSE
  )
}

read_abundance_matrix <- function(path) {
  source <- read_source_table(path)
  if (!"marker" %in% names(source)) preparation_error("Abundance heatmap matrix lacks marker: ", path)
  matrix <- as.matrix(source[setdiff(names(source), "marker")])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- as.character(source$marker)
  matrix
}

read_status_file <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "=", fixed = TRUE)
  fields <- fields[lengths(fields) == 2L]
  stats::setNames(
    lapply(fields, function(field) field[[2L]]),
    vapply(fields, `[[`, character(1), 1L)
  )
}

load_abundance_heatmaps <- function(paths, donor, abundance_samples) {
  manifest <- read_source_table(paths[[paste(donor, "abundance_heatmap_manifest", sep = "_")]])
  sample_counts <- read_source_table(paths[[paste(donor, "abundance_heatmap_sample_counts", sep = "_")]])
  require_source_columns(manifest, c("contrast", "celltype", "matrix_file", "color_min", "color_max"), "abundance heatmap manifest")
  matrices <- list()
  for (i in seq_len(nrow(manifest))) {
    contrast <- as.character(manifest$contrast[[i]])
    celltype <- as.character(manifest$celltype[[i]])
    if (is.null(matrices[[contrast]])) matrices[[contrast]] <- list()
    source_path <- as.character(manifest$matrix_file[[i]])
    if (!file.exists(source_path)) {
      source_path <- file.path(
        dirname(paths[[paste(donor, "abundance_heatmap_manifest", sep = "_")]]),
        "run_daa_heatmaps", basename(source_path)
      )
    }
    matrices[[contrast]][[celltype]] <- read_abundance_matrix(source_path)
  }
  first_matrix <- matrices[[1L]][[1L]]
  sample_order <- colnames(first_matrix)
  annotations <- unique(abundance_samples[c("sample_id", "donor", "condition")])
  annotations <- annotations[match(sample_order, annotations$sample_id), , drop = FALSE]
  if (anyNA(annotations$sample_id)) preparation_error("Heatmap sample columns do not match the input RDS.")
  names(annotations)[names(annotations) == "sample_id"] <- "sample"
  list(manifest = manifest, sample_counts = sample_counts, samples = annotations, matrices = matrices)
}

load_report <- function(paths, abundance_samples) {
  source_results <- eligibility <- profiles <- settings <- abundance_heatmaps <-
    setNames(vector("list", length(DEMO_DONORS)), DEMO_DONORS)
  for (donor in DEMO_DONORS) {
    source_results[[donor]] <- eligibility[[donor]] <- profiles[[donor]] <- list()
    for (readout in DEMO_READOUTS) {
      source_results[[donor]][[readout]] <- read_source_table(paths[[paste(donor, readout, "results", sep = "_")]])
      eligibility[[donor]][[readout]] <- read_source_table(paths[[paste(donor, readout, "eligibility", sep = "_")]])
      profiles[[donor]][[readout]] <- if (readout == "abundance") {
        abundance_samples[abundance_samples$donor == donor, , drop = FALSE]
      } else {
        read_source_table(paths[[paste(donor, readout, "profiles", sep = "_")]])
      }
    }
    settings[[donor]] <- yaml::read_yaml(paths[[paste(donor, "settings", sep = "_")]])
    abundance_heatmaps[[donor]] <- load_abundance_heatmaps(paths, donor, abundance_samples)
  }
  list(
    alignment = readRDS(paths[["report_alignment"]]), source_results = source_results,
    eligibility = eligibility, profiles = profiles, settings = settings,
    abundance_heatmaps = abundance_heatmaps,
    status = list(
      analysis = read_status_file(paths[["analysis_status"]]),
      validation = read_status_file(paths[["validation_status"]]),
      results = read_source_table(paths[["result_summary"]])
    )
  )
}

build_precomputed_rds <- function(input_rds, analysis_dir, output_rds,
    reduction = "umap.cca", pca_reduction = "pca", celltype_column = "celltype_broad",
    cluster_column = "cluster", donor_column = "donor", condition_column = "condition",
    sample_column = "sample_alias", assay = "PNA", abundance_layer = "data") {
  defaults <- list(reduction = reduction, pca_reduction = pca_reduction,
    celltype_column = celltype_column, cluster_column = cluster_column, donor_column = donor_column,
    condition_column = condition_column, sample_column = sample_column, assay = assay,
    abundance_layer = abundance_layer)
  if (any(!vapply(defaults, function(x) is.character(x) && length(x) == 1L &&
      !is.na(x) && nzchar(trimws(x)), logical(1)))) {
    preparation_error("Reduction, metadata, assay, and layer settings must be nonempty strings.")
  }
  arguments <- c(input_rds, analysis_dir, output_rds)
  if (length(arguments) != 3L || anyNA(arguments) || any(!nzchar(arguments))) {
    preparation_error("input_rds, analysis_dir, and output_rds must be nonempty paths.")
  }
  input_rds <- normalizePath(input_rds, mustWork = TRUE)
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  if (dir.exists(input_rds)) preparation_error("input_rds must be a file.")
  if (tolower(tools::file_ext(output_rds)) != "rds") preparation_error("output_rds must end in .rds.")
  dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
  output_rds <- file.path(normalizePath(dirname(output_rds), mustWork = TRUE), basename(output_rds))
  if (identical(input_rds, normalizePath(output_rds, mustWork = FALSE))) {
    preparation_error("output_rds must not overwrite input_rds.")
  }
  demo_source_root <- normalizePath(dirname(analysis_dir), mustWork = TRUE)
  if (identical(output_rds, demo_source_root) ||
      startsWith(output_rds, paste0(demo_source_root, .Platform$file.sep))) {
    preparation_error("output_rds must be outside the demo_analysis source tree.")
  }
  for (package in c("SeuratObject", "Matrix", "jsonlite", "yaml")) {
    if (!requireNamespace(package, quietly = TRUE)) preparation_error("Install ", package, ".")
  }
  if (requireNamespace("pixelatorR", quietly = TRUE)) loadNamespace("pixelatorR")

  paths <- demo_paths(analysis_dir)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) preparation_error("Missing demo analysis sources: ", paste(basename(missing), collapse = ", "))
  object <- readRDS(input_rds)
  if (!inherits(object, "Seurat")) preparation_error("input_rds must contain a Seurat object.")
  required_metadata <- unlist(defaults[c("celltype_column", "cluster_column", "donor_column",
                                         "condition_column", "sample_column")], use.names = FALSE)
  if (!all(required_metadata %in% names(object@meta.data))) {
    preparation_error("Input RDS is missing metadata: ",
                      paste(setdiff(required_metadata, names(object@meta.data)), collapse = ", "))
  }
  if (!all(c(defaults$reduction, defaults$pca_reduction) %in% names(object@reductions))) {
    preparation_error("Input RDS is missing the selected UMAP or PCA reduction.")
  }
  if (!defaults$assay %in% names(object@assays)) preparation_error("Input RDS is missing the selected assay.")
  if (!defaults$abundance_layer %in% SeuratObject::Layers(object[[defaults$assay]])) {
    preparation_error("Input RDS is missing the selected abundance layer.")
  }
  features <- rownames(object[[defaults$assay]])

  input_provenance <- read_source_table(paths[["input_provenance"]])
  require_source_columns(input_provenance, c("source_md5"), "input_provenance.tsv")
  input_md5 <- unname(tools::md5sum(input_rds))
  if (!identical(input_md5, as.character(input_provenance$source_md5[[1L]]))) {
    preparation_error("input_rds MD5 does not match the frozen demo source MD5.")
  }

  cell_counts <- make_cell_counts(paths[["cell_counts"]])
  abundance_samples <- make_abundance_samples(object, defaults)
  report <- load_report(paths, abundance_samples)
  spatial_samples <- lapply(c("self_proximity", "colocalization"), function(readout) {
    do.call(rbind, lapply(DEMO_DONORS, function(donor) {
      make_spatial_samples(report$profiles[[donor]][[readout]], donor, readout, features)
    }))
  })
  names(spatial_samples) <- c("self_proximity", "colocalization")
  sample_summaries <- c(list(abundance = abundance_samples), spatial_samples)

  all_donor_contrasts <- make_donor_contrasts(
    read_source_table(paths[["donor_results_long"]]), features, abundance_samples
  )
  donor_contrasts <- lapply(DEMO_READOUTS, function(readout) {
    all_donor_contrasts[startsWith(all_donor_contrasts$analysis_id, paste0(readout, ":")), , drop = FALSE]
  })
  names(donor_contrasts) <- DEMO_READOUTS
  all_results <- make_combined_results(read_source_table(paths[["cross_donor_alignment"]]), features)
  differential_results <- lapply(DEMO_READOUTS, function(readout) {
    all_results[startsWith(all_results$analysis_id, paste0(readout, ":")), , drop = FALSE]
  })
  names(differential_results) <- DEMO_READOUTS

  pxl_inventory <- read_source_table(paths[["pxl_inventory"]])
  require_source_columns(pxl_inventory, c("sample_alias", "pxl_file", "size_bytes", "pxl_md5"), "pxl_provenance.tsv")
  if (anyDuplicated(pxl_inventory$sample_alias)) preparation_error("PXL provenance has duplicate samples.")
  pxl_by_sample <- stats::setNames(as.character(pxl_inventory$pxl_file), as.character(pxl_inventory$sample_alias))
  sample_ids <- unique(as.character(object@meta.data[[defaults$sample_column]]))
  if (!setequal(names(pxl_by_sample), sample_ids)) preparation_error("PXL provenance does not match input samples.")

  contract_environment <- load_contract()
  analysis_root <- demo_source_root
  manifest <- do.call(rbind, lapply(names(paths), function(role) source_fingerprint(paths[[role]], role, analysis_root)))
  rownames(manifest) <- NULL
  input_fingerprint <- source_fingerprint(input_rds, "input_rds", analysis_root)
  validator_fingerprint <- source_fingerprint(contract_file(), "contract")
  builder_fingerprint <- if (!is.na(.preparation_script_file)) {
    source_fingerprint(.preparation_script_file, "builder")
  } else {
    data.frame(role = "builder", relative_path = NA_character_, source_path = NA_character_,
               size_bytes = NA_real_, modified_utc = NA_character_,
               md5 = paste0("version-", PRECOMPUTED_BUILDER_VERSION), stringsAsFactors = FALSE)
  }
  identity_rows <- rbind(input_fingerprint, manifest, validator_fingerprint, builder_fingerprint)
  identity_rows <- identity_rows[order(identity_rows$role), , drop = FALSE]
  dataset_id <- paste0(
    "pixelgen-demo-v", contract_environment$PROXIOMEVIS_SCHEMA_VERSION, "-",
    substr(md5_text(c(
      paste0("builder_version=", PRECOMPUTED_BUILDER_VERSION),
      paste(names(defaults), unlist(defaults, use.names = FALSE), sep = "="),
      paste(identity_rows$role, identity_rows$size_bytes, identity_rows$md5, sep = "=")
    )), 1L, 16L)
  )

  for (assay_name in names(object@assays)) {
    assay <- object@assays[[assay_name]]
    if ("proximity" %in% methods::slotNames(assay)) {
      proximity <- methods::slot(assay, "proximity")
      methods::slot(assay, "proximity") <- proximity[0, , drop = FALSE]
      object@assays[[assay_name]] <- assay
    }
  }
  color_dir <- file.path(dirname(analysis_dir), "cca_analysis", "colors")
  object@misc$proxiomevis <- list(
    schema_version = as.integer(contract_environment$PROXIOMEVIS_SCHEMA_VERSION),
    dataset_id = dataset_id,
    title = "Pixelgen PBMC demo: PHA versus unstimulated",
    defaults = defaults,
    sample_summaries = sample_summaries,
    donor_contrasts = donor_contrasts,
    differential_results = differential_results,
    qc = list(cell_counts = make_qc_stages(object), celltype_counts = cell_counts),
    colors = list(
      celltype = read_color_map(file.path(color_dir, "celltype_color_map.json")),
      donor = read_color_map(file.path(color_dir, "donor_color_map.json")),
      condition = read_color_map(file.path(color_dir, "condition_color_map.json")),
      cluster = read_color_map(file.path(color_dir, "cluster_color_map.json")),
      sample = read_color_map(file.path(color_dir, "sample_color_map.json"))
    ),
    proximity_source = list(
      representation = "authoritative server-mounted Pixelator PXL files",
      pxl_by_sample = pxl_by_sample, pxl_manifest = pxl_inventory,
      layout_candidates = c("coarsened_pmds_3d", "wpmds_3d"),
      note = "Cell-level proximity is absent from this RDS; one selected component may be queried from its PXL file."
    ),
    report = report,
    provenance = list(
      created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      builder_version = PRECOMPUTED_BUILDER_VERSION, input_rds = input_fingerprint,
      source_table_manifest = manifest, contract = validator_fingerprint, builder = builder_fingerprint,
      analysis_parameters = object@misc$analysis_parameters,
      support_semantics = list(
        abundance_samples = "Mean normalized CLR across all cells in each sample/cell-type stratum; no statistical test recomputed.",
        spatial_samples = paste(
          "Source population mean log2-ratio includes missing observations as zero; n_total_cells and",
          "n_evaluable_cells are source n_cells, while n_scored_cells is source n_cells_detected."
        ),
        donor_contrasts = paste(
          "Original component-level p and BH-adjusted p values are preserved unchanged.",
          "Abundance value_a/value_b are frozen mean CLR values; spatial group medians were not persisted",
          "so spatial value_a/value_b are NA and effect_size is the authoritative source effect."
        ),
        combined_results = "Mean/min/max summarize available donor effects; p values are NA and inference is descriptive."
      ),
      transformation = "Existing report results were reshaped; only descriptive sample means and cross-donor summaries were computed.",
      software = list(
        R = R.version.string, SeuratObject = as.character(utils::packageVersion("SeuratObject")),
        pixelatorR = if (requireNamespace("pixelatorR", quietly = TRUE)) {
          as.character(utils::packageVersion("pixelatorR"))
        } else NA_character_
      )
    )
  )

  contract_environment$validate_precomputed_object(object)
  temporary_output <- tempfile(paste0(".", basename(output_rds), "."), tmpdir = dirname(output_rds))
  on.exit(if (file.exists(temporary_output)) unlink(temporary_output), add = TRUE)
  saveRDS(object, temporary_output, compress = "gzip")
  if (!file.rename(temporary_output, output_rds)) preparation_error("Could not atomically replace output: ", output_rds)
  invisible(output_rds)
}

if (sys.nframe() == 0L) {
  source(file.path(dirname(.preparation_script_file), "..", "R", "preparation_cli.R"))
  option_names <- c(umap = "reduction", pca = "pca_reduction", "celltype-column" = "celltype_column",
    "cluster-column" = "cluster_column", "donor-column" = "donor_column",
    "condition-column" = "condition_column", "sample-column" = "sample_column",
    assay = "assay", "abundance-layer" = "abundance_layer")
  usage <- paste("Usage: Rscript scripts/prepare_precomputed_rds.R INPUT_RDS ANALYSIS_DIR OUTPUT_RDS",
    paste(paste0("[--", names(option_names), " VALUE]"), collapse = " "))
  args <- commandArgs(trailingOnly = TRUE)
  if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0L) }
  parsed <- parse_preparation_cli(args, setNames(vector("list", length(option_names)), names(option_names)), 3L, usage)
  settings <- parsed$options[!vapply(parsed$options, is.null, logical(1))]
  names(settings) <- unname(option_names[names(settings)])
  output <- do.call(build_precomputed_rds, c(as.list(parsed$positional), settings))
  message("Wrote validated precomputed RDS: ", output)
}
