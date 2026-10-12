HYBRID_LIMITS <- list(cells = 2000L, markers = 160L, pairs = 1000L, rows = 1000000L)
HYBRID_ANALYSIS_VERSION <- "1"

validate_hybrid_request <- function(data, request, limits = HYBRID_LIMITS) {
  if (!is.list(request)) stop("An analysis request must be a list.", call. = FALSE)
  if (!is.null(request$dataset_id) && !identical(request$dataset_id, data$contract$dataset_id)) {
    stop("The dataset changed after this request was selected; select cells again.", call. = FALSE)
  }
  if (!is.null(request$bundle_md5) && !identical(request$bundle_md5, data$bundle_md5)) {
    stop("The application bundle changed after this request was selected; reload the dataset.", call. = FALSE)
  }
  for (field in c("readout", "donor", "celltype", "target", "reference")) {
    value <- request[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("Choose one ", field, ".", call. = FALSE)
    }
  }
  if (!request$readout %in% c("abundance", "self_proximity", "colocalization")) {
    stop("Unknown analysis readout.", call. = FALSE)
  }
  if (request$target == request$reference) stop("Choose two different conditions.", call. = FALSE)
  metadata <- data$metadata
  scope <- metadata$donor == request$donor & metadata$celltype == request$celltype
  if (!any(scope)) stop("No cells match this donor and cell type.", call. = FALSE)
  select_cells <- function(field, condition) {
    available <- metadata$component[scope & metadata$condition == condition]
    cells <- request[[field]]
    if (is.null(cells)) cells <- available
    if (!is.character(cells) || anyNA(cells) || anyDuplicated(cells) ||
        !all(cells %in% available)) {
      stop("Selected cells must be unique and belong to the selected donor, cell type, and condition.", call. = FALSE)
    }
    sort(cells, method = "radix")
  }
  group_a_cells <- select_cells("group_a_cells", request$target)
  group_b_cells <- select_cells("group_b_cells", request$reference)
  min_cells <- switch(request$readout, abundance = 3L, self_proximity = 20L, colocalization = 40L)
  if (min(length(group_a_cells), length(group_b_cells)) < min_cells) {
    stop("This readout requires at least ", min_cells, " selected cells in each condition.", call. = FALSE)
  }
  if (length(c(group_a_cells, group_b_cells)) > limits$cells) {
    stop("Select at most ", limits$cells, " cells in total.", call. = FALSE)
  }
  markers <- request$markers
  if (is.null(markers) || identical(markers, character())) markers <- rownames(data$abundance)
  if (!is.character(markers) || anyNA(markers) || !length(markers) ||
      !all(markers %in% rownames(data$abundance))) stop("Choose known markers.", call. = FALSE)
  markers <- sort(unique(markers), method = "radix")
  if (length(markers) > limits$markers) stop("Too many selected markers.", call. = FALSE)
  if (request$readout == "colocalization" &&
      (length(markers) < 2L || choose(length(markers), 2) > limits$pairs)) {
    stop("Select at least two markers and at most ", limits$pairs, " marker pairs (45 markers by default).", call. = FALSE)
  }
  list(dataset_id = data$contract$dataset_id, bundle_md5 = data$bundle_md5,
       readout = request$readout, donor = request$donor, celltype = request$celltype,
       target = request$target, reference = request$reference, markers = markers,
       group_a_cells = group_a_cells, group_b_cells = group_b_cells,
       min_cells = min_cells)
}

hybrid_feature_table <- function(request) {
  if (request$readout == "colocalization") {
    pairs <- utils::combn(request$markers, 2L)
    result <- data.frame(marker_1 = pairs[1L, ], marker_2 = pairs[2L, ])
  } else {
    result <- data.frame(marker_1 = request$markers, marker_2 = request$markers)
  }
  result$feature_id <- ifelse(result$marker_1 == result$marker_2, result$marker_1,
                              paste(result$marker_1, result$marker_2, sep = " × "))
  result
}

run_hybrid_analysis <- function(bundle_path, request, asset_root = NULL) {
  data <- read_precomputed_rds(bundle_path)
  request <- validate_hybrid_request(data, request)
  if (!requireNamespace("pixelatorR", quietly = TRUE)) {
    stop("Install pixelatorR in the analysis worker to run differential tests.", call. = FALSE)
  }
  cells <- c(request$group_a_cells, request$group_b_cells)
  results <- hybrid_feature_table(request)
  results$n_cells_a <- length(request$group_a_cells)
  results$n_cells_b <- length(request$group_b_cells)
  results$n_scored_a <- 0L
  results$n_scored_b <- 0L
  results$effect_size <- results$p_value <- results$p_adj <- NA_real_
  filter_provenance <- NULL
  if (request$readout == "abundance") {
    if (is.null(data$counts)) stop("This bundle does not contain counts for abundance testing.", call. = FALSE)
    metadata <- data.frame(condition = rep(c(request$target, request$reference),
                                           c(length(request$group_a_cells), length(request$group_b_cells))),
                           row.names = cells)
    object <- SeuratObject::CreateSeuratObject(data$counts[request$markers, cells, drop = FALSE],
                                               assay = "PNA", meta.data = metadata)
    SeuratObject::LayerData(object, assay = "PNA", layer = "data") <-
      data$abundance[request$markers, cells, drop = FALSE]
    tested <- pixelatorR::RunDAA(object, contrast_column = "condition", targets = request$target,
                                reference = request$reference, assay = "PNA", mean_fxn = Matrix::rowMeans,
                                fc_name = "difference", p_adjust_method = "BH", slot = "data", fc.slot = "data",
                                min.pct = 0, test.use = "wilcox", verbose = FALSE)
    index <- match(tested$marker, results$marker_1)
    results$effect_size[index] <- tested$difference
    results$p_value[index] <- tested$p
    results$p_adj[index] <- tested$p_adj
    results$n_scored_a <- results$n_cells_a
    results$n_scored_b <- results$n_cells_b
    method <- "pixelatorR::RunDAA / Seurat Wilcoxon; difference in mean stored normalized abundance"
  } else {
    if (is.null(asset_root)) {
      asset_root <- resolve_hybrid_asset(bundle_path, data$assets$proximity)
    }
    scores <- read_proximity_subset(data, asset_root, request)
    filter_provenance <- attr(scores, "provenance")
    scores$condition <- ifelse(scores$component %in% request$group_a_cells, request$target, request$reference)
    score_keys <- paste(scores$marker_1, scores$marker_2, sep = "\r")
    result_keys <- paste(results$marker_1, results$marker_2, sep = "\r")
    for (group in c("a", "b")) {
      counts <- table(score_keys[scores$condition == if (group == "a") request$target else request$reference])
      values <- as.integer(counts[result_keys])
      values[is.na(values)] <- 0L
      results[[paste0("n_scored_", group)]] <- values
    }
    eligible <- results$n_scored_a >= request$min_cells & results$n_scored_b >= request$min_cells
    if (any(eligible)) {
      tested <- pixelatorR::DifferentialProximityAnalysis(
        scores, contrast_column = "condition", targets = request$target, reference = request$reference,
        proximity_metric = "log2_ratio", metric_type = if (request$readout == "self_proximity") "self" else "co",
        p_adjust_method = "BH", min_cells_per_group = request$min_cells, verbose = FALSE)
      index <- match(paste(tested$marker_1, tested$marker_2, sep = "\r"), result_keys)
      results$effect_size[index] <- tested$diff_median
      results$p_value[index] <- tested$p
      results$p_adj[index] <- tested$p_adj
    }
    method <- "pixelatorR::DifferentialProximityAnalysis; two-sided tie-corrected rank test; target minus reference median log2 ratio"
  }
  results$status <- ifelse(results$n_scored_a < request$min_cells | results$n_scored_b < request$min_cells,
                            "insufficient_observed_scores", ifelse(is.finite(results$p_value), "tested", "no_finite_test_result"))
  results$readout <- request$readout
  results$donor <- request$donor
  results$celltype <- request$celltype
  results$group_a <- request$target
  results$group_b <- request$reference
  results$inference <- "exploratory_component_level"
  list(results = results, provenance = list(
    dataset_id = data$contract$dataset_id, bundle_md5 = data$bundle_md5,
    analysis_version = HYBRID_ANALYSIS_VERSION, request = request,
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    method = method, p_adjust_method = "BH", p_adjust_scope = "eligible tested features within this selected-cell request",
    abundance_layer = data$contract$defaults$abundance_layer,
    requested_features = nrow(results), tested_features = sum(results$status == "tested"),
    missing_score_policy = "Spatial tests use retained observed finite scores only; absent scores are not zeros.",
    inference = "Exploratory cell-level comparison within one donor; cells are not biological replicates.",
    filters = filter_provenance, software = list(R = R.version.string,
      pixelatorR = as.character(utils::packageVersion("pixelatorR")),
      Seurat = as.character(utils::packageVersion("Seurat"))))
  )
}
