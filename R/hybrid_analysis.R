HYBRID_LIMITS <- list(cells = 2000L, markers = 160L, pairs = 1000L, rows = 1000000L,
                      donors = 40L, total_cells = 80000L)
HYBRID_ANALYSIS_VERSION <- "3"

validate_hybrid_request <- function(data, request, limits = HYBRID_LIMITS) {
  if (!is.list(request)) stop("An analysis request must be a list.", call. = FALSE)
  if (!is.null(request$dataset_id) && !identical(request$dataset_id, data$contract$dataset_id)) {
    stop("The dataset changed after this request was selected; select cells again.", call. = FALSE)
  }
  if (!is.null(request$bundle_md5) && !identical(request$bundle_md5, data$bundle_md5)) {
    stop("The application bundle changed after this request was selected; reload the dataset.", call. = FALSE)
  }
  analysis_scope <- if (is.null(request$analysis_scope)) "within_donor" else request$analysis_scope
  axis <- if (is.null(request$comparison_axis)) "condition" else request$comparison_axis
  design <- if (is.null(request$design)) {
    if (identical(analysis_scope, "within_donor")) "within_donor" else "paired"
  } else request$design
  if (!identical(analysis_scope, "within_donor") && !identical(analysis_scope, "across_donors")) {
    stop("Choose within-donor exploration or an across-donor comparison.", call. = FALSE)
  }
  if (!identical(axis, "condition") && !identical(axis, "celltype")) stop("Choose a comparison axis.", call. = FALSE)
  if ((analysis_scope == "within_donor" && !identical(design, "within_donor")) ||
      (analysis_scope == "across_donors" && !identical(design, "paired") && !identical(design, "independent"))) {
    stop("Choose a donor design that matches the analysis scope.", call. = FALSE)
  }
  if (axis == "celltype" && design == "independent") {
    stop("Cell-type comparisons require paired donors.", call. = FALSE)
  }
  required <- c("readout", "celltype", "target", "reference",
                if (analysis_scope == "within_donor") "donor", if (axis == "celltype") "condition")
  for (field in required) {
    value <- request[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("Choose one ", field, ".", call. = FALSE)
    }
  }
  if (!request$readout %in% c("abundance", "self_proximity", "colocalization")) {
    stop("Unknown analysis readout.", call. = FALSE)
  }
  if (request$target == request$reference) stop("Choose two different ",
    if (axis == "condition") "conditions." else "cell types.", call. = FALSE)
  metadata <- data$metadata
  scope <- if (axis == "condition") metadata$celltype == request$celltype else metadata$condition == request$condition
  groups <- if (axis == "condition") metadata$condition else metadata$celltype
  if (analysis_scope == "within_donor") {
    scope <- scope & metadata$donor == request$donor
    if (!any(scope)) stop("No cells match this donor and cell type.", call. = FALSE)
  } else {
    donors <- request$donors
    if (is.null(donors)) donors <- sort(unique(metadata$donor[scope & groups %in% c(request$target, request$reference)]))
    if (!is.character(donors) || !length(donors) || anyNA(donors) || any(!nzchar(donors)) ||
        anyDuplicated(donors) || !all(donors %in% metadata$donor)) stop("Choose known, unique donors.", call. = FALSE)
    donors <- sort(donors, method = "radix")
    if (length(donors) > limits$donors) stop("Choose at most ", limits$donors, " donors.", call. = FALSE)
    scope <- scope & metadata$donor %in% donors
    donors_a <- sort(unique(metadata$donor[scope & groups == request$target]), method = "radix")
    donors_b <- sort(unique(metadata$donor[scope & groups == request$reference]), method = "radix")
    if (design == "independent" && length(intersect(donors_a, donors_b))) {
      stop("Independent groups cannot share donors. Use paired donors or select disjoint donor groups.", call. = FALSE)
    }
    used <- if (design == "paired") intersect(donors_a, donors_b) else union(donors_a, donors_b)
    excluded <- setdiff(donors, used)
    if (!length(used) || (design == "independent" && (!length(donors_a) || !length(donors_b)))) {
      stop("No eligible donor groups remain for this comparison.", call. = FALSE)
    }
    scope <- scope & metadata$donor %in% used
    donors_a <- intersect(donors_a, used)
    donors_b <- intersect(donors_b, used)
  }
  select_cells <- function(field, group) {
    available <- metadata$component[scope & groups == group]
    cells <- request[[field]]
    if (is.null(cells)) cells <- available
    if (!is.character(cells) || anyNA(cells) || anyDuplicated(cells) ||
        !all(cells %in% available)) {
      stop("Selected cells must be unique and belong to the selected donor, cell type, and condition.", call. = FALSE)
    }
    if (analysis_scope == "across_donors" && !setequal(cells, available)) {
      stop("Across-donor comparisons use all cells in the selected donor groups.", call. = FALSE)
    }
    sort(cells, method = "radix")
  }
  group_a_cells <- select_cells("group_a_cells", request$target)
  group_b_cells <- select_cells("group_b_cells", request$reference)
  min_cells <- switch(request$readout, abundance = 3L, self_proximity = 20L, colocalization = 40L)
  if (request$readout != "abundance" && !is.null(request$min_cells_per_group)) {
    min_cells <- request$min_cells_per_group
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1L || !is.finite(min_cells) ||
      min_cells < 1 || min_cells > limits$cells / 2 || min_cells != floor(min_cells)) {
    stop("Minimum observed cells per condition must be a whole number between 1 and ",
      floor(limits$cells / 2), ".", call. = FALSE)
  }
  min_cells <- as.integer(min_cells)
  if (analysis_scope == "within_donor" && min(length(group_a_cells), length(group_b_cells)) < min_cells) {
    stop("This readout requires at least ", min_cells, " selected cells in each condition.", call. = FALSE)
  }
  cells <- c(group_a_cells, group_b_cells)
  donor_cells <- table(metadata$donor[match(cells, metadata$component)])
  if (any(donor_cells > limits$cells)) {
    stop(if (analysis_scope == "across_donors") "Each donor must contain at most " else "The comparison must contain at most ",
      limits$cells, " cells in total. Use precomputed results or run this comparison offline.", call. = FALSE)
  }
  if (analysis_scope == "across_donors" && length(cells) > limits$total_cells) {
    stop("The comparison must contain at most ", limits$total_cells, " cells across donors.", call. = FALSE)
  }
  markers <- request$markers
  if (is.null(markers) || identical(markers, character())) markers <- rownames(data$abundance)
  if (!is.character(markers) || anyNA(markers) || !length(markers) ||
      !all(markers %in% rownames(data$abundance))) stop("Choose known markers.", call. = FALSE)
  markers <- sort(unique(markers), method = "radix")
  if (length(markers) > limits$markers) stop("Too many selected markers.", call. = FALSE)
  if (request$readout == "colocalization" && length(markers) < 2L) {
    stop("Colocalization requires at least two dataset markers.", call. = FALSE)
  }
  normalized <- list(dataset_id = data$contract$dataset_id, bundle_md5 = data$bundle_md5,
       readout = request$readout, donor = request$donor, celltype = request$celltype,
       target = request$target, reference = request$reference, markers = markers,
       group_a_cells = group_a_cells, group_b_cells = group_b_cells,
       min_cells = min_cells, min_cells_per_group = min_cells)
  if (analysis_scope == "across_donors" || any(c("analysis_scope", "comparison_axis", "design") %in% names(request))) {
    normalized$analysis_scope <- analysis_scope
    normalized$comparison_axis <- axis
    normalized$design <- design
    normalized$condition <- if (axis == "celltype") request$condition else NULL
    if (axis == "celltype") normalized$celltype <- request$target
  }
  if (analysis_scope == "across_donors") {
    normalized$donors <- donors
    normalized$donors_a <- donors_a
    normalized$donors_b <- donors_b
    normalized$excluded_donors <- data.frame(donor = excluded,
      reason = rep(if (design == "paired") "incomplete_pair" else "no_comparison_cells", length(excluded)))
  }
  normalized
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

run_donor_analysis <- function(data, request, asset_root = NULL, score_reader = read_proximity_subset) {
  request <- validate_hybrid_request(data, request)
  if (!identical(request$analysis_scope, "across_donors")) stop("Choose an across-donor request.", call. = FALSE)
  features <- hybrid_feature_table(request)
  keys <- paste(features$marker_1, features$marker_2, sep = "\r")
  axis <- request$comparison_axis
  donors <- sort(union(request$donors_a, request$donors_b), method = "radix")
  values <- list()
  filters <- list()
  for (donor in donors) {
    donor_cells <- data$metadata$component[data$metadata$donor == donor]
    donor_request <- request
    donor_request$donor <- donor
    groups <- list(a = intersect(request$group_a_cells, donor_cells), b = intersect(request$group_b_cells, donor_cells))
    donor_values <- lapply(names(groups), function(group) {
      label <- if (group == "a") request$target else request$reference
      rows <- features
      rows$donor <- donor
      rows$group <- label
      rows$value <- NA_real_
      rows$n_total_cells <- rows$n_evaluable_cells <- length(groups[[group]])
      rows$n_scored_cells <- 0L
      rows$eligible <- FALSE
      rows$source_celltype <- if (axis == "condition") request$celltype else label
      rows$condition <- if (axis == "condition") label else request$condition
      rows
    })
    names(donor_values) <- names(groups)
    if (request$readout == "abundance") {
      for (group in names(groups)) {
        cells <- groups[[group]]
        if (!length(cells)) next
        expression <- as.matrix(data$abundance[request$markers, cells, drop = FALSE])
        scored <- rowSums(is.finite(expression))
        expression[!is.finite(expression)] <- NA_real_
        eligible <- scored >= request$min_cells
        donor_values[[group]]$n_scored_cells <- as.integer(scored)
        donor_values[[group]]$eligible <- eligible
        donor_values[[group]]$value[eligible] <- rowMeans(expression[eligible, , drop = FALSE], na.rm = TRUE)
      }
    } else {
      cell_count <- sum(lengths(groups))
      batch_size <- max(1L, floor(min(HYBRID_LIMITS$pairs, HYBRID_LIMITS$rows / cell_count) /
        if (request$readout == "colocalization") length(request$markers) else 1L))
      batches <- split(request$markers, ceiling(seq_along(request$markers) / batch_size))
      for (batch in batches) {
        scores <- score_reader(data, asset_root, donor_request, marker_1 = batch)
        if (!all(c("component", "marker_1", "marker_2", "log2_ratio") %in% names(scores)) ||
            any(!scores$component %in% unlist(groups, use.names = FALSE)) ||
            anyDuplicated(scores[c("component", "marker_1", "marker_2")])) {
          stop("Spatial queries must return unique observations for the requested donor.", call. = FALSE)
        }
        provenance <- attr(scores, "provenance")
        if (!is.null(provenance)) {
          if (is.null(filters[[donor]])) filters[[donor]] <- provenance else {
            filters[[donor]]$observed_rows <- filters[[donor]]$observed_rows + provenance$observed_rows
          }
        }
        scores <- scores[is.finite(scores$log2_ratio), , drop = FALSE]
        score_keys <- paste(scores$marker_1, scores$marker_2, sep = "\r")
        for (group in names(groups)) {
          index <- scores$component %in% groups[[group]]
          observed <- split(scores$log2_ratio[index], score_keys[index])
          counts <- lengths(observed)
          positions <- match(names(counts), keys)
          if (anyNA(positions)) stop("Spatial queries returned an unexpected marker pair.", call. = FALSE)
          donor_values[[group]]$n_scored_cells[positions] <- as.integer(counts)
          eligible <- counts >= request$min_cells
          donor_values[[group]]$eligible[positions] <- eligible
          donor_values[[group]]$value[positions[eligible]] <- vapply(observed[eligible], stats::median, numeric(1))
        }
      }
    }
    values <- c(values, donor_values)
  }
  values <- do.call(rbind, values)
  rownames(values) <- NULL
  results <- features
  results$effect_size <- results$min_effect <- results$max_effect <- results$p_value <- results$p_adj <- NA_real_
  results$n_pairs <- results$n_donors_a <- results$n_donors_b <- results$n_cells_a <- results$n_cells_b <- 0L
  results$direction_consistent <- NA
  results$status <- "insufficient_donors"
  results$inference <- "descriptive"
  results$method <- if (request$design == "paired") "Paired donor t-test" else "Welch independent donor t-test"
  results$replicate_unit <- "donor"
  contrasts <- list()
  for (i in seq_len(nrow(features))) {
    rows <- values[values$feature_id == features$feature_id[[i]], , drop = FALSE]
    a <- rows[rows$group == request$target & rows$eligible, , drop = FALSE]
    b <- rows[rows$group == request$reference & rows$eligible, , drop = FALSE]
    results$n_donors_a[[i]] <- nrow(a)
    results$n_donors_b[[i]] <- nrow(b)
    results$n_cells_a[[i]] <- sum(a$n_total_cells)
    results$n_cells_b[[i]] <- sum(b$n_total_cells)
    if (request$design == "paired") {
      paired <- intersect(a$donor, b$donor)
      a <- a[match(paired, a$donor), , drop = FALSE]
      b <- b[match(paired, b$donor), , drop = FALSE]
      differences <- a$value - b$value
      results$n_pairs[[i]] <- length(paired)
      if (!length(paired)) next
      contrast <- features[rep(i, length(paired)), , drop = FALSE]
      contrast$donor <- paired
      contrast$value_a <- a$value
      contrast$value_b <- b$value
      contrast$effect_size <- differences
      contrast$n_cells_a <- a$n_total_cells
      contrast$n_cells_b <- b$n_total_cells
      contrasts[[length(contrasts) + 1L]] <- contrast
      results$effect_size[[i]] <- mean(differences)
      results$min_effect[[i]] <- min(differences)
      results$max_effect[[i]] <- max(differences)
      results$direction_consistent[[i]] <- length(unique(sign(differences))) == 1L
      enough <- length(paired) >= 2L
      tested <- if (enough) tryCatch(stats::t.test(differences, mu = 0), error = function(e) NULL) else NULL
    } else {
      if (!nrow(a) || !nrow(b)) next
      results$effect_size[[i]] <- mean(a$value) - mean(b$value)
      enough <- min(nrow(a), nrow(b)) >= 2L
      tested <- if (enough) tryCatch(stats::t.test(a$value, b$value), error = function(e) NULL) else NULL
    }
    if (enough) results$status[[i]] <- "no_finite_test_result"
    if (!is.null(tested) && is.finite(tested$p.value)) {
      results$p_value[[i]] <- tested$p.value
      results$status[[i]] <- "tested"
      results$inference[[i]] <- "donor_level"
    }
  }
  tested <- is.finite(results$p_value)
  results$p_adj[tested] <- stats::p.adjust(results$p_value[tested], method = "BH")
  if (length(contrasts)) {
    contrasts <- do.call(rbind, contrasts)
    rownames(contrasts) <- NULL
  } else {
    contrasts <- features[FALSE, , drop = FALSE]
    contrasts$donor <- character()
    contrasts$value_a <- contrasts$value_b <- contrasts$effect_size <- numeric()
    contrasts$n_cells_a <- contrasts$n_cells_b <- integer()
  }
  context <- function(rows) {
    rows$readout <- rep(request$readout, nrow(rows))
    rows$celltype <- rep(request$celltype, nrow(rows))
    rows$group_a <- rep(request$target, nrow(rows))
    rows$group_b <- rep(request$reference, nrow(rows))
    rows$comparison_axis <- rep(axis, nrow(rows))
    rows$design <- rep(request$design, nrow(rows))
    rows$fixed_condition <- rep(if (axis == "celltype") request$condition else NA_character_, nrow(rows))
    rows
  }
  list(results = context(results), donor_values = context(values), donor_contrasts = context(contrasts),
       provenance = list(dataset_id = data$contract$dataset_id, bundle_md5 = data$bundle_md5,
         analysis_version = HYBRID_ANALYSIS_VERSION, request = request,
         created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
         method = unique(results$method), p_adjust_method = "BH",
         p_adjust_scope = "eligible tested features within this donor-level request",
         abundance_layer = data$contract$defaults$abundance_layer,
         requested_features = nrow(results), tested_features = sum(tested),
         summary = if (request$readout == "abundance") "Mean stored normalized abundance per donor and group" else
           "Median retained observed finite log2 ratio per donor and group",
         sample_aggregation = "All cells from samples of the same donor and group are pooled before summarization; each donor has equal weight.",
         missing_score_policy = "Absent or filtered scores are not zeros; feature-level incomplete donor pairs are excluded from paired tests.",
         inference = "Donors are biological replicates; paired tests use complete eligible donor pairs only.",
         exclusions = request$excluded_donors, filters = filters,
         software = list(R = R.version.string)))
}

run_hybrid_analysis <- function(bundle_path, request, asset_root = NULL) {
  data <- read_precomputed_rds(bundle_path)
  request <- validate_hybrid_request(data, request)
  if (identical(request$analysis_scope, "across_donors")) {
    if (request$readout != "abundance" && is.null(asset_root)) {
      asset_root <- resolve_hybrid_asset(bundle_path, data$assets$proximity)
    }
    return(run_donor_analysis(data, request, asset_root, score_reader = read_proximity_subset))
  }
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
    # Partition by the first canonical marker; each pair belongs to exactly one batch.
    batch_size <- max(1L, floor(min(HYBRID_LIMITS$pairs, HYBRID_LIMITS$rows / length(cells)) /
      if (request$readout == "colocalization") length(request$markers) else 1L))
    batches <- split(request$markers, ceiling(seq_along(request$markers) / batch_size))
    result_keys <- paste(results$marker_1, results$marker_2, sep = "\r")
    for (batch in batches) {
      scores <- read_proximity_subset(data, asset_root, request, marker_1 = batch)
      provenance <- attr(scores, "provenance")
      if (is.null(filter_provenance)) {
        filter_provenance <- provenance
      } else {
        filter_provenance$observed_rows <- filter_provenance$observed_rows + provenance$observed_rows
      }
      scores$condition <- ifelse(scores$component %in% request$group_a_cells, request$target, request$reference)
      score_keys <- paste(scores$marker_1, scores$marker_2, sep = "\r")
      in_batch <- results$marker_1 %in% batch
      for (group in c("a", "b")) {
        counts <- table(score_keys[scores$condition == if (group == "a") request$target else request$reference])
        values <- as.integer(counts[result_keys])
        values[is.na(values)] <- 0L
        results[[paste0("n_scored_", group)]][in_batch] <- values[in_batch]
      }
      eligible <- in_batch & results$n_scored_a >= request$min_cells & results$n_scored_b >= request$min_cells
      if (any(eligible)) {
        tested <- pixelatorR::DifferentialProximityAnalysis(
          scores, contrast_column = "condition", targets = request$target, reference = request$reference,
          proximity_metric = "log2_ratio", metric_type = if (request$readout == "self_proximity") "self" else "co",
          p_adjust_method = "BH", min_cells_per_group = request$min_cells, verbose = FALSE)
        index <- match(paste(tested$marker_1, tested$marker_2, sep = "\r"), result_keys)
        results$effect_size[index] <- tested$diff_median
        results$p_value[index] <- tested$p
      }
    }
    tested <- is.finite(results$p_value)
    results$p_adj[tested] <- stats::p.adjust(results$p_value[tested], method = "BH")
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
      Seurat = if (requireNamespace("Seurat", quietly = TRUE)) {
        as.character(utils::packageVersion("Seurat"))
      } else {
        NA_character_
      }))
  )
}
