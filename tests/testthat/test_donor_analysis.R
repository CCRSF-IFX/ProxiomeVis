donor_analysis_fixture <- function() {
  metadata <- expand.grid(donor = paste0("D", 1:4), condition = c("NC", "TRT"),
                          celltype = c("B", "T"), cell = 1:4, stringsAsFactors = FALSE)
  metadata$component <- paste0("donor-cell-", seq_len(nrow(metadata)))
  metadata$sample_id <- paste(metadata$donor, metadata$condition, sep = "_")
  donor <- match(metadata$donor, paste0("D", 1:4))
  value <- donor * 10 + (metadata$condition == "TRT") * c(1, 2, 4, 7)[donor] +
    (metadata$celltype == "T") * c(2, 4, 5, 9)[donor]
  abundance <- rbind(A = value, B = value / 2)
  colnames(abundance) <- metadata$component
  list(metadata = metadata, abundance = abundance,
       contract = list(dataset_id = "donor-fixture", defaults = list(abundance_layer = "data")),
       bundle_md5 = "donor-fixture-hash")
}

donor_analysis_request <- function(design = "paired", axis = "condition", readout = "abundance") {
  list(analysis_scope = "across_donors", comparison_axis = axis, design = design,
       donors = paste0("D", 1:4), readout = readout, markers = c("B", "A"),
       celltype = "B", condition = "NC",
       target = if (axis == "condition") "TRT" else "T",
       reference = if (axis == "condition") "NC" else "B", min_cells_per_group = 3L)
}

test_that("paired condition comparisons give equal weight to donors rather than cells", {
  data <- donor_analysis_fixture()
  request <- donor_analysis_request()
  # Extra technical-sample cells must not become extra biological replicates.
  extra <- data$metadata[data$metadata$donor == "D4", , drop = FALSE]
  extra$component <- paste0(extra$component, "-technical")
  extra$sample_id <- paste0(extra$sample_id, "-technical")
  expression <- data$abundance[, sub("-technical$", "", extra$component), drop = FALSE]
  colnames(expression) <- extra$component
  data$metadata <- rbind(data$metadata, extra)
  data$abundance <- cbind(data$abundance, expression)
  result <- run_donor_analysis(data, request)
  expect_equal(result$results$effect_size, c(mean(c(1, 2, 4, 7)), mean(c(1, 2, 4, 7)) / 2))
  expect_equal(result$results$p_value, rep(t.test(c(1, 2, 4, 7))$p.value, 2))
  expect_equal(result$results$p_adj, p.adjust(result$results$p_value, "BH"))
  expect_identical(result$results$n_pairs, c(4L, 4L))
  expect_identical(result$results$n_donors_a, c(4L, 4L))
  expect_true(all(result$results$replicate_unit == "donor"))
  expect_true(all(result$results$inference == "donor_level"))
  expect_equal(result$results$min_effect, c(1, .5))
  expect_equal(result$results$max_effect, c(7, 3.5))
  expect_equal(nrow(result$donor_values), 16)
  expect_equal(nrow(result$donor_contrasts), 8)
  expect_match(result$provenance$sample_aggregation, "equal weight")
  expect_identical(result$provenance$request, validate_hybrid_request(data, result$provenance$request))
})

test_that("cell-type comparisons hold condition fixed and compare matching donors", {
  data <- donor_analysis_fixture()
  request <- donor_analysis_request(axis = "celltype")
  result <- run_donor_analysis(data, request)
  expect_equal(result$results$effect_size, c(5, 2.5))
  expect_equal(result$results$p_value, rep(t.test(c(2, 4, 5, 9))$p.value, 2))
  expect_true(all(result$donor_values$condition == "NC"))
  expect_setequal(result$donor_values$source_celltype, c("B", "T"))
  expect_true(all(result$donor_values$celltype == "T"))
  expect_true(all(result$results$fixed_condition == "NC"))
  request$design <- "independent"
  expect_error(validate_hybrid_request(data, request), "paired donors")
  within <- donor_analysis_request(axis = "celltype")
  within$analysis_scope <- within$design <- "within_donor"
  within$donor <- "D1"
  selected <- validate_hybrid_request(data, within)
  expect_length(selected$group_a_cells, 4L)
  expect_length(selected$group_b_cells, 4L)
  expect_true(all(data$metadata$condition[match(c(selected$group_a_cells, selected$group_b_cells),
                                               data$metadata$component)] == "NC"))
})

test_that("independent treatment groups use donor means and Welch tests", {
  data <- donor_analysis_fixture()
  data$metadata <- data$metadata[(data$metadata$donor %in% c("D1", "D2") & data$metadata$condition == "NC") |
    (data$metadata$donor %in% c("D3", "D4") & data$metadata$condition == "TRT"), ]
  result <- run_donor_analysis(data, donor_analysis_request("independent"))
  expected <- t.test(c(34, 47), c(10, 20))
  expect_equal(result$results$effect_size, c(25.5, 12.75))
  expect_equal(result$results$p_value, rep(expected$p.value, 2))
  expect_true(all(result$results$n_pairs == 0L))
  expect_true(all(result$results$n_donors_a == 2L & result$results$n_donors_b == 2L))
  expect_true(all(is.na(result$results$min_effect) & is.na(result$results$max_effect)))
  expect_equal(nrow(result$donor_contrasts), 0L)
  absent <- result$donor_values[result$donor_values$n_total_cells == 0, ]
  expect_false(any(absent$eligible))
  expect_true(all(is.na(absent$value)))
  expect_error(validate_hybrid_request(donor_analysis_fixture(), donor_analysis_request("independent")), "cannot share donors")
})

test_that("incomplete donor pairs are explicit exclusions and selection is pinned", {
  data <- donor_analysis_fixture()
  data$metadata <- data$metadata[!(data$metadata$donor == "D4" & data$metadata$condition == "TRT"), ]
  request <- validate_hybrid_request(data, donor_analysis_request())
  expect_identical(request$excluded_donors, data.frame(donor = "D4", reason = "incomplete_pair"))
  expect_identical(request$donors_a, paste0("D", 1:3))
  expect_identical(request$donors_b, paste0("D", 1:3))
  expect_false(any(data$metadata$donor[match(c(request$group_a_cells, request$group_b_cells),
                                           data$metadata$component)] == "D4"))
  expect_identical(validate_hybrid_request(data, request), request)
  request$group_a_cells <- request$group_a_cells[-1]
  expect_error(validate_hybrid_request(data, request), "all cells")
  request <- donor_analysis_request()
  request$dataset_id <- "changed"
  expect_error(validate_hybrid_request(data, request), "dataset changed")
  request <- donor_analysis_request()
  request$donors <- c("D1", "D1")
  expect_error(validate_hybrid_request(data, request), "unique donors")
  expect_error(validate_hybrid_request(data, donor_analysis_request(), modifyList(HYBRID_LIMITS, list(donors = 2))), "at most 2")
})

test_that("missing feature support and insufficient donors give descriptive effects only", {
  data <- donor_analysis_fixture()
  request <- donor_analysis_request()
  request$donors <- "D1"
  one <- run_donor_analysis(data, request)
  expect_equal(one$results$effect_size, c(1, .5))
  expect_true(all(is.na(one$results$p_value) & is.na(one$results$p_adj)))
  expect_true(all(one$results$inference == "descriptive"))
  expect_true(all(one$results$status == "insufficient_donors"))
  # A missing measurement is not a zero and cannot supply a complete pair.
  missing <- data$metadata$condition == "TRT" & data$metadata$donor != "D1"
  data$abundance["B", data$metadata$component[missing]] <- NA_real_
  result <- run_donor_analysis(data, donor_analysis_request())
  expect_identical(result$results$n_pairs, c(4L, 1L))
  expect_true(is.finite(result$results$p_value[1]))
  expect_true(is.na(result$results$p_value[2]))
  expect_equal(result$results$effect_size[2], .5)
  expect_true(all(is.na(result$donor_values$value[!result$donor_values$eligible])))
  # Degenerate donor differences must not crash the whole comparison.
  data$abundance[,] <- 1
  constant <- run_donor_analysis(data, donor_analysis_request())
  expect_true(all(is.na(constant$results$p_value)))
  expect_true(all(constant$results$status == "no_finite_test_result"))
})

test_that("spatial donor summaries stream one donor and retain observed-score missingness", {
  data <- donor_analysis_fixture()
  request <- donor_analysis_request(readout = "self_proximity")
  grid <- merge(data$metadata[c("component", "donor", "condition", "celltype")],
                data.frame(marker_1 = c("A", "B")), by = NULL)
  grid$marker_2 <- grid$marker_1
  grid$log2_ratio <- data$abundance[cbind(match(grid$marker_1, rownames(data$abundance)),
                                          match(grid$component, colnames(data$abundance)))]
  grid <- grid[!(grid$marker_1 == "B" & grid$donor == "D2" & grid$condition == "TRT"), ]
  calls <- character()
  reader <- function(data, asset_root, request, marker_1) {
    calls <<- c(calls, request$donor)
    cells <- data$metadata$component[data$metadata$donor == request$donor]
    cells <- intersect(cells, c(request$group_a_cells, request$group_b_cells))
    rows <- grid[grid$component %in% cells & grid$marker_1 %in% marker_1,
                 c("component", "marker_1", "marker_2", "log2_ratio")]
    attr(rows, "provenance") <- list(observed_rows = nrow(rows), background_threshold_pct = .001)
    rows
  }
  result <- run_donor_analysis(data, request, score_reader = reader)
  expect_identical(calls, paste0("D", 1:4))
  expect_equal(result$results$effect_size, c(3.5, mean(c(.5, 2, 3.5))))
  expect_identical(result$results$n_pairs, c(4L, 3L))
  expect_equal(result$results$p_value[2], t.test(c(.5, 2, 3.5))$p.value)
  missing <- subset(result$donor_values, donor == "D2" & group == "TRT" & feature_id == "B")
  expect_false(missing$eligible)
  expect_true(is.na(missing$value))
  expect_equal(missing$n_total_cells, 4L)
  expect_equal(missing$n_scored_cells, 0L)
  expect_length(result$provenance$filters, 4L)
  expect_match(result$provenance$summary, "Median retained")
  bad <- function(data, asset_root, request, marker_1) grid[, c("component", "marker_1", "marker_2", "log2_ratio")]
  expect_error(run_donor_analysis(data, request, score_reader = bad), "requested donor")
})

test_that("donor summaries are independent of metadata and marker input order", {
  data <- donor_analysis_fixture()
  expected <- run_donor_analysis(data, donor_analysis_request())
  data$metadata <- data$metadata[rev(seq_len(nrow(data$metadata))), ]
  request <- donor_analysis_request()
  request$markers <- rev(request$markers)
  request$donors <- rev(request$donors)
  result <- run_donor_analysis(data, request)
  expect_equal(result$results, expected$results)
  expect_equal(result$donor_values, expected$donor_values)
  expect_equal(result$donor_contrasts, expected$donor_contrasts)
})

test_that("the application worker uses the same donor engine as offline preparation", {
  data <- donor_analysis_fixture()
  path <- tempfile(fileext = ".rds")
  saveRDS(data, path)
  runner <- run_hybrid_analysis
  environment(runner) <- list2env(list(read_precomputed_rds = function(path) readRDS(path)),
                                  parent = environment(run_hybrid_analysis))
  result <- runner(path, donor_analysis_request())
  expected <- run_donor_analysis(data, donor_analysis_request())
  expect_equal(result$results, expected$results)
  expect_equal(result$donor_values, expected$donor_values)
  expect_equal(result$donor_contrasts, expected$donor_contrasts)
})

test_that("within-donor cell-type testing does not replace the fixed source condition", {
  skip_if_not_installed("pixelatorR")
  data <- donor_analysis_fixture()
  data$counts <- Matrix::Matrix(round(data$abundance), sparse = TRUE)
  original <- data$metadata
  path <- tempfile(fileext = ".rds")
  saveRDS(data, path)
  runner <- run_hybrid_analysis
  environment(runner) <- list2env(list(read_precomputed_rds = function(path) readRDS(path)),
                                  parent = environment(run_hybrid_analysis))
  request <- donor_analysis_request(axis = "celltype")
  request$analysis_scope <- request$design <- "within_donor"
  request$donor <- "D1"
  result <- runner(path, request)
  expect_equal(result$results$effect_size, c(2, 1))
  expect_equal(result$results$n_cells_a, c(4L, 4L))
  expect_equal(result$results$n_cells_b, c(4L, 4L))
  expect_true(all(result$results$inference == "exploratory_component_level"))
  expect_equal(result$provenance$request$condition, "NC")
  expect_identical(readRDS(path)$metadata, original)
  request$readout <- "self_proximity"
  request$min_cells_per_group <- 5L
  expect_error(runner(path, request), "at least 5")
})

test_that("real spatial partitions keep donor cutoffs and exact cell IDs for both comparison axes", {
  skip_if_not_installed("arrow")
  for (axis in c("condition", "celltype")) {
    data <- donor_analysis_fixture()
    selected <- data$metadata$donor %in% c("D1", "D2")
    if (axis == "condition") {
      selected <- selected & data$metadata$celltype == "B" &
        ((data$metadata$donor == "D1" & data$metadata$condition == "NC") |
         (data$metadata$donor == "D2" & data$metadata$condition == "TRT"))
    } else selected <- selected & data$metadata$condition == "NC"
    data$metadata <- data$metadata[selected, , drop = FALSE]
    data$abundance <- data$abundance[, data$metadata$component, drop = FALSE]
    data$cell_sources <- data.frame(component = data$metadata$component,
      original_id = rev(data$metadata$component), sample_id = data$metadata$sample_id)
    root <- withr::local_tempdir()
    raw <- merge(data$metadata, data.frame(marker_1 = "A", marker_2 = c("A", "B")), by = NULL)
    raw$log2_ratio <- data$abundance["A", raw$component] + (raw$marker_2 == "B") / 10
    raw$min_marker_fraction <- c(.02, .05, .08, .10)[raw$cell]
    raw <- raw[c("component", "marker_1", "marker_2", "log2_ratio", "min_marker_fraction")]
    samples <- sort(unique(data$metadata$sample_id))
    partitions <- lapply(seq_along(samples), function(i) {
      sample <- samples[[i]]
      file <- paste0("sample-", i, ".parquet")
      cells <- data$metadata$component[data$metadata$sample_id == sample]
      rows <- raw[raw$component %in% cells, , drop = FALSE]
      path <- file.path(root, file)
      arrow::write_parquet(rows, path, chunk_size = 2L)
      data.frame(sample_id = sample, file = file, rows = nrow(rows), md5 = unname(tools::md5sum(path)))
    })
    manifest <- list(schema_version = 1L, dataset_id = data$contract$dataset_id,
      dataset_fingerprint = hybrid_dataset_fingerprint(data), partitions = do.call(rbind, partitions),
      donor_cutoffs = c(D1 = .05, D2 = .08))
    saveRDS(manifest, file.path(root, "manifest.rds"))
    expect_length(list.files(root, pattern = "[.]parquet$"), 2L)
    for (readout in c("self_proximity", "colocalization")) {
      request <- donor_analysis_request(design = if (axis == "condition") "independent" else "paired",
        axis = axis, readout = readout)
      request$donors <- c("D1", "D2")
      request$min_cells_per_group <- 2L
      validated <- validate_hybrid_request(data, request)
      expect_error(read_proximity_subset(data, root, validated), "one selected donor")
      for (donor in request$donors) {
        query <- validated
        query$donor <- donor
        scores <- read_proximity_subset(data, root, query)
        wanted <- data$metadata$component[data$metadata$donor == donor]
        expected <- raw[raw$component %in% wanted &
          raw$min_marker_fraction >= manifest$donor_cutoffs[[donor]] &
          if (readout == "self_proximity") raw$marker_1 == raw$marker_2 else raw$marker_1 != raw$marker_2, ]
        expect_setequal(scores$component, expected$component)
        expected_key <- paste(expected$component, expected$marker_1, expected$marker_2)
        actual_key <- paste(scores$component, scores$marker_1, scores$marker_2)
        expect_equal(scores$log2_ratio, expected$log2_ratio[match(actual_key, expected_key)])
        expect_equal(attr(scores, "provenance")$background_threshold_pct, manifest$donor_cutoffs[[donor]])
        expect_true(all(data$metadata$donor[match(scores$component, data$metadata$component)] == donor))
      }
      result <- run_donor_analysis(data, request, asset_root = root)
      expect_equal(result$provenance$filters$D1$background_threshold_pct, .05)
      expect_equal(result$provenance$filters$D2$background_threshold_pct, .08)
      measured <- result$donor_values[result$donor_values$marker_1 == "A" &
        result$donor_values$n_total_cells > 0, ]
      expect_equal(measured$n_scored_cells, ifelse(measured$donor == "D1", 3L, 2L))
      if (axis == "celltype") {
        expect_equal(result$results$n_pairs[[1]], 2L)
        expect_equal(result$results$effect_size[[1]], 3)
        request$min_cells_per_group <- 3L
        insufficient <- run_donor_analysis(data, request, asset_root = root)
        expect_equal(insufficient$results$n_pairs[[1]], 1L)
        expect_true(is.na(insufficient$results$p_value[[1]]))
        expect_equal(insufficient$results$effect_size[[1]], 2)
      }
    }
  }
})
