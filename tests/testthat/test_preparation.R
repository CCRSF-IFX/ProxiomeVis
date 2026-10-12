builder_path <- testthat::test_path("..", "..", "scripts", "prepare_precomputed_rds.R")
source(builder_path)

expect_all_source_rows_preserved <- function(contract, donor, readout) {
  source <- contract$report$source_results[[donor]][[readout]]
  stored <- contract$donor_contrasts[[readout]]
  stored <- stored[stored$donor == donor, , drop = FALSE]

  if (readout == "abundance") {
    marker_1 <- marker_2 <- as.character(source$marker)
    effect <- as.numeric(source$difference)
    n_cells_a <- as.numeric(source$n_cells_target)
    n_cells_b <- as.numeric(source$n_cells_reference)
  } else {
    marker_1 <- pmin(as.character(source$marker_1), as.character(source$marker_2))
    marker_2 <- pmax(as.character(source$marker_1), as.character(source$marker_2))
    effect <- as.numeric(source$diff_median)
    n_cells_a <- as.numeric(source$n_tgt)
    n_cells_b <- as.numeric(source$n_ref)
  }
  feature_id <- ifelse(marker_1 == marker_2, marker_1, paste(marker_1, marker_2, sep = " × "))
  analysis_id <- paste(readout, source$contrast, sep = ":")
  source_key <- paste(
    analysis_id, donor, source$celltype, feature_id, source$target, source$reference,
    sep = "\r"
  )
  stored_key <- paste(
    stored$analysis_id, stored$donor, stored$celltype, stored$feature_id,
    stored$group_a, stored$group_b, sep = "\r"
  )
  matched <- match(source_key, stored_key)

  expect_identical(nrow(stored), nrow(source), info = paste(donor, readout, "row count"))
  expect_identical(anyDuplicated(source_key), 0L, info = paste(donor, readout, "source keys"))
  expect_identical(anyDuplicated(stored_key), 0L, info = paste(donor, readout, "stored keys"))
  expect_false(anyNA(matched), info = paste(donor, readout, "key coverage"))
  stored <- stored[matched, , drop = FALSE]

  expect_identical(stored$marker_1, marker_1, info = paste(donor, readout, "marker_1"))
  expect_identical(stored$marker_2, marker_2, info = paste(donor, readout, "marker_2"))
  expect_identical(stored$feature_id, feature_id, info = paste(donor, readout, "feature_id"))
  expect_equal(stored$effect_size, effect, tolerance = 0, info = paste(donor, readout, "effect"))
  expect_equal(stored$p_value, as.numeric(source$p), tolerance = 0, info = paste(donor, readout, "p"))
  expect_equal(stored$p_adj, as.numeric(source$p_adj), tolerance = 0, info = paste(donor, readout, "p_adj"))
  expect_equal(stored$n_cells_a, n_cells_a, tolerance = 0, info = paste(donor, readout, "n_cells_a"))
  expect_equal(stored$n_cells_b, n_cells_b, tolerance = 0, info = paste(donor, readout, "n_cells_b"))
  expect_identical(stored$group_a, as.character(source$target), info = paste(donor, readout, "group_a"))
  expect_identical(stored$group_b, as.character(source$reference), info = paste(donor, readout, "group_b"))
  expect_identical(
    stored$adjustment_family,
    paste(donor, readout, source$contrast, source$celltype, "BH", sep = ":"),
    info = paste(donor, readout, "adjustment family")
  )
}

test_that("demo mappings preserve source effects, p values and pair identity", {
  abundance_samples <- data.frame(
    donor = "d1", condition = c("PHA", "unstim"), celltype = "T", feature_id = "A",
    value = c(3, 1), stringsAsFactors = FALSE
  )
  source <- data.frame(
    donor = c("d1", "d1"), layer = c("abundance", "colocalization"),
    contrast = "PHA_vs_unstim", target = "PHA", reference = "unstim", celltype = "T",
    # Frozen cross-donor co tables store marker 1/2/pair across feature/marker_1/marker_2.
    feature = c("A", "A"), marker_1 = c("A", "B"), marker_2 = c("", "A/B"),
    effect = c(2, 0.5), p = c(0.01, 0.02), p_adj = c(0.03, 0.04),
    n_target = c(10, 8), n_reference = c(9, 7), stringsAsFactors = FALSE
  )
  mapped <- make_donor_contrasts(source, c("A", "B"), abundance_samples)

  expect_equal(mapped$value_a[1] - mapped$value_b[1], mapped$effect_size[1])
  expect_true(all(is.na(mapped[2, c("value_a", "value_b")])))
  expect_identical(mapped$feature_id[2], "A × B")
  expect_equal(mapped$p_value, source$p)
  expect_equal(mapped$p_adj, source$p_adj)
  expect_identical(mapped$adjustment_family[2], "d1:colocalization:PHA_vs_unstim:T:BH")
  expect_identical(mapped$replicate_unit, rep("component", 2))
})

test_that("spatial sample summaries preserve the missing-as-zero population mean policy", {
  source <- data.frame(
    celltype = "T", condition = "PHA", marker_1 = "B", marker_2 = "A",
    n_cells_detected = 7L, n_cells = 10L, mean_log2_ratio = 0.25,
    stringsAsFactors = FALSE
  )
  mapped <- make_spatial_samples(source, "d1", "colocalization", c("A", "B"))

  expect_identical(mapped$feature_id, "A × B")
  expect_equal(mapped$n_total_cells, 10)
  expect_equal(mapped$n_evaluable_cells, 10)
  expect_equal(mapped$n_scored_cells, 7)
  expect_identical(mapped$metric, "population_mean_log2_ratio_missing_zero")
})

test_that("combined results remain donor-descriptive", {
  source <- data.frame(
    layer = "colocalization", contrast = "PHA_vs_unstim", target = "PHA", reference = "unstim",
    celltype = "T", feature = "A", marker_1 = "B", marker_2 = "A/B",
    effect_d1 = 0.2, effect_d2 = 0.4, tested_d1 = TRUE, tested_d2 = TRUE,
    same_direction = TRUE, support_status = "tested_both", stringsAsFactors = FALSE
  )
  mapped <- make_combined_results(source, c("A", "B"))

  expect_identical(mapped$feature_id, "A × B")
  expect_equal(mapped$effect_size, 0.3)
  expect_equal(mapped$min_effect, 0.2)
  expect_equal(mapped$max_effect, 0.4)
  expect_equal(mapped$n_pairs, 2)
  expect_true(mapped$direction_consistent)
  expect_true(is.na(mapped$p_value) && is.na(mapped$p_adj))
  expect_identical(mapped$replicate_unit, "donor")
  expect_identical(mapped$inference, "descriptive")
})

test_that("the builder refuses every output beneath the source analysis tree", {
  source_root <- tempfile("demo-source-")
  dir.create(file.path(source_root, "tertiary_analysis"), recursive = TRUE)
  input <- tempfile(fileext = ".rds")
  saveRDS(list(), input)
  on.exit(unlink(c(source_root, input), recursive = TRUE), add = TRUE)

  expect_error(
    build_precomputed_rds(input, file.path(source_root, "tertiary_analysis"),
                          file.path(source_root, "would-overwrite-source.rds")),
    "outside the demo_analysis source tree"
  )
})

test_that("the real demo builds atomically and preserves every report statistic", {
  repository <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  demo <- file.path(repository, "..", "pixelgen_demo_data", "demo_analysis")
  input <- file.path(demo, "cca_analysis", "results", "demo_cca_clustered.rds")
  analysis <- file.path(demo, "tertiary_analysis")
  testthat::skip_if_not(file.exists(input) && dir.exists(analysis), "demo bundle is not mounted")

  output <- tempfile(fileext = ".rds")
  on.exit(unlink(output), add = TRUE)
  input_md5 <- unname(tools::md5sum(input))
  expect_identical(build_precomputed_rds(input, analysis, output), normalizePath(output))
  expect_identical(unname(tools::md5sum(input)), input_md5)

  object <- readRDS(output)
  expect_true(validate_precomputed_object(object))
  contract <- object@misc$proxiomevis
  expect_match(contract$dataset_id, "^pixelgen-demo-v1-[0-9a-f]{16}$")
  expect_identical(names(contract$report$source_results), c("d1", "d2"))
  expect_identical(names(contract$report$source_results$d1), DEMO_READOUTS)
  expect_identical(contract$report$status$analysis$status, "PASS")
  expect_identical(contract$report$status$validation$status, "PASS")
  expect_true(all(vapply(object@assays, function(assay) {
    !"proximity" %in% methods::slotNames(assay) || NROW(methods::slot(assay, "proximity")) == 0L
  }, logical(1))))

  for (donor in DEMO_DONORS) {
    for (readout in DEMO_READOUTS) {
      expect_all_source_rows_preserved(contract, donor, readout)
    }
  }
  expect_true(all(is.na(contract$donor_contrasts$self_proximity$value_a)))
  expect_true(all(is.na(contract$donor_contrasts$self_proximity$value_b)))
  expect_true(all(is.na(contract$donor_contrasts$colocalization$value_a)))
  expect_true(all(is.na(contract$donor_contrasts$colocalization$value_b)))

  expect_true(all(vapply(contract$differential_results, function(x) {
    all(is.na(x$p_value)) && all(is.na(x$p_adj)) && all(x$inference == "descriptive")
  }, logical(1))))
  expect_true(all(vapply(contract$donor_contrasts, function(x) {
    all(is.finite(x$p_value)) && all(is.finite(x$p_adj))
  }, logical(1))))

  heatmaps <- contract$report$abundance_heatmaps$d1
  expect_identical(names(heatmaps$matrices), "PHA_vs_unstim")
  expect_identical(colnames(heatmaps$matrices$PHA_vs_unstim$B), heatmaps$samples$sample)
  expect_equal(
    range(heatmaps$manifest[c("color_min", "color_max")]),
    c(unique(heatmaps$manifest$color_min), unique(heatmaps$manifest$color_max))
  )
  expect_true(any(grepl("abundance_heatmap_matrix", contract$provenance$source_table_manifest$role)))
  expect_true("input_provenance" %in% contract$provenance$source_table_manifest$role)
})
