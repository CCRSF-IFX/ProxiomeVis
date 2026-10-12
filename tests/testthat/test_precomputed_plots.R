test_that("donor views preserve stored exploratory statistics", {
  contract <- precomputed_fixture()@misc$proxiomevis
  rows <- contract$donor_contrasts$abundance
  rows$p_value <- seq_len(nrow(rows)) / 100
  rows$p_adj <- seq_len(nrow(rows)) / 50
  rows$adjustment_family <- "stored-family"
  rows$method <- "stored component test"
  rows$replicate_unit <- "component"
  rows$inference <- "exploratory_component_level"
  contract$donor_contrasts$abundance <- rows

  selected <- precomputed_result_rows_for_donor(contract, "abundance", "D1")
  expect_true(nrow(selected) > 0L)
  expect_true(all(selected$donor == "D1"))
  expect_equal(selected$p_value, rows$p_value[rows$donor == "D1"])
  expect_equal(selected$p_adj, rows$p_adj[rows$donor == "D1"])
  expect_identical(unique(selected$replicate_unit), "component")
  expect_identical(unique(selected$inference), "exploratory_component_level")
})

test_that("heatmap completion keeps unavailable values missing and uses symmetric limits", {
  results <- precomputed_fixture()@misc$proxiomevis$differential_results$abundance
  second <- results
  second$analysis_id <- "abundance:ALT_vs_NC"
  second$comparison <- "ALT vs NC"
  second$group_a <- "TRT"
  second <- second[-1L, , drop = FALSE]
  combined <- rbind(results, second)

  prepared <- precomputed_prepare_heatmap(combined, top_n = 20L)
  expect_equal(nrow(prepared), length(unique(combined$feature_id)) * 4L)
  expect_true(anyNA(prepared$effect_size))

  plot <- plot_precomputed_heatmap(combined)
  expect_s3_class(plot, "ggplot")
  limit <- max(abs(combined$effect_size), na.rm = TRUE)
  expect_equal(plot$scales$scales[[1L]]$limits, c(-limit, limit))
})

test_that("embedding helpers align cells and support sparse marker abundance", {
  app_data <- precomputed_app_data(precomputed_fixture())
  prepared <- precomputed_embedding_data(
    app_data$metadata,
    app_data$embeddings$umap,
    color_by = "abundance",
    abundance = app_data$abundance,
    marker = "A"
  )
  expect_identical(prepared$component, app_data$metadata$component)
  expect_true(is.numeric(prepared$color_value))
  expect_s3_class(
    plot_precomputed_embedding(app_data$metadata, app_data$embeddings$umap, "celltype"),
    "ggplot"
  )
})

test_that("QC retention and unavailable paired source values have explicit plots", {
  retention <- data.frame(
    stage = rep(c("00 loaded", "01 retained"), 2),
    sample_id = rep(c("D1_unstim", "D1_PHA"), each = 2),
    donor = "D1",
    condition = rep(c("unstim", "PHA"), each = 2),
    n_cells_display = c(100, 80, 120, 90),
    celltype = "All",
    stringsAsFactors = FALSE
  )
  expect_s3_class(plot_precomputed_qc_counts(retention), "ggplot")

  contrast <- precomputed_fixture()@misc$proxiomevis$donor_contrasts$colocalization[1, ]
  contrast$value_a <- NA_real_
  contrast$value_b <- NA_real_
  plot <- plot_precomputed_condition_detail(contrast)
  expect_s3_class(plot, "ggplot")
  expect_match(ggplot2::ggplot_build(plot)$data[[1L]]$label, "not stored")
})

test_that("the report UI exposes frozen report adapters", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  ui <- as.character(precomputed_ui("report"))
  expect_match(ui, "Report volcano", fixed = TRUE)
  expect_match(ui, "Spatial report", fixed = TRUE)
  expect_match(ui, "Exploratory component-level test", fixed = TRUE)
})

test_that("report cell-type choices retain source scopes that were not tested", {
  object <- precomputed_fixture()
  contract <- object@misc$proxiomevis
  contract$report <- list(
    eligibility = list(
      D1 = list(abundance = data.frame(celltype = c("B", "Not tested")))
    )
  )
  donor_rows <- contract$donor_contrasts$abundance[
    contract$donor_contrasts$abundance$donor == "D1",
    ,
    drop = FALSE
  ]
  choices <- precomputed_report_celltypes(contract, object@meta.data, "D1", "abundance", donor_rows)
  expect_identical(choices, c("B", "Not tested"))
})
