
test_that("prepared across-donor results use donor measurements in the shared viewer", {
  for (axis in c("condition", "celltype")) for (independent in c(FALSE, TRUE)) {
    if (axis == "celltype" && independent) next
    data <- comparison_view_fixture(independent)
    entry <- comparison_view_entry(data, axis, independent)
    data$contract <- hybrid_result_contract(data, entry)
    original <- data$contract
    current <- shiny::reactiveVal(data)
    shiny::testServer(precomputed_server, args = list(data = current), {
      session$setInputs(readout = "abundance", donor = PRECOMPUTED_ALL, celltype = entry$request$celltype,
        comparison = paste(entry$request$target, "vs", entry$request$reference), results_view = "Volcano")
      expect_true(scope()$across)
      expect_equal(selected_results()$effect_size, entry$result$results$effect_size)
      expect_equal(nrow(selected_samples()), 24L)
      select_feature('["A","A"]')
      plot <- detail_ggplot()
      expect_equal(sum(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomLine"), logical(1))),
        if (independent) 0L else 1L)
      expect_equal(sum(is.finite(plot$data$value)), if (independent) 4L else 8L)
      expect_equal(nrow(report_volcano_ggplot()$data), 3L)
      expect_match(output$volcano_notice$html, "Across-donor comparison", fixed = TRUE)
      expect_match(output$detail_notice$html, "donor summaries", fixed = TRUE)
      expect_match(output$methods_ui$html, "across-donor results", fixed = TRUE)
      expect_equal(nrow(read.csv(output$download_results_csv)), 3L)
      expect_equal(nrow(read.csv(output$download_samples_csv)), 24L)
      expect_s3_class(plotly::plotly_build(plotly::ggplotly(plot, tooltip = "text")), "plotly")
      expect_identical(current()$contract, original)
    })
  }
})

test_that("paired donor detail explains when only unpaired measurements remain", {
  data <- comparison_view_fixture()
  entry <- comparison_view_entry(data)
  rows <- hybrid_result_contract(data, entry)$donor_values$abundance
  rows <- rows[rows$feature_id == "A", ]
  rows$eligible <- (rows$donor %in% c("D1", "D2") & rows$group == "TRT") |
    (rows$donor %in% c("D3", "D4") & rows$group == "NC")
  rows$value[!rows$eligible] <- NA_real_
  plot <- plot_precomputed_donor_values(rows, paired = TRUE)
  expect_s3_class(plot, "ggplot")
  expect_match(ggplot2::ggplot_build(plot)$data[[1L]]$label, "No complete eligible donor pairs", fixed = TRUE)
})

test_that("new spatial summaries retain common protein order for paired and independent donor views", {
  for (axis in c("condition", "celltype")) for (independent in c(FALSE, TRUE)) {
    if (axis == "celltype" && independent) next
    data <- comparison_view_fixture(independent)
    entry <- comparison_view_entry(data, axis, independent, "colocalization")
    data$contract <- hybrid_result_contract(data, entry)
    current <- shiny::reactiveVal(data)
    shiny::testServer(precomputed_server, args = list(data = current), {
      session$setInputs(readout = "colocalization", donor = PRECOMPUTED_ALL, celltype = entry$request$celltype,
        comparison = paste(entry$request$target, "vs", entry$request$reference), results_view = "Spatial patterns",
        spatial_kind = "differences", spatial_layout = "compare")
      session$setInputs(spatial_set = "manual", spatial_markers = c("A", "B", "C"),
        spatial_order = "custom", spatial_custom_order = "C,B,A", spatial_apply = 1L)
      plots <- spatial_comparison_plots()
      for (plot in plots) {
        expect_identical(plot$scales$get_scales("x")$limits, c("C", "B", "A"))
        expect_identical(plot$scales$get_scales("y")$limits, c("A", "B", "C"))
        expect_equal(sum(plot$data$status == "supported"), 6L)
        expect_true(all(is.na(plot$data$value[plot$data$marker_1 == plot$data$marker_2])))
        expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
      }
      expect_identical(plots$reference$scales$get_scales("fill")$limits,
        plots$target$scales$get_scales("fill")$limits)
      expect_identical(plots$target$labels$fill, "Mean donor median log2 ratio")
      expect_identical(plots$difference$labels$fill, "Mean donor difference")
      expect_s3_class(plotly::plotly_build(report_spatial_plotly(plots$difference)), "plotly")
      exported <- read.csv(output$download_spatial_compare_csv)
      expect_equal(nrow(exported), 27L)
      expect_identical(unique(exported$marker_1), c("C", "B", "A"))
      expect_true(all(is.na(exported$value[exported$marker_1 == exported$marker_2])))
      network <- network_profiles()
      expect_setequal(as.character(network$condition), c(entry$request$target, entry$request$reference))
      expect_true(all(network$n_donors == if (independent) 2L else 4L))
    })
  }
})

test_that("comparison labels cannot silently mix analyses with different donor sets", {
  data <- comparison_view_fixture()
  first <- comparison_view_entry(data, id = "first")
  second <- first
  second$id <- "second"
  a <- hybrid_result_contract(data, first)
  b <- hybrid_result_contract(data, second)
  a$differential_results$abundance <- rbind(a$differential_results$abundance, b$differential_results$abundance)
  data$contract <- a
  selected <- precomputed_explorer_scope(data, "abundance", PRECOMPUTED_ALL, "B", "TRT vs NC [second]")
  expect_identical(selected$analysis_id, "second")
  expect_true(all(selected$rows$analysis_id == "second"))
  expect_length(selected$comparisons, 2L)
})

test_that("spatial panels cannot mix same-label analyses from another readout", {
  data <- comparison_view_fixture()
  self <- comparison_view_entry(data, readout = "self_proximity", id = "self")
  first <- comparison_view_entry(data, readout = "colocalization", id = "co-1")
  second <- first
  second$id <- "co-2"
  second$result$donor_values$value <- second$result$donor_values$value + 100
  second$result$donor_contrasts$value_a <- second$result$donor_contrasts$value_a + 100
  second$result$donor_contrasts$value_b <- second$result$donor_contrasts$value_b + 100
  a <- hybrid_result_contract(data, self)
  b <- hybrid_result_contract(data, first)
  c <- hybrid_result_contract(data, second)
  for (section in c("donor_contrasts", "differential_results", "donor_values"))
    a[[section]]$colocalization <- rbind(b[[section]]$colocalization, c[[section]]$colocalization)
  data$contract <- a
  selected <- precomputed_explorer_scope(data, "self_proximity", PRECOMPUTED_ALL, "B")
  scoped <- precomputed_scope_contract(data$contract, selected)
  expect_equal(nrow(scoped$donor_values$colocalization), 0L)
  alignment <- report_comparison_alignment(scoped, PRECOMPUTED_ALL, "B", selected$comparison)
  expect_true(all(alignment$profiles$marker_1 == alignment$profiles$marker_2))
  expect_equal(alignment$profiles$value[alignment$profiles$marker_1 == "A" & alignment$profiles$panel == "TRT"], 2.85)
})

test_that("comparison choices follow the selected cell type", {
  data <- comparison_view_fixture()
  first <- comparison_view_entry(data)
  second <- comparison_view_entry(data, axis = "celltype", id = "celltype")
  a <- hybrid_result_contract(data, first)
  b <- hybrid_result_contract(data, second)
  a$differential_results$abundance <- rbind(a$differential_results$abundance, b$differential_results$abundance)
  data$contract <- a
  selected <- precomputed_explorer_scope(data, "abundance", PRECOMPUTED_ALL, "T")
  expect_identical(selected$comparison, "T vs B")
  expect_true(all(selected$rows$celltype == "T"))
  expect_equal(nrow(selected$rows), 3L)
})

test_that("within-donor cell-type exploration uses the fixed condition in every view", {
  data <- comparison_view_fixture()
  profiles <- precomputed_app_data(precomputed_fixture())$contract$sample_summaries
  for (readout in c("self_proximity", "colocalization")) {
    profiles[[readout]]$metric <- "population_mean_log2_ratio_missing_zero"
    profiles[[readout]]$value <- ifelse(profiles[[readout]]$celltype == "T", 2, 1)
  }
  data$contract$sample_summaries <- profiles
  for (readout in c("abundance", "colocalization")) {
    request <- validate_hybrid_request(data, list(analysis_scope = "within_donor", comparison_axis = "celltype",
      design = "within_donor", donor = "D1", condition = "NC", celltype = "T", target = "T", reference = "B",
      readout = readout, markers = c("A", "B", "C"), min_cells_per_group = 3L))
    rows <- hybrid_feature_table(request)
    rows$effect_size <- .5
    rows$p_value <- .01
    rows$p_adj <- .03
    rows$n_cells_a <- rows$n_cells_b <- 4L
    entry <- list(id = "within-celltype", label = "Within donor cell types", request = request,
      result = list(results = rows, provenance = list(dataset_id = request$dataset_id,
        bundle_md5 = request$bundle_md5, request = request)))
    prepared <- data
    prepared$contract <- hybrid_result_contract(data, entry)
    shiny::testServer(precomputed_server, args = list(data = shiny::reactiveVal(prepared)), {
      session$setInputs(readout = readout, donor = "D1", celltype = "T", comparison = "T vs B", results_view = "Volcano")
      expect_false(scope()$across)
      expect_identical(scope()$fixed_condition, "NC")
      expect_match(output$volcano_notice$html, "Exploratory component-level test", fixed = TRUE)
      expect_equal(nrow(report_volcano_ggplot()$data), 3L)
      if (readout == "abundance") {
        select_feature('["A","A"]')
        expect_equal(nrow(detail_distribution()), 8L)
        expect_setequal(as.character(detail_distribution()$condition), c("B", "T"))
        expect_identical(detail_ggplot()$labels$x, "Cell type")
      } else {
        session$setInputs(results_view = "Spatial patterns", spatial_layout = "compare")
        expect_setequal(spatial_scope()$conditions, c("B", "T"))
        plots <- spatial_comparison_plots()
        expect_true(all(plots$reference$data$value[plots$reference$data$status == "supported"] == 1))
        expect_true(all(plots$target$data$value[plots$target$data$status == "supported"] == 2))
        expect_true(all(plots$difference$data$value[plots$difference$data$status == "supported"] == .5))
      }
      session$setInputs(explore_3d = 1L)
      expect_identical(graph_selection()$conditions, "NC")
    })
  }
})
