network_fixture <- function() {
  data <- precomputed_app_data(precomputed_fixture())
  rows <- data$contract$sample_summaries$colocalization
  rows$metric <- "population_mean_log2_ratio_missing_zero"
  rows$value <- ifelse(rows$marker_2 == "B", ifelse(rows$condition == "TRT", 1, 0.5),
    ifelse(rows$condition == "TRT", 0.1, -0.4))
  rows$n_scored_cells <- ifelse(rows$condition == "NC" & rows$marker_2 == "B", 1, 2)
  data$contract$sample_summaries$colocalization <- rows
  data
}

test_that("network profiles use only the chosen scope and enforce mean-score semantics", {
  data <- network_fixture()
  rows <- report_network_profiles(data$contract, "D1", "B", "TRT vs NC")
  expect_equal(nrow(rows), 4L)
  expect_identical(unique(rows$donor), "D1")
  expect_identical(unique(rows$celltype), "B")
  expect_identical(levels(rows$condition), c("NC", "TRT"))
  expect_equal(rows$pct_detected, rows$n_scored_cells / rows$n_total_cells)
  expect_error(report_network_profiles(data$contract, "D1", PRECOMPUTED_ALL, "TRT vs NC"), "one cell type")
  expect_error(report_network_profiles(data$contract, "D1", "B", "unknown"), "comparison")
  bad <- data$contract
  bad$sample_summaries$colocalization$metric <- "median_z"
  expect_error(report_network_profiles(bad, "D1", "B", "TRT vs NC"), "population mean")
  bad <- data$contract
  rows <- bad$sample_summaries$colocalization
  reversed <- rows
  reversed$marker_1 <- rows$marker_2
  reversed$marker_2 <- rows$marker_1
  bad$sample_summaries$colocalization <- rbind(rows, reversed)
  expect_error(report_network_profiles(bad, "D1", "B", "TRT vs NC"), "Multiple summaries")
})

test_that("network and matrix profiles share cell weighting and incomplete-sample protection", {
  contract <- network_fixture()$contract
  rows <- contract$sample_summaries$colocalization
  at <- rows$donor == "D1" & rows$celltype == "B" & rows$condition == "NC"
  rows$n_total_cells[at] <- rows$n_scored_cells[at] <- 10
  rows$value[at] <- 1
  zero <- rows[at, ]
  zero$sample_id <- "D1_NC_zero"
  zero$n_total_cells <- 30
  zero$n_scored_cells <- zero$value <- 0
  contract$sample_summaries$colocalization <- rbind(rows, zero)
  design <- contract$sample_summaries$abundance
  extra <- design[design$donor == "D1" & design$celltype == "B" & design$condition == "NC", ]
  extra$sample_id <- "D1_NC_zero"
  contract$sample_summaries$abundance <- rbind(design, extra)

  network <- report_network_profiles(contract, "D1", "B", "TRT vs NC")
  pooled <- network[network$condition == "NC", ]
  expect_equal(pooled$value, rep(.25, 2))
  expect_equal(pooled$pct_detected, rep(.25, 2))
  expect_equal(pooled$n_total_cells, rep(40, 2))
  expect_equal(pooled$n_scored_cells, rep(10, 2))
  matrix <- report_spatial_alignment(contract)$profiles
  expect_equal(matrix[donor == "D1" & celltype == "B" & panel == "NC", value], pooled$value)

  contract$sample_summaries$colocalization <- rbind(rows, zero[zero$marker_2 != "B", ])
  incomplete <- report_network_profiles(contract, "D1", "B", "TRT vs NC")
  absent <- incomplete[incomplete$condition == "NC" & incomplete$marker_2 == "B", ]
  expect_identical(absent$status, "unavailable")
  expect_true(all(is.na(absent[c("value", "pct_detected", "n_total_cells", "n_scored_cells")])))
  expect_equal(incomplete$value[incomplete$condition == "NC" & incomplete$marker_2 == "C"], .25)
  plot <- plot_report_network(incomplete, c("A", "B", "C"), min_detection = 0, min_mean = 0)
  expect_false(any(plot$data$condition == "NC" & plot$data$marker_2 == "B"))
  matrix <- report_spatial_alignment(contract)$profiles
  expect_true(is.na(matrix[donor == "D1" & celltype == "B" & panel == "NC" & marker_2 == "B", value]))
})

test_that("profile-only networks display scoped conditions without a differential test", {
  data <- network_fixture()
  rows <- data$contract$sample_summaries$colocalization
  extra <- rows[rows$donor == "D1" & rows$celltype == "B" & rows$condition == "TRT", ]
  extra$condition <- "Third"
  extra$sample_id <- "D1_Third"
  data$contract$sample_summaries$colocalization <- rbind(rows, extra)
  selected <- report_network_profiles(data$contract, "D1", "B", "TRT vs NC")
  expect_identical(levels(selected$condition), c("NC", "TRT"))
  data$contract$donor_contrasts$colocalization <- data$contract$donor_contrasts$colocalization[FALSE, ]
  original <- data$contract
  profiles <- report_network_profiles(data$contract, "D1", "B")
  expect_identical(levels(profiles$condition), c("NC", "TRT", "Third"))
  expect_equal(nrow(profiles), 6L)
  expect_identical(report_network_profiles(data$contract, "D1", "B", NA_character_), profiles)
  expect_error(report_network_profiles(data$contract, "D1", "B", "unknown"), "comparison")
  plot <- plot_report_network(profiles, c("A", "B", "C"))
  expect_setequal(as.character(plot$layers[[3L]]$data$condition), c("NC", "TRT", "Third"))
  expect_warning(plotly::plotly_build(report_network_plotly(plot)), NA)
  expect_identical(data$contract, original)

  current <- shiny::reactiveVal(data)
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(readout = "colocalization", donor = "D1", celltype = "B",
      results_view = "Spatial patterns", spatial_view = "network", network_markers = c("A", "B", "C"),
      network_detection = 0.5, network_mean = 0.1)
    expect_equal(nrow(selected_results()), 0L)
    expect_equal(nrow(network_profiles()), 6L)
    expect_true(nrow(network_ggplot()$data) > 0L)
    expect_equal(read.csv(output$download_network_csv)$value, network_ggplot()$data$value)
  })
})

test_that("networks filter each condition and retain shared positions and scales", {
  rows <- report_network_profiles(network_fixture()$contract, "D1", "B", "TRT vs NC")
  plot <- plot_report_network(rows, c("C", "A", "B"))
  edges <- plot$data
  expect_setequal(paste(edges$condition, edges$marker_2), c("TRT B", "NC C"))
  expect_setequal(edges$value, c(1, -0.4))
  nodes <- plot$layers[[3L]]$data
  expect_equal(unname(as.matrix(nodes[nodes$condition == "NC", c("x", "y")])),
    unname(as.matrix(nodes[nodes$condition == "TRT", c("x", "y")])))
  empty <- plot_report_network(rows, c("A", "B", "C"), min_mean = 2)
  expect_equal(nrow(empty$data), 0L)
  expect_equal(empty$layers[[3L]]$data, nodes)
  expect_equal(empty$scales$get_scales("colour")$limits, plot$scales$get_scales("colour")$limits)
  expect_setequal(as.character(empty$layers[[4L]]$data$condition), c("NC", "TRT"))
  two <- plot_report_network(rows, c("A", "B"))
  expect_true(all(is.finite(two$layers[[3L]]$data$x)))
  expect_true(all(is.finite(two$layers[[3L]]$data$y)))
  expect_error(plot_report_network(rows, "A"), "2 to 20")
  expect_error(plot_report_network(rows, c("A", "unknown")), "available")
  expect_error(plot_report_network(rows, c("A", "B"), min_detection = 1.1), "between 0 and 1")
  expect_error(plot_report_network(rows, c("A", "B"), min_mean = NA_real_), "finite")
  expect_warning(built <- plotly::plotly_build(report_network_plotly(plot, source = "test-network")), NA)
  keys <- unlist(lapply(built$x$data, function(trace) trace$key), use.names = FALSE)
  expect_true(all(report_feature_key(edges$marker_1, edges$marker_2) %in% keys))
  expect_true("plotly_click" %in% built$x$shinyEvents)
  hover <- unlist(lapply(built$x$data, function(trace) trace$text), use.names = FALSE)
  expect_true(any(grepl("Detected cells", hover, fixed = TRUE)))
  expect_warning(plotly::plotly_build(report_network_plotly(empty)), NA)
})

test_that("network clicks link to other views and exports use the displayed edges", {
  current <- shiny::reactiveVal(network_fixture())
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(readout = "colocalization", donor = "D1", celltype = "B", comparison = "TRT vs NC",
      results_view = "Spatial patterns", spatial_view = "network", network_markers = c("A", "B", "C"),
      network_detection = 0.5, network_mean = 0.1)
    expect_equal(nrow(network_ggplot()$data), 2L)
    csv <- utils::read.csv(output$download_network_csv)
    expect_equal(csv$value, network_ggplot()$data$value)
    expect_true(all(csv$minimum_detection_fraction == 0.5))
    expect_gt(file.info(output$download_network_svg)$size, 0)
    event_id <- paste0("plotly_click-", session$ns("network"))
    session$userData$plotlyShinyEventIDs <- event_id
    event <- jsonlite::toJSON(data.frame(key = report_feature_key("C", "A")), auto_unbox = TRUE)
    do.call(session$rootScope()$setInputs, stats::setNames(list(event), event_id))
    expect_identical(selected_key(), report_feature_key("A", "C"))
    session$setInputs(results_view = "Volcano")
    expect_identical(selected_key(), report_feature_key("A", "C"))
    session$setInputs(network_mean = 2)
    expect_equal(nrow(utils::read.csv(output$download_network_csv)), 0L)
    session$setInputs(readout = "self_proximity")
    expect_false(grepl('value="network"', as.character(output$spatial_view_control$html), fixed = TRUE))
    expect_error(network_profiles(), "Colocalization")
  })
})
