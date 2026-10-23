test_that("automatic proteins combine donors, comparisons, readouts and explicit session updates", {
  data <- precomputed_app_data(precomputed_fixture())
  contract <- data$contract
  for (section in c("donor_contrasts", "differential_results"))
    contract[[section]] <- lapply(contract[[section]], function(x) x[FALSE, ])
  self <- data$contract$donor_contrasts$self_proximity[1, ]
  self$marker_1 <- self$marker_2 <- "A"
  self$celltype <- "B"
  self$p_adj <- .8
  co <- self
  co$marker_1 <- "C"
  co$marker_2 <- "B"
  co$donor <- "D2"
  co$comparison <- "Third vs NC"
  contract$donor_contrasts$self_proximity <- self
  contract$donor_contrasts$colocalization <- co
  alignment <- report_spatial_alignment(contract)
  result <- report_spatial_protein_set(contract, "B", alignment)
  expect_identical(result$markers, c("A", "B", "C"))
  expect_length(report_spatial_protein_set(contract, "T", alignment)$markers, 0L)
  history <- list(new = list(request = list(celltype = "B", readout = "colocalization"),
    result = list(results = data.frame(marker_1 = c("D", "Missing"), marker_2 = "A", effect_size = c(.2, NA)))))
  expect_identical(report_spatial_protein_set(contract, "B", alignment, history)$markers, c("A", "B", "C", "D"))
  history$new$request$readout <- "abundance"
  expect_identical(report_spatial_protein_set(contract, "B", alignment, history)$markers, result$markers)
  original <- contract
  contract$donor_contrasts <- lapply(contract$donor_contrasts, function(x) x[FALSE, ])
  profiles <- data.frame(donor = "D1", celltype = "B", panel = "NC", marker_1 = c("A", "X"),
    marker_2 = c("B", "Y"), value = c(0, 1), status = "supported", pct_detected = c(.2, 0))
  alignment$profiles <- profiles
  alignment$differences <- alignment$differences[FALSE, ]
  fallback <- report_spatial_protein_set(contract, "B", alignment)
  expect_identical(fallback$markers, c("A", "B"))
  expect_match(fallback$source, "profiles")
  expect_equal(original$donor_contrasts$self_proximity$p_adj, .8)
})

test_that("reference ordering groups shared profiles deterministically without imputing missing scores", {
  markers <- c("A", "B", "C", "D")
  values <- matrix(c(1, .9, .1, .1, .9, 1, .1, .1, .1, .1, 1, .9, .1, .1, .9, 1), 4)
  profiles <- expand.grid(marker_1 = markers, marker_2 = markers, stringsAsFactors = FALSE)
  profiles$value <- as.vector(values)
  profiles$donor <- "D1"
  profiles$panel <- "NC"
  profiles$celltype <- "B"
  profiles$status <- "supported"
  profiles$pct_detected <- 1
  profiles$p_adj <- NA_real_
  alignment <- list(cfg = list(), markers = data.frame(), profiles = profiles)
  ordered <- report_spatial_order(alignment, "B", c("E", rev(markers)), "D1", "NC")
  expect_match(ordered$method, "Ward.D2")
  expect_equal(abs(diff(match(c("A", "B"), ordered$markers))), 1)
  expect_equal(abs(diff(match(c("C", "D"), ordered$markers))), 1)
  expect_identical(tail(ordered$markers, 1), "E")
  alignment$profiles <- profiles[rev(seq_len(nrow(profiles))), ]
  expect_identical(report_spatial_order(alignment, "B", c(markers, "E"), "D1", "NC"), ordered)
  alignment$profiles <- profiles[profiles$marker_1 == profiles$marker_2, ]
  sparse <- alignment
  fallback <- report_spatial_order(alignment, "B", rev(markers), "D1", "NC")
  expect_identical(fallback$markers, markers)
  expect_match(fallback$method, "insufficient")
  expect_identical(alignment, sparse)
  alignment$markers <- data.frame(celltype = "B", marker = c("C", "A", "D", "B"), marker_rank = c(3, 2, 4, 1))
  alignment$cfg <- list(ordering_reference_donor = "D1", ordering_reference_condition = "NC")
  expect_identical(report_spatial_order(alignment, "B", markers, "D1", "NC")$markers, c("B", "A", "C", "D"))
  expect_identical(report_spatial_order(alignment, "B", markers, "D2", "NC")$markers, markers)
})

test_that("automatic and manual axes are independent and new comparisons need an explicit refresh", {
  fixture <- precomputed_app_data(precomputed_fixture())
  for (section in c("donor_contrasts", "differential_results")) for (readout in c("self_proximity", "colocalization")) {
    rows <- fixture$contract[[section]][[readout]]
    fixture$contract[[section]][[readout]] <- rows[rows$marker_1 != "C" & rows$marker_2 != "C", ]
  }
  pending_history <- shiny::reactiveVal(list())
  server <- precomputed_server
  environment(server) <- list2env(list(hybrid_analysis_server = function(...) list(history = pending_history)),
    parent = environment(precomputed_server))
  current <- shiny::reactiveVal(fixture)
  shiny::testServer(server, args = list(data = current), {
    session$setInputs(readout = "colocalization", donor = "D1", celltype = "B")
    expect_setequal(spatial_markers(), c("A", "B"))
    before <- spatial_markers()
    session$setInputs(donor = "D2")
    expect_identical(spatial_markers(), before)
    pending_history(list(new = list(label = "New comparison", request = list(celltype = "B", readout = "colocalization"),
      result = list(results = data.frame(marker_1 = "A", marker_2 = "C", effect_size = .5)))))
    session$flushReact()
    expect_identical(spatial_markers(), before)
    expect_match(output$spatial_auto_status, "New session results")
    session$setInputs(spatial_refresh = 1L)
    expect_setequal(spatial_markers(), c("A", "B", "C"))
    session$setInputs(spatial_set = "manual", spatial_markers = c("C", "A"), spatial_order = "alphabetical", spatial_apply = 1L)
    expect_identical(spatial_markers(), c("A", "C"))
    session$setInputs(spatial_order = "custom", spatial_custom_order = "C,A", spatial_apply = 2L)
    expect_identical(spatial_markers(), c("C", "A"))
    session$setInputs(spatial_refresh = 2L, donor = "D1", celltype = "T")
    expect_identical(spatial_markers(), c("C", "A"))
    session$setInputs(spatial_custom_order = "C", spatial_apply = 3L)
    expect_match(output$spatial_selection_status, "Not applied")
    expect_identical(spatial_markers(), c("C", "A"))
    session$setInputs(celltype = "B", spatial_reset = 1L)
    expect_identical(spatial_applied_mode(), "automatic")
    expect_setequal(spatial_markers(), c("A", "B", "C"))
    expect_identical(current()$contract, fixture$contract)
  })
})
