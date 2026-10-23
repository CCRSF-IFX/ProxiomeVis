hybrid_module_fixture <- function() {
  data <- hybrid_fixture()
  data$bundle_path <- "/trusted/analysis.rds"
  data$bundle_md5 <- paste(rep("a", 32), collapse = "")
  data$assets$proximity <- list(path = "proximity")
  data$embeddings <- lapply(data$embeddings, function(coordinates) {
    matrix(seq_len(2L * nrow(data$metadata)), ncol = 2L,
           dimnames = list(data$metadata$component, colnames(coordinates)))
  })
  data
}

hybrid_module_probe <- function() {
  state <- new.env(parent = emptyenv())
  state$jobs <- list()
  module <- hybrid_analysis_server
  environment(module) <- list2env(list(start_hybrid_job = function(bundle_path, request, code_dir, manager) {
    pending <- new.env(parent = emptyenv())
    pending$request <- request
    pending$cancelled <- FALSE
    promise <- promises::promise(function(resolve, reject) {
      pending$resolve <- resolve
      pending$reject <- reject
    })
    state$jobs[[length(state$jobs) + 1L]] <- pending
    list(promise = promise, cancel = function() pending$cancelled <- TRUE)
  }), parent = environment(hybrid_analysis_server))
  list(module = module, state = state)
}

hybrid_module_result <- function(request) {
  rows <- hybrid_feature_table(request)
  rows$effect_size <- 0.2
  rows$p_value <- 0.01
  rows$p_adj <- 0.02
  rows$n_cells_a <- length(request$group_a_cells)
  rows$n_cells_b <- length(request$group_b_cells)
  list(results = rows,
       provenance = list(dataset_id = request$dataset_id, bundle_md5 = request$bundle_md5, request = request,
                         filters = list(donor_isotype_fraction_cutoff = 0.001)))
}

settle_module_task <- function(session, task) {
  deadline <- Sys.time() + 5
  while (identical(task$status(), "running") && Sys.time() < deadline) {
    later::run_now(0.01)
    session$flushReact()
  }
  expect_false(identical(task$status(), "running"))
  for (i in seq_len(3L)) {
    later::run_now(0)
    session$flushReact()
  }
}

test_that("completed analyses retain their submitted scope when controls change", {
  current <- shiny::reactiveVal(hybrid_module_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
                                             code_dir = tempdir()), {
    session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
                      readout = "abundance", acknowledge = TRUE)
    session$setInputs(run = 1L)
    expect_length(probe$state$jobs, 1L)
    submitted <- probe$state$jobs[[1L]]$request
    expect_identical(submitted$bundle_md5, current()$bundle_md5)
    session$setInputs(target = "NC", reference = "TRT")
    probe$state$jobs[[1L]]$resolve(hybrid_module_result(submitted))
    settle_module_task(session, task)
    expect_identical(result()$provenance$request, submitted)
    expect_identical(last_request(), submitted)
    expect_match(output$result_scope, "TRT minus NC", fixed = TRUE)
    expect_identical(result()$results$feature_id, c("A", "B", "C"))
    session$setInputs(acknowledge = FALSE, run = 2L)
    expect_length(probe$state$jobs, 1L)
    expect_match(output$status, "Confirm the exploratory")
    expect_identical(last_request(), submitted)
    expect_identical(result()$provenance$request, submitted)
  })
})

test_that("dataset replacement cancels the worker and hides a late result even for the same dataset ID", {
  current <- shiny::reactiveVal(hybrid_module_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
                                             code_dir = tempdir()), {
    session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
                      readout = "abundance", acknowledge = TRUE)
    session$setInputs(run = 1L)
    submitted <- probe$state$jobs[[1L]]$request
    replacement <- current()
    replacement$bundle_md5 <- paste(rep("b", 32), collapse = "")
    current(replacement)
    session$flushReact()
    expect_true(probe$state$jobs[[1L]]$cancelled)
    expect_null(last_request())
    probe$state$jobs[[1L]]$resolve(hybrid_module_result(submitted))
    settle_module_task(session, task)
    expect_error(result(), class = "shiny.silent.error")
  })
})



test_that("all readouts run without marker or cell dropdown inputs", {
  for (readout in c("abundance", "self_proximity", "colocalization")) {
    current <- shiny::reactiveVal(hybrid_module_fixture())
    probe <- hybrid_module_probe()
    shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
                                               code_dir = tempdir()), {
      session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
                        readout = readout, acknowledge = TRUE)
      session$setInputs(min_cells_per_group = 35L)
      session$setInputs(cells = "obsolete-cell-selection")
      session$setInputs(run = 1L)
      expect_length(probe$state$jobs, 1L)
      submitted <- probe$state$jobs[[1L]]$request
      expect_identical(submitted$markers, c("A", "B", "C"))
      expect_identical(submitted$min_cells_per_group, if (readout == "abundance") 3L else 35L)
      expect_length(submitted$group_a_cells, 42L)
      expect_length(submitted$group_b_cells, 42L)
      expect_match(output$selection_summary, "42 target cells; 42 reference cells", fixed = TRUE)
      probe$state$jobs[[1L]]$resolve(hybrid_module_result(submitted))
      settle_module_task(session, task)
      if (readout != "abundance") expect_match(output$result_scope, "Minimum observed cells per condition: 35")
    })
  }
})

test_that("late worker completion can settle after its module session is destroyed", {
  for (outcome in c("success", "error")) {
    current <- shiny::reactiveVal(hybrid_module_fixture())
    probe <- hybrid_module_probe()
    saved <- new.env(parent = emptyenv())
    shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
                                               code_dir = tempdir()), {
      session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
                        readout = "abundance", acknowledge = TRUE, run = 1L)
      saved$task <- task
      expect_identical(task$status(), "running")
    })
    expect_true(probe$state$jobs[[1L]]$cancelled)
    if (outcome == "success") {
      probe$state$jobs[[1L]]$resolve(hybrid_module_result(probe$state$jobs[[1L]]$request))
      expect_warning(for (i in seq_len(10L)) later::run_now(0.01), NA)
    } else {
      probe$state$jobs[[1L]]$reject(simpleError("Analysis cancelled."))
      expect_warning(for (i in seq_len(10L)) later::run_now(0.01), "Analysis cancelled")
    }
    expect_identical(shiny::isolate(saved$task$status()), outcome)
  }
})

test_that("worker failures are actionable and cancellation covers buttons and closed sessions", {
  current <- shiny::reactiveVal(hybrid_module_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
                                             code_dir = tempdir()), {
    session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
                      readout = "abundance", acknowledge = TRUE)
    session$setInputs(run = 1L)
    probe$state$jobs[[1L]]$reject(simpleError("Install pixelatorR in the analysis worker to run differential tests."))
    expect_warning(settle_module_task(session, task), "Install pixelatorR")
    expect_match(output$status, "Install pixelatorR")
    expect_error(result(), class = "shiny.silent.error")
    session$setInputs(run = 2L)
    expect_length(probe$state$jobs, 2L)
    session$setInputs(cancel = 1L)
    expect_true(probe$state$jobs[[2L]]$cancelled)
    probe$state$jobs[[2L]]$reject(simpleError("Analysis cancelled."))
    expect_warning(settle_module_task(session, task), "Analysis cancelled")
    expect_match(output$status, "Analysis cancelled")
    expect_error(result(), class = "shiny.silent.error")
    session$setInputs(run = 3L)
    expect_length(probe$state$jobs, 3L)
    session$close()
    expect_true(probe$state$jobs[[3L]]$cancelled)
    probe$state$jobs[[3L]]$reject(simpleError("Analysis cancelled."))
    expect_warning(for (i in seq_len(3L)) later::run_now(0), "Analysis cancelled")
  })
})

test_that("completed comparisons stay in history and open only on an explicit request", {
  current <- shiny::reactiveVal(hybrid_module_fixture())
  probe <- hybrid_module_probe()
  opened <- character()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
    code_dir = tempdir(), on_open = function(id) opened <<- c(opened, id)), {
    session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
      readout = "abundance", acknowledge = TRUE, run = 1L)
    first <- probe$state$jobs[[1L]]
    first$resolve(hybrid_module_result(first$request))
    settle_module_task(session, task)
    first_id <- selected_id()
    expect_length(history(), 1L)
    expect_length(opened, 0L)
    expect_null(open_result())
    session$setInputs(target = "NC", reference = "TRT", run = 2L)
    expect_identical(result()$provenance$request, first$request)
    second <- probe$state$jobs[[2L]]
    second$resolve(hybrid_module_result(second$request))
    settle_module_task(session, task)
    expect_length(history(), 2L)
    expect_identical(selected_entry()$request, second$request)
    session$setInputs(history_id = first_id, open_result = 1L)
    expect_identical(opened, first_id)
    expect_identical(open_result()$request, first$request)
    expect_identical(result()$provenance$request, first$request)
    expect_match(output$result_scope, "TRT minus NC", fixed = TRUE)
    package <- readRDS(output$result_rds)
    expect_identical(package$result_version, 1L)
    expect_identical(package$entry, history()[[first_id]])
    expect_named(package, c("result_version", "entry"))
    exported <- utils::read.csv(output$results_csv)
    expect_equal(exported$effect_size, history()[[first_id]]$result$results$effect_size)
    replacement <- current()
    replacement$bundle_md5 <- paste(rep("b", 32), collapse = "")
    current(replacement)
    session$flushReact()
    expect_length(history(), 0L)
    expect_null(open_result())
    expect_error(result(), class = "shiny.silent.error")
  })
})

test_that("preflight checks missing inputs and invalid scopes without starting a worker", {
  current <- shiny::reactiveVal(hybrid_module_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
    code_dir = tempdir()), {
    session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
      readout = "abundance", acknowledge = TRUE)
    expect_null(prepared_request()$error)
    expect_match(output$readiness$html, "42 target cells", fixed = TRUE)
    replacement <- current()
    replacement$counts <- NULL
    current(replacement)
    session$flushReact()
    expect_match(prepared_request()$error, "Raw counts are not stored")
    session$setInputs(run = 1L)
    expect_length(probe$state$jobs, 0L)
    replacement$assets$proximity <- NULL
    current(replacement)
    session$setInputs(readout = "self_proximity", min_cells_per_group = 20L, run = 2L)
    expect_match(prepared_request()$error, "proximity sidecar is not configured")
    expect_length(probe$state$jobs, 0L)
    session$setInputs(target = "NC")
    expect_match(prepared_request()$error, "two different conditions")
    session$setInputs(target = "TRT", min_cells_per_group = 43L)
    expect_match(prepared_request()$error, "at least 43")
  })
})

test_that("a late successful completion after cancellation does not enter history", {
  current <- shiny::reactiveVal(hybrid_module_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
    code_dir = tempdir()), {
    session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
      readout = "abundance", acknowledge = TRUE, run = 1L)
    session$setInputs(cancel = 1L)
    pending <- probe$state$jobs[[1L]]
    pending$resolve(hybrid_module_result(pending$request))
    settle_module_task(session, task)
    expect_length(history(), 0L)
    expect_match(output$status, "Analysis cancelled")
    expect_error(result(), class = "shiny.silent.error")
  })
})

test_that("session adapters preserve results and profiles without reusing saved differential statistics", {
  data <- hybrid_module_fixture()
  data$contract$report <- list(source_results = list(old = "old test"),
    eligibility = list(old = "old eligibility"), abundance_heatmaps = list(old = "old matrix"),
    alignment = list(profiles = data.frame(value = 1.3), differences = data.frame(value = 99),
      contrasts = "old contrast", limits = list(profiles = 2, differences = 99)))
  original <- data$contract
  request <- validate_hybrid_request(data, hybrid_request("colocalization"))
  value <- hybrid_module_result(request)
  value$results$n_scored_a <- 40L
  value$results$n_scored_b <- 39L
  value$results$status <- "tested"
  entry <- list(id = "session-1", label = "New result", request = request, result = value)
  contract <- hybrid_result_contract(data, entry)
  rows <- contract$donor_contrasts$colocalization
  for (field in names(value$results)) expect_identical(rows[[field]], value$results[[field]])
  expect_identical(rows$group_a, rep(request$target, nrow(rows)))
  expect_identical(rows$group_b, rep(request$reference, nrow(rows)))
  expect_identical(rows$metric, rep("delta_median_log2_ratio", nrow(rows)))
  expect_true(all(is.na(rows$value_a)))
  expect_true(all(is.na(rows$value_b)))
  expect_equal(nrow(contract$donor_contrasts$abundance), 0L)
  expect_true(all(vapply(contract$differential_results, nrow, integer(1)) == 0L))
  expect_identical(contract$sample_summaries, original$sample_summaries)
  expect_identical(contract$report$alignment$profiles, original$report$alignment$profiles)
  expect_null(contract$report$alignment$differences)
  expect_null(contract$report$alignment$contrasts)
  expect_null(contract$report$source_results)
  expect_null(contract$report$eligibility)
  expect_identical(data$contract, original)
  invalid <- entry
  invalid$result$provenance$dataset_id <- "another dataset"
  expect_error(hybrid_result_contract(data, invalid), "provenance does not match")
  invalid <- entry
  invalid$result$results$marker_2 <- "unknown marker"
  expect_error(hybrid_result_contract(data, invalid), "features do not match")
  entry$request$bundle_md5 <- "another bundle"
  expect_error(hybrid_result_contract(data, entry), "does not belong")
})
