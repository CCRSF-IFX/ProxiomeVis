hybrid_module_fixture <- function() {
  data <- hybrid_fixture()
  data$bundle_path <- "/trusted/analysis.rds"
  data$bundle_md5 <- paste(rep("a", 32), collapse = "")
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
  list(results = data.frame(feature_id = request$markers, effect_size = 0.2, p_value = 0.01,
                            p_adj = 0.02, n_cells_a = length(request$group_a_cells),
                            n_cells_b = length(request$group_b_cells)),
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
                                             code_dir = tempdir(), active = shiny::reactive(FALSE)), {
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
                                             code_dir = tempdir(), active = shiny::reactive(FALSE)), {
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

test_that("invalid selected cell IDs are rejected rather than silently discarded", {
  current <- shiny::reactiveVal(hybrid_module_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
                                             code_dir = tempdir(), active = shiny::reactive(FALSE)), {
    session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
                      readout = "abundance", acknowledge = TRUE)
    valid <- current()$metadata$component[c(1:3, 43:45)]
    selected_cells(c(valid, "unknown-cell"))
    session$setInputs(run = 1L)
    expect_length(probe$state$jobs, 0L)
    expect_match(output$status, "selected|selection|scope|unique", ignore.case = TRUE)
    selected_cells(c(valid, valid[1L]))
    session$setInputs(run = 2L)
    expect_length(probe$state$jobs, 0L)
    expect_match(output$status, "selected|selection|scope|unique", ignore.case = TRUE)
  })
})

test_that("all readouts run without marker or cell dropdown inputs", {
  for (readout in c("abundance", "self_proximity", "colocalization")) {
    current <- shiny::reactiveVal(hybrid_module_fixture())
    probe <- hybrid_module_probe()
    shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
                                               code_dir = tempdir(), active = shiny::reactive(FALSE)), {
      session$setInputs(donor = "D1", celltype = "B", target = "TRT", reference = "NC",
                        readout = readout, acknowledge = TRUE)
      session$setInputs(min_cells_per_group = 35L)
      selected_cells(current()$metadata$component[c(1:40, 43:82)])
      session$setInputs(run = 1L)
      expect_length(probe$state$jobs, 1L)
      submitted <- probe$state$jobs[[1L]]$request
      expect_identical(submitted$markers, c("A", "B", "C"))
      expect_identical(submitted$min_cells_per_group, if (readout == "abundance") 3L else 35L)
      expect_length(submitted$group_a_cells, 40L)
      expect_length(submitted$group_b_cells, 40L)
      session$setInputs(clear_selection = 1L)
      expect_length(selected_metadata()$component, 84L)
      selected_cells(current()$metadata$component[1:3])
      session$setInputs(donor = "D2")
      expect_length(selected_cells(), 0L)
      probe$state$jobs[[1L]]$resolve(hybrid_module_result(submitted))
      settle_module_task(session, task)
      if (readout != "abundance") expect_match(output$result_scope, "Minimum observed cells per condition: 35")
    })
  }
})

test_that("worker failures are actionable and cancellation covers buttons and closed sessions", {
  current <- shiny::reactiveVal(hybrid_module_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
                                             code_dir = tempdir(), active = shiny::reactive(FALSE)), {
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
