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

hybrid_module_donor_fixture <- function(independent = FALSE, incomplete = FALSE) {
  data <- hybrid_module_fixture()
  original <- data$metadata
  meta <- do.call(rbind, lapply(c("D1", "D2", "D3"), function(donor) {
    do.call(rbind, lapply(c("B", "T"), function(celltype) {
      rows <- original
      rows$donor <- donor
      rows$celltype <- celltype
      rows$component <- paste(donor, celltype, original$component, sep = "_")
      rows$sample_id <- paste(donor, rows$condition, sep = "_")
      rows
    }))
  }))
  if (independent) meta <- meta[meta$condition == ifelse(meta$donor == "D2", "NC", "TRT"), ]
  if (incomplete) meta <- meta[!(meta$donor == "D3" & meta$condition == "NC"), ]
  rownames(meta) <- meta$component
  index <- match(sub("^D[123]_[BT]_", "", meta$component), original$component)
  data$metadata <- meta
  data$counts <- data$counts[, index, drop = FALSE]
  colnames(data$counts) <- meta$component
  data$abundance <- log1p(data$counts)
  data$cell_sources <- data.frame(component = meta$component, original_id = meta$component, sample_id = meta$sample_id)
  data$embeddings <- lapply(data$embeddings, function(coordinates) matrix(seq_len(2L * nrow(meta)), ncol = 2L,
    dimnames = list(meta$component, colnames(coordinates))))
  data$contract$sample_summaries <- lapply(data$contract$sample_summaries,
    function(rows) rows[rows$sample_id %in% meta$sample_id, , drop = FALSE])
  data
}

test_that("explicit within-donor exploration supports both comparison axes", {
  current <- shiny::reactiveVal(hybrid_module_donor_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
    code_dir = tempdir()), {
    session$setInputs(analysis_scope = "within_donor", comparison_axis = "condition", donor = "D1",
      celltype = "B", target = "TRT", reference = "NC", readout = "abundance", acknowledge = TRUE, run = 1L)
    expect_length(probe$state$jobs, 1L)
    submitted <- probe$state$jobs[[1L]]$request
    expect_identical(submitted$analysis_scope, "within_donor")
    expect_identical(submitted$design, "within_donor")
    probe$state$jobs[[1L]]$resolve(hybrid_module_result(submitted))
    settle_module_task(session, task)
    session$setInputs(comparison_axis = "celltype", condition = "NC", target = "T", reference = "B", run = 2L)
    expect_length(probe$state$jobs, 2L)
    submitted <- probe$state$jobs[[2L]]$request
    expect_identical(submitted$comparison_axis, "celltype")
    expect_identical(submitted$condition, "NC")
    expect_identical(submitted$celltype, "T")
    expect_setequal(current()$metadata$celltype[match(submitted$group_a_cells, current()$metadata$component)], "T")
    expect_setequal(current()$metadata$celltype[match(submitted$group_b_cells, current()$metadata$component)], "B")
    probe$state$jobs[[2L]]$resolve(hybrid_module_result(submitted))
    settle_module_task(session, task)
    expect_match(output$result_scope, "Condition: NC", fixed = TRUE)
    contract <- hybrid_result_contract(current(), selected_entry())
    expect_identical(contract$schema_version, 2L)
    expect_true(all(contract$donor_contrasts$abundance$comparison_axis == "celltype"))
    expect_true(all(contract$donor_contrasts$abundance$fixed_condition == "NC"))
  })
})

test_that("paired controls show complete donors and never require the exploratory acknowledgement", {
  current <- shiny::reactiveVal(hybrid_module_donor_fixture(incomplete = TRUE))
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
    code_dir = tempdir()), {
    session$setInputs(analysis_scope = "across_donors", comparison_axis = "condition", design = "paired",
      donors = c("D1", "D2", "D3"), celltype = "B", target = "TRT", reference = "NC",
      readout = "abundance", acknowledge = FALSE)
    expect_null(prepared_request()$error)
    expect_match(output$selection_summary, "84 target cells; 84 reference cells retained", fixed = TRUE)
    expect_match(output$selection_summary, "2 complete donor pairs", fixed = TRUE)
    expect_match(output$selection_summary, "D3 (incomplete pair)", fixed = TRUE)
    session$setInputs(run = 1L)
    expect_length(probe$state$jobs, 1L)
    submitted <- probe$state$jobs[[1L]]$request
    expect_identical(submitted$design, "paired")
    expect_setequal(submitted$donors_a, c("D1", "D2"))
    expect_length(submitted$group_a_cells, 84L)
    session$setInputs(cancel = 1L)
    probe$state$jobs[[1L]]$reject(simpleError("Analysis cancelled."))
    expect_warning(settle_module_task(session, task), "Analysis cancelled")
    session$setInputs(donors = character())
    expect_match(prepared_request()$error, "Choose at least one donor")
    session$setInputs(run = 2L)
    expect_length(probe$state$jobs, 1L)
  })
})

test_that("independent controls keep donor groups separate and celltype mode uses a fixed condition", {
  current <- shiny::reactiveVal(hybrid_module_donor_fixture(independent = TRUE))
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
    code_dir = tempdir()), {
    session$setInputs(analysis_scope = "across_donors", comparison_axis = "condition", design = "independent",
      donors = c("D1", "D2", "D3"), celltype = "B", target = "TRT", reference = "NC", readout = "abundance")
    expect_null(prepared_request()$error)
    expect_match(output$selection_summary, "2 target donors; 1 reference donors", fixed = TRUE)
    expect_setequal(prepared_request()$request$donors_a, c("D1", "D3"))
    expect_identical(prepared_request()$request$donors_b, "D2")
    session$setInputs(comparison_axis = "celltype", condition = "TRT", target = "T", reference = "B", design = "paired")
    expect_null(prepared_request()$error)
    expect_match(output$selection_summary, "2 complete donor pairs", fixed = TRUE)
    request <- prepared_request()$request
    expect_true(all(current()$metadata$condition[match(c(request$group_a_cells, request$group_b_cells),
      current()$metadata$component)] == "TRT"))
    expect_identical(request$celltype, "T")
  })
})

test_that("across-donor adapters retain donor measurements and validate as prepared results", {
  for (scenario in c("paired_conditions", "independent_conditions", "paired_celltypes")) {
    data <- hybrid_module_donor_fixture(independent = scenario == "independent_conditions")
    celltypes <- scenario == "paired_celltypes"
    request <- validate_hybrid_request(data, list(analysis_scope = "across_donors",
      comparison_axis = if (celltypes) "celltype" else "condition",
      design = if (scenario == "independent_conditions") "independent" else "paired",
      donors = c("D1", "D2", "D3"), celltype = if (celltypes) "T" else "B", condition = "NC",
      target = if (celltypes) "T" else "TRT", reference = if (celltypes) "B" else "NC",
      readout = "abundance", markers = c("A", "B")))
    result <- run_donor_analysis(data, request)
    entry <- list(id = "donor-comparison", label = scenario, request = request, result = result)
    contract <- hybrid_result_contract(data, entry)
    expect_identical(contract$schema_version, 2L)
    expect_identical(contract$differential_results$abundance$effect_size, result$results$effect_size)
    expect_identical(contract$donor_values$abundance$value, result$donor_values$value)
    expect_identical(contract$donor_values$abundance$source_celltype, result$donor_values$source_celltype)
    expect_identical(contract$donor_contrasts$abundance$effect_size, result$donor_contrasts$effect_size)
    expect_true(all(contract$differential_results$abundance$replicate_unit == "donor"))
    expect_true(all(contract$donor_values$abundance$comparison_axis == request$comparison_axis))
    expect_identical(data$contract$sample_summaries, contract$sample_summaries)
    adapted <- data
    adapted$contract <- contract
    expect_true(validate_precomputed_data(adapted))
    invalid <- entry
    invalid$result$donor_values$donor[[1L]] <- "outside scope"
    expect_error(hybrid_result_contract(data, invalid), "measurements do not match")
    invalid <- entry
    invalid$result$donor_values$group[[1L]] <- "outside group"
    expect_error(hybrid_result_contract(data, invalid), "invalid groups")
    invalid <- entry
    invalid$result$donor_values <- rbind(invalid$result$donor_values, invalid$result$donor_values[1L, ])
    expect_error(hybrid_result_contract(data, invalid), "duplicate records")
    invalid <- entry
    invalid$result$donor_values <- invalid$result$donor_values[FALSE, ]
    invalid$result$donor_contrasts <- invalid$result$donor_contrasts[FALSE, ]
    expect_error(hybrid_result_contract(data, invalid), "cover every requested feature")
    invalid <- entry
    invalid$result$donor_values <- invalid$result$donor_values[-1L, ]
    expect_error(hybrid_result_contract(data, invalid), "cover every requested feature")
    invalid <- entry
    invalid$result$donor_values <- invalid$result$donor_values[invalid$result$donor_values$feature_id != "A", ]
    expect_error(hybrid_result_contract(data, invalid), "cover every requested feature")
    invalid <- entry
    invalid$result$donor_values$n_total_cells[[1L]] <- 41L
    expect_error(hybrid_result_contract(data, invalid), "inconsistent cell support")
    invalid <- entry
    invalid$result$donor_values$n_evaluable_cells[[1L]] <- 41L
    expect_error(hybrid_result_contract(data, invalid), "inconsistent cell support")
    invalid <- entry
    invalid$result$donor_values$n_scored_cells[[1L]] <- 43L
    expect_error(hybrid_result_contract(data, invalid), "inconsistent cell support")
    invalid <- entry
    invalid$result$donor_values$eligible[[1L]] <- !invalid$result$donor_values$eligible[[1L]]
    expect_error(hybrid_result_contract(data, invalid), "inconsistent cell support")
    invalid <- entry
    invalid$result$results$effect_size[[1L]] <- invalid$result$results$effect_size[[1L]] + 1
    expect_error(hybrid_result_contract(data, invalid), "effects or support do not match")
    invalid <- entry
    invalid$result$results$n_cells_a[[1L]] <- invalid$result$results$n_cells_a[[1L]] + 1L
    expect_error(hybrid_result_contract(data, invalid), "effects or support do not match")
    invalid <- entry
    invalid$result$results$n_donors_a[[1L]] <- invalid$result$results$n_donors_a[[1L]] + 1L
    expect_error(hybrid_result_contract(data, invalid), "effects or support do not match")
    invalid <- entry
    invalid$result$results$effect_size[[1L]] <- NA_real_
    expect_error(hybrid_result_contract(data, invalid), "effects or support do not match")
    invalid <- entry
    invalid$result$results$inference[[1L]] <- "donor_level"
    invalid$result$results$p_value[[1L]] <- .01
    invalid$result$results$p_adj[[1L]] <- .02
    invalid$result$results$status[[1L]] <- "no_finite_test_result"
    expect_error(hybrid_result_contract(data, invalid), "inference and test status")
    if (request$design == "paired") {
      invalid <- entry
      invalid$result$donor_contrasts <- invalid$result$donor_contrasts[-1L, ]
      expect_error(hybrid_result_contract(data, invalid), "every complete eligible donor pair")
      invalid <- entry
      invalid$result$donor_contrasts$value_a[[1L]] <- 99
      expect_error(hybrid_result_contract(data, invalid), "do not match their donor measurements")
      invalid <- entry
      invalid$result$donor_contrasts$n_cells_a[[1L]] <- 41L
      expect_error(hybrid_result_contract(data, invalid), "do not match their donor measurements")
    } else {
      invalid <- entry
      invalid$result$donor_values <- invalid$result$donor_values[invalid$result$donor_values$n_total_cells > 0, ]
      expect_error(hybrid_result_contract(data, invalid), "cover every requested feature")
    }
  }
})

test_that("engine exports cannot claim inference below the donor test minimum", {
  for (design in c("paired", "independent")) {
    data <- hybrid_module_donor_fixture(independent = design == "independent")
    request <- validate_hybrid_request(data, list(analysis_scope = "across_donors", comparison_axis = "condition",
      design = design, donors = if (design == "paired") "D1" else c("D1", "D2"), celltype = "B",
      target = "TRT", reference = "NC", readout = "abundance"))
    entry <- list(id = "small-comparison", label = "Insufficient donors", request = request,
      result = run_donor_analysis(data, request))
    expect_silent(hybrid_result_contract(data, entry))
    entry$result$results$inference[[1L]] <- "donor_level"
    entry$result$results$status[[1L]] <- "tested"
    entry$result$results$p_value[[1L]] <- .01
    entry$result$results$p_adj[[1L]] <- .02
    expect_error(hybrid_result_contract(data, entry), "inference and test status")
  }
})

test_that("across-donor history and exports retain the submitted design after controls change", {
  current <- shiny::reactiveVal(hybrid_module_donor_fixture())
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
    code_dir = tempdir()), {
    session$setInputs(analysis_scope = "across_donors", comparison_axis = "condition", design = "paired",
      donors = c("D1", "D2"), celltype = "B", target = "TRT", reference = "NC", readout = "abundance", run = 1L)
    submitted <- probe$state$jobs[[1L]]$request
    session$setInputs(donors = "D3", target = "NC", reference = "TRT", design = "independent")
    probe$state$jobs[[1L]]$resolve(run_donor_analysis(current(), submitted))
    settle_module_task(session, task)
    expect_length(history(), 1L)
    expect_identical(selected_entry()$request, submitted)
    expect_match(selected_entry()$label, "2 donors paired", fixed = TRUE)
    expect_match(output$result_scope, "Across donors: D1, D2 paired", fixed = TRUE)
    package <- readRDS(output$result_rds)
    expect_identical(package$entry$request, submitted)
    expect_identical(package$entry$result$donor_values, result()$donor_values)
    expect_match(paste(readLines(output$provenance_json), collapse = ""), '"analysis_scope": "across_donors"', fixed = TRUE)
    expect_equal(nrow(utils::read.csv(output$results_csv)), 3L)
  })
})

test_that("across-donor abundance uses normalized values without requiring raw counts", {
  data <- hybrid_module_donor_fixture()
  data$counts <- NULL
  current <- shiny::reactiveVal(data)
  probe <- hybrid_module_probe()
  shiny::testServer(probe$module, args = list(data = current, manager = new_hybrid_job_manager(),
    code_dir = tempdir()), {
    session$setInputs(analysis_scope = "across_donors", comparison_axis = "condition", design = "paired",
      donors = c("D1", "D2"), celltype = "B", target = "TRT", reference = "NC", readout = "abundance")
    expect_null(prepared_request()$error)
    expect_match(output$readiness$html, "does not require counts", fixed = TRUE)
    session$setInputs(run = 1L)
    expect_length(probe$state$jobs, 1L)
    submitted <- probe$state$jobs[[1L]]$request
    probe$state$jobs[[1L]]$resolve(run_donor_analysis(current(), submitted))
    settle_module_task(session, task)
    contract <- hybrid_result_contract(current(), selected_entry())
    expect_true(all(contract$differential_results$abundance$metric == "delta_mean_donor_mean_normalized_abundance"))
    session$setInputs(analysis_scope = "within_donor", donor = "D1")
    expect_match(prepared_request()$error, "Raw counts are not stored")
  })
})

test_that("donor result validation uses resolved aliases without changing source mappings", {
  source <- hybrid_module_donor_fixture()
  source$bundle_version <- PROXIOMEVIS_BUNDLE_VERSION
  source$metadata$subject <- source$metadata$donor
  source$metadata$donor <- source$metadata$celltype
  source$metadata$sample_source <- source$metadata$sample_id
  source$metadata$sample_id <- source$metadata$condition
  source$metadata$class <- source$metadata$cluster
  mappings <- list(donor_column = "subject", celltype_column = "donor", condition_column = "sample_id",
    sample_column = "sample_source", cluster_column = "class")
  source$contract$defaults[names(mappings)] <- mappings
  data <- precomputed_app_data(source)
  request <- validate_hybrid_request(data, list(analysis_scope = "across_donors", comparison_axis = "condition",
    design = "paired", donors = c("D1", "D2"), celltype = "B", target = "TRT", reference = "NC", readout = "abundance"))
  entry <- list(id = "mapped-comparison", label = "Mapped source", request = request,
    result = run_donor_analysis(data, request))
  contract <- hybrid_result_contract(data, entry)
  expect_identical(contract$defaults, source$contract$defaults)
  expect_true(all(contract$donor_values$abundance$source_celltype == "B"))
  expect_setequal(contract$donor_values$abundance$donor, c("D1", "D2"))
})

test_that("documentation covers all comparison designs and distinguishes cells from donors", {
  html <- as.character(documentation_ui())
  for (text in c("Within-donor exploration, Conditions:", "Within-donor exploration, Cell types:",
    "Across donors, Conditions, Paired:", "Across donors, Cell types, Paired:",
    "Across donors, Conditions, Independent:", "Biological replicates are not necessary for this exploratory view.",
    "Version 1 results keep their original methods and values.", "Version 2 adds comparison design and donor measurements.",
    "one point per donor and group", "a violin and boxplot of available cell values", "paired t-test", "Welch t-test",
    "at least two complete eligible donor pairs", "at least two eligible donors in each group",
    "do not show that a study has sufficient statistical power", "keep P values as NA", "approximately normal donor differences",
    "Raw counts are not necessary for this mode.", "mean stored normalized abundance", "median retained observed log2 ratio")) {
    expect_match(html, text, fixed = TRUE)
  }
  expect_false(grepl("<p>P values describe cells within the selected donor.", html, fixed = TRUE))
})
