source(file.path("..", "..", "R", "hybrid_jobs.R"))

hybrid_worker_fixture <- function() {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  bundle <- file.path(root, "bundle.rds")
  saveRDS(list(dataset_id = "worker-fixture"), bundle)
  for (file in c("precomputed_contract.R", "proximity_store.R")) writeLines("", file.path(root, file))
  writeLines(c(
    "run_hybrid_analysis <- function(bundle_path, request) {",
    "  if (identical(request$mode, 'error')) stop('Worker fixture error.', call. = FALSE)",
    "  if (!is.null(request$delay)) Sys.sleep(request$delay)",
    "  list(pid = Sys.getpid(), bundle_path = bundle_path, request = request)",
    "}"), file.path(root, "hybrid_analysis.R"))
  list(bundle = bundle, code_dir = root)
}

await_hybrid_job <- function(job, timeout = 30) {
  state <- new.env(parent = emptyenv())
  state$done <- FALSE
  promises::then(job$promise, onFulfilled = function(value) {
    state$value <- value
    state$done <- TRUE
  }, onRejected = function(error) {
    state$error <- error
    state$done <- TRUE
  })
  deadline <- Sys.time() + timeout
  while (!state$done && Sys.time() < deadline) later::run_now(0.05)
  if (!state$done) {
    job$cancel()
    stop("The worker test did not settle within its deadline.")
  }
  as.list(state)
}

test_that("analysis runs in a separate R process while the event loop remains responsive", {
  fixture <- hybrid_worker_fixture()
  manager <- new_hybrid_job_manager(timeout = 20)
  job <- start_hybrid_job(fixture$bundle, list(delay = 0.5), fixture$code_dir, manager)
  withr::defer(job$cancel())
  expect_identical(manager$active, 1L)
  expect_true(job$process$is_alive())
  callback_while_running <- FALSE
  later::later(function() callback_while_running <<- job$process$is_alive(), delay = 0.01)
  completed <- await_hybrid_job(job)
  expect_null(completed$error)
  expect_true(callback_while_running)
  expect_false(identical(completed$value$pid, Sys.getpid()))
  expect_identical(completed$value$bundle_path, normalizePath(fixture$bundle))
  expect_identical(manager$active, 0L)
  expect_false(job$process$is_alive())
})

test_that("worker errors reject their promise and release capacity", {
  fixture <- hybrid_worker_fixture()
  manager <- new_hybrid_job_manager(timeout = 20)
  job <- start_hybrid_job(fixture$bundle, list(mode = "error"), fixture$code_dir, manager)
  withr::defer(job$cancel())
  completed <- await_hybrid_job(job)
  expect_match(conditionMessage(completed$error), "Worker fixture error")
  expect_identical(manager$active, 0L)
  expect_false(job$process$is_alive())
  expect_error(start_hybrid_job(file.path(fixture$code_dir, "absent.rds"), list(), fixture$code_dir, manager))
  expect_identical(manager$active, 0L)
})

test_that("shared concurrency limits and repeated cancellation do not leak workers or slots", {
  fixture <- hybrid_worker_fixture()
  manager <- new_hybrid_job_manager(max_jobs = 1L, timeout = 20)
  job <- start_hybrid_job(fixture$bundle, list(delay = 5), fixture$code_dir, manager)
  withr::defer(job$cancel())
  expect_error(start_hybrid_job(fixture$bundle, list(), fixture$code_dir, manager), "worker is busy")
  expect_identical(manager$active, 1L)
  job$cancel()
  job$cancel()
  completed <- await_hybrid_job(job)
  expect_match(conditionMessage(completed$error), "cancelled")
  expect_false(job$process$is_alive())
  expect_identical(manager$active, 0L)
  replacement <- start_hybrid_job(fixture$bundle, list(), fixture$code_dir, manager)
  withr::defer(replacement$cancel())
  expect_null(await_hybrid_job(replacement)$error)
  job$cancel()
  later::run_now(0.25)
  expect_identical(manager$active, 0L)
  expect_false(replacement$process$is_alive())
})

test_that("timed-out workers are terminated and their promise rejects", {
  fixture <- hybrid_worker_fixture()
  manager <- new_hybrid_job_manager(timeout = 0.1)
  job <- start_hybrid_job(fixture$bundle, list(delay = 5), fixture$code_dir, manager)
  withr::defer(job$cancel())
  completed <- await_hybrid_job(job)
  expect_match(conditionMessage(completed$error), "time limit")
  expect_identical(manager$active, 0L)
  expect_false(job$process$is_alive())
})
