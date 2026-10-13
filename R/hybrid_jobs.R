new_hybrid_job_manager <- function(max_jobs = 1L, timeout = 120) {
  stopifnot(length(max_jobs) == 1L, is.finite(max_jobs), max_jobs >= 1L,
            length(timeout) == 1L, is.finite(timeout), timeout > 0)
  manager <- new.env(parent = emptyenv())
  manager$active <- 0L
  manager$max_jobs <- as.integer(max_jobs)
  manager$timeout <- timeout
  manager
}

# Only paths and a bounded request cross the process boundary, never the Seurat object.
start_hybrid_job <- function(bundle_path, request, code_dir, manager) {
  if (manager$active >= manager$max_jobs) {
    stop("The analysis worker is busy. Please retry after the current job finishes.", call. = FALSE)
  }
  bundle_path <- normalizePath(bundle_path, mustWork = TRUE)
  code_dir <- normalizePath(code_dir, mustWork = TRUE)
  process <- callr::r_bg(
    function(bundle_path, request, code_dir) {
      for (file in c("precomputed_contract.R", "proximity_store.R", "hybrid_analysis.R")) {
        sys.source(file.path(code_dir, file), envir = globalenv())
      }
      run_hybrid_analysis(bundle_path, request)
    },
    args = list(bundle_path = bundle_path, request = request, code_dir = code_dir),
    env = c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", ARROW_NUM_THREADS = "1"),
    supervise = TRUE
  )
  manager$active <- manager$active + 1L
  started <- Sys.time()
  cancelled <- FALSE
  result <- promises::promise(function(resolve, reject) {
    poll <- function() {
      tryCatch({
        expired <- as.numeric(difftime(Sys.time(), started, units = "secs")) > manager$timeout
        if (cancelled || expired) {
          if (process$is_alive()) process$kill()
          stop(if (cancelled) "Analysis cancelled." else
            "Analysis exceeded the time limit. Use precomputed results or run this comparison offline.", call. = FALSE)
        }
        if (process$is_alive()) {
          later::later(poll, delay = 0.2)
        } else {
          value <- process$get_result()
          manager$active <- manager$active - 1L
          resolve(value)
        }
      }, error = function(error) {
        if (process$is_alive()) process$kill()
        manager$active <- manager$active - 1L
        reject(error)
      })
    }
    later::later(poll, delay = 0.1)
  })
  list(promise = result, process = process, cancel = function() {
    cancelled <<- TRUE
    if (process$is_alive()) process$kill()
    invisible(NULL)
  })
}
