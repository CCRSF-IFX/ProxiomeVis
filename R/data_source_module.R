PROXIOMEVIS_APP_NAME <- ".ProxiomeVis"

ensure_writable_directory <- function(path, label, create = TRUE) {
  path <- path.expand(as.character(path %||% ""))
  if (!nzchar(path)) {
    stop(label, " path is empty.", call. = FALSE)
  }

  if (isTRUE(create)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path)) {
    stop(label, " does not exist: ", path, call. = FALSE)
  }

  path <- normalizePath(path, mustWork = TRUE)
  if (file.access(path, 2) != 0) {
    stop(label, " is not writable: ", path, call. = FALSE)
  }

  path
}

proxiomevis_home_dir <- function(create = TRUE) {
  writable_dir <- Sys.getenv("PROXIOMEVIS_HOME", unset = "")
  if (!nzchar(writable_dir)) {
    home_dir <- Sys.getenv("HOME", unset = "")
    if (!nzchar(home_dir)) {
      home_dir <- path.expand("~")
    }
    if (!nzchar(home_dir) || identical(home_dir, "~")) {
      stop("Could not resolve HOME for the ProxiomeVis writable directory.", call. = FALSE)
    }
    writable_dir <- file.path(home_dir, PROXIOMEVIS_APP_NAME)
  }

  writable_dir <- ensure_writable_directory(writable_dir, "ProxiomeVis writable directory", create = create)
  Sys.setenv(PROXIOMEVIS_HOME = writable_dir)
  writable_dir
}

proxiomevis_runtime_dir <- function(create = TRUE) {
  runtime_dir <- Sys.getenv("PROXIOME_RUNTIME_DIR", unset = "")
  if (!nzchar(runtime_dir)) {
    runtime_dir <- file.path(proxiomevis_home_dir(create = create), "runtime")
  }

  runtime_dir <- ensure_writable_directory(runtime_dir, "ProxiomeVis runtime directory", create = create)
  Sys.setenv(PROXIOME_RUNTIME_DIR = runtime_dir)
  runtime_dir
}

proxiomevis_cache_dir <- function(create = TRUE) {
  ensure_writable_directory(
    file.path(proxiomevis_home_dir(create = create), "cache"),
    "ProxiomeVis cache directory",
    create = create
  )
}

default_demo_cache_path <- function(app_dir = APP_DIR) {
  configured_cache <- Sys.getenv("PROXIOME_DEMO_CACHE_PATH", unset = "")
  if (nzchar(configured_cache)) {
    dir.create(dirname(configured_cache), recursive = TRUE, showWarnings = FALSE)
    return(normalizePath(configured_cache, mustWork = FALSE))
  }

  bundled_cache <- file.path(app_dir, "cache", "demo_proxiome_data.rds")
  if (file.exists(bundled_cache)) {
    return(bundled_cache)
  }

  file.path(proxiomevis_cache_dir(), "demo_proxiome_data.rds")
}

user_rds_source_signature <- function(rds_path) {
  rds_path <- normalizePath(trimws(as.character(rds_path %||% "")), mustWork = FALSE)
  info <- file.info(rds_path)

  list(
    path = rds_path,
    size = unname(info$size[1]),
    mtime = unname(as.numeric(info$mtime[1]))
  )
}

user_rds_cache_path <- function(rds_path, cache_dir = proxiomevis_cache_dir()) {
  require_namespace("digest")

  cache_dir <- ensure_writable_directory(cache_dir, "ProxiomeVis cache directory")
  signature <- user_rds_source_signature(rds_path)
  key <- digest::digest(signature, algo = "xxhash64")

  file.path(cache_dir, paste0("user-rds-", key, ".rds"))
}

user_rds_path_loading_enabled <- function(platform = APP_PLATFORM) {
  disabled <- tolower(Sys.getenv("PROXIOME_DISABLE_USER_RDS", unset = "false")) %in% c("1", "true", "yes")
  platform %in% c("ccrsf_hpc", "biowulf_hpc", "portable") && !disabled
}

data_source_controls <- function(id, platform = APP_PLATFORM) {
  if (!user_rds_path_loading_enabled(platform)) {
    return(tagList())
  }

  ns <- NS(id)
  popover(
    actionButton(
      ns("data_source_menu"),
      "Data",
      class = "btn btn-outline-secondary btn-sm data-source-button"
    ),
    div(
      class = "data-source-popover",
      textInput(
        ns("rds_server_path"),
        "RDS path",
        placeholder = "/path/to/sample.rds"
      ),
      actionButton(
        ns("validate_rds_path"),
        "Validate RDS",
        class = "btn btn-outline-secondary btn-sm"
      ),
      actionButton(
        ns("load_rds_path"),
        "Load Data",
        class = "btn btn-outline-primary btn-sm"
      ),
      actionButton(
        ns("use_demo_data"),
        "Use demo data",
        class = "btn btn-outline-secondary btn-sm"
      ),
      actionButton(
        ns("use_raji_demo_data"),
        "Use Raji/CAR-T data",
        class = "btn btn-outline-secondary btn-sm"
      ),
      uiOutput(ns("rds_schema_report")),
      div(
        id = ns("rds_load_status"),
        class = "rds-load-status idle",
        "Enter an RDS path, then click Load Data."
      ),
      div(
        id = ns("rds_load_progress"),
        class = "rds-load-progress idle",
        role = "progressbar",
        `aria-hidden` = "true",
        `aria-valuemin` = "0",
        `aria-valuemax` = "100",
        div(
          class = "rds-load-progress-track",
          div(id = ns("rds_load_progress_bar"), class = "rds-load-progress-bar", `aria-valuenow` = "0")
        ),
        div(id = ns("rds_load_elapsed"), class = "rds-load-elapsed")
      )
    ),
    title = "Data source",
    placement = "bottom",
    options = list(customClass = "data-source-popover-shell")
  )
}

data_source_module_ui <- function(id, platform = APP_PLATFORM) {
  ns <- NS(id)
  controls <- if (user_rds_path_loading_enabled(platform)) {
    data_source_controls(id, platform = platform)
  } else {
    tagList()
  }

  nav_item(
    controls,
    tags$span(class = "source-chip", textOutput(ns("source_summary"), inline = TRUE))
  )
}

load_default_proxiome_data <- function(app_dir = APP_DIR) {
  data <- load_demo_proxiome_data(
    cache_path = default_demo_cache_path(app_dir)
  )
  data$source$source_type <- "demo"
  data$source$display_name <- basename(data$source$rds_path)
  data
}

validate_rds_file_path <- function(rds_path) {
  rds_path <- trimws(as.character(rds_path %||% ""))
  if (!nzchar(rds_path)) {
    stop("Enter an RDS path.", call. = FALSE)
  }
  if (!grepl("\\.rds$", rds_path, ignore.case = TRUE)) {
    stop("RDS path must be an .rds file.", call. = FALSE)
  }
  if (!file.exists(rds_path)) {
    stop("RDS path does not exist.", call. = FALSE)
  }

  TRUE
}

format_elapsed_time <- function(seconds) {
  seconds <- suppressWarnings(as.numeric(seconds[1]))
  if (!is.finite(seconds) || seconds < 0) {
    seconds <- 0
  }

  seconds <- floor(seconds)
  minutes <- seconds %/% 60
  remaining_seconds <- seconds %% 60
  hours <- minutes %/% 60
  remaining_minutes <- minutes %% 60

  if (hours > 0) {
    return(sprintf("%dh %02dm %02ds", hours, remaining_minutes, remaining_seconds))
  }
  if (minutes > 0) {
    return(sprintf("%dm %02ds", minutes, remaining_seconds))
  }
  paste0(seconds, "s")
}

rds_load_progress_file <- function(session_token = NULL) {
  token <- gsub("[^A-Za-z0-9_-]+", "-", as.character(session_token %||% "session"))
  file.path(
    proxiomevis_runtime_dir(),
    paste0("rds-load-progress-", token, "-", format(Sys.time(), "%Y%m%d%H%M%S"), ".rds")
  )
}

write_rds_load_progress <- function(
  progress_path,
  state = "running",
  stage = "",
  message = "",
  value = 0,
  started_at = Sys.time(),
  finished_at = NULL
) {
  if (is.null(progress_path) || !nzchar(progress_path)) {
    return(invisible(NULL))
  }

  started_at <- as.POSIXct(started_at, origin = "1970-01-01")
  finished_at <- if (is.null(finished_at)) NULL else as.POSIXct(finished_at, origin = "1970-01-01")
  now <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(finished_at %||% now, started_at, units = "secs"))
  value <- suppressWarnings(as.numeric(value[1]))
  if (!is.finite(value)) {
    value <- 0
  }
  value <- max(0, min(1, value))

  progress <- list(
    state = state,
    stage = stage,
    message = message,
    value = value,
    percent = round(value * 100),
    started_at = started_at,
    finished_at = finished_at,
    elapsed_seconds = elapsed_seconds,
    updated_at = now
  )

  dir.create(dirname(progress_path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- tempfile(pattern = basename(progress_path), tmpdir = dirname(progress_path), fileext = ".tmp")
  saveRDS(progress, tmp_path)
  file.rename(tmp_path, progress_path)
  invisible(progress)
}

read_rds_load_progress <- function(progress_path) {
  if (is.null(progress_path) || !nzchar(progress_path) || !file.exists(progress_path)) {
    return(NULL)
  }

  tryCatch(
    readRDS(progress_path),
    error = function(error) NULL
  )
}

rds_load_elapsed_seconds <- function(progress, now = Sys.time()) {
  if (is.null(progress)) {
    return(0)
  }

  started_at <- progress$started_at %||% NULL
  finished_at <- progress$finished_at %||% NULL
  if (!is.null(started_at)) {
    started_at <- as.POSIXct(started_at, origin = "1970-01-01")
  }
  if (!is.null(finished_at)) {
    finished_at <- as.POSIXct(finished_at, origin = "1970-01-01")
  }

  if (!is.null(started_at) && identical(progress$state, "running")) {
    return(as.numeric(difftime(as.POSIXct(now, origin = "1970-01-01"), started_at, units = "secs")))
  }
  if (!is.null(started_at) && !is.null(finished_at)) {
    return(as.numeric(difftime(finished_at, started_at, units = "secs")))
  }

  progress$elapsed_seconds %||% 0
}

rds_load_progress_message <- function(progress, now = Sys.time()) {
  if (is.null(progress)) {
    return("")
  }

  elapsed_label <- if (identical(progress$state, "success")) "Loaded in" else "Elapsed"
  paste0(
    progress$message %||% "",
    " ",
    elapsed_label,
    ": ",
    format_elapsed_time(rds_load_elapsed_seconds(progress, now = now))
  )
}

rds_load_elapsed_label <- function(progress, now = Sys.time()) {
  if (is.null(progress)) {
    return("")
  }

  label <- if (identical(progress$state, "success")) "Loaded in" else "Elapsed"
  paste0(label, ": ", format_elapsed_time(rds_load_elapsed_seconds(progress, now = now)))
}

make_rds_load_progress_callback <- function(progress_path, started_at = Sys.time()) {
  force(progress_path)
  force(started_at)

  function(stage, message, value, state = "running") {
    write_rds_load_progress(
      progress_path,
      state = state,
      stage = stage,
      message = message,
      value = value,
      started_at = started_at
    )
  }
}

load_user_proxiome_data <- function(rds_path, app_dir = APP_DIR, progress_path = NULL) {
  started_at <- Sys.time()
  progress_callback <- make_rds_load_progress_callback(progress_path, started_at = started_at)
  progress_callback("validate", "Validating RDS path...", 0.03)
  validate_rds_file_path(rds_path)
  progress_callback("library", "Preparing R library environment...", 0.08)
  configure_app_r_libs(app_dir)
  rds_path <- normalizePath(trimws(as.character(rds_path)), mustWork = TRUE)
  cache_path <- user_rds_cache_path(rds_path)
  if (file.exists(cache_path)) {
    progress_callback("cache", "Reading processed app data cache...", 0.12)
  } else {
    progress_callback("cache", "Preparing processed app data cache...", 0.12)
  }

  data <- load_demo_proxiome_data(
    rds_path = rds_path,
    marker_selection = "all",
    cache_path = cache_path,
    force = FALSE,
    progress_callback = progress_callback
  )
  data$source$source_type <- "user_rds"
  data$source$display_name <- basename(rds_path)
  data$source$user_cache_path <- cache_path
  finished_at <- Sys.time()
  data$source$load_elapsed_seconds <- as.numeric(difftime(finished_at, started_at, units = "secs"))
  write_rds_load_progress(
    progress_path,
    state = "success",
    stage = "complete",
    message = paste0("Loaded ", basename(rds_path), "."),
    value = 1,
    started_at = started_at,
    finished_at = finished_at
  )
  data
}

data_source_summary <- function(data) {
  if (is.null(data)) {
    return("Loading demo RDS...")
  }

  source <- data$source %||% list()
  source_name <- source$display_name %||% basename(source$rds_path %||% "RDS")
  n_cells <- source$n_cells %||% NA_integer_
  n_markers <- source$n_markers %||% NA_integer_

  paste0(
    source_name,
    " | ",
    format(n_cells, big.mark = ","),
    " cells | ",
    format(n_markers, big.mark = ","),
    " markers"
  )
}

configure_user_rds_future_plan <- function(workers = Sys.getenv("PROXIOME_USER_RDS_WORKERS", unset = "1")) {
  workers <- suppressWarnings(as.integer(workers[1]))
  if (!is.finite(workers) || workers < 1) {
    workers <- 1L
  }

  future::plan(future::multisession, workers = workers)
  invisible(workers)
}

create_user_rds_load_task <- function(app_dir = APP_DIR) {
  app_dir <- normalizePath(app_dir, mustWork = TRUE)
  configure_user_rds_future_plan()

  ExtendedTask$new(function(rds_path, progress_path = NULL) {
    promises::future_promise({
      configure_app_r_libs(app_dir)
      load_user_proxiome_data(rds_path, app_dir = app_dir, progress_path = progress_path)
    })
  })
}

data_source_module_server <- function(id, app_dir = APP_DIR) {
  moduleServer(id, function(input, output, session) {
    app_dir <- normalizePath(app_dir, mustWork = TRUE)
    demo_data <- reactiveVal(NULL)
    user_rds_load_task <- if (user_rds_path_loading_enabled()) {
      create_user_rds_load_task(app_dir)
    } else {
      NULL
    }
    rds_load_state <- reactiveVal("idle")
    rds_load_message <- reactiveVal("Enter an RDS path, then click Load Data.")
    rds_load_path_label <- reactiveVal("")
    rds_load_progress_path <- reactiveVal(NULL)
    rds_schema <- reactiveVal(NULL)

    send_rds_load_state <- function(
      state,
      message,
      button_label = "Load Data",
      progress = NULL,
      elapsed_label = NULL
    ) {
      state <- state %||% "idle"
      message <- message %||% ""
      is_running <- identical(state, "running")

      rds_load_state(state)
      rds_load_message(message)
      session$sendCustomMessage("proxiome-rds-load-state", list(
        state = state,
        message = message,
        disabled = is_running,
        buttonLabel = if (is_running) "Loading..." else button_label,
        progress = progress %||% if (identical(state, "success")) 100 else 0,
        elapsedLabel = elapsed_label %||% ""
      ))
    }

    load_demo_into_app <- function(show_loaded_notification = FALSE) {
      withProgress(message = "Loading Pixelator v4.1.1 demo data", value = 0.1, {
        data <- load_default_proxiome_data(app_dir)
        incProgress(0.8)
        demo_data(data)
      })

      if (isTRUE(show_loaded_notification)) {
        showNotification("Demo data loaded.", type = "message")
      }
    }

    load_raji_demo_into_app <- function(show_loaded_notification = FALSE) {
      withProgress(message = "Loading Raji/CAR-T demo data", value = 0.1, {
        data <- load_raji_demo_proxiome_data(app_dir)
        incProgress(0.8)
        demo_data(data)
      })

      if (isTRUE(show_loaded_notification)) {
        showNotification("Raji/CAR-T demo data loaded.", type = "message")
      }
    }

    observeEvent(TRUE, {
      load_demo_into_app()
    }, once = TRUE)

    observeEvent(input$use_demo_data, {
      req(user_rds_path_loading_enabled())
      load_demo_into_app(show_loaded_notification = TRUE)
      send_rds_load_state("success", "Demo data is active.", progress = 100)
      rds_load_path_label("")
      rds_load_progress_path(NULL)
      rds_schema(NULL)
    }, ignoreInit = TRUE)

    observeEvent(input$use_raji_demo_data, {
      req(user_rds_path_loading_enabled())
      tryCatch(
        {
          load_raji_demo_into_app(show_loaded_notification = TRUE)
          send_rds_load_state("success", "Raji/CAR-T demo data is active.", progress = 100)
          rds_load_path_label("")
          rds_load_progress_path(NULL)
          rds_schema(NULL)
        },
        error = function(error) {
          send_rds_load_state("error", paste("Could not load Raji/CAR-T demo data:", conditionMessage(error)), progress = 0)
          showNotification(
            paste("Could not load Raji/CAR-T demo data:", conditionMessage(error)),
            type = "error",
            duration = NULL
          )
        }
      )
    }, ignoreInit = TRUE)

    observeEvent(input$validate_rds_path, {
      req(user_rds_path_loading_enabled())

      tryCatch(
        {
          schema <- inspect_user_rds_schema(input$rds_server_path)
          rds_schema(schema)
          send_rds_load_state(
            "success",
            paste0("Validated ", basename(schema$path), ". Review the schema summary before loading."),
            button_label = "Load Data",
            progress = 0
          )
          showNotification("RDS schema validated.", type = "message", duration = 6)
        },
        error = function(error) {
          rds_schema(NULL)
          send_rds_load_state("error", paste("Could not validate RDS:", conditionMessage(error)), progress = 0)
          showNotification(
            paste("Could not validate RDS:", conditionMessage(error)),
            type = "error",
            duration = NULL
          )
        }
      )
    }, ignoreInit = TRUE)

    observeEvent(input$load_rds_path, {
      req(user_rds_path_loading_enabled())
      req(!is.null(user_rds_load_task))

      tryCatch(
        {
          validate_rds_file_path(input$rds_server_path)
          rds_schema(inspect_user_rds_schema(input$rds_server_path))
          current_status <- isolate(user_rds_load_task$status())
          current_state <- isolate(rds_load_state())
          if (identical(current_state, "running") || identical(current_status, "running")) {
            send_rds_load_state("running", paste(
              "RDS load is already running.",
              "Large files can take several minutes; wait for the current load to finish."
            ))
            showNotification(
              "RDS load is already running. Wait for it to finish before clicking Load Data again.",
              type = "warning",
              duration = 8
            )
            return(invisible(NULL))
          }

          rds_path <- normalizePath(trimws(as.character(input$rds_server_path)), mustWork = TRUE)
          progress_path <- rds_load_progress_file(session$token %||% paste0("session-", as.integer(Sys.time())))
          rds_load_progress_path(progress_path)
          write_rds_load_progress(
            progress_path,
            state = "running",
            stage = "start",
            message = paste0("Starting RDS load for ", basename(rds_path), "."),
            value = 0.03,
            started_at = Sys.time()
          )
          rds_load_path_label(basename(rds_path))
          send_rds_load_state("running", paste0(
            "Loading ",
            basename(rds_path),
            " in the background. Large files can take several minutes."
          ), progress = 3, elapsed_label = "Elapsed: 0s")
          user_rds_load_task$invoke(input$rds_server_path, progress_path)
          showNotification(
            paste("Loading", basename(rds_path), "in the background."),
            type = "message",
            duration = 8
          )
        },
        error = function(error) {
          send_rds_load_state("error", paste("Could not start RDS load:", conditionMessage(error)), progress = 0)
          showNotification(
            paste("Could not load RDS path:", conditionMessage(error)),
            type = "error",
            duration = NULL
          )
        }
      )
    }, ignoreInit = TRUE)

    observe({
      req(identical(rds_load_state(), "running"))
      invalidateLater(1000, session)

      progress <- read_rds_load_progress(rds_load_progress_path())
      req(!is.null(progress))
      display_state <- progress$state %||% "running"
      if (identical(display_state, "success") && !identical(user_rds_load_task$status(), "success")) {
        display_state <- "running"
      }
      send_rds_load_state(
        state = display_state,
        message = rds_load_progress_message(progress),
        progress = progress$percent %||% round((progress$value %||% 0) * 100),
        elapsed_label = rds_load_elapsed_label(progress)
      )
    })

    observe({
      req(!is.null(user_rds_load_task))
      status <- user_rds_load_task$status()

      if (identical(status, "error")) {
        tryCatch(
          user_rds_load_task$result(),
          error = function(error) {
            label <- rds_load_path_label()
            if (!nzchar(label)) {
              label <- "RDS path"
            }
            progress <- read_rds_load_progress(rds_load_progress_path())
            write_rds_load_progress(
              rds_load_progress_path(),
              state = "error",
              stage = "error",
              message = paste0("Could not load ", label, ": ", conditionMessage(error)),
              value = progress$value %||% 0,
              started_at = progress$started_at %||% Sys.time(),
              finished_at = Sys.time()
            )
            error_progress <- read_rds_load_progress(rds_load_progress_path())
            send_rds_load_state(
              "error",
              rds_load_progress_message(error_progress),
              progress = error_progress$percent %||% 0,
              elapsed_label = rds_load_elapsed_label(error_progress)
            )
            showNotification(
              paste("Could not load RDS path:", conditionMessage(error)),
              type = "error",
              duration = NULL
            )
          }
        )
        return(invisible(NULL))
      }

      req(identical(status, "success"))
      data <- user_rds_load_task$result()
      demo_data(data)
      label <- data$source$display_name %||% rds_load_path_label()
      if (!nzchar(label)) {
        label <- "RDS path"
      }
      progress <- read_rds_load_progress(rds_load_progress_path())
      success_message <- rds_load_progress_message(progress)
      if (!nzchar(success_message)) {
        success_message <- paste0("Loaded ", label, ".")
      }
      send_rds_load_state(
        "success",
        success_message,
        progress = progress$percent %||% 100,
        elapsed_label = rds_load_elapsed_label(progress)
      )
      showNotification("RDS path loaded.", type = "message")
    })

    source_summary_text <- reactive({
      if (identical(rds_load_state(), "running")) {
        label <- rds_load_path_label()
        if (!nzchar(label)) {
          label <- "RDS path"
        }
        return(paste("Loading", label, "..."))
      }

      data_source_summary(demo_data())
    })

    output$source_summary <- renderText({
      source_summary_text()
    })

    output$rds_schema_report <- renderUI({
      schema <- rds_schema()
      if (is.null(schema)) {
        return(tags$div(
          class = "rds-schema-report idle",
          "Validate RDS to preview assay, metadata, embeddings, proximity, dimensions, and cache size."
        ))
      }

      tags$pre(
        class = "rds-schema-report",
        format_user_rds_schema_report(schema)
      )
    })

    list(
      data = reactive(demo_data()),
      state = reactive(rds_load_state()),
      summary = source_summary_text
    )
  })
}
