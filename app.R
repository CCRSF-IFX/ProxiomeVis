resolve_app_dir <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) {
      NA_character_
    } else {
      as.character(frame$ofile)
    }
  }, character(1))
  frame_dirs <- dirname(normalizePath(frame_files[!is.na(frame_files)], mustWork = FALSE))

  candidates <- unique(c(
    getwd(),
    frame_dirs,
    file.path(getwd(), "shiny", "proxiome_demo")
  ))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "R", "data_adapter.R"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop("Could not resolve the Shiny app directory.", call. = FALSE)
}

APP_DIR <- resolve_app_dir()
APP_RENV_PROJECT <- APP_DIR

detect_proxiome_platform <- function(
  hostname = Sys.info()[["nodename"]] %||% "",
  env_platform = Sys.getenv("PROXIOME_PLATFORM", unset = "")
) {
  env_platform <- tolower(trimws(env_platform))
  if (nzchar(env_platform)) {
    return(env_platform)
  }

  shinyapps_markers <- c("SHINYAPPS_ACCOUNT", "SHINYAPPS_APP_NAME", "SHINYAPPS_APPLICATION_ID")
  if (any(nzchar(Sys.getenv(shinyapps_markers, unset = "")))) {
    return("shinyapps")
  }

  hostname <- tolower(trimws(as.character(hostname %||% "")))
  if (identical(hostname, "ncifcrf.gov") || grepl("\\.ncifcrf\\.gov$", hostname)) {
    return("ccrsf_hpc")
  }
  if (identical(hostname, "biowulf.nih.gov") || grepl("\\.biowulf\\.nih\\.gov$", hostname)) {
    return("biowulf_hpc")
  }

  "portable"
}

proxiome_platform_requires_shared_renv <- function(platform = APP_PLATFORM) {
  allow_system_libs <- tolower(Sys.getenv("PROXIOME_ALLOW_SYSTEM_LIBS", unset = "false")) %in% c("1", "true", "yes")
  platform %in% c("ccrsf_hpc", "biowulf_hpc") && !allow_system_libs
}

APP_PLATFORM <- detect_proxiome_platform()

required_rlang_version <- package_version("1.1.7")

needs_clean_renv_relaunch <- function() {
  if (!proxiome_platform_requires_shared_renv(APP_PLATFORM)) {
    return(FALSE)
  }
  if (identical(Sys.getenv("PROXIOME_REEXEC", unset = ""), "1")) {
    return(FALSE)
  }
  if (!"rlang" %in% loadedNamespaces()) {
    return(FALSE)
  }

  utils::packageVersion("rlang") < required_rlang_version
}

build_clean_renv_launch_expression <- function(app_dir, shiny_port, launch_browser) {
  app_dir_expr <- paste(capture.output(dput(normalizePath(app_dir, mustWork = TRUE))), collapse = "\n")
  port_expr <- paste(capture.output(dput(shiny_port)), collapse = "\n")
  launch_browser_expr <- paste(capture.output(dput(isTRUE(launch_browser))), collapse = "\n")

  paste(
    sprintf("app_dir <- %s", app_dir_expr),
    "Sys.setenv(PROXIOME_REEXEC = '1')",
    "Sys.setenv(RENV_PROJECT = app_dir)",
    "Sys.setenv(RENV_PATHS_LIBRARY = file.path(app_dir, 'renv', 'library'))",
    "Sys.setenv(RENV_CONFIG_CACHE_ENABLED = 'FALSE')",
    "Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE')",
    "Sys.setenv(RENV_DEPLOYMENT_STRICT = 'TRUE')",
    "source(file.path(app_dir, 'renv', 'activate.R'))",
    "required_packages <- c('renv', 'shiny', 'bslib', 'ggplot2', 'plotly', 'future', 'promises', 'rlang')",
    "library_root <- file.path(app_dir, 'renv', 'library')",
    "candidate_libraries <- unique(c(Sys.getenv('RENV_PATHS_LIBRARY', unset = ''), library_root, Sys.glob(file.path(library_root, 'R-*', '*')), Sys.glob(file.path(library_root, '*', 'R-*', '*'))))",
    "candidate_libraries <- candidate_libraries[nzchar(candidate_libraries) & dir.exists(candidate_libraries)]",
    "complete_libraries <- candidate_libraries[vapply(candidate_libraries, function(library_root) all(dir.exists(file.path(library_root, required_packages))), logical(1))]",
    "if (length(complete_libraries) == 0) stop('No complete app-local renv library was found before launching Shiny.', call. = FALSE)",
    ".libPaths(normalizePath(c(complete_libraries, .libPaths()), mustWork = FALSE))",
    "Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))",
    sprintf(".shiny_port <- %s", port_expr),
    sprintf("options(shiny.port = .shiny_port, shiny.launch.browser = %s)", launch_browser_expr),
    "shiny::runApp(app_dir)",
    sep = "; "
  )
}

relaunch_with_clean_renv_if_needed <- function(app_dir = APP_DIR) {
  if (!needs_clean_renv_relaunch()) {
    return(invisible(FALSE))
  }

  r_bin <- file.path(R.home("bin"), "R")
  if (!file.exists(r_bin)) {
    r_bin <- "R"
  }
  child_expr <- build_clean_renv_launch_expression(
    app_dir = app_dir,
    shiny_port = getOption("shiny.port"),
    launch_browser = getOption("shiny.launch.browser", FALSE)
  )
  env <- c(
    "PROXIOME_REEXEC=1",
    sprintf("RENV_PROJECT=%s", normalizePath(app_dir, mustWork = TRUE)),
    sprintf("RENV_PATHS_LIBRARY=%s", file.path(normalizePath(app_dir, mustWork = TRUE), "renv", "library")),
    "RENV_CONFIG_CACHE_ENABLED=FALSE",
    "RENV_CONFIG_AUTOLOADER_ENABLED=FALSE",
    "RENV_DEPLOYMENT_STRICT=TRUE"
  )

  message(
    sprintf(
      "Detected preloaded rlang %s from the launcher; relaunching with the app-local renv before loading app packages.",
      as.character(utils::packageVersion("rlang"))
    )
  )
  launch_file <- tempfile("proxiome-clean-renv-launch-", fileext = ".R")
  writeLines(child_expr, launch_file, useBytes = TRUE)
  status <- system2(r_bin, c("--no-save", "--slave", "-f", launch_file), env = env)
  if (is.null(status)) {
    status <- 0L
  }
  quit(save = "no", status = status, runLast = FALSE)
}

relaunch_with_clean_renv_if_needed(APP_DIR)

split_path_env <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(character(0))
  }
  unlist(strsplit(value, .Platform$path.sep, fixed = TRUE), use.names = FALSE)
}

is_disallowed_r_lib_path <- function(paths) {
  paths <- normalizePath(paths, mustWork = FALSE)
  grepl("/.pixi/", paths, fixed = TRUE)
}

normalize_existing_paths <- function(paths) {
  paths <- trimws(paths)
  paths <- paths[nzchar(paths)]
  paths <- paths[dir.exists(paths)]
  paths <- paths[!is_disallowed_r_lib_path(paths)]
  unique(normalizePath(paths, mustWork = TRUE))
}

renv_library_candidates <- function(project_root = APP_RENV_PROJECT) {
  project_root <- normalizePath(project_root, mustWork = FALSE)
  configured_library <- Sys.getenv("RENV_PATHS_LIBRARY", unset = "")
  project_library <- file.path(project_root, "renv", "library")
  versioned_libraries <- c(
    Sys.glob(file.path(project_library, "R-*", "*")),
    Sys.glob(file.path(project_library, "*", "R-*", "*"))
  )

  unique(c(
    split_path_env(configured_library),
    project_library,
    versioned_libraries
  ))
}

required_app_packages <- function(include_renv = FALSE) {
  packages <- c("shiny", "bslib", "ggplot2", "plotly", "future", "promises", "rlang")
  if (isTRUE(include_renv)) {
    packages <- c("renv", packages)
  }
  packages
}

library_directly_contains_package <- function(library_root, package) {
  dir.exists(file.path(library_root, package))
}

library_directly_contains_packages <- function(library_root, packages) {
  dir.exists(library_root) && all(vapply(packages, function(package) {
    library_directly_contains_package(library_root, package)
  }, logical(1)))
}

concrete_library_candidates_from_roots <- function(library_roots) {
  library_roots <- normalize_existing_paths(library_roots)
  nested_libraries <- unlist(lapply(library_roots, function(library_root) {
    c(
      Sys.glob(file.path(library_root, "R-*", "*")),
      Sys.glob(file.path(library_root, "*", "R-*", "*"))
    )
  }), use.names = FALSE)

  normalize_existing_paths(unique(c(library_roots, nested_libraries)))
}

current_r_minor_version <- function(r_version = getRversion()) {
  version <- as.character(r_version)
  paste(strsplit(version, ".", fixed = TRUE)[[1]][seq_len(2)], collapse = ".")
}

library_r_minor_version <- function(library_root) {
  matches <- regmatches(
    library_root,
    gregexpr("R-[0-9]+\\.[0-9]+", library_root, perl = TRUE)
  )[[1]]
  if (length(matches) == 0 || identical(matches, character(0))) {
    return(NA_character_)
  }

  sub("^R-", "", tail(matches, 1))
}

library_matches_current_r <- function(library_root, r_version = getRversion()) {
  library_minor <- library_r_minor_version(normalizePath(library_root, mustWork = FALSE))
  is.na(library_minor) || identical(library_minor, current_r_minor_version(r_version))
}

renv_concrete_library_candidates <- function(project_root = APP_RENV_PROJECT) {
  candidates <- concrete_library_candidates_from_roots(renv_library_candidates(project_root))
  candidates[vapply(candidates, library_matches_current_r, logical(1))]
}

restored_renv_library_paths <- function(
  project_root = APP_RENV_PROJECT,
  required_packages = required_app_packages(include_renv = TRUE)
) {
  library_roots <- renv_concrete_library_candidates(project_root)
  library_roots[vapply(library_roots, function(library_root) {
    library_directly_contains_packages(library_root, required_packages)
  }, logical(1))]
}

renv_library_is_restored <- function(library_roots) {
  library_roots <- concrete_library_candidates_from_roots(library_roots)
  required_packages <- required_app_packages(include_renv = TRUE)

  any(vapply(library_roots, function(library_root) {
    library_directly_contains_packages(library_root, required_packages)
  }, logical(1)))
}

required_startup_packages_are_available <- function(lib_paths = .libPaths()) {
  required_packages <- required_app_packages()
  all(vapply(required_packages, requireNamespace, logical(1), quietly = TRUE))
}

activate_project_renv <- function(project_root = APP_RENV_PROJECT, platform = APP_PLATFORM) {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  activate_path <- file.path(project_root, "renv", "activate.R")
  has_activate <- file.exists(activate_path)

  Sys.setenv(RENV_CONFIG_CACHE_ENABLED = "FALSE")
  if (!nzchar(Sys.getenv("RENV_PATHS_LIBRARY", unset = ""))) {
    Sys.setenv(RENV_PATHS_LIBRARY = file.path(project_root, "renv", "library"))
  }

  restored_library_paths <- if (has_activate) {
    restored_renv_library_paths(project_root)
  } else {
    character(0)
  }
  is_restored <- length(restored_library_paths) > 0

  if (!is_restored && proxiome_platform_requires_shared_renv(platform)) {
    stop(
      paste(
        "renv project library is not restored with the app packages.",
        sprintf("Project: %s", project_root),
        sprintf("Platform: %s", platform),
        "Maintainer action: restore the renv project library during deployment.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  if (!is_restored) {
    if (!required_startup_packages_are_available()) {
      stop(
        paste(
          "Required R packages are not available from the active R libraries.",
          sprintf("Project: %s", project_root),
          "Restore the project dependencies or deploy through a service that installs renv.lock dependencies.",
          sep = "\n"
        ),
        call. = FALSE
      )
    }
    .libPaths(normalize_existing_paths(.libPaths()))
    Sys.setenv(RENV_PROJECT = project_root)
    Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    return(invisible(.libPaths()))
  }

  Sys.setenv(RENV_PROJECT = project_root)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_root)
  source(activate_path, local = TRUE)
  .libPaths(normalize_existing_paths(c(restored_library_paths, .libPaths())))
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

  invisible(.libPaths())
}

resolve_app_r_libs <- function(app_dir = APP_DIR) {
  app_dir <- normalizePath(app_dir, mustWork = FALSE)
  candidates <- c(
    renv_library_candidates(app_dir),
    split_path_env(Sys.getenv("R_LIBS", unset = "")),
    split_path_env(Sys.getenv("R_LIBS_USER", unset = "")),
    split_path_env(Sys.getenv("R_LIBS_SITE", unset = "")),
    .libPaths()
  )

  libs <- normalize_existing_paths(candidates)
  libs[vapply(libs, library_matches_current_r, logical(1))]
}

configure_app_r_libs <- function(app_dir = APP_DIR) {
  app_dir <- normalizePath(app_dir, mustWork = FALSE)
  libs <- resolve_app_r_libs(app_dir)
  if (length(libs) > 0) {
    .libPaths(normalize_existing_paths(c(libs, .libPaths())))
  }

  Sys.setenv(RENV_PROJECT = app_dir)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))

  invisible(.libPaths())
}

APP_R_LIBS <- activate_project_renv(APP_RENV_PROJECT)

suppressPackageStartupMessages({
  library(bslib)
  library(ggplot2)
  library(plotly)
  library(shiny)
})

source(file.path(APP_DIR, "R", "data_adapter.R"))
source(file.path(APP_DIR, "R", "spatial_metrics.R"))
source(file.path(APP_DIR, "R", "plot_layout.R"))
source(file.path(APP_DIR, "R", "environment_module.R"), local = TRUE)

app_css <- function() {
  tags$style(
    "
    :root {
      --ink: #192124;
      --muted: #5a676d;
      --panel: #ffffff;
      --line: #d8e0df;
      --teal: #176d73;
      --coral: #c7503e;
      --amber: #c58a20;
      --mint: #dcebe7;
    }

    body {
      background: #f6f8f7;
      color: var(--ink);
      letter-spacing: 0;
    }

    .source-chip {
      max-width: 420px;
      padding: 6px 10px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fbfcfc;
      color: var(--muted);
      font-size: 0.82rem;
      word-break: break-word;
    }

    .navbar .source-chip {
      display: inline-block;
      margin-left: 12px;
    }

    .data-source-button {
      margin-left: 12px;
    }

    .data-source-popover-shell {
      width: min(360px, calc(100vw - 32px));
      max-width: min(360px, calc(100vw - 32px)) !important;
    }

    .data-source-popover-shell .popover-body {
      width: 100%;
      max-width: 100%;
      overflow: hidden;
    }

    .data-source-popover {
      width: 100%;
      max-width: 100%;
      min-width: 0;
    }

    .data-source-popover .form-group,
    .data-source-popover .form-label,
    .data-source-popover .progress {
      margin-bottom: 10px;
    }

    .data-source-popover .shiny-input-container {
      width: 100% !important;
      max-width: 100%;
    }

    .data-source-popover .input-group {
      width: 100%;
      max-width: 100%;
      min-width: 0;
      flex-wrap: nowrap;
    }

    .data-source-popover .input-group-btn,
    .data-source-popover .input-group-prepend {
      flex: 0 0 auto;
    }

    .data-source-popover .input-group > .form-control {
      flex: 1 1 auto;
      width: 1%;
      min-width: 0;
      overflow: hidden;
      text-overflow: ellipsis;
    }

    .rds-load-status {
      margin-top: 10px;
      padding: 8px 10px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fbfcfc;
      color: var(--muted);
      font-size: 0.82rem;
      line-height: 1.35;
    }

    .rds-load-status.running {
      border-color: #d8b35c;
      background: #fff8e6;
      color: #5f4600;
    }

    .rds-load-status.success {
      border-color: #9bc7bb;
      background: #edf7f3;
      color: #14565b;
    }

    .rds-load-status.error {
      border-color: #df9c91;
      background: #fff0ed;
      color: #8a2d20;
    }

    .rds-load-progress {
      display: none;
      margin-top: 10px;
    }

    .rds-load-progress.running,
    .rds-load-progress.success,
    .rds-load-progress.error {
      display: block;
    }

    .rds-load-progress-track {
      width: 100%;
      height: 8px;
      overflow: hidden;
      border-radius: 999px;
      background: #e8eeee;
    }

    .rds-load-progress-bar {
      width: 0%;
      height: 100%;
      border-radius: inherit;
      background: #176d73;
      transition: width 180ms ease-out;
    }

    .rds-load-progress.running .rds-load-progress-bar {
      background: #c58a1b;
    }

    .rds-load-progress.error .rds-load-progress-bar {
      background: #b94332;
    }

    .rds-load-progress.success .rds-load-progress-bar {
      background: #176d73;
    }

    .rds-load-elapsed {
      margin-top: 6px;
      color: var(--muted);
      font-size: 0.78rem;
      line-height: 1.25;
    }

    .bslib-sidebar-layout {
      padding: 18px 22px 20px;
    }

    .summary-box-row {
      margin-bottom: 14px;
    }

    .plot-pane {
      width: 100%;
      max-width: 100%;
      margin-inline: auto;
      min-height: 468px;
      padding-top: 12px;
      padding-bottom: 26px;
    }

    .plot-pane-compact {
      width: 100%;
      max-width: 820px;
      margin-inline: auto;
    }

    .plot-pane-standard {
      width: 100%;
      max-width: 960px;
      margin-inline: auto;
    }

    .plot-pane-wide {
      width: 100%;
      max-width: 1120px;
      margin-inline: auto;
    }

    .plot-pane-scroll {
      width: 100%;
      max-width: 100%;
      margin-inline: 0;
      overflow-x: auto;
    }

    .plot-pane > .html-widget,
    .plot-pane > .shiny-plot-output,
    .detail-grid > .html-widget,
    .detail-grid > .shiny-plot-output {
      width: 100% !important;
      max-width: 100%;
    }

    .coloc-heatmap-pane {
      overflow-x: auto;
      overflow-y: visible;
    }

    .coloc-heatmap-pane .html-widget,
    .coloc-heatmap-pane .shiny-plot-output {
      max-width: none;
    }

    .distribution-plot-pane {
      overflow-x: auto;
      overflow-y: visible;
    }

    .distribution-plot-shell {
      max-width: none;
    }

    .table-pane {
      margin-top: 24px;
      overflow-x: auto;
    }

    .detail-grid {
      display: grid;
      grid-template-columns: minmax(0, 1fr);
      gap: 16px;
      margin-top: 16px;
    }

    .differential-plot-row {
      display: grid;
      grid-template-columns: max-content max-content;
      gap: 22px;
      align-items: start;
      width: 100%;
      max-width: 100%;
      margin-inline: 0;
      overflow-x: auto;
      overflow-y: visible;
      padding-bottom: 8px;
    }

    .differential-plot-row .plot-pane {
      width: max-content;
      max-width: none;
      min-height: 0;
      margin-inline: 0;
      padding-bottom: 18px;
    }

    .differential-plot-row .plot-pane-compact,
    .differential-plot-row .plot-pane-standard {
      max-width: none;
    }

    .differential-plot-row .html-widget {
      width: auto !important;
      max-width: none;
      height: auto !important;
    }

    @media (max-width: 1100px) {
      .differential-plot-row {
        grid-template-columns: max-content;
      }
    }

    .tab-content,
    .tab-pane,
    .card,
    .card-body {
      overflow: visible;
    }

    .loading-state {
      display: flex;
      align-items: center;
      justify-content: center;
      min-height: 460px;
      color: var(--muted);
      border: 1px dashed var(--line);
      border-radius: 8px;
      background: #fbfcfc;
    }

    @media (max-width: 900px) {
      .source-chip {
        margin-top: 12px;
        max-width: none;
      }
    }
    "
  )
}

app_js <- function() {
  tags$script(HTML(
    "
    (function() {
      function setRdsLoadState(message) {
        message = message || {};
        var state = message.state || 'idle';
        var loading = state === 'running' || message.disabled === true;
        var button = document.getElementById('load_rds_path');
        var status = document.getElementById('rds_load_status');
        var progress = document.getElementById('rds_load_progress');
        var progressBar = document.getElementById('rds_load_progress_bar');
        var elapsed = document.getElementById('rds_load_elapsed');

        if (button) {
          button.disabled = loading;
          button.classList.toggle('disabled', loading);
          button.setAttribute('aria-disabled', loading ? 'true' : 'false');
          button.textContent = message.buttonLabel || (loading ? 'Loading...' : 'Load Data');
        }

        if (status) {
          status.className = 'rds-load-status ' + state;
          status.textContent = message.message || '';
        }

        if (progress) {
          progress.className = 'rds-load-progress ' + state;
          progress.setAttribute('aria-hidden', state === 'idle' ? 'true' : 'false');
        }

        if (progressBar) {
          var value = Number(message.progress || 0);
          if (!isFinite(value)) {
            value = 0;
          }
          value = Math.max(0, Math.min(100, value));
          progressBar.style.width = value + '%';
          progressBar.setAttribute('aria-valuenow', String(Math.round(value)));
        }

        if (elapsed) {
          elapsed.textContent = message.elapsedLabel || '';
        }
      }

      function registerRdsLoadHandler() {
        if (!window.Shiny || !Shiny.addCustomMessageHandler || window.proxiomeRdsLoadHandlerRegistered) {
          return;
        }
        Shiny.addCustomMessageHandler('proxiome-rds-load-state', setRdsLoadState);
        window.proxiomeRdsLoadHandlerRegistered = true;
      }

      registerRdsLoadHandler();
      document.addEventListener('shiny:connected', registerRdsLoadHandler);
      document.addEventListener('click', function(event) {
        var target = event.target;
        if (!target || !target.closest) {
          return;
        }

        var button = target.closest('#load_rds_path');
        if (!button || button.disabled) {
          return;
        }

        window.setTimeout(function() {
          var pathInput = document.getElementById('rds_server_path');
          var label = 'RDS file';
          if (pathInput && pathInput.value) {
            var parts = pathInput.value.split(/[\\\\/]/);
            label = parts[parts.length - 1] || label;
          }
          setRdsLoadState({
            state: 'running',
            disabled: true,
            buttonLabel: 'Loading...',
            message: 'Starting RDS load for ' + label + '.',
            progress: 3,
            elapsedLabel: 'Elapsed: 0s'
          });
        }, 0);
      }, true);
    })();
    "
  ))
}

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

qc_sidebar <- function() {
  sidebar(
    title = "QC controls",
    width = 300,
    accordion(
      open = c("Filters", "Cutoffs", "Display"),
      accordion_panel(
        "Filters",
        selectizeInput("qc_sample_filter", "Sample", choices = character(0), multiple = TRUE),
        selectInput(
          "qc_metadata_source",
          "Metadata",
          choices = c("Original cells" = "origin", "Filtered cells" = "filtered"),
          selected = "origin"
        )
      ),
      accordion_panel(
        "Cutoffs",
        numericInput("qc_n_umi_cutoff", "n_umi cutoff", value = 10000, min = 0, step = 500),
        numericInput("qc_isotype_cutoff", "Isotype fraction cutoff", value = 0.001, min = 0, max = 1, step = 0.0005)
      ),
      accordion_panel(
        "Display",
        selectInput(
          "qc_filter_y",
          "Filter count y-axis",
          choices = c("Number of cells" = "count", "Fraction of loaded cells" = "fraction_loaded"),
          selected = "count"
        ),
        checkboxInput("qc_filter_include_total", "Include TOTAL trajectory", value = TRUE),
        selectInput("qc_metric", "Distribution metric", choices = character(0))
      )
    )
  )
}

abundance_sidebar <- function() {
  sidebar(
    title = "Abundance controls",
    width = 300,
    conditionalPanel(
      condition = "input.abundance_mode == 'Observed'",
      accordion(
        open = c("Display", "Filters"),
        accordion_panel(
          "Display",
          selectInput("abundance_embedding", "Embedding", choices = character(0)),
          selectInput("abundance_color_by", "Color UMAP by", choices = character(0)),
          conditionalPanel(
            condition = "input.abundance_color_by == 'abundance'",
            selectInput("abundance_marker", "Marker", choices = character(0))
          ),
          selectInput("abundance_split_by", "Split UMAP by", choices = character(0)),
          sliderInput("abundance_point_size", "Dot size", min = 0.5, max = 5, value = 1.9, step = 0.1)
        ),
        accordion_panel(
          "Filters",
          selectizeInput("abundance_condition_filter", "Condition", choices = character(0), multiple = TRUE),
          selectizeInput("abundance_celltype_filter", "Cell type", choices = character(0), multiple = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.abundance_mode == 'Marker Distributions'",
      accordion(
        open = c("Display"),
        accordion_panel(
          "Display",
          selectInput("abundance_distribution_marker", "Marker", choices = character(0)),
          numericInput("abundance_distribution_columns", "Facet columns", value = 3, min = 1, max = 12, step = 1),
          numericInput("abundance_distribution_width", "Plot width (px)", value = 832, min = 420, max = 2600, step = 50),
          numericInput("abundance_distribution_height", "Plot height (px)", value = 678, min = 320, max = 2600, step = 50),
          checkboxInput("abundance_distribution_show_jitter", "Show jitter dots", value = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.abundance_mode == 'Differential'",
      accordion(
        open = c("Contrast", "Thresholds"),
        accordion_panel(
          "Contrast",
          selectInput("abundance_diff_group_a", "Group A", choices = character(0)),
          selectInput("abundance_diff_group_b", "Group B (reference)", choices = character(0)),
          selectizeInput("abundance_diff_celltype_filter", "Cell type", choices = character(0), multiple = TRUE),
          checkboxInput("abundance_diff_stratify_celltype", "Stratify by cell type", value = FALSE),
          actionButton("abundance_run_differential", "Run differential analysis", class = "btn-primary w-100")
        ),
        accordion_panel(
          "Thresholds",
          numericInput("abundance_diff_fdr", "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput("abundance_diff_effect", "Minimum effect", value = 0.25, min = 0, step = 0.05),
          numericInput("abundance_diff_min_cells", "Minimum cells per group", value = 3, min = 1, step = 1)
        ),
        accordion_panel(
          "Detail",
          selectInput("abundance_diff_marker", "Detail marker", choices = character(0))
        )
      )
    )
  )
}

clustering_sidebar <- function() {
  sidebar(
    title = "Clustering controls",
    width = 300,
    conditionalPanel(
      condition = "input.clustering_mode == 'Observed' || input.clustering_mode == 'Per Marker'",
      accordion(
        open = c("Display", "Filters"),
        accordion_panel(
          "Display",
          selectInput("clustering_marker", "Marker", choices = character(0))
        ),
        accordion_panel(
          "Filters",
          selectizeInput("clustering_condition_filter", "Condition", choices = character(0), multiple = TRUE),
          selectizeInput("clustering_celltype_filter", "Cell type", choices = character(0), multiple = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.clustering_mode == 'Summary Heatmap'",
      accordion(
        open = c("Display", "Filters"),
        accordion_panel(
          "Display",
          numericInput("clustering_heatmap_marker_count", "Top markers", value = 20, min = 2, max = 40, step = 1)
        ),
        accordion_panel(
          "Filters",
          selectizeInput("clustering_heatmap_condition_filter", "Condition", choices = character(0), multiple = TRUE),
          selectizeInput("clustering_heatmap_celltype_filter", "Cell type", choices = character(0), multiple = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.clustering_mode == 'Differential'",
      accordion(
        open = c("Contrast", "Thresholds"),
        accordion_panel(
          "Contrast",
          selectInput("clustering_diff_group_a", "Group A", choices = character(0)),
          selectInput("clustering_diff_group_b", "Group B (reference)", choices = character(0)),
          selectizeInput("clustering_diff_celltype_filter", "Cell type", choices = character(0), multiple = TRUE),
          checkboxInput("clustering_diff_stratify_celltype", "Stratify by cell type", value = FALSE),
          actionButton("clustering_run_differential", "Run differential analysis", class = "btn-primary w-100")
        ),
        accordion_panel(
          "Thresholds",
          numericInput("clustering_diff_fdr", "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput("clustering_diff_effect", "Minimum effect", value = 0.25, min = 0, step = 0.05),
          numericInput("clustering_diff_min_cells", "Minimum cells per group", value = 3, min = 1, step = 1)
        ),
        accordion_panel(
          "Detail",
          selectInput("clustering_diff_marker", "Detail marker", choices = character(0))
        )
      )
    )
  )
}

colocalization_sidebar <- function() {
  sidebar(
    title = "Colocalization controls",
    width = 300,
    conditionalPanel(
      condition = "input.colocalization_mode == 'Observed'",
      accordion(
        open = c("Display", "Filters"),
        accordion_panel(
          "Display",
          radioButtons(
            "colocalization_heatmap_display",
            "Plot style",
            choices = c("Interactive" = "interactive", "Original R plot" = "original"),
            selected = "interactive"
          ),
          selectInput(
            "colocalization_heatmap_preset",
            "Heatmap preset",
            choices = c("Custom" = "custom", "Report style" = "report"),
            selected = "custom"
          ),
          selectInput(
            "spatial_coloc_scope",
            "Heatmap scope",
            choices = c(
              "Condition summary" = "condition",
              "Sample summary" = "sample",
              "Cell type focus" = "celltype"
            ),
            selected = "condition"
          ),
          conditionalPanel(
            condition = "input.spatial_coloc_scope == 'celltype'",
            selectInput("spatial_celltype_focus", "Cell type focus", choices = character(0))
          ),
          selectInput(
            "spatial_marker_selection_mode",
            "Marker set",
            choices = c(
              "Variable detected markers" = "auto",
              "Selected markers" = "manual"
            ),
            selected = "auto"
          ),
          selectizeInput("colocalization_heatmap_markers", "Heatmap markers", choices = character(0), multiple = TRUE),
          conditionalPanel(
            condition = "input.spatial_marker_selection_mode == 'auto'",
            numericInput("spatial_top_marker_count", "Top markers", value = 20, min = 2, max = 40, step = 1),
            numericInput("spatial_min_pct_detected", "Minimum fraction detected", value = 0.25, min = 0, max = 1, step = 0.05),
            numericInput("spatial_min_log2_range", "Minimum log2 range", value = 0.2, min = 0, step = 0.05)
          ),
          selectInput("colocalization_reference_condition", "Reference condition", choices = character(0)),
          selectInput(
            "colocalization_clustering_method",
            "Marker ordering",
            choices = c("Ward D2" = "ward.D2", "Complete" = "complete", "Average" = "average", "Single" = "single"),
            selected = "ward.D2"
          ),
          numericInput("colocalization_legend_min", "Legend minimum", value = -1, step = 0.1),
          numericInput("colocalization_legend_max", "Legend maximum", value = 1, step = 0.1)
        ),
        accordion_panel(
          "Filters",
          selectizeInput("colocalization_condition_filter", "Condition", choices = character(0), multiple = TRUE),
          selectizeInput("colocalization_celltype_filter", "Cell type", choices = character(0), multiple = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.colocalization_mode == 'Differential'",
      accordion(
        open = c("Contrast", "Thresholds"),
        accordion_panel(
          "Contrast",
          selectInput("colocalization_diff_group_a", "Group A", choices = character(0)),
          selectInput("colocalization_diff_group_b", "Group B (reference)", choices = character(0)),
          selectizeInput("colocalization_diff_celltype_filter", "Cell type", choices = character(0), multiple = TRUE),
          checkboxInput("colocalization_diff_stratify_celltype", "Stratify by cell type", value = FALSE),
          actionButton("colocalization_run_differential", "Run differential analysis", class = "btn-primary w-100")
        ),
        accordion_panel(
          "Thresholds",
          numericInput("colocalization_diff_fdr", "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput("colocalization_diff_effect", "Minimum effect", value = 0.25, min = 0, step = 0.05),
          numericInput("colocalization_diff_min_cells", "Minimum cells per group", value = 3, min = 1, step = 1)
        ),
        accordion_panel(
          "Detail",
          selectInput("colocalization_diff_anchor_marker", "Anchor marker", choices = character(0)),
          selectInput("colocalization_diff_pair", "Detail pair", choices = character(0))
        )
      )
    )
  )
}

user_rds_path_loading_enabled <- function(platform = APP_PLATFORM) {
  disabled <- tolower(Sys.getenv("PROXIOME_DISABLE_USER_RDS", unset = "false")) %in% c("1", "true", "yes")
  platform %in% c("ccrsf_hpc", "biowulf_hpc", "portable") && !disabled
}

plot_pane <- function(..., size = c("standard", "compact", "wide", "scroll"), extra_class = NULL) {
  size <- match.arg(size)
  div(
    class = paste(c("plot-pane", paste0("plot-pane-", size), extra_class), collapse = " "),
    ...
  )
}

detail_pane <- function(..., size = c("standard", "compact", "wide")) {
  size <- match.arg(size)
  div(
    class = paste(c("detail-grid", paste0("plot-pane-", size)), collapse = " "),
    ...
  )
}

differential_plot_row <- function(volcano_output_id, detail_output_id) {
  div(
    class = "differential-plot-row",
    plot_pane(
      size = "compact",
      extra_class = "differential-plot-pane",
      plotlyOutput(volcano_output_id, width = "auto", height = "auto")
    ),
    plot_pane(
      size = "standard",
      extra_class = "differential-plot-pane",
      plotlyOutput(detail_output_id, width = "auto", height = "auto")
    )
  )
}

data_source_controls <- function(platform = APP_PLATFORM) {
  if (!user_rds_path_loading_enabled(platform)) {
    return(tagList())
  }

  popover(
    actionButton(
      "data_source_menu",
      "Data",
      class = "btn btn-outline-secondary btn-sm data-source-button"
    ),
    div(
      class = "data-source-popover",
      textInput(
        "rds_server_path",
        "RDS path",
        placeholder = "/path/to/sample.rds"
      ),
      actionButton(
        "load_rds_path",
        "Load Data",
        class = "btn btn-outline-primary btn-sm"
      ),
      actionButton(
        "use_demo_data",
        "Use demo data",
        class = "btn btn-outline-secondary btn-sm"
      ),
      div(
        id = "rds_load_status",
        class = "rds-load-status idle",
        "Enter an RDS path, then click Load Data."
      ),
      div(
        id = "rds_load_progress",
        class = "rds-load-progress idle",
        role = "progressbar",
        `aria-hidden` = "true",
        `aria-valuemin` = "0",
        `aria-valuemax` = "100",
        div(
          class = "rds-load-progress-track",
          div(id = "rds_load_progress_bar", class = "rds-load-progress-bar", `aria-valuenow` = "0")
        ),
        div(id = "rds_load_elapsed", class = "rds-load-elapsed")
      )
    ),
    title = "Data source",
    placement = "bottom",
    options = list(customClass = "data-source-popover-shell")
  )
}

ui <- page_navbar(
  title = "Pixelgen Proxiome Explorer",
  id = "readout_tab",
  fillable = c("QC", "Abundance", "Spatial Metrics", "Environment"),
  theme = bs_theme(
    version = 5,
    bg = "#f6f8f7",
    fg = "#192124",
    primary = "#176d73",
    secondary = "#c7503e",
    base_font = "system-ui"
  ),
  header = tagList(app_css(), app_js()),
  nav_panel(
    "QC",
    layout_sidebar(
      sidebar = qc_sidebar(),
      navset_card_underline(
        id = "qc_mode",
        title = "QC",
        full_screen = TRUE,
        nav_panel(
          "Filtering",
          uiOutput("qc_metric_row"),
          plot_pane(size = "compact", plotlyOutput("qc_filter_plot", height = proxiome_plot_height())),
          div(class = "table-pane", tableOutput("qc_filter_table"))
        ),
        nav_panel(
          "Cell Calling",
          plot_pane(size = "wide", plotlyOutput("qc_molecule_rank_plot", height = proxiome_plot_height()))
        ),
        nav_panel(
          "Distributions",
          plot_pane(size = "compact", plotlyOutput("qc_distribution_plot", height = proxiome_plot_height()))
        ),
        nav_panel(
          "Metadata",
          div(class = "table-pane", tableOutput("qc_origin_metadata_table"))
        )
      )
    )
  ),
  nav_panel(
    "Abundance",
    layout_sidebar(
      sidebar = abundance_sidebar(),
      navset_card_underline(
        id = "abundance_mode",
        title = "Abundance",
        full_screen = TRUE,
        nav_panel(
          "Observed",
          uiOutput("metric_row"),
          plot_pane(size = "standard", plotlyOutput("abundance_umap", height = proxiome_plot_height())),
          conditionalPanel(
            condition = "input.abundance_color_by == 'abundance'",
            div(class = "table-pane", tableOutput("abundance_table"))
          )
        ),
        nav_panel(
          "Marker Distributions",
          plot_pane(
            size = "scroll",
            extra_class = "distribution-plot-pane",
            uiOutput("abundance_marker_distribution_plot_ui")
          ),
          div(class = "table-pane", tableOutput("abundance_marker_distribution_table"))
        ),
        nav_panel(
          "Cell Annotation",
          plot_pane(size = "compact", plotlyOutput("abundance_celltype_composition_plot", height = proxiome_plot_height())),
          plot_pane(size = "wide", plotlyOutput("abundance_annotation_heatmap", height = proxiome_plot_height())),
          div(class = "table-pane", tableOutput("abundance_celltype_composition_table"))
        ),
        nav_panel(
          "Differential",
          uiOutput("abundance_diff_summary"),
          differential_plot_row("abundance_diff_volcano", "abundance_diff_detail"),
          div(class = "table-pane", tableOutput("abundance_diff_table"))
        )
      )
    )
  ),
  nav_panel(
    "Spatial Metrics",
    navset_tab(
      id = "spatial_metric_readout",
      nav_panel(
        "Clustering",
        layout_sidebar(
          sidebar = clustering_sidebar(),
          navset_card_underline(
            id = "clustering_mode",
            title = "Clustering",
            full_screen = TRUE,
            nav_panel(
              "Observed",
              plot_pane(size = "standard", plotlyOutput("clustering_plot", height = proxiome_plot_height())),
              div(class = "table-pane", tableOutput("clustering_table"))
            ),
            nav_panel(
              "Per Marker",
              plot_pane(size = "compact", plotlyOutput("clustering_per_marker_plot", height = proxiome_plot_height())),
              div(class = "table-pane", tableOutput("clustering_per_marker_table"))
            ),
            nav_panel(
              "Summary Heatmap",
              plot_pane(size = "scroll", plotlyOutput("clustering_summary_heatmap", height = proxiome_plot_height())),
              div(class = "table-pane", tableOutput("clustering_summary_heatmap_table"))
            ),
            nav_panel(
              "Differential",
              uiOutput("clustering_diff_summary"),
              differential_plot_row("clustering_diff_volcano", "clustering_diff_detail"),
              div(class = "table-pane", tableOutput("clustering_diff_table"))
            )
          )
        )
      ),
      nav_panel(
        "Colocalization",
        layout_sidebar(
          sidebar = colocalization_sidebar(),
          navset_card_underline(
            id = "colocalization_mode",
            title = "Colocalization",
            full_screen = TRUE,
            nav_panel(
              "Observed",
              plot_pane(
                size = "scroll",
                extra_class = "coloc-heatmap-pane",
                conditionalPanel(
                  condition = "input.colocalization_heatmap_display == 'interactive'",
                  plotlyOutput("colocalization_heatmap_interactive", height = coloc_heatmap_output_height())
                ),
                conditionalPanel(
                  condition = "input.colocalization_heatmap_display == 'original'",
                  plotOutput("colocalization_heatmap_original", height = coloc_heatmap_output_height())
                )
              ),
              div(class = "table-pane", tableOutput("colocalization_table"))
            ),
            nav_panel(
              "Differential",
              uiOutput("colocalization_diff_summary"),
              differential_plot_row("colocalization_diff_volcano", "colocalization_diff_detail"),
              div(class = "table-pane", tableOutput("colocalization_diff_table"))
            )
          )
        )
      )
    )
  ),
  nav_panel(
    "Environment",
    environment_module_ui("environment")
  ),
  nav_spacer(),
  nav_item(data_source_controls()),
  nav_item(tags$span(class = "source-chip", textOutput("source_summary", inline = TRUE)))
)

load_default_proxiome_data <- function() {
  data <- load_demo_proxiome_data(
    cache_path = default_demo_cache_path(APP_DIR)
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

make_differential_config <- function(
  group_a,
  group_b,
  celltype_filter,
  stratify_by_celltype = FALSE,
  min_cells = 3,
  fdr_cutoff = 0.05,
  effect_cutoff = 0.25,
  anchor_marker = NULL
) {
  list(
    group_a = group_a,
    group_b = group_b,
    celltype_filter = celltype_filter,
    stratify_by_celltype = isTRUE(stratify_by_celltype),
    min_cells = min_cells,
    fdr_cutoff = fdr_cutoff,
    effect_cutoff = effect_cutoff,
    anchor_marker = anchor_marker
  )
}

default_differential_config <- function(
  conditions,
  cell_types,
  anchor_marker = NULL,
  min_cells = 3,
  fdr_cutoff = 0.05,
  effect_cutoff = 0.25
) {
  make_differential_config(
    group_a = conditions[1],
    group_b = conditions[min(2, length(conditions))],
    celltype_filter = cell_types,
    stratify_by_celltype = FALSE,
    min_cells = min_cells,
    fdr_cutoff = fdr_cutoff,
    effect_cutoff = effect_cutoff,
    anchor_marker = anchor_marker
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

server <- function(input, output, session) {
  demo_data <- reactiveVal(NULL)
  abundance_diff_config <- reactiveVal(NULL)
  clustering_diff_config <- reactiveVal(NULL)
  colocalization_diff_config <- reactiveVal(NULL)
  user_rds_load_task <- if (user_rds_path_loading_enabled()) {
    create_user_rds_load_task(APP_DIR)
  } else {
    NULL
  }
  rds_load_state <- reactiveVal("idle")
  rds_load_message <- reactiveVal("Enter an RDS path, then click Load Data.")
  rds_load_path_label <- reactiveVal("")
  rds_load_progress_path <- reactiveVal(NULL)

  environment_module_server("environment", app_dir = APP_DIR)

  differential_config_from_inputs <- function(prefix, data, anchor_marker = NULL) {
    make_differential_config(
      group_a = input[[paste0(prefix, "_diff_group_a")]],
      group_b = input[[paste0(prefix, "_diff_group_b")]],
      celltype_filter = selected_or_all(
        input[[paste0(prefix, "_diff_celltype_filter")]],
        unique(data$metadata$celltype_manual)
      ),
      stratify_by_celltype = input[[paste0(prefix, "_diff_stratify_celltype")]],
      min_cells = numeric_input_value(input[[paste0(prefix, "_diff_min_cells")]], 3),
      fdr_cutoff = numeric_input_value(input[[paste0(prefix, "_diff_fdr")]], 0.05),
      effect_cutoff = numeric_input_value(input[[paste0(prefix, "_diff_effect")]], 0.25),
      anchor_marker = anchor_marker
    )
  }

  apply_differential_config <- function(prefix, target, anchor_marker = NULL) {
    data <- demo_data()
    req(data)
    target(differential_config_from_inputs(prefix, data, anchor_marker = anchor_marker))
  }

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
      data <- load_default_proxiome_data()
      incProgress(0.8)
      demo_data(data)
    })

    if (isTRUE(show_loaded_notification)) {
      showNotification("Demo data loaded.", type = "message")
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
  }, ignoreInit = TRUE)

  observeEvent(input$load_rds_path, {
    req(user_rds_path_loading_enabled())
    req(!is.null(user_rds_load_task))

    tryCatch(
      {
        validate_rds_file_path(input$rds_server_path)
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

  observe({
    data <- demo_data()
    req(data)

    embedding_choices <- available_embedding_choices(data$metadata)
    conditions <- sort(unique(data$metadata$condition))
    cell_types <- sort(unique(data$metadata$celltype_manual))
    default_group_a <- conditions[1]
    default_group_b <- conditions[min(2, length(conditions))]
    colocalization_pairs <- sort(unique(data$colocalization$marker_pair))
    colocalization_markers <- available_colocalization_marker_choices(data$colocalization)
    default_heatmap_markers <- head(colocalization_markers, min(20L, length(colocalization_markers)))
    default_colocalization_reference <- if ("CD3CD28" %in% conditions) "CD3CD28" else conditions[1]
    qc_samples <- qc_sample_choices(data$qc)
    qc_metric_choices <- available_qc_distribution_choices(data$qc$origin_metadata)

    updateSelectizeInput(session, "qc_sample_filter", choices = qc_samples, selected = qc_samples)
    updateSelectInput(session, "qc_metric", choices = qc_metric_choices, selected = qc_metric_choices[1])
    updateSelectInput(session, "abundance_embedding", choices = embedding_choices, selected = embedding_choices[1])
    updateSelectInput(session, "abundance_color_by", choices = available_abundance_color_choices(data$metadata), selected = "abundance")
    updateSelectInput(session, "abundance_marker", choices = data$marker_options, selected = data$marker_options[1])
    updateSelectInput(session, "abundance_distribution_marker", choices = data$marker_options, selected = data$marker_options[1])
    updateSelectInput(session, "abundance_split_by", choices = available_abundance_split_choices(data$metadata), selected = "")
    updateSelectInput(session, "clustering_marker", choices = data$marker_options, selected = data$marker_options[1])
    updateSelectizeInput(session, "colocalization_heatmap_markers", choices = colocalization_markers, selected = default_heatmap_markers)
    updateSelectInput(session, "spatial_celltype_focus", choices = cell_types, selected = cell_types[1])
    updateSelectInput(
      session,
      "colocalization_reference_condition",
      choices = conditions,
      selected = default_colocalization_reference
    )

    updateSelectizeInput(session, "abundance_condition_filter", choices = conditions, selected = conditions)
    updateSelectizeInput(session, "abundance_celltype_filter", choices = cell_types, selected = cell_types)
    updateSelectizeInput(session, "clustering_condition_filter", choices = conditions, selected = conditions)
    updateSelectizeInput(session, "clustering_celltype_filter", choices = cell_types, selected = cell_types)
    updateSelectizeInput(session, "clustering_heatmap_condition_filter", choices = conditions, selected = conditions)
    updateSelectizeInput(session, "clustering_heatmap_celltype_filter", choices = cell_types, selected = cell_types)
    updateSelectizeInput(session, "colocalization_condition_filter", choices = conditions, selected = conditions)
    updateSelectizeInput(session, "colocalization_celltype_filter", choices = cell_types, selected = cell_types)

    updateSelectInput(session, "abundance_diff_group_a", choices = conditions, selected = default_group_a)
    updateSelectInput(session, "abundance_diff_group_b", choices = conditions, selected = default_group_b)
    updateSelectizeInput(session, "abundance_diff_celltype_filter", choices = cell_types, selected = cell_types)
    updateSelectInput(session, "abundance_diff_marker", choices = data$marker_options, selected = data$marker_options[1])

    updateSelectInput(session, "clustering_diff_group_a", choices = conditions, selected = default_group_a)
    updateSelectInput(session, "clustering_diff_group_b", choices = conditions, selected = default_group_b)
    updateSelectizeInput(session, "clustering_diff_celltype_filter", choices = cell_types, selected = cell_types)
    updateSelectInput(session, "clustering_diff_marker", choices = data$marker_options, selected = data$marker_options[1])

    updateSelectInput(session, "colocalization_diff_group_a", choices = conditions, selected = default_group_a)
    updateSelectInput(session, "colocalization_diff_group_b", choices = conditions, selected = default_group_b)
    updateSelectizeInput(session, "colocalization_diff_celltype_filter", choices = cell_types, selected = cell_types)
    updateSelectInput(session, "colocalization_diff_anchor_marker", choices = data$marker_options, selected = data$marker_options[1])
    updateSelectInput(session, "colocalization_diff_pair", choices = colocalization_pairs, selected = colocalization_pairs[1])

    abundance_diff_config(default_differential_config(conditions, cell_types))
    clustering_diff_config(default_differential_config(conditions, cell_types))
    colocalization_diff_config(default_differential_config(
      conditions,
      cell_types,
      anchor_marker = data$marker_options[1]
    ))
  })

  observeEvent(input$abundance_run_differential, {
    apply_differential_config("abundance", abundance_diff_config)
  }, ignoreInit = TRUE)

  observeEvent(input$clustering_run_differential, {
    apply_differential_config("clustering", clustering_diff_config)
  }, ignoreInit = TRUE)

  observeEvent(input$colocalization_run_differential, {
    apply_differential_config(
      "colocalization",
      colocalization_diff_config,
      anchor_marker = input$colocalization_diff_anchor_marker
    )
  }, ignoreInit = TRUE)

  observeEvent(input$colocalization_heatmap_preset, {
    if (!identical(input$colocalization_heatmap_preset, "report")) {
      return(invisible(NULL))
    }

    updateSelectInput(session, "spatial_marker_selection_mode", selected = "auto")
    updateNumericInput(session, "spatial_top_marker_count", value = 40)
    updateNumericInput(session, "spatial_min_pct_detected", value = 0.25)
    updateNumericInput(session, "spatial_min_log2_range", value = 0.2)
    updateNumericInput(session, "colocalization_legend_min", value = -0.75)
    updateNumericInput(session, "colocalization_legend_max", value = 0.75)
    updateSelectInput(session, "colocalization_clustering_method", selected = "ward.D2")
  }, ignoreInit = TRUE)

  observe({
    data <- demo_data()
    config <- colocalization_diff_config()
    req(data, config$anchor_marker)

    choices <- sort(unique(data$colocalization$marker_pair[
      data$colocalization$marker_1 == config$anchor_marker |
        data$colocalization$marker_2 == config$anchor_marker
    ]))

    if (length(choices) == 0) {
      choices <- sort(unique(data$colocalization$marker_pair))
    }

    updateSelectInput(session, "colocalization_diff_pair", choices = choices, selected = choices[1])
  })

  output$source_summary <- renderText({
    if (identical(rds_load_state(), "running")) {
      label <- rds_load_path_label()
      if (!nzchar(label)) {
        label <- "RDS path"
      }
      return(paste("Loading", label, "..."))
    }

    data_source_summary(demo_data())
  })

  filtered_metadata_for <- function(condition_filter, celltype_filter) {
    data <- demo_data()
    req(data)

    metadata <- data$metadata
    conditions <- selected_or_all(condition_filter, unique(metadata$condition))
    cell_types <- selected_or_all(celltype_filter, unique(metadata$celltype_manual))

    metadata[
      metadata$condition %in% conditions &
        metadata$celltype_manual %in% cell_types,
      ,
      drop = FALSE
    ]
  }

  abundance_metadata <- reactive({
    filtered_metadata_for(input$abundance_condition_filter, input$abundance_celltype_filter)
  })

  clustering_metadata <- reactive({
    filtered_metadata_for(input$clustering_condition_filter, input$clustering_celltype_filter)
  })

  colocalization_metadata <- reactive({
    filtered_metadata_for(input$colocalization_condition_filter, input$colocalization_celltype_filter)
  })

  abundance_readout_with_metadata <- reactive({
    data <- demo_data()
    req(data)

    merge(
      data$abundance,
      data$metadata[, intersect(c("component", "condition", "celltype_manual"), names(data$metadata)), drop = FALSE],
      by = "component",
      all.x = FALSE,
      sort = FALSE
    )
  })

  abundance_diff_results <- reactive({
    data <- demo_data()
    config <- abundance_diff_config()
    req(data, config, config$group_a, config$group_b)

    calculate_differential_readout(
      abundance_readout_with_metadata(),
      feature_cols = "marker",
      value_col = "abundance",
      group_a = config$group_a,
      group_b = config$group_b,
      celltype_filter = config$celltype_filter,
      stratify_by_celltype = config$stratify_by_celltype,
      min_cells = config$min_cells,
      fdr_cutoff = config$fdr_cutoff
    )
  })

  clustering_diff_results <- reactive({
    data <- demo_data()
    config <- clustering_diff_config()
    req(data, config, config$group_a, config$group_b)

    calculate_differential_readout(
      data$clustering,
      feature_cols = "marker",
      value_col = "log2_ratio",
      group_a = config$group_a,
      group_b = config$group_b,
      celltype_filter = config$celltype_filter,
      stratify_by_celltype = config$stratify_by_celltype,
      min_cells = config$min_cells,
      fdr_cutoff = config$fdr_cutoff
    )
  })

  colocalization_diff_results <- reactive({
    data <- demo_data()
    config <- colocalization_diff_config()
    req(data, config, config$group_a, config$group_b)

    calculate_differential_readout(
      data$colocalization,
      feature_cols = c("marker_pair", "marker_1", "marker_2"),
      value_col = "log2_ratio",
      group_a = config$group_a,
      group_b = config$group_b,
      celltype_filter = config$celltype_filter,
      stratify_by_celltype = config$stratify_by_celltype,
      min_cells = config$min_cells,
      fdr_cutoff = config$fdr_cutoff
    )
  })

  colocalization_diff_anchor_results <- reactive({
    result <- colocalization_diff_results()
    config <- colocalization_diff_config()
    req(config$anchor_marker)

    result[
      result$marker_1 == config$anchor_marker |
        result$marker_2 == config$anchor_marker,
      ,
      drop = FALSE
    ]
  })

  qc_metadata <- reactive({
    data <- demo_data()
    req(data)

    metadata <- selected_qc_metadata(data$qc, input$qc_metadata_source)
    filter_qc_metadata_by_sample(metadata, input$qc_sample_filter)
  })

  qc_origin_metadata <- reactive({
    data <- demo_data()
    req(data)

    filter_qc_metadata_by_sample(data$qc$origin_metadata, input$qc_sample_filter)
  })

  qc_filtered_metadata <- reactive({
    data <- demo_data()
    req(data)

    filter_qc_metadata_by_sample(data$qc$filtered_metadata, input$qc_sample_filter)
  })

  qc_filter_counts <- reactive({
    data <- demo_data()
    req(data)

    qc_filter_counts_for_samples(
      data$qc$filter_counts,
      input$qc_sample_filter,
      include_total = isTRUE(input$qc_filter_include_total)
    )
  })

  output$qc_metric_row <- renderUI({
    origin_metadata <- qc_origin_metadata()
    filtered_metadata <- qc_filtered_metadata()

    loaded <- nrow(origin_metadata)
    final <- nrow(filtered_metadata)
    retained <- if (loaded > 0) final / loaded else NA_real_
    samples <- if ("sample" %in% names(origin_metadata)) length(unique(origin_metadata$sample)) else 1L

    metric_row(
      metric_tile("Loaded Cells", format(loaded, big.mark = ",")),
      metric_tile("Final Cells", format(final, big.mark = ",")),
      metric_tile("Retained", format_percent(retained)),
      metric_tile("Samples", format(samples, big.mark = ","))
    )
  })

  output$qc_filter_plot <- renderPlotly({
    counts <- qc_filter_counts()
    validate(need(nrow(counts) > 0, "No QC filter counts are available."))

    p <- plot_filter_cell_counts(
      counts,
      include_total = isTRUE(input$qc_filter_include_total),
      y = input$qc_filter_y %||% "count"
    )

    ggplotly(p, tooltip = "text") |>
      apply_proxiome_plot_frame()
  })

  output$qc_filter_table <- renderTable({
    format_qc_filter_table(qc_filter_counts())
  }, striped = TRUE, bordered = FALSE, width = "100%")

  output$qc_molecule_rank_plot <- renderPlotly({
    metadata <- qc_origin_metadata()
    validate(need("n_umi" %in% names(metadata), "The original metadata does not include n_umi."))
    validate(need(any(is.finite(metadata$n_umi) & metadata$n_umi > 0), "No positive n_umi values are available."))

    qc_molecule_rank_plotly(
      metadata,
      cutoff = numeric_input_value(input$qc_n_umi_cutoff, 10000)
    )
  })

  output$qc_distribution_plot <- renderPlotly({
    metadata <- qc_metadata()
    req(input$qc_metric)
    validate(need(input$qc_metric %in% names(metadata), "Selected QC metric is not available."))

    qc_distribution_plotly(
      metadata,
      metric = input$qc_metric,
      isotype_cutoff = numeric_input_value(input$qc_isotype_cutoff, 0.001)
    )
  })

  output$qc_origin_metadata_table <- renderTable({
    head(qc_origin_metadata(), 200)
  }, striped = TRUE, bordered = FALSE, width = "100%")

  output$metric_row <- renderUI({
    data <- demo_data()
    req(data)
    metadata <- abundance_metadata()

    metric_row(
      metric_tile("Cells", format(nrow(metadata), big.mark = ",")),
      metric_tile("Markers", format(length(data$marker_options), big.mark = ",")),
      metric_tile("Conditions", format(length(unique(metadata$condition)), big.mark = ",")),
      metric_tile("Cell Types", format(length(unique(metadata$celltype_manual)), big.mark = ","))
    )
  })

  abundance_points <- reactive({
    data <- demo_data()
    req(data, input$abundance_embedding)

    metadata <- abundance_metadata()
    color_by <- selected_abundance_color_by(input$abundance_color_by, metadata)

    if (identical(color_by, "abundance")) {
      req(input$abundance_marker)
      abundance <- data$abundance[data$abundance$marker == input$abundance_marker, , drop = FALSE]
      plot_data <- merge(metadata, abundance, by = "component", all.x = FALSE, sort = FALSE)
    } else {
      plot_data <- metadata
    }

    plot_data$abundance_color_by <- color_by
    add_embedding_columns(plot_data, input$abundance_embedding)
  })

  output$abundance_umap <- renderPlotly({
    plot_data <- abundance_points()
    validate(need(nrow(plot_data) > 0, "No cells match the selected filters."))

    color_by <- unique(plot_data$abundance_color_by)
    color_by <- color_by[1] %||% "abundance"
    point_size <- numeric_input_value(input$abundance_point_size, 1.9)

    if (identical(color_by, "abundance")) {
      plot_data$hover <- paste0(
        "Cell: ", plot_data$component,
        "<br>Condition: ", plot_data$condition,
        hover_sample_text(plot_data),
        "<br>Cell type: ", plot_data$celltype_manual,
        "<br>", input$abundance_marker, " abundance: ", round(plot_data$abundance, 3)
      )

      p <- ggplot(plot_data, aes(embedding_x, embedding_y, color = abundance, text = hover)) +
        geom_point(size = point_size, alpha = 0.82) +
        scale_color_gradientn(
          colors = c("#edf7f4", "#78aeb2", "#f0b45b", "#c7503e"),
          name = input$abundance_marker
        ) +
        labs(x = "Embedding 1", y = "Embedding 2") +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank())

      colorbar_title <- paste(input$abundance_marker, "abundance")
    } else {
      color_label <- abundance_color_label(color_by)
      plot_data$color_group <- as.factor(plot_data[[color_by]])
      plot_data$hover <- paste0(
        "Cell: ", plot_data$component,
        "<br>Condition: ", plot_data$condition,
        hover_sample_text(plot_data),
        "<br>Cell type: ", plot_data$celltype_manual,
        "<br>", color_label, ": ", plot_data$color_group
      )

      p <- ggplot(plot_data, aes(embedding_x, embedding_y, color = color_group, text = hover)) +
        geom_point(size = point_size, alpha = 0.82) +
        labs(x = "Embedding 1", y = "Embedding 2", color = color_label) +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank())

      colorbar_title <- NULL
    }

    split_col <- selected_split_column(input$abundance_split_by, plot_data)
    if (!is.null(split_col)) {
      plot_data$split_group <- plot_data[[split_col]]
      p <- p %+% plot_data +
        facet_wrap(~split_group)
    }

    ggplotly(p, tooltip = "text") |>
      apply_proxiome_plot_frame(colorbar_title = colorbar_title)
  })

  output$abundance_table <- renderTable({
    data <- demo_data()
    req(data, input$abundance_marker)

    summary <- data$abundance_summary[
      data$abundance_summary$marker == input$abundance_marker &
        data$abundance_summary$condition %in% selected_or_all(input$abundance_condition_filter, unique(data$metadata$condition)) &
        data$abundance_summary$celltype_manual %in% selected_or_all(input$abundance_celltype_filter, unique(data$metadata$celltype_manual)),
      ,
      drop = FALSE
    ]
    format_summary_table(summary, value_label = "mean_abundance")
  }, striped = TRUE, bordered = FALSE, width = "100%")

  abundance_distribution_data <- reactive({
    data <- demo_data()
    req(data, input$abundance_distribution_marker)

    plot_data <- data$abundance[data$abundance$marker == input$abundance_distribution_marker, , drop = FALSE]
    merge(
      plot_data,
      data$metadata[, intersect(c("component", "sample_alias", "sample", "condition", "celltype_manual"), names(data$metadata)), drop = FALSE],
      by = "component",
      all.x = FALSE,
      sort = FALSE
    )
  })

  update_abundance_distribution_size_controls <- function(plot_data, facet_cols = NULL) {
    dimensions <- abundance_distribution_widget_dimensions(plot_data, facet_cols = facet_cols)
    updateNumericInput(session, "abundance_distribution_columns", value = dimensions$facet_cols)
    updateNumericInput(session, "abundance_distribution_width", value = dimensions$width)
    updateNumericInput(session, "abundance_distribution_height", value = dimensions$height)
  }

  observeEvent(abundance_distribution_data(), {
    update_abundance_distribution_size_controls(abundance_distribution_data())
  }, ignoreInit = FALSE)

  observeEvent(input$abundance_distribution_columns, {
    update_abundance_distribution_size_controls(
      abundance_distribution_data(),
      facet_cols = input$abundance_distribution_columns
    )
  }, ignoreInit = TRUE)

  abundance_distribution_dimensions <- reactive({
    abundance_distribution_widget_dimensions(
      abundance_distribution_data(),
      facet_cols = input$abundance_distribution_columns,
      width_px = input$abundance_distribution_width,
      height_px = input$abundance_distribution_height
    )
  })

  output$abundance_marker_distribution_plot_ui <- renderUI({
    dimensions <- abundance_distribution_dimensions()
    div(
      class = "distribution-plot-shell",
      style = paste0("width:", dimensions$width, "px;height:", dimensions$height, "px;"),
      plotlyOutput("abundance_marker_distribution_plot", width = "100%", height = "100%")
    )
  })

  output$abundance_marker_distribution_plot <- renderPlotly({
    plot_data <- abundance_distribution_data()
    validate(need(nrow(plot_data) > 0, "No abundance rows are available for the selected marker."))

    dimensions <- abundance_distribution_dimensions()
    ggplotly(
      plot_abundance_marker_distribution(
        plot_data,
        input$abundance_distribution_marker,
        facet_cols = dimensions$facet_cols,
        show_jitter = !identical(input$abundance_distribution_show_jitter, FALSE)
      ),
      tooltip = c("x", "y", "fill"),
      width = dimensions$width,
      height = dimensions$height
    ) |>
      apply_proxiome_plot_frame()
  })

  output$abundance_marker_distribution_table <- renderTable({
    plot_data <- abundance_distribution_data()
    validate(need(nrow(plot_data) > 0, "No abundance rows are available for the selected marker."))

    summary <- aggregate_numeric_readout(
      plot_data,
      group_cols = c("marker", "sample_alias", "condition", "celltype_manual"),
      value_col = "abundance"
    )
    format_summary_table(summary, value_label = "mean_abundance")
  }, striped = TRUE, bordered = FALSE, width = "100%")

  output$abundance_celltype_composition_plot <- renderPlotly({
    data <- demo_data()
    req(data)
    plot_data <- celltype_composition_data(data$metadata)
    validate(need(nrow(plot_data) > 0, "No cell annotation rows are available."))

    ggplotly(plot_celltype_composition(plot_data), tooltip = "text") |>
      apply_proxiome_plot_frame()
  })

  output$abundance_annotation_heatmap <- renderPlotly({
    data <- demo_data()
    req(data)
    plot_data <- annotation_heatmap_data(data$abundance, data$metadata)
    validate(need(nrow(plot_data) > 0, "No abundance rows are available for annotation heatmap."))

    ggplotly(plot_annotation_heatmap(plot_data), tooltip = "text") |>
      apply_proxiome_plot_frame(colorbar_title = "Median abundance")
  })

  output$abundance_celltype_composition_table <- renderTable({
    data <- demo_data()
    req(data)
    composition <- celltype_composition_data(data$metadata)
    validate(need(nrow(composition) > 0, "No cell annotation rows are available."))
    format_numeric_table(composition)
  }, striped = TRUE, bordered = FALSE, width = "100%")

  output$abundance_diff_summary <- renderUI({
    config <- abundance_diff_config()
    req(config)
    result <- abundance_diff_results()
    differential_summary_row(
      result,
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff
    )
  })

  output$abundance_diff_volcano <- renderPlotly({
    config <- abundance_diff_config()
    req(config)
    result <- abundance_diff_results()
    validate(need(nrow(result) > 0, "Choose two different groups with enough abundance data."))

    x_label <- paste("Abundance effect:", config$group_a, "minus", config$group_b, "(reference)")
    differential_volcano_plot(
      result,
      label_col = "marker",
      x_label = x_label,
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff,
      source = "abundance_diff",
      dimensions = differential_volcano_dimensions(x_label)
    )
  })

  observeEvent(plotly::event_data("plotly_click", source = "abundance_diff"), {
    event <- plotly::event_data("plotly_click", source = "abundance_diff")
    if (!is.null(event$key) && nzchar(event$key)) {
      updateSelectInput(session, "abundance_diff_marker", selected = event$key)
    }
  })

  output$abundance_diff_detail <- renderPlotly({
    data <- demo_data()
    config <- abundance_diff_config()
    req(data, config, input$abundance_diff_marker, config$group_a, config$group_b)

    plot_data <- abundance_readout_with_metadata()
    plot_data <- plot_data[
      plot_data$marker == input$abundance_diff_marker &
        plot_data$condition %in% c(config$group_a, config$group_b) &
        plot_data$celltype_manual %in% config$celltype_filter,
      ,
      drop = FALSE
    ]
    validate(need(nrow(plot_data) > 0, "No abundance values are available for the selected marker and contrast."))

    plot_data$hover <- paste0(
      "Cell: ", plot_data$component,
      "<br>Condition: ", plot_data$condition,
      "<br>Cell type: ", plot_data$celltype_manual,
      "<br>Abundance: ", round(plot_data$abundance, 3)
    )

    y_label <- paste(input$abundance_diff_marker, "abundance")
    p <- ggplot(plot_data, aes(condition, abundance, color = condition, text = hover)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.18, linewidth = 0.5) +
      geom_jitter(width = 0.18, height = 0, alpha = 0.5, size = 1.4) +
      labs(x = NULL, y = y_label) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(), legend.position = "none")

    if (isTRUE(config$stratify_by_celltype)) {
      p <- p + facet_wrap(~celltype_manual, scales = "free_y")
    }

    dimensions <- differential_detail_dimensions(
      plot_data,
      stratify_by_celltype = isTRUE(config$stratify_by_celltype),
      y_label = y_label
    )
    ggplotly(p, tooltip = "text", width = dimensions$width, height = dimensions$height) |>
      apply_proxiome_plot_frame(dimensions = dimensions)
  })

  output$abundance_diff_table <- renderTable({
    config <- abundance_diff_config()
    req(config)
    result <- filter_differential_hits(
      abundance_diff_results(),
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff
    )
    validate(need(nrow(result) > 0, "No abundance markers pass the selected differential thresholds."))

    format_differential_table(result, effect_label = "abundance_effect_vs_reference")
  }, striped = TRUE, bordered = FALSE, width = "100%")

  clustering_points <- reactive({
    data <- demo_data()
    req(data, input$clustering_marker)

    metadata <- clustering_metadata()
    clustering <- data$clustering[data$clustering$marker == input$clustering_marker, , drop = FALSE]
    clustering <- clustering[clustering$component %in% metadata$component, , drop = FALSE]
    abundance <- data$abundance[data$abundance$marker == input$clustering_marker, , drop = FALSE]
    merge(clustering, abundance, by = c("component", "marker"), all.x = TRUE, sort = FALSE)
  })

  output$clustering_plot <- renderPlotly({
    plot_data <- clustering_points()
    validate(need(nrow(plot_data) > 0, "No self-clustering scores are available for the selected marker and filters."))

    plot_data$hover <- paste0(
      "Cell: ", plot_data$component,
      "<br>Condition: ", plot_data$condition,
      "<br>Cell type: ", plot_data$celltype_manual,
      "<br>Abundance: ", round(plot_data$abundance, 3),
      "<br>Self-clustering log2 ratio: ", round(plot_data$log2_ratio, 3)
    )

    p <- ggplot(plot_data, aes(abundance, log2_ratio, color = condition, text = hover)) +
      geom_hline(yintercept = 0, color = "#8a9699", linewidth = 0.6) +
      geom_point(size = 2, alpha = 0.76) +
      labs(x = paste(input$clustering_marker, "abundance"), y = "Self-clustering log2 ratio") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())

    ggplotly(p, tooltip = "text") |>
      apply_proxiome_plot_frame()
  })

  output$clustering_table <- renderTable({
    plot_data <- clustering_points()
    validate(need(nrow(plot_data) > 0, "No clustering rows to summarize."))

    summary <- aggregate_numeric_readout(
      plot_data,
      group_cols = c("marker", "condition", "celltype_manual"),
      value_col = "log2_ratio"
    )
    format_summary_table(summary, value_label = "mean_log2_ratio")
  }, striped = TRUE, bordered = FALSE, width = "100%")

  output$clustering_per_marker_plot <- renderPlotly({
    plot_data <- clustering_points()
    validate(need(nrow(plot_data) > 0, "No self-clustering scores are available for the selected marker and filters."))

    ggplotly(plot_clustering_per_marker(plot_data, input$clustering_marker), tooltip = c("x", "y", "fill")) |>
      apply_proxiome_plot_frame()
  })

  output$clustering_per_marker_table <- renderTable({
    plot_data <- clustering_points()
    validate(need(nrow(plot_data) > 0, "No clustering rows to summarize."))

    summary <- aggregate_numeric_readout(
      plot_data,
      group_cols = c("marker", "sample_alias", "condition", "celltype_manual"),
      value_col = "log2_ratio"
    )
    format_summary_table(summary, value_label = "mean_log2_ratio")
  }, striped = TRUE, bordered = FALSE, width = "100%")

  clustering_heatmap_summary <- reactive({
    data <- demo_data()
    req(data)

    clustering <- data$clustering
    conditions <- selected_or_all(input$clustering_heatmap_condition_filter, unique(clustering$condition))
    cell_types <- selected_or_all(input$clustering_heatmap_celltype_filter, unique(clustering$celltype_manual))
    clustering <- clustering[
      clustering$condition %in% conditions &
        clustering$celltype_manual %in% cell_types,
      ,
      drop = FALSE
    ]

    full_summary <- summarize_clustering_heatmap(
      clustering,
      selected_markers = unique(clustering$marker)
    )
    selected_markers <- select_clustering_heatmap_markers(
      full_summary,
      n_markers = numeric_input_value(input$clustering_heatmap_marker_count, 20)
    )
    summarize_clustering_heatmap(clustering, selected_markers = selected_markers)
  })

  output$clustering_summary_heatmap <- renderPlotly({
    summary <- clustering_heatmap_summary()
    validate(need(nrow(summary) > 0, "No clustering rows are available for the selected heatmap filters."))

    ggplotly(plot_clustering_summary_heatmap(summary), tooltip = "text") |>
      apply_proxiome_plot_frame(colorbar_title = "Mean log2 ratio")
  })

  output$clustering_summary_heatmap_table <- renderTable({
    summary <- clustering_heatmap_summary()
    validate(need(nrow(summary) > 0, "No clustering heatmap rows are available."))
    format_numeric_table(summary)
  }, striped = TRUE, bordered = FALSE, width = "100%")

  output$clustering_diff_summary <- renderUI({
    config <- clustering_diff_config()
    req(config)
    result <- clustering_diff_results()
    differential_summary_row(
      result,
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff
    )
  })

  output$clustering_diff_volcano <- renderPlotly({
    config <- clustering_diff_config()
    req(config)
    result <- clustering_diff_results()
    validate(need(nrow(result) > 0, "Choose two different groups with enough self-clustering data."))

    x_label <- paste("Difference in medians:", config$group_a, "minus", config$group_b, "(reference)")
    differential_volcano_plot(
      result,
      label_col = "marker",
      x_label = x_label,
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff,
      source = "clustering_diff",
      dimensions = differential_volcano_dimensions(x_label)
    )
  })

  observeEvent(plotly::event_data("plotly_click", source = "clustering_diff"), {
    event <- plotly::event_data("plotly_click", source = "clustering_diff")
    if (!is.null(event$key) && nzchar(event$key)) {
      updateSelectInput(session, "clustering_diff_marker", selected = event$key)
    }
  })

  output$clustering_diff_detail <- renderPlotly({
    data <- demo_data()
    config <- clustering_diff_config()
    req(data, config, input$clustering_diff_marker, config$group_a, config$group_b)

    plot_data <- data$clustering[
      data$clustering$marker == input$clustering_diff_marker &
        data$clustering$condition %in% c(config$group_a, config$group_b) &
        data$clustering$celltype_manual %in% config$celltype_filter,
      ,
      drop = FALSE
    ]
    validate(need(nrow(plot_data) > 0, "No self-clustering values are available for the selected marker and contrast."))

    plot_data$hover <- paste0(
      "Cell: ", plot_data$component,
      "<br>Condition: ", plot_data$condition,
      "<br>Cell type: ", plot_data$celltype_manual,
      "<br>Self-clustering log2 ratio: ", round(plot_data$log2_ratio, 3)
    )

    y_label <- paste(input$clustering_diff_marker, "self-clustering log2 ratio")
    p <- ggplot(plot_data, aes(condition, log2_ratio, color = condition, text = hover)) +
      geom_hline(yintercept = 0, color = "#8a9699", linewidth = 0.5) +
      geom_boxplot(outlier.shape = NA, alpha = 0.18, linewidth = 0.5) +
      geom_jitter(width = 0.18, height = 0, alpha = 0.5, size = 1.4) +
      labs(x = NULL, y = y_label) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(), legend.position = "none")

    if (isTRUE(config$stratify_by_celltype)) {
      p <- p + facet_wrap(~celltype_manual, scales = "free_y")
    }

    dimensions <- differential_detail_dimensions(
      plot_data,
      stratify_by_celltype = isTRUE(config$stratify_by_celltype),
      y_label = y_label
    )
    ggplotly(p, tooltip = "text", width = dimensions$width, height = dimensions$height) |>
      apply_proxiome_plot_frame(dimensions = dimensions)
  })

  output$clustering_diff_table <- renderTable({
    config <- clustering_diff_config()
    req(config)
    result <- filter_differential_hits(
      clustering_diff_results(),
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff
    )
    validate(need(nrow(result) > 0, "No clustering markers pass the selected differential thresholds."))

    format_differential_table(result, effect_label = "diff_median_vs_reference")
  }, striped = TRUE, bordered = FALSE, width = "100%")

  colocalization_heatmap_rows <- reactive({
    data <- demo_data()
    req(data)

    metadata <- colocalization_metadata()
    data$colocalization[data$colocalization$component %in% metadata$component, , drop = FALSE]
  })

  colocalization_heatmap_result <- reactive({
    data <- demo_data()
    req(data)

    colocalization <- colocalization_heatmap_rows()
    validate(need(nrow(colocalization) > 0, "No colocalization scores are available for the selected filters."))
    if (!"sample_alias" %in% names(colocalization)) {
      colocalization$sample_alias <- "sample"
    }

    scope <- input$spatial_coloc_scope %||% "condition"
    if (length(scope) == 0) {
      scope <- "condition"
    }
    if (!scope %in% c("condition", "sample", "celltype")) {
      scope <- "condition"
    }
    marker_selection_mode <- input$spatial_marker_selection_mode %||% "auto"
    if (length(marker_selection_mode) == 0) {
      marker_selection_mode <- "auto"
    }

    selected_conditions <- selected_or_all(
      input$colocalization_condition_filter,
      sort(unique(as.character(data$metadata$condition)))
    )
    selected_conditions <- intersect(selected_conditions, unique(as.character(colocalization$condition)))
    validate(need(length(selected_conditions) > 0, "Select at least one condition for the colocalization heatmap."))

    colocalization <- colocalization[colocalization$condition %in% selected_conditions, , drop = FALSE]

    if (identical(scope, "celltype")) {
      focus_celltype <- input$spatial_celltype_focus
      if (is.null(focus_celltype) || length(focus_celltype) == 0 || is.na(focus_celltype) || !nzchar(focus_celltype)) {
        focus_celltype <- sort(unique(as.character(colocalization$celltype_manual)))[1]
      }
      colocalization <- colocalization[colocalization$celltype_manual == focus_celltype, , drop = FALSE]
      validate(need(nrow(colocalization) > 0, "No colocalization scores are available for the selected cell type focus."))
    }

    available_markers <- available_colocalization_marker_choices(colocalization)
    requested_markers <- selected_or_all(input$colocalization_heatmap_markers, available_markers)
    requested_markers <- intersect(requested_markers, available_markers)
    candidate_markers <- if (identical(marker_selection_mode, "auto")) available_markers else requested_markers
    validate(need(length(candidate_markers) >= 2, "Select at least two markers for the spatial heatmap."))

    marker_summary <- spatial_heatmap_summary_for_scope(
      colocalization = colocalization,
      selected_markers = candidate_markers,
      scope = scope,
      selected_conditions = selected_conditions
    )
    selected_markers <- spatial_heatmap_selected_markers(
      summary = marker_summary,
      available_markers = candidate_markers,
      requested_markers = requested_markers,
      marker_selection_mode = marker_selection_mode,
      n_markers = input$spatial_top_marker_count,
      min_pct_detected = input$spatial_min_pct_detected,
      min_range = input$spatial_min_log2_range
    )
    validate(need(length(selected_markers) >= 2, "Select at least two markers for the spatial heatmap."))
    validate(need(length(selected_markers) <= 40, "Use 40 or fewer markers for an interpretable spatial heatmap."))

    summary <- spatial_heatmap_summary_for_scope(
      colocalization = colocalization,
      selected_markers = selected_markers,
      scope = scope,
      selected_conditions = selected_conditions
    )
    validate(need(nrow(summary) > 0, "No colocalization scores are available for the selected markers."))

    group_cols <- spatial_heatmap_group_cols(scope)
    summary <- complete_spatial_marker_pairs(
      summary = summary,
      selected_markers = selected_markers,
      group_cols = group_cols
    )
    condition_col <- if (identical(scope, "condition")) "condition" else "sample_alias"
    plot_groups <- unique(as.character(summary[[condition_col]]))
    plot_groups <- plot_groups[nzchar(plot_groups)]
    validate(need(length(plot_groups) > 0, "No spatial heatmap groups are available for the selected filters."))

    result <- make_coloc_heatmaps(
      data = summary,
      selected_markers = selected_markers,
      cell_label = spatial_heatmap_cell_label(
        scope,
        selected_celltypes = input$colocalization_celltype_filter,
        focus_celltype = input$spatial_celltype_focus
      ),
      conditions = plot_groups,
      reference_condition = if (identical(scope, "condition")) {
        input$colocalization_reference_condition %||% selected_conditions[1]
      } else {
        plot_groups[1]
      },
      condition_col = condition_col,
      clustering_method = input$colocalization_clustering_method %||% "ward.D2",
      legend_range = colocalization_legend_range(
        input$colocalization_legend_min,
        input$colocalization_legend_max
      )
    )
    result$summary <- summary
    result$scope <- scope
    result$condition_col <- condition_col
    result
  })

  output$colocalization_heatmap_interactive <- renderPlotly({
    coloc_heatmap_plotly(colocalization_heatmap_result())
  })

  output$colocalization_heatmap_original <- renderPlot({
    print(colocalization_heatmap_result()$plot)
  }, width = function() {
    coloc_heatmap_widget_dimensions(colocalization_heatmap_result()$plot_data)$width
  }, height = function() {
    coloc_heatmap_widget_dimensions(colocalization_heatmap_result()$plot_data)$height
  })

  output$colocalization_table <- renderTable({
    summary <- colocalization_heatmap_result()$summary
    validate(need(nrow(summary) > 0, "No spatial metric rows to summarize."))

    format_spatial_heatmap_table(summary)
  }, striped = TRUE, bordered = FALSE, width = "100%")

  output$colocalization_diff_summary <- renderUI({
    config <- colocalization_diff_config()
    req(config)
    result <- colocalization_diff_anchor_results()
    differential_summary_row(
      result,
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff
    )
  })

  output$colocalization_diff_volcano <- renderPlotly({
    config <- colocalization_diff_config()
    req(config)
    result <- colocalization_diff_anchor_results()
    validate(need(nrow(result) > 0, "Choose two different groups with enough colocalization data."))

    x_label <- paste("Difference in medians:", config$group_a, "minus", config$group_b, "(reference)")
    differential_volcano_plot(
      result,
      label_col = "marker_pair",
      x_label = x_label,
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff,
      source = "colocalization_diff",
      dimensions = differential_volcano_dimensions(x_label)
    )
  })

  observeEvent(plotly::event_data("plotly_click", source = "colocalization_diff"), {
    event <- plotly::event_data("plotly_click", source = "colocalization_diff")
    if (!is.null(event$key) && nzchar(event$key)) {
      updateSelectInput(session, "colocalization_diff_pair", selected = event$key)
    }
  })

  output$colocalization_diff_detail <- renderPlotly({
    data <- demo_data()
    config <- colocalization_diff_config()
    req(data, config, input$colocalization_diff_pair, config$group_a, config$group_b)

    plot_data <- data$colocalization[
      data$colocalization$marker_pair == input$colocalization_diff_pair &
        data$colocalization$condition %in% c(config$group_a, config$group_b) &
        data$colocalization$celltype_manual %in% config$celltype_filter,
      ,
      drop = FALSE
    ]
    validate(need(nrow(plot_data) > 0, "No colocalization values are available for the selected pair and contrast."))

    plot_data$hover <- paste0(
      "Cell: ", plot_data$component,
      "<br>Condition: ", plot_data$condition,
      "<br>Cell type: ", plot_data$celltype_manual,
      "<br>Colocalization log2 ratio: ", round(plot_data$log2_ratio, 3)
    )

    y_label <- paste(input$colocalization_diff_pair, "colocalization log2 ratio")
    p <- ggplot(plot_data, aes(condition, log2_ratio, color = condition, text = hover)) +
      geom_hline(yintercept = 0, color = "#8a9699", linewidth = 0.5) +
      geom_boxplot(outlier.shape = NA, alpha = 0.18, linewidth = 0.5) +
      geom_jitter(width = 0.18, height = 0, alpha = 0.5, size = 1.4) +
      labs(x = NULL, y = y_label) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(), legend.position = "none")

    if (isTRUE(config$stratify_by_celltype)) {
      p <- p + facet_wrap(~celltype_manual, scales = "free_y")
    }

    dimensions <- differential_detail_dimensions(
      plot_data,
      stratify_by_celltype = isTRUE(config$stratify_by_celltype),
      y_label = y_label
    )
    ggplotly(p, tooltip = "text", width = dimensions$width, height = dimensions$height) |>
      apply_proxiome_plot_frame(dimensions = dimensions)
  })

  output$colocalization_diff_table <- renderTable({
    config <- colocalization_diff_config()
    req(config)
    result <- filter_differential_hits(
      colocalization_diff_anchor_results(),
      fdr_cutoff = config$fdr_cutoff,
      effect_cutoff = config$effect_cutoff
    )
    validate(need(nrow(result) > 0, "No colocalization pairs pass the selected differential thresholds."))

    format_differential_table(result, effect_label = "diff_median_vs_reference")
  }, striped = TRUE, bordered = FALSE, width = "100%")
}

available_embedding_choices <- function(metadata) {
  candidates <- c(
    "Harmony UMAP" = "harmony_umap",
    "UMAP" = "umap"
  )

  available <- candidates[
    paste0(candidates, "_1") %in% names(metadata) &
      paste0(candidates, "_2") %in% names(metadata)
  ]

  if (length(available) == 0) {
    stop("No two-dimensional embedding columns are available.", call. = FALSE)
  }

  available
}

available_abundance_split_choices <- function(metadata) {
  candidates <- c(
    "None" = "",
    "Condition" = "condition",
    "Sample" = "sample"
  )

  candidates[candidates == "" | candidates %in% names(metadata)]
}

available_abundance_color_choices <- function(metadata) {
  candidates <- c(
    "Marker abundance" = "abundance",
    "Cell type" = "celltype_manual",
    "Condition" = "condition",
    "Sample" = "sample"
  )

  candidates[candidates == "abundance" | candidates %in% names(metadata)]
}

selected_abundance_color_by <- function(selected, data) {
  if (is.null(selected) || length(selected) == 0 || identical(selected, "")) {
    return("abundance")
  }
  if (identical(selected, "abundance")) {
    return("abundance")
  }
  if (!selected %in% names(data)) {
    return("abundance")
  }
  selected
}

abundance_color_label <- function(color_by) {
  labels <- c(
    abundance = "Marker abundance",
    celltype_manual = "Cell type",
    condition = "Condition",
    sample = "Sample"
  )

  label <- labels[color_by]
  if (is.na(label)) {
    return(color_by)
  }
  unname(label)
}

hover_sample_text <- function(data) {
  if (!"sample" %in% names(data)) {
    return("")
  }

  paste0("<br>Sample: ", data$sample)
}

selected_split_column <- function(selected, data) {
  if (is.null(selected) || length(selected) == 0 || identical(selected, "")) {
    return(NULL)
  }
  if (!selected %in% names(data)) {
    return(NULL)
  }
  selected
}

qc_sample_choices <- function(qc) {
  metadata <- qc$origin_metadata
  if (!"sample" %in% names(metadata)) {
    return(character(0))
  }
  sort(unique(as.character(metadata$sample)))
}

selected_qc_metadata <- function(qc, source) {
  if (identical(source, "filtered")) {
    return(qc$filtered_metadata)
  }
  qc$origin_metadata
}

filter_qc_metadata_by_sample <- function(metadata, selected_samples) {
  if (!"sample" %in% names(metadata) || is.null(selected_samples) || length(selected_samples) == 0) {
    return(metadata)
  }

  metadata[metadata$sample %in% selected_samples, , drop = FALSE]
}

qc_filter_counts_for_samples <- function(filter_counts, selected_samples, include_total = TRUE) {
  if (nrow(filter_counts) == 0 || !"sample" %in% names(filter_counts)) {
    return(filter_counts)
  }

  available_samples <- setdiff(unique(as.character(filter_counts$sample)), "TOTAL")
  selected_samples <- selected_or_all(selected_samples, available_samples)
  selected_samples <- intersect(selected_samples, available_samples)

  rows <- filter_counts[as.character(filter_counts$sample) %in% selected_samples, , drop = FALSE]
  if (nrow(rows) == 0) {
    return(filter_counts[0, , drop = FALSE])
  }

  if (isTRUE(include_total)) {
    total_rows <- if (length(selected_samples) == length(available_samples)) {
      filter_counts[as.character(filter_counts$sample) == "TOTAL", , drop = FALSE]
    } else {
      build_qc_filter_total_rows(rows, filter_counts)
    }
    rows <- rbind(rows, total_rows)
  }

  order_qc_filter_counts(rows)
}

build_qc_filter_total_rows <- function(rows, template) {
  group_cols <- intersect(c("step", "step_label"), names(rows))
  totals <- stats::aggregate(
    rows$n_cells,
    rows[, group_cols, drop = FALSE],
    sum
  )
  names(totals)[ncol(totals)] <- "n_cells"

  totals$sample <- "TOTAL"
  if ("condition" %in% names(template)) {
    totals$condition <- "TOTAL"
  }

  totals$fraction_loaded <- ifelse(totals$n_cells[1] > 0, totals$n_cells / totals$n_cells[1], NA_real_)
  align_qc_filter_count_columns(totals, template)
}

align_qc_filter_count_columns <- function(rows, template) {
  for (column in setdiff(names(template), names(rows))) {
    rows[[column]] <- NA
  }
  rows[, names(template), drop = FALSE]
}

order_qc_filter_counts <- function(filter_counts) {
  if (nrow(filter_counts) == 0 || !"step" %in% names(filter_counts)) {
    return(filter_counts)
  }

  step_order <- unique(as.character(filter_counts$step))
  filter_counts$step <- factor(filter_counts$step, levels = step_order)

  if ("sample" %in% names(filter_counts)) {
    sample_order <- unique(as.character(filter_counts$sample))
    total_index <- which(sample_order == "TOTAL")
    if (length(total_index) > 0) {
      sample_order <- c(sample_order[-total_index], "TOTAL")
    }
    filter_counts$sample <- factor(filter_counts$sample, levels = sample_order)
    filter_counts <- filter_counts[order(filter_counts$step, filter_counts$sample), , drop = FALSE]
    filter_counts$sample <- as.character(filter_counts$sample)
  } else {
    filter_counts <- filter_counts[order(filter_counts$step), , drop = FALSE]
  }

  filter_counts$step <- as.character(filter_counts$step)
  rownames(filter_counts) <- NULL
  filter_counts
}

plot_filter_cell_counts <- function(
  filter_cell_counts,
  include_total = TRUE,
  y = c("count", "fraction_loaded")
) {
  y <- match.arg(y)
  required_cols <- c("step", "n_cells")
  missing_cols <- setdiff(required_cols, names(filter_cell_counts))
  if (length(missing_cols) > 0) {
    stop("Missing columns for filter cell count plot: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  plot_df <- filter_cell_counts
  if (!"sample" %in% names(plot_df)) {
    plot_df$sample <- "All cells"
  }
  if (!"condition" %in% names(plot_df)) {
    plot_df$condition <- "All cells"
  }
  if (!"step_label" %in% names(plot_df)) {
    plot_df$step_label <- format_qc_step_label(plot_df$step)
  }

  plot_df$step <- factor(plot_df$step, levels = unique(as.character(plot_df$step)))
  plot_df$sample <- as.character(plot_df$sample)
  plot_df$condition <- as.character(plot_df$condition)

  if (!include_total) {
    plot_df <- plot_df[plot_df$sample != "TOTAL", , drop = FALSE]
  }

  group_key <- interaction(plot_df$sample, plot_df$condition, drop = TRUE, lex.order = TRUE)
  loaded_cells <- ave(plot_df$n_cells, group_key, FUN = function(values) values[1])
  plot_df$fraction_loaded_for_hover <- ifelse(loaded_cells > 0, plot_df$n_cells / loaded_cells, NA_real_)

  if (identical(y, "fraction_loaded")) {
    plot_df$value <- plot_df$fraction_loaded_for_hover
    y_label <- "Fraction of loaded cells"
    y_scale <- scale_y_continuous(labels = qc_percent_axis_labels)
  } else {
    plot_df$value <- plot_df$n_cells
    y_label <- "Number of cells"
    y_scale <- scale_y_continuous(labels = qc_count_axis_labels)
  }

  plot_df$hover <- paste0(
    "Step: ", plot_df$step_label,
    "<br>Sample: ", plot_df$sample,
    "<br>Condition: ", plot_df$condition,
    "<br>Cells: ", format(plot_df$n_cells, big.mark = ","),
    "<br>Fraction loaded: ", format_percent(plot_df$fraction_loaded_for_hover)
  )

  ggplot(plot_df, aes(x = step, y = value, group = sample, color = sample, text = hover)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    y_scale +
    labs(
      x = "Filtering step",
      y = y_label,
      color = "Sample"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

qc_percent_axis_labels <- function(values) {
  paste0(round(100 * values), "%")
}

qc_count_axis_labels <- function(values) {
  format(values, big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_qc_filter_table <- function(filter_counts) {
  if (nrow(filter_counts) == 0) {
    return(filter_counts)
  }

  filter_counts$fraction_loaded <- round(filter_counts$fraction_loaded, 3)
  filter_counts
}

format_percent <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NA_character_)
  }

  result <- rep(NA_character_, length(value))
  keep <- !is.na(value)
  result[keep] <- paste0(round(100 * value[keep], 1), "%")
  result
}

available_colocalization_marker_choices <- function(colocalization) {
  sort(unique(stats::na.omit(as.character(c(colocalization$marker_1, colocalization$marker_2)))))
}

colocalization_cell_label <- function(selected_celltypes) {
  if (is.null(selected_celltypes) || length(selected_celltypes) == 0) {
    return("selected cells")
  }
  if (length(selected_celltypes) == 1) {
    return(selected_celltypes)
  }
  paste(length(selected_celltypes), "cell types")
}

colocalization_legend_range <- function(min_value, max_value) {
  min_value <- numeric_input_value(min_value, -1)
  max_value <- numeric_input_value(max_value, 1)

  if (!is.finite(min_value) || !is.finite(max_value) || min_value == max_value) {
    return(c(-1, 1))
  }
  sort(c(min_value, max_value))
}

spatial_heatmap_group_cols <- function(scope) {
  switch(
    scope,
    sample = c("sample_alias", "condition"),
    celltype = c("sample_alias", "condition", "celltype_manual"),
    "condition"
  )
}

spatial_heatmap_summary_for_scope <- function(
  colocalization,
  selected_markers,
  scope,
  selected_conditions
) {
  if (identical(scope, "sample")) {
    return(summarize_spatial_heatmap_by_sample(
      proximity = colocalization,
      selected_markers = selected_markers
    ))
  }

  if (identical(scope, "celltype")) {
    return(summarize_spatial_heatmap_by_celltype(
      proximity = colocalization,
      selected_markers = selected_markers
    ))
  }

  summarize_colocalization_heatmap(
    colocalization,
    selected_markers = selected_markers,
    conditions = selected_conditions
  )
}

spatial_heatmap_selected_markers <- function(
  summary,
  available_markers,
  requested_markers,
  marker_selection_mode = "manual",
  n_markers = 20,
  min_pct_detected = 0.25,
  min_range = 0.2,
  max_markers = 40L
) {
  available_markers <- unique(as.character(available_markers))
  requested_markers <- intersect(unique(as.character(requested_markers)), available_markers)
  max_markers <- max(2L, as.integer(max_markers))

  if (identical(marker_selection_mode, "auto")) {
    n_markers <- numeric_input_value(n_markers, 20)
    n_markers <- max(2L, min(as.integer(n_markers), max_markers, length(available_markers)))
    return(select_spatial_heatmap_markers(
      summary = summary,
      available_markers = available_markers,
      n_markers = n_markers,
      min_pct_detected = numeric_input_value(min_pct_detected, 0.25),
      min_range = numeric_input_value(min_range, 0.2)
    ))
  }

  head(requested_markers, max_markers)
}

spatial_heatmap_cell_label <- function(scope, selected_celltypes = NULL, focus_celltype = NULL) {
  if (identical(scope, "sample")) {
    return(paste(colocalization_cell_label(selected_celltypes), "by sample"))
  }
  if (identical(scope, "celltype")) {
    if (is.null(focus_celltype) || length(focus_celltype) == 0 || is.na(focus_celltype) || !nzchar(focus_celltype)) {
      focus_celltype <- "selected cell type"
    }
    return(paste(focus_celltype, "by sample"))
  }

  colocalization_cell_label(selected_celltypes)
}

summarize_colocalization_heatmap <- function(
  colocalization,
  selected_markers,
  conditions = NULL,
  condition_col = "condition",
  marker1_col = "marker_1",
  marker2_col = "marker_2",
  value_col = "log2_ratio",
  component_col = "component"
) {
  required_cols <- c(condition_col, marker1_col, marker2_col, value_col)
  missing_cols <- setdiff(required_cols, names(colocalization))
  if (length(missing_cols) > 0) {
    stop("Missing columns for colocalization heatmap: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (nrow(colocalization) == 0 || length(selected_markers) == 0) {
    return(empty_colocalization_heatmap_summary(condition_col, marker1_col, marker2_col))
  }

  plot_data <- colocalization
  plot_data[[condition_col]] <- as.character(plot_data[[condition_col]])
  plot_data[[marker1_col]] <- as.character(plot_data[[marker1_col]])
  plot_data[[marker2_col]] <- as.character(plot_data[[marker2_col]])
  plot_data[[value_col]] <- as.numeric(plot_data[[value_col]])

  selected_markers <- as.character(selected_markers)
  if (!is.null(conditions) && length(conditions) > 0) {
    plot_data <- plot_data[plot_data[[condition_col]] %in% as.character(conditions), , drop = FALSE]
  }

  plot_data <- plot_data[
    plot_data[[marker1_col]] %in% selected_markers &
      plot_data[[marker2_col]] %in% selected_markers &
      plot_data[[marker1_col]] != plot_data[[marker2_col]] &
      is.finite(plot_data[[value_col]]),
    ,
    drop = FALSE
  ]

  if (nrow(plot_data) == 0) {
    return(empty_colocalization_heatmap_summary(condition_col, marker1_col, marker2_col))
  }

  group_cols <- c(condition_col, marker1_col, marker2_col)
  means <- aggregate(
    plot_data[[value_col]],
    plot_data[group_cols],
    function(values) mean(values, na.rm = TRUE)
  )
  names(means)[ncol(means)] <- "mean_log2_ratio"

  if (component_col %in% names(plot_data)) {
    detected <- aggregate(
      plot_data[[component_col]],
      plot_data[group_cols],
      function(values) length(unique(values))
    )
    names(detected)[ncol(detected)] <- "n_detected"

    totals <- aggregate(
      plot_data[[component_col]],
      plot_data[condition_col],
      function(values) length(unique(values))
    )
    names(totals)[ncol(totals)] <- "n_total"
  } else {
    detected <- aggregate(
      plot_data[[value_col]],
      plot_data[group_cols],
      length
    )
    names(detected)[ncol(detected)] <- "n_detected"

    totals <- aggregate(
      plot_data[[value_col]],
      plot_data[condition_col],
      length
    )
    names(totals)[ncol(totals)] <- "n_total"
  }

  summary <- merge(means, detected, by = group_cols, all.x = TRUE, sort = FALSE)
  summary <- merge(summary, totals, by = condition_col, all.x = TRUE, sort = FALSE)
  summary$pct_detected <- ifelse(summary$n_total > 0, summary$n_detected / summary$n_total, NA_real_)
  summary[, c(group_cols, "mean_log2_ratio", "pct_detected", "n_detected", "n_total"), drop = FALSE]
}

empty_colocalization_heatmap_summary <- function(condition_col, marker1_col, marker2_col) {
  result <- data.frame(
    condition = character(),
    marker_1 = character(),
    marker_2 = character(),
    mean_log2_ratio = numeric(),
    pct_detected = numeric(),
    n_detected = integer(),
    n_total = integer(),
    stringsAsFactors = FALSE
  )
  names(result)[1:3] <- c(condition_col, marker1_col, marker2_col)
  result
}

get_coloc_data <- function(
  data,
  condition_name,
  selected_markers,
  condition_col = "condition",
  marker1_col = "marker_1",
  marker2_col = "marker_2"
) {
  data[
    as.character(data[[condition_col]]) == condition_name &
      as.character(data[[marker1_col]]) %in% selected_markers &
      as.character(data[[marker2_col]]) %in% selected_markers,
    ,
    drop = FALSE
  ]
}

symmetrise_coloc_pairs <- function(data, condition_col = "condition", marker1_col = "marker_1", marker2_col = "marker_2") {
  if (nrow(data) == 0) {
    return(data)
  }

  data[[condition_col]] <- as.character(data[[condition_col]])
  data[[marker1_col]] <- as.character(data[[marker1_col]])
  data[[marker2_col]] <- as.character(data[[marker2_col]])

  reversed <- data
  reversed[[marker1_col]] <- data[[marker2_col]]
  reversed[[marker2_col]] <- data[[marker1_col]]

  combined <- rbind(data, reversed)
  key_cols <- c(condition_col, marker1_col, marker2_col)
  combined <- combined[!duplicated(combined[key_cols]), , drop = FALSE]
  rownames(combined) <- NULL
  combined
}

get_reference_marker_order <- function(
  data,
  selected_markers,
  reference_condition = "CD3CD28",
  condition_col = "condition",
  marker1_col = "marker_1",
  marker2_col = "marker_2",
  value_col = "mean_log2_ratio",
  size_col = "pct_detected",
  clustering_method = "ward.D2",
  legend_range = c(-1, 1)
) {
  selected_markers <- as.character(selected_markers)
  if (length(selected_markers) < 2) {
    return(selected_markers)
  }

  reference_data <- get_coloc_data(
    data = data,
    condition_name = reference_condition,
    selected_markers = selected_markers,
    condition_col = condition_col,
    marker1_col = marker1_col,
    marker2_col = marker2_col
  )
  if (nrow(reference_data) == 0) {
    reference_data <- data[
      as.character(data[[marker1_col]]) %in% selected_markers &
        as.character(data[[marker2_col]]) %in% selected_markers,
      ,
      drop = FALSE
    ]
  }
  if (nrow(reference_data) == 0) {
    return(selected_markers)
  }
  reference_data <- symmetrise_coloc_pairs(
    reference_data,
    condition_col = condition_col,
    marker1_col = marker1_col,
    marker2_col = marker2_col
  )

  matrix_values <- matrix(0, nrow = length(selected_markers), ncol = length(selected_markers))
  rownames(matrix_values) <- selected_markers
  colnames(matrix_values) <- selected_markers

  for (row_index in seq_len(nrow(reference_data))) {
    marker_1 <- as.character(reference_data[[marker1_col]][row_index])
    marker_2 <- as.character(reference_data[[marker2_col]][row_index])
    value <- as.numeric(reference_data[[value_col]][row_index])
    if (marker_1 %in% selected_markers && marker_2 %in% selected_markers && is.finite(value)) {
      matrix_values[marker_1, marker_2] <- value
    }
  }

  matrix_values[!is.finite(matrix_values)] <- 0
  clustering_method <- selected_clustering_method(clustering_method)

  tryCatch(
    {
      selected_markers[stats::hclust(stats::dist(matrix_values), method = clustering_method)$order]
    },
    error = function(error) selected_markers
  )
}

selected_clustering_method <- function(clustering_method) {
  valid_methods <- c("ward.D2", "complete", "average", "single", "ward.D", "mcquitty", "median", "centroid")
  if (is.null(clustering_method) || length(clustering_method) == 0 || !clustering_method %in% valid_methods) {
    return("ward.D2")
  }
  clustering_method
}

plot_coloc_fixed_order <- function(
  data,
  condition_name,
  marker_order,
  cell_label,
  selected_markers = marker_order,
  condition_col = "condition",
  marker1_col = "marker_1",
  marker2_col = "marker_2",
  value_col = "mean_log2_ratio",
  size_col = "pct_detected",
  legend_range = c(-1, 1)
) {
  plot_data <- prepare_coloc_heatmap_plot_data(
    data = data,
    marker_order = marker_order,
    conditions = condition_name,
    selected_markers = selected_markers,
    condition_col = condition_col,
    marker1_col = marker1_col,
    marker2_col = marker2_col,
    value_col = value_col,
    size_col = size_col
  )

  build_coloc_heatmap_plot(
    plot_data,
    cell_label = cell_label,
    legend_range = legend_range,
    facet = FALSE,
    title = paste0("Colocalization in ", cell_label, " (", condition_name, ")")
  )
}

make_coloc_heatmaps <- function(
  data,
  selected_markers,
  cell_label,
  conditions = c("CD3CD28", "PHA", "UNT"),
  reference_condition = "CD3CD28",
  condition_col = "condition",
  marker1_col = "marker_1",
  marker2_col = "marker_2",
  value_col = "mean_log2_ratio",
  size_col = "pct_detected",
  clustering_method = "ward.D2",
  legend_range = c(-1, 1)
) {
  selected_markers <- as.character(selected_markers)
  conditions <- as.character(conditions)
  conditions <- conditions[conditions %in% unique(as.character(data[[condition_col]]))]
  if (length(conditions) == 0) {
    conditions <- unique(as.character(data[[condition_col]]))
  }
  if (!reference_condition %in% conditions) {
    reference_condition <- conditions[1]
  }

  marker_order <- get_reference_marker_order(
    data = data,
    selected_markers = selected_markers,
    reference_condition = reference_condition,
    condition_col = condition_col,
    marker1_col = marker1_col,
    marker2_col = marker2_col,
    value_col = value_col,
    size_col = size_col,
    clustering_method = clustering_method,
    legend_range = legend_range
  )

  plot_data <- prepare_coloc_heatmap_plot_data(
    data = data,
    marker_order = marker_order,
    conditions = conditions,
    selected_markers = selected_markers,
    condition_col = condition_col,
    marker1_col = marker1_col,
    marker2_col = marker2_col,
    value_col = value_col,
    size_col = size_col
  )

  plots <- lapply(
    conditions,
    function(condition_name) {
      plot_coloc_fixed_order(
        data = data,
        condition_name = condition_name,
        marker_order = marker_order,
        cell_label = cell_label,
        selected_markers = selected_markers,
        condition_col = condition_col,
        marker1_col = marker1_col,
        marker2_col = marker2_col,
        value_col = value_col,
        size_col = size_col,
        legend_range = legend_range
      )
    }
  )
  names(plots) <- conditions

  list(
    marker_order = marker_order,
    plots = plots,
    plot_data = plot_data,
    plot = build_coloc_heatmap_plot(
      plot_data,
      cell_label = cell_label,
      legend_range = legend_range,
      facet = TRUE
    )
  )
}

coloc_heatmap_widget_dimensions <- function(plot_data) {
  condition_count <- length(unique(stats::na.omit(as.character(plot_data$plot_condition))))
  if (condition_count == 0) {
    condition_count <- 1L
  }

  margins <- coloc_heatmap_plot_margins()
  panel_px <- coloc_heatmap_panel_px()

  list(
    width = panel_px * condition_count + margins$l + margins$r,
    height = panel_px + margins$t + margins$b,
    margin = margins,
    panel_px = panel_px,
    condition_count = condition_count
  )
}

apply_coloc_heatmap_square_layout <- function(widget, dimensions) {
  widget$x$layout$autosize <- FALSE
  widget$x$layout$width <- dimensions$width
  widget$x$layout$height <- dimensions$height
  widget$x$layout$margin <- dimensions$margin

  axis_names <- names(widget$x$layout)
  xaxis_names <- axis_names[grepl("^xaxis[0-9]*$", axis_names)]
  yaxis_names <- axis_names[grepl("^yaxis[0-9]*$", axis_names)]

  for (axis_name in xaxis_names) {
    widget$x$layout[[axis_name]] <- modifyList(
      widget$x$layout[[axis_name]] %||% list(),
      list(automargin = TRUE, constrain = "domain")
    )
  }

  for (axis_name in yaxis_names) {
    suffix <- sub("^yaxis", "", axis_name)
    scaleanchor <- paste0("x", suffix)
    widget$x$layout[[axis_name]] <- modifyList(
      widget$x$layout[[axis_name]] %||% list(),
      list(
        automargin = TRUE,
        constrain = "domain",
        scaleanchor = scaleanchor,
        scaleratio = 1
      )
    )
  }

  widget
}

coloc_heatmap_plotly <- function(coloc_result, colorbar_title = "Mean log2 ratio") {
  dimensions <- coloc_heatmap_widget_dimensions(coloc_result$plot_data)
  widget <- ggplotly(
    coloc_result$plot,
    tooltip = "text",
    width = dimensions$width,
    height = dimensions$height
  )
  widget <- add_pct_detected_size_legend(widget, coloc_result$plot_data)
  widget <- apply_proxiome_plot_frame(widget, colorbar_title = colorbar_title)
  widget$x$layout$showlegend <- TRUE
  widget$x$layout$legend <- modifyList(
    widget$x$layout$legend %||% list(),
    list(title = list(text = "pct_detected"))
  )
  widget <- apply_coloc_heatmap_square_layout(widget, dimensions)
  widget
}

add_pct_detected_size_legend <- function(widget, plot_data) {
  legend_values <- pct_detected_legend_values(plot_data$plot_size)
  if (length(legend_values) == 0) {
    return(widget)
  }

  for (value in legend_values) {
    widget <- plotly::add_trace(
      widget,
      x = 0,
      y = 0,
      type = "scatter",
        mode = "markers",
        name = paste0("pct_detected ", format_percent(value)),
        marker = list(
          size = pct_detected_marker_size(value, plot_data),
          color = "rgba(80, 96, 100, 0.45)",
          line = list(color = "#263238", width = 1)
        ),
      hoverinfo = "skip",
      showlegend = TRUE,
      visible = "legendonly",
      inherit = FALSE
    )
  }

  widget
}

pct_detected_legend_values <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values) & values > 0]
  if (length(values) == 0) {
    return(numeric(0))
  }

  breaks <- c(0.25, 0.5, 0.75, 1)
  breaks[breaks <= max(values, na.rm = TRUE)]
}

pct_detected_marker_size <- function(value, plot_data = NULL) {
  value <- sqrt(pmin(pmax(value, 0), 1))
  size_range <- if (is.null(plot_data)) c(1.5, 7) else coloc_heatmap_size_range(plot_data)
  ggplot_to_plotly_size_ratio <- 3.78
  ggplot_to_plotly_size_ratio * (size_range[1] + diff(size_range) * value)
}

prepare_coloc_heatmap_plot_data <- function(
  data,
  marker_order,
  conditions,
  selected_markers = marker_order,
  condition_col = "condition",
  marker1_col = "marker_1",
  marker2_col = "marker_2",
  value_col = "mean_log2_ratio",
  size_col = "pct_detected"
) {
  grid <- expand.grid(
    condition = conditions,
    marker_1 = marker_order,
    marker_2 = marker_order,
    stringsAsFactors = FALSE
  )
  names(grid) <- c(condition_col, marker1_col, marker2_col)
  grid <- grid[grid[[marker1_col]] != grid[[marker2_col]], , drop = FALSE]

  data <- data[
    as.character(data[[condition_col]]) %in% conditions &
      as.character(data[[marker1_col]]) %in% selected_markers &
      as.character(data[[marker2_col]]) %in% selected_markers,
    ,
    drop = FALSE
  ]
  data <- symmetrise_coloc_pairs(
    data,
    condition_col = condition_col,
    marker1_col = marker1_col,
    marker2_col = marker2_col
  )

  merge_cols <- intersect(
    unique(c(condition_col, marker1_col, marker2_col, value_col, size_col, "n_detected", "n_total")),
    names(data)
  )

  plot_data <- merge(
    grid,
    data[, merge_cols, drop = FALSE],
    by = c(condition_col, marker1_col, marker2_col),
    all.x = TRUE,
    sort = FALSE
  )
  if (!"n_detected" %in% names(plot_data)) {
    plot_data$n_detected <- NA_integer_
  }
  if (!"n_total" %in% names(plot_data)) {
    plot_data$n_total <- NA_integer_
  }

  plot_data$plot_condition <- factor(plot_data[[condition_col]], levels = conditions)
  plot_data$plot_marker_1 <- factor(plot_data[[marker1_col]], levels = marker_order)
  plot_data$plot_marker_2 <- factor(plot_data[[marker2_col]], levels = rev(marker_order))
  plot_data$plot_value <- as.numeric(plot_data[[value_col]])
  plot_data$plot_size <- as.numeric(plot_data[[size_col]])
  condition_label <- if (identical(condition_col, "sample_alias")) "Sample" else "Condition"
  plot_data$hover <- paste0(
    condition_label, ": ", plot_data[[condition_col]],
    "<br>Marker 1: ", plot_data[[marker1_col]],
    "<br>Marker 2: ", plot_data[[marker2_col]],
    "<br>Mean log2 ratio: ", ifelse(is.na(plot_data$plot_value), "No data", round(plot_data$plot_value, 3)),
    "<br>Detected cells: ", plot_data$n_detected,
    "<br>Fraction detected: ", format_percent(plot_data$plot_size)
  )

  plot_data
}

coloc_heatmap_size_range <- function(plot_data) {
  row_count <- length(unique(stats::na.omit(as.character(plot_data$plot_marker_2))))
  if (row_count == 0) {
    return(c(1.5, 7))
  }

  height_px <- suppressWarnings(as.numeric(sub("px$", "", proxiome_plot_height())))
  if (!is.finite(height_px)) {
    height_px <- 430
  }
  margins <- proxiome_plot_margins()
  panel_height_px <- max(160, height_px - margins$t - margins$b)
  row_step_px <- panel_height_px / row_count
  ggplot_to_plotly_size_ratio <- 3.78
  upper <- min(7, max(0.8, (row_step_px * 0.75) / ggplot_to_plotly_size_ratio))
  lower <- min(1.5, upper * 0.4)

  c(lower, upper)
}

build_coloc_heatmap_plot <- function(plot_data, cell_label, legend_range = c(-1, 1), facet = TRUE, title = NULL) {
  size_range <- coloc_heatmap_size_range(plot_data)

  p <- ggplot(plot_data, aes(plot_marker_1, plot_marker_2, fill = plot_value, size = plot_size, text = hover)) +
    geom_point(shape = 21, color = "#263238", stroke = 0.18, alpha = 0.9, na.rm = TRUE) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    scale_fill_gradient2(
      low = "#176d73",
      mid = "#f7f8f7",
      high = "#c7503e",
      midpoint = 0,
      limits = legend_range,
      oob = squish_to_limits,
      na.value = "#e3e8e7",
      name = "Mean log2 ratio"
    ) +
    scale_size_continuous(
      range = size_range,
      limits = c(0, 1),
      breaks = c(0.25, 0.5, 0.75, 1),
      labels = qc_percent_axis_labels,
      name = "Detected"
    ) +
    coord_fixed(ratio = 1) +
    labs(
      title = title %||% paste("Colocalization in", cell_label),
      x = "Marker 1",
      y = "Marker 2"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#263238", fill = NA, linewidth = 0.6),
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )

  if (isTRUE(facet)) {
    p <- p + facet_wrap(~plot_condition, nrow = 1)
  }

  p
}

squish_to_limits <- function(x, range) {
  pmin(pmax(x, range[1]), range[2])
}

rank_qc_metadata <- function(metadata, value_col) {
  if (!"sample" %in% names(metadata)) {
    metadata$sample <- "All cells"
  }

  ranked <- lapply(split(metadata, metadata$sample), function(sample_data) {
    sample_data <- sample_data[order(sample_data[[value_col]], decreasing = TRUE), , drop = FALSE]
    sample_data$rank <- seq_len(nrow(sample_data))
    sample_data
  })

  result <- do.call(rbind, ranked)
  rownames(result) <- NULL
  result
}

qc_molecule_rank_plotly <- function(metadata, cutoff = 10000) {
  if (!"n_umi" %in% names(metadata)) {
    stop("The metadata does not include n_umi.", call. = FALSE)
  }

  ranked <- rank_qc_metadata(metadata, value_col = "n_umi")
  ranked <- ranked[is.finite(ranked$rank) & is.finite(ranked$n_umi) & ranked$rank > 0 & ranked$n_umi > 0, , drop = FALSE]
  if (nrow(ranked) == 0) {
    return(apply_proxiome_plot_frame(plotly::plot_ly()))
  }

  ranked$hover <- paste0(
    "Cell: ", ranked$component,
    "<br>Sample: ", ranked$sample,
    "<br>Rank: ", format(ranked$rank, big.mark = ","),
    "<br>n_umi: ", format(round(ranked$n_umi, 1), big.mark = ",")
  )

  widget <- plotly::plot_ly()
  for (sample_name in unique(ranked$sample)) {
    sample_data <- ranked[ranked$sample == sample_name, , drop = FALSE]
    widget <- plotly::add_trace(
      widget,
      data = sample_data,
      x = ~rank,
      y = ~n_umi,
      text = ~hover,
      name = sample_name,
      type = "scatter",
      mode = "lines+markers",
      hoverinfo = "text",
      line = list(width = 1.1),
      marker = list(size = 4.5, opacity = 0.7),
      inherit = FALSE
    )
  }

  cutoff <- numeric_input_value(cutoff, 10000)
  if (is.finite(cutoff) && cutoff > 0) {
    rank_range <- range(ranked$rank, na.rm = TRUE)
    widget <- plotly::add_trace(
      widget,
      x = rank_range,
      y = rep(cutoff, 2),
      text = paste0("n_umi cutoff: ", format(cutoff, big.mark = ",")),
      name = "n_umi cutoff",
      type = "scatter",
      mode = "lines",
      hoverinfo = "text",
      line = list(color = "#c7503e", dash = "dash", width = 1.2),
      inherit = FALSE
    )
  }

  plotly::layout(
    widget,
    xaxis = list(title = "Cell rank", type = "log"),
    yaxis = list(title = "n_umi", type = "log"),
    legend = list(title = list(text = "Sample"))
  ) |>
    apply_proxiome_plot_frame()
}

qc_distribution_plotly <- function(metadata, metric, isotype_cutoff = 0.001) {
  metric <- as.character(metric)[1]
  plot_data <- qc_distribution_plot_data(metadata, metric)
  metric_label <- qc_metric_label(metric)

  widget <- plotly::plot_ly()
  if (nrow(plot_data) == 0) {
    return(
      plotly::layout(
        widget,
        xaxis = list(title = ""),
        yaxis = list(title = metric_label)
      ) |>
        apply_proxiome_plot_frame()
    )
  }

  for (sample_name in unique(plot_data$sample)) {
    sample_data <- plot_data[plot_data$sample == sample_name, , drop = FALSE]
    widget <- plotly::add_trace(
      widget,
      x = rep(sample_name, nrow(sample_data)),
      y = sample_data$metric_value,
      text = sample_data$hover,
      name = sample_name,
      type = "violin",
      hoverinfo = "text",
      points = "all",
      jitter = 0.28,
      pointpos = 0,
      box = list(visible = TRUE),
      meanline = list(visible = TRUE),
      marker = list(size = 3.5, opacity = 0.35),
      line = list(width = 1),
      spanmode = "hard",
      inherit = FALSE
    )
  }

  yaxis <- list(title = metric_label)
  if (is_qc_log_metric(metric)) {
    yaxis$type <- "log"
  }

  shapes <- list()
  if (identical(metric, "isotype_fraction") && is.finite(isotype_cutoff)) {
    shapes <- list(list(
      type = "line",
      xref = "paper",
      x0 = 0,
      x1 = 1,
      y0 = isotype_cutoff,
      y1 = isotype_cutoff,
      line = list(color = "#c7503e", dash = "dash", width = 1.2)
    ))
  }

  plotly::layout(
    widget,
    xaxis = list(title = "", categoryorder = "array", categoryarray = unique(plot_data$sample)),
    yaxis = yaxis,
    showlegend = FALSE,
    shapes = shapes
  ) |>
    apply_proxiome_plot_frame()
}

qc_distribution_plot_data <- function(metadata, metric) {
  if (!metric %in% names(metadata)) {
    stop("Selected QC metric is not available: ", metric, call. = FALSE)
  }

  plot_data <- metadata
  if (!"component" %in% names(plot_data)) {
    plot_data$component <- seq_len(nrow(plot_data))
  }
  if (!"sample" %in% names(plot_data)) {
    plot_data$sample <- "All cells"
  }

  metric_values <- numeric_metric_vector(plot_data[[metric]], expected_length = nrow(plot_data), metric = metric)
  keep <- is.finite(metric_values)
  if (is_qc_log_metric(metric)) {
    keep <- keep & metric_values > 0
  }

  plot_data <- plot_data[keep, , drop = FALSE]
  plot_data$metric_value <- metric_values[keep]
  plot_data$sample <- as.character(plot_data$sample)

  metric_label <- qc_metric_label(metric)
  plot_data$hover <- paste0(
    "Cell: ", plot_data$component,
    "<br>Sample: ", plot_data$sample,
    "<br>", metric_label, ": ", signif(plot_data$metric_value, 4)
  )

  rownames(plot_data) <- NULL
  plot_data
}

numeric_metric_vector <- function(values, expected_length, metric) {
  if (is.data.frame(values)) {
    if (ncol(values) != 1) {
      stop("QC metric column is not one-dimensional: ", metric, call. = FALSE)
    }
    values <- values[[1]]
  }

  if (is.factor(values)) {
    values <- as.character(values)
  }
  if (is.list(values) && !is.atomic(values)) {
    values <- unlist(values, use.names = FALSE)
  }

  values <- as.vector(values)
  if (length(values) != expected_length) {
    stop("QC metric length does not match metadata rows: ", metric, call. = FALSE)
  }

  suppressWarnings(as.numeric(values))
}

qc_metric_label <- function(metric) {
  choices <- c(
    n_umi = "UMIs",
    n_edges = "Edges",
    reads_in_component = "Reads in component",
    isotype_fraction = "Isotype fraction",
    tau = "Tau"
  )

  label <- choices[metric]
  if (is.na(label)) {
    return(metric)
  }
  unname(label)
}

is_qc_log_metric <- function(metric) {
  metric %in% c("n_umi", "n_edges", "reads_in_component")
}

format_numeric_table <- function(data, digits = 3) {
  if (nrow(data) == 0) {
    return(data)
  }
  numeric_cols <- vapply(data, is.numeric, logical(1))
  data[numeric_cols] <- lapply(data[numeric_cols], function(values) round(values, digits))
  data
}

plot_dimension_override <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) == 0 || !is.finite(value[1]) || is.na(value[1]) || value[1] <= 0) {
    return(NULL)
  }

  round(value[1])
}

count_observed_values <- function(values) {
  count <- length(unique(stats::na.omit(as.character(values))))
  if (count == 0) {
    return(1L)
  }

  count
}

bounded_integer <- function(value, lower, upper) {
  value <- suppressWarnings(as.integer(round(as.numeric(value))))
  if (length(value) == 0 || !is.finite(value[1]) || is.na(value[1])) {
    return(lower)
  }

  min(max(value[1], lower), upper)
}

differential_volcano_dimensions <- function(x_label, legend = TRUE) {
  label_width <- nchar(as.character(x_label %||% ""), type = "width")
  width <- bounded_integer(360 + label_width * 7, lower = 620, upper = 1100)
  bottom_margin <- bounded_integer(86 + label_width * 1.2, lower = 130, upper = 180)
  right_margin <- if (isTRUE(legend)) 160 else 96

  list(
    width = width,
    height = 500,
    margin = list(l = 92, r = right_margin, t = 56, b = bottom_margin)
  )
}

differential_detail_dimensions <- function(plot_data, stratify_by_celltype = FALSE, y_label = NULL) {
  group_count <- if ("condition" %in% names(plot_data)) count_observed_values(plot_data$condition) else 2L
  facet_count <- if (isTRUE(stratify_by_celltype) && "celltype_manual" %in% names(plot_data)) {
    count_observed_values(plot_data$celltype_manual)
  } else {
    1L
  }
  facet_cols <- if (facet_count <= 1) 1L else ceiling(sqrt(facet_count))
  facet_rows <- ceiling(facet_count / facet_cols)
  y_label_width <- nchar(as.character(y_label %||% ""), type = "width")
  left_margin <- bounded_integer(70 + y_label_width * 1.2, lower = 90, upper = 180)
  panel_width <- max(260, group_count * 150)
  panel_height <- 250

  list(
    width = max(560, panel_width * facet_cols + left_margin + 72),
    height = max(500, panel_height * facet_rows + 170),
    margin = list(l = left_margin, r = 72, t = 56, b = 110),
    facet_count = facet_count,
    facet_cols = facet_cols,
    facet_rows = facet_rows,
    group_count = group_count
  )
}

facet_column_override <- function(value, facet_count) {
  column_count <- plot_dimension_override(value)
  if (is.null(column_count)) {
    return(NULL)
  }

  max(1L, min(as.integer(column_count), as.integer(max(1L, facet_count))))
}

abundance_distribution_widget_dimensions <- function(plot_data, facet_cols = NULL, width_px = NULL, height_px = NULL) {
  if ("sample_alias" %in% names(plot_data)) {
    sample_values <- plot_data$sample_alias
  } else if ("sample" %in% names(plot_data)) {
    sample_values <- plot_data$sample
  } else {
    sample_values <- "sample"
  }

  celltype_values <- if ("celltype_manual" %in% names(plot_data)) plot_data$celltype_manual else "cell type"
  sample_count <- count_observed_values(sample_values)
  celltype_count <- count_observed_values(celltype_values)
  auto_facet_cols <- if (celltype_count <= 1) 1L else ceiling(sqrt(celltype_count))
  facet_cols <- facet_column_override(facet_cols, celltype_count) %||% auto_facet_cols
  facet_rows <- ceiling(celltype_count / facet_cols)
  margins <- proxiome_plot_margins()
  panel_width <- max(220, sample_count * 70)
  panel_height <- 170
  auto_width <- max(520, panel_width * facet_cols + margins$l + margins$r)
  auto_height <- max(430, panel_height * facet_rows + margins$t + margins$b)
  width_override <- plot_dimension_override(width_px)
  height_override <- plot_dimension_override(height_px)

  list(
    width = width_override %||% auto_width,
    height = height_override %||% auto_height,
    margin = margins,
    sample_count = sample_count,
    celltype_count = celltype_count,
    facet_cols = facet_cols,
    facet_rows = facet_rows
  )
}

plot_abundance_marker_distribution <- function(plot_data, marker, facet_cols = NULL, show_jitter = TRUE) {
  if (!"sample_alias" %in% names(plot_data)) {
    plot_data$sample_alias <- if ("sample" %in% names(plot_data)) plot_data$sample else "sample"
  }
  p <- ggplot(plot_data, aes(sample_alias, abundance, fill = condition)) +
    geom_violin(scale = "width", color = NA, alpha = 0.85, na.rm = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.22, na.rm = TRUE) +
    labs(x = "Sample", y = paste(marker, "normalized abundance"), fill = "Condition") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid = element_blank()
    )

  if (isTRUE(show_jitter)) {
    p <- p + geom_jitter(width = 0.14, height = 0, alpha = 0.28, size = 0.5, na.rm = TRUE)
  }

  celltype_count <- count_observed_values(plot_data$celltype_manual)
  if (celltype_count > 1) {
    facet_ncol <- facet_column_override(facet_cols, celltype_count)
    if (is.null(facet_ncol)) {
      p <- p + facet_wrap(~celltype_manual, scales = "free_y")
    } else {
      p <- p + facet_wrap(~celltype_manual, scales = "free_y", ncol = facet_ncol)
    }
  }
  p
}

celltype_composition_data <- function(metadata) {
  required_cols <- c("component", "condition", "celltype_manual")
  if (!all(required_cols %in% names(metadata))) {
    return(data.frame())
  }
  plot_data <- metadata
  if (!"sample_alias" %in% names(plot_data)) {
    plot_data$sample_alias <- if ("sample" %in% names(plot_data)) plot_data$sample else "sample"
  }

  counts <- aggregate(
    plot_data$component,
    plot_data[c("sample_alias", "condition", "celltype_manual")],
    function(values) length(unique(values))
  )
  names(counts)[ncol(counts)] <- "n"
  totals <- aggregate(counts$n, counts[c("sample_alias", "condition")], sum)
  names(totals)[ncol(totals)] <- "total_cells"
  out <- merge(counts, totals, by = c("sample_alias", "condition"), all.x = TRUE, sort = FALSE)
  out$frac <- ifelse(out$total_cells > 0, out$n / out$total_cells, NA_real_)
  out$hover <- paste0(
    "Sample: ", out$sample_alias,
    "<br>Condition: ", out$condition,
    "<br>Cell type: ", out$celltype_manual,
    "<br>Cells: ", out$n,
    "<br>Fraction: ", format_percent(out$frac)
  )
  out
}

plot_celltype_composition <- function(composition) {
  ggplot(composition, aes(sample_alias, frac, fill = celltype_manual, text = hover)) +
    geom_col(position = "stack") +
    scale_y_continuous(labels = qc_percent_axis_labels, expand = expansion(0)) +
    labs(x = "Sample", y = "% cells", fill = "Cell type") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid = element_blank()
    )
}

annotation_heatmap_data <- function(abundance, metadata, max_markers = 40L) {
  if (nrow(abundance) == 0 || !"celltype_manual" %in% names(metadata)) {
    return(data.frame())
  }
  plot_data <- merge(
    abundance,
    metadata[, intersect(c("component", "celltype_manual"), names(metadata)), drop = FALSE],
    by = "component",
    all.x = FALSE,
    sort = FALSE
  )
  if (nrow(plot_data) == 0) {
    return(data.frame())
  }

  summary <- aggregate(
    plot_data$abundance,
    plot_data[c("celltype_manual", "marker")],
    function(values) stats::median(values, na.rm = TRUE)
  )
  names(summary)[ncol(summary)] <- "median_abundance"

  marker_scores <- aggregate(
    summary$median_abundance,
    summary["marker"],
    function(values) if (length(values) <= 1) 0 else stats::sd(values, na.rm = TRUE)
  )
  names(marker_scores)[ncol(marker_scores)] <- "sd"
  marker_scores <- marker_scores[order(marker_scores$sd, decreasing = TRUE), , drop = FALSE]
  selected_markers <- head(marker_scores$marker, min(max_markers, nrow(marker_scores)))
  summary <- summary[summary$marker %in% selected_markers, , drop = FALSE]
  summary$marker <- factor(summary$marker, levels = selected_markers)
  summary$hover <- paste0(
    "Cell type: ", summary$celltype_manual,
    "<br>Marker: ", summary$marker,
    "<br>Median abundance: ", round(summary$median_abundance, 3)
  )
  summary
}

plot_annotation_heatmap <- function(plot_data) {
  ggplot(plot_data, aes(marker, celltype_manual, fill = median_abundance, text = hover)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradientn(colors = c("#f7f8f7", "#78aeb2", "#f0b45b", "#c7503e"), na.value = "#e3e8e7") +
    labs(x = "Marker", y = "Cell type", fill = "Median abundance") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank()
    )
}

plot_clustering_per_marker <- function(plot_data, marker) {
  p <- ggplot(plot_data, aes(sample_alias, log2_ratio, fill = condition)) +
    geom_hline(yintercept = 0, color = "#8a9699", linewidth = 0.5) +
    geom_violin(scale = "width", color = NA, alpha = 0.85, na.rm = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.22, na.rm = TRUE) +
    geom_jitter(width = 0.14, height = 0, alpha = 0.28, size = 0.5, na.rm = TRUE) +
    labs(x = "Sample", y = paste(marker, "self-clustering log2 ratio"), fill = "Condition") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid = element_blank()
    )

  if (length(unique(stats::na.omit(plot_data$celltype_manual))) > 1) {
    p <- p + facet_wrap(~celltype_manual, scales = "free_y")
  }
  p
}

plot_clustering_summary_heatmap <- function(summary) {
  limit <- max(abs(summary$mean_log2_ratio), na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) {
    limit <- 1
  }
  summary$hover <- paste0(
    "Condition: ", summary$condition,
    "<br>Cell type: ", summary$celltype_manual,
    "<br>Marker: ", summary$marker,
    "<br>Mean log2 ratio: ", round(summary$mean_log2_ratio, 3),
    "<br>Cells: ", summary$n_cells
  )

  ggplot(summary, aes(marker, celltype_manual, fill = mean_log2_ratio, text = hover)) +
    geom_tile(color = "white", linewidth = 0.35) +
    facet_grid(condition ~ ., scales = "free", space = "free") +
    scale_fill_gradient2(
      low = "#176d73",
      mid = "#f7f8f7",
      high = "#c7503e",
      midpoint = 0,
      limits = c(-limit, limit),
      oob = squish_to_limits,
      na.value = "#e3e8e7"
    ) +
    labs(x = "Marker", y = "Cell type", fill = "Mean log2 ratio") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

selected_or_all <- function(selected, all_values) {
  if (is.null(selected) || length(selected) == 0) {
    return(all_values)
  }
  selected
}

add_embedding_columns <- function(data, embedding) {
  x_col <- paste0(embedding, "_1")
  y_col <- paste0(embedding, "_2")

  if (!x_col %in% names(data) || !y_col %in% names(data)) {
    stop("Selected embedding is not available: ", embedding, call. = FALSE)
  }

  data$embedding_x <- data[[x_col]]
  data$embedding_y <- data[[y_col]]
  data
}

metric_row <- function(...) {
  layout_column_wrap(
    width = "180px",
    fill = FALSE,
    class = "summary-box-row",
    ...
  )
}

metric_tile <- function(label, value) {
  value_box(
    title = label,
    value = value,
    theme = "text-primary",
    fill = FALSE
  )
}

format_summary_table <- function(summary, value_label) {
  if (nrow(summary) == 0) {
    return(summary)
  }

  summary$mean_value <- round(summary$mean_value, 3)
  summary$median_value <- round(summary$median_value, 3)
  names(summary)[names(summary) == "mean_value"] <- value_label
  names(summary)[names(summary) == "median_value"] <- sub("^mean", "median", value_label)
  summary
}

format_spatial_heatmap_table <- function(summary, max_rows = 30) {
  if (nrow(summary) == 0) {
    return(summary)
  }

  summary <- summary[order(abs(summary$mean_log2_ratio), summary$mean_log2_ratio, decreasing = TRUE), , drop = FALSE]
  summary <- head(summary, max_rows)

  numeric_cols <- intersect(c("mean_log2_ratio", "pct_detected"), names(summary))
  for (col in numeric_cols) {
    summary[[col]] <- round(as.numeric(summary[[col]]), 3)
  }
  count_cols <- intersect(c("n_detected", "n_total"), names(summary))
  for (col in count_cols) {
    summary[[col]] <- as.integer(summary[[col]])
  }

  display_cols <- intersect(
    c(
      "sample_alias",
      "condition",
      "celltype_manual",
      "marker_1",
      "marker_2",
      "mean_log2_ratio",
      "pct_detected",
      "n_detected",
      "n_total"
    ),
    names(summary)
  )
  summary[, display_cols, drop = FALSE]
}

numeric_input_value <- function(value, default) {
  if (is.null(value) || length(value) == 0 || is.na(value)) {
    return(default)
  }
  value
}

differential_summary_row <- function(result, fdr_cutoff, effect_cutoff) {
  hits <- filter_differential_hits(result, fdr_cutoff = fdr_cutoff, effect_cutoff = effect_cutoff)
  up <- sum(hits$effect_size > 0, na.rm = TRUE)
  down <- sum(hits$effect_size < 0, na.rm = TRUE)
  eligible <- sum(!is.na(result$p_adj))

  metric_row(
    metric_tile("Tested", format(nrow(result), big.mark = ",")),
    metric_tile("With FDR", format(eligible, big.mark = ",")),
    metric_tile("Higher in A", format(up, big.mark = ",")),
    metric_tile("Higher in B", format(down, big.mark = ","))
  )
}

filter_differential_hits <- function(result, fdr_cutoff, effect_cutoff) {
  if (nrow(result) == 0) {
    return(result)
  }

  keep <- !is.na(result$p_adj) &
    result$p_adj <= fdr_cutoff &
    !is.na(result$effect_size) &
    abs(result$effect_size) >= effect_cutoff

  result <- result[keep, , drop = FALSE]
  order_differential_results(result)
}

order_differential_results <- function(result) {
  if (nrow(result) == 0) {
    return(result)
  }

  result[order(is.na(result$p_adj), result$p_adj, -abs(result$effect_size)), , drop = FALSE]
}

differential_volcano_plot <- function(
  result,
  label_col,
  x_label,
  fdr_cutoff,
  effect_cutoff,
  source = NULL,
  dimensions = differential_volcano_dimensions(x_label)
) {
  plot_data <- prepare_differential_plot_data(result, label_col, fdr_cutoff, effect_cutoff)
  fdr_line <- -log10(max(fdr_cutoff, .Machine$double.xmin))

  p <- ggplot(plot_data, aes(effect_size, neg_log10_fdr, color = threshold_status, text = hover, key = plot_key)) +
    geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), color = "#b3bdbf", linetype = "dashed", linewidth = 0.45) +
    geom_hline(yintercept = fdr_line, color = "#b3bdbf", linetype = "dashed", linewidth = 0.45) +
    geom_point(size = 2.3, alpha = 0.78) +
    scale_color_manual(
      values = c(
        "Higher in A" = "#c7503e",
        "Higher in B" = "#176d73",
        "Not significant" = "#9aa5a8"
      ),
      name = NULL
    ) +
    labs(x = x_label, y = "-log10(FDR)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  ggplotly(p, tooltip = "text", source = source, width = dimensions$width, height = dimensions$height) |>
    apply_differential_plot_frame(dimensions = dimensions)
}

prepare_differential_plot_data <- function(result, label_col, fdr_cutoff, effect_cutoff) {
  plot_data <- result
  plot_data$p_adj_for_plot <- ifelse(is.na(plot_data$p_adj), 1, pmax(plot_data$p_adj, .Machine$double.xmin))
  plot_data$neg_log10_fdr <- -log10(plot_data$p_adj_for_plot)
  plot_data$threshold_status <- differential_threshold_status(plot_data, fdr_cutoff, effect_cutoff)
  plot_data$plot_label <- paste(plot_data[[label_col]], plot_data$celltype_manual, sep = " | ")
  plot_data$plot_key <- plot_data[[label_col]]
  plot_data$hover <- paste0(
    plot_data[[label_col]],
    "<br>Cell type: ", plot_data$celltype_manual,
    "<br>Effect: ", round(plot_data$effect_size, 3),
    "<br>FDR: ", format_p_value(plot_data$p_adj),
    "<br>", plot_data$group_a, " mean: ", round(plot_data$mean_a, 3),
    "<br>", plot_data$group_b, " reference mean: ", round(plot_data$mean_b, 3),
    "<br>", plot_data$direction
  )
  plot_data
}

differential_threshold_status <- function(result, fdr_cutoff, effect_cutoff) {
  is_hit <- !is.na(result$p_adj) &
    result$p_adj <= fdr_cutoff &
    !is.na(result$effect_size) &
    abs(result$effect_size) >= effect_cutoff

  ifelse(
    is_hit,
    ifelse(result$effect_size > 0, "Higher in A", "Higher in B"),
    "Not significant"
  )
}

format_differential_table <- function(result, effect_label, max_rows = 30) {
  result <- order_differential_results(result)
  result <- head(result, max_rows)

  numeric_cols <- intersect(c("mean_a", "mean_b", "median_a", "median_b", "effect_size"), names(result))
  for (col in numeric_cols) {
    result[[col]] <- round(result[[col]], 3)
  }
  result$p_value <- format_p_value(result$p_value)
  result$p_adj <- format_p_value(result$p_adj)
  names(result)[names(result) == "effect_size"] <- effect_label
  result
}

format_p_value <- function(values) {
  ifelse(
    is.na(values),
    NA_character_,
    formatC(values, format = "e", digits = 2)
  )
}

shinyApp(ui, server)
