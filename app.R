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
  packages <- c("shiny", "bslib", "ggplot2", "plotly", "future", "promises", "rlang", "data.table")
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
source(file.path(APP_DIR, "R", "differential_helpers.R"), local = TRUE)
source(file.path(APP_DIR, "R", "environment_module.R"), local = TRUE)
source(file.path(APP_DIR, "R", "data_source_module.R"), local = TRUE)
source(file.path(APP_DIR, "R", "qc_module.R"), local = TRUE)
source(file.path(APP_DIR, "R", "abundance_module.R"), local = TRUE)
source(file.path(APP_DIR, "R", "clustering_module.R"), local = TRUE)
source(file.path(APP_DIR, "R", "colocalization_module.R"), local = TRUE)

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

    .plot-download-controls {
      display: flex;
      justify-content: flex-end;
      gap: 8px;
      margin-bottom: 8px;
    }

    .plot-download-button {
      min-width: 58px;
      padding: 3px 10px;
      font-size: 0.78rem;
      line-height: 1.35;
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

    .umap-plot-pane {
      overflow-x: auto;
      overflow-y: visible;
    }

    .umap-plot-shell {
      max-width: none;
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
      function byIdOrSuffix(id) {
        return document.getElementById(id) || document.querySelector('[id$=\"' + id + '\"]');
      }

      function setRdsLoadState(message) {
        message = message || {};
        var state = message.state || 'idle';
        var loading = state === 'running' || message.disabled === true;
        var button = byIdOrSuffix('load_rds_path');
        var status = byIdOrSuffix('rds_load_status');
        var progress = byIdOrSuffix('rds_load_progress');
        var progressBar = byIdOrSuffix('rds_load_progress_bar');
        var elapsed = byIdOrSuffix('rds_load_elapsed');

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

        var button = target.closest('#load_rds_path, [id$=\"load_rds_path\"]');
        if (!button || button.disabled) {
          return;
        }

        window.setTimeout(function() {
          var pathInput = byIdOrSuffix('rds_server_path');
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

user_rds_path_loading_enabled <- function(platform = APP_PLATFORM) {
  disabled <- tolower(Sys.getenv("PROXIOME_DISABLE_USER_RDS", unset = "false")) %in% c("1", "true", "yes")
  platform %in% c("ccrsf_hpc", "biowulf_hpc", "portable") && !disabled
}

plot_pane <- function(
  ...,
  size = c("standard", "compact", "wide", "scroll"),
  extra_class = NULL,
  download_id = NULL,
  ns = identity
) {
  size <- match.arg(size)
  children <- list(...)
  if (!is.null(download_id) && nzchar(download_id)) {
    children <- c(list(plot_download_controls(ns, download_id)), children)
  }

  do.call(
    div,
    c(
      list(class = paste(c("plot-pane", paste0("plot-pane-", size), extra_class), collapse = " ")),
      children
    )
  )
}

detail_pane <- function(..., size = c("standard", "compact", "wide")) {
  size <- match.arg(size)
  div(
    class = paste(c("detail-grid", paste0("plot-pane-", size)), collapse = " "),
    ...
  )
}

build_app_ui <- function(show_environment = environment_diagnostics_enabled()) {
  show_environment <- isTRUE(show_environment)
  fillable_panels <- c("QC", "Abundance", "Spatial Metrics")
  nav_items <- list(
    qc_module_ui("qc"),
    abundance_module_ui("abundance"),
    nav_panel(
      "Spatial Metrics",
      navset_tab(
        id = "spatial_metric_readout",
        clustering_module_ui("clustering"),
        colocalization_module_ui("colocalization")
      )
    )
  )

  if (show_environment) {
    fillable_panels <- c(fillable_panels, "Environment")
    nav_items <- c(nav_items, list(
      nav_panel(
        "Environment",
        environment_module_ui("environment")
      )
    ))
  }

  do.call(
    page_navbar,
    c(
      list(
        title = "ProxiomeVis",
        id = "readout_tab",
        fillable = fillable_panels,
        theme = bs_theme(
          version = 5,
          bg = "#f6f8f7",
          fg = "#192124",
          primary = "#176d73",
          secondary = "#c7503e",
          base_font = "system-ui"
        ),
        header = tagList(app_css(), app_js())
      ),
      nav_items,
      list(
        nav_spacer(),
        data_source_module_ui("data_source")
      )
    )
  )
}

ui <- build_app_ui()

server <- function(input, output, session) {
  data_source <- data_source_module_server("data_source", app_dir = APP_DIR)
  demo_data <- data_source$data
  qc_module_server("qc", data = demo_data)
  abundance_module_server("abundance", data = demo_data)
  clustering_module_server("clustering", data = demo_data)
  colocalization_module_server("colocalization", data = demo_data)
  if (environment_diagnostics_enabled()) {
    environment_module_server("environment", app_dir = APP_DIR)
  }
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

facet_column_override <- function(value, facet_count) {
  column_count <- plot_dimension_override(value)
  if (is.null(column_count)) {
    return(NULL)
  }

  max(1L, min(as.integer(column_count), as.integer(max(1L, facet_count))))
}

selected_or_all <- function(selected, all_values) {
  if (is.null(selected) || length(selected) == 0) {
    return(all_values)
  }
  selected
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

shinyApp(ui, server)
