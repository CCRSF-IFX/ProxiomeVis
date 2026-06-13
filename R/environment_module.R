environment_runtime_paths <- function(app_dir) {
  runtime_rows <- data.frame(
    item = c(
      "R",
      "Rscript",
      "R.home",
      "R version",
      "APP_DIR",
      "PROXIOMEVIS_HOME",
      "PROXIOME_RUNTIME_DIR",
      "RENV_PROJECT",
      "RENV_PATHS_LIBRARY",
      "R_LIBS",
      "R_LIBS_USER",
      "R_LIBS_SITE"
    ),
    value = c(
      unname(Sys.which("R")),
      unname(Sys.which("Rscript")),
      R.home(),
      R.version.string,
      normalizePath(app_dir, mustWork = FALSE),
      Sys.getenv("PROXIOMEVIS_HOME", unset = ""),
      Sys.getenv("PROXIOME_RUNTIME_DIR", unset = ""),
      Sys.getenv("RENV_PROJECT", unset = ""),
      Sys.getenv("RENV_PATHS_LIBRARY", unset = ""),
      Sys.getenv("R_LIBS", unset = ""),
      Sys.getenv("R_LIBS_USER", unset = ""),
      Sys.getenv("R_LIBS_SITE", unset = "")
    ),
    stringsAsFactors = FALSE
  )

  lib_rows <- data.frame(
    item = rep(".libPaths()", length(.libPaths())),
    value = .libPaths(),
    stringsAsFactors = FALSE
  )

  resolved_libs <- resolve_app_r_libs(app_dir)
  resolved_rows <- data.frame(
    item = rep("resolved app R library", length(resolved_libs)),
    value = resolved_libs,
    stringsAsFactors = FALSE
  )

  rows <- rbind(runtime_rows, lib_rows, resolved_rows)
  rows[nzchar(rows$value), , drop = FALSE]
}

environment_package_paths <- function(package_filter = character(0), app_dir = NULL) {
  lib_paths <- .libPaths()
  if (!is.null(app_dir)) {
    lib_paths <- normalize_existing_paths(c(resolve_app_r_libs(app_dir), lib_paths))
  }

  packages <- as.data.frame(
    installed.packages(lib.loc = lib_paths)[, c("Package", "Version", "LibPath"), drop = FALSE],
    stringsAsFactors = FALSE
  )
  rownames(packages) <- NULL

  if (length(package_filter) > 0) {
    packages <- packages[packages$Package %in% package_filter, , drop = FALSE]
  }

  packages[order(packages$Package, packages$LibPath), , drop = FALSE]
}

environment_module_ui <- function(id) {
  ns <- NS(id)

  layout_sidebar(
    sidebar = sidebar(
      title = "Environment",
      width = 300,
      actionButton(ns("refresh"), "Refresh", class = "btn btn-outline-secondary btn-sm"),
      selectizeInput(ns("package_filter"), "Package", choices = character(0), multiple = TRUE)
    ),
    navset_card_underline(
      id = ns("environment_view"),
      title = "Environment",
      full_screen = TRUE,
      nav_panel(
        "R Paths",
        div(class = "table-pane", tableOutput(ns("r_paths")))
      ),
      nav_panel(
        "Packages",
        div(class = "table-pane", tableOutput(ns("package_paths")))
      )
    )
  )
}

environment_module_server <- function(id, app_dir) {
  moduleServer(id, function(input, output, session) {
    runtime_paths <- reactiveVal(environment_runtime_paths(app_dir))
    package_paths <- reactiveVal(environment_package_paths(app_dir = app_dir))

    observeEvent(input$refresh, {
      runtime_paths(environment_runtime_paths(app_dir))
      package_paths(environment_package_paths(app_dir = app_dir))
    }, ignoreInit = TRUE)

    observe({
      packages <- sort(unique(package_paths()$Package))
      updateSelectizeInput(session, "package_filter", choices = packages, selected = character(0))
    })

    output$r_paths <- renderTable({
      runtime_paths()
    }, striped = TRUE, bordered = FALSE, width = "100%")

    output$package_paths <- renderTable({
      environment_package_paths(input$package_filter, app_dir = app_dir)
    }, striped = TRUE, bordered = FALSE, width = "100%")
  })
}
