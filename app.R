library(shiny)
library(bslib)
library(Seurat)

source(file.path("R", "precomputed_contract.R"), local = TRUE)
source(file.path("R", "report_plot_sources.R"), local = TRUE)
source(file.path("R", "report_plots.R"), local = TRUE)
source(file.path("R", "proximity_store.R"), local = TRUE)
source(file.path("R", "hybrid_analysis.R"), local = TRUE)
source(file.path("R", "hybrid_jobs.R"), local = TRUE)
source(file.path("R", "hybrid_module.R"), local = TRUE)
source(file.path("R", "cell_graph_assets.R"), local = TRUE)
source(file.path("R", "cell_graph_export.R"), local = TRUE)
source(file.path("R", "cell_graph_module.R"), local = TRUE)
source(file.path("R", "precomputed_module.R"), local = TRUE)
source(file.path("R", "precomputed_ui.R"), local = TRUE)

demo_bundle <- file.path("output", "hybrid_demo", "analysis.rds")
if (!file.exists(demo_bundle)) demo_bundle <- file.path("output", "proxiomevis_demo_v1.rds")
default_path <- Sys.getenv("PROXIOME_PRECOMPUTED_RDS", unset = demo_bundle)
allow_local_load <- identical(Sys.getenv("PROXIOME_ALLOW_LOCAL_LOAD"), "true")
job_manager <- new_hybrid_job_manager()
analysis_code_dir <- normalizePath("R", mustWork = TRUE)

# A configured dataset is loaded once per R process and shared read-only by sessions.
initial_data <- NULL
initial_error <- NULL
if (file.exists(default_path)) {
  tryCatch(initial_data <- read_precomputed_rds(default_path),
           error = function(error) initial_error <<- conditionMessage(error))
}

ui <- page_fluid(
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#245C72",
    secondary = "#59696B", info = "#246A8E", warning = "#825600",
    success = "#18765E", danger = "#B73225", "link-color" = "#245C72",
    "navbar-light-bg" = "#f8f9fa", "navbar-light-color" = "#245C72",
    "navbar-light-hover-color" = "#212529", "navbar-light-active-color" = "#212529",
    "navbar-light-brand-color" = "#245C72", "navbar-light-brand-hover-color" = "#212529",
    "navbar-light-icon-color" = "#245C72", "navbar-light-toggler-border-color" = "#59696B"),
  title = "ProxiomeVis — Explore and compare",
  tags$head(tags$style(HTML("body {max-width: 1700px; margin: auto; padding: 0 16px 20px;} .dataTables_wrapper {overflow-x:auto;} #report-provenance_text {max-height: 22rem; overflow: auto; white-space: pre-wrap; overflow-wrap: anywhere;} .nav {flex-wrap: wrap;} .dataset-button {max-width: min(20rem, calc(100vw - 5rem)); overflow: hidden; text-overflow: ellipsis; white-space: nowrap;} #load_status {margin: .75rem 0;} @media (max-width: 576px) {body {padding: 0 10px 10px;}}"))),
  precomputed_ui("report",
    dataset_control = actionButton("dataset", "Dataset", class = "btn-outline-primary dataset-button"),
    dataset_status = uiOutput("load_status"))
)

server <- function(input, output, session) {
  data <- reactiveVal(initial_data)
  error <- reactiveVal(initial_error)
  loaded_path <- reactiveVal(if (!is.null(initial_data)) default_path else NULL)
  observeEvent(input$dataset, {
    showModal(modalDialog(
      title = "Dataset", easyClose = TRUE, size = "l",
      if (allow_local_load) tagList(
        textInput("rds_path", "Precomputed RDS path on this server",
          value = if (is.null(loaded_path())) default_path else loaded_path(), width = "100%"),
        actionButton("load_rds", "Load results", class = "btn-primary")
      ),
      tags$p(class = "text-muted mt-2", "Use a trusted application bundle. Prepared views need no PXL files; optional sidecars enable spatial tests and cell-level 3D."),
      uiOutput("dataset_details"),
      footer = modalButton("Close")
    ))
  })
  observeEvent(input$load_rds, {
    req(allow_local_load)
    error(NULL)
    tryCatch({
      result <- withProgress(message = "Loading precomputed results", value = 0.2, {
        read_precomputed_rds(trimws(input$rds_path))
      })
      data(result)
      loaded_path(trimws(input$rds_path))
    }, error = function(e) error(conditionMessage(e)))
  })
  observe({
    title <- if (is.null(data())) "Dataset: not loaded" else paste("Dataset:", data()$contract$title)
    updateActionButton(session, "dataset", label = title)
  })
  output$load_status <- renderUI({
    if (!is.null(error())) return(tags$div(class = "alert alert-danger mt-2", error(),
                                         if (!is.null(data())) tags$p("The previously loaded dataset remains active.")))
    if (is.null(data())) return(tags$p("Configure PROXIOME_PRECOMPUTED_RDS with a trusted ProxiomeVis application bundle or legacy prepared RDS to begin."))
    tags$div(class = "small text-muted", role = "status",
      paste(data()$contract$title, "·", format(nrow(data()$metadata), big.mark = ","), "cells · Prepared results; no PXL upload required"))
  })
  output$dataset_details <- renderUI({
    if (is.null(data())) return(tagList(tags$p("No dataset is loaded."),
      if (!is.null(error())) tags$div(class = "alert alert-danger", error())))
    value <- data()
    tagList(
      tags$h5(value$contract$title),
      tags$p(paste(format(nrow(value$metadata), big.mark = ","), "cells ·", length(unique(value$metadata$donor)), "donor(s)")),
      tags$p(class = "small text-break", loaded_path()),
      tags$ul(
        tags$li(if (is.null(value$assets$proximity)) "Proximity sidecar: not configured" else "Proximity sidecar: configured; validated when used"),
        tags$li(if (is.null(value$assets$graphs)) "3D graph sidecar: not configured" else "3D graph sidecar: configured; validated when used")
      ),
      if (!is.null(error())) tags$div(class = "alert alert-danger", error(), " The previously loaded dataset remains active.")
    )
  })
  precomputed_server("report", data, job_manager, analysis_code_dir)
}

shinyApp(ui, server)
