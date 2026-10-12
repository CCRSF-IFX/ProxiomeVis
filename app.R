library(shiny)
library(bslib)

source(file.path("R", "precomputed_contract.R"), local = TRUE)
source(file.path("R", "report_plot_sources.R"), local = TRUE)
source(file.path("R", "report_plots.R"), local = TRUE)
source(file.path("R", "proximity_store.R"), local = TRUE)
source(file.path("R", "hybrid_analysis.R"), local = TRUE)
source(file.path("R", "hybrid_jobs.R"), local = TRUE)
source(file.path("R", "hybrid_module.R"), local = TRUE)
source(file.path("R", "cell_graph_assets.R"), local = TRUE)
source(file.path("R", "cell_graph_module.R"), local = TRUE)
source(file.path("R", "precomputed_module.R"), local = TRUE)

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
    success = "#18765E", danger = "#B73225", "link-color" = "#245C72"),
  title = "ProxiomeVis — Explore and compare",
  tags$head(tags$style(HTML("body {max-width: 1700px; margin: auto; padding: 20px;} .dataTables_wrapper {overflow-x:auto;} #report-provenance_text {max-height: 22rem; overflow: auto; white-space: pre-wrap; overflow-wrap: anywhere;} .nav {flex-wrap: wrap;} @media (max-width: 576px) {body {padding: 10px;}}"))),
  tags$h2("ProxiomeVis"),
  tags$p("Explore prepared results, run bounded exploratory comparisons, and inspect individual cell graphs."),
  accordion(accordion_panel(
    "Dataset",
    if (allow_local_load) tagList(
      textInput("rds_path", "Precomputed RDS path on this server", value = default_path, width = "100%"),
      actionButton("load_rds", "Load results", class = "btn-primary")
    ),
    tags$p(class = "text-muted mt-2", "Use a trusted application bundle. Prepared views need no PXL files; optional sidecars enable spatial tests and cell-level 3D."),
    uiOutput("load_status")
  ), open = if (is.null(initial_data)) TRUE else FALSE),
  conditionalPanel("output.has_data", precomputed_ui("report"))
)

server <- function(input, output, session) {
  data <- reactiveVal(initial_data)
  error <- reactiveVal(initial_error)
  loaded_path <- reactiveVal(if (!is.null(initial_data)) default_path else NULL)
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
  output$has_data <- reactive(!is.null(data()))
  outputOptions(output, "has_data", suspendWhenHidden = FALSE)
  output$load_status <- renderUI({
    if (!is.null(error())) return(tags$div(class = "alert alert-danger mt-2", error(),
                                         if (!is.null(data())) tags$p("The previously loaded dataset remains active.")))
    if (is.null(data())) return(tags$p("Configure PROXIOME_PRECOMPUTED_RDS with a trusted ProxiomeVis application bundle or legacy prepared RDS to begin."))
    tags$div(class = "alert alert-success mt-2",
             paste(data()$contract$title, "—", format(nrow(data()$metadata), big.mark = ","), "cells"),
             tags$div(class = "small", loaded_path()))
  })
  precomputed_server("report", data, job_manager, analysis_code_dir)
}

shinyApp(ui, server)
