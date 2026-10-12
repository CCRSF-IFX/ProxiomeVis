hybrid_analysis_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "New exploratory comparison", width = 310,
      shiny::selectInput(ns("readout"), "Readout", choices = PRECOMPUTED_READOUTS),
      shiny::selectInput(ns("donor"), "One donor", choices = character()),
      shiny::selectInput(ns("celltype"), "One cell type", choices = character()),
      shiny::selectInput(ns("target"), "Target condition", choices = character()),
      shiny::selectInput(ns("reference"), "Reference condition", choices = character()),
      shiny::selectizeInput(ns("markers"), "Markers to test", choices = character(), multiple = TRUE),
      shiny::selectizeInput(ns("cells"), "Cell subset (empty = all in scope)",
                           choices = character(), multiple = TRUE),
      shiny::actionButton(ns("clear_selection"), "Use all cells in scope", class = "btn-sm"),
      shiny::textOutput(ns("selection_summary")),
      shiny::checkboxInput(ns("acknowledge"), "I understand this is exploratory within-donor analysis.", FALSE),
      bslib::input_task_button(ns("run"), "Run comparison"),
      shiny::actionButton(ns("cancel"), "Cancel running analysis", class = "btn-sm btn-outline-secondary"),
      shiny::tags$p(class = "small text-muted",
        "Limits: 2,000 cells, 160 markers, 1,000 pairs, 1 million score rows, 120 seconds. Minimum cells per group: abundance 3, self-proximity 20, colocalization 40."),
      shiny::tags$div(class = "alert alert-warning small",
        "Cells from one donor are not independent biological replicates. Selecting cells using the same marker or embedding can bias tests. These results do not establish population-level significance.")
    ),
    shiny::uiOutput(ns("availability")),
    shiny::tags$div(role = "status", `aria-live` = "polite", shiny::textOutput(ns("status"))),
    bslib::card(full_screen = TRUE,
      bslib::card_header("Select cells in the prepared embedding"),
      bslib::card_body(fill = FALSE,
        "Drag a lasso or box to select cells; condition determines the two comparison groups. Nothing is tested until you click Run comparison."),
      bslib::card_body(plotly::plotlyOutput(ns("cells_plot"), height = "430px"))
    ),
    bslib::card(full_screen = TRUE,
      bslib::card_header("New comparison — not the frozen report"),
      bslib::card_body(fill = FALSE, shiny::textOutput(ns("result_scope")),
        precomputed_download_group(
          shiny::downloadButton(ns("results_csv"), "Results CSV"),
          shiny::downloadButton(ns("provenance_json"), "Analysis provenance JSON")
        )),
      bslib::card_body(plotly::plotlyOutput(ns("volcano"), height = "500px"))
    ),
    bslib::card(full_screen = TRUE,
      bslib::card_header("All requested features and support counts"),
      DT::DTOutput(ns("results_table"))
    )
  )
}

plot_hybrid_volcano <- function(result) {
  rows <- result$results
  rows <- rows[is.finite(rows$effect_size) & is.finite(rows$p_adj), , drop = FALSE]
  if (!nrow(rows)) return(plotly::ggplotly(precomputed_empty_plot(
    "No tested features met the required support. See the results table.", "Exploratory comparison")))
  rows$significance <- -log10(pmax(rows$p_adj, .Machine$double.xmin))
  rows$hover <- paste0("Feature: ", htmltools::htmlEscape(rows$feature_id),
    "<br>Effect (target − reference): ", signif(rows$effect_size, 5),
    "<br>P value: ", signif(rows$p_value, 5), "<br>Adjusted P value: ", signif(rows$p_adj, 5),
    "<br>Target cells: ", rows$n_cells_a, "<br>Reference cells: ", rows$n_cells_b)
  plotly::plot_ly(rows, x = ~effect_size, y = ~significance, type = "scatter", mode = "markers",
    text = ~hover, hoverinfo = "text", marker = list(color = "#3B6EA8", size = 7)) |>
    plotly::layout(xaxis = list(title = "Effect (target − reference)"),
      yaxis = list(title = "-log10(adjusted P)"), dragmode = "zoom", hovermode = "closest") |>
    plotly::config(scrollZoom = TRUE, displayModeBar = TRUE, displaylogo = FALSE, doubleClick = "reset")
}

hybrid_analysis_server <- function(id, data, manager, code_dir, active = shiny::reactive(TRUE)) {
  shiny::moduleServer(id, function(input, output, session) {
    job <- NULL
    message <- shiny::reactiveVal("Choose a scope and markers, then run a comparison.")
    validation_error <- shiny::reactiveVal(NULL)
    last_request <- shiny::reactiveVal(NULL)
    task <- shiny::ExtendedTask$new(function(path, request) {
      job <<- start_hybrid_job(path, request, code_dir, manager)
      job$promise
    }) |> bslib::bind_task_button("run")
    session$onSessionEnded(function() if (!is.null(job)) job$cancel())

    shiny::observeEvent(data(), {
      value <- data()
      if (!is.null(job)) job$cancel()
      last_request(NULL)
      shiny::updateSelectInput(session, "donor", choices = sort(unique(value$metadata$donor)))
      markers <- rownames(value$abundance)
      preferred <- intersect(c("CD3e", "CD4", "CD8", "CD19", "CD14", "CD28", "CD69", "HLA-DR"), markers)
      if (!length(preferred)) preferred <- head(markers, 6L)
      shiny::updateSelectizeInput(session, "markers", choices = markers, selected = preferred, server = TRUE)
    })
    shiny::observeEvent(list(data(), input$donor), {
      meta <- data()$metadata
      choices <- sort(unique(meta$celltype[meta$donor == input$donor]))
      shiny::updateSelectInput(session, "celltype", choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$celltype), choices, "CD4 T"))
    })
    scope <- shiny::reactive({
      shiny::req(input$donor, input$celltype)
      meta <- data()$metadata
      meta[meta$donor == input$donor & meta$celltype == input$celltype, , drop = FALSE]
    })
    shiny::observeEvent(scope(), {
      conditions <- sort(unique(scope()$condition))
      reference <- if ("unstim" %in% conditions) "unstim" else head(conditions, 1L)
      target <- head(setdiff(conditions, reference), 1L)
      shiny::updateSelectInput(session, "reference", choices = conditions, selected = reference)
      shiny::updateSelectInput(session, "target", choices = conditions, selected = target)
      shiny::updateSelectizeInput(session, "cells", choices = scope()$component, selected = character(), server = TRUE)
    })
    shiny::observeEvent(input$clear_selection, {
      shiny::updateSelectizeInput(session, "cells", selected = character())
    })
    output$availability <- shiny::renderUI({
      value <- data()
      if (is.null(value$bundle_path)) return(shiny::tags$p(class = "alert alert-info", "Save a prepared application bundle to enable background analysis."))
      if (!identical(input$readout, "abundance") && is.null(value$assets$proximity)) {
        shiny::tags$p(class = "alert alert-info", "Spatial reanalysis needs the prepared proximity sidecar. Frozen report views remain available without it.")
      }
    })
    output$cells_plot <- plotly::renderPlotly({
      meta <- scope()
      embedding <- data()$embeddings[[data()$contract$defaults$reduction]][meta$component, , drop = FALSE]
      rows <- data.frame(meta, x = embedding[, 1L], y = embedding[, 2L])
      rows$hover <- paste0("Cell: ", htmltools::htmlEscape(rows$component), "<br>Condition: ", htmltools::htmlEscape(rows$condition))
      plotly::plot_ly(rows, x = ~x, y = ~y, color = ~condition, key = ~component,
        source = session$ns("cells"), type = "scatter", mode = "markers", text = ~hover,
        hoverinfo = "text", marker = list(size = 6)) |>
        plotly::layout(dragmode = "lasso", xaxis = list(title = colnames(embedding)[1L]),
                       yaxis = list(title = colnames(embedding)[2L])) |>
        plotly::config(displaylogo = FALSE) |>
        plotly::event_register("plotly_selected")
    })
    selection_event <- shiny::reactive({
      shiny::req(active())
      scope()
      plotly::event_data("plotly_selected", source = session$ns("cells"))
    })
    shiny::observeEvent(selection_event(), {
      selection <- plotly::event_data("plotly_selected", source = session$ns("cells"))
      ids <- intersect(as.character(selection$key), scope()$component)
      if (length(ids)) shiny::updateSelectizeInput(session, "cells", selected = ids)
    })
    selected_metadata <- shiny::reactive({
      meta <- scope()
      if (length(input$cells)) meta <- meta[meta$component %in% input$cells, , drop = FALSE]
      meta
    })
    output$selection_summary <- shiny::renderText({
      meta <- selected_metadata()
      paste(sum(meta$condition == input$target), "target cells;",
            sum(meta$condition == input$reference), "reference cells.")
    })
    shiny::observeEvent(input$run, {
      if (identical(task$status(), "running")) return()
      validation_error(NULL)
      tryCatch({
        if (!isTRUE(input$acknowledge)) stop("Confirm the exploratory-analysis notice before running.", call. = FALSE)
        value <- data()
        if (is.null(value$bundle_path)) stop("No saved application bundle is available.", call. = FALSE)
        if (!length(input$markers)) stop("Choose at least one marker to test.", call. = FALSE)
        if (length(input$cells) && (!is.character(input$cells) || anyNA(input$cells) ||
            anyDuplicated(input$cells) || !all(input$cells %in% scope()$component))) {
          stop("The cell selection is stale or invalid. Select unique cells from the current scope.", call. = FALSE)
        }
        meta <- selected_metadata()
        request <- validate_hybrid_request(value, list(
          dataset_id = value$contract$dataset_id, bundle_md5 = value$bundle_md5, readout = input$readout,
          donor = input$donor, celltype = input$celltype, target = input$target, reference = input$reference,
          markers = input$markers, group_a_cells = meta$component[meta$condition == input$target],
          group_b_cells = meta$component[meta$condition == input$reference]
        ))
        last_request(request)
        message("Running in a separate R process. You can continue exploring the reports.")
        task$invoke(value$bundle_path, request)
      }, error = function(error) validation_error(conditionMessage(error)))
    })
    shiny::observeEvent(input$cancel, {
      if (identical(task$status(), "running") && !is.null(job)) job$cancel()
    })
    output$status <- shiny::renderText({
      if (!is.null(validation_error())) return(validation_error())
      status <- task$status()
      if (identical(status, "error")) return(tryCatch(task$result(), error = function(e) conditionMessage(e)))
      if (identical(status, "success")) return("Completed. The displayed scope is recorded below; changing controls does not rerun the analysis.")
      message()
    })
    result <- shiny::reactive({
      shiny::req(identical(task$status(), "success"))
      value <- task$result()
      shiny::req(last_request(), identical(value$provenance$dataset_id, data()$contract$dataset_id))
      value
    })
    output$result_scope <- shiny::renderText({
      result()
      request <- last_request()
      paste(request$donor, request$celltype, request$readout,
        paste(request$target, "minus", request$reference),
        paste(length(request$group_a_cells), "vs", length(request$group_b_cells), "selected cells"),
        "BH adjustment applies to this requested feature family; frozen report results are unchanged.", sep = " | ")
    })
    output$volcano <- plotly::renderPlotly(plot_hybrid_volcano(result()))
    output$results_table <- DT::renderDT({
      DT::datatable(result()$results, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
    })
    output$results_csv <- shiny::downloadHandler(
      filename = function() "proxiomevis_exploratory_comparison.csv",
      content = function(file) utils::write.csv(result()$results, file, row.names = FALSE)
    )
    output$provenance_json <- shiny::downloadHandler(
      filename = function() "proxiomevis_analysis_provenance.json",
      content = function(file) jsonlite::write_json(list(request = last_request(), provenance = result()$provenance),
        file, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
    )
    invisible(list(task = task, result = result, selected = selected_metadata))
  })
}
