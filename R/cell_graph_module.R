cell_graph_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(title = "Prepared cell graph", width = 300,
      shiny::selectInput(ns("donor"), "Donor", choices = character()),
      shiny::selectInput(ns("celltype"), "Cell type", choices = character()),
      shiny::selectInput(ns("condition"), "Condition", choices = character()),
      shiny::selectizeInput(ns("cell"), "Cell", choices = character()),
      shiny::selectizeInput(ns("markers"), "Highlight up to three markers", choices = character(), multiple = TRUE,
        options = list(maxItems = 3)),
      shiny::actionButton(ns("refresh"), "Refresh prepared cells", class = "btn-sm"),
      shiny::tags$p(class = "small text-muted",
        "Only the chosen cell is loaded. No PXL upload is needed. Rotate by dragging; scroll to zoom. Unhighlighted nodes remain visible."),
      shiny::tags$div(class = "alert alert-info small",
        "These are abstract graph-layout coordinates, not measured physical positions or distances. Shared-node marker highlights are not a statistical colocalization test.")
    ),
    shiny::uiOutput(ns("coverage")),
    bslib::card(full_screen = TRUE,
      bslib::card_header("Cell-level 3D view"),
      bslib::card_body(fill = FALSE, shiny::textOutput(ns("summary")),
        precomputed_download_group(
          shiny::downloadButton(ns("nodes_csv"), "Node layout and markers CSV"),
          shiny::downloadButton(ns("provenance_json"), "Graph provenance JSON")
        )),
      bslib::card_body(plotly::plotlyOutput(ns("plot"), height = "700px"))
    )
  )
}

cell_graph_server <- function(id, data, active = shiny::reactive(TRUE)) {
  shiny::moduleServer(id, function(input, output, session) {
    index <- shiny::reactive({
      shiny::req(active())
      input$refresh
      value <- data()
      shiny::validate(shiny::need(!is.null(value$bundle_path) && !is.null(value$assets$graphs),
        "This dataset has no prepared 3D assets. Prepare the graph sidecar to enable this view."))
      tryCatch({
        directory <- resolve_hybrid_asset(value$bundle_path, value$assets$graphs)
        prepared <- read_cell_graph_index(directory, value$contract$dataset_id, value$metadata$component)
        if (!is.null(value$cell_sources)) {
          sources <- value$cell_sources[match(prepared$cells$cell_id, value$cell_sources$component), , drop = FALSE]
          if (!identical(prepared$cells$source_cell_id, sources$original_id) ||
              !identical(prepared$cells$sample, sources$sample_id)) {
            stop("The graph index does not match the bundle's source cell mapping.", call. = FALSE)
          }
        }
        inventory <- value$contract$proximity_source$pxl_manifest
        if (!is.null(inventory)) {
          if (!is.data.frame(inventory) || !all(c("sample_alias", "pxl_md5") %in% names(inventory)) ||
              anyNA(inventory$sample_alias) || anyDuplicated(inventory$sample_alias) ||
              any(!nzchar(as.character(inventory$sample_alias)))) {
            stop("The bundle's frozen PXL source inventory is invalid.", call. = FALSE)
          }
          expected <- as.character(inventory$pxl_md5[match(prepared$samples$sample, inventory$sample_alias)])
          if (anyNA(expected) || !identical(prepared$samples$source_pxl_md5, expected)) {
            stop("The graph index does not match the bundle's frozen PXL source checksums.", call. = FALSE)
          }
        }
        list(directory = directory, index = prepared)
      }, error = function(error) shiny::validate(shiny::need(FALSE, conditionMessage(error))))
    })
    available <- shiny::reactive({
      entries <- index()$index$cells
      meta <- data()$metadata
      meta[match(entries$cell_id, meta$component), , drop = FALSE]
    })
    output$coverage <- shiny::renderUI({
      value <- index()$index
      shiny::tags$p(class = "alert alert-info small",
        paste(format(nrow(value$cells), big.mark = ","), "of", format(nrow(data()$metadata), big.mark = ","),
          "dataset cells have prepared 3D assets. All nodes are preserved; oversized cells are reported, never silently truncated."))
    })
    shiny::observeEvent(available(), {
      choices <- sort(unique(available()$donor))
      shiny::updateSelectInput(session, "donor", choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$donor), choices))
    })
    shiny::observeEvent(list(available(), input$donor), {
      meta <- available()
      choices <- sort(unique(meta$celltype[meta$donor == input$donor]))
      shiny::updateSelectInput(session, "celltype", choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$celltype), choices, "CD4 T"))
    })
    shiny::observeEvent(list(available(), input$donor, input$celltype), {
      meta <- available()
      choices <- sort(unique(meta$condition[meta$donor == input$donor & meta$celltype == input$celltype]))
      shiny::updateSelectInput(session, "condition", choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$condition), choices))
    })
    shiny::observeEvent(list(available(), input$donor, input$celltype, input$condition), {
      meta <- available()
      choices <- meta$component[meta$donor == input$donor & meta$celltype == input$celltype & meta$condition == input$condition]
      shiny::updateSelectizeInput(session, "cell", choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$cell), choices), server = TRUE)
    })
    asset <- shiny::reactive({
      shiny::req(input$cell)
      prepared <- index()
      tryCatch(load_cell_graph_asset(prepared$index, prepared$directory, input$cell),
        error = function(error) shiny::validate(shiny::need(FALSE, conditionMessage(error))))
    })
    shiny::observeEvent(asset(), {
      choices <- sort(unique(as.character(unlist(asset()$nodes$markers, use.names = FALSE))))
      selected <- intersect(shiny::isolate(input$markers), choices)
      if (!length(selected)) selected <- head(intersect(c("CD3e", "CD4", "CD8", "CD19", "CD14"), choices), 2L)
      shiny::updateSelectizeInput(session, "markers", choices = choices, selected = selected, server = TRUE)
    })
    output$summary <- shiny::renderText({
      value <- asset()
      paste(format(nrow(value$nodes), big.mark = ","), "nodes;", value$layout,
        "— prepared coordinates and marker identities; nodes only.")
    })
    output$plot <- plotly::renderPlotly({
      tryCatch(plot_cell_graph_asset(asset(), input$markers),
        error = function(error) shiny::validate(shiny::need(FALSE, conditionMessage(error))))
    })
    output$nodes_csv <- shiny::downloadHandler(
      filename = function() "proxiomevis_cell_coordinates.csv",
      content = function(file) {
        nodes <- asset()$nodes
        nodes$markers <- vapply(nodes$markers, function(markers) {
          as.character(jsonlite::toJSON(markers, auto_unbox = FALSE))
        }, character(1))
        utils::write.csv(nodes, file, row.names = FALSE)
      }
    )
    output$provenance_json <- shiny::downloadHandler(
      filename = function() "proxiomevis_cell_graph_provenance.json",
      content = function(file) {
        value <- asset()
        jsonlite::write_json(value[setdiff(names(value), "nodes")], file,
          auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
      }
    )
  })
}
