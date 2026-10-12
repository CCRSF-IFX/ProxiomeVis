cell_graph_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(title = "Prepared cell graph", width = 300,
      shiny::selectInput(ns("donor"), "Donor", choices = character()),
      shiny::selectInput(ns("celltype"), "Cell type", choices = character()),
      shiny::selectInput(ns("condition"), "Condition", choices = character()),
      shiny::selectizeInput(ns("pair_1"), "Ranking marker 1", choices = character()),
      shiny::selectizeInput(ns("pair_2"), "Ranking marker 2", choices = character()),
      shiny::selectInput(ns("rank_view"), "Ranked cells to show", choices = c(
        "Top 10" = "top", "Bottom 10" = "bottom", "Top and bottom 10" = "both", "All scored cells" = "all")),
      shiny::selectizeInput(ns("cell"), "Cell (rank / log2_ratio / ID)", choices = character(),
        options = list(sortField = list(list(field = "$order", direction = "asc")))),
      shiny::selectizeInput(ns("markers"), "Highlight marker (one only)", choices = character()),
      shiny::actionButton(ns("display"), "Display", class = "btn-primary"),
      shiny::actionButton(ns("refresh"), "Refresh prepared cells", class = "btn-sm"),
      shiny::tags$p(class = "small text-muted",
        "Select the same ranking marker twice for self-clustering. Pair order does not matter. Clear both ranking markers to browse without ranking."),
      shiny::tags$p(class = "small text-muted",
        "Only Display loads a cell and updates the plot. Changing filters or highlights leaves the last displayed view unchanged. No PXL upload is needed."),
      shiny::tags$div(class = "alert alert-info small",
        "These are abstract graph-layout coordinates, not measured physical positions or distances. Shared-node marker highlights are not a statistical colocalization test.")
    ),
    shiny::singleton(shiny::tags$script(src = "cell-graph-gif.js")),
    shiny::singleton(shiny::tags$script(src = "cell-graph.js")),
    shiny::uiOutput(ns("coverage")),
    shiny::textOutput(ns("ranking_status")),
    shiny::downloadButton(ns("ranking_csv"), "Current ranking CSV"),
    bslib::card(full_screen = TRUE,
      bslib::card_header("Cell-level 3D view"),
      bslib::card_body(fill = FALSE, shiny::textOutput(ns("summary")),
        precomputed_download_group(
          shiny::downloadButton(ns("nodes_csv"), "Node layout and markers CSV"),
          shiny::downloadButton(ns("provenance_json"), "Graph provenance JSON"),
          shiny::tags$button(id = ns("plot_reset"), type = "button", class = "btn btn-outline-primary", "Reset view"),
          shiny::tags$button(id = ns("plot_png"), type = "button", class = "btn btn-outline-primary", "Save PNG")
        ),
        shiny::tags$details(
          shiny::tags$summary("Rotating GIF settings"),
          shiny::tags$div(class = "d-flex flex-wrap gap-3 mt-2",
            shiny::numericInput(ns("plot_gif_seconds"), "Duration (seconds)", 6, min = 1, max = 20, step = 1, width = "160px"),
            shiny::numericInput(ns("plot_gif_fps"), "Frames / second", 15, min = 1, max = 30, step = 1, width = "160px"),
            shiny::numericInput(ns("plot_gif_pixels"), "Width = height (pixels)", 600, min = 128, max = 1024, step = 1, width = "180px")
          ),
          shiny::tags$p(class = "small text-muted",
            "Exports the displayed cell and highlight through one full turn. Maximum 300 frames and 120 million total pixels; reduce duration, FPS, or size for large exports.")
        ),
        shiny::tags$div(class = "d-flex flex-wrap gap-2 align-items-center",
          shiny::tags$button(id = ns("plot_gif"), type = "button", disabled = "disabled",
            class = "btn btn-outline-primary", "Download rotating GIF"),
          shiny::tags$button(id = ns("plot_gif_cancel"), type = "button", hidden = NA,
            class = "btn btn-outline-secondary", "Cancel GIF"),
          shiny::tags$progress(id = ns("plot_gif_progress"), max = 1, value = 0, hidden = NA,
            `aria-label` = "GIF export progress")
        ),
        shiny::tags$p(id = ns("plot_gif_status"), role = "status", `aria-live` = "polite",
          class = "small text-muted", "Display a cell before exporting.")
      ),
      bslib::card_body(
        shiny::tags$p(class = "small text-muted", "Red: highlighted marker. Grey: other nodes. Drag to rotate; scroll to zoom; hover for marker identities."),
        threejs::scatterplotThreeOutput(ns("plot"), height = "700px"))
    )
  )
}

cell_graph_server <- function(id, data, active = shiny::reactive(TRUE), selection = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {
    handoff <- shiny::reactiveVal(NULL)
    handoff_pending <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(selection(), {
      request <- selection()
      if (is.null(request)) {
        handoff(NULL)
        handoff_pending(FALSE)
        return()
      }
      fields <- c("donor", "celltype", "marker_1", "marker_2")
      shiny::req(is.list(request), all(vapply(request[fields], function(value) {
        is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
      }, logical(1))))
      value <- data()
      meta <- value$metadata
      cohort <- meta[meta$donor == request$donor & meta$celltype == request$celltype, , drop = FALSE]
      markers <- rownames(value$abundance)
      shiny::req(nrow(cohort), request$marker_1 %in% markers, request$marker_2 %in% markers,
        is.character(request$conditions), length(request$conditions), !anyNA(request$conditions))
      request$conditions <- intersect(unique(request$conditions), unique(cohort$condition))
      handoff_pending(TRUE)
      handoff(request)
      shiny::updateSelectInput(session, "donor", choices = sort(unique(meta$donor)), selected = request$donor)
      shiny::updateSelectInput(session, "celltype", choices = sort(unique(meta$celltype[meta$donor == request$donor])),
        selected = request$celltype)
      shiny::updateSelectInput(session, "condition", choices = c("Choose a condition" = "",
        stats::setNames(request$conditions, request$conditions)), selected = "")
      choices <- c("Choose a marker" = "", stats::setNames(markers, markers))
      shiny::updateSelectizeInput(session, "pair_1", choices = choices, selected = request$marker_1, server = TRUE)
      shiny::updateSelectizeInput(session, "pair_2", choices = choices, selected = request$marker_2, server = TRUE)
      shiny::updateSelectizeInput(session, "markers", choices = choices, selected = "", server = TRUE)
      shiny::updateSelectizeInput(session, "cell", choices = character(), selected = character(), server = TRUE)
    }, ignoreNULL = FALSE, priority = 100)
    shiny::observeEvent(list(input$donor, input$celltype, input$pair_1, input$pair_2,
                            input$condition, input$markers), {
      request <- handoff()
      if (handoff_pending() && identical(input$donor, request$donor) &&
          identical(input$celltype, request$celltype) && identical(input$pair_1, request$marker_1) &&
          identical(input$pair_2, request$marker_2) &&
          (!length(input$condition) || identical(input$condition, "")) &&
          (!length(input$markers) || identical(input$markers, ""))) handoff_pending(FALSE)
    })
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
    shiny::observeEvent(data(), {
      if (handoff_pending()) return()
      choices <- c("Choose a marker" = "", stats::setNames(rownames(data()$abundance), rownames(data()$abundance)))
      for (id in c("pair_1", "pair_2", "markers")) {
        selected <- shiny::isolate(input[[id]])
        if (length(selected) != 1L || !selected %in% choices) selected <- ""
        shiny::updateSelectizeInput(session, id, choices = choices, selected = selected, server = TRUE)
      }
    })
    ranking_requested <- shiny::reactive({
      any(nzchar(c(input$pair_1, input$pair_2)))
    })
    ranked <- shiny::reactive({
      shiny::req(active(), !handoff_pending(), ranking_requested(), input$donor, input$celltype, input$condition)
      shiny::validate(shiny::need(length(input$pair_1) == 1L && nzchar(input$pair_1) &&
        length(input$pair_2) == 1L && nzchar(input$pair_2), "Select one marker in each ranking dropdown."))
      value <- data()
      shiny::validate(shiny::need(!is.null(value$assets$proximity),
        "Prepare a proximity sidecar to rank cells. Clear both ranking markers to browse prepared graphs."))
      tryCatch(shiny::withProgress(message = "Ranking cells", value = 0.1,
        detail = "Validating sidecars and reading the selected pair…", {
          root <- resolve_hybrid_asset(value$bundle_path, value$assets$proximity)
          result <- rank_cell_proximity(value, root, input$donor, input$celltype, input$condition,
                                        input$pair_1, input$pair_2)
          shiny::incProgress(0.8, detail = "Matching ranked cells to prepared graphs…")
          entries <- index()$index$cells
          result$graph_prepared <- result$component %in% entries$cell_id
          result$n_nodes <- entries$n_nodes[match(result$component, entries$cell_id)]
          result
        }), error = function(error) shiny::validate(shiny::need(FALSE, conditionMessage(error))))
    })
    shown_ranking <- shiny::reactive({
      result <- ranked()
      view <- if (is.null(input$rank_view)) "top" else input$rank_view
      shiny::validate(shiny::need(length(view) == 1L && view %in% c("top", "bottom", "both", "all"),
        "Choose a ranking view."))
      rows <- seq_len(nrow(result))
      chosen <- switch(view, top = head(rows, 10L), bottom = rev(tail(rows, 10L)),
                       both = unique(c(head(rows, 10L), rev(tail(rows, 10L)))), all = rows)
      result[chosen, , drop = FALSE]
    })
    cell_choices <- shiny::reactive({
      shiny::req(!handoff_pending())
      if (ranking_requested()) {
        result <- shown_ranking()
        if (!nrow(result)) return(character())
        labels <- paste0("#", result$rank, " / ", formatC(result$log2_ratio, digits = 4, format = "f"),
          " / ", result$component, ifelse(result$graph_prepared, "", " [graph not prepared]"),
          ifelse(!is.na(result$n_nodes) & result$n_nodes > 200000L, " [over viewer limit]", ""))
        return(stats::setNames(result$component, labels))
      }
      meta <- available()
      choices <- meta$component[meta$donor == input$donor & meta$celltype == input$celltype & meta$condition == input$condition]
      stats::setNames(choices, choices)
    })
    output$ranking_status <- shiny::renderText({
      if (handoff_pending() || (!is.null(handoff()) && !length(input$condition[nzchar(input$condition)]))) {
        return("Choose a condition to rank cells, then select a cell and highlight marker before clicking Display.")
      }
      if (!ranking_requested()) return("Select two ranking markers to rank cells, or browse prepared cells without ranking.")
      result <- ranked()
      total <- attr(result, "provenance")$cohort_cells
      if (!nrow(result)) return("Ranking complete: no finite observed score for this pair in the selected group. Missing scores are not replaced with zero.")
      paste0("Ranking complete: ", nrow(result), " of ", total, " cells have a finite observed score; showing ",
        nrow(shown_ranking()), ". Rank 1 is highest; ties use cell ID. ",
        sum(!result$graph_prepared), " scored cells lack prepared graphs. No abundance or differential-test cutoff is applied.")
    })
    output$ranking_csv <- shiny::downloadHandler(
      filename = function() "proxiomevis_ranked_cells.csv",
      content = function(file) utils::write.csv(shown_ranking(), file, row.names = FALSE)
    )
    output$coverage <- shiny::renderUI({
      value <- index()$index
      shiny::tags$p(class = "alert alert-info small",
        paste(format(nrow(value$cells), big.mark = ","), "of", format(nrow(data()$metadata), big.mark = ","),
          "dataset cells have prepared 3D assets. All nodes are preserved; oversized cells are reported, never silently truncated."))
    })
    shiny::observeEvent(available(), {
      if (handoff_pending()) return()
      meta <- if (is.null(handoff())) available() else data()$metadata
      choices <- sort(unique(meta$donor))
      shiny::updateSelectInput(session, "donor", choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$donor), choices))
    })
    shiny::observeEvent(list(available(), input$donor), {
      if (handoff_pending()) return()
      meta <- if (is.null(handoff())) available() else data()$metadata
      choices <- sort(unique(meta$celltype[meta$donor == input$donor]))
      shiny::updateSelectInput(session, "celltype", choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$celltype), choices, "CD4 T"))
    })
    shiny::observeEvent(list(available(), input$donor, input$celltype), {
      if (handoff_pending()) return()
      meta <- if (is.null(handoff())) available() else data()$metadata
      choices <- sort(unique(meta$condition[meta$donor == input$donor & meta$celltype == input$celltype]))
      request <- handoff()
      if (!is.null(request) && identical(input$donor, request$donor) && identical(input$celltype, request$celltype)) {
        choices <- intersect(request$conditions, choices)
        choices <- c("Choose a condition" = "", stats::setNames(choices, choices))
      }
      shiny::updateSelectInput(session, "condition", choices = choices,
        selected = precomputed_keep_selection(shiny::isolate(input$condition), choices))
    })
    shiny::observe({
      if (handoff_pending()) return()
      choices <- tryCatch(cell_choices(), error = function(error) character())
      shiny::updateSelectizeInput(session, "cell", choices = choices,
        selected = if (ranking_requested()) unname(head(choices, 1L)) else
          precomputed_keep_selection(shiny::isolate(input$cell), unname(choices)), server = TRUE)
    })
    displayed <- shiny::eventReactive(input$display, {
      shiny::req(active())
      shiny::validate(shiny::need(!handoff_pending() && length(input$condition) == 1L && nzchar(input$condition),
        "Choose a condition before clicking Display."))
      choices <- cell_choices()
      shiny::validate(shiny::need(length(input$cell) == 1L && input$cell %in% choices,
        "Select a cell from the current list before clicking Display."))
      shiny::validate(shiny::need(length(input$markers) == 1L && nzchar(input$markers),
        "Select exactly one highlight marker, then click Display."))
      prepared <- index()
      tryCatch(shiny::withProgress(message = "Loading selected cell", value = 0.1, {
        asset <- load_cell_graph_asset(prepared$index, prepared$directory, input$cell)
        shiny::validate(shiny::need(input$markers %in% unlist(asset$nodes$markers, use.names = FALSE),
          "The selected highlight marker has no nodes in this cell. Choose another marker and click Display."))
        shiny::incProgress(0.8, detail = "Preparing the 3D view…")
        ranking <- if (ranking_requested()) {
          list(cell = ranked()[ranked()$component == input$cell, , drop = FALSE],
               provenance = attr(ranked(), "provenance"))
        } else NULL
        list(asset = asset, marker = input$markers, ranking = ranking)
      }),
        error = function(error) shiny::validate(shiny::need(FALSE, conditionMessage(error))))
    }, ignoreInit = TRUE)
    asset <- shiny::reactive(displayed()$asset)
    output$summary <- shiny::renderText({
      value <- asset()
      paste(format(nrow(value$nodes), big.mark = ","), "nodes;", value$layout,
        "— prepared coordinates and marker identities; nodes only. Highlight:", displayed()$marker,
        if (!is.null(displayed()$ranking)) paste0(" | Rank ", displayed()$ranking$cell$rank,
          " | log2_ratio ", signif(displayed()$ranking$cell$log2_ratio, 5L)) else "")
    })
    output$plot <- threejs::renderScatterplotThree({
      tryCatch(htmlwidgets::onRender(plot_cell_graph_asset(asset(), displayed()$marker),
        "function(el, x) { window.proxiomeCellGraph(el, x); }"),
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
        provenance <- value[setdiff(names(value), "nodes")]
        provenance$highlight_marker <- displayed()$marker
        provenance$ranking <- displayed()$ranking
        jsonlite::write_json(provenance, file,
          auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
      }
    )
  })
}
