patch_sidebar <- function(id) {
  ns <- NS(id)

  sidebar(
    title = "Patch controls",
    width = 300,
    accordion(
      open = c("Markers"),
      accordion_panel(
        "Markers",
        selectizeInput(ns("patch_label_filter"), "Marker class", choices = character(0), multiple = TRUE)
      )
    )
  )
}

patch_module_ui <- function(id) {
  ns <- NS(id)

  nav_panel(
    "Patch Analysis",
    layout_sidebar(
      sidebar = patch_sidebar(id),
      navset_card_underline(
        id = ns("patch_mode"),
        title = "Patch Analysis",
        full_screen = TRUE,
        nav_panel(
          "Markers",
          uiOutput(ns("patch_metric_row")),
          plot_pane(
            size = "compact",
            download_id = "patch_marker_unmixing_plot",
            ns = ns,
            plotlyOutput(ns("patch_marker_unmixing_plot"), height = "auto")
          ),
          div(class = "table-pane", tableOutput(ns("patch_marker_table")))
        ),
        nav_panel(
          "Raji Signal",
          plot_pane(
            size = "compact",
            download_id = "patch_raji_proximity_plot",
            ns = ns,
            plotlyOutput(ns("patch_raji_proximity_plot"), height = "auto")
          ),
          div(class = "table-pane", tableOutput(ns("patch_raji_abundance_table")))
        ),
        nav_panel(
          "Patch Burden",
          div(class = "table-pane", tableOutput(ns("patch_burden_table")))
        )
      )
    )
  )
}

patch_module_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    observe({
      current_data <- data()
      labels <- character(0)
      if (!is.null(current_data$patch$marker_unmixing) && "label" %in% names(current_data$patch$marker_unmixing)) {
        labels <- sort(unique(as.character(current_data$patch$marker_unmixing$label)))
      }
      updateSelectizeInput(session, "patch_label_filter", choices = labels, selected = labels, server = TRUE)
    })

    patch_payload <- reactive({
      current_data <- data()
      req(current_data)
      validate(need(!is.null(current_data$patch), "Patch analysis is available for the Raji/CAR-T demo data source."))
      current_data$patch
    })

    filtered_marker_unmixing <- reactive({
      table <- patch_payload()$marker_unmixing
      validate(need(!is.null(table) && nrow(table) > 0, "No patch marker table is available."))
      if (!is.null(input$patch_label_filter) && length(input$patch_label_filter) > 0 && "label" %in% names(table)) {
        table <- table[table$label %in% input$patch_label_filter, , drop = FALSE]
      }
      table
    })

    output$patch_metric_row <- renderUI({
      plan <- patch_payload()$run_plan
      detection <- as.logical(patch_scalar(plan, "run_patch_detection", FALSE))

      metric_row(
        metric_tile("Patch Detection", if (isTRUE(detection)) "Run" else "Skipped"),
        metric_tile("Cells Selected", format(patch_scalar(plan, c("n_cart_cells_selected", "n_cd8t_cells_selected")), big.mark = ",")),
        metric_tile("Patch Markers", format(patch_scalar(plan, "n_patch_markers"), big.mark = ",")),
        metric_tile("Receiver Markers", format(patch_scalar(plan, "n_receiver_markers"), big.mark = ","))
      )
    })

    patch_marker_unmixing_ggplot <- reactive({
      table <- filtered_marker_unmixing()
      required <- c("receiver_freq", "target_freq", "label", "marker")
      validate(need(all(required %in% names(table)), "Patch marker table is missing frequency columns."))

      table$text <- paste0(
        table$marker,
        "<br>Class: ", table$label,
        "<br>Receiver frequency: ", signif(table$receiver_freq, 3),
        "<br>Target frequency: ", signif(table$target_freq, 3)
      )

      ggplot(table, aes(x = receiver_freq, y = target_freq, color = label, text = text)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#7b8588") +
        geom_point(size = 2.4, alpha = 0.85) +
        labs(x = "Receiver frequency", y = "Target frequency", color = "Class") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
    })

    patch_marker_unmixing_dimensions <- reactive({
      plot_options_input_dimensions(input, "patch_marker_unmixing_plot")
    })
    patch_marker_unmixing_export_dimensions <- reactive({
      plot_options_export_dimensions(input, "patch_marker_unmixing_plot")
    })

    output$patch_marker_unmixing_plot <- renderPlotly({
      dimensions <- patch_marker_unmixing_dimensions()
      display_dimensions <- plotly_display_dimensions(dimensions)
      ggplotly(
        patch_marker_unmixing_ggplot(),
        tooltip = "text",
        width = display_dimensions$width,
        height = display_dimensions$height
      ) |>
        apply_proxiome_plot_frame(dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "patch_marker_unmixing_plot",
      patch_marker_unmixing_ggplot,
      filename_prefix = "patch-marker-unmixing",
      width = function() plot_download_size_from_dimensions(patch_marker_unmixing_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(patch_marker_unmixing_export_dimensions())$height
    )

    output$patch_marker_table <- renderTable({
      format_patch_table(filtered_marker_unmixing(), max_rows = 40)
    }, striped = TRUE, bordered = FALSE, width = "100%")

    patch_raji_proximity_ggplot <- reactive({
      table <- patch_payload()$raji_marker_proximity
      validate(need(!is.null(table) && nrow(table) > 0, "No Raji-marker proximity table is available."))
      validate(need(all(c("raji_marker_count", "log2_ratio", "celltype_condition") %in% names(table)), "Raji-marker proximity table is missing required columns."))

      table <- table[is.finite(table$raji_marker_count) & is.finite(table$log2_ratio), , drop = FALSE]
      table <- patch_plot_rows(table)
      table$text <- paste0(
        table$celltype_condition,
        "<br>Raji marker count: ", table$raji_marker_count,
        "<br>log2 ratio: ", signif(table$log2_ratio, 3)
      )

      ggplot(table, aes(x = raji_marker_count, y = log2_ratio, color = celltype_condition, text = text)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "#7b8588") +
        geom_point(size = 1.5, alpha = 0.65) +
        labs(x = "Raji marker count", y = "Joint proximity log2 ratio", color = "Cell state") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
    })

    patch_raji_proximity_dimensions <- reactive({
      plot_options_input_dimensions(input, "patch_raji_proximity_plot")
    })
    patch_raji_proximity_export_dimensions <- reactive({
      plot_options_export_dimensions(input, "patch_raji_proximity_plot")
    })

    output$patch_raji_proximity_plot <- renderPlotly({
      dimensions <- patch_raji_proximity_dimensions()
      display_dimensions <- plotly_display_dimensions(dimensions)
      ggplotly(
        patch_raji_proximity_ggplot(),
        tooltip = "text",
        width = display_dimensions$width,
        height = display_dimensions$height
      ) |>
        apply_proxiome_plot_frame(dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "patch_raji_proximity_plot",
      patch_raji_proximity_ggplot,
      filename_prefix = "raji-marker-proximity",
      width = function() plot_download_size_from_dimensions(patch_raji_proximity_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(patch_raji_proximity_export_dimensions())$height
    )

    output$patch_raji_abundance_table <- renderTable({
      abundance <- patch_payload()$raji_marker_abundance
      validate(need(!is.null(abundance) && nrow(abundance) > 0, "No Raji-marker abundance table is available."))
      format_patch_table(abundance, max_rows = 40)
    }, striped = TRUE, bordered = FALSE, width = "100%")

    output$patch_burden_table <- renderTable({
      burden <- patch_payload()$patch_burden
      validate(need(!is.null(burden) && nrow(burden) > 0, "Graph-level patch detection was not run for this demo."))
      format_patch_table(burden, max_rows = 40)
    }, striped = TRUE, bordered = FALSE, width = "100%")
  })
}

patch_scalar <- function(table, columns, default = NA_integer_) {
  if (is.null(table) || !is.data.frame(table) || nrow(table) == 0) {
    return(default)
  }
  column <- columns[columns %in% names(table)][1]
  if (is.na(column)) {
    return(default)
  }

  table[[column]][1] %||% default
}

patch_plot_rows <- function(table, max_rows = 4000L) {
  if (nrow(table) <= max_rows) {
    return(table)
  }

  table[unique(round(seq(1, nrow(table), length.out = max_rows))), , drop = FALSE]
}

format_patch_table <- function(table, max_rows = 40L) {
  table <- head(as.data.frame(table), max_rows)
  for (column in names(table)) {
    if (is.numeric(table[[column]])) {
      table[[column]] <- round(table[[column]], 4)
    }
  }

  table
}
