abundance_sidebar <- function(id) {
  ns <- NS(id)

  sidebar(
    title = "Abundance controls",
    width = 300,
    conditionalPanel(
      condition = "input.abundance_mode == 'Observed'",
      ns = ns,
      accordion(
        open = c("Display", "Filters"),
        accordion_panel(
          "Display",
          selectInput(ns("abundance_embedding"), "Embedding", choices = character(0)),
          selectInput(ns("abundance_color_by"), "Color UMAP by", choices = character(0)),
          conditionalPanel(
            condition = "input.abundance_color_by == 'abundance'",
            ns = ns,
            selectInput(ns("abundance_marker"), "Marker", choices = character(0))
          ),
          selectInput(ns("abundance_split_by"), "Split UMAP by", choices = character(0)),
          sliderInput(ns("abundance_point_size"), "Dot size", min = 0.5, max = 5, value = 0.6, step = 0.1)
        ),
        accordion_panel(
          "Filters",
          selectizeInput(ns("abundance_condition_filter"), "Condition", choices = character(0), multiple = TRUE),
          selectizeInput(ns("abundance_celltype_filter"), "Cell type", choices = character(0), multiple = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.abundance_mode == 'Marker Distributions'",
      ns = ns,
      accordion(
        open = c("Display"),
        accordion_panel(
          "Display",
          selectInput(ns("abundance_distribution_marker"), "Marker", choices = character(0)),
          numericInput(ns("abundance_distribution_columns"), "Facet columns", value = 3, min = 1, max = 12, step = 1),
          checkboxInput(ns("abundance_distribution_show_jitter"), "Show jitter dots", value = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.abundance_mode == 'Differential'",
      ns = ns,
      accordion(
        open = c("Contrast", "Thresholds"),
        accordion_panel(
          "Contrast",
          selectInput(ns("abundance_diff_group_a"), "Group A", choices = character(0)),
          selectInput(ns("abundance_diff_group_b"), "Group B (reference)", choices = character(0)),
          selectizeInput(ns("abundance_diff_celltype_filter"), "Cell type", choices = character(0), multiple = TRUE),
          checkboxInput(ns("abundance_diff_stratify_celltype"), "Stratify by cell type", value = FALSE),
          actionButton(ns("abundance_run_differential"), "Run differential analysis", class = "btn-primary w-100")
        ),
        accordion_panel(
          "Thresholds",
          numericInput(ns("abundance_diff_fdr"), "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput(ns("abundance_diff_effect"), "Minimum effect", value = 0.25, min = 0, step = 0.05),
          numericInput(ns("abundance_diff_min_cells"), "Minimum cells per group", value = 3, min = 1, step = 1)
        ),
        accordion_panel(
          "Detail",
          selectInput(ns("abundance_diff_marker"), "Detail marker", choices = character(0))
        )
      )
    )
  )
}

abundance_module_ui <- function(id) {
  ns <- NS(id)

  nav_panel(
    "Abundance",
    layout_sidebar(
      sidebar = abundance_sidebar(id),
      navset_card_underline(
        id = ns("abundance_mode"),
        title = "Abundance",
        full_screen = TRUE,
        nav_panel(
          "Observed",
          uiOutput(ns("metric_row")),
          plot_pane(
            size = "standard",
            extra_class = "umap-plot-pane",
            download_id = "abundance_umap",
            ns = ns,
            controls = plot_resize_controls(ns, "abundance_umap_width", "abundance_umap_height", width_value = 832, height_value = 520),
            uiOutput(ns("abundance_umap_ui"))
          ),
          conditionalPanel(
            condition = "input.abundance_color_by == 'abundance'",
            ns = ns,
            div(class = "table-pane", tableOutput(ns("abundance_table")))
          )
        ),
        nav_panel(
          "Marker Distributions",
          plot_pane(
            size = "scroll",
            extra_class = "distribution-plot-pane",
            download_id = "abundance_marker_distribution_plot",
            ns = ns,
            controls = plot_resize_controls(ns, "abundance_distribution_width", "abundance_distribution_height", width_value = 832, height_value = 678),
            uiOutput(ns("abundance_marker_distribution_plot_ui"))
          ),
          div(class = "table-pane", tableOutput(ns("abundance_marker_distribution_table")))
        ),
        nav_panel(
          "Cell Annotation",
          plot_pane(
            size = "compact",
            download_id = "abundance_celltype_composition_plot",
            ns = ns,
            plotlyOutput(ns("abundance_celltype_composition_plot"), height = proxiome_plot_height())
          ),
          plot_pane(
            size = "wide",
            download_id = "abundance_annotation_heatmap",
            ns = ns,
            plotlyOutput(ns("abundance_annotation_heatmap"), height = proxiome_plot_height())
          ),
          div(class = "table-pane", tableOutput(ns("abundance_celltype_composition_table")))
        ),
        nav_panel(
          "Differential",
          uiOutput(ns("abundance_diff_summary")),
          differential_plot_row(ns("abundance_diff_volcano"), ns("abundance_diff_detail")),
          div(class = "table-pane", tableOutput(ns("abundance_diff_table")))
        )
      )
    )
  )
}

abundance_module_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    abundance_diff_config <- reactiveVal(NULL)

    abundance_differential_config_from_inputs <- function(current_data) {
      make_differential_config(
        group_a = input$abundance_diff_group_a,
        group_b = input$abundance_diff_group_b,
        celltype_filter = selected_or_all(
          input$abundance_diff_celltype_filter,
          unique(current_data$metadata$celltype_manual)
        ),
        stratify_by_celltype = input$abundance_diff_stratify_celltype,
        min_cells = numeric_input_value(input$abundance_diff_min_cells, 3),
        fdr_cutoff = numeric_input_value(input$abundance_diff_fdr, 0.05),
        effect_cutoff = numeric_input_value(input$abundance_diff_effect, 0.25)
      )
    }

    observe({
      current_data <- data()
      req(current_data)

      embedding_choices <- available_embedding_choices(current_data$metadata)
      conditions <- sort(unique(current_data$metadata$condition))
      cell_types <- sort(unique(current_data$metadata$celltype_manual))
      default_group_a <- conditions[1]
      default_group_b <- conditions[min(2, length(conditions))]

      updateSelectInput(session, "abundance_embedding", choices = embedding_choices, selected = embedding_choices[1])
      updateSelectInput(
        session,
        "abundance_color_by",
        choices = available_abundance_color_choices(current_data$metadata),
        selected = "abundance"
      )
      updateSelectInput(session, "abundance_marker", choices = current_data$marker_options, selected = current_data$marker_options[1])
      updateSelectInput(
        session,
        "abundance_distribution_marker",
        choices = current_data$marker_options,
        selected = current_data$marker_options[1]
      )
      updateSelectInput(
        session,
        "abundance_split_by",
        choices = available_abundance_split_choices(current_data$metadata),
        selected = ""
      )
      updateSelectizeInput(session, "abundance_condition_filter", choices = conditions, selected = conditions)
      updateSelectizeInput(session, "abundance_celltype_filter", choices = cell_types, selected = cell_types)
      updateSelectInput(session, "abundance_diff_group_a", choices = conditions, selected = default_group_a)
      updateSelectInput(session, "abundance_diff_group_b", choices = conditions, selected = default_group_b)
      updateSelectizeInput(session, "abundance_diff_celltype_filter", choices = cell_types, selected = cell_types)
      updateSelectInput(session, "abundance_diff_marker", choices = current_data$marker_options, selected = current_data$marker_options[1])

      abundance_diff_config(default_differential_config(conditions, cell_types))
    })

    observeEvent(input$abundance_run_differential, {
      current_data <- data()
      req(current_data)
      abundance_diff_config(abundance_differential_config_from_inputs(current_data))
    }, ignoreInit = TRUE)

    filtered_metadata_for <- function(condition_filter, celltype_filter) {
      current_data <- data()
      req(current_data)

      metadata <- current_data$metadata
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

    abundance_readout_with_metadata <- reactive({
      current_data <- data()
      req(current_data)

      merge(
        current_data$abundance,
        current_data$metadata[, intersect(c("component", "condition", "celltype_manual"), names(current_data$metadata)), drop = FALSE],
        by = "component",
        all.x = FALSE,
        sort = FALSE
      )
    })

    abundance_diff_results <- reactive({
      current_data <- data()
      config <- abundance_diff_config()
      req(current_data, config, config$group_a, config$group_b)

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

    output$metric_row <- renderUI({
      current_data <- data()
      req(current_data)
      metadata <- abundance_metadata()

      metric_row(
        metric_tile("Cells", format(nrow(metadata), big.mark = ",")),
        metric_tile("Markers", format(length(current_data$marker_options), big.mark = ",")),
        metric_tile("Conditions", format(length(unique(metadata$condition)), big.mark = ",")),
        metric_tile("Cell Types", format(length(unique(metadata$celltype_manual)), big.mark = ","))
      )
    })

    abundance_points <- reactive({
      current_data <- data()
      req(current_data, input$abundance_embedding)

      metadata <- abundance_metadata()
      color_by <- selected_abundance_color_by(input$abundance_color_by, metadata)

      if (identical(color_by, "abundance")) {
        req(input$abundance_marker)
        abundance <- current_data$abundance[current_data$abundance$marker == input$abundance_marker, , drop = FALSE]
        plot_data <- merge(metadata, abundance, by = "component", all.x = FALSE, sort = FALSE)
      } else {
        plot_data <- metadata
      }

      plot_data$abundance_color_by <- color_by
      add_embedding_columns(plot_data, input$abundance_embedding)
    })

    abundance_umap_ggplot <- reactive({
      plot_data <- abundance_points()
      validate(need(nrow(plot_data) > 0, "No cells match the selected filters."))

      color_by <- unique(plot_data$abundance_color_by)
      color_by <- color_by[1] %||% "abundance"
      point_size <- numeric_input_value(input$abundance_point_size, 0.6)

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
      }

      split_col <- selected_split_column(input$abundance_split_by, plot_data)
      if (!is.null(split_col)) {
        plot_data$split_group <- plot_data[[split_col]]
        p <- p %+% plot_data +
          facet_wrap(~split_group)
      }

      p
    })

    abundance_umap_dimensions <- reactive({
      abundance_umap_widget_dimensions(
        width_px = input$abundance_umap_width,
        height_px = input$abundance_umap_height
      )
    })

    output$abundance_umap_ui <- renderUI({
      dimensions <- abundance_umap_dimensions()
      div(
        class = "umap-plot-shell",
        style = paste0("width:", dimensions$width, "px;height:", dimensions$height, "px;"),
        plotlyOutput(session$ns("abundance_umap"), width = "100%", height = "100%")
      )
    })

    output$abundance_umap <- renderPlotly({
      plot_data <- abundance_points()
      color_by <- unique(plot_data$abundance_color_by)
      color_by <- color_by[1] %||% "abundance"
      colorbar_title <- if (identical(color_by, "abundance")) paste(input$abundance_marker, "abundance") else NULL
      dimensions <- abundance_umap_dimensions()

      ggplotly(
        abundance_umap_ggplot(),
        tooltip = "text",
        width = dimensions$width,
        height = dimensions$height
      ) |>
        apply_proxiome_plot_frame(colorbar_title = colorbar_title)
    })
    register_ggplot_downloads(
      output,
      "abundance_umap",
      abundance_umap_ggplot,
      filename_prefix = function() paste("abundance-umap", input$abundance_color_by %||% "abundance", input$abundance_marker %||% "", sep = "-"),
      width = function() plot_download_size_from_dimensions(abundance_umap_dimensions())$width,
      height = function() plot_download_size_from_dimensions(abundance_umap_dimensions())$height
    )

    output$abundance_table <- renderTable({
      current_data <- data()
      req(current_data, input$abundance_marker)

      summary <- current_data$abundance_summary[
        current_data$abundance_summary$marker == input$abundance_marker &
          current_data$abundance_summary$condition %in% selected_or_all(
            input$abundance_condition_filter,
            unique(current_data$metadata$condition)
          ) &
          current_data$abundance_summary$celltype_manual %in% selected_or_all(
            input$abundance_celltype_filter,
            unique(current_data$metadata$celltype_manual)
          ),
        ,
        drop = FALSE
      ]
      format_summary_table(summary, value_label = "mean_abundance")
    }, striped = TRUE, bordered = FALSE, width = "100%")

    abundance_distribution_data <- reactive({
      current_data <- data()
      req(current_data, input$abundance_distribution_marker)

      plot_data <- current_data$abundance[current_data$abundance$marker == input$abundance_distribution_marker, , drop = FALSE]
      merge(
        plot_data,
        current_data$metadata[, intersect(c("component", "sample_alias", "sample", "condition", "celltype_manual"), names(current_data$metadata)), drop = FALSE],
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
        plotlyOutput(session$ns("abundance_marker_distribution_plot"), width = "100%", height = "100%")
      )
    })

    abundance_marker_distribution_ggplot <- reactive({
      plot_data <- abundance_distribution_data()
      validate(need(nrow(plot_data) > 0, "No abundance rows are available for the selected marker."))

      dimensions <- abundance_distribution_dimensions()
      plot_abundance_marker_distribution(
        plot_data,
        input$abundance_distribution_marker,
        facet_cols = dimensions$facet_cols,
        show_jitter = !identical(input$abundance_distribution_show_jitter, FALSE)
      )
    })

    output$abundance_marker_distribution_plot <- renderPlotly({
      dimensions <- abundance_distribution_dimensions()
      ggplotly(
        abundance_marker_distribution_ggplot(),
        tooltip = c("x", "y", "fill"),
        width = dimensions$width,
        height = dimensions$height
      ) |>
        apply_proxiome_plot_frame()
    })
    register_ggplot_downloads(
      output,
      "abundance_marker_distribution_plot",
      abundance_marker_distribution_ggplot,
      filename_prefix = function() paste("abundance-marker-distribution", input$abundance_distribution_marker %||% "marker", sep = "-"),
      width = function() plot_download_size_from_dimensions(abundance_distribution_dimensions())$width,
      height = function() plot_download_size_from_dimensions(abundance_distribution_dimensions())$height
    )

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

    abundance_celltype_composition_ggplot <- reactive({
      current_data <- data()
      req(current_data)
      plot_data <- celltype_composition_data(current_data$metadata)
      validate(need(nrow(plot_data) > 0, "No cell annotation rows are available."))

      plot_celltype_composition(plot_data)
    })

    output$abundance_celltype_composition_plot <- renderPlotly({
      ggplotly(abundance_celltype_composition_ggplot(), tooltip = "text") |>
        apply_proxiome_plot_frame()
    })
    register_ggplot_downloads(
      output,
      "abundance_celltype_composition_plot",
      abundance_celltype_composition_ggplot,
      filename_prefix = "abundance-celltype-composition",
      width = 7,
      height = 5
    )

    abundance_annotation_heatmap_ggplot <- reactive({
      current_data <- data()
      req(current_data)
      plot_data <- annotation_heatmap_data(current_data$abundance, current_data$metadata)
      validate(need(nrow(plot_data) > 0, "No abundance rows are available for annotation heatmap."))

      plot_annotation_heatmap(plot_data)
    })

    output$abundance_annotation_heatmap <- renderPlotly({
      ggplotly(abundance_annotation_heatmap_ggplot(), tooltip = "text") |>
        apply_proxiome_plot_frame(colorbar_title = "Median abundance")
    })
    register_ggplot_downloads(
      output,
      "abundance_annotation_heatmap",
      abundance_annotation_heatmap_ggplot,
      filename_prefix = "abundance-annotation-heatmap",
      width = 10,
      height = 5
    )

    output$abundance_celltype_composition_table <- renderTable({
      current_data <- data()
      req(current_data)
      composition <- celltype_composition_data(current_data$metadata)
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

    abundance_diff_volcano_x_label <- reactive({
      config <- abundance_diff_config()
      req(config)
      paste("Abundance effect:", config$group_a, "minus", config$group_b, "(reference)")
    })

    abundance_diff_volcano_dimensions <- reactive({
      differential_volcano_dimensions(abundance_diff_volcano_x_label())
    })

    abundance_diff_volcano_ggplot <- reactive({
      config <- abundance_diff_config()
      req(config)
      result <- abundance_diff_results()
      validate(need(nrow(result) > 0, "Choose two different groups with enough abundance data."))

      differential_volcano_ggplot(
        result,
        label_col = "marker",
        x_label = abundance_diff_volcano_x_label(),
        fdr_cutoff = config$fdr_cutoff,
        effect_cutoff = config$effect_cutoff
      )
    })

    output$abundance_diff_volcano <- renderPlotly({
      dimensions <- abundance_diff_volcano_dimensions()
      ggplotly(
        abundance_diff_volcano_ggplot(),
        tooltip = "text",
        source = "abundance_diff",
        width = dimensions$width,
        height = dimensions$height
      ) |>
        apply_differential_plot_frame(dimensions = dimensions)
    })
    register_ggplot_downloads(
      output,
      "abundance_diff_volcano",
      abundance_diff_volcano_ggplot,
      filename_prefix = function() paste("abundance-differential-volcano", abundance_diff_volcano_x_label(), sep = "-"),
      width = function() plot_download_size_from_dimensions(abundance_diff_volcano_dimensions())$width,
      height = function() plot_download_size_from_dimensions(abundance_diff_volcano_dimensions())$height
    )

    observeEvent(plotly::event_data("plotly_click", source = "abundance_diff"), {
      event <- plotly::event_data("plotly_click", source = "abundance_diff")
      if (!is.null(event$key) && nzchar(event$key)) {
        updateSelectInput(session, "abundance_diff_marker", selected = event$key)
      }
    })

    abundance_diff_detail_data <- reactive({
      config <- abundance_diff_config()
      req(config, input$abundance_diff_marker, config$group_a, config$group_b)

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
      dimensions <- differential_detail_dimensions(
        plot_data,
        stratify_by_celltype = isTRUE(config$stratify_by_celltype),
        y_label = y_label
      )

      list(
        config = config,
        plot_data = plot_data,
        y_label = y_label,
        dimensions = dimensions
      )
    })

    abundance_diff_detail_ggplot <- reactive({
      detail <- abundance_diff_detail_data()
      plot_data <- detail$plot_data
      config <- detail$config

      p <- ggplot(plot_data, aes(condition, abundance, color = condition, text = hover)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.18, linewidth = 0.5) +
        geom_jitter(width = 0.18, height = 0, alpha = 0.5, size = 1.4) +
        labs(x = NULL, y = detail$y_label) +
        theme_minimal(base_size = 12) +
        theme(panel.grid.minor = element_blank(), legend.position = "none")

      if (isTRUE(config$stratify_by_celltype)) {
        p <- p + facet_wrap(~celltype_manual, scales = "free_y")
      }

      p
    })

    output$abundance_diff_detail <- renderPlotly({
      dimensions <- abundance_diff_detail_data()$dimensions
      ggplotly(abundance_diff_detail_ggplot(), tooltip = "text", width = dimensions$width, height = dimensions$height) |>
        apply_proxiome_plot_frame(dimensions = dimensions)
    })
    register_ggplot_downloads(
      output,
      "abundance_diff_detail",
      abundance_diff_detail_ggplot,
      filename_prefix = function() paste("abundance-differential-detail", input$abundance_diff_marker %||% "marker", sep = "-"),
      width = function() plot_download_size_from_dimensions(abundance_diff_detail_data()$dimensions)$width,
      height = function() plot_download_size_from_dimensions(abundance_diff_detail_data()$dimensions)$height
    )

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
  })
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

abundance_umap_widget_dimensions <- function(width_px = NULL, height_px = NULL) {
  width_override <- plot_dimension_override(width_px)
  height_override <- plot_dimension_override(height_px)

  list(
    width = width_override %||% 832,
    height = height_override %||% 520,
    margin = proxiome_plot_margins()
  )
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
