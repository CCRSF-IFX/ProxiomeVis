clustering_sidebar <- function(id) {
  ns <- NS(id)

  sidebar(
    title = "Clustering controls",
    width = 300,
    conditionalPanel(
      condition = "input.clustering_mode == 'Observed' || input.clustering_mode == 'Per Marker'",
      ns = ns,
      accordion(
        open = c("Display", "Filters"),
        accordion_panel(
          "Display",
          selectInput(ns("clustering_marker"), "Marker", choices = character(0))
        ),
        accordion_panel(
          "Filters",
          selectizeInput(ns("clustering_condition_filter"), "Condition", choices = character(0), multiple = TRUE),
          selectizeInput(ns("clustering_celltype_filter"), "Cell type", choices = character(0), multiple = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.clustering_mode == 'Summary Heatmap'",
      ns = ns,
      accordion(
        open = c("Display", "Filters"),
        accordion_panel(
          "Display",
          numericInput(ns("clustering_heatmap_marker_count"), "Top markers", value = 20, min = 2, max = 40, step = 1)
        ),
        accordion_panel(
          "Filters",
          selectizeInput(ns("clustering_heatmap_condition_filter"), "Condition", choices = character(0), multiple = TRUE),
          selectizeInput(ns("clustering_heatmap_celltype_filter"), "Cell type", choices = character(0), multiple = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.clustering_mode == 'Differential'",
      ns = ns,
      accordion(
        open = c("Contrast", "Thresholds"),
        accordion_panel(
          "Contrast",
          selectInput(ns("clustering_diff_group_a"), "Group A", choices = character(0)),
          selectInput(ns("clustering_diff_group_b"), "Group B (reference)", choices = character(0)),
          selectizeInput(ns("clustering_diff_celltype_filter"), "Cell type", choices = character(0), multiple = TRUE),
          checkboxInput(ns("clustering_diff_stratify_celltype"), "Stratify by cell type", value = FALSE),
          actionButton(ns("clustering_run_differential"), "Run differential analysis", class = "btn-primary w-100")
        ),
        accordion_panel(
          "Thresholds",
          numericInput(ns("clustering_diff_fdr"), "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput(ns("clustering_diff_effect"), "Minimum effect", value = 0.25, min = 0, step = 0.05),
          numericInput(ns("clustering_diff_min_cells"), "Minimum cells per group", value = 3, min = 1, step = 1)
        ),
        accordion_panel(
          "Detail",
          selectInput(ns("clustering_diff_marker"), "Detail marker", choices = character(0))
        )
      )
    )
  )
}

clustering_module_ui <- function(id) {
  ns <- NS(id)

  nav_panel(
    "Clustering",
    layout_sidebar(
      sidebar = clustering_sidebar(id),
      navset_card_underline(
        id = ns("clustering_mode"),
        title = "Clustering",
        full_screen = TRUE,
        nav_panel(
          "Observed",
          plot_pane(
            size = "standard",
            download_id = "clustering_plot",
            ns = ns,
            plotlyOutput(ns("clustering_plot"), height = "auto")
          ),
          div(class = "table-pane", tableOutput(ns("clustering_table")))
        ),
        nav_panel(
          "Per Marker",
          plot_pane(
            size = "compact",
            download_id = "clustering_per_marker_plot",
            ns = ns,
            plotlyOutput(ns("clustering_per_marker_plot"), height = "auto")
          ),
          div(class = "table-pane", tableOutput(ns("clustering_per_marker_table")))
        ),
        nav_panel(
          "Summary Heatmap",
          plot_pane(
            size = "scroll",
            download_id = "clustering_summary_heatmap",
            ns = ns,
            plotlyOutput(ns("clustering_summary_heatmap"), height = "auto")
          ),
          div(class = "table-pane", tableOutput(ns("clustering_summary_heatmap_table")))
        ),
        nav_panel(
          "Differential",
          uiOutput(ns("clustering_diff_summary")),
          differential_plot_row(ns("clustering_diff_volcano"), ns("clustering_diff_detail")),
          div(class = "table-pane", tableOutput(ns("clustering_diff_table")))
        )
      )
    )
  )
}

clustering_module_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    clustering_diff_config <- reactiveVal(NULL)

    clustering_differential_config_from_inputs <- function(current_data) {
      make_differential_config(
        group_a = input$clustering_diff_group_a,
        group_b = input$clustering_diff_group_b,
        celltype_filter = selected_or_all(
          input$clustering_diff_celltype_filter,
          unique(current_data$metadata$celltype_manual)
        ),
        stratify_by_celltype = input$clustering_diff_stratify_celltype,
        min_cells = numeric_input_value(input$clustering_diff_min_cells, 3),
        fdr_cutoff = numeric_input_value(input$clustering_diff_fdr, 0.05),
        effect_cutoff = numeric_input_value(input$clustering_diff_effect, 0.25)
      )
    }

    observe({
      current_data <- data()
      req(current_data)

      conditions <- sort(unique(current_data$metadata$condition))
      cell_types <- sort(unique(current_data$metadata$celltype_manual))
      default_group_a <- conditions[1]
      default_group_b <- conditions[min(2, length(conditions))]

      updateSelectInput(session, "clustering_marker", choices = current_data$marker_options, selected = current_data$marker_options[1])
      updateSelectizeInput(session, "clustering_condition_filter", choices = conditions, selected = conditions)
      updateSelectizeInput(session, "clustering_celltype_filter", choices = cell_types, selected = cell_types)
      updateSelectizeInput(session, "clustering_heatmap_condition_filter", choices = conditions, selected = conditions)
      updateSelectizeInput(session, "clustering_heatmap_celltype_filter", choices = cell_types, selected = cell_types)
      updateSelectInput(session, "clustering_diff_group_a", choices = conditions, selected = default_group_a)
      updateSelectInput(session, "clustering_diff_group_b", choices = conditions, selected = default_group_b)
      updateSelectizeInput(session, "clustering_diff_celltype_filter", choices = cell_types, selected = cell_types)
      updateSelectInput(session, "clustering_diff_marker", choices = current_data$marker_options, selected = current_data$marker_options[1])

      clustering_diff_config(default_differential_config(conditions, cell_types))
    })

    observeEvent(input$clustering_run_differential, {
      current_data <- data()
      req(current_data)
      clustering_diff_config(clustering_differential_config_from_inputs(current_data))
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

    clustering_metadata <- reactive({
      filtered_metadata_for(input$clustering_condition_filter, input$clustering_celltype_filter)
    })

    clustering_diff_results <- reactive({
      current_data <- data()
      config <- clustering_diff_config()
      req(current_data, config, config$group_a, config$group_b)

      calculate_differential_readout(
        current_data$clustering,
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

    clustering_points <- reactive({
      current_data <- data()
      req(current_data, input$clustering_marker)

      metadata <- clustering_metadata()
      clustering <- current_data$clustering[current_data$clustering$marker == input$clustering_marker, , drop = FALSE]
      clustering <- clustering[clustering$component %in% metadata$component, , drop = FALSE]
      abundance <- current_data$abundance[current_data$abundance$marker == input$clustering_marker, , drop = FALSE]
      merge(clustering, abundance, by = c("component", "marker"), all.x = TRUE, sort = FALSE)
    })

    clustering_plot_ggplot <- reactive({
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

      p
    })

    clustering_plot_dimensions <- reactive({
      plot_options_input_dimensions(input, "clustering_plot")
    })
    clustering_plot_export_dimensions <- reactive({
      plot_options_export_dimensions(input, "clustering_plot")
    })

    output$clustering_plot <- renderPlotly({
      dimensions <- clustering_plot_dimensions()
      display_dimensions <- plotly_display_dimensions(dimensions)
      ggplotly(clustering_plot_ggplot(), tooltip = "text", width = display_dimensions$width, height = display_dimensions$height) |>
        apply_proxiome_plot_frame(dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "clustering_plot",
      clustering_plot_ggplot,
      filename_prefix = function() paste("clustering-observed", input$clustering_marker %||% "marker", sep = "-"),
      width = function() plot_download_size_from_dimensions(clustering_plot_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(clustering_plot_export_dimensions())$height
    )

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

    clustering_per_marker_ggplot <- reactive({
      plot_data <- clustering_points()
      validate(need(nrow(plot_data) > 0, "No self-clustering scores are available for the selected marker and filters."))

      plot_clustering_per_marker(plot_data, input$clustering_marker)
    })

    clustering_per_marker_plot_dimensions <- reactive({
      plot_options_input_dimensions(input, "clustering_per_marker_plot")
    })
    clustering_per_marker_plot_export_dimensions <- reactive({
      plot_options_export_dimensions(input, "clustering_per_marker_plot")
    })

    output$clustering_per_marker_plot <- renderPlotly({
      dimensions <- clustering_per_marker_plot_dimensions()
      display_dimensions <- plotly_display_dimensions(dimensions)
      ggplotly(clustering_per_marker_ggplot(), tooltip = c("x", "y", "fill"), width = display_dimensions$width, height = display_dimensions$height) |>
        apply_proxiome_plot_frame(dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "clustering_per_marker_plot",
      clustering_per_marker_ggplot,
      filename_prefix = function() paste("clustering-per-marker", input$clustering_marker %||% "marker", sep = "-"),
      width = function() plot_download_size_from_dimensions(clustering_per_marker_plot_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(clustering_per_marker_plot_export_dimensions())$height
    )

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
      current_data <- data()
      req(current_data)

      clustering <- current_data$clustering
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

    clustering_summary_heatmap_ggplot <- reactive({
      summary <- clustering_heatmap_summary()
      validate(need(nrow(summary) > 0, "No clustering rows are available for the selected heatmap filters."))

      plot_clustering_summary_heatmap(summary)
    })

    clustering_summary_heatmap_dimensions <- reactive({
      plot_options_input_dimensions(input, "clustering_summary_heatmap")
    })
    clustering_summary_heatmap_export_dimensions <- reactive({
      plot_options_export_dimensions(input, "clustering_summary_heatmap")
    })

    output$clustering_summary_heatmap <- renderPlotly({
      dimensions <- clustering_summary_heatmap_dimensions()
      display_dimensions <- plotly_display_dimensions(dimensions)
      ggplotly(clustering_summary_heatmap_ggplot(), tooltip = "text", width = display_dimensions$width, height = display_dimensions$height) |>
        apply_proxiome_plot_frame(colorbar_title = "Mean log2 ratio", dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "clustering_summary_heatmap",
      clustering_summary_heatmap_ggplot,
      filename_prefix = "clustering-summary-heatmap",
      width = function() plot_download_size_from_dimensions(clustering_summary_heatmap_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(clustering_summary_heatmap_export_dimensions())$height
    )

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

    clustering_diff_volcano_x_label <- reactive({
      config <- clustering_diff_config()
      req(config)
      paste("Difference in medians:", config$group_a, "minus", config$group_b, "(reference)")
    })

    clustering_diff_volcano_dimensions <- reactive({
      plot_options_view_overrides(input, "clustering_diff_volcano", differential_volcano_dimensions(clustering_diff_volcano_x_label()))
    })
    clustering_diff_volcano_export_dimensions <- reactive({
      plot_options_export_overrides(input, "clustering_diff_volcano", differential_volcano_dimensions(clustering_diff_volcano_x_label()))
    })

    clustering_diff_volcano_ggplot <- reactive({
      config <- clustering_diff_config()
      req(config)
      result <- clustering_diff_results()
      validate(need(nrow(result) > 0, "Choose two different groups with enough self-clustering data."))

      differential_volcano_ggplot(
        result,
        label_col = "marker",
        x_label = clustering_diff_volcano_x_label(),
        fdr_cutoff = config$fdr_cutoff,
        effect_cutoff = config$effect_cutoff
      )
    })

    output$clustering_diff_volcano <- renderPlotly({
      dimensions <- clustering_diff_volcano_dimensions()
      display_dimensions <- responsive_plotly_dimensions(dimensions)
      ggplotly(
        clustering_diff_volcano_ggplot(),
        tooltip = "text",
        source = "clustering_diff",
        width = display_dimensions$width,
        height = display_dimensions$height
      ) |>
        apply_differential_plot_frame(dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "clustering_diff_volcano",
      clustering_diff_volcano_ggplot,
      filename_prefix = function() paste("clustering-differential-volcano", clustering_diff_volcano_x_label(), sep = "-"),
      width = function() plot_download_size_from_dimensions(clustering_diff_volcano_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(clustering_diff_volcano_export_dimensions())$height
    )

    observeEvent(plotly::event_data("plotly_click", source = "clustering_diff"), {
      event <- plotly::event_data("plotly_click", source = "clustering_diff")
      if (!is.null(event$key) && nzchar(event$key)) {
        updateSelectInput(session, "clustering_diff_marker", selected = event$key)
      }
    })

    clustering_diff_detail_data <- reactive({
      current_data <- data()
      config <- clustering_diff_config()
      req(current_data, config, input$clustering_diff_marker, config$group_a, config$group_b)

      plot_data <- current_data$clustering[
        current_data$clustering$marker == input$clustering_diff_marker &
          current_data$clustering$condition %in% c(config$group_a, config$group_b) &
          current_data$clustering$celltype_manual %in% config$celltype_filter,
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
      base_dimensions <- differential_detail_dimensions(
        plot_data,
        stratify_by_celltype = isTRUE(config$stratify_by_celltype),
        y_label = y_label
      )
      dimensions <- plot_options_view_overrides(input, "clustering_diff_detail", base_dimensions)
      export_dimensions <- plot_options_export_overrides(input, "clustering_diff_detail", base_dimensions)

      list(
        config = config,
        plot_data = plot_data,
        y_label = y_label,
        dimensions = dimensions,
        export_dimensions = export_dimensions
      )
    })

    clustering_diff_detail_ggplot <- reactive({
      detail <- clustering_diff_detail_data()
      plot_data <- detail$plot_data
      config <- detail$config

      p <- ggplot(plot_data, aes(condition, log2_ratio, color = condition, text = hover)) +
        geom_hline(yintercept = 0, color = "#8a9699", linewidth = 0.5) +
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

    output$clustering_diff_detail <- renderPlotly({
      dimensions <- clustering_diff_detail_data()$dimensions
      display_dimensions <- responsive_plotly_dimensions(dimensions)
      ggplotly(clustering_diff_detail_ggplot(), tooltip = "text", width = display_dimensions$width, height = display_dimensions$height) |>
        apply_proxiome_plot_frame(dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "clustering_diff_detail",
      clustering_diff_detail_ggplot,
      filename_prefix = function() paste("clustering-differential-detail", input$clustering_diff_marker %||% "marker", sep = "-"),
      width = function() plot_download_size_from_dimensions(clustering_diff_detail_data()$export_dimensions)$width,
      height = function() plot_download_size_from_dimensions(clustering_diff_detail_data()$export_dimensions)$height
    )

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
  })
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
