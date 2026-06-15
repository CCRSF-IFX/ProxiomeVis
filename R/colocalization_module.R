colocalization_sidebar <- function(id) {
  ns <- NS(id)

  sidebar(
    title = "Colocalization controls",
    width = 300,
    conditionalPanel(
      condition = "input.colocalization_mode == 'Observed'",
      ns = ns,
      accordion(
        open = c("Display", "Filters"),
        accordion_panel(
          "Display",
          radioButtons(
            ns("colocalization_heatmap_display"),
            "Plot style",
            choices = c("Interactive" = "interactive", "Original R plot" = "original"),
            selected = "interactive"
          ),
          selectInput(
            ns("colocalization_heatmap_preset"),
            "Heatmap preset",
            choices = c("Custom" = "custom", "Report style" = "report"),
            selected = "custom"
          ),
          selectInput(
            ns("spatial_coloc_scope"),
            "Heatmap scope",
            choices = c(
              "Condition summary" = "condition",
              "Sample summary" = "sample",
              "Cell type focus" = "celltype"
            ),
            selected = "condition"
          ),
          conditionalPanel(
            condition = "input.spatial_coloc_scope == 'celltype'",
            ns = ns,
            selectInput(ns("spatial_celltype_focus"), "Cell type focus", choices = character(0))
          ),
          selectInput(
            ns("spatial_marker_selection_mode"),
            "Marker set",
            choices = c(
              "Variable detected markers" = "auto",
              "Selected markers" = "manual"
            ),
            selected = "auto"
          ),
          selectizeInput(ns("colocalization_heatmap_markers"), "Heatmap markers", choices = character(0), multiple = TRUE),
          conditionalPanel(
            condition = "input.spatial_marker_selection_mode == 'auto'",
            ns = ns,
            numericInput(ns("spatial_top_marker_count"), "Top markers", value = 20, min = 2, max = 40, step = 1),
            numericInput(ns("spatial_min_pct_detected"), "Minimum fraction detected", value = 0.25, min = 0, max = 1, step = 0.05),
            numericInput(ns("spatial_min_log2_range"), "Minimum log2 range", value = 0.2, min = 0, step = 0.05)
          ),
          selectInput(ns("colocalization_reference_condition"), "Reference condition", choices = character(0)),
          selectInput(
            ns("colocalization_clustering_method"),
            "Marker ordering",
            choices = c("Ward D2" = "ward.D2", "Complete" = "complete", "Average" = "average", "Single" = "single"),
            selected = "ward.D2"
          ),
          numericInput(ns("colocalization_legend_min"), "Legend minimum", value = -1, step = 0.1),
          numericInput(ns("colocalization_legend_max"), "Legend maximum", value = 1, step = 0.1)
        ),
        accordion_panel(
          "Filters",
          selectizeInput(ns("colocalization_condition_filter"), "Condition", choices = character(0), multiple = TRUE),
          selectizeInput(ns("colocalization_celltype_filter"), "Cell type", choices = character(0), multiple = TRUE)
        )
      )
    ),
    conditionalPanel(
      condition = "input.colocalization_mode == 'Differential'",
      ns = ns,
      accordion(
        open = c("Contrast", "Thresholds"),
        accordion_panel(
          "Contrast",
          selectInput(ns("colocalization_diff_group_a"), "Group A", choices = character(0)),
          selectInput(ns("colocalization_diff_group_b"), "Group B (reference)", choices = character(0)),
          selectizeInput(ns("colocalization_diff_celltype_filter"), "Cell type", choices = character(0), multiple = TRUE),
          checkboxInput(ns("colocalization_diff_stratify_celltype"), "Stratify by cell type", value = FALSE),
          actionButton(ns("colocalization_run_differential"), "Run differential analysis", class = "btn-primary w-100")
        ),
        accordion_panel(
          "Thresholds",
          numericInput(ns("colocalization_diff_fdr"), "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput(ns("colocalization_diff_effect"), "Minimum effect", value = 0.25, min = 0, step = 0.05),
          numericInput(ns("colocalization_diff_min_cells"), "Minimum cells per group", value = 3, min = 1, step = 1)
        ),
        accordion_panel(
          "Detail",
          selectInput(ns("colocalization_diff_anchor_marker"), "Anchor marker", choices = character(0)),
          selectInput(ns("colocalization_diff_pair"), "Detail pair", choices = character(0))
        )
      )
    ),
    conditionalPanel(
      condition = "input.colocalization_mode == '3D Layout'",
      ns = ns,
      accordion(
        open = c("Cell", "Markers"),
        accordion_panel(
          "Cell",
          selectInput(ns("colocalization_3d_sample"), "Sample", choices = character(0)),
          selectizeInput(ns("colocalization_3d_celltype_filter"), "Cell type", choices = character(0), multiple = TRUE),
          selectInput(ns("colocalization_3d_component"), "Cell/component", choices = character(0)),
          numericInput(ns("colocalization_3d_max_background"), "Max background nodes", value = 7000, min = 0, max = 50000, step = 500)
        ),
        accordion_panel(
          "Markers",
          selectizeInput(ns("colocalization_3d_markers"), "Highlighted markers", choices = character(0), multiple = TRUE)
        )
      )
    )
  )
}

colocalization_module_ui <- function(id) {
  ns <- NS(id)

  nav_panel(
    "Colocalization",
    layout_sidebar(
      sidebar = colocalization_sidebar(id),
      navset_card_underline(
        id = ns("colocalization_mode"),
        title = "Colocalization",
        full_screen = TRUE,
        nav_panel(
          "Observed",
          plot_pane(
            size = "scroll",
            extra_class = "coloc-heatmap-pane",
            download_id = "colocalization_heatmap",
            ns = ns,
            conditionalPanel(
              condition = "input.colocalization_heatmap_display == 'interactive'",
              ns = ns,
              plotlyOutput(ns("colocalization_heatmap_interactive"), height = "auto")
            ),
            conditionalPanel(
              condition = "input.colocalization_heatmap_display == 'original'",
              ns = ns,
              plotOutput(ns("colocalization_heatmap_original"), height = "auto")
            )
          ),
          div(class = "table-pane", tableOutput(ns("colocalization_table")))
        ),
        nav_panel(
          "Differential",
          uiOutput(ns("colocalization_diff_summary")),
          differential_plot_row(ns("colocalization_diff_volcano"), ns("colocalization_diff_detail")),
          div(class = "table-pane", tableOutput(ns("colocalization_diff_table")))
        ),
        nav_panel(
          "3D Layout",
          plot_pane(
            size = "wide",
            controls = plot_options_controls(
              ns,
              "colocalization_3d_layout_width",
              "colocalization_3d_layout_height",
              width_value = 832,
              height_value = 620,
              min_height = 420
            ),
            plotlyOutput(ns("colocalization_3d_layout"), height = "auto")
          ),
          div(class = "table-pane", tableOutput(ns("colocalization_3d_component_table")))
        )
      )
    )
  )
}

colocalization_module_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    colocalization_diff_config <- reactiveVal(NULL)

    colocalization_differential_config_from_inputs <- function(current_data, anchor_marker = NULL) {
      make_differential_config(
        group_a = input$colocalization_diff_group_a,
        group_b = input$colocalization_diff_group_b,
        celltype_filter = selected_or_all(
          input$colocalization_diff_celltype_filter,
          unique(current_data$metadata$celltype_manual)
        ),
        stratify_by_celltype = input$colocalization_diff_stratify_celltype,
        min_cells = numeric_input_value(input$colocalization_diff_min_cells, 3),
        fdr_cutoff = numeric_input_value(input$colocalization_diff_fdr, 0.05),
        effect_cutoff = numeric_input_value(input$colocalization_diff_effect, 0.25),
        anchor_marker = anchor_marker
      )
    }

    observe({
      current_data <- data()
      req(current_data)

      conditions <- sort(unique(current_data$metadata$condition))
      cell_types <- sort(unique(current_data$metadata$celltype_manual))
      default_group_a <- conditions[1]
      default_group_b <- conditions[min(2, length(conditions))]
      colocalization_pairs <- sort(unique(current_data$colocalization$marker_pair))
      colocalization_markers <- available_colocalization_marker_choices(current_data$colocalization)
      default_heatmap_markers <- head(colocalization_markers, min(20L, length(colocalization_markers)))
      default_colocalization_reference <- if ("CD3CD28" %in% conditions) "CD3CD28" else conditions[1]

      updateSelectizeInput(session, "colocalization_heatmap_markers", choices = colocalization_markers, selected = default_heatmap_markers)
      updateSelectInput(session, "spatial_celltype_focus", choices = cell_types, selected = cell_types[1])
      updateSelectInput(
        session,
        "colocalization_reference_condition",
        choices = conditions,
        selected = default_colocalization_reference
      )

      updateSelectizeInput(session, "colocalization_condition_filter", choices = conditions, selected = conditions)
      updateSelectizeInput(session, "colocalization_celltype_filter", choices = cell_types, selected = cell_types)

      updateSelectInput(session, "colocalization_diff_group_a", choices = conditions, selected = default_group_a)
      updateSelectInput(session, "colocalization_diff_group_b", choices = conditions, selected = default_group_b)
      updateSelectizeInput(session, "colocalization_diff_celltype_filter", choices = cell_types, selected = cell_types)
      updateSelectInput(session, "colocalization_diff_anchor_marker", choices = current_data$marker_options, selected = current_data$marker_options[1])
      updateSelectInput(session, "colocalization_diff_pair", choices = colocalization_pairs, selected = colocalization_pairs[1])

      sample_col <- colocalization_3d_sample_column(current_data$metadata)
      samples <- sort(unique(as.character(current_data$metadata[[sample_col]])))
      default_sample <- if ("3_CD3CD28" %in% samples) "3_CD3CD28" else samples[1]
      default_3d_markers <- intersect(c("ICAM-1", "CD54", "CD40", "CD8", "CD3e", "CD81", "CD82"), current_data$marker_options)
      if (length(default_3d_markers) == 0) {
        default_3d_markers <- head(current_data$marker_options, min(4L, length(current_data$marker_options)))
      }
      updateSelectInput(session, "colocalization_3d_sample", choices = samples, selected = default_sample)
      updateSelectizeInput(
        session,
        "colocalization_3d_celltype_filter",
        choices = cell_types,
        selected = if ("CD8 T" %in% cell_types) "CD8 T" else cell_types[1]
      )
      updateSelectizeInput(session, "colocalization_3d_markers", choices = current_data$marker_options, selected = default_3d_markers)

      colocalization_diff_config(default_differential_config(
        conditions,
        cell_types,
        anchor_marker = current_data$marker_options[1]
      ))
    })

    observeEvent(input$colocalization_run_differential, {
      current_data <- data()
      req(current_data)
      colocalization_diff_config(colocalization_differential_config_from_inputs(
        current_data,
        anchor_marker = input$colocalization_diff_anchor_marker
      ))
    }, ignoreInit = TRUE)

    observeEvent(input$colocalization_heatmap_preset, {
      if (!identical(input$colocalization_heatmap_preset, "report")) {
        return(invisible(NULL))
      }

      updateSelectInput(session, "spatial_marker_selection_mode", selected = "auto")
      updateNumericInput(session, "spatial_top_marker_count", value = 40)
      updateNumericInput(session, "spatial_min_pct_detected", value = 0.25)
      updateNumericInput(session, "spatial_min_log2_range", value = 0.2)
      updateNumericInput(session, "colocalization_legend_min", value = -0.75)
      updateNumericInput(session, "colocalization_legend_max", value = 0.75)
      updateSelectInput(session, "colocalization_clustering_method", selected = "ward.D2")
    }, ignoreInit = TRUE)

    observe({
      current_data <- data()
      config <- colocalization_diff_config()
      req(current_data, config$anchor_marker)

      choices <- sort(unique(current_data$colocalization$marker_pair[
        current_data$colocalization$marker_1 == config$anchor_marker |
          current_data$colocalization$marker_2 == config$anchor_marker
      ]))

      if (length(choices) == 0) {
        choices <- sort(unique(current_data$colocalization$marker_pair))
      }

      updateSelectInput(session, "colocalization_diff_pair", choices = choices, selected = choices[1])
    })

    observe({
      current_data <- data()
      req(current_data)

      choices <- colocalization_3d_component_choices(
        current_data$metadata,
        sample = input$colocalization_3d_sample,
        cell_types = input$colocalization_3d_celltype_filter
      )
      updateSelectInput(session, "colocalization_3d_component", choices = choices, selected = unname(choices[1]))
    })

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

    colocalization_metadata <- reactive({
      filtered_metadata_for(input$colocalization_condition_filter, input$colocalization_celltype_filter)
    })

    colocalization_diff_results <- reactive({
      current_data <- data()
      config <- colocalization_diff_config()
      req(current_data, config, config$group_a, config$group_b)

      calculate_differential_readout(
        current_data$colocalization,
        feature_cols = c("marker_pair", "marker_1", "marker_2"),
        value_col = "log2_ratio",
        group_a = config$group_a,
        group_b = config$group_b,
        celltype_filter = config$celltype_filter,
        stratify_by_celltype = config$stratify_by_celltype,
        min_cells = config$min_cells,
        fdr_cutoff = config$fdr_cutoff
      )
    })

    colocalization_diff_anchor_results <- reactive({
      result <- colocalization_diff_results()
      config <- colocalization_diff_config()
      req(config$anchor_marker)

      result[
        result$marker_1 == config$anchor_marker |
          result$marker_2 == config$anchor_marker,
        ,
        drop = FALSE
      ]
    })

    colocalization_heatmap_rows <- reactive({
      current_data <- data()
      req(current_data)

      metadata <- colocalization_metadata()
      current_data$colocalization[current_data$colocalization$component %in% metadata$component, , drop = FALSE]
    })

    colocalization_heatmap_result <- reactive({
      current_data <- data()
      req(current_data)

      colocalization <- colocalization_heatmap_rows()
      validate(need(nrow(colocalization) > 0, "No colocalization scores are available for the selected filters."))
      if (!"sample_alias" %in% names(colocalization)) {
        colocalization$sample_alias <- "sample"
      }

      scope <- input$spatial_coloc_scope %||% "condition"
      if (length(scope) == 0) {
        scope <- "condition"
      }
      if (!scope %in% c("condition", "sample", "celltype")) {
        scope <- "condition"
      }
      marker_selection_mode <- input$spatial_marker_selection_mode %||% "auto"
      if (length(marker_selection_mode) == 0) {
        marker_selection_mode <- "auto"
      }

      selected_conditions <- selected_or_all(
        input$colocalization_condition_filter,
        sort(unique(as.character(current_data$metadata$condition)))
      )
      selected_conditions <- intersect(selected_conditions, unique(as.character(colocalization$condition)))
      validate(need(length(selected_conditions) > 0, "Select at least one condition for the colocalization heatmap."))

      colocalization <- colocalization[colocalization$condition %in% selected_conditions, , drop = FALSE]

      if (identical(scope, "celltype")) {
        focus_celltype <- input$spatial_celltype_focus
        if (is.null(focus_celltype) || length(focus_celltype) == 0 || is.na(focus_celltype) || !nzchar(focus_celltype)) {
          focus_celltype <- sort(unique(as.character(colocalization$celltype_manual)))[1]
        }
        colocalization <- colocalization[colocalization$celltype_manual == focus_celltype, , drop = FALSE]
        validate(need(nrow(colocalization) > 0, "No colocalization scores are available for the selected cell type focus."))
      }

      available_markers <- available_colocalization_marker_choices(colocalization)
      requested_markers <- selected_or_all(input$colocalization_heatmap_markers, available_markers)
      requested_markers <- intersect(requested_markers, available_markers)
      candidate_markers <- if (identical(marker_selection_mode, "auto")) available_markers else requested_markers
      validate(need(length(candidate_markers) >= 2, "Select at least two markers for the spatial heatmap."))

      marker_summary <- spatial_heatmap_summary_for_scope(
        colocalization = colocalization,
        selected_markers = candidate_markers,
        scope = scope,
        selected_conditions = selected_conditions
      )
      selected_markers <- spatial_heatmap_selected_markers(
        summary = marker_summary,
        available_markers = candidate_markers,
        requested_markers = requested_markers,
        marker_selection_mode = marker_selection_mode,
        n_markers = input$spatial_top_marker_count,
        min_pct_detected = input$spatial_min_pct_detected,
        min_range = input$spatial_min_log2_range
      )
      validate(need(length(selected_markers) >= 2, "Select at least two markers for the spatial heatmap."))
      validate(need(length(selected_markers) <= 40, "Use 40 or fewer markers for an interpretable spatial heatmap."))

      summary <- spatial_heatmap_summary_for_scope(
        colocalization = colocalization,
        selected_markers = selected_markers,
        scope = scope,
        selected_conditions = selected_conditions
      )
      validate(need(nrow(summary) > 0, "No colocalization scores are available for the selected markers."))

      group_cols <- spatial_heatmap_group_cols(scope)
      summary <- complete_spatial_marker_pairs(
        summary = summary,
        selected_markers = selected_markers,
        group_cols = group_cols
      )
      condition_col <- if (identical(scope, "condition")) "condition" else "sample_alias"
      plot_groups <- unique(as.character(summary[[condition_col]]))
      plot_groups <- plot_groups[nzchar(plot_groups)]
      validate(need(length(plot_groups) > 0, "No spatial heatmap groups are available for the selected filters."))

      result <- make_coloc_heatmaps(
        data = summary,
        selected_markers = selected_markers,
        cell_label = spatial_heatmap_cell_label(
          scope,
          selected_celltypes = input$colocalization_celltype_filter,
          focus_celltype = input$spatial_celltype_focus
        ),
        conditions = plot_groups,
        reference_condition = if (identical(scope, "condition")) {
          input$colocalization_reference_condition %||% selected_conditions[1]
        } else {
          plot_groups[1]
        },
        condition_col = condition_col,
        clustering_method = input$colocalization_clustering_method %||% "ward.D2",
        legend_range = colocalization_legend_range(
          input$colocalization_legend_min,
          input$colocalization_legend_max
        )
      )
      result$summary <- summary
      result$scope <- scope
      result$condition_col <- condition_col
      result
    })

    output$colocalization_heatmap_interactive <- renderPlotly({
      coloc_heatmap_plotly(colocalization_heatmap_result(), dimensions = colocalization_heatmap_dimensions())
    })

    output$colocalization_heatmap_original <- renderPlot({
      print(colocalization_heatmap_result()$plot)
    }, width = function() {
      colocalization_heatmap_dimensions()$width
    }, height = function() {
      colocalization_heatmap_dimensions()$height
    })

    colocalization_heatmap_ggplot <- reactive({
      colocalization_heatmap_result()$plot
    })
    colocalization_heatmap_dimensions <- reactive({
      apply_plot_options_overrides(
        coloc_heatmap_widget_dimensions(colocalization_heatmap_result()$plot_data),
        width_px = input$colocalization_heatmap_width,
        height_px = input$colocalization_heatmap_height
      )
    })
    register_ggplot_downloads(
      output,
      "colocalization_heatmap",
      colocalization_heatmap_ggplot,
      filename_prefix = "colocalization-heatmap",
      width = function() plot_download_size_from_dimensions(colocalization_heatmap_dimensions())$width,
      height = function() plot_download_size_from_dimensions(colocalization_heatmap_dimensions())$height
    )

    output$colocalization_table <- renderTable({
      summary <- colocalization_heatmap_result()$summary
      validate(need(nrow(summary) > 0, "No spatial metric rows to summarize."))

      format_spatial_heatmap_table(summary)
    }, striped = TRUE, bordered = FALSE, width = "100%")

    colocalization_3d_dimensions <- reactive({
      plot_options_input_dimensions(
        input,
        "colocalization_3d_layout",
        default_width = 832,
        default_height = 620,
        margin = list(l = 0, r = 0, t = 40, b = 0)
      )
    })

    colocalization_3d_layout_raw <- reactive({
      current_data <- data()
      req(current_data, input$colocalization_3d_sample, input$colocalization_3d_component)

      layout_path <- pixelator_layout_pxl_path(input$colocalization_3d_sample, source = current_data$source)
      validate(need(nzchar(layout_path), paste("No Pixelator 3D layout file found for", input$colocalization_3d_sample)))

      raw_component <- pixelator_raw_component_id(input$colocalization_3d_component, input$colocalization_3d_sample)
      layout <- read_pixelator_3d_layout(layout_path, raw_component)
      validate(need(nrow(layout) > 0, "No 3D layout nodes are available for the selected component."))
      layout
    })

    colocalization_3d_markers <- reactive({
      selected <- input$colocalization_3d_markers
      selected <- as.character(selected)
      selected <- selected[!is.na(selected) & nzchar(selected)]
      if (length(selected) > 0) {
        return(selected)
      }

      layout <- colocalization_3d_layout_raw()
      head(sort(unique(layout$marker[layout$marker != "unlabeled"])), 4L)
    })

    colocalization_3d_nodes <- reactive({
      prepare_pixelator_3d_layout(
        colocalization_3d_layout_raw(),
        highlighted_markers = colocalization_3d_markers(),
        max_background_nodes = numeric_input_value(input$colocalization_3d_max_background, 7000)
      )
    })

    output$colocalization_3d_layout <- renderPlotly({
      nodes <- colocalization_3d_nodes()
      validate(need(nrow(nodes) > 0, "No 3D layout nodes are available for the selected component."))

      pixelator_3d_layout_plot(
        nodes,
        highlighted_markers = colocalization_3d_markers(),
        title = paste("3D layout:", input$colocalization_3d_component),
        dimensions = colocalization_3d_dimensions()
      )
    })

    output$colocalization_3d_component_table <- renderTable({
      current_data <- data()
      req(current_data, input$colocalization_3d_component)

      metadata <- current_data$metadata[current_data$metadata$component == input$colocalization_3d_component, , drop = FALSE]
      validate(need(nrow(metadata) > 0, "No metadata are available for the selected component."))
      cols <- intersect(c("sample", "sample_alias", "condition", "celltype_manual", "component", "n_umi", "n_edges"), names(metadata))
      metadata[1, cols, drop = FALSE]
    }, striped = TRUE, bordered = FALSE, width = "100%")

    output$colocalization_diff_summary <- renderUI({
      config <- colocalization_diff_config()
      req(config)
      result <- colocalization_diff_anchor_results()
      differential_summary_row(
        result,
        fdr_cutoff = config$fdr_cutoff,
        effect_cutoff = config$effect_cutoff
      )
    })

    colocalization_diff_volcano_x_label <- reactive({
      config <- colocalization_diff_config()
      req(config)
      paste("Difference in medians:", config$group_a, "minus", config$group_b, "(reference)")
    })

    colocalization_diff_volcano_dimensions <- reactive({
      apply_plot_options_overrides(
        differential_volcano_dimensions(colocalization_diff_volcano_x_label()),
        width_px = input$colocalization_diff_volcano_width,
        height_px = input$colocalization_diff_volcano_height
      )
    })

    colocalization_diff_volcano_ggplot <- reactive({
      config <- colocalization_diff_config()
      req(config)
      result <- colocalization_diff_anchor_results()
      validate(need(nrow(result) > 0, "Choose two different groups with enough colocalization data."))

      differential_volcano_ggplot(
        result,
        label_col = "marker_pair",
        x_label = colocalization_diff_volcano_x_label(),
        fdr_cutoff = config$fdr_cutoff,
        effect_cutoff = config$effect_cutoff
      )
    })

    output$colocalization_diff_volcano <- renderPlotly({
      dimensions <- colocalization_diff_volcano_dimensions()
      ggplotly(
        colocalization_diff_volcano_ggplot(),
        tooltip = "text",
        source = "colocalization_diff",
        width = dimensions$width,
        height = dimensions$height
      ) |>
        apply_differential_plot_frame(dimensions = dimensions)
    })
    register_ggplot_downloads(
      output,
      "colocalization_diff_volcano",
      colocalization_diff_volcano_ggplot,
      filename_prefix = function() paste("colocalization-differential-volcano", colocalization_diff_volcano_x_label(), sep = "-"),
      width = function() plot_download_size_from_dimensions(colocalization_diff_volcano_dimensions())$width,
      height = function() plot_download_size_from_dimensions(colocalization_diff_volcano_dimensions())$height
    )

    observeEvent(plotly::event_data("plotly_click", source = "colocalization_diff"), {
      event <- plotly::event_data("plotly_click", source = "colocalization_diff")
      if (!is.null(event$key) && nzchar(event$key)) {
        updateSelectInput(session, "colocalization_diff_pair", selected = event$key)
      }
    })

    colocalization_diff_detail_data <- reactive({
      current_data <- data()
      config <- colocalization_diff_config()
      req(current_data, config, input$colocalization_diff_pair, config$group_a, config$group_b)

      plot_data <- current_data$colocalization[
        current_data$colocalization$marker_pair == input$colocalization_diff_pair &
          current_data$colocalization$condition %in% c(config$group_a, config$group_b) &
          current_data$colocalization$celltype_manual %in% config$celltype_filter,
        ,
        drop = FALSE
      ]
      validate(need(nrow(plot_data) > 0, "No colocalization values are available for the selected pair and contrast."))

      plot_data$hover <- paste0(
        "Cell: ", plot_data$component,
        "<br>Condition: ", plot_data$condition,
        "<br>Cell type: ", plot_data$celltype_manual,
        "<br>Colocalization log2 ratio: ", round(plot_data$log2_ratio, 3)
      )

      y_label <- paste(input$colocalization_diff_pair, "colocalization log2 ratio")
      dimensions <- differential_detail_dimensions(
        plot_data,
        stratify_by_celltype = isTRUE(config$stratify_by_celltype),
        y_label = y_label
      )
      dimensions <- apply_plot_options_overrides(
        dimensions,
        width_px = input$colocalization_diff_detail_width,
        height_px = input$colocalization_diff_detail_height
      )

      list(
        config = config,
        plot_data = plot_data,
        y_label = y_label,
        dimensions = dimensions
      )
    })

    colocalization_diff_detail_ggplot <- reactive({
      detail <- colocalization_diff_detail_data()
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

    output$colocalization_diff_detail <- renderPlotly({
      dimensions <- colocalization_diff_detail_data()$dimensions
      ggplotly(colocalization_diff_detail_ggplot(), tooltip = "text", width = dimensions$width, height = dimensions$height) |>
        apply_proxiome_plot_frame(dimensions = dimensions)
    })
    register_ggplot_downloads(
      output,
      "colocalization_diff_detail",
      colocalization_diff_detail_ggplot,
      filename_prefix = function() paste("colocalization-differential-detail", input$colocalization_diff_pair %||% "pair", sep = "-"),
      width = function() plot_download_size_from_dimensions(colocalization_diff_detail_data()$dimensions)$width,
      height = function() plot_download_size_from_dimensions(colocalization_diff_detail_data()$dimensions)$height
    )

    output$colocalization_diff_table <- renderTable({
      config <- colocalization_diff_config()
      req(config)
      result <- filter_differential_hits(
        colocalization_diff_anchor_results(),
        fdr_cutoff = config$fdr_cutoff,
        effect_cutoff = config$effect_cutoff
      )
      validate(need(nrow(result) > 0, "No colocalization pairs pass the selected differential thresholds."))

      format_differential_table(result, effect_label = "diff_median_vs_reference")
    }, striped = TRUE, bordered = FALSE, width = "100%")
  })
}

colocalization_3d_sample_column <- function(metadata) {
  for (column in c("sample", "sample_id", "sample_alias")) {
    if (column %in% names(metadata)) {
      return(column)
    }
  }
  "component"
}

colocalization_3d_component_choices <- function(metadata, sample, cell_types) {
  if (!"component" %in% names(metadata) || nrow(metadata) == 0) {
    return(character(0))
  }

  sample_col <- colocalization_3d_sample_column(metadata)
  rows <- metadata
  sample <- as.character(sample)[1]
  if (!is.na(sample) && nzchar(sample) && sample_col %in% names(rows)) {
    rows <- rows[as.character(rows[[sample_col]]) == sample, , drop = FALSE]
  }

  cell_types <- as.character(cell_types)
  cell_types <- cell_types[!is.na(cell_types) & nzchar(cell_types)]
  if (length(cell_types) > 0 && "celltype_manual" %in% names(rows)) {
    rows <- rows[as.character(rows$celltype_manual) %in% cell_types, , drop = FALSE]
  }
  if (nrow(rows) == 0) {
    return(character(0))
  }

  score_col <- if ("n_umi" %in% names(rows)) "n_umi" else if ("n_edges" %in% names(rows)) "n_edges" else NULL
  if (!is.null(score_col)) {
    rows <- rows[order(abs(as.numeric(rows[[score_col]]) - 10000)), , drop = FALSE]
  }

  labels <- as.character(rows$component)
  if ("celltype_manual" %in% names(rows)) {
    labels <- paste(rows$celltype_manual, labels, sep = " | ")
  }
  if (!is.null(score_col)) {
    labels <- paste0(labels, " | ", score_col, ": ", format(as.numeric(rows[[score_col]]), big.mark = ",", trim = TRUE))
  }

  stats::setNames(as.character(rows$component), labels)
}
