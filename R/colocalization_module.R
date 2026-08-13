colocalization_sidebar <- function(id) {
  ns <- NS(id)

  sidebar(
    title = "Colocalization controls",
    width = 300,
    conditionalPanel(
      condition = "input.colocalization_mode == 'Observed'",
      ns = ns,
      accordion(
        open = c("Cell population", "Display"),
        accordion_panel(
          "Cell population",
          selectizeInput(ns("colocalization_celltype_filter"), "Cell population", choices = character(0), multiple = TRUE)
        ),
        accordion_panel(
          "Display",
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
              "Analysis group summary" = "condition",
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
            ns("colocalization_heatmap_view"),
            "Group display",
            choices = c("Focused group" = "focused", "Compare groups" = "compare"),
            selected = "focused"
          ),
          conditionalPanel(
            condition = "input.colocalization_heatmap_view == 'focused'",
            ns = ns,
            selectInput(ns("colocalization_heatmap_focus_group"), "Displayed group", choices = character(0))
          ),
          conditionalPanel(
            condition = "input.colocalization_heatmap_view == 'compare'",
            ns = ns,
            selectizeInput(ns("colocalization_heatmap_compare_groups"), "Groups to compare", choices = character(0), multiple = TRUE),
            helpText("Comparison uses a two-column grid. For the clearest labels, use 15 or fewer markers.")
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
          selectInput(ns("colocalization_reference_condition"), "Reference group", choices = character(0)),
          selectInput(
            ns("colocalization_clustering_method"),
            "Marker ordering",
            choices = c("Ward D2" = "ward.D2", "Complete" = "complete", "Average" = "average", "Single" = "single"),
            selected = "ward.D2"
          ),
          numericInput(ns("colocalization_legend_min"), "Legend minimum", value = -1, step = 0.1),
          numericInput(ns("colocalization_legend_max"), "Legend maximum", value = 1, step = 0.1),
          helpText("Report style applies 15 variable markers, a comparison grid, Ward D2 ordering, and a -0.75 to 0.75 legend."),
          conditionalPanel(
            condition = "input.colocalization_heatmap_preset == 'custom'",
            ns = ns,
            uiOutput(ns("colocalization_heatmap_settings_status")),
            actionButton(ns("apply_colocalization_heatmap"), "Apply heatmap settings", class = "btn-primary w-100")
          )
        ),
        accordion_panel(
          "Filters",
          selectizeInput(ns("colocalization_condition_filter"), "Analysis group", choices = character(0), multiple = TRUE),
          checkboxInput(ns("colocalization_pixelator_filter_enabled"), "Apply Pixelator proximity filters", value = FALSE),
          conditionalPanel(
            condition = "input.colocalization_pixelator_filter_enabled",
            ns = ns,
            numericInput(ns("colocalization_background_threshold_pct"), "Minimum marker fraction", value = 0.001, min = 0, max = 1, step = 0.0005),
            numericInput(ns("colocalization_background_threshold_count"), "Minimum marker count", value = 0, min = 0, step = 1),
            numericInput(ns("colocalization_min_cells_count"), "Minimum cells per pair", value = 1, min = 1, step = 1),
            helpText("Compatible with pixelatorR FilterProximityScores: both markers must meet the fraction and count thresholds; the pair must remain in the requested number of cells. Set a threshold to 0 to disable it.")
          )
        ),
        accordion_panel(
          "Advanced",
          selectInput(
            ns("colocalization_mean_type"),
            "Heatmap mean",
            choices = c(
              "Population mean" = "population",
              "Detected-cell mean" = "detected"
            ),
            selected = "population"
          ),
          helpText("Population mean includes cells without the pair as zero. Detected-cell mean averages only cells where the pair is recorded. Dot size always shows the detected fraction.")
        ),
        accordion_panel(
          "Interpretation",
          div(
            class = "small",
            p(class = "mb-2", tags$strong("Positive:"), " closer than expected by chance."),
            p(class = "mb-2", tags$strong("Zero:"), " approximately random spatial organization."),
            p(class = "mb-0", tags$strong("Negative:"), " spatial segregation.")
          )
        )
      )
    ),
    conditionalPanel(
      condition = "input.colocalization_mode == 'Differential'",
      ns = ns,
      accordion(
        open = c("Cell population", "Contrast", "Pair display"),
        accordion_panel(
          "Cell population",
          selectInput(ns("colocalization_diff_celltype_filter"), "Cell population", choices = character(0)),
          helpText("Differential colocalization compares one cell population at a time to avoid mixing population composition with spatial effects.")
        ),
        accordion_panel(
          "Contrast",
          selectInput(ns("colocalization_diff_group_a"), "Group A", choices = character(0)),
          selectInput(ns("colocalization_diff_group_b"), "Group B (reference)", choices = character(0)),
          selectInput(
            ns("colocalization_diff_mean_type"),
            "Sample statistic",
            choices = c(
              "Population mean" = "population",
              "Detected-cell mean" = "detected"
            ),
            selected = "population"
          ),
          numericInput(ns("colocalization_diff_min_samples"), "Minimum samples per group", value = 2, min = 2, step = 1),
          actionButton(ns("colocalization_run_differential"), "Run differential analysis", class = "btn-primary w-100")
        ),
        accordion_panel(
          "Pair display",
          selectInput(
            ns("colocalization_diff_pair_scope"),
            "Pairs shown",
            choices = c("All marker pairs" = "all", "Pairs containing one marker" = "anchor"),
            selected = "all"
          ),
          conditionalPanel(
            condition = "input.colocalization_diff_pair_scope == 'anchor'",
            ns = ns,
            selectInput(ns("colocalization_diff_anchor_marker"), "Marker", choices = character(0))
          ),
          selectizeInput(ns("colocalization_diff_pair"), "Detail pair", choices = character(0))
        ),
        accordion_panel(
          "Thresholds",
          numericInput(ns("colocalization_diff_fdr"), "FDR threshold", value = 0.05, min = 0, max = 1, step = 0.01),
          numericInput(ns("colocalization_diff_effect"), "Minimum median-sample difference", value = 0.25, min = 0, step = 0.05)
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

make_colocalization_heatmap_config <- function(
  markers,
  conditions,
  cell_types,
  preset = "custom"
) {
  reference_condition <- if ("CD3CD28" %in% conditions) "CD3CD28" else conditions[1]
  list(
    preset = preset,
    scope = "condition",
    focus_celltype = cell_types[1],
    view_mode = "focused",
    focus_group = reference_condition,
    comparison_groups = head(conditions, 6L),
    marker_selection_mode = "auto",
    markers = head(markers, min(20L, length(markers))),
    top_marker_count = 20L,
    min_pct_detected = 0.25,
    min_log2_range = 0.2,
    mean_type = "population",
    reference_condition = reference_condition,
    clustering_method = "ward.D2",
    legend_min = -1,
    legend_max = 1
  )
}

report_colocalization_heatmap_config <- function(config) {
  config$preset <- "report"
  config$view_mode <- "compare"
  config$marker_selection_mode <- "auto"
  config$top_marker_count <- 15L
  config$min_pct_detected <- 0.25
  config$min_log2_range <- 0.2
  config$mean_type <- "population"
  config$clustering_method <- "ward.D2"
  config$legend_min <- -0.75
  config$legend_max <- 0.75
  config
}

colocalization_heatmap_view_plan <- function(
  view_mode,
  marker_count,
  available_groups,
  focus_group = NULL,
  comparison_groups = character(),
  reference_group = NULL
) {
  available_groups <- unique(as.character(available_groups))
  available_groups <- available_groups[!is.na(available_groups) & nzchar(available_groups)]
  view_mode <- as.character(view_mode[1])
  if (length(view_mode) == 0 || is.na(view_mode) || !view_mode %in% c("focused", "compare")) {
    view_mode <- "focused"
  }

  marker_count <- suppressWarnings(as.integer(marker_count[1]))
  if (length(marker_count) == 0 || is.na(marker_count)) {
    marker_count <- 0L
  }
  density_notice <- NULL
  if (identical(view_mode, "compare") && marker_count > 20L) {
    view_mode <- "focused"
    density_notice <- "Comparison view supports up to 20 markers. Showing one focused group for readability."
  } else if (identical(view_mode, "compare") && marker_count > 15L) {
    density_notice <- "For clearer comparison labels, use 15 or fewer markers."
  }

  if (identical(view_mode, "compare")) {
    plot_groups <- intersect(as.character(comparison_groups), available_groups)
    if (length(plot_groups) == 0) {
      plot_groups <- head(available_groups, 6L)
    }
    facet_columns <- min(2L, length(plot_groups))
  } else {
    focus_group <- as.character(focus_group[1])
    if (length(focus_group) == 0 || is.na(focus_group) || !focus_group %in% available_groups) {
      reference_group <- as.character(reference_group[1])
      focus_group <- if (
        length(reference_group) > 0 &&
          !is.na(reference_group) &&
          reference_group %in% available_groups
      ) {
        reference_group
      } else {
        available_groups[1]
      }
    }
    plot_groups <- focus_group
    facet_columns <- 1L
  }

  list(
    view_mode = view_mode,
    plot_groups = plot_groups,
    facet_columns = facet_columns,
    density_notice = density_notice
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
          uiOutput(ns("colocalization_heatmap_notice")),
          plot_pane(
            size = "scroll",
            extra_class = "coloc-heatmap-pane",
            download_id = "colocalization_heatmap",
            ns = ns,
            controls = plot_options_controls(
              ns,
              "colocalization_heatmap_width",
              "colocalization_heatmap_height",
              width_value = 1000,
              height_value = 1200,
              max_value = 5000,
              show_view_dimensions = FALSE
            ),
            plotlyOutput(ns("colocalization_heatmap"), height = "auto")
          ),
          card(
            class = "mt-3",
            card_header("Pair detail"),
            card_body(
              p(class = "text-muted mb-2", "Click a heatmap point to inspect that marker pair by sample and cell population."),
              uiOutput(ns("colocalization_pair_detail_title")),
              uiOutput(ns("colocalization_pair_detail_metrics")),
              plot_pane(
                size = "wide",
                download_id = "colocalization_pair_detail",
                ns = ns,
                controls = plot_options_controls(
                  ns,
                  "colocalization_pair_detail_width",
                  "colocalization_pair_detail_height",
                  width_value = 900,
                  height_value = 520,
                  max_value = 5000,
                  show_view_dimensions = FALSE
                ),
                plotlyOutput(ns("colocalization_pair_detail"), height = "auto")
              ),
              div(class = "table-pane", tableOutput(ns("colocalization_pair_detail_table")))
            )
          ),
          div(class = "table-pane", tableOutput(ns("colocalization_table")))
        ),
        nav_panel(
          "Differential",
          uiOutput(ns("colocalization_diff_method")),
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
    colocalization_heatmap_config <- reactiveVal(NULL)
    custom_colocalization_heatmap_config <- reactiveVal(NULL)
    colocalization_heatmap_summary_config <- reactiveVal(NULL)
    colocalization_pair_selection <- reactiveVal(NULL)
    colocalization_pixelator_filter_config <- debounce(reactive(list(
      enabled = isTRUE(input$colocalization_pixelator_filter_enabled),
      background_threshold_pct = numeric_input_value(input$colocalization_background_threshold_pct, 0.001),
      background_threshold_count = numeric_input_value(input$colocalization_background_threshold_count, 0),
      min_cells_count = numeric_input_value(input$colocalization_min_cells_count, 1)
    )), 500)

    set_colocalization_heatmap_config <- function(config) {
      summary_config <- list(
        scope = config$scope,
        focus_celltype = if (identical(config$scope, "celltype")) config$focus_celltype else NULL,
        mean_type = config$mean_type %||% "population"
      )
      if (!identical(isolate(colocalization_heatmap_summary_config()), summary_config)) {
        colocalization_heatmap_summary_config(summary_config)
      }
      if (!identical(isolate(colocalization_heatmap_config()), config)) {
        colocalization_heatmap_config(config)
      }
      invisible(config)
    }

    update_colocalization_heatmap_inputs <- function(config) {
      updateSelectInput(session, "colocalization_heatmap_preset", selected = config$preset)
      updateSelectInput(session, "spatial_coloc_scope", selected = config$scope)
      updateSelectInput(session, "spatial_celltype_focus", selected = config$focus_celltype)
      updateSelectInput(session, "colocalization_heatmap_view", selected = config$view_mode)
      updateSelectInput(session, "colocalization_heatmap_focus_group", selected = config$focus_group)
      updateSelectizeInput(session, "colocalization_heatmap_compare_groups", selected = config$comparison_groups)
      updateSelectInput(session, "spatial_marker_selection_mode", selected = config$marker_selection_mode)
      updateSelectizeInput(session, "colocalization_heatmap_markers", selected = config$markers)
      updateNumericInput(session, "spatial_top_marker_count", value = config$top_marker_count)
      updateNumericInput(session, "spatial_min_pct_detected", value = config$min_pct_detected)
      updateNumericInput(session, "spatial_min_log2_range", value = config$min_log2_range)
      updateSelectInput(session, "colocalization_mean_type", selected = config$mean_type %||% "population")
      updateSelectInput(session, "colocalization_reference_condition", selected = config$reference_condition)
      updateSelectInput(session, "colocalization_clustering_method", selected = config$clustering_method)
      updateNumericInput(session, "colocalization_legend_min", value = config$legend_min)
      updateNumericInput(session, "colocalization_legend_max", value = config$legend_max)
      invisible(config)
    }

    colocalization_heatmap_config_from_inputs <- function() {
      current_config <- isolate(colocalization_heatmap_config())
      scope <- input$spatial_coloc_scope
      if (is.null(scope) || length(scope) != 1 || !scope %in% c("condition", "sample", "celltype")) {
        scope <- current_config$scope %||% "condition"
      }
      focus_celltype <- input$spatial_celltype_focus
      if (is.null(focus_celltype) || length(focus_celltype) != 1 || is.na(focus_celltype) || !nzchar(focus_celltype)) {
        focus_celltype <- current_config$focus_celltype
      }
      reference_condition <- input$colocalization_reference_condition
      if (is.null(reference_condition) || length(reference_condition) != 1 || is.na(reference_condition) || !nzchar(reference_condition)) {
        reference_condition <- current_config$reference_condition
      }
      view_mode <- input$colocalization_heatmap_view
      if (is.null(view_mode) || length(view_mode) != 1 || !view_mode %in% c("focused", "compare")) {
        view_mode <- current_config$view_mode %||% "focused"
      }
      mean_type <- input$colocalization_mean_type
      if (is.null(mean_type) || length(mean_type) != 1 || !mean_type %in% c("population", "detected")) {
        mean_type <- current_config$mean_type %||% "population"
      }

      list(
        preset = "custom",
        scope = scope,
        focus_celltype = focus_celltype,
        view_mode = view_mode,
        focus_group = input$colocalization_heatmap_focus_group %||% current_config$focus_group,
        comparison_groups = as.character(input$colocalization_heatmap_compare_groups %||% current_config$comparison_groups),
        marker_selection_mode = input$spatial_marker_selection_mode %||% "auto",
        markers = as.character(input$colocalization_heatmap_markers %||% character()),
        top_marker_count = numeric_input_value(input$spatial_top_marker_count, 20),
        min_pct_detected = numeric_input_value(input$spatial_min_pct_detected, 0.25),
        min_log2_range = numeric_input_value(input$spatial_min_log2_range, 0.2),
        mean_type = mean_type,
        reference_condition = reference_condition,
        clustering_method = input$colocalization_clustering_method %||% "ward.D2",
        legend_min = numeric_input_value(input$colocalization_legend_min, -1),
        legend_max = numeric_input_value(input$colocalization_legend_max, 1)
      )
    }

    colocalization_differential_config_from_inputs <- function() {
      config <- make_differential_config(
        group_a = input$colocalization_diff_group_a,
        group_b = input$colocalization_diff_group_b,
        celltype_filter = input$colocalization_diff_celltype_filter,
        stratify_by_celltype = FALSE,
        min_cells = numeric_input_value(input$colocalization_diff_min_samples, 2)
      )
      config$mean_type <- if (identical(input$colocalization_diff_mean_type, "detected")) "detected" else "population"
      config
    }

    observe({
      current_data <- data()
      req(current_data)

      conditions <- sort(unique(current_data$metadata$condition))
      cell_types <- sort(unique(current_data$metadata$celltype_manual))
      default_group_a <- conditions[1]
      default_group_b <- conditions[min(2, length(conditions))]
      pair_source <- current_data$colocalization_summary %||% current_data$colocalization
      colocalization_pairs <- sort(unique(pair_source$marker_pair))
      marker_source <- current_data$colocalization_sample_summary %||% pair_source
      colocalization_markers <- available_colocalization_marker_choices(marker_source)
      default_heatmap_markers <- head(colocalization_markers, min(20L, length(colocalization_markers)))
      default_colocalization_reference <- if ("CD3CD28" %in% conditions) "CD3CD28" else conditions[1]

      heatmap_config <- isolate(colocalization_heatmap_config())
      if (is.null(heatmap_config)) {
        heatmap_config <- make_colocalization_heatmap_config(colocalization_markers, conditions, cell_types)
      } else {
        heatmap_config$markers <- intersect(heatmap_config$markers, colocalization_markers)
        if (length(heatmap_config$markers) < 2) {
          heatmap_config$markers <- default_heatmap_markers
        }
        if (
          is.null(heatmap_config$focus_celltype) ||
            length(heatmap_config$focus_celltype) != 1 ||
            !heatmap_config$focus_celltype %in% cell_types
        ) {
          heatmap_config$focus_celltype <- cell_types[1]
        }
        if (
          is.null(heatmap_config$reference_condition) ||
            length(heatmap_config$reference_condition) != 1 ||
            !heatmap_config$reference_condition %in% conditions
        ) {
          heatmap_config$reference_condition <- default_colocalization_reference
        }
        heatmap_config$view_mode <- heatmap_config$view_mode %||% "focused"
        heatmap_config$focus_group <- heatmap_config$focus_group %||% default_colocalization_reference
        heatmap_config$comparison_groups <- heatmap_config$comparison_groups %||% head(conditions, 6L)
        heatmap_config$mean_type <- heatmap_config$mean_type %||% "population"
      }

      updateSelectizeInput(session, "colocalization_heatmap_markers", choices = colocalization_markers, selected = heatmap_config$markers)
      updateSelectInput(session, "spatial_celltype_focus", choices = cell_types, selected = heatmap_config$focus_celltype)
      updateSelectInput(
        session,
        "colocalization_reference_condition",
        choices = conditions,
        selected = heatmap_config$reference_condition
      )
      update_colocalization_heatmap_inputs(heatmap_config)
      set_colocalization_heatmap_config(heatmap_config)
      if (identical(heatmap_config$preset, "custom")) {
        custom_colocalization_heatmap_config(heatmap_config)
      } else if (is.null(isolate(custom_colocalization_heatmap_config()))) {
        custom_colocalization_heatmap_config(
          make_colocalization_heatmap_config(colocalization_markers, conditions, cell_types)
        )
      }

      updateSelectizeInput(session, "colocalization_condition_filter", choices = conditions, selected = conditions)
      updateSelectizeInput(session, "colocalization_celltype_filter", choices = cell_types, selected = cell_types)

      updateSelectInput(session, "colocalization_diff_group_a", choices = conditions, selected = default_group_a)
      updateSelectInput(session, "colocalization_diff_group_b", choices = conditions, selected = default_group_b)
      default_diff_celltype <- if ("CD8 T" %in% cell_types) "CD8 T" else cell_types[1]
      updateSelectInput(session, "colocalization_diff_celltype_filter", choices = cell_types, selected = default_diff_celltype)
      updateSelectInput(session, "colocalization_diff_anchor_marker", choices = current_data$marker_options, selected = current_data$marker_options[1])
      updateSelectizeInput(session, "colocalization_diff_pair", choices = colocalization_pairs, selected = colocalization_pairs[1], server = TRUE)

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

      default_diff_config <- default_differential_config(
        conditions,
        default_diff_celltype,
        min_cells = 2
      )
      default_diff_config$mean_type <- "population"
      colocalization_diff_config(default_diff_config)
    })

    observeEvent(input$colocalization_run_differential, {
      current_data <- data()
      req(current_data)
      colocalization_diff_config(colocalization_differential_config_from_inputs())
    }, ignoreInit = TRUE)

    observeEvent(input$colocalization_heatmap_preset, {
      current_config <- isolate(colocalization_heatmap_config())
      req(current_config)

      if (identical(input$colocalization_heatmap_preset, "report")) {
        if (identical(current_config$preset, "custom")) {
          custom_colocalization_heatmap_config(current_config)
        }
        new_config <- report_colocalization_heatmap_config(current_config)
      } else {
        new_config <- isolate(custom_colocalization_heatmap_config()) %||% current_config
        new_config$preset <- "custom"
      }

      update_colocalization_heatmap_inputs(new_config)
      set_colocalization_heatmap_config(new_config)
    }, ignoreInit = TRUE)

    observeEvent(input$apply_colocalization_heatmap, {
      new_config <- colocalization_heatmap_config_from_inputs()
      custom_colocalization_heatmap_config(new_config)
      updateSelectInput(session, "colocalization_heatmap_preset", selected = "custom")
      set_colocalization_heatmap_config(new_config)
    }, ignoreInit = TRUE)

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

    observe({
      current_data <- data()
      req(current_data)

      scope <- input$spatial_coloc_scope %||% "condition"
      metadata <- colocalization_metadata()
      if (identical(scope, "celltype") && !is.null(input$spatial_celltype_focus)) {
        metadata <- metadata[metadata$celltype_manual == input$spatial_celltype_focus, , drop = FALSE]
      }
      groups <- spatial_heatmap_group_values(metadata, scope)
      req(length(groups) > 0)

      active_config <- isolate(colocalization_heatmap_config())
      focus_group <- input$colocalization_heatmap_focus_group
      if (is.null(focus_group) || !focus_group %in% groups) {
        focus_group <- active_config$focus_group
      }
      if (is.null(focus_group) || !focus_group %in% groups) {
        focus_group <- groups[1]
      }

      comparison_groups <- intersect(input$colocalization_heatmap_compare_groups %||% character(), groups)
      if (length(comparison_groups) == 0) {
        comparison_groups <- intersect(active_config$comparison_groups %||% character(), groups)
      }
      if (length(comparison_groups) == 0) {
        comparison_groups <- head(groups, 6L)
      }

      updateSelectInput(session, "colocalization_heatmap_focus_group", choices = groups, selected = focus_group)
      updateSelectizeInput(
        session,
        "colocalization_heatmap_compare_groups",
        choices = groups,
        selected = comparison_groups
      )
    })

    output$colocalization_heatmap_settings_status <- renderUI({
      req(identical(input$colocalization_heatmap_preset, "custom"))
      active_config <- colocalization_heatmap_config()
      input_config <- colocalization_heatmap_config_from_inputs()
      fields <- c(
        "scope", "focus_celltype", "view_mode", "focus_group", "comparison_groups",
        "marker_selection_mode", "markers", "top_marker_count", "min_pct_detected",
        "min_log2_range", "mean_type", "reference_condition", "clustering_method", "legend_min", "legend_max"
      )
      pending <- !identical(active_config[fields], input_config[fields])

      if (pending) {
        return(div(class = "alert alert-warning py-2 mb-2", "Settings changed — click Apply to update the heatmap."))
      }
      div(class = "small text-muted mb-2", "Heatmap settings are applied.")
    })

    colocalization_diff_sample_data <- reactive({
      current_data <- data()
      config <- colocalization_diff_config()
      req(current_data, config, config$group_a, config$group_b)

      sample_summary <- current_data$colocalization_sample_summary %||%
        summarize_colocalization_by_sample(current_data$colocalization)
      prepare_colocalization_differential_samples(
        sample_summary,
        current_data$metadata,
        conditions = c(config$group_a, config$group_b),
        celltype = config$celltype_filter[1],
        mean_type = config$mean_type %||% "population"
      )
    })

    colocalization_diff_results <- reactive({
      config <- colocalization_diff_config()
      req(config, config$group_a, config$group_b)

      calculate_differential_readout(
        colocalization_diff_sample_data(),
        feature_cols = c("marker_pair", "marker_1", "marker_2"),
        value_col = "mean_log2_ratio",
        group_a = config$group_a,
        group_b = config$group_b,
        celltype_filter = config$celltype_filter,
        stratify_by_celltype = FALSE,
        min_cells = config$min_cells,
        fdr_cutoff = config$fdr_cutoff
      )
    })

    colocalization_diff_display_results <- reactive({
      result <- colocalization_diff_results()
      if (!identical(input$colocalization_diff_pair_scope, "anchor")) {
        return(result)
      }

      anchor <- input$colocalization_diff_anchor_marker
      result[result$marker_1 == anchor | result$marker_2 == anchor, , drop = FALSE]
    })

    colocalization_diff_thresholds <- reactive(list(
      fdr = numeric_input_value(input$colocalization_diff_fdr, 0.05),
      effect = numeric_input_value(input$colocalization_diff_effect, 0.25)
    ))

    observe({
      choices <- sort(unique(as.character(colocalization_diff_display_results()$marker_pair)))
      selected <- input$colocalization_diff_pair
      if (is.null(selected) || !selected %in% choices) {
        selected <- choices[1]
      }
      updateSelectizeInput(session, "colocalization_diff_pair", choices = choices, selected = selected, server = TRUE)
    })

    colocalization_all_marker_summary <- reactive({
      current_data <- data()
      summary_config <- colocalization_heatmap_summary_config()
      pixelator_filter_config <- colocalization_pixelator_filter_config()
      req(current_data, summary_config, pixelator_filter_config)

      metadata <- colocalization_metadata()
      validate(need(nrow(metadata) > 0, "No cells are available for the selected filters."))

      scope <- summary_config$scope
      if (!scope %in% c("condition", "sample", "celltype")) {
        scope <- "condition"
      }

      if (identical(scope, "celltype")) {
        focus_celltype <- summary_config$focus_celltype
        if (is.null(focus_celltype) || length(focus_celltype) == 0 || is.na(focus_celltype) || !nzchar(focus_celltype)) {
          focus_celltype <- sort(unique(as.character(metadata$celltype_manual)))[1]
        }
        metadata <- metadata[metadata$celltype_manual == focus_celltype, , drop = FALSE]
        validate(need(nrow(metadata) > 0, "No cells are available for the selected cell type focus."))
      }

      sample_summary <- current_data$colocalization_sample_summary
      if (isTRUE(pixelator_filter_config$enabled)) {
        filtered_colocalization <- tryCatch(
          filter_pixelator_proximity_scores(
            current_data$colocalization,
            metadata,
            current_data$abundance,
            background_threshold_pct = pixelator_filter_config$background_threshold_pct,
            background_threshold_count = pixelator_filter_config$background_threshold_count,
            min_cells_count = pixelator_filter_config$min_cells_count
          ),
          error = function(error) validate(need(FALSE, conditionMessage(error)))
        )
        validate(need(nrow(filtered_colocalization) > 0, "No colocalization scores pass the Pixelator filters."))
        sample_summary <- summarize_colocalization_by_sample(filtered_colocalization)
      } else if (is.null(sample_summary)) {
        sample_summary <- summarize_colocalization_by_sample(current_data$colocalization)
      }
      available_markers <- available_colocalization_marker_choices(sample_summary)
      validate(need(length(available_markers) >= 2, "No colocalization scores are available for the selected filters."))

      spatial_heatmap_summary_for_scope(
        selected_markers = available_markers,
        scope = scope,
        sample_summary = sample_summary,
        metadata = metadata,
        include_missing_obs = !identical(summary_config$mean_type, "detected")
      )
    })

    colocalization_heatmap_result <- reactive({
      config <- colocalization_heatmap_config()
      marker_summary <- colocalization_all_marker_summary()
      req(config)
      validate(need(nrow(marker_summary) > 0, "No colocalization scores are available for the selected filters."))

      scope <- config$scope
      marker_selection_mode <- config$marker_selection_mode
      available_markers <- available_colocalization_marker_choices(marker_summary)
      requested_markers <- selected_or_all(config$markers, available_markers)
      requested_markers <- intersect(requested_markers, available_markers)
      candidate_markers <- if (identical(marker_selection_mode, "auto")) available_markers else requested_markers
      validate(need(length(candidate_markers) >= 2, "Select at least two markers for the spatial heatmap."))

      candidate_summary <- marker_summary[
        marker_summary$marker_1 %in% candidate_markers &
          marker_summary$marker_2 %in% candidate_markers,
        ,
        drop = FALSE
      ]
      selected_markers <- spatial_heatmap_selected_markers(
        summary = candidate_summary,
        available_markers = candidate_markers,
        requested_markers = requested_markers,
        marker_selection_mode = marker_selection_mode,
        n_markers = config$top_marker_count,
        min_pct_detected = config$min_pct_detected,
        min_range = config$min_log2_range
      )
      validate(need(length(selected_markers) >= 2, "Select at least two markers for the spatial heatmap."))
      validate(need(length(selected_markers) <= 40, "Use 40 or fewer markers for an interpretable spatial heatmap."))

      summary <- marker_summary[
        marker_summary$marker_1 %in% selected_markers &
          marker_summary$marker_2 %in% selected_markers,
        ,
        drop = FALSE
      ]
      validate(need(nrow(summary) > 0, "No colocalization scores are available for the selected markers."))

      group_cols <- spatial_heatmap_group_cols(scope)
      summary <- complete_spatial_marker_pairs(
        summary = summary,
        selected_markers = selected_markers,
        group_cols = group_cols
      )
      condition_col <- if (identical(scope, "condition")) "condition" else "sample_alias"
      available_plot_groups <- unique(as.character(summary[[condition_col]]))
      available_plot_groups <- available_plot_groups[nzchar(available_plot_groups)]
      validate(need(length(available_plot_groups) > 0, "No spatial heatmap groups are available for the selected filters."))
      selected_conditions <- unique(as.character(summary$condition))
      reference_condition <- config$reference_condition
      if (
        is.null(reference_condition) ||
          length(reference_condition) != 1 ||
          !reference_condition %in% selected_conditions
      ) {
        reference_condition <- selected_conditions[1]
      }

      view_plan <- colocalization_heatmap_view_plan(
        view_mode = config$view_mode %||% "focused",
        marker_count = length(selected_markers),
        available_groups = available_plot_groups,
        focus_group = config$focus_group,
        comparison_groups = config$comparison_groups %||% character(),
        reference_group = if (identical(scope, "condition")) reference_condition else NULL
      )
      plot_groups <- view_plan$plot_groups

      result <- make_coloc_heatmaps(
        data = summary,
        selected_markers = selected_markers,
        cell_label = spatial_heatmap_cell_label(
          scope,
          selected_celltypes = input$colocalization_celltype_filter,
          focus_celltype = config$focus_celltype
        ),
        conditions = plot_groups,
        reference_condition = if (identical(scope, "condition")) {
          reference_condition
        } else {
          plot_groups[1]
        },
        condition_col = condition_col,
        clustering_method = config$clustering_method,
        facet_columns = view_plan$facet_columns,
        value_label = colocalization_mean_label(config$mean_type),
        legend_range = colocalization_legend_range(
          config$legend_min,
          config$legend_max
        )
      )
      result$summary <- summary[summary[[condition_col]] %in% plot_groups, , drop = FALSE]
      result$scope <- scope
      result$condition_col <- condition_col
      result$view_mode <- view_plan$view_mode
      result$density_notice <- view_plan$density_notice
      result
    })

    output$colocalization_heatmap_notice <- renderUI({
      notice <- colocalization_heatmap_result()$density_notice
      if (is.null(notice) || !nzchar(notice)) {
        return(NULL)
      }
      div(class = "alert alert-warning py-2 mb-2", notice)
    })

    output$colocalization_heatmap <- renderPlotly({
      coloc_heatmap_plotly(
        colocalization_heatmap_result(),
        dimensions = plotly_display_dimensions(colocalization_heatmap_dimensions()),
        source = "colocalization_heatmap"
      )
    })

    observeEvent(plotly::event_data("plotly_click", source = "colocalization_heatmap"), {
      event <- plotly::event_data("plotly_click", source = "colocalization_heatmap")
      result <- isolate(colocalization_heatmap_result())
      clicked <- result$plot_data[
        as.character(result$plot_data$pair_key) == as.character(event$key),
        ,
        drop = FALSE
      ]
      if (nrow(clicked) > 0) {
        colocalization_pair_selection(list(
          marker_1 = as.character(clicked$marker_1[1]),
          marker_2 = as.character(clicked$marker_2[1]),
          group = as.character(clicked[[result$condition_col]][1])
        ))
      }
    }, ignoreInit = TRUE)

    selected_colocalization_pair <- reactive({
      result <- colocalization_heatmap_result()
      resolve_colocalization_pair_selection(
        result$plot_data,
        selection = colocalization_pair_selection(),
        condition_col = result$condition_col
      )
    })

    colocalization_pair_detail_data <- reactive({
      current_data <- data()
      result <- colocalization_heatmap_result()
      selection <- selected_colocalization_pair()
      filter_config <- colocalization_pixelator_filter_config()
      req(current_data, selection, filter_config)

      metadata <- colocalization_metadata()
      if (identical(result$scope, "celltype")) {
        focus_celltype <- colocalization_heatmap_config()$focus_celltype
        metadata <- metadata[metadata$celltype_manual == focus_celltype, , drop = FALSE]
      }
      displayed_groups <- unique(as.character(result$summary[[result$condition_col]]))
      metadata <- metadata[
        as.character(metadata[[result$condition_col]]) %in% displayed_groups,
        ,
        drop = FALSE
      ]
      validate(need(nrow(metadata) > 0, "No cells are available for the selected pair detail."))

      sample_summary <- current_data$colocalization_sample_summary
      if (isTRUE(filter_config$enabled) || is.null(sample_summary)) {
        pair_scores <- filter_colocalization_marker_pair(
          current_data$colocalization,
          selection$marker_1,
          selection$marker_2
        )
        if (isTRUE(filter_config$enabled) && nrow(pair_scores) > 0) {
          pair_scores <- tryCatch(
            filter_pixelator_proximity_scores(
              pair_scores,
              metadata,
              current_data$abundance,
              background_threshold_pct = filter_config$background_threshold_pct,
              background_threshold_count = filter_config$background_threshold_count,
              min_cells_count = filter_config$min_cells_count
            ),
            error = function(error) validate(need(FALSE, conditionMessage(error)))
          )
        }
        sample_summary <- summarize_colocalization_by_sample(pair_scores)
      }

      detail <- summarize_colocalization_pair_detail(
        sample_summary,
        metadata,
        marker_1 = selection$marker_1,
        marker_2 = selection$marker_2,
        mean_type = colocalization_heatmap_config()$mean_type
      )
      validate(need(nrow(detail) > 0, "No sample or cell-population detail is available for this pair."))

      list(
        summary = detail,
        selection = selection,
        value_label = result$value_label,
        legend_range = colocalization_legend_range(
          colocalization_heatmap_config()$legend_min,
          colocalization_heatmap_config()$legend_max
        )
      )
    })

    output$colocalization_pair_detail_title <- renderUI({
      selection <- selected_colocalization_pair()
      tags$h4(class = "mb-3", paste(selection$marker_1, "↔", selection$marker_2))
    })

    output$colocalization_pair_detail_metrics <- renderUI({
      selection <- selected_colocalization_pair()
      result <- colocalization_heatmap_result()
      group_label <- if (identical(result$condition_col, "sample_alias")) "Selected sample" else "Selected analysis group"
      score <- if (isTRUE(selection$pair_observed)) {
        format(round(selection$plot_value, 3), trim = TRUE)
      } else {
        "Not available"
      }

      metric_row(
        metric_tile(group_label, selection$group),
        metric_tile(result$value_label, score),
        metric_tile("Detected cells", paste0(format(selection$n_detected, big.mark = ","), " / ", format(selection$n_total, big.mark = ","))),
        metric_tile("Detected fraction", format_percent(selection$plot_size))
      )
    })

    colocalization_pair_detail_ggplot <- reactive({
      detail <- colocalization_pair_detail_data()
      plot_colocalization_pair_detail(
        detail$summary,
        marker_1 = detail$selection$marker_1,
        marker_2 = detail$selection$marker_2,
        value_label = detail$value_label,
        legend_range = detail$legend_range
      )
    })
    colocalization_pair_detail_dimensions <- reactive({
      plot_options_view_overrides(
        input,
        "colocalization_pair_detail",
        colocalization_pair_detail_base_dimensions(colocalization_pair_detail_data()$summary)
      )
    })
    colocalization_pair_detail_export_dimensions <- reactive({
      plot_options_export_overrides(
        input,
        "colocalization_pair_detail",
        colocalization_pair_detail_base_dimensions(colocalization_pair_detail_data()$summary)
      )
    })

    output$colocalization_pair_detail <- renderPlotly({
      dimensions <- plotly_display_dimensions(colocalization_pair_detail_dimensions())
      ggplotly(
        colocalization_pair_detail_ggplot(),
        tooltip = "text",
        width = dimensions$width,
        height = dimensions$height
      ) |>
        apply_proxiome_plot_frame(
          colorbar_title = colocalization_pair_detail_data()$value_label,
          dimensions = dimensions
        )
    })
    register_ggplot_downloads(
      output,
      "colocalization_pair_detail",
      colocalization_pair_detail_ggplot,
      filename_prefix = function() {
        detail <- colocalization_pair_detail_data()
        paste(
          "colocalization-pair-detail",
          detail$selection$marker_1,
          detail$selection$marker_2,
          detail$value_label,
          sep = "-"
        )
      },
      width = function() plot_download_size_from_dimensions(colocalization_pair_detail_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(colocalization_pair_detail_export_dimensions())$height
    )

    output$colocalization_pair_detail_table <- renderTable({
      detail <- colocalization_pair_detail_data()
      format_colocalization_pair_detail_table(detail$summary, detail$value_label)
    }, striped = TRUE, bordered = FALSE, width = "100%")

    colocalization_heatmap_ggplot <- reactive({
      colocalization_heatmap_result()$plot
    })
    colocalization_heatmap_dimensions <- reactive({
      result <- colocalization_heatmap_result()
      plot_options_view_overrides(
        input,
        "colocalization_heatmap",
        coloc_heatmap_widget_dimensions(result$plot_data, facet_columns = result$facet_columns)
      )
    })
    colocalization_heatmap_export_dimensions <- reactive({
      result <- colocalization_heatmap_result()
      plot_options_export_overrides(
        input,
        "colocalization_heatmap",
        coloc_heatmap_widget_dimensions(result$plot_data, facet_columns = result$facet_columns)
      )
    })
    register_ggplot_downloads(
      output,
      "colocalization_heatmap",
      colocalization_heatmap_ggplot,
      filename_prefix = function() paste("colocalization-heatmap", colocalization_heatmap_result()$value_label, sep = "-"),
      width = function() plot_download_size_from_dimensions(colocalization_heatmap_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(colocalization_heatmap_export_dimensions())$height
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
      dimensions <- plotly_display_dimensions(colocalization_3d_dimensions())
      validate(need(nrow(nodes) > 0, "No 3D layout nodes are available for the selected component."))

      pixelator_3d_layout_plot(
        nodes,
        highlighted_markers = colocalization_3d_markers(),
        title = paste("3D layout:", input$colocalization_3d_component),
        dimensions = dimensions
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
      thresholds <- colocalization_diff_thresholds()
      result <- colocalization_diff_display_results()
      differential_summary_row(
        result,
        fdr_cutoff = thresholds$fdr,
        effect_cutoff = thresholds$effect
      )
    })

    output$colocalization_diff_method <- renderUI({
      config <- colocalization_diff_config()
      sample_data <- colocalization_diff_sample_data()
      req(config)
      sample_counts <- vapply(c(config$group_a, config$group_b), function(group) {
        length(unique(sample_data$sample_alias[sample_data$condition == group]))
      }, integer(1))
      replicated <- all(sample_counts >= config$min_cells)

      div(
        class = paste("alert py-2 mb-3", if (replicated) "alert-info" else "alert-warning"),
        tags$strong("Sample-aware analysis. "),
        paste0(
          colocalization_mean_label(config$mean_type),
          " is calculated within each sample; effects compare median sample values. ",
          config$group_a, ": ", sample_counts[1], " sample(s); ",
          config$group_b, ": ", sample_counts[2], " sample(s)."
        ),
        if (!replicated) {
          paste0(" At least ", config$min_cells, " samples per group are required for p-values and FDR; effects remain descriptive.")
        }
      )
    })

    colocalization_diff_volcano_x_label <- reactive({
      config <- colocalization_diff_config()
      req(config)
      paste(
        "Median difference in sample-level",
        paste0(tolower(colocalization_mean_label(config$mean_type)), ":"),
        config$group_a, "minus", config$group_b, "(reference)"
      )
    })

    colocalization_diff_volcano_dimensions <- reactive({
      plot_options_view_overrides(input, "colocalization_diff_volcano", differential_volcano_dimensions(colocalization_diff_volcano_x_label()))
    })
    colocalization_diff_volcano_export_dimensions <- reactive({
      plot_options_export_overrides(input, "colocalization_diff_volcano", differential_volcano_dimensions(colocalization_diff_volcano_x_label()))
    })

    colocalization_diff_volcano_ggplot <- reactive({
      config <- colocalization_diff_config()
      thresholds <- colocalization_diff_thresholds()
      req(config)
      result <- colocalization_diff_display_results()
      validate(need(nrow(result) > 0, "Choose two different groups with sample-level colocalization data."))

      differential_volcano_ggplot(
        result,
        label_col = "marker_pair",
        x_label = colocalization_diff_volcano_x_label(),
        fdr_cutoff = thresholds$fdr,
        effect_cutoff = thresholds$effect
      )
    })

    output$colocalization_diff_volcano <- renderPlotly({
      dimensions <- colocalization_diff_volcano_dimensions()
      display_dimensions <- responsive_plotly_dimensions(dimensions)
      ggplotly(
        colocalization_diff_volcano_ggplot(),
        tooltip = "text",
        source = "colocalization_diff",
        width = display_dimensions$width,
        height = display_dimensions$height
      ) |>
        apply_differential_plot_frame(dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "colocalization_diff_volcano",
      colocalization_diff_volcano_ggplot,
      filename_prefix = function() paste("colocalization-differential-volcano", colocalization_diff_volcano_x_label(), sep = "-"),
      width = function() plot_download_size_from_dimensions(colocalization_diff_volcano_export_dimensions())$width,
      height = function() plot_download_size_from_dimensions(colocalization_diff_volcano_export_dimensions())$height
    )

    observeEvent(plotly::event_data("plotly_click", source = "colocalization_diff"), {
      event <- plotly::event_data("plotly_click", source = "colocalization_diff")
      if (!is.null(event$key) && nzchar(event$key)) {
        updateSelectizeInput(session, "colocalization_diff_pair", selected = event$key, server = TRUE)
      }
    })

    colocalization_diff_detail_data <- reactive({
      config <- colocalization_diff_config()
      req(config, input$colocalization_diff_pair, config$group_a, config$group_b)

      plot_data <- colocalization_diff_sample_data()
      plot_data <- plot_data[
        plot_data$marker_pair == input$colocalization_diff_pair,
        ,
        drop = FALSE
      ]
      validate(need(nrow(plot_data) > 0, "No sample-level values are available for the selected pair and contrast."))

      plot_data$hover <- paste0(
        "Sample: ", plot_data$sample_alias,
        "<br>Analysis group: ", plot_data$condition,
        "<br>Cell population: ", plot_data$celltype_manual,
        "<br>", colocalization_mean_label(config$mean_type), ": ", round(plot_data$mean_log2_ratio, 3),
        "<br>Pair status: ", ifelse(plot_data$pair_observed, "Detected pair", "No detected pair"),
        "<br>Detected cells: ", plot_data$n_detected, " / ", plot_data$n_total,
        "<br>Detected fraction: ", format_percent(plot_data$pct_detected)
      )
      plot_data$display_value <- ifelse(plot_data$pair_observed, plot_data$mean_log2_ratio, 0)

      y_label <- paste(input$colocalization_diff_pair, colocalization_mean_label(config$mean_type))
      base_dimensions <- differential_detail_dimensions(
        plot_data,
        stratify_by_celltype = FALSE,
        y_label = y_label
      )
      dimensions <- plot_options_view_overrides(input, "colocalization_diff_detail", base_dimensions)
      export_dimensions <- plot_options_export_overrides(input, "colocalization_diff_detail", base_dimensions)

      list(
        config = config,
        plot_data = plot_data,
        y_label = y_label,
        dimensions = dimensions,
        export_dimensions = export_dimensions
      )
    })

    colocalization_diff_detail_ggplot <- reactive({
      detail <- colocalization_diff_detail_data()
      plot_data <- detail$plot_data
      observed <- plot_data[plot_data$pair_observed, , drop = FALSE]
      missing <- plot_data[!plot_data$pair_observed, , drop = FALSE]

      ggplot(plot_data, aes(condition, color = condition, text = hover)) +
        geom_hline(yintercept = 0, color = "#8a9699", linewidth = 0.5) +
        geom_boxplot(aes(y = mean_log2_ratio), outlier.shape = NA, alpha = 0.18, linewidth = 0.5, na.rm = TRUE) +
        geom_jitter(data = observed, aes(y = mean_log2_ratio), width = 0.1, height = 0, alpha = 0.8, size = 2.8) +
        geom_point(data = missing, aes(y = display_value), shape = 4, color = "#8a9493", size = 2.8, stroke = 0.9) +
        labs(x = NULL, y = detail$y_label, caption = "Each point = one sample; × = no detected pair") +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          plot.caption = element_text(color = "#5f6b6a", hjust = 0),
          legend.position = "none"
        )
    })

    output$colocalization_diff_detail <- renderPlotly({
      dimensions <- colocalization_diff_detail_data()$dimensions
      display_dimensions <- responsive_plotly_dimensions(dimensions)
      ggplotly(colocalization_diff_detail_ggplot(), tooltip = "text", width = display_dimensions$width, height = display_dimensions$height) |>
        apply_proxiome_plot_frame(dimensions = display_dimensions)
    })
    register_ggplot_downloads(
      output,
      "colocalization_diff_detail",
      colocalization_diff_detail_ggplot,
      filename_prefix = function() paste(
        "colocalization-differential-detail",
        input$colocalization_diff_pair %||% "pair",
        colocalization_mean_label(colocalization_diff_config()$mean_type),
        sep = "-"
      ),
      width = function() plot_download_size_from_dimensions(colocalization_diff_detail_data()$export_dimensions)$width,
      height = function() plot_download_size_from_dimensions(colocalization_diff_detail_data()$export_dimensions)$height
    )

    output$colocalization_diff_table <- renderTable({
      config <- colocalization_diff_config()
      thresholds <- colocalization_diff_thresholds()
      req(config)
      result <- colocalization_diff_display_results()
      validate(need(nrow(result) > 0, "No marker pairs are available for the selected differential comparison."))
      result$threshold_status <- differential_threshold_status(
        result,
        fdr_cutoff = thresholds$fdr,
        effect_cutoff = thresholds$effect
      )
      result$is_significant <- NULL

      result <- format_differential_table(
        result,
        effect_label = paste0("median_sample_", config$mean_type, "_difference")
      )
      names(result)[names(result) == "n_a"] <- "samples_a"
      names(result)[names(result) == "n_b"] <- "samples_b"
      result
    }, striped = TRUE, bordered = FALSE, width = "100%")
  })
}

prepare_colocalization_differential_samples <- function(
  sample_summary,
  metadata,
  conditions,
  celltype,
  mean_type = "population"
) {
  metadata_cols <- c("component", "sample_alias", "condition", "celltype_manual")
  summary_cols <- c(
    "sample_alias", "celltype_manual", "marker_1", "marker_2",
    "sum_log2_ratio", "n_detected"
  )
  missing_cols <- unique(c(
    setdiff(metadata_cols, names(metadata)),
    setdiff(summary_cols, names(sample_summary))
  ))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns for sample-aware differential colocalization: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  require_namespace("data.table")
  celltype <- as.character(celltype)[1]
  conditions <- unique(as.character(conditions))
  metadata_dt <- unique(data.table::as.data.table(metadata)[, ..metadata_cols], by = "component")
  metadata_dt <- metadata_dt[
    as.character(condition) %in% conditions &
      as.character(celltype_manual) == celltype
  ]
  totals <- metadata_dt[, .(n_total = data.table::uniqueN(component)), keyby = .(sample_alias, condition, celltype_manual)]

  empty_result <- data.frame(
    sample_alias = character(),
    condition = character(),
    celltype_manual = character(),
    marker_pair = character(),
    marker_1 = character(),
    marker_2 = character(),
    mean_log2_ratio = numeric(),
    pct_detected = numeric(),
    n_detected = integer(),
    n_total = integer(),
    pair_observed = logical(),
    stringsAsFactors = FALSE
  )
  if (nrow(totals) == 0) {
    return(empty_result)
  }

  summary_dt <- data.table::as.data.table(sample_summary[, summary_cols, drop = FALSE])[
    as.character(sample_alias) %in% totals$sample_alias &
      as.character(celltype_manual) == celltype &
      as.character(marker_1) != as.character(marker_2)
  ]
  if (nrow(summary_dt) == 0) {
    return(empty_result)
  }

  summary_dt[, original_marker_1 := as.character(marker_1)]
  summary_dt[, original_marker_2 := as.character(marker_2)]
  summary_dt[, marker_1 := ifelse(original_marker_1 <= original_marker_2, original_marker_1, original_marker_2)]
  summary_dt[, marker_2 := ifelse(original_marker_1 <= original_marker_2, original_marker_2, original_marker_1)]
  summary_dt[, preferred_orientation := original_marker_1 == marker_1]
  data.table::setorderv(
    summary_dt,
    c("sample_alias", "celltype_manual", "marker_1", "marker_2", "preferred_orientation"),
    c(1L, 1L, 1L, 1L, -1L)
  )
  summary_dt <- unique(
    summary_dt,
    by = c("sample_alias", "celltype_manual", "marker_1", "marker_2")
  )
  summary_dt[, c("original_marker_1", "original_marker_2", "preferred_orientation") := NULL]

  pairs <- unique(summary_dt[, .(marker_1, marker_2)])
  totals[, cross_key := 1L]
  pairs[, cross_key := 1L]
  grid <- pairs[totals, on = "cross_key", allow.cartesian = TRUE]
  grid[, cross_key := NULL]
  detail <- summary_dt[grid, on = c("sample_alias", "celltype_manual", "marker_1", "marker_2")]
  detail[is.na(sum_log2_ratio), sum_log2_ratio := 0]
  detail[is.na(n_detected), n_detected := 0L]
  detail[, n_detected := as.integer(n_detected)]
  detail[, pair_observed := n_detected > 0]
  denominator <- if (identical(mean_type, "detected")) detail$n_detected else detail$n_total
  detail[, mean_log2_ratio := ifelse(denominator > 0, sum_log2_ratio / denominator, NA_real_)]
  detail[, pct_detected := ifelse(n_total > 0, n_detected / n_total, NA_real_)]
  detail[, marker_pair := paste(marker_1, marker_2, sep = " / ")]
  data.table::setorderv(detail, c("marker_pair", "condition", "sample_alias"))

  as.data.frame(detail[, .(
    sample_alias,
    condition,
    celltype_manual,
    marker_pair,
    marker_1,
    marker_2,
    mean_log2_ratio,
    pct_detected,
    n_detected,
    n_total,
    pair_observed
  )])
}

resolve_colocalization_pair_selection <- function(plot_data, selection = NULL, condition_col = "condition") {
  required_cols <- c(
    condition_col, "marker_1", "marker_2", "plot_value", "plot_size",
    "pair_observed", "n_detected", "n_total"
  )
  if (!is.data.frame(plot_data) || nrow(plot_data) == 0 || !all(required_cols %in% names(plot_data))) {
    return(NULL)
  }

  selected_rows <- plot_data[FALSE, , drop = FALSE]
  if (!is.null(selection) && all(c("marker_1", "marker_2") %in% names(selection))) {
    selected_rows <- plot_data[
      as.character(plot_data$marker_1) == as.character(selection$marker_1)[1] &
        as.character(plot_data$marker_2) == as.character(selection$marker_2)[1],
      ,
      drop = FALSE
    ]
    if (nrow(selected_rows) > 0 && "group" %in% names(selection)) {
      group_rows <- selected_rows[
        as.character(selected_rows[[condition_col]]) == as.character(selection$group)[1],
        ,
        drop = FALSE
      ]
      if (nrow(group_rows) > 0) {
        selected_rows <- group_rows
      }
    }
  }

  if (nrow(selected_rows) == 0) {
    selected_rows <- plot_data[plot_data$pair_observed, , drop = FALSE]
    if (nrow(selected_rows) == 0) {
      selected_rows <- plot_data
    }
    selected_rows <- selected_rows[
      order(abs(selected_rows$plot_value), selected_rows$plot_size, decreasing = TRUE, na.last = TRUE),
      ,
      drop = FALSE
    ]
  }

  selected <- as.list(selected_rows[1, , drop = FALSE])
  selected$marker_1 <- as.character(selected$marker_1)
  selected$marker_2 <- as.character(selected$marker_2)
  selected$group <- as.character(selected_rows[[condition_col]][1])
  selected
}

filter_colocalization_marker_pair <- function(data, marker_1, marker_2) {
  if (!is.data.frame(data) || nrow(data) == 0 || !all(c("marker_1", "marker_2") %in% names(data))) {
    return(data)
  }

  marker_1_values <- as.character(data$marker_1)
  marker_2_values <- as.character(data$marker_2)
  keep <- (marker_1_values == marker_1 & marker_2_values == marker_2) |
    (marker_1_values == marker_2 & marker_2_values == marker_1)
  data[keep, , drop = FALSE]
}

summarize_colocalization_pair_detail <- function(
  sample_summary,
  metadata,
  marker_1,
  marker_2,
  mean_type = "population"
) {
  metadata_cols <- c("component", "sample_alias", "condition", "celltype_manual")
  summary_cols <- c(
    "sample_alias", "celltype_manual", "marker_1", "marker_2",
    "sum_log2_ratio", "n_detected"
  )
  missing_metadata_cols <- setdiff(metadata_cols, names(metadata))
  missing_summary_cols <- setdiff(summary_cols, names(sample_summary))
  if (length(missing_metadata_cols) > 0 || length(missing_summary_cols) > 0) {
    stop(
      "Missing columns for colocalization pair detail: ",
      paste(unique(c(missing_metadata_cols, missing_summary_cols)), collapse = ", "),
      call. = FALSE
    )
  }

  require_namespace("data.table")
  metadata_dt <- unique(data.table::as.data.table(metadata)[, ..metadata_cols], by = "component")
  totals <- metadata_dt[, .(n_total = data.table::uniqueN(component)), keyby = .(condition, sample_alias, celltype_manual)]
  if (nrow(totals) == 0) {
    return(data.frame())
  }

  pair_rows <- filter_colocalization_marker_pair(sample_summary, marker_1, marker_2)
  pair_dt <- data.table::as.data.table(pair_rows[, summary_cols, drop = FALSE])
  if (nrow(pair_dt) > 0) {
    pair_dt[, preferred_orientation := marker_1 == ..marker_1 & marker_2 == ..marker_2]
    data.table::setorderv(
      pair_dt,
      c("sample_alias", "celltype_manual", "preferred_orientation"),
      c(1L, 1L, -1L)
    )
    pair_dt <- unique(pair_dt, by = c("sample_alias", "celltype_manual"))
    pair_dt[, preferred_orientation := NULL]
  }

  detail <- merge(
    totals,
    pair_dt[, .(sample_alias, celltype_manual, sum_log2_ratio, n_detected)],
    by = c("sample_alias", "celltype_manual"),
    all.x = TRUE,
    sort = FALSE
  )
  detail[is.na(sum_log2_ratio), sum_log2_ratio := 0]
  detail[is.na(n_detected), n_detected := 0L]
  detail[, n_detected := as.integer(n_detected)]
  detail[, pair_observed := n_detected > 0]
  denominator <- if (identical(mean_type, "detected")) detail$n_detected else detail$n_total
  detail[, mean_log2_ratio := ifelse(denominator > 0, sum_log2_ratio / denominator, NA_real_)]
  detail[, pct_detected := ifelse(n_total > 0, n_detected / n_total, NA_real_)]
  data.table::setorderv(detail, c("condition", "sample_alias", "celltype_manual"))

  as.data.frame(detail[, .(
    condition,
    sample_alias,
    celltype_manual,
    mean_log2_ratio,
    pct_detected,
    n_detected,
    n_total,
    pair_observed
  )])
}

colocalization_pair_detail_base_dimensions <- function(summary) {
  group_count <- max(1L, length(unique(as.character(summary$condition))))
  facet_columns <- min(2L, group_count)
  facet_rows <- ceiling(group_count / facet_columns)
  samples_per_group <- table(as.character(summary$condition), as.character(summary$sample_alias)) > 0
  max_samples <- if (length(samples_per_group) == 0) 1L else max(rowSums(samples_per_group))
  population_count <- max(1L, length(unique(as.character(summary$celltype_manual))))

  list(
    width = min(1800, max(760, facet_columns * max(340, max_samples * 82) + 210)),
    height = min(3000, max(480, facet_rows * max(230, population_count * 42) + 150)),
    margin = list(l = 120, r = 180, t = 64, b = 120)
  )
}

plot_colocalization_pair_detail <- function(
  summary,
  marker_1,
  marker_2,
  value_label = "Population mean log2 ratio",
  legend_range = c(-1, 1)
) {
  plot_data <- summary
  plot_data$plot_condition <- factor(plot_data$condition, levels = unique(plot_data$condition))
  plot_data$plot_sample <- factor(plot_data$sample_alias, levels = unique(plot_data$sample_alias))
  plot_data$plot_celltype <- factor(plot_data$celltype_manual, levels = rev(unique(plot_data$celltype_manual)))
  plot_data$plot_value <- plot_data$mean_log2_ratio
  plot_data$plot_value[!plot_data$pair_observed] <- NA_real_
  plot_data$hover <- paste0(
    "Analysis group: ", plot_data$condition,
    "<br>Sample: ", plot_data$sample_alias,
    "<br>Cell population: ", plot_data$celltype_manual,
    "<br>", value_label, ": ", ifelse(plot_data$pair_observed, round(plot_data$plot_value, 3), "Not available"),
    "<br>Pair status: ", ifelse(plot_data$pair_observed, "Detected pair", "No detected pair"),
    "<br>Detected cells: ", plot_data$n_detected, " / ", plot_data$n_total,
    "<br>Detected fraction: ", format_percent(plot_data$pct_detected)
  )

  missing_pairs <- plot_data[!plot_data$pair_observed, , drop = FALSE]
  observed_pairs <- plot_data[plot_data$pair_observed, , drop = FALSE]
  ggplot(plot_data, aes(plot_sample, plot_celltype, text = hover)) +
    geom_point(data = missing_pairs, shape = 4, color = "#8a9493", size = 2.8, stroke = 0.9, show.legend = FALSE) +
    geom_point(data = observed_pairs, aes(fill = plot_value, size = pct_detected), shape = 21, color = "#263238", stroke = 0.25) +
    scale_fill_gradient2(
      low = "#176d73",
      mid = "#f7f8f7",
      high = "#c7503e",
      midpoint = 0,
      limits = legend_range,
      oob = squish_to_limits,
      name = value_label
    ) +
    scale_size_continuous(
      range = c(2.5, 9),
      limits = c(0, 1),
      breaks = c(0.25, 0.5, 0.75, 1),
      labels = qc_percent_axis_labels,
      name = "Detected fraction"
    ) +
    facet_wrap(~plot_condition, ncol = min(2L, length(unique(plot_data$plot_condition))), scales = "free_x") +
    labs(
      title = paste(marker_1, "↔", marker_2, "by sample and cell population"),
      x = "Sample",
      y = "Cell population",
      caption = "× = no detected pair; dot size = detected fraction"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(color = "#5f6b6a", hjust = 0),
      legend.position = "right"
    )
}

format_colocalization_pair_detail_table <- function(summary, value_label, max_rows = 100L) {
  result <- head(summary, max_rows)
  result$mean_log2_ratio <- ifelse(
    result$pair_observed,
    round(result$mean_log2_ratio, 3),
    NA_real_
  )
  result$pct_detected <- format_percent(result$pct_detected)
  result$pair_status <- ifelse(result$pair_observed, "Detected pair", "No detected pair")
  result <- result[, c(
    "condition", "sample_alias", "celltype_manual", "mean_log2_ratio",
    "pct_detected", "n_detected", "n_total", "pair_status"
  ), drop = FALSE]
  names(result) <- c(
    "analysis_group", "sample", "cell_population", value_label,
    "detected_fraction", "detected_cells", "total_cells", "pair_status"
  )
  result
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
