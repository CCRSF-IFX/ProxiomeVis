qc_sidebar <- function(id) {
  ns <- NS(id)

  sidebar(
    title = "QC controls",
    width = 300,
    accordion(
      open = c("Filters", "Cutoffs", "Display"),
      accordion_panel(
        "Filters",
        selectizeInput(ns("qc_sample_filter"), "Sample", choices = character(0), multiple = TRUE),
        selectInput(
          ns("qc_metadata_source"),
          "Metadata",
          choices = c("Original cells" = "origin", "Filtered cells" = "filtered"),
          selected = "origin"
        )
      ),
      accordion_panel(
        "Cutoffs",
        numericInput(ns("qc_n_umi_cutoff"), "n_umi cutoff", value = 10000, min = 0, step = 500),
        numericInput(ns("qc_isotype_cutoff"), "Isotype fraction cutoff", value = 0.001, min = 0, max = 1, step = 0.0005)
      ),
      accordion_panel(
        "Display",
        selectInput(
          ns("qc_filter_y"),
          "Filter count y-axis",
          choices = c("Number of cells" = "count", "Fraction of loaded cells" = "fraction_loaded"),
          selected = "count"
        ),
        checkboxInput(ns("qc_filter_include_total"), "Include TOTAL trajectory", value = TRUE),
        selectInput(ns("qc_metric"), "Distribution metric", choices = character(0))
      )
    )
  )
}

qc_module_ui <- function(id) {
  ns <- NS(id)

  nav_panel(
    "QC",
    layout_sidebar(
      sidebar = qc_sidebar(id),
      navset_card_underline(
        id = ns("qc_mode"),
        title = "QC",
        full_screen = TRUE,
        nav_panel(
          "Filtering",
          uiOutput(ns("qc_metric_row")),
          plot_pane(size = "compact", plotlyOutput(ns("qc_filter_plot"), height = proxiome_plot_height())),
          div(class = "table-pane", tableOutput(ns("qc_filter_table")))
        ),
        nav_panel(
          "Cell Calling",
          plot_pane(size = "wide", plotlyOutput(ns("qc_molecule_rank_plot"), height = proxiome_plot_height()))
        ),
        nav_panel(
          "Distributions",
          plot_pane(size = "compact", plotlyOutput(ns("qc_distribution_plot"), height = proxiome_plot_height()))
        ),
        nav_panel(
          "Metadata",
          div(class = "table-pane", tableOutput(ns("qc_origin_metadata_table")))
        )
      )
    )
  )
}

qc_module_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    observe({
      app_data <- data()
      req(app_data)

      qc_samples <- qc_sample_choices(app_data$qc)
      qc_metric_choices <- available_qc_distribution_choices(app_data$qc$origin_metadata)

      updateSelectizeInput(session, "qc_sample_filter", choices = qc_samples, selected = qc_samples)
      updateSelectInput(session, "qc_metric", choices = qc_metric_choices, selected = qc_metric_choices[1])
    })

    qc_metadata <- reactive({
      app_data <- data()
      req(app_data)

      metadata <- selected_qc_metadata(app_data$qc, input$qc_metadata_source)
      filter_qc_metadata_by_sample(metadata, input$qc_sample_filter)
    })

    qc_origin_metadata <- reactive({
      app_data <- data()
      req(app_data)

      filter_qc_metadata_by_sample(app_data$qc$origin_metadata, input$qc_sample_filter)
    })

    qc_filtered_metadata <- reactive({
      app_data <- data()
      req(app_data)

      filter_qc_metadata_by_sample(app_data$qc$filtered_metadata, input$qc_sample_filter)
    })

    qc_filter_counts <- reactive({
      app_data <- data()
      req(app_data)

      qc_filter_counts_for_samples(
        app_data$qc$filter_counts,
        input$qc_sample_filter,
        include_total = isTRUE(input$qc_filter_include_total)
      )
    })

    output$qc_metric_row <- renderUI({
      origin_metadata <- qc_origin_metadata()
      filtered_metadata <- qc_filtered_metadata()

      loaded <- nrow(origin_metadata)
      final <- nrow(filtered_metadata)
      retained <- if (loaded > 0) final / loaded else NA_real_
      samples <- if ("sample" %in% names(origin_metadata)) length(unique(origin_metadata$sample)) else 1L

      metric_row(
        metric_tile("Loaded Cells", format(loaded, big.mark = ",")),
        metric_tile("Final Cells", format(final, big.mark = ",")),
        metric_tile("Retained", format_percent(retained)),
        metric_tile("Samples", format(samples, big.mark = ","))
      )
    })

    output$qc_filter_plot <- renderPlotly({
      counts <- qc_filter_counts()
      validate(need(nrow(counts) > 0, "No QC filter counts are available."))

      p <- plot_filter_cell_counts(
        counts,
        include_total = isTRUE(input$qc_filter_include_total),
        y = input$qc_filter_y %||% "count"
      )

      ggplotly(p, tooltip = "text") |>
        apply_proxiome_plot_frame()
    })

    output$qc_filter_table <- renderTable({
      format_qc_filter_table(qc_filter_counts())
    }, striped = TRUE, bordered = FALSE, width = "100%")

    output$qc_molecule_rank_plot <- renderPlotly({
      metadata <- qc_origin_metadata()
      validate(need("n_umi" %in% names(metadata), "The original metadata does not include n_umi."))
      validate(need(any(is.finite(metadata$n_umi) & metadata$n_umi > 0), "No positive n_umi values are available."))

      qc_molecule_rank_plotly(
        metadata,
        cutoff = numeric_input_value(input$qc_n_umi_cutoff, 10000)
      )
    })

    output$qc_distribution_plot <- renderPlotly({
      metadata <- qc_metadata()
      req(input$qc_metric)
      validate(need(input$qc_metric %in% names(metadata), "Selected QC metric is not available."))

      qc_distribution_plotly(
        metadata,
        metric = input$qc_metric,
        isotype_cutoff = numeric_input_value(input$qc_isotype_cutoff, 0.001)
      )
    })

    output$qc_origin_metadata_table <- renderTable({
      head(qc_origin_metadata(), 200)
    }, striped = TRUE, bordered = FALSE, width = "100%")
  })
}

qc_sample_choices <- function(qc) {
  metadata <- qc$origin_metadata
  if (!"sample" %in% names(metadata)) {
    return(character(0))
  }
  sort(unique(as.character(metadata$sample)))
}

selected_qc_metadata <- function(qc, source) {
  if (identical(source, "filtered")) {
    return(qc$filtered_metadata)
  }
  qc$origin_metadata
}

filter_qc_metadata_by_sample <- function(metadata, selected_samples) {
  if (!"sample" %in% names(metadata) || is.null(selected_samples) || length(selected_samples) == 0) {
    return(metadata)
  }

  metadata[metadata$sample %in% selected_samples, , drop = FALSE]
}

qc_filter_counts_for_samples <- function(filter_counts, selected_samples, include_total = TRUE) {
  if (nrow(filter_counts) == 0 || !"sample" %in% names(filter_counts)) {
    return(filter_counts)
  }

  available_samples <- setdiff(unique(as.character(filter_counts$sample)), "TOTAL")
  selected_samples <- selected_or_all(selected_samples, available_samples)
  selected_samples <- intersect(selected_samples, available_samples)

  rows <- filter_counts[as.character(filter_counts$sample) %in% selected_samples, , drop = FALSE]
  if (nrow(rows) == 0) {
    return(filter_counts[0, , drop = FALSE])
  }

  if (isTRUE(include_total)) {
    total_rows <- if (length(selected_samples) == length(available_samples)) {
      filter_counts[as.character(filter_counts$sample) == "TOTAL", , drop = FALSE]
    } else {
      build_qc_filter_total_rows(rows, filter_counts)
    }
    rows <- rbind(rows, total_rows)
  }

  order_qc_filter_counts(rows)
}

build_qc_filter_total_rows <- function(rows, template) {
  group_cols <- intersect(c("step", "step_label"), names(rows))
  totals <- stats::aggregate(
    rows$n_cells,
    rows[, group_cols, drop = FALSE],
    sum
  )
  names(totals)[ncol(totals)] <- "n_cells"

  totals$sample <- "TOTAL"
  if ("condition" %in% names(template)) {
    totals$condition <- "TOTAL"
  }

  totals$fraction_loaded <- ifelse(totals$n_cells[1] > 0, totals$n_cells / totals$n_cells[1], NA_real_)
  align_qc_filter_count_columns(totals, template)
}

align_qc_filter_count_columns <- function(rows, template) {
  for (column in setdiff(names(template), names(rows))) {
    rows[[column]] <- NA
  }
  rows[, names(template), drop = FALSE]
}

order_qc_filter_counts <- function(filter_counts) {
  if (nrow(filter_counts) == 0 || !"step" %in% names(filter_counts)) {
    return(filter_counts)
  }

  step_order <- unique(as.character(filter_counts$step))
  filter_counts$step <- factor(filter_counts$step, levels = step_order)

  if ("sample" %in% names(filter_counts)) {
    sample_order <- unique(as.character(filter_counts$sample))
    total_index <- which(sample_order == "TOTAL")
    if (length(total_index) > 0) {
      sample_order <- c(sample_order[-total_index], "TOTAL")
    }
    filter_counts$sample <- factor(filter_counts$sample, levels = sample_order)
    filter_counts <- filter_counts[order(filter_counts$step, filter_counts$sample), , drop = FALSE]
    filter_counts$sample <- as.character(filter_counts$sample)
  } else {
    filter_counts <- filter_counts[order(filter_counts$step), , drop = FALSE]
  }

  filter_counts$step <- as.character(filter_counts$step)
  rownames(filter_counts) <- NULL
  filter_counts
}

plot_filter_cell_counts <- function(
  filter_cell_counts,
  include_total = TRUE,
  y = c("count", "fraction_loaded")
) {
  y <- match.arg(y)
  required_cols <- c("step", "n_cells")
  missing_cols <- setdiff(required_cols, names(filter_cell_counts))
  if (length(missing_cols) > 0) {
    stop("Missing columns for filter cell count plot: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  plot_df <- filter_cell_counts
  if (!"sample" %in% names(plot_df)) {
    plot_df$sample <- "All cells"
  }
  if (!"condition" %in% names(plot_df)) {
    plot_df$condition <- "All cells"
  }
  if (!"step_label" %in% names(plot_df)) {
    plot_df$step_label <- format_qc_step_label(plot_df$step)
  }

  plot_df$step <- factor(plot_df$step, levels = unique(as.character(plot_df$step)))
  plot_df$sample <- as.character(plot_df$sample)
  plot_df$condition <- as.character(plot_df$condition)

  if (!include_total) {
    plot_df <- plot_df[plot_df$sample != "TOTAL", , drop = FALSE]
  }

  group_key <- interaction(plot_df$sample, plot_df$condition, drop = TRUE, lex.order = TRUE)
  loaded_cells <- ave(plot_df$n_cells, group_key, FUN = function(values) values[1])
  plot_df$fraction_loaded_for_hover <- ifelse(loaded_cells > 0, plot_df$n_cells / loaded_cells, NA_real_)

  if (identical(y, "fraction_loaded")) {
    plot_df$value <- plot_df$fraction_loaded_for_hover
    y_label <- "Fraction of loaded cells"
    y_scale <- scale_y_continuous(labels = qc_percent_axis_labels)
  } else {
    plot_df$value <- plot_df$n_cells
    y_label <- "Number of cells"
    y_scale <- scale_y_continuous(labels = qc_count_axis_labels)
  }

  plot_df$hover <- paste0(
    "Step: ", plot_df$step_label,
    "<br>Sample: ", plot_df$sample,
    "<br>Condition: ", plot_df$condition,
    "<br>Cells: ", format(plot_df$n_cells, big.mark = ","),
    "<br>Fraction loaded: ", format_percent(plot_df$fraction_loaded_for_hover)
  )

  ggplot(plot_df, aes(x = step, y = value, group = sample, color = sample, text = hover)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    y_scale +
    labs(
      x = "Filtering step",
      y = y_label,
      color = "Sample"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

qc_percent_axis_labels <- function(values) {
  paste0(round(100 * values), "%")
}

qc_count_axis_labels <- function(values) {
  format(values, big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_qc_filter_table <- function(filter_counts) {
  if (nrow(filter_counts) == 0) {
    return(filter_counts)
  }

  filter_counts$fraction_loaded <- round(filter_counts$fraction_loaded, 3)
  filter_counts
}

rank_qc_metadata <- function(metadata, value_col) {
  if (!"sample" %in% names(metadata)) {
    metadata$sample <- "All cells"
  }

  ranked <- lapply(split(metadata, metadata$sample), function(sample_data) {
    sample_data <- sample_data[order(sample_data[[value_col]], decreasing = TRUE), , drop = FALSE]
    sample_data$rank <- seq_len(nrow(sample_data))
    sample_data
  })

  result <- do.call(rbind, ranked)
  rownames(result) <- NULL
  result
}

qc_molecule_rank_plotly <- function(metadata, cutoff = 10000) {
  if (!"n_umi" %in% names(metadata)) {
    stop("The metadata does not include n_umi.", call. = FALSE)
  }

  ranked <- rank_qc_metadata(metadata, value_col = "n_umi")
  ranked <- ranked[is.finite(ranked$rank) & is.finite(ranked$n_umi) & ranked$rank > 0 & ranked$n_umi > 0, , drop = FALSE]
  if (nrow(ranked) == 0) {
    return(apply_proxiome_plot_frame(plotly::plot_ly()))
  }

  ranked$hover <- paste0(
    "Cell: ", ranked$component,
    "<br>Sample: ", ranked$sample,
    "<br>Rank: ", format(ranked$rank, big.mark = ","),
    "<br>n_umi: ", format(round(ranked$n_umi, 1), big.mark = ",")
  )

  widget <- plotly::plot_ly()
  for (sample_name in unique(ranked$sample)) {
    sample_data <- ranked[ranked$sample == sample_name, , drop = FALSE]
    widget <- plotly::add_trace(
      widget,
      data = sample_data,
      x = ~rank,
      y = ~n_umi,
      text = ~hover,
      name = sample_name,
      type = "scatter",
      mode = "lines+markers",
      hoverinfo = "text",
      line = list(width = 1.1),
      marker = list(size = 4.5, opacity = 0.7),
      inherit = FALSE
    )
  }

  cutoff <- numeric_input_value(cutoff, 10000)
  if (is.finite(cutoff) && cutoff > 0) {
    rank_range <- range(ranked$rank, na.rm = TRUE)
    widget <- plotly::add_trace(
      widget,
      x = rank_range,
      y = rep(cutoff, 2),
      text = paste0("n_umi cutoff: ", format(cutoff, big.mark = ",")),
      name = "n_umi cutoff",
      type = "scatter",
      mode = "lines",
      hoverinfo = "text",
      line = list(color = "#c7503e", dash = "dash", width = 1.2),
      inherit = FALSE
    )
  }

  plotly::layout(
    widget,
    xaxis = list(title = "Cell rank", type = "log"),
    yaxis = list(title = "n_umi", type = "log"),
    legend = list(title = list(text = "Sample"))
  ) |>
    apply_proxiome_plot_frame()
}

qc_distribution_plotly <- function(metadata, metric, isotype_cutoff = 0.001) {
  metric <- as.character(metric)[1]
  plot_data <- qc_distribution_plot_data(metadata, metric)
  metric_label <- qc_metric_label(metric)

  widget <- plotly::plot_ly()
  if (nrow(plot_data) == 0) {
    return(
      plotly::layout(
        widget,
        xaxis = list(title = ""),
        yaxis = list(title = metric_label)
      ) |>
        apply_proxiome_plot_frame()
    )
  }

  for (sample_name in unique(plot_data$sample)) {
    sample_data <- plot_data[plot_data$sample == sample_name, , drop = FALSE]
    widget <- plotly::add_trace(
      widget,
      x = rep(sample_name, nrow(sample_data)),
      y = sample_data$metric_value,
      text = sample_data$hover,
      name = sample_name,
      type = "violin",
      hoverinfo = "text",
      points = "all",
      jitter = 0.28,
      pointpos = 0,
      box = list(visible = TRUE),
      meanline = list(visible = TRUE),
      marker = list(size = 3.5, opacity = 0.35),
      line = list(width = 1),
      spanmode = "hard",
      inherit = FALSE
    )
  }

  yaxis <- list(title = metric_label)
  if (is_qc_log_metric(metric)) {
    yaxis$type <- "log"
  }

  shapes <- list()
  if (identical(metric, "isotype_fraction") && is.finite(isotype_cutoff)) {
    shapes <- list(list(
      type = "line",
      xref = "paper",
      x0 = 0,
      x1 = 1,
      y0 = isotype_cutoff,
      y1 = isotype_cutoff,
      line = list(color = "#c7503e", dash = "dash", width = 1.2)
    ))
  }

  plotly::layout(
    widget,
    xaxis = list(title = "", categoryorder = "array", categoryarray = unique(plot_data$sample)),
    yaxis = yaxis,
    showlegend = FALSE,
    shapes = shapes
  ) |>
    apply_proxiome_plot_frame()
}

qc_distribution_plot_data <- function(metadata, metric) {
  if (!metric %in% names(metadata)) {
    stop("Selected QC metric is not available: ", metric, call. = FALSE)
  }

  plot_data <- metadata
  if (!"component" %in% names(plot_data)) {
    plot_data$component <- seq_len(nrow(plot_data))
  }
  if (!"sample" %in% names(plot_data)) {
    plot_data$sample <- "All cells"
  }

  metric_values <- numeric_metric_vector(plot_data[[metric]], expected_length = nrow(plot_data), metric = metric)
  keep <- is.finite(metric_values)
  if (is_qc_log_metric(metric)) {
    keep <- keep & metric_values > 0
  }

  plot_data <- plot_data[keep, , drop = FALSE]
  plot_data$metric_value <- metric_values[keep]
  plot_data$sample <- as.character(plot_data$sample)

  metric_label <- qc_metric_label(metric)
  plot_data$hover <- paste0(
    "Cell: ", plot_data$component,
    "<br>Sample: ", plot_data$sample,
    "<br>", metric_label, ": ", signif(plot_data$metric_value, 4)
  )

  rownames(plot_data) <- NULL
  plot_data
}

numeric_metric_vector <- function(values, expected_length, metric) {
  if (is.data.frame(values)) {
    if (ncol(values) != 1) {
      stop("QC metric column is not one-dimensional: ", metric, call. = FALSE)
    }
    values <- values[[1]]
  }

  if (is.factor(values)) {
    values <- as.character(values)
  }
  if (is.list(values) && !is.atomic(values)) {
    values <- unlist(values, use.names = FALSE)
  }

  values <- as.vector(values)
  if (length(values) != expected_length) {
    stop("QC metric length does not match metadata rows: ", metric, call. = FALSE)
  }

  suppressWarnings(as.numeric(values))
}

qc_metric_label <- function(metric) {
  choices <- c(
    n_umi = "UMIs",
    n_edges = "Edges",
    reads_in_component = "Reads in component",
    isotype_fraction = "Isotype fraction",
    tau = "Tau"
  )

  label <- choices[metric]
  if (is.na(label)) {
    return(metric)
  }
  unname(label)
}

is_qc_log_metric <- function(metric) {
  metric %in% c("n_umi", "n_edges", "reads_in_component")
}
