hybrid_analysis_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "New comparison", width = 310,
      shiny::selectInput(ns("analysis_scope"), "Analysis scope", choices = c(
        "Within-donor exploration" = "within_donor", "Across donors" = "across_donors")),
      shiny::selectInput(ns("comparison_axis"), "Compare", choices = c(
        "Conditions" = "condition", "Cell types" = "celltype")),
      shiny::selectInput(ns("readout"), "Readout", choices = PRECOMPUTED_READOUTS),
      shiny::conditionalPanel("input.analysis_scope !== 'across_donors'", ns = ns,
        shiny::selectInput(ns("donor"), "One donor", choices = character())),
      shiny::conditionalPanel("input.analysis_scope === 'across_donors'", ns = ns,
        shiny::selectizeInput(ns("donors"), "Donors", choices = character(), multiple = TRUE),
        shiny::selectInput(ns("design"), "Donor design", choices = c(
          "Paired — same donors in both groups" = "paired",
          "Independent — different donors in each group" = "independent"))),
      shiny::conditionalPanel("input.comparison_axis !== 'celltype'", ns = ns,
        shiny::selectInput(ns("celltype"), "One cell type", choices = character())),
      shiny::conditionalPanel("input.comparison_axis === 'celltype'", ns = ns,
        shiny::selectInput(ns("condition"), "One condition", choices = character())),
      shiny::selectInput(ns("target"), "Target condition", choices = character()),
      shiny::selectInput(ns("reference"), "Reference condition", choices = character()),
      shiny::conditionalPanel("input.readout !== 'abundance'", ns = ns,
        shiny::numericInput(ns("min_cells_per_group"), "Minimum observed cells per donor and group",
          value = 20L, min = 1L, max = 1000L, step = 1L),
        shiny::tags$p(class = "small text-muted",
          "Isotype-background filtering is enabled using the prepared donor cutoff. The minimum applies to each marker or pair after filtering (min_cells_per_group).")
      ),
      shiny::tags$p(class = "small text-muted", "All dataset markers are tested. Spatial tests use available cell-level scores."),
      shiny::tags$p(class = "small text-muted", "All cells in the selected scope are included. Technical samples do not count as independent donors."),
      shiny::textOutput(ns("selection_summary")),
      shiny::conditionalPanel("input.analysis_scope !== 'across_donors'", ns = ns,
        shiny::checkboxInput(ns("acknowledge"), "I understand this is exploratory within-donor analysis.", FALSE)),
      bslib::input_task_button(ns("run"), "Run comparison"),
      shiny::actionButton(ns("cancel"), "Cancel running analysis", class = "btn-sm btn-outline-secondary"),
      shiny::tags$p(class = "small text-muted",
        "Limits: 2,000 cells per donor, 40 donors, 80,000 cells in total, 160 markers, 120 seconds. Spatial tests run in batches of at most 1,000 pairs and 1 million score rows. Abundance requires 3 cells per donor and group. Spatial minimum defaults: self-clustering 20, colocalization 40."),
      shiny::uiOutput(ns("inference_notice"))
    ),
    shiny::uiOutput(ns("availability")),
    shiny::uiOutput(ns("readiness")),
    shiny::tags$div(role = "status", `aria-live` = "polite", shiny::textOutput(ns("status"))),
    bslib::card(full_screen = TRUE, fill = FALSE,
      bslib::card_header("Completed comparisons in this session"),
      bslib::card_body(fill = FALSE,
        shiny::tags$p("Choose a completed comparison, then open it in Results. Saved dataset results stay unchanged."),
        shiny::selectInput(ns("history_id"), "Completed comparison", choices = character()),
        shiny::uiOutput(ns("open_control")),
        shiny::textOutput(ns("result_scope")),
        precomputed_download_group(
          shiny::downloadButton(ns("results_csv"), "Results CSV"),
          shiny::downloadButton(ns("provenance_json"), "Analysis provenance JSON"),
          shiny::downloadButton(ns("result_rds"), "Comparison RDS")
        ),
        shiny::tags$p(class = "small text-muted", "History is cleared when you load another dataset or close this session. Download a comparison to keep it. The RDS contains results and settings, not the cell-level dataset."))
    )
  )
}

# Use the same canonical tables as saved analyses; never reinterpret a profile mean as a median.
hybrid_result_contract <- function(data, entry) {
  request <- entry$request
  result <- entry$result
  if (!is.list(request) || !is.list(result) || !is.data.frame(result$results) ||
      !identical(request$dataset_id, data$contract$dataset_id) ||
      !identical(request$bundle_md5, data$bundle_md5)) {
    stop("This comparison does not belong to the loaded application bundle.", call. = FALSE)
  }
  for (field in c("id", "label")) {
    if (!is.character(entry[[field]]) || length(entry[[field]]) != 1L ||
        is.na(entry[[field]]) || !nzchar(entry[[field]])) stop("Comparison ", field, " must be a nonempty string.", call. = FALSE)
  }
  if (!identical(result$provenance$dataset_id, request$dataset_id) ||
      !identical(result$provenance$bundle_md5, request$bundle_md5) ||
      !identical(result$provenance$request, request)) {
    stop("Comparison provenance does not match its submitted request.", call. = FALSE)
  }
  request <- validate_hybrid_request(data, request)
  readout <- request$readout
  rows <- as.data.frame(result$results)
  required <- c("feature_id", "marker_1", "marker_2", "effect_size", "p_value", "p_adj", "n_cells_a", "n_cells_b")
  if (!all(required %in% names(rows))) stop("Comparison results are missing required feature or statistic columns.", call. = FALSE)
  expected <- hybrid_feature_table(request)
  keys <- paste(rows$marker_1, rows$marker_2, sep = "\r")
  index <- match(keys, paste(expected$marker_1, expected$marker_2, sep = "\r"))
  if (nrow(rows) != nrow(expected) || anyDuplicated(keys) || anyNA(index) ||
      !identical(as.character(rows$feature_id), expected$feature_id[index])) {
    stop("Comparison features do not match its submitted marker set.", call. = FALSE)
  }
  abundance <- data$contract$sample_summaries$abundance
  controls <- unique(c(data$contract$control_markers, abundance$marker_1[abundance$is_control %in% TRUE]))
  across <- identical(request$analysis_scope, "across_donors")
  decorate <- function(table, values = FALSE, aggregate = FALSE) {
    table <- as.data.frame(table)
    table$analysis_id <- rep(entry$id, nrow(table))
    table$celltype <- rep(request$celltype, nrow(table))
    table$group_a <- rep(request$target, nrow(table))
    table$group_b <- rep(request$reference, nrow(table))
    table$comparison <- rep(paste(request$target, "vs", request$reference), nrow(table))
    table$comparison_axis <- rep(if (is.null(request$comparison_axis)) "condition" else request$comparison_axis, nrow(table))
    table$design <- rep(if (is.null(request$design)) "within_donor" else request$design, nrow(table))
    table$fixed_condition <- rep(if (identical(request$comparison_axis, "celltype")) request$condition else NA_character_, nrow(table))
    table$metric <- rep(if (aggregate) {
      if (readout == "abundance") "delta_mean_donor_mean_normalized_abundance" else "delta_mean_donor_median_log2_ratio"
    } else if (readout == "abundance") {
      if (values) "mean_normalized_abundance" else "delta_mean_normalized_abundance"
    } else if (values) "median_log2_ratio" else "delta_median_log2_ratio", nrow(table))
    table$is_control <- table$marker_1 %in% controls | table$marker_2 %in% controls
    table
  }
  rows <- decorate(rows, aggregate = across)
  if (across) rows$adjustment_family <- paste(entry$id, readout, request$celltype, sep = ":")
  contract <- data$contract
  for (section in c("donor_contrasts", "differential_results")) {
    contract[[section]] <- lapply(contract[[section]], function(table) table[FALSE, , drop = FALSE])
  }
  if (across) {
    if (!is.data.frame(result$donor_values) || !is.data.frame(result$donor_contrasts)) {
      stop("Across-donor results must include donor measurements and paired contrasts.", call. = FALSE)
    }
    value_columns <- c("feature_id", "marker_1", "marker_2", "donor", "group", "value",
      "n_total_cells", "n_evaluable_cells", "n_scored_cells", "eligible", "source_celltype", "condition")
    contrast_columns <- c("feature_id", "marker_1", "marker_2", "donor", "value_a", "value_b",
      "effect_size", "n_cells_a", "n_cells_b")
    if (!all(value_columns %in% names(result$donor_values)) ||
        !all(contrast_columns %in% names(result$donor_contrasts))) {
      stop("Donor measurements are missing required feature, group, or support columns.", call. = FALSE)
    }
    values <- decorate(result$donor_values, values = TRUE)
    contrasts <- decorate(result$donor_contrasts)
    for (table in list(values, contrasts)) {
      matched <- match(paste(table$marker_1, table$marker_2, sep = "\r"), keys)
      if (anyNA(matched) || !identical(as.character(table$feature_id), as.character(rows$feature_id[matched])) ||
          !all(table$donor %in% union(request$donors_a, request$donors_b))) {
        stop("Donor measurements do not match the submitted scope.", call. = FALSE)
      }
    }
    if (!all(values$group %in% c(request$target, request$reference)) ||
        anyDuplicated(values[c("feature_id", "donor", "group")]) ||
        anyDuplicated(contrasts[c("feature_id", "donor")]) ||
        (identical(request$design, "independent") && nrow(contrasts))) {
      stop("Donor measurements contain invalid groups or duplicate records.", call. = FALSE)
    }
    retained <- union(request$donors_a, request$donors_b)
    if (nrow(values) != nrow(rows) * length(retained) * 2L) {
      stop("Donor measurements must cover every requested feature, retained donor, and group.", call. = FALSE)
    }
    total_cells <- integer(nrow(values))
    for (group in c("a", "b")) {
      selected <- values$group == if (group == "a") request$target else request$reference
      counts <- table(data$metadata$donor[match(request[[paste0("group_", group, "_cells")]], data$metadata$component)])
      support <- as.integer(counts[values$donor[selected]])
      support[is.na(support)] <- 0L
      total_cells[selected] <- support
    }
    support <- values[c("n_total_cells", "n_evaluable_cells", "n_scored_cells")]
    if (!all(vapply(support, is.numeric, logical(1))) || anyNA(support) ||
        any(as.matrix(support) < 0 | as.matrix(support) != floor(as.matrix(support))) ||
        any(values$n_total_cells != total_cells | values$n_evaluable_cells != total_cells |
          values$n_scored_cells > total_cells) || !is.logical(values$eligible) || anyNA(values$eligible) ||
        any(values$eligible != (values$n_scored_cells >= request$min_cells)) ||
        !is.numeric(values$value) || any(values$eligible & !is.finite(values$value)) ||
        any(!values$eligible & !is.na(values$value))) {
      stop("Donor measurements have inconsistent cell support, eligibility, or values.", call. = FALSE)
    }
    aggregate_columns <- c("n_donors_a", "n_donors_b", "n_pairs", "n_cells_a", "n_cells_b",
      "effect_size", "min_effect", "max_effect", "direction_consistent")
    if (!all(c(aggregate_columns, "inference", "status") %in% names(rows))) {
      stop("Across-donor results are missing donor support or effect columns.", call. = FALSE)
    }
    if (anyNA(rows[c("inference", "status")]) || !all(rows$inference %in% c("donor_level", "descriptive")) ||
        !all(rows$status %in% c("tested", "insufficient_donors", "no_finite_test_result"))) {
      stop("Across-donor results have invalid inference or test status.", call. = FALSE)
    }
    same <- function(actual, expected) isTRUE(all.equal(actual, expected, check.attributes = FALSE, tolerance = 1e-8))
    for (i in seq_len(nrow(rows))) {
      feature <- values[values$feature_id == rows$feature_id[[i]], , drop = FALSE]
      a <- feature[feature$group == request$target & feature$eligible, , drop = FALSE]
      b <- feature[feature$group == request$reference & feature$eligible, , drop = FALSE]
      expected <- list(n_donors_a = nrow(a), n_donors_b = nrow(b), n_pairs = 0L,
        n_cells_a = sum(a$n_total_cells), n_cells_b = sum(b$n_total_cells),
        effect_size = NA_real_, min_effect = NA_real_, max_effect = NA_real_, direction_consistent = NA)
      if (request$design == "paired") {
        pairs <- intersect(a$donor, b$donor)
        a <- a[match(pairs, a$donor), , drop = FALSE]
        b <- b[match(pairs, b$donor), , drop = FALSE]
        contrast <- contrasts[contrasts$feature_id == rows$feature_id[[i]], , drop = FALSE]
        if (!setequal(contrast$donor, pairs)) stop("Paired contrasts must cover every complete eligible donor pair.", call. = FALSE)
        contrast <- contrast[match(pairs, contrast$donor), , drop = FALSE]
        differences <- a$value - b$value
        required_contrast <- list(value_a = a$value, value_b = b$value, effect_size = differences,
          n_cells_a = a$n_total_cells, n_cells_b = b$n_total_cells)
        if (!all(vapply(names(required_contrast), function(field) same(contrast[[field]], required_contrast[[field]]), logical(1)))) {
          stop("Paired contrasts do not match their donor measurements or cell support.", call. = FALSE)
        }
        expected$n_pairs <- length(pairs)
        if (length(pairs)) {
          expected$effect_size <- mean(differences)
          expected$min_effect <- min(differences)
          expected$max_effect <- max(differences)
          expected$direction_consistent <- length(unique(sign(differences))) == 1L
        }
      } else if (nrow(a) && nrow(b)) expected$effect_size <- mean(a$value) - mean(b$value)
      if (!all(vapply(names(expected), function(field) same(rows[[field]][[i]], expected[[field]]), logical(1)))) {
        stop("Across-donor effects or support do not match their donor measurements.", call. = FALSE)
      }
      enough <- if (request$design == "paired") expected$n_pairs >= 2L else min(expected$n_donors_a, expected$n_donors_b) >= 2L
      tested <- rows$status[[i]] == "tested"
      untested_status <- if (enough) "no_finite_test_result" else "insufficient_donors"
      if ((tested && (!enough || rows$inference[[i]] != "donor_level" || !is.finite(rows$effect_size[[i]]) ||
          !is.finite(rows$p_value[[i]]) || !is.finite(rows$p_adj[[i]]))) ||
          (!tested && (rows$inference[[i]] != "descriptive" || !is.na(rows$p_value[[i]]) || !is.na(rows$p_adj[[i]]) ||
            rows$status[[i]] != untested_status))) {
        stop("Across-donor inference and test status do not match the eligible donor support.", call. = FALSE)
      }
    }
    contract$schema_version <- 2L
    contract$donor_values <- stats::setNames(lapply(PROXIOMEVIS_READOUTS,
      function(name) values[FALSE, , drop = FALSE]), PROXIOMEVIS_READOUTS)
    contract$donor_values[[readout]] <- values
    contract$donor_contrasts[[readout]] <- contrasts
    contract$differential_results[[readout]] <- rows
  } else {
    rows$donor <- request$donor
    if (!"value_a" %in% names(rows)) rows$value_a <- NA_real_
    if (!"value_b" %in% names(rows)) rows$value_b <- NA_real_
    contract$donor_contrasts[[readout]] <- rows
    if (identical(request$comparison_axis, "celltype")) contract$schema_version <- 2L
    if (is.list(contract$donor_values)) {
      contract$donor_values <- lapply(contract$donor_values, function(table) table[FALSE, , drop = FALSE])
    }
  }
  contract$session_result <- list(id = entry$id, label = entry$label, request = request, provenance = result$provenance)
  if (is.list(contract$report)) {
    contract$report$source_results <- NULL
    contract$report$eligibility <- NULL
    contract$report$abundance_heatmaps <- NULL
    if (is.list(contract$report$alignment)) {
      contract$report$alignment$differences <- NULL
      contract$report$alignment$contrasts <- NULL
    }
  }
  if (across) {
    checked <- data
    checked$contract <- contract
    mappings <- list(sample_column = "sample_id", donor_column = "donor", condition_column = "condition",
      celltype_column = "celltype", cluster_column = "cluster")
    checked$contract$defaults[names(mappings)] <- mappings
    validate_precomputed_data(checked)
  }
  contract
}

hybrid_analysis_server <- function(id, data, manager, code_dir, on_open = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    job <- NULL
    message <- shiny::reactiveVal("Choose a scope, then run a comparison across all dataset markers.")
    validation_error <- shiny::reactiveVal(NULL)
    last_request <- shiny::reactiveVal(NULL)
    history <- shiny::reactiveVal(list())
    selected_id <- shiny::reactiveVal(NULL)
    open_result <- shiny::reactiveVal(NULL)
    generation <- 0L
    run_number <- 0L
    active_run <- NULL
    # Cancellation settles asynchronously, after the module's reactive state may be destroyed.
    task <- shiny::withReactiveDomain(NULL, shiny::ExtendedTask$new(function(path, request, token, dataset_generation) {
      job <<- start_hybrid_job(path, request, code_dir, manager)
      promises::then(job$promise, function(value) list(result = value, request = request,
        token = token, generation = dataset_generation))
    })) |> bslib::bind_task_button("run")
    session$onSessionEnded(function() if (!is.null(job)) job$cancel())

    shiny::observeEvent(input$readout, {
      shiny::updateNumericInput(session, "min_cells_per_group",
        value = if (identical(input$readout, "colocalization")) 40L else 20L)
    })

    analysis_scope <- shiny::reactive(if (is.null(input$analysis_scope)) "within_donor" else input$analysis_scope)
    comparison_axis <- shiny::reactive(if (is.null(input$comparison_axis)) "condition" else input$comparison_axis)
    selected_donors <- shiny::reactive({
      if (identical(analysis_scope(), "across_donors")) {
        if (is.null(input$donors)) sort(unique(data()$metadata$donor)) else input$donors
      } else input$donor
    })

    shiny::observeEvent(data(), {
      value <- data()
      if (!is.null(job)) job$cancel()
      generation <<- generation + 1L
      active_run <<- NULL
      last_request(NULL)
      history(list())
      selected_id(NULL)
      open_result(NULL)
      validation_error(NULL)
      message("Choose a scope, then run a comparison across all dataset markers.")
      shiny::updateSelectInput(session, "history_id", choices = character())
      donors <- sort(unique(value$metadata$donor))
      shiny::updateSelectInput(session, "donor", choices = donors)
      shiny::updateSelectizeInput(session, "donors", choices = donors, selected = donors)
    }, priority = 100)
    shiny::observeEvent(list(data(), selected_donors()), {
      meta <- data()$metadata
      meta <- meta[meta$donor %in% selected_donors(), , drop = FALSE]
      choices <- sort(unique(meta$celltype))
      shiny::updateSelectInput(session, "celltype", choices = choices,
        selected = if (length(choices)) precomputed_keep_selection(shiny::isolate(input$celltype), choices, "CD4 T") else character())
      conditions <- sort(unique(meta$condition))
      shiny::updateSelectInput(session, "condition", choices = conditions,
        selected = if (length(conditions)) precomputed_keep_selection(shiny::isolate(input$condition), conditions, "unstim") else character())
    })
    scope <- shiny::reactive({
      if (analysis_scope() == "within_donor") shiny::req(length(selected_donors()))
      meta <- data()$metadata
      fixed <- if (comparison_axis() == "celltype") {
        shiny::req(input$condition)
        meta$condition == input$condition
      } else {
        shiny::req(input$celltype)
        meta$celltype == input$celltype
      }
      meta[meta$donor %in% selected_donors() & fixed, , drop = FALSE]
    })
    shiny::observeEvent(list(scope(), comparison_axis()), {
      groups <- sort(unique(scope()[[comparison_axis()]]))
      default_reference <- if ("unstim" %in% groups) "unstim" else head(groups, 1L)
      reference <- if (length(groups)) precomputed_keep_selection(shiny::isolate(input$reference), groups, default_reference) else character()
      targets <- setdiff(groups, reference)
      target <- if (length(targets)) precomputed_keep_selection(shiny::isolate(input$target), targets, head(targets, 1L)) else character()
      group_name <- if (comparison_axis() == "celltype") "cell type" else "condition"
      shiny::updateSelectInput(session, "reference", label = paste("Reference", group_name), choices = groups, selected = reference)
      shiny::updateSelectInput(session, "target", label = paste("Target", group_name), choices = groups, selected = target)
    })
    output$inference_notice <- shiny::renderUI({
      shiny::tags$div(class = "alert alert-warning small", if (analysis_scope() == "across_donors") {
        paste("Donors are the biological replicates, with equal weight. Each donor contributes a mean normalized abundance or a median observed proximity score per group.",
          "Paired designs test complete donor differences with a paired t-test; independent designs use a Welch t-test on different donors.",
          "These tests assume independent donors and approximately normal donor differences (paired) or donor summaries (independent). Few donors limit assessment of these assumptions. Limited support or zero variance can give descriptive results without P values.")
      } else {
        "Cells from one donor are not independent biological replicates. These exploratory cell-level results do not establish population-level significance."
      })
    })
    output$availability <- shiny::renderUI({
      value <- data()
      if (is.null(value$bundle_path)) return(shiny::tags$p(class = "alert alert-info", "Save a prepared application bundle to enable background analysis."))
      if (!identical(input$readout, "abundance") && is.null(value$assets$proximity)) {
        shiny::tags$p(class = "alert alert-info", "Spatial reanalysis needs the prepared proximity sidecar. Saved results remain available without it.")
      }
    })
    prepared_request <- shiny::reactive({
      value <- data()
      meta <- scope()
      tryCatch({
        if (analysis_scope() == "across_donors" && !length(selected_donors())) {
          stop("Choose at least one donor.", call. = FALSE)
        }
        selected <- list(
          dataset_id = value$contract$dataset_id, bundle_md5 = value$bundle_md5, readout = input$readout,
          donor = input$donor,
          celltype = if (comparison_axis() == "celltype") input$target else input$celltype,
          target = input$target, reference = input$reference,
          min_cells_per_group = input$min_cells_per_group,
          markers = rownames(value$abundance))
        # Omitted scenario controls retain the original serialized request for older clients.
        if (!is.null(input$analysis_scope) || !is.null(input$comparison_axis)) {
          selected$analysis_scope <- analysis_scope()
          selected$comparison_axis <- comparison_axis()
          selected$design <- if (analysis_scope() == "across_donors") input$design else "within_donor"
          if (is.null(selected$design)) selected$design <- "paired"
          if (comparison_axis() == "celltype") selected$condition <- input$condition
          if (analysis_scope() == "across_donors") selected$donors <- selected_donors()
        }
        if (analysis_scope() == "within_donor") {
          selected$group_a_cells <- meta$component[meta[[comparison_axis()]] == input$target]
          selected$group_b_cells <- meta$component[meta[[comparison_axis()]] == input$reference]
        }
        request <- validate_hybrid_request(value, selected)
        if (is.null(value$bundle_path)) stop("Save a prepared application bundle to run background analysis.", call. = FALSE)
        if (request$readout == "abundance" && analysis_scope() == "within_donor" && is.null(value$counts)) {
          stop("Raw counts are not stored. Add counts to the bundle to run abundance comparisons.", call. = FALSE)
        }
        if (request$readout != "abundance" && is.null(value$assets$proximity)) {
          stop("The proximity sidecar is not configured. Add it to run spatial comparisons.", call. = FALSE)
        }
        list(request = request, error = NULL)
      }, error = function(error) list(request = NULL, error = conditionMessage(error)))
    })
    output$readiness <- shiny::renderUI({
      ready <- prepared_request()
      value <- data()
      shiny::tags$div(class = paste("alert", if (is.null(ready$error)) "alert-success" else "alert-warning"),
        shiny::tags$strong(if (is.null(ready$error)) "Scope is ready" else "Scope is not ready"),
        shiny::tags$p(selection_summary()),
        shiny::tags$p(paste(nrow(value$abundance), "markers.")),
        shiny::tags$p(paste("Raw counts:", if (is.null(value$counts)) "not stored." else "stored.",
          if (analysis_scope() == "across_donors") "Across-donor abundance uses stored normalized values and does not require counts.",
          "Proximity sidecar:", if (is.null(value$assets$proximity)) "not configured." else "configured (checked by the worker).")),
        if (!is.null(ready$error)) shiny::tags$p(ready$error),
        shiny::tags$p(class = "small mb-0", "Limits are checked before a worker starts. Score coverage and sidecar integrity are checked in the worker. For a larger scope, prepare results offline and load them in Results."))
    })
    selection_summary <- shiny::reactive({
      meta <- scope()
      cells <- paste(sum(meta[[comparison_axis()]] == input$target), "target cells;",
        sum(meta[[comparison_axis()]] == input$reference), "reference cells.")
      if (analysis_scope() != "across_donors") return(cells)
      request <- prepared_request()$request
      if (is.null(request)) return(paste(cells, "Choose valid donor groups to check eligibility."))
      cells <- paste(length(request$group_a_cells), "target cells;", length(request$group_b_cells), "reference cells retained.")
      donors <- if (identical(request$design, "paired")) {
        paste(length(intersect(request$donors_a, request$donors_b)), "complete donor pairs.")
      } else paste(length(request$donors_a), "target donors;", length(request$donors_b), "reference donors.")
      excluded <- request$excluded_donors
      exclusions <- if (is.data.frame(excluded) && nrow(excluded)) {
        paste("Excluded:", paste(paste0(excluded$donor, " (", gsub("_", " ", excluded$reason), ")"), collapse = ", "))
      } else "No donors excluded at scope selection."
      paste(cells, donors, exclusions, "Feature-specific score coverage is checked in the worker.")
    })
    output$selection_summary <- shiny::renderText(selection_summary())
    shiny::observeEvent(input$run, {
      if (identical(task$status(), "running")) return()
      validation_error(NULL)
      tryCatch({
        if (analysis_scope() == "within_donor" && !isTRUE(input$acknowledge)) {
          stop("Confirm the exploratory-analysis notice before running.", call. = FALSE)
        }
        value <- data()
        if (is.null(value$bundle_path)) stop("No saved application bundle is available.", call. = FALSE)
        ready <- prepared_request()
        if (!is.null(ready$error)) stop(ready$error, call. = FALSE)
        request <- ready$request
        run_number <<- run_number + 1L
        active_run <<- paste0("session-", generation, "-", run_number)
        last_request(request)
        message("Running in a separate R process. You can continue exploring the reports.")
        task$invoke(value$bundle_path, request, active_run, generation)
      }, error = function(error) validation_error(conditionMessage(error)))
    })
    shiny::observeEvent(input$cancel, {
      if (identical(task$status(), "running") && !is.null(job)) {
        active_run <<- NULL
        job$cancel()
        message("Analysis cancelled. Previous completed comparisons remain available.")
      }
    })
    shiny::observeEvent(task$status(), {
      if (!identical(task$status(), "success")) return()
      completed <- task$result()
      if (!identical(completed$generation, generation) || !identical(completed$token, active_run)) return()
      request <- completed$request
      if (!identical(request$dataset_id, data()$contract$dataset_id) ||
          !identical(request$bundle_md5, data()$bundle_md5)) return()
      if (!identical(completed$result$provenance$dataset_id, request$dataset_id) ||
          !identical(completed$result$provenance$bundle_md5, request$bundle_md5)) {
        active_run <<- NULL
        validation_error("The worker returned a result for a different bundle. Reload the dataset and run again.")
        return()
      }
      id <- completed$token
      donor_label <- if (identical(request$analysis_scope, "across_donors")) {
        paste(length(unique(c(request$donors_a, request$donors_b))), "donors", request$design)
      } else paste(request$donor, "exploratory")
      label <- paste(paste0("#", run_number), donor_label,
        if (identical(request$comparison_axis, "celltype")) request$condition else request$celltype, request$readout,
        paste(request$target, "vs", request$reference), sep = " · ")
      entries <- history()
      entries[[id]] <- list(id = id, label = label, request = request, result = completed$result)
      history(entries)
      selected_id(id)
      active_run <<- NULL
      shiny::updateSelectInput(session, "history_id", choices = stats::setNames(names(entries),
        vapply(entries, `[[`, character(1), "label")), selected = id)
      message("Completed. Open the comparison in Results. Changing controls does not rerun an analysis.")
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$history_id, {
      if (input$history_id %in% names(history())) selected_id(input$history_id)
    })
    selected_entry <- shiny::reactive({
      shiny::req(selected_id(), selected_id() %in% names(history()))
      history()[[selected_id()]]
    })
    output$open_control <- shiny::renderUI({
      shiny::actionButton(session$ns("open_result"), "Open in Results", class = "btn-primary",
        disabled = if (!length(history())) "disabled" else NULL)
    })
    shiny::observeEvent(input$open_result, {
      entry <- selected_entry()
      open_result(entry)
      if (is.function(on_open)) on_open(entry$id)
    })
    output$status <- shiny::renderText({
      if (!is.null(validation_error())) return(validation_error())
      status <- task$status()
      if (identical(status, "error") && !is.null(last_request())) {
        return(tryCatch(task$result(), error = function(e) conditionMessage(e)))
      }
      message()
    })
    result <- shiny::reactive({
      selected_entry()$result
    })
    output$result_scope <- shiny::renderText({
      result()
      request <- selected_entry()$request
      across <- identical(request$analysis_scope, "across_donors")
      cutoff <- result()$provenance$filters$donor_isotype_fraction_cutoff
      paste(c(if (across) paste("Across donors:", paste(request$donors, collapse = ", "), request$design) else {
          paste("Within-donor exploration:", request$donor)
        }, if (identical(request$comparison_axis, "celltype")) paste("Condition:", request$condition) else request$celltype, request$readout,
        paste(request$target, "minus", request$reference),
        paste(length(request$group_a_cells), "vs", length(request$group_b_cells), "cells in scope"),
        if (request$readout != "abundance") paste0(if (across || identical(request$comparison_axis, "celltype")) {
            "Minimum observed cells per donor and group: "
          } else "Minimum observed cells per condition: ", request$min_cells,
          if (is.numeric(cutoff) && length(cutoff) == 1L) paste0("; isotype-background cutoff: ", signif(cutoff, 4))),
        "BH adjustment applies to this requested feature family; frozen report results are unchanged."), collapse = " | ")
    })
    output$results_csv <- shiny::downloadHandler(
      filename = function() "proxiomevis_comparison.csv",
      content = function(file) utils::write.csv(result()$results, file, row.names = FALSE)
    )
    output$provenance_json <- shiny::downloadHandler(
      filename = function() "proxiomevis_analysis_provenance.json",
      content = function(file) jsonlite::write_json(list(request = selected_entry()$request, provenance = result()$provenance),
        file, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
    )
    output$result_rds <- shiny::downloadHandler(
      filename = function() "proxiomevis_comparison.rds",
      content = function(file) saveRDS(list(result_version = 1L, entry = selected_entry()), file, compress = "xz")
    )
    invisible(list(task = task, result = result, history = history, selected_id = selected_id,
      selected_entry = selected_entry, open_result = open_result))
  })
}
