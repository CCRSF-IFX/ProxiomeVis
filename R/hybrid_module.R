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
      shiny::conditionalPanel("input.readout !== 'abundance'", ns = ns,
        shiny::numericInput(ns("min_cells_per_group"), "Minimum observed cells per condition",
          value = 20L, min = 1L, max = 1000L, step = 1L),
        shiny::tags$p(class = "small text-muted",
          "Isotype-background filtering is enabled using the prepared donor cutoff. The minimum applies to each marker or pair after filtering (min_cells_per_group).")
      ),
      shiny::tags$p(class = "small text-muted", "All dataset markers are tested. Spatial tests use available cell-level scores."),
      shiny::tags$p(class = "small text-muted", "All cells in the selected donor, cell type, and conditions are included."),
      shiny::textOutput(ns("selection_summary")),
      shiny::checkboxInput(ns("acknowledge"), "I understand this is exploratory within-donor analysis.", FALSE),
      bslib::input_task_button(ns("run"), "Run comparison"),
      shiny::actionButton(ns("cancel"), "Cancel running analysis", class = "btn-sm btn-outline-secondary"),
      shiny::tags$p(class = "small text-muted",
        "Limits: 2,000 cells, 160 markers, 120 seconds. Spatial tests run in batches of at most 1,000 pairs and 1 million score rows. Abundance requires 3 cells per condition. Spatial minimum defaults: self-clustering 20, colocalization 40."),
      shiny::tags$div(class = "alert alert-warning small",
        "Cells from one donor are not independent biological replicates. These results do not establish population-level significance.")
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
  rows$analysis_id <- entry$id
  rows$donor <- request$donor
  rows$celltype <- request$celltype
  rows$group_a <- request$target
  rows$group_b <- request$reference
  rows$comparison <- paste(request$target, "vs", request$reference)
  rows$metric <- if (readout == "abundance") "delta_mean_normalized_abundance" else "delta_median_log2_ratio"
  rows$value_a <- rows$value_b <- NA_real_
  abundance <- data$contract$sample_summaries$abundance
  controls <- unique(c(data$contract$control_markers, abundance$marker_1[abundance$is_control %in% TRUE]))
  rows$is_control <- rows$marker_1 %in% controls | rows$marker_2 %in% controls
  contract <- data$contract
  for (section in c("donor_contrasts", "differential_results")) {
    contract[[section]] <- lapply(contract[[section]], function(table) table[FALSE, , drop = FALSE])
  }
  contract$donor_contrasts[[readout]] <- rows
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
      shiny::updateSelectInput(session, "donor", choices = sort(unique(value$metadata$donor)))
    }, priority = 100)
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
        request <- validate_hybrid_request(value, list(
          dataset_id = value$contract$dataset_id, bundle_md5 = value$bundle_md5, readout = input$readout,
          donor = input$donor, celltype = input$celltype, target = input$target, reference = input$reference,
          min_cells_per_group = input$min_cells_per_group,
          markers = rownames(value$abundance), group_a_cells = meta$component[meta$condition == input$target],
          group_b_cells = meta$component[meta$condition == input$reference]
        ))
        if (is.null(value$bundle_path)) stop("Save a prepared application bundle to run background analysis.", call. = FALSE)
        if (request$readout == "abundance" && is.null(value$counts)) {
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
      meta <- scope()
      shiny::tags$div(class = paste("alert", if (is.null(ready$error)) "alert-success" else "alert-warning"),
        shiny::tags$strong(if (is.null(ready$error)) "Scope is ready" else "Scope is not ready"),
        shiny::tags$p(paste(sum(meta$condition == input$target), "target cells;",
          sum(meta$condition == input$reference), "reference cells;", nrow(value$abundance), "markers.")),
        shiny::tags$p(paste("Raw counts:", if (is.null(value$counts)) "not stored." else "stored.",
          "Proximity sidecar:", if (is.null(value$assets$proximity)) "not configured." else "configured (checked by the worker).")),
        if (!is.null(ready$error)) shiny::tags$p(ready$error),
        shiny::tags$p(class = "small mb-0", "Limits are checked before a worker starts. Score coverage and sidecar integrity are checked in the worker. For a larger scope, prepare results offline and load them in Results."))
    })
    output$selection_summary <- shiny::renderText({
      meta <- scope()
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
      label <- paste(paste0("#", run_number), request$donor, request$celltype, request$readout,
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
      paste(c(request$donor, request$celltype, request$readout,
        paste(request$target, "minus", request$reference),
        paste(length(request$group_a_cells), "vs", length(request$group_b_cells), "cells in scope"),
        if (request$readout != "abundance") paste0("Minimum observed cells per condition: ", request$min_cells,
          "; isotype-background cutoff: ", signif(result()$provenance$filters$donor_isotype_fraction_cutoff, 4)),
        "BH adjustment applies to this requested feature family; frozen report results are unchanged."), collapse = " | ")
    })
    output$results_csv <- shiny::downloadHandler(
      filename = function() "proxiomevis_exploratory_comparison.csv",
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
