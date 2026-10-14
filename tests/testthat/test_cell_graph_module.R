testthat::test_that("GIF export controls expose bounded settings and accessible progress", {
  html <- as.character(cell_graph_ui("graph"))
  testthat::expect_match(html, 'id="graph-plot_gif_seconds"', fixed = TRUE)
  testthat::expect_match(html, 'id="graph-plot_gif_fps"', fixed = TRUE)
  testthat::expect_match(html, 'id="graph-plot_gif_pixels"', fixed = TRUE)
  testthat::expect_match(html, "Download rotating GIF", fixed = TRUE)
  testthat::expect_match(html, 'aria-label="GIF export progress"', fixed = TRUE)
  testthat::expect_match(html, 'aria-live="polite"', fixed = TRUE)
  testthat::expect_match(html, "cell-graph-gif.js", fixed = TRUE)
})

cell_graph_module_fixture <- function() {
  skip_if_not_installed("arrow")
  skip_if_not_installed("dplyr")
  root <- withr::local_tempdir(.local_envir = parent.frame())
  directory <- file.path(root, "graphs")
  dir.create(directory)
  data <- precomputed_app_data(precomputed_fixture())
  data$bundle_path <- file.path(root, "analysis.rds")
  data$assets <- list(graphs = list(path = "graphs", index = "index.rds"))
  data$cell_sources <- data.frame(component = data$metadata$component,
                                  original_id = paste0("source-", seq_len(nrow(data$metadata))),
                                  sample_id = data$metadata$sample_id)
  samples <- lapply(1:2, function(i) {
    mapping <- data$cell_sources[i, ]
    nodes <- data.frame(cell_id = mapping$component,
                        x = c(0, 1), y = c(1, 0), z = c(0, 0))
    nodes$markers <- list(character(), "A")
    temporary <- tempfile(tmpdir = directory)
    arrow::write_parquet(nodes, temporary)
    checksum <- unname(tools::md5sum(temporary))
    filename <- paste0(checksum, ".parquet")
    stopifnot(file.rename(temporary, file.path(directory, filename)))
    data.frame(sample = mapping$sample_id, file = filename, md5 = checksum,
               size_bytes = file.info(file.path(directory, filename))$size,
               source_pxl = paste0(mapping$sample_id, ".pxl"),
               source_pxl_md5 = paste(rep("a", 32), collapse = ""), layout = "wpmds_3d",
               n_nodes = 2L, n_cells = 1L)
  })
  index <- list(schema_version = 3L, dataset_id = data$contract$dataset_id,
                n_total_cells = nrow(data$metadata), coverage = "selected_cells",
                cells = data.frame(cell_id = data$cell_sources$component[1:2],
                  source_cell_id = data$cell_sources$original_id[1:2],
                  sample = data$cell_sources$sample_id[1:2], n_nodes = 2L),
                samples = do.call(rbind, samples))
  index_path <- file.path(directory, "index.rds")
  saveRDS(index, index_path)
  list(data = data, index = index, index_path = index_path)
}

flush_graph_outputs <- function(session) {
  session$flushReact()
  for (i in seq_len(3L)) {
    later::run_now(0)
    session$flushReact()
  }
}

test_that("a report-only dataset and missing graph directory give validation messages without ending the session", {
  fixture <- cell_graph_module_fixture()
  report_only <- fixture$data
  report_only$assets <- NULL
  current <- shiny::reactiveVal(report_only)
  shiny::testServer(cell_graph_server, args = list(data = current), {
    session$setInputs(donor = "D1", celltype = "B", condition = "NC", cell = "cell1", markers = "A")
    expect_error(index(), "no prepared 3D assets", class = "shiny.silent.error")
    expect_error(output$coverage, "no prepared 3D assets", class = "shiny.silent.error")
    expect_false(session$isClosed())
    missing <- fixture$data
    missing$assets$graphs$path <- "absent-graphs"
    current(missing)
    session$flushReact()
    expect_error(index(), class = "shiny.silent.error")
    expect_false(session$isClosed())
    current(fixture$data)
    flush_graph_outputs(session)
    session$setInputs(display = 1L)
    expect_identical(asset()$cell_id, "cell1")
    expect_match(output$summary, "2 nodes; wpmds_3d", fixed = TRUE)
    expect_match(output$summary, "marker identities; nodes only", fixed = TRUE)
    expect_no_error(output$plot)
    exported <- utils::read.csv(output$nodes_csv, stringsAsFactors = FALSE)
    expect_identical(names(exported), c("cell_id", "x", "y", "z", "markers"))
    expect_length(jsonlite::fromJSON(exported$markers[1L]), 0L)
    expect_identical(jsonlite::fromJSON(exported$markers[2L]), "A")
    provenance <- jsonlite::read_json(output$provenance_json)
    expect_equal(provenance$schema_version, 3L)
    expect_false(any(c("nodes", "counts", "edges") %in% names(provenance)))
    expect_false(session$isClosed())
    model <- output$model_glb
    expect_equal(rawToChar(readBin(model, "raw", n = 4)), "glTF")
    original_size <- file.info(model)$size
    expect_equal(powerpoint_appearance()$background_opacity, 1)
    expect_equal(powerpoint_appearance()$background_radius, 0.005)
    expect_equal(powerpoint_appearance()$highlight_radius, 0.007)
    session$setInputs(background_color = "#123456", highlight_color = "#00FF00",
      background_opacity = 0.4, highlight_opacity = 1, background_radius = 0.005, highlight_radius = 0.007,
      glb_background_opacity = 0.8, glb_highlight_opacity = 1, glb_background_radius = 0.008, glb_highlight_radius = 0.009,
      glb_background_retention = 50)
    expect_equal(glb_retention(), 25)
    expect_equal(appearance(), cell_graph_appearance())
    expect_equal(powerpoint_appearance()$background_radius, 0.005)
    session$setInputs(apply_appearance = 1L)
    expect_equal(appearance()$background_color, "#123456")
    expect_equal(glb_retention(), 50)
    expect_equal(appearance()$highlight_radius, 0.007)
    expect_equal(powerpoint_appearance()$background_color, appearance()$background_color)
    expect_equal(powerpoint_appearance()$background_opacity, 0.8)
    expect_equal(powerpoint_appearance()$highlight_radius, 0.009)
    session$setInputs(glb_background_radius = -1, apply_appearance = 2L)
    expect_equal(powerpoint_appearance()$background_radius, 0.008)
    expect_equal(appearance()$background_opacity, 0.4)
    session$setInputs(background_color = "invalid", apply_appearance = 3L)
    expect_equal(appearance()$background_color, "#123456")
    expect_gt(original_size, 192)
  })
})

test_that("a forged graph-to-bundle source mapping is rejected and refresh recovers", {
  fixture <- cell_graph_module_fixture()
  forged <- fixture$index
  forged$cells$source_cell_id[1L] <- "different-original-cell"
  saveRDS(forged, fixture$index_path)
  current <- shiny::reactiveVal(fixture$data)
  shiny::testServer(cell_graph_server, args = list(data = current), {
    session$setInputs(donor = "D1", celltype = "B", condition = "NC", cell = "cell1", markers = "A")
    expect_error(index(), "bundle's source cell mapping", class = "shiny.silent.error")
    expect_false(session$isClosed())
    saveRDS(fixture$index, fixture$index_path)
    session$setInputs(refresh = 1L)
    flush_graph_outputs(session)
    session$setInputs(display = 1L)
    expect_identical(asset()$source_cell_id, "source-1")
    expect_no_error(output$plot)
    expect_false(session$isClosed())
  })
})

test_that("oversized graph selection is explained and choosing a smaller cell restores the view", {
  fixture <- cell_graph_module_fixture()
  oversized <- fixture$index
  oversized$cells$n_nodes[1L] <- 200001L
  oversized$samples$n_nodes[1L] <- 200001L
  saveRDS(oversized, fixture$index_path)
  current <- shiny::reactiveVal(fixture$data)
  shiny::testServer(cell_graph_server, args = list(data = current), {
    session$setInputs(donor = "D1", celltype = "B", condition = "NC", cell = "cell1", markers = "A")
    session$setInputs(display = 1L)
    expect_error(asset(), "size limit", class = "shiny.silent.error")
    expect_error(output$summary, "no nodes were silently removed", class = "shiny.silent.error")
    expect_false(session$isClosed())
    session$setInputs(donor = "D2", cell = "cell2")
    flush_graph_outputs(session)
    session$setInputs(display = 2L)
    expect_identical(asset()$cell_id, "cell2")
    expect_equal(nrow(asset()$nodes), 2L)
    expect_match(output$summary, "2 nodes; wpmds_3d", fixed = TRUE)
    expect_no_error(output$plot)
    expect_false(session$isClosed())
  })
})

test_that("graph source hashes match the frozen inventory without reading any PXL file", {
  fixture <- cell_graph_module_fixture()
  fixture$data$contract$proximity_source$pxl_manifest <- data.frame(
    sample_alias = rev(fixture$index$samples$sample),
    pxl_md5 = rev(fixture$index$samples$source_pxl_md5),
    pxl_file = "/unmounted/no-runtime-access.pxl")
  current <- shiny::reactiveVal(fixture$data)
  shiny::testServer(cell_graph_server, args = list(data = current), {
    session$setInputs(donor = "D1", celltype = "B", condition = "NC", cell = "cell1", markers = "A")
    session$setInputs(display = 1L)
    expect_identical(asset()$cell_id, "cell1")
    expect_no_error(index())
    changed <- fixture$data
    changed$contract$proximity_source$pxl_manifest$pxl_md5[1L] <- paste(rep("b", 32), collapse = "")
    current(changed)
    session$flushReact()
    session$setInputs(display = 2L)
    expect_error(index(), "frozen PXL source checksums", class = "shiny.silent.error")
    expect_error(asset(), "frozen PXL source checksums", class = "shiny.silent.error")
    expect_false(session$isClosed())
    current(fixture$data)
    flush_graph_outputs(session)
    session$setInputs(display = 3L)
    expect_identical(asset()$cell_id, "cell1")
    expect_no_error(output$plot)
  })
})

test_that("only Display loads a graph and snapshots exactly one highlight", {
  fixture <- cell_graph_module_fixture()
  current <- shiny::reactiveVal(fixture$data)
  calls <- 0L
  loader <- load_cell_graph_asset
  server <- cell_graph_server
  environment(server) <- list2env(list(load_cell_graph_asset = function(...) {
    calls <<- calls + 1L
    loader(...)
  }), parent = environment(cell_graph_server))
  shiny::testServer(server, args = list(data = current), {
    session$setInputs(donor = "D1", celltype = "B", condition = "NC", cell = "cell1", markers = "A")
    expect_equal(calls, 0L)
    session$setInputs(display = 1L)
    expect_identical(asset()$cell_id, "cell1")
    expect_identical(displayed()$marker, "A")
    expect_equal(calls, 1L)
    session$setInputs(donor = "D2", cell = "cell2", markers = "B")
    expect_identical(asset()$cell_id, "cell1")
    expect_identical(displayed()$marker, "A")
    expect_equal(calls, 1L)
    expect_identical(jsonlite::read_json(output$provenance_json)$highlight_marker, "A")
    session$setInputs(markers = "A", display = 2L)
    expect_identical(asset()$cell_id, "cell2")
    expect_equal(calls, 2L)
    session$setInputs(markers = c("A", "B"), display = 3L)
    expect_error(asset(), "exactly one highlight", class = "shiny.silent.error")
    expect_equal(calls, 2L)
    session$setInputs(markers = "", display = 4L)
    expect_error(asset(), "exactly one highlight", class = "shiny.silent.error")
    expect_equal(calls, 2L)
  })
})

test_that("ranking controls preserve pair symmetry, describe missing graphs and wait for Display", {
  fixture <- cell_graph_module_fixture()
  fixture$data$assets$proximity <- list(path = "proximity")
  ranking_fixture(fixture$data, file.path(dirname(fixture$data$bundle_path), "proximity"))
  current <- shiny::reactiveVal(fixture$data)
  shiny::testServer(cell_graph_server, args = list(data = current), {
    session$setInputs(donor = "D1", celltype = "B", condition = "NC", pair_1 = "A", pair_2 = "B", rank_view = "top")
    expect_identical(ranked()$component, c("cell9", "cell1"))
    expect_identical(ranked()$graph_prepared, c(FALSE, TRUE))
    expect_match(names(cell_choices())[1L], "graph not prepared")
    expect_match(output$ranking_status, "Ranking complete: 2 of 2")
    ab <- ranked()
    session$setInputs(pair_1 = "B", pair_2 = "A")
    expect_identical(ranked(), ab)
    session$setInputs(rank_view = "bottom")
    expect_identical(shown_ranking()$component, c("cell1", "cell9"))
    session$setInputs(rank_view = "both")
    expect_identical(shown_ranking()$component, c("cell9", "cell1"))
    expect_equal(nrow(read.csv(output$ranking_csv)), 2L)
    session$setInputs(cell = "cell1", markers = "A", display = 1L)
    expect_identical(displayed()$ranking$cell$component, "cell1")
    expect_equal(displayed()$ranking$cell$rank, 2L)
    session$setInputs(pair_1 = "A", pair_2 = "A")
    expect_true(all(ranked()$marker_1 == ranked()$marker_2))
    expect_identical(displayed()$ranking$cell$marker_2, "B")
    session$setInputs(pair_1 = "C", pair_2 = "C")
    expect_length(cell_choices(), 0L)
    expect_match(output$ranking_status, "no finite observed score")
    expect_identical(asset()$cell_id, "cell1")
    session$setInputs(display = 2L)
    expect_error(asset(), "current list", class = "shiny.silent.error")
    session$setInputs(pair_1 = "A", pair_2 = "B", cell = "cell9", display = 3L)
    expect_error(asset(), "No prepared 3D asset", class = "shiny.silent.error")
    session$setInputs(pair_1 = "", pair_2 = "", cell = "cell1", display = 4L)
    expect_identical(asset()$cell_id, "cell1")
    expect_null(displayed()$ranking)
  })
})

test_that("missing or duplicate expected graph sources fail closed", {
  fixture <- cell_graph_module_fixture()
  inventory <- data.frame(sample_alias = fixture$index$samples$sample,
                          pxl_md5 = fixture$index$samples$source_pxl_md5)
  current <- shiny::reactiveVal(fixture$data)
  shiny::testServer(cell_graph_server, args = list(data = current), {
    session$setInputs(cell = "cell1", markers = "A")
    for (invalid in list(inventory[-1L, ], rbind(inventory, inventory[1L, ]),
                         inventory["sample_alias"], transform(inventory, pxl_md5 = NA_character_))) {
      changed <- fixture$data
      changed$contract$proximity_source$pxl_manifest <- invalid
      current(changed)
      session$flushReact()
      expect_error(index(), "frozen PXL source", class = "shiny.silent.error")
      expect_false(session$isClosed())
    }
  })
})

test_that("report handoff prefills controls without loading or redrawing until Display", {
  fixture <- cell_graph_module_fixture()
  fixture$data$assets$proximity <- list(path = "proximity")
  ranking_fixture(fixture$data, file.path(dirname(fixture$data$bundle_path), "proximity"))
  current <- shiny::reactiveVal(fixture$data)
  selection <- shiny::reactiveVal(NULL)
  loads <- plots <- 0L
  loader <- load_cell_graph_asset
  plotter <- plot_cell_graph_asset
  server <- cell_graph_server
  environment(server) <- list2env(list(
    load_cell_graph_asset = function(...) { loads <<- loads + 1L; loader(...) },
    plot_cell_graph_asset = function(...) { plots <<- plots + 1L; plotter(...) }
  ), parent = environment(cell_graph_server))
  messages <- list()
  mock <- shiny::MockShinySession$new()
  mock$sendInputMessage <- function(id, message) messages[[sub("^.*-", "", id)]] <<- message
  shiny::testServer(server, args = list(data = current, selection = selection), session = mock, {
    session$setInputs(donor = "D1", celltype = "B", condition = "NC", pair_1 = "", pair_2 = "",
                      cell = "cell1", markers = "A", display = 0L)
    session$setInputs(display = 1L)
    expect_identical(asset()$cell_id, "cell1")
    expect_no_error(output$plot)
    original_plots <- plots
    request <- list(donor = "D2", celltype = "B", marker_1 = "B", marker_2 = "A",
                    conditions = c("TRT", "NC"), token = 1L)
    selection(request)
    flush_graph_outputs(session)
    expect_true(handoff_pending())
    expect_identical(messages$donor$value, "D2")
    expect_identical(messages$celltype$value, "B")
    expect_identical(messages$pair_1$value, "B")
    expect_identical(messages$pair_2$value, "A")
    expect_identical(messages$condition$value, "")
    expect_match(as.character(messages$condition$options), "TRT")
    expect_match(as.character(messages$condition$options), "NC")
    expect_identical(messages$markers$value, "")
    expect_identical(asset()$cell_id, "cell1")
    expect_identical(displayed()$marker, "A")
    expect_no_error(output$plot)
    expect_equal(loads, 1L)
    expect_equal(plots, original_plots)
    expect_match(output$ranking_status, "Choose a condition")

    # Client updates arrive in separate flushes; old condition/highlight must not win.
    session$setInputs(donor = "D2")
    expect_true(handoff_pending())
    expect_identical(messages$condition$value, "")
    expect_error(ranked(), class = "shiny.silent.error")
    session$setInputs(celltype = "B", pair_1 = "B", pair_2 = "A", condition = "", markers = "")
    flush_graph_outputs(session)
    expect_false(handoff_pending())
    expect_identical(messages$condition$value, "")
    expect_match(as.character(messages$condition$options), "TRT")
    expect_equal(loads, 1L)
    expect_equal(plots, original_plots)
    session$setInputs(display = 2L)
    expect_error(asset(), "Choose a condition", class = "shiny.silent.error")
    expect_equal(loads, 1L)

    session$setInputs(condition = "NC", cell = "cell2", markers = "A")
    expect_identical(ranked()$component, c("cell10", "cell2"))
    expect_equal(loads, 1L)
    session$setInputs(display = 3L)
    expect_identical(asset()$cell_id, "cell2")
    expect_identical(displayed()$marker, "A")
    expect_equal(loads, 2L)
    expect_no_error(output$plot)
    updated_plots <- plots
    request$token <- 2L
    selection(request)
    flush_graph_outputs(session)
    expect_true(handoff_pending())
    expect_identical(messages$condition$value, "")
    expect_identical(messages$markers$value, "")
    expect_identical(asset()$cell_id, "cell2")
    expect_equal(loads, 2L)
    expect_equal(plots, updated_plots)
  })
})

test_that("handoff without sidecars remains friendly and does not load a graph", {
  fixture <- cell_graph_module_fixture()
  fixture$data$assets <- NULL
  current <- shiny::reactiveVal(fixture$data)
  selection <- shiny::reactiveVal(NULL)
  loads <- 0L
  server <- cell_graph_server
  environment(server) <- list2env(list(load_cell_graph_asset = function(...) {
    loads <<- loads + 1L
    stop("Unexpected graph read")
  }), parent = environment(cell_graph_server))
  messages <- list()
  mock <- shiny::MockShinySession$new()
  mock$sendInputMessage <- function(id, message) messages[[sub("^.*-", "", id)]] <<- message
  shiny::testServer(server, args = list(data = current, selection = selection), session = mock, {
    selection(list(donor = "D2", celltype = "B", marker_1 = "A", marker_2 = "A",
                   conditions = c("TRT", "NC"), token = 1L))
    flush_graph_outputs(session)
    expect_identical(messages$donor$value, "D2")
    expect_identical(messages$pair_1$value, "A")
    expect_identical(messages$pair_2$value, "A")
    expect_identical(messages$condition$value, "")
    expect_identical(messages$markers$value, "")
    expect_error(output$coverage, "no prepared 3D assets", class = "shiny.silent.error")
    session$setInputs(donor = "D2", celltype = "B", pair_1 = "A", pair_2 = "A", condition = NULL, markers = NULL)
    expect_false(handoff_pending())
    session$setInputs(condition = "NC")
    expect_error(ranked(), "Prepare a proximity sidecar", class = "shiny.silent.error")
    expect_equal(loads, 0L)
    expect_false(session$isClosed())
  })
})
