Sys.setenv(PROXIOME_ALLOW_SYSTEM_LIBS = "true")
source("../../app.R")

repo_root <- normalizePath(file.path(APP_DIR, "..", ".."), mustWork = TRUE)
app_renv_root <- normalizePath(APP_DIR, mustWork = TRUE)

expect_output_wrapped_with <- function(html, output_id, wrapper_class) {
  output_marker <- paste0('id="', output_id, '"')
  output_pos <- regexpr(output_marker, html, fixed = TRUE)[[1]]
  expect_true(output_pos > 0)

  preceding_html <- substr(html, max(1, output_pos - 2400), output_pos)
  wrapper_pattern <- paste0('class="[^"]*\\b', wrapper_class, '\\b[^"]*"')
  expect_true(grepl(wrapper_pattern, preceding_html, perl = TRUE))
}

test_that("readout description cards are not shown in the app chrome", {
  html <- htmltools::renderTags(ui)$html

  expect_false(grepl("readout-strip", html, fixed = TRUE))
  expect_false(grepl("Per-cell marker abundance from the PNA assay", html, fixed = TRUE))
  expect_false(grepl("Self-proximity for a marker within each cell graph", html, fixed = TRUE))
  expect_false(grepl("Marker-pair proximity between two different proteins", html, fixed = TRUE))
})

test_that("app name is consistent across UI and README", {
  html <- htmltools::renderTags(ui)$html
  readme_lines <- readLines(file.path(APP_DIR, "README.md"), warn = FALSE)
  readme <- paste(readme_lines, collapse = "\n")

  expect_true(grepl("ProxiomeVis", html, fixed = TRUE))
  expect_true(startsWith(readme, "# ProxiomeVis"))
  expect_false(grepl("Pixelgen Proxiome Explorer", html, fixed = TRUE))
  expect_false(grepl("# Pixelgen Proxiome Shiny Demo", readme, fixed = TRUE))
  expect_false(any(grepl("# ProxiomeVis", readme_lines[-1], fixed = TRUE)))
})

test_that("each readout tab has only the controls it needs", {
  html <- htmltools::renderTags(ui)$html

  expect_true(grepl('data-value="QC"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_mode"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_sample_filter"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_metric"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_n_umi_cutoff"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_isotype_cutoff"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_filter_include_total"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_filter_y"', html, fixed = TRUE))
  expect_false(grepl('id="pixelator_output_dir"', html, fixed = TRUE))
  expect_false(grepl('data-value="Sequencing"', html, fixed = TRUE))
  expect_false(grepl('data-value="Cell Recovery"', html, fixed = TRUE))
  expect_false(grepl('data-value="Graph Metrics"', html, fixed = TRUE))
  expect_false(grepl('data-value="Control Markers"', html, fixed = TRUE))
  expect_false(grepl('id="qc_sequencing_plot"', html, fixed = TRUE))
  expect_false(grepl('id="qc_sequencing_saturation_plot"', html, fixed = TRUE))
  expect_false(grepl('id="qc_cell_recovery_plot"', html, fixed = TRUE))
  expect_false(grepl('id="qc_graph_metrics_plot"', html, fixed = TRUE))
  expect_false(grepl('id="qc_control_markers_plot"', html, fixed = TRUE))

  expect_true(grepl('id="abundance-abundance_embedding"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_color_by"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_marker"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_distribution_marker"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_distribution_columns"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_distribution_width"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_distribution_height"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_distribution_show_jitter"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_split_by"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_point_size"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_marker"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_heatmap_marker_count"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_per_marker_plot"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_summary_heatmap"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_heatmap_markers"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_heatmap_display"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_heatmap_preset"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-spatial_coloc_scope"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-spatial_celltype_focus"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-spatial_marker_selection_mode"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-spatial_top_marker_count"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-spatial_min_pct_detected"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-spatial_min_log2_range"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_reference_condition"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_clustering_method"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_legend_min"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_legend_max"', html, fixed = TRUE))

  expect_false(grepl('id="anchor_marker"', html, fixed = TRUE))
  expect_false(grepl('id="colocalization_anchor_marker"', html, fixed = TRUE))
  expect_false(grepl('id="embedding"', html, fixed = TRUE))
  expect_false(grepl('id="marker"', html, fixed = TRUE))

  abundance_section <- sub(".*Abundance", "", html)
  abundance_section <- sub("Clustering.*", "", abundance_section)
  clustering_section <- sub(".*Clustering", "", html)
  clustering_section <- sub("Spatial Metrics.*", "", clustering_section)

  expect_false(grepl("Colocalization anchor", abundance_section, fixed = TRUE))
  expect_false(grepl("Colocalization anchor", clustering_section, fixed = TRUE))
})

test_that("spatial metrics tab wires PixelatorES-style heatmap controls", {
  html <- htmltools::renderTags(ui)$html
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_true(grepl('data-value="Spatial Metrics"', html, fixed = TRUE))
  expect_false(grepl('fillable = c("QC", "Abundance", "Clustering", "Spatial Metrics", "Environment")', app_source, fixed = TRUE))
  expect_true(grepl('id = "spatial_metric_readout"', app_source, fixed = TRUE))
  clustering_module_source <- paste(readLines(file.path(APP_DIR, "R", "clustering_module.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('nav_panel\\(\\s*"Clustering"', clustering_module_source, perl = TRUE))
  colocalization_module_source <- paste(readLines(file.path(APP_DIR, "R", "colocalization_module.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('nav_panel\\(\\s*"Colocalization"', colocalization_module_source, perl = TRUE))
  expect_match(app_source, 'source\\(file\\.path\\(APP_DIR, "R", "spatial_metrics\\.R"\\)\\)')
  expect_true(grepl("summarize_spatial_heatmap_by_sample", app_source, fixed = TRUE))
  expect_true(grepl("summarize_spatial_heatmap_by_celltype", app_source, fixed = TRUE))
  expect_true(grepl("select_spatial_heatmap_markers", app_source, fixed = TRUE))
  expect_true(grepl("complete_spatial_marker_pairs", colocalization_module_source, fixed = TRUE))
  expect_true(grepl('else "sample_alias"', colocalization_module_source, fixed = TRUE))
})

test_that("each readout exposes observed and differential modes", {
  html <- htmltools::renderTags(ui)$html

  expect_true(grepl('id="abundance-abundance_mode"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_mode"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_mode"', html, fixed = TRUE))
  expect_true(grepl('data-value="Marker Distributions"', html, fixed = TRUE))
  expect_true(grepl('data-value="Cell Annotation"', html, fixed = TRUE))
  expect_true(grepl('data-value="Per Marker"', html, fixed = TRUE))
  expect_true(grepl('data-value="Summary Heatmap"', html, fixed = TRUE))

  expect_true(grepl('id="abundance-abundance_diff_group_a"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_diff_group_b"', html, fixed = TRUE))
  expect_true(grepl("Group B (reference)", html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_diff_fdr"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_diff_effect"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_diff_marker"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_diff_stratify_celltype"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_run_differential"', html, fixed = TRUE))

  expect_true(grepl('id="clustering-clustering_diff_group_a"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_diff_group_b"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_diff_fdr"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_diff_effect"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_diff_marker"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_diff_stratify_celltype"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_run_differential"', html, fixed = TRUE))

  expect_true(grepl('id="colocalization-colocalization_diff_group_a"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_diff_group_b"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_diff_fdr"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_diff_effect"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_diff_anchor_marker"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_diff_pair"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_diff_stratify_celltype"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_run_differential"', html, fixed = TRUE))
})

test_that("differential views include plot, detail, and table outputs", {
  html <- htmltools::renderTags(ui)$html

  expect_true(grepl('id="qc-qc_filter_plot"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_molecule_rank_plot"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_distribution_plot"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_filter_table"', html, fixed = TRUE))
  expect_true(grepl('id="qc-qc_origin_metadata_table"', html, fixed = TRUE))

  expect_true(grepl('id="abundance-abundance_marker_distribution_plot_ui"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_celltype_composition_plot"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_annotation_heatmap"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_diff_volcano"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_diff_detail"', html, fixed = TRUE))
  expect_true(grepl('id="abundance-abundance_diff_table"', html, fixed = TRUE))

  expect_true(grepl('id="clustering-clustering_per_marker_plot"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_summary_heatmap"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_diff_volcano"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_diff_detail"', html, fixed = TRUE))
  expect_true(grepl('id="clustering-clustering_diff_table"', html, fixed = TRUE))
  expect_false(grepl('id="clustering_diff_rank"', html, fixed = TRUE))

  expect_true(grepl('id="colocalization-colocalization_diff_volcano"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_diff_detail"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_diff_table"', html, fixed = TRUE))
  expect_false(grepl('id="colocalization_diff_heatmap"', html, fixed = TRUE))
})

test_that("differential volcano and box plots share a side-by-side row", {
  html <- htmltools::renderTags(ui)$html
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_true(grepl("differential-plot-row", html, fixed = TRUE))
  expect_true(grepl(".differential-plot-row", app_source, fixed = TRUE))
  expect_true(grepl("grid-template-columns: max-content max-content;", app_source, fixed = TRUE))
  expect_true(grepl("overflow-x: auto;", app_source, fixed = TRUE))
  expect_true(grepl("@media (max-width: 1100px)", app_source, fixed = TRUE))
  abundance_module_source <- paste(readLines(file.path(APP_DIR, "R", "abundance_module.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('differential_plot_row(ns("abundance_diff_volcano"), ns("abundance_diff_detail"))', abundance_module_source, fixed = TRUE))
  clustering_module_source <- paste(readLines(file.path(APP_DIR, "R", "clustering_module.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('differential_plot_row(ns("clustering_diff_volcano"), ns("clustering_diff_detail"))', clustering_module_source, fixed = TRUE))
  colocalization_module_source <- paste(readLines(file.path(APP_DIR, "R", "colocalization_module.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl('differential_plot_row(ns("colocalization_diff_volcano"), ns("colocalization_diff_detail"))', colocalization_module_source, fixed = TRUE))
})

test_that("plot outputs are wrapped in readable width classes", {
  html <- htmltools::renderTags(ui)$html

  compact_outputs <- c(
    "qc-qc_filter_plot",
    "qc-qc_distribution_plot",
    "abundance-abundance_diff_volcano",
    "clustering-clustering_diff_volcano",
    "colocalization-colocalization_diff_volcano"
  )
  standard_outputs <- c(
    "abundance-abundance_umap",
    "clustering-clustering_plot",
    "abundance-abundance_diff_detail",
    "clustering-clustering_diff_detail",
    "colocalization-colocalization_diff_detail"
  )
  wide_outputs <- c("qc-qc_molecule_rank_plot")
  scroll_outputs <- c("abundance-abundance_marker_distribution_plot_ui", "colocalization-colocalization_heatmap_interactive", "colocalization-colocalization_heatmap_original")

  for (output_id in compact_outputs) {
    expect_output_wrapped_with(html, output_id, "plot-pane-compact")
  }
  for (output_id in standard_outputs) {
    expect_output_wrapped_with(html, output_id, "plot-pane-standard")
  }
  for (output_id in wide_outputs) {
    expect_output_wrapped_with(html, output_id, "plot-pane-wide")
  }
  for (output_id in scroll_outputs) {
    expect_output_wrapped_with(html, output_id, "plot-pane-scroll")
  }
})

test_that("figure panes expose PNG and SVG downloads", {
  html <- htmltools::renderTags(ui)$html
  downloadable_outputs <- c(
    "qc-qc_filter_plot",
    "qc-qc_molecule_rank_plot",
    "qc-qc_distribution_plot",
    "abundance-abundance_umap",
    "abundance-abundance_marker_distribution_plot",
    "abundance-abundance_celltype_composition_plot",
    "abundance-abundance_annotation_heatmap",
    "abundance-abundance_diff_volcano",
    "abundance-abundance_diff_detail",
    "clustering-clustering_plot",
    "clustering-clustering_per_marker_plot",
    "clustering-clustering_summary_heatmap",
    "clustering-clustering_diff_volcano",
    "clustering-clustering_diff_detail",
    "colocalization-colocalization_heatmap",
    "colocalization-colocalization_diff_volcano",
    "colocalization-colocalization_diff_detail"
  )

  for (output_id in downloadable_outputs) {
    expect_true(grepl(paste0('id="', output_id, '_download_png"'), html, fixed = TRUE))
    expect_true(grepl(paste0('id="', output_id, '_download_svg"'), html, fixed = TRUE))
  }
})

test_that("plot download helper writes PNG and SVG files from ggplot objects", {
  plot <- ggplot(data.frame(x = 1:3, y = c(1, 4, 2)), aes(x, y)) +
    geom_line()
  png_path <- tempfile(fileext = ".png")
  svg_path <- tempfile(fileext = ".svg")

  save_ggplot_download(plot, png_path, format = "png", width = 3, height = 2, dpi = 96)
  save_ggplot_download(plot, svg_path, format = "svg", width = 3, height = 2, dpi = 96)

  expect_true(file.exists(png_path))
  expect_true(file.exists(svg_path))
  expect_gt(file.info(png_path)$size, 100)
  expect_match(paste(readLines(svg_path, warn = FALSE), collapse = "\n"), "<svg", fixed = TRUE)
})

test_that("app uses a bslib navbar page instead of custom top chrome", {
  html <- htmltools::renderTags(ui)$html

  expect_true(grepl("navbar", html, fixed = TRUE))
  expect_true(grepl("navbar-brand", html, fixed = TRUE))
  expect_false(grepl("app-topbar", html, fixed = TRUE))
})

test_that("data loading is isolated behind the Data Source module", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_true(file.exists(file.path(APP_DIR, "R", "data_source_module.R")))
  expect_match(app_source, 'source\\(file\\.path\\(APP_DIR, "R", "data_source_module\\.R"\\), local = TRUE\\)')
  expect_true(grepl('data_source_module_ui("data_source")', app_source, fixed = TRUE))
  expect_true(grepl('data_source <- data_source_module_server("data_source", app_dir = APP_DIR)', app_source, fixed = TRUE))
  expect_true(grepl("demo_data <- data_source$data", app_source, fixed = TRUE))
  expect_false(grepl("send_rds_load_state <- function", app_source, fixed = TRUE))
  expect_false(grepl("load_demo_into_app <- function", app_source, fixed = TRUE))
})

test_that("QC readout is isolated behind the QC module", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_true(file.exists(file.path(APP_DIR, "R", "qc_module.R")))
  expect_match(app_source, 'source\\(file\\.path\\(APP_DIR, "R", "qc_module\\.R"\\), local = TRUE\\)')
  expect_true(grepl('qc_module_ui("qc")', app_source, fixed = TRUE))
  expect_true(grepl('qc_module_server("qc", data = demo_data)', app_source, fixed = TRUE))
  expect_false(grepl("qc_sidebar <- function", app_source, fixed = TRUE))
  expect_false(grepl("output$qc_", app_source, fixed = TRUE))
  expect_false(grepl("qc_sample_choices <- function", app_source, fixed = TRUE))
  expect_false(grepl("qc_molecule_rank_plotly <- function", app_source, fixed = TRUE))
})

test_that("Abundance readout is isolated behind the Abundance module", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_true(file.exists(file.path(APP_DIR, "R", "abundance_module.R")))
  expect_match(app_source, 'source\\(file\\.path\\(APP_DIR, "R", "abundance_module\\.R"\\), local = TRUE\\)')
  expect_true(grepl('abundance_module_ui("abundance")', app_source, fixed = TRUE))
  expect_true(grepl('abundance_module_server("abundance", data = demo_data)', app_source, fixed = TRUE))
  expect_false(grepl("abundance_sidebar <- function", app_source, fixed = TRUE))
  expect_false(grepl("abundance_diff_config <- reactiveVal", app_source, fixed = TRUE))
  expect_false(grepl("output$abundance_", app_source, fixed = TRUE))
  expect_false(grepl("output$metric_row", app_source, fixed = TRUE))
  expect_false(grepl("abundance_points <- reactive", app_source, fixed = TRUE))
  expect_false(grepl("abundance_distribution_data <- reactive", app_source, fixed = TRUE))
  expect_false(grepl("plot_abundance_marker_distribution <- function", app_source, fixed = TRUE))
})

test_that("Clustering readout is isolated behind the Clustering module", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_true(file.exists(file.path(APP_DIR, "R", "clustering_module.R")))
  expect_match(app_source, 'source\\(file\\.path\\(APP_DIR, "R", "clustering_module\\.R"\\), local = TRUE\\)')
  expect_true(grepl('clustering_module_ui("clustering")', app_source, fixed = TRUE))
  expect_true(grepl('clustering_module_server("clustering", data = demo_data)', app_source, fixed = TRUE))
  expect_false(grepl("clustering_sidebar <- function", app_source, fixed = TRUE))
  expect_false(grepl("clustering_diff_config <- reactiveVal", app_source, fixed = TRUE))
  expect_false(grepl("output$clustering_", app_source, fixed = TRUE))
  expect_false(grepl("clustering_points <- reactive", app_source, fixed = TRUE))
  expect_false(grepl("clustering_heatmap_summary <- reactive", app_source, fixed = TRUE))
  expect_false(grepl("plot_clustering_per_marker <- function", app_source, fixed = TRUE))
  expect_false(grepl("plot_clustering_summary_heatmap <- function", app_source, fixed = TRUE))
})

test_that("Colocalization readout is isolated behind the Colocalization module", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_true(file.exists(file.path(APP_DIR, "R", "colocalization_module.R")))
  expect_match(app_source, 'source\\(file\\.path\\(APP_DIR, "R", "colocalization_module\\.R"\\), local = TRUE\\)')
  expect_true(grepl('colocalization_module_ui("colocalization")', app_source, fixed = TRUE))
  expect_true(grepl('colocalization_module_server("colocalization", data = demo_data)', app_source, fixed = TRUE))
  expect_false(grepl("colocalization_sidebar <- function", app_source, fixed = TRUE))
  expect_false(grepl("colocalization_diff_config <- reactiveVal", app_source, fixed = TRUE))
  expect_false(grepl("output$colocalization_", app_source, fixed = TRUE))
  expect_false(grepl("colocalization_metadata <- reactive", app_source, fixed = TRUE))
  expect_false(grepl("colocalization_heatmap_result <- reactive", app_source, fixed = TRUE))
  expect_false(grepl("colocalization_diff_results <- reactive", app_source, fixed = TRUE))
})

test_that("Differential helpers are shared outside app.R", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")
  helper_source <- paste(readLines(file.path(APP_DIR, "R", "differential_helpers.R"), warn = FALSE), collapse = "\n")

  expect_true(file.exists(file.path(APP_DIR, "R", "differential_helpers.R")))
  expect_match(app_source, 'source\\(file\\.path\\(APP_DIR, "R", "differential_helpers\\.R"\\), local = TRUE\\)')
  expect_true(grepl("differential_plot_row <- function", helper_source, fixed = TRUE))
  expect_true(grepl("make_differential_config <- function", helper_source, fixed = TRUE))
  expect_true(grepl("default_differential_config <- function", helper_source, fixed = TRUE))
  expect_true(grepl("differential_volcano_plot <- function", helper_source, fixed = TRUE))
  expect_true(grepl("format_differential_table <- function", helper_source, fixed = TRUE))
  expect_false(grepl("differential_plot_row <- function", app_source, fixed = TRUE))
  expect_false(grepl("make_differential_config <- function", app_source, fixed = TRUE))
  expect_false(grepl("default_differential_config <- function", app_source, fixed = TRUE))
  expect_false(grepl("differential_volcano_plot <- function", app_source, fixed = TRUE))
  expect_false(grepl("format_differential_table <- function", app_source, fixed = TRUE))
})

test_that("project includes renv deployment infrastructure", {
  lock_path <- file.path(app_renv_root, "renv.lock")
  activate_path <- file.path(app_renv_root, "renv", "activate.R")
  settings_path <- file.path(app_renv_root, "renv", "settings.json")
  profile_path <- file.path(app_renv_root, ".Rprofile")
  renvignore_path <- file.path(app_renv_root, ".renvignore")
  gitignore_path <- file.path(app_renv_root, ".gitignore")

  expect_true(file.exists(lock_path))
  expect_true(file.exists(activate_path))
  expect_true(file.exists(settings_path))
  expect_true(file.exists(profile_path))
  expect_true(file.exists(renvignore_path))
  expect_false(file.exists(file.path(repo_root, "renv.lock")))

  gitignore_entries <- trimws(readLines(gitignore_path, warn = FALSE))
  expect_false("renv/" %in% gitignore_entries)
  expect_true("renv/library/" %in% gitignore_entries)
  expect_true("renv/local/" %in% gitignore_entries)
  expect_true("renv/staging/" %in% gitignore_entries)

  lock <- jsonlite::fromJSON(lock_path, simplifyVector = FALSE)
  expect_true("Packages" %in% names(lock))
  expect_true(all(c("shiny", "bslib", "ggplot2", "plotly", "future", "promises", "rlang", "data.table", "Seurat", "SeuratObject", "pixelatorR", "callr") %in% names(lock$Packages)))
  expect_true("renv" %in% names(lock$Packages))
})

test_that("app activates shared renv before loading packages", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_match(app_source, "activate_project_renv")
  expect_lt(
    regexpr("activate_project_renv", app_source, fixed = TRUE)[[1]],
    regexpr("library\\(bslib\\)", app_source)[[1]]
  )
  expect_false(grepl("renv::restore", app_source, fixed = TRUE))
  expect_false(grepl("install.packages", app_source, fixed = TRUE))
})

test_that("platform detection recognizes current HPC, Biowulf, and portable runtimes", {
  expect_equal(detect_proxiome_platform("ncifcrf.gov"), "ccrsf_hpc")
  expect_equal(detect_proxiome_platform("fsitgl-head01p.ncifcrf.gov"), "ccrsf_hpc")
  expect_equal(detect_proxiome_platform("biowulf.nih.gov"), "biowulf_hpc")
  expect_equal(detect_proxiome_platform("laptop.local"), "portable")
  expect_equal(detect_proxiome_platform("random-host", env_platform = "shinyapps"), "shinyapps")
  withr::local_envvar(c(SHINYAPPS_ACCOUNT = "demo"))
  expect_equal(detect_proxiome_platform("random-host"), "shinyapps")

  withr::local_envvar(c(PROXIOME_ALLOW_SYSTEM_LIBS = NA))
  expect_true(proxiome_platform_requires_shared_renv("ccrsf_hpc"))
  expect_true(proxiome_platform_requires_shared_renv("biowulf_hpc"))
  expect_false(proxiome_platform_requires_shared_renv("portable"))
  expect_false(proxiome_platform_requires_shared_renv("shinyapps"))
})

test_that("portable runtimes ignore restored renv libraries for a different R minor version", {
  project_root <- tempfile("mismatched-renv-")
  wrong_minor <- if (identical(current_r_minor_version(), "9.9")) "9.8" else "9.9"
  library_root <- file.path(project_root, "renv", "library", "linux-test", paste0("R-", wrong_minor), "x86_64-test")
  for (package in required_app_packages(include_renv = TRUE)) {
    dir.create(file.path(library_root, package), recursive = TRUE, showWarnings = FALSE)
  }

  expect_false(library_matches_current_r(library_root))
  expect_equal(restored_renv_library_paths(project_root), character(0))
})

test_that("portable runtimes can use installed packages when app renv is absent", {
  project_root <- tempfile("portable-proxiome-")
  dir.create(project_root, recursive = TRUE)

  libs <- activate_project_renv(project_root, platform = "portable")

  expect_true(length(libs) > 0)
  expect_true(all(c("shiny", "bslib", "ggplot2", "plotly", "future", "promises", "rlang", "data.table") %in% rownames(installed.packages(lib.loc = libs))))
})

test_that("renv restored library detection supports platform library layouts", {
  library_root <- file.path(tempdir(), "renv-library-layout")
  package_root <- file.path(library_root, "linux-ol-8.10", paste0("R-", current_r_minor_version()), "x86_64-redhat-linux-gnu")
  packages <- c("renv", "shiny", "bslib", "ggplot2", "plotly", "future", "promises", "rlang", "data.table")
  vapply(packages, function(package) {
    dir.create(file.path(package_root, package), recursive = TRUE, showWarnings = FALSE)
  }, logical(1))

  expect_true(renv_library_is_restored(library_root))
})

test_that("renv restored library resolution ignores incomplete platform libraries", {
  project_root <- tempfile("mixed-renv-library-layout-")
  library_root <- file.path(project_root, "renv", "library")
  complete_root <- file.path(library_root, "linux-ol-8.10", paste0("R-", current_r_minor_version()), "x86_64-pc-linux-gnu")
  incomplete_root <- file.path(library_root, "linux-ol-8.10", paste0("R-", current_r_minor_version()), "x86_64-redhat-linux-gnu")
  packages <- c("renv", "shiny", "bslib", "ggplot2", "plotly", "future", "promises", "rlang", "data.table")

  vapply(packages, function(package) {
    dir.create(file.path(complete_root, package), recursive = TRUE, showWarnings = FALSE)
  }, logical(1))
  vapply(setdiff(packages, "promises"), function(package) {
    dir.create(file.path(incomplete_root, package), recursive = TRUE, showWarnings = FALSE)
  }, logical(1))

  withr::local_envvar(c(RENV_PATHS_LIBRARY = library_root))
  restored_libraries <- restored_renv_library_paths(project_root)

  expect_equal(restored_libraries, normalizePath(complete_root, mustWork = TRUE))
  expect_false(incomplete_root %in% restored_libraries)
})

expect_environment_nav_visible <- function(html) {
  expect_true(grepl('data-value="Environment"', html, fixed = TRUE))
  expect_true(grepl('id="environment-r_paths"', html, fixed = TRUE))
  expect_true(grepl('id="environment-package_paths"', html, fixed = TRUE))
  expect_true(grepl('id="environment-package_filter"', html, fixed = TRUE))
  expect_true(grepl("R Paths", html, fixed = TRUE))
  expect_true(grepl("Packages", html, fixed = TRUE))
}

expect_environment_nav_hidden <- function(html) {
  expect_false(grepl('data-value="Environment"', html, fixed = TRUE))
  expect_false(grepl('id="environment-r_paths"', html, fixed = TRUE))
  expect_false(grepl('id="environment-package_paths"', html, fixed = TRUE))
  expect_false(grepl('id="environment-package_filter"', html, fixed = TRUE))
  expect_false(grepl("R Paths", html, fixed = TRUE))
  expect_false(grepl("Packages", html, fixed = TRUE))
}

test_that("environment diagnostics are hidden unless the debug env var is enabled", {
  html <- htmltools::renderTags(ui)$html
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")
  environment_source <- paste(readLines(file.path(APP_DIR, "R", "environment_module.R"), warn = FALSE), collapse = "\n")

  expect_false(environment_diagnostics_enabled())
  expect_environment_nav_hidden(html)
  expect_true(grepl('PROXIOME_SHOW_ENVIRONMENT', environment_source, fixed = TRUE))
  expect_true(grepl('if (environment_diagnostics_enabled())', app_source, fixed = TRUE))

  withr::local_envvar(c(PROXIOME_SHOW_ENVIRONMENT = "true"))
  expect_true(environment_diagnostics_enabled())
  expect_environment_nav_visible(htmltools::renderTags(build_app_ui())$html)
})

test_that("environment diagnostics list R runtime and package paths", {
  runtime <- environment_runtime_paths(APP_DIR)

  expect_s3_class(runtime, "data.frame")
  expect_named(runtime, c("item", "value"))
  expect_true(all(c("R", "Rscript", "R.home", ".libPaths()", "RENV_PROJECT", "RENV_PATHS_LIBRARY") %in% runtime$item))
  expect_true(app_renv_root %in% runtime$value)
  expect_true(file.path(app_renv_root, "renv", "library") %in% runtime$value)
  expect_true(any(.libPaths() %in% runtime$value))

  packages <- environment_package_paths()

  expect_s3_class(packages, "data.frame")
  expect_true(all(c("Package", "Version", "LibPath") %in% names(packages)))
  expect_true("shiny" %in% packages$Package)
  expect_true(any(packages$LibPath %in% .libPaths()))
})

test_that("environment module receives app_dir explicitly", {
  runtime_default <- paste(deparse(formals(environment_runtime_paths)$app_dir), collapse = "")
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")

  expect_true("app_dir" %in% names(formals(environment_runtime_paths)))
  expect_true("app_dir" %in% names(formals(environment_module_server)))
  expect_false(grepl("APP_DIR", runtime_default, fixed = TRUE))
  expect_true(grepl('source(file.path(APP_DIR, "R", "environment_module.R"), local = TRUE)', app_source, fixed = TRUE))
  expect_true(grepl('environment_module_server("environment", app_dir = APP_DIR)', app_source, fixed = TRUE))
})

test_that("environment package diagnostics include the app renv library explicitly", {
  module_library <- "/mnt/nasapps/production/R/4.5.0/lib64/R/library"
  old_lib_paths <- .libPaths()
  on.exit(.libPaths(old_lib_paths), add = TRUE)

  if (dir.exists(module_library)) {
    .libPaths(module_library)
  }

  packages <- environment_package_paths(app_dir = APP_DIR)
  app_package_rows <- packages[
    packages$Package %in% c("renv", "pixelatorR", "Seurat", "plotly", "bslib", "shiny"),
    ,
    drop = FALSE
  ]

  expect_true(any(startsWith(app_package_rows$LibPath, file.path(app_renv_root, "renv", "library"))))
  expect_true("renv" %in% packages$Package)
})

test_that("navbar exposes server-side RDS path loading on HPC and desktop only", {
  html <- htmltools::renderTags(ui)$html
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")
  data_source_source <- paste(readLines(file.path(APP_DIR, "R", "data_source_module.R"), warn = FALSE), collapse = "\n")

  expect_true(user_rds_path_loading_enabled("ccrsf_hpc"))
  expect_true(user_rds_path_loading_enabled("biowulf_hpc"))
  expect_true(user_rds_path_loading_enabled("portable"))
  expect_false(user_rds_path_loading_enabled("shinyapps"))

  expect_true(grepl('id="data_source-data_source_menu"', html, fixed = TRUE))
  expect_false(grepl('id="upload_rds"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-rds_server_path"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-validate_rds_path"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-load_rds_path"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-use_demo_data"', html, fixed = TRUE))
  expect_false(grepl("Upload RDS", html, fixed = TRUE))
  expect_true(grepl("RDS path", html, fixed = TRUE))
  expect_true(grepl("Validate RDS", html, fixed = TRUE))
  expect_true(grepl("Load Data", html, fixed = TRUE))
  expect_true(grepl("Use demo data", html, fixed = TRUE))
  expect_true(grepl('id="data_source-source_summary"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-rds_schema_report"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-rds_load_status"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-rds_load_progress"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-rds_load_progress_bar"', html, fixed = TRUE))
  expect_true(grepl('id="data_source-rds_load_elapsed"', html, fixed = TRUE))
  expect_true(grepl("rds-load-status idle", html, fixed = TRUE))
  expect_true(grepl("rds-load-progress idle", html, fixed = TRUE))
  expect_true(grepl("Enter an RDS path, then click Load Data.", html, fixed = TRUE))

  expect_false(grepl("observeEvent(input$upload_rds", app_source, fixed = TRUE))
  expect_false(grepl("observeEvent(input$load_rds_path", app_source, fixed = TRUE))
  expect_true(grepl("observeEvent(input$validate_rds_path", data_source_source, fixed = TRUE))
  expect_true(grepl("observeEvent(input$load_rds_path", data_source_source, fixed = TRUE))
  expect_true(grepl("inspect_user_rds_schema(input$rds_server_path)", data_source_source, fixed = TRUE))
  expect_true(grepl("format_user_rds_schema_report", data_source_source, fixed = TRUE))
  expect_false(grepl("observeEvent(input$upload_rds_too_large", app_source, fixed = TRUE))
})

test_that("server RDS path loading can be disabled for hosted runtimes", {
  html <- htmltools::renderTags(data_source_controls("data_source", platform = "shinyapps"))$html

  expect_false(grepl('id="upload_rds"', html, fixed = TRUE))
  expect_false(grepl('id="rds_server_path"', html, fixed = TRUE))
  expect_false(grepl('id="load_rds_path"', html, fixed = TRUE))
})

test_that("server RDS path validation requires an existing RDS file", {
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(list(ok = TRUE), rds_path)

  expect_true(validate_rds_file_path(rds_path))
  expect_error(validate_rds_file_path(tempfile(fileext = ".txt")), "must be an .rds file")
  expect_error(validate_rds_file_path(tempfile(fileext = ".rds")), "does not exist")
  expect_error(validate_rds_file_path(""), "Enter an RDS path")
})

test_that("RDS load progress records stage, percent, and elapsed time", {
  progress_path <- tempfile(fileext = ".rds")
  started_at <- Sys.time() - 75

  write_rds_load_progress(
    progress_path,
    state = "running",
    stage = "proximity",
    message = "Reading stored proximity scores...",
    value = 0.65,
    started_at = started_at
  )
  progress <- read_rds_load_progress(progress_path)

  expect_equal(progress$state, "running")
  expect_equal(progress$stage, "proximity")
  expect_equal(progress$value, 0.65)
  expect_equal(progress$percent, 65)
  expect_match(format_elapsed_time(progress$elapsed_seconds), "1m")

  message <- rds_load_progress_message(progress)
  expect_match(message, "Reading stored proximity scores", fixed = TRUE)
  expect_match(message, "elapsed", ignore.case = TRUE)
})

test_that("running RDS load elapsed display advances between progress writes", {
  started_at <- as.POSIXct("2026-06-12 10:00:00", tz = "UTC")
  progress <- list(
    state = "running",
    stage = "read_rds",
    message = "Reading RDS file...",
    value = 0.18,
    percent = 18,
    started_at = started_at,
    elapsed_seconds = 0
  )

  expect_equal(
    rds_load_elapsed_label(progress, now = started_at + 65),
    "Elapsed: 1m 05s"
  )
  expect_match(
    rds_load_progress_message(progress, now = started_at + 65),
    "Elapsed: 1m 05s",
    fixed = TRUE
  )
})

test_that("user RDS cache path is keyed by source file metadata", {
  home_dir <- tempfile("proxiomevis-home-")
  dir.create(home_dir, recursive = TRUE)
  withr::local_envvar(c(
    PROXIOMEVIS_HOME = home_dir,
    PROXIOME_RUNTIME_DIR = NA
  ))

  rds_path <- tempfile(fileext = ".rds")
  saveRDS(list(version = 1), rds_path)
  first_cache_path <- user_rds_cache_path(rds_path)

  Sys.sleep(1.1)
  saveRDS(list(version = 2), rds_path)
  second_cache_path <- user_rds_cache_path(rds_path)

  expect_true(startsWith(basename(first_cache_path), "user-rds-"))
  expect_equal(dirname(first_cache_path), proxiomevis_cache_dir())
  expect_false(identical(first_cache_path, second_cache_path))
})

test_that("server RDS path loading refreshes app-local libraries before reading", {
  rds_path <- tempfile(fileext = ".rds")
  saveRDS(list(ok = TRUE), rds_path)
  wrong_library <- tempfile("r-lib-")
  dir.create(wrong_library, recursive = TRUE)

  old_lib_paths <- .libPaths()
  on.exit(.libPaths(old_lib_paths), add = TRUE)

  original_loader <- load_demo_proxiome_data
  on.exit(assign("load_demo_proxiome_data", original_loader, envir = globalenv()), add = TRUE)

  captured_lib_paths <- character()
  captured_loader_args <- list()
  captured_cache_path <- NULL
  captured_force <- NULL
  assign(
    "load_demo_proxiome_data",
    function(rds_path, cache_path, force, ...) {
      captured_lib_paths <<- .libPaths()
      captured_loader_args <<- list(...)
      captured_cache_path <<- cache_path
      captured_force <<- force
      list(source = list(rds_path = rds_path, n_cells = 1L, n_markers = 1L))
    },
    envir = globalenv()
  )

  .libPaths(wrong_library)
  loaded <- load_user_proxiome_data(rds_path, app_dir = APP_DIR)
  seurat_path <- find.package("Seurat", quiet = TRUE)

  expect_equal(loaded$source$source_type, "user_rds")
  expect_equal(captured_loader_args$marker_selection, "all")
  expect_true(startsWith(basename(captured_cache_path), "user-rds-"))
  expect_false(captured_force)
  expect_true(any(startsWith(
    normalizePath(captured_lib_paths, mustWork = FALSE),
    file.path(app_renv_root, "renv", "library")
  )))
  expect_length(seurat_path, 1)
  expect_false(startsWith(normalizePath(seurat_path, mustWork = TRUE), normalizePath(wrong_library, mustWork = TRUE)))
})

test_that("server RDS path loading is delegated to an async task", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")
  data_source_source <- paste(readLines(file.path(APP_DIR, "R", "data_source_module.R"), warn = FALSE), collapse = "\n")
  load_observer_start <- regexpr("observeEvent(input$load_rds_path", data_source_source, fixed = TRUE)[[1]]
  expect_gt(load_observer_start, 0)
  load_observer <- substr(data_source_source, load_observer_start, load_observer_start + 2600)

  expect_true(grepl("create_user_rds_load_task", data_source_source, fixed = TRUE))
  expect_true(grepl("ExtendedTask$new", data_source_source, fixed = TRUE))
  expect_true(grepl("promises::future_promise", data_source_source, fixed = TRUE))
  expect_false(grepl("future::future", data_source_source, fixed = TRUE))
  expect_true(grepl("future::plan(future::multisession", data_source_source, fixed = TRUE))
  expect_true(grepl("rds_load_message <- reactiveVal", data_source_source, fixed = TRUE))
  expect_true(grepl('rds_load_state <- reactiveVal("idle")', data_source_source, fixed = TRUE))
  expect_true(grepl("rds_load_progress_path <- reactiveVal", data_source_source, fixed = TRUE))
  expect_true(grepl("send_rds_load_state <- function", data_source_source, fixed = TRUE))
  expect_true(grepl("write_rds_load_progress", data_source_source, fixed = TRUE))
  expect_true(grepl("read_rds_load_progress", data_source_source, fixed = TRUE))
  expect_true(grepl("invalidateLater(1000", data_source_source, fixed = TRUE))
  expect_true(grepl('session$sendCustomMessage("proxiome-rds-load-state"', data_source_source, fixed = TRUE))
  expect_true(grepl("Shiny.addCustomMessageHandler('proxiome-rds-load-state'", app_source, fixed = TRUE))
  expect_true(grepl("rds_load_progress_bar", app_source, fixed = TRUE))
  expect_true(grepl("current_status <- isolate(user_rds_load_task$status())", load_observer, fixed = TRUE))
  expect_true(grepl('identical(current_status, "running")', load_observer, fixed = TRUE))
  expect_true(grepl('identical(current_state, "running")', load_observer, fixed = TRUE))
  expect_true(grepl("RDS load is already running", load_observer, fixed = TRUE))
  expect_false(grepl("validate_optional_pixelator_dir", load_observer, fixed = TRUE))
  expect_true(grepl("user_rds_load_task$invoke(input$rds_server_path, progress_path)", load_observer, fixed = TRUE))
  expect_false(grepl("load_user_proxiome_data(", load_observer, fixed = TRUE))
  expect_true(grepl("user_rds_load_task$result()", data_source_source, fixed = TRUE))
  expect_true(grepl('identical(rds_load_state(), "running")', data_source_source, fixed = TRUE))
})

test_that("ProxiomeVis writable directory defaults under user home", {
  home_dir <- tempfile("proxiomevis-home-")
  dir.create(home_dir, recursive = TRUE)
  withr::local_envvar(c(
    HOME = home_dir,
    PROXIOMEVIS_HOME = NA,
    PROXIOME_RUNTIME_DIR = NA
  ))

  writable_dir <- proxiomevis_home_dir()
  runtime_dir <- proxiomevis_runtime_dir()
  cache_dir <- proxiomevis_cache_dir()

  expect_equal(writable_dir, normalizePath(file.path(home_dir, ".ProxiomeVis"), mustWork = TRUE))
  expect_equal(runtime_dir, normalizePath(file.path(home_dir, ".ProxiomeVis", "runtime"), mustWork = TRUE))
  expect_equal(cache_dir, normalizePath(file.path(home_dir, ".ProxiomeVis", "cache"), mustWork = TRUE))
  expect_equal(Sys.getenv("PROXIOMEVIS_HOME", unset = ""), writable_dir)
  expect_equal(Sys.getenv("PROXIOME_RUNTIME_DIR", unset = ""), runtime_dir)
})

test_that("app resolves active R library paths without pixi libraries", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")
  resolved_libs <- resolve_app_r_libs(APP_DIR)
  renv_settings <- jsonlite::fromJSON(file.path(APP_DIR, "renv", "settings.json"))
  active_non_pixi_libs <- normalizePath(.libPaths(), mustWork = TRUE)
  active_non_pixi_libs <- active_non_pixi_libs[!grepl("/.pixi/", active_non_pixi_libs, fixed = TRUE)]

  expect_true(all(active_non_pixi_libs %in% resolved_libs))
  expect_equal(Sys.getenv("RENV_PATHS_LIBRARY", unset = ""), file.path(app_renv_root, "renv", "library"))
  expect_false(isTRUE(renv_settings$use.cache))
  expect_false(any(grepl("/.pixi/", resolved_libs, fixed = TRUE)))
  expect_false(grepl(".pixi/envs/r/lib/R/library", app_source, fixed = TRUE))
})

test_that("generic Open OnDemand launcher is redirected into a clean renv process", {
  app_source <- paste(readLines(file.path(APP_DIR, "app.R"), warn = FALSE), collapse = "\n")
  launch_expr <- build_clean_renv_launch_expression(APP_DIR, structure("/tmp/app.sock", mask = 63L), FALSE)

  expect_true(grepl("needs_clean_renv_relaunch", app_source, fixed = TRUE))
  expect_true(grepl("PROXIOME_REEXEC", app_source, fixed = TRUE))
  expect_true(grepl("RENV_CONFIG_CACHE_ENABLED = 'FALSE'", launch_expr, fixed = TRUE))
  expect_true(grepl("source(file.path(app_dir, 'renv', 'activate.R'))", launch_expr, fixed = TRUE))
  expect_true(grepl("shiny::runApp(app_dir)", launch_expr, fixed = TRUE))
  expect_true(grepl("complete_libraries", launch_expr, fixed = TRUE))
  expect_true(grepl(".libPaths(normalizePath(c(complete_libraries, .libPaths())", launch_expr, fixed = TRUE))
  expect_lt(
    regexpr(".libPaths(normalizePath(c(complete_libraries, .libPaths())", launch_expr, fixed = TRUE)[[1]],
    regexpr("shiny::runApp(app_dir)", launch_expr, fixed = TRUE)[[1]]
  )
  expect_true(grepl(".shiny_port <- structure(\"/tmp/app.sock\", mask = 63L)", launch_expr, fixed = TRUE))
})

test_that("data source summary uses source display names when present", {
  data <- list(
    source = list(
      rds_path = "/shared/proxiome/demo.rds",
      display_name = "patient_sample.rds",
      n_cells = 1234,
      n_markers = 56
    )
  )

  expect_equal(
    data_source_summary(data),
    "patient_sample.rds | 1,234 cells | 56 markers"
  )
  expect_equal(data_source_summary(NULL), "Loading demo RDS...")
})

test_that("readout modes are full-screen bslib card navigation", {
  html <- htmltools::renderTags(ui)$html

  expect_true(grepl("nav-underline", html, fixed = TRUE))
  expect_true(grepl("bslib-full-screen-enter", html, fixed = TRUE))
  expect_false(grepl("navset_card_tab", html, fixed = TRUE))
})

test_that("summary metrics use bslib value boxes", {
  html <- htmltools::renderTags(metric_row(metric_tile("Cells", "10")))$html

  expect_true(grepl("value-box", html, fixed = TRUE))
  expect_false(grepl('class="metric"', html, fixed = TRUE))
})

test_that("percent formatting supports QC hover vectors", {
  expect_equal(format_percent(c(1, 0.625, NA)), c("100%", "62.5%", NA))
})

test_that("filter cell count plot follows the notebook helper interface", {
  filter_cell_counts <- data.frame(
    step = rep(c("00_loaded", "01_after_n_umi_filter"), each = 3),
    sample = rep(c("s1", "s2", "TOTAL"), times = 2),
    condition = rep(c("UNT", "PHA", "TOTAL"), times = 2),
    n_cells = c(10, 20, 30, 6, 15, 21),
    stringsAsFactors = FALSE
  )

  p <- plot_filter_cell_counts(filter_cell_counts, include_total = FALSE, y = "fraction_loaded")

  expect_s3_class(p, "ggplot")
  expect_equal(unique(as.character(p$data$sample)), c("s1", "s2"))
  expect_equal(p$data$value, c(1, 1, 0.6, 0.75))
  expect_equal(p$labels$y, "Fraction of loaded cells")
})

test_that("cell calling rank plot emits visible Plotly traces", {
  metadata <- data.frame(
    component = paste0("cell-", 1:6),
    sample = rep(c("s1", "s2"), each = 3),
    n_umi = c(120, 90, 60, 200, 160, 80),
    stringsAsFactors = FALSE
  )

  widget <- qc_molecule_rank_plotly(metadata, cutoff = 100)

  expect_s3_class(widget, "plotly")
  expect_true(length(widget$x$data) >= 3)
  expect_true(all(vapply(widget$x$data[1:2], function(trace) length(trace$x) > 0, logical(1))))
  expect_true(any(vapply(widget$x$data, function(trace) identical(trace$mode, "lines+markers"), logical(1))))
})

test_that("QC distribution plot emits visible Plotly traces", {
  metadata <- data.frame(
    component = paste0("cell-", 1:6),
    sample = rep(c("s1", "s2"), each = 3),
    n_umi = c(120, 90, 60, 200, 160, 80),
    stringsAsFactors = FALSE
  )

  widget <- qc_distribution_plotly(metadata, metric = "n_umi")

  expect_s3_class(widget, "plotly")
  expect_equal(length(widget$x$data), 2L)
  expect_true(all(vapply(widget$x$data, function(trace) identical(trace$type, "violin"), logical(1))))
  expect_true(all(vapply(widget$x$data, function(trace) length(trace$y) > 0, logical(1))))
})

test_that("colocalization heatmap helper follows the notebook shared-order contract", {
  coloc_summary <- expand.grid(
    condition = c("CD3CD28", "UNT"),
    marker_1 = c("CD3e", "CD4", "CD8"),
    marker_2 = c("CD3e", "CD4", "CD8"),
    stringsAsFactors = FALSE
  )
  coloc_summary <- coloc_summary[coloc_summary$marker_1 != coloc_summary$marker_2, , drop = FALSE]
  coloc_summary$mean_log2_ratio <- seq(-0.6, 0.6, length.out = nrow(coloc_summary))
  coloc_summary$pct_detected <- seq(0.2, 0.9, length.out = nrow(coloc_summary))

  result <- make_coloc_heatmaps(
    data = coloc_summary,
    selected_markers = c("CD3e", "CD4", "CD8"),
    cell_label = "CD8 T",
    conditions = c("CD3CD28", "UNT"),
    reference_condition = "CD3CD28",
    clustering_method = "ward.D2",
    legend_range = c(-1, 1)
  )

  expect_true(all(c("marker_order", "plots", "plot_data", "plot") %in% names(result)))
  expect_setequal(result$marker_order, c("CD3e", "CD4", "CD8"))
  expect_named(result$plots, c("CD3CD28", "UNT"))
  expect_s3_class(result$plot, "ggplot")
  expect_equal(result$plot$coordinates$ratio, 1)
  expect_true(inherits(result$plot$theme$panel.border, "element_rect"))
  expect_true(all(c("marker_1", "marker_2", "mean_log2_ratio", "pct_detected", "hover") %in% names(result$plot_data)))
})

test_that("colocalization heatmap mirrors one-direction marker pairs", {
  coloc_summary <- data.frame(
    condition = "CD3CD28",
    marker_1 = c("CD3e", "CD3e", "CD4"),
    marker_2 = c("CD4", "CD8", "CD8"),
    mean_log2_ratio = c(0.5, -0.2, 0.7),
    pct_detected = c(0.6, 0.8, 0.4),
    n_detected = c(6L, 8L, 4L),
    n_total = 10L,
    stringsAsFactors = FALSE
  )

  result <- make_coloc_heatmaps(
    data = coloc_summary,
    selected_markers = c("CD3e", "CD4", "CD8"),
    cell_label = "CD8 T",
    conditions = "CD3CD28",
    reference_condition = "CD3CD28",
    clustering_method = "ward.D2",
    legend_range = c(-1, 1)
  )

  find_pair <- function(marker_1, marker_2) {
    result$plot_data[
      as.character(result$plot_data$marker_1) == marker_1 &
        as.character(result$plot_data$marker_2) == marker_2,
      ,
      drop = FALSE
    ]
  }

  forward <- find_pair("CD3e", "CD4")
  reverse <- find_pair("CD4", "CD3e")

  expect_equal(nrow(forward), 1L)
  expect_equal(nrow(reverse), 1L)
  expect_equal(reverse$plot_value, forward$plot_value)
  expect_equal(reverse$plot_size, forward$plot_size)
  expect_equal(reverse$n_detected, forward$n_detected)
  expect_equal(reverse$n_total, forward$n_total)
  expect_false(any(is.na(result$plot_data$plot_value)))
})

test_that("colocalization heatmap Plotly output displays a pct_detected size legend", {
  coloc_summary <- expand.grid(
    condition = c("CD3CD28", "UNT"),
    marker_1 = c("CD3e", "CD4", "CD8"),
    marker_2 = c("CD3e", "CD4", "CD8"),
    stringsAsFactors = FALSE
  )
  coloc_summary <- coloc_summary[coloc_summary$marker_1 != coloc_summary$marker_2, , drop = FALSE]
  coloc_summary$mean_log2_ratio <- seq(-0.6, 0.6, length.out = nrow(coloc_summary))
  coloc_summary$pct_detected <- seq(0.2, 0.9, length.out = nrow(coloc_summary))

  result <- make_coloc_heatmaps(
    data = coloc_summary,
    selected_markers = c("CD3e", "CD4", "CD8"),
    cell_label = "CD8 T",
    conditions = c("CD3CD28", "UNT"),
    reference_condition = "CD3CD28"
  )
  widget <- coloc_heatmap_plotly(result)

  expect_s3_class(widget, "plotly")
  expect_true(isTRUE(widget$x$layout$showlegend))
  expect_equal(widget$x$layout$legend$title$text, "pct_detected")
  expect_true(any(grepl("pct_detected", vapply(widget$x$data, function(trace) trace$name %||% "", character(1)), fixed = TRUE)))
})

test_that("colocalization observed heatmap has interactive and original R plot renderers", {
  html <- htmltools::renderTags(ui)$html
  colocalization_module_source <- paste(readLines(file.path(APP_DIR, "R", "colocalization_module.R"), warn = FALSE), collapse = "\n")
  result_start <- regexpr("colocalization_heatmap_result <- reactive", colocalization_module_source, fixed = TRUE)[[1]]
  interactive_start <- regexpr("output$colocalization_heatmap_interactive <- renderPlotly", colocalization_module_source, fixed = TRUE)[[1]]
  original_start <- regexpr("output$colocalization_heatmap_original <- renderPlot", colocalization_module_source, fixed = TRUE)[[1]]

  expect_true(grepl('id="colocalization-colocalization_heatmap_display"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_heatmap_interactive"', html, fixed = TRUE))
  expect_true(grepl('id="colocalization-colocalization_heatmap_original"', html, fixed = TRUE))
  expect_gt(result_start, 0)
  expect_gt(interactive_start, 0)
  expect_gt(original_start, 0)
  expect_lt(result_start, interactive_start)
  expect_lt(result_start, original_start)
  expect_true(grepl("colocalization_heatmap_result()", colocalization_module_source, fixed = TRUE))
  expect_true(grepl("coloc_heatmap_plotly(colocalization_heatmap_result())", colocalization_module_source, fixed = TRUE))
  expect_true(grepl("print(colocalization_heatmap_result()$plot)", colocalization_module_source, fixed = TRUE))
})

test_that("colocalization heatmap Plotly output preserves square heatmap panels", {
  coloc_summary <- expand.grid(
    condition = "CD3CD28",
    marker_1 = c("CD3e", "CD4", "CD8"),
    marker_2 = c("CD3e", "CD4", "CD8"),
    stringsAsFactors = FALSE
  )
  coloc_summary <- coloc_summary[coloc_summary$marker_1 != coloc_summary$marker_2, , drop = FALSE]
  coloc_summary$mean_log2_ratio <- seq(-0.6, 0.6, length.out = nrow(coloc_summary))
  coloc_summary$pct_detected <- seq(0.2, 0.9, length.out = nrow(coloc_summary))

  result <- make_coloc_heatmaps(
    data = coloc_summary,
    selected_markers = c("CD3e", "CD4", "CD8"),
    cell_label = "CD8 T",
    conditions = "CD3CD28",
    reference_condition = "CD3CD28"
  )
  widget <- coloc_heatmap_plotly(result)
  dimensions <- coloc_heatmap_widget_dimensions(result$plot_data)

  expect_false(isTRUE(widget$x$layout$autosize))
  expect_equal(widget$x$layout$width, dimensions$width)
  expect_equal(widget$x$layout$height, dimensions$height)
  expect_equal(widget$x$layout$margin, dimensions$margin)
  expect_equal(widget$x$layout$yaxis$scaleanchor, "x")
  expect_equal(widget$x$layout$yaxis$scaleratio, 1)
})

test_that("colocalization heatmap caps marker size for dense marker panels", {
  markers <- paste0("M", seq_len(15))
  coloc_summary <- expand.grid(
    condition = c("CD3CD28", "PHA", "UNT"),
    marker_1 = markers,
    marker_2 = markers,
    stringsAsFactors = FALSE
  )
  coloc_summary <- coloc_summary[coloc_summary$marker_1 != coloc_summary$marker_2, , drop = FALSE]
  coloc_summary$mean_log2_ratio <- 0.4
  coloc_summary$pct_detected <- 1

  result <- make_coloc_heatmaps(
    data = coloc_summary,
    selected_markers = markers,
    cell_label = "selected cells",
    conditions = c("CD3CD28", "PHA", "UNT"),
    reference_condition = "CD3CD28"
  )
  widget <- coloc_heatmap_plotly(result)
  heatmap_sizes <- unlist(lapply(widget$x$data, function(trace) {
    if (identical(trace$visible, "legendonly")) {
      return(numeric(0))
    }
    if ((trace$mode %||% "") != "markers" || grepl("^pct_detected", trace$name %||% "")) {
      return(numeric(0))
    }
    sizes <- trace$marker$size %||% numeric(0)
    if (length(sizes) == 0) {
      return(numeric(0))
    }
    sizes
  }))

  expect_gt(length(heatmap_sizes), 0)
  expect_true(max(heatmap_sizes, na.rm = TRUE) <= 14)
})

test_that("differential plot dimensions scale with labels and facets", {
  short_volcano <- differential_volcano_dimensions("Effect")
  long_volcano <- differential_volcano_dimensions("Abundance effect: CD3CD28 minus UNT (reference)")
  small_detail <- differential_detail_dimensions(
    data.frame(
      condition = rep(c("UNT", "CD3CD28"), each = 4),
      celltype_manual = "CD8 T",
      stringsAsFactors = FALSE
    ),
    stratify_by_celltype = FALSE,
    y_label = "CD25 abundance"
  )
  faceted_detail <- differential_detail_dimensions(
    expand.grid(
      condition = c("UNT", "CD3CD28"),
      celltype_manual = paste0("Cell type ", seq_len(6)),
      stringsAsFactors = FALSE
    ),
    stratify_by_celltype = TRUE,
    y_label = "CD25 abundance"
  )

  expect_gt(long_volcano$width, short_volcano$width)
  expect_true(long_volcano$margin$b >= 130)
  expect_true(long_volcano$margin$r >= 150)
  expect_equal(small_detail$facet_rows, 1)
  expect_gt(faceted_detail$height, small_detail$height)
  expect_gt(faceted_detail$width, small_detail$width)
})

test_that("differential volcano plots use content-aware Plotly dimensions and margins", {
  result <- data.frame(
    marker_pair = c("CD3e / CD4", "CD3e / CD8", "CD4 / CD8"),
    celltype_manual = "combined",
    effect_size = c(-0.4, 0.15, 0.6),
    p_adj = c(0.01, 0.8, 0.03),
    group_a = "PHA",
    group_b = "UNT",
    mean_a = c(0.3, 0.2, 0.7),
    mean_b = c(0.7, 0.05, 0.1),
    direction = c("Higher in B", "Not significant", "Higher in A"),
    stringsAsFactors = FALSE
  )
  dimensions <- differential_volcano_dimensions("Difference in medians: PHA minus UNT (reference)")

  widget <- differential_volcano_plot(
    result,
    label_col = "marker_pair",
    x_label = "Difference in medians: PHA minus UNT (reference)",
    fdr_cutoff = 0.05,
    effect_cutoff = 0.25,
    source = "colocalization_diff",
    dimensions = dimensions
  )

  expect_equal(widget$x$layout$width, dimensions$width)
  expect_equal(widget$x$layout$height, dimensions$height)
  expect_equal(widget$x$layout$margin, dimensions$margin)
  expect_true(isTRUE(widget$x$layout$xaxis$automargin))
  expect_true(isTRUE(widget$x$layout$yaxis$automargin))
  expect_equal(widget$x$layout$xaxis$domain, c(0, 1))
  expect_equal(widget$x$layout$yaxis$domain, c(0, 1))
})

test_that("Open OnDemand template launches the app with the shared renv environment", {
  script_path <- file.path(APP_DIR, "template", "script.sh.erb")

  expect_true(file.exists(script_path))

  script <- paste(readLines(script_path), collapse = "\n")
  expect_true(grepl("renv/activate.R", script, fixed = TRUE))
  expect_true(grepl('HOSTNAME_FQDN="${HOSTNAME:-$(hostname -f 2>/dev/null || hostname 2>/dev/null || echo unknown)}"', script, fixed = TRUE))
  expect_true(grepl('ncifcrf.gov|*.ncifcrf.gov)', script, fixed = TRUE))
  expect_true(grepl('biowulf.nih.gov|*.biowulf.nih.gov)', script, fixed = TRUE))
  expect_true(grepl('DEFAULT_APP_DIR="/mnt/ccrsf-static/illumina/RnD_pixelgen_CAR-T_datasets/shiny/proxiome_demo"', script, fixed = TRUE))
  expect_true(grepl('APP_DIR="${PROXIOME_APP_DIR:-${DEFAULT_APP_DIR}}"', script, fixed = TRUE))
  expect_true(grepl('PROXIOME_R_MODULE="${PROXIOME_R_MODULE:-${DEFAULT_R_MODULE}}"', script, fixed = TRUE))
  expect_true(grepl('if [[ -n "${PROXIOME_R_MODULE:-}" ]]; then', script, fixed = TRUE))
  expect_true(grepl('module load "${PROXIOME_R_MODULE}"', script, fixed = TRUE))
  expect_true(grepl("RENV_CONFIG_AUTOLOADER_ENABLED", script, fixed = TRUE))
  expect_true(grepl("RENV_CONFIG_CACHE_ENABLED", script, fixed = TRUE))
  expect_true(grepl("RENV_DEPLOYMENT_STRICT", script, fixed = TRUE))
  expect_true(grepl("RENV_PATHS_LIBRARY", script, fixed = TRUE))
  expect_true(grepl("PROXIOME_DEMO_RDS", script, fixed = TRUE))
  expect_true(grepl('PROXIOMEVIS_HOME="${PROXIOMEVIS_HOME:-${HOME}/.ProxiomeVis}"', script, fixed = TRUE))
  expect_true(grepl('PROXIOME_RUNTIME_DIR="${PROXIOME_RUNTIME_DIR:-${PROXIOMEVIS_HOME}/runtime}"', script, fixed = TRUE))
  expect_true(grepl('mkdir -p "${PROXIOMEVIS_HOME}" "${PROXIOMEVIS_HOME}/cache" "${PROXIOME_RUNTIME_DIR}"', script, fixed = TRUE))
  expect_true(grepl('RENV_PROJECT="${APP_DIR}"', script, fixed = TRUE))
  expect_true(grepl('${APP_DIR}/renv/library', script, fixed = TRUE))
  expect_true(grepl('R_BIN="${R_BIN:-$(command -v Rscript)}"', script, fixed = TRUE))
  expect_true(grepl("${port:-${PORT:-${PORT1:-6323}}}", script, fixed = TRUE))
  expect_true(grepl("shiny::runApp", script, fixed = TRUE))
  expect_true(grepl("complete_libraries", script, fixed = TRUE))
  expect_true(grepl(".libPaths(normalizePath(c(complete_libraries, .libPaths())", script, fixed = TRUE))
  expect_lt(
    regexpr(".libPaths(normalizePath(c(complete_libraries, .libPaths())", script, fixed = TRUE)[[1]],
    regexpr("shiny::runApp", script, fixed = TRUE)[[1]]
  )
  expect_true(grepl('SHINY_HOST="${SHINY_HOST:-0.0.0.0}"', script, fixed = TRUE))
  expect_true(grepl("for required_package in renv shiny bslib ggplot2 plotly future promises rlang data.table", script, fixed = TRUE))
  expect_false(grepl(".pixi/envs/r/bin/Rscript", script, fixed = TRUE))
  expect_false(grepl(".pixi/envs/r/lib/R/library", script, fixed = TRUE))
  expect_false(grepl('APP_ROOT="/mnt/ccrsf-static/illumina/RnD_pixelgen_CAR-T_datasets"', script, fixed = TRUE))
  expect_false(grepl("renv::restore", script, fixed = TRUE))
  expect_false(grepl("install.packages", script, fixed = TRUE))
})

test_that("sidebars organize controls with accordions", {
  html <- htmltools::renderTags(ui)$html

  expect_true(grepl("accordion", html, fixed = TRUE))
  expect_true(grepl("Display", html, fixed = TRUE))
  expect_true(grepl("Filters", html, fixed = TRUE))
  expect_true(grepl("Contrast", html, fixed = TRUE))
  expect_true(grepl("Thresholds", html, fixed = TRUE))
  expect_true(grepl("Detail", html, fixed = TRUE))
})

test_that("abundance UMAP can be split by condition or sample", {
  choices <- available_abundance_split_choices(
    data.frame(
      condition = c("UNT", "PHA"),
      sample = c("s1", "s2"),
      celltype_manual = c("B", "T")
    )
  )

  expect_named(choices, c("None", "Condition", "Sample"))
  expect_equal(unname(choices), c("", "condition", "sample"))
})

test_that("abundance UMAP can be colored by marker, cell type, condition, or sample", {
  choices <- available_abundance_color_choices(
    data.frame(
      condition = c("UNT", "PHA"),
      sample = c("s1", "s2"),
      celltype_manual = c("B", "T")
    )
  )

  expect_named(choices, c("Marker abundance", "Cell type", "Condition", "Sample"))
  expect_equal(unname(choices), c("abundance", "celltype_manual", "condition", "sample"))
})

test_that("abundance marker distribution plot keeps violin groups valid for Plotly", {
  plot_data <- data.frame(
    component = paste0("cell-", seq_len(12)),
    marker = "CD3e",
    abundance = c(0.2, 0.4, 0.7, 1.0, 0.3, 0.6, 0.8, 1.2, 0.5, 0.9, 1.1, 1.4),
    sample_alias = rep(c("S1", "S2"), each = 6),
    condition = rep(c("UNT", "PHA"), each = 6),
    celltype_manual = rep(c("CD4 T", "CD8 T"), times = 6),
    stringsAsFactors = FALSE
  )

  plot <- plot_abundance_marker_distribution(plot_data, "CD3e", facet_cols = 1, show_jitter = TRUE)
  plot_without_jitter <- plot_abundance_marker_distribution(plot_data, "CD3e", facet_cols = 1, show_jitter = FALSE)
  built <- ggplot2::ggplot_build(plot)
  built_without_jitter <- ggplot2::ggplot_build(plot_without_jitter)

  expect_gt(nrow(built$data[[1]]), 0)
  expect_equal(max(built$layout$layout$COL), 1)
  expect_length(built$data, 3)
  expect_length(built_without_jitter$data, 2)
  expect_s3_class(plotly::ggplotly(plot, tooltip = c("x", "y", "fill")), "plotly")
})

test_that("abundance marker distribution controls include visible defaults", {
  abundance_module_source <- paste(readLines(file.path(APP_DIR, "R", "abundance_module.R"), warn = FALSE), collapse = "\n")

  expect_true(grepl('numericInput(ns("abundance_distribution_columns"), "Facet columns", value = 3', abundance_module_source, fixed = TRUE))
  expect_true(grepl('numericInput(ns("abundance_distribution_width"), "Plot width (px)", value = 832', abundance_module_source, fixed = TRUE))
  expect_true(grepl('numericInput(ns("abundance_distribution_height"), "Plot height (px)", value = 678', abundance_module_source, fixed = TRUE))
  expect_true(grepl("update_abundance_distribution_size_controls", abundance_module_source, fixed = TRUE))
  expect_true(grepl('checkboxInput(ns("abundance_distribution_show_jitter"), "Show jitter dots", value = TRUE)', abundance_module_source, fixed = TRUE))
})

test_that("abundance marker distribution dimensions scale and accept user overrides", {
  small_data <- data.frame(
    sample_alias = rep(c("S1", "S2"), each = 4),
    celltype_manual = rep(c("CD4 T", "CD8 T"), times = 4),
    stringsAsFactors = FALSE
  )
  large_data <- expand.grid(
    sample_alias = paste0("S", seq_len(5)),
    celltype_manual = paste0("Cell type ", seq_len(7)),
    stringsAsFactors = FALSE
  )

  small_dimensions <- abundance_distribution_widget_dimensions(small_data)
  large_dimensions <- abundance_distribution_widget_dimensions(large_data)
  override_dimensions <- abundance_distribution_widget_dimensions(large_data, width_px = 900, height_px = 720)
  single_column_dimensions <- abundance_distribution_widget_dimensions(large_data, facet_cols = 1)
  four_column_dimensions <- abundance_distribution_widget_dimensions(large_data, facet_cols = 4)

  expect_gt(large_dimensions$width, small_dimensions$width)
  expect_gt(large_dimensions$height, small_dimensions$height)
  expect_equal(override_dimensions$width, 900)
  expect_equal(override_dimensions$height, 720)
  expect_equal(single_column_dimensions$facet_cols, 1)
  expect_equal(single_column_dimensions$facet_rows, 7)
  expect_equal(four_column_dimensions$facet_cols, 4)
  expect_equal(four_column_dimensions$facet_rows, 2)
  expect_gt(four_column_dimensions$width, single_column_dimensions$width)
  expect_lt(four_column_dimensions$height, single_column_dimensions$height)
})

test_that("clustering per-marker distribution plot keeps violin groups valid for Plotly", {
  plot_data <- data.frame(
    component = paste0("cell-", seq_len(12)),
    marker = "CD3e",
    log2_ratio = c(-0.4, -0.1, 0.2, 0.5, -0.2, 0.1, 0.4, 0.8, -0.3, 0.0, 0.3, 0.7),
    sample_alias = rep(c("S1", "S2"), each = 6),
    condition = rep(c("UNT", "PHA"), each = 6),
    celltype_manual = rep(c("CD4 T", "CD8 T"), times = 6),
    stringsAsFactors = FALSE
  )

  plot <- plot_clustering_per_marker(plot_data, "CD3e")
  built <- ggplot2::ggplot_build(plot)

  expect_gt(nrow(built$data[[2]]), 0)
  expect_s3_class(plotly::ggplotly(plot, tooltip = c("x", "y", "fill")), "plotly")
})
