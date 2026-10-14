precomputed_ui <- function(id, dataset_control = NULL, dataset_status = NULL) {
  ns <- shiny::NS(id)
  bslib::navset_bar(
    id = ns("report_page"), title = "ProxiomeVis", fillable = FALSE,
    navbar_options = bslib::navbar_options(bg = "#f8f9fa", theme = "light"),
    header = shiny::tagList(dataset_status,
      shiny::tags$head(shiny::tags$script(src = "precomputed.js"),
        shiny::tags$style(shiny::HTML(
          ".spatial-plot:has(.js-plotly-plot):not(:has(.recalculating)) .spatial-loading, .spatial-plot:has(.shiny-output-error) .spatial-loading {display:none;}
           .spatial-matrix {width:100%;max-width:900px;}
           .spatial-matrix .plotly {aspect-ratio:6 / 5;min-height:400px;}")))),
    bslib::nav_panel(
      "Overview / QC",
      shiny::uiOutput(ns("dataset_heading")),
      bslib::layout_column_wrap(
        width = "180px", fill = FALSE,
        precomputed_metric_card("Cells", ns("metric_cells")),
        precomputed_metric_card("Samples", ns("metric_samples")),
        precomputed_metric_card("Donors", ns("metric_donors")),
        precomputed_metric_card("Cell types", ns("metric_celltypes"))
      ),
      bslib::layout_columns(
        col_widths = c(7, 5), fill = FALSE,
        bslib::card(full_screen = TRUE,
          bslib::card_header("QC cell counts"),
          plotly::plotlyOutput(ns("qc_plot"), height = "480px")
        ),
        bslib::card(full_screen = TRUE,
          bslib::card_header("QC table"), DT::DTOutput(ns("qc_table"))
        )
      ),
      bslib::card(
        class = "mt-3", bslib::card_header("Analysis provenance"),
        precomputed_download_group(
          shiny::downloadButton(ns("download_provenance"), "Provenance JSON", class = "btn-sm btn-outline-primary")
        ),
        shiny::tags$pre(class = "small", shiny::textOutput(ns("provenance_text")))
      )
    ),
    bslib::nav_panel(
      "Embeddings",
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          title = "Embedding controls", width = 300,
          shiny::selectInput(ns("embedding_reduction"), "Embedding", choices = character()),
          shiny::selectInput(ns("embedding_split"), "Split by", choices = c(
            "None" = "none", "Donor" = "donor", "Sample alias" = "sample_id")),
          shiny::selectInput(ns("embedding_color"), "Color cells by", choices = c(
            "Cell type" = "celltype", "Cluster" = "cluster", "Condition" = "condition",
            "Donor" = "donor", "Marker abundance" = "abundance"), selected = "celltype"),
          shiny::conditionalPanel("input.embedding_color == 'abundance'", ns = ns,
            shiny::selectizeInput(ns("embedding_marker"), "Abundance marker", choices = character())
          ),
          shiny::tags$p(class = "text-muted small",
            "Embedding colors are annotations or normalized abundance from the loaded RDS; no spatial file is read.")
        ),
          bslib::card(full_screen = TRUE,
            bslib::card_header("Embedding"),
            precomputed_download_group(
              shiny::downloadButton(ns("download_embedding_svg"), "SVG", class = "btn-sm btn-outline-primary")
            ),
            plotly::plotlyOutput(ns("embedding_plot"), height = "620px")
          )
      )
    ),
    bslib::nav_panel(
      "Precomputed results",
      bslib::layout_column_wrap(
        width = "220px", fill = FALSE,
        shiny::uiOutput(ns("donor_control")),
        shiny::selectInput(ns("celltype"), "Cell type", choices = character()),
        shiny::selectInput(ns("comparison"), "Comparison (target − reference)", choices = character()),
        shiny::selectInput(ns("readout"), "Readout", choices = PRECOMPUTED_READOUTS)
      ),
      shiny::uiOutput(ns("result_summary")),
      bslib::navset_underline(
        id = ns("results_view"),
        bslib::nav_panel(
          "Tables",
          shiny::uiOutput(ns("fdr_controls")),
          shiny::tags$p(class = "small text-muted",
            "Table filters and searches affect only tables and filtered CSV downloads. All stored rows exports ignore those filters. Click a result row to select its feature."),
          bslib::navset_card_tab(
            full_screen = TRUE, id = ns("result_tables"),
            bslib::nav_panel(
              "Results",
              precomputed_download_group(
                shiny::downloadButton(ns("download_filtered_results_csv"), "Filtered table CSV", class = "btn-sm btn-outline-primary"),
                shiny::downloadButton(ns("download_results_csv"), "All stored rows CSV", class = "btn-sm btn-outline-primary")
              ),
              DT::DTOutput(ns("results_table"))
            ),
            bslib::nav_panel(
              "Sample support",
              precomputed_download_group(
                shiny::downloadButton(ns("download_samples_csv"), "Sample support CSV", class = "btn-sm btn-outline-primary")
              ),
              DT::DTOutput(ns("samples_table"))
            )
          ),
          bslib::accordion(bslib::accordion_panel("Statistical methods and original thresholds",
            shiny::uiOutput(ns("methods_ui"))), open = FALSE, class = "mt-3")
        ),
        bslib::nav_panel(
          "Volcano",
          shiny::tags$div(class = "alert alert-warning small mt-2",
            shiny::tags$strong("Exploratory component-level test. "),
            "Stored P values describe cells within the selected donor, not a population-level treatment effect. Table filters do not change the original report thresholds or plot."
          ),
          bslib::layout_columns(
            col_widths = c(6, 6), fill = FALSE,
            bslib::card(full_screen = TRUE,
              bslib::card_header("Stored report volcano"),
              precomputed_download_group(
                shiny::downloadButton(ns("download_report_volcano_svg"), "Original report SVG", class = "btn-sm btn-outline-primary"),
                shiny::downloadButton(ns("download_report_volcano_csv"), "All stored volcano rows CSV", class = "btn-sm btn-outline-primary")
              ),
              shiny::tags$p(class = "small text-muted",
                "Click a point to update the cell distribution alongside. Hover for details; scroll or drag a box to zoom; double-click to reset."),
              plotly::plotlyOutput(ns("report_volcano_plot"), height = "650px")
            ),
            bslib::card(full_screen = TRUE,
              bslib::card_header("Selected feature — cell distribution"),
              shiny::tags$p(class = "small text-muted",
                "Violin and boxplot show available cell values, not a new test. Plotted cells can differ from the original test cohort. Spatial distributions require the proximity sidecar."),
              plotly::plotlyOutput(ns("volcano_distribution"), height = "650px")
            )
          )
        ),
        bslib::nav_panel(
          "Spatial patterns",
          shiny::uiOutput(ns("spatial_view_control")),
          shiny::conditionalPanel("!input.spatial_view || input.spatial_view == 'matrix'", ns = ns,
          bslib::layout_column_wrap(
            width = "280px", fill = FALSE,
            shiny::radioButtons(ns("spatial_kind"), "Spatial quantity", choices = c(
              "Differential effect (difference of medians)" = "differences",
              "Condition profile (mean proximity)" = "profiles"), selected = "differences"),
            shiny::conditionalPanel("input.spatial_kind == 'profiles'", ns = ns,
              shiny::selectInput(ns("spatial_panel"), "Condition", choices = character())
            )
          ),
          shiny::tags$p(class = "small text-muted",
            "Available for spatial readouts; combines self-proximity and colocalization. A mean condition profile is not the quantity used for a difference-of-medians effect. No proximity is recalculated."),
          bslib::card(full_screen = TRUE,
            bslib::card_header("Aligned spatial matrix"),
            precomputed_download_group(
              shiny::downloadButton(ns("download_spatial_report_svg"), "Original report SVG", class = "btn-sm btn-outline-primary"),
              shiny::downloadButton(ns("download_spatial_csv"), "Displayed spatial rows CSV", class = "btn-sm btn-outline-primary")
            ),
            shiny::tags$p(class = "small text-muted",
              "Hover for values and availability; click a marker pair to select it for inspection in Tables or Volcano. Marker order and scales follow the frozen report alignment."),
            shiny::div(class = "spatial-plot spatial-matrix",
              shiny::div(id = ns("spatial_loading"), class = "spatial-loading alert alert-info",
                role = "status", `aria-live` = "polite",
                shiny::span(class = "spinner-border spinner-border-sm me-2", `aria-hidden` = "true"),
                "Loading spatial matrix…"),
              plotly::plotlyOutput(ns("spatial_report_plot"), height = "auto")
            )
          )),
          shiny::conditionalPanel("input.spatial_view == 'network' && input.readout == 'colocalization'", ns = ns,
            bslib::layout_column_wrap(width = "280px", fill = FALSE,
              shiny::selectizeInput(ns("network_markers"), "Network markers (2–20)", choices = character(),
                multiple = TRUE, options = list(maxItems = 20)),
              shiny::numericInput(ns("network_detection"), "Minimum detection fraction (greater than)",
                value = 0.5, min = 0, max = 1, step = 0.05),
              shiny::numericInput(ns("network_mean"), "Minimum absolute mean proximity (greater than)",
                value = 0.1, min = 0, step = 0.05)
            ),
            shiny::tags$p(class = "small text-muted",
              "Each condition is filtered separately. The initial marker set uses up to 10 markers from the strongest stored pairs. Color shows mean log2 ratio; line width shows its absolute magnitude. Positions and scales are shared across conditions. These display thresholds are not significance tests. Stored population means include missing scores as zero."),
            bslib::card(full_screen = TRUE,
              bslib::card_header("Colocalization network — mean proximity by condition"),
              precomputed_download_group(
                shiny::downloadButton(ns("download_network_svg"), "Network SVG", class = "btn-sm btn-outline-primary"),
                shiny::downloadButton(ns("download_network_csv"), "Displayed edges CSV", class = "btn-sm btn-outline-primary")
              ),
              shiny::tags$p(class = "small text-muted",
                "Hover over an edge for its score and support. Click an edge to select its pair for Tables or Volcano. No proximity sidecar is required."),
              shiny::textOutput(ns("network_selected_pair")),
              shiny::div(class = "spatial-plot",
                shiny::div(class = "spatial-loading alert alert-info", role = "status", `aria-live` = "polite",
                  shiny::span(class = "spinner-border spinner-border-sm me-2", `aria-hidden` = "true"),
                  "Loading colocalization network…"),
                plotly::plotlyOutput(ns("network_plot"), height = "650px")
              )
            )
          )
        )
      ),
      shiny::conditionalPanel("input.results_view != 'Spatial patterns'", ns = ns,
      bslib::card(full_screen = TRUE,
        class = "mt-3", bslib::card_header("Selected feature — shared across views"),
        bslib::layout_columns(
          col_widths = c(4, 8), fill = FALSE,
          shiny::tagList(
            shiny::selectizeInput(ns("detail_feature"), "Marker / pair", choices = character()),
            shiny::conditionalPanel("input.results_view == 'Spatial patterns'", ns = ns,
              shiny::radioButtons(ns("detail_display"), "Detail display", choices = c(
                "Stored effects" = "effects", "Condition values, when stored" = "conditions"), selected = "effects")
            ),
            shiny::uiOutput(ns("feature_context")),
            shiny::actionButton(ns("explore_3d"), "Explore cells in 3D", class = "btn-outline-primary"),
            shiny::tags$p(class = "small text-muted mt-2", "Transfers the selected scope and pair. Choose a condition and one highlight marker, then click Display to load a graph."),
            shiny::downloadButton(ns("download_detail_svg"), "Feature detail SVG", class = "btn-sm btn-outline-primary")
          ),
          shiny::conditionalPanel("input.results_view != 'Volcano'", ns = ns,
            plotly::plotlyOutput(ns("effect_detail"), height = "390px"))
        )
      ))
    ),
    bslib::nav_panel("Custom analysis", hybrid_analysis_ui(ns("analysis"))),
    bslib::nav_panel("Cell 3D", cell_graph_ui(ns("graph"))),
    bslib::nav_panel("Documentation", documentation_ui()),
    bslib::nav_spacer(),
    bslib::nav_item(dataset_control)
  )
}

documentation_ui <- function() {
  tags <- shiny::tags
  tags$article(style = "max-width: 75ch; margin: 1.5rem auto;",
    tags$h1("ProxiomeVis documentation"),
    tags$p("ProxiomeVis shows prepared Pixelgen results. It also calculates custom comparisons and displays 3D cell graphs."),
    tags$p(tags$a(href = "https://github.com/CCRSF-IFX/ProxiomeVis",
      target = "_blank", rel = "noopener noreferrer", "GitHub repository (opens in a new tab)"),
      ". The repository contains the source code, preparation instructions, and installation requirements."),
    bslib::accordion(open = "Start with a dataset",
      bslib::accordion_panel("Start with a dataset",
        tags$p("The administrator supplies an approved RDS file. The public application does not accept file uploads."),
        tags$ol(tags$li("Click Dataset at the top right."),
          tags$li("Check the dataset name and available files.")),
        tags$p("A sidecar is a separate file that supplies additional data. The RDS file identifies the sidecar locations."),
        tags$ul(tags$li("Prepared tables, embeddings, and spatial summaries need only the RDS file."),
          tags$li("Cell-level spatial distributions, cell ranking, and custom spatial tests need the proximity sidecar."),
          tags$li("The Cell 3D view needs the graph sidecar. PXL files are not necessary during application use.")),
        tags$p("For local use, the administrator can enable server-path entry. Use RDS files only from approved sources.")
      ),
      bslib::accordion_panel("Overview / QC and Embeddings",
        tags$p("QC means quality control. Overview / QC shows cell counts before and after each filter."),
        tags$p("The QC table identifies the sample, donor, and condition. Provenance records the data sources and analysis settings."),
        tags$ol(tags$li("Open Embeddings."), tags$li("Select an Embedding."),
          tags$li("Select Color cells by."), tags$li("For Marker abundance, select an Abundance marker."),
          tags$li("Use Split by to make separate panels for each donor or sample alias.")),
        tags$p("An embedding places cells in a reduced coordinate space. Split panels use the same axes and color scale.")
      ),
      bslib::accordion_panel("Precomputed results",
        tags$p("Precomputed results come from an analysis completed before application use. Display controls do not calculate new tests."),
        tags$ol(tags$li("Open Precomputed results."), tags$li("If more than one donor is available, select one donor."),
          tags$li("Select a Cell type."), tags$li("Select a Comparison (target − reference)."), tags$li("Select a Readout.")),
        tags$p("Abundance describes marker amounts. Self-clustering describes proximity for one marker. Colocalization describes proximity between different markers."),
        tags$h3("Tables and Volcano"),
        tags$ol(tags$li("Click a result row in Tables to select a marker or pair."),
          tags$li("Open Volcano to see effects and statistical evidence."),
          tags$li("Click a volcano point to change the adjacent violin and boxplot.")),
        tags$p("The distribution shows available cell values by condition. Missing proximity scores are excluded, not replaced with zero."),
        tags$p("The displayed cells can differ from the original test group. The distribution does not calculate a new test."),
        tags$h3("Statistical limits"),
        tags$p("A positive effect means that the target value exceeds the reference value. Spatial differential effects are differences between medians."),
        tags$p("P values describe cells within the selected donor. Cells are not independent biological replicates."),
        tags$p("These exploratory results do not establish a treatment effect across a population of donors."),
        tags$p("FDR means false discovery rate. Table FDR controls do not change the original report thresholds."),
        tags$p("The Statistical methods and original thresholds section gives the stored method settings.")
      ),
      bslib::accordion_panel("Spatial patterns",
        tags$p("Spatial patterns is unavailable for abundance. Network is available only for colocalization."),
        tags$ol(tags$li("Select a spatial readout under Precomputed results."), tags$li("Open Spatial patterns."),
          tags$li("Select Matrix or Network.")),
        tags$p("Matrix shows differential effects or condition profiles. Condition profiles contain means, not differences between medians."),
        tags$ol(tags$li("For Network, select 2 to 20 Network markers."),
          tags$li("Set Minimum detection fraction if necessary."), tags$li("Set Minimum absolute mean proximity if necessary."),
          tags$li("Point to an edge to read its score and cell counts."),
          tags$li("Click an edge to select its pair for Tables or Volcano.")),
        tags$p("The network uses the same marker positions and scales for both conditions. Each condition uses the display filters separately."),
        tags$p("Edge color shows the mean log2 ratio. Edge width shows its absolute magnitude."),
        tags$p("Network means include missing scores as zero. Network thresholds control the display; they do not test statistical significance.")
      ),
      bslib::accordion_panel("Custom analysis",
        tags$ol(tags$li("Open Custom analysis."), tags$li("Select one donor."), tags$li("Select one cell type."),
          tags$li("Select the target and reference conditions."), tags$li("Select the readout."),
          tags$li("For spatial tests, set Minimum observed cells per condition if necessary."),
          tags$li("Read the exploratory-analysis notice."), tags$li("Select its confirmation checkbox."),
          tags$li("Click Run comparison.")),
        tags$p("The comparison uses all cells in the selected scope. It tests all dataset markers or different-marker pairs."),
        tags$p("Spatial tests apply isotype-background filtering before the minimum observation count. The analysis adjusts P values across the tested features."),
        tags$p("One job permits at most 2,000 cells, 160 markers, and 120 seconds. Larger comparisons require analysis outside the application."),
        tags$p("You can use other tabs during the analysis. Custom results do not replace precomputed results.")
      ),
      bslib::accordion_panel("Cell 3D, GIF, and PowerPoint model",
        tags$ol(tags$li("Open Cell 3D."), tags$li("Select a donor, cell type, and condition."),
          tags$li("Select the two ranking markers."), tags$li("Select a rank group."), tags$li("Select a cell with prepared graph data."),
          tags$li("Select one highlight marker."), tags$li("Click Display.")),
        tags$p("For a self-pair, use the same marker twice. Marker order does not change the rank."),
        tags$p("Ranking uses observed log2 ratios. Some ranked cells can have no prepared graph data."),
        tags$p("Cell and marker changes do not update the graph until you click Display."),
        tags$p("Open Appearance and export to set shared colors and the view/GIF opacity and dot sizes. Open PowerPoint export preset (GLB only) for separate model opacity and sizes. Click Apply appearance to save both sets of settings. Pending edits are not exported."),
        tags$p("Click Download 3D model (GLB) to save a model for PowerPoint. The file contains tetrahedra, colors, and opacity. No companion file is necessary. Lighting and transparency can differ between viewers."),
        tags$p("The coordinates describe an abstract graph layout. They are not measured physical positions or distances."),
        tags$p("A shared highlight is not a statistical test of colocalization."),
        tags$ol(tags$li("To export a GIF, open Rotating GIF settings."), tags$li("Set the duration in seconds."),
          tags$li("Set the frame rate in frames per second."), tags$li("Set the image size in pixels."),
          tags$li("Click Download rotating GIF.")),
        tags$p("The GIF shows one rotation of the displayed cell. Larger images and more frames increase export work.")
      ),
      bslib::accordion_panel("Plot controls and exports",
        tags$ol(tags$li("Point to a Plotly data point to see its values."),
          tags$li("Use the plot toolbar to zoom or restore the view.")),
        tags$p("Filtered table CSV uses the table filters and search. All stored rows CSV ignores these filters."),
        tags$p("Table filters do not change the volcano or spatial report thresholds."),
        tags$p("SVG is a vector image format. CSV stores table data. JSON stores structured information, including provenance."),
        tags$p("Original report SVG preserves the report design. Network SVG and Displayed edges CSV use the current network selections."),
        tags$p("A missing-sidecar message means that the required file is unavailable. Ask the administrator to check the prepared files.")
      )
    )
  )
}
