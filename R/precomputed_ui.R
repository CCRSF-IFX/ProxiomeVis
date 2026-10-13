precomputed_ui <- function(id, dataset_control = NULL, dataset_status = NULL) {
  ns <- shiny::NS(id)
  bslib::navset_bar(
    id = ns("report_page"), title = "ProxiomeVis", fillable = FALSE,
    navbar_options = bslib::navbar_options(bg = "#f8f9fa", theme = "light"),
    header = shiny::tagList(dataset_status,
      shiny::tags$head(shiny::tags$script(src = "precomputed.js"),
        shiny::tags$style(shiny::HTML(
          ".spatial-plot:has(.js-plotly-plot):not(:has(.recalculating)) .spatial-loading, .spatial-plot:has(.shiny-output-error) .spatial-loading {display:none;}")))),
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
              "Volcano rows",
              precomputed_download_group(
                shiny::downloadButton(ns("download_filtered_report_csv"), "Filtered volcano table CSV", class = "btn-sm btn-outline-primary")
              ),
              DT::DTOutput(ns("report_volcano_table"))
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
          "Spatial matrix",
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
              "Hover for values and availability; click a marker pair to inspect it below. Marker order and scales follow the frozen report alignment."),
            shiny::div(class = "spatial-plot",
              shiny::div(id = ns("spatial_loading"), class = "spatial-loading alert alert-info",
                role = "status", `aria-live` = "polite",
                shiny::span(class = "spinner-border spinner-border-sm me-2", `aria-hidden` = "true"),
                "Loading spatial matrix…"),
              plotly::plotlyOutput(ns("spatial_report_plot"), height = "750px")
            )
          )
        )
      ),
      bslib::card(full_screen = TRUE,
        class = "mt-3", bslib::card_header("Selected feature — shared across views"),
        bslib::layout_columns(
          col_widths = c(4, 8), fill = FALSE,
          shiny::tagList(
            shiny::selectizeInput(ns("detail_feature"), "Marker / pair", choices = character()),
            shiny::conditionalPanel("input.results_view == 'Spatial matrix'", ns = ns,
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
      )
    ),
    bslib::nav_panel("Custom analysis", hybrid_analysis_ui(ns("analysis"))),
    bslib::nav_panel("Cell 3D", cell_graph_ui(ns("graph"))),
    bslib::nav_spacer(),
    bslib::nav_item(dataset_control)
  )
}
