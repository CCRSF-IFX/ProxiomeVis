"""ProxiomeVis: the Python Shiny implementation."""

from __future__ import annotations

import asyncio
import math
from datetime import datetime
from itertools import combinations
from pathlib import Path
from queue import Empty, SimpleQueue
from time import perf_counter

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_plotly

from proxiome import (
    AppData,
    SpatialRetrieval,
    analysis_grouping_summary,
    apply_analysis_grouping,
    build_spatial_retrieval,
    calculate_differential,
    default_h5ad_path,
    default_pxl_spec,
    filter_pixelator_proximity,
    load_h5ad_data,
    load_pxl_proximity,
    mapping_for_column,
    new_analysis_grouping_config,
    read_component_layout,
    resolve_h5ad_path,
    sample_pxl_colocalization,
    select_colocalization_heatmap_markers,
    select_proximity_profile_markers,
    summarize_numeric,
    summarize_sample_colocalization,
    summarize_spatial,
    update_analysis_grouping_config,
)


APP_DIR = Path(__file__).resolve().parent
TEAL = "#176d73"
CORAL = "#c7503e"
PLOT_COLORS = ["#176d73", "#c7503e", "#c58a20", "#62879a", "#7b6aa2", "#56875d"]


def selectize(id: str, label: str, *, multiple: bool = False):
    return ui.input_selectize(id, label, [], multiple=multiple, remove_button=multiple)


def plot_pane(output_id: str, *, height: str = "520px", class_: str = "plot-pane-standard"):
    return ui.div(
        ui.div(
            ui.download_button(f"{output_id}_png", "PNG", class_="btn-sm btn-outline-secondary"),
            ui.download_button(f"{output_id}_pdf", "PDF", class_="btn-sm btn-outline-secondary"),
            class_="plot-pane-controls",
        ),
        output_widget(output_id, height=height),
        class_=f"plot-pane {class_}",
    )


def split_umap_height(facet_count: int, facet_columns: int) -> int:
    columns = max(1, min(facet_columns, max(1, facet_count)))
    rows = math.ceil(max(1, facet_count) / columns)
    return max(620, 120 + 320 * rows)


def label_embedding_axes(figure: go.Figure) -> go.Figure:
    for axis in figure.select_xaxes():
        axis.update(title_text="Embedding 1" if axis.title.text else None, showgrid=False)
    for axis in figure.select_yaxes():
        axis.update(title_text="Embedding 2" if axis.title.text else None, showgrid=False)
    return figure


def table_pane(output_id: str):
    return ui.div(ui.output_data_frame(output_id), class_="table-pane")


def differential_panel(prefix: str):
    return ui.div(
        ui.output_ui(f"{prefix}_summary"),
        ui.layout_columns(
            plot_pane(f"{prefix}_volcano", height="500px"),
            plot_pane(f"{prefix}_detail", height="500px"),
            col_widths=(6, 6),
        ),
        table_pane(f"{prefix}_table"),
    )


def qc_ui():
    sidebar = ui.sidebar(
        ui.accordion(
            ui.accordion_panel(
                "Filters",
                selectize("qc_sample_filter", "Sample", multiple=True),
            ),
            ui.accordion_panel(
                "Cutoffs",
                ui.input_numeric("qc_n_umi_cutoff", "n_umi cutoff", 10000, min=0, step=500),
                ui.input_numeric("qc_isotype_cutoff", "Isotype fraction cutoff", 0.001, min=0, max=1, step=0.0005),
            ),
            ui.accordion_panel(
                "Display",
                ui.input_select("qc_filter_y", "Filter count y-axis", {"count": "Number of cells", "fraction_loaded": "Fraction of loaded cells"}),
                ui.input_select("qc_metric", "Distribution metric", []),
            ),
            open=["Filters", "Cutoffs", "Display"],
        ),
        title="QC controls",
        width=300,
    )
    return ui.nav_panel(
        "QC",
        ui.layout_sidebar(
            sidebar,
            ui.navset_card_underline(
                ui.nav_panel("Filtering", ui.output_ui("qc_metric_row"), plot_pane("qc_filter_plot"), table_pane("qc_filter_table")),
                ui.nav_panel("Cell Calling", plot_pane("qc_molecule_rank_plot", height="560px")),
                ui.nav_panel("Distributions", plot_pane("qc_distribution_plot")),
                ui.nav_panel("Metadata", table_pane("qc_metadata_table")),
                id="qc_mode",
                title="QC",
                full_screen=True,
            ),
        ),
    )


def abundance_ui():
    sidebar = ui.sidebar(
        ui.panel_conditional(
            "input.abundance_mode === 'Observed'",
            ui.accordion(
                ui.accordion_panel(
                    "Display",
                    ui.input_select("abundance_embedding", "Embedding", []),
                    ui.input_select("abundance_color_by", "Color UMAP by", {
                        "abundance": "Marker abundance", "celltype_manual": "Cell type", "condition": "Analysis group", "sample_alias": "Sample"
                    }),
                    ui.panel_conditional(
                        "input.abundance_color_by === 'abundance'",
                        selectize("abundance_marker", "Marker"),
                    ),
                    ui.input_select("abundance_split_by", "Split UMAP by", {"": "None", "condition": "Analysis group", "sample_alias": "Sample"}),
                    ui.input_numeric("abundance_split_columns", "Split columns", 2, min=1, max=12),
                    ui.input_slider("abundance_point_size", "Dot size", 1, 8, 3, step=0.5),
                ),
                ui.accordion_panel(
                    "Filters",
                    selectize("abundance_condition_filter", "Analysis group", multiple=True),
                    selectize("abundance_celltype_filter", "Cell type", multiple=True),
                ),
                open=["Display", "Filters"],
            ),
        ),
        ui.panel_conditional(
            "input.abundance_mode === 'Marker Distributions'",
            ui.accordion(
                ui.accordion_panel(
                    "Display",
                    selectize("abundance_distribution_marker", "Marker"),
                    ui.input_numeric("abundance_distribution_columns", "Facet columns", 3, min=1, max=12),
                    ui.input_checkbox("abundance_distribution_show_jitter", "Show jitter dots", True),
                ),
                open="Display",
            ),
        ),
        ui.panel_conditional(
            "input.abundance_mode === 'Differential'",
            differential_controls("abundance", "marker"),
        ),
        title="Abundance controls",
        width=300,
    )
    return ui.nav_panel(
        "Abundance",
        ui.layout_sidebar(
            sidebar,
            ui.navset_card_underline(
                ui.nav_panel("Observed", ui.output_ui("abundance_metric_row"), plot_pane("abundance_umap", height="auto"), table_pane("abundance_table")),
                ui.nav_panel("Marker Distributions", plot_pane("abundance_distribution", height="680px"), table_pane("abundance_distribution_table")),
                ui.nav_panel(
                    "Cell Annotation",
                    plot_pane("abundance_composition"),
                    plot_pane("abundance_annotation_heatmap", height="680px"),
                    table_pane("abundance_composition_table"),
                ),
                ui.nav_panel("Differential", differential_panel("abundance_diff")),
                id="abundance_mode",
                title="Abundance",
                full_screen=True,
            ),
        ),
    )


def differential_controls(prefix: str, detail_label: str, *, show_run_button: bool = True):
    contrast_controls = [
        ui.input_select(f"{prefix}_diff_group_a", "Group A", []),
        ui.input_select(f"{prefix}_diff_group_b", "Group B (reference)", []),
        selectize(f"{prefix}_diff_celltype_filter", "Cell type", multiple=True),
        ui.input_checkbox(f"{prefix}_diff_stratify", "Stratify by cell type", False),
    ]
    if show_run_button:
        contrast_controls.append(
            ui.input_action_button(
                f"{prefix}_run_differential",
                "Run differential analysis",
                class_="btn-primary w-100",
            )
        )
    return ui.accordion(
        ui.accordion_panel(
            "Contrast",
            *contrast_controls,
        ),
        ui.accordion_panel(
            "Thresholds",
            ui.input_numeric(f"{prefix}_diff_fdr", "FDR threshold", 0.05, min=0, max=1, step=0.01),
            ui.input_numeric(f"{prefix}_diff_effect", "Minimum effect", 0.25, min=0, step=0.05),
            ui.input_numeric(f"{prefix}_diff_min", "Minimum observations per group", 3, min=1),
        ),
        ui.accordion_panel("Detail", selectize(f"{prefix}_diff_feature", f"Detail {detail_label}")),
        open=["Contrast", "Thresholds"],
    )


def spatial_retrieval_ui():
    sidebar = ui.sidebar(
        ui.accordion(
            ui.accordion_panel(
                "Population",
                selectize("spatial_retrieval_conditions", "Analysis group", multiple=True),
                selectize("spatial_retrieval_samples", "Sample", multiple=True),
                selectize("spatial_retrieval_celltypes", "Cell type", multiple=True),
            ),
            ui.accordion_panel(
                "Markers",
                ui.input_select(
                    "spatial_retrieval_marker_mode",
                    "Number of markers",
                    {
                        "all": "All markers",
                        "top": "Top abundance markers",
                        "manual": "Selected markers",
                    },
                    selected="all",
                ),
                ui.panel_conditional(
                    "input.spatial_retrieval_marker_mode === 'top'",
                    ui.input_numeric("spatial_retrieval_marker_count", "Top marker count", 40, min=1),
                ),
                ui.panel_conditional(
                    "input.spatial_retrieval_marker_mode === 'manual'",
                    selectize("spatial_retrieval_markers", "Markers", multiple=True),
                ),
            ),
            open=["Population", "Markers"],
        ),
        ui.input_task_button("retrieve_spatial_data", "Retrieve Spatial Data", class_="btn-primary w-100"),
        title="Retrieval controls",
        width=320,
    )
    return ui.nav_panel(
        "Retrieve Data",
        ui.layout_sidebar(
            sidebar,
            ui.div(
                ui.output_ui("spatial_retrieval_status"),
                ui.card(
                    ui.card_header("Active retrieval"),
                    ui.card_body(ui.output_ui("spatial_retrieval_summary")),
                ),
                ui.card(
                    ui.card_header("Retrieved cells"),
                    ui.card_body(ui.output_data_frame("spatial_retrieval_table")),
                ),
                class_="p-3",
            ),
        ),
    )


def clustering_ui():
    sidebar = ui.sidebar(
        ui.panel_conditional(
            "input.clustering_mode === 'Observed' || input.clustering_mode === 'Per Marker'",
            ui.accordion(
                ui.accordion_panel("Display", selectize("clustering_marker", "Marker")),
                ui.accordion_panel(
                    "Filters",
                    selectize("clustering_condition_filter", "Analysis group", multiple=True),
                    selectize("clustering_celltype_filter", "Cell type", multiple=True),
                ),
                open=["Display", "Filters"],
            ),
        ),
        ui.panel_conditional(
            "input.clustering_mode === 'Summary Heatmap'",
            ui.accordion(
                ui.accordion_panel("Display", ui.input_numeric("clustering_heatmap_marker_count", "Top markers", 20, min=2, max=40)),
                ui.accordion_panel(
                    "Filters",
                    selectize("clustering_heatmap_condition_filter", "Analysis group", multiple=True),
                    selectize("clustering_heatmap_celltype_filter", "Cell type", multiple=True),
                ),
                open=["Display", "Filters"],
            ),
        ),
        ui.panel_conditional(
            "input.clustering_mode === 'Differential'",
            differential_controls("clustering", "marker", show_run_button=False),
        ),
        title="Clustering controls",
        width=300,
    )
    return ui.nav_panel(
        "Clustering",
        ui.layout_sidebar(
            sidebar,
            ui.div(
                ui.output_ui("clustering_retrieval_notice"),
                ui.navset_card_underline(
                    ui.nav_panel("Observed", plot_pane("clustering_plot"), table_pane("clustering_table")),
                    ui.nav_panel("Per Marker", plot_pane("clustering_per_marker"), table_pane("clustering_per_marker_table")),
                    ui.nav_panel("Summary Heatmap", plot_pane("clustering_summary_heatmap", height="680px"), table_pane("clustering_summary_table")),
                    ui.nav_panel("Differential", differential_panel("clustering_diff")),
                    id="clustering_mode",
                    title="Clustering",
                    full_screen=True,
                ),
            ),
        ),
    )


def colocalization_ui():
    sidebar = ui.sidebar(
        ui.panel_conditional(
            "input.colocalization_mode === 'Observed'",
            ui.accordion(
                ui.accordion_panel("Cell population", selectize("coloc_celltype_filter", "Cell population", multiple=True)),
                ui.accordion_panel(
                    "View",
                    ui.input_select("coloc_preset", "Heatmap preset", {"report": "Notebook-compatible", "custom": "Custom"}),
                    ui.input_select("coloc_scope", "Heatmap scope", {"condition": "Analysis group summary", "sample_alias": "Sample summary", "celltype": "Cell type focus"}),
                    ui.panel_conditional("input.coloc_scope === 'celltype'", selectize("coloc_celltype_focus", "Cell type focus")),
                    ui.input_select("coloc_view", "Group display", {"focused": "Focused group", "compare": "Compare groups"}),
                    ui.panel_conditional("input.coloc_view === 'focused'", selectize("coloc_focus_group", "Displayed group")),
                    ui.panel_conditional(
                        "input.coloc_view === 'compare'",
                        selectize("coloc_compare_groups", "Groups to compare", multiple=True),
                        ui.input_select("coloc_reference", "Clustering reference", []),
                    ),
                    ui.input_select("coloc_ordering", "Marker ordering", {"ward": "Ward", "complete": "Complete", "average": "Average", "single": "Single"}),
                ),
                ui.accordion_panel(
                    "Marker selection",
                    ui.input_select(
                        "coloc_marker_mode",
                        "Marker set",
                        {
                            "profile": "Notebook proximity profile",
                            "abundance": "Top abundance markers",
                            "manual": "Selected markers",
                        },
                    ),
                    ui.panel_conditional(
                        "input.coloc_marker_mode === 'profile'",
                        ui.input_numeric("coloc_top_pairs", "Strongest proximity pairs", 60, min=1, max=500),
                        ui.help_text("Markers are taken from pairs detected in more than 50% of cells."),
                    ),
                    ui.panel_conditional(
                        "input.coloc_marker_mode === 'abundance'",
                        ui.input_numeric("coloc_top_markers", "Top abundance markers", 40, min=2, max=40),
                    ),
                    ui.panel_conditional(
                        "input.coloc_marker_mode === 'manual'",
                        selectize("coloc_markers", "Heatmap markers", multiple=True),
                    ),
                ),
                ui.accordion_panel(
                    "Appearance",
                    ui.input_numeric("coloc_legend_min", "Legend minimum", -1, step=0.1),
                    ui.input_numeric("coloc_legend_max", "Legend maximum", 1, step=0.1),
                    selectize("coloc_detail_pair", "Pair detail"),
                ),
                ui.accordion_panel(
                    "Filters",
                    selectize("coloc_condition_filter", "Analysis group", multiple=True),
                    ui.input_checkbox("coloc_pixelator_filter", "Apply Pixelator proximity filters", False),
                    ui.input_numeric("coloc_min_fraction", "Minimum marker fraction", 0.001, min=0, max=1, step=0.0005),
                    ui.input_numeric("coloc_min_count", "Minimum marker count", 0, min=0),
                    ui.input_numeric("coloc_min_cells", "Minimum cells per pair", 1, min=1),
                ),
                ui.accordion_panel(
                    "Advanced",
                    ui.input_select("coloc_mean_type", "Heatmap mean", {"population": "Population mean", "detected": "Detected-cell mean"}),
                ),
                ui.accordion_panel(
                    "Interpretation",
                    ui.p(ui.strong("Positive:"), " closer than expected by chance."),
                    ui.p(ui.strong("Zero:"), " approximately random spatial organization."),
                    ui.p(ui.strong("Negative:"), " spatial segregation."),
                ),
                open=["Cell population", "View", "Marker selection"],
            ),
        ),
        ui.panel_conditional(
            "input.colocalization_mode === 'Differential'",
            ui.accordion(
                ui.accordion_panel(
                    "Cell population",
                    selectize("coloc_diff_celltype_filter", "Cell population", multiple=True),
                    ui.input_select("coloc_diff_mean", "Sample summary", {"population": "Population mean", "detected": "Detected-cell mean"}),
                ),
                ui.accordion_panel(
                    "Contrast",
                    ui.input_select("coloc_diff_group_a", "Group A", []),
                    ui.input_select("coloc_diff_group_b", "Group B (reference)", []),
                    ui.input_numeric("coloc_diff_min_samples", "Minimum samples per group", 2, min=1),
                ),
                ui.accordion_panel(
                    "Pair display",
                    ui.input_select("coloc_diff_pair_scope", "Pairs shown", {"all": "All marker pairs", "anchor": "Pairs containing one marker"}),
                    selectize("coloc_diff_anchor", "Marker"),
                    selectize("coloc_diff_pair", "Detail pair"),
                ),
                ui.accordion_panel(
                    "Thresholds",
                    ui.input_numeric("coloc_diff_fdr", "FDR threshold", 0.05, min=0, max=1, step=0.01),
                    ui.input_numeric("coloc_diff_effect", "Minimum median-sample difference", 0.25, min=0, step=0.05),
                ),
                open=["Cell population", "Contrast", "Pair display"],
            ),
        ),
        ui.panel_conditional(
            "input.colocalization_mode === '3D Layout'",
            ui.accordion(
                ui.accordion_panel(
                    "Cell",
                    selectize("coloc_3d_sample", "Sample"),
                    selectize("coloc_3d_celltypes", "Cell type", multiple=True),
                    selectize("coloc_3d_component", "Cell/component"),
                    ui.input_numeric("coloc_3d_max_background", "Max background nodes", 7000, min=0, max=50000, step=500),
                ),
                ui.accordion_panel("Markers", selectize("coloc_3d_markers", "Highlighted markers", multiple=True)),
                open=["Cell", "Markers"],
            ),
        ),
        title="Colocalization controls",
        width=300,
    )
    return ui.nav_panel(
        "Colocalization",
        ui.layout_sidebar(
            sidebar,
            ui.div(
                ui.output_ui("colocalization_retrieval_notice"),
                ui.navset_card_underline(
                    ui.nav_panel(
                        "Observed",
                        ui.output_ui("coloc_notice"),
                        plot_pane(
                            "coloc_heatmap",
                            height="auto",
                            class_="plot-pane-scroll coloc-heatmap-pane",
                        ),
                        ui.card(
                            ui.card_header("Pair detail"),
                            ui.card_body(ui.output_ui("coloc_pair_metrics"), plot_pane("coloc_pair_detail"), table_pane("coloc_pair_table")),
                        ),
                        table_pane("coloc_table"),
                    ),
                    ui.nav_panel("Differential", ui.output_ui("coloc_diff_method"), differential_panel("coloc_diff")),
                    ui.nav_panel("3D Layout", plot_pane("coloc_3d_layout", height="640px"), table_pane("coloc_3d_table")),
                    id="colocalization_mode",
                    title="Colocalization",
                    full_screen=True,
                ),
            ),
        ),
    )


def patch_ui():
    sidebar = ui.sidebar(
        ui.accordion(ui.accordion_panel("Markers", selectize("patch_label_filter", "Marker class", multiple=True)), open="Markers"),
        title="Patch controls",
        width=300,
    )
    return ui.nav_panel(
        "Patch Analysis",
        ui.layout_sidebar(
            sidebar,
            ui.navset_card_underline(
                ui.nav_panel("Markers", ui.output_ui("patch_metric_row"), plot_pane("patch_marker_plot"), table_pane("patch_marker_table")),
                ui.nav_panel("Raji Signal", plot_pane("patch_raji_plot"), table_pane("patch_raji_table")),
                ui.nav_panel("Patch Burden", table_pane("patch_burden_table")),
                id="patch_mode",
                title="Patch Analysis",
                full_screen=True,
            ),
        ),
    )


def activity_log_ui():
    return ui.nav_panel(
        "Activity Log",
        ui.div(
            ui.div(
                ui.div(
                    ui.h2("Activity Log", class_="mb-1"),
                    ui.p("Live session activity for data loading and PXL queries.", class_="text-muted mb-0"),
                ),
                ui.input_action_button("clear_activity_log", "Clear log", class_="btn-outline-secondary"),
                class_="d-flex justify-content-between align-items-center mb-3",
            ),
            ui.p(ui.strong("Latest: "), ui.output_text("activity_status", inline=True)),
            ui.output_data_frame("activity_log_table"),
            class_="container-fluid py-3",
        ),
    )



def data_popover():
    return ui.nav_control(
        ui.popover(
            ui.input_action_button("data_button", "Data", class_="btn btn-outline-light data-source-button"),
            ui.div(
                ui.input_text_area("h5ad_path", "Processed .h5ad path", default_h5ad_path(), rows=4),
                ui.input_text_area(
                    "pxl_path",
                    "PXL path(s) for proximity and cellgraph data",
                    default_pxl_spec(),
                    rows=3,
                ),
                ui.div(
                    ui.input_action_button("inspect_h5ad", "Inspect", class_="btn-outline-secondary"),
                    ui.input_task_button("load_h5ad", "Load Data"),
                    class_="data-source-actions",
                ),
                ui.output_ui("load_status"),
                ui.hr(),
                ui.input_action_button(
                    "configure_analysis_grouping",
                    "Analysis grouping…",
                    class_="btn-outline-primary w-100",
                ),
                ui.output_text("analysis_grouping_summary", inline=True),
                ui.output_ui("source_summary"),
                class_="data-source-popover",
            ),
            title="Data source",
            placement="bottom",
            options={"html": True},
            class_="data-source-popover-shell",
        )
    )


app_ui = ui.page_navbar(
    qc_ui(),
    abundance_ui(),
    ui.nav_panel(
        "Spatial Metrics",
        ui.navset_tab(spatial_retrieval_ui(), clustering_ui(), colocalization_ui(), id="spatial_metric_readout"),
    ),
    patch_ui(),
    activity_log_ui(),
    ui.nav_spacer(),
    data_popover(),
    title="ProxiomeVis",
    id="readout_tab",
    fillable=["QC", "Abundance", "Spatial Metrics", "Patch Analysis", "Activity Log"],
    header=ui.include_css(APP_DIR / "www" / "proixome.css"),
    navbar_options=ui.navbar_options(bg=TEAL, theme="dark", underline=True),
)


def empty_figure(message: str = "Load a processed H5AD from the Data menu.") -> go.Figure:
    figure = go.Figure()
    figure.add_annotation(text=message, x=0.5, y=0.5, xref="paper", yref="paper", showarrow=False)
    figure.update_xaxes(visible=False)
    figure.update_yaxes(visible=False)
    return style_figure(figure)


def style_figure(figure: go.Figure, *, height: int | None = None) -> go.Figure:
    figure.update_layout(
        template="plotly_white",
        colorway=PLOT_COLORS,
        margin=dict(l=70, r=50, t=55, b=80),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        height=height,
    )
    return figure


def clickable_volcano(figure: go.Figure, input_id: str, session: Session):
    widget = go.FigureWidget(figure)

    def select_point(trace, points, _state):
        if points.point_inds:
            ui.update_selectize(
                input_id,
                selected=str(trace.customdata[points.point_inds[0]][0]),
                session=session,
            )

    for trace in widget.data:
        trace.on_click(select_point)
    return widget


def selected(value, choices: pd.Series | list[str]) -> list[str]:
    all_values = [str(item) for item in pd.Series(choices).dropna().unique()]
    if value is None or value == "" or value == () or value == []:
        return all_values
    return [str(item) for item in (value if isinstance(value, (list, tuple)) else [value])]


def metric_boxes(items: list[tuple[str, str]]):
    return ui.layout_columns(
        *[ui.value_box(label, value, theme="teal", height="110px") for label, value in items],
        col_widths=tuple(12 // len(items) for _ in items),
    )


def embedding_columns(metadata: pd.DataFrame) -> dict[str, tuple[str, str]]:
    result = {}
    allowed = ("umap", "pca", "harmony", "tsne")
    for column in metadata.columns:
        match = column.rsplit("_", 1)
        if (
            len(match) == 2
            and match[1] == "1"
            and f"{match[0]}_2" in metadata
            and any(name in match[0].lower() for name in allowed)
        ):
            result[match[0]] = (column, f"{match[0]}_2")
    return result


def server(input: Inputs, output: Outputs, session: Session):
    data_state: reactive.Value[AppData | None] = reactive.Value(None)
    grouping_state: reactive.Value[dict | None] = reactive.Value(None)
    spatial_retrieval_state: reactive.Value[SpatialRetrieval | None] = reactive.Value(None)
    inspect_message = reactive.Value("No data loaded.")
    activity_queue: SimpleQueue[dict] = SimpleQueue()
    activity_state: reactive.Value[tuple[dict, ...]] = reactive.Value(())

    def log_activity(operation: str, status: str, details: str = "", elapsed: float | None = None):
        activity_queue.put({
            "time": datetime.now().astimezone().strftime("%H:%M:%S"),
            "operation": operation,
            "status": status,
            "details": details,
            "seconds": round(elapsed, 2) if elapsed is not None else None,
        })

    log_activity("Session", "Ready", "Waiting for data.")

    @reactive.effect
    def _drain_activity_queue():
        reactive.invalidate_later(0.5)
        pending = []
        while True:
            try:
                pending.append(activity_queue.get_nowait())
            except Empty:
                break
        if pending:
            activity_state.set(tuple([*activity_state.get(), *pending][-500:]))

    @ui.bind_task_button(button_id="load_h5ad")
    @reactive.extended_task
    async def load_task(path: str, pxl_path: str) -> AppData:
        started = perf_counter()
        log_activity("Load data", "Started", f"Reading {Path(path).name} and resolving PXL files.")
        try:
            loaded = await asyncio.to_thread(load_h5ad_data, path, pxl_spec=pxl_path or None)
        except Exception as error:
            log_activity("Load data", "Failed", str(error), perf_counter() - started)
            raise
        log_activity(
            "Load data",
            "Completed",
            f"{loaded.source['n_cells']:,} cells, {len(loaded.marker_options):,} markers, "
            f"{len(loaded.pxl_files):,} PXL files.",
            perf_counter() - started,
        )
        return loaded

    @ui.bind_task_button(button_id="retrieve_spatial_data")
    @reactive.extended_task
    async def spatial_retrieval_task(data: AppData, request: tuple) -> tuple[int, SpatialRetrieval]:
        conditions, samples, celltypes, marker_mode, marker_count, markers = request
        started = perf_counter()
        log_activity(
            "Spatial retrieval",
            "Started",
            f"Resolving {len(celltypes):,} cell type(s) and "
            f"{'all' if marker_mode == 'all' else marker_count if marker_mode == 'top' else len(markers)} markers.",
        )
        try:
            retrieval = await asyncio.to_thread(
                build_spatial_retrieval,
                data,
                conditions=conditions,
                samples=samples,
                celltypes=celltypes,
                marker_mode=marker_mode,
                n_markers=marker_count or len(data.marker_options),
                plot_markers=markers,
                retrieved_at=datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"),
                request=request,
            )
        except Exception as error:
            log_activity("Spatial retrieval", "Failed", str(error), perf_counter() - started)
            raise
        log_activity(
            "Spatial retrieval",
            "Completed",
            f"{len(retrieval.metadata):,} cells and {len(retrieval.markers):,} markers.",
            perf_counter() - started,
        )
        return id(data), retrieval

    @reactive.effect
    @reactive.event(input.inspect_h5ad)
    def _inspect_h5ad():
        try:
            path = resolve_h5ad_path(input.h5ad_path())
            inspect_message.set(f"Ready: {path.name}, {path.stat().st_size / 1024**2:.1f} MiB.")
            log_activity("Inspect H5AD", "Completed", f"{path.name}: {path.stat().st_size / 1024**2:.1f} MiB.")
        except Exception as error:
            inspect_message.set(str(error))
            log_activity("Inspect H5AD", "Failed", str(error))

    @reactive.effect
    @reactive.event(input.load_h5ad)
    def _start_load():
        inspect_message.set("Loading processed AnnData…")
        load_task.invoke(input.h5ad_path(), input.pxl_path())

    def current_spatial_request() -> tuple | None:
        data = data_state.get()
        if data is None:
            return None
        metadata = data.metadata
        marker_value = input.spatial_retrieval_markers()
        markers = tuple(
            map(str, marker_value if isinstance(marker_value, (list, tuple)) else [marker_value])
        ) if marker_value else ()
        marker_mode = str(input.spatial_retrieval_marker_mode() or "all")
        return (
            tuple(sorted(selected(input.spatial_retrieval_conditions(), metadata["condition"]))),
            tuple(sorted(selected(input.spatial_retrieval_samples(), metadata["sample_alias"]))),
            tuple(sorted(selected(input.spatial_retrieval_celltypes(), metadata["celltype_manual"]))),
            marker_mode,
            max(1, int(input.spatial_retrieval_marker_count() or 40)) if marker_mode == "top" else None,
            tuple(sorted(set(markers))) if marker_mode == "manual" else (),
        )

    @reactive.effect
    @reactive.event(input.retrieve_spatial_data)
    def _start_spatial_retrieval():
        data = data_state.get()
        request = current_spatial_request()
        if data is None:
            ui.notification_show("Load H5AD and PXL data before retrieving spatial data.", type="warning")
            return
        if not data.pxl_files:
            ui.notification_show("Assign one or more PXL files before retrieving spatial data.", type="warning")
            return
        spatial_retrieval_task.invoke(data, request)

    @reactive.effect
    def _activate_loaded_data():
        loaded = load_task.result()
        data_state.set(loaded)
        spatial_retrieval_state.set(None)
        inspect_message.set(
            f"Loaded {loaded.source['n_cells']:,} cells and {len(loaded.marker_options):,} markers; "
            f"{len(loaded.pxl_files):,} PXL file(s) assigned."
        )
        grouping_state.set(new_analysis_grouping_config(loaded.metadata))

    @reactive.effect
    def _activate_spatial_retrieval():
        source_id, retrieval = spatial_retrieval_task.result()
        current_data = data_state.get()
        if current_data is None or id(current_data) != source_id:
            log_activity("Spatial retrieval", "Discarded", "The source dataset changed before retrieval completed.")
            return
        spatial_retrieval_state.set(retrieval)
        ui.notification_show(
            f"Spatial retrieval ready: {len(retrieval.metadata):,} cells and {len(retrieval.markers):,} markers.",
            type="message",
        )

    @output
    @render.ui
    def load_status():
        return ui.div(inspect_message(), class_="rds-load-status")

    @reactive.effect
    @reactive.event(input.clear_activity_log)
    def _clear_activity_log():
        while True:
            try:
                activity_queue.get_nowait()
            except Empty:
                break
        activity_state.set(())
        log_activity("Activity log", "Cleared")

    @output
    @render.text
    def activity_status():
        events = activity_state.get()
        if not events:
            return "No activity yet."
        latest = events[-1]
        return f"{latest['operation']} — {latest['status']}"

    @output
    @render.data_frame
    def activity_log_table():
        events = activity_state.get()
        columns = ["time", "operation", "status", "details", "seconds"]
        frame = pd.DataFrame(events, columns=columns).iloc[::-1].reset_index(drop=True)
        return render.DataGrid(frame, filters=True, width="100%", height="620px")

    @output
    @render.ui
    def spatial_retrieval_status():
        data = data_state.get()
        retrieval = spatial_retrieval_state.get()
        if data is None:
            return ui.div("Load H5AD and PXL data first.", class_="alert alert-warning")
        if not data.pxl_files:
            return ui.div("Assign PXL files in Data before retrieving spatial data.", class_="alert alert-warning")
        if retrieval is None:
            return ui.div("Choose a population and retrieve spatial data to enable the spatial views.", class_="alert alert-warning")
        if retrieval.request != current_spatial_request():
            return ui.div(
                "Retrieval settings changed. The current visualizations still use the active retrieval; "
                "click Retrieve Spatial Data to replace it.",
                class_="alert alert-warning",
            )
        return ui.div("The active retrieval matches the current settings.", class_="alert alert-success")

    @output
    @render.ui
    def spatial_retrieval_summary():
        retrieval = spatial_retrieval_state.get()
        if retrieval is None:
            return ui.p("No spatial data have been retrieved.", class_="text-muted")
        metadata = retrieval.metadata
        return ui.div(
            metric_boxes([
                ("Cells", f"{len(metadata):,}"),
                ("Samples", f"{metadata['sample_alias'].nunique():,}"),
                ("Cell types", f"{metadata['celltype_manual'].nunique():,}"),
                ("Markers", f"{len(retrieval.markers):,}"),
            ]),
            ui.p(
                f"Retrieved {retrieval.retrieved_at} · marker selection: {retrieval.marker_mode}.",
                class_="text-muted small mb-0",
            ),
        )

    @output
    @render.data_frame
    def spatial_retrieval_table():
        retrieval = spatial_retrieval_state.get()
        if retrieval is None:
            return render.DataGrid(pd.DataFrame())
        columns = [
            column for column in ("component", "sample_alias", "condition", "celltype_manual")
            if column in retrieval.metadata
        ]
        return render.DataGrid(
            retrieval.metadata[columns].head(2000),
            filters=True,
            width="100%",
            height="520px",
        )

    def active_retrieval_notice():
        retrieval = spatial_retrieval_state.get()
        if retrieval is None:
            return ui.div(
                "Retrieve data in Spatial Metrics > Retrieve Data to enable this view.",
                class_="alert alert-warning m-3 mb-0",
            )
        return ui.div(
            f"Active retrieval: {len(retrieval.metadata):,} cells, "
            f"{retrieval.metadata['sample_alias'].nunique():,} samples, and "
            f"{len(retrieval.markers):,} markers. View controls can only narrow this scope.",
            class_="alert alert-info m-3 mb-0 py-2",
        )

    @output
    @render.ui
    def clustering_retrieval_notice():
        return active_retrieval_notice()

    @output
    @render.ui
    def colocalization_retrieval_notice():
        return active_retrieval_notice()

    @output(id="analysis_grouping_summary")
    @render.text
    def analysis_grouping_summary_text():
        return analysis_grouping_summary(grouping_state.get())

    @output
    @render.ui
    def source_summary():
        data = data_state.get()
        if data is None:
            return ui.p("Load a processed H5AD file to begin.", class_="text-muted small")
        return ui.div(
            ui.hr(),
            ui.p(ui.strong(data.source.get("display_name", "AnnData"))),
            ui.p(f"{data.source['n_cells']:,} cells · {data.source['n_markers']:,} markers"),
            ui.p(f"Spatial metrics: {'available from PXL' if data.source.get('has_spatial_metrics') else 'assign PXL files'}"),
            ui.p(f"Patch analysis: {'available' if data.source.get('has_patch_analysis') else 'not stored in H5AD'}"),
            ui.p(f"PXL files: {len(data.pxl_files):,}"),
            ui.p(f"Grouping: {data.source.get('analysis_group_label', 'condition')}"),
            class_="source-chip",
        )

    @reactive.effect
    def _refresh_inputs():
        data = data_state.get()
        if data is None:
            return
        metadata = data.metadata
        markers = list(data.marker_options)
        conditions = sorted(metadata["condition"].dropna().astype(str).unique())
        celltypes = sorted(metadata["celltype_manual"].dropna().astype(str).unique())
        samples = sorted(metadata["sample_alias"].dropna().astype(str).unique())
        embeddings = list(embedding_columns(metadata))
        numeric_qc = [
            column for column in ("n_umi", "n_edges", "reads_in_component", "isotype_fraction", "tau")
            if column in metadata and pd.api.types.is_numeric_dtype(metadata[column])
        ]

        ui.update_selectize("qc_sample_filter", choices=samples, selected=samples, server=True)
        ui.update_select("qc_metric", choices=numeric_qc, selected=numeric_qc[0] if numeric_qc else None)
        ui.update_select("abundance_embedding", choices=embeddings, selected=embeddings[0] if embeddings else None)
        for input_id in (
            "abundance_marker", "abundance_distribution_marker", "abundance_diff_feature",
        ):
            ui.update_selectize(input_id, choices=markers, selected=markers[0] if markers else None, server=True)
        ui.update_selectize("abundance_condition_filter", choices=conditions, selected=conditions, server=True)
        ui.update_selectize("abundance_celltype_filter", choices=celltypes, selected=celltypes, server=True)
        for input_id in ("abundance_diff_celltype_filter",):
            ui.update_selectize(input_id, choices=celltypes, selected=celltypes[:1], server=True)
        for prefix in ("abundance",):
            ui.update_select(f"{prefix}_diff_group_a", choices=conditions, selected=conditions[0] if conditions else None)
            ui.update_select(f"{prefix}_diff_group_b", choices=conditions, selected=conditions[min(1, len(conditions) - 1)] if conditions else None)
        ui.update_selectize("spatial_retrieval_conditions", choices=conditions, selected=conditions, server=True)
        ui.update_selectize("spatial_retrieval_samples", choices=samples, selected=samples, server=True)
        ui.update_selectize("spatial_retrieval_celltypes", choices=celltypes, selected=celltypes[:1], server=True)
        ui.update_selectize(
            "spatial_retrieval_markers",
            choices=markers,
            selected=markers[: min(40, len(markers))],
            server=True,
        )
        marker_table = data.patch.get("marker_unmixing")
        labels = sorted(marker_table["label"].dropna().astype(str).unique()) if marker_table is not None and "label" in marker_table else []
        ui.update_selectize("patch_label_filter", choices=labels, selected=labels, server=True)

    @reactive.effect
    def _refresh_spatial_inputs():
        retrieval = spatial_retrieval_state.get()
        if retrieval is None:
            return
        metadata = retrieval.metadata
        markers = list(retrieval.markers)
        conditions = sorted(metadata["condition"].dropna().astype(str).unique())
        celltypes = sorted(metadata["celltype_manual"].dropna().astype(str).unique())
        samples = sorted(metadata["sample_alias"].dropna().astype(str).unique())
        pairs = [f"{marker_1} / {marker_2}" for marker_1, marker_2 in combinations(sorted(markers), 2)]

        for input_id in ("clustering_marker", "clustering_diff_feature", "coloc_diff_anchor"):
            ui.update_selectize(input_id, choices=markers, selected=markers[0] if markers else None, server=True)
        for input_id in ("coloc_markers", "coloc_3d_markers"):
            ui.update_selectize(input_id, choices=markers, selected=markers[: min(15, len(markers))], server=True)
        for input_id in (
            "clustering_condition_filter", "clustering_heatmap_condition_filter", "coloc_condition_filter",
        ):
            ui.update_selectize(input_id, choices=conditions, selected=conditions, server=True)
        for input_id in (
            "clustering_celltype_filter", "clustering_heatmap_celltype_filter",
            "clustering_diff_celltype_filter", "coloc_celltype_filter",
            "coloc_diff_celltype_filter", "coloc_3d_celltypes",
        ):
            ui.update_selectize(input_id, choices=celltypes, selected=celltypes[:1], server=True)
        ui.update_select("clustering_diff_group_a", choices=conditions, selected=conditions[0] if conditions else None)
        ui.update_select(
            "clustering_diff_group_b",
            choices=conditions,
            selected=conditions[min(1, len(conditions) - 1)] if conditions else None,
        )
        ui.update_select("coloc_diff_group_a", choices=conditions, selected=conditions[0] if conditions else None)
        ui.update_select(
            "coloc_diff_group_b",
            choices=conditions,
            selected=conditions[min(1, len(conditions) - 1)] if conditions else None,
        )
        ui.update_select("coloc_reference", choices=conditions, selected=conditions[0] if conditions else None)
        ui.update_selectize("coloc_celltype_focus", choices=celltypes, selected=celltypes[0] if celltypes else None, server=True)
        ui.update_selectize("coloc_focus_group", choices=conditions, selected=conditions[0] if conditions else None, server=True)
        ui.update_selectize("coloc_compare_groups", choices=conditions, selected=conditions[:6], server=True)
        ui.update_selectize("coloc_detail_pair", choices=pairs, selected=pairs[0] if pairs else None, server=True)
        ui.update_selectize("coloc_diff_pair", choices=pairs, selected=pairs[0] if pairs else None, server=True)
        ui.update_selectize("coloc_3d_sample", choices=samples, selected=samples[0] if samples else None, server=True)

    @reactive.effect
    def _refresh_colocalization_groups():
        retrieval = spatial_retrieval_state.get()
        if retrieval is None:
            return
        scope = input.coloc_scope() or "condition"
        group_col = "sample_alias" if scope == "sample_alias" else "condition"
        groups = sorted(retrieval.metadata[group_col].dropna().astype(str).unique())
        selected_group = groups[0] if groups else None
        ui.update_selectize(
            "coloc_focus_group", choices=groups, selected=selected_group, server=True
        )
        ui.update_selectize(
            "coloc_compare_groups", choices=groups, selected=groups[:6], server=True
        )
        ui.update_select("coloc_reference", choices=groups, selected=selected_group)

    @reactive.effect
    def _update_3d_components():
        data = data_state.get()
        retrieval = spatial_retrieval_state.get()
        if data is None or retrieval is None or not input.coloc_3d_sample():
            return
        rows = retrieval.metadata[
            retrieval.metadata["sample_alias"].astype(str) == str(input.coloc_3d_sample())
        ]
        celltypes = selected(input.coloc_3d_celltypes(), rows["celltype_manual"])
        rows = rows[rows["celltype_manual"].astype(str).isin(celltypes)]
        components = rows["component"].astype(str).tolist()
        if data.component_layouts and not data.pxl_files:
            components = [component for component in components if component in data.component_layouts]
        ui.update_selectize("coloc_3d_component", choices=components, selected=components[0] if components else None, server=True)

    def grouping_editor_rows() -> tuple[str, str, pd.DataFrame]:
        config = grouping_state.get()
        if config is None:
            return "column", "condition", pd.DataFrame()
        mode = input.analysis_group_mode() or config["mode"]
        column = input.analysis_group_column() or config["column"]
        if column not in config["columns"]:
            column = config["column"]
        mapping = (
            config["mapping"]
            if mode == "custom" and config["mode"] == "custom" and column == config["column"]
            else mapping_for_column(config["source"], column)
        )
        rows = config["source"][["sample_alias", column]].copy()
        rows = rows.rename(columns={column: "source_value"})
        rows["analysis_group"] = rows["sample_alias"].astype(str).map(mapping)
        return mode, column, rows

    @reactive.effect
    @reactive.event(input.configure_analysis_grouping)
    def _configure_analysis_grouping():
        config = grouping_state.get()
        if config is None:
            ui.notification_show("Load data before changing analysis grouping.", type="warning")
            return
        ui.modal_show(ui.modal(
            ui.input_radio_buttons(
                "analysis_group_mode",
                "Grouping mode",
                {"column": "Use metadata column", "custom": "Edit sample groups"},
                selected=config["mode"],
                inline=True,
            ),
            ui.input_select(
                "analysis_group_column",
                "Sample-level metadata column",
                config["columns"],
                selected=config["column"],
            ),
            ui.output_ui("analysis_group_editor"),
            title="Analysis grouping",
            footer=ui.tags.div(
                ui.modal_button("Cancel"),
                ui.input_action_button("reset_analysis_grouping", "Reset to condition", class_="btn-outline-secondary"),
                ui.input_action_button("apply_analysis_grouping", "Apply grouping", class_="btn-primary"),
                class_="d-flex gap-2",
            ),
            easy_close=True,
            size="l",
        ))

    @output
    @render.ui
    def analysis_group_editor():
        mode, column, rows = grouping_editor_rows()
        editable = mode == "custom"
        body = []
        for index, row in rows.reset_index(drop=True).iterrows():
            value = str(row.analysis_group)
            body.append(ui.tags.tr(
                ui.tags.td(str(row.sample_alias)),
                ui.tags.td("Unassigned" if pd.isna(row.source_value) else str(row.source_value)),
                ui.tags.td(
                    ui.input_text(f"analysis_group_value_{index}", None, value=value, width="100%")
                    if editable else value
                ),
            ))
        return ui.tags.table(
            ui.tags.thead(ui.tags.tr(
                ui.tags.th("Sample"),
                ui.tags.th(column),
                ui.tags.th("Analysis group"),
            )),
            ui.tags.tbody(*body),
            class_="table table-sm align-middle",
        )

    @reactive.effect
    @reactive.event(input.apply_analysis_grouping)
    def _apply_analysis_grouping():
        config = grouping_state.get()
        data = data_state.get()
        if config is None or data is None:
            return
        try:
            mode, column, rows = grouping_editor_rows()
            custom_groups = None
            if mode == "custom":
                custom_groups = {
                    str(row.sample_alias): str(getattr(input, f"analysis_group_value_{index}")() or "")
                    for index, row in rows.reset_index(drop=True).iterrows()
                }
            config = update_analysis_grouping_config(
                config,
                mode=mode,
                column=column,
                custom_groups=custom_groups,
            )
            data_state.set(apply_analysis_grouping(data, config["mapping"], config["label"]))
            spatial_retrieval_state.set(None)
            grouping_state.set(config)
            ui.modal_remove()
            ui.notification_show(analysis_grouping_summary(config), type="message")
            log_activity("Analysis grouping", "Completed", analysis_grouping_summary(config))
        except Exception as error:
            ui.notification_show(str(error), type="error", duration=None)
            log_activity("Analysis grouping", "Failed", str(error))

    @reactive.effect
    @reactive.event(input.reset_analysis_grouping)
    def _reset_analysis_grouping():
        config = grouping_state.get()
        data = data_state.get()
        if config is None or data is None:
            return
        column = "condition" if "condition" in config["columns"] else "sample_alias"
        config = update_analysis_grouping_config(config, mode="column", column=column)
        data_state.set(apply_analysis_grouping(data, config["mapping"], config["label"]))
        spatial_retrieval_state.set(None)
        grouping_state.set(config)
        ui.modal_remove()
        ui.notification_show(analysis_grouping_summary(config), type="message")
        log_activity("Analysis grouping", "Reset", analysis_grouping_summary(config))

    def get_data() -> AppData | None:
        return data_state.get()

    def get_spatial_retrieval() -> SpatialRetrieval | None:
        return spatial_retrieval_state.get()

    def tracked_pxl_proximity(operation: str, data: AppData, metadata: pd.DataFrame, **kwargs) -> pd.DataFrame:
        started = perf_counter()
        markers = kwargs.get("markers")
        marker_count = len(markers) if markers is not None else len(data.marker_options)
        log_activity(operation, "Started", f"Querying {len(metadata):,} cells and {marker_count:,} markers from PXL.")
        try:
            rows = load_pxl_proximity(data, metadata, **kwargs)
        except Exception as error:
            log_activity(operation, "Failed", str(error), perf_counter() - started)
            raise
        log_activity(operation, "Completed", f"Returned {len(rows):,} rows.", perf_counter() - started)
        return rows

    def tracked_sample_colocalization(
        data: AppData,
        metadata: pd.DataFrame,
        *,
        operation: str = "Differential colocalization",
        **kwargs,
    ) -> pd.DataFrame:
        started = perf_counter()
        log_activity(
            operation,
            "Started",
            f"Aggregating {len(metadata):,} cells inside PXL DuckDB.",
        )
        try:
            rows = sample_pxl_colocalization(data, metadata, **kwargs)
        except Exception as error:
            log_activity(operation, "Failed", str(error), perf_counter() - started)
            raise
        log_activity(
            operation,
            "Completed",
            f"Returned {len(rows):,} sample-pair rows.",
            perf_counter() - started,
        )
        return rows

    def tracked_component_layout(data: AppData, sample: str, component: str) -> pd.DataFrame:
        started = perf_counter()
        log_activity("3D layout", "Started", f"Reading component {component} from PXL.")
        try:
            rows = read_component_layout(data, sample, component)
        except Exception as error:
            log_activity("3D layout", "Failed", str(error), perf_counter() - started)
            raise
        log_activity("3D layout", "Completed", f"Returned {len(rows):,} nodes.", perf_counter() - started)
        return rows

    def register_downloads(output_id: str, producer):
        @output(id=f"{output_id}_png")
        @render.download_button(filename=f"{output_id.replace('_', '-')}.png")
        def _png():
            yield producer().to_image(format="png", width=1200, height=800, scale=2)

        @output(id=f"{output_id}_pdf")
        @render.download_button(filename=f"{output_id.replace('_', '-')}.pdf")
        def _pdf():
            yield producer().to_image(format="pdf", width=1200, height=800)

    def grid(frame: pd.DataFrame, *, max_rows: int = 2000):
        return render.DataGrid(frame.head(max_rows), filters=True, width="100%", height="520px")

    def filtered_qc_metadata() -> pd.DataFrame:
        data = get_data()
        if data is None:
            return pd.DataFrame()
        frame = data.metadata
        samples = selected(input.qc_sample_filter(), frame["sample_alias"])
        return frame[frame["sample_alias"].astype(str).isin(samples)].copy()

    @output
    @render.ui
    def qc_metric_row():
        data = get_data()
        if data is None:
            return metric_boxes([("Loaded Cells", "—"), ("Final Cells", "—"), ("Retained", "—"), ("Samples", "—")])
        samples = selected(input.qc_sample_filter(), data.metadata["sample_alias"])
        history = data.qc_filter_counts[data.qc_filter_counts["sample"].astype(str).isin(samples)]
        loaded = history.loc[history["step"].astype(str) == "00_loaded", "n_cells"].sum()
        final = len(data.metadata[data.metadata["sample_alias"].astype(str).isin(samples)])
        loaded = int(loaded) if loaded else final
        retained = final / loaded if loaded else math.nan
        return metric_boxes([
            ("Loaded Cells", f"{loaded:,}"), ("Final Cells", f"{final:,}"),
            ("Retained", f"{retained:.1%}" if np.isfinite(retained) else "—"), ("Samples", f"{len(samples):,}"),
        ])

    def fig_qc_filter():
        data = get_data()
        if data is None:
            return empty_figure()
        rows = data.qc_filter_counts.copy()
        samples = selected(input.qc_sample_filter(), rows["sample"])
        rows = rows[rows["sample"].astype(str).isin(samples)]
        value = "fraction_loaded" if input.qc_filter_y() == "fraction_loaded" else "n_cells"
        figure = px.bar(rows, x="step_label", y=value, color="sample", barmode="group", hover_data=["condition"])
        figure.update_yaxes(title="Fraction of loaded cells" if value == "fraction_loaded" else "Number of cells")
        return style_figure(figure)

    @output
    @render_plotly
    def qc_filter_plot():
        return fig_qc_filter()

    register_downloads("qc_filter_plot", fig_qc_filter)

    @output
    @render.data_frame
    def qc_filter_table():
        data = get_data()
        if data is None:
            return grid(pd.DataFrame())
        rows = data.qc_filter_counts.copy()
        samples = selected(input.qc_sample_filter(), rows["sample"])
        rows = rows[rows["sample"].astype(str).isin(samples)]
        return grid(rows)

    def fig_qc_rank():
        frame = filtered_qc_metadata()
        if frame.empty or "n_umi" not in frame:
            return empty_figure("No n_umi values are available.")
        rows = []
        for sample, chunk in frame.groupby("sample_alias", observed=True):
            chunk = chunk.sort_values("n_umi", ascending=False).copy()
            chunk["rank"] = np.arange(1, len(chunk) + 1)
            chunk["sample_alias"] = sample
            rows.append(chunk)
        plot_data = pd.concat(rows, ignore_index=True)
        figure = px.line(plot_data, x="rank", y="n_umi", color="sample_alias", hover_data=["component"])
        figure.add_hline(y=float(input.qc_n_umi_cutoff()), line_dash="dash", line_color=CORAL)
        figure.update_xaxes(type="log", title="Cell rank")
        figure.update_yaxes(type="log", title="n_umi")
        return style_figure(figure)

    @output
    @render_plotly
    def qc_molecule_rank_plot():
        return fig_qc_rank()

    register_downloads("qc_molecule_rank_plot", fig_qc_rank)

    def fig_qc_distribution():
        frame = filtered_qc_metadata()
        metric = input.qc_metric()
        if frame.empty or not metric or metric not in frame:
            return empty_figure("No numeric QC metric is available.")
        figure = px.violin(frame, x="sample_alias", y=metric, color="condition", box=True, points=False)
        if metric == "isotype_fraction":
            figure.add_hline(y=float(input.qc_isotype_cutoff()), line_dash="dash", line_color=CORAL)
        if metric in {"n_umi", "n_edges", "reads_in_component"}:
            figure.update_yaxes(type="log")
        return style_figure(figure)

    @output
    @render_plotly
    def qc_distribution_plot():
        return fig_qc_distribution()

    register_downloads("qc_distribution_plot", fig_qc_distribution)

    @output
    @render.data_frame
    def qc_metadata_table():
        return grid(filtered_qc_metadata())

    def filtered_metadata(condition_input, celltype_input) -> pd.DataFrame:
        data = get_data()
        if data is None:
            return pd.DataFrame()
        metadata = data.metadata
        conditions = selected(condition_input, metadata["condition"])
        celltypes = selected(celltype_input, metadata["celltype_manual"])
        return metadata[
            metadata["condition"].astype(str).isin(conditions)
            & metadata["celltype_manual"].astype(str).isin(celltypes)
        ].copy()

    def abundance_with_metadata() -> pd.DataFrame:
        data = get_data()
        if data is None:
            return pd.DataFrame()
        columns = ["component", "sample_alias", "condition", "celltype_manual"]
        return data.abundance.merge(data.metadata[columns], on="component", how="inner")

    @output
    @render.ui
    def abundance_metric_row():
        data = get_data()
        if data is None:
            return metric_boxes([("Cells", "—"), ("Markers", "—"), ("Groups", "—"), ("Cell Types", "—")])
        metadata = filtered_metadata(input.abundance_condition_filter(), input.abundance_celltype_filter())
        return metric_boxes([
            ("Cells", f"{len(metadata):,}"), ("Markers", f"{len(data.marker_options):,}"),
            ("Groups", f"{metadata['condition'].nunique():,}"), ("Cell Types", f"{metadata['celltype_manual'].nunique():,}"),
        ])

    def abundance_points() -> pd.DataFrame:
        data = get_data()
        if data is None:
            return pd.DataFrame()
        metadata = filtered_metadata(input.abundance_condition_filter(), input.abundance_celltype_filter())
        color_by = input.abundance_color_by() or "abundance"
        if color_by == "abundance":
            marker = input.abundance_marker()
            rows = data.abundance[data.abundance["marker"].astype(str) == str(marker)]
            metadata = metadata.merge(rows, on="component", how="inner")
        return metadata

    def fig_abundance_umap():
        data = get_data()
        plot_data = abundance_points()
        if data is None or plot_data.empty:
            return empty_figure("No cells match the abundance filters.")
        embedding = input.abundance_embedding()
        columns = embedding_columns(data.metadata)
        if embedding not in columns:
            return empty_figure("No two-dimensional embedding is available.")
        x, y = columns[embedding]
        color_by = input.abundance_color_by() or "abundance"
        color = "abundance" if color_by == "abundance" else color_by
        split = input.abundance_split_by() or None
        facet_columns = max(1, int(input.abundance_split_columns() or 2))
        figure = px.scatter(
            plot_data,
            x=x,
            y=y,
            color=color,
            facet_col=split,
            facet_col_wrap=facet_columns if split else 0,
            facet_row_spacing=0.1 if split else None,
            hover_data=[column for column in ("component", "sample_alias", "condition", "celltype_manual") if column in plot_data],
            color_continuous_scale=["#edf7f4", "#78aeb2", "#f0b45b", CORAL] if color == "abundance" else None,
            render_mode="webgl",
        )
        figure.update_traces(marker={"size": float(input.abundance_point_size() or 3), "opacity": 0.82})
        label_embedding_axes(figure)
        facet_count = plot_data[split].nunique(dropna=False) if split else 0
        height = split_umap_height(facet_count, facet_columns) if split else 540
        return style_figure(figure, height=height)

    @output
    @render_plotly
    def abundance_umap():
        return fig_abundance_umap()

    register_downloads("abundance_umap", fig_abundance_umap)

    @output
    @render.data_frame
    def abundance_table():
        rows = abundance_with_metadata()
        marker = input.abundance_marker()
        if rows.empty or not marker:
            return grid(pd.DataFrame())
        rows = rows[rows["marker"].astype(str) == str(marker)]
        return grid(summarize_numeric(rows, ["marker", "condition", "celltype_manual"], "abundance"))

    def abundance_distribution_data() -> pd.DataFrame:
        rows = abundance_with_metadata()
        marker = input.abundance_distribution_marker()
        return rows[rows["marker"].astype(str) == str(marker)].copy() if not rows.empty and marker else pd.DataFrame()

    def fig_abundance_distribution():
        rows = abundance_distribution_data()
        if rows.empty:
            return empty_figure("No abundance values are available for this marker.")
        figure = px.violin(
            rows,
            x="celltype_manual",
            y="abundance",
            color="condition",
            facet_col="sample_alias",
            facet_col_wrap=max(1, int(input.abundance_distribution_columns() or 3)),
            box=True,
            points="all" if input.abundance_distribution_show_jitter() else False,
            hover_data=["component"],
        )
        figure.update_traces(jitter=0.25, marker={"size": 2.5, "opacity": 0.45})
        return style_figure(figure, height=max(520, 300 * math.ceil(rows["sample_alias"].nunique() / max(1, int(input.abundance_distribution_columns() or 3)))))

    @output
    @render_plotly
    def abundance_distribution():
        return fig_abundance_distribution()

    register_downloads("abundance_distribution", fig_abundance_distribution)

    @output
    @render.data_frame
    def abundance_distribution_table():
        rows = abundance_distribution_data()
        return grid(summarize_numeric(rows, ["marker", "sample_alias", "condition", "celltype_manual"], "abundance"))

    def composition_data() -> pd.DataFrame:
        data = get_data()
        if data is None:
            return pd.DataFrame()
        counts = data.metadata.groupby(["condition", "celltype_manual"], observed=True).size().rename("n_cells").reset_index()
        counts["fraction"] = counts["n_cells"] / counts.groupby("condition", observed=True)["n_cells"].transform("sum")
        return counts

    def fig_abundance_composition():
        rows = composition_data()
        if rows.empty:
            return empty_figure("No cell annotations are available.")
        figure = px.bar(rows, x="condition", y="fraction", color="celltype_manual", hover_data=["n_cells"], barmode="stack")
        figure.update_yaxes(tickformat=".0%", title="Fraction of cells")
        return style_figure(figure)

    @output
    @render_plotly
    def abundance_composition():
        return fig_abundance_composition()

    register_downloads("abundance_composition", fig_abundance_composition)

    def annotation_heatmap_data() -> pd.DataFrame:
        rows = abundance_with_metadata()
        if rows.empty:
            return rows
        return rows.groupby(["celltype_manual", "marker"], observed=True)["abundance"].median().rename("median_abundance").reset_index()

    def fig_abundance_annotation_heatmap():
        rows = annotation_heatmap_data()
        if rows.empty:
            return empty_figure("No abundance values are available for the annotation heatmap.")
        matrix = rows.pivot(index="celltype_manual", columns="marker", values="median_abundance")
        figure = px.imshow(matrix, color_continuous_scale=["#edf7f4", "#78aeb2", "#f0b45b", CORAL], aspect="auto", labels={"color": "Median abundance"})
        return style_figure(figure, height=max(480, 36 * len(matrix)))

    @output
    @render_plotly
    def abundance_annotation_heatmap():
        return fig_abundance_annotation_heatmap()

    register_downloads("abundance_annotation_heatmap", fig_abundance_annotation_heatmap)

    @output
    @render.data_frame
    def abundance_composition_table():
        return grid(composition_data())

    def differential_result(prefix: str, rows: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
        if rows.empty:
            return pd.DataFrame()
        group_a = getattr(input, f"{prefix}_diff_group_a")()
        group_b = getattr(input, f"{prefix}_diff_group_b")()
        if not group_a or not group_b or group_a == group_b:
            return pd.DataFrame()
        return calculate_differential(
            rows,
            feature_cols=feature_cols,
            value_col="abundance" if prefix == "abundance" else "log2_ratio",
            group_a=group_a,
            group_b=group_b,
            celltypes=selected(getattr(input, f"{prefix}_diff_celltype_filter")(), rows["celltype_manual"]),
            stratify=bool(getattr(input, f"{prefix}_diff_stratify")()),
            min_observations=max(1, int(getattr(input, f"{prefix}_diff_min")() or 3)),
        )

    def abundance_diff_result() -> pd.DataFrame:
        return differential_result("abundance", abundance_with_metadata(), ["marker"])

    def volcano_figure(result: pd.DataFrame, label_col: str, fdr: float, effect: float, x_label: str):
        if result.empty:
            return empty_figure("Choose two different groups with enough observations.")
        rows = result.copy()
        rows["p_value"] = pd.to_numeric(rows["p_value"], errors="coerce").fillna(1.0)
        rows["p_adj"] = pd.to_numeric(rows["p_adj"], errors="coerce").fillna(1.0)
        rows["effect_size"] = pd.to_numeric(rows["effect_size"], errors="coerce").fillna(0.0)
        rows["minus_log10_fdr"] = -np.log10(rows["p_adj"].clip(lower=np.finfo(float).tiny))
        rows["threshold"] = np.where(
            (rows["p_adj"] <= fdr) & (rows["effect_size"].abs() >= effect), "Pass", "Below threshold"
        )
        figure = px.scatter(
            rows,
            x="effect_size",
            y="minus_log10_fdr",
            color="threshold",
            color_discrete_map={"Pass": CORAL, "Below threshold": "#8c9998"},
            hover_name=label_col,
            hover_data=["p_adj", "n_a", "n_b", "direction"],
            custom_data=[label_col],
        )
        figure.add_vline(x=effect, line_dash="dash", line_color="#7b8588")
        figure.add_vline(x=-effect, line_dash="dash", line_color="#7b8588")
        figure.add_hline(y=-math.log10(max(fdr, np.finfo(float).tiny)), line_dash="dash", line_color="#7b8588")
        figure.update_xaxes(title=x_label)
        figure.update_yaxes(title="−log10 adjusted p-value")
        return style_figure(figure)

    def diff_summary(result: pd.DataFrame, fdr: float, effect: float):
        tested = len(result)
        passing = int(((result.get("p_adj", pd.Series(dtype=float)) <= fdr) & (result.get("effect_size", pd.Series(dtype=float)).abs() >= effect)).sum())
        return metric_boxes([("Features tested", f"{tested:,}"), ("Threshold hits", f"{passing:,}"), ("FDR", f"{fdr:g}"), ("Minimum effect", f"{effect:g}")])

    @output
    @render.ui
    def abundance_diff_summary():
        return diff_summary(abundance_diff_result(), float(input.abundance_diff_fdr() or 0.05), float(input.abundance_diff_effect() or 0.25))

    def fig_abundance_diff_volcano():
        return volcano_figure(
            abundance_diff_result(), "marker", float(input.abundance_diff_fdr() or 0.05),
            float(input.abundance_diff_effect() or 0.25),
            f"Abundance effect: {input.abundance_diff_group_a()} minus {input.abundance_diff_group_b()} (reference)",
        )

    @output
    @render_plotly
    def abundance_diff_volcano():
        return clickable_volcano(
            fig_abundance_diff_volcano(), "abundance_diff_feature", session
        )

    register_downloads("abundance_diff_volcano", fig_abundance_diff_volcano)

    def fig_abundance_diff_detail():
        rows = abundance_with_metadata()
        marker = input.abundance_diff_feature()
        groups = [input.abundance_diff_group_a(), input.abundance_diff_group_b()]
        if rows.empty or not marker:
            return empty_figure("Select a detail marker.")
        rows = rows[(rows["marker"].astype(str) == str(marker)) & rows["condition"].isin(groups)]
        celltypes = selected(input.abundance_diff_celltype_filter(), rows["celltype_manual"])
        rows = rows[rows["celltype_manual"].astype(str).isin(celltypes)]
        figure = px.violin(
            rows, x="condition", y="abundance", color="condition",
            facet_col="celltype_manual" if input.abundance_diff_stratify() else None,
            box=True, points="all", hover_data=["component"],
        )
        figure.update_traces(jitter=0.25, marker={"size": 3, "opacity": 0.55})
        return style_figure(figure)

    @output
    @render_plotly
    def abundance_diff_detail():
        return fig_abundance_diff_detail()

    register_downloads("abundance_diff_detail", fig_abundance_diff_detail)

    @output
    @render.data_frame
    def abundance_diff_table():
        return grid(abundance_diff_result())

    def spatial_metadata(condition_input, celltype_input) -> pd.DataFrame:
        retrieval = get_spatial_retrieval()
        if retrieval is None:
            return pd.DataFrame()
        metadata = retrieval.metadata
        conditions = selected(condition_input, metadata["condition"])
        celltypes = selected(celltype_input, metadata["celltype_manual"])
        return metadata[
            metadata["condition"].astype(str).isin(conditions)
            & metadata["celltype_manual"].astype(str).isin(celltypes)
        ].copy()

    @reactive.calc
    def clustering_marker_rows() -> pd.DataFrame:
        data = get_data()
        retrieval = get_spatial_retrieval()
        marker = input.clustering_marker()
        metadata = spatial_metadata(input.clustering_condition_filter(), input.clustering_celltype_filter())
        if data is None or retrieval is None or metadata.empty or not marker or str(marker) not in retrieval.markers:
            return pd.DataFrame()
        rows = tracked_pxl_proximity(
            "Clustering marker",
            data,
            metadata,
            markers=[str(marker)],
            pair_type="self",
        )
        rows["marker"] = rows["marker_1"].astype(str)
        return rows

    def fig_clustering_observed():
        data = get_data()
        retrieval = get_spatial_retrieval()
        rows = clustering_marker_rows()
        if data is None or retrieval is None or rows.empty:
            return empty_figure("No PXL self-proximity values match the selected filters.")
        embedding = next(iter(embedding_columns(retrieval.metadata)), None)
        if not embedding:
            return empty_figure("No two-dimensional embedding is available.")
        x, y = embedding_columns(retrieval.metadata)[embedding]
        plot_data = rows.merge(retrieval.metadata[["component", x, y]], on="component", how="inner")
        figure = px.scatter(
            plot_data, x=x, y=y, color="log2_ratio", hover_data=["component", "condition", "celltype_manual"],
            color_continuous_scale="RdBu_r", color_continuous_midpoint=0, render_mode="webgl",
        )
        figure.update_traces(marker={"size": 4, "opacity": 0.82})
        return style_figure(figure)

    @output
    @render_plotly
    def clustering_plot():
        return fig_clustering_observed()

    register_downloads("clustering_plot", fig_clustering_observed)

    @output
    @render.data_frame
    def clustering_table():
        return grid(summarize_numeric(clustering_marker_rows(), ["marker", "condition", "celltype_manual"], "log2_ratio"))

    def fig_clustering_per_marker():
        rows = clustering_marker_rows()
        if rows.empty:
            return empty_figure("No PXL self-proximity values match the selected filters.")
        figure = px.violin(rows, x="celltype_manual", y="log2_ratio", color="condition", box=True, points="all", hover_data=["component"])
        figure.add_hline(y=0, line_color="#7b8588")
        figure.update_traces(jitter=0.25, marker={"size": 3, "opacity": 0.5})
        return style_figure(figure)

    @output
    @render_plotly
    def clustering_per_marker():
        return fig_clustering_per_marker()

    register_downloads("clustering_per_marker", fig_clustering_per_marker)

    @output
    @render.data_frame
    def clustering_per_marker_table():
        return grid(summarize_numeric(clustering_marker_rows(), ["marker", "condition", "celltype_manual"], "log2_ratio"))

    @reactive.calc
    def clustering_heatmap_data() -> pd.DataFrame:
        data = get_data()
        retrieval = get_spatial_retrieval()
        metadata = spatial_metadata(
            input.clustering_heatmap_condition_filter(),
            input.clustering_heatmap_celltype_filter(),
        )
        if data is None or retrieval is None or metadata.empty:
            return pd.DataFrame()
        rows = tracked_pxl_proximity(
            "Clustering heatmap",
            data,
            metadata,
            markers=retrieval.markers,
            pair_type="self",
        )
        if rows.empty:
            return rows
        rows["marker"] = rows["marker_1"].astype(str)
        summary = rows.groupby(["marker", "condition", "celltype_manual"], observed=True)["log2_ratio"].mean().rename("mean_log2_ratio").reset_index()
        count = max(2, int(input.clustering_heatmap_marker_count() or 20))
        ranking = summary.groupby("marker", observed=True)["mean_log2_ratio"].std().fillna(0).sort_values(ascending=False)
        return summary[summary["marker"].isin(ranking.head(count).index)]

    def fig_clustering_heatmap():
        rows = clustering_heatmap_data()
        if rows.empty:
            return empty_figure("No PXL clustering summary matches the selected filters.")
        rows["population"] = rows["condition"].astype(str) + " · " + rows["celltype_manual"].astype(str)
        matrix = rows.pivot(index="population", columns="marker", values="mean_log2_ratio")
        figure = px.imshow(matrix, color_continuous_scale="RdBu_r", color_continuous_midpoint=0, aspect="auto", labels={"color": "Mean log2 ratio"})
        return style_figure(figure, height=max(520, 26 * len(matrix)))

    @output
    @render_plotly
    def clustering_summary_heatmap():
        return fig_clustering_heatmap()

    register_downloads("clustering_summary_heatmap", fig_clustering_heatmap)

    @output
    @render.data_frame
    def clustering_summary_table():
        return grid(clustering_heatmap_data())

    @reactive.calc
    def clustering_diff_rows() -> pd.DataFrame:
        data = get_data()
        retrieval = get_spatial_retrieval()
        if data is None or retrieval is None:
            return pd.DataFrame()
        groups = [input.clustering_diff_group_a(), input.clustering_diff_group_b()]
        metadata = retrieval.metadata[retrieval.metadata["condition"].isin(groups)].copy()
        celltypes = selected(input.clustering_diff_celltype_filter(), metadata["celltype_manual"])
        metadata = metadata[metadata["celltype_manual"].astype(str).isin(celltypes)]
        if metadata.empty:
            return pd.DataFrame()
        rows = tracked_pxl_proximity(
            "Differential clustering",
            data,
            metadata,
            markers=retrieval.markers,
            pair_type="self",
        )
        rows["marker"] = rows["marker_1"].astype(str)
        return rows

    @reactive.calc
    def clustering_diff_result():
        return differential_result("clustering", clustering_diff_rows(), ["marker"])

    @output
    @render.ui
    def clustering_diff_summary():
        return diff_summary(clustering_diff_result(), float(input.clustering_diff_fdr() or 0.05), float(input.clustering_diff_effect() or 0.25))

    def fig_clustering_diff_volcano():
        return volcano_figure(
            clustering_diff_result(), "marker", float(input.clustering_diff_fdr() or 0.05),
            float(input.clustering_diff_effect() or 0.25),
            f"Clustering effect: {input.clustering_diff_group_a()} minus {input.clustering_diff_group_b()} (reference)",
        )

    @output
    @render_plotly
    def clustering_diff_volcano():
        return clickable_volcano(
            fig_clustering_diff_volcano(), "clustering_diff_feature", session
        )

    register_downloads("clustering_diff_volcano", fig_clustering_diff_volcano)

    def fig_clustering_diff_detail():
        marker = input.clustering_diff_feature()
        rows = clustering_diff_rows()
        if rows.empty:
            return empty_figure("No PXL self-proximity values match this contrast.")
        rows = rows[
            (rows["marker"].astype(str) == str(marker))
            & rows["condition"].isin([input.clustering_diff_group_a(), input.clustering_diff_group_b()])
        ]
        celltypes = selected(input.clustering_diff_celltype_filter(), rows["celltype_manual"])
        rows = rows[rows["celltype_manual"].astype(str).isin(celltypes)]
        if rows.empty:
            return empty_figure("No clustering values are available for this marker.")
        figure = px.violin(
            rows, x="condition", y="log2_ratio", color="condition",
            facet_col="celltype_manual" if input.clustering_diff_stratify() else None,
            box=True, points="all", hover_data=["component"],
        )
        figure.add_hline(y=0, line_color="#7b8588")
        figure.update_traces(jitter=0.25, marker={"size": 3, "opacity": 0.55})
        return style_figure(figure)

    @output
    @render_plotly
    def clustering_diff_detail():
        return fig_clustering_diff_detail()

    register_downloads("clustering_diff_detail", fig_clustering_diff_detail)

    @output
    @render.data_frame
    def clustering_diff_table():
        return grid(clustering_diff_result())

    @reactive.effect
    def _apply_report_preset():
        if input.coloc_preset() != "report":
            return
        ui.update_select("coloc_view", selected="focused")
        ui.update_select("coloc_marker_mode", selected="profile")
        ui.update_numeric("coloc_top_pairs", value=60)
        ui.update_select("coloc_mean_type", selected="population")
        ui.update_select("coloc_ordering", selected="ward")
        ui.update_numeric("coloc_legend_min", value=-1)
        ui.update_numeric("coloc_legend_max", value=1)
        ui.update_checkbox("coloc_pixelator_filter", value=True)
        ui.update_numeric("coloc_min_fraction", value=0.001)
        ui.update_numeric("coloc_min_count", value=0)
        ui.update_numeric("coloc_min_cells", value=1)

    def colocalization_metadata() -> pd.DataFrame:
        retrieval = get_spatial_retrieval()
        if retrieval is None:
            return pd.DataFrame()
        metadata = retrieval.metadata.copy()
        conditions = selected(input.coloc_condition_filter(), metadata["condition"])
        celltypes = selected(input.coloc_celltype_filter(), metadata["celltype_manual"])
        metadata = metadata[
            metadata["condition"].astype(str).isin(conditions)
            & metadata["celltype_manual"].astype(str).isin(celltypes)
        ]
        if input.coloc_scope() == "celltype" and input.coloc_celltype_focus():
            metadata = metadata[metadata["celltype_manual"].astype(str) == str(input.coloc_celltype_focus())]
        return metadata

    def apply_colocalization_filters(scores: pd.DataFrame) -> pd.DataFrame:
        if input.coloc_pixelator_filter():
            try:
                return filter_pixelator_proximity(
                    scores,
                    min_marker_fraction=float(input.coloc_min_fraction() or 0),
                    min_marker_count=float(input.coloc_min_count() or 0),
                    min_cells=int(input.coloc_min_cells() or 1),
                )
            except ValueError as error:
                ui.notification_show(str(error), type="warning")
                return pd.DataFrame()
        return scores

    def heatmap_markers(metadata: pd.DataFrame) -> list[str]:
        data = get_data()
        retrieval = get_spatial_retrieval()
        if data is None or retrieval is None or metadata.empty:
            return []
        manual_value = input.coloc_markers()
        manual = (
            [manual_value] if isinstance(manual_value, str) else list(manual_value or [])
        ) if input.coloc_marker_mode() == "manual" else None
        return select_colocalization_heatmap_markers(
            data.abundance,
            metadata,
            retrieval.markers,
            n_markers=max(2, int(input.coloc_top_markers() or 40)),
            plot_markers=manual,
        )

    def ordered_coloc_markers(summary: pd.DataFrame, markers: list[str], group_col: str) -> list[str]:
        if len(markers) < 3 or summary.empty:
            return markers
        requested_group = (
            input.coloc_reference()
            if input.coloc_view() == "compare"
            else input.coloc_focus_group()
        )
        groups = set(summary[group_col].astype(str))
        group = requested_group if requested_group in groups else str(summary[group_col].iloc[0])
        rows = summary[summary[group_col].astype(str) == group]
        matrix = rows.pivot(index="marker_1", columns="marker_2", values="mean_log2_ratio").reindex(index=markers, columns=markers).fillna(0)
        try:
            distances = pdist(matrix.to_numpy(), metric="euclidean")
            if not distances.size or np.allclose(distances, 0):
                return markers
            tree = linkage(
                distances,
                method=input.coloc_ordering() or "ward",
                optimal_ordering=True,
            )
            return [markers[index] for index in leaves_list(tree)]
        except Exception:
            return markers

    @reactive.calc
    def colocalization_heatmap_base() -> tuple[pd.DataFrame, str, list[str]]:
        data = get_data()
        retrieval = get_spatial_retrieval()
        metadata = colocalization_metadata()
        if data is None or retrieval is None or metadata.empty:
            return pd.DataFrame(), "condition", []
        scope = input.coloc_scope() or "condition"
        group_col = "sample_alias" if scope == "sample_alias" else "condition"
        mode = input.coloc_marker_mode() or "profile"
        if mode == "profile":
            available_markers = list(retrieval.markers)
            samples = tracked_sample_colocalization(
                data,
                metadata,
                operation="Colocalization heatmap",
                markers=available_markers,
                mean_type="population",
                pair_type="all",
                min_marker_fraction=(
                    float(input.coloc_min_fraction() or 0)
                    if input.coloc_pixelator_filter()
                    else 0
                ),
                min_marker_count=(
                    float(input.coloc_min_count() or 0)
                    if input.coloc_pixelator_filter()
                    else 0
                ),
            )
            summary = summarize_sample_colocalization(
                samples,
                group_col=group_col,
                markers=available_markers,
                mean_type=input.coloc_mean_type() or "population",
            )
            minimum_cells = int(input.coloc_min_cells() or 1)
            if minimum_cells > 1 and not summary.empty:
                below_minimum = summary["n_detected"] < minimum_cells
                summary.loc[
                    below_minimum,
                    ["sum_log2_ratio", "detected_mean", "mean_log2_ratio", "n_detected", "pct_detected"],
                ] = 0
            markers = select_proximity_profile_markers(
                summary,
                n_pairs=max(1, int(input.coloc_top_pairs() or 60)),
            )
            summary = summary[
                summary["marker_1"].isin(markers) & summary["marker_2"].isin(markers)
            ].copy()
        else:
            markers = heatmap_markers(metadata)
            scores = tracked_pxl_proximity(
                "Colocalization heatmap",
                data,
                metadata,
                markers=markers,
                pair_type="all",
                add_marker_counts=bool(input.coloc_pixelator_filter()),
            )
            scores = apply_colocalization_filters(scores)
            if scores.empty:
                return pd.DataFrame(), group_col, markers
            summary = summarize_spatial(
                scores,
                metadata,
                group_col=group_col,
                markers=markers,
                mean_type=input.coloc_mean_type() or "population",
            )
        if summary.empty:
            return summary, group_col, markers
        return summary, group_col, markers

    @reactive.calc
    def colocalization_heatmap_data() -> tuple[pd.DataFrame, str, list[str]]:
        summary, group_col, markers = colocalization_heatmap_base()
        if summary.empty:
            return summary, group_col, markers
        available_groups = list(summary[group_col].dropna().astype(str).unique())
        if input.coloc_view() == "compare":
            groups = selected(input.coloc_compare_groups(), available_groups)
        else:
            groups = [str(input.coloc_focus_group() or available_groups[0])]
        summary = summary[summary[group_col].astype(str).isin(groups)].copy()
        markers = ordered_coloc_markers(summary, markers, group_col)
        summary["marker_1"] = pd.Categorical(summary["marker_1"], markers, ordered=True)
        summary["marker_2"] = pd.Categorical(summary["marker_2"], list(reversed(markers)), ordered=True)
        return summary, group_col, markers

    @reactive.effect
    def _refresh_colocalization_detail_pairs():
        _, _, markers = colocalization_heatmap_base()
        pairs = [f"{first} / {second}" for first, second in combinations(markers, 2)]
        ui.update_selectize(
            "coloc_detail_pair",
            choices=pairs,
            selected=pairs[0] if pairs else None,
            server=True,
        )

    @output
    @render.ui
    def coloc_notice():
        data = get_data()
        retrieval = get_spatial_retrieval()
        if data is not None and retrieval is None:
            return ui.div("Retrieve data in Spatial Metrics > Retrieve Data first.", class_="alert alert-warning")
        summary, group_col, markers = colocalization_heatmap_data()
        if summary.empty:
            return ui.div("No spatial metric rows match the selected settings.", class_="alert alert-warning")
        notices = [f"Showing {len(markers)} markers across {summary[group_col].nunique()} {group_col.replace('_', ' ')} group(s)."]
        if input.coloc_marker_mode() == "profile":
            notices.append(
                f"Markers come from the {int(input.coloc_top_pairs() or 60)} strongest pairs "
                "detected in more than half of cells."
            )
        if input.coloc_view() == "compare" and len(markers) > 20:
            notices.append("Large comparison matrices are stacked vertically so labels remain readable.")
        elif input.coloc_view() == "compare" and len(markers) > 15:
            notices.append("For clearer comparison labels, use 15 or fewer markers.")
        return ui.div(" ".join(notices), class_="alert alert-info py-2")

    def fig_coloc_heatmap():
        summary, group_col, markers = colocalization_heatmap_data()
        if summary.empty:
            return empty_figure("No PXL colocalization scores match the selected settings.")
        legend_min = float(input.coloc_legend_min() or -1)
        legend_max = float(input.coloc_legend_max() or 1)
        if legend_min >= legend_max:
            legend_min, legend_max = -1, 1
        figure = px.scatter(
            summary,
            x="marker_1",
            y="marker_2",
            color="mean_log2_ratio",
            size="pct_detected",
            facet_col=group_col,
            facet_col_wrap=(
                1 if len(markers) > 20 else 2
            ) if input.coloc_view() == "compare" else 1,
            color_continuous_scale="RdBu_r",
            range_color=(legend_min, legend_max),
            size_max=18,
            hover_data={
                "mean_log2_ratio": ":.3f",
                "n_detected": True,
                "n_total": True,
                "pct_detected": ":.1%",
            },
            category_orders={"marker_1": markers, "marker_2": list(reversed(markers))},
            labels={"mean_log2_ratio": "Mean log2 ratio", "pct_detected": "Detected fraction"},
        )
        figure.update_traces(marker={"sizemin": 0, "line": {"width": 0.65, "color": "#222222"}})
        panels = max(1, summary[group_col].nunique())
        facet_columns = 1 if len(markers) > 20 or input.coloc_view() != "compare" else 2
        panel_rows = math.ceil(panels / facet_columns)
        figure = style_figure(
            figure,
            height=max(720, 26 * len(markers) * panel_rows + 240),
        )
        figure.update_xaxes(
            side="top",
            tickangle=-45,
            showgrid=True,
            gridcolor="#dfe7eb",
            zeroline=False,
            title_text=None,
            tickfont={"size": 10},
        )
        figure.update_yaxes(
            showgrid=True,
            gridcolor="#dfe7eb",
            zeroline=False,
            title_text=None,
            tickfont={"size": 11},
        )
        figure.for_each_annotation(
            lambda annotation: annotation.update(
                text="" if panels == 1 else annotation.text.split("=", 1)[-1],
                font={"size": 15},
            )
        )
        groups = summary[group_col].dropna().astype(str).drop_duplicates().tolist()
        celltypes = colocalization_metadata()["celltype_manual"].dropna().astype(str).unique()
        population = str(celltypes[0]) if len(celltypes) == 1 else "selected cells"
        title = (
            f"Colocalization in {population} ({groups[0]})"
            if len(groups) == 1
            else f"Colocalization in {population}"
        )
        figure.update_layout(
            title={
                "text": title,
                "x": 0.01,
                "xanchor": "left",
            },
            margin={"l": 150, "r": 150, "t": 180, "b": 55},
            legend={
                "title": {"text": "Detected fraction"},
                "orientation": "h",
                "yanchor": "bottom",
                "y": 1.08,
                "xanchor": "right",
                "x": 0.84,
            },
            coloraxis_colorbar={
                "title": {"text": "Mean<br>log2 ratio"},
                "thickness": 22,
                "len": 0.72,
            },
        )
        return figure

    @output
    @render_plotly
    def coloc_heatmap():
        return fig_coloc_heatmap()

    register_downloads("coloc_heatmap", fig_coloc_heatmap)

    @output
    @render.data_frame
    def coloc_table():
        summary, _, _ = colocalization_heatmap_data()
        return grid(summary)

    def pair_markers(pair: str | None) -> tuple[str, str] | None:
        if not pair or " / " not in pair:
            return None
        return tuple(pair.split(" / ", 1))

    @reactive.calc
    def pair_detail_data() -> pd.DataFrame:
        data = get_data()
        retrieval = get_spatial_retrieval()
        metadata = colocalization_metadata()
        pair = pair_markers(input.coloc_detail_pair())
        if (
            data is None
            or retrieval is None
            or metadata.empty
            or pair is None
            or not set(pair).issubset(retrieval.markers)
        ):
            return pd.DataFrame()
        marker_1, marker_2 = pair
        scores = tracked_pxl_proximity(
            "Colocalization pair",
            data,
            metadata,
            markers=[marker_1, marker_2],
            pair_type="nonself",
            add_marker_counts=bool(input.coloc_pixelator_filter()),
        )
        scores = apply_colocalization_filters(scores)
        rows = scores[
            ((scores["marker_1"].astype(str) == marker_1) & (scores["marker_2"].astype(str) == marker_2))
            | ((scores["marker_1"].astype(str) == marker_2) & (scores["marker_2"].astype(str) == marker_1))
        ].copy()
        totals = metadata.groupby(["sample_alias", "condition", "celltype_manual"], observed=True)["component"].nunique().rename("n_total").reset_index()
        detected = rows.groupby(["sample_alias", "condition", "celltype_manual"], observed=True).agg(
            sum_log2_ratio=("log2_ratio", "sum"), detected_mean=("log2_ratio", "mean"), n_detected=("component", "nunique")
        ).reset_index()
        detail = totals.merge(detected, on=["sample_alias", "condition", "celltype_manual"], how="left")
        detail[["sum_log2_ratio", "n_detected"]] = detail[["sum_log2_ratio", "n_detected"]].fillna(0)
        detail["pct_detected"] = detail["n_detected"] / detail["n_total"]
        detail["mean_log2_ratio"] = np.where(
            input.coloc_mean_type() == "detected", detail["detected_mean"], detail["sum_log2_ratio"] / detail["n_total"]
        )
        detail["marker_pair"] = input.coloc_detail_pair()
        return detail

    @output
    @render.ui
    def coloc_pair_metrics():
        rows = pair_detail_data()
        if rows.empty:
            return metric_boxes([("Marker pair", "—"), ("Samples", "—"), ("Detected cells", "—"), ("Detected fraction", "—")])
        return metric_boxes([
            ("Marker pair", str(input.coloc_detail_pair())), ("Samples", f"{rows['sample_alias'].nunique():,}"),
            ("Detected cells", f"{int(rows['n_detected'].sum()):,} / {int(rows['n_total'].sum()):,}"),
            ("Detected fraction", f"{rows['n_detected'].sum() / rows['n_total'].sum():.1%}"),
        ])

    def fig_coloc_pair_detail():
        rows = pair_detail_data()
        if rows.empty:
            return empty_figure("No sample or cell-population detail is available for this pair.")
        figure = px.scatter(
            rows, x="sample_alias", y="celltype_manual", color="mean_log2_ratio", size="pct_detected",
            facet_col="condition", color_continuous_scale="RdBu_r", color_continuous_midpoint=0,
            hover_data=["n_detected", "n_total", "pct_detected"],
        )
        figure.update_traces(marker={"sizemin": 2, "line": {"width": 0.5, "color": "#556264"}})
        return style_figure(figure, height=max(500, 34 * rows["celltype_manual"].nunique()))

    @output
    @render_plotly
    def coloc_pair_detail():
        return fig_coloc_pair_detail()

    register_downloads("coloc_pair_detail", fig_coloc_pair_detail)

    @output
    @render.data_frame
    def coloc_pair_table():
        return grid(pair_detail_data())

    @reactive.calc
    def coloc_sample_data() -> pd.DataFrame:
        data = get_data()
        retrieval = get_spatial_retrieval()
        if data is None or retrieval is None:
            return pd.DataFrame()
        groups = [input.coloc_diff_group_a(), input.coloc_diff_group_b()]
        metadata = retrieval.metadata[retrieval.metadata["condition"].isin(groups)].copy()
        celltypes = selected(input.coloc_diff_celltype_filter(), metadata["celltype_manual"])
        metadata = metadata[metadata["celltype_manual"].astype(str).isin(celltypes)]
        anchor = (
            str(input.coloc_diff_anchor())
            if input.coloc_diff_pair_scope() == "anchor" and input.coloc_diff_anchor()
            else None
        )
        summary = tracked_sample_colocalization(
            data,
            metadata,
            markers=retrieval.markers,
            mean_type=input.coloc_diff_mean() or "population",
            anchor=anchor,
        )
        if summary.empty:
            return summary
        summary["celltype_manual"] = "Pooled cell types"
        summary["sample_value"] = summary["mean_log2_ratio"]
        return summary

    @reactive.calc
    def coloc_diff_result() -> pd.DataFrame:
        rows = coloc_sample_data()
        if rows.empty or input.coloc_diff_group_a() == input.coloc_diff_group_b():
            return pd.DataFrame()
        result = calculate_differential(
            rows,
            feature_cols=["marker_pair"],
            value_col="sample_value",
            group_a=input.coloc_diff_group_a(),
            group_b=input.coloc_diff_group_b(),
            min_observations=max(1, int(input.coloc_diff_min_samples() or 2)),
        )
        return result

    @output
    @render.ui
    def coloc_diff_method():
        rows = coloc_sample_data()
        if rows.empty:
            return ui.div("Assign matching PXL files for sample-aware differential colocalization.", class_="alert alert-warning")
        counts = rows.groupby("condition", observed=True)["sample_alias"].nunique()
        a, b = input.coloc_diff_group_a(), input.coloc_diff_group_b()
        minimum = int(input.coloc_diff_min_samples() or 2)
        replicated = counts.get(a, 0) >= minimum and counts.get(b, 0) >= minimum
        message = (
            f"Effects compare median sample-level {input.coloc_diff_mean()} means. "
            f"{a}: {counts.get(a, 0)} sample(s); {b}: {counts.get(b, 0)} sample(s)."
        )
        if not replicated:
            message += f" At least {minimum} samples per group are required for p-values and FDR."
        return ui.div(ui.strong("Sample-aware analysis. "), message, class_=f"alert {'alert-info' if replicated else 'alert-warning'}")

    @output
    @render.ui
    def coloc_diff_summary():
        return diff_summary(coloc_diff_result(), float(input.coloc_diff_fdr() or 0.05), float(input.coloc_diff_effect() or 0.25))

    def fig_coloc_diff_volcano():
        return volcano_figure(
            coloc_diff_result(), "marker_pair", float(input.coloc_diff_fdr() or 0.05),
            float(input.coloc_diff_effect() or 0.25),
            f"Median sample effect: {input.coloc_diff_group_a()} minus {input.coloc_diff_group_b()} (reference)",
        )

    @output
    @render_plotly
    def coloc_diff_volcano():
        return clickable_volcano(
            fig_coloc_diff_volcano(), "coloc_diff_pair", session
        )

    register_downloads("coloc_diff_volcano", fig_coloc_diff_volcano)

    def fig_coloc_diff_detail():
        rows = coloc_sample_data()
        pair = input.coloc_diff_pair()
        if rows.empty or not pair:
            return empty_figure("Select a detail pair.")
        rows = rows[
            (rows["marker_pair"].astype(str) == str(pair))
            & rows["condition"].isin([input.coloc_diff_group_a(), input.coloc_diff_group_b()])
        ]
        if rows.empty:
            return empty_figure("No sample-level values are available for this pair.")
        figure = px.box(rows, x="condition", y="sample_value", color="condition", points="all", hover_data=["sample_alias", "n_detected", "n_total"])
        figure.add_hline(y=0, line_color="#7b8588")
        figure.update_yaxes(title=f"{pair} {input.coloc_diff_mean()} mean")
        return style_figure(figure)

    @output
    @render_plotly
    def coloc_diff_detail():
        return fig_coloc_diff_detail()

    register_downloads("coloc_diff_detail", fig_coloc_diff_detail)

    @output
    @render.data_frame
    def coloc_diff_table():
        return grid(coloc_diff_result())

    def fig_coloc_3d():
        data = get_data()
        retrieval = get_spatial_retrieval()
        sample = input.coloc_3d_sample()
        component = input.coloc_3d_component()
        available_components = set(retrieval.metadata["component"].astype(str)) if retrieval is not None else set()
        if data is None or retrieval is None or not sample or not component or str(component) not in available_components:
            return empty_figure("Retrieve spatial data and select one of its components.")
        try:
            nodes = tracked_component_layout(data, str(sample), str(component))
        except Exception as error:
            return empty_figure(str(error))
        if "marker" not in nodes:
            nodes["marker"] = "unlabeled"
        highlights = [
            marker for marker in selected(input.coloc_3d_markers(), nodes["marker"])
            if marker in retrieval.markers
        ]
        foreground = nodes[nodes["marker"].astype(str).isin(highlights)].copy()
        background = nodes[~nodes["marker"].astype(str).isin(highlights)].copy()
        maximum = max(0, int(input.coloc_3d_max_background() or 7000))
        if len(background) > maximum:
            background = background.sample(maximum, random_state=42)
        background["marker_group"] = "Other"
        foreground["marker_group"] = foreground["marker"].astype(str)
        nodes = pd.concat([background, foreground], ignore_index=True)
        for axis in ("x", "y", "z"):
            values = pd.to_numeric(nodes[axis], errors="coerce")
            centered = values - values.median()
            nodes[axis] = centered / (centered.abs().max() or 1)
        figure = px.scatter_3d(
            nodes, x="x", y="y", z="z", color="marker_group",
            hover_data=[column for column in ("marker", "umi") if column in nodes],
            color_discrete_map={"Other": "#c8cecd"}, opacity=0.78,
        )
        figure.update_traces(marker={"size": 2})
        figure.update_layout(scene=dict(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False))
        return style_figure(figure, height=620)

    @output
    @render_plotly
    def coloc_3d_layout():
        return fig_coloc_3d()

    register_downloads("coloc_3d_layout", fig_coloc_3d)

    @output
    @render.data_frame
    def coloc_3d_table():
        retrieval = get_spatial_retrieval()
        if retrieval is None or not input.coloc_3d_component():
            return grid(pd.DataFrame())
        columns = [
            column for column in ("sample", "sample_alias", "condition", "celltype_manual", "component", "n_umi", "n_edges")
            if column in retrieval.metadata
        ]
        return grid(
            retrieval.metadata.loc[
                retrieval.metadata["component"].astype(str) == str(input.coloc_3d_component()),
                columns,
            ]
        )

    def patch_table(name: str) -> pd.DataFrame:
        data = get_data()
        if data is None or data.patch.get(name) is None:
            return pd.DataFrame()
        return data.patch[name].copy()

    @output
    @render.ui
    def patch_metric_row():
        plan = patch_table("run_plan")
        if plan.empty:
            return metric_boxes([("Patch Detection", "Unavailable"), ("Cells Selected", "—"), ("Patch Markers", "—"), ("Receiver Markers", "—")])
        row = plan.iloc[0]
        cells = row.get("n_cart_cells_selected", row.get("n_cd8t_cells_selected", 0))
        return metric_boxes([
            ("Patch Detection", "Run" if bool(row.get("run_patch_detection", False)) else "Skipped"),
            ("Cells Selected", f"{int(cells):,}"), ("Patch Markers", f"{int(row.get('n_patch_markers', 0)):,}"),
            ("Receiver Markers", f"{int(row.get('n_receiver_markers', 0)):,}"),
        ])

    def filtered_marker_unmixing() -> pd.DataFrame:
        rows = patch_table("marker_unmixing")
        if rows.empty:
            return rows
        labels = selected(input.patch_label_filter(), rows["label"]) if "label" in rows else []
        return rows[rows["label"].astype(str).isin(labels)] if labels else rows

    def fig_patch_marker():
        rows = filtered_marker_unmixing()
        if rows.empty or not {"receiver_freq", "target_freq"}.issubset(rows):
            return empty_figure('Patch tables are not stored in this H5AD under uns["proxiome"]["patch"].')
        figure = px.scatter(rows, x="receiver_freq", y="target_freq", color="label" if "label" in rows else None, hover_name="marker" if "marker" in rows else None)
        maximum = max(rows["receiver_freq"].max(), rows["target_freq"].max())
        figure.add_shape(type="line", x0=0, y0=0, x1=maximum, y1=maximum, line={"dash": "dash", "color": "#7b8588"})
        return style_figure(figure)

    @output
    @render_plotly
    def patch_marker_plot():
        return fig_patch_marker()

    register_downloads("patch_marker_plot", fig_patch_marker)

    @output
    @render.data_frame
    def patch_marker_table():
        return grid(filtered_marker_unmixing())

    def fig_patch_raji():
        rows = patch_table("raji_marker_proximity")
        if rows.empty or not {"raji_marker_count", "log2_ratio"}.issubset(rows):
            return empty_figure("No Raji joint-proximity table is stored in this H5AD.")
        color = "celltype_condition" if "celltype_condition" in rows else None
        figure = px.scatter(rows, x="raji_marker_count", y="log2_ratio", color=color, hover_data=[column for column in ("component", "join_count") if column in rows], opacity=0.65)
        figure.add_hline(y=0, line_color="#7b8588")
        figure.update_xaxes(type="log")
        return style_figure(figure)

    @output
    @render_plotly
    def patch_raji_plot():
        return fig_patch_raji()

    register_downloads("patch_raji_plot", fig_patch_raji)

    @output
    @render.data_frame
    def patch_raji_table():
        return grid(patch_table("raji_marker_abundance"))

    @output
    @render.data_frame
    def patch_burden_table():
        return grid(patch_table("patch_burden"))


app = App(app_ui, server, static_assets=APP_DIR / "www")
