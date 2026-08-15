"""ProxiomeVis: the Python Shiny implementation."""

from __future__ import annotations

import asyncio
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.cluster.hierarchy import leaves_list, linkage
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_plotly

from proxiome import (
    AppData,
    apply_analysis_grouping,
    calculate_differential,
    default_pxl_spec,
    filter_pixelator_proximity,
    load_pxl_data,
    mapping_for_column,
    parse_group_mapping,
    read_component_layout,
    resolve_pxl_paths,
    sample_colocalization,
    sample_level_columns,
    summarize_numeric,
    summarize_spatial,
)


APP_DIR = Path(__file__).resolve().parent
TEAL = "#176d73"
CORAL = "#c7503e"
PLOT_COLORS = ["#176d73", "#c7503e", "#c58a20", "#62879a", "#7b6aa2", "#56875d"]


def selectize(id: str, label: str, *, multiple: bool = False):
    return ui.input_selectize(id, label, [], multiple=multiple, remove_button=multiple)


def plot_pane(output_id: str, *, height: str = "520px"):
    return ui.div(
        ui.div(
            ui.download_button(f"{output_id}_png", "PNG", class_="btn-sm btn-outline-secondary"),
            ui.download_button(f"{output_id}_pdf", "PDF", class_="btn-sm btn-outline-secondary"),
            class_="plot-pane-controls",
        ),
        output_widget(output_id, height=height),
        class_="plot-pane plot-pane-standard",
    )


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
                ui.input_select("qc_metadata_source", "Metadata", {"origin": "Original cells", "filtered": "Filtered cells"}),
            ),
            ui.accordion_panel(
                "Cutoffs",
                ui.input_numeric("qc_n_umi_cutoff", "n_umi cutoff", 10000, min=0, step=500),
                ui.input_numeric("qc_isotype_cutoff", "Isotype fraction cutoff", 0.001, min=0, max=1, step=0.0005),
            ),
            ui.accordion_panel(
                "Display",
                ui.input_select("qc_filter_y", "Filter count y-axis", {"count": "Number of cells", "fraction_loaded": "Fraction of loaded cells"}),
                ui.input_checkbox("qc_filter_include_total", "Include TOTAL rows in table", True),
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
                    selectize("abundance_marker", "Marker"),
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
                ui.nav_panel("Observed", ui.output_ui("abundance_metric_row"), plot_pane("abundance_umap"), table_pane("abundance_table")),
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


def differential_controls(prefix: str, detail_label: str):
    return ui.accordion(
        ui.accordion_panel(
            "Contrast",
            ui.input_select(f"{prefix}_diff_group_a", "Group A", []),
            ui.input_select(f"{prefix}_diff_group_b", "Group B (reference)", []),
            selectize(f"{prefix}_diff_celltype_filter", "Cell type", multiple=True),
            ui.input_checkbox(f"{prefix}_diff_stratify", "Stratify by cell type", False),
            ui.input_action_button(f"{prefix}_run_differential", "Run differential analysis", class_="btn-primary w-100"),
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
        ui.panel_conditional("input.clustering_mode === 'Differential'", differential_controls("clustering", "marker")),
        title="Clustering controls",
        width=300,
    )
    return ui.nav_panel(
        "Clustering",
        ui.layout_sidebar(
            sidebar,
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
    )


def colocalization_ui():
    sidebar = ui.sidebar(
        ui.panel_conditional(
            "input.colocalization_mode === 'Observed'",
            ui.accordion(
                ui.accordion_panel("Cell population", selectize("coloc_celltype_filter", "Cell population", multiple=True)),
                ui.accordion_panel(
                    "Display",
                    ui.input_select("coloc_preset", "Heatmap preset", {"custom": "Custom", "report": "Report style"}),
                    ui.input_select("coloc_scope", "Heatmap scope", {"condition": "Analysis group summary", "sample_alias": "Sample summary", "celltype": "Cell type focus"}),
                    selectize("coloc_celltype_focus", "Cell type focus"),
                    ui.input_select("coloc_view", "Group display", {"focused": "Focused group", "compare": "Compare groups"}),
                    selectize("coloc_focus_group", "Displayed group"),
                    selectize("coloc_compare_groups", "Groups to compare", multiple=True),
                    ui.input_select("coloc_marker_mode", "Marker set", {"auto": "Variable detected markers", "manual": "Selected markers"}),
                    selectize("coloc_markers", "Heatmap markers", multiple=True),
                    ui.input_numeric("coloc_top_markers", "Top markers", 20, min=2, max=40),
                    ui.input_numeric("coloc_min_detected", "Minimum fraction detected", 0.25, min=0, max=1, step=0.05),
                    ui.input_numeric("coloc_min_range", "Minimum log2 range", 0.2, min=0, step=0.05),
                    ui.input_select("coloc_reference", "Reference group", []),
                    ui.input_select("coloc_ordering", "Marker ordering", {"ward": "Ward", "complete": "Complete", "average": "Average", "single": "Single"}),
                    ui.input_numeric("coloc_legend_min", "Legend minimum", -1, step=0.1),
                    ui.input_numeric("coloc_legend_max", "Legend maximum", 1, step=0.1),
                    selectize("coloc_detail_pair", "Pair detail"),
                    ui.input_action_button("apply_coloc", "Apply heatmap settings", class_="btn-primary w-100"),
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
                open=["Cell population", "Display"],
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
                    ui.input_action_button("coloc_run_differential", "Run differential analysis", class_="btn-primary w-100"),
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
            ui.navset_card_underline(
                ui.nav_panel(
                    "Observed",
                    ui.output_ui("coloc_notice"),
                    plot_pane("coloc_heatmap", height="900px"),
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


def data_popover():
    return ui.nav_control(
        ui.popover(
            ui.input_action_button("data_button", "Data", class_="btn btn-outline-light data-source-button"),
            ui.div(
                ui.input_text_area("pxl_spec", ".layout.pxl path(s), directory, or glob", default_pxl_spec(), rows=4),
                ui.input_text("annotation_path", "Optional cell annotation CSV", os.getenv("PROXIOME_ANNOTATIONS", "")),
                ui.input_text("patch_dir", "Optional patch-analysis directory", os.getenv("PROXIOME_PATCH_DIR", "")),
                ui.layout_columns(
                    ui.input_action_button("inspect_pxl", "Inspect", class_="btn-outline-secondary w-100"),
                    ui.input_task_button("load_pxl", "Load Data", class_="w-100"),
                    col_widths=(5, 7),
                ),
                ui.output_ui("load_status"),
                ui.hr(),
                ui.h6("Analysis grouping"),
                ui.input_select("grouping_column", "Sample-level metadata column", []),
                ui.input_text_area("custom_grouping", "Custom sample=group map", rows=5),
                ui.input_action_button("apply_grouping", "Apply grouping", class_="btn-outline-primary w-100"),
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
    ui.nav_panel("Spatial Metrics", ui.navset_tab(clustering_ui(), colocalization_ui(), id="spatial_metric_readout")),
    patch_ui(),
    ui.nav_spacer(),
    data_popover(),
    title="ProxiomeVis",
    id="readout_tab",
    fillable=["QC", "Abundance", "Spatial Metrics", "Patch Analysis"],
    header=ui.include_css(APP_DIR / "www" / "proixome.css"),
    navbar_options=ui.navbar_options(bg=TEAL, theme="dark", underline=True),
)


def empty_figure(message: str = "Load PXL data from the Data menu.") -> go.Figure:
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
    inspect_message = reactive.Value("No data loaded.")

    @ui.bind_task_button(button_id="load_pxl")
    @reactive.extended_task
    async def load_task(spec: str, annotations: str, patch_dir: str, force: bool) -> AppData:
        return await asyncio.to_thread(
            load_pxl_data,
            spec,
            annotation_path=annotations or None,
            patch_dir=patch_dir or None,
            force=force,
        )

    @reactive.effect
    @reactive.event(input.inspect_pxl)
    def _inspect_pxl():
        try:
            paths = resolve_pxl_paths(input.pxl_spec())
            size = sum(path.stat().st_size for path in paths)
            inspect_message.set(f"Ready: {len(paths)} PXL file(s), {size / 1024**3:.1f} GiB total.")
        except Exception as error:
            inspect_message.set(str(error))

    @reactive.effect
    @reactive.event(input.load_pxl)
    def _start_load():
        inspect_message.set("Loading PXL modalities and building the compact app cache…")
        load_task.invoke(input.pxl_spec(), input.annotation_path(), input.patch_dir(), False)

    @reactive.effect
    def _activate_loaded_data():
        loaded = load_task.result()
        data_state.set(loaded)
        inspect_message.set(
            f"Loaded {loaded.source['n_cells']:,} cells and {len(loaded.marker_options):,} selected markers."
        )
        samples = loaded.metadata[["sample_alias", "condition"]].drop_duplicates("sample_alias")
        mapping_text = "\n".join(f"{row.sample_alias}={row.condition}" for row in samples.itertuples())
        ui.update_text_area("custom_grouping", value=mapping_text)

    @output
    @render.ui
    def load_status():
        return ui.div(inspect_message(), class_="rds-load-status")

    @output
    @render.ui
    def source_summary():
        data = data_state.get()
        if data is None:
            return ui.p("Load one or more PXL files to begin.", class_="text-muted small")
        return ui.div(
            ui.hr(),
            ui.p(ui.strong(data.source.get("display_name", "PXL data"))),
            ui.p(f"{data.source['n_cells']:,} cells · {data.source['n_markers']:,} available markers"),
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
        pairs = sorted(data.colocalization["marker_pair"].dropna().astype(str).unique())
        numeric_qc = [
            column for column in ("n_umi", "n_edges", "reads_in_component", "isotype_fraction", "tau")
            if column in data.qc_origin and pd.api.types.is_numeric_dtype(data.qc_origin[column])
        ]

        ui.update_selectize("qc_sample_filter", choices=samples, selected=samples, server=True)
        ui.update_select("qc_metric", choices=numeric_qc, selected=numeric_qc[0] if numeric_qc else None)
        ui.update_select("abundance_embedding", choices=embeddings, selected=embeddings[0] if embeddings else None)

        single_marker_ids = [
            "abundance_marker", "abundance_distribution_marker", "abundance_diff_feature",
            "clustering_marker", "clustering_diff_feature", "coloc_diff_anchor",
        ]
        for input_id in single_marker_ids:
            ui.update_selectize(input_id, choices=markers, selected=markers[0] if markers else None, server=True)
        multiple_marker_ids = ["coloc_markers", "coloc_3d_markers"]
        for input_id in multiple_marker_ids:
            ui.update_selectize(input_id, choices=markers, selected=markers[: min(15, len(markers))], server=True)

        for input_id in (
            "abundance_condition_filter", "clustering_condition_filter", "clustering_heatmap_condition_filter",
            "coloc_condition_filter",
        ):
            ui.update_selectize(input_id, choices=conditions, selected=conditions, server=True)
        for input_id in (
            "abundance_celltype_filter", "abundance_diff_celltype_filter", "clustering_celltype_filter",
            "clustering_heatmap_celltype_filter", "clustering_diff_celltype_filter", "coloc_celltype_filter",
            "coloc_diff_celltype_filter", "coloc_3d_celltypes",
        ):
            ui.update_selectize(input_id, choices=celltypes, selected=celltypes, server=True)

        for prefix in ("abundance", "clustering"):
            ui.update_select(f"{prefix}_diff_group_a", choices=conditions, selected=conditions[0] if conditions else None)
            ui.update_select(f"{prefix}_diff_group_b", choices=conditions, selected=conditions[min(1, len(conditions) - 1)] if conditions else None)
        ui.update_select("coloc_diff_group_a", choices=conditions, selected=conditions[0] if conditions else None)
        ui.update_select("coloc_diff_group_b", choices=conditions, selected=conditions[min(1, len(conditions) - 1)] if conditions else None)
        ui.update_select("coloc_reference", choices=conditions, selected=conditions[0] if conditions else None)
        ui.update_selectize("coloc_celltype_focus", choices=celltypes, selected=celltypes[0] if celltypes else None, server=True)
        ui.update_selectize("coloc_focus_group", choices=conditions, selected=conditions[0] if conditions else None, server=True)
        ui.update_selectize("coloc_compare_groups", choices=conditions, selected=conditions[:6], server=True)
        ui.update_selectize("coloc_detail_pair", choices=pairs, selected=pairs[0] if pairs else None, server=True)
        ui.update_selectize("coloc_diff_pair", choices=pairs, selected=pairs[0] if pairs else None, server=True)
        ui.update_selectize("coloc_3d_sample", choices=samples, selected=samples[0] if samples else None, server=True)
        grouping_columns = sample_level_columns(metadata)
        ui.update_select("grouping_column", choices=grouping_columns, selected="condition" if "condition" in grouping_columns else (grouping_columns[0] if grouping_columns else None))

        labels = []
        marker_table = data.patch.get("marker_unmixing")
        if marker_table is not None and "label" in marker_table:
            labels = sorted(marker_table["label"].dropna().astype(str).unique())
        ui.update_selectize("patch_label_filter", choices=labels, selected=labels, server=True)

    @reactive.effect
    @reactive.event(input.apply_grouping)
    def _apply_grouping():
        data = data_state.get()
        if data is None:
            ui.notification_show("Load data before changing analysis grouping.", type="warning")
            return
        try:
            custom = input.custom_grouping().strip()
            if custom:
                mapping = parse_group_mapping(custom)
                label = "Custom sample groups"
            else:
                column = input.grouping_column()
                mapping = mapping_for_column(data.metadata, column)
                label = column
            data_state.set(apply_analysis_grouping(data, mapping, label))
            ui.notification_show("Analysis grouping applied.", type="message")
        except Exception as error:
            ui.notification_show(str(error), type="error", duration=None)

    @reactive.effect
    def _update_3d_components():
        data = data_state.get()
        if data is None:
            return
        sample = input.coloc_3d_sample()
        if not sample:
            return
        rows = data.metadata[data.metadata["sample_alias"].astype(str) == str(sample)]
        celltypes = selected(input.coloc_3d_celltypes(), rows["celltype_manual"])
        rows = rows[rows["celltype_manual"].astype(str).isin(celltypes)]
        components = rows["component"].astype(str).tolist()
        ui.update_selectize("coloc_3d_component", choices=components, selected=components[0] if components else None, server=True)

    def get_data() -> AppData | None:
        return data_state.get()

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
        frame = data.qc_origin if input.qc_metadata_source() == "origin" else data.qc_filtered
        samples = selected(input.qc_sample_filter(), frame["sample_alias"])
        return frame[frame["sample_alias"].astype(str).isin(samples)].copy()

    @output
    @render.ui
    def qc_metric_row():
        data = get_data()
        if data is None:
            return metric_boxes([("Loaded Cells", "—"), ("Final Cells", "—"), ("Retained", "—"), ("Samples", "—")])
        samples = selected(input.qc_sample_filter(), data.qc_origin["sample_alias"])
        origin = data.qc_origin[data.qc_origin["sample_alias"].astype(str).isin(samples)]
        filtered = data.qc_filtered[data.qc_filtered["sample_alias"].astype(str).isin(samples)]
        retained = len(filtered) / len(origin) if len(origin) else math.nan
        return metric_boxes([
            ("Loaded Cells", f"{len(origin):,}"), ("Final Cells", f"{len(filtered):,}"),
            ("Retained", f"{retained:.1%}" if np.isfinite(retained) else "—"), ("Samples", f"{origin['sample_alias'].nunique():,}"),
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
        figure = px.scatter(
            plot_data,
            x=x,
            y=y,
            color=color,
            facet_col=split,
            facet_col_wrap=max(1, int(input.abundance_split_columns() or 2)) if split else 0,
            hover_data=[column for column in ("component", "sample_alias", "condition", "celltype_manual") if column in plot_data],
            color_continuous_scale=["#edf7f4", "#78aeb2", "#f0b45b", CORAL] if color == "abundance" else None,
            render_mode="webgl",
        )
        figure.update_traces(marker={"size": float(input.abundance_point_size() or 3), "opacity": 0.82})
        figure.update_xaxes(title="Embedding 1", showgrid=False)
        figure.update_yaxes(title="Embedding 2", showgrid=False)
        return style_figure(figure, height=620 if split else 540)

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
        return fig_abundance_diff_volcano()

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

    def filtered_clustering(condition_input, celltype_input) -> pd.DataFrame:
        data = get_data()
        if data is None:
            return pd.DataFrame()
        conditions = selected(condition_input, data.clustering["condition"])
        celltypes = selected(celltype_input, data.clustering["celltype_manual"])
        return data.clustering[
            data.clustering["condition"].astype(str).isin(conditions)
            & data.clustering["celltype_manual"].astype(str).isin(celltypes)
        ].copy()

    def clustering_marker_rows() -> pd.DataFrame:
        rows = filtered_clustering(input.clustering_condition_filter(), input.clustering_celltype_filter())
        marker = input.clustering_marker()
        return rows[rows["marker"].astype(str) == str(marker)] if not rows.empty and marker else pd.DataFrame()

    def fig_clustering_observed():
        data = get_data()
        rows = clustering_marker_rows()
        if data is None or rows.empty:
            return empty_figure("No self-proximity values match the selected filters.")
        embedding = next(iter(embedding_columns(data.metadata)), None)
        if not embedding:
            return empty_figure("No two-dimensional embedding is available.")
        x, y = embedding_columns(data.metadata)[embedding]
        plot_data = rows.merge(data.metadata[["component", x, y]], on="component", how="inner")
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
            return empty_figure("No self-proximity values match the selected filters.")
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

    def clustering_heatmap_data() -> pd.DataFrame:
        rows = filtered_clustering(input.clustering_heatmap_condition_filter(), input.clustering_heatmap_celltype_filter())
        if rows.empty:
            return rows
        summary = rows.groupby(["marker", "condition", "celltype_manual"], observed=True)["log2_ratio"].mean().rename("mean_log2_ratio").reset_index()
        count = max(2, int(input.clustering_heatmap_marker_count() or 20))
        ranking = summary.groupby("marker", observed=True)["mean_log2_ratio"].std().fillna(0).sort_values(ascending=False)
        return summary[summary["marker"].isin(ranking.head(count).index)]

    def fig_clustering_heatmap():
        rows = clustering_heatmap_data()
        if rows.empty:
            return empty_figure("No clustering summary is available.")
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

    def clustering_diff_result():
        data = get_data()
        return differential_result("clustering", data.clustering if data is not None else pd.DataFrame(), ["marker"])

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
        return fig_clustering_diff_volcano()

    register_downloads("clustering_diff_volcano", fig_clustering_diff_volcano)

    def fig_clustering_diff_detail():
        data = get_data()
        if data is None:
            return empty_figure()
        marker = input.clustering_diff_feature()
        rows = data.clustering[
            (data.clustering["marker"].astype(str) == str(marker))
            & data.clustering["condition"].isin([input.clustering_diff_group_a(), input.clustering_diff_group_b()])
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
        ui.update_select("coloc_view", selected="compare")
        ui.update_select("coloc_marker_mode", selected="auto")
        ui.update_numeric("coloc_top_markers", value=15)
        ui.update_numeric("coloc_min_detected", value=0.25)
        ui.update_numeric("coloc_min_range", value=0.2)
        ui.update_select("coloc_mean_type", selected="population")
        ui.update_select("coloc_ordering", selected="ward")
        ui.update_numeric("coloc_legend_min", value=-0.75)
        ui.update_numeric("coloc_legend_max", value=0.75)

    def colocalization_scores() -> tuple[pd.DataFrame, pd.DataFrame]:
        data = get_data()
        if data is None:
            return pd.DataFrame(), pd.DataFrame()
        metadata = data.metadata.copy()
        conditions = selected(input.coloc_condition_filter(), metadata["condition"])
        celltypes = selected(input.coloc_celltype_filter(), metadata["celltype_manual"])
        metadata = metadata[
            metadata["condition"].astype(str).isin(conditions)
            & metadata["celltype_manual"].astype(str).isin(celltypes)
        ]
        if input.coloc_scope() == "celltype" and input.coloc_celltype_focus():
            metadata = metadata[metadata["celltype_manual"].astype(str) == str(input.coloc_celltype_focus())]
        scores = data.proximity[data.proximity["component"].isin(metadata["component"])].copy()
        if input.coloc_pixelator_filter():
            scores = filter_pixelator_proximity(
                scores,
                min_marker_fraction=float(input.coloc_min_fraction() or 0),
                min_marker_count=float(input.coloc_min_count() or 0),
                min_cells=int(input.coloc_min_cells() or 1),
            )
        return scores, metadata

    def auto_coloc_markers(scores: pd.DataFrame, metadata: pd.DataFrame, group_col: str) -> list[str]:
        data = get_data()
        if data is None or scores.empty:
            return []
        candidates = list(data.marker_options)
        pair_summary = scores.groupby([group_col, "marker_1", "marker_2"], observed=True).agg(
            mean=("log2_ratio", "mean"), detected=("component", "nunique")
        ).reset_index()
        totals = metadata.groupby(group_col, observed=True)["component"].nunique()
        pair_summary["pct"] = pair_summary.apply(lambda row: row.detected / max(1, totals.get(row[group_col], 1)), axis=1)
        ranges = pair_summary.groupby(["marker_1", "marker_2"], observed=True).agg(
            value_range=("mean", lambda values: values.max() - values.min()), max_pct=("pct", "max")
        ).reset_index()
        ranges = ranges[
            (ranges["max_pct"] >= float(input.coloc_min_detected() or 0))
            & (ranges["value_range"] >= float(input.coloc_min_range() or 0))
        ]
        marker_scores = pd.concat([
            ranges[["marker_1", "value_range"]].rename(columns={"marker_1": "marker"}),
            ranges[["marker_2", "value_range"]].rename(columns={"marker_2": "marker"}),
        ]).groupby("marker", observed=True)["value_range"].max().sort_values(ascending=False)
        selected_markers = [marker for marker in marker_scores.index if marker in candidates]
        count = max(2, int(input.coloc_top_markers() or 20))
        return (selected_markers + [marker for marker in candidates if marker not in selected_markers])[:count]

    def ordered_coloc_markers(summary: pd.DataFrame, markers: list[str], group_col: str) -> list[str]:
        if len(markers) < 3 or summary.empty:
            return markers
        reference = input.coloc_reference()
        group = reference if reference in set(summary[group_col].astype(str)) else str(summary[group_col].iloc[0])
        rows = summary[summary[group_col].astype(str) == group]
        matrix = rows.pivot(index="marker_1", columns="marker_2", values="mean_log2_ratio").reindex(index=markers, columns=markers).fillna(0)
        method = input.coloc_ordering() or "ward"
        try:
            return [markers[index] for index in leaves_list(linkage(matrix.to_numpy(), method=method))]
        except Exception:
            return markers

    def colocalization_heatmap_data() -> tuple[pd.DataFrame, str, list[str]]:
        scores, metadata = colocalization_scores()
        if scores.empty or metadata.empty:
            return pd.DataFrame(), "condition", []
        scope = input.coloc_scope() or "condition"
        group_col = "sample_alias" if scope == "sample_alias" else "condition"
        manual = selected(input.coloc_markers(), scores["marker_1"])
        markers = manual if input.coloc_marker_mode() == "manual" else auto_coloc_markers(scores, metadata, group_col)
        markers = [marker for marker in markers if marker in set(scores["marker_1"]).union(scores["marker_2"])]
        summary = summarize_spatial(
            scores,
            metadata,
            group_col=group_col,
            markers=markers,
            mean_type=input.coloc_mean_type() or "population",
        )
        if summary.empty:
            return summary, group_col, markers
        available_groups = list(summary[group_col].dropna().astype(str).unique())
        if input.coloc_view() == "compare":
            groups = selected(input.coloc_compare_groups(), available_groups)
            if len(markers) > 20:
                groups = [str(input.coloc_focus_group() or available_groups[0])]
        else:
            groups = [str(input.coloc_focus_group() or available_groups[0])]
        summary = summary[summary[group_col].astype(str).isin(groups)]
        markers = ordered_coloc_markers(summary, markers, group_col)
        summary["marker_1"] = pd.Categorical(summary["marker_1"], markers, ordered=True)
        summary["marker_2"] = pd.Categorical(summary["marker_2"], list(reversed(markers)), ordered=True)
        return summary, group_col, markers

    @output
    @render.ui
    def coloc_notice():
        summary, group_col, markers = colocalization_heatmap_data()
        if summary.empty:
            return ui.div("No spatial metric rows match the selected settings.", class_="alert alert-warning")
        notices = [f"Showing {len(markers)} markers across {summary[group_col].nunique()} {group_col.replace('_', ' ')} group(s)."]
        if input.coloc_view() == "compare" and len(markers) > 20:
            notices.append("Comparison view supports up to 20 markers; one focused group is shown.")
        elif input.coloc_view() == "compare" and len(markers) > 15:
            notices.append("For clearer comparison labels, use 15 or fewer markers.")
        return ui.div(" ".join(notices), class_="alert alert-info py-2")

    def fig_coloc_heatmap():
        summary, group_col, markers = colocalization_heatmap_data()
        if summary.empty:
            return empty_figure("No colocalization rows match the selected settings.")
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
            facet_col_wrap=2 if input.coloc_view() == "compare" else 1,
            color_continuous_scale="RdBu_r",
            range_color=(legend_min, legend_max),
            hover_data=["n_detected", "n_total", "pct_detected"],
            category_orders={"marker_1": markers, "marker_2": list(reversed(markers))},
            labels={"mean_log2_ratio": "Mean log2 ratio", "pct_detected": "Detected fraction"},
        )
        figure.update_traces(marker={"sizemin": 2, "line": {"width": 0.4, "color": "#556264"}})
        panels = max(1, summary[group_col].nunique())
        return style_figure(figure, height=max(620, 34 * len(markers) * math.ceil(panels / 2)))

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
        marker_1, marker_2 = pair.split(" / ", 1)
        return marker_1, marker_2

    def pair_detail_data() -> pd.DataFrame:
        scores, metadata = colocalization_scores()
        pair = pair_markers(input.coloc_detail_pair())
        if scores.empty or metadata.empty or pair is None:
            return pd.DataFrame()
        marker_1, marker_2 = pair
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
            input.coloc_mean_type() == "detected",
            detail["detected_mean"],
            detail["sum_log2_ratio"] / detail["n_total"],
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

    def coloc_sample_data() -> pd.DataFrame:
        data = get_data()
        if data is None:
            return pd.DataFrame()
        celltypes = selected(input.coloc_diff_celltype_filter(), data.metadata["celltype_manual"])
        summary = sample_colocalization(
            data.colocalization,
            data.metadata,
            celltypes=celltypes,
            mean_type=input.coloc_diff_mean() or "population",
        )
        summary = summary[
            summary["marker_1"].astype(str) < summary["marker_2"].astype(str)
        ].copy()
        summary["celltype_manual"] = "Pooled cell types"
        summary["sample_value"] = summary["mean_log2_ratio"]
        return summary

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
        if input.coloc_diff_pair_scope() == "anchor" and input.coloc_diff_anchor():
            anchor = str(input.coloc_diff_anchor())
            result = result[result["marker_pair"].str.split(" / ").apply(lambda pair: anchor in pair)]
        return result

    @output
    @render.ui
    def coloc_diff_method():
        rows = coloc_sample_data()
        if rows.empty:
            return ui.div("Load data to run sample-aware differential colocalization.", class_="alert alert-warning")
        counts = rows.groupby("condition", observed=True)["sample_alias"].nunique()
        a, b = input.coloc_diff_group_a(), input.coloc_diff_group_b()
        minimum = int(input.coloc_diff_min_samples() or 2)
        replicated = counts.get(a, 0) >= minimum and counts.get(b, 0) >= minimum
        message = (
            f"Sample-aware analysis. Effects compare median sample-level {input.coloc_diff_mean()} means. "
            f"{a}: {counts.get(a, 0)} sample(s); {b}: {counts.get(b, 0)} sample(s)."
        )
        if not replicated:
            message += f" At least {minimum} samples per group are required for p-values and FDR."
        return ui.div(ui.strong("Sample-aware analysis. "), message.removeprefix("Sample-aware analysis. "), class_=f"alert {'alert-info' if replicated else 'alert-warning'}")

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
        return fig_coloc_diff_volcano()

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
        sample = input.coloc_3d_sample()
        component = input.coloc_3d_component()
        if data is None or not sample or not component:
            return empty_figure("Choose a sample and component with a stored 3D layout.")
        try:
            nodes = read_component_layout(data, str(sample), str(component))
        except Exception as error:
            return empty_figure(str(error))
        if "marker" not in nodes:
            nodes["marker"] = "unlabeled"
        highlights = selected(input.coloc_3d_markers(), nodes["marker"])
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
            scale = centered.abs().max() or 1
            nodes[axis] = centered / scale
        figure = px.scatter_3d(
            nodes, x="x", y="y", z="z", color="marker_group", hover_data=[column for column in ("marker", "umi") if column in nodes],
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
        data = get_data()
        if data is None or not input.coloc_3d_component():
            return grid(pd.DataFrame())
        columns = [column for column in ("sample", "sample_alias", "condition", "celltype_manual", "component", "n_umi", "n_edges") if column in data.metadata]
        return grid(data.metadata.loc[data.metadata["component"].astype(str) == str(input.coloc_3d_component()), columns])

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
            return empty_figure("Patch analysis is available when PROXIOME_PATCH_DIR contains the patch tables.")
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
            return empty_figure("No Raji joint-proximity table is available.")
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
