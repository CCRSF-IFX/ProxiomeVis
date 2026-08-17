"""ProxiomeVis: the Python Shiny implementation."""

from __future__ import annotations

import asyncio
import math
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_plotly

from proxiome import (
    AppData,
    apply_analysis_grouping,
    calculate_differential,
    default_h5ad_path,
    load_h5ad_data,
    mapping_for_column,
    parse_group_mapping,
    resolve_h5ad_path,
    sample_level_columns,
    summarize_numeric,
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



def data_popover():
    return ui.nav_control(
        ui.popover(
            ui.input_action_button("data_button", "Data", class_="btn btn-outline-light data-source-button"),
            ui.div(
                ui.input_text_area("h5ad_path", "Processed .h5ad path", default_h5ad_path(), rows=4),
                ui.layout_columns(
                    ui.input_action_button("inspect_h5ad", "Inspect", class_="btn-outline-secondary w-100"),
                    ui.input_task_button("load_h5ad", "Load Data", class_="w-100"),
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
    ui.nav_spacer(),
    data_popover(),
    title="ProxiomeVis",
    id="readout_tab",
    fillable=["QC", "Abundance"],
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

    @ui.bind_task_button(button_id="load_h5ad")
    @reactive.extended_task
    async def load_task(path: str) -> AppData:
        return await asyncio.to_thread(load_h5ad_data, path)

    @reactive.effect
    @reactive.event(input.inspect_h5ad)
    def _inspect_h5ad():
        try:
            path = resolve_h5ad_path(input.h5ad_path())
            inspect_message.set(f"Ready: {path.name}, {path.stat().st_size / 1024**2:.1f} MiB.")
        except Exception as error:
            inspect_message.set(str(error))

    @reactive.effect
    @reactive.event(input.load_h5ad)
    def _start_load():
        inspect_message.set("Loading processed AnnData…")
        load_task.invoke(input.h5ad_path())

    @reactive.effect
    def _activate_loaded_data():
        loaded = load_task.result()
        data_state.set(loaded)
        inspect_message.set(
            f"Loaded {loaded.source['n_cells']:,} cells and {len(loaded.marker_options):,} markers."
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
            return ui.p("Load a processed H5AD file to begin.", class_="text-muted small")
        return ui.div(
            ui.hr(),
            ui.p(ui.strong(data.source.get("display_name", "AnnData"))),
            ui.p(f"{data.source['n_cells']:,} cells · {data.source['n_markers']:,} markers"),
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
        for input_id in ("abundance_marker", "abundance_distribution_marker", "abundance_diff_feature"):
            ui.update_selectize(input_id, choices=markers, selected=markers[0] if markers else None, server=True)
        ui.update_selectize("abundance_condition_filter", choices=conditions, selected=conditions, server=True)
        for input_id in ("abundance_celltype_filter", "abundance_diff_celltype_filter"):
            ui.update_selectize(input_id, choices=celltypes, selected=celltypes, server=True)
        ui.update_select("abundance_diff_group_a", choices=conditions, selected=conditions[0] if conditions else None)
        ui.update_select(
            "abundance_diff_group_b",
            choices=conditions,
            selected=conditions[min(1, len(conditions) - 1)] if conditions else None,
        )
        grouping_columns = sample_level_columns(metadata)
        selected_grouping = "condition" if "condition" in grouping_columns else (grouping_columns[0] if grouping_columns else None)
        ui.update_select("grouping_column", choices=grouping_columns, selected=selected_grouping)

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


app = App(app_ui, server, static_assets=APP_DIR / "www")
