from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from proxiome import (
    AppData,
    adjust_bh,
    analysis_grouping_summary,
    apply_analysis_grouping,
    build_spatial_retrieval,
    calculate_differential,
    clear_pxl_cache,
    load_h5ad_data,
    load_pxl_proximity,
    new_analysis_grouping_config,
    read_component_layout,
    resolve_pxl_paths,
    resolve_h5ad_path,
    sample_pxl_colocalization,
    sample_level_columns,
    select_colocalization_heatmap_markers,
    summarize_spatial,
    update_analysis_grouping_config,
)


def make_data() -> AppData:
    metadata = pd.DataFrame(
        {
            "component": ["c1", "c2", "c3", "c4"],
            "sample": ["s1", "s1", "s2", "s2"],
            "sample_alias": ["s1", "s1", "s2", "s2"],
            "condition": ["A", "A", "B", "B"],
            "celltype_manual": ["T", "T", "T", "T"],
            "batch": ["x", "x", "y", "y"],
            "umap_1": [0, 1, 2, 3],
            "umap_2": [0, 1, 0, 1],
        }
    )
    qc_counts = pd.DataFrame(
        {
            "sample": ["s1", "s2", "TOTAL"],
            "condition": ["A", "B", "TOTAL"],
            "step": ["loaded", "loaded", "loaded"],
            "n_cells": [2, 2, 4],
        }
    )
    return AppData(
        source={},
        marker_options=("A", "B"),
        metadata=metadata,
        abundance=pd.DataFrame(),
        qc_filter_counts=qc_counts,
    )


def test_h5ad_analysis_payload_and_optional_pxl_layout(tmp_path: Path, monkeypatch):
    from anndata import AnnData

    path = tmp_path / "processed.h5ad"
    adata = AnnData(
        X=np.array([[1, 2], [3, 4]], dtype=np.uint32),
        obs=pd.DataFrame(
            {
                "sample": ["s1", "s2"],
                "sample_alias": ["s1", "s2"],
                "condition": ["A", "B"],
                "celltype_manual": ["T", "B"],
                "n_umi": [10, 20],
            },
            index=["c1", "c2"],
        ),
        var=pd.DataFrame(index=["CD3", "CD19"]),
    )
    adata.layers["clr"] = np.log1p(adata.X).astype(np.float32)
    adata.obsm["X_umap"] = np.array([[0, 0], [1, 1]], dtype=float)
    adata.uns["qc_cell_counts_by_step"] = {
        "sample": ["s1", "s2", "s1", "s2"],
        "step": ["00_loaded", "00_loaded", "03_after_isotype_filter", "03_after_isotype_filter"],
        "n_cells": [2, 2, 1, 1],
        "fraction_loaded": [1.0, 1.0, 0.5, 0.5],
    }
    adata.uns["proxiome"] = {
        "proximity": pd.DataFrame(
            {
                "component": ["c1", "c1", "c2"],
                "marker_1": ["CD3", "CD3", "CD19"],
                "marker_2": ["CD3", "CD19", "CD19"],
                "log2_ratio": [0.4, -0.2, 0.6],
            },
            index=["p1", "p2", "p3"],
        ),
        "patch": {
            "run_plan": pd.DataFrame(
                {"run_patch_detection": [True], "n_patch_markers": [1]}, index=["plan"]
            ),
            "marker_unmixing": pd.DataFrame(
                {"marker": ["CD19"], "label": ["patch"], "receiver_freq": [0.1], "target_freq": [0.8]},
                index=["CD19"],
            ),
        },
        "component_layouts": {
            "c1": pd.DataFrame(
                {"x": [0.0], "y": [0.0], "z": [0.0], "marker": ["CD3"]}, index=["node1"]
            )
        },
    }
    adata.write_h5ad(path)
    pxl_path = tmp_path / "s1.layout.pxl"
    pxl_path.write_bytes(b"")

    data = load_h5ad_data(path, pxl_spec=pxl_path)
    assert resolve_h5ad_path(path) == path.resolve()
    assert data.source["source_type"] == "h5ad"
    assert data.marker_options == ("CD3", "CD19")
    assert set(data.metadata) >= {"component", "umap_1", "umap_2"}
    assert data.qc_filter_counts["n_cells"].sum() == 6
    assert data.source["has_spatial_metrics"] is True
    assert not hasattr(data, "proximity")
    assert data.patch["marker_unmixing"].iloc[0]["marker"] == "CD19"
    assert data.component_layouts["c1"].iloc[0]["marker"] == "CD3"
    assert data.pxl_files == (str(pxl_path.resolve()),)
    assert resolve_pxl_paths("") == []
    layout = pd.DataFrame({"x": [1.0], "y": [2.0], "z": [3.0], "marker": ["CD3"]})

    class FakeDataset:
        def components(self):
            return {"c1"}

        def filter(self, *, components):
            assert components == ["c1"]
            return self

        def precomputed_layouts(self, *, add_marker_counts):
            assert add_marker_counts is False
            return SimpleNamespace(to_df=lambda: layout)

    clear_pxl_cache()
    monkeypatch.setitem(__import__("sys").modules, "pixelator", SimpleNamespace(read_pna=lambda _: FakeDataset()))
    pxl_data = replace(data, component_layouts={})
    assert read_component_layout(pxl_data, "s1", "c1").equals(layout)
    with pytest.raises(ValueError, match="Expected an .h5ad"):
        resolve_h5ad_path(tmp_path / "input.pxl")


def test_proximity_is_queried_from_pxl_and_joined_to_h5ad_metadata(tmp_path: Path, monkeypatch):
    pxl_path = tmp_path / "sample.layout.pxl"
    pxl_path.write_bytes(b"")
    data = replace(make_data(), pxl_files=(str(pxl_path),))
    proximity = pd.DataFrame(
        {
            "component": ["c1", "c1"],
            "marker_1": ["A", "A"],
            "marker_2": ["A", "B"],
            "log2_ratio": [0.4, -0.2],
            "marker_1_freq": [0.2, 0.2],
            "marker_2_freq": [0.2, 0.3],
            "min_count": [4, 4],
        }
    )

    class FakeDataset:
        def components(self):
            return {"c1", "c2"}

        def filter(self, *, components, markers):
            assert components == ["c1", "c2"]
            assert markers == ["A", "B"]
            return self

        def proximity(self, *, add_marker_counts, add_logratio, calculate_from_edgelist):
            assert add_marker_counts and add_logratio and not calculate_from_edgelist
            return SimpleNamespace(to_df=lambda: proximity)

    clear_pxl_cache()
    monkeypatch.setitem(__import__("sys").modules, "pixelator", SimpleNamespace(read_pna=lambda _: FakeDataset()))
    rows = load_pxl_proximity(
        data,
        data.metadata.iloc[:2],
        markers=["A", "B"],
        pair_type="nonself",
        add_marker_counts=True,
    )
    assert rows[["marker_1", "marker_2"]].values.tolist() == [["A", "B"]]
    assert rows.iloc[0]["condition"] == "A"
    assert rows.iloc[0]["sample_alias"] == "s1"
    clear_pxl_cache()


def test_analysis_grouping_edits_and_resets_conditions_from_preserved_source():
    data = make_data()
    assert set(sample_level_columns(data.metadata)) >= {"sample_alias", "sample", "condition", "batch"}
    config = new_analysis_grouping_config(data.metadata)
    assert config["column"] == "condition"
    by_sample = update_analysis_grouping_config(config, mode="column", column="sample_alias")
    assert by_sample["mapping"] == {"s1": "s1", "s2": "s2"}
    config = update_analysis_grouping_config(
        config,
        mode="custom",
        column="condition",
        custom_groups={"s1": "baseline", "s2": "stimulated"},
    )
    grouped = apply_analysis_grouping(data, config["mapping"], config["label"])
    assert grouped.metadata["condition"].tolist() == ["baseline", "baseline", "stimulated", "stimulated"]
    assert grouped.qc_filter_counts["condition"].tolist() == ["baseline", "stimulated", "TOTAL"]
    assert config["source"]["condition"].tolist() == ["A", "B"]
    assert analysis_grouping_summary(config) == "Analysis grouping: Custom sample groups · 2 groups"
    with pytest.raises(ValueError, match="cannot be blank"):
        update_analysis_grouping_config(
            config,
            mode="custom",
            column="condition",
            custom_groups={"s1": "", "s2": "stimulated"},
        )
    reset = update_analysis_grouping_config(config, mode="column", column="condition")
    restored = apply_analysis_grouping(grouped, reset["mapping"], reset["label"])
    assert restored.metadata["condition"].tolist() == ["A", "A", "B", "B"]


def test_differential_uses_median_effect_and_bh_adjustment():
    rows = pd.DataFrame(
        {
            "marker": ["M"] * 8,
            "condition": ["A"] * 4 + ["B"] * 4,
            "celltype_manual": ["T"] * 8,
            "abundance": [5, 6, 7, 8, 1, 2, 3, 4],
        }
    )
    result = calculate_differential(
        rows,
        feature_cols=["marker"],
        value_col="abundance",
        group_a="A",
        group_b="B",
        min_observations=3,
    )
    assert result.iloc[0].effect_size == 4
    assert result.iloc[0].n_a == result.iloc[0].n_b == 4
    assert np.isfinite(result.iloc[0].p_adj)
    assert np.allclose(adjust_bh([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])


def test_heatmap_markers_use_equal_sample_weighting_and_manual_override():
    metadata = pd.DataFrame(
        {
            "component": ["s1a", "s1b", "s2a", "s1c", "s2b"],
            "sample_alias": ["s1", "s1", "s2", "s1", "s2"],
            "celltype_manual": ["T", "T", "T", "B", "B"],
        }
    )
    values = {
        "s1a": {"A": 10, "B": 0, "PD-L1": 1},
        "s1b": {"A": 10, "B": 0, "PD-L1": 1},
        "s2a": {"A": 0, "B": 12, "PD-L1": 1},
        "s1c": {"A": 20, "B": 0, "PD-L1": 1},
        "s2b": {"A": 20, "B": 0, "PD-L1": 1},
    }
    abundance = pd.DataFrame(
        [
            {"component": component, "marker": marker, "abundance": value}
            for component, markers in values.items()
            for marker, value in markers.items()
        ]
    )
    available = ["A", "B", "PD-L1"]

    t_cells = metadata[metadata["celltype_manual"] == "T"]
    assert select_colocalization_heatmap_markers(
        abundance, t_cells, available, n_markers=1
    ) == ["B"]
    b_cells = metadata[metadata["celltype_manual"] == "B"]
    assert select_colocalization_heatmap_markers(
        abundance, b_cells, available, n_markers=1
    ) == ["A"]
    assert select_colocalization_heatmap_markers(
        abundance,
        metadata,
        available,
        plot_markers=["PD-L1", "B"],
    ) == ["PD-L1", "B"]


def test_spatial_summary_requires_both_markers_in_selected_set():
    metadata = pd.DataFrame(
        {
            "component": ["c1", "c2"],
            "condition": ["control", "control"],
            "celltype_manual": ["T", "T"],
        }
    )
    scores = pd.DataFrame(
        {
            "component": ["c1", "c1", "c2"],
            "marker_1": ["A", "B", "A"],
            "marker_2": ["B", "C", "A"],
            "log2_ratio": [0.5, 0.7, 0.2],
            "condition": ["control", "control", "control"],
        }
    )
    summary = summarize_spatial(
        scores,
        metadata,
        group_col="condition",
        markers=["A", "B"],
    )
    assert not ((summary["marker_1"] == "C") | (summary["marker_2"] == "C")).any()


def test_spatial_retrieval_freezes_population_and_defaults_to_all_markers(monkeypatch):
    data = replace(make_data(), pxl_files=("fake.pxl",))
    monkeypatch.setattr(
        "proxiome.filter_pxl_metadata",
        lambda _data, metadata=None: _data.metadata if metadata is None else metadata,
    )
    request = (("A",), ("s1",), ("T",), "all", None, ())
    retrieval = build_spatial_retrieval(
        data,
        conditions=["A"],
        samples=["s1"],
        celltypes=["T"],
        retrieved_at="now",
        request=request,
    )

    assert retrieval.request == request
    assert retrieval.metadata["component"].tolist() == ["c1", "c2"]
    assert retrieval.markers == ("A", "B")
    data.metadata.loc[data.metadata["component"] == "c1", "condition"] = "changed"
    assert retrieval.metadata["condition"].tolist() == ["A", "A"]


def test_sample_colocalization_query_is_restricted_to_retrieved_markers(monkeypatch):
    data = replace(make_data(), pxl_files=("fake.pxl",))
    captured = {}

    class FakeConnection:
        def execute(self, query, parameters):
            captured["query"] = query
            captured["parameters"] = parameters
            return self

        def fetchdf(self):
            return pd.DataFrame(
                columns=[
                    "sample",
                    "marker_1",
                    "marker_2",
                    "sum_log2_ratio",
                    "detected_mean",
                    "n_detected",
                ]
            )

    class FakeSession:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def get_connection(self):
            return FakeConnection()

    fake_dataset = SimpleNamespace(view=SimpleNamespace(open=FakeSession))
    monkeypatch.setattr(
        "proxiome.filter_pxl_metadata",
        lambda _data, metadata=None: _data.metadata if metadata is None else metadata,
    )
    monkeypatch.setattr("proxiome._read_pxl_dataset", lambda _paths: fake_dataset)

    rows = sample_pxl_colocalization(data, data.metadata, markers=["A", "B"])

    assert "marker_1 IN $markers AND marker_2 IN $markers" in captured["query"]
    assert captured["parameters"]["markers"] == ["A", "B"]
    assert set(rows["marker_pair"]) == {"A / B"}


def test_python_app_exposes_h5ad_and_pxl_spatial_modules():
    import app
    import plotly.express as px

    html = str(app.app_ui)
    source = Path("app.py").read_text()
    for label in (
        "QC", "Filtering", "Cell Calling", "Distributions", "Metadata",
            "Abundance", "Observed", "Marker Distributions", "Cell Annotation", "Differential",
            "Spatial Metrics", "Retrieve Data", "Retrieve Spatial Data", "Number of markers",
            "All markers", "Clustering", "Colocalization", "3D Layout", "Patch Analysis",
            "Activity Log", "Clear log",
            "Processed .h5ad path", "PXL path(s) for proximity and cellgraph data",
            "Data source", "Analysis grouping",
        ):
            assert label in html
    for label in ("Optional patch-analysis directory", "Analyze PXL"):
        assert label not in html
    assert 'id="configure_analysis_grouping"' in html
    assert 'id="analysis_grouping_summary"' in html
    assert 'id="custom_grouping"' not in html
    assert 'id="apply_coloc"' not in html
    assert 'id="coloc_run_differential"' not in html
    assert 'id="clustering_run_differential"' not in html
    assert 'class="data-source-actions"' in html
    assert "grid-template-columns: repeat(2, minmax(0, 1fr))" in Path("www/proixome.css").read_text()
    for token in ("Use metadata column", "Edit sample groups", "Reset to condition", "analysis_group_editor"):
        assert token in source
    assert "input.abundance_color_by === 'abundance'" in source
    assert (
        'id="abundance_umap" class="shiny-ipywidget-output '
        'shiny-report-size shiny-report-theme" style="height:auto;"'
    ) in html

    facets = pd.DataFrame(
        {"x": list(range(6)), "y": list(range(6)), "sample": list("ABCDEF")}
    )
    figure = px.scatter(facets, x="x", y="y", facet_col="sample", facet_col_wrap=3)
    app.label_embedding_axes(figure)
    assert [axis.title.text for axis in figure.select_xaxes()] == ["Embedding 1"] * 3 + [None] * 3
    assert app.split_umap_height(6, 3) == 760
    assert app.split_umap_height(6, 1) == 2040

    embeddings = app.embedding_columns(
        pd.DataFrame(columns=["umap_1", "umap_2", "k_core_1", "k_core_2"])
    )
    assert embeddings == {"umap": ("umap_1", "umap_2")}
