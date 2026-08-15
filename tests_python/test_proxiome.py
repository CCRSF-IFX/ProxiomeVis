from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from proxiome import (
    AppData,
    adjust_bh,
    apply_analysis_grouping,
    calculate_differential,
    filter_pixelator_proximity,
    parse_group_mapping,
    resolve_pxl_paths,
    sample_level_columns,
    select_markers,
    summarize_spatial,
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
    proximity = pd.DataFrame(
        {
            "component": ["c1", "c3"],
            "marker_1": ["A", "A"],
            "marker_2": ["B", "B"],
            "log2_ratio": [2.0, 4.0],
            "marker_1_freq": [0.2, 0.2],
            "marker_2_freq": [0.2, 0.01],
            "min_count": [20, 3],
            "sample_alias": ["s1", "s2"],
            "condition": ["A", "B"],
            "celltype_manual": ["T", "T"],
        }
    )
    qc_counts = pd.DataFrame(
        {
            "sample": ["s1", "s2"],
            "condition": ["A", "B"],
            "step": ["loaded", "loaded"],
            "n_cells": [2, 2],
        }
    )
    return AppData(
        source={},
        marker_options=("A", "B"),
        metadata=metadata,
        abundance=pd.DataFrame(),
        proximity=proximity,
        clustering=pd.DataFrame(),
        colocalization=proximity,
        qc_origin=metadata,
        qc_filtered=metadata,
        qc_filter_counts=qc_counts,
        patch={},
        pxl_files=(),
    )


def test_resolve_pxl_paths_and_marker_selection(tmp_path: Path):
    first = tmp_path / "a.layout.pxl"
    second = tmp_path / "b.layout.pxl"
    first.touch()
    second.touch()
    assert resolve_pxl_paths(tmp_path) == [first.resolve(), second.resolve()]
    assert resolve_pxl_paths(f"{first},{second}") == [first.resolve(), second.resolve()]
    assert select_markers(["CD8", "CD3e", "X"]) == ["CD3e", "CD8"]
    with pytest.raises(ValueError, match="No .pxl"):
        resolve_pxl_paths(tmp_path / "missing*.pxl")


def test_analysis_grouping_updates_every_sample_backed_table():
    data = make_data()
    assert set(sample_level_columns(data.metadata)) >= {"condition", "batch"}
    grouped = apply_analysis_grouping(data, {"s1": "baseline", "s2": "treated"}, "Custom")
    assert grouped.metadata["condition"].tolist() == ["baseline", "baseline", "treated", "treated"]
    assert grouped.proximity["condition"].tolist() == ["baseline", "treated"]
    assert grouped.qc_filter_counts["condition"].tolist() == ["baseline", "treated"]
    assert parse_group_mapping("s1=baseline\ns2=treated") == {"s1": "baseline", "s2": "treated"}


def test_population_spatial_mean_includes_undetected_cells_and_symmetrizes():
    data = make_data()
    summary = summarize_spatial(
        data.proximity,
        data.metadata,
        group_col="condition",
        markers=["A", "B"],
        mean_type="population",
    )
    a_to_b = summary[(summary.condition == "A") & (summary.marker_1 == "A") & (summary.marker_2 == "B")].iloc[0]
    b_to_a = summary[(summary.condition == "A") & (summary.marker_1 == "B") & (summary.marker_2 == "A")].iloc[0]
    assert a_to_b.mean_log2_ratio == 1.0
    assert a_to_b.pct_detected == 0.5
    assert b_to_a.mean_log2_ratio == a_to_b.mean_log2_ratio


def test_pixelator_filters_use_native_fraction_count_and_cell_thresholds():
    rows = make_data().proximity
    filtered = filter_pixelator_proximity(rows, min_marker_fraction=0.05, min_marker_count=10)
    assert filtered.component.tolist() == ["c1"]
    assert filter_pixelator_proximity(rows, min_cells=3).empty


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


def test_python_app_contains_all_r_app_sections():
    import app

    html = str(app.app_ui)
    for label in (
        "QC", "Filtering", "Cell Calling", "Distributions", "Metadata",
        "Abundance", "Observed", "Marker Distributions", "Cell Annotation", "Differential",
        "Spatial Metrics", "Clustering", "Summary Heatmap", "Colocalization", "Pair detail", "3D Layout",
        "Patch Analysis", "Markers", "Raji Signal", "Patch Burden", "Data source", "Analysis grouping",
    ):
        assert label in html

    embeddings = app.embedding_columns(
        pd.DataFrame(columns=["umap_1", "umap_2", "k_core_1", "k_core_2"])
    )
    assert embeddings == {"umap": ("umap_1", "umap_2")}
