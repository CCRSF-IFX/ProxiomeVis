from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from proxiome import (
    AppData,
    adjust_bh,
    apply_analysis_grouping,
    calculate_differential,
    load_h5ad_data,
    parse_group_mapping,
    resolve_h5ad_path,
    sample_level_columns,
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
        qc_filter_counts=qc_counts,
    )


def test_h5ad_is_the_only_supported_input(tmp_path: Path):
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
    adata.write_h5ad(path)

    data = load_h5ad_data(path)
    assert resolve_h5ad_path(path) == path.resolve()
    assert data.source["source_type"] == "h5ad"
    assert data.marker_options == ("CD3", "CD19")
    assert set(data.metadata) >= {"component", "umap_1", "umap_2"}
    assert data.qc_filter_counts["n_cells"].sum() == 6
    with pytest.raises(ValueError, match="Expected an .h5ad"):
        resolve_h5ad_path(tmp_path / "input.pxl")


def test_analysis_grouping_updates_every_sample_backed_table():
    data = make_data()
    assert set(sample_level_columns(data.metadata)) >= {"condition", "batch"}
    grouped = apply_analysis_grouping(data, {"s1": "baseline", "s2": "treated"}, "Custom")
    assert grouped.metadata["condition"].tolist() == ["baseline", "baseline", "treated", "treated"]
    assert grouped.qc_filter_counts["condition"].tolist() == ["baseline", "treated"]
    assert parse_group_mapping("s1=baseline\ns2=treated") == {"s1": "baseline", "s2": "treated"}


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


def test_python_app_exposes_only_h5ad_qc_and_abundance():
    import app

    html = str(app.app_ui)
    for label in (
        "QC", "Filtering", "Cell Calling", "Distributions", "Metadata",
        "Abundance", "Observed", "Marker Distributions", "Cell Annotation", "Differential",
        "Processed .h5ad path", "Data source", "Analysis grouping",
    ):
        assert label in html
    for label in (".layout.pxl", "Spatial Metrics", "Patch Analysis", "3D Layout"):
        assert label not in html

    embeddings = app.embedding_columns(
        pd.DataFrame(columns=["umap_1", "umap_2", "k_core_1", "k_core_2"])
    )
    assert embeddings == {"umap": ("umap_1", "umap_2")}
