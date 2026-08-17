"""Data and analysis helpers for the Python ProxiomeVis app."""

from __future__ import annotations

import os
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


REFERENCE_H5AD = Path(
    "/Volumes/ccrsf-static/illumina/CCRSFIFX-23_MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/"
    "python_results/pg_data_combined_filtered_annotated.h5ad"
)



@dataclass(frozen=True)
class AppData:
    source: dict
    marker_options: tuple[str, ...]
    metadata: pd.DataFrame
    abundance: pd.DataFrame
    qc_filter_counts: pd.DataFrame


def resolve_h5ad_path(spec: str | Path) -> Path:
    """Validate one processed AnnData file supplied by server-visible path."""
    path = Path(os.path.expandvars(os.path.expanduser(str(spec).strip()))).resolve()
    if path.suffix.lower() != ".h5ad":
        raise ValueError(f"Expected an .h5ad file: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"H5AD file not found: {path}")
    return path


def default_h5ad_path() -> str:
    return os.getenv("PROXIOME_H5AD", str(REFERENCE_H5AD)).strip()


def load_h5ad_data(spec: str | Path) -> AppData:
    """Load a processed H5AD without rerunning Pixelator analysis."""
    from anndata import read_h5ad

    path = resolve_h5ad_path(spec)
    adata = read_h5ad(path)
    if not adata.n_obs or not adata.n_vars:
        raise ValueError("H5AD must contain at least one observation and one marker.")

    metadata = build_metadata(adata)
    markers = tuple(map(str, adata.var_names))
    return AppData(
        source={
            "display_name": path.name,
            "source_type": "h5ad",
            "n_cells": int(adata.n_obs),
            "n_markers": int(adata.n_vars),
            "h5ad_path": str(path),
            "analysis_group_label": "condition",
        },
        marker_options=markers,
        metadata=metadata,
        abundance=build_abundance(adata),
        qc_filter_counts=build_h5ad_qc_filter_counts(adata, metadata),
    )




def build_metadata(adata) -> pd.DataFrame:
    metadata = adata.obs.copy()
    metadata.insert(0, "component", metadata.index.astype(str))
    metadata = metadata.reset_index(drop=True)
    if "sample" not in metadata:
        metadata["sample"] = "sample"
    metadata["sample"] = metadata["sample"].astype("string").fillna("sample").astype(str)
    if "sample_alias" not in metadata:
        metadata["sample_alias"] = metadata["sample"]
    if "condition" not in metadata:
        metadata["condition"] = metadata["sample"]
    if "celltype_manual" not in metadata:
        metadata["celltype_manual"] = "unannotated"

    metadata["condition"] = metadata["condition"].astype("string").fillna(metadata["sample"]).astype(str)
    metadata["celltype_manual"] = metadata["celltype_manual"].astype("string").fillna("unannotated").astype(str)
    metadata["sample_alias"] = metadata["sample_alias"].astype("string").fillna(metadata["sample"]).astype(str)

    embedding_names = ("umap", "pca", "harmony", "tsne")
    embeddings = {
        key.removeprefix("X_"): np.asarray(value)
        for key, value in adata.obsm.items()
        if np.asarray(value).ndim == 2
        and np.asarray(value).shape[1] >= 2
        and any(name in key.lower() for name in embedding_names)
    }
    for name, values in embeddings.items():
        for index in range(min(values.shape[1], 30)):
            metadata[f"{name}_{index + 1}"] = values[:, index]
    return metadata

def build_abundance(adata) -> pd.DataFrame:
    counts = adata.to_df()
    abundance = adata.to_df("clr") if "clr" in adata.layers else np.log1p(counts)
    counts.index = counts.index.astype(str)
    abundance.index = abundance.index.astype(str)
    counts_long = counts.rename_axis("component").stack().rename("count")
    abundance_long = abundance.rename_axis("component").stack().rename("abundance")
    result = pd.concat([abundance_long, counts_long], axis=1).reset_index()
    result = result.rename(columns={result.columns[1]: "marker"})
    result["count"] = pd.to_numeric(result["count"], errors="coerce").fillna(0)
    result["abundance"] = pd.to_numeric(result["abundance"], errors="coerce")
    return result



def build_qc_filter_counts(metadata: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample, values in metadata.groupby("sample_alias", observed=True):
        condition = values["condition"].iloc[0]
        count = len(values)
        rows.extend(
            (
                {"sample": sample, "condition": condition, "step": "00_loaded", "step_label": "Loaded", "n_cells": count, "fraction_loaded": 1.0},
                {"sample": sample, "condition": condition, "step": "99_filtered", "step_label": "Filtered", "n_cells": count, "fraction_loaded": 1.0},
            )
        )
    return pd.DataFrame(rows)


def build_h5ad_qc_filter_counts(adata, metadata: pd.DataFrame) -> pd.DataFrame:
    """Use QC history stored by the reference notebook, with a safe fallback."""
    history = adata.uns.get("qc_cell_counts_by_step")
    if not isinstance(history, Mapping):
        return build_qc_filter_counts(metadata)

    rows = pd.DataFrame(history)
    required = {"sample", "step", "n_cells"}
    if rows.empty or not required.issubset(rows):
        return build_qc_filter_counts(metadata)

    sample_metadata = metadata[["sample", "sample_alias", "condition"]].drop_duplicates("sample")
    sample_alias = sample_metadata.set_index("sample")["sample_alias"].astype(str)
    conditions = sample_metadata.set_index("sample")["condition"].astype(str)
    rows["sample"] = rows["sample"].astype(str)
    rows["condition"] = rows["sample"].map(conditions).fillna(rows["sample"])
    rows["sample"] = rows["sample"].map(sample_alias).fillna(rows["sample"])
    rows["n_cells"] = pd.to_numeric(rows["n_cells"], errors="coerce").fillna(0).astype(int)
    labels = {
        "00_loaded": "Loaded",
        "01_after_n_umi_filter": "After UMI filter",
        "02_after_tau_filter": "After tau filter",
        "03_after_isotype_filter": "After isotype filter",
    }
    rows["step_label"] = rows["step"].astype(str).map(labels).fillna(rows["step"].astype(str))
    if "fraction_loaded" not in rows:
        loaded = rows.loc[rows["step"].astype(str) == "00_loaded"].set_index("sample")["n_cells"]
        rows["fraction_loaded"] = rows.apply(
            lambda row: row["n_cells"] / loaded.get(row["sample"], row["n_cells"] or 1), axis=1
        )
    else:
        rows["fraction_loaded"] = pd.to_numeric(rows["fraction_loaded"], errors="coerce")
    return rows[["sample", "condition", "step", "step_label", "n_cells", "fraction_loaded"]]


def sample_level_columns(metadata: pd.DataFrame) -> list[str]:
    ignored = {"component", "sample", "sample_alias", "celltype_manual"}
    result = []
    for column in metadata.columns:
        if column in ignored or column.rsplit("_", 1)[-1].isdigit():
            continue
        if metadata.groupby("sample_alias", observed=True)[column].nunique(dropna=False).max() <= 1:
            result.append(column)
    return result


def apply_analysis_grouping(data: AppData, mapping: Mapping[str, str], label: str) -> AppData:
    clean = {str(sample): str(group).strip() for sample, group in mapping.items() if str(group).strip()}
    samples = set(data.metadata["sample_alias"].astype(str))
    if samples.difference(clean):
        raise ValueError("Every sample must have a non-empty analysis group.")

    metadata = data.metadata.copy()
    metadata["condition"] = metadata["sample_alias"].astype(str).map(clean)

    def regroup(frame: pd.DataFrame) -> pd.DataFrame:
        frame = frame.copy()
        if "sample_alias" in frame:
            frame["condition"] = frame["sample_alias"].astype(str).map(clean)
        elif "sample" in frame and set(frame["sample"].astype(str)).issubset(clean):
            frame["condition"] = frame["sample"].astype(str).map(clean)
        elif "component" in frame:
            conditions = metadata.set_index("component")["condition"]
            frame["condition"] = frame["component"].astype(str).map(conditions)
        return frame

    source = dict(data.source, analysis_group_label=label, analysis_group_count=len(set(clean.values())))
    return replace(
        data,
        source=source,
        metadata=metadata,
        qc_filter_counts=regroup(data.qc_filter_counts),
    )


def mapping_for_column(metadata: pd.DataFrame, column: str) -> dict[str, str]:
    if column not in sample_level_columns(metadata):
        raise ValueError(f"{column!r} is not constant within samples.")
    rows = metadata[["sample_alias", column]].drop_duplicates("sample_alias")
    return dict(zip(rows["sample_alias"].astype(str), rows[column].fillna("missing").astype(str)))


def parse_group_mapping(text: str) -> dict[str, str]:
    mapping = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        if "=" not in line:
            raise ValueError(f"Expected sample=group: {line}")
        sample, group = line.split("=", 1)
        mapping[sample.strip()] = group.strip()
    return mapping



def summarize_numeric(data: pd.DataFrame, groups: Sequence[str], value: str) -> pd.DataFrame:
    if data.empty:
        return pd.DataFrame(columns=[*groups, "mean_value", "median_value", "n_cells"])
    grouped = data.groupby(list(groups), observed=True, dropna=False)
    summary = grouped[value].agg(mean_value="mean", median_value="median").reset_index()
    if "component" in data:
        counts = grouped["component"].nunique().rename("n_cells").reset_index()
    else:
        counts = grouped.size().rename("n_cells").reset_index()
    return summary.merge(counts, on=list(groups), how="left")



def adjust_bh(values: Sequence[float]) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    adjusted = np.full(values.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    if not len(valid):
        return adjusted
    order = valid[np.argsort(values[valid])]
    ranked = values[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.minimum(ranked, 1)
    return adjusted


def calculate_differential(
    data: pd.DataFrame,
    *,
    feature_cols: Sequence[str],
    value_col: str,
    group_a: str,
    group_b: str,
    celltypes: Sequence[str] | None = None,
    stratify: bool = False,
    min_observations: int = 3,
    group_col: str = "condition",
    celltype_col: str = "celltype_manual",
) -> pd.DataFrame:
    rows = data[data[group_col].isin([group_a, group_b])].copy()
    if celltypes and celltype_col in rows:
        rows = rows[rows[celltype_col].isin(celltypes)]
    if celltype_col not in rows:
        rows[celltype_col] = "All cells"
    if not stratify:
        rows[celltype_col] = "Pooled cell types"
    groups = [*feature_cols, celltype_col]
    results = []
    for keys, chunk in rows.groupby(groups, observed=True, dropna=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        a = pd.to_numeric(chunk.loc[chunk[group_col] == group_a, value_col], errors="coerce").dropna()
        b = pd.to_numeric(chunk.loc[chunk[group_col] == group_b, value_col], errors="coerce").dropna()
        p_value = np.nan
        if len(a) >= min_observations and len(b) >= min_observations and pd.concat([a, b]).nunique() > 1:
            p_value = mannwhitneyu(b, a, alternative="two-sided").pvalue
        row = dict(zip(groups, keys))
        row.update(
            group_a=group_a,
            group_b=group_b,
            mean_a=a.mean() if len(a) else np.nan,
            mean_b=b.mean() if len(b) else np.nan,
            median_a=a.median() if len(a) else np.nan,
            median_b=b.median() if len(b) else np.nan,
            n_a=len(a),
            n_b=len(b),
            effect_size=(a.median() - b.median()) if len(a) and len(b) else np.nan,
            p_value=p_value,
        )
        results.append(row)
    result = pd.DataFrame(results)
    if result.empty:
        return result
    result["p_adj"] = adjust_bh(result["p_value"])
    result["direction"] = np.select(
        [result["effect_size"] > 0, result["effect_size"] < 0],
        [f"Higher in {group_a}", f"Higher in {group_b}"],
        default="No change",
    )
    return result.sort_values(["p_adj", "effect_size"], ascending=[True, False], na_position="last")
