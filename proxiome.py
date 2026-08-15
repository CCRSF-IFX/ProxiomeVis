"""Data and analysis helpers for the Python ProxiomeVis app."""

from __future__ import annotations

import glob
import gzip
import hashlib
import os
import pickle
import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


DEFAULT_MARKERS = (
    "CD3e", "CD4", "CD8", "CD45", "CD81", "CD82", "CD53", "CD58",
    "CD69", "CD25", "CD279", "CD28", "CD40", "CD44", "HLA-DR", "B2M",
)

PATCH_TABLES = {
    "run_plan": ("patch_detection_run_plan.csv",),
    "marker_unmixing": ("patch_marker_unmixing_table.csv",),
    "raji_marker_abundance": ("raji_marker_abundance_sanity_check.csv",),
    "raji_marker_proximity": ("raji_marker_joint_proximity_sanity_check.csv",),
    "patch_burden": ("patch_burden_by_cart_cell.csv", "patch_burden_by_cd8t_cell.csv"),
}


@dataclass(frozen=True)
class AppData:
    source: dict
    marker_options: tuple[str, ...]
    metadata: pd.DataFrame
    abundance: pd.DataFrame
    proximity: pd.DataFrame
    clustering: pd.DataFrame
    colocalization: pd.DataFrame
    qc_origin: pd.DataFrame
    qc_filtered: pd.DataFrame
    qc_filter_counts: pd.DataFrame
    patch: dict[str, pd.DataFrame | None]
    pxl_files: tuple[str, ...]


def resolve_pxl_paths(spec: str | Path | Iterable[str | Path]) -> list[Path]:
    """Resolve files, directories, globs, or newline/comma-separated PXL paths."""
    if isinstance(spec, (str, Path)):
        parts = re.split(r"[\n,]+", str(spec))
    else:
        parts = [str(value) for value in spec]

    paths: list[Path] = []
    for raw in parts:
        value = os.path.expandvars(os.path.expanduser(raw.strip()))
        if not value:
            continue
        path = Path(value)
        if path.is_dir():
            paths.extend(sorted(path.glob("*.pxl")))
        elif glob.has_magic(value):
            paths.extend(Path(match) for match in sorted(glob.glob(value)))
        else:
            paths.append(path)

    unique = list(dict.fromkeys(path.resolve() for path in paths))
    missing = [str(path) for path in unique if not path.is_file()]
    invalid = [str(path) for path in unique if path.suffix != ".pxl"]
    if not unique:
        raise ValueError("No .pxl files matched the supplied path.")
    if missing:
        raise FileNotFoundError("PXL file(s) not found: " + ", ".join(missing))
    if invalid:
        raise ValueError("Expected .pxl files: " + ", ".join(invalid))
    return unique


def default_pxl_spec() -> str:
    configured = os.getenv("PROXIOME_PXL", "").strip()
    if configured:
        return configured
    data_dir = Path(__file__).resolve().parents[2] / "data"
    return str(data_dir / "PNA065*.layout.pxl")


def cache_dir() -> Path:
    root = Path(os.getenv("PROXIOMEVIS_HOME", Path.home() / ".ProxiomeVis"))
    path = root / "python-cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _cache_path(
    paths: Sequence[Path],
    markers: Sequence[str] | None,
    annotation_path: str | Path | None,
) -> Path:
    signature = []
    for path in paths:
        stat = path.stat()
        signature.append((str(path), stat.st_size, stat.st_mtime_ns))
    if annotation_path:
        sidecar = Path(annotation_path).expanduser().resolve()
        stat = sidecar.stat()
        signature.append((str(sidecar), stat.st_size, stat.st_mtime_ns))
    signature.append(tuple(markers or DEFAULT_MARKERS))
    digest = hashlib.sha256(repr(signature).encode()).hexdigest()[:20]
    return cache_dir() / f"proxiome-{digest}.pkl.gz"


def load_pxl_data(
    spec: str | Path | Iterable[str | Path],
    *,
    markers: Sequence[str] | None = None,
    annotation_path: str | Path | None = None,
    patch_dir: str | Path | None = None,
    force: bool = False,
) -> AppData:
    """Load and aggregate one or more PXL files using Pixelator's public Python API."""
    paths = resolve_pxl_paths(spec)
    cached = _cache_path(paths, markers, annotation_path)
    if cached.exists() and not force:
        with gzip.open(cached, "rb") as handle:
            data = pickle.load(handle)
        if isinstance(data, AppData):
            return replace(data, patch=load_patch_tables(patch_dir))

    from pixelator import read_pna

    dataset = read_pna(paths)
    available = sorted(dataset.markers())
    selected = select_markers(available, markers)
    filtered = dataset.filter(markers=selected)
    adata = filtered.adata(add_log1p_transform=True, add_clr_transform=True)
    metadata = build_metadata(adata, annotation_path=annotation_path)
    abundance = build_abundance(adata)
    proximity = filtered.proximity().to_df()
    proximity = normalize_proximity(proximity, metadata)
    clustering, colocalization = split_proximity(proximity)
    qc_counts = build_qc_filter_counts(metadata)
    data = AppData(
        source={
            "display_name": f"{len(paths)} PXL file(s)",
            "source_type": "pxl",
            "n_cells": len(metadata),
            "n_markers": len(available),
            "selected_markers": len(selected),
            "pxl_paths": [str(path) for path in paths],
            "analysis_group_label": "condition",
        },
        marker_options=tuple(selected),
        metadata=metadata,
        abundance=abundance,
        proximity=proximity,
        clustering=clustering,
        colocalization=colocalization,
        qc_origin=metadata.copy(),
        qc_filtered=metadata.copy(),
        qc_filter_counts=qc_counts,
        patch=load_patch_tables(patch_dir),
        pxl_files=tuple(str(path) for path in paths),
    )
    with gzip.open(cached, "wb", compresslevel=3) as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return data


def select_markers(available: Sequence[str], requested: Sequence[str] | None = None) -> list[str]:
    available = list(dict.fromkeys(map(str, available)))
    requested = list(requested or DEFAULT_MARKERS)
    selected = [marker for marker in requested if marker in available]
    return selected or available[: min(16, len(available))]


def _read_annotations(path: str | Path | None) -> pd.DataFrame | None:
    if not path:
        return None
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Annotation file not found: {path}")
    annotations = pd.read_csv(path)
    if "component" not in annotations.columns:
        annotations = annotations.rename(columns={annotations.columns[0]: "component"})
    if "cell_type" in annotations.columns and "celltype_manual" not in annotations.columns:
        annotations = annotations.rename(columns={"cell_type": "celltype_manual"})
    return annotations.drop_duplicates("component")


def build_metadata(adata, annotation_path: str | Path | None = None) -> pd.DataFrame:
    metadata = adata.obs.copy()
    metadata.insert(0, "component", metadata.index.astype(str))
    metadata = metadata.reset_index(drop=True)
    if "sample" not in metadata:
        metadata["sample"] = "sample"
    metadata["sample"] = metadata["sample"].fillna("sample").astype(str)
    if "sample_alias" not in metadata:
        metadata["sample_alias"] = metadata["sample"]
    if "condition" not in metadata:
        metadata["condition"] = metadata["sample"]
    if "celltype_manual" not in metadata:
        metadata["celltype_manual"] = "unannotated"

    annotations = _read_annotations(annotation_path)
    if annotations is not None:
        columns = [column for column in annotations if column != "component"]
        metadata = metadata.merge(annotations, on="component", how="left", suffixes=("", "_annotation"))
        for column in columns:
            annotation_col = f"{column}_annotation"
            if annotation_col in metadata:
                metadata[column] = metadata[annotation_col].combine_first(metadata.get(column))
                metadata = metadata.drop(columns=annotation_col)
    metadata["condition"] = metadata["condition"].fillna(metadata["sample"]).astype(str)
    metadata["celltype_manual"] = metadata["celltype_manual"].fillna("unannotated").astype(str)
    metadata["sample_alias"] = metadata["sample_alias"].fillna(metadata["sample"]).astype(str)

    embedding_names = ("umap", "pca", "harmony", "tsne")
    embeddings = {
        key.removeprefix("X_"): np.asarray(value)
        for key, value in adata.obsm.items()
        if np.asarray(value).ndim == 2
        and np.asarray(value).shape[1] >= 2
        and any(name in key.lower() for name in embedding_names)
    }
    if not embeddings:
        embeddings["umap"] = compute_umap(adata)
    for name, values in embeddings.items():
        for index in range(min(values.shape[1], 30)):
            metadata[f"{name}_{index + 1}"] = values[:, index]
    return metadata


def compute_umap(adata) -> np.ndarray:
    n_cells = adata.n_obs
    if n_cells < 3:
        return np.column_stack((np.arange(n_cells), np.zeros(n_cells)))
    if adata.n_vars < 2:
        values = np.asarray(adata.layers["clr"] if "clr" in adata.layers else adata.X).reshape(n_cells, -1)
        return np.column_stack((values[:, 0], np.zeros(n_cells)))
    import scanpy as sc

    n_pcs = min(30, adata.n_vars - 1, n_cells - 1)
    values = adata.layers["clr"] if "clr" in adata.layers else adata.X
    adata.obsm["X_pca"] = sc.tl.pca(values, n_comps=n_pcs, random_state=42)
    sc.pp.neighbors(
        adata,
        n_neighbors=min(30, n_cells - 1),
        n_pcs=min(10, n_pcs),
        use_rep="X_pca",
        metric="cosine",
    )
    sc.tl.umap(adata, min_dist=0.3, random_state=42)
    return np.asarray(adata.obsm["X_umap"])


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


def normalize_proximity(proximity: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    proximity = proximity.copy()
    required = {"component", "marker_1", "marker_2", "log2_ratio"}
    missing = required.difference(proximity.columns)
    if missing:
        raise ValueError("PXL proximity data is missing: " + ", ".join(sorted(missing)))
    proximity["component"] = proximity["component"].astype(str)
    meta_cols = ["component", "condition", "celltype_manual", "sample_alias", "sample"]
    for column in meta_cols[1:]:
        proximity = proximity.drop(columns=column, errors="ignore")
    proximity = proximity.merge(metadata[meta_cols], on="component", how="left")
    proximity["log2_ratio"] = pd.to_numeric(proximity["log2_ratio"], errors="coerce")
    return proximity


def split_proximity(proximity: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    self_mask = proximity["marker_1"].astype(str) == proximity["marker_2"].astype(str)
    clustering = proximity.loc[self_mask].copy()
    clustering["marker"] = clustering["marker_1"].astype(str)
    colocalization = proximity.loc[~self_mask].copy()
    colocalization["marker_pair"] = (
        colocalization["marker_1"].astype(str) + " / " + colocalization["marker_2"].astype(str)
    )
    return clustering, colocalization


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
        proximity=regroup(data.proximity),
        clustering=regroup(data.clustering),
        colocalization=regroup(data.colocalization),
        qc_origin=regroup(data.qc_origin),
        qc_filtered=regroup(data.qc_filtered),
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


def filter_pixelator_proximity(
    proximity: pd.DataFrame,
    *,
    min_marker_fraction: float = 0,
    min_marker_count: float = 0,
    min_cells: int = 1,
) -> pd.DataFrame:
    rows = proximity.copy()
    if min_marker_fraction > 0:
        if {"marker_1_freq", "marker_2_freq"}.issubset(rows):
            rows = rows[
                (pd.to_numeric(rows["marker_1_freq"], errors="coerce") >= min_marker_fraction)
                & (pd.to_numeric(rows["marker_2_freq"], errors="coerce") >= min_marker_fraction)
            ]
        else:
            raise ValueError("Marker-fraction filtering requires marker_1_freq and marker_2_freq.")
    if min_marker_count > 0:
        if "min_count" not in rows:
            raise ValueError("Marker-count filtering requires min_count.")
        rows = rows[pd.to_numeric(rows["min_count"], errors="coerce") >= min_marker_count]
    if min_cells > 1 and len(rows):
        counts = rows.groupby(["marker_1", "marker_2"], observed=True)["component"].transform("nunique")
        rows = rows[counts >= int(min_cells)]
    return rows


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


def summarize_spatial(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    group_col: str,
    markers: Sequence[str],
    celltypes: Sequence[str] | None = None,
    mean_type: str = "population",
) -> pd.DataFrame:
    meta = metadata.copy()
    if celltypes:
        meta = meta[meta["celltype_manual"].isin(celltypes)]
    rows = scores[
        scores["component"].isin(meta["component"])
        & scores["marker_1"].isin(markers)
        & scores["marker_2"].isin(markers)
    ].copy()
    group_cols = [group_col, "marker_1", "marker_2"]
    detected = rows.groupby(group_cols, observed=True, dropna=False).agg(
        sum_log2_ratio=("log2_ratio", "sum"),
        detected_mean=("log2_ratio", "mean"),
        n_detected=("component", "nunique"),
    ).reset_index()
    totals = meta.groupby(group_col, observed=True)["component"].nunique().rename("n_total").reset_index()
    summary = detected.merge(totals, on=group_col, how="left")
    summary["pct_detected"] = summary["n_detected"] / summary["n_total"]
    summary["mean_log2_ratio"] = np.where(
        mean_type == "detected",
        summary["detected_mean"],
        summary["sum_log2_ratio"] / summary["n_total"],
    )
    return complete_spatial(summary, group_col, markers, totals)


def complete_spatial(
    summary: pd.DataFrame,
    group_col: str,
    markers: Sequence[str],
    totals: pd.DataFrame,
) -> pd.DataFrame:
    if summary.empty:
        return summary
    reverse = summary.rename(columns={"marker_1": "marker_2", "marker_2": "marker_1"})
    summary = pd.concat([summary, reverse], ignore_index=True).drop_duplicates(
        [group_col, "marker_1", "marker_2"], keep="first"
    )
    grid = pd.MultiIndex.from_product(
        [totals[group_col].astype(str), markers, markers], names=[group_col, "marker_1", "marker_2"]
    ).to_frame(index=False)
    result = grid.merge(summary, on=[group_col, "marker_1", "marker_2"], how="left")
    result = result.merge(totals, on=group_col, how="left", suffixes=("", "_grid"))
    result["n_total"] = result["n_total"].fillna(result.pop("n_total_grid"))
    for column in ("sum_log2_ratio", "n_detected", "pct_detected"):
        result[column] = result[column].fillna(0)
    return result


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


def sample_colocalization(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    celltypes: Sequence[str] | None = None,
    mean_type: str = "population",
) -> pd.DataFrame:
    summary = summarize_spatial(
        scores,
        metadata,
        group_col="sample_alias",
        markers=sorted(set(scores["marker_1"]).union(scores["marker_2"])),
        celltypes=celltypes,
        mean_type=mean_type,
    )
    sample_groups = metadata[["sample_alias", "condition"]].drop_duplicates("sample_alias")
    summary = summary.merge(sample_groups, on="sample_alias", how="left")
    summary["marker_pair"] = summary["marker_1"].astype(str) + " / " + summary["marker_2"].astype(str)
    summary["pair_observed"] = summary["n_detected"] > 0
    return summary


def load_patch_tables(patch_dir: str | Path | None) -> dict[str, pd.DataFrame | None]:
    configured = str(patch_dir or os.getenv("PROXIOME_PATCH_DIR", "")).strip()
    result: dict[str, pd.DataFrame | None] = {key: None for key in PATCH_TABLES}
    if not configured:
        return result
    root = Path(configured).expanduser().resolve()
    tables = root / "tables" if (root / "tables").is_dir() else root
    for key, filenames in PATCH_TABLES.items():
        for filename in filenames:
            path = tables / filename
            if path.is_file():
                result[key] = pd.read_csv(path)
                break
    return result


def read_component_layout(data: AppData, sample: str, component: str) -> pd.DataFrame:
    from pixelator import read_pna

    candidates = [Path(path) for path in data.pxl_files if sample in Path(path).stem]
    if not candidates:
        raise FileNotFoundError(f"No PXL file is associated with sample {sample}.")
    dataset = read_pna(candidates[0])
    available = dataset.components()
    raw_component = component
    prefix = f"{sample}_"
    if raw_component not in available and raw_component.startswith(prefix):
        raw_component = raw_component[len(prefix):]
    if raw_component not in available:
        raise ValueError(f"Component {component} is not present in {candidates[0].name}.")
    layout = dataset.filter(components=[raw_component]).precomputed_layouts(add_marker_counts=False).to_df()
    if layout.empty:
        raise ValueError("No precomputed 3D layout is stored for this component.")
    return layout
