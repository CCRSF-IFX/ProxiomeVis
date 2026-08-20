"""Data and analysis helpers for the Python ProxiomeVis app."""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field, replace
from functools import lru_cache
from io import BytesIO
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from typing import Literal, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


REFERENCE_H5AD = Path(
    "/Volumes/ccrsf-static/illumina/CCRSFIFX-23_MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/"
    "python_results/pg_data_combined_filtered_annotated.h5ad"
)
REFERENCE_PXL = (
    "/Volumes/ccrsf-static/singlecell_projects/MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/"
    "Analysis_2nd_combo/Analysis/pixelator/*.pxl"
)

PATCH_TABLES = (
    "run_plan",
    "marker_unmixing",
    "target_marker_abundance",
    "target_marker_proximity",
    "patch_sizes",
    "raji_marker_abundance",
    "raji_marker_proximity",
    "patch_burden",
    "patch_composition",
)


def empty_patch() -> dict[str, pd.DataFrame | None]:
    return {name: None for name in PATCH_TABLES}



@dataclass(frozen=True)
class AppData:
    source: dict
    marker_options: tuple[str, ...]
    metadata: pd.DataFrame
    abundance: pd.DataFrame
    qc_filter_counts: pd.DataFrame
    patch: dict[str, pd.DataFrame | None] = field(default_factory=empty_patch)
    component_layouts: dict[str, pd.DataFrame] = field(default_factory=dict)
    pxl_files: tuple[str, ...] = ()


@dataclass(frozen=True)
class SpatialRetrieval:
    retrieved_at: str
    request: tuple
    metadata: pd.DataFrame
    markers: tuple[str, ...]
    marker_mode: str
    self_proximity: pd.DataFrame


@dataclass(frozen=True)
class PatchRetrieval:
    prepared_at: str
    request: tuple
    metadata: pd.DataFrame
    candidate_data: pd.DataFrame
    target_markers: tuple[str, ...]
    receiver_markers: tuple[str, ...]


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


def default_pxl_spec() -> str:
    return os.getenv("PROXIOME_PXL", REFERENCE_PXL).strip()


def dataframe_export_bytes(frame: pd.DataFrame, file_type: Literal["csv", "xlsx"]) -> bytes:
    """Serialize a complete table for browser download."""
    if file_type == "csv":
        return frame.to_csv(index=False).encode("utf-8-sig")
    if file_type != "xlsx":
        raise ValueError("Table export type must be 'csv' or 'xlsx'.")
    if len(frame) > 1_048_575 or len(frame.columns) > 16_384:
        raise ValueError("This table exceeds Excel's worksheet limits; export CSV instead.")
    payload = BytesIO()
    with pd.ExcelWriter(payload, engine="openpyxl") as writer:
        frame.to_excel(writer, index=False, sheet_name="Data")
    return payload.getvalue()


def resolve_pxl_paths(spec: str | Path | None) -> list[Path]:
    """Resolve PXL files used for proximity queries and component layouts."""
    if not spec or not str(spec).strip():
        return []
    paths = []
    for raw in str(spec).replace(",", "\n").splitlines():
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
    paths = list(dict.fromkeys(path.resolve() for path in paths))
    missing = [str(path) for path in paths if not path.is_file()]
    invalid = [str(path) for path in paths if path.suffix.lower() != ".pxl"]
    if missing:
        raise FileNotFoundError("PXL file(s) not found: " + ", ".join(missing))
    if invalid:
        raise ValueError("Expected .pxl files: " + ", ".join(invalid))
    return paths


def load_h5ad_data(spec: str | Path, *, pxl_spec: str | Path | None = None) -> AppData:
    """Load cell data from H5AD and assign PXL files for spatial queries."""
    from anndata import read_h5ad

    path = resolve_h5ad_path(spec)
    pxl_files = tuple(map(str, resolve_pxl_paths(pxl_spec)))
    adata = read_h5ad(path)
    if not adata.n_obs or not adata.n_vars:
        raise ValueError("H5AD must contain at least one observation and one marker.")

    # Proximity is intentionally sourced from PXL. Dropping an embedded copy
    # here prevents later helpers from accidentally expanding it in memory.
    proxiome_payload = adata.uns.get("proxiome")
    if isinstance(proxiome_payload, dict):
        proxiome_payload.pop("proximity", None)

    metadata = build_metadata(adata)
    markers = tuple(map(str, adata.var_names))
    patch = load_h5ad_patch(adata)
    component_layouts = load_h5ad_component_layouts(adata)
    return AppData(
        source={
            "display_name": path.name,
            "source_type": "h5ad",
            "n_cells": int(adata.n_obs),
            "n_markers": int(adata.n_vars),
            "h5ad_path": str(path),
            "analysis_group_label": "condition",
            "has_spatial_metrics": bool(pxl_files),
            "has_patch_analysis": any(
                patch.get(name) is not None and not patch[name].empty
                for name in ("patch_sizes", "patch_burden", "patch_composition")
            ),
            "has_component_layouts": bool(component_layouts),
            "pxl_files": len(pxl_files),
        },
        marker_options=markers,
        metadata=metadata,
        abundance=build_abundance(adata),
        qc_filter_counts=build_h5ad_qc_filter_counts(adata, metadata),
        patch=patch,
        component_layouts=component_layouts,
        pxl_files=pxl_files,
    )


def h5ad_proxiome_payload(adata) -> Mapping:
    payload = adata.uns.get("proxiome")
    return payload if isinstance(payload, Mapping) else adata.uns


def payload_frame(payload: Mapping, name: str) -> pd.DataFrame | None:
    value = payload.get(name)
    if value is None:
        return None
    if isinstance(value, pd.DataFrame):
        return value.copy()
    if isinstance(value, Mapping):
        return pd.DataFrame(value)
    raise ValueError(f'H5AD uns["proxiome"]["{name}"] must be a table.')


@lru_cache(maxsize=4)
def _read_pxl_dataset(pxl_files: tuple[str, ...]):
    try:
        from pixelator import read_pna
    except ImportError as error:
        raise RuntimeError("PXL proximity queries require the pixelgen-pixelator package.") from error
    paths = [Path(path) for path in pxl_files]
    return read_pna(paths[0] if len(paths) == 1 else paths)


@lru_cache(maxsize=4)
def _pxl_component_ids(pxl_files: tuple[str, ...]) -> frozenset[str]:
    return frozenset(map(str, _read_pxl_dataset(pxl_files).components()))


def clear_pxl_cache() -> None:
    """Clear cached PXL readers (primarily useful after replacing a file)."""
    _pxl_component_ids.cache_clear()
    _read_pxl_dataset.cache_clear()


def filter_pxl_metadata(data: AppData, metadata: pd.DataFrame | None = None) -> pd.DataFrame:
    """Keep H5AD observations that exist in the assigned PXL files."""
    metadata = data.metadata if metadata is None else metadata
    if not data.pxl_files or metadata.empty:
        return metadata.iloc[0:0].copy()
    components = _pxl_component_ids(data.pxl_files)
    return metadata[metadata["component"].astype(str).isin(components)].copy()


def build_spatial_retrieval(
    data: AppData,
    *,
    conditions: Sequence[str],
    samples: Sequence[str],
    celltypes: Sequence[str],
    marker_mode: Literal["all", "top", "manual"] = "all",
    n_markers: int = 40,
    plot_markers: Sequence[str] | None = None,
    retrieved_at: str,
    request: tuple,
) -> SpatialRetrieval:
    """Freeze the component and marker scope used by every spatial view."""
    if marker_mode not in {"all", "top", "manual"}:
        raise ValueError("Marker mode must be 'all', 'top', or 'manual'.")
    conditions = set(map(str, conditions))
    samples = set(map(str, samples))
    celltypes = set(map(str, celltypes))
    metadata = data.metadata[
        data.metadata["condition"].astype(str).isin(conditions)
        & data.metadata["sample_alias"].astype(str).isin(samples)
        & data.metadata["celltype_manual"].astype(str).isin(celltypes)
    ].copy()
    metadata = filter_pxl_metadata(data, metadata)
    if metadata.empty:
        raise ValueError("No H5AD components in the selected population are present in the assigned PXL files.")

    if marker_mode == "all":
        markers = list(map(str, data.marker_options))
    elif marker_mode == "manual":
        requested = set(map(str, plot_markers or ()))
        markers = [str(marker) for marker in data.marker_options if str(marker) in requested]
    else:
        markers = select_colocalization_heatmap_markers(
            data.abundance,
            metadata,
            data.marker_options,
            n_markers=n_markers,
        )
    if not markers:
        raise ValueError("Select at least one marker for spatial retrieval.")
    self_proximity = load_pxl_proximity(
        data,
        metadata,
        markers=markers,
        pair_type="self",
    )
    if self_proximity.empty:
        raise ValueError("No PXL self-proximity rows match the selected retrieval scope.")
    self_proximity["marker"] = self_proximity["marker_1"].astype(str)
    self_proximity = self_proximity.drop(columns=["marker_1", "marker_2"])
    for column in ("marker", "condition", "celltype_manual", "sample_alias", "sample"):
        self_proximity[column] = self_proximity[column].astype("category")
    return SpatialRetrieval(
        retrieved_at=retrieved_at,
        request=request,
        metadata=metadata,
        markers=tuple(markers),
        marker_mode=marker_mode,
        self_proximity=self_proximity,
    )


def load_pxl_proximity(
    data: AppData,
    metadata: pd.DataFrame | None = None,
    *,
    markers: Sequence[str] | None = None,
    pair_type: Literal["all", "self", "nonself"] = "all",
    anchor: str | None = None,
    add_marker_counts: bool = False,
) -> pd.DataFrame:
    """Query selected proximity rows from assigned PXL files."""
    if pair_type not in {"all", "self", "nonself"}:
        raise ValueError("pair_type must be 'all', 'self', or 'nonself'.")
    if not data.pxl_files:
        return pd.DataFrame()

    metadata = filter_pxl_metadata(data, metadata)
    components = metadata["component"].dropna().astype(str).drop_duplicates().tolist()
    if not components:
        return pd.DataFrame()
    selected_markers = None if markers is None else list(dict.fromkeys(map(str, markers)))
    if selected_markers == []:
        return pd.DataFrame()

    dataset = _read_pxl_dataset(data.pxl_files)
    if add_marker_counts:
        filtered = dataset.filter(components=components, markers=selected_markers)
        proximity = filtered.proximity(
            add_marker_counts=True,
            add_logratio=True,
            calculate_from_edgelist=False,
        ).to_df()
        if pair_type == "self":
            proximity = proximity[proximity["marker_1"] == proximity["marker_2"]]
        elif pair_type == "nonself":
            proximity = proximity[proximity["marker_1"] != proximity["marker_2"]]
        if anchor:
            proximity = proximity[
                proximity["marker_1"].eq(str(anchor)) | proximity["marker_2"].eq(str(anchor))
            ]
    else:
        conditions = ["component IN $components"]
        parameters: dict[str, object] = {"components": components}
        if selected_markers is not None:
            conditions.append("marker_1 IN $markers AND marker_2 IN $markers")
            parameters["markers"] = selected_markers
        if pair_type == "self":
            conditions.append("marker_1 = marker_2")
        elif pair_type == "nonself":
            conditions.append("marker_1 != marker_2")
        if anchor:
            conditions.append("(marker_1 = $anchor OR marker_2 = $anchor)")
            parameters["anchor"] = str(anchor)
        query = f"""
            SELECT component, marker_1, marker_2,
                   log2(
                       greatest(CAST(join_count AS DOUBLE), 1) /
                       greatest(join_count_expected_mean, 1)
                   ) AS log2_ratio
            FROM proximity
            WHERE {' AND '.join(conditions)}
        """
        with dataset.view.open() as session:
            proximity = session.get_connection().execute(query, parameters).fetchdf()

    required = {"component", "marker_1", "marker_2", "log2_ratio"}
    if missing := required.difference(proximity.columns):
        raise ValueError("PXL proximity table is missing: " + ", ".join(sorted(missing)))
    proximity = proximity.copy()
    for column in ("component", "marker_1", "marker_2"):
        proximity[column] = proximity[column].astype(str)
    proximity["log2_ratio"] = pd.to_numeric(proximity["log2_ratio"], errors="coerce")
    meta_columns = ["component", "condition", "celltype_manual", "sample_alias", "sample"]
    proximity = proximity.drop(columns=meta_columns[1:], errors="ignore")
    return proximity.merge(
        metadata[meta_columns].drop_duplicates("component"),
        on="component",
        how="inner",
        validate="many_to_one",
    )


def sample_pxl_colocalization(
    data: AppData,
    metadata: pd.DataFrame,
    *,
    markers: Sequence[str] | None = None,
    mean_type: str = "population",
    anchor: str | None = None,
    pair_type: Literal["all", "self", "nonself"] = "nonself",
    min_marker_fraction: float = 0,
    min_marker_count: float = 0,
) -> pd.DataFrame:
    """Aggregate PXL colocalization in DuckDB before returning sample rows."""
    metadata = filter_pxl_metadata(data, metadata)
    if metadata.empty:
        return pd.DataFrame()
    components = metadata["component"].dropna().astype(str).drop_duplicates().tolist()
    if pair_type not in {"all", "self", "nonself"}:
        raise ValueError("pair_type must be 'all', 'self', or 'nonself'.")
    conditions = ["proximity.component IN $components"]
    if pair_type == "self":
        conditions.append("marker_1 = marker_2")
    elif pair_type == "nonself":
        conditions.append("marker_1 != marker_2")
    parameters: dict[str, object] = {"components": components}
    selected_markers = list(dict.fromkeys(map(str, data.marker_options if markers is None else markers)))
    if not selected_markers:
        return pd.DataFrame()
    conditions.append("marker_1 IN $markers AND marker_2 IN $markers")
    parameters["markers"] = selected_markers
    if anchor:
        conditions.append("(marker_1 = $anchor OR marker_2 = $anchor)")
        parameters["anchor"] = str(anchor)
    eligibility_joins = ""
    eligibility = None
    if min_marker_fraction > 0 or min_marker_count > 0:
        counts = data.abundance[["component", "marker", "count"]].copy()
        counts["component"] = counts["component"].astype(str)
        counts["marker"] = counts["marker"].astype(str)
        counts["count"] = pd.to_numeric(counts["count"], errors="coerce").fillna(0)
        counts = counts[counts["component"].isin(components)]
        totals_by_cell = counts.groupby("component", observed=True)["count"].transform("sum")
        counts["marker_fraction"] = counts["count"] / totals_by_cell.replace(0, np.nan)
        eligibility = counts[
            counts["marker"].isin(selected_markers)
            & (counts["count"] >= float(min_marker_count))
            & (counts["marker_fraction"] >= float(min_marker_fraction))
        ][["component", "marker"]].drop_duplicates()
        eligibility_joins = """
            INNER JOIN marker_eligibility AS marker_1_ok
                ON proximity.component = marker_1_ok.component
               AND proximity.marker_1 = marker_1_ok.marker
            INNER JOIN marker_eligibility AS marker_2_ok
                ON proximity.component = marker_2_ok.component
               AND proximity.marker_2 = marker_2_ok.marker
        """
    query = f"""
        WITH selected AS (
            SELECT sample, marker_1, marker_2, proximity.component,
                   log2(
                       greatest(CAST(join_count AS DOUBLE), 1) /
                       greatest(join_count_expected_mean, 1)
                   ) AS log2_ratio
            FROM proximity
            {eligibility_joins}
            WHERE {' AND '.join(conditions)}
        )
        SELECT sample, marker_1, marker_2,
               sum(log2_ratio) AS sum_log2_ratio,
               avg(log2_ratio) AS detected_mean,
               count(DISTINCT component) AS n_detected
        FROM selected
        GROUP BY sample, marker_1, marker_2
    """
    dataset = _read_pxl_dataset(data.pxl_files)
    with dataset.view.open() as session:
        connection = session.get_connection()
        if eligibility is not None:
            connection.register("marker_eligibility", eligibility)
        detected = connection.execute(query, parameters).fetchdf()

    sample_metadata = metadata[["sample", "sample_alias", "condition"]].drop_duplicates("sample").copy()
    sample_metadata["sample"] = sample_metadata["sample"].astype(str)
    totals = (
        metadata.groupby(["sample", "sample_alias", "condition"], observed=True)["component"]
        .nunique()
        .rename("n_total")
        .reset_index()
    )
    totals["sample"] = totals["sample"].astype(str)
    pair_builder = combinations_with_replacement if pair_type == "all" else combinations
    marker_pairs = [
        pair
        for pair in pair_builder(sorted(selected_markers), 2)
        if (pair_type != "self" or pair[0] == pair[1]) and (not anchor or anchor in pair)
    ]
    if not marker_pairs or totals.empty:
        return pd.DataFrame()
    pair_grid = pd.DataFrame(marker_pairs, columns=["marker_1", "marker_2"])
    grid = totals.merge(pair_grid, how="cross")

    if not detected.empty:
        detected["sample"] = detected["sample"].astype(str)
        detected = detected.merge(sample_metadata, on="sample", how="inner")
        first = detected["marker_1"].astype(str)
        second = detected["marker_2"].astype(str)
        detected["marker_1"] = first.where(first <= second, second)
        detected["marker_2"] = second.where(first <= second, first)
        detected = detected.groupby(
            ["sample", "sample_alias", "condition", "marker_1", "marker_2"],
            observed=True,
            as_index=False,
        ).agg(
            sum_log2_ratio=("sum_log2_ratio", "sum"),
            n_detected=("n_detected", "sum"),
        )
        # PXL stores one orientation per unordered marker pair.
        detected["detected_mean"] = detected["sum_log2_ratio"] / detected["n_detected"]
    else:
        detected = pd.DataFrame(
            columns=[
                "sample",
                "sample_alias",
                "condition",
                "marker_1",
                "marker_2",
                "sum_log2_ratio",
                "detected_mean",
                "n_detected",
            ]
        )
    result = grid.merge(
        detected.drop(columns=["n_total"], errors="ignore") if not detected.empty else detected,
        on=["sample", "sample_alias", "condition", "marker_1", "marker_2"],
        how="left",
    )
    for column in ("sum_log2_ratio", "n_detected"):
        result[column] = pd.to_numeric(result[column], errors="coerce").fillna(0)
    result["pct_detected"] = result["n_detected"] / result["n_total"]
    result["mean_log2_ratio"] = np.where(
        mean_type == "detected",
        result["detected_mean"],
        result["sum_log2_ratio"] / result["n_total"],
    )
    result["marker_pair"] = result["marker_1"] + " / " + result["marker_2"]
    result["pair_observed"] = result["n_detected"] > 0
    return result


def joint_marker_proximity(
    data: AppData,
    metadata: pd.DataFrame,
    markers: Sequence[str],
) -> pd.DataFrame:
    """Calculate the joint proximity score for one marker set per component."""
    metadata = filter_pxl_metadata(data, metadata)
    selected_markers = list(dict.fromkeys(map(str, markers)))
    if metadata.empty or not selected_markers:
        return pd.DataFrame()
    components = metadata["component"].dropna().astype(str).drop_duplicates().tolist()
    query = """
        SELECT component,
               sum(CAST(join_count AS DOUBLE)) AS join_count,
               sum(join_count_expected_mean) AS join_count_expected_mean,
               log2(
                   greatest(sum(CAST(join_count AS DOUBLE)), 1) /
                   greatest(sum(join_count_expected_mean), 1)
               ) AS log2_ratio
        FROM proximity
        WHERE component IN $components
          AND marker_1 IN $markers
          AND marker_2 IN $markers
        GROUP BY component
    """
    dataset = _read_pxl_dataset(data.pxl_files)
    with dataset.view.open() as session:
        rows = session.get_connection().execute(
            query,
            {"components": components, "markers": selected_markers},
        ).fetchdf()
    if rows.empty:
        rows = pd.DataFrame(
            columns=["component", "join_count", "join_count_expected_mean", "log2_ratio"]
        )
    else:
        rows["component"] = rows["component"].astype(str)
    meta_columns = [
        column
        for column in ("component", "sample", "sample_alias", "condition", "celltype_manual")
        if column in metadata
    ]
    return metadata[meta_columns].drop_duplicates("component").merge(
        rows,
        on="component",
        how="left",
        validate="one_to_one",
    )


def build_patch_retrieval(
    data: AppData,
    request: tuple,
    *,
    prepared_at: str,
) -> PatchRetrieval:
    """Freeze one receiver/target patch scope and its candidate signal."""
    (
        conditions,
        population_col,
        receiver,
        target,
        target_markers,
        receiver_markers,
    ) = request
    if population_col not in data.metadata:
        raise ValueError("Select a valid population metadata field.")
    if not receiver or not target or receiver == target:
        raise ValueError("Select two different receiver and target populations.")
    if len(target_markers) < 2:
        raise ValueError("Select at least two target markers.")

    metadata = data.metadata[
        data.metadata["condition"].astype(str).isin(conditions)
    ].copy()
    receiver_metadata = metadata[
        metadata[population_col].astype(str) == str(receiver)
    ].copy()
    target_metadata = metadata[
        metadata[population_col].astype(str) == str(target)
    ]
    if receiver_metadata.empty or target_metadata.empty:
        raise ValueError("The active analysis groups must contain both populations.")

    if data.pxl_files:
        candidate_data = joint_marker_proximity(data, receiver_metadata, target_markers)
        counts = data.abundance[
            data.abundance["component"].astype(str).isin(
                receiver_metadata["component"].astype(str)
            )
            & data.abundance["marker"].astype(str).isin(target_markers)
        ].copy()
        counts["count"] = pd.to_numeric(counts["count"], errors="coerce").fillna(0)
        counts = (
            counts.groupby("component", observed=True)["count"]
            .sum()
            .rename("target_marker_count")
        )
        candidate_data = candidate_data.merge(counts, on="component", how="left")
        candidate_data["source"] = "Live PXL query"
    else:
        candidate_data = data.patch.get("target_marker_proximity")
        candidate_data = candidate_data.copy() if candidate_data is not None else pd.DataFrame()
        if not candidate_data.empty and "component" in candidate_data:
            candidate_data = candidate_data.rename(
                columns={"raji_marker_count": "target_marker_count"}
            )
            if "target_marker_count" not in candidate_data:
                abundance = data.patch.get("target_marker_abundance")
                abundance = abundance.copy() if abundance is not None else pd.DataFrame()
                count_col = next(
                    (
                        column
                        for column in ("target_marker_count", "raji_marker_count", "count")
                        if column in abundance
                    ),
                    None,
                )
                if "component" in abundance and count_col:
                    abundance = abundance[["component", count_col]].rename(
                        columns={count_col: "target_marker_count"}
                    )
                    abundance["component"] = abundance["component"].astype(str)
                    candidate_data["component"] = candidate_data["component"].astype(str)
                    candidate_data = candidate_data.merge(
                        abundance.drop_duplicates("component"),
                        on="component",
                        how="left",
                    )
            columns = [
                column
                for column in (
                    "component",
                    "sample",
                    "sample_alias",
                    "condition",
                    "celltype_manual",
                )
                if column in receiver_metadata
            ]
            candidate_data["component"] = candidate_data["component"].astype(str)
            candidate_data = candidate_data.drop(
                columns=columns[1:], errors="ignore"
            ).merge(
                receiver_metadata[columns].drop_duplicates("component"),
                on="component",
                how="inner",
            )
            candidate_data["source"] = "Stored H5AD table"
        else:
            candidate_data = pd.DataFrame()

    return PatchRetrieval(
        prepared_at=prepared_at,
        request=request,
        metadata=metadata,
        candidate_data=candidate_data,
        target_markers=tuple(target_markers),
        receiver_markers=tuple(receiver_markers),
    )


def summarize_sample_colocalization(
    samples: pd.DataFrame,
    *,
    group_col: str,
    markers: Sequence[str],
    mean_type: str = "population",
) -> pd.DataFrame:
    """Pool sample-level PXL summaries and return a complete symmetric matrix."""
    if samples.empty or group_col not in samples:
        return pd.DataFrame()
    grouped = samples.groupby(
        [group_col, "marker_1", "marker_2"], observed=True, dropna=False
    ).agg(
        sum_log2_ratio=("sum_log2_ratio", "sum"),
        n_detected=("n_detected", "sum"),
        n_total=("n_total", "sum"),
    ).reset_index()
    grouped["detected_mean"] = grouped["sum_log2_ratio"] / grouped["n_detected"].replace(0, np.nan)
    grouped["pct_detected"] = grouped["n_detected"] / grouped["n_total"].replace(0, np.nan)
    grouped["mean_log2_ratio"] = np.where(
        mean_type == "detected",
        grouped["detected_mean"],
        grouped["sum_log2_ratio"] / grouped["n_total"].replace(0, np.nan),
    )
    if group_col == "sample_alias":
        totals = samples[[group_col, "n_total"]].drop_duplicates(group_col)
    else:
        sample_totals = samples[["sample_alias", group_col, "n_total"]].drop_duplicates("sample_alias")
        totals = sample_totals.groupby(group_col, observed=True, as_index=False)["n_total"].sum()
    return complete_spatial(grouped, group_col, markers, totals)


def select_proximity_profile_markers(
    summary: pd.DataFrame,
    *,
    n_pairs: int = 60,
    min_detected_fraction: float = 0.5,
) -> list[str]:
    """Select markers from the strongest detected pairs, matching the notebook."""
    if summary.empty:
        return []
    rows = summary[summary["pct_detected"] > float(min_detected_fraction)].copy()
    if rows.empty:
        return []
    first = rows["marker_1"].astype(str)
    second = rows["marker_2"].astype(str)
    rows["pair_1"] = first.where(first <= second, second)
    rows["pair_2"] = second.where(first <= second, first)
    rows["absolute_mean"] = pd.to_numeric(rows["mean_log2_ratio"], errors="coerce").abs()
    rows = (
        rows.sort_values("absolute_mean", ascending=False, kind="stable")
        .drop_duplicates(["pair_1", "pair_2"])
        .head(max(1, int(n_pairs)))
    )
    return list(dict.fromkeys([*rows["pair_1"].tolist(), *rows["pair_2"].tolist()]))


def patch_marker_profiles(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    population_col: str,
    receiver_population: str,
    target_population: str,
    min_fraction: float = 0.01,
    min_enrichment: float = 3,
) -> pd.DataFrame:
    """Screen raw marker frequencies for receiver- and target-specific markers."""
    required_abundance = {"component", "marker", "count"}
    required_metadata = {"component", population_col}
    if missing := required_abundance.difference(abundance.columns):
        raise ValueError("Abundance table is missing: " + ", ".join(sorted(missing)))
    if missing := required_metadata.difference(metadata.columns):
        raise ValueError("Metadata table is missing: " + ", ".join(sorted(missing)))
    if receiver_population == target_population:
        raise ValueError("Receiver and target populations must be different.")

    population = metadata[["component", population_col]].drop_duplicates("component").copy()
    population["component"] = population["component"].astype(str)
    population[population_col] = population[population_col].astype(str)
    rows = abundance[["component", "marker", "count"]].copy()
    rows["component"] = rows["component"].astype(str)
    rows["marker"] = rows["marker"].astype(str)
    rows["count"] = pd.to_numeric(rows["count"], errors="coerce").fillna(0)
    rows = rows.merge(population, on="component", how="inner")
    rows = rows[
        rows[population_col].isin([str(receiver_population), str(target_population)])
    ]
    if rows.empty:
        return pd.DataFrame()

    counts = rows.groupby([population_col, "marker"], observed=True)["count"].sum()
    totals = counts.groupby(level=0).transform("sum").replace(0, np.nan)
    frequencies = (counts / totals).rename("frequency").reset_index()
    profiles = frequencies.pivot(index="marker", columns=population_col, values="frequency")
    profiles = profiles.reindex(
        columns=[str(receiver_population), str(target_population)], fill_value=0
    ).fillna(0)
    profiles.columns = ["receiver_freq", "target_freq"]
    profiles = profiles.reset_index()
    receiver = profiles["receiver_freq"]
    target = profiles["target_freq"]
    enrichment = max(float(min_enrichment), 1)
    threshold = max(float(min_fraction), 0)
    profiles["marker_class"] = np.select(
        [
            (target >= threshold) & (target >= receiver * enrichment),
            (receiver >= threshold) & (receiver >= target * enrichment),
        ],
        ["Target marker", "Receiver marker"],
        default="Shared / low signal",
    )
    profiles["target_enrichment"] = target / receiver.replace(0, np.nan)
    profiles["receiver_enrichment"] = receiver / target.replace(0, np.nan)
    return profiles.sort_values(
        ["marker_class", "target_freq", "receiver_freq", "marker"],
        ascending=[True, False, False, True],
        kind="stable",
    ).reset_index(drop=True)


def load_h5ad_patch(adata) -> dict[str, pd.DataFrame | None]:
    payload = h5ad_proxiome_payload(adata).get("patch", {})
    if payload is None:
        return empty_patch()
    if not isinstance(payload, Mapping):
        raise ValueError('H5AD uns["proxiome"]["patch"] must be a mapping of tables.')
    result = {
        name: payload_frame(payload, name)
        for name in PATCH_TABLES
        if name in payload
    }
    aliases = {
        "target_marker_abundance": "raji_marker_abundance",
        "target_marker_proximity": "raji_marker_proximity",
    }
    for generic, legacy in aliases.items():
        if generic not in result and legacy in result:
            result[generic] = result[legacy].copy()
    return result


def load_h5ad_component_layouts(adata) -> dict[str, pd.DataFrame]:
    payload = h5ad_proxiome_payload(adata).get("component_layouts", {})
    if payload is None:
        return {}
    if not isinstance(payload, Mapping):
        raise ValueError('H5AD uns["proxiome"]["component_layouts"] must be a mapping of tables.')
    layouts = {}
    for component in payload:
        frame = payload_frame(payload, component)
        if frame is not None:
            layouts[str(component)] = frame
    return layouts




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



def empty_qc_filter_counts() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "sample",
            "condition",
            "step",
            "step_label",
            "n_cells",
            "fraction_loaded",
        ]
    )


def build_h5ad_qc_filter_counts(adata, metadata: pd.DataFrame) -> pd.DataFrame:
    """Use QC history stored by the reference notebook when available."""
    history = adata.uns.get("qc_cell_counts_by_step")
    if not isinstance(history, Mapping):
        return empty_qc_filter_counts()

    rows = pd.DataFrame(history)
    required = {"sample", "step", "n_cells"}
    if rows.empty or not required.issubset(rows):
        return empty_qc_filter_counts()

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
    ignored = {"component", "celltype_manual"}
    result = []
    for column in metadata.columns:
        if column in ignored or column.rsplit("_", 1)[-1].isdigit():
            continue
        if metadata.groupby("sample_alias", observed=True)[column].nunique(dropna=False).max() <= 1:
            result.append(column)
    return result


def apply_analysis_grouping(data: AppData, mapping: Mapping[str, str], label: str) -> AppData:
    clean = {str(sample): str(group).strip() for sample, group in mapping.items()}
    if any(not group for group in clean.values()):
        raise ValueError("Analysis group labels cannot be blank.")
    samples = set(data.metadata["sample_alias"].astype(str))
    if samples.difference(clean):
        raise ValueError("Every sample must have a non-empty analysis group.")

    metadata = data.metadata.copy()
    metadata["condition"] = metadata["sample_alias"].astype(str).map(clean)

    def regroup(frame: pd.DataFrame) -> pd.DataFrame:
        frame = frame.copy()
        if "sample_alias" in frame:
            groups = frame["sample_alias"].astype(str).map(clean)
        elif "sample" in frame:
            groups = frame["sample"].astype(str).map(clean)
        elif "component" in frame:
            conditions = metadata.set_index("component")["condition"]
            groups = frame["component"].astype(str).map(conditions)
        else:
            return frame
        if "condition" in frame:
            frame["condition"] = groups.where(groups.notna(), frame["condition"])
        else:
            frame["condition"] = groups
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
    if column == "sample_alias":
        samples = metadata["sample_alias"].astype(str).drop_duplicates()
        return dict(zip(samples, samples))
    rows = metadata[["sample_alias", column]].drop_duplicates("sample_alias")
    values = rows[column].astype("string").fillna("Unassigned").str.strip().replace("", "Unassigned")
    return dict(zip(rows["sample_alias"].astype(str), values.astype(str)))


def new_analysis_grouping_config(metadata: pd.DataFrame) -> dict:
    columns = sample_level_columns(metadata)
    if not columns:
        raise ValueError("No sample-level metadata columns are available for analysis grouping.")
    source_columns = list(dict.fromkeys(["sample_alias", *columns]))
    source = metadata[source_columns].drop_duplicates("sample_alias").copy()
    source["sample_alias"] = source["sample_alias"].astype(str)
    source = source.sort_values("sample_alias").reset_index(drop=True)
    column = "condition" if "condition" in columns else "sample_alias"
    return {
        "mode": "column",
        "column": column,
        "label": column,
        "columns": columns,
        "source": source,
        "mapping": mapping_for_column(source, column),
    }


def update_analysis_grouping_config(
    config: Mapping,
    *,
    mode: str,
    column: str,
    custom_groups: Mapping[str, str] | None = None,
) -> dict:
    if mode not in {"column", "custom"}:
        raise ValueError("Grouping mode must be 'column' or 'custom'.")
    source = config["source"]
    mapping = mapping_for_column(source, column)
    if mode == "custom":
        custom_groups = {str(sample): str(group).strip() for sample, group in (custom_groups or {}).items()}
        missing = set(mapping).difference(custom_groups)
        if missing:
            raise ValueError("Enter one analysis group for every sample.")
        mapping = {sample: custom_groups[sample] for sample in mapping}
        if any(not group for group in mapping.values()):
            raise ValueError("Analysis group labels cannot be blank.")
    return {
        **config,
        "mode": mode,
        "column": column,
        "label": "Custom sample groups" if mode == "custom" else column,
        "mapping": mapping,
    }


def analysis_grouping_summary(config: Mapping | None) -> str:
    if not config:
        return "Analysis grouping unavailable"
    return f"Analysis grouping: {config['label']} · {len(set(config['mapping'].values()))} groups"


def filter_pixelator_proximity(
    proximity: pd.DataFrame,
    *,
    min_marker_fraction: float = 0,
    min_marker_count: float = 0,
    min_cells: int = 1,
) -> pd.DataFrame:
    rows = proximity.copy()
    if min_marker_fraction > 0:
        required = {"marker_1_freq", "marker_2_freq"}
        if not required.issubset(rows):
            raise ValueError("Marker-fraction filtering requires marker_1_freq and marker_2_freq.")
        rows = rows[
            (pd.to_numeric(rows["marker_1_freq"], errors="coerce") >= min_marker_fraction)
            & (pd.to_numeric(rows["marker_2_freq"], errors="coerce") >= min_marker_fraction)
        ]
    if min_marker_count > 0:
        if "min_count" not in rows:
            raise ValueError("Marker-count filtering requires min_count.")
        rows = rows[pd.to_numeric(rows["min_count"], errors="coerce") >= min_marker_count]
    if min_cells > 1 and not rows.empty:
        counts = rows.groupby(["marker_1", "marker_2"], observed=True)["component"].transform("nunique")
        rows = rows[counts >= int(min_cells)]
    return rows


def select_colocalization_heatmap_markers(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    available_markers: Sequence[str],
    *,
    n_markers: int = 40,
    plot_markers: Sequence[str] | None = None,
) -> list[str]:
    """Select PixelatorES heatmap markers by sample-weighted mean abundance."""
    available = list(dict.fromkeys(map(str, available_markers)))
    available_set = set(available)
    if plot_markers is not None:
        return [
            marker
            for marker in dict.fromkeys(map(str, plot_markers))
            if marker in available_set
        ]
    if abundance.empty or metadata.empty or not available:
        return []

    required_abundance = {"component", "marker", "abundance"}
    required_metadata = {"component", "sample_alias"}
    if missing := required_abundance.difference(abundance.columns):
        raise ValueError("Abundance table is missing: " + ", ".join(sorted(missing)))
    if missing := required_metadata.difference(metadata.columns):
        raise ValueError("Metadata table is missing: " + ", ".join(sorted(missing)))

    sample_by_component = (
        metadata[["component", "sample_alias"]]
        .drop_duplicates("component")
        .assign(component=lambda frame: frame["component"].astype(str))
        .set_index("component")["sample_alias"]
    )
    rows = abundance.loc[
        abundance["component"].astype(str).isin(sample_by_component.index)
        & abundance["marker"].astype(str).isin(available_set),
        ["component", "marker", "abundance"],
    ].copy()
    rows["component"] = rows["component"].astype(str)
    rows["marker"] = rows["marker"].astype(str)
    rows["sample_alias"] = rows["component"].map(sample_by_component)
    rows["abundance"] = pd.to_numeric(rows["abundance"], errors="coerce")
    rows = rows.dropna(subset=["sample_alias", "abundance"])
    if rows.empty:
        return []

    sample_means = (
        rows.groupby(["sample_alias", "marker"], observed=True, dropna=False)["abundance"]
        .mean()
        .rename("sample_mean_abundance")
        .reset_index()
    )
    marker_means = (
        sample_means.groupby("marker", observed=True)["sample_mean_abundance"]
        .mean()
        .rename("mean_abundance")
        .reset_index()
        .sort_values(["mean_abundance", "marker"], ascending=[False, True], kind="stable")
    )
    count = min(max(1, int(n_markers)), len(marker_means))
    return marker_means["marker"].head(count).tolist()


def select_clustering_heatmap_markers(
    summary: pd.DataFrame,
    available_markers: Sequence[str],
    *,
    mode: Literal["top", "manual"] = "top",
    n_markers: int = 20,
    selected_markers: Sequence[str] | None = None,
) -> list[str]:
    """Select variable or explicit proteins for the clustering heatmap."""
    available = set(map(str, available_markers))
    if mode == "manual":
        return [
            marker
            for marker in dict.fromkeys(map(str, selected_markers or ()))
            if marker in available
        ]
    if mode != "top":
        raise ValueError("Clustering heatmap marker mode must be 'top' or 'manual'.")
    if summary.empty:
        return []
    ranking = (
        summary[summary["marker"].astype(str).isin(available)]
        .groupby("marker", observed=True)["mean_log2_ratio"]
        .std()
        .fillna(0)
        .sort_values(ascending=False, kind="stable")
    )
    return ranking.head(max(1, int(n_markers))).index.astype(str).tolist()



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
    if detected.empty or totals.empty:
        return pd.DataFrame()
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
    summary = summary.copy()
    summary[group_col] = summary[group_col].astype(str)
    totals = totals.copy()
    totals[group_col] = totals[group_col].astype(str)
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
    for column in ("sum_log2_ratio", "detected_mean", "mean_log2_ratio", "n_detected", "pct_detected"):
        if column not in result:
            result[column] = 0.0
        result[column] = result[column].fillna(0)
    return result


def sample_colocalization(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    celltypes: Sequence[str] | None = None,
    mean_type: str = "population",
) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame()
    summary = summarize_spatial(
        scores,
        metadata,
        group_col="sample_alias",
        markers=sorted(set(scores["marker_1"]).union(scores["marker_2"])),
        celltypes=celltypes,
        mean_type=mean_type,
    )
    if summary.empty:
        return summary
    sample_groups = metadata[["sample_alias", "condition"]].drop_duplicates("sample_alias")
    summary = summary.merge(sample_groups, on="sample_alias", how="left")
    summary["marker_pair"] = summary["marker_1"].astype(str) + " / " + summary["marker_2"].astype(str)
    summary["pair_observed"] = summary["n_detected"] > 0
    return summary


def read_component_layout(data: AppData, sample: str, component: str) -> pd.DataFrame:
    if component in data.component_layouts:
        layout = data.component_layouts[component].copy()
    else:
        if not data.pxl_files:
            raise ValueError("Assign a .layout.pxl path in the Data menu to display this cellgraph.")
        try:
            from pixelator import read_pna
        except ImportError as error:
            raise RuntimeError("Cellgraph display requires the pixelgen-pixelator package.") from error

        metadata = data.metadata[data.metadata["component"].astype(str) == str(component)]
        sample_ids = {str(sample)}
        if not metadata.empty:
            sample_ids.update(str(metadata.iloc[0][column]) for column in ("sample", "sample_alias"))
        candidates = [
            Path(path) for path in data.pxl_files
            if any(sample_id and sample_id in Path(path).stem for sample_id in sample_ids)
        ]
        if not candidates and len(data.pxl_files) == 1:
            candidates = [Path(data.pxl_files[0])]
        if not candidates:
            raise FileNotFoundError(f"No assigned PXL filename matches sample {sample}.")

        dataset = read_pna(candidates[0])
        available = set(map(str, dataset.components()))
        raw_component = str(component)
        for sample_id in sorted(sample_ids, key=len, reverse=True):
            prefix = f"{sample_id}_"
            if raw_component not in available and raw_component.startswith(prefix):
                raw_component = raw_component[len(prefix):]
        if raw_component not in available:
            raise ValueError(f"Component {component} is not present in {candidates[0].name}.")
        layout = dataset.filter(components=[raw_component]).precomputed_layouts(add_marker_counts=False).to_df()
        if layout.empty:
            raise ValueError("No precomputed 3D layout is stored for this component.")
    if not {"x", "y", "z"}.issubset(layout):
        raise ValueError("Stored component layout must contain x, y, and z columns.")
    return layout



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
