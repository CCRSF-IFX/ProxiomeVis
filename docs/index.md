# ProxiomeVis R

ProxiomeVis is an R Shiny explorer for processed Pixelgen data. It combines
prepared cell annotations and embeddings with frozen abundance,
self-proximity, and colocalization results from the tertiary analysis.

The app provides:

- QC cell counts and stored cell annotations;
- interactive embeddings colored by annotation or marker abundance;
- per-sample summaries, donor-specific contrasts, and stored cross-donor
  descriptive effects for the three readouts;
- downloadable tables for the currently displayed stored results;
- bounded background comparisons on selected cells within one donor/cell type;
- on-demand, node-only 3D cell views with marker-identity highlights.

The app does not open PXL files or recalculate proximity scores. A preparation
workflow builds a compact plain-R bundle and optional disk-backed proximity and
graph sidecars. Sessions share the startup bundle read-only within each R process;
analysis workers load only a bounded selection of observed scores. The graph
sidecar uses one Parquet file per sample and a version 3 index. Arrow selects one
cell's nodes for viewing; no edges or marker-count matrices are stored there.

## Interpretation boundary

The demo contains two donors. Combined donor-matched effects are descriptive,
contain no p-values, and are not population-level inference. Donor-specific
report volcanoes reproduce stored exploratory component-level p-values from the
tertiary analysis; Shiny does not recalculate them. Cells are not independent
biological replicates.

Start with [Getting started](getting-started.md), then review the
[precomputed-results contract](precomputed-results.md) and
[hybrid analysis guide](hybrid-analysis.md) before preparing a new dataset.
