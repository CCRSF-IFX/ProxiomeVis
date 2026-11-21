# ProxiomeVis

ProxiomeVis is a Shiny app for processed Pixelgen data. It combines prepared
cell annotations and embeddings with protein abundance, self-proximity, and
colocalization results. Use your own prepared data or the supplied demo.

The app provides:

- QC cell counts and stored cell annotations;
- interactive embeddings colored by annotation or marker abundance;
- one Results viewer for prepared results and completed session comparisons;
- sample profiles, within-donor exploration, and across-donor comparisons;
- downloadable tables for the currently displayed stored results;
- background comparisons of conditions or cell types, with explicit donor groups;
- a shared protein list and order for spatial matrices;
- on-demand, node-only 3D cell views with marker-identity highlights.

The app does not open PXL files or recalculate proximity scores. A preparation
workflow builds a compact plain-R bundle and optional disk-backed proximity and
graph sidecars. Sessions share the startup bundle read-only within each R process;
analysis workers load only a bounded selection of observed scores. The graph
sidecar uses one Parquet file per sample and a version 3 index. Arrow selects one
cell's nodes for viewing; no edges or marker-count matrices are stored there.

## Select the analysis design

Keep **Within-donor exploration** when biological replicates are not available.
It compares conditions in one cell type, or cell types in one condition, within
one donor. Both prepared and new results support this mode. Its cell-level
P values are exploratory, not evidence of reproducibility across donors.

Use **Across donors** when donors are the biological replicates. Paired analysis
compares conditions or cell types within matching donors. Independent analysis
compares conditions with different donors in each group. Samples from the same
donor do not become additional replicates. Donors have equal weight. Missing
measurements are not zeros. Two donors per group or two pairs are only a
computational test minimum, not evidence of adequate statistical power.

Prepare a comparison offline with `scripts/run_comparison.R`, or calculate it in
**New comparison**. Both use the same analysis code. Saved comparisons can be
appended to a new application bundle without removing earlier results. Analysis
schema version 2 stores donor measurements and comparison design. Schema
version 1 datasets remain readable with their original statistics.

## Demo interpretation

The demo contains two donors. Combined donor-matched effects are descriptive,
contain no p-values, and are not population-level inference. Donor-specific
report volcanoes reproduce stored exploratory component-level p-values from the
tertiary analysis; Shiny does not recalculate them. Cells are not independent
biological replicates.

Start with [Getting started](getting-started.md), then review the
[precomputed-results contract](precomputed-results.md) and
[hybrid analysis guide](hybrid-analysis.md) before preparing a new dataset.
