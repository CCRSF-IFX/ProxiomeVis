# Patch Analysis

Patch Analysis follows Pixelgen's receiver/target workflow without assuming a
particular experiment or cell type. A patch represents a connected subgraph on
a receiver cell enriched for markers from exactly one target population.

## Prepare Patch Data

Population, analysis-group, and marker controls are a draft until you click
**Prepare Patch Data**. Preparation freezes the selected scope, runs the
filtered candidate PXL query, and records the preparation time. The active
scope summary reports receiver cells, target cells, analysis groups, and target
markers. If a retrieval control changes, existing candidate and stored-result
views remain pinned to the active scope and show an unapplied-changes warning.
Candidate count and score thresholds are display controls and do not require
preparing the data again.

## Marker Selection

Choose a categorical H5AD metadata field, one receiver population, and one
target population. The app calculates each marker's raw fraction of total
counts in both populations and suggests markers using configurable minimum
fraction and fold-enrichment thresholds. Suggested target and receiver markers
remain editable. The selected analysis groups scope marker selection, candidate
screening, and component-level stored results consistently.

This raw-frequency screen is not the abundance-unmixing method implemented by
PixelatorR. The UI reports the combined target-marker fraction because the
Pixelgen tutorial recommends that selected markers represent at least 20–30%
of the target population's signal.

## Candidate Signal

For receiver cells in the selected analysis groups, the app queries PXL for the
joint target-marker proximity score and combines it with target-marker UMI
counts from H5AD. Adjustable count and log2-ratio reference lines identify
candidate cells for review. They do not perform connected-subgraph patch
detection.

## Detected Patches

Patch detection is an experimental, graph-intensive PixelatorR workflow. The
Python app does not silently approximate it from proximity scores. Run patch
detection upstream in bounded batches and store its results under
`uns["proxiome"]["patch"]`.

The preferred generic tables are:

- `run_plan`: detection status, receiver/target definitions, marker counts, and
  algorithm settings.
- `patch_sizes`: one row per component and patch, including `patch`, `count`,
  and optionally graph fraction `p` and patch size `type`.
- `patch_burden`: one row per receiver component with numeric patch count or
  coverage metrics.
- `patch_composition`: long-form `component`, `patch`, `marker`, and `count`
  rows, optionally with a `source` column. A wide marker-count table is also
  accepted.
- `marker_unmixing`: optional upstream marker-selection provenance.
- `target_marker_abundance` and `target_marker_proximity`: optional upstream
  candidate-screen tables used when live PXL querying is unavailable.

Legacy `raji_marker_abundance` and `raji_marker_proximity` tables are mapped to
their generic target-marker names when loaded.

Include `component` identifiers in candidate, burden, and composition tables so
the app can join them to H5AD metadata and apply the active receiver and analysis
group scope. When `run_plan` includes `receiver_population` and
`target_population`, those values initialize the population selectors.

Design reference: [Pixelgen Patch analysis tutorial](https://software.pixelgen.com/pna-analysis/R/tutorials/patch_analysis/).
