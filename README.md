# ProxiomeVis

ProxiomeVis is an R Shiny app for exploring Pixelgen proxiome data. It combines
frozen report views with bounded, background differential analysis on selected
cells and an on-demand 3D cell graph viewer. Expensive proximity scoring and
graph preparation happen upstream; no runtime view or test opens PXL files.

The default viewer loads a compact, versioned plain-R bundle. Optional sidecar
files supply observed cell-level scores and per-sample 3D node tables without
loading those whole datasets into each Shiny session. New exploratory results
never replace the frozen report results.

## Requirements

- Viewer: R with `Matrix`, `shiny`, `bslib`, `ggplot2`,
  `plotly`, `htmltools`, `DT`, `jsonlite`, `data.table`, `ggrepel`, `scales`, `tibble`,
  `pheatmap`, `yaml`, `callr`, `promises`, `later`, and `digest`
- Analysis worker: additionally `SeuratObject`, `Seurat`, `pixelatorR`, `arrow`,
  and `dplyr`; preparation also needs `DBI` and `duckdb`
- Cell 3D: additionally `arrow` and `dplyr` to query prepared node tables
- A trusted application bundle; legacy Seurat RDS inputs remain supported

Plain application bundles do not require Seurat or pixelatorR to open.
Legacy PNA/Seurat RDS inputs still need their class-defining packages.

## Prepare the RDS

Build the runtime file from the annotated Seurat RDS and the matching analysis
output directory:

```bash
Rscript --vanilla scripts/prepare_precomputed_rds.R \
  ../pixelgen_demo_data/demo_analysis/cca_analysis/results/demo_cca_clustered.rds \
  ../pixelgen_demo_data/demo_analysis/tertiary_analysis \
  output/proxiomevis_demo_v1.rds
```

Then prepare the compact application bundle and its sidecars:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 \
  Rscript --vanilla scripts/prepare_hybrid_bundle.R \
  output/proxiomevis_demo_v1.rds output/hybrid_demo
```

The result is `output/hybrid_demo/analysis.rds`, `proximity/` (partitioned
Parquet), and `graphs/` (one node-only Parquet file per prepared sample plus a
version 2 `index.rds`). Keep these together
when deploying. The optional third argument overrides the source PXL directory;
a fourth argument lists selected graph cells, one exact ID per line. This limits
only 3D preparation, not the reports or proximity data. See
[Getting started](docs/getting-started.md#2-prepare-the-hybrid-application-bundle).
Rerunning validates completed sample files and preserves their prepared cells.
Adding cells rewrites the affected sample file atomically with the union of
previously prepared and newly requested cells. Do not run two preparers against
the same output directory concurrently. Version 1 per-cell graph exports must
be rebuilt in a new directory; they are not automatically migrated or deleted.
The legacy precomputed RDS alone still supports all existing report views.

`output/` is intentionally gitignored. Build the artifact on the deployment
host or copy it there through the project's controlled data-delivery process;
do not commit the demo RDS to this source repository.

The builder copies the frozen `abundance`, `self_proximity`, and
`colocalization` readouts into the versioned `object@misc$proxiomevis`
payload. See [Precomputed results](docs/precomputed-results.md) for the input
files, stored fields, and scientific limits.

The current builder targets this demo analysis layout and its `d1`/`d2`
outputs. The versioned runtime contract and Shiny reader remain separate from
that demo-specific preparation code.

`../pixelgen_demo_data` is a read-only source dataset. The builder reads the CCA
object and tertiary results there and writes only the requested `OUTPUT_RDS`.

## Run the app

By default, the app reads `output/hybrid_demo/analysis.rds` when present,
otherwise `output/proxiomevis_demo_v1.rds`. To use a
different trusted server-local file, set `PROXIOME_PRECOMPUTED_RDS`:

```bash
export PROXIOME_PRECOMPUTED_RDS=/absolute/path/to/hybrid_demo/analysis.rds
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

The repository workspace also provides an isolated R executable. Limit BLAS
threads when using it on a shared host:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 \
  ../.pixi/envs/r/bin/Rscript --vanilla \
  -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

Then open `http://127.0.0.1:3838` from the same machine, or use the secure
proxy provided by the deployment environment.

In the report tabs, donor and cell-type controls filter results already computed.
They do not recompute a multi-donor aggregate for an arbitrary subset. This
dataset contains two donors, so the combined donor-matched results are
descriptive and have no p-values. Donor-specific report volcanoes retain the
tertiary analysis's stored component-level p-values and adjusted p-values; they
are labeled exploratory and are not population-level donor inference. The app
never treats cells as independent biological replicates.

The **Report volcano** view supports marker hover, scroll-wheel zoom, and
drag-to-box zoom for all three readouts. Double-click to reset, or use the
Plotly toolbar to pan/reset. Hover replaces repelled labels on screen; the
Volcano SVG download retains the original labeled report plot. Spatial report
and abundance heatmaps remain static.

**Custom analysis:** choose one donor, one cell type, two conditions, and markers.
Optionally select cells using the embedding lasso/box or searchable cell selector.
Confirm the exploratory-analysis notice, then run. The separate R worker allows
continued navigation and can be cancelled. Results and JSON provenance record
the submitted cell groups, feature family, bundle checksum, methods, and software.
Changing controls does not change the completed result until another run.

**Cell 3D:** choose donor, cell type, and condition, then two ranking markers
(including self-pairs). Ranking uses observed log2 ratios, regardless of pair order,
with top 10, bottom 10, their union, or all scored cells available.
Choose a cell and one highlight marker, then click **Display**.
Changing controls leaves the last displayed view unchanged until Display is clicked. All nodes
remain visible. The node-only sidecar stores coordinates and marker identities,
not edges or marker counts. Rotate, zoom, hover, and export the selected cell's
nodes and provenance. These are abstract graph-layout coordinates, not measured
physical distances; shared-node highlights are not a statistical colocalization
test. The viewer refuses cells above 200,000 nodes without silently downsampling.
Prepared coverage is shown explicitly. The first access to a sample may checksum
its whole file; later cell selections reuse validation and query matching
Parquet row groups rather than load the sample into R.

See [Hybrid analysis](docs/hybrid-analysis.md) for methods, data schemas, and limits.

## Data handling

Use RDS paths that are trusted and already visible to the R process. `readRDS()`
can deserialize executable R objects, so do not expose unrestricted public RDS
uploads. For hosted use, place prepared files in an administrator-controlled
directory and pass the selected path through `PROXIOME_PRECOMPUTED_RDS`.

The in-app server-path field is disabled by default. A trusted local operator
can enable it before startup with `PROXIOME_ALLOW_LOCAL_LOAD=true`. Never enable
that flag on a public deployment: it lets a visitor request any RDS path that
the Shiny process can read.

The deployed bundle and sidecars need no PXL mount. Keep them read-only and
replace datasets atomically; requests are pinned to the bundle loaded by the UI.
Each Shiny R process permits one analysis worker and rejects overlapping work
with a retry message. Jobs are limited to 120 seconds, 2,000 cells, 160 markers,
1,000 pairs, and 1 million collected score rows. These are work limits, not an
OS memory guarantee; set deployment-level CPU/RAM quotas too. Multiple Shiny R
processes each have their own worker limit. Results are session-local; there is
no shared result cache or durable multi-server job queue.

## Tests

```bash
Rscript --vanilla tests/testthat.R
```

Tests use synthetic Seurat inputs, small locally generated PXL-like databases,
and actual isolated R subprocesses; no source demo download is required. Install
the worker/preparation dependencies to exercise the full suite.

With the prepared demo app and Chromium CDP on port 9223:

```bash
node tests/browser_smoke.mjs
node tests/browser_hybrid_smoke.mjs
```

The 3D browser check needs WebGL (software rendering is sufficient).

## Python implementation

The Python implementation remains preserved and maintained on the
`py_version` branch. The R implementation is developed independently; its
history does not reset or overwrite the Python branch.
