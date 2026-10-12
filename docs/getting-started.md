# Getting started

Run commands from the `ProxiomeVis` directory.

## 1. Prepare the demo RDS

The repository expects the demo source at `../pixelgen_demo_data`. Build the
runtime file from its CCA object and tertiary results:

```bash
Rscript --vanilla scripts/prepare_precomputed_rds.R \
  ../pixelgen_demo_data/demo_analysis/cca_analysis/results/demo_cca_clustered.rds \
  ../pixelgen_demo_data/demo_analysis/tertiary_analysis \
  output/proxiomevis_demo_v1.rds
```

The source tree is read-only. The script creates or replaces only the explicit
output path after validating the source object and frozen tables.

## 2. Prepare the hybrid application bundle

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 \
  Rscript --vanilla scripts/prepare_hybrid_bundle.R \
  output/proxiomevis_demo_v1.rds output/hybrid_demo
```

This exports a plain-R `analysis.rds`, observed proximity rows in `proximity/`,
and node-only 3D tables in `graphs/`: one Parquet file per prepared sample and a
version 2 `index.rds`. Each row contains a cell ID, node ID, 3D coordinates, and
marker identities; no edges or marker counts are exported. Rerunning reuses
validated completed sample files.
Supply the source PXL directory as the optional third argument if paths moved.
Keep the bundle and sidecar directories together when copying to a server.

To prepare 3D assets for only one or a few cells, supply a fourth argument: a
text file containing one exact cell ID per line, with no header, blank lines,
or duplicates. The third argument must be the PXL directory:

```bash
Rscript --vanilla scripts/prepare_hybrid_bundle.R \
  output/proxiomevis_demo_v1.rds output/hybrid_selected \
  ../pixelgen_demo_data/data/nf_pixelator_outdir/pixelator graph_cells.txt
```

Only graph preparation is restricted by this list; the analysis bundle and
proximity sidecar still cover the dataset. **Cell 3D** shows prepared-cell
coverage and offers only available assets. Rerunning with additional cells
atomically rewrites each affected sample file with the union of old and new
cells and updates the graph index. Previously prepared cells are retained. Run
preparers sequentially, not concurrently against the same directory. Omit the
fourth argument to prepare all cells, or call
`build_hybrid_bundle(..., graphs = FALSE)` from R to omit graph preparation.

Old version 1 per-cell graph RDS exports are not compatible with this reader.
Preserve that directory and rebuild into a new output directory, then deploy the
new bundle and sidecars together. No automatic destructive migration is performed.

## 3. Start Shiny

The default input is `output/hybrid_demo/analysis.rds` when it exists, falling
back to `output/proxiomevis_demo_v1.rds`:

```bash
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

Open `http://127.0.0.1:3838`. To use another trusted prepared file:

```bash
export PROXIOME_PRECOMPUTED_RDS=/absolute/path/to/results.rds
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

Plain bundles require no Seurat/PNA class packages for viewing. Legacy inputs
still require those packages; new statistical tests require the analysis-worker
dependencies listed in the README. **Cell 3D** requires `arrow` and `dplyr` in
the Shiny R environment. No runtime PXL mount is required.

The app loads the configured RDS once and does not expose a visitor-selectable
server path. For trusted local administration only, set
`PROXIOME_ALLOW_LOCAL_LOAD=true` before startup to show the path field. Never
enable that setting on a public deployment.

## 4. Explore or run a new comparison

Donor, comparison, cell-type, and marker controls filter stored rows. They do
not trigger a new statistical test or rebuild an aggregate in the report tabs.
The combined two-donor tables are descriptive and have no p-values.
Donor-specific volcanoes reproduce stored exploratory component-level tests;
they do not provide population-level donor inference.

Use **Custom analysis** for a new selected-cell comparison: choose one donor,
one cell type, two conditions and markers; optionally lasso cells in the embedding;
acknowledge the exploratory scope and click **Run comparison**. You can navigate
while the separate R process works or cancel it. Export results and their JSON
provenance. Report tabs remain unchanged. Missing spatial sidecars disable only
spatial reanalysis, not the report viewer.

Use **Cell 3D** to rank cells within donor, cell type, and condition by two markers'
observed log2 ratio (self-pairs allowed; pair order does not matter).
Choose a ranked cell and one highlight marker, then click **Display** to load it.
Changing controls does not update the displayed graph until Display is clicked. Rotate,
zoom, hover, and export its nodes and provenance. The view is nodes only, and all
nodes remain visible unless the cell exceeds the 200,000-node limit. Coverage is
shown explicitly. Initial access may checksum the whole sample file; later cell
queries use Arrow row-group pruning and cached checksum validation. Coordinates
represent a graph layout, not physical microscopy positions, and marker
highlights indicate identities rather than counts or a colocalization test.

Read [Hybrid analysis](hybrid-analysis.md) before interpreting new P values.

## Tests

```bash
Rscript --vanilla tests/testthat.R
```

Synthetic fixtures include small PXL-like databases and background-worker tests;
CI does not need the original PXL files or the demo dataset.
