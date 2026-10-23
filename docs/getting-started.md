# Getting started

Run commands from the `ProxiomeVis` directory.

The Cell 3D renderer uses the CRAN `threejs` package (>= 0.3.4):
`install.packages("threejs")`. Plotly remains in use for the other interactive plots.

For your own annotated Seurat object, follow
[Prepare your data](prepare-user-data.md). The general builder does not require
the demo report folders. The following preparation steps are for the demo.
If you already have a trusted application bundle, go to **3. Start Shiny**.

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
version 3 `index.rds`. Each row contains a cell ID, 3D coordinates, and
marker identities; no edges or marker counts are exported. Rerunning reuses
validated completed sample files.
Existing version 2 graph sidecars remain readable, but upgrading to the smaller
version 3 files requires preparation in a new output directory.
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
coverage; ranked cells without an asset are labeled and cannot be displayed.
Rerunning with additional cells
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

Install the packages listed in `DESCRIPTION`'s `Imports`, including Seurat.
A plain bundle does not contain Seurat or PNA objects. Legacy inputs and new
tests need the applicable analysis packages. **Cell 3D** requires `arrow` and
`dplyr` in the Shiny environment. The running app does not need the source PXL
files.

The app loads the configured RDS once and does not expose a visitor-selectable
server path. For trusted local administration only, set
`PROXIOME_ALLOW_LOCAL_LOAD=true` before startup to show the path field inside
the upper-right **Dataset** panel. Never
enable that setting on a public deployment.

## 4. Explore or run a new comparison

Use the top navigation for **Overview / QC**, **Embeddings**, **Results**,
**New comparison**, **Cell 3D**, and **Documentation**. Open **Dataset** at the
upper right to check the active bundle and configured sidecars.

### View results

Open **Results**. Select **Prepared results** or a completed session comparison
under **Result source**. Both use the same views. For prepared results, choose
a donor, cell type, comparison, and readout. A dataset with one donor shows a
fixed donor label. A session comparison retains its submitted scope.

- **Tables:** inspect searchable results, sample support, method settings, and CSV downloads.
- **Volcano:** inspect the test beside a combined violin and boxplot. Click a volcano point to change the displayed cell distribution. Hover for details. Scroll or drag a box to zoom. Double-click to reset.
- **Spatial patterns:** inspect self-proximity and colocalization profiles or differential effects. For colocalization, select the network view to compare condition profiles with shared node positions.

The volcano and spatial views require one cell type. **Spatial patterns** is
disabled for abundance. Select self-proximity or colocalization to enable it.
Profiles need stored summaries but do not require a differential test.

Click a plot feature or result-table row to select that feature across views.
The pairs A–B and B–A identify the same feature. **Explore cells in 3D** transfers
the donor, cell type, and pair. It does not load a graph.

FDR and search controls affect tables and **Filtered table CSV** downloads.
**All stored rows CSV** ignores these filters within the chosen scope. A table
filter does not change the volcano or rerun a test. Stored summaries and test
results do not need sidecar reads. Cell-level spatial distributions require
the proximity sidecar; the app shows a message if it is unavailable.

### Align spatial matrices

In **Spatial patterns**, open **Shared proteins and order**. Select proteins
and choose the prepared order, alphabetical order, or a custom ordered list.
For a custom list, separate names with commas or new lines. For example:

```text
CD19, CD37, CD53, CD81, CD82
```

Click **Apply**. This list controls both axes of every spatial matrix.
Click **Reset** to restore the dataset defaults. The applied list stays the
same when you switch profiles, comparisons, or result sources. Loading another
dataset resets it. Selection and ordering do not rerun tests or change adjusted
P values.

Choose **Single** for one matrix. Choose **Compare profiles and difference**
to show the reference profile, target profile, and differential effect together.
All three panels use the same protein positions. The two profiles share a fixed
color scale. The differential panel uses a separate fixed scale. Gray cells
mean that a result is unavailable or untested; missing values are not zero.

Profile means are not the group medians used for differential effects.
Positive effects mean target minus reference. Diagonal entries show
self-proximity; other entries show colocalization. Use the SVG and CSV buttons
to export the applied view. **Protein order CSV** exports the ordered list.

The colocalization network has its own marker and display-threshold controls.
Its stored population means include missing scores as zero. Its display
thresholds are not significance tests.

### Run and keep a comparison

Open **New comparison**. Choose one donor, one cell type, two different
conditions, and a readout. The test uses all cells in this scope and all dataset
markers. The app has no marker-test or cell-subset selector.

Read the readiness message before you run. It shows cell and marker counts,
scope errors, and the required input status. Abundance requires stored counts.
Spatial tests require a configured proximity sidecar. The worker checks
sidecar integrity and score coverage after submission. Adjust **Minimum
observed cells per condition** for a spatial test if necessary.

Confirm the exploratory-analysis notice and click **Run comparison**. You can
use other views while the separate R process runs. Click **Cancel running
analysis** to stop it. Comparisons are limited to 2,000 cells, 160 markers,
and 120 seconds. Spatial processing uses batches of at most 1,000 pairs and
1,000,000 score rows. There is no 1,000-pair limit on the complete comparison.
For a larger scope, prepare results offline.

Each completed test appears in session history. Select a completed comparison
and click **Open in Results**. Prepared results remain unchanged. Use the CSV,
JSON, or **Comparison RDS** download to keep the selected result. The RDS stores
results and settings, not cell-level measurements. See
[Prepare your data](prepare-user-data.md) for offline reuse.

Loading another dataset cancels the active test and clears history. Closing
the session also removes its history. Download results that you need to keep.
Within-donor cell-level P values do not establish a population-level effect.

### Display a cell

Use **Cell 3D** to rank cells within donor, cell type, and condition by two markers'
observed log2 ratio (self-pairs allowed; pair order does not matter).
After a feature handoff, choose a condition from the comparison and one highlight
marker. Choose a ranked cell, then click **Display** to load it.
Changing controls does not update the displayed graph until Display is clicked. Rotate,
zoom, hover, and export its nodes and provenance. The view is nodes only, and all
nodes remain visible unless the cell exceeds the 200,000-node limit. Coverage is
shown explicitly. Initial access may checksum the whole sample file; later cell
queries use Arrow row-group pruning and cached checksum validation. Coordinates
represent a graph layout, not physical microscopy positions, and marker
highlights indicate identities rather than counts or a colocalization test.

Read [Results and new comparisons](hybrid-analysis.md) before you interpret new P values.

## Tests

```bash
Rscript --vanilla tests/testthat.R
```

Synthetic fixtures include small PXL-like databases and background-worker tests;
CI does not need the original PXL files or the demo dataset.
