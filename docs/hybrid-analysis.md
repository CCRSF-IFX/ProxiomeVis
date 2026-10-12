# Prepared reports and new exploratory analyses

The application separates three operations: viewing frozen results, submitting
a bounded new comparison to an R worker, and opening prepared cell graphs.
Viewing or filtering a report never changes its stored statistics. Runtime
analysis reads a prepared proximity sidecar, not PXL files.

## Running a new comparison

Open the exploratory-analysis tab, choose one donor, one cell type, two distinct
conditions, a readout, and the markers to test. Use all cells in that scope, or
select a subset using the embedding's lasso/box tools or the cell selector. The
selected conditions assign cells to target and reference groups; cells cannot
belong to both. The displayed counts update before submission.

Confirm the exploratory-analysis notice and click **Run comparison**. An R
process performs the test while the reports remain usable. Changing controls
does not rerun a test. The result identifies the submitted scope; downloads
include the results CSV and JSON provenance with the exact selected cell IDs.
The volcano supports hover, box/scroll zoom, and double-click reset. The results
table also includes requested features that could not be tested.

Only one background analysis runs at a time per application process. The server
limits each request to 2,000 cells, 160 markers, 1,000 unordered non-self pairs,
and 1,000,000 retained score rows, with a 120-second timeout. Colocalization tests
all combinations of the chosen markers: 45 markers produce 990 pairs; 46 exceed
the cap. Oversized selections are rejected, not silently downsampled. Cancel a
running job or narrow the cells/markers if it exceeds a limit.

## Methods and interpretation

| Readout | Effect: target minus reference | Test and minimum support |
| --- | --- | --- |
| Abundance | Mean stored normalized abundance | `pixelatorR::RunDAA`, using Seurat's two-sided Wilcoxon test; at least 3 selected cells per group |
| Self-proximity | Median observed `log2_ratio` for each marker paired with itself | `pixelatorR::DifferentialProximityAnalysis(metric_type = "self")`; at least 20 retained observations per marker in each group |
| Colocalization | Median observed `log2_ratio` for each unordered pair of different markers | `pixelatorR::DifferentialProximityAnalysis(metric_type = "co")`; at least 40 retained observations per pair in each group |

Abundance uses the bundle's existing transformed layer without renormalization:
`mean_fxn = Matrix::rowMeans`, `fc_name = "difference"`, `slot = fc.slot = "data"`
on a temporary assay containing those stored values, `min.pct = 0`,
`test.use = "wilcox"`. `RunDAA` sets `logfc.threshold = 0` and applies BH correction
to the returned marker p-values. The demo layer contains CLR abundance. Three
cells is the API's testing minimum, not evidence that a small comparison is
well powered; the original pipeline recommends 20 for abundance.

Spatial tests use the package's two-sided, tie-corrected rank test with
continuity correction. The effect is a **difference of median log2 ratios**,
not the log fold-change in mean abundance. The preparation formula is
`log2(pmax(join_count, 1) / pmax(join_count_expected_mean, 1))`.

BH adjustment is performed separately over the eligible tested features in
each submitted request. Restricting markers changes this testing family, so
adjusted p-values can differ from the frozen report even when individual
effects and raw p-values agree. The app does not recalculate or replace the
original report's adjustment family.

These are exploratory comparisons of cells within **one donor**. Cells are not
independent biological replicates, and the results do not establish an effect
across donors. Selecting cells using the same markers or an embedding derived
from them can introduce selection bias. Neither the donor restriction nor a
small p-value resolves that problem.

## Missing scores and filtering

Preparation retains observed cell–pair rows. Before a spatial test, both
markers must pass the donor's pinned background-fraction threshold, and the
score must be finite. A missing, filtered-out, or non-finite observation is
never converted to zero. A measured score of zero remains an observation.

The demo's threshold is computed from **all prepared cells in that donor** as
`max(median(isotype_fraction) + 2 * mad(isotype_fraction), 0.001)`. It is fixed
during preparation, not recomputed after lasso selection. Both demo donors have
a cutoff of 0.001. The source pipeline's median/MAD settings are read from its
stored report configuration. The current preparer supports this isotype mode
without extra marker-count or pair-count filters; incompatible configurations
are rejected instead of silently changed.

Results distinguish `n_cells_a/b` (all selected cells) from `n_scored_a/b`
(retained observations for this particular feature). Spatial eligibility uses
the latter. Unsupported features have missing statistics and an explicit
status. Stored spatial report profiles may use a different missing-observation
summary convention; those means are never used to reconstruct new tests.

## Application bundle and sidecar schemas

The version 2 application RDS is a plain R list, not a Seurat object:

- `bundle_version = 2L`; `contract` retains the version 1 frozen-analysis
  contract, report settings, defaults, QC, and provenance.
- `metadata` and `embeddings` contain prepared cell annotations and coordinates.
- `abundance` and optional `counts` are aligned marker-by-cell sparse matrices.
  Counts are required for new abundance tests and proximity preparation.
- `cell_sources`, when present, maps `component`, `original_id`, and `sample_id`
  exactly from the original FSMap; IDs are never inferred by stripping prefixes.
- `assets$proximity = list(path = "proximity", manifest = "manifest.rds")`
  locates the proximity sidecar relative to the bundle directory.
- `assets$graphs = list(path = "graphs", index = "index.rds")` locates the
  optional node-only 3D sidecar.

### Proximity sidecar

The proximity directory contains a version 1 `manifest.rds` and one Parquet
partition per sample. Each row contains `component`, canonical `marker_1` and
`marker_2`, `log2_ratio`, and `min_marker_fraction` (the smaller of the two marker
counts divided by that cell's total count). Storing the minimum is equivalent
to requiring both marker fractions to pass the same threshold.

The manifest contains the dataset ID, a fingerprint of metadata, counts,
cell-ID mapping and filter settings, donor cutoffs, row counts, partition MD5s,
and source PXL paths/checksums. Preparation verifies source checksums against
the frozen report inventory and checks that inputs did not change during
extraction. PXL connections are read-only; SQL exports selected columns/cells
directly to Parquet without collecting the full proximity table in R.

The reader validates dataset identity, partition checksums, and paths within
the bundle directory, then filters cells/markers before collecting at most the
row limit plus one. File-system traversal and escaping symlinks are rejected.
The worker also checks the load-time RDS checksum pinned to the user's request,
so replacing an RDS under the same dataset ID cannot silently run a different
analysis. RDS inputs and sidecars must still come from trusted preparation.

### Node-only 3D sidecar

The graph directory contains one Parquet file per prepared sample and a small
`index.rds` with `schema_version = 3L`. This graph-index version is separate from
the application bundle version and the frozen-report contract version.

Parquet contains exactly these columns:

| Column | Contents |
| --- | --- |
| `cell_id` | Exact current cell ID in the application bundle |
| `x`, `y`, `z` | Finite prepared 3D graph-layout coordinates |
| `markers` | List of marker identities present at that node; no counts |

Source node IDs are used only during preparation to join marker identities and
reject duplicate layout nodes, then discarded before writing Parquet. Hover
shows marker identities; node CSV exports contain the same five columns above.
Keep the source PXL for node-level traceability: exported rows no longer carry
source node IDs. Existing version 2 Parquet files remain readable; Arrow excludes
their node-ID column before loading, and the displayed asset uses version 3.

All source layout nodes are retained, including nodes without a highlighted
marker. No edges, marker-count matrices, dataset metadata, or repeated source
provenance are stored in the Parquet rows. The index holds the dataset identity,
prepared coverage, and two tables:

| Index table | Columns |
| --- | --- |
| `cells` | `cell_id`, `source_cell_id`, `sample`, `n_nodes` |
| `samples` | `sample`, `file`, `md5`, `size_bytes`, `source_pxl`, `source_pxl_md5`, `layout`, `graph_projection`, `n_nodes`, `n_cells` |

The source cell mapping is checked against the bundle's exact FSMap-derived
mapping. Source PXL checksums are checked against the frozen inventory when it
is present, without opening the original PXL files. Source/layout provenance is
stored once in the index and included in the selected cell's JSON export.
The demo's observed source layout is `coarsened_pmds_3d` with
`graph_projection = "full"`. This corrects the earlier `wpmds_3d` provenance
label; preparation preserves the source coordinates rather than recomputing them.

Runtime queries filter `cell_id` through Arrow before collecting nodes. Parquet
row-group pruning avoids collecting the whole sample. The first access to a
sample may read its entire file to validate the checksum; subsequent selections
reuse checksum validation for the unchanged file. That initial disk scan is
different from loading all sample nodes into R. The viewer rejects cells above
200,000 nodes instead of silently truncating them.

**Cell 3D** uses CRAN `threejs` for rotation, zoom, hover, one highlight marker, and
node/provenance downloads. It has no edge display. Highlights represent marker
identities, not signal counts; shared-node highlights do not constitute a
statistical colocalization test. Coordinates are abstract graph-layout positions,
not physical microscopy positions or measured distances.

The viewer keeps all nodes and uses one common coordinate scale, not separate
axis stretching. **Reset view** restores the camera; **Save PNG** captures the
canvas. Switching the highlight for the same cell preserves the camera, while
changing cells resets it. Old WebGL resources are released on each Display.
Plotly remains the renderer for embeddings and precomputed results.

**Download rotating GIF** exports the displayed cell/highlight through one full
turn around the current camera's up axis, without changing the visible camera.
Expand **Rotating GIF settings** for duration (1–20 whole seconds), frame rate
(1–30 FPS), and square image size (128–1024 pixels per side). Defaults are
6 seconds, 15 FPS, and 600 × 600 pixels. GIF frame delays have centisecond
precision; rounding is distributed so the full animation retains its duration.
The export has a white background and loops continuously. GIF's 256-color
palette may slightly change the appearance of translucent points.

Encoding runs in a browser worker using locally bundled
[gifenc 1.0.3](https://github.com/mattdesl/gifenc) (MIT; license in
`www/vendor/gifenc-LICENSE.txt`), with one captured frame in flight at a time.
No external CDN, PXL read, server-side analysis, or new R dependency is needed.
Progress and **Cancel GIF** remain available; Display is disabled until the
export ends. Replacing/unbinding the graph cancels its export. Limits of
300 frames, 120 million total pixels, and 64 MiB encoded output bound work and
memory; reduce duration, FPS, or image size if a limit is reached.

Choose donor, cell type, condition, and two ranking markers (the same marker
twice for a self-pair). Reversed pair orders are equivalent. Ranking uses finite
observed log2 ratios, descending with cell-ID tie breaking; it does not apply
abundance/differential-test cutoffs or replace missing scores with zero.
A progress indicator appears while ranking. Choose top 10, bottom 10, their
deduplicated union, or all scored cells; download the current ranking as CSV.
Cells without prepared graphs are identified in the list. Clear both ranking
markers to browse prepared cells without ranking.

Select one highlight marker and click **Display** to load the selected cell.
Changing filters, ranking, cell, or highlight does not redraw the last displayed
graph. Node and provenance downloads refer to that displayed snapshot.

Preparation can target selected cells. A later invocation retains the union of
previously prepared and newly requested cells, atomically rewrites affected
sample files, and updates the index. Unaffected, validated sample files are
reused. Run preparers sequentially against a given directory. Version 1 per-cell
RDS graph indexes are rejected with a rebuild instruction: preserve the old
directory and prepare version 3 assets in a new one. Version 2 sample Parquet
sidecars are read-only: rebuild in a new directory to remove their stored node
IDs. The preparer does not perform an automatic destructive migration or mix
schema versions within one graph directory.

## Dependencies and preparation

Install the CRAN renderer with `install.packages("threejs")` (version >= 0.3.4).

The plain-bundle viewer uses the packages in `DESCRIPTION`'s `Imports`; it
does not require Seurat/PNA class definitions merely to load data. Legacy
Seurat RDS inputs still require SeuratObject, and PNA inputs require pixelatorR.
Background differential analysis requires pixelatorR, Seurat, and SeuratObject;
spatial querying additionally requires Arrow and dplyr. **Cell 3D** also requires
Arrow and dplyr in the Shiny R environment, but not Seurat or pixelatorR for a
plain bundle. Offline graph preparation requires the source Seurat/PNA object,
pixelatorR, Arrow, and dplyr; PXL-to-Parquet preparation also uses DBI and DuckDB.
Neither the analysis worker nor the 3D viewer needs PXL mounts or permission to
access the original files.

### Renderer comparison (2026-10-12)

The saved measurements are in
[`docs/benchmarks/2026-10-12-cell-renderers/`](benchmarks/2026-10-12-cell-renderers/):
[R timings and sizes (CSV)](benchmarks/2026-10-12-cell-renderers/r-results.csv),
[browser measurements (JSON)](benchmarks/2026-10-12-cell-renderers/browser-results.json),
and [test cases and package versions (JSON)](benchmarks/2026-10-12-cell-renderers/cases.json).
These are unchanged copies of the measured results, including the optional
OrbitControls comparison. The `file` fields name generated HTML snapshots in
`output/renderer-benchmark`, not files included in this small results archive.

Compared CRAN `threejs` 0.3.4 with the existing Plotly 4.12.1 renderer from commit
`bb341b2`, using three real demo cells, all nodes, the same CD44 highlight, marker
hover labels, preserved aspect ratios, and 1000 × 700 widgets. The cells are the
smallest, median, and largest among those within the unchanged 200,000-node limit.

| Nodes | Plotly browser ready | threejs browser ready | Plotly / threejs HTML payload |
| ---: | ---: | ---: | ---: |
| 10,021 | 1,165 ms | 159 ms | 1.65 / 0.88 MB |
| 41,111 | 1,636 ms | 347 ms | 6.66 / 3.48 MB |
| 197,680 | 5,439 ms | 713 ms | 32.29 / 16.92 MB |

These are observations from one controlled browser pass, not end-to-end Shiny
latency or hardware-GPU guarantees. Browser-ready time is navigation to two
animation frames after the htmlwidget render hook. Each case ran in a fresh
Chrome 149.0.7827.55 process, after warming dependencies on the smallest cell,
with back/forward caching disabled. The server has no available hardware WebGL:
the test used SwiftShader with four workers, CPU affinity 0–3, single-threaded
V8, and device scale 1. Shared-library download sizes are excluded from the HTML
payload; MB here means decimal megabytes.

During a two-second mouse drag, animation-frame callback rates were approximately
5.8 / 60.0, 2.0 / 27.8, and 0.5 / 7.8 Hz (Plotly / threejs). This measures browser
responsiveness, not an independent GPU frame counter. The post-GC JavaScript heap
was 10.8 / 4.2, 17.4 / 9.0, and 51.7 / 33.7 MiB; it excludes native/GPU memory.
The alternative threejs OrbitControls were also tested; the app retains the
default TrackballControls for unrestricted rotation. The app disables residual
camera inertia so Reset stays at the initial view, and wakes rendering for idle
zoom events; these small interaction fixes are not part of the renderer benchmark.

R construction is not uniformly faster: medians of three warmed builds were
146 / 134, 560 / 592, and 2,687 / 3,097 ms. The largest threejs R widget was also
larger (15.85 versus 14.31 MB). The switch therefore targets browser responsiveness
and the roughly 47% smaller serialized plot payload, not R-side memory savings.
No downsampling, sidecar schema change, or node-limit increase was used.

To reproduce, run `Rscript scripts/benchmark_cell_renderers.R` from the repository
root with the demo sidecars and CRAN threejs installed. Serve
`output/renderer-benchmark` locally on port 3840, then run
`PROXIOME_CHROME=/path/to/chromium PROXIOME_BENCH_CPUS=0-3 node tests/browser_renderer_benchmark.mjs`.
Omit the CPU-affinity variable outside Linux or choose valid CPUs for the host.
The scripts save raw R and browser measurements under `output/renderer-benchmark`.
The benchmark needs the baseline commit in local Git history and Node.js with
built-in `fetch` and `WebSocket`. The app does not run these benchmarks.

After preparing the legacy frozen-results RDS, run from the repository:

```sh
Rscript scripts/prepare_hybrid_bundle.R output/proxiomevis_demo_v1.rds output/hybrid_demo
```

This writes `analysis.rds` and prepares sidecars. The optional third argument is
the source PXL directory; an optional fourth argument is a text file containing
one exact graph cell ID per line. It limits only graph preparation, not the
reports or proximity sidecar. See [Getting started](getting-started.md) for the
selected-cell command. To prepare only a proximity directory from a compatible
source bundle:

```sh
Rscript scripts/prepare_proximity_store.R output/proxiomevis_demo_v1.rds output/proximity
```

The standalone proximity preparer refuses a nonempty output directory. Preserve
the manifest with its partitions when moving a prepared dataset. The deployed
layout must keep asset paths relative to `analysis.rds`; original PXL files are
not uploaded through the application.
