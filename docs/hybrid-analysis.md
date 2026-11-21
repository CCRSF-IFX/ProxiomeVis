# Results and new comparisons

Use **Results** to view saved analyses and completed session comparisons.
Use **New comparison** to run a within-donor or an across-donor analysis.
Both result sources use the same tables, volcano, and spatial views.
Display filters do not change stored statistics. New spatial tests read the
prepared proximity files, not the source PXL files.

For your own dataset, use [Prepare your data](prepare-user-data.md).
The general preparation script does not require the demo report folders.

## Run a new comparison

1. Open **New comparison**.
2. Select **Analysis scope**, **Compare**, and a readout. Use the table below to select the groups.
3. For a spatial test, set **Minimum observed cells per donor and group** if necessary.
4. Read the readiness message. Check the target and reference cell counts, protein count, and required data.
5. For within-donor exploration, confirm the exploratory-analysis notice. Click **Run comparison**.

| Analysis scope | Compare | Select |
| --- | --- | --- |
| Within-donor exploration | Conditions | One donor, one cell type, and two conditions |
| Within-donor exploration | Cell types | One donor, one condition, and two cell types |
| Across donors, paired | Conditions | Matching donors, one cell type, and two conditions |
| Across donors, paired | Cell types | Matching donors, one condition, and two cell types |
| Across donors, independent | Conditions | Different donors in each condition, one cell type, and two conditions |

Within-donor exploration is available for prepared results and new comparisons.
It does not require biological replicates. Its cell-level tests are exploratory;
they do not establish that an effect is reproducible across donors.

Use **Paired** when the same donors provide both groups. Donors without cells in
both groups are excluded and listed. Each feature can have fewer complete pairs
after score filtering. Use **Independent** only when the two groups have no
donors in common. The app rejects overlapping donor groups for this design.
Cell-type comparisons use paired donors, or one donor for exploration.

The test uses all dataset markers and all cells in the selected scope.
There is no marker-test selector or cell-subset selector. The embedding does
not select cells for this test. Target and reference must be different groups.

Within-donor abundance tests require stored raw counts. Across-donor abundance
tests use stored normalized values and do not require counts. Spatial tests require a configured
proximity sidecar. Readiness checks the scope and input declarations before
the worker starts. The worker checks sidecar integrity and available score
coverage. A configured sidecar does not guarantee that all pairs have enough
observations.

A separate R process runs the test. You can use the other app views while it
runs. Click **Cancel running analysis** to stop it. A control change does not
start another test or change a completed result.

Only one analysis runs at a time per application process. The app permits at
most 2,000 cells per donor, 40 selected donors, 80,000 cells in total, 160 markers,
and 120 seconds per comparison. A within-donor comparison has a 2,000-cell limit
across its two groups. Spatial tests use batches
of at most 1,000 pairs and 1,000,000 score rows. These are batch limits, not a
1,000-pair limit on the complete comparison. Colocalization tests all pairs
of different markers. For example, 46 markers produce 1,035 pairs and can use
more than one batch. The app does not silently remove cells to meet a limit.
For data above these cell or protein limits, use your analysis workflow outside
the app and import validated result tables. The offline command below uses the
same cell, donor, protein, and batch limits.

## Open and keep completed comparisons

The **Completed comparisons in this session** list retains each completed test.
Select a comparison, then click **Open in Results**. You can also use
**Result source** in Results to switch between prepared results and session
comparisons. This selection does not change the loaded dataset or the Cell 3D
view. A failed or cancelled test does not remove earlier completed comparisons.

Use **Results CSV**, **Analysis provenance JSON**, or **Comparison RDS** to
keep the selected comparison. The JSON includes the submitted scope and cell
IDs. The RDS includes results, settings, and donor summaries for an across-donor
analysis. It does not include the source Seurat object or cell-level measurements. See
[Prepare your data](prepare-user-data.md) for offline reuse of a comparison.
Use only trusted RDS files.

Loading another dataset clears session history and cancels the active test.
A late result from the old dataset is not added to history. Closing the session
also removes its history. Download results that you need to keep.

## Compare spatial patterns

Open **Results → Spatial patterns**. Select a spatial readout and one cell type.
Condition profiles can be viewed without a differential test when the bundle
contains the required summaries.

1. Open **Shared proteins and order**.
2. Select **Proteins**. Choose the prepared order, alphabetical order, or a custom ordered list.
3. For a custom order, enter protein names separated by commas or new lines. This list sets both selection and order.
4. Click **Apply**. Click **Reset** to restore the dataset defaults.

All spatial matrices use the applied protein list on both axes. The selection
is retained when you change the result source, condition, comparison, or matrix
quantity. Loading a different dataset restores its defaults. Protein selection
and ordering do not rerun tests or change adjusted P values.

Choose **Single** for one matrix. Choose **Compare profiles and difference**
for the reference profile, target profile, and differential effect. The panels
use the same protein positions and matrix dimensions. On a narrow screen, the
panels are stacked.

The two profiles share a fixed color scale. Differential effects use a separate
fixed scale. New across-donor spatial profiles use the mean of eligible donor
medians. Each donor has equal weight. Paired profiles and effects use only
complete eligible donor pairs for each feature. Group eligibility counts can
include a donor with a measurement in only one group; `n_pairs` records the
complete pairs used for the test.

Older prepared profiles can use mean proximity scores. Their stored spatial
effects can be differences in medians. Those profile means do not reconstruct
the stored effects. The app keeps these older summary methods unchanged.
Diagonal entries show self-proximity. Other entries show colocalization.
The pairs A–B and B–A identify the same feature.

Unavailable or untested pairs keep their positions and appear as gray cells.
They are not replaced with zero. Hover for values and availability. Click a
pair to select it for inspection in the other result views.

Use the matrix or comparison SVG and CSV buttons to export the applied view.
Use **Protein order CSV** to export the ordered list. The selection does not
change the underlying saved tables or their full-row exports.

## Methods and interpretation

### Within-donor exploration

| Readout | Effect: target minus reference | Test and default minimum support |
| --- | --- | --- |
| Abundance | Mean stored normalized abundance | `pixelatorR::RunDAA`, using Seurat's two-sided Wilcoxon test; at least 3 cells per group |
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
each submitted request. This family can differ from the saved analysis family.
Adjusted P values can therefore differ from saved results even when effects
and raw P values agree. Display filters and shared protein controls do not
change either testing family.

These are exploratory comparisons of cells within **one donor**. Cells are not
independent biological replicates, and the results do not establish an effect
across donors. A small P value does not remove this limitation.

### Across-donor comparisons

The analysis first creates one value for each donor, group, and feature:

- Abundance: mean stored normalized abundance.
- Self-proximity and colocalization: median retained, finite, observed `log2_ratio`.

Cells from samples of the same donor and group are combined before this summary.
Additional samples do not increase the number of biological replicates. Donors
have equal weight, even when they have different cell counts. Abundance requires
at least three cells per donor and group. Spatial minimums apply to observed
scores for each marker or pair after filtering.

For a **paired** comparison, the effect is the mean of donor target-minus-reference
differences. A two-sided one-sample t-test tests those differences against zero.
The displayed minimum and maximum effects are the range of donor differences,
not a confidence interval. Only complete, eligible donor pairs enter the test.

For an **independent** comparison, the effect is the target donor mean minus the
reference donor mean. A two-sided Welch t-test compares the donor values. There
are no paired differences or donor-effect ranges. Different donors must be
independent. Paired tests also require independent donors and approximately
normal donor differences. Welch tests assume approximately normal group values;
they do not require equal group variances.

The software needs at least two complete pairs, or two eligible donors in each
independent group, to attempt a test. This is a computational minimum, not
evidence of adequate power or reliable assumptions. Inspect the donor points.
Use caution with few donors or strong outliers. A test that lacks support or
cannot return a finite P value has descriptive results and `NA` P values.

BH correction applies to the eligible tested features in that request. Changing
the donor set requires a new comparison. Display filters do not recalculate
P values. The model does not adjust for age, race, or other covariates.

## Missing scores and filtering

Preparation retains observed cell–pair rows. Before a spatial test, both
markers must pass the donor's pinned background-fraction threshold, and the
score must be finite. A missing, filtered-out, or non-finite observation is
never converted to zero. A measured score of zero remains an observation.

The demo's threshold is computed from **all prepared cells in that donor** as
`max(median(isotype_fraction) + 2 * mad(isotype_fraction), 0.001)`. It is fixed
during preparation. Display filters do not recompute it. Both demo donors have
a cutoff of 0.001. The source pipeline's median/MAD settings are read from its
stored report configuration. The current preparer supports this isotype mode
without extra marker-count or pair-count filters; incompatible configurations
are rejected instead of silently changed.

Within-donor results distinguish `n_cells_a/b` (all cells in scope) from `n_scored_a/b`
(retained observations for this particular feature). Spatial eligibility uses
the latter. Unsupported features have missing statistics and an explicit
status. Stored spatial report profiles may use a different missing-observation
summary convention; those means are never used to reconstruct new tests.

Across-donor results also record eligible donor counts for each group and the
number of complete pairs. An unavailable donor measurement stays `NA`; it is
not a zero. A measured score of zero remains valid. The donor plots connect
paired measurements only. Independent-group points have no connecting lines.

## Prepare a comparison outside the app

Use the same calculation code to prepare a result before launch. The request
is a JSON file. Use exact donor, condition, and cell-type names from the bundle.
These examples test all dataset markers. Add `"markers": ["CD19", "CD82"]`
only when you need a restricted test family outside the app.

Paired conditions, saved as `paired_conditions.json`:

```json
{"analysis_scope":"across_donors","comparison_axis":"condition","design":"paired","readout":"abundance","celltype":"B","target":"treatment","reference":"control","donors":["D1","D2","D3"]}
```

Paired cell types in one condition:

```json
{"analysis_scope":"across_donors","comparison_axis":"celltype","design":"paired","readout":"self_proximity","celltype":"T","condition":"control","target":"T","reference":"B","donors":["D1","D2","D3"],"min_cells_per_group":20}
```

Independent conditions. The selected donors must belong to separate groups:

```json
{"analysis_scope":"across_donors","comparison_axis":"condition","design":"independent","readout":"colocalization","celltype":"B","target":"treatment","reference":"control","donors":["C1","C2","T1","T2"],"min_cells_per_group":40}
```

Run the command from the repository:

```sh
Rscript scripts/run_comparison.R output/my_dataset/analysis.rds \
  paired_conditions.json paired_comparison.rds \
  --id paired-treatment-B --label "Paired treatment comparison: B cells" \
  --output-bundle output/my_dataset/with_comparison.rds
```

The three positional arguments are the input bundle, request JSON, and new
comparison RDS. `--id` and `--label` are optional. `--output-bundle` also appends
the comparison to a new application RDS. Existing results are retained. The
import validates the result and rejects a duplicate analysis ID. Existing files
are not overwritten. With sidecars, keep the new bundle beside the input bundle.
The command has the same calculation and size limits as the app. It does not
have the Shiny worker's 120-second timeout. Use only trusted input files.

## Application bundle and sidecar schemas

The version 2 application RDS is a plain R list, not a Seurat object:

- `bundle_version = 2L`; `contract` contains report settings, defaults, QC, and
  provenance. The reader accepts analysis `schema_version = 1L` or `2L`.
  New donor-level and cell-type comparisons use schema version 2. Old version 1
  result tables retain their original values and meaning.
- `metadata` and `embeddings` contain prepared cell annotations and coordinates.
- `abundance` and optional `counts` are aligned marker-by-cell sparse matrices.
  Counts are required for within-donor abundance tests and proximity preparation,
  not for across-donor abundance tests.
- `cell_sources`, when present, maps `component`, `original_id`, and `sample_id`
  exactly from the original FSMap; IDs are never inferred by stripping prefixes.
- `assets$proximity = list(path = "proximity", manifest = "manifest.rds")`
  locates the proximity sidecar relative to the bundle directory.
- `assets$graphs = list(path = "graphs", index = "index.rds")` locates the
  optional node-only 3D sidecar.

Schema version 2 adds comparison-axis and donor-design fields. Its optional
`donor_values` section contains a named table for each readout. Each row stores
one donor, group, and feature value, cell counts, and an eligibility flag.
`source_celltype` identifies the actual source cell type. The `celltype` field
anchors a cell-type comparison to its target cell type. `fixed_condition`
identifies its common condition. `donor_contrasts` stores paired donor
differences or within-donor exploratory results. `differential_results` stores
across-donor effects and tests. Independent comparisons have no paired rows;
their `n_pairs` is zero and their donor-effect range is `NA`.

### Proximity sidecar

The proximity directory contains a version 1 `manifest.rds` and one Parquet
partition per sample. Each row contains `component`, canonical `marker_1` and
`marker_2`, `log2_ratio`, and `min_marker_fraction` (the smaller of the two marker
counts divided by that cell's total count). Storing the minimum is equivalent
to requiring both marker fractions to pass the same threshold.

`component` contains the current application cell ID. Preparation converts the
source PXL IDs with the exact per-sample cell mapping. The reader preserves these
IDs and selects only the requested application cells; it does not infer the ID type.

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

Install the packages in `DESCRIPTION`'s `Imports`, including Seurat.
A plain application bundle does not contain Seurat or PNA objects. Legacy
Seurat RDS inputs require SeuratObject, and PNA inputs require pixelatorR.
Within-donor tests require pixelatorR, Seurat, and SeuratObject. Across-donor
abundance uses stored normalized values. Spatial tests also need the proximity
reader packages below.
Spatial queries and **Cell 3D** require Arrow and dplyr in the Shiny environment.
Offline graph preparation requires the source Seurat/PNA object,
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
