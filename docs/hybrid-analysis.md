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
`index.rds` with `schema_version = 2L`. This graph-index version is separate from
the application bundle version and the frozen-report contract version.

Parquet contains exactly these columns:

| Column | Contents |
| --- | --- |
| `cell_id` | Exact current cell ID in the application bundle |
| `node_id` | Source layout node ID, unique within its cell |
| `x`, `y`, `z` | Finite prepared 3D graph-layout coordinates |
| `markers` | List of marker identities present at that node; no counts |

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

**Cell 3D** supports rotation, zoom, hover, up to three marker highlights, and
node/provenance downloads. It has no edge display. Highlights represent marker
identities, not signal counts; shared-node highlights do not constitute a
statistical colocalization test. Coordinates are abstract graph-layout positions,
not physical microscopy positions or measured distances.

Preparation can target selected cells. A later invocation retains the union of
previously prepared and newly requested cells, atomically rewrites affected
sample files, and updates the index. Unaffected, validated sample files are
reused. Run preparers sequentially against a given directory. Version 1 per-cell
RDS graph indexes are rejected with a rebuild instruction: preserve the old
directory and prepare version 2 assets in a new one. The preparer does not
perform an automatic destructive migration.

## Dependencies and preparation

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
