# Precomputed results

The report views consume a versioned frozen-analysis contract. The recommended
application file is a compact plain-R bundle; legacy Seurat RDS files remain
supported. Filtering a report never changes its stored statistics. The separate
**Custom analysis** tab can submit a new bounded comparison to an R worker using
prepared cell-level inputs; see [Hybrid analysis](hybrid-analysis.md).

## Prepare the frozen results

Run the builder once after the upstream analysis has completed:

```bash
Rscript --vanilla scripts/prepare_precomputed_rds.R \
  ../pixelgen_demo_data/demo_analysis/cca_analysis/results/demo_cca_clustered.rds \
  ../pixelgen_demo_data/demo_analysis/tertiary_analysis \
  output/proxiomevis_demo_v1.rds
```

The arguments are:

1. `INPUT_RDS`: the annotated Seurat object used by the analysis.
2. `ANALYSIS_DIR`: the analysis directory containing the frozen result tables.
3. `OUTPUT_RDS`: the legacy frozen-results Seurat RDS; normally
   `output/proxiomevis_demo_v1.rds`.

Optional settings follow the paths as `--option VALUE`. Omit an option to use its default:

| Option | Default |
| --- | --- |
| `--umap` | `umap.cca` |
| `--pca` | `pca` |
| `--celltype-column` | `celltype_broad` |
| `--cluster-column` | `cluster` |
| `--donor-column` | `donor` |
| `--condition-column` | `condition` |
| `--sample-column` | `sample_alias` |
| `--assay` | `PNA` |
| `--abundance-layer` | `data` |

These names must exist in the input object. The selected settings are saved in
the RDS contract and included in its dataset identity. For example, append
`--umap umap --celltype-column celltype --sample-column sample_id` to the command
if your input uses those names. Run the script with `--help` to list its options.
The same settings are named R arguments to `build_precomputed_rds()`; use
`reduction`, `pca_reduction`, and underscores instead of hyphens for the column/layer arguments.

These options map input fields; they do not recompute or relabel frozen report
results. The analysis directory must still follow the demo report layout, and
its provenance must match the input RDS. Selected annotations, sample identities,
and abundance values must correspond to the saved analysis.

The `output/` directory is gitignored. Build or deliver the RDS separately on
the deployment host rather than committing it with the application source.

The input object must contain cell metadata that can be mapped to these fields:

| Field | Meaning |
| --- | --- |
| `component` | Current cell/component identifier |
| `sample_id` | Experimental sample identifier |
| `donor` | Biological donor identifier |
| `condition` | Experimental condition |
| `celltype` | Cell-type annotation |
| `cluster` | Stored cluster annotation |

The builder validates its inputs and writes a separate output file. It does not
modify the input RDS or anything below the `../pixelgen_demo_data` source tree.
This builder is intentionally specific to the demo directory layout and its
`d1`/`d2` donor outputs; the versioned application contract is the reusable
boundary for other preparation pipelines.

## Export a plain application bundle

To view frozen results without Seurat/PNA classes on the deployment host,
export the prepared file in the analysis environment:

```bash
Rscript --vanilla -e 'source("R/precomputed_contract.R"); write_precomputed_bundle(read_precomputed_rds("output/proxiomevis_demo_v1.rds"), "output/proxiomevis_app_v2.rds")'
```

Set `PROXIOME_PRECOMPUTED_RDS` to the exported file when starting Shiny. This
export needs no PXL access and does not add spatial or graph sidecars. Its sparse
counts and stored abundance can support new abundance tests when the worker
dependencies are installed.

For the complete demo bundle, including proximity and cell-graph sidecars, run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 \
  Rscript --vanilla scripts/prepare_hybrid_bundle.R \
  output/proxiomevis_demo_v1.rds output/hybrid_demo
```

This preparation needs the source PXL files. The optional third argument is
their directory if they moved. The deployed output contains `analysis.rds`,
`proximity/`, and `graphs/`; keep the relative directory layout intact. Full graph
preparation can be lengthy and is resumable, but do not run two preparers in the
same output directory. The graph sidecar is one node-only Parquet file per
prepared sample plus a version 3 index, not per-cell RDS files. Source and layout
provenance belong in that index; Parquet rows contain only cell IDs,
coordinates, and marker identities. Version 2 remains readable without loading
node IDs. Rebuild older graph exports in a new directory to obtain version 3;
do not change their version field by hand. See
[Getting started](getting-started.md) for startup and
[Hybrid analysis](hybrid-analysis.md) for sidecar schemas and dependencies.

### Export only ranked cells for 3D

Append `--graph-markers` and a comma-separated marker list to the **hybrid bundle** command to limit graph export:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 \
  Rscript --vanilla scripts/prepare_hybrid_bundle.R \
  output/proxiomevis_demo_v1.rds output/hybrid_ranked \
  --graph-markers CD19,CD37,CD53,CD81,CD82
```

Within each **donor × cell type × condition**, the exporter selects:

- The top 2 and bottom 2 cells for every different-marker pair among the specified markers (colocalization).
- The union of these cells, without duplicate graph exports.

Self-marker pairs are excluded. Both markers in each pair must appear in the
list. At least two distinct known marker names are required; names are case-sensitive.
Omitting the option keeps the existing all-cell or explicit cell-list export.
Each selected cell retains its complete layout and all node-marker identities,
including markers outside the selection list. No node filtering is performed.

Ranking uses finite observed `log2_ratio` scores from the proximity sidecar.
Reversed marker orders represent the same pair. Duplicate observations for one
cell and unordered pair are rejected. Ties use ascending exact cell IDs in each
tail; groups with fewer than two scores contribute all available cells. Missing
scores are not replaced with zero. Like the application's ranking, this does not
apply the differential-test abundance cutoff. No finite scores means export
stops with an error. The top/bottom count is fixed at two, not configurable.

Two is a **per-pair, per-tail** limit, not a total-cell limit. Five markers give
ten different-marker pairs and at most 40 cells per donor/type/condition stratum.
Overlapping selections reduce this count. A small stratum can still contribute
all its cells. Ranking runs in DuckDB with one thread and a
512 MB database memory limit, using temporary disk space for large sorts; this
is not a limit on total process memory. Only selected cell IDs and flags enter R.

The bundle root also contains:

| File | Contents |
| --- | --- |
| `graph_export_summary.tsv` | Total, selected, and exported cells per sample, including zero-export samples; separate tail counts before the union |
| `graph_selected_cells.tsv` | Exported cell IDs, sample/donor/type/condition, and flags identifying the selecting tails |
| `graph_selection.rds` | The selection, summary, ranking parameters, and source proximity-manifest checksum |

The console reports exported counts after graph preparation. Tail counts overlap
and must not be added to obtain the union count. The proximity sidecar, application
metadata, and frozen results retain all cells.

Do not combine `--graph-markers` with the fourth positional cell-list argument.
Use a new output directory when existing graphs include cells outside the new
selection. Matching partial exports can resume; the exporter does not delete
previously prepared cells. Supply the ranking option again when resuming a ranked
directory. From R, use `build_hybrid_bundle(..., graph_markers = c("CD19", "CD37", "CD53", "CD81", "CD82"))`.

The plain R list has these fields:

| Field | Contents |
| --- | --- |
| `bundle_version` | `2L`, the application file format version |
| `contract` | The unchanged frozen contract with `schema_version = 1L` |
| `metadata`, `embeddings` | Prepared cell metadata and all stored reductions |
| `abundance` | Aligned sparse marker-by-cell transformed abundance |
| `counts` | Optional aligned sparse counts; needed for new abundance tests and proximity preparation |
| `cell_sources` | Optional exact `component`, `original_id`, and `sample_id` mapping from FSMap |
| `assets` | Optional relative descriptors for proximity and graph sidecars |

No Seurat/PNA object, per-cell proximity table, or 3D node table is included in
this application bundle. The reader adds `bundle_path` and a load-time
`bundle_md5` in memory; these are runtime identity fields and are not exported
into the RDS.

## Frozen analysis contract

The first builder's output remains a Seurat object, with its contract at:

```r
object@misc$proxiomevis
```

In the plain bundle the same contract is stored at `bundle$contract`.
Its version 1 fields contain:

```text
schema_version         integer contract version
dataset_id, title      stable identifier and display title
defaults               reduction, metadata-column, assay, and layer mappings
sample_summaries       one stored value per sample and feature
donor_contrasts        stored group A minus group B effect for each donor
differential_results   descriptive effects combined across donors
qc                     frozen cell-count table
provenance             source and preparation details
report                 frozen source rows and settings used by report plots
proximity_source       optional source-PXL locations for offline preparation
```

Each of the three result sections must contain these named readouts:

- `abundance` stores marker abundance summaries and contrasts.
- `self_proximity` stores marker self-proximity, where the two marker fields
  are identical.
- `colocalization` stores proximity between two different markers.

All readout tables identify the cell type, feature, ordered marker pair,
metric, and control-marker status. Sample summaries additionally retain the
sample, donor, condition, value, and total/evaluable/scored cell counts. Donor
contrasts retain the analysis ID, comparison labels, donor, source effect, cell
counts, and original exploratory test fields. Differential results retain the
effect range, donor-pair count, direction agreement, method, replicate unit,
inference label, and multiple-testing fields.

Abundance sample values use `mean_clr`, and abundance contrasts use
`delta_mean_clr`; their frozen group means are stored in `value_a` and
`value_b`. Spatial sample values use the source population mean `log2_ratio`
with the original missing-as-zero policy. Spatial effects use the original
difference of observed medians (`delta_median_log2_ratio`). Because those two
source group medians were not persisted, spatial `value_a` and `value_b` are
`NA`; `effect_size` is the authoritative stored effect and must not be inferred
from the sample-profile means.

For this two-donor dataset, the combined `differential_results` rows are
descriptive and their `p_value` and `p_adj` fields must be `NA`. The contract
validator rejects a descriptive row carrying inferential probabilities.

The optional `report` bundle is deliberately separate. It retains the original
donor-specific tertiary outputs needed to reproduce the established figures:
source result rows, eligibility tables, spatial profiles, resolved settings,
the cross-donor alignment bundle, and frozen abundance heatmap matrices. Those
source result rows can contain the tertiary analysis's component-level
`p`/`p_adj` values. The app displays them unchanged as exploratory
component-level tests; it does not reinterpret them as donor-level inference or
combine them into a new test.

In a legacy input, the Seurat object supplies per-cell metadata, embeddings, and
abundance values. The plain bundle stores those extracted values directly.
`defaults` retains the prepared reductions, metadata mappings, assay name, and
abundance-layer identity:

```text
reduction, pca_reduction
celltype_column, cluster_column, donor_column, condition_column, sample_column
assay, abundance_layer
```

The reductions must contain finite coordinates aligned to the cells. A legacy
input's configured assay must contain `counts` and the configured abundance
layer, and its cell-level proximity slots must be empty. Plain bundles validate
the extracted matrices and metadata without constructing a Seurat assay. Large
per-cell proximity tables belong in optional disk-backed sidecars, not the RDS.

Consumers should validate the contract version and required columns before
rendering any view. A newer incompatible layout should use a new version rather
than silently changing the meaning of an existing field.

## Explore one donor's stored results

The top navigation groups stored views under **Precomputed results**. Its shared
selection bar contains donor, cell type, comparison, and readout. A single-donor
dataset displays the donor as a fixed label; multi-donor datasets offer a donor
dropdown, but every result view still uses one donor's original stored contrast.
The descriptive cross-donor rows remain preserved in the contract; they are not
required or substituted for donor-specific effects by this interface.

The four subtabs are:

| Subtab | Contents |
| --- | --- |
| Summary | Differential effects for the selected cell type, or across all cell types; optional original abundance heatmap |
| Volcano | One donor and one cell type's original stored test with hover, zoom, and feature selection |
| Spatial matrix | Aligned self-proximity and colocalization, either condition profiles or differential effects |
| Tables | Searchable result and volcano rows, sample support, exports, method settings, and a combined cell-level violin/boxplot for the selected feature |

Choosing **All cell types (summary)** is supported by the heatmap and tables.
Volcano and spatial views then ask for one cell type rather than silently
choosing a different scope. Comparison labels retain target/reference direction:
a positive effect means target exceeds reference. Results come from stored
rows; these controls do not reopen PXL, recalculate proximity, rerun eligibility
filters, or fit a new test.

Click a heatmap tile, volcano point, spatial marker pair, or result-table row
to select the same feature in the detail panel shared by all four views.
Pair keys are order-independent. The panel shows the selected scope's stored
effect, available probabilities, and cell support; condition values are shown
only if present in the contract. A feature without a stored test is labeled as
such instead of being converted to a zero effect.

**Explore cells in 3D** transfers that donor, cell type, and pair. The user must
choose a condition from the comparison and one highlight marker, then click
**Display**. The handoff does not load a graph or alter the previously displayed
3D view. It can initiate ranking only after a condition is chosen and the
required proximity sidecar is available.

### Statistical meaning and spatial exploration

Original donor-specific volcanoes may display the stored component-level P
values described above. These are exploratory within-donor tests, not evidence
of a population-level treatment effect. Cells from one donor are not independent
biological replicates. The same limitation applies to a dataset containing only
one donor; no combined or pooled test is invented.

The spatial matrix deliberately combines self-proximity and colocalization,
independent of the readout control used by the other subtabs. Hover exposes the
pair, stored value, and availability; click selects the pair for detail. The
marker order and aligned scales follow the frozen report. Condition profiles
show stored population means with the upstream missing-as-zero policy, while
differential effects show differences of observed medians. The former must not
be subtracted to reconstruct the latter. Untested or unavailable entries are
distinguished from measured values near zero.

**Custom analysis** instead submits an explicit selected-cell comparison to a
background R process with its own result and provenance. Those new P values
remain exploratory, and selection using the same data can bias testing.

### Filters, exports, and rendering

FDR controls are labeled table filters; the app does not apply them to the
original report plots. When no adjusted P values were stored, these controls
are unavailable rather than fabricated. Table searches and FDR filters affect
the corresponding **Filtered table CSV**, including all matching rows rather
than just the current page. **All stored rows CSV** ignores both filters and
exports the selected donor/cell-type/comparison/readout scope. Heatmap feature
counts and control-marker display options affect only the summary display.

**Original report SVG** preserves the upstream report plotting strategy and
thresholds. The interactive spatial view also offers its displayed rows as CSV.
The methods section exposes stored effect metrics, original settings, test
methods, and multiple-testing information where available; it does not infer
missing method details.

Plot outputs are suspended while their subtab is hidden; shared selections and
the selected feature persist while navigating. Prepared views read only the
loaded RDS. Sidecars are used for explicit cell-level ranking, graph display,
or custom spatial analysis, not merely for viewing a stored result.

## Runtime and PXL files

Both supported RDS formats contain all three frozen readouts, so those views
need no mounted PXL file. New spatial tests require the prepared proximity
sidecar; **Cell 3D** requires prepared graph assets. These runtime features also
need no PXL access. Missing sidecars do not prevent frozen report viewing.
Cell 3D additionally requires Arrow and dplyr to select nodes from its sample's
Parquet file. Graph data is not collected for all cells when opening the RDS.

Only offline preparation needs the original PXL files. Plain-bundle loading
requires Matrix but not Seurat or pixelatorR. Legacy Seurat inputs still require
SeuratObject, and PNA inputs require pixelatorR's class definitions. The
`read_precomputed_rds()` helper handles supported input formats; background
testing has additional package requirements documented in
[Hybrid analysis](hybrid-analysis.md#dependencies-and-preparation).

## Security

RDS is an R serialization format, not a safe public upload format. Only load
files produced by a trusted analysis pipeline. Prefer an administrator-managed
server directory and `PROXIOME_PRECOMPUTED_RDS`; do not accept arbitrary
internet uploads and pass them to `readRDS()`.

The optional `PROXIOME_ALLOW_LOCAL_LOAD=true` flag exposes a server-path field
for a trusted local operator. Leave it unset on public deployments because it
allows the visitor to select any readable server-side RDS path.
