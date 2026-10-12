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
prepared sample plus a version 2 index, not per-cell RDS files. Source and layout
provenance belong in that index; Parquet rows contain only cell/node IDs,
coordinates, and marker identities. Rebuild old version 1 graph exports in a
new directory rather than changing their version field. See
[Getting started](getting-started.md) for startup and
[Hybrid analysis](hybrid-analysis.md) for sidecar schemas and dependencies.

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

## What donor and cell-type selection means

In the report tabs, selections narrow the stored rows. They do not reopen PXL, recalculate a
proximity score, rerun marker eligibility filters, or recompute a donor-matched
aggregate for the selected subset. Consequently, a result can only be shown if
that donor/cell-type/marker combination was present when the frozen tables were
built.

This analysis contains two donors. Donor points and treatment-minus-reference
effects are useful for exploratory comparison, but two donors do not support
population-level significance claims. The app therefore presents descriptive
combined effects without p-values. Donor-specific report volcanoes may display
the already stored component-level p-values described above. Report controls
never run a new test. **Custom analysis** instead submits an explicit selected-cell
comparison to a background R process, with its own result and provenance. Those
new p-values remain exploratory: cells within a donor are not independent
biological replicates, and selection using the same data can bias testing.

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
