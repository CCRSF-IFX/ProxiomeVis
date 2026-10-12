# Precomputed results

The R app consumes a compact, versioned payload stored inside a Seurat RDS.
This keeps expensive PXL access and proximity calculations in the analysis
pipeline, while the Shiny process performs only filtering, plotting, and table
export.

## Build the runtime file

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
3. `OUTPUT_RDS`: the new runtime RDS; normally
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

## Stored contract

The output remains a Seurat object. Its application payload is stored at:

```r
object@misc$proxiomevis
```

The version 1 payload contains:

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

The Seurat object itself remains the source of per-cell metadata, embeddings,
and abundance values. `defaults` maps the app to:

```text
reduction, pca_reduction
celltype_column, cluster_column, donor_column, condition_column, sample_column
assay, abundance_layer
```

The named reductions must contain finite coordinates aligned to the cells. The
configured assay must contain `counts` and the configured abundance layer.
Cell-level assay proximity slots must be empty; only the frozen result tables
belong in the runtime RDS.

Consumers should validate the contract version and required columns before
rendering any view. A newer incompatible layout should use a new version rather
than silently changing the meaning of an existing field.

## What donor and cell-type selection means

Selections narrow the stored rows. They do not reopen PXL, recalculate a
proximity score, rerun marker eligibility filters, or recompute a donor-matched
aggregate for the selected subset. Consequently, a result can only be shown if
that donor/cell-type/marker combination was present when the frozen tables were
built.

This analysis contains two donors. Donor points and treatment-minus-reference
effects are useful for exploratory comparison, but two donors do not support
population-level significance claims. The app therefore presents descriptive
combined effects without p-values. Donor-specific report volcanoes may display
the already stored component-level p-values described above. No p-value is
calculated in Shiny, and cells within a donor are not independent biological
replicates.

## Runtime and PXL files

The runtime RDS contains all three readouts, so the app needs no mounted PXL
file. PXL files remain inputs to the upstream analysis pipeline and can stay
outside the hosted Shiny environment.

The Seurat object contains pixelatorR PNA S4 classes. Load `pixelatorR` before
calling `readRDS()` so those class definitions are registered, even though no
PXL query is performed.

## Security

RDS is an R serialization format, not a safe public upload format. Only load
files produced by a trusted analysis pipeline. Prefer an administrator-managed
server directory and `PROXIOME_PRECOMPUTED_RDS`; do not accept arbitrary
internet uploads and pass them to `readRDS()`.

The optional `PROXIOME_ALLOW_LOCAL_LOAD=true` flag exposes a server-path field
for a trusted local operator. Leave it unset on public deployments because it
allows the visitor to select any readable server-side RDS path.
