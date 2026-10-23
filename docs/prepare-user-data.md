# Prepare your data

Use this procedure for your own annotated Seurat object. The demo report folders are not required.

Run preparation on a trusted local computer. Use only RDS files from a trusted source. An RDS file is not a safe public-upload format.

## Required input

The Seurat object must contain:

- A protein assay with a normalized abundance layer.
- Unique cell IDs and protein names.
- Metadata for sample, donor, condition, cell type, and cluster. Each sample must have one donor and one condition. These values must not be missing.
- A default embedding and a PCA embedding. Each embedding must have at least two finite coordinates for each cell.

A counts layer is optional for viewing. Include counts to run abundance comparisons or prepare proximity files. Counts must be nonnegative integers. Counts and normalized abundance must use the same proteins and cells.

The script does not normalize values, assign cell types, or run differential tests. Complete these steps in your analysis first. Join split assay layers before preparation. Do not use a raw counts layer as the normalized abundance layer.

## Prepare the application folder

From the ProxiomeVis repository, run:

```sh
Rscript scripts/prepare_user_bundle.R my_analysis.rds output/my_dataset \
  --title "My PNA dataset" \
  --assay PNA \
  --abundance-layer data \
  --reduction umap \
  --pca-reduction pca \
  --celltype-column celltype \
  --cluster-column seurat_clusters \
  --donor-column donor \
  --condition-column condition \
  --sample-column sample_id
```

Change the names to match your object. Run `Rscript scripts/prepare_user_bundle.R --help` to see all options. The script also works from another folder if you use its full path.

Use a new output folder. The script does not overwrite an existing folder or the input RDS. It prepares files in a temporary folder. It moves the complete folder to the output path only after validation succeeds.

The output contains `analysis.rds`. This version 2 application bundle contains metadata, embeddings, abundance, optional counts, QC counts, and result tables. It does not contain the source Seurat object or its cell-level proximity table.

Start the app with:

```sh
PROXIOME_PRECOMPUTED_RDS=output/my_dataset/analysis.rds Rscript -e 'shiny::runApp()'
```

Without saved comparisons, the bundle has empty differential-result tables. Abundance profiles use the mean stored normalized value for each sample and cell type. Their metric is `mean_normalized_abundance`; the script does not assume that your values use CLR normalization. Spatial profiles need imported summaries or PXL preparation. New comparisons need the applicable counts or proximity files.

## QC counts

By default, `00_loaded` is the number of cells in the input object for each sample. It is not the number of cells before your earlier analysis. Counts for earlier filter steps are unavailable, not zero.

To include known filter counts, use `--qc qc_counts.tsv`. The TSV must contain:

| Column | Content |
|---|---|
| `sample_id` | Sample ID from the mapped metadata column |
| `donor` | Donor for that sample |
| `condition` | Condition for that sample |
| `stage` | One of the four stage names below |
| `n_cells` | Nonnegative integer cell count |

Stage names are `00_loaded`, `01_after_n_umi_filter`, `02_after_tau_filter`, and `03_after_isotype_filter`. Include each sample and stage at most once. Do not add filter counts that you did not measure.

## Import saved results

Use `--saved-results saved_results.rds` to import canonical ProxiomeVis result sections. This file is an R list. It can contain `sample_summaries`, `donor_contrasts`, and `differential_results`. Each section is a named list of data frames for `abundance`, `self_proximity`, and `colocalization`. You can omit a section or readout. Omitted differential tables stay empty. The default abundance profile is retained unless you supply an abundance summary table.

Use the table templates to see the required columns and R types:

```r
source("R/precomputed_contract.R")
source("R/bundle_builder.R")
saved <- empty_user_results()
str(saved$sample_summaries$colocalization)
str(saved$donor_contrasts$colocalization)
str(saved$differential_results$colocalization)
```

Common fields are `celltype`, `feature_id`, `marker_1`, `marker_2`, `metric`, and logical `is_control`. Use the same feature ID for a pair in all sections. Store the pair in canonical order: `marker_1 <= marker_2`. Self-proximity and abundance use the same marker twice. Colocalization uses two different markers.

- `sample_summaries` contains the sample, donor, condition, value, and cell-support counts.
- `donor_contrasts` contains the analysis ID, donor, comparison groups, group values, effect, and group cell counts. `group_a` is the target. `group_b` is the reference. If both group values are available, the effect must equal `value_a - value_b`.
- `differential_results` contains the analysis ID, groups, effect range, number of pairs, P values, and method fields. A descriptive result must use `NA_real_` for `p_value` and `p_adj`.

Use real analysis results for these fields. Do not fill missing test results with zero. The contract validator checks the required columns, pair types, sample design, probability bounds, and support counts.

You can copy canonical result sections from an existing bundle:

```r
old <- readRDS("existing/analysis.rds")
saved <- old$contract[c("sample_summaries", "donor_contrasts", "differential_results")]
saveRDS(saved, "saved_results.rds")
```

Only use results that match the new object's cells, sample design, and analysis. The import can also contain `report` and `provenance` entries from an existing ProxiomeVis contract. A native pixelatorR result table alone is not this format; map its columns to the templates first. The builder does not infer the meaning of arbitrary table columns.

## Prepare proximity and 3D files

To prepare proximity files, add `--pxl-manifest pxl_files.tsv`. This TSV has two columns:

```text
sample_id	pxl_file
sample_A	data/sample_A.pxl
sample_B	data/sample_B.pxl
```

Relative PXL paths are relative to the manifest file. Include one PXL file for every sample in the object. The PXL files must contain the original proximity measurements.

Cell IDs must have an exact mapping to the source PXL files. The builder uses the protein assay's stored FSMap if it is available. Otherwise, supply `--cell-map cell_map.tsv` with these columns:

```text
component	original_id	sample_id
current_cell_A	original_pxl_cell_A	sample_A
current_cell_B	original_pxl_cell_B	sample_B
```

Include every cell in the input object. `component` is the current Seurat cell ID. `original_id` is the original PXL cell ID. A current ID must be unique. An original ID must be unique within its sample. The builder does not guess cell mappings.

The counts layer and an isotype-fraction metadata column are required. Use `--isotype-column NAME` if the column is not named `isotype_fraction`. With no imported filter settings, the donor cutoff is the larger of `median + 2 × MAD` and `0.001`. The prepared file records the minimum marker fraction for each pair. The app applies the cutoff when it reads scores. Existing supported report filter settings are retained.

PXL preparation also creates spatial profiles for each sample and cell type. The profile is the sum of finite scores that pass the cutoff, divided by the number of cells in that sample and cell type. Missing scores contribute zero to this population mean. The summary records the number of observed scores.

If a pair has a passing score in a donor and cell type, the output includes that pair for every sample in the same donor and cell type, across conditions. A sample with no passing scores has a zero mean, zero observed scores, and its full cell count. This keeps its cells in the denominator when the app combines samples. For example, 10 cells with a score sum of 10 and 30 cells with no scores give a combined mean of `10 / 40 = 0.25`, not `1`. This zero is from the stated missing-score policy; it is not evidence that a score was measured. A pair with no passing scores anywhere in that donor and cell type stays unavailable. The builder does not create every possible pair.

The metric is `population_mean_log2_ratio_missing_zero`. These profiles do not require a differential test. Imported nonempty spatial summaries take priority over these generated profiles. Include the complete sample contributions when you import summaries with this metric; an omitted sample row does not establish its denominator or prove a zero mean.

Add `--graphs true` to prepare 3D layout files. This option also requires pixelatorR, a source object with an FSMap, and all PXL files in one source directory. It exports layouts for all input cells. It can take substantial disk space and preparation time. A cell-map TSV alone does not replace the FSMap for graph preparation.

Example:

```sh
Rscript scripts/prepare_user_bundle.R my_analysis.rds output/my_dataset \
  --cluster-column cluster \
  --pxl-manifest pxl_files.tsv \
  --cell-map cell_map.tsv \
  --isotype-column isotype_fraction
```

The bundle stores relative paths to `proximity/` and, when requested, `graphs/`. Keep these folders next to `analysis.rds` when you deploy or move the dataset. The app reads the prepared files; users do not need to upload the original PXL files.

## Use the R function

```r
source("R/precomputed_contract.R")
source("R/proximity_store.R")
source("R/cell_graph_assets.R")
source("R/bundle_builder.R")

path <- prepare_user_bundle(
  input_rds = "my_analysis.rds",
  output_dir = "output/my_dataset",
  title = "My PNA dataset",
  celltype_column = "celltype",
  cluster_column = "cluster"
)
```

For an in-memory object, call `user_bundle_data(object, ...)`. It returns validated application data. `empty_user_results()` returns the import templates. Optional `control_markers` lists known isotype controls; the default is an empty list. The builder does not identify control proteins by name.

## Reuse a comparison exported from the app

The comparison-history download is a separate format from the canonical result-section file above. To save it as an application bundle, run this trusted offline import:

```r
source("R/precomputed_contract.R")
source("R/hybrid_analysis.R")
source("R/hybrid_module.R")
source("R/bundle_builder.R")

import_user_comparison(
  bundle_path = "output/my_dataset/analysis.rds",
  comparison_path = "downloaded_comparison.rds",
  output_rds = "output/my_dataset/with_comparison.rds"
)
```

Use the exact input bundle from that session. The import checks the dataset ID, bundle checksum, submitted scope, and result fields. The new file contains the selected comparison in place of the earlier differential comparisons. It retains metadata, abundance, profiles, and sidecar paths. The original file is unchanged. Save the new RDS beside the original when the dataset has sidecars. Use the new RDS path to launch the app.
