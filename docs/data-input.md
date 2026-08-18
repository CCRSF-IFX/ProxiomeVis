# Data Input

## Python app: H5AD analysis data

The Python Shiny app accepts one server-visible `.h5ad` file for all analysis
data. An optional PXL path is used only to display precomputed 3D component
layouts; the app does not rerun Pixelator analysis from PXL.

The reference dataset is:

```text
/Volumes/ccrsf-static/illumina/CCRSFIFX-23_MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/python_results/pg_data_combined_filtered_annotated.h5ad
```

Override the default path before launch when needed:

```bash
export PROXIOME_H5AD=/path/to/processed_data.h5ad
export PROXIOME_PXL='/path/to/layouts/*.layout.pxl'
```

## Expected AnnData content

- `X`: component-by-marker counts
- `obs`: sample, condition, cell annotation, and QC metadata
- `var_names`: antibody marker names
- `layers["clr"]`: normalized abundance when available
- `obsm`: stored two-dimensional or higher embeddings such as UMAP, PCA, or Harmony
- `uns["qc_cell_counts_by_step"]`: optional notebook-generated QC retention history
- `uns["proxiome"]["proximity"]`: optional table with `component`, `marker_1`,
  `marker_2`, and `log2_ratio`; optional Pixelator filter columns are
  `marker_1_freq`, `marker_2_freq`, and `min_count`
- `uns["proxiome"]["patch"]`: optional mapping containing `run_plan`,
  `marker_unmixing`, `raji_marker_abundance`, `raji_marker_proximity`, and
  `patch_burden` tables
- `uns["proxiome"]["component_layouts"]`: optional mapping from component ID
  to a node table containing `x`, `y`, and `z`

Missing `sample_alias`, `condition`, or `celltype_manual` fields receive safe
fallback values. Missing optional tables leave the corresponding module visible
with an unavailable-data message; the app never substitutes abundance
correlations for spatial proximity.

## Optional PXL cellgraph path

In **Data**, assign a `.layout.pxl` file, directory, glob, or comma/newline-
separated paths. When a cell is selected in **Spatial Metrics >
Colocalization > 3D Layout**, the app matches its sample to a PXL filename and
reads that component's stored layout on demand. PXL is not used for QC,
abundance, proximity, colocalization, or patch analysis.

## Analysis grouping

Use **Data > Analysis grouping** after loading a dataset. Select **Use metadata
column** to group by a column that is constant within each sample, or **Edit
sample groups** to assign a group to every sample. **Reset to condition**
restores the original `condition` values imported from the H5AD, even after
custom groups have been applied. Grouping changes are session-local and do not
modify the H5AD file.
