# Data Input

## Python app: processed H5AD only

The Python Shiny app accepts one server-visible `.h5ad` file. It does not read
Pixelator `.pxl` files or run Pixelator analysis.

The reference dataset is:

```text
/Volumes/ccrsf-static/illumina/CCRSFIFX-23_MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/python_results/pg_data_combined_filtered_annotated.h5ad
```

Override the default path before launch when needed:

```bash
export PROXIOME_H5AD=/path/to/processed_data.h5ad
```

## Expected AnnData content

- `X`: component-by-marker counts
- `obs`: sample, condition, cell annotation, and QC metadata
- `var_names`: antibody marker names
- `layers["clr"]`: normalized abundance when available
- `obsm`: stored two-dimensional or higher embeddings such as UMAP, PCA, or Harmony
- `uns["qc_cell_counts_by_step"]`: optional notebook-generated QC retention history

Missing `sample_alias`, `condition`, or `celltype_manual` fields receive safe
fallback values. The app uses stored processed values and does not calculate
PXL proximity, colocalization, patch, or 3D layout data.

## Analysis grouping

Use **Data > Analysis grouping** after loading a dataset to choose a metadata
column that is constant within each sample, or assign custom sample groups.
Grouping changes are session-local and do not modify the H5AD file.
