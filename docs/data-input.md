# Data Input

## Python app: H5AD and PXL data

The Python Shiny app accepts one server-visible `.h5ad` file for cell data and
matching `.pxl` files for proximity and 3D component layouts. It reads the
precomputed PXL proximity table; it does not recalculate proximity from edges.

The reference dataset is:

```text
/Volumes/ccrsf-static/illumina/CCRSFIFX-23_MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/python_results/pg_data_combined_filtered_annotated.h5ad
```

Its matching PXL files are prefilled with:

```text
/Volumes/ccrsf-static/singlecell_projects/MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/Analysis_2nd_combo/Analysis/pixelator/*.pxl
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
- `uns["proxiome"]["patch"]`: optional mapping containing `run_plan`,
  `marker_unmixing`, `raji_marker_abundance`, `raji_marker_proximity`, and
  `patch_burden` tables
- `uns["proxiome"]["component_layouts"]`: optional mapping from component ID
  to a node table containing `x`, `y`, and `z`

Missing `sample_alias`, `condition`, or `celltype_manual` fields receive safe
fallback values. An embedded `uns["proxiome"]["proximity"]` table is ignored.
Missing optional patch/layout tables leave their corresponding views visible
with an unavailable-data message.

## PXL paths

In **Data**, assign a `.layout.pxl` file, directory, glob, or comma/newline-
separated paths. Clustering and colocalization query the PXL proximity table by
the H5AD component IDs and selected markers. Differential colocalization is
aggregated inside the PXL DuckDB database before results are returned. The 3D
view reads the selected component's stored layout on demand. PXL is not used
for QC, abundance, annotations, or patch analysis.

## Analysis grouping

Use **Data > Analysis grouping** after loading a dataset. Select **Use metadata
column** to group by a column that is constant within each sample, or **Edit
sample groups** to assign a group to every sample. **Reset to condition**
restores the original `condition` values imported from the H5AD, even after
custom groups have been applied. Grouping changes are session-local and do not
modify the H5AD file.
