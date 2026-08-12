# Data Input

## RDS requirements

ProxiomeVis expects a Pixelator-compatible Seurat object saved as an RDS file.
The object should include:

- cell metadata
- embeddings such as UMAP
- PNA assay abundance data
- stored Pixelator proximity outputs for clustering and colocalization

RDS loading reads stored proximity values. It does **not** call
`pixelatorR::ProximityScores()` again.

## Demo data

The CCRSF deployment is configured around a Pixelator v4.1.1 demo RDS:

```text
RnD_CS041188_BaoTran_XiaolinWu_3_Pixelgen_042126/notebooks/r/pg_data_combined_fil.pixelator_v4.1.1.rds
```

Deployments can override the demo RDS with:

```bash
export PROXIOME_DEMO_RDS=/path/to/demo.rds
```

## 3D layout files

The 3D layout view uses Pixelator `.layout.pxl` files. These files are separate
from the RDS and are usually found under a Pixelator results directory, for
example:

```text
results/run_pixelator-4.1.1_merged_pixelator_v0.27.2/pixelator/3_CD3CD28.layout.pxl
```

The app searches for `<sample>.layout.pxl` near the loaded RDS. If files are in
a separate directory, set:

```bash
export PROXIOME_LAYOUT_DIR=/path/to/pixelator/layout/files
```

If a matching `.layout.pxl` is not found, the 3D Layout tab shows a message
instead of a plot.

## User cache

Processed app data are cached under:

```text
$HOME/.ProxiomeVis/cache
```

The cache avoids repeated conversion work for the same RDS. It does not replace
the source RDS.

## Analysis grouping

Use **Data > Analysis grouping…** after loading a dataset to choose the metadata
column used as the global analysis group. Only columns with one value per
sample are offered. Choose **Edit sample groups** to assign a custom group to
each sample instead.

The selection applies to abundance, clustering, colocalization, QC labels,
filters, heatmaps, and differential comparisons. Changes are session-local:
the source RDS and processed cache are not modified. Loading another dataset
resets grouping to its `condition` column when available.
