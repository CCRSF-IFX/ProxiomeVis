# Getting Started

## Launch the Python app

From this directory:

```bash
uv sync --extra test
uv run shiny run --reload app.py
```

## Load data

Open **Data**, enter the full path to one readable processed `.h5ad` file, and
assign the matching PXL paths before clicking **Load Data**. The reference H5AD
is prefilled unless `PROXIOME_H5AD` overrides it.

## Navigate the app

The Python app contains four analysis tabs:

- **QC** for filtering history, molecule-rank curves, distributions, and metadata
- **Abundance** for embeddings, marker distributions, annotations, and differential abundance
- **Spatial Metrics** for PXL-backed self-clustering, colocalization,
  differential spatial scores, and 3D cellgraphs
- **Patch Analysis** for stored patch marker, Raji signal, and burden tables

Spatial plots require assigned PXL files with matching H5AD component IDs.
Patch plots require their optional precomputed H5AD payload.

## Recommended first checks

1. Open **QC > Filtering** and confirm the expected samples and retained counts.
2. Open **QC > Cell Calling** and inspect the stored `n_umi` values.
3. Open **Abundance > Observed** and inspect marker signal on a stored embedding.
4. Open **Abundance > Cell Annotation** to verify the cell populations.
