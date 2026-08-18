# Getting Started

## Launch the Python app

From this directory:

```bash
uv sync --extra test
uv run shiny run --reload app.py
```

## Load data

Open **Data**, enter the full path to one readable processed `.h5ad` file, and
click **Load Data**. Optionally assign PXL paths for 3D cellgraphs. The
reference dataset is prefilled unless `PROXIOME_H5AD` overrides it.

## Navigate the app

The Python app contains four analysis tabs:

- **QC** for filtering history, molecule-rank curves, distributions, and metadata
- **Abundance** for embeddings, marker distributions, annotations, and differential abundance
- **Spatial Metrics** for stored self-clustering, colocalization, differential
  spatial scores, and optional PXL-backed 3D cellgraphs
- **Patch Analysis** for stored patch marker, Raji signal, and burden tables

Spatial and patch plots require their precomputed H5AD payloads. The bundled
reference H5AD does not currently include those payloads.

## Recommended first checks

1. Open **QC > Filtering** and confirm the expected samples and retained counts.
2. Open **QC > Cell Calling** and inspect the stored `n_umi` values.
3. Open **Abundance > Observed** and inspect marker signal on a stored embedding.
4. Open **Abundance > Cell Annotation** to verify the cell populations.
