# Getting Started

## Launch the Python app

From this directory:

```bash
uv sync --extra test
uv run shiny run --reload app.py
```

## Load data

Open **Data**, enter the full path to one readable processed `.h5ad` file, and
click **Load Data**. The reference dataset is prefilled unless
`PROXIOME_H5AD` overrides it.

## Navigate the app

The Python app contains two analysis tabs:

- **QC** for filtering history, molecule-rank curves, distributions, and metadata
- **Abundance** for embeddings, marker distributions, annotations, and differential abundance

It intentionally has no PXL, proximity, colocalization, patch-analysis, or 3D
layout workflow.

## Recommended first checks

1. Open **QC > Filtering** and confirm the expected samples and retained counts.
2. Open **QC > Cell Calling** and inspect the stored `n_umi` values.
3. Open **Abundance > Observed** and inspect marker signal on a stored embedding.
4. Open **Abundance > Cell Annotation** to verify the cell populations.
