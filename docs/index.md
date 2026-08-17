# ProxiomeVis Python

ProxiomeVis Python is an interactive Shiny application for reviewing processed
Pixelgen AnnData. It accepts one `.h5ad` file and uses the stored counts,
normalized abundance, QC history, annotations, and embeddings.

## What the app shows

- **QC**: filtering summaries, molecule-rank plots, QC distributions, and metadata
- **Abundance**: marker abundance on embeddings, distributions, cell-type composition, annotation heatmaps, and differential abundance

The Python app does not read `.pxl` files and does not expose proximity,
clustering, colocalization, patch-analysis, or 3D-layout views.

## Typical workflow

1. Open the app.
2. Load a processed H5AD path from **Data**.
3. Review QC.
4. Explore abundance and annotation.
5. Optionally change the session-local analysis grouping.
6. Download figures or tables as needed.
