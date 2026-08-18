# ProxiomeVis Python

ProxiomeVis Python is an interactive Shiny application for reviewing processed
Pixelgen AnnData. It accepts one `.h5ad` file and uses the stored counts,
normalized abundance, QC history, annotations, and embeddings.

## What the app shows

- **QC**: filtering summaries, molecule-rank plots, QC distributions, and metadata
- **Abundance**: marker abundance on embeddings, distributions, cell-type composition, annotation heatmaps, and differential abundance
- **Spatial Metrics**: H5AD-stored self-clustering and marker-pair proximity,
  plus optional PXL-backed 3D cellgraphs
- **Patch Analysis**: H5AD-stored patch marker, Raji signal, and burden tables

The Python app uses H5AD for every analysis table. It reads assigned PXL files
only when displaying a selected component's precomputed 3D layout.

## Typical workflow

1. Open the app.
2. Load a processed H5AD path from **Data**.
3. Review QC.
4. Explore abundance and annotation.
5. Review spatial or patch modules when the H5AD contains those tables.
6. Optionally assign PXL paths and display a selected cellgraph.
7. Optionally change the session-local analysis grouping.
8. Download figures or tables as needed.
