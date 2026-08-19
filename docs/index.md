# ProxiomeVis Python

ProxiomeVis Python is an interactive Shiny application for reviewing processed
Pixelgen AnnData. It accepts one `.h5ad` file and uses the stored counts,
normalized abundance, QC history, annotations, and embeddings.

## What the app shows

- **QC**: filtering summaries, molecule-rank plots, QC distributions, and metadata
- **Abundance**: marker abundance on embeddings, distributions, cell-type composition, annotation heatmaps, and differential abundance
- **Spatial Metrics**: PXL-backed self-clustering, marker-pair proximity, and
  3D cellgraphs
- **Patch Analysis**: H5AD-stored patch marker, Raji signal, and burden tables

The Python app uses H5AD for cell metadata, abundance, embeddings, and optional
patch tables. It queries proximity and layouts from assigned PXL files.

## Typical workflow

1. Open the app.
2. Load a processed H5AD path from **Data**.
3. Review QC.
4. Explore abundance and annotation.
5. Assign matching PXL paths and review spatial metrics or a cellgraph.
6. Review patch modules when the H5AD contains those tables.
7. Optionally change the session-local analysis grouping.
8. Download figures or tables as needed.
