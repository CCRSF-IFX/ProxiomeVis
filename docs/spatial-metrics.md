# Spatial Metrics

Spatial Metrics contains two readouts:

- **Clustering**: marker self-proximity within each cell graph.
- **Colocalization**: marker-pair proximity between two different proteins.

These readouts query the precomputed `proximity` table in the assigned PXL
files. The app uses H5AD component IDs to push the active cell and marker
filters into PXL and does not recalculate proximity from graph edges. An
embedded H5AD proximity table is ignored.

## Clustering

Use **Clustering > Observed** to inspect a selected marker's self-proximity
across cells.

Use **Clustering > Per Marker** to compare one marker across conditions or cell
types.

Use **Clustering > Summary Heatmap** for a marker-level heatmap across selected
conditions and cell types.

Use **Clustering > Differential** to compare marker self-proximity between two
groups.

## Colocalization

Use **Colocalization > Observed** to view marker-pair proximity heatmaps. The
heatmap can summarize by condition, sample, or focused cell type.

Use **Colocalization > Differential** to compare marker-pair proximity between
two groups within one cell population. The analysis first calculates the
selected population or detected-cell mean for every sample, then compares the
median sample values between groups. Samples—not cells—are the statistical
replicates. If either group has fewer than the requested number of samples, the
app reports descriptive effects but leaves p-values and FDR unavailable.

The volcano shows every marker pair by default. Use **Pairs shown** to focus on
pairs containing one marker, and click a volcano point to update the sample-level
detail plot. Missing pairs are included as zero only for the population mean.

Use **Colocalization > 3D Layout** to inspect one selected Pixelator cell graph
in 3D. Assign one or more `.layout.pxl` paths in **Data**. The app reads that
component's stored layout on demand and renders an interactive Plotly scatter
plot. The same assigned PXL files provide the spatial analysis tables.

Controls for the 3D layout:

- **Sample**: chooses the `.layout.pxl` file.
- **Cell type**: filters available cell/component choices.
- **Cell/component**: chooses the graph component to render.
- **Max background nodes**: caps non-highlighted nodes for browser performance.
- **Highlighted markers**: markers drawn larger and in color.
- **Options**: plot canvas width and height.

If the plot is slow, reduce **Max background nodes** or select a smaller
component.
