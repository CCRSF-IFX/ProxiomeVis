# Spatial Metrics

Spatial Metrics contains two readouts:

- **Clustering**: marker self-proximity within each cell graph.
- **Colocalization**: marker-pair proximity between two different proteins.

These readouts use stored Pixelator proximity values from the loaded data. The
app does not recompute proximity scores during interactive use.

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
two groups.

Use **Colocalization > 3D Layout** to inspect one selected Pixelator cell graph
in 3D. The app reads stored `wpmds_3d` coordinates from `.layout.pxl`, labels
nodes with marker names from the edge list, and renders an interactive Plotly
scatter plot.

Controls for the 3D layout:

- **Sample**: chooses the `.layout.pxl` file.
- **Cell type**: filters available cell/component choices.
- **Cell/component**: chooses the graph component to render.
- **Max background nodes**: caps non-highlighted nodes for browser performance.
- **Highlighted markers**: markers drawn larger and in color.
- **Options**: plot canvas width and height.

If the plot is slow, reduce **Max background nodes** or select a smaller
component.
