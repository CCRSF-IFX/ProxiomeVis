# Spatial Metrics

Spatial Metrics contains a shared retrieval step and four downstream views:

- **Retrieve Data**: defines the cell population and marker scope.
- **Proximity Profile**: profiles marker-pair proximity scores and pair details.
- **Clustering**: marker self-proximity within each cell graph.
- **Colocalization**: differential marker-pair tests.
- **3D Layout**: interactive views of individual Pixelator cell graphs.

These readouts query the precomputed `proximity` table in the assigned PXL
files. The app uses H5AD component IDs to push the active cell and marker
filters into PXL and does not recalculate proximity from graph edges. An
embedded H5AD proximity table is ignored.

## Retrieve Data

Choose the analysis groups, samples, cell types, and number of markers, then
click **Retrieve Spatial Data**. Marker selection defaults to **All markers**;
you can instead retrieve the top abundance-ranked markers or a manual set.

The completed retrieval becomes the active spatial scope. Clustering,
colocalization, differential analyses, pair details, and 3D component choices
are all restricted to its cells and markers. Controls inside those views can
only narrow the active scope. If you edit retrieval settings, the existing
visualizations continue to use the prior retrieval until you click **Retrieve
Spatial Data** again. Loading another H5AD or changing analysis grouping clears
the active retrieval.

The app freezes component IDs and marker choices rather than copying the full
PXL proximity table into browser memory. Each visualization pushes its smaller
query into the PXL-backed DuckDB view. This preserves a consistent retrieval
boundary without materializing every marker pair.

## Proximity Profile

Use **Proximity Profile** immediately after retrieval to view marker-pair
proximity heatmaps. The heatmap can summarize by analysis group or sample and
can focus on one group or compare groups.

The default retrieval includes all markers, while the default heatmap uses the
PixelatorES proximity profile for readability and query performance. A manual
heatmap marker set must be a subset of the retrieved markers. The applied
population and proximity definition are reused by differential colocalization.

## Clustering

Use **Clustering > Per Marker** to compare one marker across conditions or cell
types.

Use **Clustering > Summary Heatmap** for a marker-level heatmap across selected
conditions and cell types. Keep **Protein set** at **Top variable proteins** to
rank proteins by variation in mean self-proximity, or choose **Custom proteins**
to display an explicit ordered list from the active spatial retrieval.

Use **Clustering > Differential** to compare marker self-proximity between two
groups.

Clustering controls update their plots reactively inside the active retrieval;
there is no separate run button.

## Colocalization

Use **Colocalization** to compare marker-pair proximity between
two groups within one cell population. The analysis first calculates the
selected population or detected-cell mean for every sample, then compares the
median sample values between groups. Samples—not cells—are the statistical
replicates. If either group has fewer than the requested number of samples, the
app reports descriptive effects but leaves p-values and FDR unavailable.

The volcano shows every marker pair by default. Use **Pairs shown** to focus on
pairs containing one marker, and click a volcano point to update the sample-level
detail plot. Missing pairs are included as zero only for the population mean.

Use **3D Layout** to inspect one selected Pixelator cell graph
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
