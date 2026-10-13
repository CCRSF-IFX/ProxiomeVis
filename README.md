# ProxiomeVis

[Open the demo application](https://01a1327a-475a-181b-4bfa-2118f4c243fb.share.connect.posit.cloud/), hosted on Posit Connect Cloud.

ProxiomeVis is an R Shiny application for Pixelgen proxiome data.
The application shows precomputed results, calculates custom analyses, and shows 3D cell graphs.
The source analysis calculates proximity scores before application use.
Preparation scripts write scores and graph data to sidecar files.
The application does not open PXL files during use.

The application reads a small RDS file with a version number.
Sidecar files contain proximity scores and 3D graph data.
The application reads selected data from these files.
It does not load all sidecar data into each Shiny session.
Custom analyses do not replace precomputed results.

## Technical terms

| Term | Meaning |
| --- | --- |
| Application bundle | An RDS file that contains prepared data, precomputed results, and references to sidecar files. |
| Sidecar | A file outside the application bundle that supplies data to the application. |
| Precomputed results | Results that an analysis calculated before application use. The application does not change these results. |
| Custom analysis | An exploratory analysis that the user starts in the application. |
| Analysis worker | An isolated R process that calculates a custom analysis. |
| Provenance | A record of data sources, analysis settings, methods, and software. |

Package names, file formats, scientific terms, and interface labels keep their technical meanings.

## Requirements

Use R version 4.3.0 or later.
Install the packages for each necessary function.

| Function | R packages |
| --- | --- |
| Application | `Matrix`, `shiny`, `bslib`, `ggplot2`, `plotly`, `threejs`, `htmlwidgets`, `htmltools`, `DT`, `jsonlite`, `data.table`, `ggrepel`, `scales`, `tibble`, `pheatmap`, `yaml`, `callr`, `promises`, `later`, `digest` |
| Analysis worker | The application packages, plus `SeuratObject`, `Seurat`, `pixelatorR`, `arrow`, and `dplyr` |
| Data preparation | The analysis worker packages, plus `DBI` and `duckdb` |
| Cell 3D | The application packages, plus `arrow` and `dplyr` |

Use `threejs` version 0.3.4 or later from CRAN.
The [DESCRIPTION](DESCRIPTION) file gives the package requirements.

Seurat and pixelatorR are not necessary to load a plain application bundle.
For older PNA or Seurat RDS files, install the packages that define their object classes.
Use RDS files only from approved sources.

## Prepare the data

Run all commands from the `ProxiomeVis` directory.

### 1. Prepare the precomputed RDS file

The preparation script uses an annotated Seurat RDS file and the related analysis directory.
Run this command for the demo data:

```bash
Rscript --vanilla scripts/prepare_precomputed_rds.R \
  ../pixelgen_demo_data/demo_analysis/cca_analysis/results/demo_cca_clustered.rds \
  ../pixelgen_demo_data/demo_analysis/tertiary_analysis \
  output/proxiomevis_demo_v1.rds
```

The script copies the `abundance`, `self_proximity`, and `colocalization` results into `object@misc$proxiomevis`.
This structure includes a version number.
[Precomputed results](docs/precomputed-results.md) gives the input files, stored fields, and scientific limits.

The script uses the demo directory structure and its `d1` and `d2` results.
The application reader uses a versioned data format, not the demo directory structure.

Optional arguments select the UMAP, PCA, metadata columns, assay, and abundance layer.
For example, append `--umap umap --celltype-column celltype` when these names match your input.
The [preparation options](docs/precomputed-results.md#prepare-the-frozen-results) table gives all nine arguments and their defaults.
These arguments do not calculate new report results.

The script reads the CCA object and tertiary results from `../pixelgen_demo_data`.
It does not change this source dataset.
It writes only the specified `OUTPUT_RDS` file.

### 2. Prepare the application bundle

Run this command:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 \
  Rscript --vanilla scripts/prepare_hybrid_bundle.R \
  output/proxiomevis_demo_v1.rds output/hybrid_demo
```

The script writes these files under `output/hybrid_demo/`:

| Path | Contents |
| --- | --- |
| `analysis.rds` | The application bundle |
| `proximity/` | Proximity scores in partitioned Parquet files |
| `graphs/` | One Parquet file for each prepared sample, plus a version 3 `index.rds` file |

The graph files contain cell IDs, 3D coordinates, and marker identities.
They do not contain node IDs, edges, or marker counts.

An optional third argument specifies the source PXL directory.
An optional fourth argument specifies a text file with the selected cell IDs.
Each line must contain one cell ID that matches an ID in the dataset.
This selection changes only graph preparation.
It does not change the precomputed results or proximity data.

Alternatively, append `--graph-markers CD19,CD37,CD53,CD81,CD82` to select cells by their observed proximity scores.
The script ranks cells separately within each donor, cell type, and condition.
It selects the top 2 and bottom 2 cells for each different-marker pair among the specified markers.
It excludes self-marker pairs.
It exports the union of these cells without duplicates.
Each selected cell retains its complete 3D layout, including nodes for other markers.
Five markers give ten different-marker pairs and at most 40 selected cells per donor, cell type, and condition.
Do not combine this option with a cell-ID text file.

Use a new output directory if existing graphs contain cells outside the selection.
The script writes per-sample counts to `graph_export_summary.tsv` and reports them in the console.
It also saves the selected cell IDs and ranking settings.
See [ranked graph export](docs/precomputed-results.md#export-only-ranked-cells-for-3d) for details.

[Getting started](docs/getting-started.md#2-prepare-the-hybrid-application-bundle) gives the argument details.

### 3. Run preparation again or add cells

The preparation script validates completed sample files before it uses them again.
It keeps the cells that these files contain.
For more cells, the script writes a replacement sample file.
The replacement contains the existing cells and the new cells.
The script completes the replacement before it makes the file available.

Do not run two preparation processes in the same output directory at the same time.

The viewer can read version 2 graph sidecars.
It cannot read version 1 graph sidecars.
To replace version 1 or version 2 sidecars, prepare version 3 files in a new directory.
The preparation script does not automatically convert or delete the older files.

The precomputed RDS file alone is sufficient for the precomputed report views.

### 4. Put the prepared data on the server

Keep `analysis.rds`, `proximity/`, and `graphs/` together.
Git ignores the `output/` directory.
Do not commit the demo RDS file to this repository.

Use one of these methods:

- Prepare the files on the application server.
- Copy the files to the server through the project's controlled data-transfer process.

## Start the application

The default file is `output/hybrid_demo/analysis.rds`.
If this file is missing, the application uses `output/proxiomevis_demo_v1.rds`.

To use a different approved file on the server, set `PROXIOME_PRECOMPUTED_RDS`.
Start the application with this command:

```bash
export PROXIOME_PRECOMPUTED_RDS=/absolute/path/to/hybrid_demo/analysis.rds
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

This workspace also contains an isolated R executable.
The command sets the BLAS, OpenMP, and Arrow thread limits on a shared server.
Use this command if the workspace R executable is available:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 \
  ../.pixi/envs/r/bin/Rscript --vanilla \
  -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

On the same computer, open `http://127.0.0.1:3838`.
For remote access, use the secure proxy that your server environment supplies.

## Use the application

The navigation bar contains these tabs:

- **Overview / QC**
- **Embeddings**
- **Precomputed results**
- **Custom analysis**
- **Cell 3D**

The **Dataset** button is at the top right.
This button opens dataset information and, when enabled, controls to load an approved file from the server.
The sidecar status identifies the configured files.
The application validates each sidecar when it uses the file.

### Precomputed results

The **Tables**, **Summary**, **Volcano**, and **Spatial matrix** views use the same selections.
Under **Tables**, select a marker or pair to show its cell-level violin and boxplot by condition.
The plot uses normalized abundance or observed proximity scores for the selected donor and cell type.
Missing proximity scores are excluded; they are not replaced with zero.
The plot does not calculate a new statistical test.
These selections are the donor, cell type, comparison, and readout.
If the dataset has one donor, the application shows the donor name as a fixed label.
If the dataset has more than one donor, select one donor.

The results come from the selected donor.
They do not include results from other donors.
The summary can show all cell types.
One cell type is necessary for the volcano and spatial views.

A positive effect means that the target value is larger than the reference value.
The effect is the target value minus the reference value.
Stored component-level P values give exploratory evidence within one donor.
These P values do not give population-level evidence of treatment effects.
Cells are not independent biological replicates.

To select a feature, click a summary tile, volcano point, spatial pair, or result-table row.
The shared detail panel shows the selected feature.
The order of markers in a pair does not change the selection.

The volcano and spatial plots have hover and zoom functions.
The report SVG files keep the precomputed plot design and statistics.
The precomputed abundance heatmap is static.
The **Summary** view keeps this heatmap in a closed section.

Spatial condition profiles contain means.
Differential effects contain differences between medians.
These values describe different calculations.

False discovery rate (FDR) controls and table searches change only the tables and **Filtered table CSV** exports.
The **All stored rows CSV** export does not use these filters.
Heatmap controls and table filters do not change the precomputed report thresholds.
The application does not calculate hidden plots.
These views do not load proximity or graph sidecars.

### Custom analysis

1. Select one donor.
2. Select one cell type.
3. Select two conditions.
4. Select the markers.
5. If necessary, select cells with the embedding lasso, embedding box, or searchable cell selector.
6. Select the checkbox for the exploratory-analysis notice.
7. Start the analysis.

The analysis worker runs in an isolated R process.
You can use other tabs while the analysis runs.
You can cancel the analysis.

The results and JSON provenance contain the submitted cell groups, feature family, bundle checksum, methods, and software.
Different control settings do not change a completed result.
A new analysis uses the new settings.

### Cell 3D

1. Select a donor.
2. Select a cell type.
3. Select a condition.
4. Select the first ranking marker.
5. Select the second ranking marker.
6. Select a rank group.
7. Select a cell.
8. Select one highlight marker.
9. Click **Display**.

For a self-pair, select the same marker in both ranking fields.
The application ranks cells by observed log2 ratios.
The order of the two markers does not change the rank.
The rank groups contain the top 10 cells, bottom 10 cells, their union, or all cells with scores.
The union contains no duplicate cells.

Different control settings do not change the displayed graph until you click **Display**.
The graph contains all nodes for the displayed cell.
You can turn the graph, zoom, and use hover to read marker identities.
You can download the node data and provenance.
**Reset view** sets the camera to its initial view.
**Save PNG** downloads an image of the graph.

The coordinates give an abstract graph layout.
They do not show measured physical positions or distances.
Shared-node highlights are not a statistical colocalization test.
The viewer rejects cells with more than 200,000 nodes.
It does not remove nodes to stay below this limit.

The application shows how many cells have prepared graph data.
On first access, the application can read the full sample file to validate its checksum.
For subsequent cell selections, it uses the stored validation if the file is unchanged.
It reads only the matching Parquet row groups into R.

The selected-feature panel also contains **Explore cells in 3D**.
This button transfers the donor, cell type, and marker pair to **Cell 3D**.
It does not load a graph.
After this transfer, select a comparison condition.
Then select one highlight marker.
Click **Display** to load the graph.

### Download a rotating GIF

The GIF shows one full rotation of the displayed cell and highlight.
The export does not change the displayed camera view.

1. Open **Rotating GIF settings**.
2. Set the duration in seconds.
3. Set the frame rate in frames per second (FPS).
4. Set the image size in pixels.
5. Click **Download rotating GIF**.

| Setting | Default | Permitted values |
| --- | --- | --- |
| Duration | 6 seconds | Whole numbers from 1 to 20 |
| Frame rate | 15 FPS | Whole numbers from 1 to 30 |
| Width and height | 600 × 600 pixels | Whole numbers from 128 to 1024 for each side |

The width and height are equal.
Each export has a maximum of 300 frames and 120 million total pixels.
The encoded file limit is 64 MiB.
The application shows export progress.
To stop an export, click **Cancel GIF**.

The browser makes the GIF without PXL access or a new analysis.
The GIF plays continuously.
Each GIF frame has a maximum of 256 colors.
This color palette can change the appearance of transparent points.

[Hybrid analysis](docs/hybrid-analysis.md) gives the methods, data formats, and limits.

## Keep the data secure

### RDS files and server access

Use RDS files only from approved sources.
Make sure that the R process can read the selected file.
The `readRDS()` function can load R objects that contain executable code.
For public access, disable RDS uploads from unapproved sources.

For a hosted application, put the prepared files in a directory that an administrator controls.
Set `PROXIOME_PRECOMPUTED_RDS` to the selected file.

The **Dataset** panel disables server-path entry by default.
For trusted local use, set `PROXIOME_ALLOW_LOCAL_LOAD=true` before you start the application.

**SECURITY NOTICE:** Do not set `PROXIOME_ALLOW_LOCAL_LOAD=true` on a public server.
This setting lets a visitor request any RDS file that the Shiny process can read.

### Prepared data and analysis limits

Access to the source PXL files is not necessary for the application bundle and sidecars.
Keep the prepared files read-only.
Use atomic replacement when you change a dataset.
Each request uses the bundle that the interface loaded for that request.

Each Shiny R process can run one analysis worker at a time.
If a worker is active, the application rejects a second job with a retry message.

| Limit | Maximum for one job |
| --- | --- |
| Time | 120 seconds |
| Cells | 2,000 |
| Markers | 160 |
| Pairs | 1,000 |
| Collected score rows | 1 million |

These limits do not guarantee a maximum for total memory use.
Set CPU and memory limits in the server environment.
Each Shiny R process has its own worker limit.

Results stay in the user's session.
The application does not share a result cache between sessions.
It does not have a persistent job queue across servers.

## Run the tests

Install the packages for the analysis worker and data preparation before you run all tests.

Run the R tests:

```bash
Rscript --vanilla tests/testthat.R
```

The tests use synthetic Seurat inputs, small PXL-like databases, and separate R processes.
A download of the source demo data is not necessary.

For the browser tests, start the application with the prepared demo data.
Start Chromium with the Chrome DevTools Protocol (CDP) on port 9223.
Run these commands:

```bash
node tests/browser_smoke.mjs
node tests/browser_hybrid_smoke.mjs
```

WebGL is necessary for the 3D browser test.
Software rendering is sufficient.
The GIF test also uses the browser's image decoder.
