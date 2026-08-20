# ProxiomeVis

ProxiomeVis is an interactive Shiny application for exploring Pixelator-derived
single-cell protein spatial data. It focuses on marker abundance, cell
annotation, differential readouts, self-clustering, and marker-pair
colocalization from Pixelator v4.1.1 Seurat objects.

The repository contains both the original R application (`app.R`) and a Python
Shiny port (`app.py`). The Python app accepts one processed `.h5ad` file and
one or more matching `.pxl` files. Cell-level data come from H5AD; proximity
scores and component layouts are queried from PXL on demand.

## Features

- **QC**: cell filtering summaries, cell-calling rank plots, QC metric
  distributions, and original metadata inspection.
- **Abundance**: UMAP marker abundance views, marker distribution plots,
  cell-type composition, annotation heatmaps, and differential abundance.
- **Spatial Metrics**: one shared retrieval feeds proximity profiling,
  clustering, differential colocalization, and 3D layouts.
- **Patch Analysis**: marker unmixing, Raji signal, and patch-burden views when
  their precomputed tables are stored in the H5AD.
- **On-demand PXL queries**: read only the selected components and markers for
  spatial tables and stored 3D layouts.
- **Activity Log**: watch severity-tagged H5AD loading, spatial retrieval, PXL
  queries, and 3D layout reads for the current session, then download a
  sanitized diagnostics bundle for support.
- **Global analysis grouping**: choose a sample-level metadata column or assign
  custom sample groups for filters, summaries, heatmaps, and comparisons, with
  one-click reset to the H5AD's original `condition` values.

## Data Model

The Python app expects processed AnnData with `obs`, marker names, abundance in
`X` or `layers["clr"]`, and stored embeddings in `obsm`.

The original R app expects a Pixelator-compatible Seurat object with metadata, embeddings,
PNA assay abundance layers, and a stored assay `proximity` slot. RDS loading
reads proximity values from that stored slot. It does **not** rerun
`pixelatorR::ProximityScores()`.

The default demo data path used in the CCRSF deployment is:

```text
RnD_CS041188_BaoTran_XiaolinWu_3_Pixelgen_042126/notebooks/r/pg_data_combined_fil.pixelator_v4.1.1.rds
```

## Run The App

### Python

The Python app supports Python 3.10 and 3.11. From this directory, create the
environment and start the Python Shiny app with:

```bash
uv sync --extra test
uv run shiny run --reload app.py
```

Open the app's **Data** panel and enter the path to one processed `.h5ad` file.
Assign the matching `.pxl` files to enable spatial metrics and 3D cellgraphs.
Alternatively, set the default before starting the app:

```bash
export PROXIOME_H5AD='/path/to/processed_data.h5ad'
export PROXIOME_PXL='/path/to/layouts/*.layout.pxl'
export PROXIOMEVIS_HOME="$HOME/.ProxiomeVis"
uv run shiny run app.py
```

The Python server writes one structured JSONL log per browser session under
`$PROXIOMEVIS_HOME/runtime` (default: `$HOME/.ProxiomeVis/runtime`). Each event
includes the session ID, severity, app version, Git commit, and traceback for
failures. Set `PROXIOMEVIS_VERSION` and `PROXIOMEVIS_COMMIT` in packaged
deployments where Git metadata is unavailable. Browser errors and downloaded
diagnostics redact full server paths.

The default reference is:

```text
/Volumes/ccrsf-static/illumina/CCRSFIFX-23_MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/python_results/pg_data_combined_filtered_annotated.h5ad
```

The matching default PXL glob is:

```text
/Volumes/ccrsf-static/singlecell_projects/MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/Analysis_2nd_combo/Analysis/pixelator/*.pxl
```

The Python app reads observations, marker names, abundance layers, QC history,
cell annotations, embeddings, and optional patch tables from H5AD. It ignores
any embedded H5AD proximity table and reads precomputed spatial scores from the
assigned PXL files. PXL proximity queries are filtered and aggregated before
their results are returned to the app.

Before opening a spatial visualization, use **Spatial Metrics > Retrieve Data**
to choose analysis groups, samples, cell types, and markers. **All markers** is
the default. Clicking **Retrieve Spatial Data** replaces the active retrieval;
editing the controls alone does not change it. Every downstream spatial view
can narrow, but never expand, that active population and marker scope.

### R

From the parent analysis repository:

```bash
pixi run -e r serve-shiny-proxiome
```

From this Shiny app directory with a restored R environment:

```bash
Rscript -e "shiny::runApp('.')"
```

The first load of a new RDS can take several minutes because the app builds
compact tables for interactive use. Processed app data are cached under:

```text
$HOME/.ProxiomeVis/cache
```

If a bundled `cache/demo_proxiome_data.rds` exists, the app uses it for the
demo dataset. Otherwise it writes a user-local cache under `$HOME/.ProxiomeVis`.

## Tests

Run the Python tests with:

```bash
uv run --extra test pytest -q
```

Run the R tests from the parent analysis repository:

```bash
pixi run -e r test-shiny-proxiome
```

From this Shiny app directory with dependencies restored:

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Open OnDemand Deployment

The Open OnDemand deployment uses this app folder as a shared `renv` project.
The app directory can be copied to a shared location and restored there; it is
not copied into user home directories at launch time.

Common deployment overrides:

```bash
export PROXIOME_APP_DIR=/path/to/shared/proxiome_demo
export PROXIOME_R_MODULE=R/4.5.2
export PROXIOME_DEMO_RDS=/path/to/shared/data/demo.rds
export PROXIOMEVIS_HOME=$HOME/.ProxiomeVis
```

Runtime cache and diagnostics use `$HOME/.ProxiomeVis`, not the shared
application directory. Browser file upload is disabled because of Open
OnDemand proxy limits. Users can load their own data by entering an `.rds`
path that is visible on the HPC or desktop filesystem.

Git tracks `renv/activate.R`, `renv/settings.json`, and `renv.lock`; generated
package libraries under `renv/library/` stay out of Git. Deployers restore the
project library during deployment or app updates:

```bash
cd /path/to/shared/proxiome_demo
Rscript -e 'renv::restore(prompt = FALSE)'
```

On Biowulf, run the restore with the same R module used by Open OnDemand:

```bash
module load R/4.5.2
cd /path/to/shared/proxiome_demo
Rscript -e 'renv::restore(prompt = FALSE)'
```

The app startup only activates the restored project library. It does not run
`renv::restore()` or install packages for end users.
