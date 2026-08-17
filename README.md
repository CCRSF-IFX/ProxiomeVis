# ProxiomeVis

ProxiomeVis is an interactive Shiny application for exploring Pixelator-derived
single-cell protein spatial data. It focuses on marker abundance, cell
annotation, differential readouts, self-clustering, and marker-pair
colocalization from Pixelator v4.1.1 Seurat objects.

The repository contains both the original R application (`app.R`) and a Python
Shiny port (`app.py`). The Python app accepts one processed `.h5ad` file and
does not load or analyze Pixelator `.pxl` files.

## Features

- **QC**: cell filtering summaries, cell-calling rank plots, QC metric
  distributions, and original metadata inspection.
- **Abundance**: UMAP marker abundance views, marker distribution plots,
  cell-type composition, annotation heatmaps, and differential abundance.
- **H5AD-only Python input**: load processed AnnData by server-visible path,
  with no PXL ingestion or proximity analysis.
- **Global analysis grouping**: choose a sample-level metadata column or assign
  custom sample groups for filters, summaries, heatmaps, and comparisons.

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
Alternatively, set the default before starting the app:

```bash
export PROXIOME_H5AD='/path/to/processed_data.h5ad'
uv run shiny run app.py
```

The default reference is:

```text
/Volumes/ccrsf-static/illumina/CCRSFIFX-23_MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/python_results/pg_data_combined_filtered_annotated.h5ad
```

The Python app reads stored observations, marker names, abundance layers, QC
history, cell annotations, and embeddings. It does not rerun Pixelator or
provide proximity, colocalization, patch-analysis, or 3D-layout views.

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
