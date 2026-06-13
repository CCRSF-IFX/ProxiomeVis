# ProxiomeVis

ProxiomeVis is an interactive Shiny application for exploring Pixelator-derived
single-cell protein spatial data. It focuses on marker abundance, cell
annotation, differential readouts, self-clustering, and marker-pair
colocalization from Pixelator v4.1.1 Seurat objects.

## Features

- **QC**: cell filtering summaries, cell-calling rank plots, QC metric
  distributions, and original metadata inspection.
- **Abundance**: UMAP marker abundance views, marker distribution plots,
  cell-type composition, annotation heatmaps, and differential abundance.
- **Spatial Metrics**: clustering and colocalization readouts from stored
  proximity scores, including observed views, summary heatmaps, and
  differential analyses.
- **PixelatorES-style heatmaps**: condition, sample, and cell-type focused
  colocalization heatmaps with marker selection and legend controls.
- **Server-side RDS loading**: users on supported desktop or HPC runtimes can
  load an `.rds` file by path, with background progress reporting and
  processed app-data caching.

## Data Model

The app expects a Pixelator-compatible Seurat object with metadata, embeddings,
PNA assay abundance layers, and a stored assay `proximity` slot. RDS loading
reads proximity values from that stored slot. It does **not** rerun
`pixelatorR::ProximityScores()`.

The default demo data path used in the CCRSF deployment is:

```text
RnD_CS041188_BaoTran_XiaolinWu_3_Pixelgen_042126/notebooks/r/pg_data_combined_fil.pixelator_v4.1.1.rds
```

## Run The App

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

From the parent analysis repository:

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

Maintainers restore the project library during deployment or app updates:

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
