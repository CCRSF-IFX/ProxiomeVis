# Pixelgen Proxiome Shiny Demo

This Shiny app uses the Pixelator v4.1.1 Seurat object at:

`RnD_CS041188_BaoTran_XiaolinWu_3_Pixelgen_042126/notebooks/r/pg_data_combined_fil.pixelator_v4.1.1.rds`

The app organizes the demo around Abundance and Spatial Metrics:

- `Abundance`: per-cell PNA marker abundance on the Seurat UMAP embeddings.
- `Spatial Metrics`: proximity readouts, with `Clustering` for marker self-proximity (`marker_1 == marker_2`) and `Colocalization` for marker-pair proximity (`marker_1 != marker_2`).

RDS loading reads proximity values from the stored assay proximity slot. It does
not rerun `pixelatorR::ProximityScores()`.

The Spatial Metrics > Colocalization observed view uses a PixelatorES-style heatmap strategy:

- `Condition summary` keeps the original condition-level colocalization heatmap.
- `Sample summary` groups mean marker-pair `log2_ratio` by sample alias and condition.
- `Cell type focus` filters to one cell type, then groups sample-level marker-pair summaries.
- `Variable detected markers` ranks markers by variable mean `log2_ratio` among detectable pairs, capped at 40 markers for readable heatmaps.

Run from the repository root:

```bash
pixi run -e r serve-shiny-proxiome
```

The first launch loads the full RDS and writes a compact cache. If a bundled
`cache/demo_proxiome_data.rds` exists, the app uses it. Otherwise it writes the
cache under the user's hidden writable directory at
`$HOME/.ProxiomeVis/cache/demo_proxiome_data.rds`.

Run the app tests:

```bash
pixi run -e r test-shiny-proxiome
```

## Shared Open OnDemand Deployment

The Open OnDemand deployment uses this app folder as a shared renv project. The
`proxiome_demo` directory can be copied to a shared app location and restored
there; it is not copied into user home directories at launch time.

The launcher keeps the current CCRSF path as the default, but supports
site-specific overrides for other HPC systems such as Biowulf:

```bash
export PROXIOME_APP_DIR=/path/to/shared/proxiome_demo
export PROXIOME_R_MODULE=R/4.5.2
export PROXIOME_DEMO_RDS=/path/to/shared/data/demo.rds
export PROXIOMEVIS_HOME=$HOME/.ProxiomeVis
```

Runtime cache and diagnostics use `$HOME/.ProxiomeVis`, not the shared
application directory. Browser file upload is disabled because of Open OnDemand
proxy limits. Users can load their own data by entering an `.rds` path that is
already visible on the HPC filesystem or local desktop filesystem.

Developers update dependencies in the shared project lockfile when package
requirements change:

```r
renv::snapshot()
```

Maintainers restore the project library during deployment or application
updates:

```bash
cd /mnt/ccrsf-static/illumina/RnD_pixelgen_CAR-T_datasets
cd shiny/proxiome_demo
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
# ProxiomeVis
