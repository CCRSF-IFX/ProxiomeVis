# Getting Started

## Launch the app

From the parent analysis repository, launch ProxiomeVis with:

```bash
pixi run -e r serve-shiny-proxiome
```

From this Shiny app directory, if the R environment is already restored:

```bash
Rscript -e "shiny::runApp('.')"
```

The app opens in a browser. On Open OnDemand, use the URL provided by the
interactive session.

## Load data

The app starts with the configured demo dataset when available. To load another
dataset, use the **Data source** controls and enter the full path to a readable
`.rds` file.

The first load of a new RDS can take several minutes. ProxiomeVis builds compact
tables for interactive plotting and stores them in a user-local cache. Loading
the same RDS again should be faster if the cache is still valid.

## Navigate the app

The main tabs are:

- **QC**
- **Abundance**
- **Spatial Metrics**

Each tab has a left sidebar for filters and data-changing controls. Plot size
and download controls are placed next to plots in the **Options** popover.

## Recommended first checks

1. Open **QC > Filtering** and confirm the expected samples are present.
2. Open **QC > Cell Calling** and check cell-rank curves.
3. Open **Abundance > Observed** and inspect marker signal on UMAP.
4. Open **Spatial Metrics > Colocalization > Observed** for marker-pair spatial
   patterns.
