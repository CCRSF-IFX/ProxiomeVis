# Getting started

Run commands from the `ProxiomeVis` directory.

## 1. Prepare the demo RDS

The repository expects the demo source at `../pixelgen_demo_data`. Build the
runtime file from its CCA object and tertiary results:

```bash
Rscript --vanilla scripts/prepare_precomputed_rds.R \
  ../pixelgen_demo_data/demo_analysis/cca_analysis/results/demo_cca_clustered.rds \
  ../pixelgen_demo_data/demo_analysis/tertiary_analysis \
  output/proxiomevis_demo_v1.rds
```

The source tree is read-only. The script creates or replaces only the explicit
output path after validating the source object and frozen tables.

## 2. Start Shiny

The default input is `output/proxiomevis_demo_v1.rds`:

```bash
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

Open `http://127.0.0.1:3838`. To use another trusted prepared file:

```bash
export PROXIOME_PRECOMPUTED_RDS=/absolute/path/to/results.rds
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

`pixelatorR` must be installed when the real demo RDS contains its custom PNA
classes. No PXL mount is required after the runtime RDS has been built.

The app loads the configured RDS once and does not expose a visitor-selectable
server path. For trusted local administration only, set
`PROXIOME_ALLOW_LOCAL_LOAD=true` before startup to show the path field. Never
enable that setting on a public deployment.

## 3. Use the stored views

Donor, comparison, cell-type, and marker controls filter stored rows. They do
not trigger a new statistical test or rebuild an aggregate. If a combination
was not computed by the tertiary analysis, the app cannot create it at runtime.
The combined two-donor tables are descriptive and have no p-values.
Donor-specific volcanoes reproduce stored exploratory component-level tests;
they do not provide population-level donor inference.

## Tests

```bash
Rscript --vanilla tests/testthat.R
```

The test fixture is synthetic and uses a standard `SeuratObject` assay, so CI
does not need a real PXL file or the demo dataset.
