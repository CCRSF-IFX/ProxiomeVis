# ProxiomeVis

ProxiomeVis is an R Shiny app for viewing frozen Pixelgen proxiome results.
The expensive proximity analysis is completed before the app starts; the app
filters and plots compact tables stored in a prepared RDS file.

This mode is intentionally read-only. It does not recalculate proximity scores,
run tests on selected cells, or require the original PXL files at runtime.

## Requirements

- R with `SeuratObject`, `Matrix`, `pixelatorR`, `shiny`, `bslib`, `ggplot2`,
  `plotly`, `DT`, `jsonlite`, `data.table`, `ggrepel`, `scales`, `tibble`,
  `pheatmap`, and `yaml`
- a prepared RDS following the [precomputed-results contract](docs/precomputed-results.md)

`pixelatorR` must be installed even though the app does not query PXL files:
the saved Seurat object contains custom PNA classes defined by that package.

## Prepare the RDS

Build the runtime file from the annotated Seurat RDS and the matching analysis
output directory:

```bash
Rscript --vanilla scripts/prepare_precomputed_rds.R \
  ../pixelgen_demo_data/demo_analysis/cca_analysis/results/demo_cca_clustered.rds \
  ../pixelgen_demo_data/demo_analysis/tertiary_analysis \
  output/proxiomevis_demo_v1.rds
```

For this repository, the conventional output path is:

```text
output/proxiomevis_demo_v1.rds
```

`output/` is intentionally gitignored. Build the artifact on the deployment
host or copy it there through the project's controlled data-delivery process;
do not commit the demo RDS to this source repository.

The builder copies the frozen `abundance`, `self_proximity`, and
`colocalization` readouts into the versioned `object@misc$proxiomevis`
payload. See [Precomputed results](docs/precomputed-results.md) for the input
files, stored fields, and scientific limits.

The current builder targets this demo analysis layout and its `d1`/`d2`
outputs. The versioned runtime contract and Shiny reader remain separate from
that demo-specific preparation code.

`../pixelgen_demo_data` is a read-only source dataset. The builder reads the CCA
object and tertiary results there and writes only the requested `OUTPUT_RDS`.

## Run the app

By default, the app reads `output/proxiomevis_demo_v1.rds`. To use a
different trusted server-local file, set `PROXIOME_PRECOMPUTED_RDS`:

```bash
export PROXIOME_PRECOMPUTED_RDS=/absolute/path/to/proxiomevis_demo_v1.rds
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

The repository workspace also provides an isolated R executable. Limit BLAS
threads when using it on a shared host:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../.pixi/envs/r/bin/Rscript --vanilla \
  -e 'shiny::runApp(".", host="127.0.0.1", port=3838)'
```

Then open `http://127.0.0.1:3838` from the same machine, or use the secure
proxy provided by the deployment environment.

The donor and cell-type controls filter results that were already computed.
They do not recompute a multi-donor aggregate for an arbitrary subset. This
dataset contains two donors, so the combined donor-matched results are
descriptive and have no p-values. Donor-specific report volcanoes retain the
tertiary analysis's stored component-level p-values and adjusted p-values; they
are labeled exploratory and are not population-level donor inference. The app
never treats cells as independent biological replicates.

## Data handling

Use RDS paths that are trusted and already visible to the R process. `readRDS()`
can deserialize executable R objects, so do not expose unrestricted public RDS
uploads. For hosted use, place prepared files in an administrator-controlled
directory and pass the selected path through `PROXIOME_PRECOMPUTED_RDS`.

The in-app server-path field is disabled by default. A trusted local operator
can enable it before startup with `PROXIOME_ALLOW_LOCAL_LOAD=true`. Never enable
that flag on a public deployment: it lets a visitor request any RDS path that
the Shiny process can read.

Because the prepared RDS contains all runtime readouts, the server does not
need the large PXL files mounted in this mode.

## Tests

```bash
Rscript --vanilla tests/testthat.R
```

The unit-test fixture uses `SeuratObject` only and does not download or open a
real PXL file.

## Python implementation

The Python implementation remains preserved and maintained on the
`py_version` branch. The R implementation is developed independently; its
history does not reset or overwrite the Python branch.
