# ProxiomeVis

ProxiomeVis is a Python Shiny application for exploring processed Pixelgen
single-cell protein spatial data. The retired R implementation is preserved on
the `legacy/r-version` branch.

The app accepts one processed `.h5ad` file and matching `.pxl` files. Cell
metadata, abundance, annotations, embeddings, and optional patch results come
from H5AD. Proximity scores and component layouts are queried from PXL on
demand.

## Features

- **QC**: stored filtering history, cell-calling rank plots, QC distributions,
  and metadata inspection.
- **Abundance**: embedding views, marker distributions, cell-type composition,
  annotation heatmaps, and differential abundance.
- **Spatial Metrics**: one prepared retrieval feeds proximity profiling,
  clustering, differential colocalization, and 3D layouts.
- **Patch Analysis**: receiver/target marker selection, prepared candidate
  screening from PXL, and stored patch burden/composition views.
- **Documentation**: an in-app production workflow, parameter reference,
  interpretation guidance, and reporting checklist shared with the MkDocs site.
- **Activity Log**: structured session events and downloadable sanitized
  diagnostics.
- **Analysis grouping**: metadata-derived or custom sample groups shared across
  the application.

## Requirements

- Python 3.10 or 3.11
- One processed AnnData `.h5ad` file
- Matching `.pxl` files for spatial metrics, patch candidate screening, and 3D
  cellgraphs

## Run the app

From this directory:

```bash
uv sync --extra test
uv run shiny run --reload app.py
```

Open **Data** and enter the server-visible H5AD and PXL paths. Defaults can be
set before launch:

```bash
export PROXIOME_H5AD='/path/to/processed_data.h5ad'
export PROXIOME_PXL='/path/to/layouts/*.layout.pxl'
export PROXIOMEVIS_HOME="$HOME/.ProxiomeVis"
uv run shiny run app.py
```

The reference H5AD is:

```text
/Volumes/ccrsf-static/illumina/CCRSFIFX-23_MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/python_results/pg_data_combined_filtered_annotated.h5ad
```

The matching default PXL glob is:

```text
/Volumes/ccrsf-static/singlecell_projects/MarinaDobrovolskaia_CS041374_6_Pixelgen_062226/Analysis_2nd_combo/Analysis/pixelator/*.pxl
```

## Data and retrieval model

H5AD supplies observations, markers, abundance, QC history, annotations,
embeddings, and optional patch tables. Embedded H5AD proximity tables are
ignored; the app filters and aggregates precomputed PXL proximity data before
returning results.

Use **Spatial Metrics > Retrieve Data** to prepare a population and marker
scope. Use **Prepare Patch Data** to freeze receiver, target, analysis-group,
and marker selections for Patch Analysis. Editing retrieval controls does not
replace active results until the corresponding preparation button is clicked.

## Tests

```bash
uv run --extra test pytest -q
```

## Diagnostics

The server writes one structured JSONL log per browser session under
`$PROXIOMEVIS_HOME/runtime`. Events include session ID, severity, app version,
Git commit, and server traceback. Browser messages and downloaded diagnostics
redact full server paths.

Set `PROXIOMEVIS_VERSION` and `PROXIOMEVIS_COMMIT` when Git metadata is not
available in a packaged deployment.

## Open OnDemand deployment

The launcher in `template/script.sh.erb` expects a deployed virtual environment
at `${PROXIOME_APP_DIR}/.venv` by default. Restore it during deployment:

```bash
cd /path/to/shared/proxiome_demo
uv sync --frozen
```

Common overrides are:

```bash
export PROXIOME_APP_DIR=/path/to/shared/proxiome_demo
export PROXIOME_PYTHON=/path/to/python
export PROXIOME_H5AD=/path/to/processed_data.h5ad
export PROXIOME_PXL='/path/to/layouts/*.layout.pxl'
export PROXIOMEVIS_HOME=$HOME/.ProxiomeVis
```

Runtime cache and diagnostics use `$PROXIOMEVIS_HOME`, not the shared
application directory. The launcher validates the Python interpreter and does
not install dependencies at user startup.
