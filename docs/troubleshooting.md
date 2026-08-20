# Troubleshooting

## The H5AD load takes a long time

Large H5AD files can take time to materialize metadata and abundance tables.
Use the Activity Log to distinguish H5AD loading from later PXL queries.

## The app says the H5AD path is not readable

Check that the path is absolute and visible from the machine running Shiny. On
Open OnDemand, the path must be readable from the compute session, not only from
your local computer. The input must be a processed `.h5ad` file.

## A plot is empty

Check the sidebar filters. Empty plots usually mean the selected condition,
cell type, marker, or contrast has no matching cells.

For differential views, verify that group A and group B both contain enough
cells after filtering.

## The 3D Layout tab cannot find a layout file

In the Python app, assign a matching `.layout.pxl` file, directory, or glob in
the **Data** menu, or set:

```bash
export PROXIOME_PXL='/path/to/pixelator/layout/files/*.layout.pxl'
```

Then reload the H5AD. The PXL filename must contain the selected sample name
unless only one PXL file is assigned.

## The 3D Layout plot is slow

Large Pixelator graph components can contain many nodes. Reduce **Max background
nodes** or choose a smaller component.

Highlighted marker nodes are kept; the cap mainly reduces non-highlighted
background nodes.

## SVG download fails

Open **Activity Log** and check the latest operation. SVG export uses Plotly's
static-image engine and requires the application's Kaleido runtime. If the
problem persists, download the diagnostics bundle and provide it to the
application administrator.

## An operation reports a reference ID

Open **Activity Log** and click **Download diagnostics**. The ZIP contains
sanitized runtime metadata and structured events for the current browser
session, including app version, Git commit, severity, and error traceback. It
does not include H5AD/PXL contents or full server paths. Server-side JSONL logs
are stored under `$PROXIOMEVIS_HOME/runtime` (default:
`$HOME/.ProxiomeVis/runtime`) for administrators.

## Build the documentation locally

From the app directory:

```bash
mkdocs build --strict
```

To preview locally:

```bash
mkdocs serve
```
