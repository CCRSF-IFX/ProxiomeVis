# Troubleshooting

## The RDS load takes a long time

The first load of a new RDS can take several minutes because ProxiomeVis builds
compact tables for interactive plotting. Loading the same RDS again should be
faster if the cache under `$HOME/.ProxiomeVis/cache` is valid.

## The app says the RDS path is not readable

Check that the path is absolute and visible from the machine running Shiny. On
Open OnDemand, the path must be readable from the compute session, not only from
your local computer.

## A plot is empty

Check the sidebar filters. Empty plots usually mean the selected condition,
cell type, marker, or contrast has no matching cells.

For differential views, verify that group A and group B both contain enough
cells after filtering.

## The 3D Layout tab cannot find a layout file

The app needs a matching `<sample>.layout.pxl` file. Either place the Pixelator
layout files near the loaded RDS in the expected results directory, or set:

```bash
export PROXIOME_LAYOUT_DIR=/path/to/pixelator/layout/files
```

Then restart the app.

## The 3D Layout plot is slow

Large Pixelator graph components can contain many nodes. Reduce **Max background
nodes** or choose a smaller component.

Highlighted marker nodes are kept; the cap mainly reduces non-highlighted
background nodes.

## PNG or SVG download fails

Try a smaller plot width and height first. If SVG fails for a complex plot, use
PNG. The 3D Layout view is interactive and does not currently use the app-level
PNG/SVG download buttons.

## Build the documentation locally

From the app directory:

```bash
mkdocs build --strict
```

To preview locally:

```bash
mkdocs serve
```
