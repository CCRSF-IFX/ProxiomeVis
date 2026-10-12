# Troubleshooting

## The app cannot find its RDS

Build `output/proxiomevis_demo_v1.rds`, or set
`PROXIOME_PRECOMPUTED_RDS` to an absolute path readable by the Shiny process.
The configured file must exist and end in `.rds`.

## The RDS fails contract validation

Rebuild it with `scripts/prepare_precomputed_rds.R`. The app accepts schema
version 1 and checks metadata mappings, reductions, assay layers, result-table
columns, marker-pair types, support counts, and descriptive-inference rules.
It deliberately rejects an old or partially assembled payload.

## R reports an unknown PNA class

Install a `pixelatorR` version compatible with the source analysis. The app
does not query PXL, but R must register pixelatorR's S4 classes before restoring
the saved Seurat object.

## A filter produces no rows

The controls only narrow the frozen tertiary results. Choose a donor,
comparison, cell type, and feature present in the prepared tables. The app
cannot compute a missing combination from individual cells.

## PXL files are unavailable

This is expected at runtime. PXL files are upstream analysis inputs and are not
read by precomputed mode.

## The port is unavailable

Use another local port, for example:

```bash
Rscript --vanilla -e 'shiny::runApp(".", host="127.0.0.1", port=3839)'
```

For a shared deployment, keep Shiny behind the environment's authenticated
proxy. Do not expose the server-local RDS path control or unrestricted RDS
uploads to untrusted users.

## The server-path field is missing

This is the secure default. Configure `PROXIOME_PRECOMPUTED_RDS` before startup.
A trusted local operator can set `PROXIOME_ALLOW_LOCAL_LOAD=true` to expose the
field temporarily. Do not use that flag on a public app.
