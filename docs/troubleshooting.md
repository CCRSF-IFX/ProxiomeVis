# Troubleshooting

## The app cannot find its RDS

The default is `output/hybrid_demo/analysis.rds` when present, otherwise
`output/proxiomevis_demo_v1.rds`. Follow [Getting started](getting-started.md) to
prepare these, or set
`PROXIOME_PRECOMPUTED_RDS` to an absolute path readable by the Shiny process.
The configured file must exist and end in `.rds`.

## The RDS fails contract validation

The app accepts a plain application bundle with `bundle_version = 2L` or a
legacy Seurat RDS. Both retain the frozen contract's `schema_version = 1L`;
these version numbers describe different structures and should not be changed
to make an incompatible file load.

Validation checks metadata mappings, coordinates, matrix/cell alignment,
result-table columns, marker pairs, support counts, and descriptive-inference
rules. Legacy inputs additionally need their configured assay layers and empty
cell-level proximity slots. Rebuild the frozen demo contract with
`scripts/prepare_precomputed_rds.R`, then export a plain bundle as described in
[Precomputed results](precomputed-results.md). Do not patch invalid results
merely to bypass validation.

## R reports an unknown PNA class

This applies to legacy analysis files, not plain application bundles. Install
SeuratObject and a pixelatorR version compatible with the source analysis in the
R environment that runs the app, or export a plain bundle in the analysis
environment and deploy that instead. Plain bundles need no Seurat/PNA classes
for viewing; they still need the viewer packages in `DESCRIPTION`.

Check which R executable starts Shiny if packages appear installed but cannot
be found. This workspace's project environment is
`../.pixi/envs/r/bin/Rscript`; it can differ from the system `Rscript`.

## A filter produces no rows

Report controls only narrow the frozen tertiary results. Choose a donor,
comparison, cell type, and feature present in those tables. A missing stored
combination is not computed automatically.

To run a new comparison, use **Custom analysis** and explicitly submit a valid
cell selection. Its result is separate from the report and requires the
appropriate cell-level inputs and worker dependencies.

## PXL files are unavailable

PXL files are preparation inputs, not runtime requirements. Report views use
the RDS, spatial tests use the prepared proximity sidecar, and **Cell 3D** loads
the chosen cell's nodes from its sample's prepared Parquet file. If preparation
itself fails because PXL paths moved, pass their directory as the third argument to
`scripts/prepare_hybrid_bundle.R`.

## Spatial analysis or Cell 3D says assets are missing

A report-only RDS does not contain the optional sidecars. Prepare them upstream,
then deploy `analysis.rds` together with its `proximity/` and `graphs/`
directories. Asset descriptors are relative to the bundle directory; absolute
paths, traversal, and symlinks escaping that directory are rejected.

Graph coverage is shown explicitly. If preparation is still running, click
**Refresh prepared cells** to see newly completed cells. Rerun the preparation
command after an interruption to reuse validated completed sample files. Adding
cells atomically rewrites affected files with the union of old and new cells;
it does not remove previously prepared cells from coverage. Do not run two
preparers concurrently in the same output directory. A failed partial
proximity export without a manifest needs a fresh output directory; the
standalone proximity preparer refuses to overwrite a nonempty directory.

## Cell 3D reports an unsupported graph index

The current graph sidecar uses `schema_version = 3L`: one node-only Parquet file
per sample plus `index.rds`, with no stored node IDs. Version 2 sample Parquet
sidecars remain readable but cannot be extended by the preparer; rebuild in a
new directory to get smaller version 3 files. Old version 1 per-cell RDS assets
are not compatible.
Do not change the version field by hand or delete the old directory to bypass
the check. Preserve it and run `scripts/prepare_hybrid_bundle.R` with a new
output directory, then deploy the new bundle and its sidecars together. There
is no automatic destructive migration.

## Cell 3D is missing Arrow/dplyr or is slow on first load

Install `arrow` and `dplyr` in the R environment running Shiny. Frozen reports
still work without these optional graph-query dependencies. No pixelatorR or
source PXL access is needed for viewing a plain bundle's prepared nodes.

The first selection from a sample may checksum its whole Parquet file. Later
selections reuse validation for the unchanged file and use Arrow to filter the
cell and prune unrelated row groups before collecting nodes. A checksum scan
does not collect the full sample into R. Keep deployed assets read-only; do not
disable validation to remove the initial disk scan.

## Graph preparation stops with exit 137 and no R error

Check the host's resource limits before assuming the graph data is corrupt.
Use `ulimit -t` or scheduler diagnostics to check the CPU-time allowance. A full
export can hit a time limit even with bounded memory. Do not change shared-host
limits. Resume completed samples or submit smaller, sequential cell batches
through the existing preparer:

```bash
Rscript --vanilla scripts/prepare_cell_graph_assets.R \
  output/proxiomevis_demo_v1.rds /path/to/pxl_directory \
  output/hybrid_demo/graphs pixelgen-demo-v1-d20537e8d79c4e12 cells.txt
```

`cells.txt` contains one exact prepared cell ID per line. Use the dataset ID
from your source contract for non-demo data. Each invocation validates and
retains previously published cells, atomically rewrites affected sample files
with the combined selection, and updates the same index. Run batches
sequentially, never against the same directory
concurrently. Exit 137 can also indicate a memory limit or an external kill;
check the scheduler/system diagnostics rather than repeatedly restarting a
healthy process.

## A worker cannot start or reports a missing package

The analysis subprocess uses the app's R installation and library paths. New
tests require pixelatorR, Seurat, and SeuratObject; spatial querying also needs
Arrow and dplyr. The report-only viewer does not install these optional packages.
See [Hybrid analysis](hybrid-analysis.md#dependencies-and-preparation).

Missing worker dependencies or sidecars do not invalidate the frozen reports.
Install the required packages in the correct environment, then retry the
explicit comparison. The app does not substitute a different test.

## The worker is busy, times out, or rejects a large selection

Only one background analysis runs at a time **per Shiny R process**, shared by
that process's sessions. Busy requests are rejected with a retry message, not
queued. Retry after the current job finishes or is cancelled. Closing a session
cancels its running analysis.

Each job has a 120-second timeout and limits of 2,000 selected cells, 160 markers,
1,000 unordered colocalization pairs, and 1 million retained score rows. Select
fewer cells or markers to reduce work; colocalization is limited to 45 markers
because 46 produce more than 1,000 pairs. The app never silently downsamples a
request. Deployment-level RAM/CPU quotas are still needed: these workload limits
are not a memory guarantee, and multiple Shiny processes each have their own
worker limit.

## A new result has missing statistics or differs from the report

Abundance requires at least 3 selected cells per condition. Spatial features
require at least 20 retained observations per group for self-proximity, or 40
for colocalization. Inspect `n_cells_a/b`, `n_scored_a/b`, and `status` in the
result table. Missing or filtered scores are not imputed as zeros, so selecting
enough cells does not guarantee every feature is testable.

New tests use the submitted cells and eligible selected-marker testing family.
Restricting markers can change BH-adjusted p-values even when raw p-values agree
with the report. Changing controls does not change the completed result until
another run; check its recorded scope and downloaded provenance. Neither result
supports population-level inference from cells in a single donor.

## The bundle or sidecar checksum does not match

Reload the trusted dataset if its RDS changed after the app loaded it. Requests
are pinned to that loaded file's checksum, not just its dataset name. Sidecars
must match the bundle's data and prepared source inputs; do not mix files from
different preparation runs or edit checksums to suppress the error. Prepare a
new matching bundle in a separate directory when source data changes, and
deploy it read-only with its sidecars.

## A graph is too large or the 3D panel is blank

The viewer refuses cells above 200,000 nodes. Select a smaller cell to continue;
no nodes are silently dropped. The sidecar and viewer are nodes only: there are
no stored edges, edge controls, or marker-count matrices.

Plotly 3D requires a browser with working WebGL. Check browser GPU/WebGL support
when graph coverage and node totals load but the canvas stays blank. The coordinates
are abstract graph-layout positions, not physical microscopy coordinates.

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
