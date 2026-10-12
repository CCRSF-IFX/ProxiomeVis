# ProxiomeVis R

ProxiomeVis is an R Shiny report for processed Pixelgen data. It combines the
annotated Seurat object from the demo CCA analysis with frozen abundance,
self-proximity, and colocalization results from the tertiary analysis.

The app provides:

- QC cell counts and stored cell annotations;
- interactive embeddings colored by annotation or marker abundance;
- per-sample summaries, donor-specific contrasts, and stored cross-donor
  descriptive effects for the three readouts;
- downloadable tables for the currently displayed stored results.

The app does not open PXL files or recalculate proximity scores. A preparation
script builds one versioned runtime RDS, and every Shiny session reads that
prepared file.

## Interpretation boundary

The demo contains two donors. Combined donor-matched effects are descriptive,
contain no p-values, and are not population-level inference. Donor-specific
report volcanoes reproduce stored exploratory component-level p-values from the
tertiary analysis; Shiny does not recalculate them. Cells are not independent
biological replicates.

Start with [Getting started](getting-started.md), then review the
[precomputed-results contract](precomputed-results.md) before preparing a new
dataset.
