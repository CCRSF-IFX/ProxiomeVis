# QC

Use the QC tab to check whether the processed H5AD looks reasonable before
interpreting abundance.

## Filtering

The **Filtering** view summarizes how many cells remain after each QC step. Use
the sidebar to choose samples and switch the y-axis between number of cells and
fraction of loaded cells.

## Cell Calling

The **Cell Calling** view shows molecule rank curves. Use the `n_umi cutoff`
control to inspect how the selected threshold relates to the ranked cell
distribution.

## Distributions

The **Distributions** view shows selected QC metrics across samples. Use this
view to check whether one sample has an unusual distribution before comparing
biological readouts.

## Metadata

The **Metadata** view shows the processed `obs` rows for inspection. This is useful
when validating sample labels, conditions, and cell annotations.
