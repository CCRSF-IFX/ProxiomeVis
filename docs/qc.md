# QC

Use the QC tab to check whether the processed H5AD looks reasonable before
interpreting abundance.

## Filtering

The **Filtering** view summarizes how many cells remain after each QC step. Use
the sidebar to choose samples and switch the y-axis between number of cells and
fraction of loaded cells. If the H5AD does not contain notebook-generated QC
history, the app reports that the history is unavailable rather than inferring
filtering counts from the final cells.

## Cell Calling

The **Cell Calling** view shows molecule rank curves. Use the `n_umi reference
line` control to inspect how a candidate threshold relates to the ranked cell
distribution. Reference lines are visual guides and never filter or modify the
loaded cells.

## Distributions

The **Distributions** view shows selected QC metrics across samples. Use this
view to check whether one sample has an unusual distribution before comparing
biological readouts.

## Metadata

The **Metadata** view shows the processed `obs` rows for inspection. This is useful
when validating sample labels, conditions, and cell annotations.
