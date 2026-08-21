# User Guide and Parameter Reference

This guide applies to ProxiomeVis Python 0.1.0. It explains the intended
workflow, when controls take effect, and how to interpret each analysis. The
app is exploratory: confirm important findings with the study design, upstream
Pixelgen outputs, and an appropriate statistical analysis.

## Before You Begin

ProxiomeVis combines two data sources:

| Source | Used for |
| --- | --- |
| Processed H5AD | Cell metadata, raw marker counts, normalized abundance, embeddings, optional QC history, annotations, and optional stored patch results |
| Matching PXL file(s) | Precomputed proximity rows, patch candidate signal, and precomputed 3D component layouts |

The H5AD and PXL files must describe the same experiment and share component
identifiers. The app ignores any proximity table embedded in H5AD and queries
the assigned PXL files instead. All changes are session-local; the source files
are never modified.

## Quick Start

1. Open **Data** in the top navigation.
2. Enter one server-visible processed `.h5ad` path and the matching `.pxl`
   path, directory, or glob.
3. Click **Inspect** to validate the H5AD path, then click **Load Data**.
4. Confirm the loaded cell, marker, and PXL counts in the data-source summary.
5. Review **QC**, then verify annotations and normalized signal in
   **Abundance**.
6. In **Spatial Metrics > Retrieve Data**, choose the population and marker
   scope and click **Retrieve Spatial Data**.
7. Continue in order through **Proximity Profile**, **Clustering**,
   **Colocalization**, and **3D Layout**.
8. For patch screening, select the receiver, target, groups, and markers in
   **Patch Analysis**, then click **Prepare Patch Data**.
9. Use **Activity Log** to monitor long operations. Download diagnostics if an
   error displays a reference ID.

## Understand When Controls Take Effect

The app uses three control behaviors. Check the status message whenever a plot
does not reflect a recent change.

| Control behavior | Examples | When results change |
| --- | --- | --- |
| Live display or filter | QC reference lines, plot facets, group display, legend limits, differential thresholds | Immediately |
| Applied analysis settings | Proximity population, mean calculation, marker filters, heatmap marker set | After **Apply analysis settings** |
| Prepared retrieval scope | Spatial population and markers; patch receiver, target, groups, and markers | After **Retrieve Spatial Data** or **Prepare Patch Data** |

An **Unapplied changes** warning means the visible results still use the last
applied or prepared settings. Changing a loaded dataset or analysis grouping
clears prepared spatial and patch scopes.

## Data and Analysis Grouping

### Data controls

| Parameter | Meaning | Recommended use |
| --- | --- | --- |
| **Processed .h5ad path** | Full path visible to the Shiny server | Use one processed AnnData file; browser uploads are not supported |
| **PXL path(s) for proximity and cellgraph data** | One PXL file, directory, glob, or comma/newline-separated list | Assign the PXL files produced for the same samples as the H5AD |
| **Inspect** | Checks that the H5AD path is readable and reports file size | Use before loading a new dataset |
| **Load Data** | Reads H5AD and resolves PXL files | Required after changing either path |

### Analysis grouping

**Analysis grouping** defines the `condition` used throughout plots and
contrasts. **Use metadata column** is appropriate when a metadata value is
constant within each sample. **Edit sample groups** allows an explicit group
for every sample. **Reset to condition** restores the original H5AD condition.

Use biologically meaningful sample-level groups. Do not encode individual
cells as independent treatment groups. Grouping changes do not modify the
H5AD, but they reset active spatial and patch preparations.

## QC

QC describes the cells already present in the processed H5AD. It does not
rerun filtering or cell calling.

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Sample** | All samples | Restricts QC plots and tables to selected samples |
| **n_umi reference line** | 10,000 | Visual guide on molecule-rank and QC plots; it does not remove cells |
| **Isotype fraction reference line** | 0.001 | Visual guide for isotype fraction; it does not remove cells |
| **Filter count y-axis** | Number of cells | Switches stored QC history between cell counts and fraction of loaded cells |
| **Distribution metric** | First available QC field | Selects the H5AD observation field summarized across samples |

**Filtering** displays `uns["qc_cell_counts_by_step"]` when available. If the
upstream history is missing, the app reports **QC history unavailable** and
does not reconstruct an apparently real history from retained cells.

**Cell Calling** shows ranked `n_umi` values. **Distributions** helps identify
sample-specific shifts in QC fields. **Metadata** exposes processed H5AD
observation rows for label and annotation checks.

## Abundance

The app uses normalized abundance from `layers["clr"]` when present; otherwise,
it uses the available H5AD expression matrix. Raw counts remain available for
count-based marker filters and patch screening.

### Observed

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Embedding** | First available embedding | Selects a stored two-dimensional UMAP, PCA, Harmony, or t-SNE representation |
| **Color UMAP by** | Marker abundance | Colors cells by normalized marker abundance, cell type, analysis group, or sample |
| **Marker** | First marker | Marker used when coloring by abundance |
| **Split UMAP by** | None | Facets the same embedding by analysis group or sample |
| **Split columns** | 2 | Controls facet layout only |
| **Dot size** | 3 | Changes point rendering only |
| **Analysis group** | All groups | Live cell filter |
| **Cell type** | All cell types | Live cell filter |

### Marker Distributions and Cell Annotation

**Marker Distributions** shows one marker by cell type, colored by analysis
group and faceted by sample. **Facet columns** defaults to 3; **Show jitter
dots** defaults to on. Jittered dots are cells, not biological replicates.

**Cell Annotation** reports cell-type fractions within each analysis group and
the median normalized abundance of every marker within each cell type. Large
differences in cell counts can affect pooled summaries.

### Abundance Differential

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Group A** | Data-dependent | Numerator group in the reported effect |
| **Group B (reference)** | `NC` when available | Reference group; effect is median A minus median B |
| **Cell type** | One cell type initially | Cells included in the comparison |
| **Stratify by cell type** | Off | Off pools selected cell types; on runs a separate test for each cell type |
| **FDR threshold** | 0.05 | Classifies volcano hits using Benjamini-Hochberg adjusted p-values |
| **Minimum effect** | 0.25 | Requires the absolute median normalized-abundance difference to reach this value |
| **Minimum observations per group** | 3 | Minimum cells required in each group for a Mann-Whitney U test |
| **Detail marker** | Data-dependent | Marker shown in the companion violin plot; clicking a volcano point also selects it |

Abundance differential uses cells as observations. It does not model sample
pairing or repeated measures. Treat its p-values as exploratory when the study
contains biological replicates, and confirm findings with a sample-aware
method appropriate to the design.

## Spatial Metrics

Follow this sequence:

**Retrieve Data → Proximity Profile → Clustering → Colocalization → 3D Layout**

All downstream tabs are bounded by one active retrieval. Downstream controls
can narrow that scope but cannot add cells or markers that were not retrieved.

### Retrieve Data

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Analysis group** | All groups | Groups included in the active spatial scope |
| **Sample** | All samples | Samples included in the active spatial scope |
| **Cell type** | First cell type | Cell populations included; start with one population for interpretable summaries |
| **Number of markers: All markers** | Selected | Makes every H5AD marker available downstream |
| **Top abundance markers** | 40 | Retrieves markers ranked by sample-weighted mean normalized abundance |
| **Selected markers** | None | Retrieves only the explicit marker set |

Click **Retrieve Spatial Data** after any retrieval change. The active summary
reports the frozen population and marker scope. The server extracts every PXL
proximity row matching those cells and markers into the active retrieval. Plots
receive summarized results rather than the full table in browser memory. For
large panels, reduce the cell or marker scope if retrieval time or server memory
is constrained.

### Interpret Proximity Scores

PXL proximity is represented as a log2 observed-to-expected ratio:

| Value | Interpretation |
| --- | --- |
| Greater than 0 | Markers are closer than expected by the PXL null model |
| Near 0 | Approximately expected spatial organization |
| Less than 0 | Markers are spatially segregated relative to expectation |

The magnitude is model-relative and should not be interpreted as physical
distance. A missing marker pair can mean it was not detected or did not pass
eligibility filters; it is not automatically evidence of segregation.

### Proximity Profile: Applied Analysis Settings

The following controls remain drafts until **Apply analysis settings** is
clicked.

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Cell population** | First population in retrieval | Exactly one cell type used for proximity summaries and downstream colocalization |
| **Mean calculation: Population mean** | Selected | Sum of detected log2 ratios divided by all cells in the group; undetected pairs contribute zero |
| **Detected-cell mean** | Not selected | Mean among cells with a proximity row; emphasizes score strength conditional on detection |
| **Apply Pixelator marker filters** | On | Requires both markers to meet per-cell count and fraction thresholds |
| **Minimum marker fraction** | 0.001 | Minimum raw count fraction for each marker within a cell |
| **Minimum marker count** | 0 | Minimum raw count for each marker within a cell |
| **Summarize by** | Analysis group | Aggregates the heatmap by analysis group or by sample |
| **Marker set: PixelatorES proximity profile** | Selected | Takes markers from the 60 strongest absolute-mean pairs detected in more than 50% of cells |
| **Strongest proximity pairs** | 60 | Number of qualifying pairs used to construct the profile marker set |
| **Top abundance markers** | 40 | Alternative sample-weighted abundance display set |
| **Selected markers** | None | Explicit marker set; at least two markers are required |
| **Minimum detected cells per summarized group** | 1 | Entries below this support are shown as zero |

**Load PixelatorES defaults** resets these values but does not apply them.
Click **Apply analysis settings** afterward. Apply filters and summarizes the
active retrieval in memory; it does not query the PXL files again.

For top abundance markers, the app first calculates each marker's mean
normalized abundance within each sample and then averages those sample means.
Each sample therefore has equal weight in marker ranking. This top-marker rule
is a display filter only; it does not determine whether a PXL proximity score
was calculated.

### Proximity Profile: Live Display Controls

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Group display** | Focused group | Shows one summarized group or compares multiple groups |
| **Displayed group / Groups to compare** | Data-dependent | Chooses visible panels without changing applied analysis |
| **Clustering reference group** | `NC` when available | Defines the group used to calculate one shared marker order for all panels |
| **Marker ordering** | Ward | Hierarchical-linkage method used only for row and column ordering |
| **Legend minimum / maximum** | -1 / 1 | Color scale limits; values outside the range are clipped visually, not filtered |
| **Pair detail** | First available pair | Selects the sample-level detail plot; clicking a heatmap dot also selects the pair |

The same reference-derived order is used when **Displayed group** changes,
which makes panels comparable. The heatmap column direction is reversed to
match PixelatorES. Dot color encodes mean log2 ratio and dot size encodes the
fraction of cells with the pair.

Use **Settings JSON** to preserve the active retrieval, applied analysis, and
display configuration with exported results.

### Clustering

Clustering is marker self-proximity (`marker_1 = marker_2`), not unsupervised
cell clustering.

| Control | Interpretation |
| --- | --- |
| **Marker** | Marker whose self-proximity is shown in the violin plot |
| **Analysis group / Cell type** | Live filters within the active retrieval |
| **Protein set** | Uses the top variable proteins by default or an explicit custom list from the active spatial retrieval |
| **Top proteins** | Number ranked by variation in mean self-proximity; default 20, maximum 40 |
| **Proteins** | Ordered custom list used when **Custom proteins** is selected |
| **Differential contrast** | Median self-proximity in group A minus group B |

Clustering differential uses a two-sided Mann-Whitney U test on detected
cell-level self-proximity rows with Benjamini-Hochberg correction. Its default
FDR, minimum effect, and minimum observations are 0.05, 0.25, and 3. Missing
self-proximity rows are not inserted as zeros. As with abundance differential,
these are cell-level exploratory statistics rather than sample-aware tests.

### Colocalization Differential

Colocalization reuses the cell population, mean type, and marker-filter
settings last applied in **Proximity Profile**.

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Group A** | Non-reference group when available | Numerator group |
| **Group B (reference)** | `NC` when available | Reference group |
| **Minimum samples per group** | 2 | Minimum biological samples required in each group for p-values and FDR |
| **Pairs shown** | All marker pairs | Tests all retrieved pairs or only pairs containing one anchor marker |
| **Marker** | Data-dependent | Anchor marker when the restricted pair scope is selected |
| **Detail pair** | First available pair | Pair shown in the sample-level box plot; clicking a volcano point selects it |
| **FDR threshold** | 0.05 | Volcano hit threshold after Benjamini-Hochberg correction |
| **Minimum median-sample difference** | 0.25 | Required absolute difference between group medians of sample-level means |

This is the sample-aware spatial comparison: one summarized value per sample
and marker pair is the observation. The effect is median sample value in group
A minus group B. If either group has too few samples, descriptive effects may
remain available, but p-values and FDR are not reported.

### 3D Layout

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Sample** | First available sample | Chooses the matching assigned PXL file |
| **Cell type** | One cell type | Narrows component choices within the active retrieval |
| **Cell/component** | First available component | Cell graph rendered from the stored precomputed layout |
| **Max background nodes** | 7,000 | Caps non-highlighted nodes for browser performance; highlighted nodes remain visible |
| **Highlighted markers** | Data-dependent | Draws selected marker nodes larger and in color |

The 3D view is descriptive. Node positions are the stored PXL layout, not
microscopy coordinates or measured physical distances.

## Patch Analysis

Patch Analysis is an experimental receiver/target workflow. A candidate signal
is not a detected patch. Connected-subgraph patch detection must be run
upstream and stored in the H5AD for **Detected Patches** and **Composition**.

### Prepared Patch Scope

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Population metadata** | Preferred annotation field | Categorical H5AD field defining receiver and target populations |
| **Receiver population** | Choose a population | Cells on which candidate or detected patches are evaluated |
| **Target population** | Choose a different population | Population whose markers define the putative patch signal |
| **Analysis groups** | All groups | Groups included in marker screening, PXL candidates, and stored results |
| **Minimum population fraction** | 0.01 | Raw marker fraction required before a marker can be suggested |
| **Minimum fold enrichment** | 3 | Target-to-receiver or receiver-to-target enrichment required for suggestion |
| **Target/patch markers** | Suggested, editable | Markers expected from the target population; at least two are required |
| **Receiver/blocking markers** | Suggested, editable | Receiver-specific markers retained as provenance for the analysis scope |

Marker suggestions use pooled raw count fractions, not PixelatorR abundance
unmixing. Aim for selected target markers that collectively cover at least
20–30% of target-population counts. Click **Use suggested markers** to restore
the current suggestions, then click **Prepare Patch Data** to freeze the scope
and query receiver-cell candidate signal.

### Candidate and Stored-Result Controls

| Parameter | Default | Interpretation |
| --- | --- | --- |
| **Minimum target-marker count** | 100 | Display threshold for summed target-marker raw counts in a receiver cell |
| **Minimum joint log2 ratio** | 0.3 | Display threshold for joint target-marker proximity |
| **Patch metric** | Preferred stored numeric field | Metric plotted from stored patch burden or derived patch-size summaries |
| **Group results by** | Preferred stored categorical field | Category used for the burden comparison |

Candidate thresholds update live and do not require preparation again. A cell
must meet both thresholds to be labeled **Candidate**. This label is a review
screen only and does not call a connected graph patch.

Stored patch results are restricted to the active receiver and analysis-group
scope. If no patch tables are stored under `uns["proxiome"]["patch"]`, the app
reports them as unavailable instead of estimating patches from proximity.

## Tables, Downloads, and Reproducibility

- Plot **SVG** exports are vector images suitable for reports, slides, and
  downstream editing.
- Every table provides **CSV** and **Excel** downloads for its complete
  underlying result. Browser grids remain capped at 2,000 rows for performance.
- Interactive tables can be sorted and filtered without changing the analysis.
- Record the H5AD name, PXL names, analysis grouping, active retrieval summary,
  applied proximity settings, contrasts, and app commit with every report.
- Use **Settings JSON** for proximity-profile provenance.

## Activity Log and Diagnostics

**Activity Log** reports session-local operations with timestamp, severity,
status, elapsed time, app version, and Git commit. It is the first place to
check when a PXL query appears slow.

If the app reports a reference ID, click **Download diagnostics** and provide
the ZIP to the application administrator. It contains sanitized runtime
metadata, structured activity, and tracebacks for the current browser session.
It does not include H5AD/PXL contents or full server paths.

## Interpretation Checklist

Before reporting a result, confirm all of the following:

- The H5AD and PXL files are from the same processed experiment.
- Sample labels, analysis grouping, and cell annotations are correct.
- The visible spatial result uses the intended active retrieval and applied
  proximity settings; no unapplied-changes warning is present.
- The chosen mean type is stated: population mean or detected-cell mean.
- Marker filters, marker-selection method, and clustering reference are stated.
- Sample-level inference is used when biological samples are the replicates.
- Candidate patch signal is not described as detected patches.
- Figures are accompanied by their underlying table or settings record.
