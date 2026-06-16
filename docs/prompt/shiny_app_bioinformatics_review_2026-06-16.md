# ProxiomeVis Shiny App Bioinformatics Review

Date: 2026-06-16

Reviewer perspective: senior bioinformatician reviewing the Shiny app with the Raji/CAR-T demo as the default dataset.

## Evidence Used

- Live Shiny app: `http://127.0.0.1:3853`
- Default loaded source: `Raji/CAR-T patch demo | 4,470 cells | 159 markers`
- Browser screenshots: `screenshots/ui-audit-2026-06-16/`
- Automated R tests after the default-data change: `PASS 883`, `FAIL 0`
- Raji cache metadata checked from `cache/raji_cart_demo_data.rds`

The browser pass visited QC, Abundance, Spatial Metrics/Clustering, Spatial Metrics/Colocalization, Patch Analysis, the Data popover, Options popovers, and visible plot/download controls. Download links were checked for presence and hrefs; full PNG/SVG file fetching was skipped in the broad pass because it made the scan too slow.

## Highest Priority Issues

1. **The Raji demo is now the startup dataset, but the Data popover status is stale.**
   - Observed: the source chip shows Raji/CAR-T correctly, but the Data popover still says `Enter an RDS path, then click Load Data.`
   - Fix: after startup `load_raji_demo_into_app()`, also call `send_rds_load_state("success", "Raji/CAR-T demo data is active.", progress = 100)`.
   - File: `R/data_source_module.R`

2. **The app hides the useful Raji biology because it keys many controls to `celltype_manual`, which is all `unannotated`.**
   - Raji cache has richer columns: `celltype`, `celltype_condition`, `cd8_state`, `contains_raji`, `timepoint_hr`, `ratio`, `donor`, and marker-derived scores.
   - `celltype_manual` has 4,470/4,470 cells as `unannotated`.
   - `celltype` separates `Other`, `Raji`, `CD8T`, `CD4T`, and `T_other`.
   - `celltype_condition` separates `CD8T_alone`, `CD8T_coculture`, `Raji_coculture`, etc.
   - Fix: add an active annotation field helper, prefer `celltype_manual` when informative, otherwise fall back to `celltype` or expose a "Cell annotation" selector.
   - Files: `R/data_adapter.R`, `R/abundance_module.R`, `R/clustering_module.R`, `R/colocalization_module.R`

3. **Colocalization Observed does not show the heatmap in the browser pass.**
   - Observed: table rows render, but the heatmap area is blank in `screenshots/ui-audit-2026-06-16/colocalization-observed.png`.
   - The Plot style radio also appears unselected in the screenshot, despite `Interactive` being the intended default.
   - Fix: verify `colocalization_heatmap_display` initializes inside the conditional sidebar, add an explicit fallback to interactive rendering, and add a small UI test that the observed heatmap output is visible after selecting the Colocalization tab.
   - File: `R/colocalization_module.R`

4. **Differential tabs initially show empty plot panes with active PNG/SVG/Options controls.**
   - Observed in Abundance, Clustering, and Colocalization Differential screenshots.
   - This is confusing because users see export controls before any result exists.
   - Server warnings also show `plotly_click` event sources are used without `event_register()`.
   - Fix: show a clear empty state before running differential analysis, hide or disable download controls until plots exist, and register click events on volcano plots.
   - Files: `R/differential_helpers.R`, `R/abundance_module.R`, `R/clustering_module.R`, `R/colocalization_module.R`

5. **3D Layout defaults to a sample with no layout file.**
   - Observed: `No Pixelator 3D layout file found for PNA065_CARTcells_Raji_24hrs_1to1_S05`.
   - Server warning: `colocalization_3d_component` has a large number of options and should use server-side selectize.
   - Fix: detect samples with layout files and select one by default; if none exist, show a tab-level empty state and hide Options/export sizing. Convert component selection to server-side selectize.
   - File: `R/colocalization_module.R`

## Tab-by-Tab Bioinformatics Enhancements

### QC

- Add a per-sample QC summary table with loaded cells, retained cells, median UMIs, median antibodies, isotype fraction, and intracellular fraction.
- Show filter losses by sample, donor, condition, timepoint, and ratio. This is important for the Raji/CAR-T demo because donor and coculture condition can confound downstream comparisons.
- Add recommended threshold bands or visual cutoff lines to the molecule rank and distribution plots.
- Add a "QC outlier samples" callout based on median UMIs, retained fraction, and isotype fraction.

### Abundance

- Make `celltype` or `celltype_condition` the default annotation for the Raji demo instead of `celltype_manual`.
- Add marker program views for T cell, CD8, CD4, Raji/B cell, CAR/FMC63, activation, and exhaustion scores.
- Add condition-aware abundance panels: cart alone vs coculture, 4 h vs 24 h, 1:1 vs 5:1, donor PNA065 vs PNA066.
- For Marker Distributions, add percent-positive summaries alongside abundance. A thresholded detection metric is often easier to interpret than abundance alone.
- For Differential, prefer a replicate-aware or sample-level summary when possible. Per-cell tests can overstate significance when cells are treated as independent replicates.
- Add labeled marker callouts on volcano plots for key markers: `CD3e`, `CD4`, `CD8`, `FMC63`, `CD19`, `CD20`, `CD40`, `CD54`, `CD25`, `CD279`, `HLA-DR-DP-DQ`.

### Spatial Metrics - Clustering

- Add a short interpretation label for self-clustering log2 ratio: positive means same-marker proximity is enriched, negative means depleted.
- Add replicate-aware summaries by donor/sample rather than only pooled cell-level summaries.
- Add a paired marker ranking table: marker, condition, median self-clustering, effect vs reference, adjusted p value, n cells.
- Add confidence intervals or sample-level spread for top markers.
- Default the cell-type dimension to the informative Raji annotation field, not `celltype_manual`.

### Spatial Metrics - Colocalization

- Fix the observed heatmap default rendering first.
- Keep the heatmap symmetric for observed pairwise marker results and make the symmetry rule explicit in the plot caption or tooltip.
- Add pair-level detail views for selected marker pairs: abundance of marker 1, abundance of marker 2, observed proximity, expected/null proximity, and log2 enrichment.
- Add a marker-pair network view for the strongest condition-specific colocalizations.
- Add top gained/lost pair tables by condition and by timepoint/ratio.
- Add an option to focus on Raji interaction biology: CAR/T markers vs B-cell/Raji markers.

### Spatial Metrics - 3D Layout

- Select only samples/components that have a readable layout file.
- Add component search by sample, condition, cell type, UMI count, and marker signal.
- Add a compact legend explaining highlighted markers and background nodes.
- Add presets for Raji contact biology: CAR-T markers, Raji/B-cell markers, activation markers, and exhaustion markers.

### Patch Analysis

- Patch Markers: add threshold lines and a compact classification summary for receiver-enriched, target-enriched, and unspecific markers.
- Raji Signal: summarize log2 proximity by condition, timepoint, ratio, donor, and inferred cell type.
- Patch Burden: add plots. The current tab is a raw table; users need a burden distribution plot and an aggregate condition summary.
- Add a clear method box showing how patch markers were selected, how patch nodes were called, and what `patch_node_fraction` means.
- Add a top-cell/component drilldown for high patch burden cells.

## Global Enhancements

- Add a dataset provenance panel: cache path, source type, cells, markers, samples, conditions, donors, annotation field used, and cache schema version.
- Add a parameter summary near exported plots or as a sidecar CSV/JSON: selected markers, filters, thresholds, annotation field, and display/export dimensions.
- Add a small automated browser smoke test for the main tabs, Options popovers, and one PNG/SVG download.
- Add empty states everywhere a plot depends on a button click or an unavailable source file.
- Add a "Biology view" preset for Raji/CAR-T that sets annotation to `celltype_condition`, selects Raji/CAR-T markers, and uses condition/time/ratio as the primary grouping variables.

## Suggested Implementation Order

1. Fix startup Data popover status and active annotation-field selection for the Raji demo.
2. Fix Colocalization Observed heatmap rendering and add its UI regression check.
3. Add empty states and disabled downloads for Differential tabs before analysis is run.
4. Fix 3D Layout defaults and convert component selection to server-side selectize.
5. Add Raji/CAR-T biology presets and Patch Burden plots.
6. Add a lightweight browser smoke test for tab navigation, Options popovers, and one PNG/SVG export.

