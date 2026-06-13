# ProxiomeVis Full Repository Review

Review date: 2026-06-11

Reviewer perspective: senior R/Shiny engineer, with focus on the Shiny app in `shiny/proxiome_demo` and its deployment inside the broader Pixelgen CAR-T/Raji reproduction repository.

## Review Scope And Evidence

Inspected source and configuration:

- Shiny app: `shiny/proxiome_demo/app.R`, `R/data_adapter.R`, `R/environment_module.R`, `R/plot_layout.R`
- Shiny deployment: `shiny/proxiome_demo/template/script.sh.erb`, `.Rprofile`, `.renvignore`, `renv.lock`, `renv/settings.json`, app README
- Shiny tests: `shiny/proxiome_demo/tests/testthat/*.R`
- Repository workflow: `README.md`, `pixi.toml`, `configs/analysis_config.yaml`, `scripts/*.py`, `scripts/*.R`, `env/*.R`, `.github/workflows/ci_mkdocs.yaml`
- Generated notebooks, figures, tables, and large data/results directories were inventoried but not audited line-by-line as maintainable source.

Verification commands run:

```bash
module load R/4.5.0 && R_BIN=$(command -v Rscript) && cd shiny/proxiome_demo && \
  /usr/bin/env RENV_PROJECT=$PWD RENV_PATHS_LIBRARY=$PWD/renv/library \
  RENV_CONFIG_CACHE_ENABLED=FALSE "$R_BIN" \
  -e 'source("renv/activate.R"); testthat::test_dir("tests/testthat", reporter = "summary")'
```

Result: passed. Output ended with `DONE` after `app_ui`, `data_adapter`, and `plot_layout` tests.

```bash
bash -n shiny/proxiome_demo/template/script.sh.erb
Rscript --vanilla -e 'parse("shiny/proxiome_demo/app.R"); parse("shiny/proxiome_demo/R/data_adapter.R"); parse("shiny/proxiome_demo/R/environment_module.R"); parse("shiny/proxiome_demo/R/plot_layout.R"); cat("parse ok\n")'
```

Result: passed. Output included `parse ok`.

## Architecture Overview

### Repository

The repository has two related but distinct surfaces:

1. The reproducibility workflow at the repository root.
   - Entry points are `pixi.toml` tasks, especially `analysis`, `notebooks`, `reproduce`, and `build-docs`.
   - `scripts/pixelgen_repro.py` is the central Python pipeline. It reads Pixelator `.layout.pxl` files through DuckDB, annotates cells, computes abundance/proximity summaries, writes figures/tables, and copies outputs into `docs/assets/generated`.
   - `scripts/build_notebooks.py` generates and optionally executes notebooks from reusable script snippets.
   - `scripts/compare_r_outputs.R` validates R/Python-derived outputs.

2. The Shiny app in `shiny/proxiome_demo`.
   - This is a separate R application with its own `renv.lock`, `renv/activate.R`, `.Rprofile`, tests, and Open OnDemand launcher template.
   - It is not packaged as a formal R package or golem app. The executable entry point is `app.R`.

### Shiny Entry Points

- `app.R`
  - Resolves `APP_DIR`, detects platform, activates the app-local renv, loads `bslib`, `ggplot2`, `plotly`, and `shiny`, sources helpers, defines UI, defines server, and calls `shinyApp(ui, server)`.
  - It also contains most nontrivial plotting and reactive helper functions below the server.

- `R/data_adapter.R`
  - Loads demo/user RDS data, extracts Seurat metadata/embeddings, computes abundance, clustering, colocalization, QC payloads, and differential readout statistics.

- `R/environment_module.R`
  - Defines the Environment tab. It reports R executable paths, library paths, selected environment variables, and installed package locations.

- `R/plot_layout.R`
  - Defines shared Plotly sizing/margin/colorbar helpers.

- `template/script.sh.erb`
  - Open OnDemand launcher. It detects CCRSF/Biowulf-like hostnames, exports renv variables, checks required packages in the restored app library, then runs `shiny::runApp()`.

### Shiny Data Flow

1. Startup activates renv in `app.R:204-282`.
2. The server creates `demo_data <- reactiveVal(NULL)` in `app.R:921-922`.
3. Startup observer calls `load_default_proxiome_data()` in `app.R:926-940`.
4. Demo/user data loading calls `load_demo_proxiome_data()` in `R/data_adapter.R:79-132`.
5. `build_demo_proxiome_data()` creates a list with:
   - `source`
   - `qc`
   - `metadata`
   - `abundance`
   - `abundance_summary`
   - `clustering`
   - `clustering_summary`
   - `colocalization`
   - `colocalization_summary`
6. A broad observer updates all input choices from the loaded data in `app.R:969-1022`.
7. Output reactives filter subsets and render Plotly plots/tables for QC, Abundance, Clustering, Colocalization, Differential, and Environment tabs.

### Reactive Flow

- Primary app state is centralized in `demo_data`.
- Filters are per-readout: QC, abundance, clustering, and colocalization each have independent inputs.
- Differential reactives call `calculate_differential_readout()` for abundance, clustering, and colocalization in `app.R:1085-1130`.
- Plotly click events update detail selectors for differential views in `app.R:1368-1372`, `app.R:1488-1492`, and `app.R:1636-1640`.
- Most plots are synchronous `renderPlotly()` blocks.

## Prompt Coverage Matrix

| Prompt area | Covered by |
| --- | --- |
| Bugs and potential runtime errors | Findings 8, 13, and 14 |
| Reactive programming issues | Findings 1, 5, and 11 |
| Performance bottlenecks | Findings 1, 3, 4, and 5 |
| Memory inefficiencies | Findings 2 and 3 |
| Security concerns | Finding 10 |
| Maintainability issues | Findings 6, 7, and 11 |
| Scalability limitations | Findings 1, 2, 3, 4, and 13 |
| UX/UI problems | Findings 4, 8, and 9 |
| Shiny best practices | Shiny Best-Practice Evaluation section |
| Technical debt | Findings 6, 7, 12, and 13 |
| Deployment readiness | Findings 10, 12, 13, and 15 |

## Major Findings

### 1. Large user RDS loading is synchronous and blocks the Shiny session

Severity: High

Evidence:

- Startup data loading runs in a normal observer with `withProgress()` in `app.R:926-940`.
- User RDS loading also runs synchronously in the click observer in `app.R:947-967`.
- The load path calls `readRDS()` and expensive Seurat/Pixelator transformations in `R/data_adapter.R:123-124`.
- `pixelatorR::ProximityScores()` is called during load in `R/data_adapter.R:147-151`.

Why this is a problem:

For large Seurat objects, `readRDS()`, abundance extraction, proximity scoring, and QC construction can take long enough to freeze the Shiny R process. On Open OnDemand this can look like a hung app or upstream timeout, especially when users load their own data.

Specific fix:

Move user data loading into a background task and keep the UI responsive. In modern Shiny, use `ExtendedTask` or `future`/`promises` with a bounded worker plan. Keep demo cache synchronous only if it is guaranteed small and already precomputed.

Example direction:

```r
load_user_task <- ExtendedTask$new(function(rds_path, app_dir) {
  future::future({
    load_user_proxiome_data(rds_path, app_dir = app_dir)
  })
})

observeEvent(input$load_rds_path, {
  req(user_rds_path_loading_enabled())
  load_user_task$invoke(input$rds_server_path, APP_DIR)
})

observe({
  result <- load_user_task$result()
  req(result)
  demo_data(result)
})
```

Estimated impact: large improvement in perceived reliability and responsiveness for real user datasets.

### 2. Sparse Seurat assay layers are densified, creating avoidable memory pressure

Severity: High

Evidence:

- `build_abundance_long()` pulls Seurat assay layers in `R/data_adapter.R:336-337`.
- `matrix_to_long()` calls `as.matrix()` in `R/data_adapter.R:352-353`.

Why this is a problem:

Seurat assay layers are commonly sparse matrices. `as.matrix()` converts the selected layer into a dense R matrix. Even with marker subsetting, this can allocate much more memory than needed and can fail for larger objects or broader marker panels.

Specific fix:

Keep sparse matrices sparse and convert only nonzero entries, then explicitly join zeros only where the UI needs them. If a full marker-by-cell table is required, build it in chunks and cache the result.

Example direction:

```r
matrix_to_long_sparse <- function(x, value_name) {
  summary_x <- Matrix::summary(x)
  data.frame(
    marker = rownames(x)[summary_x$i],
    component = colnames(x)[summary_x$j],
    value = summary_x$x,
    stringsAsFactors = FALSE
  ) |>
    stats::setNames(c("marker", "component", value_name))
}
```

Estimated impact: lower peak memory, fewer failed user loads, and better scalability.

### 3. The app computes all Pixelator proximity scores before marker filtering

Severity: High

Evidence:

- `pixelatorR::ProximityScores(object, assay = assay, meta_data_columns = ...)` runs on the whole object in `R/data_adapter.R:147-151`.
- Marker filtering happens only after that call in the following proximity subset.

Why this is a problem:

Proximity scores are the largest readout surface in the app. Computing all marker-pair scores before filtering creates unnecessary CPU and memory cost when the UI only needs a selected marker subset.

Specific fix:

If `pixelatorR::ProximityScores()` supports marker restriction, pass the selected markers into the call. If not, precompute/cache proximity summaries during deployment or create a slim app-specific RDS/cache format that stores only `metadata`, `abundance`, `qc`, and selected proximity summaries.

Estimated impact: shorter load time and lower memory use, especially for user-loaded RDS files.

### 4. Colocalization defaults select all markers and produce O(marker squared) plot data

Severity: High

Evidence:

- Available colocalization markers are derived from all marker pairs in `app.R:979`.
- The update observer selects all of them by default in `app.R:992`.
- The heatmap renderer summarizes and plots all selected marker pairs in `app.R:1559-1598`.
- `prepare_coloc_heatmap_plot_data()` creates a full condition by marker by marker grid in `app.R:2370-2377`.

Why this is a problem:

Marker-pair visualizations scale quadratically. A 30-marker panel creates 870 off-diagonal cells per condition; a 100-marker panel creates 9,900 per condition. Plotly rendering, hover text, and size legends become slow and cluttered.

Specific fix:

Default to a bounded marker subset, such as 12-15 markers, and add explicit controls for "select all" or "top markers by detection". Enforce a hard warning threshold before rendering very large heatmaps.

Example direction:

```r
default_heatmap_markers <- function(markers, max_markers = 15) {
  head(markers, max_markers)
}

updateSelectizeInput(
  session,
  "colocalization_heatmap_markers",
  choices = colocalization_markers,
  selected = default_heatmap_markers(colocalization_markers)
)
```

Estimated impact: faster default rendering and a clearer user experience.

### 5. Differential statistics recompute when only display thresholds change

Severity: Medium

Evidence:

- Differential reactives pass `input$*_diff_fdr` into `calculate_differential_readout()` in `app.R:1085-1130`.
- The function uses `fdr_cutoff` only to set `is_significant` in `R/data_adapter.R:544`.
- Plot/table filtering separately applies thresholds in `app.R:2771-2856`.

Why this is a problem:

Changing the FDR threshold should not rerun Wilcoxon tests. In the current reactive graph, threshold changes invalidate the expensive differential calculation for all readout types.

Specific fix:

Make `calculate_differential_readout()` independent of FDR and effect-size UI thresholds. Compute p-values and adjusted p-values from group/celltype/min-cell inputs only, then apply thresholds downstream in plot and table helpers.

Estimated impact: faster differential interaction and less CPU churn.

### 6. `app.R` is too large and mixes startup, UI, server, statistics, plotting, and deployment helpers

Severity: Medium

Evidence:

- `app.R` is 2,885 lines.
- It contains platform detection, renv activation, CSS, sidebars, UI, server, QC plotting, colocalization heatmap logic, differential plotting, and utility functions.
- Only three helper files exist: `data_adapter.R`, `environment_module.R`, and `plot_layout.R`.

Why this is a problem:

The current file is hard to review and risky to change. Many functions below the server are pure helpers that could be tested independently, but they are coupled to global app startup because tests source `app.R`.

Specific fix:

Split by responsibility:

- `R/runtime.R`: app dir, platform, renv, writable dirs.
- `R/ui_controls.R`: sidebars and navbar controls.
- `R/qc_module.R`: QC UI/server/helpers.
- `R/abundance_module.R`: abundance observed and differential views.
- `R/clustering_module.R`: clustering views.
- `R/colocalization_module.R`: heatmaps and differential colocalization.
- `R/differential.R`: shared statistics and volcano plot helpers.

Estimated impact: much better maintainability and lower regression risk.

### 7. Module boundaries are incomplete because `environment_module.R` depends on globals from `app.R`

Severity: Medium

Evidence:

- `environment_module.R` calls `resolve_app_r_libs()` and `normalize_existing_paths()` in `R/environment_module.R:40` and `R/environment_module.R:54`.
- Those functions are defined in `app.R:256-280`, not in the module file.
- Tests explicitly check that the module is sourced with `local = TRUE` and receives `APP_DIR` in `tests/testthat/test_app_ui.R:210-218`.

Why this is a problem:

The module is not reusable or independently sourceable. It works because of a specific sourcing order in `app.R`, which is fragile during refactoring.

Specific fix:

Move runtime path helpers into `R/runtime.R` and source that before all modules. Then make `environment_module.R` depend only on functions from `runtime.R` and explicit `app_dir`.

Estimated impact: clearer module contract and easier testing.

### 8. User RDS validation is minimal and errors are not diagnostic enough

Severity: Medium

Evidence:

- `validate_rds_file_path()` checks only nonempty path, `.rds` suffix, and existence in `app.R:871-883`.
- `require_namespace()` reports only the package name in `R/data_adapter.R:587-591`.
- RDS loading then calls `readRDS()` without schema validation in `R/data_adapter.R:123`.

Why this is a problem:

Users can select an RDS that is not a compatible Seurat object, lacks the expected assay layers, lacks required reductions, or lacks Pixelator proximity support. The resulting error may surface late and not tell the user what is missing.

Specific fix:

Add a `validate_proxiome_object()` step immediately after `readRDS()`:

```r
validate_proxiome_object <- function(object) {
  require_namespace("Seurat")
  if (!inherits(object, "Seurat")) {
    stop("RDS must contain a Seurat object.", call. = FALSE)
  }
  reductions <- Seurat::Reductions(object)
  if (!any(reductions %in% c("umap", "harmony_umap"))) {
    stop("Seurat object must include a umap or harmony_umap reduction.", call. = FALSE)
  }
  invisible(TRUE)
}
```

Also include `.libPaths()` in package-missing diagnostics for deployment support.

Estimated impact: fewer confusing upload/load failures and faster support.

### 9. Table outputs are static and capped rather than interactive

Severity: Medium

Evidence:

- Metadata table uses `head(qc_origin_metadata(), 200)` in `app.R:1236-1238`.
- Summary/differential tables use `renderTable()` throughout the server.
- The renv lock includes Plotly/Shiny stack but not a user-facing table package such as DT in the app code.

Why this is a problem:

Users cannot search, sort, export, or page through QC metadata and differential hits. For large datasets, static tables are both less useful and more memory-heavy in the browser.

Specific fix:

Use `DT::renderDT()` or `reactable` for all larger tables, with server-side processing for metadata and row limits for very large readouts.

Estimated impact: better UX and lower browser rendering pressure.

### 10. Server-side RDS paths and RDS deserialization need a documented trust boundary

Severity: Medium

Evidence:

- The user data control accepts a server-visible path in `app.R:703-735`.
- Validation checks suffix and existence, but not allowed root directories or object trust level, in `app.R:871-883`.
- The app calls `readRDS()` on the path in `R/data_adapter.R:123`.

Why this is a problem:

On Open OnDemand the Shiny process normally runs as the launching user, which limits filesystem access to that user's permissions. That is a reasonable model, but it should be explicit. Without an allowlist or clear documentation, users can point the app at any readable `.rds` path, and support staff may not know whether path access is intended. R deserialization is also not a safe boundary for untrusted files.

Specific fix:

Document that RDS path loading is same-user trusted input. For shared deployments, optionally enforce an allowlist through `PROXIOME_ALLOWED_RDS_ROOTS`, and reject paths outside those roots before `readRDS()`.

Example direction:

```r
validate_rds_allowed_root <- function(rds_path, allowed_roots = Sys.getenv("PROXIOME_ALLOWED_RDS_ROOTS", "")) {
  roots <- strsplit(allowed_roots, .Platform$path.sep, fixed = TRUE)[[1]]
  roots <- normalizePath(roots[nzchar(roots)], mustWork = FALSE)
  if (length(roots) == 0) {
    return(invisible(TRUE))
  }
  resolved <- normalizePath(rds_path, mustWork = TRUE)
  if (!any(startsWith(resolved, paste0(roots, .Platform$file.sep)))) {
    stop("RDS path is outside the allowed data roots.", call. = FALSE)
  }
  invisible(TRUE)
}
```

Estimated impact: clearer security posture and fewer accidental reads of unsupported files.

### 11. Tests cover many helpers but not full Shiny server behavior

Severity: Medium

Evidence:

- Existing tests verify UI text/controls, helper outputs, heatmap behavior, renv files, and launcher strings.
- There is no `testServer()` coverage for reactive behavior after data load, input updates, Plotly click events, or error notifications.
- There is no visual regression or browser-level test for the UI.

Why this is a problem:

Many prior issues in this app were UI/runtime issues. Static HTML grep tests catch structure, but they do not prove that the server updates controls correctly, plots render after changing inputs, or user data load errors are handled.

Specific fix:

Add `testServer()` tests using small synthetic app data for:

- Demo load and source summary.
- User RDS load success and failure.
- Differential group changes without recomputing thresholds.
- Colocalization anchor click updates.

Add a small number of `shinytest2` or Playwright checks for the main tabs.

Estimated impact: stronger regression protection for user-facing behavior.

### 12. CI does not run Shiny tests

Severity: Medium

Evidence:

- `.github/workflows/ci_mkdocs.yaml` only runs `pixi run build-docs`.
- The workflow paths do not include `shiny/proxiome_demo/**`.
- Pixi defines `test-shiny-proxiome` in `pixi.toml:91-96`, but CI does not call it.

Why this is a problem:

Shiny regressions will not be caught in PRs or pushes. The app has a real test suite, but it is currently local-only.

Specific fix:

Add a CI job or workflow that runs `pixi run -e r test-shiny-proxiome` for `shiny/proxiome_demo/**`, `pixi.toml`, and `pixi.lock` changes. If app-local renv is the source of truth, add an R 4.5 job that restores `shiny/proxiome_demo/renv.lock` and runs the same tests.

Estimated impact: earlier detection of app regressions.

### 13. R version strategy is split between Pixi and the Shiny renv/Open OnDemand deployment

Severity: Medium

Evidence:

- Pixi R feature pins `r-base = ">=4.4,<4.5"` in `pixi.toml:42-43`.
- Open OnDemand launcher defaults to R 4.5.2 for Biowulf in `template/script.sh.erb:14-17`.
- Current HPC verification used R 4.5.0 and the app-local renv.

Why this is a problem:

Developers can run Shiny tests through Pixi under R 4.4.x while production OnDemand runs R 4.5.x. Compiled packages such as Seurat, SeuratObject, Matrix, and pixelatorR are sensitive to R minor versions and library paths.

Specific fix:

Choose one Shiny development path:

- Either update Pixi R to match production R 4.5.x and make `pixi run -e r test-shiny-proxiome` use the app-local renv, or
- Document that Pixi R is only for the reproduction workflow and Shiny development must use `shiny/proxiome_demo/renv.lock` under R 4.5.x.

Estimated impact: fewer local/OnDemand discrepancies.

### 14. Top-level reproduction workflow has hard-coded biological/sample assumptions

Severity: Low

Evidence:

- `sample_metadata_from_name()` infers condition, time, and ratio from filename substrings in `scripts/pixelgen_repro.py:98-132`.
- `annotate_cells()` uses fixed marker heuristics and quantiles in `scripts/pixelgen_repro.py:208-277`.
- `configs/analysis_config.yaml` includes `max_umap_cells: 100000`, but `add_umap()` uses all cells in `scripts/pixelgen_repro.py:280-299`.

Why this is a problem:

This is acceptable for a reproduction repository, but fragile if reused as a general Pixelgen workflow. The unused `max_umap_cells` setting suggests intended scalability behavior that is not implemented.

Specific fix:

Document the workflow as dataset-specific, or move sample parsing and thresholds fully into config. Either implement `max_umap_cells` in `add_umap()` or remove the unused config key.

Estimated impact: clearer reproducibility boundaries and fewer surprises when adapting the workflow.

### 15. Logging is limited to user notifications and process output

Severity: Low

Evidence:

- User-visible errors are mostly `showNotification()` in `app.R:947-967`.
- Startup/deployment errors print to process logs.
- No app-level structured logging or session ID correlation exists.

Why this is a problem:

Support for Open OnDemand sessions depends on reading raw output logs. When a user reports a failure, there is no structured record of data source, app version, package paths, or failing step.

Specific fix:

Add a lightweight logger that writes per-session events under `$PROXIOMEVIS_HOME/runtime`, with sensitive paths redacted or basename-only where appropriate.

Estimated impact: easier support without exposing sensitive data.

## Shiny Best-Practice Evaluation

Strengths:

- Uses `bslib` navigation, sidebars, accordions, and value boxes consistently.
- Has separate controls for readouts and avoids unnecessary colocalization controls in abundance/clustering panels.
- Uses `validate(need(...))` in many render functions for empty states.
- Uses an app-local renv strategy and avoids package installation at startup.
- Includes meaningful unit tests for many helper functions and UI structure.

Gaps:

- App is not modularized into Shiny modules by readout.
- Most expensive operations are synchronous.
- Large tables use static `renderTable()`.
- Reactive graph recomputes statistics more often than needed.
- Environment module depends on globals from `app.R`.
- Testing is helper-heavy and lacks server/browser-level coverage.
- Logging is insufficient for production support.

## Prioritized Recommendations

### Quick Wins (<1 day)

1. Cap default colocalization heatmap markers to 12-15 and require explicit user opt-in for larger marker sets.
   - Impact: high UI responsiveness gain, low implementation risk.

2. Remove FDR threshold from differential statistic reactives.
   - Impact: faster threshold interaction, small code change.

3. Improve user RDS validation and package-missing errors.
   - Impact: fewer confusing failures, easier support.

4. Add CI coverage for `shiny/proxiome_demo/tests/testthat`.
   - Impact: catches app regressions before deployment.

5. Document the R version split or align Pixi and Shiny R versions.
   - Impact: fewer package/library-path incidents.

### Medium Improvements (1-3 days)

1. Replace static tables with DT/reactable equivalents.
   - Impact: better UX for metadata and differential results.

2. Extract runtime and environment helpers into `R/runtime.R`.
   - Impact: cleaner module contracts and easier testing.

3. Add `testServer()` coverage for data load, input updates, and key reactive outputs.
   - Impact: stronger regression coverage for actual Shiny behavior.

4. Make user RDS loading schema-aware.
   - Impact: robust support for user data and clearer failure messages.

5. Implement sparse abundance extraction.
   - Impact: lower peak memory for larger datasets.

### Major Refactors (>3 days)

1. Modularize the app by readout.
   - Impact: substantial maintainability improvement and lower future feature cost.

2. Move expensive user-data processing to async/background tasks.
   - Impact: prevents blocked sessions and improves Open OnDemand reliability.

3. Introduce a slim app cache/data contract.
   - Impact: faster startup and less dependence on full Seurat object recomputation.

4. Add browser-level regression testing.
   - Impact: catches layout, modal, Plotly, and table regressions that unit tests miss.

## Refactoring Roadmap

### Phase 1: Stabilize Deployment And Fast Interactions

- Add CI for Shiny tests.
- Align/document R version strategy.
- Bound colocalization default markers.
- Decouple threshold controls from differential test recomputation.
- Improve RDS schema/package diagnostics.

Expected outcome: fewer deployment surprises and better responsiveness with minimal architecture change.

### Phase 2: Reduce Memory And Runtime Cost

- Replace `as.matrix()` abundance extraction with sparse-aware conversion.
- Cache expensive proximity summaries per source RDS.
- Limit or precompute colocalization summaries.
- Add object schema validation before heavy computations.

Expected outcome: user-loaded RDS files become more practical and less likely to fail.

### Phase 3: Modularize The Shiny App

- Create `R/runtime.R`.
- Move QC, abundance, clustering, colocalization, and differential logic into separate files or Shiny modules.
- Keep pure functions testable without sourcing the whole app.
- Add `testServer()` tests for each module.

Expected outcome: easier development and safer changes.

### Phase 4: Productionize For Shared HPC Use

- Add async loading for user RDS paths.
- Add session-level structured logging under `$HOME/.ProxiomeVis/runtime`.
- Add browser-level smoke tests for Open OnDemand-compatible workflows.
- Create a deployment checklist for CCRSF, Biowulf, and portable desktop use.

Expected outcome: robust multi-HPC operation with better supportability.

## Estimated Impact Summary

| Recommendation | Severity Addressed | Effort | Expected Impact |
| --- | --- | --- | --- |
| Async user RDS loading | High | Major | Prevents blocked sessions and proxy-like failures for large data |
| Sparse abundance conversion | High | Medium | Reduces peak memory and supports larger Seurat objects |
| Pre-filter/cache proximity scores | High | Medium/Major | Shortens load time and reduces memory |
| Cap heatmap markers by default | High | Quick | Improves first-render speed and readability |
| Decouple differential thresholds | Medium | Quick | Faster volcano/table threshold changes |
| Split `app.R` into modules | Medium | Major | Improves maintainability and testability |
| Move runtime helpers out of `app.R` | Medium | Medium | Makes modules sourceable and contracts explicit |
| Add RDS schema validation | Medium | Quick/Medium | Clearer user errors and less support time |
| Document/enforce RDS path trust boundary | Medium | Quick/Medium | Clarifies security posture for shared deployments |
| Use interactive tables | Medium | Medium | Better UX for metadata and result exploration |
| Add Shiny CI | Medium | Quick | Prevents untested app regressions |
| Align R versions | Medium | Quick/Medium | Reduces package/library mismatch failures |
| Add session logging | Low | Medium | Better production support |

## Suggested Next Implementation Order

1. Add CI for Shiny tests and R parse checks.
2. Cap default colocalization heatmap markers.
3. Decouple differential thresholds from statistical recomputation.
4. Add user RDS object validation and richer package diagnostics.
5. Extract runtime helpers into `R/runtime.R`.
6. Convert abundance extraction to sparse-aware code.
7. Add async user RDS loading.
8. Split the app into readout modules.

This order prioritizes lower-risk changes that reduce the current support burden before larger architectural refactors.
