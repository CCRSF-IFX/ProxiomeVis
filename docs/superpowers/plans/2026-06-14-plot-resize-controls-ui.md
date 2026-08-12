# Plot Resize Controls UI Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make per-plot resize controls feel like plot tools, not raw sidebar controls moved into the chart area.

**Architecture:** Keep resize controls inside `plot_pane()`, next to downloads, but collapse width and height into a compact `Size` popover. Reuse the existing input IDs so current reactive plot sizing and download sizing do not change.

**Tech Stack:** Shiny, bslib `popover()`, existing CSS in `app.R`, testthat UI source tests.

---

### Task 1: Capture The Current Problem

**Files:**
- Reference screenshots:
  - `screenshots/proxiomevis-abundance-observed-resize.png`
  - `screenshots/proxiomevis-abundance-marker-distributions-resize.png`

- [ ] **Step 1: Review screenshots**

Confirm the current issue:
- `Plot width` and `Plot height` look like ordinary form fields floating beside the plot.
- Labels create visual noise in the plot header.
- Marker Distributions has the controls sitting awkwardly above the plot, competing with downloads.

### Task 2: Add A Failing UI Test

**Files:**
- Modify: `tests/testthat/test_app_ui.R`

- [ ] **Step 1: Replace the current resize-control placement expectation**

Update `resizable plots keep size controls with the plot pane controls` so it expects a compact popover:

```r
expect_true(grepl("plot_size_controls <- function", plot_layout_source, fixed = TRUE))
expect_true(grepl("plot-size-button", app_source, fixed = TRUE))
expect_true(grepl("plot-size-popover", app_source, fixed = TRUE))
expect_true(grepl('controls = plot_size_controls(ns, "abundance_umap_width", "abundance_umap_height"', abundance_module_source, fixed = TRUE))
expect_true(grepl('controls = plot_size_controls(ns, "abundance_distribution_width", "abundance_distribution_height"', abundance_module_source, fixed = TRUE))
expect_false(grepl('class = "plot-resize-controls"', plot_layout_source, fixed = TRUE))
```

- [ ] **Step 2: Run the focused UI test and confirm it fails**

Run:

```bash
/mnt/ccrsf-static/illumina/RnD_pixelgen_CAR-T_datasets/.pixi/envs/r/bin/Rscript -e "Sys.setenv(PROXIOMEVIS_HOME='/tmp/xies4/proxiomevis-test', PROXIOME_RUNTIME_DIR='/tmp/xies4/proxiomevis-test/runtime', XDG_CACHE_HOME='/tmp/xies4/.cache'); setwd('shiny/proxiome_demo'); testthat::test_file('tests/testthat/test_app_ui.R')"
```

Expected: FAIL because `plot_size_controls()` and the popover classes do not exist yet.

### Task 3: Replace Raw Inline Inputs With A Size Popover

**Files:**
- Modify: `R/plot_layout.R`
- Modify: `R/abundance_module.R`

- [ ] **Step 1: Rename and restyle the helper**

Replace `plot_resize_controls()` with `plot_size_controls()`:

```r
plot_size_controls <- function(
  ns,
  width_id,
  height_id,
  width_value,
  height_value,
  min_width = 420,
  min_height = 320,
  max_value = 2600,
  step = 50
) {
  if (is.null(ns)) {
    ns <- identity
  }

  bslib::popover(
    actionButton(
      ns(paste0(width_id, "_size")),
      "Size",
      class = "btn btn-outline-secondary btn-sm plot-size-button"
    ),
    div(
      class = "plot-size-popover",
      div(
        class = "plot-size-field",
        numericInput(ns(width_id), "Width", value = width_value, min = min_width, max = max_value, step = step),
        span("px", class = "plot-size-unit")
      ),
      div(
        class = "plot-size-field",
        numericInput(ns(height_id), "Height", value = height_value, min = min_height, max = max_value, step = step),
        span("px", class = "plot-size-unit")
      )
    ),
    title = "Plot size",
    placement = "bottom"
  )
}
```

- [ ] **Step 2: Update Abundance plot panes**

Change both Abundance calls:

```r
controls = plot_size_controls(ns, "abundance_umap_width", "abundance_umap_height", width_value = 832, height_value = 520)
```

```r
controls = plot_size_controls(ns, "abundance_distribution_width", "abundance_distribution_height", width_value = 832, height_value = 678)
```

### Task 4: Add Minimal CSS Polish

**Files:**
- Modify: `app.R`

- [ ] **Step 1: Replace resize-input CSS with popover CSS**

Remove `.plot-resize-controls` rules and add:

```css
.plot-size-button {
  min-width: 58px;
  padding: 3px 10px;
  font-size: 0.78rem;
  line-height: 1.35;
}

.plot-size-popover {
  display: grid;
  gap: 10px;
  min-width: 180px;
}

.plot-size-field {
  display: grid;
  grid-template-columns: minmax(0, 1fr) auto;
  align-items: end;
  gap: 6px;
}

.plot-size-field .shiny-input-container {
  width: 100%;
  margin-bottom: 0;
}

.plot-size-field label {
  margin-bottom: 2px;
  color: var(--muted);
  font-size: 0.75rem;
}

.plot-size-field input.form-control {
  min-height: 32px;
  padding: 3px 8px;
  font-size: 0.82rem;
}

.plot-size-unit {
  padding-bottom: 6px;
  color: var(--muted);
  font-size: 0.75rem;
}
```

### Task 5: Verify

**Files:**
- Test: `tests/testthat/test_app_ui.R`

- [ ] **Step 1: Run focused UI tests**

Run:

```bash
/mnt/ccrsf-static/illumina/RnD_pixelgen_CAR-T_datasets/.pixi/envs/r/bin/Rscript -e "Sys.setenv(PROXIOMEVIS_HOME='/tmp/xies4/proxiomevis-test', PROXIOME_RUNTIME_DIR='/tmp/xies4/proxiomevis-test/runtime', XDG_CACHE_HOME='/tmp/xies4/.cache'); setwd('shiny/proxiome_demo'); testthat::test_file('tests/testthat/test_app_ui.R')"
```

Expected: PASS.

- [ ] **Step 2: Run full tests**

Run:

```bash
/mnt/ccrsf-static/illumina/RnD_pixelgen_CAR-T_datasets/.pixi/envs/r/bin/Rscript -e "Sys.setenv(PROXIOMEVIS_HOME='/tmp/xies4/proxiomevis-test', PROXIOME_RUNTIME_DIR='/tmp/xies4/proxiomevis-test/runtime', XDG_CACHE_HOME='/tmp/xies4/.cache'); setwd('shiny/proxiome_demo'); testthat::test_dir('tests/testthat')"
```

Expected: PASS with no failures or warnings.

### Deferred

Skip `Auto / Custom` for this pass. It needs real behavior: auto means computed dimensions ignore manual values, custom means user values override. Add it only if users need to reset plots frequently after changing marker, facet, sample, or cell type counts.
