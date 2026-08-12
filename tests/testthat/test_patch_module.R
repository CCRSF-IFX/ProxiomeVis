source("../../R/data_adapter.R")

test_that("Raji demo cache path is separate from the default demo cache", {
  app_dir <- normalizePath("../..", mustWork = TRUE)

  expect_equal(basename(raji_demo_cache_path(app_dir)), "raji_cart_demo_data.rds")
  expect_false(identical(
    normalizePath(default_demo_cache_path(app_dir), mustWork = FALSE),
    normalizePath(raji_demo_cache_path(app_dir), mustWork = FALSE)
  ))
})

test_that("Raji patch demo payload keeps marker, plan, and sanity-check tables", {
  patch_dir <- tempfile("patch-demo-")
  dir.create(file.path(patch_dir, "tables"), recursive = TRUE)
  write.csv(
    data.frame(run_patch_detection = FALSE, n_patch_markers = 2, stringsAsFactors = FALSE),
    file.path(patch_dir, "tables", "patch_detection_run_plan.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(marker = c("CD19", "CD8"), label = c("marker for Raji_coculture", "marker for CD8T_coculture")),
    file.path(patch_dir, "tables", "patch_marker_unmixing_table.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(condition = "coculture", celltype_condition = "CD8T_coculture", n_cells = 1, median_raji_score = 1.2),
    file.path(patch_dir, "tables", "raji_marker_abundance_sanity_check.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(component = "cell-a", celltype_condition = "CD8T_coculture", raji_marker_count = 10, log2_ratio = 0.8),
    file.path(patch_dir, "tables", "raji_marker_joint_proximity_sanity_check.csv"),
    row.names = FALSE
  )

  patch <- build_raji_patch_demo_data(patch_dir)

  expect_named(patch, c("run_plan", "marker_unmixing", "raji_marker_abundance", "raji_marker_proximity", "patch_burden"))
  expect_equal(patch$run_plan$n_patch_markers, 2)
  expect_equal(patch$marker_unmixing$marker, c("CD19", "CD8"))
  expect_equal(patch$raji_marker_proximity$log2_ratio, 0.8)
  expect_null(patch$patch_burden)
})
