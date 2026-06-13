source("../../R/spatial_metrics.R")

test_that("sample-level spatial heatmap summary uses stored log2 ratios", {
  proximity <- data.frame(
    component = c("cell-a", "cell-b", "cell-c", "cell-a", "cell-c"),
    sample_alias = c("S1", "S1", "S2", "S1", "S2"),
    condition = c("UNT", "UNT", "PHA", "UNT", "PHA"),
    marker_1 = c("CD3e", "CD3e", "CD3e", "CD4", "CD4"),
    marker_2 = c("CD4", "CD4", "CD4", "CD8", "CD8"),
    log2_ratio = c(1, 3, -1, 0.5, 1.5),
    stringsAsFactors = FALSE
  )

  result <- summarize_spatial_heatmap_by_sample(
    proximity,
    selected_markers = c("CD3e", "CD4", "CD8")
  )

  cd3_cd4_s1 <- result[
    result$sample_alias == "S1" &
      result$marker_1 == "CD3e" &
      result$marker_2 == "CD4",
    ,
    drop = FALSE
  ]

  expect_equal(cd3_cd4_s1$mean_log2_ratio, 2)
  expect_equal(cd3_cd4_s1$n_detected, 2L)
  expect_equal(cd3_cd4_s1$n_total, 2L)
  expect_equal(cd3_cd4_s1$pct_detected, 1)
})

test_that("per-cell-type spatial summary groups by sample and cell type", {
  proximity <- data.frame(
    component = c("cell-a", "cell-b", "cell-c", "cell-d"),
    sample_alias = c("S1", "S1", "S1", "S2"),
    condition = c("UNT", "UNT", "UNT", "PHA"),
    celltype_manual = c("CD8 T", "CD8 T", "Raji", "CD8 T"),
    marker_1 = c("CD3e", "CD3e", "CD19", "CD3e"),
    marker_2 = c("CD4", "CD4", "CD20", "CD4"),
    log2_ratio = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  result <- summarize_spatial_heatmap_by_celltype(
    proximity,
    selected_markers = c("CD3e", "CD4", "CD19", "CD20")
  )

  expect_true(all(c("sample_alias", "condition", "celltype_manual") %in% names(result)))
  expect_equal(
    result$mean_log2_ratio[result$sample_alias == "S1" & result$celltype_manual == "CD8 T"],
    1.5
  )
})

test_that("spatial heatmap completion mirrors observed marker pairs before filling missing pairs", {
  summary <- data.frame(
    sample_alias = "S1",
    condition = "UNT",
    marker_1 = "CD3e",
    marker_2 = "CD4",
    mean_log2_ratio = 1.2,
    pct_detected = 0.5,
    n_detected = 3L,
    n_total = 6L,
    stringsAsFactors = FALSE
  )

  completed <- complete_spatial_marker_pairs(
    summary,
    selected_markers = c("CD3e", "CD4", "CD8"),
    group_cols = c("sample_alias", "condition")
  )

  reverse_pair <- completed[completed$marker_1 == "CD4" & completed$marker_2 == "CD3e", , drop = FALSE]
  expect_equal(nrow(reverse_pair), 1L)
  expect_equal(reverse_pair$mean_log2_ratio, 1.2)
  expect_equal(reverse_pair$pct_detected, 0.5)
  expect_equal(reverse_pair$n_detected, 3L)
  expect_equal(reverse_pair$n_total, 6L)
  expect_true(any(completed$marker_1 == "CD3e" & completed$marker_2 == "CD8"))
  expect_equal(
    completed$mean_log2_ratio[completed$marker_1 == "CD3e" & completed$marker_2 == "CD8"],
    0
  )
})

test_that("spatial marker selection ranks variable detected markers", {
  summary <- data.frame(
    sample_alias = rep(c("S1", "S2"), each = 4),
    condition = rep(c("UNT", "PHA"), each = 4),
    marker_1 = rep(c("CD3e", "CD3e", "CD19", "CD19"), times = 2),
    marker_2 = rep(c("CD4", "CD8", "CD20", "CD21"), times = 2),
    mean_log2_ratio = c(0.1, 0.2, 1.5, 1.6, 0.1, 0.3, -1.2, -1.1),
    pct_detected = c(0.9, 0.9, 0.8, 0.8, 0.9, 0.9, 0.8, 0.8),
    stringsAsFactors = FALSE
  )

  selected <- select_spatial_heatmap_markers(
    summary,
    available_markers = c("CD3e", "CD4", "CD8", "CD19", "CD20", "CD21"),
    n_markers = 3,
    min_pct_detected = 0.25,
    min_range = 0.2
  )

  expect_true("CD19" %in% selected)
  expect_lte(length(selected), 3)
})

test_that("clustering heatmap summary uses stored self-proximity scores", {
  clustering <- data.frame(
    component = c("cell-a", "cell-b", "cell-c", "cell-d", "cell-e"),
    marker = c("CD3e", "CD3e", "CD4", "CD4", "CD8"),
    condition = c("UNT", "UNT", "PHA", "PHA", "PHA"),
    celltype_manual = c("CD8 T", "CD8 T", "CD4 T", "CD4 T", "CD8 T"),
    log2_ratio = c(1, 3, -1, 1, 0.5),
    stringsAsFactors = FALSE
  )

  summary <- summarize_clustering_heatmap(
    clustering,
    selected_markers = c("CD3e", "CD4", "CD8")
  )

  expect_true(all(c("marker", "condition", "celltype_manual", "mean_log2_ratio", "n_cells") %in% names(summary)))
  expect_equal(
    summary$mean_log2_ratio[summary$marker == "CD3e" & summary$condition == "UNT"],
    2
  )
  expect_equal(
    summary$n_cells[summary$marker == "CD4" & summary$condition == "PHA"],
    2L
  )

  selected <- select_clustering_heatmap_markers(summary, n_markers = 2)
  expect_lte(length(selected), 2)
  expect_true("CD3e" %in% selected)
})

test_that("spatial heatmap summaries avoid base R split lapply and merge hot paths", {
  spatial_source <- paste(deparse(body(summarize_spatial_heatmap)), collapse = "\n")
  clustering_source <- paste(deparse(body(summarize_clustering_heatmap)), collapse = "\n")

  for (source_text in c(spatial_source, clustering_source)) {
    expect_false(grepl("\\bsplit\\s*\\(", source_text))
    expect_false(grepl("\\blapply\\s*\\(", source_text))
    expect_false(grepl("\\bmerge\\s*\\(", source_text))
  }
})
