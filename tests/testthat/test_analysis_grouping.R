source("../../R/data_source_module.R")

grouping_test_metadata <- function() {
  data.frame(
    component = paste0("cell-", 1:4),
    sample_alias = rep(c("sample-a", "sample-b"), each = 2),
    sample = rep(c("raw-a", "raw-b"), each = 2),
    condition = rep(c("control", "treated"), each = 2),
    donor = rep(c("donor-1", "donor-2"), each = 2),
    celltype_manual = c("T", "B", "T", "B"),
    qc_score = 1:4,
    stringsAsFactors = FALSE
  )
}

test_that("analysis grouping offers metadata columns that are constant within samples", {
  metadata <- grouping_test_metadata()

  expect_setequal(
    sample_level_grouping_columns(metadata),
    c("sample_alias", "sample", "condition", "donor")
  )
  expect_false("celltype_manual" %in% sample_level_grouping_columns(metadata))
  expect_false("qc_score" %in% sample_level_grouping_columns(metadata))

  config <- new_analysis_grouping_config(metadata)
  expect_equal(config$column, "condition")
  expect_equal(config$map$analysis_group, c("control", "treated"))
})

test_that("custom analysis grouping updates sample-backed app data", {
  metadata <- grouping_test_metadata()
  config <- new_analysis_grouping_config(metadata)
  config <- update_analysis_grouping_config(
    config,
    mode = "custom",
    custom_groups = c("sample-a" = "baseline", "sample-b" = "stimulated")
  )
  data <- list(
    metadata = metadata,
    clustering = data.frame(
      component = c("cell-1", "cell-3"),
      sample_alias = c("sample-a", "sample-b"),
      condition = c("control", "treated"),
      stringsAsFactors = FALSE
    ),
    colocalization = data.frame(
      component = c("cell-2", "cell-4"),
      sample_alias = c("sample-a", "sample-b"),
      condition = c("control", "treated"),
      stringsAsFactors = FALSE
    ),
    qc = list(filter_counts = data.frame(
      sample = c("raw-a", "raw-b", "TOTAL"),
      condition = c("control", "treated", "TOTAL"),
      stringsAsFactors = FALSE
    )),
    source = list()
  )

  grouped <- apply_analysis_grouping(data, config)

  expect_equal(grouped$metadata$condition, c("baseline", "baseline", "stimulated", "stimulated"))
  expect_equal(grouped$clustering$condition, c("baseline", "stimulated"))
  expect_equal(grouped$colocalization$condition, c("baseline", "stimulated"))
  expect_equal(grouped$qc$filter_counts$condition, c("baseline", "stimulated", "TOTAL"))
  expect_equal(grouped$source$analysis_group_label, "Custom sample groups")
  expect_equal(grouped$source$analysis_group_count, 2L)
  expect_equal(config$source$condition, c("control", "treated"))
})

test_that("analysis grouping validates custom labels and can reset to source condition", {
  config <- new_analysis_grouping_config(grouping_test_metadata())

  expect_error(
    update_analysis_grouping_config(config, mode = "custom", custom_groups = c("", "treated")),
    "cannot be blank"
  )

  custom <- update_analysis_grouping_config(config, mode = "custom", custom_groups = c("A", "B"))
  reset <- update_analysis_grouping_config(custom, mode = "column", column = "condition")
  expect_equal(reset$map$analysis_group, c("control", "treated"))
})
