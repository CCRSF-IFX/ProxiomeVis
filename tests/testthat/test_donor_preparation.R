source(file.path("..", "..", "R", "bundle_builder.R"))

test_that("import appends donor results and keeps existing within-donor statistics", {
  root <- tempfile("donor-import-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  data <- comparison_view_fixture()
  data$contract <- precomputed_app_data(precomputed_fixture())$contract
  data$contract$proximity_source <- NULL
  data$contract$donor_contrasts$abundance$p_value <- .02
  data$contract$donor_contrasts$abundance$p_adj <- .04
  source <- file.path(root, "analysis.rds")
  write_precomputed_bundle(data, source)
  data <- read_precomputed_rds(source)
  entry <- comparison_view_entry(data)
  exported <- file.path(root, "comparison.rds")
  saveRDS(list(result_version = 1L, entry = entry), exported)
  output <- import_user_comparison(source, exported, file.path(root, "with_results.rds"))
  imported <- read_precomputed_rds(output)
  expect_equal(imported$contract$schema_version, 2L)
  expect_identical(imported$contract$sample_summaries, data$contract$sample_summaries)
  old <- imported$contract$donor_contrasts$abundance
  old <- old[old$analysis_id != entry$id, , drop = FALSE]
  expect_equal(old$effect_size, data$contract$donor_contrasts$abundance$effect_size)
  expect_equal(old$p_adj, data$contract$donor_contrasts$abundance$p_adj)
  expect_true(all(old$design == "within_donor"))
  expect_true(entry$id %in% imported$contract$differential_results$abundance$analysis_id)
  expect_identical(imported$contract$comparison_provenance[[entry$id]]$request, entry$request)
  scope <- precomputed_explorer_scope(imported, "abundance", PRECOMPUTED_ALL, "B", "TRT vs NC [new-1]")
  expect_identical(scope$analysis_id, entry$id)
  expect_identical(scope$rows$p_adj, entry$result$results$p_adj)
  within <- precomputed_explorer_scope(imported, "abundance", "D1", "B")
  expect_true(all(within$rows$p_adj == .04))
  second <- comparison_view_entry(imported, id = "new-2")
  saveRDS(list(result_version = 1L, entry = second), exported)
  next_output <- import_user_comparison(output, exported, file.path(root, "two_results.rds"))
  again <- read_precomputed_rds(next_output)
  expect_setequal(names(again$contract$comparison_provenance), c("new-1", "new-2"))
  expect_true(all(c("new-1", "new-2") %in% again$contract$differential_results$abundance$analysis_id))
  duplicate <- comparison_view_entry(again, id = "new-2")
  saveRDS(list(result_version = 1L, entry = duplicate), exported)
  expect_error(import_user_comparison(next_output, exported, file.path(root, "duplicate.rds")), "already stored")
  expect_false(file.exists(file.path(root, "duplicate.rds")))
})

test_that("the offline command uses the same engine for every donor comparison design", {
  root <- tempfile("donor-cli-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  script <- normalizePath(file.path("..", "..", "scripts", "run_comparison.R"))
  previous <- setwd(root)
  on.exit(setwd(previous), add = TRUE)
  for (scenario in c("paired-condition", "paired-celltype", "independent-condition")) {
    independent <- startsWith(scenario, "independent")
    axis <- if (endsWith(scenario, "celltype")) "celltype" else "condition"
    data <- comparison_view_fixture(independent)
    data$counts <- NULL
    bundle <- file.path(root, paste0(scenario, "-analysis.rds"))
    write_precomputed_bundle(data, bundle)
    data <- read_precomputed_rds(bundle)
    expected <- comparison_view_entry(data, axis, independent, id = scenario)
    request <- expected$request
    request[c("group_a_cells", "group_b_cells", "dataset_id", "bundle_md5")] <- NULL
    input <- file.path(root, paste0(scenario, ".json"))
    jsonlite::write_json(request, input, auto_unbox = TRUE)
    comparison <- file.path(root, paste0(scenario, ".rds"))
    prepared <- file.path(root, paste0(scenario, "-prepared.rds"))
    command <- c(script, bundle, input, comparison, "--id", scenario, "--output-bundle", prepared)
    output <- system2(file.path(R.home("bin"), "Rscript"), shQuote(command), stdout = TRUE, stderr = TRUE)
    expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
    expect_true(file.exists(comparison))
    expect_true(file.exists(prepared))
    saved <- readRDS(comparison)
    expect_identical(saved$result_version, 1L)
    expect_equal(saved$entry$result$results, expected$result$results)
    loaded <- read_precomputed_rds(prepared)
    expect_identical(loaded$contract$comparison_provenance[[scenario]]$id, scenario)
    expect_identical(loaded$contract$differential_results$abundance$p_adj, expected$result$results$p_adj)
  }
})
