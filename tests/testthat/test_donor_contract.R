donor_contract_fixture <- function(design = "paired", axis = "condition") {
  data <- precomputed_app_data(precomputed_fixture())
  data$contract$schema_version <- 2L
  data$contract$donor_values <- list()
  for (readout in PROXIOMEVIS_READOUTS) {
    rows <- data$contract$sample_summaries[[readout]]
    rows <- rows[if (axis == "condition") rows$celltype == "B" else rows$condition == "NC", , drop = FALSE]
    rows$source_celltype <- rows$celltype
    rows$celltype <- "B"
    rows$comparison_axis <- axis
    rows$design <- design
    rows$fixed_condition <- if (axis == "celltype") "NC" else NA_character_
    rows$group <- if (axis == "celltype") rows$source_celltype else rows$condition
    rows$group_a <- if (axis == "celltype") "B" else "TRT"
    rows$group_b <- if (axis == "celltype") "T" else "NC"
    rows$comparison <- paste(rows$group_a, "vs", rows$group_b)
    rows$analysis_id <- paste(readout, axis, design, sep = ":")
    rows$eligible <- TRUE
    if (design == "independent") rows <- rows[(rows$donor == "D1" & rows$group == rows$group_a) |
      (rows$donor == "D2" & rows$group == rows$group_b), , drop = FALSE]
    data$contract$donor_values[[readout]] <- rows
    contrasts <- data$contract$donor_contrasts[[readout]]
    contrasts <- contrasts[contrasts$celltype == "B", , drop = FALSE]
    contrasts$comparison_axis <- axis
    contrasts$design <- design
    contrasts$fixed_condition <- if (axis == "celltype") "NC" else NA_character_
    contrasts$group_a <- rows$group_a[1]
    contrasts$group_b <- rows$group_b[1]
    contrasts$comparison <- rows$comparison[1]
    contrasts$analysis_id <- rows$analysis_id[1]
    key <- paste(rows$donor, rows$feature_id, rows$group)
    at_a <- match(paste(contrasts$donor, contrasts$feature_id, contrasts$group_a), key)
    at_b <- match(paste(contrasts$donor, contrasts$feature_id, contrasts$group_b), key)
    contrasts$value_a <- rows$value[at_a]
    contrasts$value_b <- rows$value[at_b]
    contrasts$effect_size <- contrasts$value_a - contrasts$value_b
    if (design == "independent") contrasts <- contrasts[FALSE, ]
    data$contract$donor_contrasts[[readout]] <- contrasts
    result <- data$contract$differential_results[[readout]]
    result <- result[result$celltype == "B", , drop = FALSE]
    result$comparison_axis <- axis
    result$design <- design
    result$fixed_condition <- rows$fixed_condition[1]
    result$group_a <- rows$group_a[1]
    result$group_b <- rows$group_b[1]
    result$comparison <- rows$comparison[1]
    result$analysis_id <- rows$analysis_id[1]
    result$effect_size <- if (axis == "celltype") 0 else if (design == "independent") 1 else 1.5
    result$min_effect <- if (axis == "celltype") 0 else if (design == "independent") NA_real_ else 1
    result$max_effect <- if (axis == "celltype") 0 else if (design == "independent") NA_real_ else 2
    result$n_pairs <- if (design == "independent") 0L else 2L
    result$n_donors_a <- result$n_donors_b <- if (design == "independent") 1L else 2L
    data$contract$differential_results[[readout]] <- result
  }
  data
}

test_that("schema 1 is unchanged and schema 2 supports all donor designs", {
  old <- precomputed_app_data(precomputed_fixture())
  expect_true(validate_precomputed_data(old))
  expect_identical(old$contract$schema_version, 1L)
  upgraded <- old
  upgraded$contract$schema_version <- 2L
  expect_true(validate_precomputed_data(upgraded))
  expect_identical(upgraded$contract$donor_contrasts, old$contract$donor_contrasts)
  expect_true(validate_precomputed_data(donor_contract_fixture()))
  expect_true(validate_precomputed_data(donor_contract_fixture(axis = "celltype")))
  expect_true(validate_precomputed_data(donor_contract_fixture(design = "independent")))
})

test_that("new comparisons validate their axis and fixed population", {
  data <- donor_contract_fixture(axis = "celltype")
  data$contract$donor_values$abundance$fixed_condition <- NA_character_
  expect_error(validate_precomputed_data(data), "fixed condition")
  data <- donor_contract_fixture(axis = "celltype")
  data$contract$donor_values$abundance$source_celltype[1] <- "absent"
  expect_error(validate_precomputed_data(data), "comparison groups")
  data <- donor_contract_fixture()
  data$contract$differential_results$abundance$group_a <- "absent"
  expect_error(validate_precomputed_data(data), "invalid contrast")
  data <- donor_contract_fixture()
  data$contract$schema_version <- 1L
  expect_error(validate_precomputed_data(data), "schema_version = 2")
  data <- donor_contract_fixture(design = "independent", axis = "celltype")
  expect_error(validate_precomputed_data(data), "require paired donors")
})

test_that("donor measurements retain missing values and enforce support", {
  data <- donor_contract_fixture()
  rows <- data$contract$donor_values$abundance
  data$contract$donor_values$abundance <- rbind(rows, rows[1, ])
  expect_error(validate_precomputed_data(data), "duplicate keys")
  data$contract$donor_values$abundance <- rows
  data$contract$donor_values$abundance$n_total_cells[1] <- 100L
  expect_error(validate_precomputed_data(data), "support exceeds")
  data$contract$donor_values$abundance <- rows
  data$contract$donor_values$abundance$n_scored_cells[1] <- 0L
  expect_error(validate_precomputed_data(data), "eligibility")
  data$contract$donor_values$abundance$eligible[1] <- FALSE
  data$contract$donor_values$abundance$value[1] <- 0
  expect_error(validate_precomputed_data(data), "missing measurements must stay NA")
  data$contract$donor_values$abundance$value[1] <- NA_real_
  row <- data$contract$donor_values$abundance[1, ]
  result <- data$contract$differential_results$abundance
  at <- match(row$feature_id, result$feature_id)
  field <- if (row$group == row$group_a) "n_donors_a" else "n_donors_b"
  result[[field]][at] <- result[[field]][at] - 1L
  result$n_pairs[at] <- result$n_pairs[at] - 1L
  data$contract$differential_results$abundance <- result
  paired <- data$contract$donor_contrasts$abundance
  data$contract$donor_contrasts$abundance <- paired[!(paired$donor == row$donor &
    paired$feature_id == row$feature_id), , drop = FALSE]
  expect_true(validate_precomputed_data(data))
})

test_that("independent comparisons cannot contain pairing or repeated donors", {
  data <- donor_contract_fixture(design = "independent")
  data$contract$differential_results$abundance$n_pairs <- 1L
  expect_error(validate_precomputed_data(data), "donor support")
  data <- donor_contract_fixture(design = "independent")
  data$contract$differential_results$abundance$min_effect <- 0
  expect_error(validate_precomputed_data(data), "donor-effect ranges")
  data <- donor_contract_fixture(design = "independent")
  data$contract$donor_values$abundance$donor <- "D1"
  expect_error(validate_precomputed_data(data), "cannot share donors")
  data <- donor_contract_fixture()
  data$contract$differential_results$abundance$n_donors_b <- 1L
  expect_error(validate_precomputed_data(data), "donor support")
  data <- donor_contract_fixture()
  data$contract$differential_results$abundance$p_value <- 2
  expect_error(validate_precomputed_data(data), "invalid probabilities")
  data <- donor_contract_fixture()
  data$contract$differential_results$abundance$replicate_unit <- "cell"
  expect_error(validate_precomputed_data(data), "donor replication")
})

test_that("result context and donor support agree across stored sections", {
  data <- donor_contract_fixture()
  data$contract$differential_results$abundance$n_donors_b <- 3L
  expect_error(validate_precomputed_data(data), "do not match eligible")
  data <- donor_contract_fixture()
  data$contract$donor_contrasts$abundance$value_a[1] <- 10
  data$contract$donor_contrasts$abundance$effect_size[1] <- 9
  expect_error(validate_precomputed_data(data), "paired contrasts do not match")
  data <- donor_contract_fixture()
  data$contract$donor_contrasts$abundance$group_a <- "NC"
  data$contract$donor_contrasts$abundance$group_b <- "TRT"
  expect_error(validate_precomputed_data(data), "context across result tables")
  data <- donor_contract_fixture()
  data$contract$differential_results$abundance$inference <- "donor_level"
  expect_error(validate_precomputed_data(data), "valid inference labels")
})

test_that("within-donor cell-type exploration is accepted without donor-level inference", {
  data <- donor_contract_fixture(axis = "celltype")
  data$contract$donor_values <- NULL
  for (readout in PROXIOMEVIS_READOUTS) {
    data$contract$donor_contrasts[[readout]]$design <- "within_donor"
    data$contract$donor_contrasts[[readout]]$inference <- "exploratory_component_level"
    data$contract$donor_contrasts[[readout]]$p_adj <- 0.01
    data$contract$differential_results[[readout]] <- data$contract$differential_results[[readout]][FALSE, ]
  }
  expect_true(validate_precomputed_data(data))
})

test_that("analysis, session adaptation and prepared RDS loading share the same scenario contract", {
  skip_if_not_installed("pixelatorR")
  for (scenario in c("within_condition", "within_celltype", "paired_condition", "paired_celltype", "independent_condition")) {
    data <- precomputed_app_data(precomputed_fixture())
    source <- rep(data$metadata$component, each = 2L)
    meta <- data$metadata[match(source, data$metadata$component), , drop = FALSE]
    meta$component <- rownames(meta) <- paste0("scenario-cell-", seq_len(nrow(meta)))
    donor <- match(meta$donor, c("D1", "D2"))
    value <- donor * 10 + (meta$condition == "TRT") * c(1, 3)[donor] +
      (meta$celltype == "T") * c(2, 5)[donor]
    data$abundance <- Matrix::Matrix(rbind(A = value, B = value / 2, C = value * 2), sparse = TRUE)
    colnames(data$abundance) <- meta$component
    data$counts <- round(data$abundance * 2)
    data$embeddings <- lapply(data$embeddings, function(x) {
      x <- x[source, , drop = FALSE]
      rownames(x) <- meta$component
      x
    })
    if (startsWith(scenario, "independent")) {
      meta$donor <- meta$sample_id
      for (readout in PROXIOMEVIS_READOUTS) {
        data$contract$sample_summaries[[readout]]$donor <- data$contract$sample_summaries[[readout]]$sample_id
        data$contract$donor_contrasts[[readout]] <- data$contract$donor_contrasts[[readout]][FALSE, ]
        data$contract$differential_results[[readout]] <- data$contract$differential_results[[readout]][FALSE, ]
      }
    }
    data$metadata <- meta
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)
    write_precomputed_bundle(data, path)
    loaded <- read_precomputed_rds(path)
    within <- startsWith(scenario, "within")
    celltype <- endsWith(scenario, "celltype")
    request <- validate_hybrid_request(loaded, list(readout = "abundance", markers = c("A", "B"),
      analysis_scope = if (within) "within_donor" else "across_donors",
      comparison_axis = if (celltype) "celltype" else "condition",
      design = if (within) "within_donor" else if (startsWith(scenario, "paired")) "paired" else "independent",
      donor = if (within) "D1" else NULL, donors = unique(meta$donor),
      condition = if (celltype) "NC" else NULL, celltype = if (celltype) "T" else "B",
      target = if (celltype) "T" else "TRT", reference = if (celltype) "B" else "NC"))
    result <- run_hybrid_analysis(path, request)
    view <- loaded
    view$contract <- hybrid_result_contract(loaded, list(id = scenario, label = scenario,
      request = request, result = result))
    expect_true(validate_precomputed_data(view), info = scenario)
    expect_identical(view$contract$sample_summaries, loaded$contract$sample_summaries, info = scenario)
    output <- tempfile(fileext = ".rds")
    on.exit(unlink(output), add = TRUE)
    write_precomputed_bundle(view, output)
    restored <- read_precomputed_rds(output)
    expect_identical(restored$contract, view$contract, info = scenario)
    if (within) {
      expect_true(all(view$contract$donor_contrasts$abundance$inference == "exploratory_component_level"))
    } else {
      expect_true(all(view$contract$differential_results$abundance$replicate_unit == "donor"))
      expect_true(all(is.finite(view$contract$differential_results$abundance$p_value)))
    }
  }
})
