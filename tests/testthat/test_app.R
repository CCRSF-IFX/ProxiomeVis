test_that("documentation is available without a dataset", {
  html <- as.character(precomputed_ui("report"))
  expect_match(html, "Documentation", fixed = TRUE)
  expect_match(html, "ProxiomeVis documentation", fixed = TRUE)
  expect_match(html, 'href="https://github.com/CCRSF-IFX/ProxiomeVis"', fixed = TRUE)
  expect_match(html, 'rel="noopener noreferrer"', fixed = TRUE)
  expect_match(html, "Cells are not independent biological replicates.", fixed = TRUE)
  expect_match(html, "Download rotating GIF", fixed = TRUE)
})

test_that("the hosted app rejects visitor-selected server paths", {
  previous_directory <- setwd(file.path("..", ".."))
  on.exit(setwd(previous_directory), add = TRUE)
  fixture_path <- tempfile(fileext = ".rds")
  saveRDS(precomputed_fixture(), fixture_path)
  on.exit(unlink(fixture_path), add = TRUE)
  withr::local_envvar(c(PROXIOME_PRECOMPUTED_RDS = fixture_path, PROXIOME_ALLOW_LOCAL_LOAD = "false"))
  environment <- new.env(parent = globalenv())
  sys.source("app.R", envir = environment)
  expect_false(grepl('id="rds_path"', as.character(environment$ui), fixed = TRUE))
  shiny::testServer(environment$server, {
    modal <- NULL
    session$sendModal <- function(type, message) modal <<- message$html
    session$setInputs(dataset = 1)
    expect_match(modal, "Dataset", fixed = TRUE)
    expect_false(grepl('id="rds_path"', modal, fixed = TRUE))
    session$setInputs(rds_path = "/untrusted/not-present.rds", load_rds = 1)
    expect_identical(data()$contract$dataset_id, "fixture-v1")
    expect_null(error())
  })
})

test_that("failed trusted-local reload keeps the active dataset", {
  previous_directory <- setwd(file.path("..", ".."))
  on.exit(setwd(previous_directory), add = TRUE)
  fixture_path <- tempfile(fileext = ".rds")
  saveRDS(precomputed_fixture(), fixture_path)
  on.exit(unlink(fixture_path), add = TRUE)
  withr::local_envvar(c(PROXIOME_PRECOMPUTED_RDS = fixture_path, PROXIOME_ALLOW_LOCAL_LOAD = "true"))
  environment <- new.env(parent = globalenv())
  sys.source("app.R", envir = environment)
  expect_false(grepl('id="rds_path"', as.character(environment$ui), fixed = TRUE))
  shiny::testServer(environment$server, {
    modal <- NULL
    session$sendModal <- function(type, message) modal <<- message$html
    session$setInputs(dataset = 1)
    expect_true(grepl('id="rds_path"', modal, fixed = TRUE))
    session$setInputs(rds_path = "/untrusted/not-present.rds", load_rds = 1)
    expect_identical(data()$contract$dataset_id, "fixture-v1")
    expect_match(error(), "existing precomputed")
  })
})
