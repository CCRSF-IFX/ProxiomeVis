ranking_fixture <- function(data = precomputed_app_data(precomputed_fixture()), root = NULL) {
  skip_if_not_installed("arrow")
  skip_if_not_installed("dplyr")
  if (is.null(root)) root <- withr::local_tempdir(.local_envir = parent.frame())
  dir.create(root, showWarnings = FALSE)
  if (is.null(data$cell_sources)) data$cell_sources <- data.frame(component = data$metadata$component,
    original_id = paste0("raw", seq_len(nrow(data$metadata))), sample_id = data$metadata$sample_id)
  samples <- unique(data$metadata$sample_id)
  partitions <- lapply(seq_along(samples), function(i) {
    cells <- data$metadata$component[data$metadata$sample_id == samples[i]]
    scores <- data.frame(component = rep(cells, 2L), marker_1 = rep(c("A", "B"), each = length(cells)),
      marker_2 = "A", log2_ratio = rep(seq_along(cells), 2L), min_marker_fraction = 0)
    file <- paste0("sample-", i, ".parquet")
    arrow::write_parquet(scores, file.path(root, file))
    data.frame(sample_id = samples[i], file = file, md5 = unname(tools::md5sum(file.path(root, file))),
      source_pxl_md5 = paste(rep("a", 32), collapse = ""))
  })
  manifest <- list(schema_version = 1L, dataset_id = data$contract$dataset_id,
    dataset_fingerprint = hybrid_dataset_fingerprint(data), partitions = do.call(rbind, partitions),
    donor_cutoffs = c(D1 = .001, D2 = .001))
  saveRDS(manifest, file.path(root, "manifest.rds"))
  list(data = data, root = root, manifest = manifest)
}
