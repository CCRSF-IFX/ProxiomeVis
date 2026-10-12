source(file.path("..", "..", "R", "hybrid_analysis.R"))
source(file.path("..", "..", "R", "proximity_store.R"))

hybrid_fixture <- function() {
  data <- precomputed_app_data(precomputed_fixture())
  metadata <- data$metadata[rep(1L, 84L), , drop = FALSE]
  metadata$component <- rownames(metadata) <- paste0("selected", seq_len(nrow(metadata)))
  metadata$condition <- rep(c("TRT", "NC"), each = 42L)
  metadata$sample_id <- paste(metadata$donor, metadata$condition, sep = "_")
  metadata$isotype_fraction <- .0001
  counts <- matrix(rep(seq_len(84L), each = 3L), 3L, dimnames = list(c("A", "B", "C"), metadata$component))
  data$metadata <- metadata
  data$counts <- Matrix::Matrix(counts, sparse = TRUE)
  data$abundance <- log1p(data$counts)
  data$cell_sources <- data.frame(component = metadata$component, original_id = paste0("raw", seq_len(84L)),
                                  sample_id = metadata$sample_id)
  data
}

hybrid_request <- function(readout = "abundance") {
  list(readout = readout, donor = "D1", celltype = "B", target = "TRT", reference = "NC", markers = c("A", "B"))
}
