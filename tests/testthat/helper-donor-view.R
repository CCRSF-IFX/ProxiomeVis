comparison_view_fixture <- function(independent = FALSE) {
  data <- precomputed_app_data(precomputed_fixture())
  meta <- expand.grid(donor = paste0("D", 1:4), condition = c("NC", "TRT"),
    celltype = c("B", "T"), cell = 1:4, stringsAsFactors = FALSE)
  if (independent) meta <- meta[(meta$donor %in% c("D1", "D2") & meta$condition == "NC") |
    (meta$donor %in% c("D3", "D4") & meta$condition == "TRT"), ]
  meta$component <- rownames(meta) <- paste0("view-cell-", seq_len(nrow(meta)))
  meta$sample_id <- paste(meta$donor, meta$condition, sep = "_")
  meta$cluster <- ifelse(meta$celltype == "B", "1", "2")
  at <- match(meta$donor, paste0("D", 1:4))
  value <- at * 10 + (meta$condition == "TRT") * c(1, 2, 4, 7)[at] +
    (meta$celltype == "T") * c(2, 4, 5, 9)[at]
  data$metadata <- meta
  data$abundance <- Matrix::Matrix(rbind(A = value, B = value / 2, C = value / 4), sparse = TRUE)
  colnames(data$abundance) <- meta$component
  data$counts <- Matrix::Matrix(matrix(1, 3, nrow(meta), dimnames = dimnames(data$abundance)), sparse = TRUE)
  data$embeddings <- lapply(data$embeddings, function(x) matrix(seq_len(nrow(meta) * 2L), ncol = 2,
    dimnames = list(meta$component, colnames(x))))
  for (section in c("sample_summaries", "donor_contrasts", "differential_results"))
    data$contract[[section]] <- lapply(data$contract[[section]], function(rows) rows[FALSE, , drop = FALSE])
  data$contract$proximity_source <- NULL
  data$bundle_md5 <- "fixture-bundle-md5"
  data$bundle_path <- "/trusted/analysis.rds"
  data
}

comparison_view_entry <- function(data, axis = "condition", independent = FALSE, readout = "abundance", id = "new-1") {
  request <- validate_hybrid_request(data, list(analysis_scope = "across_donors", comparison_axis = axis,
    design = if (independent) "independent" else "paired", donors = paste0("D", 1:4),
    readout = readout, celltype = "B", condition = "NC", markers = c("A", "B", "C"),
    target = if (axis == "celltype") "T" else "TRT", reference = if (axis == "celltype") "B" else "NC",
    min_cells_per_group = 3L))
  reader <- function(data, asset_root, request, marker_1) {
    cells <- intersect(c(request$group_a_cells, request$group_b_cells),
      data$metadata$component[data$metadata$donor == request$donor])
    features <- hybrid_feature_table(request)
    features <- features[features$marker_1 %in% marker_1, , drop = FALSE]
    rows <- merge(data.frame(component = cells), features[c("marker_1", "marker_2")], by = NULL)
    rows$log2_ratio <- as.numeric(data$abundance["A", rows$component]) / 10
    rows
  }
  result <- run_donor_analysis(data, request, score_reader = reader)
  list(id = id, label = "Fixture comparison", request = request, result = result)
}
