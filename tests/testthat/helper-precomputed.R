source(file.path("..", "..", "R", "precomputed_contract.R"))
source(file.path("..", "..", "R", "report_plot_sources.R"))
source(file.path("..", "..", "R", "report_plots.R"))
if (file.exists(file.path("..", "..", "R", "precomputed_module.R"))) {
  source(file.path("..", "..", "R", "precomputed_module.R"))
}

precomputed_fixture <- function() {
  meta <- expand.grid(donor = c("D1", "D2"), condition = c("NC", "TRT"),
                      celltype = c("B", "T"), replicate = 1:2, stringsAsFactors = FALSE)
  meta$sample_id <- paste(meta$donor, meta$condition, sep = "_")
  meta$cluster <- ifelse(meta$celltype == "B", "1", "2")
  rownames(meta) <- paste0("cell", seq_len(nrow(meta)))
  counts <- matrix(rep(1:6, length.out = 3 * nrow(meta)), nrow = 3,
                   dimnames = list(c("A", "B", "C"), rownames(meta)))
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts, assay = "PNA", meta.data = meta)
  SeuratObject::LayerData(object, assay = "PNA", layer = "data") <- log1p(counts)
  for (name in c("umap", "pca")) {
    coords <- cbind(seq_len(nrow(meta)), seq_len(nrow(meta)) / 2)
    rownames(coords) <- rownames(meta)
    colnames(coords) <- paste0(toupper(name), "_", 1:2)
    object[[name]] <- SeuratObject::CreateDimReducObject(coords, key = paste0(toupper(name), "_"), assay = "PNA")
  }
  contract <- list(
    schema_version = 1L, dataset_id = "fixture-v1", title = "Precomputed fixture",
    defaults = list(reduction = "umap", pca_reduction = "pca", celltype_column = "celltype",
                    cluster_column = "cluster", donor_column = "donor", condition_column = "condition",
                    sample_column = "sample_id", assay = "PNA", abundance_layer = "data"),
    sample_summaries = list(), donor_contrasts = list(), differential_results = list(),
    qc = list(cell_counts = data.frame(stage = c("loaded", "retained"), n_cells = c(20L, 16L))),
    provenance = list(inference = "descriptive", missing_score_policy = "missing_not_zero", source = "synthetic fixture"),
    proximity_source = list(pxl_by_sample = c(D1_NC = "/unmounted/not-needed.pxl"))
  )
  for (readout in PROXIOMEVIS_READOUTS) {
    feature <- if (readout == "colocalization") c("A × B", "A × C") else c("A", "B")
    grid <- merge(unique(meta[c("sample_id", "donor", "condition", "celltype")]),
                  data.frame(feature_id = feature), by = NULL)
    grid$marker_1 <- if (readout == "colocalization") "A" else grid$feature_id
    grid$marker_2 <- if (readout == "colocalization") sub("A × ", "", grid$feature_id) else grid$feature_id
    grid$value <- 1 + (grid$condition == "TRT") * ifelse(grid$donor == "D1", 1, 2)
    grid$n_total_cells <- 2L
    grid$n_evaluable_cells <- 2L
    grid$n_scored_cells <- 2L
    grid$metric <- if (readout == "abundance") "median_clr" else "median_z"
    grid$is_control <- FALSE
    contract$sample_summaries[[readout]] <- grid
    contrast <- grid[grid$condition == "TRT", c("donor", "celltype", "feature_id", "marker_1", "marker_2", "metric", "is_control")]
    contrast$analysis_id <- paste(readout, "TRT_vs_NC", sep = ":")
    contrast$group_a <- "TRT"
    contrast$group_b <- "NC"
    contrast$comparison <- "TRT vs NC"
    contrast$value_a <- ifelse(contrast$donor == "D1", 2, 3)
    contrast$value_b <- 1
    contrast$effect_size <- contrast$value_a - contrast$value_b
    contrast$n_cells_a <- 2L
    contrast$n_cells_b <- 2L
    contract$donor_contrasts[[readout]] <- contrast
    result <- contrast[contrast$donor == "D1", setdiff(names(contrast), c("donor", "value_a", "value_b", "n_cells_a", "n_cells_b"))]
    result$effect_size <- 1.5
    result$min_effect <- 1
    result$max_effect <- 2
    result$n_pairs <- 2L
    result$direction_consistent <- TRUE
    result$p_value <- NA_real_
    result$p_adj <- NA_real_
    result$adjustment_family <- paste(readout, "TRT_vs_NC", result$celltype, sep = ":")
    result$method <- "mean of paired donor differences"
    result$replicate_unit <- "donor"
    result$inference <- "descriptive"
    contract$differential_results[[readout]] <- result
  }
  object@misc$proxiomevis <- contract
  object
}
