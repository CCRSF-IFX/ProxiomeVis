#!/usr/bin/env Rscript
# Input RDS files must be trusted. This script does not run differential tests.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
for (file in c("precomputed_contract.R", "preparation_cli.R", "proximity_store.R", "cell_graph_assets.R", "bundle_builder.R")) {
  source(file.path(root, "R", file))
}
usage <- paste(
  "Usage: Rscript scripts/prepare_user_bundle.R INPUT_RDS NEW_OUTPUT_DIR [options]",
  "  --title TEXT                 Dataset title (default: My proximity dataset)",
  "  --dataset-id TEXT            Optional stable dataset identifier",
  "  --assay NAME                 Protein assay (default: PNA)",
  "  --abundance-layer NAME       Existing normalized layer (default: data)",
  "  --reduction NAME             Default embedding (default: umap)",
  "  --pca-reduction NAME         PCA embedding (default: pca)",
  "  --celltype-column NAME       Cell-type annotation (default: celltype)",
  "  --cluster-column NAME        Cluster annotation (default: seurat_clusters)",
  "  --donor-column NAME          Donor annotation (default: donor)",
  "  --condition-column NAME      Condition annotation (default: condition)",
  "  --sample-column NAME         Sample annotation (default: sample_id)",
  "  --saved-results RDS          Optional canonical result sections",
  "  --qc TSV                     Optional QC counts by sample and stage",
  "  --cell-map TSV               Exact component, original_id, sample_id mapping",
  "  --pxl-manifest TSV           sample_id and pxl_file; prepares proximity assets",
  "  --isotype-column NAME        Background fraction (default: isotype_fraction)",
  "  --graphs true|false          Also prepare 3D assets (default: false; needs FSMap)",
  "  --control-markers A,B        Optional control-marker list",
  sep = "\n")
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0L) }
options <- list(title = "My proximity dataset", "dataset-id" = NULL, assay = "PNA", "abundance-layer" = "data",
  reduction = "umap", "pca-reduction" = "pca", "celltype-column" = "celltype", "cluster-column" = "seurat_clusters",
  "donor-column" = "donor", "condition-column" = "condition", "sample-column" = "sample_id", "saved-results" = NULL,
  qc = NULL, "cell-map" = NULL, "pxl-manifest" = NULL, "isotype-column" = "isotype_fraction", graphs = "false",
  "control-markers" = NULL)
parsed <- parse_preparation_cli(args, options, 2L, usage)
opt <- parsed$options
if (!opt$graphs %in% c("true", "false")) stop("--graphs must be true or false.", call. = FALSE)
read_table <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
pxl <- NULL
if (!is.null(opt[["pxl-manifest"]])) {
  manifest_path <- normalizePath(opt[["pxl-manifest"]], mustWork = TRUE)
  manifest <- read_table(manifest_path)
  if (!all(c("sample_id", "pxl_file") %in% names(manifest))) {
    stop("The PXL manifest requires sample_id and pxl_file columns.", call. = FALSE)
  }
  paths <- as.character(manifest$pxl_file)
  relative <- !grepl("^/|^[A-Za-z]:[/\\\\]", paths)
  paths[relative] <- file.path(dirname(manifest_path), paths[relative])
  pxl <- stats::setNames(paths, as.character(manifest$sample_id))
}
controls <- opt[["control-markers"]]
if (!is.null(controls) && grepl(",\\s*$", controls)) stop("A control marker is missing after the final comma.", call. = FALSE)
path <- prepare_user_bundle(parsed$positional[[1L]], parsed$positional[[2L]],
  title = opt$title, dataset_id = opt[["dataset-id"]], assay = opt$assay, abundance_layer = opt[["abundance-layer"]],
  reduction = opt$reduction, pca_reduction = opt[["pca-reduction"]], celltype_column = opt[["celltype-column"]],
  cluster_column = opt[["cluster-column"]], donor_column = opt[["donor-column"]], condition_column = opt[["condition-column"]],
  sample_column = opt[["sample-column"]],
  saved_results = if (is.null(opt[["saved-results"]])) NULL else readRDS(opt[["saved-results"]]),
  qc = if (is.null(opt$qc)) NULL else read_table(opt$qc),
  cell_sources = if (is.null(opt[["cell-map"]])) NULL else read_table(opt[["cell-map"]]),
  pxl_by_sample = pxl, prepare_graphs = identical(opt$graphs, "true"), isotype_column = opt[["isotype-column"]],
  control_markers = if (is.null(controls)) character() else trimws(strsplit(controls, ",", fixed = TRUE)[[1L]]))
message("Prepared application bundle: ", path)
