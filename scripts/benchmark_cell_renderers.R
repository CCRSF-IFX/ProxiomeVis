# Optional benchmark dependency: install.packages("threejs"). Run from the repository root.
# R_LIBS_USER may point to an isolated library containing CRAN threejs.
source("R/cell_graph_assets.R")
baseline <- new.env(parent = globalenv())
eval(parse(text = system2("git", c("show", "bb341b2:R/cell_graph_assets.R"), stdout = TRUE)), baseline)
args <- commandArgs(trailingOnly = TRUE)
directory <- if (length(args)) args[1] else "output/hybrid_demo/graphs"
output <- if (length(args) > 1) args[2] else "output/renderer-benchmark"
stopifnot(requireNamespace("threejs", quietly = TRUE))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
saved <- readRDS(file.path(directory, "index.rds"))
index <- read_cell_graph_index(directory, saved$dataset_id)
cells <- index$cells[index$cells$n_nodes <= 200000L, ]
cells <- cells[order(cells$n_nodes), ]
selected <- cells[unique(c(1L, ceiling(nrow(cells) / 2), nrow(cells))), ]

candidate <- function(asset, marker, orbit = FALSE) {
  nodes <- asset$nodes
  coordinates <- as.matrix(nodes[c("x", "y", "z")])
  ranges <- apply(coordinates, 2, range)
  span <- max(ranges[2, ] - ranges[1, ])
  coordinates <- sweep(coordinates, 2, colMeans(ranges)) / if (span > 0) span else 1
  colors <- ifelse(vapply(nodes$markers, function(x) marker %in% x, logical(1)), "#D73027D9", "#DEDEDED9")
  labels <- vapply(nodes$markers, function(x) paste0("Markers: ",
    if (length(x)) paste(htmltools::htmlEscape(x), collapse = "<br>") else "None"), character(1))
  threejs::scatterplot3js(coordinates, color = colors, labels = labels,
    axis = FALSE, grid = FALSE, size = 0.1, center = TRUE,
    xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5), zlim = c(-0.5, 0.5),
    use.orbitcontrols = orbit, width = 1000, height = 700)
}

results <- list()
for (i in seq_len(nrow(selected))) {
  asset <- load_cell_graph_asset(index, directory, selected$cell_id[i])
  marker <- names(sort(table(unlist(asset$nodes$markers)), decreasing = TRUE))[1]
  for (renderer in c("plotly", "threejs", "threejs_orbit")) {
    build <- function() {
      if (renderer == "plotly") {
        plot <- baseline$plot_cell_graph_asset(asset, marker)
        plot$width <- 1000
        plot$height <- 700
        plotly::plotly_build(plot)
      } else candidate(asset, marker, renderer == "threejs_orbit")
    }
    invisible(build())
    elapsed <- replicate(3, system.time(build())[["elapsed"]])
    widget <- build()
    widget$elementId <- "bench"
    widget <- htmlwidgets::onRender(widget, "function(el,x){requestAnimationFrame(function(){requestAnimationFrame(function(){window.benchReady=performance.now();});});}")
    filename <- paste0(renderer, "-", i, ".html")
    htmlwidgets::saveWidget(widget, file.path(output, filename), selfcontained = FALSE,
      libdir = "lib", title = paste(renderer, nrow(asset$nodes), "nodes"))
    results[[length(results) + 1L]] <- data.frame(renderer, file = filename,
      cell_id = asset$cell_id, nodes = nrow(asset$nodes), marker,
      r_build_ms = median(elapsed) * 1000,
      r_widget_bytes = as.numeric(object.size(widget)), html_bytes = file.info(file.path(output, filename))$size)
    message(renderer, ": ", nrow(asset$nodes), " nodes, ", median(elapsed), " s")
  }
}
results <- do.call(rbind, results)
write.csv(results, file.path(output, "r-results.csv"), row.names = FALSE)
jsonlite::write_json(list(versions = list(R = as.character(getRversion()),
  plotly = as.character(packageVersion("plotly")), threejs = as.character(packageVersion("threejs"))),
  cases = results), file.path(output, "cases.json"), auto_unbox = TRUE, pretty = TRUE)
print(results)
