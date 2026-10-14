graph_asset_fixture <- function() {
  nodes <- data.frame(cell_id = "cell1", x = as.double(1:4), y = c(0,1,0,1), z = c(1,1,2,2))
  nodes$markers <- list(character(), "CD81", "CD82", c("CD81", "CD82"))
  list(schema_version = 3L, cell_id = "cell1", nodes = nodes)
}

read_graph_glb <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  header <- readBin(con, "integer", n = 5, size = 4, endian = "little")
  expect_equal(header[c(1, 2, 3, 5)], c(0x46546C67, 2, file.info(path)$size, 0x4E4F534A))
  doc <- jsonlite::fromJSON(rawToChar(readBin(con, "raw", n = header[4])), simplifyVector = FALSE)
  bin <- readBin(con, "integer", n = 2, size = 4, endian = "little")
  expect_equal(bin, c(doc$buffers[[1]]$byteLength, 0x004E4942))
  list(doc = doc, bytes = readBin(con, "raw", n = bin[1]))
}

test_that("GLB preserves nodes, outward tetrahedra, selected marker and applied appearance", {
  asset <- graph_asset_fixture()
  style <- cell_graph_appearance("#123456", "#00FF00", 0.4, 0.9, 0.005, 0.007)
  path <- tempfile(fileext = ".glb")
  on.exit(unlink(path))
  write_cell_graph_glb(asset, "CD82", path, style, background_retention = 100)
  model <- read_graph_glb(path)
  expect_length(model$bytes, nrow(asset$nodes) * 96)
  expect_equal(model$doc$extras$highlight_marker, "CD82")
  expect_equal(model$doc$materials[[1]]$alphaMode, "BLEND")
  expect_equal(model$doc$materials[[2]]$pbrMetallicRoughness$baseColorFactor, as.list(c(0, 1, 0, 0.9)))
  for (g in 1:2) {
    pv <- model$doc$bufferViews[[g * 2 - 1]]
    iv <- model$doc$bufferViews[[g * 2]]
    positions <- matrix(readBin(model$bytes[seq.int(pv$byteOffset + 1, length.out = pv$byteLength)],
      "double", n = pv$byteLength / 4, size = 4, endian = "little"), ncol = 3, byrow = TRUE)
    indices <- readBin(model$bytes[seq.int(iv$byteOffset + 1, length.out = iv$byteLength)],
      "integer", n = iv$byteLength / 4, size = 4, endian = "little")
    expect_true(all(indices >= 0 & indices < nrow(positions)))
    rows <- if (g == 1) 1:2 else 3:4
    radius <- if (g == 1) 0.005 else 0.007
    for (j in 1:2) {
      pts <- positions[(j - 1) * 4 + 1:4, ]
      center <- as.numeric(asset$nodes[rows[j], c("x", "y", "z")])
      expect_equal(colMeans(pts), center, tolerance = 1e-6)
      expect_equal(sqrt(rowSums(sweep(pts, 2, center)^2)), rep(radius * 3, 4), tolerance = 1e-6)
      for (face in split(indices[(j - 1) * 12 + 1:12] %% 4 + 1, rep(1:4, each = 3))) {
        a <- pts[face[2], ] - pts[face[1], ]; b <- pts[face[3], ] - pts[face[1], ]
        normal <- c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
        expect_gt(sum(normal * (pts[face[1], ] - center)), 0)
      }
    }
  }
  asset$nodes$markers <- rep(list("CD82"), 4)
  write_cell_graph_glb(asset, "CD82", path)
  expect_length(read_graph_glb(path)$doc$meshes[[1]]$primitives, 1)
  expect_error(write_cell_graph_glb(asset, "absent", path), "highlight marker")
})

test_that("appearance is bounded and default preview retains the threejs style", {
  expect_error(cell_graph_appearance(background_color = "red"), "hex")
  expect_error(cell_graph_appearance(highlight_opacity = 2), "Opacity")
  expect_error(cell_graph_appearance(background_radius = NA_real_), "Opacity")
  asset <- graph_asset_fixture()
  plain <- plot_cell_graph_asset(asset, "CD82")
  styled <- plot_cell_graph_asset(asset, "CD82", cell_graph_appearance())
  expect_equal(plain$x$color, styled$x$color)
  expect_equal(plain$x$alpha, styled$x$alpha)
  expect_true(all(unlist(styled$x$size) == 0.1))
})

test_that("GLB retention is reproducible, preserves highlights and scale, and leaves RNG untouched", {
  asset <- graph_asset_fixture()
  asset$nodes <- asset$nodes[rep(1:4, 5), ]
  asset$nodes$x <- as.double(1:20)
  full <- tempfile(fileext = ".glb")
  subset <- tempfile(fileext = ".glb")
  on.exit(unlink(c(full, subset)))
  write_cell_graph_glb(asset, "CD82", full, background_retention = 100)
  withr::local_seed(123)
  state <- .Random.seed
  write_cell_graph_glb(asset, "CD82", subset)
  expect_identical(.Random.seed, state)
  original <- read_graph_glb(full)
  sampled <- read_graph_glb(subset)
  expect_equal(sampled$doc$extras$node_count, 12)
  expect_equal(sampled$doc$extras$original_node_count, 20)
  expect_equal(sampled$doc$extras$highlighted_node_count, 10)
  expect_equal(sampled$doc$extras$retained_background_nodes, 2)
  expect_true(sampled$doc$extras$subsampled)
  expect_length(sampled$bytes, 12 * 96)
  positions <- function(model, group) {
    view <- model$doc$bufferViews[[group * 2 - 1]]
    model$bytes[seq.int(view$byteOffset + 1, length.out = view$byteLength)]
  }
  expect_identical(positions(original, 2), positions(sampled, 2))
  full_background <- matrix(readBin(positions(original, 1), "double", n = 120, size = 4, endian = "little"), ncol = 12, byrow = TRUE)
  sampled_background <- matrix(readBin(positions(sampled, 1), "double", n = 24, size = 4, endian = "little"), ncol = 12, byrow = TRUE)
  expect_true(all(apply(sampled_background, 1, function(row) any(apply(full_background, 1, function(other) identical(row, other))))))
  checksum <- unname(tools::md5sum(subset))
  set.seed(999)
  write_cell_graph_glb(asset, "CD82", subset)
  expect_equal(unname(tools::md5sum(subset)), checksum)
  write_cell_graph_glb(asset, "CD82", subset, background_retention = 0)
  expect_equal(read_graph_glb(subset)$doc$extras$node_count, 10)
  for (bad in list(NA_real_, -1, 101, "25", numeric()))
    expect_error(write_cell_graph_glb(asset, "CD82", subset, background_retention = bad), "retention")
})
