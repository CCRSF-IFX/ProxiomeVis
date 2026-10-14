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
  write_cell_graph_glb(asset, "CD82", path, style)
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
