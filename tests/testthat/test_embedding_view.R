test_that("embedding splits align metadata and preserve shared scales", {
  data <- precomputed_app_data(precomputed_fixture())
  meta <- data$metadata[nrow(data$metadata):1L, ]
  for (split in c("donor", "sample_id")) {
    for (color in c("celltype", "abundance")) {
      plot <- plot_precomputed_embedding(meta, data$embeddings$umap, color,
        abundance = data$abundance, marker = "A", split_by = split)
      expect_identical(plot$data$split_group, data$metadata[[split]])
      built <- ggplot2::ggplot_build(plot)
      expect_equal(nrow(built$layout$layout), length(unique(meta[[split]])))
      expect_equal(unique(built$layout$layout$SCALE_X), 1L)
      expect_equal(unique(built$layout$layout$SCALE_Y), 1L)
      expect_s3_class(plotly::plotly_build(plotly::ggplotly(plot, tooltip = "text")), "plotly")
    }
  }
  expect_error(plot_precomputed_embedding(meta, data$embeddings$umap,
    split_by = "invalid"), "Choose donor or sample alias")
})

test_that("single embedding view follows reduction and split selection and exports SVG", {
  current <- shiny::reactiveVal(precomputed_app_data(precomputed_fixture()))
  shiny::testServer(precomputed_server, args = list(data = current), {
    session$setInputs(embedding_color = "celltype", embedding_split = "none")
    expect_equal(embedding_ggplot()$labels$x, "UMAP_1")
    expect_s3_class(embedding_ggplot()$facet, "FacetNull")
    session$setInputs(embedding_reduction = "pca", embedding_split = "sample_id")
    expect_equal(embedding_ggplot()$labels$x, "PCA_1")
    expect_equal(embedding_ggplot()$data$split_group, current()$metadata$sample_id)
    expect_gt(file.info(output$download_embedding_svg)$size, 0)
    session$setInputs(embedding_split = "donor", embedding_color = "abundance", embedding_marker = "A")
    expect_equal(embedding_ggplot()$data$color_value, as.numeric(current()$abundance["A", ]))
    expect_equal(embedding_ggplot()$data$split_group, current()$metadata$donor)
  })
  html <- as.character(precomputed_ui("report"))
  expect_match(html, "report-embedding_plot")
  expect_false(grepl("report-(umap|pca)_plot", html))
})
