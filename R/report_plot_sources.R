# Frozen plotting functions from the Pixelgen tertiary-analysis workflow.
#
# These functions are copied into isolated environments so the deployed app
# does not source or execute the analysis scripts.  Update a copy only after
# reviewing the corresponding source and recording its new checksum here.
REPORT_PLOT_SOURCE_PROVENANCE <- data.frame(
  component = c(
    "abundance_volcano", "self_proximity_volcano",
    "colocalization_volcano", "spatial_layout", "comparable_spatial"
  ),
  source_path = c(
    "tertiary_analysis/scripts/differential_protein_abundance.R",
    "tertiary_analysis/scripts/differential_clustering.R",
    "tertiary_analysis/scripts/differential_colocalization.R",
    "tertiary_analysis/scripts/lib/spatial_plot_utils.R",
    "tertiary_analysis/scripts/plot_comparable_colocalization.R"
  ),
  md5 = c(
    "c257fa8cb870cdb6ebdbdced5e7b14f2",
    "aeb462d48e9cc3f36b16c4807bf522bb",
    "f902b4c4b846eb10c65bbb8157193444",
    "6600dc96a94d37a0d7f54bc9d719abf5",
    "76726d896b877253bb0907ade80863ae"
  ),
  stringsAsFactors = FALSE
)

REPORT_GGPLOT_BINDING_NAMES <- c(
  "aes", "annotate", "coord_cartesian", "coord_fixed", "element_blank",
  "element_rect", "element_text", "facet_wrap", "geom_blank", "geom_hline",
  "geom_point", "geom_text", "geom_tile", "geom_vline", "ggplot", "labs",
  "margin", "scale_color_gradient2", "scale_color_manual",
  "scale_fill_gradient2", "scale_size_area", "scale_size_continuous",
  "scale_x_discrete", "scale_y_discrete", "theme", "theme_bw", "theme_void"
)

new_report_plot_source_environment <- function(bindings = list()) {
  environment <- new.env(parent = baseenv())
  ggplot_bindings <- setNames(
    lapply(REPORT_GGPLOT_BINDING_NAMES, utils::getFromNamespace, ns = "ggplot2"),
    REPORT_GGPLOT_BINDING_NAMES
  )
  list2env(c(ggplot_bindings, bindings), envir = environment)
  environment
}

.report_spatial_source <- new_report_plot_source_environment()
evalq({
  spatial_reference_lines <- function(config, xintercept = NULL, yintercept = NULL,
                                      color = "grey35") {
    # Do not add invisible layers: hidden guides must not train the plot scales.
    if (!isTRUE(config$show_reference_lines)) return(list())
    list(
      if (!is.null(xintercept)) ggplot2::geom_vline(
        xintercept = xintercept, linetype = "dashed", color = color,
        alpha = config$reference_line_alpha
      ),
      if (!is.null(yintercept)) ggplot2::geom_hline(
        yintercept = yintercept, linetype = "dashed", color = color,
        alpha = config$reference_line_alpha
      )
    )
  }

  # Keep every cell type visible, including when no test survives the cell cutoff.
  spatial_volcano_layout <- function(celltypes, tested_celltypes, minimum_cells, unit) {
    scaffold <- data.frame(celltype = factor(celltypes, levels = celltypes))
    untested <- scaffold[!celltypes %in% tested_celltypes, , drop = FALSE]
    limits <- if (length(minimum_cells) == 1L && is.null(names(minimum_cells))) {
      rep(minimum_cells, nrow(untested))
    } else minimum_cells[as.character(untested$celltype)]
    stopifnot(length(limits) == nrow(untested), !anyNA(limits))
    untested$label <- paste0(
      "Not tested\nNo supported ", unit, " tests\nMinimum: ", limits,
      " cells per condition"
    )[seq_len(nrow(untested))]
    list(
      ggplot2::geom_blank(data = scaffold, ggplot2::aes(x = 0, y = 0), inherit.aes = FALSE),
      if (length(tested_celltypes)) ggplot2::geom_vline(
        data = scaffold[celltypes %in% tested_celltypes, , drop = FALSE],
        ggplot2::aes(xintercept = 0), linetype = "dashed", color = "grey55", linewidth = 0.35
      ),
      # A stable 0-1.5 range avoids tiny, negative axes in empty facets.
      if (nrow(untested)) ggplot2::geom_blank(data = untested, ggplot2::aes(x = 0, y = 1.5), inherit.aes = FALSE),
      if (nrow(untested)) ggplot2::geom_text(data = untested, ggplot2::aes(x = 0, y = 0.75, label = label),
                         inherit.aes = FALSE, color = "black", size = 3, lineheight = 1.2),
      ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.08))),
      ggplot2::theme_bw(),
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(color = "black", fill = NA))
    )
  }
}, envir = .report_spatial_source)

.report_abundance_source <- new_report_plot_source_environment(list(
  copy = data.table::copy,
  data.table = data.table::data.table,
  fifelse = data.table::fifelse,
  setNames = stats::setNames
))
evalq({
  make_volcano_plot <- function(
      daa,
      support,
      contrast_name,
      numerator,
      denominator,
      celltypes,
      facet_columns) {
    plot_data <- copy(daa)
    # Keep every cell type in the same order, including unavailable tests.
    facet_metadata <- copy(support)[match(celltypes, celltype)]
    facet_metadata[, facet_label := celltype]
    facet_metadata[
      status == "tested" & !meets_configured_minimum,
      facet_label := paste0(
        celltype, "\n[below recommended n: ", n_cells_target, " vs ",
        n_cells_reference, "]"
      )
    ]
    facet_levels <- facet_metadata$facet_label
    facet_lookup <- setNames(facet_metadata$facet_label, facet_metadata$celltype)
    plot_data[, facet_label := factor(facet_lookup[celltype], levels = facet_levels)]
    scaffold <- data.table(
      facet_label = factor(facet_levels, levels = facet_levels),
      scaffold_x = 0,
      scaffold_y = fifelse(facet_metadata$status == "tested", 0, 1)
    )
    unsupported <- facet_metadata[status != "tested"]
    if (nrow(unsupported)) {
      unsupported[, annotation := mapply(
        function(reason, n_target, n_reference) {
          wrapped_reason <- paste(strwrap(reason, width = 45), collapse = "\n")
          paste0(
            "Test unavailable\n", wrapped_reason,
            "\nTarget n = ", n_target, "; reference n = ", n_reference
          )
        },
        reason,
        n_cells_target,
        n_cells_reference,
        USE.NAMES = FALSE
      )]
      unsupported[, facet_label := factor(facet_label, levels = facet_levels)]
    }
    finite_p <- plot_data$p_adj[is.finite(plot_data$p_adj) & plot_data$p_adj > 0]
    p_floor <- if (length(finite_p)) min(finite_p) / 10 else .Machine$double.xmin
    p_floor <- max(p_floor, .Machine$double.xmin)
    plot_data[, plot_p_adj := fifelse(
      is.finite(p_adj) & p_adj >= 0,
      pmax(p_adj, p_floor),
      NA_real_
    )]

    finite_difference <- plot_data$difference[is.finite(plot_data$difference)]
    # The floor of 1 keeps the reference lines visible, even for empty results.
    x_extent <- max(1, abs(finite_difference))
    x_limits <- c(-1.1 * x_extent, 1.1 * x_extent)

    plot <- ggplot() +
      geom_blank(
        data = scaffold,
        aes(x = scaffold_x, y = scaffold_y),
        inherit.aes = FALSE
      ) +
      geom_point(
        data = plot_data,
        aes(x = difference, y = -log10(plot_p_adj), color = difference),
        na.rm = TRUE
      ) +
      scale_color_gradient2(
        low = "#0066FFFF",
        mid = "#74D055FF",
        high = "#FF0000FF",
        midpoint = 0
      ) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_vline(xintercept = c(-1, 1), color = "grey") +
      coord_cartesian(xlim = x_limits, clip = "off") +
      theme_bw() +
      geom_text(
        data = plot_data,
        aes(
          x = difference,
          y = -log10(plot_p_adj),
          label = ifelse(abs(difference) > 1, as.character(marker), "")
        ),
        hjust = -0.2,
        vjust = 0,
        size = 2.5,
        check_overlap = TRUE,
        na.rm = TRUE
      ) +
      facet_wrap(~facet_label, scales = "free_y", ncol = facet_columns, drop = FALSE) +
      labs(
        title = paste0(contrast_name, ": ", numerator, " vs ", denominator),
        subtitle = paste0(
          "RunDAA performed independently within each cell type; unavailable tests are annotated; ",
          "testable groups below the recommended cell count are marked in facet labels"
        ),
        x = "Difference in mean normalized abundance",
        y = expression(-log[10](adjusted~italic(P))),
        color = "Difference"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(5.5, 24, 5.5, 5.5)
      )
    if (nrow(unsupported)) {
      plot <- plot + geom_text(
        data = unsupported,
        aes(x = 0, y = 0.5, label = annotation),
        inherit.aes = FALSE,
        size = 3.2,
        color = "grey30",
        lineheight = 1.05
      )
    }
    plot
  }
}, envir = .report_abundance_source)

.report_self_proximity_source <- new_report_plot_source_environment(list(
  copy = data.table::copy,
  data.table = data.table::data.table,
  fifelse = data.table::fifelse,
  spatial_reference_lines = .report_spatial_source$spatial_reference_lines,
  spatial_volcano_layout = .report_spatial_source$spatial_volcano_layout
))
evalq({
  make_volcano_plot <- function(results, profile, contrast_name, numerator,
                                denominator, celltypes, profile_value_column,
                                clustering_cfg, volcano_x_limit,
                                minimum_cells = clustering_cfg$min_cells_per_group) {
    reference_lines <- spatial_reference_lines(
      clustering_cfg,
      xintercept = c(-clustering_cfg$volcano_difference_threshold,
                     clustering_cfg$volcano_difference_threshold),
      yintercept = -log10(clustering_cfg$fdr_threshold)
    )
    plot_data <- copy(results)[contrast == contrast_name]
    target_profile <- if (nrow(profile) && profile_value_column %in% names(profile)) {
      copy(profile)[
        condition == numerator & marker_1 == marker_2,
        .(
          celltype,
          marker_1,
          marker_2,
          target_profile_value = get(profile_value_column)
        )
      ]
    } else {
      data.table(
        celltype = character(), marker_1 = character(), marker_2 = character(),
        target_profile_value = numeric()
      )
    }
    plot_data <- merge(
      plot_data,
      target_profile,
      by = c("celltype", "marker_1", "marker_2"),
      all.x = TRUE,
      sort = FALSE
    )
    finite_p <- plot_data$p_adj[is.finite(plot_data$p_adj) & plot_data$p_adj > 0]
    p_floor <- if (length(finite_p)) min(finite_p) / 10 else .Machine$double.xmin
    p_floor <- max(p_floor, .Machine$double.xmin)
    plot_data[, plot_p_adj := fifelse(is.finite(p_adj), pmax(p_adj, p_floor), NA_real_)]
    plot_data[, celltype := factor(celltype, levels = celltypes)]
    title <- paste0(contrast_name, ": differential protein clustering")
    subtitle <- paste0(
      "DifferentialProximityAnalysis per cell type; ", numerator, " - ", denominator,
      "; color is target-condition ", gsub("_", " ", profile_value_column),
      "; cutoffs: |difference| >= ",
      signif(clustering_cfg$volcano_difference_threshold, 4),
      " and adjusted P <= ", signif(clustering_cfg$fdr_threshold, 4),
      "; shared x-axis: [", signif(-volcano_x_limit, 4),
      ", ", signif(volcano_x_limit, 4), "]"
    )

    plot <- ggplot(
      plot_data,
      aes(x = diff_median, y = -log10(plot_p_adj), color = target_profile_value)
    ) +
      spatial_volcano_layout(celltypes, unique(plot_data$celltype),
                             minimum_cells, "marker") +
      geom_point(na.rm = TRUE) +
      ggrepel::geom_text_repel(
        aes(label = marker_1),
        color = "black",
        max.overlaps = clustering_cfg$volcano_label_max_overlaps,
        size = 2.6,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      scale_color_gradient2(
        low = "#0066FFFF",
        mid = "#F7F7F7",
        high = "#FF0000FF",
        midpoint = 0,
        na.value = "grey70"
      ) +
      reference_lines +
      coord_cartesian(
        xlim = c(-volcano_x_limit, volcano_x_limit),
        clip = "off"
      ) +
      facet_wrap(~celltype, scales = "free_y", ncol = min(3L, length(celltypes))) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      ) +
      labs(
        title = title,
        subtitle = subtitle,
        x = "Difference in medians",
        y = expression(-log[10] ~ (adj. ~ p - value)),
        color = paste("Target", gsub("_", " ", profile_value_column))
      )
    plot
  }
}, envir = .report_self_proximity_source)

.report_colocalization_source <- new_report_plot_source_environment(list(
  copy = data.table::copy,
  head = utils::head,
  tail = utils::tail,
  spatial_reference_lines = .report_spatial_source$spatial_reference_lines,
  spatial_volcano_layout = .report_spatial_source$spatial_volcano_layout
))
evalq({
  safe_p_floor <- function(p_values) {
    finite_positive <- p_values[is.finite(p_values) & p_values > 0]
    floor_value <- if (length(finite_positive)) min(finite_positive) / 10 else .Machine$double.xmin
    max(floor_value, .Machine$double.xmin)
  }

  make_volcano_plot <- function(results, contrast_name, numerator, denominator,
                                label_extremes, difference_threshold,
                                fdr_threshold, x_limits, line_cfg, celltypes,
                                minimum_cells_per_group) {
    reference_lines <- spatial_reference_lines(
      line_cfg, xintercept = c(-difference_threshold, difference_threshold),
      yintercept = -log10(fdr_threshold), color = "grey40"
    )
    plot_data <- copy(results)
    plot_data[, celltype := factor(celltype, levels = celltypes)]
    p_floor <- safe_p_floor(plot_data$p_adj)
    plot_data[, plot_p_adj := ifelse(
      is.finite(p_adj), pmax(p_adj, p_floor), NA_real_
    )]
    plot_data[, pair_label := NA_character_]
    if (label_extremes > 0L) {
      plot_data[, pair_label := {
        labels <- rep(NA_character_, .N)
        finite_rows <- which(is.finite(diff_median))
        if (length(finite_rows)) {
          ordered <- finite_rows[order(diff_median[finite_rows])]
          selected <- unique(c(head(ordered, label_extremes), tail(ordered, label_extremes)))
          labels[selected] <- pair[selected]
        }
        labels
      }, by = celltype]
    }

    ggplot(
      plot_data,
      aes(x = diff_median, y = -log10(plot_p_adj), label = pair_label)
    ) +
      spatial_volcano_layout(celltypes, unique(plot_data$celltype),
                             minimum_cells_per_group, "marker-pair") +
      reference_lines +
      geom_point(alpha = 0.7, color = "#3B6EA8") +
      ggrepel::geom_text_repel(
        color = "black",
        na.rm = TRUE,
        max.overlaps = max(10L, 2L * label_extremes),
        size = 2.5,
        min.segment.length = 0
      ) +
      facet_wrap(~celltype, scales = "free_y", ncol = 3) +
      coord_cartesian(xlim = x_limits, clip = "off") +
      labs(
        title = contrast_name,
        subtitle = paste0(
          "DifferentialProximityAnalysis for each cell type: ",
          numerator, " - ", denominator,
          "; cutoffs: |difference| = ", signif(difference_threshold, 4),
          ", adjusted P = ", signif(fdr_threshold, 4)
        ),
        x = "Difference in medians",
        y = expression(-log[10] ~ (adj. ~ p - value))
      ) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold")
      )
  }
}, envir = .report_colocalization_source)

.report_comparable_source <- new_report_plot_source_environment(list(copy = data.table::copy))
evalq({
  axis_markers <- function(markers, kind) {
    stopifnot(kind %in% c("profiles", "differences"))
    if (kind == "profiles") markers else markers[tested_clustering | tested_colocalization]
  }

  comparison_plot <- function(bundle, celltype_value, kind, donor_value, panel_value) {
    cfg <- bundle$cfg
    panels <- if (kind == "profiles") bundle$conditions else bundle$contrasts
    stopifnot(kind %in% c("profiles", "differences"), donor_value %in% bundle$donors,
              panel_value %in% panels)
    order <- axis_markers(bundle$markers, kind)[celltype == celltype_value, marker]
    # Subset only the displayed observations; retain the cross-donor marker order/scales.
    data <- copy(bundle[[kind]][celltype == celltype_value & donor == donor_value & panel == panel_value])
    title <- paste(celltype_value, donor_value, panel_value, sep = " | ")
    if (!length(order)) {
      label <- if (kind == "profiles") "No supported proximity profiles after filtering" else
        "No tested clustering or colocalization pairs"
      return(ggplot() + annotate("text", x = 0, y = 0, label = label) +
               theme_void() + labs(title = title))
    }
    data[, `:=`(marker_1 = factor(marker_1, levels = order), marker_2 = factor(marker_2, levels = order))]
    available <- data[status == "supported"]
    descriptive_only <- kind == "profiles" && !any(bundle$markers[
      celltype == celltype_value, tested_clustering | tested_colocalization])
    plot <- ggplot(data, aes(marker_1, marker_2)) +
      geom_tile(data = data[status == "unavailable"], fill = "grey82") +
      scale_x_discrete(limits = order, drop = FALSE, position = "top") +
      scale_y_discrete(limits = rev(order), drop = FALSE) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                           limits = c(-1, 1) * bundle$limits[[kind]], oob = scales::squish) +
      coord_fixed() + theme_bw(base_size = 10) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
            axis.text.y = element_text(size = 7), panel.border = element_rect(color = "black", fill = NA)) +
      labs(x = NULL, y = NULL, title = title,
           subtitle = paste(if (descriptive_only) "Descriptive profiles (no differential spatial tests)"
                            else if (kind == "profiles") "Profiles for proteins tested by self or co analysis"
                            else "Differential clustering and colocalization",
                            "- shared marker order and scales"),
           caption = paste("Grey: unavailable / not tested.",
                           if (kind == "differences") "Diagonal: clustering. Off-diagonal: colocalization."
                           else "Profile diagonal: supported self-pairs included.",
                           "White dots: near zero. Colors saturate outside the displayed range.",
                           if (kind == "differences") "Effects: target minus reference.\nSelf/co adjusted P values retained separately."
                           else if (descriptive_only) "Descriptive only: these patterns do not establish a treatment effect."
                           else "", sep = "\n"))
    if (kind == "profiles") {
      plot <- plot + geom_point(data = available, aes(fill = value, size = pct_detected), shape = 21, color = "grey30", stroke = 0.2) +
        scale_size_area(limits = c(0, 1), max_size = 3) + labs(fill = bundle$metric, size = "Fraction detected")
    } else {
      # Display-only encodings: never export a second copy of P-derived quantities.
      available[, `:=`(
        point_size = pmin(-log10(pmax(p_adj, .Machine$double.xmin)), cfg$significance_size_cap),
        significance = ifelse(p_adj <= cfg$significance_threshold, "Significant", "Not significant")
      )]
      plot <- plot + geom_point(data = available, aes(fill = value, size = point_size, color = significance), shape = 21, stroke = 0.35) +
        scale_size_continuous(limits = c(0, cfg$significance_size_cap), range = c(0.6, 3)) +
        scale_color_manual(values = c(Significant = "black", `Not significant` = "grey65"),
                           limits = c("Significant", "Not significant"), drop = FALSE) +
        labs(fill = "Difference in medians", color = paste0("Adjusted P <= ", cfg$significance_threshold),
             size = paste0("-log10(adj. P)\n(capped at ", cfg$significance_size_cap, ")"))
    }
    if (!nrow(available)) plot <- plot + annotate("label", x = order[ceiling(length(order)/2)],
      y = order[ceiling(length(order)/2)], size = 3, fill = "white",
      label = if (kind == "profiles") "No supported pairs" else "No tested pairs\nin shared marker set")
    plot
  }
}, envir = .report_comparable_source)
