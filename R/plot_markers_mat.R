get_gene_colors <- function(genes,
                            lateral_genes = NULL,
                            noisy_genes = NULL,
                            disjoined_genes = NULL) {
    gene_colors <- tibble(gene = genes) %>%
        mutate(color = case_when(
            gene %in% lateral_genes & gene %in% noisy_genes ~ "purple",
            gene %in% lateral_genes ~ "blue",
            gene %in% noisy_genes ~ "red",
            gene %in% disjoined_genes ~ "darkgray",
            TRUE ~ "black"
        )) %>%
        deframe()
    return(gene_colors)
}


# The plot of tgutil::tgplot_heatmap(mat, col_names, col_names_orient = "slanted", interleave) + fill_scale, with
# the matrix drawn as a single raster image instead of a geom_raster layer with a row per cell, whose build
# is most of the render time of large heatmaps. The axes and ranges are the same, so shiny's hover, click
# and brush coordinates are unchanged.
tgplot_heatmap_raster <- function(mat, fill_scale, col_names = FALSE, interleave = FALSE) {
    nr <- nrow(mat)
    nc <- ncol(mat)
    row_names <- rownames(mat)
    colors <- matrix(fill_scale$map(as.vector(mat)), nrow = nr)

    # geom_blank trains the x, y and fill scales on the raster extent, as geom_raster did. show.legend = TRUE
    # keeps the fill colorbar, which ggplot otherwise drops since geom_blank has no fill aesthetic.
    extent <- data.frame(x = c(0.5, nc + 0.5), y = c(0.5, nr + 0.5), val = fill_scale$get_limits())
    p <- ggplot(extent, aes(x = x, y = y, fill = val)) +
        geom_blank(show.legend = TRUE) +
        annotation_raster(colors[nr:1, , drop = FALSE], xmin = 0.5, xmax = nc + 0.5, ymin = 0.5, ymax = nr + 0.5) +
        fill_scale +
        theme_bw() +
        theme(axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

    if (isFALSE(col_names)) {
        p <- p + scale_x_continuous(expand = c(0, 0)) + theme(axis.text.x = element_blank())
    } else {
        p <- p +
            scale_x_continuous(breaks = 1:nc, labels = col_names, expand = c(0, 0), sec.axis = dup_axis()) +
            theme(
                axis.text.x.top = element_text(color = "black", angle = -45, hjust = 1, vjust = 0),
                axis.text.x.bottom = element_text(color = "black", angle = -45, hjust = 0, vjust = 1)
            )
    }

    if (interleave) {
        p <- p + scale_y_continuous(
            breaks = seq(1, nr, 2), labels = row_names[seq(1, nr, 2)], expand = c(0, 0),
            sec.axis = dup_axis(breaks = seq(2, nr, 2), labels = row_names[seq(2, nr, 2)])
        )
    } else {
        p <- p + scale_y_continuous(breaks = 1:nr, labels = row_names, expand = c(0, 0), sec.axis = dup_axis())
    }

    return(p)
}

plot_markers_mat <- function(mat,
                             metacell_types,
                             cell_type_colors,
                             dataset,
                             low_color = "blue",
                             high_color = "red",
                             mid_color = "white",
                             midpoint = 0,
                             min_lfp = NULL,
                             max_lfp = NULL,
                             plot_legend = TRUE,
                             top_cell_type_bar = TRUE,
                             metadata = NULL,
                             gene_colors = NULL,
                             col_names = FALSE,
                             interleave = TRUE,
                             vertial_gridlines = FALSE,
                             separate_gtable = FALSE) {
    min_lfp <- min_lfp %||% -3
    max_lfp <- max_lfp %||% 3

    gene_colors <- gene_colors %||% rep("black", nrow(mat))

    if (col_names) {
        col_names <- colnames(mat)
    }

    p_mat <- tgplot_heatmap_raster(
        clip_vals(mat, min_lfp, max_lfp),
        scale_fill_gradient2(name = "", low = low_color, high = high_color, mid = mid_color, midpoint = midpoint, limits = c(min_lfp, max_lfp)),
        col_names = col_names,
        interleave = interleave
    )

    if (interleave) {
        p_mat <- p_mat +
            theme(
                axis.text.y.right = ggtext::element_markdown(color = gene_colors[seq(2, length(gene_colors), 2)]),
                axis.text.y.left =  ggtext::element_markdown(color = gene_colors[seq(1, length(gene_colors), 2)])
            )
    } else {
        p_mat <- p_mat +
            theme(axis.text.y = ggtext::element_markdown(color = gene_colors))
    }

    if (vertial_gridlines) {
        p_mat <- p_mat + geom_hline(yintercept = 1:nrow(mat) - 0.5, color = "gray", size = 0.1)
    }

    cell_type_colors <- cell_type_colors %>% select(cell_type, color)
    mc_types <- tibble(metacell = colnames(mat)) %>%
        left_join(metacell_types %>% select("metacell", "cell_type"), by = "metacell") %>%
        left_join(cell_type_colors, by = "cell_type")

    p_mat <- p_mat + theme(legend.position = "top")

    # build the heatmap once: the gtable is made from it and print.gt_custom reuses it for the plot coordinates
    p_build <- ggplot_build(p_mat)
    p_gtable <- ggplot_gtable(p_build)

    if (plot_legend) {
        legend_point_size <- max(1, min(10, 250 / nrow(cell_type_colors)))
        legend <- cowplot::get_legend(cell_type_colors %>%
            ggplot(aes(x = cell_type, color = cell_type, y = 1)) +
            geom_point() +
            scale_color_manual("", values = deframe(cell_type_colors)) +
            guides(color = guide_legend(override.aes = list(size = legend_point_size), ncol = 1)))

        p <- add_markers_colorbars(p_gtable, mc_types, dataset, top_cell_type_bar, metadata)

        if (separate_gtable) {
            return(list(p = p_mat, build = p_build, gtable = p, legend = legend))
        }

        cowplot::ggdraw(cowplot::plot_grid(p, legend, nrow = 1, rel_widths = c(0.8, 0.15)))
    } else {
        p <- add_markers_colorbars(p_gtable, mc_types, dataset, top_cell_type_bar, metadata)

        if (separate_gtable) {
            return(list(p = p_mat, build = p_build, gtable = p))
        }

        cowplot::ggdraw(p)
    }
}

add_markers_colorbars <- function(p, mc_types, dataset, top_cell_type_bar = TRUE, metadata = NULL) {
    if (!is.null(metadata)) {
        metadata <- mc_types %>%
            select(metacell) %>%
            left_join(metadata %>% mutate(metacell = as.character(metacell)), by = "metacell")

        for (md in rev(colnames(metadata)[-1])) {
            md_colors <- get_metadata_colors(dataset, md, metadata = metadata)
            if (is_numeric_field(metadata, md)) {
                palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
                p <- p %>%
                    tgplot_add_axis_annotation(palette(metadata[[md]]), position = "bottom", label = md)
            } else {
                p <- p %>%
                    tgplot_add_axis_annotation(md_colors[as.character(metadata[[md]])], position = "bottom", label = md)
            }
        }
    }
    p <- p %>% tgplot_add_axis_annotation(mc_types$color, label = "Cell type", plot_left = FALSE)
    if (top_cell_type_bar) {
        p <- p %>% tgplot_add_axis_annotation(mc_types$color, position = "top", label = "Cell type", plot_right = FALSE)
    }
    return(p)
}
