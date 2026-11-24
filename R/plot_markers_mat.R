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
    colors <- matrix(map_fill_colors(as.vector(mat), fill_scale), nrow = nr)

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

# fill_scale$map(values) through a table of 1024 colors along the limits of the scale: mapping every cell of
# a large heatmap through the scale took ~0.3 s per render. The colors differ from the exact ones by at most
# 2/255 in a channel. Values outside of the limits and NAs get the scale's na.value, as in the scale.
map_fill_colors <- function(values, fill_scale, n = 1024) {
    limits <- fill_scale$get_limits()
    if (length(values) <= n || !all(is.finite(limits)) || diff(limits) <= 0) {
        return(fill_scale$map(values))
    }
    lut <- fill_scale$map(seq(limits[1], limits[2], length.out = n))
    idx <- round((values - limits[1]) / diff(limits) * (n - 1)) + 1L
    idx[!is.finite(idx) | idx < 1 | idx > n] <- NA
    colors <- lut[idx]
    colors[is.na(colors)] <- fill_scale$na.value
    colors
}

# The labels of the genes, each in its color. Plain text elements with a vector of colors: ggplot2 warns that
# vectors are not officially supported, but ggtext::element_markdown(), which was used for it, took ~1 s per
# render of 100 genes (it parses every label as markdown) and is still used when the plain elements fail.
gene_labels_theme <- function(gene_colors, interleave, element = element_text) {
    if (interleave) {
        theme(
            axis.text.y.right = element(color = gene_colors[seq(2, length(gene_colors), 2)]),
            axis.text.y.left = element(color = gene_colors[seq(1, length(gene_colors), 2)])
        )
    } else {
        theme(axis.text.y = element(color = gene_colors))
    }
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

    gene_colors <- unname(gene_colors)

    if (vertial_gridlines) {
        p_mat <- p_mat + geom_hline(yintercept = 1:nrow(mat) - 0.5, color = "gray", size = 0.1)
    }

    cell_type_colors <- cell_type_colors %>% select(cell_type, color)
    mc_types <- tibble(metacell = colnames(mat)) %>%
        left_join(metacell_types %>% select("metacell", "cell_type"), by = "metacell") %>%
        left_join(cell_type_colors, by = "cell_type")

    p_mat <- p_mat + theme(legend.position = "top")
    p_plain <- p_mat

    # build the heatmap once: the gtable is made from it and print.gt_custom reuses it for the plot coordinates
    built <- tryCatch(
        suppressWarnings({
            p_mat <- p_plain + gene_labels_theme(gene_colors, interleave)
            p_build <- ggplot_build(p_mat)
            list(p = p_mat, build = p_build, gtable = ggplot_gtable(p_build))
        }),
        error = function(e) NULL
    )
    if (is.null(built)) {
        p_mat <- p_plain + gene_labels_theme(gene_colors, interleave, element = ggtext::element_markdown)
        p_build <- ggplot_build(p_mat)
        built <- list(p = p_mat, build = p_build, gtable = ggplot_gtable(p_build))
    }
    p_mat <- built$p
    p_build <- built$build
    p_gtable <- built$gtable

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
                    add_axis_annotation_bar(palette(metadata[[md]]), position = "bottom", label = md)
            } else {
                p <- p %>%
                    add_axis_annotation_bar(md_colors[as.character(metadata[[md]])], position = "bottom", label = md)
            }
        }
    }
    p <- p %>% add_axis_annotation_bar(mc_types$color, label = "Cell type", plot_left = FALSE)
    if (top_cell_type_bar) {
        p <- p %>% add_axis_annotation_bar(mc_types$color, position = "top", label = "Cell type", plot_right = FALSE)
    }
    return(p)
}

# tgutil::tgplot_add_axis_annotation() for a bar of colors above or below the heatmap. tgutil builds a ggplot
# for the bar and two more for its label, ~0.5 s per bar, on every render. Here the bar is a raster grob with
# the same border, and the grobs of a label, which depend only on its text, are made once per process.
add_axis_annotation_bar <- function(heatmap, colors, position = "bottom", size = 0.02, label = NULL, plot_left = TRUE, plot_right = TRUE) {
    colors <- unname(as.character(colors))
    colors[is.na(colors)] <- "gray50" # as geom_raster drew NA fills
    bar <- grid::grobTree(
        grid::rasterGrob(matrix(colors, nrow = 1), width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = FALSE),
        grid::rectGrob(gp = grid::gpar(fill = NA, col = "black", lwd = 0.2 * ggplot2::.pt))
    )
    heatmap <- cowplot::insert_xaxis_grob(heatmap, bar, grid::unit(size, "null"), position = position)

    if (!is.null(label)) {
        grobs <- annotation_label_grobs(label)
        left <- heatmap$layout[heatmap$layout$name == "axis-l", ]
        right <- heatmap$layout[heatmap$layout$name == "axis-r", ]
        suffix <- if (position == "top") "t" else "b"
        row <- function(pp) if (position == "top") pp$t - 1 else pp$b + 1
        if (plot_left) {
            heatmap <- gtable::gtable_add_grob(heatmap, grobs$left, row(left), left$l, row(left), left$r, name = paste0("xaxis-label-grob-", suffix))
        }
        if (plot_right) {
            heatmap <- gtable::gtable_add_grob(heatmap, grobs$right, row(right), right$l, row(right), right$r, name = paste0("xaxis-label-grob-", suffix, "r"))
        }
    }
    heatmap
}

annotation_label_cache <- new.env(parent = emptyenv())

# The axis grobs that tgutil draws next to an annotation bar as its label, on the left and on the right
annotation_label_grobs <- function(label) {
    label <- as.character(label)
    if (is.null(annotation_label_cache[[label]])) {
        gg <- ggplot(data.frame(), aes(y = 1, x = 1)) +
            theme_bw() +
            theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(color = "black")) +
            scale_y_continuous(breaks = 1, labels = label, sec.axis = dup_axis())
        annotation_label_cache[[label]] <- list(
            left = cowplot::get_plot_component(gg, "axis-l"),
            right = cowplot::get_plot_component(gg, "axis-r")
        )
    }
    annotation_label_cache[[label]]
}
