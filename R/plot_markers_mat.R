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
                             vertical_gridlines = FALSE,
                             separate_gtable = FALSE) {
    min_lfp <- min_lfp %||% -3
    max_lfp <- max_lfp %||% 3

    gene_colors <- gene_colors %||% rep("black", nrow(mat))

    if (col_names) {
        col_names <- colnames(mat)
    }

    p_mat <- mcview_heatmap_geom(
        clip_vals(mat, min_lfp, max_lfp),
        col_names = col_names,
        interleave = interleave
    ) +
        scale_fill_gradient2(name = "", low = low_color, high = high_color, mid = mid_color, midpoint = midpoint, limits = c(min_lfp, max_lfp))

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

    if (vertical_gridlines) {
        p_mat <- p_mat + geom_hline(yintercept = 1:nrow(mat) - 0.5, color = "gray", size = 0.1)
    }

    cell_type_colors <- cell_type_colors %>% select(cell_type, color)
    mc_types <- tibble(metacell = colnames(mat)) %>%
        left_join(metacell_types %>% select("metacell", "cell_type"), by = "metacell") %>%
        left_join(cell_type_colors, by = "cell_type")

    p_mat <- p_mat + theme(legend.position = "top")

    if (plot_legend) {
        legend_point_size <- max(1, min(10, 250 / nrow(cell_type_colors)))
        legend <- cowplot::get_legend(cell_type_colors %>%
            ggplot(aes(x = cell_type, color = cell_type, y = 1)) +
            geom_point() +
            scale_color_manual("", values = deframe(cell_type_colors)) +
            guides(color = guide_legend(override.aes = list(size = legend_point_size), ncol = 1)))

        p <- add_markers_colorbars(p_mat, mc_types, dataset, top_cell_type_bar, metadata)

        if (separate_gtable) {
            return(list(p = p_mat, gtable = p, legend = legend))
        }

        cowplot::ggdraw(cowplot::plot_grid(p, legend, nrow = 1, rel_widths = c(0.8, 0.15)))
    } else {
        p <- add_markers_colorbars(p_mat, mc_types, dataset, top_cell_type_bar, metadata)

        if (separate_gtable) {
            return(list(p = p_mat, gtable = p))
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


#' In-package replacement for tgutil::tgplot_heatmap (slanted-labels variant)
#'
#' Reproduces the ggplot built by tgutil::tgplot_heatmap with
#' col_names_orient = "slanted" for the parameters plot_markers_mat uses
#' (col_names = FALSE or character vector; row_names always taken from
#' rownames(mtrx); interleave splits odd/even rows between left/right axes).
#' Other tgutil entry points (clip_vals, tgplot_add_axis_annotation, fread)
#' still go through the package; this function only inlines the heatmap path
#' so the markers panel's plot construction lives entirely in MCView.
#'
#' @param mtrx Numeric matrix to render as a raster heatmap.
#' @param col_names FALSE / NULL to suppress column labels, or a character
#'   vector to use as column labels (slanted, on both top and bottom axes).
#' @param interleave If TRUE, odd-row labels go on the left and even-row
#'   labels go on the right (matches tgplot_heatmap's interleave behaviour).
#' @return A ggplot object - the caller adds scale_fill_gradient2 etc.
#' @noRd
mcview_heatmap_geom <- function(mtrx, col_names = NULL, interleave = FALSE) {
    if (is.null(col_names)) {
        col_names <- colnames(mtrx)
    } else if (length(col_names) == 1 && is.logical(col_names) && !col_names) {
        col_names <- NULL
    }
    row_names <- rownames(mtrx)

    df <- expand.grid(
        y = seq_len(nrow(mtrx)),
        x = seq_len(ncol(mtrx)),
        KEEP.OUT.ATTRS = FALSE
    )
    df$val <- as.numeric(mtrx)

    ggp <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$val)) +
        ggplot2::geom_raster() +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.ticks = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank()
        )

    if (!is.null(col_names)) {
        ggp <- ggp +
            ggplot2::scale_x_continuous(
                breaks = seq_len(ncol(mtrx)),
                labels = col_names,
                expand = c(0, 0),
                sec.axis = ggplot2::dup_axis()
            ) +
            ggplot2::theme(
                axis.text.x.top    = ggplot2::element_text(angle = -45, hjust = 1, vjust = 0),
                axis.text.x.bottom = ggplot2::element_text(angle = -45, hjust = 0, vjust = 1)
            )
    } else {
        ggp <- ggp +
            ggplot2::scale_x_continuous(expand = c(0, 0)) +
            ggplot2::theme(axis.text.x = ggplot2::element_blank())
    }

    if (!is.null(row_names)) {
        if (interleave) {
            ggp <- ggp + ggplot2::scale_y_continuous(
                breaks = seq(1, nrow(mtrx), 2),
                labels = row_names[seq(1, nrow(mtrx), 2)],
                expand = c(0, 0),
                sec.axis = ggplot2::dup_axis(
                    breaks = seq(2, nrow(mtrx), 2),
                    labels = row_names[seq(2, nrow(mtrx), 2)]
                )
            )
        } else {
            ggp <- ggp + ggplot2::scale_y_continuous(
                breaks = seq_len(nrow(mtrx)),
                labels = row_names,
                expand = c(0, 0),
                sec.axis = ggplot2::dup_axis()
            )
        }
    } else {
        ggp <- ggp +
            ggplot2::scale_y_continuous(expand = c(0, 0)) +
            ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }

    ggp
}
