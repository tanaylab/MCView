plot_markers_mat <- function(mc_fp,
                             metacell_types,
                             cell_type_colors,
                             dataset,
                             colors = c("darkblue", "blue", "lightblue", "white", "red", "darkred"),
                             mid_color = 4,
                             min_lfp = NULL,
                             max_lfp = NULL,
                             plot_legend = TRUE,
                             top_cell_type_bar = TRUE,
                             metadata = NULL) {
    min_lfp <- min_lfp %||% -3
    max_lfp <- max_lfp %||% 3

    values <- c(
        seq(min_lfp, 0, length.out = mid_color),
        seq(0, max_lfp, length.out = length(colors) - mid_color + 1)[-1]
    )

    gene_ord <- order(apply(mc_fp, 1, which.max))
    mat <- log2(mc_fp[gene_ord, ])

    p_mat <- tgutil::tgplot_heatmap(
        clip_vals(mat, min_lfp, max_lfp),
        col_names = FALSE,
        interleave = TRUE
    ) +
        scale_fill_gradientn(name = "Fold change", colors = colors, values = scales::rescale(values))

    cell_type_colors <- cell_type_colors %>% select(cell_type, color)
    mc_types <- tibble(metacell = colnames(mat)) %>%
        left_join(metacell_types %>% select("metacell", "cell_type"), by = "metacell") %>%
        left_join(cell_type_colors, by = "cell_type")

    if (plot_legend) {
        legend <- cowplot::get_legend(cell_type_colors %>%
            ggplot(aes(x = cell_type, color = cell_type, y = 1)) +
            geom_point() +
            scale_color_manual("", values = deframe(cell_type_colors)) +
            guides(color = guide_legend(ncol = 1)))
        p_mat <- p_mat + theme(legend.position = "top")

        p <- add_markers_colorbars(p_mat, mc_types, dataset, top_cell_type_bar, metadata)

        cowplot::ggdraw(cowplot::plot_grid(p, legend, nrow = 1, rel_widths = c(0.8, 0.15)))
    } else {
        p <- add_markers_colorbars(p_mat, mc_types, dataset, top_cell_type_bar, metadata)

        cowplot::ggdraw(p)
    }
}

add_markers_colorbars <- function(p, mc_types, dataset, top_cell_type_bar = TRUE, metadata = NULL) {
    if (!is.null(metadata)) {
        metadata <- metadata %>%
            mutate(metacell = as.character(metacell)) %>%
            right_join(mc_types %>% select(metacell), by = "metacell")

        for (md in rev(colnames(metadata)[-1])) {
            md_colors <- get_metadata_colors(dataset, md, metadata = metadata)
            palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
            p <- p %>% tgplot_add_axis_annotation(palette(metadata[[md]]), position = "bottom", label = md)
        }
    }
    p <- p %>% tgplot_add_axis_annotation(mc_types$color, label = "Cell type", plot_left = FALSE)
    if (top_cell_type_bar) {
        p <- p %>% tgplot_add_axis_annotation(mc_types$color, position = "top", label = "Cell type", plot_right = FALSE)
    }
    return(p)
}