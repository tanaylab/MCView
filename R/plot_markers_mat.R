plot_markers_mat <- function(mc_fp,
                             metacell_types,
                             cell_type_colors,
                             colors = c("darkblue", "blue", "lightblue", "white", "red", "darkred"),
                             mid_color = 4,
                             min_lfp = NULL,
                             max_lfp = NULL,
                             plot_legend = TRUE) {
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

    mc_types <- tibble(metacell = colnames(mat)) %>%
        left_join(metacell_types %>% select("metacell", "cell_type"), by = "metacell") %>%
        left_join(cell_type_colors %>% select(cell_type, color), by = "cell_type")

    if (plot_legend) {
        legend <- cowplot::get_legend(cell_type_colors %>%
            ggplot(aes(x = cell_type, color = cell_type, y = 1)) +
            geom_point() +
            scale_color_manual("", values = deframe(cell_type_colors[, 1:2])) +
            guides(color = guide_legend(ncol = 1)))
        p_mat <- p_mat + theme(legend.position = "top")
        cowplot::ggdraw(cowplot::plot_grid(p_mat %>% tgplot_add_axis_annotation(mc_types$color), legend, nrow = 1, rel_widths = c(0.8, 0.15)))
    } else {
        cowplot::ggdraw(p_mat %>% tgplot_add_axis_annotation(mc_types$color))
    }
}
