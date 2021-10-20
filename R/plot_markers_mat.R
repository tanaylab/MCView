plot_markers_mat <- function(mc_fp,
                             metacell_types,
                             cell_type_colors,
                             colors = c("darkblue", "blue", "lightblue", "white", "red", "darkred"),
                             values = c(-3, -2, -1, 0, 1.5, 3)) {
    min_lfp <- min(values)
    max_lfp <- max(values)

    gene_ord <- order(apply(mc_fp, 1, which.max))
    mat <- log2(mc_fp[gene_ord, ])

    p_mat <- tgutil::tgplot_heatmap(clip_vals(mat, min_lfp, max_lfp), col_names = FALSE) +
        scale_fill_gradientn(name = "Fold change", colors = colors, values = scales::rescale(values))

    mc_types <- tibble(metacell = colnames(mat)) %>%
        left_join(metacell_types %>% select("metacell", "cell_type"), by = "metacell") %>%
        left_join(cell_type_colors %>% select(cell_type, color), by = "cell_type")

    cowplot::ggdraw(p_mat %>% tgplot_add_axis_annotation(mc_types$color))
}
