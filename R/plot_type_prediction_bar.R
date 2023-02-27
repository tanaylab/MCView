plot_type_predictions_bar <- function(dataset, metacell_types, cell_type_colors) {
    renderPlot(
        {
            df_fracs <- get_mc_data(dataset(), "query_atlas_cell_type_fracs")
            req(!is.null(df_fracs))
            req(has_atlas(dataset()))
            atlas_colors <- get_mc_data(dataset(), "cell_type_colors", atlas = TRUE)
            atlas_colors <- atlas_colors %>%
                select(cell_type, color) %>%
                deframe()

            fracs_mat <- df_fracs %>%
                spread(type, fraction) %>%
                column_to_rownames("metacell") %>%
                as.matrix()

            ord <- slanter::slanted_orders(fracs_mat)

            df_fracs <- df_fracs %>%
                mutate(metacell = factor(metacell, levels = rownames(fracs_mat)[ord$rows])) %>%
                mutate(type = factor(type, levels = colnames(fracs_mat)[ord$cols]))

            p_fracs <- df_fracs %>%
                ggplot(aes(x = metacell, y = fraction, fill = type, customdata = metacell)) +
                geom_col() +
                scale_fill_manual(name = "", values = atlas_colors) +
                theme(
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank()
                ) +
                scale_y_continuous(expand = c(0, 0)) +
                ylab("Fraction") +
                xlab("Metacell") +
                guides(fill = "none")

            cell_type_colors <- cell_type_colors() %>%
                select(cell_type, color) %>%
                deframe()

            p_annot <- metacell_types() %>%
                mutate(metacell = factor(metacell, levels = rownames(fracs_mat)[ord$rows])) %>%
                ggplot(aes(x = metacell, y = 1, fill = cell_type)) +
                geom_raster() +
                scale_fill_manual(name = "", values = cell_type_colors) +
                theme_void() +
                ggplot2::theme(panel.border = ggplot2::element_rect(size = 0.2, color = "black", fill = NA))

            p <- cowplot::insert_xaxis_grob(p_fracs, p_annot, grid::unit(0.05, "null"), position = "bottom")
            cowplot::ggdraw(p)
        },
        res = 92
    )
}
