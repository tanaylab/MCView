plot_markers_mat <- function(mc_fp,
                             metacell_types,
                             cell_type_colors,
                             min_lfp = NULL,
                             max_lfp = NULL) {
    min_lfp <- min_lfp %||% -3
    max_lfp <- max_lfp %||% 3

    hc_genes <- tgs_dist(mc_fp) %>%
        stats::hclust(method = "ward.D2")

    hc_metacells <- tgs_dist(t(mc_fp)) %>%
        stats::hclust(method = "ward.D2")

    mat <- log2(mc_fp[hc_genes$ord, hc_metacells$ord])

    # text_df <- mat %>%
    #     as.data.frame() %>%
    #     rownames_to_column("gene") %>%
    #     gather("metacell", "value", -gene) %>%
    #     left_join(metacell_types, by = "metacell") %>%
    #     mutate(text = paste(
    #         glue("Metacell: {metacell}"),
    #         glue("Gene: {gene}"),
    #         glue("Footprint (log2): {value}"),
    #         glue("Cell type: {cell_type}"),
    #         sep = "<br>"
    #     ))

    # heatmap_text <- matrix(
    #     unlist(text_df$text),
    #     ncol = ncol(mat),
    #     byrow = FALSE)


    x_axis_options <- list(
        title = "Metacell",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        ticks = ""
    )

    y_axis_options <- list(
        title = "Gene",
        tickfont = list(size = 8)
    )

    heatmap <- plotly::plot_ly(
        y = rownames(mat),
        z = mat,
        zmin = min_lfp,
        zmax = max_lfp,
        zmid = 0, 
        colorscale = list(c(min_lfp, 0, max_lfp), c("#053061", "white", "#67001F")),        
        type = "heatmap" # no rownames are shown in webGL version so we use a regular heatmap
        # hoverinfo = "text", # rendering text was very slow
        # text = heatmap_text
    ) %>% plotly::layout(xaxis = x_axis_options, yaxis = y_axis_options)

    mc_types <- tibble(metacell = colnames(mat)) %>%
        left_join(metacell_types %>% select("metacell", "cell_type"), by = "metacell") %>%
        left_join(cell_type_colors %>% select(cell_type, color), by = "cell_type") %>%
        mutate(Metacell = paste(
            glue("{metacell}"),
            glue("Cell type: {cell_type}"),
            sep = "<br>"
        ))

    # colors_heatmap <- mc_types %>%
    #     ggplot(aes(x=1:nrow(mc_types), y=1, fill=color, label = Metacell)) +
    #     scale_fill_identity() +
    #     coord_cartesian(expand = FALSE) +
    #     geom_raster() +
    #     theme_void() +
    #     guides(fill = "none", color="none") +
    #     theme(panel.border = element_rect(size = 0.2, color = "black", fill = NA))

    # colors_heatmap <- colors_heatmap %>% plotly::ggplotly() %>% plotly::layout(xaxis = x_axis_options, yaxis = x_axis_options) %>% plotly::hide_legend()

    empty_axis_options <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        ticks = ""
    )

    colors_heatmap <- plotly::plot_ly(
        # x = 1:nrow(mc_types),
        z = matrix(as.numeric(mc_types$metacell), nrow = 1, ncol = nrow(mc_types)),
        type = "heatmap",
        hoverinfo = 'text',
        opacity = 1,
        xgap = 0.5,
        text = matrix(
            glue("Metacell: {mc_types$metacell}<br>Cell type:{mc_types$cell_type}"),
            nrow=1,
            ncol=nrow(mc_types)),
        colors = mc_types$color,
        showscale = FALSE
        ) %>% plotly::layout(xaxis = empty_axis_options, yaxis = empty_axis_options)

    fig <- plotly::subplot(
        heatmap,
        colors_heatmap,
        nrows = 2,
        shareX = FALSE,
        shareY = FALSE,
        heights = c(0.95,0.05))

    return(fig)
}
