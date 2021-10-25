#' Plot 2d projection of mc2d colored by metadata
#'
#' @param dataset name of metacell object
#' @param md name of the metadata field
#'
#' @noRd
mc2d_plot_metadata_ggp <- function(dataset,
                                   md,
                                   colors = c("white", "white", "white", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F", "black"),
                                   color_breaks = NULL,
                                   point_size = initial_proj_point_size(dataset),
                                   min_d = min_edge_length(dataset),
                                   stroke = initial_proj_stroke(dataset),
                                   graph_color = "black",
                                   graph_width = 0.1,
                                   id = NULL,
                                   scale_edges = FALSE) {
    mc2d <- get_mc_data(dataset, "mc2d")
    metadata <- get_mc_data(dataset, "metadata") %>% mutate(metacell = as.character(metacell))
    metadata_colors <- get_mc_data(dataset, "metadata_colors")
    metacell_types <- get_mc_data(dataset, "metacell_types")

    mc2d_df <- mc2d_to_df(mc2d) %>%
        left_join(metacell_types, by = "metacell") %>%
        left_join(metadata %>% select(metacell, !!md), by = "metacell") %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        rename(
            `Cell type` = cell_type
        )

    if (has_name(df, "mc_age")) {
        mc2d_df <- mc2d_df %>% rename(`Age` = mc_age)
    }

    graph <- mc2d_to_graph_df(mc2d, min_d = min_d)

    if (is.null(id)) {
        mc2d_df <- mc2d_df %>% mutate(id = metacell)
    } else {
        mc2d_df <- mc2d_df %>% mutate(id = paste(id, metacell, sep = "\t"))
    }

    mc2d_df <- mc2d_df %>%
        mutate(
            Metacell = paste(
                glue("{metacell}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                paste0(md, ": ", round(.[[md]], digits = 3)),
                ifelse(has_name(df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    min_val <- min(mc2d_df[[md]], na.rm = TRUE)
    max_val <- max(mc2d_df[[md]], na.rm = TRUE)

    if (min_val == max_val) {
        min_val <- min_val - 1e-5
    }

    color_breaks <- color_breaks %||% seq(min_val, max_val, length.out = length(colors))

    palette <- circlize::colorRamp2(colors = colors, breaks = color_breaks)

    mc2d_df <- mc2d_df %>%
        mutate(col_x = palette(!!sym(md)))

    p <- mc2d_df %>%
        ggplot(aes_string(x = "x", y = "y", label = "metacell", fill = "col_x", color = md, tooltip_text = "Metacell", customdata = "id"))

    if (nrow(graph) > 0) {
        if (scale_edges) {
            p <- p +
                geom_segment(data = graph, inherit.aes = FALSE, aes(x = x_mc1, y = y_mc1, xend = x_mc2, yend = y_mc2, size = d_norm), color = graph_color) +
                scale_size_continuous(range = c(0, graph_width)) +
                guides(size = "none")
        } else {
            p <- p + geom_segment(data = graph, inherit.aes = FALSE, aes(x = x_mc1, y = y_mc1, xend = x_mc2, yend = y_mc2), color = graph_color, size = graph_width)
        }
    }

    # Due to https://github.com/ropensci/plotly/issues/1234 we need to plot geom_point twice
    p <- p +
        geom_point(size = point_size) +
        geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
        theme_void() +
        guides(fill = "none")

    p <- p +
        scale_color_gradientn(name = md, colors = colors, values = color_breaks) +
        scale_fill_identity()

    return(p)
}
