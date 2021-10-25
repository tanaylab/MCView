get_md_attribute <- function(dataset, md, attr, default) {
    metadata_colors <- get_mc_data(dataset, "metadata_colors")
    if (has_name(metadata_colors, md)) {
        md_attr <- metadata_colors[[md]][[attr]]
        if (!is.null(md_attr)) {
            return(md_attr)
        }
    }

    return(default)
}


get_metadata_colors <- function(dataset, md, colors = NULL, color_breaks = NULL, metadata = NULL, default_colors = c("white", "white", "white", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F", "black")) {
    colors <- colors %||% get_md_attribute(dataset, md, "colors", default_colors)
    color_breaks <- color_breaks %||% get_md_attribute(dataset, md, "breaks", NULL)

    if (is.null(color_breaks)) {
        metadata <- metadata %||% get_mc_data(dataset, "metadata")
        min_val <- min(metadata[[md]], na.rm = TRUE)
        max_val <- max(metadata[[md]], na.rm = TRUE)

        if (min_val == max_val) {
            min_val <- min_val - 1e-5
        }

        color_breaks <- seq(min_val, max_val, length.out = length(colors))
    }

    return(
        list(colors = colors, breaks = color_breaks)
    )
}


#' Plot 2d projection of mc2d colored by metadata
#'
#' @param dataset name of metacell object
#' @param md name of the metadata field
#'
#' @noRd
mc2d_plot_metadata_ggp <- function(dataset,
                                   md,
                                   colors = NULL,
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

    md_colors <- get_metadata_colors(dataset, md, colors = colors, color_breaks = color_breaks, metadata = metadata)
    palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)


    mc2d_df <- mc2d_df %>%
        mutate(col_x = palette(.[[md]]))

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
        scale_color_gradientn(name = md, colors = md_colors$colors, values = scales::rescale(md_colors$breaks)) +
        scale_fill_identity()

    return(p)
}


#' Plot gene expression vs time scatter of metacells
#'
#' @param dataset name of metacell object
#' @param md1 name of the first metadata field
#' @param md2 name of the second metadata field
#' @param color_by_md name of the metadata field to color by
#'
#' @noRd
plot_md_md_mc <- function(dataset, md1, md2, color_by_md = NULL, colors = NULL, color_breaks = NULL, metacell_types = get_mc_data(dataset, "metacell_types"), cell_type_colors = get_mc_data(dataset, "cell_type_colors"), point_size = initial_scatters_point_size(dataset), stroke = initial_scatters_stroke(dataset), plot_text = TRUE) {
    metadata <- get_mc_data(dataset, "metadata") %>% mutate(metacell = as.character(metacell))
    metadata_colors <- get_mc_data(dataset, "metadata_colors")

    if (!is.null(color_by_md)) {
        metadata <- metadata %>% select(metacell, !!md1, !!md2, !!color_by_md)
    } else {
        metadata <- metadata %>% select(metacell, !!md1, !!md2)
    }

    df <- metacell_types %>%
        left_join(metadata, by = "metacell") %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        mutate(color = ifelse(!is.null(color_by_md), color_by_md, cell_type)) %>%
        rename(
            `Cell type` = cell_type
        )

    if (!is.null(color_by_md)) {
        md_colors <- get_metadata_colors(dataset, color_by_md, colors = colors, color_breaks = color_breaks, metadata = metadata)
        palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
        df$color <- palette(df[[color_by_md]])
        df$color_values <- df[[color_by_md]]
    } else {
        df$color <- df[["Cell type"]]
        df$color_values <- df$color
    }

    # We call the text field "Metacell" in order for plotly to show "Metacell:" in the tooltip
    md1_name <- md1
    md2_name <- md2

    if (!is.null(color_by_md)) {
        color_by_str <- glue("{color_by_md}: {values}", values = round(df[[color_by_md]], digits = 3))
    } else {
        color_by_str <- ""
    }

    df <- df %>%
        mutate(
            Metacell = paste(
                glue("{metacell}"),
                glue("{md1_name}: {md1_values}", md1_values = round(!!sym(md1), digits = 3)),
                glue("{md2_name}: {md2_values}", md2_values = round(!!sym(md2), digits = 3)),
                color_by_str,
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                ifelse(has_name(df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    p <- ggplot(
        data = df,
        aes(
            x = !!sym(md1),
            y = !!sym(md2),
            fill = color,
            color = color_values,
            label = metacell,
            customdata = metacell,
            tooltip_text = Metacell
        )
    ) +
        xlab(md1) +
        ylab(md2)

    if (is.null(color_by_md)) {
        col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)
        p <- p +
            geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
            scale_fill_manual(values = col_to_ct) +
            guides(color = "none")
    } else {
        p <- p +
            geom_point(size = point_size) +
            geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
            guides(fill = "none")

        p <- p +
            scale_color_gradientn(name = color_by_md, colors = md_colors$colors, values = scales::rescale(md_colors$breaks)) +
            scale_fill_identity()
    }

    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    return(p)
}
