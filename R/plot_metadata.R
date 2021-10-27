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



plot_mc_scatter <- function(dataset,
                            x_var,
                            y_var,
                            color_var = NULL,
                            x_type = "Metadata",
                            y_type = "Metadata",
                            color_type = NULL,
                            colors = NULL,
                            color_breaks = NULL,
                            metacell_types = get_mc_data(dataset, "metacell_types"),
                            cell_type_colors = get_mc_data(dataset, "cell_type_colors"),
                            point_size = initial_scatters_point_size(dataset),
                            stroke = initial_scatters_stroke(dataset),
                            expr_colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
                            plot_text = TRUE) {
    metadata <- get_mc_data(dataset, "metadata") %>% mutate(metacell = as.character(metacell))
    metadata_colors <- get_mc_data(dataset, "metadata_colors")

    df <- metacell_types %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        mutate(`Cell type` = cell_type)

    # set x variable
    x_name <- x_var
    if (x_type == "Metadata") {
        df <- df %>%
            select(-any_of(x_var)) %>%
            left_join(metadata %>% select(metacell, !!x_var), by = "metacell") %>%
            mutate(x_str = glue("{x_name}: {x_values}", x_values = round(!!sym(x_var), digits = 3)))
    } else {
        egc_x <- get_gene_egc(x_var, dataset) + egc_epsilon
        df <- df %>%
            mutate(!!x_var := egc_x) %>%
            mutate(x_str = glue("{x_name} expression: {expr_text}", expr_text = scales::scientific(!!sym(x_var))))
    }

    # set y variable
    y_name <- y_var
    if (y_type == "Metadata") {
        df <- df %>%
            select(-any_of(y_var)) %>%
            left_join(metadata %>% select(metacell, !!y_var), by = "metacell") %>%
            mutate(y_str = glue("{y_name}: {y_values}", y_values = round(!!sym(x_var), digits = 3)))
    } else {
        egc_y <- get_gene_egc(y_var, dataset) + egc_epsilon
        df <- df %>%
            mutate(!!y_var := egc_y) %>%
            mutate(y_str = glue("{y_name} expression: {expr_text}", expr_text = scales::scientific(!!sym(y_var))))
    }

    # set color variable
    color_name <- color_var
    if (is.null(color_var)) {
        df <- df %>%
            mutate(color = cell_type, color_values = cell_type) %>%
            mutate(color_str = glue("Cell type: {`Cell type`}"))
    } else if (color_type == "Metadata") {
        df <- df %>%
            select(-any_of(color_var)) %>%
            left_join(metadata %>% select(metacell, !!color_var), by = "metacell")
        md_colors <- get_metadata_colors(dataset, color_var, colors = colors, color_breaks = color_breaks, metadata = metadata)
        palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
        df$color <- palette(df[[color_var]])
        df$color_values <- df[[color_var]]
        df <- df %>%
            mutate(color_str = glue("{color_name}: {color_values}\nCell type: {`Cell type`}", color_values = round(!!sym(color_var), digits = 3)))
    } else if (color_type == "Gene") {
        egc_color <- get_gene_egc(color_var, dataset) + egc_epsilon
        df <- df %>%
            mutate(expression = log2(egc_color[df$metacell]))
        min_expr <- min(df$expression, na.rm = TRUE)
        max_expr <- max(df$expression, na.rm = TRUE)

        color_breaks <- seq(min_expr, max_expr, length.out = length(expr_colors))
        md_colors <- list(colors = expr_colors, breaks = color_breaks)
        palette <- circlize::colorRamp2(colors = expr_colors, breaks = color_breaks)
        df$color <- palette(df$expression)
        df$color_values <- df$expression

        df <- df %>%
            mutate(color_str = glue("{color_name}: {color_values}\nCell type: {`Cell type`}", color_values = round(expression, digits = 3)))
    }

    # set tooltip
    df <- df %>%
        mutate(
            Metacell = paste0(
                glue("{metacell}\n{x_str}\n{y_str}\n{color_str}\n"),
                glue("Top genes: {`Top genes`}\n"),
                ifelse(has_name(df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), "")
            )
        )


    p <- ggplot(
        data = df,
        aes(
            x = !!sym(x_var),
            y = !!sym(y_var),
            fill = color,
            color = color_values,
            label = metacell,
            customdata = metacell,
            tooltip_text = Metacell
        )
    ) +
        xlab(x_var) +
        ylab(y_var)

    # set color plotting
    if (is.null(color_var)) {
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
            scale_color_gradientn(name = color_var, colors = md_colors$colors, values = scales::rescale(md_colors$breaks)) +
            scale_fill_identity()
    }

    # arrange axis for gene expression
    xylims <- c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1)

    if (x_type == "Gene") {
        xmax <- min(c(1:length(xylims))[xylims >= max(egc_x)])
        xmin <- max(c(1:length(xylims))[xylims <= min(egc_x)])
        p <- p +
            scale_x_continuous(limits = c(xylims[xmin], xylims[xmax]), trans = "log2", breaks = xylims[xmin:xmax], labels = scales::scientific(xylims[xmin:xmax])) +
            xlab(glue("{x_var} Expression")) +
            theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))
    }

    if (y_type == "Gene") {
        ymax <- min(c(1:length(xylims))[xylims >= max(egc_y)])
        ymin <- max(c(1:length(xylims))[xylims <= min(egc_y)])
        p <- p +
            scale_y_continuous(limits = c(xylims[ymin], xylims[ymax]), trans = "log2", breaks = xylims[ymin:ymax], labels = scales::scientific(xylims[ymin:ymax])) +
            ylab(glue("{y_var} Expression"))
    }


    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    return(p)
}
