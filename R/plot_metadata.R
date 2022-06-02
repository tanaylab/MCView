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
                                   scale_edges = FALSE,
                                   metacell_types = NULL,
                                   atlas = FALSE,
                                   metadata = NULL) {
    mc2d <- get_mc_data(dataset, "mc2d", atlas = atlas)
    metadata <- metadata %||% get_mc_data(dataset, "metadata", atlas = atlas)

    metadata <- metadata %>% mutate(metacell = as.character(metacell))
    metacell_types <- metacell_types %||% get_mc_data(dataset, "metacell_types")

    metacell_types <- metacell_types %>%
        select(metacell, cell_type, top1_gene, top2_gene, top1_lfp, top2_lfp, mc_col)

    mc2d_df <- mc2d_to_df(mc2d) %>%
        left_join(metacell_types, by = "metacell") %>%
        left_join(metadata %>% select(metacell, !!md), by = "metacell") %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        )

    if (md != "Cell type") {
        mc2d_df <- mc2d_df %>% rename(
            `Cell type` = cell_type
        )
    }

    if (has_name(df, "mc_age")) {
        mc2d_df <- mc2d_df %>% rename(`Age` = mc_age)
    }

    graph <- mc2d_to_graph_df(mc2d, min_d = min_d)

    if (is.null(id)) {
        mc2d_df <- mc2d_df %>% mutate(id = metacell)
    } else {
        mc2d_df <- mc2d_df %>% mutate(id = paste(id, metacell, sep = "\t"))
    }

    if (is_numeric_field(mc2d_df, md)) {
        p <- mc2d_plot_metadata_ggp_numeric(mc2d_df, graph, dataset, metadata, md, colors, color_breaks, point_size, min_d, stroke, graph_color, graph_width, id, scale_edges)
    } else {
        p <- mc2d_plot_metadata_ggp_categorical(mc2d_df, graph, dataset, metadata, md, point_size, min_d, stroke, graph_color, graph_width, id, scale_edges, colors) %>%
            plotly::ggplotly(tooltip = "tooltip_text", source = source) %>%
            rm_plotly_grid()
    }

    return(p)
}

mc2d_plot_metadata_ggp_categorical <- function(mc2d_df,
                                               graph,
                                               dataset,
                                               metadata,
                                               md,
                                               point_size,
                                               min_d,
                                               stroke,
                                               graph_color,
                                               graph_width,
                                               id,
                                               scale_edges,
                                               colors = NULL) {
    mc2d_df <- mc2d_df %>%
        mutate(
            Metacell = paste(
                glue("{metacell}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                ifelse(md != "Cell type", paste0(md, ": ", mc2d_df[[md]]), ""),
                ifelse(has_name(mc2d_df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    if (is.null(colors)) {
        colors <- colors %||% get_metadata_colors(dataset, md, metadata = metadata)
    }

    p <- mc2d_df %>%
        ggplot(aes(x = x, y = y, label = metacell, fill = !!sym(md), tooltip_text = Metacell, customdata = id))

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

    p <- p +
        geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
        theme_void() +
        guides(fill = "none")

    p <- p +
        scale_fill_manual(name = md, values = colors)

    return(p)
}

mc2d_plot_metadata_ggp_numeric <- function(mc2d_df,
                                           graph,
                                           dataset,
                                           metadata,
                                           md,
                                           colors,
                                           color_breaks,
                                           point_size,
                                           min_d,
                                           stroke,
                                           graph_color,
                                           graph_width,
                                           id,
                                           scale_edges) {
    mc2d_df <- mc2d_df %>%
        mutate(
            Metacell = paste(
                glue("{metacell}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                paste0(md, ": ", round(mc2d_df[[md]], digits = 3)),
                ifelse(has_name(mc2d_df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    md_colors <- get_metadata_colors(dataset, md, colors = colors, color_breaks = color_breaks, metadata = metadata)
    palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)

    mc2d_df <- mc2d_df %>%
        mutate(col_x = palette(.[[md]])) %>%
        arrange(desc(!!sym(md)))

    p <- mc2d_df %>%
        ggplot(aes(x = x, y = y, label = metacell, fill = col_x, color = !!sym(md), tooltip_text = Metacell, customdata = id))

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
        scale_color_gradientn(name = md, colors = md_colors$colors, values = scales::rescale(md_colors$breaks, c(0, 1)), breaks = round(md_colors$breaks, digits = 2)) +
        scale_fill_identity()

    return(p)
}


plot_mc_scatter <- function(dataset,
                            x_var,
                            y_var,
                            color_var = NULL,
                            gene_modules = NULL,
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
                            plot_text = TRUE,
                            atlas = FALSE,
                            metadata = get_mc_data(dataset, "metadata", atlas = atlas),
                            x_limits = NULL,
                            y_limits = NULL,
                            fixed_limits = FALSE,
                            xyline = FALSE,
                            metacell_filter = NULL) {
    if (!is.null(metadata)) {
        metadata <- metadata %>% mutate(metacell = as.character(metacell))
    }
    metadata_colors <- get_mc_data(dataset, "metadata_colors", atlas = atlas)

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
        req(metadata)
        df <- df %>%
            select(-any_of(x_var)) %>%
            left_join(metadata %>% select(metacell, !!x_var), by = "metacell") %>%
            mutate(x_str = glue("{x_name}: {x_values}", x_values = round(!!sym(x_var), digits = 3)))
    } else {
        if (x_type == "Gene module") {
            egc_x <- get_gene_module_egc(x_var, dataset, gene_modules, atlas = atlas) + egc_epsilon
        } else {
            egc_x <- get_gene_egc(x_var, dataset, atlas = atlas) + egc_epsilon
        }
        df <- df %>%
            mutate(!!x_var := egc_x[metacell]) %>%
            mutate(x_str = glue("{x_name} expression: {expr_text}", expr_text = scales::scientific(!!sym(x_var))))
    }

    # set y variable
    y_name <- y_var
    if (y_type == "Metadata") {
        req(metadata)
        df <- df %>%
            select(-any_of(y_var)) %>%
            left_join(metadata %>% select(metacell, !!y_var), by = "metacell") %>%
            mutate(y_str = glue("{y_name}: {y_values}", y_values = round(!!sym(y_var), digits = 3)))
    } else {
        if (y_type == "Gene module") {
            egc_y <- get_gene_module_egc(y_var, dataset, gene_modules, atlas = atlas) + egc_epsilon
        } else {
            egc_y <- get_gene_egc(y_var, dataset, atlas = atlas) + egc_epsilon
        }

        df <- df %>%
            mutate(!!y_var := egc_y[metacell]) %>%
            mutate(y_str = glue("{y_name} expression: {expr_text}", expr_text = scales::scientific(!!sym(y_var))))
    }

    # set color variable
    color_name <- color_var
    categorical_md <- FALSE
    if (is.null(color_var)) {
        df <- df %>%
            mutate(color = cell_type, color_values = cell_type) %>%
            mutate(color_str = glue("Cell type: {`Cell type`}"))
    } else if (color_type == "Metadata") {
        req(metadata)
        df <- df %>%
            select(-any_of(color_var)) %>%
            left_join(metadata %>% select(metacell, !!color_var), by = "metacell")
        md_colors <- get_metadata_colors(dataset, color_var, colors = colors, color_breaks = color_breaks, metadata = metadata)
        if (is_numeric_field(metadata, color_var)) {
            palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
            df$color <- palette(df[[color_var]])
            df$color_values <- df[[color_var]]
            df <- df %>%
                mutate(color_str = glue("{color_name}: {color_values}\nCell type: {`Cell type`}", color_values = round(!!sym(color_var), digits = 3)))
        } else {
            df <- df %>%
                mutate(color = !!sym(color_var), color_values = !!sym(color_var)) %>%
                mutate(color_str = glue("{color_name}: {color_values}"))
            categorical_md <- TRUE
        }
    } else if (color_type %in% c("Gene", "Gene module")) {
        if (color_type == "Gene module") {
            egc_color <- get_gene_module_egc(color_var, dataset, gene_modules, atlas = atlas) + egc_epsilon
        } else {
            egc_color <- get_gene_egc(color_var, dataset, atlas = atlas) + egc_epsilon
        }
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
            mutate(color_str = glue("{color_name}: {color_values}\nCell type: {`Cell type`}\n", color_values = round(expression, digits = 3)))
    }

    # set tooltip
    df <- df %>%
        mutate(
            Metacell = paste0(
                glue("{metacell}\n{x_str}\n{y_str}\n{color_str}\nTop genes: {`Top genes`}\n"),
                ifelse(has_name(df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), "")
            )
        )

    if (!is.null(metacell_filter)) {
        df <- df %>%
            filter(metacell %in% metacell_filter)
    }

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

    if (xyline) {
        p <- p + geom_abline(linetype = "dashed")
    }

    # set color plotting
    if (is.null(color_var)) {
        if (atlas) {
            col_to_ct <- get_cell_type_colors(dataset, NULL, atlas = TRUE)
        } else {
            col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)
        }

        p <- p +
            geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
            scale_fill_manual("", values = col_to_ct) +
            guides(color = "none")
    } else if (categorical_md) {
        p <- p +
            geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
            scale_fill_manual(name = color_var, values = md_colors) +
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
    xylims <- expr_breaks

    if (fixed_limits && x_type %in% c("Gene", "Gene module") && y_type %in% c("Gene", "Gene module")) {
        x_limits <- x_limits %||% c(min(egc_x), max(egc_x))
        y_limits <- y_limits %||% c(min(egc_y), max(egc_y))
        x_limits <- c(min(c(x_limits[1], y_limits[1])), max(c(x_limits[2], y_limits[2])))
        y_limits <- x_limits
    }

    if (x_type %in% c("Gene", "Gene module")) {
        x_limits <- x_limits %||% c(min(egc_x), max(egc_x))
        xmax <- min(c(1:length(xylims))[xylims >= x_limits[2]])
        xmin <- max(c(1:length(xylims))[xylims <= x_limits[1]])
        p <- p +
            scale_x_continuous(limits = c(xylims[xmin], xylims[xmax]), trans = "log2", breaks = xylims[xmin:xmax], labels = scales::scientific(xylims[xmin:xmax])) +
            xlab(glue("{x_var} Expression")) +
            theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))
    }

    if (y_type %in% c("Gene", "Gene module")) {
        y_limits <- y_limits %||% c(min(egc_y), max(egc_y))
        ymax <- min(c(1:length(xylims))[xylims >= y_limits[2]])
        ymin <- max(c(1:length(xylims))[xylims <= y_limits[1]])
        p <- p +
            scale_y_continuous(limits = c(xylims[ymin], xylims[ymax]), trans = "log2", breaks = xylims[ymin:ymax], labels = scales::scientific(xylims[ymin:ymax])) +
            ylab(glue("{y_var} Expression"))
    }


    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    return(p)
}

plot_sample_scatter <- function(dataset,
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
                                cell_types = NULL,
                                point_size = initial_scatters_point_size(dataset),
                                stroke = initial_scatters_stroke(dataset),
                                expr_colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
                                plot_text = TRUE) {
    metadata <- get_mc_data(dataset, "cell_metadata")
    metadata_colors <- get_mc_data(dataset, "metadata_colors")

    req(metadata)
    req(metadata$samp_id)
    req(metadata$metacell)

    if (!is.null(cell_types)) {
        selected_mc <- metacell_types %>%
            select(metacell, cell_type) %>%
            filter(cell_type %in% cell_types) %>%
            pull(metacell)
    } else {
        selected_mc <- metacell_types$metacell
    }

    df <- metadata %>%
        distinct(samp_id, .keep_all = TRUE)


    if (any(c(x_type, y_type, color_type) == "Cell type")) {
        cell_type_fracs <- metadata %>%
            mutate(metacell = as.character(metacell)) %>%
            left_join(metacell_types %>% select(metacell, cell_type), by = "metacell") %>%
            count(samp_id, cell_type) %>%
            group_by(samp_id) %>%
            mutate(frac = n / sum(n)) %>%
            ungroup()
    }

    if (any(c(x_type, y_type, color_type) == "Metadata")) {
        samp_md <- get_samp_metadata(dataset)
    }

    get_var_md <- function(var, var_name, str_name = NULL) {
        if (var %in% colnames(samp_md)) {
            res <- samp_md %>%
                select(samp_id, !!var)
        } else {
            res <- metadata %>%
                select(samp_id, !!var) %>%
                group_by(samp_id) %>%
                summarise(!!var := mean(!!sym(var)))
        }
        if (!is.null(str_name)) {
            res <- res %>%
                mutate(!!str_name := glue("{var_name}: {values}", values = round(!!sym(var), digits = 3)))
        }
        return(res)
    }

    # set x variable
    x_name <- x_var
    if (x_type == "Metadata") {
        req(x_var %in% colnames(metadata))
        df <- df %>%
            select(-any_of(x_var)) %>%
            left_join(get_var_md(x_var, x_name, "x_str"), by = "samp_id")
    } else if (x_type == "Gene") {
        req(x_var %in% gene_names(dataset))
        egc_x <- get_samples_gene_egc(x_var, dataset, selected_mc) + egc_epsilon
        df <- df %>%
            mutate(!!x_var := egc_x[df$samp_id]) %>%
            mutate(x_str = glue("{x_name} expression (log2): {expr_text}", expr_text = round(log2(!!sym(x_var)), digits = 2)))
    } else {
        req(x_var %in% cell_type_colors$cell_type)
        df <- cell_type_fracs %>%
            filter(cell_type == !!x_var) %>%
            mutate(!!x_var := frac) %>%
            mutate(x_str = glue("{x_name}: {x_values}", x_values = scales::percent(!!sym(x_var))))
    }

    # set y variable
    y_name <- y_var
    if (y_type == "Metadata") {
        req(y_var %in% colnames(metadata))
        df <- df %>%
            select(-any_of(y_var)) %>%
            left_join(get_var_md(y_var, y_name, "y_str"), by = "samp_id")
    } else if (y_type == "Gene") {
        req(y_var %in% gene_names(dataset))
        egc_y <- get_samples_gene_egc(y_var, dataset, selected_mc) + egc_epsilon
        df <- df %>%
            mutate(!!y_var := egc_y[df$samp_id]) %>%
            mutate(y_str = glue("{y_name} expression (log2): {expr_text}", expr_text = round(log2(!!sym(y_var)), digits = 2)))
    } else {
        req(y_var %in% cell_type_colors$cell_type)
        y_df <- cell_type_fracs %>%
            filter(cell_type == !!y_var) %>%
            mutate(!!y_var := frac)

        df <- df %>%
            select(-any_of(y_var)) %>%
            left_join(y_df %>% select(samp_id, !!y_var), by = "samp_id") %>%
            mutate(y_str = glue("{y_name}: {y_values}", y_values = scales::percent(!!sym(y_var))))
    }


    # set color variable
    color_name <- color_var
    color_var_type <- "cont"
    if (!is.null(color_var) && color_var != "None") {
        if (color_type == "Metadata") {
            req(color_var %in% colnames(metadata))
            df <- df %>%
                select(-any_of(color_var)) %>%
                left_join(get_var_md(color_var, color_name, NULL), by = "samp_id")

            if (is_numeric_field(df, color_var)) {
                md_colors <- get_metadata_colors(dataset, color_var, colors = colors, color_breaks = color_breaks, metadata = metadata)
                palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
                df$color <- palette(df[[color_var]])
                df$color_values <- df[[color_var]]
                df <- df %>%
                    mutate(color_str = glue("{color_name}: {color_values}", color_values = round(!!sym(color_var), digits = 3)))
            } else {
                metadata_colors <- get_mc_data(dataset, "metadata_colors")
                if (is.null(metadata_colors[[color_var]])) {
                    categories <- unique(df[[color_var]])
                    colors <- chameleon::distinct_colors(length(categories))$name
                    names(colors) <- categories
                } else {
                    colors <- metadata_colors[[color_var]]
                }

                df <- df %>%
                    mutate(color_str = paste0(color_name, ": ", !!sym(color_var)))
                color_var_type <- "discrete"
            }
        } else if (color_type == "Gene") {
            req(color_var %in% gene_names(dataset))
            egc_color <- get_samples_gene_egc(color_var, dataset, selected_mc) + egc_epsilon
            df <- df %>%
                mutate(expression = log2(egc_color[df$samp_id]))
            min_expr <- min(df$expression, na.rm = TRUE)
            max_expr <- max(df$expression, na.rm = TRUE)
            if (min_expr == max_expr) {
                min_val <- min_val - 1e-5
            }

            color_breaks <- seq(min_expr, max_expr, length.out = length(expr_colors))
            md_colors <- list(colors = expr_colors, breaks = color_breaks)
            palette <- circlize::colorRamp2(colors = expr_colors, breaks = color_breaks)
            df$color <- palette(df$expression)
            df$color_values <- df$expression

            df <- df %>%
                mutate(color_str = glue("{color_name} (log2): {color_values}\n", color_values = round(expression, digits = 3)))
        } else {
            req(color_var %in% cell_type_colors$cell_type)
            color_var_df <- cell_type_fracs %>%
                filter(cell_type == !!color_var) %>%
                mutate(!!color_var := frac * 100)

            df <- df %>%
                select(-any_of(color_var)) %>%
                left_join(color_var_df %>% select(samp_id, !!color_var), by = "samp_id")

            md_colors <- get_metadata_colors(dataset, color_var, colors = colors, color_breaks = color_breaks, metadata = df)
            palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
            df$color <- palette(df[[color_var]])
            df$color_values <- df[[color_var]]

            df <- df %>%
                mutate(color_str = glue("{color_name}: {color_values}", color_values = scales::percent(!!sym(color_var))))
        }
    } else {
        df <- df %>%
            mutate(color_str = "")
    }


    # set tooltip
    df <- df %>%
        mutate(
            Sample = paste0(
                glue("{samp_id}\n{x_str}\n{y_str}\n{color_str}\n")
            )
        )


    # set color plotting
    if (!is.null(color_var) && color_var != "None") {
        if (color_var_type == "cont") {
            p <- ggplot(
                data = df,
                aes(
                    x = !!sym(x_var),
                    y = !!sym(y_var),
                    fill = color,
                    color = color_values,
                    label = samp_id,
                    customdata = samp_id,
                    tooltip_text = Sample
                )
            ) +
                geom_point(size = point_size) +
                geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
                guides(fill = "none")
            p <- p +
                scale_color_gradientn(name = color_var, colors = md_colors$colors, values = scales::rescale(md_colors$breaks)) +
                scale_fill_identity()
        } else {
            p <- ggplot(
                data = df,
                aes(
                    x = !!sym(x_var),
                    y = !!sym(y_var),
                    fill = !!sym(color_var),
                    label = samp_id,
                    customdata = samp_id,
                    tooltip_text = Sample
                )
            ) +
                geom_point(size = point_size) +
                scale_fill_manual(values = colors)
        }
    } else {
        p <- ggplot(
            data = df,
            aes(
                x = !!sym(x_var),
                y = !!sym(y_var),
                label = samp_id,
                customdata = samp_id,
                tooltip_text = Sample
            )
        ) +
            geom_point(size = point_size, color = "black")
    }

    p <- p +
        xlab(x_var) +
        ylab(y_var)

    # arrange axis for gene expression
    xylims <- expr_breaks

    if (x_type %in% c("Gene", "Gene module")) {
        xmax <- min(c(1:length(xylims))[xylims >= max(egc_x)])
        xmin <- max(c(1:length(xylims))[xylims <= min(egc_x)])
        p <- p +
            scale_x_continuous(limits = c(xylims[xmin], xylims[xmax]), trans = "log2", breaks = xylims[xmin:xmax], labels = scales::scientific(xylims[xmin:xmax])) +
            xlab(glue("{x_var} Expression")) +
            theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))
    } else if (x_type == "Cell type") {
        p <- p +
            scale_x_continuous(labels = scales::percent)
    }

    if (y_type %in% c("Gene", "Gene module")) {
        ymax <- min(c(1:length(xylims))[xylims >= max(egc_y)])
        ymin <- max(c(1:length(xylims))[xylims <= min(egc_y)])
        p <- p +
            scale_y_continuous(limits = c(xylims[ymin], xylims[ymax]), trans = "log2", breaks = xylims[ymin:ymax], labels = scales::scientific(xylims[ymin:ymax])) +
            ylab(glue("{y_var} Expression"))
    } else if (y_type == "Cell type") {
        p <- p +
            scale_y_continuous(labels = scales::percent)
    }


    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    return(p)
}


plot_obs_proj_scatter <- function(dataset,
                                  axis_var,
                                  color_var = NULL,
                                  axis_type = "Metadata",
                                  color_type = NULL,
                                  colors = NULL,
                                  color_breaks = NULL,
                                  metacell_types = get_mc_data(dataset, "metacell_types"),
                                  cell_type_colors = get_mc_data(dataset, "cell_type_colors"),
                                  cell_types = NULL,
                                  point_size = initial_scatters_point_size(dataset),
                                  stroke = initial_scatters_stroke(dataset),
                                  expr_colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
                                  plot_text = TRUE) {
    atlas_metadata <- get_mc_data(dataset, "metadata", atlas = TRUE)
    query_metadata <- get_mc_data(dataset, "metadata", atlas = FALSE)
    if (!is.null(atlas_metadata)) {
        atlas_metadata <- atlas_metadata %>% mutate(metacell = as.character(metacell))
    }
    if (!is.null(query_metadata)) {
        query_metadata <- query_metadata %>% mutate(metacell = as.character(metacell))
    }

    df <- metacell_types %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        mutate(`Cell type` = cell_type)

    # set axis variables
    axis_name <- axis_var
    if (axis_type == "Metadata") {
        req(atlas_metadata)
        proj_w <- get_mc_data(dataset, "proj_weights")
        req(proj_w)
        # not implemented yet
        req(FALSE)

        # df <- df %>%
        #     select(-any_of(axis_var)) %>%
        #     left_join(metadata %>% select(metacell, !!x_var), by = "metacell") %>%
        #     mutate(x_str = glue("{x_name}: {x_values}", x_values = round(!!sym(x_var), digits = 3)))
    } else {
        egc_obs <- get_gene_egc(axis_var, dataset) + egc_epsilon
        egc_proj <- get_gene_egc(axis_var, dataset, projected = TRUE) + egc_epsilon
        x_var <- glue("{axis_var} - observed")
        y_var <- glue("{axis_var} - projected")
        df <- df %>%
            mutate(!!x_var := egc_obs[metacell], !!y_var := egc_proj[metacell]) %>%
            mutate(x_str = glue("{axis_name} obs: {expr_text}", expr_text = scales::scientific(!!sym(x_var)))) %>%
            mutate(y_str = glue("{axis_name} proj: {expr_text}", expr_text = scales::scientific(!!sym(x_var))))
    }

    categorical_md <- FALSE
    color_name <- color_var
    if (is.null(color_var)) {
        df <- df %>%
            mutate(color = cell_type, color_values = cell_type) %>%
            mutate(color_str = glue("Cell type: {`Cell type`}"))
    } else if (color_type == "Metadata") {
        if (grepl("_atlas$", color_var) && !is.null(atlas_metadata) && has_name(atlas_metadata, sub("_atlas$", "", color_var))) {
            req(atlas_metadata)
            proj_w <- get_mc_data(dataset, "proj_weights")
            req(proj_w)
            color_var <- sub("_atlas$", "", color_var)
            color_name <- color_var
            proj_md <- proj_w %>%
                left_join(
                    atlas_metadata %>%
                        select(atlas = metacell, !!color_var),
                    by = "atlas"
                ) %>%
                group_by(query) %>%
                summarise(!!color_var := sum(weight * !!sym(color_var))) %>%
                rename(metacell = query)

            df <- df %>%
                select(-any_of(color_var)) %>%
                left_join(proj_md, by = "metacell")

            md_colors <- get_metadata_colors(dataset, color_var, colors = colors, color_breaks = color_breaks, metadata = atlas_metadata, atlas = TRUE)
        } else {
            req(query_metadata)
            df <- df %>%
                select(-any_of(color_var)) %>%
                left_join(query_metadata %>% select(metacell, !!color_var), by = "metacell")
            md_colors <- get_metadata_colors(dataset, color_var, colors = colors, color_breaks = color_breaks, metadata = query_metadata)
            categorical_md <- !is_numeric_field(query_metadata, color_var)
        }

        if (categorical_md) {
            df <- df %>%
                mutate(color = !!sym(color_var), color_values = !!sym(color_var)) %>%
                mutate(color_str = glue("{color_name}: {color_values}"))
        } else {
            palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
            df$color <- palette(df[[color_var]])
            df$color_values <- df[[color_var]]
            df <- df %>%
                mutate(color_str = glue("{color_name}: {color_values}\nCell type: {`Cell type`}", color_values = round(!!sym(color_var), digits = 3)))
        }
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
            mutate(color_str = glue("{color_name}: {color_values}\nCell type: {`Cell type`}\n", color_values = round(expression, digits = 3)))
    }

    # set tooltip
    df <- df %>%
        mutate(
            Metacell = paste0(
                glue("{metacell}\n{x_str}\n{y_str}\n{color_str}\nTop genes: {`Top genes`}\n"),
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
        ylab(y_var) +
        geom_abline(linetype = "dashed")

    # set color plotting
    if (is.null(color_var)) {
        col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)
        p <- p +
            geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
            scale_fill_manual(values = col_to_ct) +
            guides(color = "none")
    } else if (categorical_md) {
        p <- p +
            geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
            scale_fill_manual(name = color_var, values = md_colors) +
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
    xylims <- expr_breaks

    if (axis_type %in% c("Gene", "Gene module")) {
        xmax <- min(c(1:length(xylims))[xylims >= max(egc_obs)])
        xmin <- max(c(1:length(xylims))[xylims <= min(egc_obs)])
        ymax <- min(c(1:length(xylims))[xylims >= max(egc_proj)])
        ymin <- max(c(1:length(xylims))[xylims <= min(egc_proj)])
        p <- p +
            scale_x_continuous(limits = c(xylims[xmin], xylims[xmax]), trans = "log2", breaks = xylims[xmin:xmax], labels = scales::scientific(xylims[xmin:xmax])) +
            xlab(glue("{x_var} Expression")) +
            scale_y_continuous(limits = c(xylims[ymin], xylims[ymax]), trans = "log2", breaks = xylims[ymin:ymax], labels = scales::scientific(xylims[ymin:ymax])) +
            ylab(glue("{y_var} Expression")) +
            theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))
    }

    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    return(p)
}
