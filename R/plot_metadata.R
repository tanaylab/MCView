# ---------------------------------------------------------------------------
# Shared helpers for scatter plot functions
# ---------------------------------------------------------------------------

#' Apply gene-expression log2 axis scale to a ggplot (x or y).
#' @noRd
apply_gene_axis_scale <- function(p,
                                  axis = c("x", "y"),
                                  var_name,
                                  egc_values,
                                  xylims,
                                  limits = NULL,
                                  log_labels = FALSE,
                                  corrected = FALSE,
                                  tolerance = 1e-10,
                                  rotate_x = TRUE) {
    axis <- match.arg(axis)
    limits <- limits %||% c(min(egc_values), max(egc_values))

    idx_max <- min(which(xylims >= limits[2] - tolerance))
    idx_min <- max(which(xylims <= limits[1] + tolerance))

    lab <- glue("{var_name} Expression")

    if (log_labels) {
        labels <- log2(xylims[idx_min:idx_max])
        lab <- glue("{lab} (log2)")
    } else {
        labels <- scales::scientific(xylims[idx_min:idx_max])
    }

    if (corrected) {
        lab <- glue("{lab} (corrected)")
    }

    scale_fn <- if (axis == "x") scale_x_continuous else scale_y_continuous
    lab_fn   <- if (axis == "x") xlab else ylab

    p <- p +
        scale_fn(
            limits = c(xylims[idx_min], xylims[idx_max]),
            trans  = "log2",
            breaks = xylims[idx_min:idx_max],
            labels = labels
        ) +
        lab_fn(lab)

    if (axis == "x" && rotate_x) {
        p <- p + theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))
    }

    return(p)
}


#' Apply color/fill geom_point + scale layers to a scatter ggplot.
#' Dispatches on NULL (cell-type), categorical, or continuous color_var.
#' @noRd
apply_scatter_color_layer <- function(p,
                                      color_var,
                                      categorical_md,
                                      point_size,
                                      stroke,
                                      col_to_ct = NULL,
                                      md_colors = NULL,
                                      fill_name = "") {
    if (is.null(color_var)) {
        p <- p +
            geom_point(size = point_size, shape = 21, stroke = stroke, color = "black") +
            scale_fill_manual(fill_name, values = col_to_ct) +
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

    return(p)
}


#' Resolve color for gene expression: log2 palette + color/color_values/color_str columns.
#' Returns list(df, md_colors).
#' @noRd
resolve_gene_color <- function(df,
                               egc_values,
                               id_col,
                               color_name,
                               expr_colors,
                               color_str_suffix = "",
                               guard_equal_range = FALSE) {
    df <- df %>%
        mutate(expression = log2(egc_values[df[[id_col]]]))
    min_expr <- min(df$expression, na.rm = TRUE)
    max_expr <- max(df$expression, na.rm = TRUE)
    if (guard_equal_range && min_expr == max_expr) {
        min_expr <- min_expr - 1e-5
    }

    color_breaks <- seq(min_expr, max_expr, length.out = length(expr_colors))
    md_colors <- list(colors = expr_colors, breaks = color_breaks)
    palette <- circlize::colorRamp2(colors = expr_colors, breaks = color_breaks)
    df$color <- palette(df$expression)
    df$color_values <- df$expression

    template <- paste0("{color_name}: {color_values}", color_str_suffix)
    df <- df %>%
        mutate(color_str = glue(template,
            color_values = round(expression, digits = 3)
        ))

    return(list(df = df, md_colors = md_colors))
}


#' Resolve color for numeric metadata: colorRamp2 palette + color/color_values/color_str columns.
#' @noRd
resolve_numeric_md_color <- function(df, color_var, color_name, md_colors, color_str_suffix = "") {
    palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
    df$color <- palette(df[[color_var]])
    df$color_values <- df[[color_var]]
    template <- paste0("{color_name}: {color_values}", color_str_suffix)
    df <- df %>%
        mutate(color_str = glue(template,
            color_values = round(!!sym(color_var), digits = 3)
        ))
    return(df)
}


# ---------------------------------------------------------------------------
# mc2d metadata projection plots (plotly-based, not part of scatter unification)
# ---------------------------------------------------------------------------

#' Apply the standard plotly layout for mc2d projection plots.
#'
#' This sets hidden axes and zero margins, which is shared across categorical,
#' numeric, and gene-expression projection plots.
#'
#' @param fig A plotly figure
#' @param legend_title Optional legend title. When provided as a list
#'   (e.g. `list(text = "...")`) it is passed to \code{legend}; when a scalar
#'   string it is passed to \code{plotly::colorbar}.
#' @param use_colorbar If TRUE, add a colorbar with \code{legend_title} instead
#'   of a discrete legend title.
#' @return The plotly figure with layout applied.
#' @noRd
mc2d_plotly_proj_layout <- function(fig, legend_title = NULL, use_colorbar = FALSE) {
    layout_args <- list(
        xaxis = list(showgrid = FALSE, zeroline = FALSE, visible = FALSE),
        yaxis = list(showgrid = FALSE, zeroline = FALSE, visible = FALSE),
        margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0)
    )
    if (!is.null(legend_title) && !use_colorbar) {
        layout_args$legend <- list(title = list(text = legend_title))
    }
    fig <- do.call(plotly::layout, c(list(p = fig), layout_args))
    if (!is.null(legend_title) && use_colorbar) {
        fig <- fig %>% plotly::colorbar(title = legend_title)
    }
    return(fig)
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
                                   scale_edges = FALSE,
                                   metacell_types = NULL,
                                   atlas = FALSE,
                                   metadata = NULL,
                                   graph_name = NULL,
                                   mc2d = NULL,
                                   selected_cell_types = NULL) {
    mc2d <- mc2d %||% get_mc_data(dataset, "mc2d", atlas = atlas)
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
    } else {
        if (!is.null(colors)) {
            mc2d_df <- mc2d_df %>%
                mutate(`Cell type` = factor(`Cell type`, levels = sort(names(colors))))
        }
    }

    if (has_name(mc2d_df, "mc_age")) {
        mc2d_df <- mc2d_df %>% rename(`Age` = mc_age)
    }

    if (!is.null(graph_name) && graph_name != "metacell") {
        graph <- get_mc_data(dataset, "metacell_graphs")[[graph_name]]
    } else {
        graph <- NULL
    }

    graph <- mc2d_to_graph_df(mc2d, min_d = min_d, graph = graph)

    if (is.null(id)) {
        mc2d_df <- mc2d_df %>% mutate(id = metacell)
    } else {
        mc2d_df <- mc2d_df %>% mutate(id = paste(id, metacell, sep = "\t"))
    }

    if (!is.null(selected_cell_types)) {
        mc2d_df <- mc2d_df %>%
            filter(`Cell type` %in% selected_cell_types())
    }

    if (is_numeric_field(mc2d_df, md)) {
        p <- mc2d_plot_metadata_ggp_numeric(mc2d_df, graph, dataset, metadata, md, colors, color_breaks, point_size, min_d, stroke, graph_color, graph_width, id, scale_edges)
    } else {
        p <- mc2d_plot_metadata_ggp_categorical(mc2d_df, graph, dataset, metadata, md, point_size, min_d, stroke, graph_color, graph_width, id, scale_edges, colors)
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
            text = paste(
                glue("Metacell: {metacell}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                sep = "\n"
            )
        )

    if (md != "Cell type") {
        mc2d_df <- mc2d_df %>%
            mutate(text = paste0(text, "\n", md, ": ", mc2d_df[[md]]))
    }
    if (has_name(mc2d_df, "Age")) {
        mc2d_df <- mc2d_df %>%
            mutate(text = paste0(text, "\n", glue("Metacell age (E[t]): {round(Age, digits=2)}")))
    }

    if (is.null(colors)) {
        colors <- get_metadata_colors(dataset, md, metadata = metadata)
    }

    mc2d_df <- mc2d_df %>%
        arrange(desc(!!sym(md))) %>%
        mutate(value = !!sym(md))

    legend_title <- md

    add_scatter_layer <- function(x, showlegend = FALSE) {
        plotly::add_trace(x,
            data = mc2d_df,
            x = ~x,
            y = ~y,
            color = ~value,
            split = ~value,
            text = ~text,
            customdata = ~id,
            legendgroup = ~value,
            hoverinfo = "text",
            type = "scatter",
            mode = "markers",
            colors = colors,
            marker = list(
                size = point_size * 4,
                line = list(
                    color = "black",
                    width = stroke %||% 0.2
                )
            ),
            showlegend = showlegend
        )
    }
    # Draw edges first so scatter points appear on top
    fig <- plotly::plot_ly() %>%
        mc2d_add_graph_edges(graph, graph_color, graph_width) %>%
        add_scatter_layer(showlegend = TRUE)

    fig <- fig %>% mc2d_plotly_proj_layout(legend_title = legend_title)

    return(fig)
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
            text = paste(
                glue("Metacell: {metacell}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                paste0(md, ": ", round(mc2d_df[[md]], digits = 3)),
                ifelse(has_name(mc2d_df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    md_colors <- get_metadata_colors(dataset, md, colors = colors, color_breaks = color_breaks, metadata = metadata)
    palette <- circlize::colorRamp2(colors = md_colors$colors, breaks = md_colors$breaks)
    colors <- palette(seq(min(md_colors$breaks), max(md_colors$breaks), length.out = 100))

    mc2d_df <- mc2d_df %>%
        arrange(desc(!!sym(md))) %>%
        mutate(value = !!sym(md))

    legend_title <- md

    add_scatter_layer <- function(x, showlegend = FALSE) {
        plotly::add_trace(x,
            data = mc2d_df,
            x = ~x,
            y = ~y,
            color = ~value,
            text = ~text,
            customdata = ~id,
            hoverinfo = "text",
            type = "scatter",
            mode = "markers",
            colors = colors,
            marker = list(
                size = point_size * 4,
                line = list(
                    color = "black",
                    width = stroke %||% 0.2
                )
            ),
            showlegend = showlegend
        )
    }
    # Draw edges first so scatter points appear on top
    fig <- plotly::plot_ly() %>%
        mc2d_add_graph_edges(graph, graph_color, graph_width) %>%
        add_scatter_layer(showlegend = TRUE)

    fig <- fig %>% mc2d_plotly_proj_layout(legend_title = legend_title, use_colorbar = TRUE)

    return(fig)
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
                            expr_colors = mcview_palette$expression,
                            plot_text = TRUE,
                            atlas = FALSE,
                            metadata = get_mc_data(dataset, "metadata", atlas = atlas),
                            x_limits = NULL,
                            y_limits = NULL,
                            fixed_limits = FALSE,
                            xyline = FALSE,
                            metacell_filter = NULL,
                            show_correlation = TRUE,
                            correlation_type = "pearson",
                            corrected = FALSE,
                            log_labels = default_scatters_log_labels(dataset),
                            xylims = NULL) {
    if (!is.null(metadata)) {
        metadata <- metadata %>% mutate(metacell = as.character(metacell))
    }
    df <- metacell_types %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_na_value_to_level(cell_type, "(Missing)")) %>%
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
            egc_x <- get_gene_module_egc(x_var, dataset, gene_modules, atlas = atlas) + mcv_get("egc_epsilon")
        } else {
            egc_x <- get_gene_egc(x_var, dataset, atlas = atlas, corrected = corrected) + mcv_get("egc_epsilon")
        }
        df <- df %>%
            mutate(!!x_var := egc_x[metacell]) %>%
            mutate(x_str = glue("{x_name} expression: {expr_text} ({expr_text_log2})", expr_text = scales::scientific(!!sym(x_var)), expr_text_log2 = round(log2(!!sym(x_var)), digits = 2)))
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
            egc_y <- get_gene_module_egc(y_var, dataset, gene_modules, atlas = atlas) + mcv_get("egc_epsilon")
        } else {
            egc_y <- get_gene_egc(y_var, dataset, atlas = atlas, corrected = corrected) + mcv_get("egc_epsilon")
        }

        df <- df %>%
            mutate(!!y_var := egc_y[metacell]) %>%
            mutate(y_str = glue("{y_name} expression: {expr_text}, ({expr_text_log2})", expr_text = scales::scientific(!!sym(y_var)), expr_text_log2 = round(log2(!!sym(y_var)), digits = 2)))
    }

    # set color variable
    color_name <- color_var
    col_to_ct <- NULL
    categorical_md <- FALSE
    md_colors <- NULL
    if (is.null(color_var)) {
        if (atlas) {
            col_to_ct <- get_cell_type_colors(dataset, NULL, atlas = TRUE)
        } else {
            col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)
        }

        df <- df %>%
            mutate(color = cell_type, color_values = cell_type) %>%
            mutate(color_str = glue("Cell type: {`Cell type`}")) %>%
            mutate(color = factor(color, levels = sort(names(col_to_ct))))
    } else if (color_type == "Metadata") {
        req(metadata)
        df <- df %>%
            select(-any_of(color_var)) %>%
            left_join(metadata %>% select(metacell, !!color_var), by = "metacell")
        md_colors <- get_metadata_colors(dataset, color_var, colors = colors, color_breaks = color_breaks, metadata = metadata)
        if (is_numeric_field(metadata, color_var)) {
            df <- resolve_numeric_md_color(df, color_var, color_name, md_colors, color_str_suffix = "\nCell type: {`Cell type`}")
        } else {
            df <- df %>%
                mutate(color = !!sym(color_var), color_values = !!sym(color_var)) %>%
                mutate(color_str = glue("{color_name}: {color_values}"))
            categorical_md <- TRUE
        }
    } else if (color_type %in% c("Gene", "Gene module")) {
        if (color_type == "Gene module") {
            egc_color <- get_gene_module_egc(color_var, dataset, gene_modules, atlas = atlas) + mcv_get("egc_epsilon")
        } else {
            egc_color <- get_gene_egc(color_var, dataset, atlas = atlas, corrected = corrected) + mcv_get("egc_epsilon")
        }
        res <- resolve_gene_color(df, egc_color, "metacell", color_name, expr_colors, color_str_suffix = "\nCell type: {`Cell type`}\n")
        df <- res$df
        md_colors <- res$md_colors
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

    if (show_correlation) {
        x_vals <- df[[x_var]][!is.na(df[[x_var]]) & !is.na(df[[y_var]])]
        y_vals <- df[[y_var]][!is.na(df[[x_var]]) & !is.na(df[[y_var]])]
        if (length(unique(x_vals)) > 1 && length(unique(y_vals)) > 1) {
            correlation <- cor(df[[x_var]], df[[y_var]], method = correlation_type, use = "pairwise.complete.obs")
            correlation_text <- glue("Correlation: {round(correlation, 3)} ({correlation_type})")
        } else {
            correlation <- NA
            correlation_text <- "Correlation: N/A (zero variance)"
        }
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

    if (show_correlation) {
        p <- p + labs(title = correlation_text)
    }

    p <- apply_scatter_color_layer(p, color_var, categorical_md, point_size, stroke, col_to_ct = col_to_ct, md_colors = md_colors, fill_name = "")

    if (fixed_limits && x_type %in% c("Gene", "Gene module") && y_type %in% c("Gene", "Gene module")) {
        x_limits <- x_limits %||% c(min(egc_x), max(egc_x))
        y_limits <- y_limits %||% c(min(egc_y), max(egc_y))
        x_limits <- c(min(c(x_limits[1], y_limits[1])), max(c(x_limits[2], y_limits[2])))
        y_limits <- x_limits
    }

    if (is.null(xylims)) {
        if (log_labels) {
            xylims <- 2^seq(-17, 0, by = 1)
        } else {
            xylims <- mcv_get("expr_breaks")
        }
    }

    if (x_type %in% c("Gene", "Gene module")) {
        x_limits <- x_limits %||% c(min(egc_x), max(egc_x))
        p <- apply_gene_axis_scale(p, "x", x_var, egc_x, xylims, limits = x_limits, log_labels = log_labels, corrected = corrected)
    }

    if (y_type %in% c("Gene", "Gene module")) {
        y_limits <- y_limits %||% c(min(egc_y), max(egc_y))
        p <- apply_gene_axis_scale(p, "y", y_var, egc_y, xylims, limits = y_limits, log_labels = log_labels, corrected = corrected, rotate_x = FALSE)
    }

    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    p <- p + theme_mcview()

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
                                expr_colors = mcview_palette$expression,
                                plot_text = TRUE) {
    metadata <- get_mc_data(dataset, "cell_metadata")

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
        egc_x <- get_samples_gene_egc(x_var, dataset, selected_mc) + mcv_get("egc_epsilon")
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
        egc_y <- get_samples_gene_egc(y_var, dataset, selected_mc) + mcv_get("egc_epsilon")
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
    md_colors <- NULL
    if (!is.null(color_var) && color_var != "None") {
        if (color_type == "Metadata") {
            req(color_var %in% colnames(metadata))
            df <- df %>%
                select(-any_of(color_var)) %>%
                left_join(get_var_md(color_var, color_name, NULL), by = "samp_id")

            if (is_numeric_field(df, color_var)) {
                md_colors <- get_metadata_colors(dataset, color_var, colors = colors, color_breaks = color_breaks, metadata = metadata)
                df <- resolve_numeric_md_color(df, color_var, color_name, md_colors)
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
            egc_color <- get_samples_gene_egc(color_var, dataset, selected_mc) + mcv_get("egc_epsilon")
            res <- resolve_gene_color(df, egc_color, "samp_id", color_name, expr_colors, color_str_suffix = "\n", guard_equal_range = TRUE)
            df <- res$df
            md_colors <- res$md_colors
            # Override color_str to match original format with " (log2): " prefix
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
            df <- resolve_numeric_md_color(df, color_var, color_name, md_colors)
            # Override color_str for percent formatting
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

    # build ggplot + color layers
    base_aes <- aes(x = !!sym(x_var), y = !!sym(y_var), label = samp_id, customdata = samp_id, tooltip_text = Sample)
    if (!is.null(color_var) && color_var != "None" && color_var_type == "cont") {
        p <- ggplot(data = df, mapping = utils::modifyList(base_aes, aes(fill = color, color = color_values)))
        p <- apply_scatter_color_layer(p, color_var, FALSE, point_size, stroke, md_colors = md_colors)
    } else if (!is.null(color_var) && color_var != "None" && color_var_type == "discrete") {
        p <- ggplot(data = df, mapping = utils::modifyList(base_aes, aes(fill = !!sym(color_var)))) +
            geom_point(size = point_size) +
            scale_fill_manual(values = colors)
    } else {
        p <- ggplot(data = df, mapping = base_aes) +
            geom_point(size = point_size, color = "black")
    }

    p <- p +
        xlab(x_var) +
        ylab(y_var)

    # arrange axis for gene expression
    xylims <- mcv_get("expr_breaks")

    if (x_type %in% c("Gene", "Gene module")) {
        p <- apply_gene_axis_scale(p, "x", x_var, egc_x, xylims, log_labels = FALSE, corrected = FALSE, tolerance = 0)
    } else if (x_type == "Cell type") {
        p <- p +
            scale_x_continuous(labels = scales::percent)
    }

    if (y_type %in% c("Gene", "Gene module")) {
        p <- apply_gene_axis_scale(p, "y", y_var, egc_y, xylims, log_labels = FALSE, corrected = FALSE, tolerance = 0, rotate_x = FALSE)
    } else if (y_type == "Cell type") {
        p <- p +
            scale_y_continuous(labels = scales::percent)
    }

    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    p <- p + theme_mcview()

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
                                  expr_colors = mcview_palette$expression,
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
        mutate(cell_type = forcats::fct_na_value_to_level(cell_type, "(Missing)")) %>%
        mutate(`Cell type` = cell_type)

    # set axis variables
    correction_factor <- NULL
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
        egc_obs <- get_gene_egc(axis_var, dataset, corrected = TRUE) + mcv_get("egc_epsilon")
        egc_proj <- get_gene_egc(axis_var, dataset, projected = TRUE) + mcv_get("egc_epsilon")

        x_var <- glue("{axis_var} - observed (corrected)")
        y_var <- glue("{axis_var} - projected")
        df <- df %>%
            mutate(!!x_var := egc_obs[metacell], !!y_var := egc_proj[metacell]) %>%
            mutate(x_str = glue("{axis_name} obs: {expr_text}", expr_text = scales::scientific(!!sym(x_var)))) %>%
            mutate(y_str = glue("{axis_name} proj: {expr_text}", expr_text = scales::scientific(!!sym(y_var))))

        # get correction factor if exists
        gene_qc <- get_gene_qc(dataset)
        if (!is.null(gene_qc) && has_name(gene_qc, "correction_factor") && axis_var %in% gene_qc$gene) {
            correction_factor <- gene_qc$correction_factor[gene_qc$gene == axis_var]
        }
    }

    categorical_md <- FALSE
    col_to_ct <- NULL
    md_colors <- NULL
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
            df <- resolve_numeric_md_color(df, color_var, color_name, md_colors, color_str_suffix = "\nCell type: {`Cell type`}")
        }
    } else if (color_type == "Gene") {
        egc_color <- get_gene_egc(color_var, dataset) + mcv_get("egc_epsilon")
        res <- resolve_gene_color(df, egc_color, "metacell", color_name, expr_colors, color_str_suffix = "\nCell type: {`Cell type`}\n")
        df <- res$df
        md_colors <- res$md_colors
    }

    # set tooltip
    df <- df %>%
        mutate(
            Metacell = paste(
                glue("{metacell}\n{x_str}\n{y_str}\n{color_str}\nTop genes: {`Top genes`}"),
                ifelse(has_name(df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                ifelse(!is.null(correction_factor), glue("Correction factor: {round(correction_factor, 3)}"), ""),
                sep = "\n"
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

    if (!is.null(correction_factor)) {
        p <- p +
            geom_abline(intercept = -correction_factor, slope = 1, linetype = "dotted", color = "red")
    }

    # set color plotting
    if (is.null(color_var)) {
        col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)
    }
    p <- apply_scatter_color_layer(p, color_var, categorical_md, point_size, stroke, col_to_ct = col_to_ct, md_colors = md_colors)

    # arrange axis for gene expression
    xylims <- mcv_get("expr_breaks")

    if (axis_type %in% c("Gene", "Gene module")) {
        p <- apply_gene_axis_scale(p, "x", x_var, egc_obs, xylims)
        p <- apply_gene_axis_scale(p, "y", y_var, egc_proj, xylims, rotate_x = FALSE)
    }

    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    p <- p + theme_mcview()

    return(p)
}
