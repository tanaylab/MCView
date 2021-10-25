
#' Plot 2d projection of mc2d colored by cell types
#'
#' @param dataset name of metacell object
#' @param highlight data.frame with 'metacell',"label" and 'color'
#'
#' @noRd
mc2d_plot_ggp <- function(dataset, highlight = NULL, point_size = initial_proj_point_size(dataset), min_d = min_edge_length(dataset), stroke = NULL, graph_color = "black", graph_width = 0.1, scale_edges = FALSE, id = NULL, metacell_types = get_mc_data(dataset, "metacell_types"), cell_type_colors = get_mc_data(dataset, "cell_type_colors")) {
    mc2d <- get_mc_data(dataset, "mc2d")

    mc2d_df <- mc2d_to_df(mc2d) %>%
        left_join(metacell_types, by = "metacell") %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        rename(
            `Cell type` = cell_type,
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
                ifelse(has_name(mc2d_df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    p <- mc2d_df %>%
        ggplot(aes(x = x, y = y, label = metacell, fill = `Cell type`, tooltip_text = Metacell, customdata = id)) +
        scale_fill_manual(name = "", values = get_cell_type_colors(dataset, cell_type_colors))


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

    stroke <- stroke %||% min(0.2, 5e5 / nrow(mc2d_df)^2)

    p <- p + geom_point(size = point_size, shape = 21, stroke = stroke) +
        theme_void() +
        guides(size = "none")

    if (!is.null(highlight)) {
        cols <- highlight %>%
            select(label, color) %>%
            tibble::deframe()
        cols <- c(cols, "Other" = "black")

        p <- p +
            geom_point(
                inherit.aes = TRUE,
                data = mc2d_df %>% inner_join(highlight, by = "metacell"),
                shape = 21,
                size = point_size * 2.5,
                stroke = 1,
                fill = NA,
                aes(color = label)
            ) + scale_color_manual(values = cols)

        p <- p + guides(color = "none")
    }


    return(p)
}


#' Plot 2d projection of mc2d colored by gene
#'
#' @param dataset name of metacell object
#'
#' @noRd
mc2d_plot_gene_ggp <- function(dataset, gene, point_size = initial_proj_point_size(dataset), min_d = min_edge_length(dataset), stroke = initial_proj_stroke(dataset), graph_color = "black", graph_width = 0.1, id = NULL, max_lfp = NULL, min_lfp = NULL, max_expr = NULL, min_expr = NULL, scale_edges = FALSE, stat = "expression") {
    mc2d <- get_mc_data(dataset, "mc2d")
    metacell_types <- get_mc_data(dataset, "metacell_types")
    min_lfp <- min_lfp %||% -3
    max_lfp <- max_lfp %||% 3

    mc_fp <- get_gene_fp(gene, dataset)

    lfp <- get_gene_egc(gene, dataset) + egc_epsilon

    metacell_names <- mc2d_to_df(mc2d)$metacell

    mc2d_df <- mc2d_to_df(mc2d) %>%
        left_join(metacell_types, by = "metacell") %>%
        mutate(
            !!gene := lfp[metacell_names],
            expression = log2(lfp[metacell_names]),
            enrich := pmin(pmax(log2(mc_fp[metacell_names]), min_lfp), max_lfp) - min_lfp,
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        rename(
            `Cell type` = cell_type
        )

    min_expr <- min_expr %||% min(mc2d_df$expression, na.rm = TRUE)
    max_expr <- max_expr %||% max(mc2d_df$expression, na.rm = TRUE)

    if (has_name(df, "mc_age")) {
        mc2d_df <- mc2d_df %>% rename(`Age` = mc_age)
    }

    graph <- mc2d_to_graph_df(mc2d, min_d = min_d)

    if (is.null(id)) {
        mc2d_df <- mc2d_df %>% mutate(id = metacell)
    } else {
        mc2d_df <- mc2d_df %>% mutate(id = paste(id, metacell, sep = "\t"))
    }

    # define colors
    colspec <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

    mc2d_df <- mc2d_df %>%
        mutate(
            expr_clipped = tgutil::clip_vals(expression, min_val = min_expr, max_val = max_expr),
            expr_trans = abs(expr_clipped) - abs(max_expr),
            expr_text = scales::scientific(!!sym(gene)),
            Metacell = paste(
                glue("{metacell}"),
                glue("{gene} expression: {expr_text}"),
                glue("{gene} enrichment: {round(log2(mc_fp[metacell_names]), digits=2)}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                ifelse(has_name(df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    if (stat == "enrichment") {
        shades <- grDevices::colorRampPalette(colspec)(100 * (max_lfp - min_lfp) + 1)
        p <- mc2d_df %>%
            mutate(col_x = shades[round(100 * enrich) + 1]) %>%
            ggplot(aes(x = x, y = y, label = metacell, fill = col_x, color = enrich + min_lfp, tooltip_text = Metacell, customdata = id))
        legend_title <- glue("{gene}\nEnrichment.\n(log2)")
        shades_subset <- shades[seq(round(100 * min(mc2d_df$enrich)), round(100 * max(mc2d_df$enrich)), 1) + 1]
    } else { # expression
        shades <- rev(grDevices::colorRampPalette(colspec)(100 * abs(max_expr - min_expr) + 1))
        p <- mc2d_df %>%
            mutate(
                col_x = shades[round(100 * expr_trans) + 1]
            ) %>%
            ggplot(aes(x = x, y = y, label = metacell, fill = col_x, color = expr_clipped, tooltip_text = Metacell, customdata = id))
        legend_title <- glue("{gene}\nExpression.\n(log2)")

        shades_subset <- rev(shades[abs(seq(round(100 * min(mc2d_df$expr_trans)), round(100 * max(mc2d_df$expr_trans)), 1) + 1)])
    }

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
        scale_color_gradientn(name = legend_title, colors = shades_subset) +
        scale_fill_identity()

    return(p)
}

render_2d_plotly <- function(input, output, session, dataset, values, metacell_types, cell_type_colors, source, buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = NULL, refresh_on_gene_change = FALSE) {
    plotly::renderPlotly({
        req(input$color_proj)
        req(input$point_size)
        req(input$stroke)
        req(input$min_edge_size)

        show_selected_metacells <- !is.null(input$show_selected_metacells) && input$show_selected_metacells

        if (show_selected_metacells && !is.null(input$metacell1) && !is.null(input$metacell2) && input$color_proj == "Cell type") {
            highlight <- tibble::tibble(
                metacell = c(input$metacell1, input$metacell2),
                label = c("metacell1", "metacell2"),
                color = c("darkred", "darkblue")
            )
        } else {
            highlight <- NULL
        }

        plot_2d_gene <- function(gene) {
            req(input$proj_stat)
            if (input$proj_stat == "enrichment") {
                req(input$lfp)
            }
            if (input$set_range) {
                min_expr <- input$expr_range[1]
                max_expr <- input$expr_range[2]
            } else {
                min_expr <- NULL
                max_expr <- NULL
            }

            fig <- mc2d_plot_gene_ggp(
                dataset(),
                gene,
                min_lfp = input$lfp[1],
                max_lfp = input$lfp[2],
                min_expr = min_expr,
                max_expr = max_expr,
                point_size = input$point_size,
                min_d = input$min_edge_size,
                stat = input$proj_stat
            ) %>%
                plotly::ggplotly(tooltip = "tooltip_text", source = source)

            # This ugly hack is due to https://github.com/ropensci/plotly/issues/1234
            # We need to remove the legend generated by scale_color_identity
            fig$x$data <- fig$x$data %>% purrr::map(~ {
                .x$showlegend <- FALSE
                .x
            })

            return(fig)
        }

        plot_2d_metadata <- function(md) {
            fig <- mc2d_plot_metadata_ggp(
                dataset(),
                md,
                point_size = input$point_size,
                min_d = input$min_edge_size
            ) %>%
                plotly::ggplotly(tooltip = "tooltip_text", source = source)

            # This ugly hack is due to https://github.com/ropensci/plotly/issues/1234
            # We need to remove the legend generated by scale_color_identity
            fig$x$data <- fig$x$data %>% purrr::map(~ {
                .x$showlegend <- FALSE
                .x
            })

            return(fig)
        }

        if (input$color_proj == "Cell type") {
            if (refresh_on_gene_change) {
                req(values$gene1)
                req(values$gene2)
            }
            fig <- mc2d_plot_ggp(
                dataset(),
                metacell_types = metacell_types(),
                cell_type_colors = cell_type_colors(),
                point_size = input$point_size,
                stroke = input$stroke,
                min_d = input$min_edge_size,
                highlight = highlight
            ) %>%
                plotly::ggplotly(tooltip = "tooltip_text", source = source)
            if (show_selected_metacells) {
                fig <- fig %>% plotly::hide_legend()
            }
        } else if (input$color_proj == "Gene A") {
            req(values$gene1)
            fig <- plot_2d_gene(values$gene1)
        } else if (input$color_proj == "Gene B") {
            req(values$gene2)
            fig <- plot_2d_gene(values$gene2)
        } else {
            fig <- plot_2d_metadata(input$color_proj)
        }

        fig <- fig %>% plotly::event_register("plotly_restyle")

        fig <- fig %>% sanitize_plotly_buttons(buttons = buttons)

        fig$x$source <- source

        if (!is.null(dragmode)) {
            fig <- fig %>% plotly::layout(dragmode = dragmode)
        }

        fig <- fig %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            arrange_2d_proj_tooltip()

        return(fig)
    })
}

# TODO: find a better heuristic that takes into account the plot size
initial_proj_point_size <- function(dataset) {
    if (!is.null(config$datasets[[dataset]]$projection_point_size)) {
        return(config$datasets[[dataset]]$projection_point_size)
    }
    n_metacells <- length(get_mc_data(dataset, "mc_sum"))
    return(max(1, min(1.5, 3e4 / n_metacells)))
}

initial_proj_stroke <- function(dataset) {
    if (!is.null(config$datasets[[dataset]]$projection_stroke)) {
        return(config$datasets[[dataset]]$projection_stroke)
    }
    n_metacells <- length(get_mc_data(dataset, "mc_sum"))
    return(min(0.1, 5e5 / n_metacells^2))
}

min_edge_length <- function(dataset) {
    config$datasets[[dataset]]$min_d %||% 0.05
}
