
#' Plot 2d projection of mc2d colored by cell types
#'
#' @param dataset name of metacell object
#' @param highlight data.frame with 'metacell',"label" and 'color'
#'
#' @noRd
mc2d_plot_ggp <- function(dataset, highlight = NULL, point_size = initial_proj_point_size(dataset), min_d = min_edge_length(dataset), stroke = NULL, graph_color = "black", graph_width = 0.1, scale_edges = FALSE, id = NULL, atlas = FALSE, metacell_types = get_mc_data(dataset, "metacell_types", atlas = atlas), cell_type_colors = get_mc_data(dataset, "cell_type_colors", atlas = atlas)) {
    mc2d <- get_mc_data(dataset, "mc2d", atlas = atlas)

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
        scale_fill_manual(name = "", values = get_cell_type_colors(dataset, cell_type_colors, atlas = atlas))


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
mc2d_plot_gene_ggp <- function(dataset, gene, point_size = initial_proj_point_size(dataset), min_d = min_edge_length(dataset), stroke = initial_proj_stroke(dataset), graph_color = "black", graph_width = 0.1, id = NULL, max_lfp = NULL, min_lfp = NULL, max_expr = NULL, min_expr = NULL, scale_edges = FALSE, stat = "expression", atlas = FALSE, gene_name = NULL) {
    mc2d <- get_mc_data(dataset, "mc2d", atlas = atlas)
    metacell_types <- get_mc_data(dataset, "metacell_types", atlas = atlas)
    min_lfp <- min_lfp %||% -3
    max_lfp <- max_lfp %||% 3
    if (length(gene) > 1) {
        lfp <- colSums(get_mc_egc(dataset, genes = gene, atlas = atlas), na.rm = TRUE) + egc_epsilon
        mc_fp <- lfp / median(lfp, na.rm = TRUE)
        gene <- gene_name %||% gene[1]
    } else if (length(gene) == 1) {
        mc_fp <- get_gene_fp(gene, dataset, atlas = atlas)
        lfp <- get_gene_egc(gene, dataset, atlas = atlas) + egc_epsilon
    } else {
        stop("gene paramater should have at least one gene")
    }

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

render_2d_plotly <- function(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, source, buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = NULL, refresh_on_gene_change = FALSE, atlas = FALSE, query_types = NULL, group = NULL, groupA = NULL, groupB = NULL, selected_metacell_types = NULL) {
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

        plot_2d_gene <- function(gene, gene_name = NULL) {
            req(input$proj_stat)
            if (input$proj_stat == "enrichment") {
                req(input$lfp)
            }
            if (!is.null(input$set_range) && input$set_range) {
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
                stat = input$proj_stat,
                atlas = atlas,
                gene_name = gene_name
            ) %>%
                plotly::ggplotly(tooltip = "tooltip_text", source = source) %>%
                rm_plotly_grid()

            # This ugly hack is due to https://github.com/ropensci/plotly/issues/1234
            # We need to remove the legend generated by scale_color_identity
            fig$x$data <- fig$x$data %>% purrr::map(~ {
                .x$showlegend <- FALSE
                .x
            })

            return(fig)
        }

        plot_2d_metadata <- function(md, metadata = NULL, colors = NULL, color_breaks = NULL) {
            fig <- mc2d_plot_metadata_ggp(
                dataset(),
                md,
                point_size = input$point_size,
                min_d = input$min_edge_size,
                metacell_types = metacell_types(),
                atlas = atlas,
                metadata = metadata,
                colors = colors,
                color_breaks = color_breaks
            ) %>%
                plotly::ggplotly(tooltip = "tooltip_text", source = source) %>%
                rm_plotly_grid()

            metadata <- metadata %||% get_mc_data(dataset(), "metadata", atlas = atlas)
            if (!is.null(metadata) && is_numeric_field(metadata, md)) {
                # This ugly hack is due to https://github.com/ropensci/plotly/issues/1234
                # We need to remove the legend generated by scale_color_identity
                fig$x$data <- fig$x$data %>% purrr::map(~ {
                    .x$showlegend <- FALSE
                    .x
                })
            }

            return(fig)
        }

        plot_2d_atlas_proj <- function(type) {
            req(input$color_by_scale)
            all_mc_w <- get_mc_data(dataset(), "proj_weights")
            req(all_mc_w)

            if (type == "Query cell type") {
                req(input$selected_cell_types)
                req(query_types)
                metacells <- query_types() %>%
                    filter(cell_type %in% input$selected_cell_types) %>%
                    pull(metacell)
                mc_proj_w <- all_mc_w %>%
                    filter(query %in% metacells)
            } else {
                req(input$selected_metacell)
                mc_proj_w <- all_mc_w %>%
                    filter(query == input$selected_metacell)
            }

            mc_proj_w <- mc_proj_w %>%
                select(metacell = atlas, Weight = weight) %>%
                mutate(query = "query")
            metadata <- tibble(metacell = colnames(get_mc_data(dataset(), "mc_mat", atlas = TRUE))) %>%
                left_join(mc_proj_w, by = "metacell") %>%
                tidyr::replace_na(replace = list(Weight = 0, query = "other"))

            req(input$query_threshold)
            if (input$color_by_scale == "Discrete") {
                metadata <- metadata %>%
                    mutate(query = ifelse(Weight <= input$query_threshold, "other", query))
                fig <- plot_2d_metadata("query", metadata = metadata, colors = c("query" = "darkred", "other" = "white"))
            } else if (input$color_by_scale == "Continuous") {
                fig <- plot_2d_metadata("Weight", metadata = metadata, colors = c("white", viridis::viridis_pal()(6)), color_breaks = c(0, seq(input$query_threshold, 1, length.out = 6)))
            } else {
                metadata <- metadata %>%
                    mutate(query = ifelse(Weight <= input$query_threshold, "other", query)) %>%
                    left_join(metacell_types() %>% select(metacell, cell_type), by = "metacell") %>%
                    mutate(query = ifelse(query != "other", cell_type, query))
                fig <- plot_2d_metadata("query", metadata = metadata, colors = c("other" = "white", get_cell_type_colors(dataset(), cell_type_colors())))
            }

            return(fig)
        }

        if (input$color_proj == "Cell type") {
            req(metacell_types())
            req(cell_type_colors())

            if (refresh_on_gene_change) {
                req(input$gene1)
                req(input$gene2)
            }

            fig <- mc2d_plot_ggp(
                dataset(),
                metacell_types = metacell_types(),
                cell_type_colors = cell_type_colors(),
                point_size = input$point_size,
                stroke = input$stroke,
                min_d = input$min_edge_size,
                highlight = highlight,
                atlas = atlas
            ) %>%
                plotly::ggplotly(tooltip = "tooltip_text", source = source)
            if (show_selected_metacells) {
                fig <- fig %>% plotly::hide_legend()
            }
        } else if (input$color_proj == "Gene A") {
            req(input$gene1)
            fig <- plot_2d_gene(input$gene1)
        } else if (input$color_proj == "Gene B") {
            req(input$gene2)
            fig <- plot_2d_gene(input$gene2)
        } else if (input$color_proj == "Metadata") {
            req(input$color_proj_metadata)
            metadata <- get_mc_data(dataset(), "metadata") %>%
                mutate(Clipboard = ifelse(metacell %in% globals$clipboard, "selected", "not selected"))
            fig <- plot_2d_metadata(input$color_proj_metadata, metadata = metadata)
        } else if (input$color_proj == "Gene") {
            req(input$color_proj_gene)
            fig <- plot_2d_gene(input$color_proj_gene)
        } else if (input$color_proj == "Gene module") {
            req(input$color_proj_gene_module)
            genes <- get_module_genes(input$color_proj_gene_module, gene_modules())
            fig <- plot_2d_gene(genes, gene_name = input$color_proj_gene_module)
        } else if (input$color_proj == "Sample") {
            req(input$samp1)
            fig <- plot_2d_metadata(paste0("samp_id: ", input$samp1))
        } else if (input$color_proj == "Similarity") {
            fig <- plot_2d_metadata("similar", colors = c("similar" = "white", "dissimilar" = "darkred"))
        } else if (input$color_proj %in% c("Query cell type", "Query metacell")) {
            fig <- plot_2d_atlas_proj(input$color_proj)
        } else if (input$color_proj == "Query Metadata") {
            req(input$color_proj_query_metadata)
            fig <- plot_2d_metadata(input$color_proj_query_metadata)
        } else if (input$color_proj == "Atlas Metadata") {
            req(input$color_proj_atlas_metadata)
            atlas_metadata <- get_mc_data(dataset(), "metadata", atlas = TRUE)
            proj_w <- get_mc_data(dataset(), "proj_weights")
            req(proj_w)
            metadata <- proj_w %>%
                mutate(atlas = as.character(atlas)) %>%
                left_join(
                    atlas_metadata %>%
                        select(atlas = metacell, !!input$color_proj_atlas_metadata) %>%
                        mutate(atlas = as.character(atlas)),
                    by = "atlas"
                ) %>%
                group_by(query) %>%
                summarise(!!input$color_proj_atlas_metadata := sum(weight * !!sym(input$color_proj_atlas_metadata))) %>%
                rename(metacell = query)
            fig <- plot_2d_metadata(input$color_proj_atlas_metadata, metadata = metadata)
        } else if (input$color_proj == "Selected") {
            if (!is.null(input$mode) && input$mode == "Groups") {
                req(groupA)
                req(groupB)
                selected_metacells1 <- groupA()
                selected_metacells2 <- groupB()
                metadata <- metacell_types() %>%
                    select(metacell) %>%
                    mutate(grp = case_when(
                        metacell %in% selected_metacells1 ~ "Group A",
                        metacell %in% selected_metacells2 ~ "Group B",
                        TRUE ~ "other"
                    ))
                fig <- plot_2d_metadata("grp", metadata = metadata, colors = c("Group A" = "red", "Group B" = "blue", "other" = "gray"))
            } else {
                if (!is.null(selected_metacell_types)) {
                    selected_metacells <- selected_metacell_types()$metacell
                } else {
                    req(input$mode)
                    if (input$mode == "MC") {
                        req(input$metacell1)
                        selected_metacells <- input$metacell1
                    } else if (input$mode == "Type") {
                        req(input$metacell1)
                        selected_metacells <- metacell_types() %>%
                            filter(cell_type %in% input$metacell1) %>%
                            pull(metacell)
                    } else if (input$mode == "Group") {
                        selected_metacells <- group()
                    } else {
                        req(FALSE)
                    }
                }
                metadata <- metacell_types() %>%
                    select(metacell) %>%
                    mutate(grp = ifelse(metacell %in% selected_metacells, "selected",
                        "other"
                    ))
                fig <- plot_2d_metadata("grp", metadata = metadata, colors = c("selected" = "red", "other" = "gray"))
            }
        }

        fig <- fig %>% plotly::event_register("plotly_restyle")

        fig <- fig %>% sanitize_plotly_buttons(buttons = buttons)

        fig$x$source <- source

        if (!is.null(dragmode)) {
            fig <- fig %>% plotly::layout(dragmode = dragmode)
        } else if (!is.null(input$mode) && input$mode %in% c("Groups", "Group")) {
            fig <- fig %>% plotly::layout(dragmode = "select")
        }

        fig <- fig %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            arrange_2d_proj_tooltip() %>%
            rm_plotly_grid()

        return(fig)
    })
}



initial_proj_point_size <- function(dataset, screen_width = NULL, screen_height = NULL, weight = 1, atlas = FALSE) {
    if (!is.null(config$datasets[[dataset]]$projection_point_size)) {
        return(config$datasets[[dataset]]$projection_point_size * weight)
    }
    n_metacells <- length(get_mc_data(dataset, "mc_sum", atlas = atlas))
    screen_width <- screen_width %||% 1920
    screen_height <- screen_height %||% 1080
    point_size <- screen_width * screen_height / (n_metacells * 1400) * weight

    return(max(1, min(3, point_size)))
}

initial_proj_stroke <- function(dataset) {
    if (!is.null(config$datasets[[dataset]]$projection_stroke)) {
        return(config$datasets[[dataset]]$projection_stroke)
    }
    return(0.1)
}

min_edge_length <- function(dataset) {
    config$datasets[[dataset]]$min_d %||% 0.025
}
