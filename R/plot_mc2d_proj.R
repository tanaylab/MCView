#' Plot 2d projection of mc2d colored by gene
#'
#' @param dataset name of metacell object
#'
#' @noRd
mc2d_plot_gene_ggp <- function(dataset, gene, point_size = initial_proj_point_size(dataset), min_d = min_edge_length(dataset), stroke = initial_proj_stroke(dataset), graph_color = "black", graph_width = 0.1, id = NULL, max_lfp = NULL, min_lfp = NULL, max_expr = NULL, min_expr = NULL, stat = "expression", atlas = FALSE, gene_name = NULL, graph_name = NULL, mc2d = NULL, selected_cell_types = NULL, metacell_types = NULL) {
    mc2d <- mc2d %||% get_mc_data(dataset, "mc2d", atlas = atlas)
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

    # define colors
    colspec <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

    mc2d_df <- mc2d_df %>%
        mutate(
            expr_clipped = tgutil::clip_vals(expression, min_val = min_expr, max_val = max_expr),
            expr_text = round(expression, digits = 2),
            text = paste(
                glue("Metacell: {metacell}"),
                glue("{gene} expression (log2): {expr_text}"),
                glue("{gene} enrichment: {round(log2(mc_fp[metacell_names]), digits=2)}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                ifelse(has_name(df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    if (stat == "enrichment") {
        shades <- grDevices::colorRampPalette(colspec)(100 * (max_lfp - min_lfp) + 1)
        df <- mc2d_df %>%
            mutate(value = enrich)

        legend_title <- glue("{gene}\nEnrichment.\n(log2)")
    } else { # expression
        shades <- grDevices::colorRampPalette(colspec)(100 * abs(max_expr - min_expr) + 1)
        df <- mc2d_df %>%
            mutate(value = expr_clipped)

        legend_title <- glue("{gene}\nExpression.\n(log2)")
    }

    if (!is.null(selected_cell_types)) {
        df <- df %>%
            filter(`Cell type` %in% selected_cell_types())
    }

    add_scatter_layer <- function(x, showlegend = FALSE) {
        plotly::add_trace(x,
            data = df,
            x = ~x,
            y = ~y,
            color = ~value,
            text = ~text,
            customdata = ~id,
            hoverinfo = "text",
            type = "scatter",
            mode = "markers",
            colors = shades,
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

    fig <- plotly::plot_ly() %>% add_scatter_layer()

    if (nrow(graph) > 0) {
        edges_x <- c(rbind(graph$x_mc1, graph$x_mc2, NA))
        edges_y <- c(rbind(graph$y_mc1, graph$y_mc2, NA))

        fig <- fig %>%
            plotly::add_trace(
                x = edges_x,
                y = edges_y,
                type = "scatter",
                mode = "lines",
                line = list(
                    color = graph_color,
                    width = graph_width * 5
                ),
                showlegend = FALSE
            )
    }

    fig <- fig %>% add_scatter_layer(showlegend = TRUE)


    fig <- fig %>%
        plotly::layout(
            xaxis = list(
                showgrid = FALSE,
                zeroline = FALSE,
                visible = FALSE
            ),
            yaxis = list(
                showgrid = FALSE,
                zeroline = FALSE,
                visible = FALSE
            ),
            margin = list(
                l = 0,
                r = 0,
                b = 0,
                t = 0,
                pad = 0
            )
        ) %>%
        plotly::colorbar(title = legend_title)

    return(fig)
}

render_2d_plotly <- function(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, source, buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = NULL, refresh_on_gene_change = FALSE, atlas = FALSE, query_types = NULL, group = NULL, groupA = NULL, groupB = NULL, selected_metacell_types = NULL, selected_cell_types = NULL) {
    plotly::renderPlotly({
        req(input$color_proj)
        req(input$point_size)
        req(input$stroke)
        req(input$min_edge_size)

        if (atlas) {
            mc2d <- NULL
        } else {
            mc2d <- globals$mc2d %||% get_mc_data(dataset, "mc2d")
        }

        plot_2d_gene <- function(gene, gene_name = NULL) {
            req(proj_stat)
            if (proj_stat == "enrichment") {
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
                stat = proj_stat,
                atlas = atlas,
                stroke = input$stroke,
                gene_name = gene_name,
                graph_name = input$graph_name,
                mc2d = mc2d,
                metacell_types = metacell_types(),
                selected_cell_types = selected_cell_types
            )

            fig <- fig %>% rm_plotly_grid()

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
                stroke = input$stroke,
                color_breaks = color_breaks,
                graph_name = input$graph_name,
                mc2d = mc2d,
                selected_cell_types = selected_cell_types
            )

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
                fig <- plot_2d_metadata("query", stroke = input$stroke, metadata = metadata, colors = c("query" = "darkred", "other" = "white"))
            } else if (input$color_by_scale == "Continuous") {
                fig <- plot_2d_metadata("Weight", metadata = metadata, stroke = input$stroke, colors = c("white", viridis::viridis_pal()(6)), color_breaks = c(0, seq(input$query_threshold, 1, length.out = 6)))
            } else {
                metadata <- metadata %>%
                    mutate(query = ifelse(Weight <= input$query_threshold, "other", query)) %>%
                    left_join(metacell_types() %>% select(metacell, cell_type), by = "metacell") %>%
                    mutate(query = ifelse(query != "other", cell_type, query))
                fig <- plot_2d_metadata("query", metadata = metadata, stroke = input$stroke, colors = c("other" = "white", get_cell_type_colors(dataset(), cell_type_colors())))
            }

            return(fig)
        }

        color_proj <- input$color_proj
        color_proj_gene <- input$color_proj_gene
        color_proj_metadata <- input$color_proj_metadata
        color_proj_gene_module <- input$color_proj_gene_module
        proj_stat <- input$proj_stat

        if (color_proj == "Scatter Axis") {
            req(input$scatter_axis_proj)
            color_proj <- input[[glue("{input$scatter_axis_proj}_axis_type")]]
            color_proj_gene <- input[[glue("{input$scatter_axis_proj}_axis_var")]]
            color_proj_metadata <- input[[glue("{input$scatter_axis_proj}_axis_var")]]
            color_proj_gene_module <- input[[glue("{input$scatter_axis_proj}_axis_var")]]
            proj_stat <- "expression"
            req(axis_vars_ok(dataset(), input, "metadata", gene_modules, atlas = atlas))
        }


        if (color_proj == "Cell type") {
            req(metacell_types())
            req(cell_type_colors())
            fig <- mc2d_plot_metadata_ggp(
                dataset(),
                "Cell type",
                point_size = input$point_size,
                min_d = input$min_edge_size,
                metacell_types = metacell_types(),
                atlas = atlas,
                metadata = metacell_types() %>% rename(`Cell type` = cell_type),
                colors = get_cell_type_colors(dataset, cell_type_colors = cell_type_colors(), atlas = atlas),
                stroke = input$stroke,
                graph_name = input$graph_name,
                mc2d = mc2d,
                selected_cell_types = selected_cell_types
            )
        } else if (color_proj == "Gene A") {
            req(input$gene1)
            fig <- plot_2d_gene(input$gene1)
        } else if (color_proj == "Gene B") {
            req(input$gene2)
            fig <- plot_2d_gene(input$gene2)
        } else if (color_proj == "Metadata") {
            req(color_proj_metadata)
            metadata <- get_mc_data(dataset(), "metadata")
            if (is.null(metadata)) {
                metadata <- metacell_types() %>% select(metacell)
            }
            metadata <- metadata %>%
                mutate(Clipboard = ifelse(metacell %in% globals$clipboard, "selected", "not selected"))
            fig <- plot_2d_metadata(color_proj_metadata, metadata = metadata)
        } else if (color_proj == "Gene") {
            req(color_proj_gene)
            fig <- plot_2d_gene(color_proj_gene)
        } else if (color_proj == "Gene module") {
            req(color_proj_gene_module)
            genes <- get_module_genes(color_proj_gene_module, gene_modules())
            fig <- plot_2d_gene(genes, gene_name = color_proj_gene_module)
        } else if (color_proj == "Sample") {
            req(input$samp1)
            fig <- plot_2d_metadata(paste0("samp_id: ", input$samp1))
        } else if (color_proj == "Similarity") {
            fig <- plot_2d_metadata("similar", colors = c("similar" = "white", "dissimilar" = "darkred"))
        } else if (input$color_proj %in% c("Query cell type", "Query metacell")) {
            fig <- plot_2d_atlas_proj(input$color_proj)
        } else if (color_proj == "Query Metadata") {
            req(input$color_proj_query_metadata)
            fig <- plot_2d_metadata(input$color_proj_query_metadata)
        } else if (color_proj == "Atlas Metadata") {
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
        } else if (color_proj == "Selected") {
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

        fig$x$source <- source

        if (!is.null(dragmode)) {
            fig <- fig %>% plotly::layout(dragmode = dragmode)
        } else if (!is.null(input$mode) && input$mode %in% c("Groups", "Group")) {
            fig <- fig %>% plotly::layout(dragmode = "select")
            buttons <- buttons[!(buttons %in% c("select2d", "lasso2d"))]
        }

        fig <- fig %>%
            sanitize_plotly_buttons(buttons = buttons) %>%
            sanitize_plotly_download(globals)

        if (!is.null(input$legend_orientation)) {
            if (input$legend_orientation == "Horizontal") {
                orientation <- "h"
            } else if (input$legend_orientation == "Vertical") {
                orientation <- "v"
            }
            fig <- fig %>% plotly::layout(legend = list(orientation = orientation))
        }

        if (!is.null(input$show_legend_projection) && !input$show_legend_projection) {
            fig <- plotly::hide_legend(fig)
        }

        fig <- fig %>%
            sanitize_for_WebGL()
        fig <- fig %>%
            plotly::toWebGL()
        # fig <- fig %>%
        #     arrange_2d_proj_tooltip()
        fig <- fig %>%
            rm_plotly_grid()

        fig <- fig %>%
            sanitize_plotly_download(globals)

        return(fig)
    })
}



initial_proj_point_size <- function(dataset, screen_width = NULL, screen_height = NULL, weight = 1, atlas = FALSE) {
    if (!is.null(config$datasets[[dataset]]$projection_point_size)) {
        return(config$datasets[[dataset]]$projection_point_size * weight)
    } else if (!is.null(config$projection_point_size)) {
        return(config$projection_point_size * weight)
    }
    n_metacells <- length(get_mc_data(dataset, "mc_sum", atlas = atlas))
    screen_width <- screen_width %||% 1920
    screen_height <- screen_height %||% 1080
    desired_area <- (screen_width * screen_height) / (n_metacells * 100) * weight
    size_pixels <- sqrt(desired_area / pi)

    return(max(1, min(3, size_pixels)))
}

initial_proj_stroke <- function(dataset) {
    if (!is.null(config$datasets[[dataset]]$projection_stroke)) {
        return(config$datasets[[dataset]]$projection_stroke)
    }
    return(0.2)
}

min_edge_length <- function(dataset) {
    config$datasets[[dataset]]$min_d %||% 0.025
}
