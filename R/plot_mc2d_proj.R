#' Add graph edges to a plotly mc2d projection figure
#'
#' @param fig A plotly figure object
#' @param graph A data frame with graph edge coordinates (x_mc1, y_mc1, x_mc2, y_mc2)
#' @param graph_color Color for graph edges
#' @param graph_width Width multiplier for graph edges
#' @return The plotly figure with edges added
#' @noRd
mc2d_add_graph_edges <- function(fig, graph, graph_color = "black", graph_width = 0.1) {
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
    return(fig)
}

#' Plot 2d projection of mc2d colored by gene
#'
#' @param dataset name of metacell object
#'
#' @noRd
mc2d_plot_gene_ggp <- function(dataset,
                               gene,
                               point_size = initial_proj_point_size(dataset),
                               min_d = min_edge_length(dataset),
                               stroke = initial_proj_stroke(dataset),
                               graph_color = "black",
                               graph_width = 0.1,
                               id = NULL,
                               max_lfp = NULL,
                               min_lfp = NULL,
                               max_expr = NULL,
                               min_expr = NULL,
                               stat = "expression",
                               atlas = FALSE,
                               gene_name = NULL,
                               graph_name = NULL,
                               mc2d = NULL,
                               selected_cell_types = NULL,
                               metacell_types = NULL) {
    mc2d <- mc2d %||% get_mc_data(dataset, "mc2d", atlas = atlas)
    metacell_types <- get_mc_data(dataset, "metacell_types", atlas = atlas)
    min_lfp <- min_lfp %||% -3
    max_lfp <- max_lfp %||% 3
    if (length(gene) > 1) {
        lfp <- colSums(get_mc_egc(dataset, genes = gene, atlas = atlas), na.rm = TRUE) + mcv_get("egc_epsilon")
        mc_fp <- lfp / median(lfp, na.rm = TRUE)
        gene <- gene_name %||% gene[1]
    } else if (length(gene) == 1) {
        mc_fp <- get_gene_fp(gene, dataset, atlas = atlas)
        lfp <- get_gene_egc(gene, dataset, atlas = atlas) + mcv_get("egc_epsilon")
    } else {
        stop("gene paramater should have at least one gene")
    }

    mc2d_df <- mc2d_to_df(mc2d)
    metacell_names <- mc2d_df$metacell

    mc2d_df <- mc2d_df %>%
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
                ifelse(has_name(mc2d_df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
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

    # Draw edges first so scatter points appear on top
    fig <- plotly::plot_ly() %>%
        mc2d_add_graph_edges(graph, graph_color, graph_width) %>%
        add_scatter_layer(showlegend = TRUE)

    fig <- fig %>% mc2d_plotly_proj_layout(legend_title = legend_title, use_colorbar = TRUE)

    return(fig)
}

# ---------------------------------------------------------------------------
# Module-scope handler functions for render_2d_plotly dispatch
# ---------------------------------------------------------------------------

#' Create a gene projection figure for render_2d_plotly.
#'
#' Wraps mc2d_plot_gene_ggp with the standard Shiny input bindings.
#'
#' @noRd
render_2d_gene_fig <- function(gene, input, dataset, atlas, mc2d,
                               metacell_types, selected_cell_types,
                               proj_stat, gene_name = NULL) {
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

    return(fig)
}

#' Create a metadata projection figure for render_2d_plotly.
#'
#' Wraps mc2d_plot_metadata_ggp with the standard Shiny input bindings.
#'
#' @noRd
render_2d_metadata_fig <- function(md, input, dataset, atlas, mc2d,
                                   metacell_types, selected_cell_types,
                                   metadata = NULL, colors = NULL,
                                   color_breaks = NULL) {
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

# -- Individual color-mode handlers -------------------------------------------

#' @noRd
handle_2d_cell_type <- function(input, dataset, atlas, mc2d, metacell_types,
                                cell_type_colors, selected_cell_types) {
    req(metacell_types())
    req(cell_type_colors())
    mc2d_plot_metadata_ggp(
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
}

#' @noRd
handle_2d_gene_a <- function(input, dataset, atlas, mc2d, metacell_types,
                             selected_cell_types, proj_stat) {
    req(input$gene1)
    render_2d_gene_fig(input$gene1, input, dataset, atlas, mc2d,
        metacell_types, selected_cell_types, proj_stat)
}

#' @noRd
handle_2d_gene_b <- function(input, dataset, atlas, mc2d, metacell_types,
                             selected_cell_types, proj_stat) {
    req(input$gene2)
    render_2d_gene_fig(input$gene2, input, dataset, atlas, mc2d,
        metacell_types, selected_cell_types, proj_stat)
}

#' @noRd
handle_2d_metadata <- function(input, dataset, atlas, mc2d, metacell_types,
                               selected_cell_types, globals,
                               color_proj_metadata) {
    req(color_proj_metadata)
    metadata <- get_mc_data(dataset(), "metadata")
    if (is.null(metadata)) {
        metadata <- metacell_types() %>% select(metacell)
    }
    metadata <- metadata %>%
        mutate(Clipboard = ifelse(metacell %in% globals$clipboard, "selected", "not selected"))
    render_2d_metadata_fig(color_proj_metadata, input, dataset, atlas, mc2d,
        metacell_types, selected_cell_types, metadata = metadata)
}

#' @noRd
handle_2d_gene <- function(input, dataset, atlas, mc2d, metacell_types,
                           selected_cell_types, proj_stat, color_proj_gene) {
    req(color_proj_gene)
    render_2d_gene_fig(color_proj_gene, input, dataset, atlas, mc2d,
        metacell_types, selected_cell_types, proj_stat)
}

#' @noRd
handle_2d_gene_module <- function(input, dataset, atlas, mc2d, metacell_types,
                                  selected_cell_types, proj_stat,
                                  gene_modules, color_proj_gene_module) {
    req(color_proj_gene_module)
    genes <- get_module_genes(color_proj_gene_module, gene_modules())
    render_2d_gene_fig(genes, input, dataset, atlas, mc2d,
        metacell_types, selected_cell_types, proj_stat,
        gene_name = color_proj_gene_module)
}

#' @noRd
handle_2d_sample <- function(input, dataset, atlas, mc2d, metacell_types,
                             selected_cell_types) {
    req(input$samp1)
    render_2d_metadata_fig(paste0("samp_id: ", input$samp1), input, dataset,
        atlas, mc2d, metacell_types, selected_cell_types)
}

#' @noRd
handle_2d_similarity <- function(input, dataset, atlas, mc2d, metacell_types,
                                 selected_cell_types) {
    render_2d_metadata_fig("similar", input, dataset, atlas, mc2d,
        metacell_types, selected_cell_types,
        colors = c("similar" = "white", "dissimilar" = "darkred"))
}

#' @noRd
handle_2d_query_proj <- function(input, dataset, atlas, mc2d, metacell_types,
                                 selected_cell_types, cell_type_colors,
                                 query_types) {
    type <- input$color_proj
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
    atlas_daf <- get_atlas_daf()
    metadata <- tibble(metacell = dafr::axis_entries(atlas_daf, "metacell")) %>%
        left_join(mc_proj_w, by = "metacell") %>%
        tidyr::replace_na(replace = list(Weight = 0, query = "other"))

    req(input$query_threshold)
    if (input$color_by_scale == "Discrete") {
        metadata <- metadata %>%
            mutate(query = ifelse(Weight <= input$query_threshold, "other", query))
        fig <- render_2d_metadata_fig("query", input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types,
            metadata = metadata, colors = c("query" = "darkred", "other" = "white"))
    } else if (input$color_by_scale == "Continuous") {
        fig <- render_2d_metadata_fig("Weight", input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types,
            metadata = metadata,
            colors = c("white", viridis::viridis_pal()(6)),
            color_breaks = c(0, seq(input$query_threshold, 1, length.out = 6)))
    } else {
        metadata <- metadata %>%
            mutate(query = ifelse(Weight <= input$query_threshold, "other", query)) %>%
            left_join(metacell_types() %>% select(metacell, cell_type), by = "metacell") %>%
            mutate(query = ifelse(query != "other", cell_type, query))
        fig <- render_2d_metadata_fig("query", input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types,
            metadata = metadata,
            colors = c("other" = "white", get_cell_type_colors(dataset(), cell_type_colors())))
    }

    return(fig)
}

#' @noRd
handle_2d_query_metadata <- function(input, dataset, atlas, mc2d,
                                     metacell_types, selected_cell_types) {
    req(input$color_proj_query_metadata)
    render_2d_metadata_fig(input$color_proj_query_metadata, input, dataset,
        atlas, mc2d, metacell_types, selected_cell_types)
}

#' @noRd
handle_2d_atlas_metadata <- function(input, dataset, atlas, mc2d,
                                     metacell_types, selected_cell_types) {
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
    render_2d_metadata_fig(input$color_proj_atlas_metadata, input, dataset,
        atlas, mc2d, metacell_types, selected_cell_types,
        metadata = metadata)
}

#' @noRd
handle_2d_selected <- function(input, dataset, atlas, mc2d, metacell_types,
                               selected_cell_types, groupA, groupB, group,
                               selected_metacell_types) {
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
        return(render_2d_metadata_fig("grp", input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types,
            metadata = metadata,
            colors = c("Group A" = "red", "Group B" = "blue", "other" = "gray")))
    }

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
        mutate(grp = ifelse(metacell %in% selected_metacells, "selected", "other"))
    render_2d_metadata_fig("grp", input, dataset, atlas, mc2d,
        metacell_types, selected_cell_types,
        metadata = metadata,
        colors = c("selected" = "red", "other" = "gray"))
}

#' Dispatch to the appropriate handler for a given color_proj value.
#'
#' @return A plotly figure.
#' @noRd
dispatch_2d_color_proj <- function(color_proj, input, dataset, atlas, mc2d,
                                   metacell_types, cell_type_colors,
                                   gene_modules, globals, selected_cell_types,
                                   proj_stat, color_proj_gene,
                                   color_proj_metadata, color_proj_gene_module,
                                   query_types, group, groupA, groupB,
                                   selected_metacell_types) {
    switch(color_proj,
        "Cell type" = handle_2d_cell_type(input, dataset, atlas, mc2d,
            metacell_types, cell_type_colors, selected_cell_types),
        "Gene A" = handle_2d_gene_a(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types, proj_stat),
        "Gene B" = handle_2d_gene_b(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types, proj_stat),
        "Metadata" = handle_2d_metadata(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types, globals, color_proj_metadata),
        "Gene" = handle_2d_gene(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types, proj_stat, color_proj_gene),
        "Gene module" = handle_2d_gene_module(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types, proj_stat, gene_modules,
            color_proj_gene_module),
        "Sample" = handle_2d_sample(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types),
        "Similarity" = handle_2d_similarity(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types),
        "Query cell type" = ,
        "Query metacell" = handle_2d_query_proj(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types, cell_type_colors,
            query_types),
        "Query Metadata" = handle_2d_query_metadata(input, dataset, atlas,
            mc2d, metacell_types, selected_cell_types),
        "Atlas Metadata" = handle_2d_atlas_metadata(input, dataset, atlas,
            mc2d, metacell_types, selected_cell_types),
        "Selected" = handle_2d_selected(input, dataset, atlas, mc2d,
            metacell_types, selected_cell_types, groupA, groupB, group,
            selected_metacell_types)
    )
}

#' Apply common post-processing to a 2d projection plotly figure.
#'
#' Handles event registration, source assignment, drag mode, toolbar buttons,
#' legend orientation, WebGL conversion, and grid cleanup.
#'
#' @noRd
finalize_2d_plotly <- function(fig, input, globals, source, buttons, dragmode) {
    fig <- fig %>%
        plotly::event_register("plotly_restyle") %>%
        plotly::event_register("plotly_click") %>%
        plotly::event_register("plotly_selected") %>%
        plotly::event_register("plotly_deselect")
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
        sanitize_for_WebGL() %>%
        plotly::toWebGL()

    return(fig)
}

# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

render_2d_plotly <- function(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, source, buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = NULL, refresh_on_gene_change = FALSE, atlas = FALSE, query_types = NULL, group = NULL, groupA = NULL, groupB = NULL, selected_metacell_types = NULL, selected_cell_types = NULL, tab_guard = NULL) {
    plotly::renderPlotly({
        # Defer computation until the tab is active (when tab_guard is specified)
        if (!is.null(tab_guard)) {
            req(globals$current_tab == tab_guard)
        }
        req(input$color_proj)
        req(input$point_size)
        req(input$stroke)
        req(input$min_edge_size)

        if (atlas) {
            mc2d <- NULL
        } else {
            mc2d <- globals$mc2d %||% get_mc_data(dataset, "mc2d")
        }

        # Resolve color mode and variable names (handle Scatter Axis remapping)
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

        # Dispatch to the appropriate handler
        fig <- dispatch_2d_color_proj(
            color_proj, input, dataset, atlas, mc2d,
            metacell_types, cell_type_colors, gene_modules, globals,
            selected_cell_types, proj_stat, color_proj_gene,
            color_proj_metadata, color_proj_gene_module,
            query_types, group, groupA, groupB, selected_metacell_types
        )

        # Apply common post-processing
        fig <- finalize_2d_plotly(fig, input, globals, source, buttons, dragmode)

        return(fig)
    })
}


initial_proj_point_size <- function(dataset, screen_width = NULL, screen_height = NULL, weight = 1, atlas = FALSE) {
    if (!is.null(app_config("datasets")[[dataset]]$projection_point_size)) {
        return(app_config("datasets")[[dataset]]$projection_point_size * weight)
    } else if (!is.null(app_config("projection_point_size"))) {
        return(app_config("projection_point_size") * weight)
    }
    # Use DAF axis length instead of loading full mc_sum vector
    daf_obj <- get_daf_for_query(dataset, atlas)
    n_metacells <- if (!is.null(daf_obj)) dafr::axis_length(daf_obj, "metacell") else 100
    screen_width <- screen_width %||% 1920
    screen_height <- screen_height %||% 1080
    desired_area <- (screen_width * screen_height) / (n_metacells * 100) * weight
    size_pixels <- sqrt(desired_area / pi)

    return(max(1, min(3, size_pixels)))
}

initial_proj_stroke <- function(dataset) {
    if (!is.null(app_config("datasets")[[dataset]]$projection_stroke)) {
        return(app_config("datasets")[[dataset]]$projection_stroke)
    }
    return(0.2)
}

min_edge_length <- function(dataset) {
    app_config("datasets")[[dataset]]$min_d %||% 0.025
}
