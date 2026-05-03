# utils_heatmap_server.R - Main heatmap module body (reactive setup + render)
#
# Split from R/utils_heatmap.R (2026-05-01). Hosts heatmap_reactives(),
# the moduleServer that wires inputs, computes and renders the heatmap
# plot, and registers the major outputs.  Companions:
#   - R/utils_heatmap_data.R     - heatmap_matrix_reactives() (data layer)
#   - R/utils_heatmap_handlers.R - tooltip + download handlers
#   - R/utils_heatmap_ui.R       - UI builders (existing, unchanged)
#   - R/utils_heatmap_help.R     - help-modal content (existing, unchanged)

heatmap_reactives <- function(id, dataset, metacell_types, gene_modules, cell_type_colors, state, markers, lfp_range, mode, genes = NULL, highlighted_genes = NULL, height = "80vh") {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            metacell_filter <- reactiveVal()
            selected_metacells <- reactiveVal()

            # Applied parameters for batched heatmap updates
            applied_params <- reactiveVal(NULL)

            genes <- genes %||% reactiveVal()
            highlighted_genes <- highlighted_genes %||% reactiveVal()

            observeEvent(input$show_help, {
                showModal(create_heatmap_help_modal(mode))
            })

            mat <- reactive({
                req(markers())
                req(dataset())
                req(metacell_types())
                req(applied_params())
                req(is.null(applied_params()$selected_cell_types) || all(applied_params()$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

                m <- get_marker_matrix(
                    dataset(),
                    markers(),
                    applied_params()$selected_cell_types,
                    metacell_types(),
                    cell_type_colors(),
                    gene_modules(),
                    force_cell_type = applied_params()$force_cell_type %||% TRUE,
                    mode = mode,
                    notify_var_genes = TRUE,
                    metadata_order = input$metadata_order_var,
                    cell_type_metadata_order = applied_params()$metadata_order_cell_type_var,
                    recalc = input$mat_value == "Local",
                    metacells = metacell_filter()
                )


                if (input$filter_by_clipboard) {
                    if (!is.null(state$selection$clipboard) && length(state$selection$clipboard) > 0) {
                        m <- m[, intersect(colnames(m), state$selection$clipboard), drop = FALSE]
                    }
                }

                # Apply categorical filtering if enabled
                if (!is.null(applied_params()$enable_categorical_filter) && applied_params()$enable_categorical_filter &&
                    !is.null(applied_params()$categorical_filter_field) && !is.null(applied_params()$categorical_filter_values) &&
                    length(applied_params()$categorical_filter_values) > 0) {
                    metadata <- get_mc_data(dataset(), "metadata")
                    if (!is.null(metadata) && applied_params()$categorical_filter_field %in% colnames(metadata)) {
                        field_values <- metadata[[applied_params()$categorical_filter_field]]
                        valid_rows <- !is.na(field_values) & field_values != ""
                        filtered_metacells <- metadata$metacell[
                            valid_rows & field_values %in% applied_params()$categorical_filter_values
                        ]
                        m <- m[, intersect(colnames(m), filtered_metacells), drop = FALSE]
                    }
                }

                # Add genes to the matrix if exists
                if (!is.null(genes) && length(genes()) > 0 && !is.null(input$show_genes) && input$show_genes) {
                    m <- add_genes_to_marker_matrix(m, genes(), dataset())
                }

                return(m)
            }) %>% bindCache(
                id, dataset(), metacell_types(), cell_type_colors(), markers(), gene_modules(),
                # Use only mat-relevant fields from applied_params to avoid cache misses
                # when unrelated fields (e.g., selected_md) change
                applied_params()$selected_cell_types,
                applied_params()$force_cell_type,
                applied_params()$metadata_order_cell_type_var,
                applied_params()$enable_categorical_filter,
                applied_params()$categorical_filter_field,
                applied_params()$categorical_filter_values,
                genes(), input$show_genes, clipboard_changed(), mode,
                input$metadata_order_var, metacell_filter(), input$mat_value
            )

            heatmap_download_handlers(output, mat, markers, metacell_filter, dataset, input, metacell_types, state, applied_params, ns)

            heatmap_matrix_reactives(ns, input, output, session, dataset, metacell_types, cell_type_colors, state, markers, lfp_range, mode, metacell_filter, mat)

            output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors, metacell_types = metacell_types)

            output$metadata_list <- metadata_selector(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata", additional_fields = "Clipboard")

            # Initialize applied_params with default values
            observe({
                req(cell_type_colors())
                if (is.null(applied_params())) {
                    all_types <- cell_type_colors() %>%
                        filter(cell_type %in% metacell_types()$cell_type) %>%
                        pull(cell_type)
                    applied_params(list(
                        selected_cell_types = all_types,
                        selected_md = NULL,
                        force_cell_type = TRUE,
                        metadata_order_cell_type_var = "Hierarchical-Clustering",
                        enable_categorical_filter = FALSE,
                        categorical_filter_field = NULL,
                        categorical_filter_values = NULL
                    ))
                }
            })

            # Dynamic UI for apply button
            output$apply_heatmap_changes_ui <- renderUI({
                actionButton(
                    ns("apply_heatmap_changes"),
                    "Update Heatmap",
                    class = "btn-primary",
                    style = "font-weight: bold; margin-bottom: 5px; width: 85%;",
                    title = "Click to apply all filtering and display changes from the sidebar."
                )
            })

            # Observer for apply changes button
            observeEvent(input$apply_heatmap_changes, {
                applied_params(list(
                    selected_cell_types = input$selected_cell_types,
                    selected_md = input$selected_md,
                    force_cell_type = input$force_cell_type %||% TRUE,
                    metadata_order_cell_type_var = input$metadata_order_cell_type_var,
                    enable_categorical_filter = input$enable_categorical_filter,
                    categorical_filter_field = input$categorical_filter_field,
                    categorical_filter_values = input$categorical_filter_values
                ))
                showNotification("Applied changes - updating heatmap...", type = "message")
            })

            output$reset_zoom_ui <- renderUI({
                shinyWidgets::actionGroupButtons(ns("reset_zoom"), labels = "Reset zoom", size = "sm")
            })

            output$copy_metacells_ui <- renderUI({
                shinyWidgets::actionGroupButtons(ns("copy_metacells"), labels = "Copy metacells", size = "sm")
            })

            output$mat_value_ui <- renderUI({
                shinyWidgets::radioGroupButtons(
                    inputId = ns("mat_value"),
                    label = "Enrichment type:",
                    choices = c("Global", "Local"),
                    selected = "Global",
                    size = "sm",
                    justified = TRUE
                )
            })

            observe({
                shinyjs::toggle(
                    id = "reset_zoom_ui",
                    condition = !is.null(metacell_filter()) && length(metacell_filter()) > 0
                )
                shinyjs::toggle(
                    id = "mat_value_ui",
                    condition = (!is.null(metacell_filter()) && length(metacell_filter()) > 0) ||
                        (
                            !is.null(input$selected_cell_types) &&
                                length(input$selected_cell_types) < nrow(cell_type_colors())
                        ) ||
                        (
                            !is.null(input$filter_by_clipboard) &&
                                input$filter_by_clipboard &&
                                !is.null(state$selection$clipboard) &&
                                length(state$selection$clipboard) > 0
                        )
                )
                shinyjs::toggle(
                    id = "copy_metacells_ui",
                    condition = !is.null(selected_metacells()) &&
                        length(selected_metacells()) > 0 &&
                        !is.null(input$brush_action) &&
                        input$brush_action == "Select"
                )
            })

            observeEvent(input$reset_zoom, {
                metacell_filter(NULL)
            })

            observeEvent(input$copy_metacells, {
                state$selection$clipboard <- selected_metacells()
                showNotification(glue("Copied {length(selected_metacells())} metacells to clipboard"))
                selected_metacells(character(0))
            })

            observeEvent(input$highlight_genes, {
                req(input$selected_marker_genes)
                highlighted_genes(input$selected_marker_genes)

                showNotification(
                    glue::glue("Highlighted {length(input$selected_marker_genes)} gene(s)"),
                    type = "message"
                )
            })

            observeEvent(input$clear_highlights, {
                highlighted_genes(NULL)
                showNotification("Cleared gene highlights", type = "message")
            })

            output$plotting_area <- renderUI({
                heatmap <- plotOutput(
                    ns("heatmap"),
                    height = height,
                    dblclick = dblclickOpts(ns("heatmap_dblclick"), clip = TRUE),
                    hover = hoverOpts(ns("heatmap_hover"), delay = 250, delayType = "debounce"),
                    brush = brushOpts(
                        id = ns("heatmap_brush"),
                        direction = "x",
                        resetOnNew = TRUE,
                        delay = 1000
                    )
                )


                if (!is.null(input$plot_legend) && input$plot_legend) {
                    req(input$legend_width)
                    legend_column <- column(
                        width = input$legend_width,
                        plotOutput(ns("markers_legend"), height = height)
                    )
                    heatmap_column <- column(
                        width = 12 - input$legend_width,
                        heatmap
                    )
                    shinycssloaders::withSpinner(
                        fluidRow(heatmap_column, legend_column)
                    )
                } else {
                    shinycssloaders::withSpinner(
                        heatmap
                    )
                }
            })

            output$markers_legend <- renderPlot(
                {
                    req(cell_type_colors())
                    req(input$plot_legend)
                    colors <- cell_type_colors()
                    colors <- colors %>% filter(cell_type %in% metacell_types()$cell_type)
                    legend_point_size <- max(1, min(2, 250 / nrow(colors)))
                    p <- colors %>%
                        ggplot(aes(x = cell_type, fill = cell_type, y = 1))

                    if (input$plot_cell_type_legend) {
                        p <- p + geom_point(shape = 21) +
                            scale_fill_manual("Cell types", values = deframe(colors %>% select(cell_type, color))) +
                            guides(fill = guide_legend(override.aes = list(size = legend_point_size), ncol = 1))
                    }

                    if (input$plot_genes_legend) {
                        gene_colors <- data.frame(type = c("lateral+noisy", "lateral", "noisy", "disjoined", "module", "highlighted"), color = c("purple", "blue", "red", "darkgray", "#012901", input$highlight_color))
                        p <- p +
                            geom_text(data = gene_colors, inherit.aes = FALSE, x = 1, y = 1, aes(label = type, color = type)) +
                            scale_color_manual("Genes", values = deframe(gene_colors))
                    }

                    p <- p + theme(legend.position = c(0.5, 0.5))

                    legend <- cowplot::get_legend(p)

                    cowplot::ggdraw(legend)
                },
                res = 96
            ) %>% bindCache(id, dataset(), cell_type_colors(), input$plot_legend, input$plot_cell_type_legend, input$plot_genes_legend)

            # We use this reactive in order to invalidate the cache only when input$filter_by_clipboard is TRUE
            clipboard_changed <- reactive({
                if ((!is.null(input$filter_by_clipboard) && input$filter_by_clipboard) || (!is.null(applied_params()) && !is.null(applied_params()$selected_md) && "Clipboard" %in% applied_params()$selected_md)) {
                    return(c(input$filter_by_clipboard, state$selection$clipboard))
                } else {
                    return(FALSE)
                }
            })

            output$heatmap <- renderPlot(
                {
                    req(dataset())
                    req(lfp_range())
                    req(metacell_types())
                    req(cell_type_colors())
                    req(input$midpoint)
                    req(input$low_color)
                    req(input$high_color)
                    req(input$mid_color)
                    req(input$midpoint > lfp_range()[1])
                    req(input$midpoint < lfp_range()[2])

                    m <- mat()
                    req(m)

                    m <- filter_heatmap_by_metacell(m, metacell_filter())
                    metadata <- get_markers_metadata(dataset, applied_params()$selected_md, metacell_types, state)

                    req(nrow(m) > 0)
                    req(ncol(m) > 0)
                    gene_colors <- get_gene_colors(
                        rownames(m),
                        lateral_genes = get_mc_data(dataset(), "lateral_genes"),
                        disjoined_genes = get_mc_data(dataset(), "disjoined_genes_no_atlas"),
                        noisy_genes = get_mc_data(dataset(), "noisy_genes")
                    )

                    if (!is.null(genes) && length(genes()) > 0 && !is.null(input$show_genes) && input$show_genes) {
                        gene_colors <- ifelse(names(gene_colors) %in% genes(), gene_colors, "#012901")
                    }

                    if (!is.null(highlighted_genes) && length(highlighted_genes()) > 0) {
                        valid_genes <- highlighted_genes()[highlighted_genes() %in% names(gene_colors)]

                        if (length(valid_genes) > 0) {
                            gene_colors[valid_genes] <- input$highlight_color
                        }
                    }


                    res <- plot_markers_mat(
                        m,
                        metacell_types(),
                        cell_type_colors(),
                        dataset(),
                        min_lfp = lfp_range()[1],
                        max_lfp = lfp_range()[2],
                        plot_legend = FALSE,
                        high_color = input$high_color,
                        low_color = input$low_color,
                        mid_color = input$mid_color,
                        midpoint = input$midpoint,
                        metadata = metadata,
                        gene_colors = gene_colors,
                        col_names = ncol(m) <= 100,
                        top_cell_type_bar = ncol(m) <= 100,
                        interleave = nrow(m) > 80,
                        vertical_gridlines = mode %in% c("Proj"),
                        separate_gtable = TRUE
                    )

                    # we are returning the gtable and ggplot object separatly in order to allow shiny to infer positions correctly.
                    return(structure(list(p = res$p, gtable = res$gtable), class = "gt_custom"))
                },
                res = 96
            ) %>% bindCache(id, dataset(), metacell_types(), cell_type_colors(), gene_modules(), lfp_range(), metacell_filter(), applied_params(), markers(), clipboard_changed(), input$high_color, input$low_color, input$mid_color, input$midpoint, genes(), input$show_genes, highlighted_genes(), input$highlight_color, input$metadata_order_var, input$mat_value, mode)

            observeEvent(input$heatmap_brush, {
                req(input$brush_action)
                m <- filter_heatmap_by_metacell(mat(), metacell_filter())
                range <- max(1, round(input$heatmap_brush$xmin)):min(ncol(m), round(input$heatmap_brush$xmax))
                req(all(range >= 1) && all(range <= ncol(m)))
                metacells <- colnames(m)[range]
                if (input$brush_action == "Zoom") {
                    metacell_filter(metacells)
                } else if (input$brush_action == "Select") {
                    selected_metacells(metacells)
                }
            })

            observeEvent(input$heatmap_dblclick, {
                m <- filter_heatmap_by_metacell(mat(), metacell_filter())
                gene <- get_gene_by_heatmap_coord(m, input$heatmap_dblclick$y)
                metacell <- get_metacell_by_heatmap_coord(m, input$heatmap_dblclick$x)

                if (gene %in% gene_names(dataset())) {
                    if (input$gene_select == "X axis") {
                        state$selection$selected_gene_x_axis <- gene
                    } else {
                        state$selection$selected_gene_y_axis <- gene
                    }
                    state$selection$selected_query_gene <- gene
                }

                if (input$metacell_select == "Metacell A") {
                    state$selection$selected_metacellA <- metacell
                } else {
                    state$selection$selected_metacellB <- metacell
                }

                state$selection$selected_query_metacell <- metacell
            })

            heatmap_tooltip_handler(output, mat, metacell_filter, metacell_types, cell_type_colors, dataset, input, state, mode, gene_modules, genes, applied_params)
        }
    )
}


