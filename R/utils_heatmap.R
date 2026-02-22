heatmap_matrix_reactives <- function(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, mode, metacell_filter, mat) {
    observe({
        choices <- markers()
        if (!is.null(choices)) {
            names(choices) <- gene_label(choices, dataset())
            updateSelectInput(session, "selected_marker_genes", choices = choices)
        }
        shinyWidgets::updateVirtualSelect(session = session, inputId = "metadata_order_var", choices = c("Hierarchical-Clustering", dataset_metadata_fields_numeric(dataset())), selected = "Hierarchical-Clustering")
        shinyWidgets::updateVirtualSelect(session = session, inputId = "metadata_order_cell_type_var", choices = c("Hierarchical-Clustering", "Colors table", dataset_metadata_fields_numeric(dataset())), selected = "Hierarchical-Clustering")

        shinyjs::toggle(id = "metadata_order_cell_type_var", condition = input$force_cell_type)
    })

    observe({
        if (is.null(markers())) {
            config <- mcv_get("config")
            if (config$light_version) {
                max_gene_num <- 100
            } else {
                req(input$max_gene_num)
                max_gene_num <- input$max_gene_num
            }

            initial_markers <- choose_markers(get_marker_genes(dataset(), mode = mode), max_gene_num, dataset = dataset())
            markers(initial_markers)
        }

        lfp_range(input$lfp_range)
    })

    observeEvent(input$load_genes, {
        req(input$load_genes)
        req(input$load_genes$datapath)
        req(input$load_genes$datapath != "")
        new_markers <- data.table::fread(input$load_genes$datapath, header = FALSE, data.table = FALSE)[, 1]
        new_markers <- new_markers[new_markers %in% gene_names(dataset())]
        if (length(new_markers) == 0) {
            showNotification("No valid genes were loaded", type = "warning")
        } else {
            markers(new_markers)
            showNotification(glue("Loaded {length(new_markers)} genes"), type = "message")
        }
    })

    observeEvent(input$use_de_genes, {
        if (!is.null(globals$significant_genes) && length(globals$significant_genes) > 0) {
            markers(globals$significant_genes)
            showNotification(glue("Updated markers with {length(globals$significant_genes)} significant genes from current Diff. Expression comparison"), type = "message")
        } else {
            showNotification("No significant genes available in the current Diff. Expression comparison.", type = "warning")
        }
    })

    observeEvent(input$update_genes, {
        req(metacell_types())
        req(cell_type_colors())
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

        if (input$mat_value == "Local") {
            mc_egc <- get_mc_egc(dataset(), metacells = colnames(mat()))
            markers_df <- calc_marker_genes(mc_egc, genes_per_metacell = 20)
        } else {
            if (!is.null(input$selected_cell_types)) {
                markers_df <- metacell_types() %>%
                    filter(cell_type %in% input$selected_cell_types)
            } else {
                markers_df <- metacell_types()
            }

            if (!is.null(input$use_markers) && input$use_markers) {
                new_markers_df <- get_marker_genes(dataset(), mode = "Markers")
            } else {
                new_markers_df <- get_marker_genes(dataset(), mode = mode)
            }

            if (has_name(new_markers_df, "metacell") && mode != "Outliers") {
                markers_df <- markers_df %>%
                    select(metacell) %>%
                    inner_join(new_markers_df, by = "metacell")
            } else {
                markers_df <- new_markers_df
            }

            if (!is.null(input$filter_by_clipboard) && input$filter_by_clipboard) {
                if (!is.null(globals$clipboard) && length(globals$clipboard) > 0) {
                    markers_df <- markers_df %>%
                        filter(metacell %in% globals$clipboard)
                }
            }
        }

        req(input$max_gene_num)
        markers_df <- filter_genes_by_flags(
            markers_df,
            lateral_genes = get_mc_data(dataset(), "lateral_genes"),
            noisy_genes = get_mc_data(dataset(), "noisy_genes"),
            include_lateral = is.null(input$include_lateral) || input$include_lateral,
            include_noisy = is.null(input$include_noisy) || input$include_noisy
        )
        if (!is.null(input$show_only_fitted) && input$show_only_fitted) {
            gene_metadata <- get_mc_data(dataset(), "gene_metadata")

            if (!is.null(gene_metadata)) {
                fitted_genes <- gene_metadata %>%
                    filter(fitted_gene) %>%
                    pull(gene) %>%
                    unique()
                markers_df <- markers_df %>%
                    filter(gene %in% fitted_genes)
            }
        }

        new_markers <- choose_markers(markers_df, input$max_gene_num, dataset = dataset())

        # If we did not choose all the cell types
        if (!is.null(input$selected_cell_types) && length(input$selected_cell_types) != nrow(cell_type_colors())) {
            # TODO: add also genes that are distictive for the specific cell type (although not variable within it)
        }

        markers(new_markers)
    })

    observeEvent(input$remove_genes, {
        req(markers())
        new_markers <- markers()[!(markers() %in% input$selected_marker_genes)]
        markers(new_markers)
        shinyWidgets::updateVirtualSelect(session = session, inputId = "genes_to_add", selected = c())
    })

    output$add_genes_ui <- renderUI({
        if (mode %in% c("Inner", "Outliers", "Stdev")) {
            if (mode == "Inner") {
                mc_fp <- get_mc_data(dataset(), "inner_fold_mat")
            } else if (mode == "Stdev") {
                mc_fp <- get_mc_data(dataset(), "inner_stdev_mat")
            } else {
                mc_fp <- get_mc_data(dataset(), "deviant_fold_mat")
            }
            req(mc_fp)
            gene_choices <- rownames(mc_fp)[Matrix::rowSums(mc_fp) > 0]
            req(markers)
            gene_choices <- gene_choices[!(gene_choices %in% markers())]
        } else {
            gene_choices <- gene_names(dataset())
        }


        tagList(
            shinyWidgets::virtualSelectInput(ns("genes_to_add"),
                label = "Add genes",
                choices = gene_choices,
                selected = c(),
                multiple = TRUE,
                showSelectedOptionsFirst = TRUE,
                search = TRUE,
                markSearchResults = TRUE,
                searchByStartsWith = TRUE,
                disableSelectAll = TRUE
            ),
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm")
        )
    })

    observeEvent(input$add_genes, {
        new_markers <- sort(unique(c(markers(), input$genes_to_add)))
        markers(new_markers)
        shinyWidgets::updateVirtualSelect(session = session, inputId = "genes_to_add", selected = character(0))
    })
}


heatmap_reactives <- function(id, dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, mode, genes = NULL, highlighted_genes = NULL, height = "80vh") {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            metacell_filter <- reactiveVal()
            selected_metacells <- reactiveVal()

            genes <- genes %||% reactiveVal()
            highlighted_genes <- highlighted_genes %||% reactiveVal()

            observeEvent(input$show_help, {
                showModal(create_heatmap_help_modal(mode))
            })

            mat <- reactive({
                req(markers())
                req(dataset())
                req(metacell_types())
                req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

                m <- get_marker_matrix(
                    dataset(),
                    markers(),
                    input$selected_cell_types,
                    metacell_types(),
                    cell_type_colors(),
                    gene_modules(),
                    force_cell_type = input$force_cell_type %||% TRUE,
                    mode = mode,
                    notify_var_genes = TRUE,
                    metadata_order = input$metadata_order_var,
                    cell_type_metadata_order = input$metadata_order_cell_type_var,
                    recalc = input$mat_value == "Local",
                    metacells = metacell_filter()
                )


                if (input$filter_by_clipboard) {
                    if (!is.null(globals$clipboard) && length(globals$clipboard) > 0) {
                        m <- m[, intersect(colnames(m), globals$clipboard), drop = FALSE]
                    }
                }

                # Add genes to the matrix if exists
                if (!is.null(genes) && length(genes()) > 0 && !is.null(input$show_genes) && input$show_genes) {
                    m <- add_genes_to_marker_matrix(m, genes(), dataset())
                }

                return(m)
            }) %>% bindCache(id, dataset(), metacell_types(), cell_type_colors(), markers(), gene_modules(), input$selected_cell_types, input$force_cell_type, genes(), input$show_genes, clipboard_changed(), mode, input$metadata_order_var, input$metadata_order_cell_type_var, metacell_filter(), input$mat_value)

            heatmap_download_handlers(output, mat, markers, metacell_filter, dataset, input, metacell_types, globals)

            heatmap_matrix_reactives(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, mode, metacell_filter, mat)

            output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors, metacell_types = metacell_types)

            output$metadata_list <- metadata_selector(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata", additional_fields = "Clipboard")


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
                                !is.null(globals$clipboard) &&
                                length(globals$clipboard) > 0
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
                globals$clipboard <- selected_metacells()
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
                if ((!is.null(input$filter_by_clipboard) && input$filter_by_clipboard) || (!is.null(input$selected_md) && "Clipboard" %in% input$selected_md)) {
                    return(c(input$filter_by_clipboard, globals$clipboard))
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
                    metadata <- get_markers_metadata(dataset, input, metacell_types, globals)

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
                        vertical_gridlines = mode %in% c("Inner", "Proj", "Stdev"),
                        separate_gtable = TRUE
                    )

                    # we are returning the gtable and ggplot object separatly in order to allow shiny to infer positions correctly.
                    return(structure(list(p = res$p, gtable = res$gtable), class = "gt_custom"))
                },
                res = 96
            ) %>% bindCache(id, dataset(), metacell_types(), cell_type_colors(), gene_modules(), lfp_range(), metacell_filter(), input$plot_legend, input$plot_cell_type_legend, input$plot_genes_legend, input$selected_md, markers(), input$selected_cell_types, input$force_cell_type, clipboard_changed(), input$high_color, input$low_color, input$mid_color, input$midpoint, genes(), input$show_genes, highlighted_genes(), input$highlight_color, input$max_gene_num, input$metadata_order_var, input$metadata_order_cell_type_var, input$mat_value)

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
                        globals$selected_gene_x_axis <- gene
                    } else {
                        globals$selected_gene_y_axis <- gene
                    }
                    globals$selected_query_gene <- gene
                }

                if (input$metacell_select == "Metacell A") {
                    globals$selected_metacellA <- metacell
                } else {
                    globals$selected_metacellB <- metacell
                }

                globals$selected_query_metacell <- metacell
            })

            heatmap_tooltip_handler(output, mat, metacell_filter, metacell_types, cell_type_colors, dataset, input, globals, mode, gene_modules, genes)
        }
    )
}


#' Build the hover tooltip UI for the heatmap
#'
#' Extracted from heatmap_reactives to reduce function complexity.
#' Renders the hover_info output that displays gene/metacell details on hover.
#'
#' @param output Shiny output object
#' @param mat Reactive expression returning the heatmap matrix
#' @param metacell_filter Reactive value for metacell filtering
#' @param metacell_types Reactive expression returning metacell type annotations
#' @param cell_type_colors Reactive expression returning cell type color mappings
#' @param dataset Reactive expression returning the current dataset
#' @param input Shiny input object
#' @param globals Reactive values for global state
#' @param mode Character string indicating heatmap mode
#' @param gene_modules Reactive expression returning gene module data
#' @param genes Reactive expression returning selected genes (can be NULL)
#' @noRd
heatmap_tooltip_handler <- function(output, mat, metacell_filter, metacell_types, cell_type_colors, dataset, input, globals, mode, gene_modules, genes) {
    output$hover_info <- renderUI({
        m <- mat()
        req(m)
        req(input$heatmap_hover)
        req(metacell_types())
        req(cell_type_colors())

        hover <- input$heatmap_hover
        m <- filter_heatmap_by_metacell(m, metacell_filter())
        gene <- get_gene_by_heatmap_coord(m, hover$y)
        metacell <- get_metacell_by_heatmap_coord(m, hover$x)
        value <- m[gene, metacell]

        req(metacell)
        req(gene)

        # taken from https://gitlab.com/-/snippets/16220
        left_px <- hover$coords_css$x
        top_px <- hover$coords_css$y

        mcell_stats <- metacell_types() %>%
            filter(metacell == !!metacell)

        style <- glue(
            "position:absolute; z-index:100; left: {left_px + 2}px; top: {top_px + 2}px;"
        )

        top_genes <- mcell_stats %>%
            glue::glue_data("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")

        if (mode == "Outliers") {
            mcell_tooltip <- paste(
                glue("Gene: {gene_name}", gene_name = gene_label(gene, dataset(), gene_modules())),
                glue("Cell: {metacell}"),
                glue("Value: {round(value, digits=2)}"),
                glue("Most similar metacell: {mcell_stats$most_similar_metacell}"),
                glue("Cell type: {mcell_stats$cell_type}"),
                glue("Top genes: {top_genes}"),
                ifelse(has_name(mcell_stats, "mc_age"), glue("Metacell age (E[t]): {round(mcell_stats$mc_age, digits=2)}"), ""),
                sep = "<br/>"
            )
        } else {
            gene_prefix <- "Gene"
            if (mode == "Gene modules" && !(gene %in% genes())) {
                gene_prefix <- "Gene module"
            }

            gene_name <- gene_label(gene, dataset(), gene_modules())
            if (mode == "Gene modules") {
                gene_name <- gene
            }

            mcell_tooltip <- paste(
                glue("{gene_prefix}: {gene_name}"),
                glue("Metacell: {metacell}"),
                glue("Value: {round(value, digits=2)}"),
                glue("Cell type: {mcell_stats$cell_type}"),
                glue("Top genes: {top_genes}"),
                ifelse(has_name(mcell_stats, "mc_age"), glue("Metacell age (E[t]): {round(mcell_stats$mc_age, digits=2)}"), ""),
                sep = "<br/>"
            )

            metadata <- get_markers_metadata(dataset, input, metacell_types, globals)
            if (!is.null(metadata)) {
                mc_md <- metadata %>%
                    filter(metacell == !!metacell) %>%
                    select(-metacell)
                md_tooltip <- purrr::map_chr(colnames(mc_md), ~
                    glue("{.x}: {mc_md[[.x]][1]}")) %>%
                    paste(collapse = "<br/>")
                mcell_tooltip <- paste(mcell_tooltip, md_tooltip)
            }
        }

        wellPanel(
            style = style,
            p(HTML(mcell_tooltip))
        )
    })
}


#' Set up download handlers for the heatmap
#'
#' Extracted from heatmap_reactives to reduce function complexity.
#' Creates download handlers for the marker matrix and gene list.
#'
#' @param output Shiny output object
#' @param mat Reactive expression returning the heatmap matrix
#' @param markers Reactive value containing the current marker genes
#' @param metacell_filter Reactive value for metacell filtering
#' @param dataset Reactive expression returning the current dataset
#' @param input Shiny input object
#' @param metacell_types Reactive expression returning metacell type annotations
#' @param globals Reactive values for global state
#' @noRd
heatmap_download_handlers <- function(output, mat, markers, metacell_filter, dataset, input, metacell_types, globals) {
    output$download_matrix <- downloadHandler(
        filename = function() {
            paste("markers_matrix-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            m <- mat()[rev(1:nrow(mat())), , drop = FALSE]
            if (length(metacell_filter()) > 0) {
                m <- m[, intersect(colnames(m), metacell_filter()), drop = FALSE]
            }
            if (input$include_metadata) {
                metadata <- get_markers_metadata(dataset, input, metacell_types, globals)
                metadata_m <- metadata %>%
                    as.data.frame() %>%
                    column_to_rownames("metacell") %>%
                    t() %>%
                    as.matrix()
                ct_m <- metacell_types() %>%
                    select(metacell, cell_type) %>%
                    deframe()
                ct_m <- ct_m[colnames(m)]
                metadata_m <- rbind(t(as.matrix(ct_m)), metadata_m)
                rownames(metadata_m)[1] <- "Cell type"
                m <- rbind(m, metadata_m)
            }
            fwrite(
                m %>%
                    as.data.frame() %>%
                    rownames_to_column("gene"),
                file,
                row.names = FALSE
            )
        }
    )

    output$download_genes <- downloadHandler(
        filename = function() {
            paste("markers-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            fwrite(
                data.frame(gene = markers()),
                file,
                row.names = FALSE,
                col.names = FALSE
            )
        }
    )
}
