heatmap_box <- function(id,
                        title = "Heatmap",
                        fold_change_range = c(-3, 3),
                        midpoint = 0,
                        low_color = "blue",
                        mid_color = "white",
                        high_color = "red",
                        gene_select_label = "Select on double-click (gene)",
                        gene_select_choices = c("X axis", "Y axis"),
                        legend_width = 2,
                        height = "80vh") {
    ns <- NS(id)
    div(
        generic_box(
            id = ns("heatmap_box"),
            title = title,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            height = height,
            sidebar = shinydashboardPlus::boxSidebar(
                startOpen = FALSE,
                width = 25,
                id = ns("heatmap_sidebar"),
                checkboxInput(ns("force_cell_type"), "Force cell type", value = TRUE),
                checkboxInput(ns("filter_by_clipboard"), "Filter by clipboard", value = FALSE),
                shinyWidgets::numericRangeInput(ns("lfp_range"), "Fold change range", fold_change_range, width = "80%", separator = " to "),
                numericInput(ns("midpoint"), "Midpoint", midpoint),
                colourpicker::colourInput(ns("low_color"), "Low color", low_color),
                colourpicker::colourInput(ns("high_color"), "High color", high_color),
                colourpicker::colourInput(ns("mid_color"), "Mid color", mid_color),
                checkboxInput(ns("plot_legend"), "Plot legend", value = TRUE),
                numericInput(ns("legend_width"), "Legend width", min = 1, max = 11, step = 1, value = legend_width),
                shinyWidgets::prettyRadioButtons(
                    inputId = ns("gene_select"),
                    label = gene_select_label,
                    choices = gene_select_choices,
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                ),
                shinyWidgets::prettyRadioButtons(
                    inputId = ns("metacell_select"),
                    label = "Select on double-click (metacell)",
                    choices = c("Metacell A", "Metacell B"),
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                )
            ),
            uiOutput(ns("plotting_area"))
        ),
        style = "position:relative",
        uiOutput(ns("hover_info"), style = "pointer-events: none")
    )
}

heatmap_sidebar <- function(id, ...) {
    ns <- NS(id)
    if (config$light_version) {
        max_gene_num_ui <- NULL
        remove_genes_ui <- NULL
        add_genes_ui <- NULL
        update_genes_ui <- NULL
    } else {
        max_gene_num_ui <- numericInput(ns("max_gene_num"), "Maximal number of genes", value = 100)
        remove_genes_ui <- shinyWidgets::actionGroupButtons(ns("remove_genes"), labels = "Remove selected genes", size = "sm")
        add_genes_ui <- uiOutput(ns("add_genes_ui"))
        update_genes_ui <- shinyWidgets::actionGroupButtons(ns("update_genes"), labels = "Update genes", size = "sm")
    }

    list(
        uiOutput(ns("reset_zoom_ui")),
        shinyWidgets::radioGroupButtons(
            inputId = ns("brush_action"),
            label = "Brush action:",
            choices = c("Zoom", "Select"),
            selected = "Zoom",
            size = "sm",
            justified = TRUE
        ),
        uiOutput(ns("copy_metacells_ui")),
        tags$hr(),
        uiOutput(ns("cell_type_list")),
        uiOutput(ns("metadata_list")),
        ...,
        update_genes_ui,
        max_gene_num_ui,
        add_genes_ui,
        selectInput(
            ns("selected_marker_genes"),
            "Genes",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            size = 30,
            selectize = FALSE
        ),
        remove_genes_ui,
        tags$hr(),
        downloadButton(ns("download_matrix"), "Download matrix", align = "center", style = "margin: 5px 5px 5px 15px; ")
    )
}

heatmap_matrix_reactives <- function(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, mode) {
    observe({
        updateSelectInput(session, "selected_marker_genes", choices = markers())
    })

    observe({
        if (is.null(markers())) {
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

    observeEvent(input$update_genes, {
        req(metacell_types())
        req(cell_type_colors())
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

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

        req(input$max_gene_num)
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
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm"),
            shinyWidgets::virtualSelectInput(ns("genes_to_add"),
                label = "Add genes",
                choices = gene_choices,
                selected = c(),
                multiple = TRUE,
                showSelectedOptionsFirst = TRUE,
                search = TRUE
            )
        )
    })

    observeEvent(input$add_genes, {
        new_markers <- sort(unique(c(markers(), input$genes_to_add)))
        markers(new_markers)
        shinyWidgets::updateVirtualSelect(session = session, inputId = "genes_to_add", selected = character(0))
    })
}


heatmap_reactives <- function(id, dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, mode, genes = NULL, highlighted_genes = NULL, highlight_color = "red", height = "80vh") {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            metacell_filter <- reactiveVal()
            selected_metacells <- reactiveVal()

            genes <- genes %||% reactiveVal()
            highlighted_genes <- highlighted_genes %||% reactiveVal()

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
                    gene_modules(),
                    force_cell_type = input$force_cell_type,
                    mode = mode,
                    notify_var_genes = TRUE
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
            }) %>% bindCache(id, dataset(), metacell_types(), cell_type_colors(), markers(), gene_modules(), input$selected_cell_types, input$force_cell_type, genes(), input$show_genes, clipboard_changed(), mode)

            output$download_matrix <- downloadHandler(
                filename = function() {
                    paste("markers_matrix-", Sys.Date(), ".csv", sep = "")
                },
                content = function(file) {
                    m <- mat()[rev(1:nrow(mat())), , drop = FALSE]
                    if (length(metacell_filter()) > 0) {
                        m <- m[, intersect(colnames(m), metacell_filter()), drop = FALSE]
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


            heatmap_matrix_reactives(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, mode)

            output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors)

            output$metadata_list <- metadata_selector(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata", additional_fields = "Clipboard")


            output$reset_zoom_ui <- renderUI({
                shinyWidgets::actionGroupButtons(ns("reset_zoom"), labels = "Reset zoom", size = "sm")
            })

            output$copy_metacells_ui <- renderUI({
                shinyWidgets::actionGroupButtons(ns("copy_metacells"), labels = "Copy metacells", size = "sm")
            })

            observe({
                shinyjs::toggle(id = "reset_zoom_ui", condition = !is.null(metacell_filter()) && length(metacell_filter()) > 0)
                shinyjs::toggle(id = "copy_metacells_ui", condition = !is.null(selected_metacells()) && length(selected_metacells()) > 0 && !is.null(input$brush_action) && input$brush_action == "Select")
            })

            observeEvent(input$reset_zoom, {
                metacell_filter(NULL)
            })

            observeEvent(input$copy_metacells, {
                globals$clipboard <- selected_metacells()
                showNotification(glue("Copied {length(selected_metacells())} metacells to clipboard"))
                selected_metacells(character(0))
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
                    legend_point_size <- max(1, min(2, 250 / nrow(cell_type_colors())))
                    legend <- cowplot::get_legend(cell_type_colors() %>%
                        ggplot(aes(x = cell_type, color = cell_type, y = 1)) +
                        geom_point() +
                        scale_color_manual("", values = deframe(cell_type_colors() %>% select(cell_type, color))) +
                        guides(color = guide_legend(override.aes = list(size = legend_point_size), ncol = 1)) +
                        theme(legend.position = c(0.5, 0.5)))
                    cowplot::ggdraw(legend)
                },
                res = 96
            ) %>% bindCache(id, dataset(), cell_type_colors(), input$plot_legend)

            # We use this reactive in order to invalidate the cache only when input$filter_by_clipboard is TRUE
            clipboard_changed <- reactive({
                if ((!is.null(input$filter_by_clipboard) && input$filter_by_clipboard) || (!is.null(input$selected_md) && input$selected_md == "Clipboard")) {
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

                    if (!is.null(input$selected_md)) {
                        metadata <- get_mc_data(dataset(), "metadata")
                        if (is.null(metadata)) {
                            metadata <- metacell_types() %>% select(metacell)
                        }
                        metadata <- metadata %>%
                            mutate(Clipboard = ifelse(metacell %in% globals$clipboard, "selected", "not selected")) %>%
                            select(metacell, one_of(input$selected_md))
                    } else {
                        metadata <- NULL
                    }

                    req(nrow(m) > 0)
                    req(ncol(m) > 0)

                    if (!is.null(genes) && length(genes()) > 0 && !is.null(input$show_genes) && input$show_genes) {
                        # m <- add_genes_to_marker_matrix(m, genes(), dataset())
                        other_genes <- genes()

                        gene_colors <- tibble(gene = rownames(m), color = ifelse(gene %in% genes(), "blue", "black")) %>% deframe()
                    } else {
                        gene_colors <- get_gene_colors(
                            rownames(m),
                            lateral_genes = get_mc_data(dataset(), "lateral_genes"),
                            disjoined_genes = get_mc_data(dataset(), "disjoined_genes_no_atlas")
                        )
                    }

                    if (!is.null(highlighted_genes) && length(highlighted_genes()) > 0 && highlighted_genes() %in% names(gene_colors)) {
                        gene_colors[highlighted_genes()] <- highlight_color
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
                        vertial_gridlines = mode %in% c("Inner", "Proj", "Stdev"),
                        separate_gtable = TRUE
                    )

                    # we are returning the gtable and ggplot object separatly in order to allow shiny to infer positions correctly.
                    return(structure(list(p = res$p, gtable = res$gtable), class = "gt_custom"))
                },
                res = 96
            ) %>% bindCache(id, dataset(), metacell_types(), cell_type_colors(), gene_modules(), lfp_range(), metacell_filter(), input$plot_legend, input$selected_md, markers(), input$selected_cell_types, input$force_cell_type, clipboard_changed(), input$high_color, input$low_color, input$mid_color, input$midpoint, genes(), input$show_genes, highlighted_genes(), highlight_color)

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
                }

                wellPanel(
                    style = style,
                    p(HTML(mcell_tooltip))
                )
            })
        }
    )
}


#' Print method for "gt_custom" class
#' we are overriding the print function, similiar to the one at shiny/render-plot.R:
#' https://github.com/rstudio/shiny/blob/main/R/render-plot.R
#' in order to provide a separate ggplot_build and gtable objects
#'
#' @param x An object of class "gt_custom"
#'
#' @export print.gt_custom
#' @export
print.gt_custom <- function(x) {
    build <- ggplot_build(x$p)

    grid::grid.newpage()
    grid::grid.draw(x$gtable)

    structure(list(
        build = build,
        gtable = x$gtable
    ), class = "ggplot_build_gtable")
}

filter_heatmap_by_metacell <- function(m, f) {
    if (!is.null(f) && length(f) > 0) {
        f <- f[f %in% colnames(m)]
        m <- m[, f, drop = FALSE]
    }
    return(m)
}

get_gene_by_heatmap_coord <- function(m, coord) {
    y_coord <- round(coord)
    req(y_coord > 0 & y_coord <= nrow(m))
    return(rownames(m)[y_coord])
}

get_metacell_by_heatmap_coord <- function(m, coord) {
    x_coord <- round(coord)
    req(x_coord > 0 & x_coord <= ncol(m))
    return(colnames(m)[x_coord])
}
