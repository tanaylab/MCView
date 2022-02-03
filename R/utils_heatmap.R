heatmap_box <- function(ns,
                        title = "Heatmap",
                        fold_change_range = c(-3, 3),
                        midpoint = 0,
                        low_color = "blue",
                        mid_color = "white",
                        high_color = "red",
                        gene_select_label = "Select on double-click (gene)",
                        gene_select_choices = c("X axis", "Y axis"),
                        id_prefix = "") {
    div(
        shinydashboardPlus::box(
            id = ns(glue("{id_prefix}heatmap_box")),
            title = title,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            height = "80vh",
            sidebar = shinydashboardPlus::boxSidebar(
                startOpen = FALSE,
                width = 25,
                id = ns(glue("{id_prefix}heatmap_sidebar")),
                checkboxInput(ns(glue("{id_prefix}force_cell_type")), "Force cell type", value = TRUE),
                shinyWidgets::numericRangeInput(ns(glue("{id_prefix}lfp_range")), "Fold change range", fold_change_range, width = "80%", separator = " to "),
                numericInput(ns(glue("{id_prefix}midpoint")), "Midpoint", midpoint),
                colourpicker::colourInput(ns(glue("{id_prefix}low_color")), "Low color", low_color),
                colourpicker::colourInput(ns(glue("{id_prefix}high_color")), "High color", high_color),
                colourpicker::colourInput(ns(glue("{id_prefix}mid_color")), "Mid color", mid_color),
                checkboxInput(ns(glue("{id_prefix}plot_legend")), "Plot legend", value = TRUE),
                numericInput(ns(glue("{id_prefix}legend_width")), "Legend width", min = 1, max = 11, step = 1, value = 2),
                shinyWidgets::prettyRadioButtons(
                    inputId = ns(glue("{id_prefix}gene_select")),
                    label = gene_select_label,
                    choices = gene_select_choices,
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                ),
                shinyWidgets::prettyRadioButtons(
                    inputId = ns(glue("{id_prefix}metacell_select")),
                    label = "Select on double-click (metacell)",
                    choices = c("Metacell A", "Metacell B"),
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                )
            ),
            uiOutput(ns(glue("{id_prefix}plotting_area")))
        ),
        style = "position:relative",
        uiOutput(ns(glue("{id_prefix}hover_info")), style = "pointer-events: none")
    )
}

heatmap_sidebar <- function(ns) {
    list(
        uiOutput(ns("reset_zoom_ui")),
        uiOutput(ns("cell_type_list")),
        uiOutput(ns("metadata_list")),
        shinyWidgets::actionGroupButtons(ns("update_genes"), labels = "Update genes", size = "sm"),
        numericInput(ns("max_gene_num"), "Maximal number of genes", value = 100),
        uiOutput(ns("add_genes_ui")),
        uiOutput(ns("marker_genes_list"))
    )
}

heatmap_matrix_reactives <- function(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, mode) {
    output$marker_genes_list <- renderUI({
        tagList(
            selectInput(
                ns("selected_marker_genes"),
                "Genes",
                choices = markers(),
                selected = NULL,
                multiple = TRUE,
                size = 30,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(ns("remove_genes"), labels = "Remove selected genes", size = "sm")
        )
    })

    observe({
        if (is.null(markers())) {
            req(input$max_gene_num)
            initial_markers <- choose_markers(get_marker_genes(dataset(), mode = mode), input$max_gene_num, dataset = dataset(), add_systematic = mode == "Proj")
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

        new_markers_df <- get_marker_genes(dataset(), mode = mode)
        if (has_name(new_markers_df, "metacell")) {
            markers_df <- markers_df %>%
                select(metacell) %>%
                inner_join(new_markers_df, by = "metacell")
        } else {
            markers_df <- new_markers_df
        }

        req(input$max_gene_num)
        new_markers <- choose_markers(markers_df, input$max_gene_num, dataset = dataset(), add_systematic = mode == "Proj")

        # If we did not choose all the cell types
        if (!is.null(input$selected_cell_types) && length(input$selected_cell_types) != nrow(cell_type_colors())) {
            # TODO: add also genes that are distictive for the specific cell type (although not variable within it)
            # browser()
        }

        markers(new_markers)
    })

    observeEvent(input$remove_genes, {
        req(markers())
        new_markers <- markers()[!(markers() %in% input$selected_marker_genes)]
        markers(new_markers)
        shinyWidgets::updatePickerInput(session, ns("genes_to_add"), selected = c())
    })

    output$add_genes_ui <- renderUI({
        if (mode == "Inner") {
            mc_fp <- get_mc_data(dataset(), "inner_fold_mat")
            req(mc_fp)
            gene_choices <- rownames(mc_fp)[Matrix::rowSums(mc_fp) > 0]
            req(markers)
            gene_choices <- gene_choices[!(gene_choices %in% markers())]
        } else {
            gene_choices <- gene_names(dataset())
        }

        tagList(
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm"),
            shinyWidgets::pickerInput(ns("genes_to_add"),
                choices = gene_choices,
                selected = c(),
                multiple = TRUE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
            )
        )
    })

    observeEvent(input$add_genes, {
        new_markers <- sort(unique(c(markers(), input$genes_to_add)))
        markers(new_markers)
        shinyWidgets::updatePickerInput(session = session, inputId = "genes_to_add", selected = character(0))
    })
}


heatmap_reactives <- function(ns, input, output, session, dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, mode) {
    ns <- session$ns

    metacell_filter <- reactiveVal()

    mat <- reactive({
        req(markers())
        req(metacell_types())
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

        get_marker_matrix(
            dataset(),
            markers(),
            input$selected_cell_types,
            metacell_types(),
            gene_modules(),
            force_cell_type = input$force_cell_type,
            mode = mode,
            notify_var_genes = TRUE
        )
    }) %>% bindCache(dataset(), metacell_types(), cell_type_colors(), markers(), gene_modules(), input$selected_cell_types, input$force_cell_type, mode)

    heatmap_matrix_reactives(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, mode)

    output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors)

    output$metadata_list <- metadata_selector(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata")


    output$reset_zoom_ui <- renderUI({
        shinyWidgets::actionGroupButtons(ns("reset_zoom"), labels = "Reset zoom", size = "sm")
    })

    observe({
        shinyjs::toggle(id = "reset_zoom_ui", condition = !is.null(metacell_filter()) && length(metacell_filter()) > 0)
    })

    observeEvent(input$reset_zoom, {
        metacell_filter(NULL)
    })

    output$plotting_area <- renderUI({
        heatmap <- plotOutput(
            ns("heatmap"),
            height = "80vh",
            dblclick = dblclickOpts(ns("heatmap_dblclick"), clip = TRUE),
            hover = hoverOpts(ns("heatmap_hover"), delay = 500, delayType = "debounce"),
            brush = brushOpts(
                id = ns("heatmap_brush"),
                direction = "x",
                resetOnNew = TRUE
            )
        )


        if (!is.null(input$plot_legend) && input$plot_legend) {
            req(input$legend_width)
            legend_column <- column(
                width = input$legend_width,
                plotOutput(ns("markers_legend"), height = "80vh")
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

    output$markers_legend <- renderPlot({
        req(cell_type_colors())
        req(input$plot_legend)
        legend_point_size <- max(1, min(10, 250 / nrow(cell_type_colors())))
        legend <- cowplot::get_legend(cell_type_colors() %>%
            ggplot(aes(x = cell_type, color = cell_type, y = 1)) +
            geom_point() +
            scale_color_manual("", values = deframe(cell_type_colors() %>% select(cell_type, color))) +
            guides(color = guide_legend(override.aes = list(size = legend_point_size), ncol = 1)) +
            theme(legend.position = c(0.5, 0.5)))
        cowplot::ggdraw(legend)
    }) %>% bindCache(dataset(), cell_type_colors(), input$plot_legend)

    output$heatmap <- renderPlot({
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

        m <- filter_heatmap_by_metacell(mat(), metacell_filter())

        if (!is.null(input$selected_md)) {
            metadata <- get_mc_data(dataset(), "metadata") %>%
                select(metacell, one_of(input$selected_md))
        } else {
            metadata <- NULL
        }

        disjoined_genes <- get_mc_data(dataset(), "disjoined_genes_no_atlas")
        forbidden_genes <- get_mc_data(dataset(), "forbidden_genes")
        systematic_genes <- get_mc_data(dataset(), "systematic_genes")

        res <- plot_markers_mat(
            m,
            metacell_types(),
            cell_type_colors(),
            dataset(),
            min_lfp = lfp_range()[1],
            max_lfp = lfp_range()[2],
            plot_legend = FALSE,
            high_color =  input$high_color,
            low_color =  input$low_color,
            mid_color =  input$mid_color,
            midpoint = input$midpoint,
            metadata = metadata,
            forbidden_genes = forbidden_genes,
            systematic_genes = systematic_genes,
            disjoined_genes = disjoined_genes,
            col_names = ncol(m) <= 100,
            top_cell_type_bar = ncol(m) <= 100,
            interleave = nrow(m) > 80,
            vertial_gridlines = mode %in% c("Inner", "Proj"),
            separate_gtable = TRUE
        )

        # we are returning the gtable and ggplot object separatly in order to allow shiny to infer positions correctly.
        return(structure(list(p = res$p, gtable = res$gtable), class = "gt_custom"))
    }) %>% bindCache(dataset(), metacell_types(), cell_type_colors(), gene_modules(), lfp_range(), metacell_filter(), input$plot_legend, input$selected_md, markers(), input$selected_cell_types, input$force_cell_type, input$high_color, input$low_color, input$mid_color, input$midpoint)

    observeEvent(input$heatmap_brush, {
        m <- filter_heatmap_by_metacell(mat(), metacell_filter())
        range <- round(input$heatmap_brush$xmin):round(input$heatmap_brush$xmax)
        req(all(range > 0) && all(range <= ncol(m)))
        metacells <- colnames(m)[range]
        metacell_filter(metacells)
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
        req(mat())
        req(input$heatmap_hover)
        req(metacell_types())
        req(cell_type_colors())

        hover <- input$heatmap_hover
        m <- filter_heatmap_by_metacell(mat(), metacell_filter())
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

        if (mode == "Gene modules") {
            gene_prefix <- "Gene module:"
        } else {
            gene_prefix <- "Gene"
        }
        mcell_tooltip <- paste(
            glue("{gene_prefix}: {gene}"),
            glue("Metacell: {metacell}"),
            glue("Value: {round(value, digits=2)}"),
            glue("Cell type: {mcell_stats$cell_type}"),
            glue("Top genes: {top_genes}"),
            ifelse(has_name(mcell_stats, "mc_age"), glue("Metacell age (E[t]): {round(mcell_stats$mc_age, digits=2)}"), ""),
            sep = "<br/>"
        )

        wellPanel(
            style = style,
            p(HTML(mcell_tooltip))
        )
    })
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
        m <- m[, f]
        if (is.null(ncol(m))) { # a single metacell
            m <- as.matrix(m)
            colnames(m) <- f
        }
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
