heatmap_box <- function(ns,
                        title = "Heatmap",
                        fold_change_range = c(-3, 3),
                        midpoint = 0,
                        low_color = "blue",
                        mid_color = "white",
                        high_color = "red") {
    div(
        shinydashboardPlus::box(
            id = ns("markers_heatmap_box"),
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
                id = ns("markers_heatmap_sidebar"),
                checkboxInput(ns("force_cell_type"), "Force cell type", value = TRUE),
                shinyWidgets::numericRangeInput(ns("lfp_range"), "Fold change range", fold_change_range, width = "80%", separator = " to "),
                numericInput(ns("midpoint"), "Midpoint", midpoint),
                colourpicker::colourInput(ns("low_color"), "Low color", low_color),
                colourpicker::colourInput(ns("high_color"), "High color", high_color),
                colourpicker::colourInput(ns("mid_color"), "Mid color", mid_color),
                checkboxInput(ns("plot_legend"), "Plot legend", value = TRUE),
                shinyWidgets::prettyRadioButtons(
                    inputId = ns("gene_select"),
                    label = "Select on click (gene)",
                    choices = c("X axis", "Y axis"),
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                ),
                shinyWidgets::prettyRadioButtons(
                    inputId = ns("metacell_select"),
                    label = "Select on click (metacell)",
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

heatmap_sidebar <- function(ns) {
    list(
        uiOutput(ns("cell_type_list")),
        uiOutput(ns("metadata_list")),
        shinyWidgets::actionGroupButtons(ns("update_markers"), labels = "Update genes", size = "sm"),
        numericInput(ns("max_gene_num"), "Maximal number of genes", value = 100),
        uiOutput(ns("add_genes_ui")),
        uiOutput(ns("marker_genes_list"))
    )
}

heatmap_reactives <- function(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, mode) {
    ns <- session$ns

    markers <- reactiveVal()
    lfp_range <- reactiveVal()

    output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors)

    output$metadata_list <- metadata_selector(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata")

    output$marker_genes_list <- renderUI({
        tagList(
            selectInput(
                ns("selected_marker_genes"),
                "Marker genes",
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

    markers_matrix <- reactive({
        req(markers())
        req(metacell_types())
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

        get_marker_matrix(
            dataset(),
            markers(),
            input$selected_cell_types,
            metacell_types(),
            force_cell_type = input$force_cell_type,
            mode = mode,
            notify_var_genes = TRUE
        )
    }) %>% bindCache(dataset(), metacell_types(), cell_type_colors(), markers(), input$selected_cell_types, input$force_cell_type, mode)


    observeEvent(input$update_markers, {
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
            browser()
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

    output$plotting_area <- renderUI({
        heatmap <- plotOutput(
            ns("markers_heatmap"),
            height = "80vh",
            click = clickOpts(ns("markers_heatmap_click")),
            hover = hoverOpts(ns("markers_heatmap_hover"), delay = 500, delayType = "debounce")
        )


        if (!is.null(input$plot_legend) && input$plot_legend) {
            legend_column <- column(
                width = 2,
                plotOutput(ns("markers_legend"))
            )
            heatmap_column <- column(
                width = 10,
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

    output$markers_heatmap <- renderPlot({
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

        mat <- markers_matrix()
        req(mat)

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
            mat,
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
            col_names = ncol(mat) <= 100,
            top_cell_type_bar = ncol(mat) <= 100,
            interleave = nrow(mat) > 80,
            vertial_gridlines = mode %in% c("Inner", "Proj"),
            separate_gtable = TRUE
        )

        # we are returning the gtable and ggplot object separatly in order to allow shiny to infer positions correctly.
        return(structure(list(p = res$p, gtable = res$gtable), class = "gt_custom"))
    }) %>% bindCache(dataset(), metacell_types(), cell_type_colors(), lfp_range(), input$plot_legend, input$selected_md, markers(), input$selected_cell_types, input$force_cell_type, input$high_color, input$low_color, input$mid_color, input$midpoint)

    observeEvent(input$markers_heatmap_click, {
        gene <- rownames(markers_matrix())[round(input$markers_heatmap_click$y)]
        metacell <- colnames(markers_matrix())[round(input$markers_heatmap_click$x)]

        if (input$gene_select == "X axis") {
            globals$selected_gene_x_axis <- gene
        } else {
            globals$selected_gene_y_axis <- gene
        }
        globals$selected_query_gene <- gene

        if (input$metacell_select == "Metacell A") {
            globals$selected_metacellA <- metacell
        } else {
            globals$selected_metacellB <- metacell
        }

        globals$selected_query_metacell <- metacell
    })

    output$hover_info <- renderUI({
        req(markers_matrix())
        req(input$markers_heatmap_hover)
        req(metacell_types())
        req(cell_type_colors())

        hover <- input$markers_heatmap_hover

        gene <- rownames(markers_matrix())[round(hover$y)]
        metacell <- colnames(markers_matrix())[round(hover$x)]
        value <- markers_matrix()[gene, metacell]

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

        mcell_tooltip <- paste(
            glue("Gene: {gene}"),
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

# we are overriding the print function, similiar to the one at shiny/render-plot.R:
# https://github.com/rstudio/shiny/blob/main/R/render-plot.R
# in order to provide a separate ggplot_build and gtable objects
print.gt_custom <<- function(x) {
    build <- ggplot_build(x$p)

    grid::grid.newpage()
    grid::grid.draw(x$gtable)

    structure(list(
        build = build,
        gtable = x$gtable
    ), class = "ggplot_build_gtable")
}
