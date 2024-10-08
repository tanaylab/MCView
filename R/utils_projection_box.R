projection_box <- function(ns,
                           id,
                           ...,
                           color_choices = c("Cell type", "Gene", "Gene module", "Metadata"),
                           title = "2D Projection",
                           height = NULL,
                           plotly_height = "400px",
                           additional_elements = NULL,
                           collapsed_accordion = TRUE,
                           legend_orientation = "Vertical",
                           show_legend = TRUE) {
    generic_box(
        id = ns(id),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        height = height,
        sidebar = shinydashboardPlus::boxSidebar(
            startOpen = FALSE,
            width = 80,
            id = ns(glue("{id}_sidebar")),
            uiOutput(ns("graph_select_ui")),
            uiOutput(ns("proj_stat_ui")),
            uiOutput(ns("set_range_ui")),
            uiOutput(ns("expr_range_ui")),
            uiOutput(ns("enrich_range_ui")),
            uiOutput(ns("point_size_ui")),
            uiOutput(ns("stroke_ui")),
            uiOutput(ns("edge_distance_ui")),
            shinyWidgets::prettyRadioButtons(
                ns("legend_orientation"),
                label = "Legend orientation:",
                choices = c("Vertical", "Horizontal"),
                selected = legend_orientation,
                inline = TRUE,
                status = "danger",
                fill = TRUE
            )
        ),
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("plot_gene_proj_2d"), height = plotly_height)
        ),
        shinydashboardPlus::accordion(
            id = ns("proj_accordion"),
            shinydashboardPlus::accordionItem(
                title = "Modify colors",
                collapsed = collapsed_accordion,
                shinyWidgets::prettyRadioButtons(
                    ns("color_proj"),
                    label = "Color by:",
                    choices = color_choices,
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                ),
                ...,
                shinyWidgets::virtualSelectInput(
                    ns("color_proj_gene"),
                    "Gene:",
                    choices = c(),
                    multiple = FALSE,
                    search = TRUE,
                    dropboxWrapper = "body",
                    markSearchResults = TRUE,
                    searchByStartsWith = TRUE
                ),
                shinyWidgets::virtualSelectInput(
                    ns("color_proj_metadata"),
                    "Metadata:",
                    choices = c(),
                    multiple = FALSE,
                    search = TRUE,
                    dropboxWrapper = "body"
                ),
                shinyWidgets::virtualSelectInput(
                    ns("color_proj_gene_module"),
                    "Gene module:",
                    choices = c(),
                    multiple = FALSE,
                    search = TRUE,
                    dropboxWrapper = "body"
                ),
                shinyWidgets::prettyRadioButtons(
                    ns("scatter_axis_proj"),
                    label = "Axis:",
                    choices = c("x", "y"),
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                ),
                checkboxInput(ns("show_legend_projection"), "Show legend", value = show_legend)
            )
        )
    )
}

projection_selectors <- function(ns, dataset, output, input, gene_modules, globals, session, weight = 1, atlas = FALSE) {
    observe({
        shinyWidgets::updateVirtualSelect(
            session = session,
            inputId = "color_proj_gene",
            choices = gene_names_label(dataset(), atlas = atlas),
            selected = default_gene1
        )

        shinyWidgets::updateVirtualSelect(
            session = session,
            inputId = "color_proj_metadata",
            choices = c("Clipboard", dataset_metadata_fields(dataset(), atlas = atlas)),
            selected = dataset_metadata_fields(dataset(), atlas = atlas)[1]
        )

        shinyWidgets::updateVirtualSelect(
            session = session,
            inputId = "color_proj_gene_module",
            choices = levels(gene_modules()$module),
            selected = NULL
        )
    })

    picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE)

    observe({
        req(input$color_proj)
        shinyjs::toggle(id = "color_proj_gene", condition = input$color_proj == "Gene")
        shinyjs::toggle(id = "color_proj_metadata", condition = input$color_proj == "Metadata")
        shinyjs::toggle(id = "color_proj_gene_module", condition = input$color_proj == "Gene module")
        shinyjs::toggle(id = "scatter_axis_proj", condition = input$color_proj == "Scatter Axis")
    })


    output$proj_stat_ui <- renderUI({
        req(input$color_proj == "Gene" || input$color_proj == "Gene A" || input$color_proj == "Gene B" || input$color_proj == "Gene module")
        selectInput(ns("proj_stat"), label = "Statistic", choices = c("Expression" = "expression", "Enrichment" = "enrichment"), selected = "Expression", multiple = FALSE, selectize = FALSE)
    })

    output$graph_select_ui <- renderUI({
        choices <- c("metacell")
        graphs <- get_mc_data(dataset(), "metacell_graphs")
        if (!is.null(graphs)) {
            choices <- c(choices, names(graphs))
        }
        selectInput(ns("graph_name"), label = "Graph", choices = choices, selected = "metacell", multiple = FALSE, selectize = FALSE)
    })

    # Expression range
    output$set_range_ui <- renderUI({
        req(input$color_proj == "Gene" || input$color_proj == "Gene A" || input$color_proj == "Gene B" || input$color_proj == "Gene module")
        req(input$proj_stat == "expression")
        checkboxInput(ns("set_range"), "Manual range", value = FALSE)
    })

    output$expr_range_ui <- renderUI({
        req(input$color_proj == "Gene" || input$color_proj == "Gene A" || input$color_proj == "Gene B" || input$color_proj == "Gene module")
        req(input$proj_stat == "expression")
        req(input$set_range)
        shinyWidgets::numericRangeInput(ns("expr_range"), "Expression range", c(-18, -5), width = "80%", separator = " to ")
    })

    # Enrichment range
    output$enrich_range_ui <- renderUI({
        req(input$color_proj == "Gene" || input$color_proj == "Gene A" || input$color_proj == "Gene B" || input$color_proj == "Gene module")
        req(input$proj_stat == "enrichment")
        shinyWidgets::numericRangeInput(ns("lfp"), "Enrichment range", c(-3, 3), width = "80%", separator = " to ")
    })

    # Point size selectors
    output$point_size_ui <- renderUI({
        req(globals$screen_height)
        req(globals$screen_width)
        req(dataset())
        numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset(), globals$screen_width, globals$screen_height, weight = weight, atlas = atlas), min = 0.1, max = 3, step = 0.1)
    })

    output$stroke_ui <- renderUI({
        numericInput(ns("stroke"), label = "Stroke width", value = initial_proj_stroke(dataset()), min = 0, max = 3, step = 0.01)
    })

    # Minimal edge length selector
    output$edge_distance_ui <- renderUI({
        graph <- input$graph_name
        if (is.null(graph) || graph == "metacell") {
            sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
        } else {
            graph <- get_mc_data(dataset(), "metacell_graphs")[[graph]]
            sliderInput(ns("min_edge_size"), label = "Min weight", min = min(graph$weight, na.rm = TRUE), max = max(graph$weight, na.rm = TRUE), value = median(graph$weight), step = (max(graph$weight, na.rm = TRUE) - min(graph$weight, na.rm = TRUE)) / 50)
        }
    })
}
