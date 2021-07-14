#' gene_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_manifold_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 12,
                shinydashboardPlus::box(
                    id = ns("gene_projection"),
                    title = "Manifold",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "80vh",
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("gene_projection_sidebar"),
                        selectInput(ns("proj_stat"), label = "Statistic", choices = c("Expression" = "expression", "Enrichment" = "enrichment"), selected = "Expression", multiple = FALSE, selectize = FALSE),
                        uiOutput(ns("set_range_ui")),
                        uiOutput(ns("expr_range_ui")),
                        uiOutput(ns("enrich_range_ui")),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_manifold_proj_2d"), height = "80vh")
                    ),
                    shinyWidgets::prettyRadioButtons(
                        ns("color_proj"),
                        label = "Color by:",
                        choices = c("Cell type", "Gene A", "Gene B"),
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    )
                )
            )
        )
    )
}


#' manifold sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_manifold_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            uiOutput(ns("gene_selectors")),
            tags$hr(),
            uiOutput(ns("top_correlated_select_gene1")),
            uiOutput(ns("top_correlated_select_gene2")),
            tags$hr(),
            uiOutput(ns("genecards_buttons"))
        )
    )
}

#' gene_mc Server Function
#'
#' @noRd
mod_manifold_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    # gene selectors
    values <- reactiveValues(gene1 = default_gene1, gene2 = default_gene2)
    server_gene_selectors(input, output, session, values, dataset, ns)

    # Expression range
    output$set_range_ui <- renderUI({
        req(input$proj_stat == "expression")
        checkboxInput(ns("set_range"), "Manual range", value = FALSE)
    })

    output$expr_range_ui <- renderUI({
        req(input$proj_stat == "expression")
        req(input$set_range)
        shinyWidgets::numericRangeInput(ns("expr_range"), "Expression range", c(-18, -5), width = "80%", separator = " to ")
    })

    # Enrichment range
    output$enrich_range_ui <- renderUI({
        req(input$proj_stat == "enrichment")
        shinyWidgets::numericRangeInput(ns("lfp"), "Enrichment range", c(-3, 3), width = "80%", separator = " to ")
    })

    # Point size selector
    output$point_size_ui <- renderUI({
        numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset()), min = 0.1, max = 3, step = 0.1)
    })

    # Minimal edge length selector
    output$edge_distance_ui <- renderUI({
        sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
    })

    # Projection plots
    output$plot_manifold_proj_2d <- render_2d_plotly(input, output, session, dataset, values, metacell_types, cell_type_colors, source = "proj_manifold_plot")
}
