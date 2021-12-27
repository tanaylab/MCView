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
                    title = "2D Projection",
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
                        uiOutput(ns("stroke_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinyWidgets::prettyRadioButtons(
                        ns("color_proj"),
                        label = "Color by:",
                        choices = c("Cell type", "Gene A", "Gene B"),
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_manifold_proj_2d"), height = "80vh")
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
            shinycssloaders::withSpinner(uiOutput(ns("gene_selectors"))),
            shinycssloaders::withSpinner(uiOutput(ns("top_correlated_select_gene1"))),
            shinycssloaders::withSpinner(uiOutput(ns("top_correlated_select_gene2"))),
            shinycssloaders::withSpinner(uiOutput(ns("genecards_buttons")))
        )
    )
}

#' gene_mc Server Function
#'
#' @noRd
mod_manifold_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    # gene selectors
    server_gene_selectors(input, output, session, dataset, ns)

    projection_selectors(ns, dataset, output, input)

    # Projection plots
    output$plot_manifold_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, source = "proj_manifold_plot") %>%
        bindCache(dataset(), input$gene1, input$gene2, input$color_proj, metacell_types(), cell_type_colors(), input$point_size, input$stroke, input$min_edge_size, input$show_selected_metacells, input$metacell1, input$metacell2, input$proj_stat, input$expr_range, input$lfp, input$set_range)
}
