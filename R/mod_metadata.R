#' gene_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_metadata_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 12,
                shinydashboardPlus::box(
                    id = ns("metadata_projection"),
                    title = "Metadata",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "80vh",
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("metadata_projection_sidebar"),
                        uiOutput(ns("metadata_range_ui")),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_metadata_proj_2d"), height = "80vh")
                    ),
                    uiOutput(ns("manifold_metadata_select_ui"))
                )
            )
        )
    )
}


#' metadata sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_metadata_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            # uiOutput(ns("gene_selectors")),
            # tags$hr(),
            # uiOutput(ns("top_correlated_select_gene1")),
            # uiOutput(ns("top_correlated_select_gene2")),
            # tags$hr(),
            # uiOutput(ns("genecards_buttons"))
        )
    )
}

#' gene_mc Server Function
#'
#' @noRd
mod_metadata_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    # Manifold metadata select
    output$manifold_metadata_select_ui <- renderUI({
        if (!has_metadata(dataset())) {
            print(glue("Dataset \"{dataset()}\" doesn't have any metadata. Use `update_metadata` to add it to your dataset."))
        } else {
            selectInput(ns("manifold_metadata"), label = "Metadata", choices = dataset_metadata_fields(dataset()), selected = dataset_metadata_fields(dataset())[1])
        }
    })

    output$metadata_range_ui <- renderUI({
        shinyWidgets::numericRangeInput(ns("expr_range"), "Expression range", c(-18, -5), width = "80%", separator = " to ")
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
    output$plot_metadata_proj_2d <- render_2d_plotly(input, output, session, dataset, values, metacell_types, cell_type_colors, source = "proj_metadata_plot")
}
