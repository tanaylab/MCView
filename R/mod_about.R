#' about UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_about_ui <- function(id) {
    ns <- NS(id)
    tagList(
        includeRMarkdown(about_file),
        plotOutput(ns("plot_manifold_proj_2d"))
    )
}

#' about sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_about_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list()
    )
}

#' about Server Function
#'
#' @noRd
mod_about_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    ns <- session$ns

    output$plot_manifold_proj_2d <- renderPlot({
        req(metacell_types())
        req(cell_type_colors())
        mc2d_plot_ggp(
            dataset(),
            metacell_types = metacell_types(),
            cell_type_colors = cell_type_colors(),
            point_size = initial_proj_point_size(dataset(), globals$screen_width, globals$screen_height),
            stroke = initial_proj_stroke(dataset()),
            min_d = min_edge_length(dataset())
        ) + theme(aspect.ratio = 1)
    })
}
