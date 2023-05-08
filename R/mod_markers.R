#' markers UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_markers_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_column(
                width = 12,
                heatmap_box(ns("markers_heatmap"), "Markers Heatmap")
            )
        )
    )
}


#' markers sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_markers_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        heatmap_sidebar(ns("markers_heatmap"))
    )
}

#' markers Server Function
#'
#' @noRd
mod_markers_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            markers <- reactiveVal()
            lfp_range <- reactiveVal()
            heatmap_reactives("markers_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, "Markers")
        }
    )
}
