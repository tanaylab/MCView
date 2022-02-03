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
            resizable_column(
                width = 12,
                heatmap_box(ns, "Markers Heatmap")
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
        heatmap_sidebar(ns)
    )
}

#' markers Server Function
#'
#' @noRd
mod_markers_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    ns <- session$ns

    markers <- reactiveVal()
    lfp_range <- reactiveVal()

    heatmap_reactives(ns, input, output, session, dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, "Markers")
}
