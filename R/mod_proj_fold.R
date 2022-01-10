#' proj_fold UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_proj_fold_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 12,
                heatmap_box(ns, "Projected-Fold Heatmap")
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
mod_proj_fold_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        heatmap_sidebar(ns)
    )
}

#' proj_fold Server Function
#'
#' @noRd
mod_proj_fold_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals) {
    ns <- session$ns

    markers <- reactiveVal()
    lfp_range <- reactiveVal()

    heatmap_reactives(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, "Proj")
}
