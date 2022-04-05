#' inner_fold UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_inner_fold_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            resizable_column(
                width = 12,
                heatmap_box(ns("inner_fold_heatmap"), "Inner-Fold Heatmap")
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
mod_inner_fold_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        heatmap_sidebar(
            ns("inner_fold_heatmap"),
            checkboxInput(NS(ns("inner_fold_heatmap"))("use_markers"), "Use Markers", value = FALSE)
        )
    )
}

#' inner_fold Server Function
#'
#' @noRd
mod_inner_fold_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            markers <- reactiveVal()
            lfp_range <- reactiveVal()
            heatmap_reactives("inner_fold_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, "Inner")
        }
    )
}
