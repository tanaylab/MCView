#' stdev_fold UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_stdev_fold_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_column(
                width = 12,
                heatmap_box(ns("stdev_fold_heatmap"), "Stdev-fold Heatmap", fold_change_range = c(0, 1), low_color = "white", midpoint = 0.01, mid_color = "white", high_color = "red")
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
mod_stdev_fold_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        heatmap_sidebar(
            ns("stdev_fold_heatmap"),
            checkboxInput(NS(ns("stdev_fold_heatmap"))("use_markers"), "Use Markers", value = FALSE)
        )
    )
}

#' stdev_fold Server Function
#'
#' @noRd
mod_stdev_fold_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            markers <- reactiveVal()
            lfp_range <- reactiveVal()
            heatmap_reactives("stdev_fold_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, "Stdev")
        }
    )
}
