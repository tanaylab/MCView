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
            resizable_column(
                width = 12,
                uiOutput(ns("hm_box"))
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
        heatmap_sidebar(ns("proj_fold_heatmap"))
    )
}

#' proj_fold Server Function
#'
#' @noRd
mod_proj_fold_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            markers <- reactiveVal()
            lfp_range <- reactiveVal()

            output$hm_box <- renderUI({
                max_fold <- get_mc_data(dataset(), "project_max_projection_fold_factor") %||% 3
                heatmap_box(
                    ns("proj_fold_heatmap"),
                    "Projected-Fold Heatmap",
                    fold_change_range = c(-max_fold - 1, max_fold + 1),
                    midpoint = 0,
                    low_color = "blue",
                    mid_color = "white",
                    high_color = "red"
                )
            })
            heatmap_reactives("proj_fold_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, "Proj")
        }
    )



    # heatmap_reactives(ns("proj_fold"), input, output, session, dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, "Proj")
}
