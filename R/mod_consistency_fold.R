#' consistency_fold UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_consistency_fold_ui <- function(id) {
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
mod_consistency_fold_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        heatmap_sidebar(ns)
    )
}

#' consistency_fold Server Function
#'
#' @noRd
mod_consistency_fold_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals) {
    ns <- session$ns

    markers <- reactiveVal()
    lfp_range <- reactiveVal()

    output$hm_box <- renderUI({
        max_fold <- get_mc_data(dataset(), "project_max_consistency_fold_factor") %||% 4
        heatmap_box(
            ns,
            "Consistency-Fold Heatmap",
            fold_change_range = c(max_fold - 1, max_fold + 1),
            midpoint = max_fold,
            low_color = "white",
            mid_color = "red",
            high_color = "black"
        )
    })

    heatmap_reactives(ns, input, output, session, dataset, metacell_types, cell_type_colors, globals, markers, lfp_range, "Consistency")
}
