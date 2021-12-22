#' atlas UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_atlas_ui <- function(id) {
    ns <- NS(id)
    tagList()
}


#' atlas sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_atlas_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list()
    )
}

#' atlas Server Function
#'
#' @noRd
mod_atlas_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns
}
