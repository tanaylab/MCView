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
        includeRMarkdown(about_file)
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
mod_about_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
        }
    )
}
