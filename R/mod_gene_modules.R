#' gene_modules UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_gene_modules_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            resizable_column(
                width = 12,
                heatmap_box(ns, "Gene modules Heatmap")
            )
        )
    )
}


#' gene_modules sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_gene_modules_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            uiOutput(ns("reset_zoom_ui")),
            uiOutput(ns("cell_type_list")),
            uiOutput(ns("metadata_list")),
            uiOutput(ns("gene_module_list"))
        )
    )
}

#' gene_modules Server Function
#'
#' @noRd
mod_gene_modules_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    ns <- session$ns

    markers <- reactiveVal()
    lfp_range <- reactiveVal()

    output$gene_module_list <- gene_modules_selector(dataset, gene_modules, ns, id = "selected_gene_modules", selected = "all")

    observe({
        req(input$selected_gene_modules)
        markers(input$selected_gene_modules)
    })

    # mat <- reactive({}) # summarise the gene modules using tgs_matrix_tapply and then divide by the median for mc_fp.

    heatmap_reactives(ns, input, output, session, dataset, metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, "Gene modules")
}
