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
                width = 9,
                heatmap_box(ns("gene_modules_heatmap"), "Gene modules Heatmap", legend_width = 2)
            )
        ),
        fluidRow(
            resizable_column(
                width = 9,
                heatmap_box(ns("genes_heatmap"), "Genes Heatmap", legend_width = 2)
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
    ns_heatmap <- NS(ns("gene_modules_heatmap"))
    ns_genes <- NS(ns("genes_heatmap"))
    tagList(
        list(
            HTML("<h5><b><center>Genes modules heatmap</center></b></h5>"),
            uiOutput(ns_heatmap("reset_zoom_ui")),
            uiOutput(ns("shown_gene_modules_ui")),
            uiOutput(ns_heatmap("cell_type_list")),
            uiOutput(ns_heatmap("metadata_list")),
            tags$hr(),
            HTML("<h5><b><center>Genes heatmap</center></b></h5>"),
            uiOutput(ns_genes("reset_zoom_ui")),
            uiOutput(ns("selected_gene_modules_ui")),
            uiOutput(ns_genes("cell_type_list")),
            uiOutput(ns_genes("metadata_list"))
        )
    )
}

#' gene_modules Server Function
#'
#' @noRd
mod_gene_modules_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            shown_gene_modules <- reactiveVal()
            genes <- reactiveVal()
            lfp_range <- reactiveVal()

            output$shown_gene_modules_ui <- gene_modules_selector(
                dataset,
                gene_modules,
                ns,
                label = "Shown gene modules",
                id = "shown_gene_modules",
                selected = "all"
            )

            observe({
                shown_gene_modules(input$shown_gene_modules)
            })

            heatmap_reactives("gene_modules_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, shown_gene_modules, lfp_range, "Gene modules")

            # Genes of a specific gene module
            output$selected_gene_modules_ui <- renderUI({
                modules <- unique(gene_modules()$module)

                shinyWidgets::pickerInput(
                    ns("selected_gene_modules"),
                    "Selected gene modules",
                    choices = modules,
                    selected = NULL,
                    multiple = TRUE,
                    options = shinyWidgets::pickerOptions(
                        liveSearch = TRUE,
                        liveSearchNormalize = TRUE,
                        liveSearchStyle = "startsWith",
                        `dropup-auto` = FALSE,
                        `max-options` = 3,
                        `max-options-text` = "Cannot choose more than 3 gene modules"
                    )
                )
            })

            observe({
                if (is.null(input$selected_gene_modules) || length(input$selected_gene_modules) == 0) {
                    genes(character(0))
                } else {
                    req(gene_modules())
                    new_genes <- gene_modules() %>%
                        filter(module %in% input$selected_gene_modules) %>%
                        pull(gene)
                    genes(new_genes)
                }
            })

            heatmap_reactives("genes_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, genes, lfp_range, "Markers")
        }
    )
}
