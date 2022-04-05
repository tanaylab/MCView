#' outliers UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_outliers_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            resizable_column(
                width = 12,
                heatmap_box(ns("outliers_heatmap"), "Outliers Heatmap")
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
mod_outliers_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        heatmap_sidebar(ns("outliers_heatmap"))
    )
}

#' outliers Server Function
#'
#' @noRd
mod_outliers_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            markers <- reactiveVal()
            lfp_range <- reactiveVal()

            # use the 'most_similar' field in the outliers_metadata table to get the cell type for each outlier
            outliers_metacell_types <- reactive({
                outliers_metadata <- get_mc_data(dataset(), "outliers_metadata")
                req(outliers_metadata)
                req(outliers_metadata$most_similar)
                outliers_metadata %>%
                    select(cell_id, metacell = most_similar) %>%
                    mutate(metacell = as.character(metacell)) %>%
                    left_join(metacell_types(), by = "metacell") %>%
                    mutate(most_similar_metacell = metacell) %>%
                    select(-metacell) %>%
                    rename(metacell = cell_id)
            })

            heatmap_reactives("outliers_heatmap", dataset, outliers_metacell_types, gene_modules, cell_type_colors, globals, markers, lfp_range, "Outliers")
        }
    )
}
