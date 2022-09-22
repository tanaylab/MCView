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
    tagList(
        fluidRow(
            resizable_column(
                width = 12,
                shinydashboardPlus::box(
                    id = ns("metacell_projection"),
                    title = "Atlas 2D Projection",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "80vh",
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 80,
                        id = ns("gene_projection_sidebar"),
                        uiOutput(ns("proj_stat_ui")),
                        uiOutput(ns("set_range_ui")),
                        uiOutput(ns("expr_range_ui")),
                        uiOutput(ns("enrich_range_ui")),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("stroke_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_proj_2d"), height = "80vh")
                    )
                )
            )
        )
    )
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
        list(
            shinyWidgets::prettyRadioButtons(
                ns("color_proj"),
                label = "Color by:",
                choices = c("Cell type", "Gene", "Gene module", "Metadata", "Query metacell", "Query cell type"),
                inline = FALSE,
                status = "danger",
                fill = TRUE
            ),
            uiOutput(ns("gene_selector")),
            uiOutput(ns("metadata_selector")),
            uiOutput(ns("metacell_selector")),
            uiOutput(ns("cell_type_selector")),
            uiOutput(ns("proj_gene_module_selector")),
            shinyWidgets::prettyRadioButtons(
                ns("color_by_scale"),
                label = "Color scale:",
                choices = c("Cell type", "Discrete", "Continuous"),
                inline = FALSE,
                status = "danger",
                fill = TRUE
            ),
            sliderInput(ns("query_threshold"), "Threshold", 0, 1, 0.1),
            uiOutput(ns("top_correlated_select_color_proj"))
        )
    )
}

#' atlas Server Function
#'
#' @noRd
mod_atlas_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            atlas_colors <- reactive({
                req(has_atlas(dataset()))
                get_mc_data(dataset(), "cell_type_colors", atlas = TRUE)
            })

            atlas_metacell_types <- reactive({
                req(has_atlas(dataset()))
                get_mc_data(dataset(), "metacell_types", atlas = TRUE)
            })

            atlas_gene_modules <- gene_modules <- reactive({
                mods <- get_mc_data(dataset(), "gene_modules", atlas = TRUE)
                if (is.null(mods)) {
                    mods <- tibble(gene = character(0), module = factor())
                }
                return(mods)
            })

            projection_selectors(ns, dataset, output, input, atlas_gene_modules, globals, weight = 1, atlas = TRUE)
            output$top_correlated_select_color_proj <- renderUI({
                req(input$gene1)
                req(has_gg_mc_top_cor(project, dataset()))
                tagList(
                    selectInput(
                        ns("selected_top_gene"),
                        glue("Top correlated to {input$color_proj_gene}:"),
                        choices = c(get_top_cor_gene(dataset(), input$color_proj_gene, type = "pos"), rev(get_top_cor_gene(dataset(), input$color_proj_gene, type = "neg"))),
                        selected = NULL,
                        size = 10,
                        selectize = FALSE
                    )
                )
            })

            output$cell_type_selector <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", cell_type_colors = cell_type_colors, selected = "all")

            output$metacell_selector <- metacell_selector(dataset, ns, id = "selected_metacell", label = "Metacell")

            observe({
                req(input$color_proj)
                shinyjs::toggle(id = "cell_type_selector", condition = input$color_proj == "Query cell type")
                shinyjs::toggle(id = "metacell_selector", condition = input$color_proj == "Query metacell")
                shinyjs::toggle(id = "color_by_scale", condition = input$color_proj %in% c("Query metacell", "Query cell type"))
                shinyjs::toggle(id = "query_threshold", condition = input$color_proj %in% c("Query metacell", "Query cell type"))
            })


            output$plot_mc_proj_2d <- render_2d_plotly(input, output, session, dataset, atlas_metacell_types, atlas_colors, globals = globals, atlas = TRUE, query_types = metacell_types, atlas_gene_modules, source = "proj_mc_plot_proj_tab") %>% bindCache(dataset(), input$color_proj, atlas_metacell_types(), atlas_colors(), input$point_size, input$stroke, input$min_edge_size, input$set_range, input$proj_stat, input$expr_range, input$lfp, input$color_proj_gene, input$color_proj_metadata, input$selected_cell_types, input$selected_metacell, input$color_by_scale, input$query_threshold, input$color_proj_gene_module, input$graph_name)
        }
    )
}
