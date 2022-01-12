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
                choices = c("Cell type", "Gene", "Metadata", "Query metacell", "Query cell type"),
                inline = FALSE,
                status = "danger",
                fill = TRUE
            ),
            uiOutput(ns("gene_selector")),
            uiOutput(ns("metadata_selector")),
            uiOutput(ns("metacell_selector")),
            uiOutput(ns("cell_type_selector")),
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
mod_atlas_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals) {
    ns <- session$ns

    projection_selectors(ns, dataset, output, input, globals, weight = 1, atlas = TRUE)
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

    picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE)

    output$gene_selector <- renderUI({
        shinyWidgets::pickerInput(
            ns("color_proj_gene"),
            label = "Gene:",
            choices = gene_names(dataset(), atlas = TRUE),
            selected = default_gene1,
            multiple = FALSE,
            options = picker_options
        )
    })

    output$metadata_selector <- renderUI({
        if (!has_metadata(dataset(), atlas = TRUE)) {
            HTML(glue("&nbspAtlas doesn't have any metadata."))
        } else {
            shinyWidgets::pickerInput(
                ns("color_proj_metadata"),
                label = "Metadata:",
                choices = dataset_metadata_fields(dataset(), atlas = TRUE),
                selected = dataset_metadata_fields(dataset(), atlas = TRUE)[1],
                multiple = FALSE,
                options = picker_options
            )
        }
    })

    output$cell_type_selector <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", cell_type_colors = cell_type_colors, selected = "all")

    output$metacell_selector <- metacell_selector(dataset, ns, id = "selected_metacell", label = "Metacell")


    observe({
        req(input$color_proj)
        shinyjs::toggle(id = "gene_selector", condition = input$color_proj == "Gene")
        shinyjs::toggle(id = "metadata_selector", condition = input$color_proj == "Metadata")
        shinyjs::toggle(id = "cell_type_selector", condition = input$color_proj == "Query cell type")
        shinyjs::toggle(id = "metacell_selector", condition = input$color_proj == "Query metacell")
        shinyjs::toggle(id = "color_by_scale", condition = input$color_proj %in% c("Query metacell", "Query cell type"))
        shinyjs::toggle(id = "query_threshold", condition = input$color_proj %in% c("Query metacell", "Query cell type"))
    })

    atlas_colors <- reactive({
        req(has_atlas(dataset()))
        get_mc_data(dataset(), "cell_type_colors", atlas = TRUE)
    })

    atlas_metacell_types <- reactive({
        req(has_atlas(dataset()))
        get_mc_data(dataset(), "metacell_types", atlas = TRUE)
    })

    output$plot_mc_proj_2d <- render_2d_plotly(input, output, session, dataset, atlas_metacell_types, atlas_colors, atlas = TRUE, query_types = metacell_types, source = "proj_mc_plot_proj_tab") %>% bindCache(dataset(), input$color_proj, atlas_metacell_types(), atlas_colors(), input$point_size, input$stroke, input$min_edge_size, input$set_range, input$proj_stat, input$expr_range, input$lfp, input$color_proj_gene, input$color_proj_metadata, input$selected_cell_types, input$selected_metacell, input$color_by_scale, input$query_threshold)
}
