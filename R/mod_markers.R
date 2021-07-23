#' gene_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_markers_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 12,
                shinydashboardPlus::box(
                    id = ns("markers_heatmap_box"),
                    title = "Markers heatmap",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "80vh",
                    # sidebar = shinydashboardPlus::boxSidebar(
                    #     startOpen = FALSE,
                    #     width = 25,
                    #     id = ns("gene_projection_sidebar"),
                    #     selectInput(ns("proj_stat"), label = "Statistic", choices = c("Expression" = "expression", "Enrichment" = "enrichment"), selected = "Expression", multiple = FALSE, selectize = FALSE),
                    #     uiOutput(ns("set_range_ui")),
                    #     uiOutput(ns("expr_range_ui")),
                    #     uiOutput(ns("enrich_range_ui")),
                    #     uiOutput(ns("point_size_ui")),
                    #     uiOutput(ns("edge_distance_ui"))
                    # ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("markers_heatmap"), height = "80vh")
                    ) # ,
                    # shinyWidgets::prettyRadioButtons(
                    #     ns("color_proj"),
                    #     label = "Color by:",
                    #     choices = c("Cell type", "Gene A", "Gene B"),
                    #     inline = TRUE,
                    #     status = "danger",
                    #     fill = TRUE
                    # )
                )
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
mod_markers_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            uiOutput(ns("marker_genes_list")),
            uiOutput(ns("add_genes_ui"))
        )
    )
}

#' gene_mc Server Function
#'
#' @noRd
mod_markers_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    markers <- reactiveVal()

    observe({
        initial_markers <- sort(get_mc_data(dataset(), "marker_genes"))

        markers(initial_markers)
    })

    output$marker_genes_list <- renderUI({
        tagList(
            selectInput(
                ns("selected_marker_genes"),
                "Marker genes",
                choices = markers(),
                selected = NULL,
                multiple = TRUE,
                size = 30,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(ns("remove_genes"), labels = "Remove genes", size = "sm")
        )
    })

    observeEvent(input$remove_genes, {
        new_markers <- markers()[!(markers() %in% input$selected_marker_genes)]
        markers(new_markers)
        shinyWidgets::updatePickerInput(session, ns("genes_to_add"), selected = c())
    })

    output$add_genes_ui <- renderUI({
        tagList(
            shinyWidgets::pickerInput(ns("genes_to_add"), "",
                choices = gene_names,
                selected = c(),
                multiple = TRUE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
            ),
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm")
        )
    })

    observeEvent(input$add_genes, {
        new_markers <- sort(unique(c(markers(), input$genes_to_add)))
        markers(new_markers)
        shinyWidgets::updatePickerInput(session = session, inputId = "genes_to_add", selected = character(0))
    })

    output$markers_heatmap <- plotly::renderPlotly({
        req(dataset())
        req(markers())

        plot_markers_mat(
            get_mc_fp(dataset(), markers()),
            metacell_types(),
            cell_type_colors()
        )
    })
}
