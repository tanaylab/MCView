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
                    title = "Markers Heatmap",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "80vh",
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("markers_heatmap_sidebar"),
                        checkboxInput(ns("force_cell_type"), "Force cell type", value = FALSE),
                        shinyWidgets::numericRangeInput(ns("lfp_range"), "Fold change range", c(-3, 3), width = "80%", separator = " to "),
                        checkboxInput(ns("plot_legend"), "Plot legend", value = FALSE)
                    ),
                    shinycssloaders::withSpinner(
                        plotOutput(ns("markers_heatmap"), height = "80vh")
                    )
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
            shinyWidgets::actionGroupButtons(ns("apply"), labels = "Draw Heatmap"),
            uiOutput(ns("cell_type_list")),
            shinyWidgets::actionGroupButtons(ns("update_markers"), labels = "Update markers", size = "sm"),
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
    markers_matrix <- reactiveVal()
    lfp_range <- reactiveVal()

    observe({
        initial_markers <- choose_markers(get_markers(dataset()), 80)
        markers(initial_markers)

        mat <- get_marker_matrix(
            dataset(),
            initial_markers
        )
        markers_matrix(mat)

        lfp_range(c(-3, 3))
    })

    observeEvent(input$update_markers, {
        req(metacell_types())
        req(input$selected_cell_types)

        markers_df <- metacell_types() %>%
            filter(cell_type %in% input$selected_cell_types) %>%
            select(metacell) %>%
            inner_join(get_markers(dataset()), by = "metacell")
        new_markers <- choose_markers(markers_df, 100)

        markers(new_markers)
    })

    output$cell_type_list <- renderUI({
        cell_types_hex <- col2hex(cell_type_colors()$color)
        cell_types <- cell_type_colors()$cell_type
        shinyWidgets::pickerInput(ns("selected_cell_types"), "Cell types",
            choices = cell_types,
            selected = cell_types,
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `dropup-auto` = FALSE),
            choicesOpt = list(
                style = paste0("color: ", cell_types_hex, ";")
            )
        )
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
            shinyWidgets::actionGroupButtons(ns("remove_genes"), labels = "Remove selected genes", size = "sm")
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

    observeEvent(input$apply, {
        shinyjs::hide("markers_heatmap")
        selected_cell_types <- input$selected_cell_types %||% cell_type_colors()$cell_type
        force_cell_type <- input$force_cell_type %||% FALSE
        req(markers())

        mat <- get_marker_matrix(
            dataset(),
            markers(),
            input$selected_cell_types,
            metacell_types(),
            force_cell_type = input$force_cell_type
        )

        markers_matrix(mat)

        lfp_range(input$lfp_range)
        shinyjs::show("markers_heatmap")
    })

    output$markers_heatmap <- renderPlot({
        req(dataset())
        req(markers_matrix())

        plot_markers_mat(
            markers_matrix(),
            metacell_types(),
            cell_type_colors(),
            min_lfp = lfp_range()[1],
            max_lfp = lfp_range()[2],
            plot_legend = input$plot_legend %||% TRUE
        )
    })
}

get_marker_matrix <- function(dataset, markers, cell_types = NULL, metacell_types = NULL, force_cell_type = FALSE) {
    mc_fp <- get_mc_fp(dataset, markers)

    if (!is.null(cell_types)) {
        mat <- filter_mat_by_cell_types(mc_fp, cell_types, metacell_types)
    } else {
        mat <- mc_fp
    }

    mc_order <- order_mc_by_most_var_genes(mat, force_cell_type = force_cell_type, metacell_types = metacell_types)

    mat <- mat[, mc_order]

    return(mat)
}
