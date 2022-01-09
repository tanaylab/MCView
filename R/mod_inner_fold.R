#' gene_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_inner_fold_ui <- function(id) {
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
                        checkboxInput(ns("force_cell_type"), "Force cell type", value = TRUE),
                        shinyWidgets::numericRangeInput(ns("lfp_range"), "Fold change range", c(-1, 3), width = "80%", separator = " to "),
                        colourpicker::colourInput(ns("low_color"), "Low color", "blue"),
                        colourpicker::colourInput(ns("high_color"), "High color", "red"),
                        colourpicker::colourInput(ns("mid_color"), "Mid color", "white"),
                        checkboxInput(ns("plot_legend"), "Plot legend", value = TRUE)
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
mod_inner_fold_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            uiOutput(ns("cell_type_list")),
            uiOutput(ns("metadata_list")),
            shinyWidgets::actionGroupButtons(ns("update_markers"), labels = "Update markers", size = "sm"),
            numericInput(ns("max_gene_num"), "Maximal number of genes", value = 100),
            uiOutput(ns("add_genes_ui")),
            uiOutput(ns("marker_genes_list"))
        )
    )
}

#' gene_mc Server Function
#'
#' @noRd
mod_inner_fold_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals) {
    ns <- session$ns

    markers <- reactiveVal()
    lfp_range <- reactiveVal()

    output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors())

    output$metadata_list <- metadata_selector(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata")

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

    observe({
        if (is.null(markers())) {
            req(input$max_gene_num)
            initial_markers <- choose_markers(get_marker_genes(dataset(), mode = "Inner"), input$max_gene_num, dataset = dataset(), add_systematic = FALSE)
            markers(initial_markers)
        }

        lfp_range(input$lfp_range)
    })

    markers_matrix <- reactive({
        req(markers)
        req(metacell_types())
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

        get_marker_matrix(
            dataset(),
            markers(),
            input$selected_cell_types,
            metacell_types(),
            force_cell_type = input$force_cell_type,
            mode = "Inner",
            notify_var_genes = TRUE
        )
    })


    observeEvent(input$update_markers, {
        req(metacell_types())
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

        if (!is.null(input$selected_cell_types)) {
            markers_df <- metacell_types() %>%
                filter(cell_type %in% input$selected_cell_types)
        } else {
            markers_df <- metacell_types()
        }

        markers_df <- markers_df %>%
            select(metacell) %>%
            inner_join(get_marker_genes(dataset(), mode = "Inner"), by = "metacell")

        req(input$max_gene_num)
        new_markers <- choose_markers(markers_df, input$max_gene_num, dataset = dataset(), add_systematic = FALSE)

        markers(new_markers)
    })


    observeEvent(input$remove_genes, {
        req(markers)
        new_markers <- markers()[!(markers() %in% input$selected_marker_genes)]
        markers(new_markers)
        shinyWidgets::updatePickerInput(session, ns("genes_to_add"), selected = c())
    })

    output$add_genes_ui <- renderUI({
        tagList(
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm"),
            shinyWidgets::pickerInput(ns("genes_to_add"),
                choices = gene_names(dataset()),
                selected = c(),
                multiple = TRUE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
            )
        )
    })

    observeEvent(input$add_genes, {
        new_markers <- sort(unique(c(markers(), input$genes_to_add)))
        markers(new_markers)
        shinyWidgets::updatePickerInput(session = session, inputId = "genes_to_add", selected = character(0))
    })

    output$markers_heatmap <- renderPlot({
        req(dataset())
        req(lfp_range())
        req(metacell_types())
        req(cell_type_colors())

        mat <- markers_matrix()
        req(mat)

        if (!is.null(input$selected_md)) {
            metadata <- get_mc_data(dataset(), "metadata") %>%
                select(metacell, one_of(input$selected_md))
        } else {
            metadata <- NULL
        }

        marker_genes <- markers()
        req(marker_genes)


        forbidden_genes <- get_mc_data(dataset(), "forbidden_genes")
        systematic_genes <- get_mc_data(dataset(), "systematic_genes")

        plot_markers_mat(
            mat,
            metacell_types(),
            cell_type_colors(),
            dataset(),
            min_lfp = lfp_range()[1],
            max_lfp = lfp_range()[2],
            plot_legend = input$plot_legend %||% TRUE,
            high_color =  input$high_color,
            low_color =  input$low_color,
            mid_color =  input$mid_color,
            metadata = metadata,
            forbidden_genes = forbidden_genes,
            systematic_genes = systematic_genes,
            disjoined_genes = NULL
        )
    }) %>% bindCache(dataset(), metacell_types(), cell_type_colors(), lfp_range(), input$plot_legend, input$selected_md, markers(), input$selected_cell_types, input$force_cell_type, input$high_color, input$low_color, input$mid_color)
}
