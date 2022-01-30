#' cell type UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_cell_type_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            resizable_column(
                width = 12,
                shinydashboardPlus::box(
                    id = ns("boxplot_box"),
                    title = "Cell types",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "70vh",
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 50,
                        id = ns("gene_gene_sidebar"),
                        axis_selector("boxplot_axis", "Gene", ns, choices = c("Metadata", "Gene")),
                        uiOutput(ns("confusion_color_by_selector")),
                        # axis_selector("color_by", "Metadata", ns, choices = c("Metadata", "Gene")),
                        uiOutput(ns("cell_type_list"))
                        # uiOutput(ns("boxplot_point_size_ui")),
                        # uiOutput(ns("boxplot_stroke_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("cell_type_boxplot"), height = "70vh")
                    )
                )
            )
        )
    )
}


#' cell type sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_cell_type_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            div()
        )
    )
}

#' cell type Server Function
#'
#' @noRd
mod_cell_type_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals) {
    ns <- session$ns

    # scatter_selectors(ns, dataset, output, globals, "boxplot")
    top_correlated_selector("boxplot_axis_var", "boxplot_axis", "boxplot_axis_type", input, output, session, dataset, ns, button_labels = c("Axes", "Color"), ids = c("boxplot_axis", "color"))

    output$boxplot_axis_select <- render_axis_select_ui("boxplot_axis", "Data", md_choices = dataset_metadata_fields(dataset()), md_selected = dataset_metadata_fields(dataset())[1], selected_gene = default_gene1, input = input, ns = ns, dataset = dataset)

    output$cell_type_list <- cell_type_selector(dataset, ns, id = "boxplot_cell_types", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors)

    output$confusion_color_by_selector <- renderUI({
        shinyWidgets::prettyRadioButtons(
            ns("confusion_color_by"),
            label = "Normalize by:",
            choices = c("X axis", "Y Axis"),
            inline = TRUE,
            status = "danger",
            fill = TRUE
        )
    })

    observe({
        req(input$boxplot_axis_type)
        req(input$boxplot_axis_var)
        metadata <- get_mc_data(dataset(), "metadata")
        req(metadata)

        shinyjs::toggle(id = "confusion_color_by_selector", condition = input$boxplot_axis_type == "Metadata" && input$boxplot_axis_var %in% colnames(metadata) && !is_numeric_field(metadata, input$boxplot_axis_var))
    })

    output$cell_type_boxplot <- plotly::renderPlotly({
        req(input$boxplot_axis_type)
        req(dataset())
        req(metacell_types())
        req(cell_type_colors())
        req(input$boxplot_cell_types)
        req(input$boxplot_axis_var)

        if (input$boxplot_axis_type == "Gene") {
            req(input$boxplot_axis_var %in% gene_names(dataset()))
            p <- cell_type_gene_boxplot(
                input$boxplot_axis_var,
                dataset(),
                cell_types = input$boxplot_cell_types,
                metacell_types = metacell_types(),
                cell_type_colors = cell_type_colors()
            )
        } else {
            metadata <- get_mc_data(dataset(), "metadata")
            req(!is.null(metadata))
            req(input$boxplot_axis_var %in% colnames(metadata))
            if (is_numeric_field(metadata, input$boxplot_axis_var)) {
                p <- cell_type_metadata_boxplot(
                    input$boxplot_axis_var,
                    dataset(),
                    cell_types = input$boxplot_cell_types,
                    metacell_types = metacell_types(),
                    cell_type_colors = cell_type_colors()
                )
            } else {
                p <- cell_type_metadata_confusion(
                    input$boxplot_axis_var,
                    dataset(),
                    color_by = input$confusion_color_by,
                    cell_types = input$boxplot_cell_types,
                    metacell_types = metacell_types()
                )
            }
        }

        req(p)

        fig <- plotly::ggplotly(p, source = "cell_type_boxplot") %>%
            sanitize_plotly_buttons()

        return(fig)
    })
}