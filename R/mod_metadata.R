#' metadata UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_metadata_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 12,
                shinydashboardPlus::box(
                    id = ns("metadata_projection"),
                    title = "Metadata",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "80vh",
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("metadata_projection_sidebar"),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("stroke_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    uiOutput(ns("manifold_metadata_select_ui")),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_metadata_proj_2d"), height = "80vh")
                    )
                )
            )
        ),
        fluidRow(
            column(
                width = 6,
                shinydashboardPlus::box(
                    id = ns("md_md_box"),
                    title = "Comparison",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("md_md_sidebar"),
                        uiOutput(ns("md_md_point_size_ui")),
                        uiOutput(ns("md_md_stroke_ui"))
                    ),
                    uiOutput(ns("md1_select")),
                    uiOutput(ns("md2_select")),
                    uiOutput(ns("color_by_md_select")),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_md_md_mc"))
                    )
                )
            )
        )
    )
}


#' metadata sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_metadata_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list()
    )
}

#' metadata Server Function
#'
#' @noRd
mod_metadata_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    values <- reactiveValues()

    # Manifold metadata select
    output$manifold_metadata_select_ui <- renderUI({
        req(dataset())
        if (!has_metadata(dataset())) {
            print(glue("Dataset \"{dataset()}\" doesn't have any metadata. Use `update_metadata` to add it to your dataset."))
        } else {
            shinyWidgets::pickerInput(
                ns("color_proj"),
                label = "Color by:",
                choices = c("Cell type", dataset_metadata_fields(dataset())),
                selected = dataset_metadata_fields(dataset())[1],
                width = "70%",
                multiple = FALSE
            )
        }
    })

    # Point size selector
    output$point_size_ui <- renderUI({
        numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset()), min = 0.1, max = 3, step = 0.1)
    })

    output$md_md_point_size_ui <- renderUI({
        numericInput(ns("md_md_point_size"), label = "Point size", value = initial_scatters_point_size(dataset()), min = 0.05, max = 3, step = 0.1)
    })

    output$md_md_stroke_ui <- renderUI({
        numericInput(ns("md_md_stroke"), label = "Stroke width", value = initial_scatters_stroke(dataset()), min = 0, max = 3, step = 0.01)
    })

    # Minimal edge length selector
    output$edge_distance_ui <- renderUI({
        sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
    })

    output$stroke_ui <- renderUI({
        numericInput(ns("stroke"), label = "Stroke width", value = initial_proj_stroke(dataset()), min = 0, max = 3, step = 0.01)
    })

    # Projection plots
    output$plot_metadata_proj_2d <- render_2d_plotly(input, output, session, dataset, values, metacell_types, cell_type_colors, source = "proj_metadata_plot")

    # Metadata/Metadata plots
    output$md1_select <- renderUI({
        req(dataset())
        shinyWidgets::pickerInput(ns("md1"), "Metadata A",
            choices = dataset_metadata_fields(dataset()),
            selected = dataset_metadata_fields(dataset())[1], multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
        )
    })

    output$md2_select <- renderUI({
        req(dataset())
        req(length(dataset_metadata_fields(dataset())) > 1)
        shinyWidgets::pickerInput(ns("md2"), "Metadata B",
            choices = dataset_metadata_fields(dataset()),
            selected = dataset_metadata_fields(dataset())[2], multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
        )
    })

    output$color_by_md_select <- renderUI({
        req(dataset())
        shinyWidgets::pickerInput(ns("color_by_md"), "Color by",
            choices = c("Cell type", dataset_metadata_fields(dataset())),
            selected = "Cell type", multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
        )
    })



    output$plot_md_md_mc <- plotly::renderPlotly({
        req(input$md1)
        req(input$md2)
        req(input$color_by_md)
        req(input$md_md_point_size)
        req(input$md_md_stroke)

        color_by_md <- input$color_by_md
        if (input$color_by_md == "Cell type") {
            color_by_md <- NULL
        }

        fig <- plot_md_md_mc(
            dataset(),
            input$md1,
            input$md2,
            color_by_md,
            metacell_types = metacell_types(),
            cell_type_colors = cell_type_colors(),
            point_size = input$md_md_point_size,
            stroke = input$md_md_stroke,
            plot_text = FALSE
        ) %>%
            plotly::ggplotly(tooltip = "tooltip_text", source = "md_md_plot") %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons()

        if (input$color_by_md == "Cell type") {
            fig <- plotly::hide_legend(fig)
        } else {
            # This ugly hack is due to https://github.com/ropensci/plotly/issues/1234
            # We need to remove the legend generated by scale_color_identity
            fig$x$data <- fig$x$data %>% purrr::map(~ {
                .x$showlegend <- FALSE
                .x
            })
        }



        return(fig)
    })
}
