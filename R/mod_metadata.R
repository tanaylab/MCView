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
                width = 7,
                shinydashboardPlus::box(
                    id = ns("metadata_projection"),
                    title = "Metadata",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "65vh",
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
                        plotly::plotlyOutput(ns("plot_metadata_proj_2d"), height = "65vh")
                    )
                )
            ),
            column(
                width = 5,
                shinydashboardPlus::box(
                    id = ns("md_md_box"),
                    title = "Metadata Comparison",
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
                    axis_selector("x_axis", ns),
                    axis_selector("y_axis", ns),
                    axis_selector("color_by", ns),
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
    output$x_axis_select <- render_axis_select_ui("x_axis", "X axis", md_choices = dataset_metadata_fields(dataset()), md_selected = dataset_metadata_fields(dataset())[1], selected_gene = default_gene1, input = input, ns = ns, dataset = dataset)

    output$y_axis_select <- render_axis_select_ui("y_axis", "Y axis", md_choices = dataset_metadata_fields(dataset()), md_selected = dataset_metadata_fields(dataset())[2], selected_gene = default_gene2, input = input, ns = ns, dataset = dataset)

    output$color_by_select <- render_axis_select_ui("color_by", "Color", md_choices = c("Cell type", dataset_metadata_fields(dataset())), md_selected = "Cell type", selected_gene = default_gene1, input = input, ns = ns, dataset = dataset)

    axis_vars_ok <- function(dataset, input) {
        metadata <- get_mc_data(dataset, "metadata")
        vars_ok <- purrr::map_lgl(c("x_axis", "y_axis", "color_by"), function(v) {
            type <- input[[glue("{v}_type")]]
            var <- input[[glue("{v}_var")]]
            if (type == "Metadata" && (var %in% c(colnames(metadata), "Cell type"))) {
                return(TRUE)
            } else if (type == "Gene" && var %in% gene_names) {
                return(TRUE)
            } else {
                return(FALSE)
            }
        })
        return(all(vars_ok))
    }

    output$plot_md_md_mc <- plotly::renderPlotly({
        req(input$x_axis_var)
        req(input$y_axis_var)
        req(input$color_by_var)
        req(input$x_axis_type)
        req(input$y_axis_type)
        req(input$color_by_type)
        req(input$md_md_point_size)
        req(input$md_md_stroke)

        req(axis_vars_ok(dataset(), input))

        color_var <- input$color_by_var
        if (input$color_by_var == "Cell type") {
            color_var <- NULL
        }

        fig <- plot_mc_scatter(
            dataset(),
            input$x_axis_var,
            input$y_axis_var,
            color_var,
            x_type = input$x_axis_type,
            y_type = input$y_axis_type,
            color_type = input$color_by_type,
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

        if (input$color_by_var == "Cell type") {
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


axis_selector <- function(axis, ns) {
    fluidRow(
        column(
            width = 9,
            uiOutput(ns(glue("{axis}_select")))
        ),
        column(
            width = 3,
            shinyWidgets::prettyRadioButtons(
                ns(glue("{axis}_type")),
                label = "",
                choices = c("Metadata", "Gene"),
                inline = TRUE,
                status = "danger",
                fill = TRUE
            )
        )
    )
}

render_axis_select_ui <- function(axis, title, md_choices, md_selected, selected_gene, ns, input, dataset) {
    picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")

    renderUI({
        req(dataset())
        req(input[[glue("{axis}_type")]])
        if (input[[glue("{axis}_type")]] == "Metadata") {
            shinyWidgets::pickerInput(
                ns(glue("{axis}_var")),
                title,
                choices = md_choices,
                selected = md_selected,
                multiple = FALSE,
                options = picker_options
            )
        } else if (input[[glue("{axis}_type")]] == "Gene") {
            shinyWidgets::pickerInput(
                ns(glue("{axis}_var")),
                title,
                choices = gene_names,
                selected = selected_gene,
                multiple = FALSE,
                options = picker_options
            )
        }
    })
}
