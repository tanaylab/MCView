#' gene_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_gene_mc_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 7,
                shinydashboardPlus::box(
                    id = ns("gene_projection"),
                    title = "Gene projections",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("gene_projection_sidebar"),
                        uiOutput(ns("proj_stat_ui")),
                        uiOutput(ns("set_range_ui")),
                        uiOutput(ns("expr_range_ui")),
                        uiOutput(ns("enrich_range_ui")),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("stroke_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    uiOutput(ns("manifold_select_ui")),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_proj_2d"))
                    ),
                    shinyWidgets::prettyRadioButtons(
                        ns("color_proj"),
                        label = "Color by:",
                        choices = c("Cell type", "Gene", "Metadata"),
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    )
                )
            ),
            column(
                width = 5,
                shinydashboardPlus::box(
                    id = ns("gene_gene_box"),
                    title = "Gene/Gene",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("gene_gene_sidebar"),
                        uiOutput(ns("gene_gene_point_size_ui")),
                        uiOutput(ns("gene_gene_stroke_ui"))
                    ),
                    axis_selector("x_axis", "Gene", ns),
                    axis_selector("y_axis", "Gene", ns),
                    axis_selector("color_by", "Metadata", ns),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_gene_mc"))
                    )
                )
            )
        )
    )
}


#' gene_mc sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_gene_mc_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            uiOutput(ns("top_correlated_select_x_axis")),
            uiOutput(ns("top_correlated_select_y_axis")),
            uiOutput(ns("top_correlated_select_color_by")),
            uiOutput(ns("top_correlated_select_color_proj"))
        )
    )
}

#' gene_mc Server Function
#'
#' @noRd
mod_gene_mc_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    top_correlated_selectors(input, output, session, dataset, ns)
    mod_gene_mc_plotly_observers(input, session)

    # Manifold selectors
    output$manifold_select_ui <- renderUI({
        req(dataset())
        req(input$color_proj)
        picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
        if (input$color_proj == "Metadata") {
            if (!has_metadata(dataset())) {
                print(glue("Dataset \"{dataset()}\" doesn't have any metadata. Use `update_metadata` to add it to your dataset."))
            } else {
                shinyWidgets::pickerInput(
                    ns("color_proj_metadata"),
                    label = "Color by:",
                    choices = c("Cell type", dataset_metadata_fields(dataset())),
                    selected = dataset_metadata_fields(dataset())[1],
                    width = "70%",
                    multiple = FALSE,
                    options = picker_options
                )
            }
        } else if (input$color_proj == "Gene") {
            shinyWidgets::pickerInput(
                ns("color_proj_gene"),
                label = "Color by:",
                choices = gene_names,
                selected = default_gene1,
                width = "70%",
                multiple = FALSE,
                options = picker_options
            )
        }
    })

    scatter_selectors(ns, dataset, output)
    projection_selectors(ns, dataset, output, input)

    # Projection plots
    output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, values, metacell_types, cell_type_colors, source = "proj_mc_plot_gene_tab")

    # Metadata/Metadata plots
    output$x_axis_select <- render_axis_select_ui("x_axis", "X axis", md_choices = dataset_metadata_fields(dataset()), md_selected = dataset_metadata_fields(dataset())[1], selected_gene = default_gene1, input = input, ns = ns, dataset = dataset)

    output$y_axis_select <- render_axis_select_ui("y_axis", "Y axis", md_choices = dataset_metadata_fields(dataset()), md_selected = dataset_metadata_fields(dataset())[2], selected_gene = default_gene2, input = input, ns = ns, dataset = dataset)

    output$color_by_select <- render_axis_select_ui("color_by", "Color", md_choices = c("Cell type", dataset_metadata_fields(dataset())), md_selected = "Cell type", selected_gene = default_gene1, input = input, ns = ns, dataset = dataset)

    output$plot_gene_gene_mc <- plotly::renderPlotly({
        req(input$x_axis_var)
        req(input$y_axis_var)
        req(input$color_by_var)
        req(input$x_axis_type)
        req(input$y_axis_type)
        req(input$color_by_type)
        req(input$gene_gene_point_size)
        req(input$gene_gene_stroke)

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
            point_size = input$gene_gene_point_size,
            stroke = input$gene_gene_stroke,
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


mod_gene_mc_plotly_observers <- function(input, session) {
    observeEvent(plotly::event_data("plotly_click", source = "mc_mc_plot"), {
        el <- plotly::event_data("plotly_click", source = "mc_mc_plot")

        gene <- el$customdata
        if (input$x_axis_type == "Gene") {
            shinyWidgets::updatePickerInput(session, "x_axis_var", selected = gene)
        }
        showNotification(glue("Selected {gene} in \"Genes\" tab"))
    })
}
