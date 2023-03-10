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
            generic_column(
                width = 7,
                scatter_box(ns, "gene_gene_box", x_selected = "Gene", y_selected = "Gene", color_selected = "Metadata", collapsed_accordion = FALSE),
                uiOutput(ns("atlas_gene_gene_box_ui"))
            ),
            generic_column(
                width = 5,
                projection_box(ns, "gene_projection", title = "Gene projections", collapsed_accordion = FALSE, show_legend = FALSE, color_choices = c("Scatter Axis", "Cell type", "Gene", "Gene module", "Metadata"))
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
mod_gene_mc_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            top_correlated_selectors(input, output, session, dataset, ns)
            mod_gene_mc_plotly_observers(input, session)
            mod_gene_mc_globals_observers(input, session, globals, dataset)

            scatter_selectors(ns, dataset, output, globals)
            projection_selectors(ns, dataset, output, input, gene_modules, globals, session, weight = 0.6)

            clipboard_changed <- clipboard_changed_2d_reactive(input, globals)


            # Projection plots
            output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, source = "proj_mc_plot_gene_tab") %>%
                bindCache(dataset(), input$color_proj, metacell_types(), cell_type_colors(), input$point_size, input$stroke, input$min_edge_size, input$set_range, input$metacell1, input$metacell2, input$proj_stat, input$expr_range, input$lfp, input$color_proj_gene, input$color_proj_metadata, input$color_proj_gene_module, clipboard_changed(), input$graph_name, input$legend_orientation, input$show_legend_projection, input$scatter_axis_proj, {
                    if (input$color_proj == "Scatter Axis") {
                        if (input$scatter_axis_proj == "x"){
                            c(input$x_axis_var, input$x_axis_type)
                        } else {
                            c(input$y_axis_var, input$y_axis_type)
                        }
                    }
                })

            connect_gene_plots(input, output, session, ns, source = "proj_mc_plot_gene_tab")

            scatter_box_outputs(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, ns, plotly_source = "md_md_plot")

            atlas_gene_gene(input, output, session, dataset, metacell_types, cell_type_colors, globals, ns)
        }
    )
}


mod_gene_mc_plotly_observers <- function(input, session, source = "mc_mc_plot", notification_suffix = " in \"Genes\" tab") {
    observeEvent(plotly::event_data("plotly_click", source = source), {
        el <- plotly::event_data("plotly_click", source = source)

        gene <- el$customdata
        req(input$x_axis_type == "Gene")
        shinyWidgets::updatePickerInput(session, "x_axis_var", selected = gene)
        showNotification(glue("Selected {gene}{notification_suffix}"))
    })
}

mod_gene_mc_globals_observers <- function(input, session, globals, dataset, notification_suffix = " in \"Genes\" tab") {
    observe({
        req(globals$selected_gene_x_axis)
        req(input$x_axis_type == "Gene")
        req(input$x_axis_var)
        shinyWidgets::updatePickerInput(session, "x_axis_var", selected = globals$selected_gene_x_axis)

        if (has_atlas(dataset())) {
            shinyWidgets::updatePickerInput(session, "atlas_x_axis_var", selected = globals$selected_gene_x_axis)
        }

        showNotification(glue("Selected {globals$selected_gene_x_axis}{notification_suffix}"))
        globals$selected_gene_x_axis <- NULL
    })

    observe({
        req(globals$selected_gene_y_axis)
        req(input$y_axis_type == "Gene")
        req(input$y_axis_var)
        shinyWidgets::updatePickerInput(session, "y_axis_var", selected = globals$selected_gene_y_axis)

        if (has_atlas(dataset())) {
            shinyWidgets::updatePickerInput(session, "atlas_y_axis_var", selected = globals$selected_gene_y_axis)
        }

        showNotification(glue("Selected {globals$selected_gene_y_axis}{notification_suffix}"))
        globals$selected_gene_y_axis <- NULL
    })
}

atlas_gene_gene <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals, ns) {
    output$atlas_gene_gene_box_ui <- renderUI({
        req(has_atlas(dataset()))
        generic_box(
            id = ns("atlas_gene_gene_box"),
            title = "Atlas Gene/Gene",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            sidebar = shinydashboardPlus::boxSidebar(
                startOpen = FALSE,
                width = 100,
                id = ns("atlas_gene_gene_sidebar"),
                uiOutput(ns("atlas_gene_gene_xyline_ui")),
                uiOutput(ns("atlas_gene_gene_fixed_limits_ui")),
                checkboxInput(ns("use_query_limits"), label = "Use query limits", value = FALSE),
                uiOutput(ns("atlas_gene_gene_point_size_ui")),
                uiOutput(ns("atlas_gene_gene_stroke_ui"))
            ),
            shinycssloaders::withSpinner(
                plotly::plotlyOutput(ns("atlas_plot_gene_gene_mc"))
            ),
            shinydashboardPlus::accordion(
                id = ns("gene_gene_atlas_accordion"),
                shinydashboardPlus::accordionItem(
                    title = "Select axes",
                    axis_selector("atlas_x_axis", "Gene", ns),
                    axis_selector("atlas_y_axis", "Gene", ns),
                    axis_selector("atlas_color_by", "Metadata", ns),
                )
            )
        )
    })

    gene_modules <- reactive({
        mods <- get_mc_data(dataset(), "gene_modules", atlas = TRUE)
        if (is.null(mods)) {
            mods <- tibble(gene = character(0), module = factor())
        }
        return(mods)
    })

    scatter_selectors(ns, dataset, output, globals, prefix = "atlas_gene_gene")

    # Metadata/Metadata plots
    render_axis_select_ui("atlas_x_axis", "X axis", "atlas_x_axis_select", md_choices = dataset_metadata_fields_numeric(dataset(), atlas = TRUE), md_selected = dataset_metadata_fields_numeric(dataset(), atlas = TRUE)[1], selected_gene = default_gene1, input = input, output = output, ns = ns, dataset = dataset, gene_modules = gene_modules, session = session, atlas = TRUE)

    render_axis_select_ui("atlas_y_axis", "Y axis", "atlas_y_axis_select", md_choices = dataset_metadata_fields_numeric(dataset(), atlas = TRUE), md_selected = dataset_metadata_fields_numeric(dataset(), atlas = TRUE)[2], selected_gene = default_gene2, input = input, output = output, ns = ns, dataset = dataset, gene_modules = gene_modules, session = session, atlas = TRUE)

    render_axis_select_ui("atlas_color_by", "Color", "atlas_color_by_select", md_choices = c("Cell type", dataset_metadata_fields_numeric(dataset(), atlas = TRUE)), md_selected = "Cell type", selected_gene = default_gene1, input = input, output = output, ns = ns, dataset = dataset, gene_modules = gene_modules, session = session, atlas = TRUE)

    output$atlas_plot_gene_gene_mc <- plotly::renderPlotly({
        req(has_atlas(dataset()))
        req(input$atlas_x_axis_var)
        req(input$atlas_y_axis_var)
        req(input$atlas_color_by_var)
        req(input$atlas_x_axis_type)
        req(input$atlas_y_axis_type)
        req(input$atlas_color_by_type)
        req(input$atlas_gene_gene_point_size)
        req(input$atlas_gene_gene_stroke)
        req(!is.null(input$atlas_gene_gene_fixed_limits))
        req(axis_vars_ok(dataset(), input, "metadata", gene_modules, axes = c("atlas_x_axis", "atlas_y_axis", "atlas_color_by"), atlas = TRUE))

        color_var <- input$atlas_color_by_var
        if (input$atlas_color_by_var == "Cell type") {
            color_var <- NULL
        }

        x_limits <- NULL
        y_limits <- NULL
        if (input$use_query_limits) {
            if (input$atlas_x_axis_type == "Gene") {
                egc_x <- get_gene_egc(input$atlas_x_axis_var, dataset(), atlas = FALSE) + egc_epsilon
                x_limits <- c(min(egc_x), max(egc_x))
            }

            if (input$atlas_y_axis_type == "Gene") {
                egc_y <- get_gene_egc(input$atlas_y_axis_var, dataset(), atlas = FALSE) + egc_epsilon
                y_limits <- c(min(egc_y), max(egc_y))
            }
        }

        fig <- plot_mc_scatter(
            dataset(),
            input$atlas_x_axis_var,
            input$atlas_y_axis_var,
            color_var,
            gene_modules = gene_modules(),
            x_type = input$atlas_x_axis_type,
            y_type = input$atlas_y_axis_type,
            color_type = input$atlas_color_by_type,
            metacell_types = get_mc_data(dataset(), "metacell_types", atlas = TRUE),
            cell_type_colors = get_mc_data(dataset(), "cell_type_colors", atlas = TRUE),
            point_size = input$atlas_gene_gene_point_size,
            stroke = input$atlas_gene_gene_stroke,
            plot_text = FALSE,
            atlas = TRUE,
            x_limits = x_limits,
            y_limits = y_limits,
            fixed_limits = input$atlas_gene_gene_fixed_limits,
            xyline = input$atlas_gene_gene_xyline %||% FALSE
        ) %>%
            plotly::ggplotly(tooltip = "tooltip_text", source = "atlas_md_md_plot") %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons()

        if (input$atlas_color_by_var == "Cell type") {
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
    }) %>% bindCache(dataset(), input$atlas_x_axis_var, input$atlas_x_axis_type, input$atlas_y_axis_var, input$atlas_y_axis_type, input$atlas_color_by_type, input$atlas_color_by_var, input$atlas_gene_gene_point_size, input$atlas_gene_gene_stroke, input$use_query_limits, input$atlas_gene_gene_fixed_limits, input$atlas_gene_gene_xyline)
}
