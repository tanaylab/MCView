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
                        selectInput(ns("proj_stat"), label = "Statistic", choices = c("Expression" = "expression", "Enrichment" = "enrichment"), selected = "Expression", multiple = FALSE, selectize = FALSE),
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

    # Expression range
    output$set_range_ui <- renderUI({
        req(input$proj_stat == "expression")
        checkboxInput(ns("set_range"), "Manual range", value = FALSE)
    })

    output$expr_range_ui <- renderUI({
        req(input$proj_stat == "expression")
        req(input$set_range)
        shinyWidgets::numericRangeInput(ns("expr_range"), "Expression range", c(-18, -5), width = "80%", separator = " to ")
    })

    # Enrichment range
    output$enrich_range_ui <- renderUI({
        req(input$proj_stat == "enrichment")
        shinyWidgets::numericRangeInput(ns("lfp"), "Enrichment range", c(-3, 3), width = "80%", separator = " to ")
    })

    # Point size selectors
    output$point_size_ui <- renderUI({
        numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset()), min = 0.1, max = 3, step = 0.1)
    })

    output$gene_gene_point_size_ui <- renderUI({
        numericInput(ns("gene_gene_point_size"), label = "Point size", value = initial_scatters_point_size(dataset()), min = 0.05, max = 3, step = 0.1)
    })

    output$gene_gene_stroke_ui <- renderUI({
        numericInput(ns("gene_gene_stroke"), label = "Stroke width", value = initial_scatters_stroke(dataset()), min = 0, max = 3, step = 0.01)
    })

    output$stroke_ui <- renderUI({
        numericInput(ns("stroke"), label = "Stroke width", value = initial_proj_stroke(dataset()), min = 0, max = 3, step = 0.01)
    })

    # Minimal edge length selector
    output$edge_distance_ui <- renderUI({
        sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
    })

    # Projection plots
    output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, values, metacell_types, cell_type_colors, source = "proj_mc_plot_gene_tab")

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

axis_selector <- function(axis, selected, ns) {
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
                fill = TRUE,
                selected = selected
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

top_correlated_selector <- function(gene_id, id, type_id, input, output, session, dataset, ns) {
    output[[glue("top_correlated_select_{id}")]] <- renderUI({
        req(has_gg_mc_top_cor(project, dataset()))
        req(input[[type_id]] == "Gene")
        req(input[[gene_id]])
        gene <- input[[gene_id]]
        req(gene %in% gene_names)
        tagList(
            selectInput(
                ns(glue("selected_top_{id}")),
                glue("Top correlated to {gene}:"),
                choices = c(get_top_cor_gene(dataset(), gene, type = "pos"), rev(get_top_cor_gene(dataset(), gene, type = "neg"))),
                selected = NULL,
                size = 10,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(
                c(ns(glue("select_top_cor_{id}_x")), ns(glue("select_top_cor_{id}_y")), ns(glue("select_top_cor_{id}_color")), ns(glue("select_top_cor_{id}_proj2d"))),
                labels = c("X", "Y", "Color", "2D"), size = "sm", fullwidth = FALSE
            ),
            shiny::actionButton(
                inputId = ns(glue("genecards_{id}")), label = glue("GeneCards: {gene}"),
                size = "sm", onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}')")
            )
        )
    })
    observeEvent(input[[glue("select_top_cor_{id}_x")]], {
        req(input[["x_axis_type"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "x_axis_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_y")]], {
        req(input[["y_axis_type"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "y_axis_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_color")]], {
        req(input[["color_by_type"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "color_by_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_proj2d")]], {
        req(input[["color_proj"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "color_proj_gene", selected = input[[glue("selected_top_{id}")]])
    })
}

top_correlated_selectors <- function(input, output, session, dataset, ns) {
    top_correlated_selector("x_axis_var", "x_axis", "x_axis_type", input, output, session, dataset, ns)
    top_correlated_selector("y_axis_var", "y_axis", "y_axis_type", input, output, session, dataset, ns)
    top_correlated_selector("color_by_var", "color_by", "color_by_type", input, output, session, dataset, ns)
    top_correlated_selector("color_proj_gene", "color_proj", "color_proj", input, output, session, dataset, ns)
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
