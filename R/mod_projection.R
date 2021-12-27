#' projection UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_projection_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 7,
                shinydashboardPlus::box(
                    id = ns("metacell_projection"),
                    title = "2D Projection",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("gene_projection_sidebar"),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("stroke_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinyWidgets::prettyRadioButtons(
                        ns("color_proj"),
                        label = "Color by:",
                        choices = c("Cell type", "Charting", "Metadata"),
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_proj_2d"))
                    )
                )
            ),
            column(
                width = 5,
                shinydashboardPlus::box(
                    title = "Diff. Expression",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_mc_gene_scatter"))
                    ),
                    shinyWidgets::prettySwitch(inputId = ns("show_diff_expr_table"), value = FALSE, label = "Show table"),
                    DT::DTOutput(ns("diff_expr_table"))
                )
            )
        ),
        fluidRow(
            column(
                width = 7,
                shinydashboardPlus::box(
                    id = ns("metacell_projection"),
                    title = "Type predictions",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_stacked_type"))
                    )
                )
            ),
            column(
                width = 3,
                uiOutput(ns("group_box"))
            )
        )
    )
}


#' projection sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_projection_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            div(
                id = ns("sidebar_select"),
                shinyWidgets::radioGroupButtons(
                    inputId = ns("mode"),
                    label = "Compare:",
                    choices = c(
                        "MC",
                        "Type",
                        "Group"
                    ),
                    justified = TRUE
                ),
                uiOutput(ns("diff_select"))
            )
        )
    )
}

#' projection Server Function
#'
#' @noRd
mod_projection_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    group <- reactiveVal()

    projected_metacell_types <- reactive({
        get_mc_data(dataset(), "projected_metacell_types") %>%
            mutate(metacell = as.character(metacell))
    })

    atlas_colors <- reactive({
        req(has_atlas(dataset()))
        get_mc_data(dataset(), "cell_type_colors", atlas = TRUE)
    })

    metacell_names <- reactive({
        req(dataset())
        colnames(get_mc_data(dataset(), "mc_mat"))
    })

    scatter_selectors(ns, dataset, output)
    projection_selectors(ns, dataset, output, input)

    group_selectors_mod_projection(input, output, session, dataset, ns, group)
    metacell_selectors_mod_projection(input, output, session, dataset, ns, metacell_names, projected_metacell_types, atlas_colors, group)

    mc_mc_gene_scatter_df <- reactive({
        req(input$mode)
        if (input$mode == "MC") {
            req(input$metacell1)
            calc_obs_exp_mc_df(dataset(), input$metacell1)
        } else if (input$mode == "Type") {
            req(input$metacell1)
            req(input$metacell1 %in% atlas_colors()$cell_type)
            req(projected_metacell_types())
            calc_obs_exp_type_df(dataset(), input$metacell1, projected_metacell_types())
        } else if (input$mode == "Group") {
            req(group())
            req(length(group()) > 1)
            group_types_df <- tibble(metacell = group(), cell_type = "Group")
            calc_obs_exp_type_df(dataset(), "Group", group_types_df)
        }
    })


    # Differential expression
    output$plot_mc_mc_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, mc_mc_gene_scatter_df, metacell_names(), atlas_colors())

    output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, mc_mc_gene_scatter_df)

    # Point size selector
    output$point_size_ui <- renderUI({
        numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset()), min = 0.1, max = 3, step = 0.1)
    })

    # Minimal edge length selector
    output$edge_distance_ui <- renderUI({
        sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
    })

    # Projection plots
    output$plot_mc_proj_2d <- render_2d_plotly(input, output, session, dataset, projected_metacell_types, atlas_colors, source = "proj_mc_plot_proj_tab")


    output$plot_mc_stacked_type <- plot_type_predictions_bar(dataset)
}

metacell_selectors_mod_projection <- function(input, output, session, dataset, ns, metacell_names, metacell_types, cell_type_colors, group) {
    output$diff_select <- renderUI({
        req(dataset())
        req(input$mode)
        if (input$mode == "MC") {
            shinyWidgets::pickerInput(ns("metacell1"), "Metacell",
                choices = metacell_names(),
                selected = config$selected_mc1, multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
            )
        } else if (input$mode == "Type") {
            req(cell_type_colors())
            cell_types_hex <- col2hex(cell_type_colors()$color)
            cell_types <- cell_type_colors()$cell_type
            shinyWidgets::pickerInput(ns("metacell1"), "Cell type",
                choices = cell_types,
                selected = cell_types[1],
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith"),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        } else if (input$mode == "Group") {
            tagList(
                shinyWidgets::pickerInput(ns("metacell"), "Metacell",
                    choices = metacell_names(),
                    selected = config$selected_mc1, multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
                ),
                shinyWidgets::actionGroupButtons(
                    ns("add_metacell_to_group"),
                    labels = c("Add to group"),
                    size = "sm"
                )
            )
        }
    })

    # Select metacell / cell type when clicking on it
    observeEvent(plotly::event_data("plotly_click", source = "proj_mc_plot_proj_tab"), {
        el <- plotly::event_data("plotly_click", source = "proj_mc_plot_proj_tab")
        metacell <- el$customdata

        if (input$mode == "MC") {
            updateSelectInput(session, "metacell1", selected = metacell)
            showNotification(glue("Selected Metacell: #{metacell}"))
        } else if (input$mode == "Type") {
            cell_type <- metacell_types() %>%
                filter(metacell == !!metacell) %>%
                slice(1) %>%
                pull(cell_type)
            updateSelectInput(session, "metacell1", selected = cell_type)
            showNotification(glue("Selected Cell type: {cell_type}"))
        } else if (input$mode == "Group") {
            if (is.null(group())) {
                group(metacell)
            } else {
                group(unique(c(group(), metacell)))
            }
        }
    })
}

group_selectors_mod_projection <- function(input, output, session, dataset, ns, group) {
    output$group_box <- renderUI({
        req(input$mode == "Group")
        shinydashboardPlus::box(
            id = ns("group_box_1"),
            title = "Group metacells",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            actionButton(ns("reset_group"), "Reset"),
            shinycssloaders::withSpinner(
                DT::dataTableOutput(ns("group_table"))
            )
        )
    })

    output$group_table <- DT::renderDataTable(
        tibble(metacell = group()),
        escape = FALSE,
        server = FALSE,
        rownames = FALSE,
        filter = "none",
        options = list(
            dom = "t",
            paging = FALSE,
            language = list(emptyTable = "Please select metacells")
        )
    )

    observeEvent(input$add_metacell_to_group, {
        if (is.null(group())) {
            group(input$metacell)
        } else {
            group(unique(c(group(), input$metacell)))
        }
    })

    observeEvent(input$reset_group, {
        group(NULL)
    })

    observeEvent(plotly::event_data("plotly_selected", source = "proj_mc_plot_proj_tab"), {
        el <- plotly::event_data("plotly_selected", source = "proj_mc_plot_proj_tab")

        selected_metacells <- el$customdata
        req(input$mode == "Group")


        if (is.null(group())) {
            group(selected_metacells)
        } else {
            group(unique(c(group(), selected_metacells)))
        }
    })
}
