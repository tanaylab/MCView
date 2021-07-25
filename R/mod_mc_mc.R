#' mc_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_mc_mc_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 7,
                shinydashboardPlus::box(
                    id = ns("metacell_projection"),
                    title = "Metacells",
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
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_proj_2d"))
                    ),
                    shinyWidgets::prettyRadioButtons(
                        inputId = ns("proj_select_main"),
                        label = "Select on click:",
                        choices = c("Metacell A", "Metacell B"),
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    )
                ),
                uiOutput(ns("metacell_flow_box"))
            ),
            column(
                width = 5,
                shinydashboardPlus::box(
                    title = "Gene expression",
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
                ),
                uiOutput(ns("traj_genes_box"))
            )
        )
    )
}


#' mc_mc sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_mc_mc_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            div(
                id = ns("sidebar_select"),
                uiOutput(ns("metacell1_select")),
                uiOutput(ns("metacell2_select")),
                shinyWidgets::actionGroupButtons(ns("switch_metacells"), labels = c("Switch"), size = "sm")
            )
        )
    )
}

#' mc_mc Server Function
#'
#' @noRd
mod_mc_mc_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    metacell_names <- reactive({
        req(dataset())
        colnames(get_mc_data(dataset(), "mc_mat"))
    })

    output$metacell1_select <- renderUI({
        req(dataset())
        shinyWidgets::pickerInput(ns("metacell1"), "Metacell A",
            choices = metacell_names(),
            selected = config$selected_mc1, multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
        )        
    })

    output$metacell2_select <- renderUI({
        req(dataset())
        shinyWidgets::pickerInput(ns("metacell2"), "Metacell B",
            choices = metacell_names(),
            selected = config$selected_mc2, multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
        )        
    })


    mc_mc_gene_scatter_df <- reactive({
        calc_mc_mc_gene_df(dataset(), input$metacell1, input$metacell2)
    })

    observeEvent(input$switch_metacells, {
        mc1 <- input$metacell1
        mc2 <- input$metacell2
        updateSelectInput(session, "metacell1", selected = mc2)
        updateSelectInput(session, "metacell2", selected = mc1)
    })

    # MC/MC plots
    output$plot_mc_mc_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, mc_mc_gene_scatter_df)

    output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, mc_mc_gene_scatter_df)

    update_traj_genes <- function(gene, input, session) {
        if (length(input$traj_genes) == 4) {
            new_selection <- unique(c(input$traj_genes[-1], gene))
        } else {
            new_selection <- unique(c(input$traj_genes, gene))
        }

        updateSelectInput(session, "traj_genes", selected = new_selection)
        showNotification(glue("Selected {gene}"))
    }

    observeEvent(input$diff_expr_table_cell_clicked, {
        gene <- input$diff_expr_table_cell_clicked$value
        update_traj_genes(gene, input, session)
    })


    observeEvent(plotly::event_data("plotly_click", source = "mc_mc_plot"), {
        el <- plotly::event_data("plotly_click", source = "mc_mc_plot")
        gene <- el$customdata

        update_traj_genes(gene, input, session)
    })

    # Point size selector
    output$point_size_ui <- renderUI({
        numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset()), min = 0.1, max = 3, step = 0.1)
    })

    # Minimal edge length selector
    output$edge_distance_ui <- renderUI({
        sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
    })

    # Projection plots
    output$plot_mc_proj_2d <- plotly::renderPlotly({
        req(input$point_size)
        req(input$min_edge_size)
        req(input$metacell1)
        req(input$metacell2)

        highlight <- tibble::tibble(
            metacell = c(input$metacell1, input$metacell2),
            label = c("metacell1", "metacell2"),
            color = c("darkred", "darkblue")
        )

        p_proj <- mc2d_plot_ggp(dataset(), metacell_types = metacell_types(), cell_type_colors = cell_type_colors(), point_size = input$point_size, min_d = input$min_edge_size, highlight = highlight)

        if (has_time(dataset())) {
            subfig <- plotly::subplot(
                plotly::ggplotly(plot_mc_time_dist(dataset(), input$metacell1, ylab = "# of cells (A)") + theme(axis.title.y = element_text(colour = "darkred")), tooltip = "tooltip_text"),
                plotly::ggplotly(plot_mc_time_dist(dataset(), input$metacell2, ylab = "# of cells (B)") + theme(axis.title.y = element_text(colour = "darkblue")), tooltip = "tooltip_text"),
                shareX = TRUE,
                shareY = TRUE,
                titleX = TRUE,
                titleY = TRUE,
                nrows = 2
            )


            fig <- plotly::subplot(
                plotly::ggplotly(p_proj, tooltip = "tooltip_text") %>%
                    plotly::hide_legend() %>%
                    plotly::layout(annotations = get_plotly_subplot_title("Projection")),
                subfig %>%
                    plotly::layout(annotations = get_plotly_subplot_title("Time distribution")),
                shareY = FALSE,
                shareX = FALSE,
                titleX = TRUE,
                titleY = TRUE,
                margin = 0.07
            ) %>%
                sanitize_for_WebGL() %>%
                plotly::toWebGL() %>%
                sanitize_plotly_buttons() %>%
                arrange_2d_proj_tooltip()
        } else {
            fig <- plotly::ggplotly(p_proj, tooltip = "tooltip_text") %>%
                plotly::hide_legend() %>%
                plotly::layout(annotations = get_plotly_subplot_title("Projection")) %>%
                sanitize_for_WebGL() %>%
                plotly::toWebGL() %>%
                sanitize_plotly_buttons() %>%
                arrange_2d_proj_tooltip()
        }

        fig$x$source <- "proj_mc_plot"
        return(fig)
    })

    # Select metacell when clicking on it
    observeEvent(plotly::event_data("plotly_click", source = "proj_mc_plot"), {
        el <- plotly::event_data("plotly_click", source = "proj_mc_plot")

        metacell <- el$customdata

        if (input$proj_select_main == "Metacell A") {
            updateSelectInput(session, "metacell1", selected = metacell)
            showNotification(glue("Selected Metacell A #{metacell}"))
        } else {
            updateSelectInput(session, "metacell2", selected = metacell)
            showNotification(glue("Selected Metacell B #{metacell}"))
        }
    })


    output$metacell_flow_box <- renderUI({
        req(has_network(dataset()))
        shinydashboardPlus::box(
            title = "Metacell flow",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            shinyWidgets::prettyRadioButtons(
                inputId = ns("flow_metacell"),
                label = "",
                choices = c("Metacell A", "Metacell B"),
                inline = TRUE,
                status = "danger",
                fill = TRUE
            ),
            div(
                id = ns("flow_plot"),
                shinycssloaders::withSpinner(
                    plotOutput(ns("plot_metacell1_flow"), height = 600)
                )
            )
        )
    })

    # Metacell flow
    output$plot_metacell1_flow <- renderPlot({
        req(has_network(dataset()))
        if (input$flow_metacell == "Metacell A") {
            req(input$metacell1)
            mc_to_plot <- input$metacell1
            color <- "darkred"
        } else {
            req(input$metacell2)
            mc_to_plot <- input$metacell2
            color <- "darkblue"
        }

        mc_data <- metacell_types() %>% filter(metacell == mc_to_plot)
        plot_propagation_net_metacell(dataset(), mc_to_plot, metacell_types = metacell_types()) +
            ggtitle(glue("metacell #{mc_to_plot} ({mc_data$cell_type[1]})")) +
            theme(plot.title = element_text(colour = color))
    })


    observeEvent(plotly::event_data("plotly_click", source = "proj_mc_plot_gene_tab"), {
        el <- plotly::event_data("plotly_click", source = "proj_mc_plot_gene_tab")

        metacell <- el$customdata
        updateSelectInput(session, "metacell1", selected = metacell)
        showNotification(glue("Selected Metacell #{metacell}"))
    })

    observeEvent(plotly::event_data("plotly_click", source = "gene_time_mc_plot1"), {
        el <- plotly::event_data("plotly_click", source = "gene_time_mc_plot1")

        metacell <- el$customdata
        updateSelectInput(session, "metacell1", selected = metacell)
        showNotification(glue("Selected Metacell #{metacell}"))
    })

    observeEvent(plotly::event_data("plotly_click", source = "gene_time_mc_plot2"), {
        el <- plotly::event_data("plotly_click", source = "gene_time_mc_plot2")

        metacell <- el$customdata
        updateSelectInput(session, "metacell1", selected = metacell)
        showNotification(glue("Selected Metacell #{metacell}"))
    })

    # Gene trajectory plots
    output$traj_genes_box <- renderUI({
        req(has_network(dataset()))
        shinydashboardPlus::box(
            title = "Genes Trajectory",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            div(
                id = ns("traj_genes_help"),
                selectizeInput(ns("traj_genes"), "Select genes:", choices = gene_names, selected = top_var_genes(), multiple = TRUE)
            ),
            shinycssloaders::withSpinner(
                plotly::plotlyOutput(ns("plot_mc_traj"))
            ),
            shinyWidgets::prettyRadioButtons(
                inputId = ns("traj_metacell"),
                label = "",
                choices = c("Metacell A", "Metacell B"),
                inline = TRUE,
                status = "danger",
                fill = TRUE
            ),
            helpText("Click on a gene trace in order to add it to the \"Genes\" tab.")
        )
    })

    top_var_genes <- reactive({
        req(input$metacell1)
        req(has_network(dataset()))

        # Initially select the genes with highest variance along the trajectory
        get_top_var_genes(dataset(), input$metacell1)
    })

    observe({
        req(input$metacell1)
        req(has_network(dataset()))
        updateSelectizeInput(session, "traj_genes", choices = gene_names, server = TRUE, selected = top_var_genes(), options = list(maxItems = 5, maxOptions = 1e5))
    })

    output$plot_mc_traj <- plotly::renderPlotly({
        req(input$traj_genes)
        req(has_network(dataset()))
        if (input$traj_metacell == "Metacell A") {
            req(input$metacell1)
            mc_to_plot <- input$metacell1
            color <- "darkred"
        } else {
            req(input$metacell2)
            mc_to_plot <- input$metacell2
            color <- "darkblue"
        }

        plotly::ggplotly(plot_gene_trajectory(dataset(), input$traj_genes, mc_to_plot, anchor_gene = NULL) + theme(axis.title.y = element_text(color = color)), source = "traj_plot", tooltip = "tooltip_text") %>% sanitize_plotly_buttons()
    })


    # Output priorities
    outputOptions(output, "plot_mc_proj_2d", priority = 6)
    outputOptions(output, "plot_mc_mc_gene_scatter", priority = 5)
    outputOptions(output, "plot_metacell1_flow", priority = 4)
    outputOptions(output, "plot_mc_traj", priority = 3)
}
