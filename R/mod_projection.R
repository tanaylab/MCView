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
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_proj_2d"))
                    ),
                    shinyWidgets::prettyRadioButtons(
                        ns("color_proj"),
                        label = "Color by:",
                        choices = c("Cell type", "Charting", "Metadata"),
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    ),
                    uiOutput(ns("projection_selectors"))
                ),
                uiOutput(ns("metacell_flow_box"))
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
                ),
                uiOutput(ns("traj_genes_box"))
            )
        ),
        fluidRow(
            column(
                width = 3,
                uiOutput(ns("groupA_box"))
            ),
            column(
                width = 3,
                uiOutput(ns("groupB_box"))
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
                        "MCs",
                        "Types",
                        "Groups"
                    ),
                    justified = TRUE
                ),
                uiOutput(ns("metacell1_select")),
                uiOutput(ns("metacell2_select")),
                shinyWidgets::actionGroupButtons(ns("switch_metacells"), labels = c("Switch"), size = "sm")
            )
        )
    )
}

#' projection Server Function
#'
#' @noRd
mod_projection_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    groupA <- reactiveVal()
    groupB <- reactiveVal()

    metacell_names <- reactive({
        req(dataset())
        colnames(get_mc_data(dataset(), "mc_mat"))
    })

    group_selectors(input, output, session, dataset, ns, groupA, groupB)
    metacell_selectors(input, output, session, dataset, ns, metacell_names, metacell_types, cell_type_colors, groupA, groupB)

    mc_mc_gene_scatter_df <- reactive({
        req(input$mode)
        if (input$mode == "MCs") {
            calc_mc_mc_gene_df(dataset(), input$metacell1, input$metacell2)
        } else if (input$mode == "Types") {
            req(metacell_types())
            browser()
            req(input$metacell1 %in% cell_type_colors()$cell_type)
            req(input$metacell2 %in% cell_type_colors()$cell_type)
            calc_ct_ct_gene_df(dataset(), input$metacell1, input$metacell2, metacell_types())
        } else if (input$mode == "Groups") {
            req(groupA())
            req(groupB())
            group_types_df <- bind_rows(
                tibble(metacell = groupA(), cell_type = "Group A"),
                tibble(metacell = groupB(), cell_type = "Group B")
            )

            calc_ct_ct_gene_df(dataset(), "Group A", "Group B", group_types_df)
        }
    })

    observeEvent(input$switch_metacells, {
        if (input$mode == "Groups") {
            temp <- groupA()
            groupA(groupB())
            groupB(temp)
        } else {
            mc1 <- input$metacell1
            mc2 <- input$metacell2
            updateSelectInput(session, "metacell1", selected = mc2)
            updateSelectInput(session, "metacell2", selected = mc1)
        }
    })

    output$projection_selectors <- renderUI({
        req(input$mode)
        if (input$mode == "MCs") {
            choices <- c("Metacell A", "Metacell B")
            label <- "Select on click:"
        } else if (input$mode == "Types") {
            choices <- c("Cell type A", "Cell type B")
            label <- "Select on click:"
        } else {
            choices <- c("Group A", "Group B")
            label <- "Select:"
        }

        shinyWidgets::prettyRadioButtons(
            inputId = ns("proj_select_main"),
            label = label,
            choices = choices,
            inline = TRUE,
            status = "danger",
            fill = TRUE
        )
    })

    # Differential expression
    output$plot_mc_mc_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, mc_mc_gene_scatter_df, metacell_names(), cell_type_colors())

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
    output$plot_mc_proj_2d <- plotly::renderPlotly({
        req(input$point_size)
        req(input$min_edge_size)
        req(input$metacell1)
        req(input$metacell2)

        if (input$mode == "MCs") {
            highlight <- tibble::tibble(
                metacell = c(input$metacell1, input$metacell2),
                label = c("metacell1", "metacell2"),
                color = c("darkred", "darkblue")
            )
        } else {
            highlight <- NULL
        }

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
                    rm_plotly_grid() %>%
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
                rm_plotly_grid() %>%
                plotly::layout(annotations = get_plotly_subplot_title("Projection")) %>%
                sanitize_for_WebGL() %>%
                plotly::toWebGL() %>%
                sanitize_plotly_buttons() %>%
                arrange_2d_proj_tooltip()
            if (input$mode == "Groups") {
                fig <- fig %>% plotly::layout(dragmode = "select")
            }
            fig
        }

        fig$x$source <- "proj_mc_plot"
        return(fig)
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
                selectizeInput(ns("traj_genes"), "Select genes:", choices = gene_names(dataset()), selected = top_var_genes(), multiple = TRUE)
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
        updateSelectizeInput(session, "traj_genes", choices = gene_names(dataset()), server = TRUE, selected = top_var_genes(), options = list(maxItems = 5, maxOptions = 1e5))
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


    # metacell click observers
    metacell_click_observer("proj_manifold_plot", session)
    metacell_click_observer("md_md_plot", session)
    metacell_click_observer("gene_gene_plot", session)
    metacell_click_observer("proj_metadata_plot", session)
    metacell_click_observer("proj_mc_plot_gene_tab", session)
    metacell_click_observer("gene_time_mc_plot1", session)
    metacell_click_observer("gene_time_mc_plot2", session)
}
