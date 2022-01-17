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
            resizable_column(
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
                        width = 80,
                        id = ns("gene_projection_sidebar"),
                        uiOutput(ns("projection_color_selectors")),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("stroke_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_proj_2d"))
                    ),
                    uiOutput(ns("projection_selectors"))
                ),
                uiOutput(ns("metacell_flow_box"))
            ),
            resizable_column(
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

#' mc_mc Server Function
#'
#' @noRd
mod_mc_mc_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals) {
    ns <- session$ns

    groupA <- reactiveVal()
    groupB <- reactiveVal()

    metacell_names <- metacell_names_reactive(dataset)
    metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)

    projection_selectors(ns, dataset, output, input, globals, weight = 0.6)
    group_selectors(input, output, session, dataset, ns, groupA, groupB, metacell_types, cell_type_colors)
    metacell_selectors(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors, groupA, groupB)

    mc_mc_gene_scatter_df <- mc_mc_gene_scatter_df_reactive(dataset, input, output, session, metacell_types, cell_type_colors, groupA, groupB)

    diff_expr_switch_metacells(dataset, input, output, session, groupA, groupB)

    output$projection_color_selectors <- renderUI({
        req(input$mode)
        if (input$mode == "MCs") {
            color_choices <- c("Cell type")
        } else if (input$mode == "Types") {
            color_choices <- c("Cell type")
        } else if (input$mode == "Groups") {
            color_choices <- c("Cell type", "Selected")
        } else {
            req(FALSE)
        }

        shinyWidgets::prettyRadioButtons(
            ns("color_proj"),
            label = "Color by:",
            choices = color_choices,
            inline = TRUE,
            status = "danger",
            fill = TRUE
        )
    })

    output$projection_selectors <- renderUI({
        req(input$mode)
        if (input$mode == "MCs") {
            choices <- c("Metacell A", "Metacell B")
            label <- "Select on click:"
        } else if (input$mode == "Types") {
            choices <- c("Cell type A", "Cell type B")
            label <- "Select on click:"
        } else if (input$mode == "Groups") {
            choices <- c("Group A", "Group B")
            label <- "Select:"
        } else {
            req(FALSE)
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

    # Projection plots
    output$plot_mc_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, groupA = groupA, groupB = groupB, source = "proj_mc_plot")

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
        updateSelectizeInput(session, "traj_genes", choices = gene_names(dataset()), server = TRUE, selected = top_var_genes(), options = list(maxItems = 5))
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

    # Output priorities
    outputOptions(output, "plot_mc_proj_2d", priority = 6)
    outputOptions(output, "plot_mc_mc_gene_scatter", priority = 5)
    outputOptions(output, "plot_metacell1_flow", priority = 4)
    outputOptions(output, "plot_mc_traj", priority = 3)
}

metacell_selectors <- function(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors, groupA, groupB) {
    output$metacell1_select <- renderUI({
        req(dataset())
        req(input$mode)
        if (input$mode == "MCs") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())
            shinyWidgets::pickerInput(ns("metacell1"), "Metacell A",
                choices = metacell_names(),
                selected = config$selected_mc1, multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        } else if (input$mode == "Types") {
            req(cell_type_colors())
            req(metacell_types())
            types_df <- cell_type_colors() %>% filter(cell_type %in% metacell_types()$cell_type)
            cell_types_hex <- col2hex(types_df$color)
            cell_types <- types_df$cell_type
            shinyWidgets::pickerInput(ns("metacell1"), "Cell type A",
                choices = cell_types,
                selected = cell_types[1],
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        } else if (input$mode == "Groups") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())
            tagList(
                shinyWidgets::pickerInput(ns("metacell"), "Metacell",
                    choices = metacell_names(),
                    selected = config$selected_mc1, multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
                    choicesOpt = list(
                        style = paste0("color: ", cell_types_hex, ";")
                    )
                ),
                shinyWidgets::actionGroupButtons(
                    c(ns("add_metacell_to_groupA"), ns("add_metacell_to_groupB")),
                    labels = c("Add to group A", "Add to group B"),
                    size = "sm"
                )
            )
        }
    })

    output$metacell2_select <- renderUI({
        req(dataset())
        req(input$mode)
        if (input$mode == "MCs") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())
            shinyWidgets::pickerInput(ns("metacell2"), "Metacell B",
                choices = metacell_names(),
                selected = config$selected_mc2, multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith"),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        } else if (input$mode == "Types") {
            req(cell_type_colors())
            types_df <- cell_type_colors() %>% filter(cell_type %in% metacell_types()$cell_type)
            cell_types_hex <- col2hex(types_df$color)
            cell_types <- types_df$cell_type
            shinyWidgets::pickerInput(ns("metacell2"), "Cell type B",
                choices = cell_types,
                selected = cell_types[2],
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith"),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        }
    })

    # Select metacell / cell type when clicking on it
    observeEvent(plotly::event_data("plotly_click", source = "proj_mc_plot"), {
        el <- plotly::event_data("plotly_click", source = "proj_mc_plot")
        metacell <- el$customdata
        req(input$proj_select_main)

        if (input$mode == "MCs") {
            if (input$proj_select_main == "Metacell A") {
                updateSelectInput(session, "metacell1", selected = metacell)
                showNotification(glue("Selected Metacell A: #{metacell}"))
            } else {
                updateSelectInput(session, "metacell2", selected = metacell)
                showNotification(glue("Selected Metacell B: #{metacell}"))
            }
        } else if (input$mode == "Types") {
            cell_type <- metacell_types() %>%
                filter(metacell == !!metacell) %>%
                slice(1) %>%
                pull(cell_type)
            if (input$proj_select_main == "Cell type A") {
                updateSelectInput(session, "metacell1", selected = cell_type)
                showNotification(glue("Selected Cell type A: {cell_type}"))
            } else {
                updateSelectInput(session, "metacell2", selected = cell_type)
                showNotification(glue("Selected Cell type B: {cell_type}"))
            }
        } else if (input$mode == "Groups") {
            if (input$proj_select_main == "Group A") {
                if (is.null(groupA)) {
                    groupA(metacell)
                } else {
                    groupA(unique(c(groupA(), metacell)))
                }
            }

            if (input$proj_select_main == "Group B") {
                if (is.null(groupB)) {
                    groupB(metacell)
                } else {
                    groupB(unique(c(groupB(), metacell)))
                }
            }
        }
    })
}

group_selectors <- function(input, output, session, dataset, ns, groupA, groupB, metacell_types, cell_type_colors) {
    output$groupA_box <- renderUI({
        req(input$mode == "Groups")
        shinydashboardPlus::box(
            id = ns("groupA_box_1"),
            title = "Group A metacells",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            actionButton(ns("reset_groupA"), "Reset"),
            actionButton(ns("remove_groupA_metacells"), "Remove"),
            shinycssloaders::withSpinner(
                DT::dataTableOutput(ns("groupA_table"))
            )
        )
    })

    output$groupB_box <- renderUI({
        req(input$mode == "Groups")
        shinydashboardPlus::box(
            id = ns("groupB_box_1"),
            title = "Group B metacells",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            actionButton(ns("reset_groupB"), "Reset"),
            actionButton(ns("remove_groupB_metacells"), "Remove"),
            shinycssloaders::withSpinner(
                DT::dataTableOutput(ns("groupB_table"))
            )
        )
    })

    output$groupA_table <- DT::renderDataTable(
        {
            req(metacell_types())
            req(cell_type_colors())
            req(groupA())
            DT::datatable(
                tibble(metacell = groupA()) %>%
                    left_join(metacell_types() %>% select(metacell, cell_type), by = "metacell"),
                escape = FALSE,
                rownames = FALSE,
                filter = "none",
                options = list(
                    dom = "t",
                    paging = FALSE,
                    language = list(emptyTable = "Please select metacells"),
                    columnDefs = list(list(visible = FALSE, targets = c(1)))
                )
            ) %>%
                DT::formatStyle(
                    "metacell", "cell_type",
                    backgroundColor = DT::styleEqual(
                        cell_type_colors()$cell_type,
                        col2hex(cell_type_colors()$color)
                    )
                )
        },
        server = FALSE
    )

    output$groupB_table <- DT::renderDataTable(
        {
            req(metacell_types())
            req(cell_type_colors())
            req(groupB())
            DT::datatable(
                tibble(metacell = groupB()) %>%
                    left_join(metacell_types() %>% select(metacell, cell_type), by = "metacell"),
                escape = FALSE,
                rownames = FALSE,
                colnames = "",
                filter = "none",
                options = list(
                    dom = "t",
                    paging = FALSE,
                    language = list(emptyTable = "Please select metacells"),
                    columnDefs = list(list(visible = FALSE, targets = c(1)))
                )
            ) %>%
                DT::formatStyle(
                    "metacell", "cell_type",
                    backgroundColor = DT::styleEqual(
                        cell_type_colors()$cell_type,
                        col2hex(cell_type_colors()$color)
                    )
                )
        },
        server = FALSE
    )

    observeEvent(input$add_metacell_to_groupA, {
        if (is.null(groupA())) {
            groupA(input$metacell)
        } else {
            groupA(unique(c(groupA(), input$metacell)))
        }
    })

    observeEvent(input$add_metacell_to_groupB, {
        if (is.null(groupB())) {
            groupB(input$metacell)
        } else {
            groupB(unique(c(groupB(), input$metacell)))
        }
    })

    observeEvent(input$remove_groupA_metacells, {
        rows <- input$groupA_table_rows_selected
        req(rows)
        req(length(rows) > 0)

        groupA(groupA()[-rows])
    })

    observeEvent(input$remove_groupB_metacells, {
        rows <- input$groupB_table_rows_selected
        req(rows)
        req(length(rows) > 0)

        groupB(groupB()[-rows])
    })

    observeEvent(input$reset_groupA, {
        groupA(NULL)
    })

    observeEvent(input$reset_groupB, {
        groupB(NULL)
    })

    observeEvent(plotly::event_data("plotly_selected", source = "proj_mc_plot"), {
        el <- plotly::event_data("plotly_selected", source = "proj_mc_plot")

        selected_metacells <- unique(el$customdata)
        req(input$mode == "Groups")

        if (input$proj_select_main == "Group A") {
            if (is.null(groupA())) {
                groupA(selected_metacells)
            } else {
                groupA(unique(c(groupA(), selected_metacells)))
            }
        }

        if (input$proj_select_main == "Group B") {
            if (is.null(groupB())) {
                groupB(selected_metacells)
            } else {
                groupB(c(unique(groupB(), selected_metacells)))
            }
        }
    })
}
