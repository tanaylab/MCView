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
                )
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
                )
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
mod_mc_mc_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            groupA <- reactiveVal()
            groupB <- reactiveVal()

            metacell_names <- metacell_names_reactive(dataset)
            metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)

            projection_selectors(ns, dataset, output, input, gene_modules, globals, weight = 0.6)
            group_selectors(input, output, session, dataset, ns, groupA, groupB, metacell_types, cell_type_colors)
            metacell_selectors(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors, groupA, groupB)

            mc_mc_gene_scatter_df <- mc_mc_gene_scatter_df_reactive(dataset, input, output, session, metacell_types, cell_type_colors, groupA, groupB)

            diff_expr_switch_metacells(dataset, input, output, session, groupA, groupB)
            mod_mc_mc_globals_observers(input, session, globals)

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

            # Projection plots
            output$plot_mc_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, groupA = groupA, groupB = groupB, source = "proj_mc_plot")


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
        }
    )
}

mod_mc_mc_globals_observers <- function(input, session, globals, notification_suffix = " in \"Diff. Expr\" tab") {
    observe({
        req(input$mode == "MCs")
        req(globals$selected_metacellA)
        req(input$metacell1)
        shinyWidgets::updatePickerInput(session, "metacell1", selected = globals$selected_metacellA)
        showNotification(glue("Selected {globals$selected_metacellA}{notification_suffix}"))
        globals$selected_metacellA <- NULL
    })

    observe({
        req(input$mode == "MCs")
        req(globals$selected_metacellB)
        req(input$metacell2)
        shinyWidgets::updatePickerInput(session, "metacell2", selected = globals$selected_metacellB)
        showNotification(glue("Selected {globals$selected_metacellB}{notification_suffix}"))
        globals$selected_metacellB <- NULL
    })
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
                selected = config$selected_mc1 %||% metacell_names()[1], multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
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
                    selected = config$selected_mc1 %||% metacell_names()[1], multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
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
                selected = config$selected_mc2 %||% metacell_names()[2], multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
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
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE),
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
