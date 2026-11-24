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
            generic_column(
                width = 5,
                generic_box(
                    title = "Diff. Expression",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = plot_box_sidebar(
                        startOpen = FALSE,
                        id = ns("mc_mc_sidebar"),
                        div(
                            title = "Hide lateral genes (genes expressed in multiple cell types) from the differential expression analysis",
                            style = "cursor: help;",
                            checkboxInput(ns("hide_lateral"), "Hide lateral genes", value = FALSE)
                        ),
                        div(
                            title = "Hide noisy genes (genes with high variability) from the differential expression analysis",
                            style = "cursor: help;",
                            checkboxInput(ns("hide_noisy"), "Hide noisy genes", value = FALSE)
                        )
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_mc_gene_scatter"))
                    ),
                    div(
                        title = "Toggle display of the differential expression results table showing detailed statistics for each gene",
                        style = "cursor: help;",
                        shinyWidgets::prettySwitch(inputId = ns("show_diff_expr_table"), value = FALSE, label = "Show table")
                    ),
                    DT::DTOutput(ns("diff_expr_table"))
                )
            ),
            generic_column(
                width = 7,
                generic_box(
                    id = ns("metacell_projection"),
                    title = "2D Projection",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = plot_box_sidebar(
                        startOpen = FALSE,
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
                div(
                    title = "Choose comparison mode: 'MCs' compares individual metacells, 'Types' compares cell types, 'Groups' allows custom selection of multiple metacells for comparison",
                    style = "cursor: help;",
                    shinyWidgets::radioGroupButtons(
                        inputId = ns("mode"),
                        label = "Compare:",
                        choices = c(
                            "MCs",
                            "Types",
                            "Groups"
                        ),
                        selected = "Types",
                        justified = TRUE
                    )
                ),
                uiOutput(ns("metacell1_select")),
                uiOutput(ns("metacell2_select")),
                div(
                    title = "Swap the two selected items (metacells, cell types, or groups) to reverse the comparison direction",
                    style = "cursor: help;",
                    shinyWidgets::actionGroupButtons(ns("switch_metacells"), labels = c("Switch"), size = "sm")
                ),
                div(
                    title = "Filter analysis to include only metacells currently in the clipboard (copied from other views or selections)",
                    style = "cursor: help;",
                    checkboxInput(ns("filter_by_clipboard"), "Filter by clipboard", value = FALSE)
                )
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
            # the groups hold metacells of a dataset
            observeEvent(dataset(), ignoreInit = TRUE, {
                groupA(NULL)
                groupB(NULL)
            })

            metacell_names <- metacell_names_reactive(dataset)
            metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)

            projection_selectors(ns, dataset, output, input, gene_modules, globals, session, weight = 0.6)
            group_selectors(input, output, session, dataset, ns, groupA, groupB, metacell_types, cell_type_colors, globals)
            metacell_selectors(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors, groupA, groupB)

            mc_mc_gene_scatter_df <- mc_mc_gene_scatter_df_reactive(dataset, input, output, session, metacell_types, cell_type_colors, globals, groupA, groupB)

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
            output$plot_mc_mc_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, globals, gene_modules, mc_mc_gene_scatter_df, metacell_names, cell_type_colors)

            output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, mc_mc_gene_scatter_df)

            # Projection plots
            # The drag mode of the projection (select in the Groups mode), for its cache key: render_2d_plotly()
            # draws with a copy that an observer of priority 0 sets, which the cached plot doesn't depend on.
            # Set after it, so that the key changes only once the plot would be drawn with the new drag mode.
            proj_group_mode <- reactiveVal(FALSE)
            observe(proj_group_mode(isTRUE(input$mode %in% c("Groups", "Group"))), priority = -1)

            output$plot_mc_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, groupA = groupA, groupB = groupB, source = "proj_mc_plot") %>%
                bindCache(
                    "proj_mc_plot", # the projections of other tabs have the same code, other sources and buttons
                    dataset(),
                    input$color_proj,
                    metacell_types(),
                    cell_type_colors(),
                    input$point_size,
                    input$stroke,
                    input$min_edge_size,
                    input$graph_name,
                    input$legend_orientation,
                    input$show_legend_projection,
                    input$set_range,
                    input$proj_stat,
                    input$expr_range,
                    input$lfp,
                    proj_color_var(input),
                    if (isTRUE(input$color_proj %in% c("Gene module", "Scatter Axis"))) gene_modules(),
                    if (isTRUE(input$color_proj == "Metadata") && isTRUE(input$color_proj_metadata == "Clipboard")) globals$clipboard,
                    # the groups, shown when coloring by "Selected" in the Groups mode
                    if (isTRUE(input$color_proj == "Selected")) list(input$mode, groupA(), groupB()),
                    proj_group_mode(),
                    globals$mc2d,
                    globals$plotly_format,
                    globals$plotly_width,
                    globals$plotly_height,
                    globals$plotly_scale
                )

            diff_expr_auto_update_globals(mc_mc_gene_scatter_df, globals)

            # metacell click observers
            metacell_click_observer("proj_manifold_plot", session)
            metacell_click_observer("md_md_plot", session)
            metacell_click_observer("gene_gene_plot", session)
            metacell_click_observer("proj_metadata_plot", session)
            metacell_click_observer("proj_mc_plot_gene_tab", session)
            metacell_click_observer("gene_time_mc_plot1", session)
            metacell_click_observer("gene_time_mc_plot2", session)

            # Output priorities
            # the comparison needs the two pickers, also when another tab started this module (globals$load_tab)
            # while the tab is hidden
            outputOptions(output, "metacell1_select", suspendWhenHidden = FALSE)
            outputOptions(output, "metacell2_select", suspendWhenHidden = FALSE)
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

# The pickers are rendered again when the annotation changes (their colors and the list of cell types):
# keep what the user chose, when it is still one of the choices
keep_selected <- function(current, choices, default) {
    if (length(current) == 1 && !is.na(current) && current %in% choices) current else default
}

# Whether there are metacells to paste from the clipboard, telling the user when there aren't
clipboard_has_metacells <- function(globals) {
    if (length(globals$clipboard) == 0) {
        showNotification("The clipboard is empty: add metacells to it first (with the \"Clipboard\" button at the top, or by copying the metacells selected in another tab)", type = "warning")
        return(FALSE)
    }
    TRUE
}

# The choices of a metacell picker, with "(dissimilar)" after the metacells that the projection on an atlas
# marks as dissimilar. Matched by metacell: the metadata can be in another order, or miss metacells.
metacell_choices_with_similarity <- function(metacells, md) {
    choices <- metacells
    if (!is.null(md) && has_name(md, "similar")) {
        dissimilar <- md$similar[match(metacells, as.character(md$metacell))] %in% "dissimilar"
        names(choices) <- ifelse(dissimilar, paste0(metacells, " (dissimilar)"), metacells)
    }
    choices
}

metacell_selectors <- function(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors, groupA, groupB) {
    output$metacell1_select <- renderUI({
        req(dataset())
        req(input$mode)
        if (input$mode == "MCs") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())
            # add 'similar' annotation
            choices <- metacell_choices_with_similarity(metacell_names(), get_mc_data(dataset(), "metadata"))
            shinyWidgets::pickerInput(ns("metacell1"), "Metacell A",
                choices = choices,
                selected = keep_selected(isolate(input$metacell1), metacell_names(), config$selected_mc1 %||% metacell_names()[1]),
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
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
                selected = keep_selected(isolate(input$metacell1), cell_types, cell_types[1]),
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
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
                    selected = keep_selected(isolate(input$metacell), metacell_names(), config$selected_mc1 %||% metacell_names()[1]), multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
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
                selected = keep_selected(isolate(input$metacell2), metacell_names(), config$selected_mc2 %||% metacell_names()[2]), multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
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
                selected = keep_selected(isolate(input$metacell2), cell_types, cell_types[2]),
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        }
    })

    # Select metacell / cell type when clicking on it
    observeEvent(plotly::event_data("plotly_click", source = "proj_mc_plot"), {
        el <- plotly::event_data("plotly_click", source = "proj_mc_plot")
        metacell <- plotly_event_ids(el)
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

group_selectors <- function(input, output, session, dataset, ns, groupA, groupB, metacell_types, cell_type_colors, globals) {
    output$groupA_box <- renderUI({
        req(input$mode == "Groups")
        generic_box(
            id = ns("groupA_box_1"),
            title = "Group A metacells",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            actionButton(ns("reset_groupA"), "Reset"),
            actionButton(ns("remove_groupA_metacells"), "Remove"),
            actionButton(ns("paste_groupA_metacells"), "Paste"),
            shinycssloaders::withSpinner(
                DT::dataTableOutput(ns("groupA_table"))
            )
        )
    })

    output$groupB_box <- renderUI({
        req(input$mode == "Groups")
        generic_box(
            id = ns("groupB_box_1"),
            title = "Group B metacells",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            actionButton(ns("reset_groupB"), "Reset"),
            actionButton(ns("remove_groupB_metacells"), "Remove"),
            actionButton(ns("paste_groupB_metacells"), "Paste"),
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

    observeEvent(input$paste_groupA_metacells, {
        req(clipboard_has_metacells(globals))
        groupA(unique(c(groupA(), globals$clipboard)))
    })

    observeEvent(input$paste_groupB_metacells, {
        req(clipboard_has_metacells(globals))
        groupB(unique(c(groupB(), globals$clipboard)))
    })

    observeEvent(input$reset_groupA, {
        groupA(NULL)
    })

    observeEvent(input$reset_groupB, {
        groupB(NULL)
    })

    observeEvent(plotly::event_data("plotly_selected", source = "proj_mc_plot"), {
        el <- plotly::event_data("plotly_selected", source = "proj_mc_plot")

        selected_metacells <- unique(plotly_event_ids(el))
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
                groupB(unique(c(groupB(), selected_metacells)))
            }
        }
    })
}
