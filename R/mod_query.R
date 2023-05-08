#' projection UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_query_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_column(
                width = 7,
                projection_box(
                    ns,
                    "metacell_projection",
                    uiOutput(ns("query_metadata_selector")),
                    uiOutput(ns("atlas_metadata_selector")),
                    title = "Query 2D Projection",
                    color_choices = c("Cell type", "Similarity", "Query Metadata", "Atlas Metadata", "Gene", "Gene module", "Selected")
                )
            ),
            generic_column(
                width = 5,
                generic_box(
                    title = "Diff. Expression",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        id = ns("diff_expr_sidebar"),
                        startOpen = FALSE,
                        width = 60,
                        checkboxInput(ns("mark_disjoined"), "Mark disjoined genes", value = TRUE)
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_mc_gene_scatter"))
                    ),
                    shinyWidgets::prettySwitch(inputId = ns("show_diff_expr_table"), value = FALSE, label = "Show table"),
                    DT::DTOutput(ns("diff_expr_table"))
                )
            )
        ),
        fluidRow(
            generic_column(
                width = 7,
                generic_box(
                    id = ns("metacell_projection"),
                    title = "Type predictions",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotOutput(ns("plot_mc_stacked_type"))
                    )
                )
            ),
            generic_column(
                width = 5,
                generic_box(
                    id = ns("scatter_box"),
                    title = "Gene/Gene",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 100,
                        id = ns("gene_gene_sidebar"),
                        uiOutput(ns("gene_gene_point_size_ui")),
                        uiOutput(ns("gene_gene_stroke_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_gene_mc"))
                    ),
                    shinydashboardPlus::accordion(
                        id = ns("gene_gene_accordion"),
                        shinydashboardPlus::accordionItem(
                            title = "Select Gene and colors",
                            axis_selector("axis", "Gene", ns, choices = c("Gene")),
                            axis_selector("color_by", "Metadata", ns, choices = c("Metadata", "Gene")),
                        )
                    )
                )
            ),
            column(
                width = 3,
                uiOutput(ns("group_box"))
            )
        ),
        fluidRow(
            generic_column(
                width = 12,
                generic_box(
                    id = ns("gene_metadata_box"),
                    title = "Gene metadata",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    uiOutput(ns("gene_metadata_cell_type_selector")),
                    shinycssloaders::withSpinner(
                        DT::dataTableOutput(ns("gene_metadata_table"))
                    )
                )
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
mod_query_sidebar_ui <- function(id) {
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
                    selected = "Type",
                    justified = TRUE
                ),
                uiOutput(ns("diff_select"))
            ),
            uiOutput(ns("top_correlated_select_axis"))
        )
    )
}

#' projection Server Function
#'
#' @noRd
mod_query_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
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

            metacell_names <- metacell_names_reactive(dataset)
            metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)

            mod_query_globals_observers(input, session, globals, dataset)

            picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE)

            output$query_metadata_selector <- renderUI({
                if (!has_metadata(dataset())) {
                    print(glue("Query doesn't have any metadata."))
                } else {
                    shinyWidgets::pickerInput(
                        ns("color_proj_query_metadata"),
                        label = "Metadata:",
                        choices = dataset_metadata_fields(dataset()),
                        selected = dataset_metadata_fields(dataset())[1],
                        width = "70%",
                        multiple = FALSE,
                        options = picker_options
                    )
                }
            })

            output$atlas_metadata_selector <- renderUI({
                if (!has_metadata(dataset(), atlas = TRUE)) {
                    print(glue("Atlas doesn't have any metadata."))
                } else {
                    shinyWidgets::pickerInput(
                        ns("color_proj_atlas_metadata"),
                        label = "Metadata:",
                        choices = dataset_metadata_fields_numeric(dataset(), atlas = TRUE),
                        selected = dataset_metadata_fields_numeric(dataset(), atlas = TRUE)[1],
                        width = "70%",
                        multiple = FALSE,
                        options = picker_options
                    )
                }
            })

            observe({
                req(input$color_proj)
                shinyjs::toggle(id = "query_metadata_selector", condition = input$color_proj == "Query Metadata")
                shinyjs::toggle(id = "atlas_metadata_selector", condition = input$color_proj == "Atlas Metadata")
            })


            scatter_selectors(ns, dataset, output, globals)
            projection_selectors(ns, dataset, output, input, gene_modules, globals, session, weight = 0.6)
            top_correlated_selector("axis_var", "axis", "axis_type", input, output, session, dataset, ns, button_labels = c("Axes", "Color"), ids = c("axis", "color"))

            group_selectors_mod_query(input, output, session, dataset, ns, group, metacell_types, cell_type_colors, globals)
            metacell_selectors_mod_query(input, output, session, dataset, ns, metacell_names, metacell_colors, projected_metacell_types, atlas_colors, group)

            mc_mc_gene_scatter_df <- reactive({
                req(input$mode)
                req(!is.null(input$mark_disjoined))

                if (input$mode == "MC") {
                    req(input$metacell1)
                    req(input$metacell1 %in% metacell_names())
                    df <- calc_obs_exp_mc_df(dataset(), input$metacell1)
                } else if (input$mode == "Type") {
                    req(input$metacell1)
                    req(input$metacell1 %in% atlas_colors()$cell_type)
                    req(input$metacell1 %in% metacell_types()$cell_type) # we cannot show diff expression if the cell type doesn't exist in the query
                    req(projected_metacell_types())
                    df <- calc_obs_exp_type_df(dataset(), input$metacell1, projected_metacell_types())
                } else if (input$mode == "Group") {
                    req(group())
                    req(length(group()) > 1)
                    group_types_df <- tibble(metacell = group(), cell_type = "Group")
                    df <- calc_obs_exp_type_df(dataset(), "Group", group_types_df)
                }

                disjoined <- get_mc_data(dataset(), "disjoined_genes_no_atlas")

                if (input$mark_disjoined) {
                    prev_levels <- levels(df$col)
                    df <- df %>%
                        mutate(col = ifelse(gene %in% disjoined, "#00b7ff", as.character(col))) %>%
                        mutate(col = factor(col, levels = c("#00b7ff", "purple", prev_levels)))
                }

                df <- df %>%
                    mutate(Disjoined = gene %in% disjoined)

                return(df)
            })

            # Projection plots
            output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, projected_metacell_types, atlas_colors, gene_modules, globals, group = group, source = "proj_mc_plot_proj_tab")

            # connect_gene_plots(input, output, session, ns, source = "proj_mc_plot_proj_tab")

            # Differential expression
            output$plot_mc_mc_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, gene_modules, mc_mc_gene_scatter_df, metacell_names(), atlas_colors())

            # Select a gene when clicking on it
            plotly_click_observer("projection_diff_expr_plot", session, "axis_var", notification_prefix = "Selected ", update_function = shinyWidgets::updatePickerInput)

            output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, mc_mc_gene_scatter_df)

            # Scatter
            render_axis_select_ui("axis", "Data", "axis_select", md_choices = dataset_metadata_fields_numeric(dataset(), atlas = TRUE), md_selected = dataset_metadata_fields_numeric(dataset(), atlas = TRUE)[1], selected_gene = default_gene1, input = input, output = output, ns = ns, dataset = dataset, gene_modules = gene_modules, session = session)

            render_axis_select_ui(
                "color_by",
                "Color",
                "color_by_select",
                md_choices = c(
                    "Cell type",
                    paste0(dataset_metadata_fields_numeric(dataset(), atlas = TRUE), "_atlas"),
                    dataset_metadata_fields(dataset(), atlas = FALSE)
                ),
                md_selected = "Cell type",
                selected_gene = default_gene1,
                input = input,
                output = output,
                ns = ns,
                dataset = dataset,
                gene_modules = gene_modules,
                session = session
            )

            output$plot_gene_gene_mc <- plotly::renderPlotly({
                req(input$axis_var)
                req(input$color_by_var)
                req(input$axis_type)
                req(input$color_by_type)
                req(input$gene_gene_point_size)
                req(input$gene_gene_stroke)
                req(axis_vars_ok(dataset(), input, "metadata", gene_modules, axes = c("axis", "color_by"), atlas = TRUE) || axis_vars_ok(dataset(), input, "metadata", axes = c("axis", "color_by"), atlas = FALSE))

                color_var <- input$color_by_var
                if (input$color_by_var == "Cell type") {
                    color_var <- NULL
                }

                fig <- plot_obs_proj_scatter(
                    dataset(),
                    input$axis_var,
                    color_var,
                    axis_type = input$axis_type,
                    color_type = input$color_by_type,
                    metacell_types = metacell_types(),
                    cell_type_colors = cell_type_colors(),
                    point_size = input$gene_gene_point_size,
                    stroke = input$gene_gene_stroke,
                    plot_text = FALSE
                ) %>%
                    plotly::ggplotly(tooltip = "tooltip_text", source = "obs_proj_plot") %>%
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

                if (!is.null(input$mode) && input$mode %in% c("Groups", "Group")) {
                    fig <- fig %>% plotly::layout(dragmode = "select")
                }


                return(fig)
            }) %>% bindCache(dataset(), input$axis_var, input$axis_type, input$color_by_type, input$color_by_var, metacell_types(), cell_type_colors(), input$gene_gene_point_size, input$gene_gene_stroke, input$mode)

            # Point size selector
            output$point_size_ui <- renderUI({
                numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset(), globals$screen_width, globals$screen_height, weight = 0.6), min = 0.1, max = 3, step = 0.1)
            })

            # Minimal edge length selector
            output$edge_distance_ui <- renderUI({
                sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
            })

            output$plot_mc_stacked_type <- plot_type_predictions_bar(dataset, metacell_types, cell_type_colors)

            output$gene_metadata_cell_type_selector <- cell_type_selector(dataset, ns, id = "gene_metadata_cell_type", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors)

            current_gene_table <- reactiveVal()

            observe({
                req(input$gene_metadata_cell_type)
                gene_metadata <- get_mc_data(dataset(), "gene_metadata")
                req(gene_metadata)
                current_gene_table(gene_metadata %>%
                    filter(cell_type %in% input$gene_metadata_cell_type))
            })

            observe({
                req(input$gene_metadata_table_rows_selected)
                gene <- current_gene_table() %>%
                    slice(input$gene_metadata_table_rows_selected) %>%
                    pull(gene)
                showNotification(glue("selecting {gene}"))
                shinyWidgets::updatePickerInput(session, "axis_var", selected = gene)
            })

            output$gene_metadata_table <- DT::renderDataTable(
                {
                    req(current_gene_table())
                    dt <- DT::datatable(
                        current_gene_table(),
                        escape = FALSE,
                        rownames = FALSE,
                        extensions = c("Scroller", "FixedColumns"),
                        selection = "single",
                        filter = "top",
                        options = list(
                            dom = "Bfrtip",
                            deferRender = TRUE,
                            scrollY = 200,
                            scroller = TRUE,
                            scrollX = TRUE,
                            fixedColumns = list(leftColumns = 2)
                        )
                    )

                    round_columns <- c("correction_factor")
                    round_columns <- round_columns[round_columns %in% colnames(current_gene_table())]
                    if (length(round_columns) > 0) {
                        dt <- dt %>% DT::formatRound(columns = round_columns, digits = 2)
                    }

                    return(dt)
                },
                server = TRUE
            )
        }
    )
}

mod_query_globals_observers <- function(input, session, globals, dataset, notification_suffix = " in \"Query\" tab") {
    observe({
        req(globals$selected_query_gene)
        req(input$axis_type == "Gene")
        req(input$axis_var)
        shinyWidgets::updatePickerInput(session, "axis_var", selected = globals$selected_query_gene)

        showNotification(glue("Selected {globals$selected_query_gene}{notification_suffix}"))
        globals$selected_query_gene <- NULL
    })

    observe({
        req(input$mode == "MC")
        req(input$metacell1)
        req(globals$selected_query_metacell)
        shinyWidgets::updatePickerInput(session, "metacell1", selected = globals$selected_query_metacell)
        showNotification(glue("Selected {globals$selected_query_metacell}{notification_suffix}"))
        globals$selected_query_metacell <- NULL
    })
}


metacell_selectors_mod_query <- function(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors, group) {
    output$diff_select <- renderUI({
        req(dataset())
        req(input$mode)
        if (input$mode == "MC") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())
            shinyWidgets::pickerInput(ns("metacell1"), "Metacell",
                choices = metacell_names(),
                selected = config$selected_mc1 %||% metacell_names()[1], multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains"),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        } else if (input$mode == "Type") {
            req(cell_type_colors())
            req(metacell_types())
            # do not show cell types that do not exist in the query
            types_df <- cell_type_colors() %>% filter(cell_type %in% metacell_types()$cell_type)
            cell_types_hex <- col2hex(types_df$color)
            cell_types <- types_df$cell_type
            shinyWidgets::pickerInput(ns("metacell1"), "Cell type",
                choices = cell_types,
                selected = cell_types[1],
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        } else if (input$mode == "Group") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())
            tagList(
                shinyWidgets::pickerInput(ns("metacell"), "Metacell",
                    choices = metacell_names(),
                    selected = config$selected_mc1 %||% metacell_names()[1], multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                    choicesOpt = list(
                        style = paste0("color: ", cell_types_hex, ";")
                    )
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
    select_metacell_plotly_event_projection("proj_mc_plot_proj_tab", input, session, metacell_types, group)
    select_metacell_plotly_event_projection("obs_proj_plot", input, session, metacell_types, group)
    select_metacell_plotly_event_projection("type_prediction_bar", input, session, metacell_types, group)
}

select_metacell_plotly_event_projection <- function(source, input, session, metacell_types, group) {
    # Select metacell / cell type when clicking on it
    observeEvent(plotly::event_data("plotly_click", source = source), {
        el <- plotly::event_data("plotly_click", source = source)
        metacell <- el$customdata

        if (input$mode == "MC") {
            updateSelectInput(session, "metacell1", selected = metacell)
            showNotification(glue("Selected Metacell: #{metacell}"))
        } else if (input$mode == "Type") {
            cell_type <- metacell_types() %>%
                filter(metacell == !!metacell) %>%
                slice(1) %>%
                pull(cell_type)
            req(cell_type)
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

group_selectors_mod_query <- function(input, output, session, dataset, ns, group, metacell_types, cell_type_colors, globals) {
    output$group_box <- renderUI({
        req(input$mode == "Group")
        generic_box(
            id = ns("group_box_1"),
            title = "Group metacells",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            actionButton(ns("reset_group"), "Reset"),
            actionButton(ns("remove_group_metacells"), "Remove"),
            actionButton(ns("paste_group_metacells"), "Paste"),
            shinycssloaders::withSpinner(
                DT::dataTableOutput(ns("group_table"))
            )
        )
    })

    output$group_table <- DT::renderDataTable(
        {
            req(metacell_types())
            req(cell_type_colors())
            req(group())
            DT::datatable(
                tibble(metacell = group()) %>%
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

    observeEvent(input$remove_group_metacells, {
        rows <- input$group_table_rows_selected
        req(rows)
        req(length(rows) > 0)

        group(group()[-rows])
    })

    observeEvent(input$paste_group_metacells, {
        metacells <- globals$clipboard
        group(unique(c(group(), metacells)))
    })

    observeEvent(plotly::event_data("plotly_selected", source = "proj_mc_plot_proj_tab"), {
        el <- plotly::event_data("plotly_selected", source = "proj_mc_plot_proj_tab")

        selected_metacells <- unique(el$customdata)
        req(input$mode == "Group")


        if (is.null(group())) {
            group(selected_metacells)
        } else {
            group(unique(c(group(), selected_metacells)))
        }
    })

    observeEvent(plotly::event_data("plotly_selected", source = "obs_proj_plot"), {
        el <- plotly::event_data("plotly_selected", source = "obs_proj_plot")

        selected_metacells <- unique(el$customdata)
        req(input$mode == "Group")


        if (is.null(group())) {
            group(selected_metacells)
        } else {
            group(unique(c(group(), selected_metacells)))
        }
    })
}
