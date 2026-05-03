#' samples UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_samples_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_column(
                width = 12,
                generic_box(
                    id = ns("sample_types_box"),
                    title = "Sample types",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_sample_stacked_types"))
                    ),
                    shinydashboardPlus::accordion(
                        id = ns("sample_types_accordion"),
                        shinydashboardPlus::accordionItem(
                            title = "Order by",
                            collapsed = FALSE,
                            shinyWidgets::virtualSelectInput(
                                ns("sample_types_ordering"),
                                "",
                                choices = c(),
                                multiple = FALSE,
                                search = TRUE,
                                dropboxWrapper = "body"
                            )
                        )
                    )
                )
            ),
            generic_column(
                width = 5,
                generic_box(
                    id = ns("sample_sample_box"),
                    title = "Sample/Sample",
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
                    textOutput(ns("please_select_cell_types")),
                    textOutput(ns("no_samples1")),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_gene_mc"))
                    ),
                    shinydashboardPlus::accordion(
                        id = ns("gene_gene_accordion"),
                        shinydashboardPlus::accordionItem(
                            title = "Select axes",
                            collapsed = FALSE,
                            axis_selector("x_axis", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type")),
                            axis_selector("y_axis", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type")),
                            axis_selector("color_by", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type"))
                        )
                    )
                ),
                uiOutput(ns("diff_expr_box")),
                uiOutput(ns("group_comparison_box"))
            ),
            generic_column(
                width = 7,
                projection_box(
                    ns,
                    "sample_projection",
                    title = "Sample projections",
                    color_choices = c("Sample", "Cell type"),
                    additional_elements = textOutput(ns("no_samples2"))
                ),
                uiOutput(ns("sample_info_box"))
            )
        ),
        fluidRow(
            generic_column(
                width = 6,
                uiOutput(ns("sample_qc_box"))
            ),
            generic_column(
                width = 6,
                uiOutput(ns("composition_ci_box"))
            )
        )
    )
}


#' samples sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_samples_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            uiOutput(ns("group_field_ui")),
            uiOutput(ns("cell_type_list")),
            uiOutput(ns("sample_select_ui")),
            uiOutput(ns("group_comparison_select_ui")),
            uiOutput(ns("top_correlated_select_x_axis")),
            uiOutput(ns("top_correlated_select_y_axis")),
            uiOutput(ns("top_correlated_select_color_by"))
        )
    )
}

#' samples Server Function
#'
#' @noRd
mod_samples_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals, state) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            top_correlated_selectors(input, output, session, dataset, metacell_types, ns, button_labels = c("X", "Y", "Color"), gene_modules = gene_modules)

            # --- Group field selector (shown only when cell-level data available) ---
            output$group_field_ui <- renderUI({
                req(dataset())
                if (!has_cell_gene_umis(dataset())) {
                    return(NULL)
                }
                fields <- get_cell_grouping_fields(dataset())
                if (length(fields) == 0) {
                    return(NULL)
                }
                default_field <- get_default_sample_field(dataset()) %||% fields[1]
                selectInput(
                    ns("group_field"),
                    label = "Grouping field",
                    choices = fields,
                    selected = default_field
                )
            })

            # Reactive for the current grouping field with fallback
            group_field <- reactive({
                if (!is.null(input$group_field) && nchar(input$group_field) > 0) {
                    input$group_field
                } else {
                    "samp_id"
                }
            })

            output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", cell_type_colors = cell_type_colors, metacell_types = metacell_types)

            observe({
                choices <- c(dataset_cell_metadata_fields_numeric(dataset()), "Default")
                shinyWidgets::updateVirtualSelect(
                    session = session,
                    inputId = "sample_types_ordering",
                    choices = choices,
                    selected = choices[1]
                )
            })


            # Composition (stacked bars + CI panel)
            setup_samples_composition(input, output, session, ns, dataset, metacell_types, cell_type_colors, group_field, globals, state)

            scatter_selectors(ns, dataset, output, globals, state)
            projection_selectors(ns, dataset, output, input, gene_modules, globals, state, session, weight = 0.6)

            output$sample_select_ui <- renderUI({
                req(dataset())
                req(input$color_proj)
                gf <- group_field()

                # Use cell-level grouping if available, otherwise fall back to samp_id
                if (has_cell_gene_umis(dataset()) && gf != "samp_id") {
                    field_map <- get_cell_field_map(dataset(), gf)
                    samp_list <- sort(unique(field_map))
                } else {
                    samp_list <- get_samples_list(dataset())
                }

                req(samp_list)
                if (length(samp_list) > 1) {
                    selected2 <- samp_list[2]
                } else {
                    selected2 <- samp_list[1]
                }

                picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE)
                tagList(
                    shinyWidgets::pickerInput(
                        ns("samp1"),
                        label = "Sample A:",
                        choices = samp_list,
                        selected = samp_list[1],
                        width = "70%",
                        multiple = FALSE,
                        options = picker_options
                    ),
                    shinyWidgets::pickerInput(
                        ns("samp2"),
                        label = "Sample B:",
                        choices = samp_list,
                        selected = selected2,
                        width = "70%",
                        multiple = FALSE,
                        options = picker_options
                    )
                )
            })

            clipboard_changed <- clipboard_changed_2d_reactive(input, globals, state)

            # Projection plots
            output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, state, source = "proj_mc_plot_gene_tab", tab_guard = "samples") %>%
                bindCache(dataset(), input$color_proj, metacell_types(), cell_type_colors(), input$point_size, input$stroke, input$min_edge_size, input$set_range, input$metacell1, input$metacell2, input$proj_stat, input$expr_range, input$lfp, input$samp1, input$color_proj_gene_module, clipboard_changed(), input$graph_name, input$legend_orientation, input$show_legend_projection, state$manifold_state$mc2d, state$session_ui$plotly_format, state$session_ui$plotly_width, state$session_ui$plotly_height, state$session_ui$plotly_scale)

            # Info box
            output$sample_info_box <- renderUI({
                req(input$samp1)
                generic_box(
                    id = ns("sample_info_box_1"),
                    title = "Sample information",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        DT::dataTableOutput(ns("sample_info_table"))
                    )
                )
            })

            sample_info <- reactive({
                req(input$samp1)
                samp_md <- get_samp_metadata(dataset())
                req(samp_md)
                samp_md %>%
                    filter(samp_id == input$samp1) %>%
                    select(-samp_id) %>%
                    gather("variable", "value")
            })

            output$sample_info_table <- DT::renderDataTable(
                sample_info(),
                escape = FALSE,
                server = FALSE,
                rownames = FALSE,
                caption = paste0("Sample ", input$samp1),
                filter = "none",
                options = list(
                    dom = "t",
                    paging = FALSE,
                    language = list(emptyTable = "Please select metacells")
                )
            )


            # Sample/sample DE + group comparison DE
            setup_samples_de(input, output, session, ns, dataset, metacell_types, cell_type_colors, gene_modules, group_field, globals, state)

            # Sample QC metrics panel
            setup_samples_qc(input, output, session, ns, dataset, group_field, globals, state)

            # Metadata/Metadata plots
            render_axis_select_ui("x_axis", "X axis", "x_axis_select", md_choices = dataset_cell_metadata_fields_numeric(dataset()), md_selected = dataset_cell_metadata_fields_numeric(dataset())[1], selected_gene = mcv_get("default_gene1"), input = input, output = output, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE), gene_modules = gene_modules, session = session)
            render_axis_select_ui("y_axis", "Y axis", "y_axis_select", md_choices = dataset_cell_metadata_fields_numeric(dataset()), md_selected = dataset_cell_metadata_fields_numeric(dataset())[2], selected_gene = mcv_get("default_gene2"), input = input, output = output, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE), gene_modules = gene_modules, session = session)
            render_axis_select_ui("color_by", "Color", "color_by_select", md_choices = c("None", dataset_cell_metadata_fields(dataset())), md_selected = "None", selected_gene = mcv_get("default_gene1"), input = input, output = output, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE), gene_modules = gene_modules, session = session)

            output$please_select_cell_types <- renderPrint({
                if (input$x_axis_type == "Gene" || input$y_axis_type == "Gene" || input$color_by_type == "Gene") {
                    if (is.null(input$selected_cell_types) || length(input$selected_cell_types) == 0) {
                        glue("Please select at least one cell type")
                    }
                } else {
                    req(FALSE)
                }
            })

            for (out in c("no_samples1", "no_samples2")) {
                output[[out]] <- renderPrint({
                    has_samp <- has_samples(dataset())
                    has_cells <- has_cell_gene_umis(dataset()) && length(get_cell_grouping_fields(dataset())) > 0
                    if (!has_samp && !has_cells) {
                        glue("No samples were loaded to MCView.\nPlease make sure your cell metadata has a field called 'samp_id' and run 'import_cell_metadata' again.")
                    } else {
                        req(FALSE)
                    }
                })
            }

            output$plot_gene_gene_mc <- plotly::renderPlotly({
                req(state$tab_state$current_tab == "samples")
                req(input$x_axis_var)
                req(input$y_axis_var)
                req(input$color_by_var)
                req(input$x_axis_type)
                req(input$y_axis_type)
                req(input$color_by_type)
                req(input$gene_gene_point_size)
                req(input$gene_gene_stroke)
                get_samp_metadata(dataset())

                req(axis_vars_ok(dataset(), input, "cell_metadata", gene_modules))

                color_var <- input$color_by_var
                if (input$color_by_var == "Cell type") {
                    color_var <- NULL
                }

                fig <- plot_sample_scatter(
                    dataset(),
                    input$x_axis_var,
                    input$y_axis_var,
                    color_var,
                    x_type = input$x_axis_type,
                    y_type = input$y_axis_type,
                    color_type = input$color_by_type,
                    metacell_types = metacell_types(),
                    cell_type_colors = cell_type_colors(),
                    cell_types = input$selected_cell_types,
                    point_size = input$gene_gene_point_size,
                    stroke = input$gene_gene_stroke,
                    plot_text = FALSE
                ) %>%
                    plotly::ggplotly(tooltip = "tooltip_text", source = "samp_samp_plot") %>%
                    sanitize_for_WebGL() %>%
                    plotly::toWebGL() %>%
                    sanitize_plotly_buttons() %>%
                    sanitize_plotly_download(globals, state)

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
            }) %>% bindCache(dataset(), input$x_axis_var, input$x_axis_type, input$y_axis_var, input$y_axis_type, input$color_by_type, input$color_by_var, metacell_types(), cell_type_colors(), input$gene_gene_point_size, input$gene_gene_stroke, input$selected_cell_types, group_field(), state$session_ui$plotly_format, state$session_ui$plotly_width, state$session_ui$plotly_height, state$session_ui$plotly_scale)

            sample_click_observer("samp_samp_plot", session, "samp1")
            sample_click_observer("samp_types_plot", session, "samp1")
            observeEvent(plotly::event_data("plotly_click", source = "samp_samp_diff_expr_plot"), {
                req(input$x_axis_type == "Gene")
                el <- plotly::event_data("plotly_click", source = "samp_samp_diff_expr_plot")
                selected <- el$customdata
                shinyWidgets::updatePickerInput(session, "x_axis_var", selected = selected)
                showNotification(glue("Selected gene {selected}"))
            })
        }
    )
}
