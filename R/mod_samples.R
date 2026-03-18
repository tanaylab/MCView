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
mod_samples_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
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

            output$plot_sample_stacked_types <- plot_sample_stacked_types(dataset, globals, metacell_types, cell_type_colors, input, group_field)

            scatter_selectors(ns, dataset, output, globals)
            projection_selectors(ns, dataset, output, input, gene_modules, globals, session, weight = 0.6)

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

            clipboard_changed <- clipboard_changed_2d_reactive(input, globals)

            # Projection plots
            output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, source = "proj_mc_plot_gene_tab", tab_guard = "samples") %>%
                bindCache(dataset(), input$color_proj, metacell_types(), cell_type_colors(), input$point_size, input$stroke, input$min_edge_size, input$set_range, input$metacell1, input$metacell2, input$proj_stat, input$expr_range, input$lfp, input$samp1, input$color_proj_gene_module, clipboard_changed(), input$graph_name, input$legend_orientation, input$show_legend_projection, globals$mc2d, globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)

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

            # Differential expression
            output$diff_expr_box <- renderUI({
                req(input$selected_cell_types)

                generic_box(
                    title = "Diff. Expression",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_samp_samp_gene_scatter"))
                    ),
                    shinyWidgets::prettySwitch(inputId = ns("show_diff_expr_table"), value = FALSE, label = "Show table"),
                    DT::DTOutput(ns("diff_expr_table"))
                )
            })

            samp_samp_scatter_df <- reactive({
                req(input$selected_cell_types)
                req(input$samp1)
                req(input$samp2)
                gf <- group_field()

                if (has_cell_gene_umis(dataset()) && gf != "samp_id") {
                    # Cell-level pseudobulk DE
                    calc_samp_samp_gene_df(dataset(), input$samp1, input$samp2, metacell_types(),
                        cell_types = input$selected_cell_types,
                        group_field = gf
                    )
                } else {
                    # Original metacell-weighted approach
                    samp_frac <- get_samp_mc_frac(dataset())
                    req(input$samp1 %in% rownames(samp_frac))
                    req(input$samp2 %in% rownames(samp_frac))
                    req(sum(samp_frac[input$samp1, ], na.rm = TRUE) > 0)
                    req(sum(samp_frac[input$samp2, ], na.rm = TRUE) > 0)
                    calc_samp_samp_gene_df(dataset(), input$samp1, input$samp2, metacell_types(),
                        cell_types = input$selected_cell_types
                    )
                }
            }) %>% bindCache(dataset(), input$selected_cell_types, input$samp1, input$samp2, metacell_types(), group_field())

            output$plot_samp_samp_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, globals, gene_modules, samp_samp_scatter_df, metacell_names, cell_type_colors, mode = "Samples", tab_guard = "samples")

            output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, samp_samp_scatter_df)

            # --- Group Comparison (multi-select A vs B) ---
            output$group_comparison_select_ui <- renderUI({
                req(dataset())
                if (!has_cell_gene_umis(dataset())) {
                    return(NULL)
                }
                gf <- group_field()
                field_map <- get_cell_field_map(dataset(), gf)
                req(field_map)
                group_values <- sort(unique(field_map))
                req(length(group_values) >= 2)

                picker_options <- shinyWidgets::pickerOptions(
                    liveSearch = TRUE,
                    liveSearchNormalize = TRUE,
                    liveSearchStyle = "contains",
                    dropupAuto = FALSE,
                    actionsBox = TRUE
                )
                tagList(
                    tags$hr(),
                    tags$strong("Group Comparison"),
                    shinyWidgets::pickerInput(
                        ns("group_a"),
                        label = "Group A:",
                        choices = group_values,
                        selected = NULL,
                        width = "100%",
                        multiple = TRUE,
                        options = picker_options
                    ),
                    shinyWidgets::pickerInput(
                        ns("group_b"),
                        label = "Group B:",
                        choices = group_values,
                        selected = NULL,
                        width = "100%",
                        multiple = TRUE,
                        options = picker_options
                    )
                )
            })

            output$group_comparison_box <- renderUI({
                req(dataset())
                req(has_cell_gene_umis(dataset()))
                req(input$group_a)
                req(input$group_b)
                req(input$selected_cell_types)

                generic_box(
                    id = ns("group_comparison_box_inner"),
                    title = "Group Comparison",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_group_comparison_scatter"))
                    ),
                    shinyWidgets::prettySwitch(inputId = ns("show_group_comparison_table"), value = FALSE, label = "Show table"),
                    DT::DTOutput(ns("group_comparison_table"))
                )
            })

            group_comparison_df <- reactive({
                req(input$group_a)
                req(input$group_b)
                req(input$selected_cell_types)
                req(has_cell_gene_umis(dataset()))
                gf <- group_field()

                shiny::withProgress(message = "Computing group comparison DE...", {
                    calc_group_diff_expr(
                        dataset(),
                        group_field = gf,
                        group1_values = input$group_a,
                        group2_values = input$group_b,
                        cell_types = input$selected_cell_types,
                        metacell_types = metacell_types()
                    )
                })
            }) %>% bindCache(dataset(), input$group_a, input$group_b, metacell_types(), group_field(), input$selected_cell_types)

            output$plot_group_comparison_scatter <- render_mc_mc_gene_plotly(
                input, output, session, ns, dataset, globals, gene_modules,
                group_comparison_df, metacell_names = NULL, cell_type_colors = NULL,
                mode = "GroupComparison", tab_guard = "samples"
            )

            output$group_comparison_table <- render_group_comparison_diff_table(input, output, session, ns, dataset, group_comparison_df)

            observeEvent(plotly::event_data("plotly_click", source = "grp_comparison_plot"), {
                req(input$x_axis_type == "Gene")
                el <- plotly::event_data("plotly_click", source = "grp_comparison_plot")
                selected <- el$customdata
                shinyWidgets::updatePickerInput(session, "x_axis_var", selected = selected)
                showNotification(glue("Selected gene {selected}"))
            })

            # --- Sample QC Metrics Panel ---
            qc_stats <- reactive({
                req(dataset())
                req(has_cell_gene_umis(dataset()))
                gf <- group_field()
                get_group_qc_stats(dataset(), gf)
            }) %>% bindCache(dataset(), group_field())

            output$sample_qc_box <- renderUI({
                req(dataset())
                if (!has_cell_gene_umis(dataset())) {
                    return(NULL)
                }

                generic_box(
                    id = ns("sample_qc_box_inner"),
                    title = "Sample QC",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_qc_cells_per_group"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_qc_umis_per_group"))
                    ),
                    shinyWidgets::prettySwitch(inputId = ns("show_qc_table"), value = FALSE, label = "Show QC table"),
                    DT::DTOutput(ns("qc_summary_table"))
                )
            })

            output$plot_qc_cells_per_group <- plotly::renderPlotly({
                req(globals$current_tab == "samples")
                stats_df <- qc_stats()
                req(nrow(stats_df) > 0)

                # Sort by cell count descending
                stats_df <- stats_df %>%
                    dplyr::arrange(dplyr::desc(n_cells)) %>%
                    dplyr::mutate(group_id = factor(group_id, levels = group_id))

                # Highlight selected samples
                samp1 <- input$samp1
                samp2 <- input$samp2
                stats_df <- stats_df %>%
                    dplyr::mutate(
                        highlight = dplyr::case_when(
                            group_id == samp1 ~ "Sample A",
                            group_id == samp2 ~ "Sample B",
                            TRUE ~ "Other"
                        )
                    )

                color_map <- c("Sample A" = "#E41A1C", "Sample B" = "#377EB8", "Other" = "#999999")

                fig <- plotly::plot_ly(
                    stats_df,
                    x = ~group_id,
                    y = ~n_cells,
                    color = ~highlight,
                    colors = color_map,
                    type = "bar",
                    hoverinfo = "text",
                    text = ~paste0(group_id, "<br>Cells: ", scales::comma(n_cells))
                ) %>%
                    plotly::layout(
                        title = list(text = "Cells per Group", font = list(size = 14)),
                        xaxis = list(title = group_field(), tickangle = -45, showticklabels = nrow(stats_df) <= 50),
                        yaxis = list(title = "# Cells"),
                        showlegend = TRUE,
                        legend = list(orientation = "h", yanchor = "bottom", y = 1.02, xanchor = "right", x = 1)
                    ) %>%
                    sanitize_plotly_buttons() %>%
                    sanitize_plotly_download(globals)

                fig
            }) %>% bindCache(dataset(), group_field(), input$samp1, input$samp2, globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)

            output$plot_qc_umis_per_group <- plotly::renderPlotly({
                req(globals$current_tab == "samples")
                stats_df <- qc_stats()
                req(nrow(stats_df) > 0)
                req(!all(is.na(stats_df$median_umis_per_cell)))

                # Sort by median UMIs descending
                stats_df <- stats_df %>%
                    dplyr::filter(!is.na(median_umis_per_cell)) %>%
                    dplyr::arrange(dplyr::desc(median_umis_per_cell)) %>%
                    dplyr::mutate(group_id = factor(group_id, levels = group_id))

                # Highlight selected samples
                samp1 <- input$samp1
                samp2 <- input$samp2
                stats_df <- stats_df %>%
                    dplyr::mutate(
                        highlight = dplyr::case_when(
                            group_id == samp1 ~ "Sample A",
                            group_id == samp2 ~ "Sample B",
                            TRUE ~ "Other"
                        )
                    )

                color_map <- c("Sample A" = "#E41A1C", "Sample B" = "#377EB8", "Other" = "#999999")

                fig <- plotly::plot_ly(
                    stats_df,
                    x = ~group_id,
                    y = ~median_umis_per_cell,
                    color = ~highlight,
                    colors = color_map,
                    type = "bar",
                    hoverinfo = "text",
                    text = ~paste0(group_id, "<br>Median UMIs/cell: ", scales::comma(round(median_umis_per_cell)))
                ) %>%
                    plotly::layout(
                        title = list(text = "Median UMIs per Cell", font = list(size = 14)),
                        xaxis = list(title = group_field(), tickangle = -45, showticklabels = nrow(stats_df) <= 50),
                        yaxis = list(title = "Median UMIs/cell"),
                        showlegend = TRUE,
                        legend = list(orientation = "h", yanchor = "bottom", y = 1.02, xanchor = "right", x = 1)
                    ) %>%
                    sanitize_plotly_buttons() %>%
                    sanitize_plotly_download(globals)

                fig
            }) %>% bindCache(dataset(), group_field(), input$samp1, input$samp2, globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)

            output$qc_summary_table <- DT::renderDT(
                {
                    stats_df <- qc_stats()
                    req(nrow(stats_df) > 0)
                    if (!is.null(input$show_qc_table) && input$show_qc_table) {
                        display_df <- stats_df %>%
                            dplyr::rename(
                                `Group` = group_id,
                                `# Cells` = n_cells,
                                `Total UMIs` = total_umis,
                                `Median UMIs/cell` = median_umis_per_cell
                            ) %>%
                            dplyr::arrange(dplyr::desc(`# Cells`))

                        DT::datatable(
                            display_df,
                            selection = "single",
                            escape = FALSE,
                            rownames = FALSE,
                            options = list(
                                dom = "frtip",
                                pageLength = 15,
                                scrollX = TRUE,
                                language = list(emptyTable = "No QC data available")
                            )
                        ) %>%
                            DT::formatRound(c("Total UMIs", "Median UMIs/cell"), digits = 0)
                    }
                },
                server = FALSE
            )

            # --- Composition Confidence Intervals Panel ---
            composition_ci_data <- reactive({
                req(dataset())
                req(has_cell_gene_umis(dataset()))
                req(input$samp1)
                req(input$samp2)
                gf <- group_field()

                composition <- get_group_cell_type_composition(dataset(), gf, metacell_types())
                req(nrow(composition) > 0)

                # Filter to only selected samples
                selected <- c(input$samp1, input$samp2)
                comp_selected <- composition %>%
                    dplyr::filter(group_id %in% selected)

                req(nrow(comp_selected) > 0)
                comp_selected
            }) %>% bindCache(dataset(), group_field(), metacell_types(), input$samp1, input$samp2)

            output$composition_ci_box <- renderUI({
                req(dataset())
                if (!has_cell_gene_umis(dataset())) {
                    return(NULL)
                }

                generic_box(
                    id = ns("composition_ci_box_inner"),
                    title = "Composition with Confidence Intervals",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_composition_ci"))
                    )
                )
            })

            output$plot_composition_ci <- plotly::renderPlotly({
                req(globals$current_tab == "samples")
                comp_df <- composition_ci_data()
                req(nrow(comp_df) > 0)

                # Get cell type colors for consistent coloring
                ct_colors <- cell_type_colors() %>%
                    dplyr::select(cell_type, color) %>%
                    tibble::deframe()

                # Order cell types by maximum fraction across both samples
                ct_order <- comp_df %>%
                    dplyr::group_by(cell_type) %>%
                    dplyr::summarise(max_frac = max(fraction), .groups = "drop") %>%
                    dplyr::arrange(dplyr::desc(max_frac)) %>%
                    dplyr::pull(cell_type)

                comp_df <- comp_df %>%
                    dplyr::mutate(cell_type = factor(cell_type, levels = rev(ct_order)))

                samp1 <- input$samp1
                samp2 <- input$samp2

                # Build side-by-side grouped bar chart with error bars
                fig <- plotly::plot_ly() %>%
                    plotly::add_trace(
                        data = comp_df %>% dplyr::filter(group_id == samp1),
                        x = ~fraction,
                        y = ~cell_type,
                        type = "bar",
                        orientation = "h",
                        name = samp1,
                        marker = list(color = "#E41A1C"),
                        error_x = list(
                            type = "data",
                            symmetric = FALSE,
                            array = ~(fraction_upper - fraction),
                            arrayminus = ~(fraction - fraction_lower),
                            color = "#333333",
                            thickness = 1.5
                        ),
                        hoverinfo = "text",
                        text = ~paste0(
                            cell_type, " (", samp1, ")",
                            "<br>Fraction: ", round(fraction, 3),
                            "<br>95% CI: [", round(fraction_lower, 3), ", ", round(fraction_upper, 3), "]",
                            "<br>Cells: ", n_cells
                        )
                    ) %>%
                    plotly::add_trace(
                        data = comp_df %>% dplyr::filter(group_id == samp2),
                        x = ~fraction,
                        y = ~cell_type,
                        type = "bar",
                        orientation = "h",
                        name = samp2,
                        marker = list(color = "#377EB8"),
                        error_x = list(
                            type = "data",
                            symmetric = FALSE,
                            array = ~(fraction_upper - fraction),
                            arrayminus = ~(fraction - fraction_lower),
                            color = "#333333",
                            thickness = 1.5
                        ),
                        hoverinfo = "text",
                        text = ~paste0(
                            cell_type, " (", samp2, ")",
                            "<br>Fraction: ", round(fraction, 3),
                            "<br>95% CI: [", round(fraction_lower, 3), ", ", round(fraction_upper, 3), "]",
                            "<br>Cells: ", n_cells
                        )
                    ) %>%
                    plotly::layout(
                        barmode = "group",
                        title = list(text = paste0("Composition: ", samp1, " vs ", samp2), font = list(size = 14)),
                        xaxis = list(title = "Fraction", range = c(0, 1)),
                        yaxis = list(title = ""),
                        legend = list(orientation = "h", yanchor = "bottom", y = 1.02, xanchor = "right", x = 1),
                        margin = list(l = 120)
                    ) %>%
                    sanitize_plotly_buttons() %>%
                    sanitize_plotly_download(globals)

                fig
            }) %>% bindCache(dataset(), group_field(), metacell_types(), cell_type_colors(), input$samp1, input$samp2, globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)

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
                req(globals$current_tab == "samples")
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
                    sanitize_plotly_download(globals)

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
            }) %>% bindCache(dataset(), input$x_axis_var, input$x_axis_type, input$y_axis_var, input$y_axis_type, input$color_by_type, input$color_by_var, metacell_types(), cell_type_colors(), input$gene_gene_point_size, input$gene_gene_stroke, input$selected_cell_types, group_field(), globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)

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
