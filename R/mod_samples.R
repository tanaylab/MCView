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
                width = 5,
                movable_box(
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
                        axis_selector("x_axis", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type")),
                        axis_selector("y_axis", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type")),
                        axis_selector("color_by", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type")),
                        uiOutput(ns("gene_gene_point_size_ui")),
                        uiOutput(ns("gene_gene_stroke_ui"))
                    ),
                    textOutput(ns("please_select_cell_types")),
                    textOutput(ns("no_samples1")),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_gene_mc"))
                    )
                ),
                uiOutput(ns("diff_expr_box"))
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
            uiOutput(ns("cell_type_list")),
            uiOutput(ns("sample_select_ui")),
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
            top_correlated_selectors(input, output, session, dataset, ns, button_labels = c("X", "Y", "Color"))

            output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", cell_type_colors = cell_type_colors)

            scatter_selectors(ns, dataset, output, globals)
            projection_selectors(ns, dataset, output, input, gene_modules, globals, weight = 0.6)

            output$sample_select_ui <- renderUI({
                req(dataset())
                req(input$color_proj)
                samp_list <- get_samples_list(dataset())
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
            output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, source = "proj_mc_plot_gene_tab") %>%
                bindCache(dataset(), input$color_proj, metacell_types(), cell_type_colors(), input$point_size, input$stroke, input$min_edge_size, input$set_range, input$metacell1, input$metacell2, input$proj_stat, input$expr_range, input$lfp, input$samp1, input$color_proj_gene_module, clipboard_changed(), input$graph_name)

            # Info box
            output$sample_info_box <- renderUI({
                req(input$samp1)
                movable_box(
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

                movable_box(
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
                calc_samp_samp_gene_df(dataset(), input$samp1, input$samp2, metacell_types(), cell_types = input$selected_cell_types)
            }) %>% bindCache(dataset(), input$selected_cell_type, input$samp1, input$samp2, metacell_types())

            output$plot_samp_samp_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, samp_samp_scatter_df, metacell_names(), cell_type_colors(), mode = "Samples")

            output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, samp_samp_scatter_df)

            # Metadata/Metadata plots
            render_axis_select_ui("x_axis", "X axis", "x_axis_select", md_choices = dataset_cell_metadata_fields_numeric(dataset()), md_selected = dataset_cell_metadata_fields_numeric(dataset())[1], selected_gene = default_gene1, input = input, output = output, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE), gene_modules = gene_modules, session = session)
            render_axis_select_ui("y_axis", "Y axis", "y_axis_select", md_choices = dataset_cell_metadata_fields_numeric(dataset()), md_selected = dataset_cell_metadata_fields_numeric(dataset())[2], selected_gene = default_gene2, input = input, output = output, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE), gene_modules = gene_modules, session = session)
            render_axis_select_ui("color_by", "Color", "color_by_select", md_choices = c("None", dataset_cell_metadata_fields(dataset())), md_selected = "None", selected_gene = default_gene1, input = input, output = output, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE), gene_modules = gene_modules, session = session)

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
                    if (!has_samples(dataset())) {
                        glue("No samples were loaded to MCView.\nPlease make sure your cell metadata has a field called 'samp_id' and run 'import_cell_metadata' again.")
                    } else {
                        req(FALSE)
                    }
                })
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
            }) %>% bindCache(dataset(), input$x_axis_var, input$x_axis_type, input$y_axis_var, input$y_axis_type, input$color_by_type, input$color_by_var, metacell_types(), cell_type_colors(), input$gene_gene_point_size, input$gene_gene_stroke, input$selected_cell_types)

            sample_click_observer("samp_samp_plot", session, "samp1")
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
