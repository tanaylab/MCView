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
            column(
                width = 5,
                shinydashboardPlus::box(
                    id = ns("sample_sample_box"),
                    title = "Sample/Sample",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("gene_gene_sidebar"),
                        uiOutput(ns("gene_gene_point_size_ui")),
                        uiOutput(ns("gene_gene_stroke_ui"))
                    ),
                    axis_selector("x_axis", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type")),
                    axis_selector("y_axis", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type")),
                    axis_selector("color_by", "Metadata", ns, choices = c("Metadata", "Gene", "Cell type")),
                    textOutput(ns("please_select_cell_types")),
                    textOutput(ns("no_samples1")),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_gene_mc"))
                    )
                ),
                uiOutput(ns("diff_expr_box"))
            ),
            column(
                width = 7,
                shinydashboardPlus::box(
                    id = ns("sample_projection"),
                    title = "Sample projection",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("gene_projection_sidebar"),
                        uiOutput(ns("proj_stat_ui")),
                        uiOutput(ns("set_range_ui")),
                        uiOutput(ns("expr_range_ui")),
                        uiOutput(ns("enrich_range_ui")),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("stroke_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    textOutput(ns("no_samples2")),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_proj_2d"))
                    ),
                    shinyWidgets::prettyRadioButtons(
                        ns("color_proj"),
                        label = "Color by:",
                        choices = c("Sample", "Cell type"),
                        selected = "Sample",
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    )
                ),
                uiOutput(ns("sample_info_box"))
            )
        )
    )
}


#' gene_mc sidebar UI Function
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

#' gene_mc Server Function
#'
#' @noRd
mod_samples_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns
    top_correlated_selectors(input, output, session, dataset, ns, button_labels = c("X", "Y", "Color"))

    output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", cell_type_colors = cell_type_colors())

    scatter_selectors(ns, dataset, output)
    projection_selectors(ns, dataset, output, input)

    output$sample_select_ui <- renderUI({
        req(dataset())
        req(input$color_proj)
        samp_list <- get_samples_list(dataset())
        if (length(samp_list) > 1) {
            selected2 <- samp_list[2]
        } else {
            selected2 <- samp_list[1]
        }

        picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
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

    # Projection plots
    output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, values, metacell_types, cell_type_colors, source = "proj_mc_plot_gene_tab") %>%
        bindCache(dataset(), input$color_proj, metacell_types(), cell_type_colors(), input$point_size, input$stroke, input$min_edge_size, input$set_range, input$show_selected_metacells, input$metacell1, input$metacell2, input$proj_stat, input$expr_range, input$lfp, input$samp1)

    # Info box
    output$sample_info_box <- renderUI({
        req(input$samp1)
        shinydashboardPlus::box(
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

        shinydashboardPlus::box(
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
    output$x_axis_select <- render_axis_select_ui("x_axis", "X axis", md_choices = dataset_cell_metadata_fields_numeric(dataset()), md_selected = dataset_cell_metadata_fields_numeric(dataset())[1], selected_gene = default_gene1, input = input, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE)) %>% bindCache(dataset(), ns, ns("x_axis"), input$x_axis_type)

    output$y_axis_select <- render_axis_select_ui("y_axis", "Y axis", md_choices = dataset_cell_metadata_fields_numeric(dataset()), md_selected = dataset_cell_metadata_fields_numeric(dataset())[2], selected_gene = default_gene2, input = input, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE)) %>% bindCache(dataset(), ns, ns("y_axis"), input$y_axis_type)

    output$color_by_select <- render_axis_select_ui("color_by", "Color", md_choices = c("None", dataset_cell_metadata_fields(dataset())), md_selected = "None", selected_gene = default_gene1, input = input, ns = ns, dataset = dataset, cell_types = sort(names(get_cell_type_colors(dataset())), decreasing = TRUE)) %>% bindCache(dataset(), ns, ns("color_by"), input$color_by_type)

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

        req(axis_vars_ok(dataset(), input, "cell_metadata"))

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