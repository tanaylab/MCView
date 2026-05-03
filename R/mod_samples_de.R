# mod_samples_de.R - Sample/sample + group-comparison DE helpers
#
# Split from R/mod_samples.R (2026-05-01). Hosts the two DE flows on the
# Samples tab: pairwise sample-vs-sample (samp_samp_*) and the multi-select
# group A vs B (group_comparison_*). Companions: mod_samples_composition.R,
# mod_samples_qc.R.

#' Set up DE reactives + outputs for the Samples tab
#'
#' Registers diff_expr_box / plot_samp_samp_gene_scatter / diff_expr_table
#' and group_comparison_select_ui / group_comparison_box /
#' plot_group_comparison_scatter / group_comparison_table, plus the
#' plotly-click observer that surfaces the selected gene from the group
#' comparison plot.
#'
#' @noRd
setup_samples_de <- function(input, output, session, ns,
                             dataset, metacell_types, cell_type_colors,
                             gene_modules, group_field, globals, state) {
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

    output$plot_samp_samp_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, globals, state, gene_modules, samp_samp_scatter_df, metacell_names, cell_type_colors, mode = "Samples", tab_guard = "samples")

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
        input, output, session, ns, dataset, globals, state, gene_modules,
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
}
