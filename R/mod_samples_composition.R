# mod_samples_composition.R - Sample-stacked-types + composition CI helpers
#
# Split from R/mod_samples.R (2026-05-01). Owns the group_composition
# reactive and the two outputs that consume it: plot_sample_stacked_types
# (top-of-tab stacked bars) and the composition CI panel (sample-vs-sample
# Wilson-interval bar chart). Companions: mod_samples_qc.R, mod_samples_de.R.

#' Set up cell-type composition reactives + outputs for the Samples tab
#'
#' Registers `output$plot_sample_stacked_types`, `output$composition_ci_box`,
#' and `output$plot_composition_ci`. Internally creates and reuses the
#' `group_composition` cached reactive.
#'
#' @noRd
setup_samples_composition <- function(input, output, session, ns,
                                      dataset, metacell_types, cell_type_colors,
                                      group_field, state) {
    # Shared reactive for cell-type composition, used by both the stacked
    # bar plot and the composition CI panel.  Cached per (dataset, group_field,
    # metacell_types) so switching back to a previous grouping field is instant.
    group_composition <- reactive({
        req(dataset())
        req(has_cell_gene_umis(dataset()))
        gf <- group_field()
        get_group_cell_type_composition(dataset(), gf, metacell_types())
    }) %>% bindCache(dataset(), group_field(), metacell_types())

    output$plot_sample_stacked_types <- plot_sample_stacked_types(dataset, state, metacell_types, cell_type_colors, input, group_field, composition_reactive = group_composition)

    # --- Composition Confidence Intervals Panel ---
    composition_ci_data <- reactive({
        req(dataset())
        req(has_cell_gene_umis(dataset()))
        req(input$samp1)
        req(input$samp2)

        # Reuse the shared composition reactive (already cached)
        composition <- group_composition()
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
        req(state$tab_state$current_tab == "samples")
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
            sanitize_plotly_download(state)

        fig
    }) %>% bindCache(dataset(), group_field(), metacell_types(), cell_type_colors(), input$samp1, input$samp2, state$session_ui$plotly_format, state$session_ui$plotly_width, state$session_ui$plotly_height, state$session_ui$plotly_scale)
}
