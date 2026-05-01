# mod_samples_qc.R - Sample QC reactives + plots
#
# Split from R/mod_samples.R (2026-05-01). Owns the qc_stats reactive and
# the four outputs that consume it: sample_qc_box (UI shell),
# plot_qc_cells_per_group, plot_qc_umis_per_group, qc_summary_table.
# Companions: mod_samples_composition.R, mod_samples_de.R.

#' Set up Sample QC reactives + outputs for the Samples tab
#'
#' Registers `qc_stats` (per-group n_cells / total_umis / median_umis_per_cell)
#' as a cached reactive, plus the four outputs the QC panel renders. Highlights
#' the currently selected samp1 / samp2 in the bar plots.
#'
#' @noRd
setup_samples_qc <- function(input, output, session, ns,
                             dataset, group_field, globals) {
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
}
