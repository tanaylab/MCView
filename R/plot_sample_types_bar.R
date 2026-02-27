plot_sample_stacked_types <- function(dataset, globals, metacell_types, cell_type_colors, input, group_field = NULL) {
    plotly::renderPlotly({
        gf <- if (!is.null(group_field)) group_field() else "samp_id"
        use_cells <- has_cell_gene_umis(dataset()) && !is.null(group_field)

        if (use_cells) {
            # Cell-based composition using the selected grouping field
            composition <- get_group_cell_type_composition(dataset(), gf, metacell_types())

            samp_types <- composition %>%
                dplyr::rename(Sample = group_id, `Fraction` = fraction, `# of cells` = n_cells) %>%
                dplyr::mutate(
                    `Cell type` = factor(cell_type, levels = cell_type_colors()$cell_type)
                )

            if (input$sample_types_ordering == "Default") {
                # Build a fraction matrix for slanter ordering
                fracs_wide <- composition %>%
                    dplyr::select(group_id, cell_type, fraction) %>%
                    tidyr::pivot_wider(names_from = cell_type, values_from = fraction, values_fill = 0) %>%
                    tibble::column_to_rownames("group_id") %>%
                    as.matrix()
                ord <- slanter::slanted_orders(fracs_wide)
                samp_types <- samp_types %>%
                    dplyr::mutate(Sample = factor(Sample, levels = rownames(fracs_wide)[ord$rows]))
                xlab <- gf
            } else {
                # Metadata-based ordering: try samp_metadata, otherwise use alphabetical
                samp_md <- get_samp_metadata(dataset())
                if (!is.null(samp_md) && input$sample_types_ordering %in% colnames(samp_md)) {
                    samp_ord <- samp_md %>%
                        dplyr::arrange(!!sym(input$sample_types_ordering)) %>%
                        dplyr::pull(samp_id)
                    samp_types <- samp_types %>%
                        dplyr::mutate(
                            Sample = factor(Sample, levels = intersect(samp_ord, unique(as.character(Sample)))),
                            `# of cells` = paste0(`# of cells`, paste0("\n", input$sample_types_ordering, ": ", samp_md[[input$sample_types_ordering]]))
                        )
                } else {
                    samp_types <- samp_types %>%
                        dplyr::mutate(Sample = factor(Sample, levels = sort(unique(as.character(Sample)))))
                }
                xlab <- input$sample_types_ordering
            }
        } else {
            # Original metacell-weighted approach
            samp_mc_count <- get_samp_mc_count(dataset())
            mc_types <- metacell_types() %>%
                select(metacell, cell_type) %>%
                deframe()
            mc_types <- mc_types[colnames(samp_mc_count)]
            samp_mc_count <- t(tgs_matrix_tapply(samp_mc_count, mc_types, sum, na.rm = TRUE))

            samp_types <- samp_mc_count %>%
                as.data.frame() %>%
                rownames_to_column("Sample") %>%
                gather("cell_type", "n", -Sample) %>%
                group_by(Sample) %>%
                mutate(`Fraction` = n / sum(n)) %>%
                as_tibble()

            samp_types <- samp_types %>%
                mutate(
                    `Cell type` = factor(cell_type, levels = cell_type_colors()$cell_type),
                    `# of cells` = n
                )

            if (input$sample_types_ordering == "Default") {
                fracs_mat <- samp_mc_count / rowSums(samp_mc_count)
                ord <- slanter::slanted_orders(fracs_mat)
                samp_types <- samp_types %>%
                    mutate(
                        Sample = factor(Sample, levels = rownames(fracs_mat)[ord$rows])
                    )
                xlab <- "Sample"
            } else {
                samp_md <- get_samp_metadata(dataset())
                req(input$sample_types_ordering %in% colnames(samp_md))
                samp_ord <- samp_md %>%
                    arrange(!!sym(input$sample_types_ordering)) %>%
                    pull(samp_id)

                samp_types <- samp_types %>%
                    mutate(
                        Sample = factor(Sample, levels = samp_ord),
                        `# of cells` = paste0(`# of cells`, paste0("\n", input$sample_types_ordering, ": ", samp_md[[input$sample_types_ordering]]))
                    )
                xlab <- input$sample_types_ordering
            }
        }

        p <- samp_types %>%
            ggplot(aes(x = Sample, y = `Fraction`, fill = `Cell type`, customdata = Sample, tooltip_text = `# of cells`)) +
            geom_col() +
            scale_fill_manual(name = "", values = cell_type_colors() %>%
                select(cell_type, color) %>%
                deframe()) +
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            ) +
            scale_y_continuous(expand = c(0, 0)) +
            ylab("Fraction") +
            xlab(xlab) +
            guides(fill = "none")

        p <- plotly::ggplotly(p, source = "samp_types_plot") %>%
            sanitize_plotly_buttons() %>%
            sanitize_plotly_download(globals) %>%
            plotly::event_register("plotly_click")

        return(p)
    }) %>% bindCache(dataset, metacell_types, cell_type_colors, input$sample_types_ordering, group_field, globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)
}
