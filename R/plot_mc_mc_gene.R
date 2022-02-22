#' Plot mc mc scatter of gene expression
#'
#' @param df output of calc_mc_mc_gene_df
#' @param metacell1 id of the first metacell
#' @param metacell2 id of the second metacell
#'
#' @noRd
plot_mc_mc_gene <- function(df, metacell1, metacell2, highlight = NULL, label_prefix = "MC #") {
    xylims <- expr_breaks

    xmax <- min(c(1:length(xylims))[xylims >= max(df[, metacell1])])
    xmin <- max(c(1:length(xylims))[xylims <= min(df[, metacell1])])
    ymax <- min(c(1:length(xylims))[xylims >= max(df[, metacell2])])
    ymin <- max(c(1:length(xylims))[xylims <= min(df[, metacell2])])

    if (!is.null(highlight)) {
        prev_levels <- levels(df$col)
        df <- df %>%
            mutate(col = ifelse(gene %in% highlight, "green", as.character(col))) %>%
            mutate(col = factor(col, levels = c("green", prev_levels)))
    }

    df <- df %>%
        arrange(col) %>%
        mutate(
            expr_text1 = scales::scientific(!!sym(metacell1)),
            expr_text2 = scales::scientific(!!sym(metacell2)),
            pval_text = ifelse(is.na(pval), "Not computed", scales::scientific(pval, digits = 2)),
            Gene = paste(
                glue("{gene}"),
                glue("{label_prefix}{metacell1} expression: {expr_text1}"),
                glue("{label_prefix}{metacell2} expression: {expr_text2}"),
                glue("Diff (log2): {round(diff, digits=3)}"),
                glue("P-value: {pval_text}"),
                sep = "\n"
            )
        )

    p <- df %>%
        ggplot(aes(x = !!sym(metacell1), y = !!sym(metacell2), label = gene, customdata = gene, col = col, tooltip_text = Gene)) +
        geom_point(size = 1, alpha = 1) +
        scale_x_continuous(limits = c(xylims[xmin], xylims[xmax]), trans = "log2", breaks = xylims[xmin:xmax], labels = scales::scientific(xylims[xmin:xmax])) +
        scale_y_continuous(limits = c(xylims[ymin], xylims[ymax]), trans = "log2", breaks = xylims[ymin:ymax], labels = scales::scientific(xylims[ymin:ymax])) +
        xlab(glue("Expression in {label_prefix}{metacell1}")) +
        ylab(glue("Expression in {label_prefix}{metacell2}")) +
        scale_color_identity() +
        theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))

    return(p)
}

render_mc_mc_gene_plotly <- function(input, output, session, ns, dataset, mc_mc_gene_scatter_df = NULL, metacell_names = NULL, cell_type_colors = NULL, mode = NULL, source_suffix = "") {
    plotly::renderPlotly({
        req(mc_mc_gene_scatter_df)

        if (!is.null(input$diff_expr_table_rows_selected)) {
            df_sig <- mc_mc_gene_scatter_df() %>%
                filter(col != "gray") %>%
                arrange(diff)
            gene <- df_sig %>%
                slice(input$diff_expr_table_rows_selected) %>%
                pull(gene)
        } else {
            gene <- NULL
        }

        mode <- input$mode %||% mode

        if (is.null(mode) || mode == "MCs") {
            req(metacell_names)
            req(input$metacell1)
            req(input$metacell2)
            req(input$metacell1 %in% metacell_names)
            req(input$metacell2 %in% metacell_names)
            xlab <- input$metacell1
            ylab <- input$metacell2
            label_prefix <- "MC #"
            source <- glue("mc_mc_plot{source_suffix}")
        } else if (mode == "Types") {
            req(cell_type_colors)
            req(input$metacell1 %in% cell_type_colors$cell_type)
            req(input$metacell2 %in% cell_type_colors$cell_type)
            xlab <- input$metacell1
            ylab <- input$metacell2
            label_prefix <- ""
            source <- glue("ct_ct_plot{source_suffix}")
        } else if (mode == "Groups") {
            label_prefix <- ""
            xlab <- "Group A"
            ylab <- "Group B"
            source <- glue("grp_grp_plot{source_suffix}")
        } else if (mode == "Samples") {
            req(input$samp1)
            req(input$samp2)
            xlab <- input$samp1
            ylab <- input$samp2
            label_prefix <- "Sample "
            source <- glue("samp_samp_diff_expr_plot{source_suffix}")
        } else if (mode == "MC") {
            req(input$metacell1)
            xlab <- "Observed"
            ylab <- "Projected"
            label_prefix <- glue::glue("MC #{input$metacell1}: ")
            source <- glue("projection_diff_expr_plot{source_suffix}")
        } else if (mode == "Type") {
            req(input$metacell1)
            xlab <- "Observed"
            ylab <- "Projected"
            label_prefix <- glue::glue("{input$metacell1}: ")
            source <- glue("projection_diff_expr_plot{source_suffix}")
        } else if (mode == "Group") {
            label_prefix <- ""
            xlab <- "Observed"
            ylab <- "Projected"
            source <- glue("projection_diff_expr_plot{source_suffix}")
        }

        fig <- plotly::ggplotly(
            plot_mc_mc_gene(
                mc_mc_gene_scatter_df(),
                xlab,
                ylab,
                highlight = gene,
                label_prefix = label_prefix
            ) +
                theme(axis.title.y = element_text(colour = "darkblue"), axis.title.x = element_text(colour = "darkred")),
            tooltip = "tooltip_text",
            source = source
        ) %>%
            plotly::hide_legend() %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons()

        return(fig)
    })
}

render_mc_mc_gene_diff_table <- function(input, output, session, ns, dataset, mc_mc_gene_scatter_df) {
    DT::renderDT(
        {
            if (!is.null(input$mode) && input$mode == "MCs") {
                req(input$metacell1)
                req(input$metacell2)
            }

            if (input$show_diff_expr_table) {
                DT::datatable(
                    mc_mc_gene_scatter_df() %>%
                        filter(col != "gray") %>%
                        arrange(diff) %>%
                        select(Gene = gene, `Diff (log2)` = diff, `P-value` = pval, any_of(c("Systematic", "Disjoined"))) %>%
                        mutate(GeneCards = glue("<a href='{link}' target='_blank'>Open</a>", link = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", Gene))),
                    selection = "single",
                    escape = FALSE,
                    rownames = FALSE,
                    extensions = c("Buttons", "Responsive"),
                    options = list(
                        dom = "Bfrtip",
                        buttons = list(
                            list(
                                extend = "copy",
                                exportOptions = list(columns = 0:2)
                            ),
                            list(
                                extend = c("csv"),
                                exportOptions = list(columns = 0:2)
                            ),
                            list(
                                extend = "excel",
                                exportOptions = list(columns = 0:2)
                            )
                        )
                    )
                ) %>%
                    DT::formatSignif(columns = c("Diff (log2)"), digits = 3) %>%
                    DT::formatSignif(columns = c("P-value"), digits = 2)
            }
        },
        server = FALSE
    )
}
