#' Plot mc mc scatter of gene expression
#'
#' @param df output of calc_mc_mc_gene_df
#' @param metacell1 id of the first metacell
#' @param metacell2 id of the second metacell
#'
#' @noRd
plot_mc_mc_gene <- function(df, metacell1, metacell2, highlight = NULL, label_prefix = "MC #", x_label_suffix = "") {
    xylims <- expr_breaks

    xmax <- min(c(1:length(xylims))[xylims >= max(df[, metacell1]) - 1e-10])
    xmin <- max(c(1:length(xylims))[xylims <= min(df[, metacell1]) + 1e-10])
    ymax <- min(c(1:length(xylims))[xylims >= max(df[, metacell2]) - 1e-10])
    ymin <- max(c(1:length(xylims))[xylims <= min(df[, metacell2]) + 1e-10])

    if (!is.null(highlight)) {
        prev_levels <- levels(df$col)
        df <- df %>%
            mutate(col = ifelse(gene %in% highlight, "green", as.character(col))) %>%
            mutate(col = factor(col, levels = c("green", prev_levels)))
    }

    # hover texts (see the hovertemplate in render_mc_mc_gene_plotly)
    df <- df %>%
        arrange(col) %>%
        mutate(
            gene_name_suffix = substring(gene_name, nchar(gene) + 1),
            expr_text1 = scales::scientific(!!sym(metacell1)),
            expr_text2 = scales::scientific(!!sym(metacell2)),
            diff_text = round(diff, digits = 3),
            pval_text = ifelse(is.na(pval), "Not computed", scales::scientific(pval, digits = 2))
        )

    p <- df %>%
        ggplot(aes(x = !!sym(metacell1), y = !!sym(metacell2), label = gene, customdata = gene, col = col)) +
        geom_point(size = 1, alpha = 1) +
        scale_x_continuous(limits = c(xylims[xmin], xylims[xmax]), trans = "log2", breaks = xylims[xmin:xmax], labels = scales::scientific(xylims[xmin:xmax])) +
        scale_y_continuous(limits = c(xylims[ymin], xylims[ymax]), trans = "log2", breaks = xylims[ymin:ymax], labels = scales::scientific(xylims[ymin:ymax])) +
        xlab(glue("Expression in {label_prefix}{metacell1}{x_label_suffix}")) +
        ylab(glue("Expression in {label_prefix}{metacell2}")) +
        # a manual scale (rather than scale_color_identity) makes ggplotly draw a trace per
        # color, instead of sending the color of each of the ~30K genes (the legend is hidden)
        scale_color_manual(values = rlang::set_names(levels(df$col))) +
        theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))

    return(p)
}

# The axes and the plotly source of the diff. expression scatter in each comparison mode (req()s the
# inputs that the mode needs)
diff_expr_plot_axes <- function(input, mode, metacell_names = NULL, cell_type_colors = NULL, source_suffix = "") {
    x_label_suffix <- ""
    if (is.null(mode) || mode == "MCs") {
        req(!is.null(metacell_names))
        req(input$metacell1)
        req(input$metacell2)
        req(input$metacell1 %in% metacell_names())
        req(input$metacell2 %in% metacell_names())
        xlab <- input$metacell1
        ylab <- input$metacell2
        label_prefix <- "MC #"
        source <- glue("mc_mc_plot{source_suffix}")
    } else if (mode == "Types") {
        req(!is.null(cell_type_colors) && !is.null(cell_type_colors()))
        req(input$metacell1 %in% cell_type_colors()$cell_type)
        req(input$metacell2 %in% cell_type_colors()$cell_type)
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
        x_label_suffix <- " (corrected)"
        source <- glue("projection_diff_expr_plot{source_suffix}")
    } else if (mode == "Type") {
        req(input$metacell1)
        xlab <- "Observed"
        ylab <- "Projected"
        label_prefix <- glue::glue("{input$metacell1}: ")
        x_label_suffix <- " (corrected)"
        source <- glue("projection_diff_expr_plot{source_suffix}")
    } else if (mode == "Group") {
        label_prefix <- ""
        xlab <- "Observed"
        ylab <- "Projected"
        x_label_suffix <- " (corrected)"
        source <- glue("projection_diff_expr_plot{source_suffix}")
    } else {
        req(FALSE)
    }

    list(xlab = xlab, ylab = ylab, label_prefix = label_prefix, x_label_suffix = x_label_suffix, source = source)
}

# The genes of the diff. expression scatter: the comparison without the genes the user hid
diff_expr_plot_df <- function(df, input, dataset, mode, metacell_types = NULL, group = NULL) {
    if (isTRUE(input$hide_lateral)) {
        df <- df %>%
            filter(!(gene %in% get_mc_data(dataset(), "lateral_genes")))
    }

    if (isTRUE(input$hide_noisy)) {
        df <- df %>%
            filter(!(gene %in% get_mc_data(dataset(), "noisy_genes")))
    }

    if (isTRUE(input$show_only_fitted)) {
        req(metacell_types)
        if (mode == "MC" || mode == "Group") {
            types <- metacells_to_types(input$metacell1, metacell_types())
        } else if (mode == "Type") {
            types <- input$metacell1
        } else {
            req(FALSE)
        }
        req(types)

        gene_metadata <- get_mc_data(dataset(), "gene_metadata")
        req(!is.null(gene_metadata))
        req(has_name(gene_metadata, "fitted_gene"))

        fitted_genes <- gene_metadata %>%
            filter(cell_type %in% types, fitted_gene) %>%
            pull(gene)

        df <- df %>%
            filter(gene %in% fitted_genes)
    }

    return(df)
}

# The rows of the diff. expression table (input$diff_expr_table_rows_selected are indices into them)
diff_expr_table_df <- function(df, input, dataset) {
    if (isTRUE(input$hide_lateral)) {
        df <- df %>%
            filter(!(gene %in% get_mc_data(dataset(), "lateral_genes")))
    }

    if (isTRUE(input$hide_noisy)) {
        df <- df %>%
            filter(!(gene %in% get_mc_data(dataset(), "noisy_genes")))
    }

    df %>%
        filter(col != "gray") %>%
        arrange(diff)
}

# The highlight of the gene selected in the table: an empty trace that is the first of every figure, whose
# point is set with plotly.restyle() (a figure for ~19-33K genes took seconds to build and send again, only to
# color one point). Being the first trace, the selected gene's own point is drawn on it, as a ring.
diff_expr_highlight_trace <- function() {
    list(
        type = "scattergl",
        mode = "markers",
        x = list(),
        y = list(),
        customdata = list(),
        marker = list(color = "#00cc00", size = 14, line = list(color = "black", width = 1)),
        # hover and clicks go to the gene's own point, under the cursor too
        hoverinfo = "skip",
        showlegend = FALSE,
        name = "highlight"
    )
}

# Tell the server each time the plot is drawn: when shiny sends it, and when it's drawn again after its tab
# was hidden (inst/app/www/webgl_contexts.js). The figure that shiny sends has an empty highlight.
notify_plotly_render <- function(fig) {
    fig$jsHooks$render <- c(fig$jsHooks$render, list(list(
        code = "function(el) { setTimeout(function() { if (window.Shiny && Shiny.setInputValue) Shiny.setInputValue(el.id + '_rendered', Date.now(), {priority: 'event'}); }, 0); }",
        data = NULL
    )))
    fig
}

diff_expr_highlight_observers <- function(input, session, dataset, mc_mc_gene_scatter_df, output_id, point_of_gene) {
    proxy <- plotly::plotlyProxy(output_id, session)
    rendered_input <- paste0(output_id, "_rendered")

    selected_gene <- reactiveVal(NULL)
    # whether the plot in the browser shows a highlight
    shown <- FALSE

    show_highlight <- function(gene) {
        point <- NULL
        if (!is.null(gene)) {
            # NULL when the gene isn't in the plot (e.g. hidden), or the plot can't be drawn
            point <- tryCatch(point_of_gene(gene), error = function(e) NULL)
        }
        if (!is.null(point)) {
            plotly::plotlyProxyInvoke(proxy, "restyle", list(x = list(I(point$x)), y = list(I(point$y)), customdata = list(I(gene))), 0)
            shown <<- TRUE
        } else if (shown) {
            plotly::plotlyProxyInvoke(proxy, "restyle", list(x = list(list()), y = list(list()), customdata = list(list())), 0)
            shown <<- FALSE
        }
    }

    # the table is in the same box as the plot, so the plot is shown when its rows are selected
    observeEvent(list(input$diff_expr_table_rows_selected, input$show_diff_expr_table), ignoreNULL = FALSE, ignoreInit = TRUE, {
        rows <- input$diff_expr_table_rows_selected
        gene <- NULL
        if (isTRUE(input$show_diff_expr_table) && length(rows) > 0) {
            gene <- tryCatch(
                diff_expr_table_df(mc_mc_gene_scatter_df(), input, dataset)$gene[rows[1]],
                error = function(e) NULL
            )
            if (length(gene) == 0 || is.na(gene)) {
                gene <- NULL
            }
        }
        selected_gene(gene)
        show_highlight(gene)
    })

    # a new figure, or the figure drawn again: show the selected gene (when the table is re-drawn, e.g. for
    # another comparison, its selection is cleared too)
    observeEvent(input[[rendered_input]], {
        shown <<- FALSE
        show_highlight(selected_gene())
    })
}

render_mc_mc_gene_plotly <- function(input, output, session, ns, dataset, globals, gene_modules, mc_mc_gene_scatter_df = NULL, metacell_names = NULL, cell_type_colors = NULL, mode = NULL, source_suffix = "", dragmode = NULL, plotly_buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), metacell_types = NULL, group = NULL, output_id = "plot_mc_mc_gene_scatter") {
    fixed_mode <- mode

    # the coordinates of a gene in the plot (as ggplotly writes them: log2 of the expression)
    point_of_gene <- function(gene) {
        mode <- input$mode %||% fixed_mode
        axes <- diff_expr_plot_axes(input, mode, metacell_names, cell_type_colors, source_suffix)
        df <- diff_expr_plot_df(mc_mc_gene_scatter_df(), input, dataset, mode, metacell_types, group)
        i <- match(gene, df$gene)
        if (is.na(i)) {
            return(NULL)
        }
        list(x = log2(df[[axes$xlab]][i]), y = log2(df[[axes$ylab]][i]))
    }
    if (!is.null(mc_mc_gene_scatter_df)) {
        diff_expr_highlight_observers(input, session, dataset, mc_mc_gene_scatter_df, output_id, point_of_gene)
    }

    plotly::renderPlotly({
        req(mc_mc_gene_scatter_df)

        mode <- input$mode %||% fixed_mode
        axes <- diff_expr_plot_axes(input, mode, metacell_names, cell_type_colors, source_suffix)
        xlab <- axes$xlab
        ylab <- axes$ylab
        label_prefix <- axes$label_prefix

        df <- diff_expr_plot_df(mc_mc_gene_scatter_df(), input, dataset, mode, metacell_types, group)

        df <- df %>%
            mutate(gene_name = gene_label(gene, dataset(), gene_modules()))

        p <- plot_mc_mc_gene(
            df,
            xlab,
            ylab,
            label_prefix = label_prefix,
            x_label_suffix = axes$x_label_suffix
        ) +
            theme(axis.title.y = element_text(colour = "darkblue"), axis.title.x = element_text(colour = "darkred"))

        # customdata[0] stays the gene, which the click and selection observers read
        hovertemplate <- paste0(
            "Gene: %{customdata[0]}%{customdata[1]}<br />",
            label_prefix, xlab, " expression: %{customdata[2]}<br />",
            label_prefix, ylab, " expression: %{customdata[3]}<br />",
            "Diff (log2): %{customdata[4]}<br />",
            "P-value: %{customdata[5]}<extra></extra>"
        )

        # tooltip = "none": the hover texts that ggplotly would write (~40% of its time on ~30K genes)
        # are replaced by the hovertemplate anyway
        fig <- plotly::ggplotly(p, source = axes$source, tooltip = "none") %>%
            set_hovertemplate(p$data, "gene", c("gene_name_suffix", "expr_text1", "expr_text2", "diff_text", "pval_text"), hovertemplate) %>%
            plotly::hide_legend() %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons(buttons = plotly_buttons) %>%
            sanitize_plotly_download(globals)

        if (!is.null(dragmode)) {
            fig <- fig %>%
                plotly::layout(dragmode = dragmode)
        }

        fig$x$data <- c(list(diff_expr_highlight_trace()), fig$x$data)

        return(notify_plotly_render(fig))
    }) %>% bindCache(
        # the same render serves several tabs, with their own sources and buttons (and one cache)
        output_id, source_suffix, dragmode, plotly_buttons, fixed_mode,
        dataset(),
        if (!is.null(mc_mc_gene_scatter_df)) mc_mc_gene_scatter_df(),
        input$mode,
        input$metacell1,
        input$metacell2,
        input$samp1,
        input$samp2,
        if (isTRUE((input$mode %||% fixed_mode) == "Types") && !is.null(cell_type_colors)) cell_type_colors(),
        input$hide_lateral,
        input$hide_noisy,
        if (isTRUE(input$show_only_fitted)) list(TRUE, if (!is.null(metacell_types)) metacell_types(), if (!is.null(group)) group()),
        gene_modules(), # the genes' labels
        globals$plotly_format,
        globals$plotly_width,
        globals$plotly_height,
        globals$plotly_scale
    )
}

render_mc_mc_gene_diff_table <- function(input, output, session, ns, dataset, mc_mc_gene_scatter_df) {
    DT::renderDT(
        {
            # no table: NULL also clears its selection
            if (!isTRUE(input$show_diff_expr_table)) {
                return(NULL)
            }

            if (!is.null(input$mode) && input$mode == "MCs") {
                req(input$metacell1)
                req(input$metacell2)
                req(input$metacell1 %in% colnames(get_mc_data(dataset(), "mc_mat")))
                req(input$metacell2 %in% colnames(get_mc_data(dataset(), "mc_mat")))
            }

            # the plot above shows the reason when there is no comparison (e.g. empty groups)
            df <- tryCatch(mc_mc_gene_scatter_df(), validation = function(e) NULL)
            if (is.null(df)) {
                return(NULL)
            }
            df <- diff_expr_table_df(df, input, dataset)

            DT::datatable(
                df %>%
                    mutate(type = annotate_genes(gene, dataset())) %>%
                    select(Gene = gene, `Diff (log2)` = diff, `P-value` = pval, Type = type) %>%
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
        },
        server = FALSE
    )
}
