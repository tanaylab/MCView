#' Plot gene gene scatter of metacells
#'
#' @param dataset name of metacell object
#' @param g1 name of the first gene
#' @param g2 name of the second gene
#'
#'
#' @noRd
plot_gg_over_mc <- function(dataset, g1, g2, metacell_type= get_mc_data(dataset, "metacell_types"), cell_type_color= get_mc_data(dataset, "cell_type_colors"), plot_text = TRUE) {
    egc_g1 <- get_gene_egc(g1, dataset) + egc_epsilon
    egc_g2 <- get_gene_egc(g2, dataset) + egc_epsilon

    df <- metacell_type%>%
        mutate(
            !!g1 := egc_g1,
            !!g2 := egc_g2
        ) %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        rename(
            `Cell type` = cell_type
        )

    if (has_name(df, "mc_age")) {
        df <- df %>% rename(`Age` = mc_age)
    }

    xylims <- c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1)

    xmax <- min(c(1:length(xylims))[xylims >= max(egc_g1)])
    xmin <- max(c(1:length(xylims))[xylims <= min(egc_g1)])
    ymax <- min(c(1:length(xylims))[xylims >= max(egc_g2)])
    ymin <- max(c(1:length(xylims))[xylims <= min(egc_g2)])

    col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)

    df <- df %>%
        mutate(
            expr_text1 = scales::scientific(!!sym(g1)),
            expr_text2 = scales::scientific(!!sym(g2)),
            Metacell = paste(
                glue("{metacell}"),
                glue("{g1} expression: {expr_text1}"),
                glue("{g2} expression: {expr_text2}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                ifelse(has_name(df, "Age"), glue("Metacell age (E[t]): {round(Age, digits=2)}"), ""),
                sep = "\n"
            )
        )

    p <- ggplot(data = df, aes(x = !!sym(g1), y = !!sym(g2), fill = `Cell type`, label = metacell, customdata = metacell, tooltip_text = Metacell)) +
        geom_point(size = 2.5, shape = 21) +
        scale_x_continuous(limits = c(xylims[xmin], xylims[xmax]), trans = "log2", breaks = xylims[xmin:xmax], labels = scales::scientific(xylims[xmin:xmax])) +
        scale_y_continuous(limits = c(xylims[ymin], xylims[ymax]), trans = "log2", breaks = xylims[ymin:ymax], labels = scales::scientific(xylims[ymin:ymax])) +
        scale_fill_manual(values = col_to_ct) +
        xlab(glue("{g1} Expression")) +
        ylab(glue("{g2} Expression")) +
        theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))

    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    return(p)
}

#' Plot gene expression vs time scatter of metacells
#'
#' @param dataset name of metacell object
#' @param gene name of the gene
#'
#' @noRd
plot_gene_time_over_mc <- function(dataset, gene, metacell_type= get_mc_data(dataset, "metacell_types"), cell_type_color= get_mc_data(dataset, "cell_type_colors")) {
    egc_gene <- get_gene_egc(gene, dataset) + egc_epsilon

    df <- metacell_type%>%
        mutate(
            !!gene := egc_gene,
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_explicit_na(cell_type)) %>%
        rename(
            `Cell type` = cell_type,
            `Age` = mc_age
        )

    ylims <- c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1)
    ymax <- min(c(1:length(ylims))[ylims >= max(egc_gene)])
    ymin <- max(c(1:length(ylims))[ylims <= min(egc_gene)])

    xmin <- min(metacell_types$mc_age)
    xmax <- max(metacell_types$mc_age)

    col_to_ct <- get_cell_type_colors(dataset, cell_type_colors)


    # We call the text field "Metacell" in order for plotly to show "Metacell:" in the tooltip
    df <- df %>%
        mutate(
            expr_text = scales::scientific(!!sym(gene)),
            Metacell = paste(
                glue("{metacell}"),
                glue("{gene} expression: {expr_text}"),
                glue("Metacell age (E[t]): {round(Age, digits=2)}"),
                glue("Cell type: {`Cell type`}"),
                glue("Top genes: {`Top genes`}"),
                sep = "\n"
            )
        )

    p <- ggplot(data = df, aes(x = `Age`, y = !!sym(gene), fill = `Cell type`, label = metacell, customdata = metacell, tooltip_text = Metacell)) +
        geom_point(size = 2.5, shape = 21, stroke = 0.1) +
        scale_x_continuous(limits = c(xmin, xmax)) +
        scale_y_continuous(limits = c(ylims[ymin], ylims[ymax]), trans = "log2", breaks = ylims[ymin:ymax], labels = scales::scientific(ylims[ymin:ymax])) +
        geom_text(size = 1, color = "black") +
        scale_fill_manual(values = col_to_ct) +
        xlab("Metacell age (E[t])") +
        ylab(glue("{gene} Expression")) +
        guides(color = "none")

    return(p)
}

connect_gene_plots <- function(input, output, session, ns, source) {
    # Connect the legend of the gene/gene plot to the expression/time plots
    plot_gene_gene_mc_proxy <- plotly::plotlyProxy(ns("plot_gene_gene_mc"), session)
    plot_gene_age_mc1_proxy <- plotly::plotlyProxy(ns("plot_gene_age_mc1"), session)
    plot_gene_age_mc2_proxy <- plotly::plotlyProxy(ns("plot_gene_age_mc2"), session)

    observe({
        restyle_events <- plotly::event_data(source = source, event = "plotly_restyle")

        plotly::plotlyProxyInvoke(plot_gene_gene_mc_proxy, "restyle", restyle_events[[1]], restyle_events[[2]])
        plotly::plotlyProxyInvoke(plot_gene_age_mc1_proxy, "restyle", restyle_events[[1]], restyle_events[[2]])
        plotly::plotlyProxyInvoke(plot_gene_age_mc2_proxy, "restyle", restyle_events[[1]], restyle_events[[2]])
    })
}
