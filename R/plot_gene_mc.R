#' Plot gene gene scatter of metacells
#'
#' @param dataset name of metacell object
#' @param g1 name of the first gene
#' @param g2 name of the second gene
#'
#'
#' @noRd
plot_gg_over_mc <- function(dataset, g1, g2, metacell_types = get_mc_data(dataset, "metacell_types"), cell_type_colors = get_mc_data(dataset, "cell_type_colors"), plot_text = TRUE, point_size = initial_scatters_point_size(dataset), stroke = initial_scatters_stroke(dataset)) {
    egc_g1 <- get_gene_egc(g1, dataset) + egc_epsilon
    egc_g2 <- get_gene_egc(g2, dataset) + egc_epsilon

    egc_g1 <- egc_g1[metacell_types$metacell]
    egc_g2 <- egc_g2[metacell_types$metacell]

    df <- metacell_types %>%
        mutate(
            !!g1 := egc_g1,
            !!g2 := egc_g2
        ) %>%
        mutate(
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_na_value_to_level(cell_type, "(Missing)")) %>%
        rename(
            `Cell type` = cell_type
        )

    if (has_name(df, "mc_age")) {
        df <- df %>% rename(`Age` = mc_age)
    }

    xylims <- expr_breaks

    xmax <- min(c(1:length(xylims))[xylims >= max(egc_g1) - 1e-10])
    xmin <- max(c(1:length(xylims))[xylims <= min(egc_g1) + 1e-10])
    ymax <- min(c(1:length(xylims))[xylims >= max(egc_g2) - 1e-10])
    ymin <- max(c(1:length(xylims))[xylims <= min(egc_g2) + 1e-10])

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
        geom_point(size = point_size, shape = 21, stroke = stroke) +
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
plot_gene_time_over_mc <- function(dataset, gene, metacell_types = get_mc_data(dataset, "metacell_types"), cell_type_colors = get_mc_data(dataset, "cell_type_colors"), point_size = initial_scatters_point_size(dataset), stroke = initial_scatters_stroke(dataset), plot_text = TRUE) {
    egc_gene <- get_gene_egc(gene, dataset) + egc_epsilon
    egc_gene <- egc_gene[metacell_types$metacell]

    df <- metacell_types %>%
        mutate(
            !!gene := egc_gene,
            `Top genes` = glue("{top1_gene} ({round(top1_lfp, digits=2)}), {top2_gene} ({round(top2_lfp, digits=2)})")
        ) %>%
        mutate(cell_type = factor(cell_type, levels = sort(as.character(cell_type_colors$cell_type)))) %>%
        mutate(cell_type = forcats::fct_na_value_to_level(cell_type, "(Missing)")) %>%
        rename(
            `Cell type` = cell_type,
            `Age` = mc_age
        )

    ylims <- expr_breaks
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
        geom_point(size = point_size, shape = 21, stroke = stroke) +
        scale_x_continuous(limits = c(xmin, xmax)) +
        scale_y_continuous(limits = c(ylims[ymin], ylims[ymax]), trans = "log2", breaks = ylims[ymin:ymax], labels = scales::scientific(ylims[ymin:ymax])) +
        scale_fill_manual(values = col_to_ct) +
        xlab("Metacell age (E[t])") +
        ylab(glue("{gene} Expression")) +
        guides(color = "none")

    if (plot_text) {
        p <- p + geom_text(size = 1, color = "black")
    }

    return(p)
}

connect_gene_plots <- function(input, output, session, ns, source) {
    # Connect the legend of the gene/gene plot to the expression/time plots
    plot_gene_gene_mc_proxy <- plotly::plotlyProxy(ns("plot_gene_gene_mc"), session)
    plot_gene_age_mc1_proxy <- plotly::plotlyProxy(ns("plot_gene_age_mc1"), session)
    plot_gene_age_mc2_proxy <- plotly::plotlyProxy(ns("plot_gene_age_mc2"), session)

    observe({
        restyle_events <- plotly::event_data(source = source, event = "plotly_restyle")

        plotly::plotlyProxyInvoke(plot_gene_age_mc1_proxy, "restyle", restyle_events[[1]], restyle_events[[2]])
        plotly::plotlyProxyInvoke(plot_gene_age_mc2_proxy, "restyle", restyle_events[[1]], restyle_events[[2]])

        req(is.null(input$color_by_var) || input$color_by_var == "Cell type")
        plotly::plotlyProxyInvoke(plot_gene_gene_mc_proxy, "restyle", restyle_events[[1]], restyle_events[[2]])
    })
}

initial_scatters_point_size <- function(dataset, screen_width = NULL, screen_height = NULL, weight = 2, atlas = FALSE) {
    if (!is.null(config$datasets[[dataset]]$scatters_point_size)) {
        return(config$datasets[[dataset]]$scatters_point_size)
    } else if (!is.null(config$scatters_point_size)) {
        return(config$scatters_point_size)
    }
    n_metacells <- length(get_mc_data(dataset, "mc_sum", atlas = atlas))
    screen_width <- screen_width %||% 1920
    screen_height <- screen_height %||% 1080
    desired_area <- screen_width * screen_height / (n_metacells * 2000) * weight
    size_pixels <- sqrt(desired_area / pi)

    return(max(1, min(2, size_pixels * 2)))
}

initial_scatters_stroke <- function(dataset) {
    if (!is.null(config$datasets[[dataset]]$scatters_stroke)) {
        return(config$datasets[[dataset]]$scatters_stroke)
    } else if (!is.null(config$scatters_stroke)) {
        return(config$scatters_stroke)
    }
    return(0.2)
}
