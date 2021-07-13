#' calculate mc mc gene expression dataframe
#'
#' @param dataset name of metacell object
#' @param metacell1 id of the first metacell
#' @param metacell2 id of the second metacell
#'
#' @noRd
calc_mc_mc_gene_df <- function(dataset, metacell1, metacell2, diff_thresh = 1.5, pval_thresh = 0.01) {
    mc_mat <- get_mc_data(dataset, "mc_mat")

    egc <- get_metacells_egc(c(metacell1, metacell2), dataset) + egc_epsilon

    df <- egc %>%
        as.data.frame()

    df$diff <- log2(df[, 1]) - log2(df[, 2])

    f <- rownames(df)[abs(df$diff) >= diff_thresh]

    m <- mc_mat[f, c(metacell1, metacell2)]
    tots <- colSums(m)
    pvals <- apply(m, 1, function(x) suppressWarnings(chisq.test(matrix(c(x, tots), nrow = 2))$p.value))

    df$pval <- NA
    df[f, ]$pval <- pvals

    df <- df %>%
        rownames_to_column("gene") %>%
        as_tibble()

    df <- df %>%
        mutate(col = case_when(
            diff >= 1.5 & pval <= pval_thresh ~ "darkred",
            diff <= -1.5 & pval <= pval_thresh ~ "darkblue",
            TRUE ~ "gray"
        ))

    return(df)
}


#' Plot mc mc scatter of gene expression
#'
#' @param df output of calc_mc_mc_gene_df
#' @param metacell1 id of the first metacell
#' @param metacell2 id of the second metacell
#'
#' @noRd
plot_mc_mc_gene <- function(df, metacell1, metacell2, highlight = NULL) {
    xylims <- c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1)

    xmax <- min(c(1:length(xylims))[xylims >= max(df[, metacell1])])
    xmin <- max(c(1:length(xylims))[xylims <= min(df[, metacell1])])
    ymax <- min(c(1:length(xylims))[xylims >= max(df[, metacell2])])
    ymin <- max(c(1:length(xylims))[xylims <= min(df[, metacell2])])

    if (!is.null(highlight)) {
        df <- df %>% mutate(col = ifelse(gene %in% highlight, "green", col))
    }

    df <- df %>%
        mutate(
            expr_text1 = scales::scientific(!!sym(metacell1)),
            expr_text2 = scales::scientific(!!sym(metacell2)),
            pval_text = ifelse(is.na(pval), "Not computed", scales::scientific(pval, digits = 2)),
            Gene = paste(
                glue("{gene}"),
                glue("MC #{metacell1} expression: {expr_text1}"),
                glue("MC #{metacell2} expression: {expr_text2}"),
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
        xlab(glue("Expression in metacell #{metacell1}")) +
        ylab(glue("Expression in metacell #{metacell2}")) +
        scale_color_identity() +
        theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 1))

    return(p)
}

render_mc_mc_gene_plotly <- function(input, output, session, ns, dataset, mc_mc_gene_scatter_df) {
    plotly::renderPlotly({
        req(input$metacell1)
        req(input$metacell2)

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

        plotly::ggplotly(
            plot_mc_mc_gene(
                mc_mc_gene_scatter_df(),
                input$metacell1,
                input$metacell2,
                highlight = gene
            ) +
                theme(axis.title.y = element_text(colour = "darkblue"), axis.title.x = element_text(colour = "darkred")),
            tooltip = "tooltip_text",
            source = "mc_mc_plot"
        ) %>%
            plotly::hide_legend() %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons()
    })
}

render_mc_mc_gene_diff_table <- function(input, output, session, ns, dataset, mc_mc_gene_scatter_df) {
    DT::renderDT({
        req(input$metacell1)
        req(input$metacell2)
        if (input$show_diff_expr_table) {
            DT::datatable(
                mc_mc_gene_scatter_df() %>%
                    filter(col != "gray") %>%
                    arrange(diff) %>%
                    select(Gene = gene, `Diff (log2)` = diff, `P-value` = pval) %>%
                    mutate(GeneCards = glue("<a href='{link}' target='_blank'>Open</a>", link = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", Gene))),
                selection = "single",
                escape = FALSE
            ) %>%
                DT::formatSignif(columns = c("Diff (log2)"), digits = 3) %>%
                DT::formatSignif(columns = c("P-value"), digits = 2)
        }
    })
}
