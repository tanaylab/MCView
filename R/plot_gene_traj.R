calc_mc_egc_t <- function(dataset, metacell, genes = NULL) {
    mct_probs_trans <- get_mc_data(dataset, "mct_probs_trans")
    mc_egc <- get_mc_egc(dataset)

    flow_prob <- mct_probs_trans[[metacell]]$probs

    # find maximal time with non-zero colSums
    t_max <- which.max(c(1:ncol(flow_prob))[colSums(flow_prob) > 0])

    # find minimal time with non-zero colSums
    t_min <- which.min(c(1:ncol(flow_prob))[colSums(flow_prob) > 0])

    flow_prob_f <- flow_prob[, t_min:t_max]
    if (is.null(ncol(flow_prob_f))) {
        flow_prob_f <- flow_prob_f %>% as.matrix(ncol = 1)
    }
    flow_prob_n <- t(t(flow_prob_f) / colSums(flow_prob_f))

    if (is.null(genes)) {
        egc_t <- mc_egc %*% flow_prob_n
    } else {
        egc_t <- mc_egc[genes, ] %*% flow_prob_n
    }


    return(egc_t)
}


plot_gene_trajectory <- function(dataset, genes, metacell, anchor_genes = NULL) {
    metacell <- as.character(metacell)
    time_annot <- get_mc_data(dataset, "time_annot")

    mc_egc_t <- calc_mc_egc_t(dataset, metacell, genes)
    lfp <- t(mc_egc_t + egc_epsilon)

    colnames(lfp) <- genes

    df <- as.data.frame(lfp) %>%
        mutate(time_bin = 1:n()) %>%
        gather("gene", "expr", -time_bin) %>%
        left_join(time_annot, by = "time_bin") %>%
        rename(time = time_desc) %>%
        as_tibble() %>%
        mutate(gene = factor(gene))

    if (!is.null(anchor_genes) && any(anchor_genes %in% genes)) {
        anchor_genes <- anchor_genes[anchor_genes %in% genes]
        df <- df %>% mutate(gene = forcats::fct_relevel(gene, anchor_genes))
    }

    ylims <- expr_breaks
    ymax <- min(c(1:length(ylims))[ylims >= max(df$expr, na.rm = TRUE)])
    ymin <- max(c(1:length(ylims))[ylims <= min(df$expr, na.rm = TRUE)])

    df <- df %>%
        mutate(
            expr_text = scales::scientific(expr),
            Gene = paste(
                glue("{gene}"),
                glue("MC #{metacell} expression: {expr_text}"),
                glue("Time: {round(time, digits=2)}"),
                sep = "\n"
            )
        )

    p <- df %>%
        ggplot(aes(x = time, y = expr, linetype = gene, group = gene, tooltip_text = Gene, customdata = gene)) +
        geom_line() +
        scale_y_continuous(
            limits = c(ylims[ymin], ylims[ymax]),
            trans = "log2",
            breaks = ylims[ymin:ymax],
            labels = scales::scientific(ylims[ymin:ymax])
        ) +
        xlab("Time (E[t])") +
        ylab(glue("Gene expression levels")) +
        scale_linetype_discrete(name = "")

    return(p)
}
