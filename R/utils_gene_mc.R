get_cell_type_colors <- function(dataset, cell_type_color= NULL, na_color = "gray") {
    if (is.null(cell_type_colors)) {
        cell_type_color<- get_mc_data(dataset, "cell_type_colors")
    }

    res <- cell_type_color%>%
        distinct(cell_type, color) %>%
        bind_rows(tibble(cell_type = "(Missing)", color = na_color)) %>%
        deframe()
    return(res)
}

get_top_cor_gene <- function(dataset, gene, type = "pos") {
    gg_mc_top_cor <- get_mc_data(dataset, "gg_mc_top_cor")

    if (gene %in% gg_mc_top_cor$gene1) {
        df <- gg_mc_top_cor %>%
            filter(gene1 == gene) %>%
            filter(type == !!type) %>%
            filter(gene1 != gene2)
    } else {
        mc_egc <- get_mc_egc(dataset)
        lfp <- log2(mc_egc + egc_epsilon)
        cors <- tgs_cor(as.matrix(lfp[gene, ]), t(lfp), pairwise.complete.obs = TRUE, spearman = FALSE)[1, ]
        if (type == "pos") {
            cors <- head(sort(cors, decreasing = TRUE), n = 31)[-1]
        } else {
            cors <- head(sort(cors), n = 30)
        }
        df <- enframe(cors, "gene2", "cor")
    }

    res <- df %>%
        mutate(
            gene2 = as.character(gene2),
            label = glue("{gene2} ({cr})", cr = round(cor, digits = 2))
        ) %>%
        select(label, gene2) %>%
        tibble::deframe()

    return(res)
}
