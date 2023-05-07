get_cell_type_colors <- function(dataset, cell_type_colors = NULL, na_color = "gray", atlas = FALSE) {
    if (is.null(cell_type_colors)) {
        cell_type_colors <- get_mc_data(dataset, "cell_type_colors", atlas = atlas)
    }

    res <- bind_rows(
        tibble(cell_type = "(Missing)", color = na_color),
        cell_type_colors %>%
            distinct(cell_type, color) %>%
            select(cell_type, color)
    ) %>%
        deframe()

    return(res)
}

get_top_cor_gene <- function(dataset, gene, type = "pos", atlas = FALSE) {
    gg_mc_top_cor <- get_mc_data(dataset, "gg_mc_top_cor", atlas = atlas)

    if (gene %in% gg_mc_top_cor$gene1) {
        df <- gg_mc_top_cor %>%
            filter(gene1 == gene) %>%
            filter(type == !!type) %>%
            filter(gene1 != gene2)
    } else {
        mc_egc <- get_mc_egc(dataset, atlas = atlas)
        req(gene %in% rownames(mc_egc))
        lfp <- log2(mc_egc + egc_epsilon)
        cors <- tgs_cor(t(as.matrix(lfp[gene, , drop = FALSE])), t(lfp), pairwise.complete.obs = TRUE, spearman = FALSE)[1, ]
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
            gene_name = gene_label(gene2, dataset),
            label = glue("{gene_name} ({cr})", cr = round(cor, digits = 2))
        ) %>%
        select(label, gene2) %>%
        tibble::deframe()

    return(res)
}
