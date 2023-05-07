annotate_genes <- function(genes, dataset) {
    lateral_genes <- get_mc_data(dataset, "lateral_genes")
    noisy_genes <- get_mc_data(dataset, "noisy_genes")

    case_when(
        (genes %in% lateral_genes) & (genes %in% noisy_genes) ~ "lateral, noisy",
        genes %in% lateral_genes ~ "lateral",
        genes %in% noisy_genes ~ "noisy",
        TRUE ~ "other"
    )
}

gene_label <- function(genes, dataset) {
    type <- annotate_genes(genes, dataset)
    ifelse(type == "other", genes, glue("{genes} ({type})"))
}
