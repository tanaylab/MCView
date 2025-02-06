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

gene_label <- function(genes, dataset, gene_modules = NULL) {
    type <- annotate_genes(genes, dataset)
    if (!is.null(gene_modules)) {
        gm_names <- gene_modules %>%
            group_by(gene) %>%
            summarize(module = paste(module, collapse = ", ")) %>%
            select(gene, module) %>%
            deframe()
        modules <- gm_names[genes]
    }
    new_names <- ifelse(type == "other", genes, glue("{genes} ({type})"))

    if (!is.null(gene_modules)) {
        new_names <- ifelse(!is.na(modules), glue("{new_names}<br />Gene module: {modules}"), new_names)
    }

    return(new_names)
}

add_gene_modules <- function(genes, dataset, gene_modules) {
    modules <- tibble(gene = genes) %>%
        left_join(gene_modules, by = join_by(gene)) %>%
        pull(module) %>%
        as.character()

    return(ifelse(!is.na(modules), glue("{genes} ({modules})"), genes))
}
