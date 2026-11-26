annotate_genes <- function(genes, dataset) {
    lateral <- genes %in% get_mc_data(dataset, "lateral_genes")
    noisy <- genes %in% get_mc_data(dataset, "noisy_genes")

    type <- rep("other", length(genes))
    type[lateral] <- "lateral"
    type[noisy] <- "noisy"
    type[lateral & noisy] <- "lateral, noisy"
    return(type)
}

# Called on all the genes (e.g. on every render of the Diff. Expression scatter), so it builds the labels only of the
# genes that need one, instead of building every label twice with ifelse() and glue().
gene_label <- function(genes, dataset, gene_modules = NULL) {
    genes <- as.character(genes)
    type <- annotate_genes(genes, dataset)
    new_names <- genes
    annotated <- type != "other"
    new_names[annotated] <- paste0(genes[annotated], " (", type[annotated], ")")

    if (!is.null(gene_modules)) {
        modules <- as.character(gene_modules$module)[match(genes, gene_modules$gene)]
        in_module <- !is.na(modules)
        new_names[in_module] <- paste0(new_names[in_module], "<br />Gene module: ", modules[in_module])
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
