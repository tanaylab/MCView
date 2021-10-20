get_mc_egc <- function(dataset, genes = NULL) {
    mc_mat <- get_mc_data(dataset, "mc_mat")
    mc_sum <- get_mc_data(dataset, "mc_sum")

    if (!is.null(genes)) {
        mc_mat <- mc_mat[genes, ]
    }

    return(t(t(mc_mat) / mc_sum))
}

get_mc_fp <- function(dataset, genes = NULL) {
    mc_egc <- get_mc_egc(dataset, genes = genes)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / apply(mc_egc_norm, 1, median, na.rm = TRUE)

    return(mc_fp)
}

filter_mat_by_cell_types <- function(mat, cell_types, metacell_types) {
    metacells <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        pull(metacell)

    mat <- mat[, metacells]
    return(mat)
}

get_gene_egc <- function(gene, dataset) {
    mc_mat <- get_mc_data(dataset, "mc_mat")
    mc_sum <- get_mc_data(dataset, "mc_sum")

    return(mc_mat[gene, ] / mc_sum)
}

get_gene_fp <- function(gene, dataset) {
    mc_egc_norm <- get_gene_egc(gene, dataset) + 1e-5

    mc_fp <- mc_egc_norm / median(mc_egc_norm, na.rm = TRUE)
    return(mc_fp)
}

get_metacells_egc <- function(metacells, dataset) {
    mc_mat <- get_mc_data(dataset, "mc_mat")
    mc_sum <- get_mc_data(dataset, "mc_sum")

    mc_egc <- t(t(mc_mat[, metacells]) / mc_sum[metacells])
    return(mc_egc)
}

get_markers <- function(dataset) {
    marker_genes <- get_mc_data(dataset, "marker_genes")
    if (!is.null(marker_genes)) {
        return(marker_genes)
    }

    marker_genes <- calc_marker_genes(get_mc_egc(dataset), 20)
    serialize_shiny_data(marker_genes, "marker_genes", dataset = dataset, cache_dir = cache_dir)

    return(marker_genes)
}
