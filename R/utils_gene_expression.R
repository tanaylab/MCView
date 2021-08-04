get_mc_egc <- function(dataset) {
    mc_mat <- get_mc_data(dataset, "mc_mat")
    mc_sum <- get_mc_data(dataset, "mc_sum")

    return(t(t(mc_mat) / mc_sum))
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