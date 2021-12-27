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

    # cell type has only a single metacell
    if (!is.matrix(mat)) {
        mat <- as.matrix(mat)
        colnames(mat) <- metacells
    }

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

get_metacells_egc <- function(metacells, dataset, projected = FALSE) {
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat")
        mc_sum <- get_mc_data(dataset, "projected_mat_sum")
    } else {
        mc_mat <- get_mc_data(dataset, "mc_mat")
        mc_sum <- get_mc_data(dataset, "mc_sum")
    }

    mc_egc <- t(t(mc_mat[, metacells]) / mc_sum[metacells])
    return(mc_egc)
}

get_cell_types_mat <- function(cell_types, metacell_types, dataset, projected = FALSE) {
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat")
    } else {
        mc_mat <- get_mc_data(dataset, "mc_mat")
    }

    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()

    mc_mat <- mc_mat[, names(mc_types)]

    if (is.null(dim(mc_mat))) {
        ct_mat <- as.matrix(mc_mat)
        colnames(ct_mat) <- mc_types
        return(ct_mat)
    }

    ct_mat <- t(tgs_matrix_tapply(mc_mat, mc_types, sum))

    return(ct_mat)
}

get_cell_types_egc <- function(cell_types, metacell_types, dataset, mat = NULL, projected = FALSE) {
    mat <- mat %||% get_cell_types_mat(cell_types, metacell_types, dataset, projected = projected)

    ct_egc <- t(t(mat) / colSums(mat))

    return(ct_egc)
}

get_samples_mat <- function(cell_types, metacell_types, dataset) {
    mc_mat <- get_mc_data(dataset, "mc_mat")

    # samp_mc_count <- get_samp_mc_count(dataset)
    # samp_mc_count <- samp_mc_count[, colnames(samp_mc_count) != "-1"]

    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()

    mc_mat <- mc_mat[, names(mc_types)]
    samp_mc_frac <- get_samp_mc_frac(dataset)[, names(mc_types)]
    samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)

    samp_mat <- mc_mat %*% t(samp_mc_frac)

    # samp_egc <- mc_egc %*% t(samp_mc_frac)
    # samp_mc_count <- samp_mc_count[, names(mc_types)]
    # mc_samp_frac <- t(t(samp_mc_count) / colSums(samp_mc_count))

    # samp_mat <- mc_mat %*% t(mc_samp_frac)

    return(samp_mat)
}

get_samples_egc <- function(cell_types, metacell_types, dataset) {
    mc_egc <- get_mc_egc(dataset)

    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()

    mc_egc <- mc_egc[, names(mc_types)]
    samp_mc_frac <- get_samp_mc_frac(dataset)[, names(mc_types)]
    samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)

    samp_egc <- mc_egc %*% t(samp_mc_frac)

    return(samp_egc)
}

get_samples_gene_egc <- function(gene, dataset, metacells = NULL) {
    g <- get_gene_egc(gene, dataset)

    samp_mc_frac <- get_samp_mc_frac(dataset)
    if (!is.null(metacells)) {
        g <- g[metacells]
        samp_mc_frac <- samp_mc_frac[, metacells]
        samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)
        samp_mc_frac[is.na(samp_mc_frac)] <- 0
    }

    samp_egc <- samp_mc_frac %*% g
    samp_egc <- setNames(samp_egc[, 1], rownames(samp_egc))

    return(samp_egc)
}
