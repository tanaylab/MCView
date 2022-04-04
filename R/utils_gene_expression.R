get_mc_egc <- function(dataset, genes = NULL, atlas = FALSE) {
    mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
    mc_sum <- get_mc_data(dataset, "mc_sum", atlas = atlas)

    if (!is.null(genes)) {
        mc_mat <- mc_mat[genes, , drop = FALSE]
    }

    return(t(t(mc_mat) / mc_sum))
}

get_mc_fp <- function(dataset, genes = NULL, atlas = FALSE) {
    mc_egc <- get_mc_egc(dataset, genes = genes, atlas = atlas)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / apply(mc_egc_norm, 1, median, na.rm = TRUE)

    return(mc_fp)
}

get_mc_gene_modules_egc <- function(dataset, modules = NULL, gene_modules = NULL, atlas = FALSE) {
    mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
    mc_sum <- get_mc_data(dataset, "mc_sum", atlas = atlas)

    gene_modules <- gene_modules %||% get_mc_data(dataset, "gene_modules", atlas = atlas)
    if (!is.null(modules)) {
        gene_modules <- gene_modules %>%
            filter(module %in% modules) %>%
            mutate(module = forcats::fct_drop(module))
    }

    mc_mat <- mc_mat[gene_modules$gene, ]

    mc_mat <- tgs_matrix_tapply(t(mc_mat), gene_modules$module, sum)

    return(t(t(mc_mat) / mc_sum))
}

get_gene_module_egc <- function(module, dataset, gene_modules = NULL, atlas = FALSE) {
    res <- get_mc_gene_modules_egc(dataset = dataset, modules = module, gene_modules = gene_modules, atlas = atlas)
    return(res[module, ])
}

get_mc_gene_modules_fp <- function(dataset, modules = NULL, gene_modules = NULL, atlas = FALSE) {
    mc_egc <- get_mc_gene_modules_egc(dataset, modules = modules, gene_modules = gene_modules, atlas = atlas)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / apply(mc_egc_norm, 1, median, na.rm = TRUE)

    return(mc_fp)
}

filter_mat_by_cell_types <- function(mat, cell_types, metacell_types) {
    metacells <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        pull(metacell)

    mat <- mat[, metacells, drop = FALSE]

    return(mat)
}

get_gene_egc <- function(gene, dataset, projected = FALSE, atlas = FALSE) {
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat")
        mc_sum <- get_mc_data(dataset, "projected_mat_sum")
    } else {
        mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
        mc_sum <- get_mc_data(dataset, "mc_sum", atlas = atlas)
    }


    return(mc_mat[gene, ] / mc_sum)
}

get_gene_fp <- function(gene, dataset, atlas = FALSE) {
    mc_egc_norm <- get_gene_egc(gene, dataset, atlas = atlas) + 1e-5

    mc_fp <- mc_egc_norm / median(mc_egc_norm, na.rm = TRUE)
    return(mc_fp)
}

get_metacells_egc <- function(metacells, dataset, projected = FALSE, atlas = FALSE) {
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat")
        mc_sum <- get_mc_data(dataset, "projected_mat_sum")
    } else {
        mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
        mc_sum <- get_mc_data(dataset, "mc_sum", atlas = atlas)
    }

    mc_egc <- t(t(mc_mat[, metacells]) / mc_sum[metacells])
    return(mc_egc)
}

get_cell_types_mat <- function(cell_types, metacell_types, dataset, projected = FALSE, atlas = FALSE) {
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat")
    } else {
        mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
    }

    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()

    mc_mat <- mc_mat[, names(mc_types), drop = FALSE]

    ct_mat <- t(tgs_matrix_tapply(mc_mat, mc_types, sum))

    return(ct_mat)
}

get_cell_types_egc <- function(cell_types, metacell_types, dataset, mat = NULL, projected = FALSE, atlas = FALSE) {
    mat <- mat %||% get_cell_types_mat(cell_types, metacell_types, dataset, projected = projected, atlas = atlas)

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

    mc_mat <- mc_mat[, names(mc_types), drop = FALSE]
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
