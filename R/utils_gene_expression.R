get_mc_egc <- function(dataset, genes = NULL, atlas = FALSE, metacells = NULL) {
    mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
    mc_sum <- get_mc_sum(dataset, atlas = atlas)

    if (!is.null(genes)) {
        mc_mat <- mc_mat[intersect(genes, rownames(mc_mat)), , drop = FALSE]
    }

    if (!is.null(metacells)) {
        mc_mat <- mc_mat[, intersect(metacells, colnames(mc_mat)), drop = FALSE]
        mc_sum <- mc_sum[intersect(metacells, names(mc_sum))]
    }

    # same values as t(t(mc_mat) / mc_sum) without the two transposed copies
    return(mc_mat / rep(unname(mc_sum), each = nrow(mc_mat)))
}

get_mc_fp <- function(dataset, genes = NULL, atlas = FALSE, metacells = NULL) {
    mc_egc <- get_mc_egc(dataset, genes = genes, atlas = atlas, metacells = metacells)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / matrixStats::rowMedians(mc_egc_norm, na.rm = TRUE)

    return(mc_fp)
}

get_mc_gene_modules_egc <- function(dataset, modules = NULL, gene_modules = NULL, atlas = FALSE) {
    mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
    mc_sum <- get_mc_sum(dataset, atlas = atlas)

    gene_modules <- gene_modules %||% get_mc_data(dataset, "gene_modules", atlas = atlas)
    if (!is.null(modules)) {
        gene_modules <- gene_modules %>%
            filter(module %in% modules) %>%
            mutate(module = forcats::fct_drop(module))
    }

    # rowsum rather than tgs_matrix_tapply: tgstat forks worker processes, which costs seconds
    # when the app process holds a large dataset
    mc_mat <- rowsum(mc_mat[gene_modules$gene, , drop = FALSE], gene_modules$module)

    return(t(t(mc_mat) / mc_sum))
}

get_gene_module_egc <- function(module, dataset, gene_modules = NULL, atlas = FALSE) {
    res <- get_mc_gene_modules_egc(dataset = dataset, modules = module, gene_modules = gene_modules, atlas = atlas)
    return(res[module, ])
}

get_mc_gene_modules_fp <- function(dataset, modules = NULL, gene_modules = NULL, atlas = FALSE) {
    mc_egc <- get_mc_gene_modules_egc(dataset, modules = modules, gene_modules = gene_modules, atlas = atlas)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / matrixStats::rowMedians(mc_egc_norm, na.rm = TRUE)

    return(mc_fp)
}

filter_mat_by_cell_types <- function(mat, cell_types, metacell_types) {
    metacells <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        pull(metacell)

    mat <- mat[, intersect(colnames(mat), metacells), drop = FALSE]

    return(mat)
}

get_gene_egc <- function(gene, dataset, projected = FALSE, atlas = FALSE, corrected = FALSE) {
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat")
        mc_sum <- get_mc_data(dataset, "projected_mat_sum")
        names(mc_sum) <- colnames(mc_mat)
    } else if (corrected) {
        mc_mat <- get_mc_data(dataset, "mc_mat_corrected", atlas = atlas)
        if (is.null(mc_mat)) { # for older versions of MCView
            mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
        }
        mc_sum <- get_mc_sum(dataset, atlas = atlas)
        names(mc_sum) <- colnames(mc_mat)
    } else {
        mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
        mc_sum <- get_mc_sum(dataset, atlas = atlas)
    }



    return(mc_mat[gene, ] / mc_sum)
}

get_gene_fp <- function(gene, dataset, atlas = FALSE) {
    mc_egc_norm <- get_gene_egc(gene, dataset, atlas = atlas) + 1e-5

    mc_fp <- mc_egc_norm / median(mc_egc_norm, na.rm = TRUE)
    return(mc_fp)
}

get_metacells_egc <- function(metacells, dataset, projected = FALSE, atlas = FALSE, corrected = FALSE) {
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat")
        mc_sum <- get_mc_data(dataset, "projected_mat_sum")
        names(mc_sum) <- colnames(mc_mat)
    } else if (corrected) {
        mc_mat <- get_mc_data(dataset, "mc_mat_corrected", atlas = atlas)
        if (is.null(mc_mat)) { # for older versions of MCView
            mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
        }
        mc_sum <- get_mc_sum(dataset, atlas = atlas)
        names(mc_sum) <- colnames(mc_mat)
    } else {
        mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
        mc_sum <- get_mc_sum(dataset, atlas = atlas)
    }

    mc_egc <- t(t(mc_mat[, metacells]) / mc_sum[metacells])
    return(mc_egc)
}

get_cell_types_mat <- function(cell_types, metacell_types, dataset, projected = FALSE, atlas = FALSE, corrected = FALSE) {
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat")
    } else if (corrected) {
        mc_mat <- get_mc_data(dataset, "mc_mat_corrected", atlas = atlas)
        if (is.null(mc_mat)) { # for older versions of MCView
            mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
        }
    } else {
        mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
    }

    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()

    ct_mat <- sum_columns_by_group(mc_mat, names(mc_types), mc_types)

    return(ct_mat)
}

# The sums of the columns of each group: a matrix with a column per group (in the order of split()).
# rowSums per group rather than tgs_matrix_tapply (which forks worker processes, slow when the app holds a
# large dataset). A dense matrix is summed in place: subsetting its columns first copied them twice (all
# the chosen columns, then those of each group), ~3 GB and 15-20 s for two groups of half of the
# metacells of a large dataset. Sparse (projected) matrices are subset as before.
sum_columns_by_group <- function(mat, columns, groups) {
    if (is.matrix(mat)) {
        col_groups <- split(match(columns, colnames(mat)), groups)
        sums <- vapply(col_groups, function(cols) matrixStats::rowSums2(mat, cols = cols, useNames = FALSE), numeric(nrow(mat)))
        if (!is.matrix(sums)) { # a single row
            sums <- matrix(sums, nrow = nrow(mat), dimnames = list(NULL, names(col_groups)))
        }
        rownames(sums) <- rownames(mat)
        return(sums)
    }
    mat <- mat[, columns, drop = FALSE]
    sapply(split(colnames(mat), groups), function(cols) rowSums(mat[, cols, drop = FALSE]))
}

get_cell_types_egc <- function(cell_types, metacell_types, dataset, mat = NULL, projected = FALSE, atlas = FALSE, corrected = FALSE) {
    mat <- mat %||% get_cell_types_mat(cell_types, metacell_types, dataset, projected = projected, atlas = atlas, corrected = corrected)

    ct_egc <- t(t(mat) / colSums(mat))

    return(ct_egc)
}

# The weights of the metacells in each sample: a matrix of all the metacells of mc_mat x the samples, with the
# fraction of the sample's cells (of the chosen cell types) in each metacell of these types, and zero rows for
# the other metacells. Samples without cells of the cell types are left out, as are metacells without cells
# in the cell metadata (which have no column in the samples' fractions).
get_samples_mc_weights <- function(cell_types, metacell_types, dataset, samples = NULL) {
    mc_names <- colnames(get_mc_data(dataset, "mc_mat"))
    samp_mc_frac <- get_samp_mc_frac(dataset)

    metacells <- metacell_types$metacell[metacell_types$cell_type %in% cell_types]
    metacells <- unique(metacells[metacells %in% colnames(samp_mc_frac) & metacells %in% mc_names])

    samp_mc_frac <- samp_mc_frac[samples %||% TRUE, metacells, drop = FALSE]
    samp_mc_frac <- samp_mc_frac[rowSums(samp_mc_frac) > 0, , drop = FALSE]
    samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)

    weights <- matrix(0, nrow = length(mc_names), ncol = nrow(samp_mc_frac), dimnames = list(mc_names, rownames(samp_mc_frac)))
    weights[match(metacells, mc_names), ] <- t(samp_mc_frac)

    return(weights)
}

# The UMIs of the samples: each sample's metacells of the chosen cell types, weighted by its cells in them. One
# product with the whole matrix (zero weights for the other metacells): subsetting its columns first copied them.
get_samples_mat <- function(cell_types, metacell_types, dataset, samples = NULL) {
    weights <- get_samples_mc_weights(cell_types, metacell_types, dataset, samples = samples)

    samp_mat <- get_mc_data(dataset, "mc_mat") %*% weights

    return(samp_mat)
}

get_samples_egc <- function(cell_types, metacell_types, dataset, samples = NULL) {
    weights <- get_samples_mc_weights(cell_types, metacell_types, dataset, samples = samples)

    # egc %*% weights, folding the normalization into the weights instead of normalizing the matrix
    used <- which(rowSums(weights) > 0)
    weights[used, ] <- weights[used, , drop = FALSE] / get_mc_sum(dataset)[used]
    samp_egc <- get_mc_data(dataset, "mc_mat") %*% weights

    return(samp_egc)
}

get_samples_gene_egc <- function(gene, dataset, metacells = NULL) {
    g <- get_gene_egc(gene, dataset)
    samp_mc_frac <- get_samp_mc_frac(dataset)
    # by metacell: the fractions have a column only for the metacells with cells, in their own order
    metacells_with_cells <- intersect(metacells %||% names(g), intersect(names(g), colnames(samp_mc_frac)))
    g <- g[metacells_with_cells]
    samp_mc_frac <- samp_mc_frac[, metacells_with_cells, drop = FALSE]
    if (!is.null(metacells)) {
        samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)
        samp_mc_frac[is.na(samp_mc_frac)] <- 0
    }

    samp_egc <- samp_mc_frac %*% g
    samp_egc <- stats::setNames(samp_egc[, 1], rownames(samp_egc))

    return(samp_egc)
}
