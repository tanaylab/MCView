#' Get expression per gene per cell (EGC) matrix
#'
#' @param dataset Dataset name
#' @param genes Optional vector of gene names to filter
#' @param atlas Whether to use atlas data
#' @param metacells Optional vector of metacell names to filter
#' @return EGC matrix (genes x metacells)
#' @export
get_mc_egc <- function(dataset, genes = NULL, atlas = FALSE, metacells = NULL) {
    # Get DAF object
    if (atlas) {
        daf_obj <- get_atlas_daf()
    } else {
        daf_obj <- get_dataset_daf(dataset)
    }

    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Use query-based matrix access for efficiency
    mc_mat <- daf_query_mc_mat(daf_obj, genes = genes, metacells = metacells, cache = TRUE)
    mc_sum <- daf_query_mc_sum(daf_obj, metacells = metacells, cache = TRUE)

    # Filter mc_sum to match metacells in matrix
    if (!is.null(metacells)) {
        mc_sum <- mc_sum[intersect(metacells, names(mc_sum))]
    }

    return(t(t(mc_mat) / mc_sum))
}

get_mc_fp <- function(dataset, genes = NULL, atlas = FALSE, metacells = NULL) {
    mc_egc <- get_mc_egc(dataset, genes = genes, atlas = atlas, metacells = metacells)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / apply(mc_egc_norm, 1, median, na.rm = TRUE)

    return(mc_fp)
}

#' Get gene modules EGC matrix
#'
#' @param dataset Dataset name
#' @param modules Optional vector of module names to filter
#' @param gene_modules Optional gene modules data frame
#' @param atlas Whether to use atlas data
#' @return EGC matrix aggregated by gene modules
#' @export
get_mc_gene_modules_egc <- function(dataset, modules = NULL, gene_modules = NULL, atlas = FALSE) {
    # Get DAF object
    if (atlas) {
        daf_obj <- get_atlas_daf()
    } else {
        daf_obj <- get_dataset_daf(dataset)
    }

    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Try DAF query first (if DAF has gene.module property)
    mod_mat <- daf_query_module_umis(daf_obj, modules = modules)
    if (!is.null(mod_mat)) {
        mc_sum <- daf_query_mc_sum(daf_obj, cache = TRUE)
        return(t(t(mod_mat) / mc_sum))
    }

    # Fallback: use gene_modules data frame and R aggregation
    gene_modules <- gene_modules %||% get_mc_data(dataset, "gene_modules", atlas = atlas)
    if (is.null(gene_modules)) {
        return(NULL)
    }

    if (!is.null(modules)) {
        gene_modules <- gene_modules %>%
            filter(module %in% modules) %>%
            mutate(module = forcats::fct_drop(module))
    }

    # Get genes for modules
    module_genes <- gene_modules$gene

    # Use query-based matrix access for module genes
    mc_mat <- daf_query_mc_mat(daf_obj, genes = module_genes, cache = TRUE)
    mc_sum <- daf_query_mc_sum(daf_obj, cache = TRUE)

    # Aggregate by module using base R (fallback)
    unique_modules <- unique(gene_modules$module)
    mod_mat <- matrix(0, nrow = length(unique_modules), ncol = ncol(mc_mat))
    rownames(mod_mat) <- unique_modules
    colnames(mod_mat) <- colnames(mc_mat)

    for (mod in unique_modules) {
        mod_genes <- gene_modules$gene[gene_modules$module == mod]
        valid_genes <- intersect(mod_genes, rownames(mc_mat))
        if (length(valid_genes) > 1) {
            mod_mat[as.character(mod), ] <- Matrix::colSums(mc_mat[valid_genes, , drop = FALSE])
        } else if (length(valid_genes) == 1) {
            mod_mat[as.character(mod), ] <- as.numeric(mc_mat[valid_genes, ])
        }
    }

    return(t(t(mod_mat) / mc_sum))
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

    mat <- mat[, intersect(colnames(mat), metacells), drop = FALSE]

    return(mat)
}

#' Get EGC for a single gene
#'
#' @param gene Gene name
#' @param dataset Dataset name
#' @param projected Use projected data (for atlas)
#' @param atlas Whether to use atlas data
#' @param corrected Use corrected data
#' @return Named vector of EGC values per metacell
#' @export
get_gene_egc <- function(gene, dataset, projected = FALSE, atlas = FALSE, corrected = FALSE) {
    # Get DAF object
    if (atlas) {
        daf_obj <- get_atlas_daf()
    } else {
        daf_obj <- get_dataset_daf(dataset)
    }

    if (is.null(daf_obj)) {
        return(NULL)
    }

    # TODO: Handle projected and corrected modes when DAF supports it
    # For now, use standard UMIs matrix
    mc_mat <- daf_query_mc_mat(daf_obj, genes = gene, cache = TRUE)
    mc_sum <- daf_query_mc_sum(daf_obj, cache = TRUE)

    # Return as vector (single gene)
    if (nrow(mc_mat) == 0) {
        return(NULL)
    }

    result <- mc_mat[1, ] / mc_sum
    names(result) <- names(mc_sum)

    return(result)
}

get_gene_fp <- function(gene, dataset, atlas = FALSE) {
    mc_egc_norm <- get_gene_egc(gene, dataset, atlas = atlas) + 1e-5

    mc_fp <- mc_egc_norm / median(mc_egc_norm, na.rm = TRUE)
    return(mc_fp)
}

#' Get EGC for specific metacells
#'
#' @param metacells Vector of metacell names
#' @param dataset Dataset name
#' @param projected Use projected data (for atlas)
#' @param atlas Whether to use atlas data
#' @param corrected Use corrected data
#' @return EGC matrix (genes x metacells)
#' @export
get_metacells_egc <- function(metacells, dataset, projected = FALSE, atlas = FALSE, corrected = FALSE) {
    # Get DAF object
    if (atlas) {
        daf_obj <- get_atlas_daf()
    } else {
        daf_obj <- get_dataset_daf(dataset)
    }

    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Use query-based matrix access for specific metacells
    mc_mat <- daf_query_mc_mat(daf_obj, metacells = metacells, cache = TRUE)
    mc_sum <- daf_query_mc_sum(daf_obj, metacells = metacells, cache = TRUE)

    mc_egc <- t(t(mc_mat) / mc_sum)
    return(mc_egc)
}

#' Get UMI matrix aggregated by cell types
#'
#' @param cell_types Vector of cell type names
#' @param metacell_types Metacell types data frame
#' @param dataset Dataset name
#' @param projected Use projected data
#' @param atlas Whether to use atlas data
#' @param corrected Use corrected data
#' @return UMI matrix (genes x cell_types)
#' @export
get_cell_types_mat <- function(cell_types, metacell_types, dataset, projected = FALSE, atlas = FALSE, corrected = FALSE) {
    # Get DAF object
    if (atlas) {
        daf_obj <- get_atlas_daf()
    } else {
        daf_obj <- get_dataset_daf(dataset)
    }

    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Use DAF query for efficient cell type aggregation
    ct_mat <- daf_query_cell_type_umis(daf_obj, cell_types = cell_types)

    return(ct_mat)
}

get_cell_types_egc <- function(cell_types, metacell_types, dataset, mat = NULL, projected = FALSE, atlas = FALSE, corrected = FALSE) {
    mat <- mat %||% get_cell_types_mat(cell_types, metacell_types, dataset, projected = projected, atlas = atlas, corrected = corrected)

    ct_egc <- t(t(mat) / colSums(mat))

    return(ct_egc)
}

#' Get samples matrix
#'
#' @param cell_types Vector of cell type names
#' @param metacell_types Metacell types data frame
#' @param dataset Dataset name
#' @return Samples UMI matrix
#' @export
get_samples_mat <- function(cell_types, metacell_types, dataset) {
    daf_obj <- get_dataset_daf(dataset)

    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Get metacells for specified cell types
    target_metacells <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        pull(metacell)

    # Use query-based matrix access for target metacells
    mc_mat <- daf_query_mc_mat(daf_obj, metacells = target_metacells, cache = TRUE)

    # Get sample-metacell fractions
    samp_mc_frac <- get_samp_mc_frac(dataset)
    if (is.null(samp_mc_frac)) {
        return(NULL)
    }

    # Intersect with available metacells
    common_metacells <- intersect(target_metacells, colnames(samp_mc_frac))
    if (length(common_metacells) == 0) {
        return(NULL)
    }

    samp_mc_frac <- samp_mc_frac[, common_metacells, drop = FALSE]
    mc_mat <- mc_mat[, common_metacells, drop = FALSE]

    samp_sums <- rowSums(samp_mc_frac)

    # remove samples with no metacells
    samp_mc_frac <- samp_mc_frac[samp_sums > 0, , drop = FALSE]
    samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)

    samp_mat <- mc_mat %*% t(samp_mc_frac)

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
    samp_egc <- stats::setNames(samp_egc[, 1], rownames(samp_egc))

    return(samp_egc)
}
