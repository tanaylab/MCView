egc_to_fp <- function(mc_egc, epsilon = 1e-5) {
    m <- mc_egc + epsilon
    sweep(m, 1, matrixStats::rowMedians(m), "/")
}

#' Get expression per gene per cell (EGC) matrix
#'
#' @param dataset Dataset name
#' @param genes Optional vector of gene names to filter
#' @param atlas Whether to use atlas data
#' @param metacells Optional vector of metacell names to filter
#' @return EGC matrix (genes x metacells)
#' @export
get_mc_egc <- function(dataset, genes = NULL, atlas = FALSE, metacells = NULL) {
    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Full unfiltered matrix: prefer a pre-computed fraction matrix
    # (`linear_fraction` then `geomean_fraction`) over rebuilding via
    # `sweep(UMIs, 2, total_UMIs, "/")`. When a fraction matrix is present
    # in the base DAF, `convert_daf_mc_egc` returns it via mmap ALTREP —
    # ~200 ms cold instead of ~2 s for the sweep. When no fraction matrix
    # is present, `convert_daf_mc_egc` falls through to a single dafr
    # `% Fraction` query.
    if (is.null(genes) && is.null(metacells) && !atlas) {
        cached <- mcv_cache_get(dataset, "mc_egc_full")
        if (!is.null(cached)) {
            return(cached)
        }
        mc_egc <- convert_daf_mc_egc(daf_obj)
        # Skip the in-R cache when the DAF carries a stored fraction matrix:
        # convert_daf_mc_egc returned an mmap-backed view, so caching the
        # materialised dense matrix in mcview_env defeats mmap and duplicates
        # dafr's own cache. For DAFs without a stored fraction we keep the
        # cache - that path computed UMIs % Fraction and the caller pays
        # repeatedly otherwise.
        has_stored_fraction <-
            dafr::has_matrix(daf_obj, "gene", "metacell", "linear_fraction") ||
            dafr::has_matrix(daf_obj, "gene", "metacell", "geomean_fraction")
        if (!has_stored_fraction) {
            mcv_cache_set(dataset, "mc_egc_full", mc_egc)
        }
        return(mc_egc)
    }

    compute_egc_from_daf(daf_obj, genes = genes, metacells = metacells)
}

#' Get the cached log-fold-point (lfp) matrix for a dataset
#'
#' Returns `log2(mc_egc_full + egc_epsilon)` with a session-level cache keyed
#' alongside `mc_egc_full`. Filtered subsets (with `genes`/`metacells`) are
#' not cached — callers should compute directly on those.
#'
#' @param dataset Dataset name
#' @param atlas Whether to use atlas data (no cache — small / rare path)
#' @return lfp matrix (genes x metacells)
#' @export
get_mc_lfp <- function(dataset, atlas = FALSE) {
    eps <- mcv_get("egc_epsilon")
    if (atlas) {
        return(dafr::fast_log(get_mc_egc(dataset, atlas = TRUE),
                              eps = eps, base = 2))
    }
    cached <- mcv_cache_get(dataset, "lfp_full")
    if (!is.null(cached)) return(cached)

    mc_egc <- get_mc_egc(dataset)
    if (is.null(mc_egc)) return(NULL)

    # dafr::fast_log dispatches to std::log2 with OpenMP — ~18x faster
    # than base-R `log2(mc_egc + eps)` on a 28K x 2400 matrix and
    # last-ULP equal (the existing lfp_full cache still amortises the
    # one-time cost across the session).
    lfp <- dafr::fast_log(mc_egc, eps = eps, base = 2)
    mcv_cache_set(dataset, "lfp_full", lfp)
    # Any cached transpose belonged to the previous lfp_full and is now
    # stale - clear it so get_mc_lfp_t() recomputes on next demand.
    mcv_cache_set(dataset, "lfp_full_t", NULL)
    lfp
}

#' Get the cached TRANSPOSE of the lfp matrix
#'
#' calc_top_cors and similar BLAS-cor callers want metacell-rows x gene-cols
#' input. Computing `t(lfp_full)` on every gene click costs ~460 ms on OBK
#' (28K x 2.4K dgeMatrix). Cache the transpose alongside `lfp_full` so each
#' gene-correlation click reuses it; ~500 MB of additional RSS in exchange.
#'
#' The cached transpose is only valid as the pair of the current lfp_full.
#' get_mc_lfp() proactively clears lfp_full_t whenever it rebuilds lfp_full,
#' but external clears via `mcv_cache_set(dataset, "lfp_full", NULL)` (e.g.
#' from the cold-path tests) bypass that. Guard against the stale-pair
#' scenario by checking lfp_full presence before trusting the transpose.
#'
#' @param dataset Dataset name
#' @return metacell x gene matrix, or NULL if lfp is unavailable
#' @noRd
get_mc_lfp_t <- function(dataset) {
    if (is.null(mcv_cache_get(dataset, "lfp_full"))) {
        # lfp_full was cleared - any cached transpose belongs to a stale
        # matrix. Drop it now so other callers don't see it either.
        mcv_cache_set(dataset, "lfp_full_t", NULL)
    } else {
        cached <- mcv_cache_get(dataset, "lfp_full_t")
        if (!is.null(cached)) return(cached)
    }

    lfp <- get_mc_lfp(dataset)
    if (is.null(lfp)) return(NULL)
    lfp_t <- t(lfp)
    mcv_cache_set(dataset, "lfp_full_t", lfp_t)
    lfp_t
}

get_mc_fp <- function(dataset, genes = NULL, atlas = FALSE, metacells = NULL) {
    mc_egc <- get_mc_egc(dataset, genes = genes, atlas = atlas, metacells = metacells)
    egc_to_fp(mc_egc)
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
    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Try DAF query first (if DAF has gene.module property). The pushed-down
    # `% Fraction` returns per-metacell module fractions directly - one DAF
    # query instead of two (mod_mat + mc_sum) + a double-transpose / divide.
    mod_mat <- daf_query_module_fraction(daf_obj, modules = modules)
    if (!is.null(mod_mat)) {
        return(mod_mat)
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
    mc_mat <- daf_query_mc_mat(daf_obj, genes = module_genes)
    mc_sum <- daf_query_mc_sum(daf_obj)

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
    return(egc_to_fp(mc_egc))
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
    if (projected && corrected) {
        stop("projected and corrected cannot both be TRUE")
    }

    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    if (projected || corrected) {
        mat_name <- if (projected) "projected_mat" else "mc_mat_corrected"
        mc_mat <- get_mc_data(dataset, mat_name, atlas = atlas)
        if (is.null(mc_mat) || is.null(rownames(mc_mat)) || !(gene %in% rownames(mc_mat))) {
            return(NULL)
        }

        mc_sum <- daf_query_mc_sum(daf_obj)
        common_metacells <- intersect(colnames(mc_mat), names(mc_sum))
        if (length(common_metacells) == 0) {
            return(NULL)
        }

        gene_mat <- mc_mat[gene, common_metacells, drop = FALSE]
        gene_vec <- as.numeric(gene_mat[1, ])
        names(gene_vec) <- common_metacells
        return(gene_vec / mc_sum[common_metacells])
    }

    # Fast path: check session-level mc_egc_full cache (instant extraction).
    # mc_mat is no longer R-cached (see utils_cache.R static_vars), so the
    # only warm slot worth checking is the EGC.
    if (!atlas) {
        cached_egc <- mcv_cache_get(dataset, "mc_egc_full")
        if (!is.null(cached_egc) && gene %in% rownames(cached_egc)) {
            return(cached_egc[gene, ])
        }
    }

    # Cold path: targeted DAF query for just this gene. ~50-700 ms cold,
    # ~30 ms warm via dafr's own per-version cache.
    egc <- compute_egc_from_daf(daf_obj, genes = gene)
    if (nrow(egc) == 0) {
        return(NULL)
    }
    egc[1, ]
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
    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    compute_egc_from_daf(daf_obj, metacells = metacells)
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
    daf_obj <- get_daf_for_query(dataset, atlas)
    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Group by the CURRENT metacell_types mapping (reflecting session
    # re-annotations), not the DAF's stored `type` vector. Aggregating via
    # `daf_query_cell_type_umis` (which groups by stored type and ignores the
    # projected/corrected flags) made DE/obs-exp ignore annotation edits and
    # made the obs-vs-projected comparison degenerate (both returned raw UMIs).
    if (is.null(metacell_types)) {
        return(NULL)
    }
    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()
    mc_types[] <- as.character(mc_types)
    names(mc_types) <- as.character(names(mc_types))
    if (length(mc_types) == 0) {
        return(NULL)
    }

    # Select the matrix variant. projected_mat / mc_mat_corrected are
    # gene x metacell UMI-like matrices reconstructed from the stored
    # fraction matrices; the raw path pulls only the requested metacells.
    if (projected) {
        mc_mat <- get_mc_data(dataset, "projected_mat", atlas = atlas)
    } else if (corrected) {
        mc_mat <- get_mc_data(dataset, "mc_mat_corrected", atlas = atlas)
        if (is.null(mc_mat)) {
            mc_mat <- daf_query_mc_mat(daf_obj, metacells = names(mc_types))
        }
    } else {
        mc_mat <- daf_query_mc_mat(daf_obj, metacells = names(mc_types))
    }
    if (is.null(mc_mat)) {
        return(NULL)
    }

    # Restrict to the requested metacells present in the matrix, preserving
    # their cell-type labels, then sum UMIs per gene within each cell type via
    # a sparse indicator multiply (keeps mc_mat sparse; result is small).
    keep <- intersect(names(mc_types), colnames(mc_mat))
    if (length(keep) == 0) {
        return(NULL)
    }
    mc_mat <- mc_mat[, keep, drop = FALSE]
    groups <- mc_types[keep]
    type_levels <- unique(groups)
    indicator <- Matrix::sparseMatrix(
        i = seq_along(groups),
        j = match(groups, type_levels),
        x = 1,
        dims = c(length(groups), length(type_levels)),
        dimnames = list(keep, type_levels)
    )
    ct_mat <- as.matrix(mc_mat %*% indicator)
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
#' @param group_field Optional grouping field for cell-level pseudobulk.
#'   When non-NULL and cell-level data is available, uses get_group_pseudobulk_mat
#'   instead of metacell-weighted approach.
#' @return Samples UMI matrix (genes x samples/groups)
#' @export
get_samples_mat <- function(cell_types, metacell_types, dataset, group_field = NULL) {
    # Cell-level pseudobulk path
    if (!is.null(group_field) && has_cell_gene_umis(dataset)) {
        all_genes <- gene_names(dataset)
        mat <- get_group_pseudobulk_mat(
            dataset, genes = all_genes, group_field = group_field,
            cell_types = cell_types, metacell_types = metacell_types
        )
        return(mat)
    }

    # Original metacell-weighted approach
    daf_obj <- get_dataset_daf(dataset)

    if (is.null(daf_obj)) {
        return(NULL)
    }

    # Get metacells for specified cell types
    target_metacells <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        pull(metacell)

    # Use query-based matrix access for target metacells
    mc_mat <- daf_query_mc_mat(daf_obj, metacells = target_metacells)

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

get_samples_egc <- function(cell_types, metacell_types, dataset, group_field = NULL) {
    # Cell-level pseudobulk path
    if (!is.null(group_field) && has_cell_gene_umis(dataset)) {
        mat <- get_samples_mat(cell_types, metacell_types, dataset, group_field = group_field)
        if (is.null(mat)) {
            return(NULL)
        }
        # Normalize to EGC (fraction of total UMIs per sample)
        col_sums <- colSums(mat)
        col_sums[col_sums == 0] <- 1 # avoid division by zero
        egc <- t(t(mat) / col_sums)
        return(egc)
    }

    # Original metacell-weighted approach
    mc_types <- metacell_types %>%
        filter(cell_type %in% cell_types) %>%
        select(metacell, cell_type) %>%
        deframe()

    # Only load EGC for the needed metacells instead of the full matrix
    mc_egc <- get_mc_egc(dataset, metacells = names(mc_types))

    samp_mc_frac <- get_samp_mc_frac(dataset)[, names(mc_types)]
    samp_mc_frac <- samp_mc_frac / rowSums(samp_mc_frac)

    samp_egc <- mc_egc %*% t(samp_mc_frac)

    return(samp_egc)
}

get_samples_gene_egc <- function(gene, dataset, metacells = NULL, group_field = NULL, metacell_types = NULL) {
    # Cell-level pseudobulk path
    if (!is.null(group_field) && has_cell_gene_umis(dataset)) {
        # Determine cell_types from metacells filter if provided
        cell_types_filter <- NULL
        if (!is.null(metacells) && !is.null(metacell_types)) {
            cell_types_filter <- metacell_types %>%
                dplyr::filter(metacell %in% metacells) %>%
                dplyr::pull(cell_type) %>%
                unique()
        }

        gene_umis <- get_group_gene_expression(
            dataset, gene, group_field,
            cell_types = cell_types_filter,
            metacell_types = metacell_types
        )

        # Get total UMIs per group for normalization
        # Sum across all genes would be expensive; use the group QC stats if available
        qc <- tryCatch(
            get_group_qc_stats(dataset, group_field),
            error = function(e) NULL
        )

        if (!is.null(qc) && !all(is.na(qc$total_umis))) {
            total_umis <- stats::setNames(qc$total_umis, qc$group_id)
            common <- intersect(names(gene_umis), names(total_umis))
            egc <- gene_umis[common] / total_umis[common]
        } else {
            # Fallback: return raw UMIs (not ideal but functional)
            egc <- gene_umis
        }
        return(egc)
    }

    # Original metacell-weighted approach
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
