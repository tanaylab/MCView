# daf_queries.R - Query-based DAF accessors (daf_query_*)
#
# Split from R/daf_data.R (2026-05-01). See R/daf_accessors.R for the
# safe-access primitives these queries build on, and R/daf_conversion.R
# for converters that consume query results.

# ==============================================================================
# Query-Based Data Access - Shared Helpers
# ==============================================================================

#' Query a named vector from a DAF axis
#'
#' Shared helper for the common pattern: get_vector + set names + optional filter.
#' Used by daf_query_mc_sum, daf_query_cell_types, etc.
#'
#' @param daf_obj DAF object
#' @param axis Axis name (e.g. "metacell")
#' @param property Vector property name (e.g. "total_UMIs", "type")
#' @param filter Optional character vector of axis entry names to keep
#'
#' @return Named vector with axis entries as names, optionally filtered
#' @noRd
daf_query_named_vector <- function(daf_obj, axis, property, filter = NULL) {
    # dafr::get_vector returns a named vector (names from NamedArrays)
    result <- dafr::get_vector(daf_obj, axis, property)

    if (!is.null(filter)) {
        valid_entries <- intersect(filter, names(result))
        result <- result[valid_entries]
    }

    return(result)
}

#' Query a gene-level aggregation of UMIs via DAF
#'
#' Shared helper for the common pattern: aggregate UMIs over metacells per gene.
#' Used by daf_query_gene_max_umis, daf_query_gene_mean_umis, daf_query_gene_sum_umis.
#'
#' @param daf_obj DAF object
#' @param agg_op DAF aggregation operator string (e.g. "Max", "Mean", "Sum")
#'
#' @return Named vector with gene names
#' @noRd
daf_query_gene_agg <- function(daf_obj, agg_op) {
    # `>-` (ReduceToRow) reduces along the row axis, returning one value per
    # column-axis entry. For `@ metacell @ gene :: UMIs` that's per-gene.
    # `>>` would collapse the whole matrix to a scalar (matches Julia
    # DataAxesFormats.jl semantics), which is not what callers want here.
    query <- paste0("@ metacell @ gene :: UMIs >- ", agg_op)
    result <- daf_obj[query]
    if (is.null(names(result)) && length(result) > 1) {
        names(result) <- dafr::axis_entries(daf_obj, "gene")
    }
    return(result)
}

# ==============================================================================
# Query-Based Data Access
# ==============================================================================

#' Build a `[ name = X | name = Y | ... ]` mask from a name vector.
#' Returns "" for empty input so callers can paste it inline.
#' @noRd
build_name_mask <- function(names) {
    if (length(names) == 0) return("")
    parts <- vapply(names, function(n) paste0("name = ", escape_daf_value(n)),
                    character(1))
    paste0(" [ ", paste(parts, collapse = " | "), " ]")
}

# Above this many entries, the mask query string gets unwieldy and the
# full-matrix load + R-side subset is cheaper. Empirical sweet spot ~500;
# OBK (~2.4K metacells, ~28K genes) lands well inside it.
.MASK_THRESHOLD <- 500L

#' Retrieve a slice of the UMIs matrix via DAF
#'
#' For small subsets (<= MASK_THRESHOLD entries on either axis), pushes a
#' `[ name = ... | name = ... ]` mask into the DAF query so dafr returns
#' just the requested rows/cols directly. For larger subsets falls back to
#' the full matrix + R-side index. The mask form requires dafr >= 0.2.7
#' (which fixed the row-mask + GroupRowsBy alignment bug).
#'
#' @param daf_obj DAF object
#' @param genes Optional vector of gene names to retrieve
#' @param metacells Optional vector of metacell names to retrieve
#'
#' @return Matrix slice with genes as rows and metacells as columns
#' @export
daf_query_mc_mat <- function(daf_obj, genes = NULL, metacells = NULL) {
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    valid_genes <- if (is.null(genes)) NULL else intersect(as.character(genes), gene_names)
    valid_mcs   <- if (is.null(metacells)) NULL else intersect(as.character(metacells), metacell_names)

    if (!is.null(valid_genes) && length(valid_genes) == 0) {
        mat <- Matrix::Matrix(0, nrow = 0, ncol = length(metacell_names), sparse = TRUE)
        colnames(mat) <- metacell_names
        return(mat)
    }
    if (!is.null(valid_mcs) && length(valid_mcs) == 0) {
        mat <- Matrix::Matrix(0, nrow = length(gene_names), ncol = 0, sparse = TRUE)
        rownames(mat) <- gene_names
        return(mat)
    }

    use_gene_mask <- !is.null(valid_genes) && length(valid_genes) <= .MASK_THRESHOLD
    use_mc_mask   <- !is.null(valid_mcs)   && length(valid_mcs)   <= .MASK_THRESHOLD

    # Single-query fast path for any combination of small subsets.
    if (use_gene_mask || use_mc_mask) {
        gmask <- if (use_gene_mask) build_name_mask(valid_genes) else ""
        mmask <- if (use_mc_mask) build_name_mask(valid_mcs) else ""
        # dafr::get_query rather than daf_obj[...]: the `[` dispatch for
        # chained DAFs from inside this package's namespace can mis-resolve
        # and raise "comparator outside of mask".
        query <- as.character(glue::glue("@ gene{gmask} @ metacell{mmask} :: UMIs"))
        mat <- tryCatch(dafr::get_query(daf_obj, query), error = function(e) NULL)
        if (!is.null(mat)) {
            # Reorder both axes to the REQUESTED order. A DAF name-mask
            # preserves the native axis order (it does not reorder to the
            # mask's listing), and the un-masked axis still needs its R-side
            # filter applied. Reordering by name in all cases keeps this fast
            # path consistent with the full-matrix fallback below, which
            # always returns requested order - positional consumers (e.g. the
            # mc_egc %*% t(samp_mc_frac) product in get_samples_egc) depend on
            # it. The masked subset is small, so the reindex is cheap.
            if (!is.null(valid_genes)) {
                mat <- mat[valid_genes, , drop = FALSE]
            }
            if (!is.null(valid_mcs)) {
                mat <- mat[, valid_mcs, drop = FALSE]
            }
            return(if (methods::is(mat, "Matrix")) mat
                   else Matrix::Matrix(mat, sparse = TRUE))
        }
        # Fall through to full-matrix load if the masked query failed for any
        # reason (e.g. an old dafr without the GroupBy fix).
    }

    # Full-matrix retrieval. Request gene x metacell directly - DAF handles
    # relayout internally; fall back to metacell x gene + transpose if the
    # native orientation isn't stored.
    mc_mat <- tryCatch(
        {
            m <- dafr::get_matrix(daf_obj, "gene", "metacell", "UMIs")
            if (is.null(rownames(m))) rownames(m) <- gene_names
            if (is.null(colnames(m))) colnames(m) <- metacell_names
            m
        },
        error = function(e) {
            m <- dafr::get_matrix(daf_obj, "metacell", "gene", "UMIs")
            m <- Matrix::t(m)
            if (is.null(rownames(m))) rownames(m) <- gene_names
            if (is.null(colnames(m))) colnames(m) <- metacell_names
            m
        }
    )

    if (!is.null(valid_genes)) {
        mc_mat <- mc_mat[valid_genes, , drop = FALSE]
    }
    if (!is.null(valid_mcs)) {
        mc_mat <- mc_mat[, valid_mcs, drop = FALSE]
    }

    return(mc_mat)
}

#' Retrieve total UMIs per metacell via DAF
#'
#' @param daf_obj DAF object
#' @param metacells Optional vector of metacell names to retrieve
#'
#' @return Named vector of total UMIs per metacell
#' @export
daf_query_mc_sum <- function(daf_obj, metacells = NULL) {
    daf_query_named_vector(daf_obj, "metacell", "total_UMIs", filter = metacells)
}

#' Retrieve cell type assignments via DAF
#'
#' @param daf_obj DAF object
#' @param metacells Optional vector of metacell names to retrieve
#'
#' @return Named vector of cell types per metacell
#' @export
daf_query_cell_types <- function(daf_obj, metacells = NULL) {
    daf_query_named_vector(daf_obj, "metacell", "type", filter = metacells)
}

#' Retrieve 2D coordinates via DAF
#'
#' @param daf_obj DAF object
#' @param metacells Optional vector of metacell names to retrieve
#'
#' @return Data frame with metacell, x, y columns
#' @export
daf_query_2d_coords <- function(daf_obj, metacells = NULL) {
    # Get axis entries
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    # Try coordinate names in priority order: x/y, then umap_x/umap_y, then u/v
    x_coords <- try_daf_names(daf_obj, "metacell", c("x", "umap_x", "u"), required = TRUE)
    y_coords <- try_daf_names(daf_obj, "metacell", c("y", "umap_y", "v"), required = TRUE)

    # Create data frame
    result <- tibble(
        metacell = metacell_names,
        x = x_coords,
        y = y_coords
    )

    # Filter by metacells if specified
    if (!is.null(metacells)) {
        result <- result %>%
            filter(metacell %in% metacells)
    }

    return(result)
}

# ==============================================================================
# DAF Aggregation Queries
# ==============================================================================

#' Build a DAF mask alternation for cell-type filtering, e.g.
#' `[ type = A | type = B | type = C ]`. Returns empty string for NULL/empty
#' inputs so callers can paste it inline unconditionally.
#' @noRd
build_type_mask <- function(cell_types) {
    if (is.null(cell_types) || length(cell_types) == 0) return("")
    parts <- vapply(cell_types, function(t) paste0("type = ", escape_daf_value(t)),
                    character(1))
    paste0(" [ ", paste(parts, collapse = " | "), " ]")
}

#' Retrieve UMIs matrix aggregated by cell type via DAF query
#'
#' Uses DAF's GroupBy query to efficiently aggregate UMIs by cell type. When
#' `cell_types` is given, the filter is pushed into the metacell mask so the
#' GroupBy only aggregates the requested types instead of all-then-filter.
#' Requires dafr >= 0.2.6 (which fixes the row-mask + GroupRowsBy alignment
#' bug); on earlier versions the masked GroupBy mis-labels output rows.
#'
#' @param daf_obj DAF object
#' @param cell_types Optional vector of cell type names to include
#'
#' @return Matrix with genes as rows and cell types as columns
#' @export
daf_query_cell_type_umis <- function(daf_obj, cell_types = NULL) {
    mask <- build_type_mask(cell_types)
    ct_mat <- daf_obj[glue::glue("@ metacell{mask} @ gene :: UMIs -/ type >- Sum")]

    # Transpose to get genes as rows, cell types as columns
    ct_mat <- t(ct_mat)

    # Belt-and-braces filter: drop any types that survived the mask but
    # weren't requested (defensive against name-typo callers).
    if (!is.null(cell_types)) {
        valid_types <- intersect(cell_types, colnames(ct_mat))
        if (length(valid_types) > 0) {
            ct_mat <- ct_mat[, valid_types, drop = FALSE]
        }
    }

    return(ct_mat)
}

#' Retrieve total UMIs aggregated by cell type via DAF query
#'
#' Uses DAF's GroupBy query to efficiently aggregate total UMIs by cell type.
#'
#' @param daf_obj DAF object
#' @param cell_types Optional vector of cell type names to include
#'
#' @return Named vector of total UMIs per cell type
#' @export
daf_query_cell_type_sum <- function(daf_obj, cell_types = NULL) {
    mask <- build_type_mask(cell_types)
    result <- daf_obj[glue::glue("@ metacell{mask} : total_UMIs / type >> Sum")]

    if (!is.null(cell_types)) {
        valid_types <- intersect(cell_types, names(result))
        if (length(valid_types) > 0) {
            result <- result[valid_types]
        }
    }

    return(result)
}

#' Get max UMIs per gene using DAF query
#'
#' @param daf_obj DAF object
#' @return Named vector of max UMIs per gene
#' @export
daf_query_gene_max_umis <- function(daf_obj) {
    daf_query_gene_agg(daf_obj, "Max")
}

#' Get mean UMIs per gene using DAF query
#'
#' @param daf_obj DAF object
#' @return Named vector of mean UMIs per gene
#' @export
daf_query_gene_mean_umis <- function(daf_obj) {
    daf_query_gene_agg(daf_obj, "Mean")
}

#' Get sum UMIs per gene using DAF query
#'
#' @param daf_obj DAF object
#' @return Named vector of sum UMIs per gene
#' @export
daf_query_gene_sum_umis <- function(daf_obj) {
    daf_query_gene_agg(daf_obj, "Sum")
}

#' Get marker genes using DAF query
#'
#' Returns genes where is_marker is TRUE
#'
#' @param daf_obj DAF object
#' @return Character vector of marker gene names
#' @export
daf_query_marker_genes <- function(daf_obj) {
    daf_query_flagged_genes(daf_obj, "is_marker")
}

#' Get lateral genes using DAF query
#'
#' Returns genes where is_lateral is TRUE
#'
#' @param daf_obj DAF object
#' @return Character vector of lateral gene names
#' @export
daf_query_lateral_genes <- function(daf_obj) {
    daf_query_flagged_genes(daf_obj, "is_lateral")
}

#' Get noisy genes using DAF query
#'
#' Returns genes where is_noisy is TRUE
#'
#' @param daf_obj DAF object
#' @return Character vector of noisy gene names
#' @export
daf_query_noisy_genes <- function(daf_obj) {
    daf_query_flagged_genes(daf_obj, "is_noisy")
}

#' Resolve gene-axis module property name.
#'
#' Prefers the canonical `primary_module` and falls back to the legacy `module`.
#' Returns NULL when neither is present.
#' @noRd
.daf_module_prop <- function(daf_obj) {
    if (dafr::has_vector(daf_obj, "gene", "primary_module")) {
        "primary_module"
    } else if (dafr::has_vector(daf_obj, "gene", "module")) {
        "module"
    } else {
        NULL
    }
}

#' Get UMIs aggregated by gene module using DAF query
#'
#' If DAF has gene.module property, uses GroupBy for efficient aggregation.
#' Otherwise returns NULL.
#'
#' @param daf_obj DAF object
#' @param modules Optional vector of module names to filter
#' @return Matrix with modules as rows and metacells as columns, or NULL
#' @export
daf_query_module_umis <- function(daf_obj, modules = NULL) {
    module_prop <- .daf_module_prop(daf_obj)
    if (is.null(module_prop)) {
        return(NULL)
    }

    # Use DAF query for efficient grouping by module
    # GroupRowsBy module property (gene axis) and reduce columns with Sum
    query <- glue::glue("@ gene @ metacell :: UMIs -/ {module_prop} >- Sum")
    mod_mat <- tryCatch(
        {
            daf_obj[query]
        },
        error = function(e) {
            return(NULL)
        }
    )

    if (is.null(mod_mat)) {
        return(NULL)
    }

    # DAF query "@ gene @ metacell :: UMIs -/ module >- Sum" returns
    # module (rows) x metacell (columns) -- already in the desired layout.
    # No transpose needed.

    # Filter by modules if specified
    if (!is.null(modules)) {
        valid_modules <- intersect(modules, rownames(mod_mat))
        if (length(valid_modules) > 0) {
            mod_mat <- mod_mat[valid_modules, , drop = FALSE]
        }
    }

    return(mod_mat)
}

#' Get per-metacell module fractions using a single DAF query
#'
#' Returns a matrix of `sum_{g in module} UMIs[g,m] / total_UMIs[m]`. This is
#' the normalized form most callers actually want; computing it via
#' `daf_query_module_umis()` + `t(t(.) / mc_sum)` requires two queries and a
#' double-transpose. Pushing `% Fraction` into the query before GroupBy gives
#' the same answer in one pass (verified bit-equivalent on a synthetic 6x3
#' fixture).
#'
#' Returns NULL if the DAF has no gene-axis module property.
#'
#' Internal helper - if needed externally, add @export and re-run
#' devtools::document().
#'
#' @param daf_obj DAF object
#' @param modules Optional vector of module names to filter
#' @return Matrix with modules as rows and metacells as columns, or NULL
#' @noRd
daf_query_module_fraction <- function(daf_obj, modules = NULL) {
    module_prop <- .daf_module_prop(daf_obj)
    if (is.null(module_prop)) {
        return(NULL)
    }

    query <- glue::glue(
        "@ gene @ metacell :: UMIs % Fraction -/ {module_prop} >- Sum"
    )
    mod_mat <- tryCatch(daf_obj[query], error = function(e) NULL)
    if (is.null(mod_mat)) {
        return(NULL)
    }

    if (!is.null(modules)) {
        valid_modules <- intersect(modules, rownames(mod_mat))
        if (length(valid_modules) > 0) {
            mod_mat <- mod_mat[valid_modules, , drop = FALSE]
        }
    }

    return(mod_mat)
}

#' Get metacell count per type using DAF query
#'
#' @param daf_obj DAF object
#' @return Named vector of counts per type
#' @export
daf_query_type_counts <- function(daf_obj) {
    result <- daf_obj["@ metacell : type / type >> Count"]
    return(result)
}

