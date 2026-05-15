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

#' Retrieve a slice of the UMIs matrix via DAF
#'
#' Retrieves required genes/metacells from the UMIs matrix.
#' Filtering is done in R after retrieval.
#'
#' @param daf_obj DAF object
#' @param genes Optional vector of gene names to retrieve
#' @param metacells Optional vector of metacell names to retrieve
#'
#' @return Matrix slice with genes as rows and metacells as columns
#' @export
daf_query_mc_mat <- function(daf_obj, genes = NULL, metacells = NULL) {
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    # For gene subsets, query per-gene vectors from DAF to avoid loading the full matrix
    if (!is.null(genes) && length(genes) > 0) {
        gene_names <- dafr::axis_entries(daf_obj, "gene")
        valid_genes <- intersect(as.character(genes), gene_names)
        if (length(valid_genes) == 0) {
            mat <- Matrix::Matrix(0, nrow = 0, ncol = length(metacell_names), sparse = TRUE)
            colnames(mat) <- metacell_names
            return(mat)
        }

        # For small gene sets, per-gene queries beat full matrix load + subset.
        # Threshold chosen empirically on the OBK dataset.
        if (length(valid_genes) <= 50) {
            # Query each gene individually - avoids loading full gene x metacell
            # matrix. dafr's query syntax requires the axis-mask to come AFTER
            # the matrix selector, e.g. `@ metacell @ gene :: UMIs @ gene = X`.
            # The earlier `@ metacell @ gene = X :: UMIs` form errors out with
            # "comparator outside of mask" on every DAF, silently sending every
            # caller through the full-matrix fallback below.
            vecs <- lapply(valid_genes, function(gene) {
                query <- glue::glue("@ metacell @ gene :: UMIs @ gene = {escape_daf_value(gene)}")
                tryCatch(daf_obj[query], error = function(e) NULL)
            })

            # Check if per-gene queries worked
            if (!any(sapply(vecs, is.null))) {
                # Build sparse matrix directly from per-gene vectors
                n_genes <- length(valid_genes)
                n_mcs <- length(metacell_names)
                numeric_vecs <- lapply(vecs, as.numeric)
                # Find non-zero entries for sparse triplet construction
                # Use list accumulation to avoid O(n²) vector growth
                i_list <- vector("list", n_genes)
                j_list <- vector("list", n_genes)
                x_list <- vector("list", n_genes)
                for (g in seq_len(n_genes)) {
                    nz <- which(numeric_vecs[[g]] != 0)
                    if (length(nz) > 0) {
                        i_list[[g]] <- rep(g, length(nz))
                        j_list[[g]] <- nz
                        x_list[[g]] <- numeric_vecs[[g]][nz]
                    }
                }
                i_idx <- unlist(i_list)
                j_idx <- unlist(j_list)
                x_vals <- unlist(x_list)
                if (length(i_idx) > 0) {
                    mat <- Matrix::sparseMatrix(
                        i = i_idx, j = j_idx, x = x_vals,
                        dims = c(n_genes, n_mcs),
                        dimnames = list(valid_genes, metacell_names)
                    )
                } else {
                    mat <- Matrix::Matrix(0, nrow = n_genes, ncol = n_mcs, sparse = TRUE)
                    rownames(mat) <- valid_genes
                    colnames(mat) <- metacell_names
                }
                if (!is.null(metacells)) {
                    valid_mc <- intersect(metacells, metacell_names)
                    mat <- mat[, valid_mc, drop = FALSE]
                }
                return(mat)
            }
            # Fall through to full matrix if per-gene queries failed
        }
        # For >50 genes, fall through to full matrix retrieval below
    }

    # For small metacell subsets (no gene filter), query per-metacell to avoid loading full matrix
    if (is.null(genes) && !is.null(metacells) && length(metacells) <= 20) {
        gene_names <- dafr::axis_entries(daf_obj, "gene")
        valid_metacells <- intersect(as.character(metacells), metacell_names)
        if (length(valid_metacells) > 0) {
            vecs <- lapply(valid_metacells, function(mc) {
                # Mask must come after the matrix selector; see the per-gene
                # fast path above for the same bug.
                query <- glue::glue("@ gene @ metacell :: UMIs @ metacell = {escape_daf_value(mc)}")
                tryCatch(daf_obj[query], error = function(e) NULL)
            })

            if (!any(sapply(vecs, is.null))) {
                mat <- do.call(cbind, lapply(vecs, as.numeric))
                rownames(mat) <- gene_names
                colnames(mat) <- valid_metacells
                return(Matrix::Matrix(mat, sparse = TRUE))
            }
        }
        # Fall through to full matrix if per-metacell queries failed
    }

    # Full matrix retrieval (needed when genes=NULL or per-gene queries failed).
    # Request gene x metacell directly - DAF handles relayout internally.
    # dafr::get_matrix() returns a named matrix when axes exist; skip redundant dimnames.
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    mc_mat <- tryCatch(
        {
            m <- dafr::get_matrix(daf_obj, "gene", "metacell", "UMIs")
            if (is.null(rownames(m))) rownames(m) <- gene_names
            if (is.null(colnames(m))) colnames(m) <- metacell_names
            m
        },
        error = function(e) {
            # Fallback: load metacell x gene and transpose
            m <- dafr::get_matrix(daf_obj, "metacell", "gene", "UMIs")
            m <- Matrix::t(m)
            if (is.null(rownames(m))) rownames(m) <- gene_names
            if (is.null(colnames(m))) colnames(m) <- metacell_names
            m
        }
    )

    if (!is.null(genes)) {
        valid_genes <- intersect(genes, gene_names)
        mc_mat <- mc_mat[valid_genes, , drop = FALSE]
    }

    if (!is.null(metacells)) {
        valid_metacells <- intersect(metacells, metacell_names)
        mc_mat <- mc_mat[, valid_metacells, drop = FALSE]
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
    # Try canonical "primary_module" first, fall back to legacy "module"
    module_prop <- NULL
    if (dafr::has_vector(daf_obj, "gene", "primary_module")) {
        module_prop <- "primary_module"
    } else if (dafr::has_vector(daf_obj, "gene", "module")) {
        module_prop <- "module"
    }
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

#' Get metacell count per type using DAF query
#'
#' @param daf_obj DAF object
#' @return Named vector of counts per type
#' @export
daf_query_type_counts <- function(daf_obj) {
    result <- daf_obj["@ metacell : type / type >> Count"]
    return(result)
}

