# daf_data.R - Data access and conversion for MCView
# Consolidates: daf_wrappers.R, daf_conversions.R

# ==============================================================================
# Safe DAF Accessors
# ==============================================================================

#' Safe DAF vector accessor
#' @param daf_obj DAF object
#' @param axis Axis name
#' @param name Property name
#' @param required If TRUE, throws error when missing
#' @return Vector or NULL
#' @export
daf_vec <- function(daf_obj, axis, name, required = TRUE) {
    if (!dafr::has_vector(daf_obj, axis, name)) {
        if (required) {
            cli_abort("Missing required vector: {axis}.{name}")
        }
        return(NULL)
    }
    dafr::get_vector(daf_obj, axis, name)
}

#' Safe DAF matrix accessor
#' @param daf_obj DAF object
#' @param axis1 First axis name
#' @param axis2 Second axis name
#' @param name Matrix name
#' @param required If TRUE, throws error when missing
#' @return Matrix or NULL
#' @export
daf_mat <- function(daf_obj, axis1, axis2, name, required = TRUE) {
    if (!dafr::has_matrix(daf_obj, axis1, axis2, name)) {
        if (required) {
            cli_abort("Missing required matrix: {axis1},{axis2}.{name}")
        }
        return(NULL)
    }
    dafr::get_matrix(daf_obj, axis1, axis2, name)
}

#' Safe DAF scalar accessor
#' @param daf_obj DAF object
#' @param name Scalar name
#' @param default Default value if missing
#' @return Scalar value or default
#' @export
daf_scalar <- function(daf_obj, name, default = NULL) {
    if (!dafr::has_scalar(daf_obj, name)) {
        return(default)
    }
    dafr::get_scalar(daf_obj, name)
}

escape_daf_value <- function(value) {
    gsub("([^A-Za-z0-9_+\\-.])", "\\\\\\1", value, perl = TRUE)
}

# ==============================================================================
# Helper Functions for Reducing Code Duplication
# ==============================================================================

#' Query genes by a boolean flag vector
#'
#' Generic function to query genes where a boolean flag is TRUE.
#' Replaces daf_query_marker_genes, daf_query_lateral_genes, daf_query_noisy_genes.
#'
#' @param daf_obj DAF object
#' @param flag_name Name of the boolean flag vector (e.g., "is_marker", "is_lateral")
#' @return Character vector of gene names where flag is TRUE
#' @export
daf_query_flagged_genes <- function(daf_obj, flag_name) {
    if (!dafr::has_vector(daf_obj, "gene", flag_name)) {
        return(character(0))
    }
    tryCatch(
        {
            query <- glue::glue("/ gene & {flag_name} : name")
            result <- daf_obj[query]
            return(result)
        },
        error = function(e) {
            flags <- dafr::get_vector(daf_obj, "gene", flag_name)
            gene_names <- dafr::axis_entries(daf_obj, "gene")
            return(gene_names[flags])
        }
    )
}

#' Convert DAF flagged genes to character vector
#'
#' Generic function to get genes where a boolean flag is TRUE.
#' Replaces convert_daf_lateral_genes, convert_daf_noisy_genes.
#'
#' @param daf_obj DAF object
#' @param flag_name Name of the boolean flag vector
#' @return Character vector of gene names where flag is TRUE
#' @export
convert_daf_flagged_genes <- function(daf_obj, flag_name) {
    flags <- daf_vec(daf_obj, "gene", flag_name, required = FALSE)
    if (is.null(flags)) {
        return(character(0))
    }
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    gene_names[flags]
}

#' Add an optional vector to a tibble
#'
#' Helper to reduce repeated pattern of checking and adding optional vectors.
#'
#' @param tbl Tibble to add vector to
#' @param daf_obj DAF object
#' @param axis Axis name
#' @param vec_name Vector name
#' @param col_name Column name in result (defaults to vec_name)
#' @return Modified tibble with vector added if it exists
#' @export
add_optional_vec <- function(tbl, daf_obj, axis, vec_name, col_name = vec_name) {
    vec <- daf_vec(daf_obj, axis, vec_name, required = FALSE)
    if (!is.null(vec)) {
        tbl[[col_name]] <- vec
    }
    tbl
}

#' Add multiple optional vectors to a tibble
#'
#' Helper to add multiple optional vectors in one call.
#'
#' @param tbl Tibble to add vectors to
#' @param daf_obj DAF object
#' @param axis Axis name
#' @param vec_names Character vector of vector names
#' @return Modified tibble with existing vectors added
#' @export
add_optional_vecs <- function(tbl, daf_obj, axis, vec_names) {
    for (vec_name in vec_names) {
        tbl <- add_optional_vec(tbl, daf_obj, axis, vec_name)
    }
    tbl
}

#' Add optional vector with fallback name
#'
#' Tries primary name first, then fallback if primary is missing.
#'
#' @param tbl Tibble to add vector to
#' @param daf_obj DAF object
#' @param axis Axis name
#' @param primary Primary vector name (also used as column name)
#' @param fallback Fallback vector name to try if primary is missing
#' @return Modified tibble with vector added if either name exists
#' @export
add_optional_vec_with_fallback <- function(tbl, daf_obj, axis, primary, fallback) {
    vec <- daf_vec(daf_obj, axis, primary, required = FALSE)
    if (is.null(vec)) {
        vec <- daf_vec(daf_obj, axis, fallback, required = FALSE)
    }
    if (!is.null(vec)) {
        tbl[[primary]] <- vec
    }
    tbl
}

#' Get DAF object for query (dataset or atlas)
#'
#' Consolidates the repeated pattern of selecting between dataset and atlas DAF.
#'
#' @param dataset Dataset name
#' @param atlas Whether to use atlas data (default: FALSE)
#' @return DAF object or NULL if not available
#' @export
get_daf_for_query <- function(dataset, atlas = FALSE) {
    if (atlas) {
        get_atlas_daf()
    } else {
        get_dataset_daf(dataset)
    }
}

#' Filter genes by lateral and noisy flags
#'
#' Consolidates the repeated pattern of filtering a gene dataframe
#' by lateral and noisy gene flags.
#'
#' @param df Data frame with a 'gene' column
#' @param lateral_genes Vector of lateral gene names (or NULL)
#' @param noisy_genes Vector of noisy gene names (or NULL)
#' @param include_lateral Whether to include lateral genes (default: TRUE)
#' @param include_noisy Whether to include noisy genes (default: TRUE)
#' @param gene_col Name of the gene column (default: "gene")
#' @return Filtered data frame
#' @export
filter_genes_by_flags <- function(df, lateral_genes = NULL, noisy_genes = NULL,
                                  include_lateral = TRUE, include_noisy = TRUE,
                                  gene_col = "gene") {
    if (!include_lateral && !is.null(lateral_genes) && length(lateral_genes) > 0) {
        df <- df[!(df[[gene_col]] %in% lateral_genes), , drop = FALSE]
    }
    if (!include_noisy && !is.null(noisy_genes) && length(noisy_genes) > 0) {
        df <- df[!(df[[gene_col]] %in% noisy_genes), , drop = FALSE]
    }
    df
}

#' Compute EGC (expression per gene per cell) from DAF
#'
#' Retrieves UMI matrix and total UMIs, then normalizes to get EGC.
#' Consolidates the repeated pattern of:
#'   mc_mat <- daf_query_mc_mat(...)
#'   mc_sum <- daf_query_mc_sum(...)
#'   return(t(t(mc_mat) / mc_sum))
#'
#' @param daf_obj DAF object
#' @param genes Optional vector of gene names to filter
#' @param metacells Optional vector of metacell names to filter
#' @return EGC matrix (genes x metacells) with columns summing to 1
#' @export
compute_egc_from_daf <- function(daf_obj, genes = NULL, metacells = NULL) {
    # NOTE: Julia EGC matrix path disabled -- JuliaCall serialization overhead
    # for a full 28K x 2.4K dense matrix (~68M elements) is much slower than
    # the R path using per-gene DAF queries + sparse matrix construction.
    # The Julia path could be re-enabled once JuliaCall has efficient shared-
    # memory transfer or for very small gene subsets.

    mc_mat <- daf_query_mc_mat(daf_obj, genes = genes, metacells = metacells)
    mc_sum <- daf_query_mc_sum(daf_obj, metacells = metacells)

    # Filter mc_sum to match metacells in matrix
    if (!is.null(metacells)) {
        mc_sum <- mc_sum[intersect(metacells, names(mc_sum))]
    }

    t(t(mc_mat) / mc_sum)
}

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
    query <- paste0("/ metacell / gene : UMIs %> ", agg_op)
    result <- daf_obj[query]
    names(result) <- dafr::axis_entries(daf_obj, "gene")
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

        # For large gene sets (>50), per-gene JuliaCall overhead (~1-3ms each) exceeds
        # full matrix load + subset cost (~100ms). Only use per-gene path for small sets.
        if (length(valid_genes) <= 50) {
            # Query each gene individually - avoids loading full gene x metacell matrix
            vecs <- lapply(valid_genes, function(gene) {
                query <- glue::glue("/ metacell / gene = {escape_daf_value(gene)} : UMIs")
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
                query <- glue::glue("/ gene / metacell = {escape_daf_value(mc)} : UMIs")
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

    # Full matrix retrieval (needed when genes=NULL or per-gene queries failed)
    # Request gene x metacell directly - DAF handles relayout internally
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    mc_mat <- tryCatch(
        {
            m <- dafr::get_matrix(daf_obj, "gene", "metacell", "UMIs")
            rownames(m) <- gene_names
            colnames(m) <- metacell_names
            m
        },
        error = function(e) {
            # Fallback: load metacell x gene and transpose
            m <- dafr::get_matrix(daf_obj, "metacell", "gene", "UMIs")
            m <- Matrix::t(m)
            rownames(m) <- gene_names
            colnames(m) <- metacell_names
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

    # Try x,y coordinates first
    if (dafr::has_vector(daf_obj, "metacell", "x") && dafr::has_vector(daf_obj, "metacell", "y")) {
        x_coords <- dafr::get_vector(daf_obj, "metacell", "x")
        y_coords <- dafr::get_vector(daf_obj, "metacell", "y")
    } else {
        # Fall back to u,v coordinates
        x_coords <- dafr::get_vector(daf_obj, "metacell", "u")
        y_coords <- dafr::get_vector(daf_obj, "metacell", "v")
    }

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

#' Retrieve UMIs matrix aggregated by cell type via DAF query
#'
#' Uses DAF's GroupBy query to efficiently aggregate UMIs by cell type.
#'
#' @param daf_obj DAF object
#' @param cell_types Optional vector of cell type names to include
#'
#' @return Matrix with genes as rows and cell types as columns
#' @export
daf_query_cell_type_umis <- function(daf_obj, cell_types = NULL) {
    # Use DAF query for efficient grouping
    # Query: / metacell / gene : UMIs @ type %> Sum
    # This groups by type and sums the UMIs
    ct_mat <- daf_obj["/ metacell / gene : UMIs @ type %> Sum"]

    # Transpose to get genes as rows, cell types as columns
    ct_mat <- t(ct_mat)

    # Filter by cell types if specified
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
    # Use DAF query for efficient grouping
    # Query: / metacell : total_UMIs @ type %> Sum
    result <- daf_obj["/ metacell : total_UMIs @ type %> Sum"]

    # Filter by cell types if specified
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
    if (!dafr::has_vector(daf_obj, "gene", "module")) {
        return(NULL)
    }

    # Use DAF query for efficient grouping by module
    # Query: / gene / metacell : UMIs @ module %> Sum
    mod_mat <- tryCatch(
        {
            daf_obj["/ gene / metacell : UMIs @ module %> Sum"]
        },
        error = function(e) {
            return(NULL)
        }
    )

    if (is.null(mod_mat)) {
        return(NULL)
    }

    # DAF query "/ gene / metacell : UMIs @ module %> Sum" returns
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
    result <- daf_obj["/ metacell : type @ type %> Count"]
    return(result)
}

# ==============================================================================
# DAF to MCView Conversion
# ==============================================================================

#' Convert DAF data to MCView format
#' @param daf_obj DAF object
#' @param var_name Variable name to convert
#' @param atlas Whether to use atlas data
#' @return Converted data in MCView format
#' @export
convert_daf_to_mcview <- function(daf_obj, var_name, atlas = FALSE) {
    switch(var_name,
        "mc_mat" = convert_daf_mc_mat(daf_obj),
        "mc_sum" = convert_daf_mc_sum(daf_obj),
        "mc2d" = convert_daf_mc2d(daf_obj),
        "metacell_types" = convert_daf_metacell_types(daf_obj),
        "cell_type_colors" = convert_daf_cell_type_colors(daf_obj),
        "metadata" = convert_daf_metadata(daf_obj),
        "lateral_genes" = convert_daf_lateral_genes(daf_obj),
        "noisy_genes" = convert_daf_noisy_genes(daf_obj),
        "marker_genes" = convert_daf_marker_genes(daf_obj),
        "gene_modules" = convert_daf_gene_modules(daf_obj),
        "inner_fold_mat" = convert_daf_inner_fold_mat(daf_obj),
        "inner_stdev_mat" = convert_daf_inner_stdev_mat(daf_obj),
        "mc_qc_metadata" = convert_daf_mc_qc_metadata(daf_obj),
        "gene_qc" = convert_daf_gene_qc(daf_obj),
        "gg_mc_top_cor" = convert_daf_gg_mc_top_cor(daf_obj),
        "gene_zero_fold" = convert_daf_gene_zero_fold(daf_obj),
        "metacell_graphs" = convert_daf_metacell_graphs(daf_obj),
        "qc_stats" = convert_daf_qc_stats(daf_obj),
        "cell_metadata" = convert_daf_cell_metadata(daf_obj),
        # Projection-related conversions
        "projected_fold" = convert_daf_projected_fold(daf_obj),
        "mc_mat_corrected" = convert_daf_mc_mat_corrected(daf_obj),
        "projected_mat" = convert_daf_projected_mat(daf_obj),
        "proj_weights" = convert_daf_proj_weights(daf_obj),
        "query_atlas_cell_type_fracs" = convert_daf_query_cell_type_fracs(daf_obj),
        "query_md" = convert_daf_query_metadata(daf_obj),
        NULL # Return NULL for unsupported data types
    )
}

# ==============================================================================
# Core Conversion Functions
# ==============================================================================

convert_daf_mc_mat <- function(daf_obj) {
    umat <- daf_mat(daf_obj, "metacell", "gene", "UMIs")
    Matrix::t(umat) # Lazy transpose, keeps sparse
}

convert_daf_mc_sum <- function(daf_obj) {
    mc_sum <- daf_vec(daf_obj, "metacell", "total_UMIs")
    names(mc_sum) <- dafr::axis_entries(daf_obj, "metacell")
    mc_sum
}

convert_daf_mc_egc <- function(daf_obj) {
    # Get UMI matrix (gene x metacell after transpose)
    mc_mat <- convert_daf_mc_mat(daf_obj)
    mc_sum <- convert_daf_mc_sum(daf_obj)

    # Compute EGC (normalized expression) - gene x metacell
    # Divide each column by its total to get fractions summing to 1
    # Matches compute_egc_from_daf() and original MCView normalization
    t(t(mc_mat) / mc_sum)
}

convert_daf_mc2d <- function(daf_obj) {
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    # Try different coordinate systems
    x_coords <- daf_vec(daf_obj, "metacell", "x", required = FALSE)
    y_coords <- daf_vec(daf_obj, "metacell", "y", required = FALSE)

    if (is.null(x_coords) || is.null(y_coords)) {
        x_coords <- daf_vec(daf_obj, "metacell", "u")
        y_coords <- daf_vec(daf_obj, "metacell", "v")
    }

    # Return as list with named vectors (expected format by mc2d_to_df)
    mc_x <- setNames(x_coords, metacell_names)
    mc_y <- setNames(y_coords, metacell_names)

    graph <- NULL
    if (dafr::has_axis(daf_obj, "metacell_graph")) {
        graph <- convert_daf_metacell_graphs(daf_obj)
        if (!is.null(graph) && "metacell" %in% names(graph)) {
            graph <- graph[["metacell"]]
        } else {
            graph <- NULL
        }
    }

    list(
        mc_id = metacell_names,
        mc_x = mc_x,
        mc_y = mc_y,
        graph = graph
    )
}

convert_daf_metacell_types <- function(daf_obj) {
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")
    cell_types <- daf_vec(daf_obj, "metacell", "type")

    mc_types <- tibble(
        metacell = metacell_names,
        cell_type = cell_types
    )

    # Add optional cell count
    n_cell <- daf_vec(daf_obj, "metacell", "n_cell", required = FALSE)
    if (!is.null(n_cell)) {
        mc_types$n_cell <- n_cell
    }

    # Add colors from metacell vector or type axis
    mc_col <- daf_vec(daf_obj, "metacell", "mc_col", required = FALSE)
    if (!is.null(mc_col)) {
        mc_types$mc_col <- mc_col
    } else {
        type_colors <- daf_vec(daf_obj, "type", "color", required = FALSE)
        if (!is.null(type_colors)) {
            type_names <- dafr::axis_entries(daf_obj, "type")
            color_map <- setNames(type_colors, type_names)
            mc_types$mc_col <- color_map[cell_types]
        }
    }

    top1_gene <- daf_vec(daf_obj, "metacell", "top1_gene", required = FALSE)
    top2_gene <- daf_vec(daf_obj, "metacell", "top2_gene", required = FALSE)
    top1_lfp <- daf_vec(daf_obj, "metacell", "top1_lfp", required = FALSE)
    top2_lfp <- daf_vec(daf_obj, "metacell", "top2_lfp", required = FALSE)

    if (is.null(top1_gene) || is.null(top2_gene) ||
        is.null(top1_lfp) || is.null(top2_lfp)) {
        top1_gene <- top1_gene %||% daf_vec(daf_obj, "metacell", "mcview_cache_top1_gene", required = FALSE)
        top2_gene <- top2_gene %||% daf_vec(daf_obj, "metacell", "mcview_cache_top2_gene", required = FALSE)
        top1_lfp <- top1_lfp %||% daf_vec(daf_obj, "metacell", "mcview_cache_top1_lfp", required = FALSE)
        top2_lfp <- top2_lfp %||% daf_vec(daf_obj, "metacell", "mcview_cache_top2_lfp", required = FALSE)
    }

    if (!is.null(top1_gene) && !is.null(top2_gene) &&
        !is.null(top1_lfp) && !is.null(top2_lfp)) {
        mc_types$top1_gene <- as.character(top1_gene)
        mc_types$top2_gene <- as.character(top2_gene)
        mc_types$top1_lfp <- as.numeric(top1_lfp)
        mc_types$top2_lfp <- as.numeric(top2_lfp)
    } else {
        # Compute top genes from expression matrix
        mc_egc <- tryCatch(
            {
                convert_daf_mc_egc(daf_obj)
            },
            error = function(e) NULL
        )

        if (!is.null(mc_egc) && ncol(mc_egc) > 0) {
            # Vectorized top 2 genes per metacell using max.col
            gene_names_egc <- rownames(mc_egc)
            egc_t <- t(as.matrix(mc_egc))

            # Top 1: find column index of max value per metacell
            top1_idx <- max.col(egc_t, ties.method = "first")
            mc_types$top1_gene <- gene_names_egc[top1_idx]
            mc_types$top1_lfp <- log2(egc_t[cbind(seq_len(nrow(egc_t)), top1_idx)] + 1e-5)

            # Top 2: mask top1 values, find next max
            egc_t[cbind(seq_len(nrow(egc_t)), top1_idx)] <- -Inf
            top2_idx <- max.col(egc_t, ties.method = "first")
            mc_types$top2_gene <- gene_names_egc[top2_idx]
            mc_types$top2_lfp <- log2(egc_t[cbind(seq_len(nrow(egc_t)), top2_idx)] + 1e-5)
        } else {
            # Fallback: add NA columns
            mc_types$top1_gene <- NA_character_
            mc_types$top2_gene <- NA_character_
            mc_types$top1_lfp <- NA_real_
            mc_types$top2_lfp <- NA_real_
        }
    }

    mc_types
}

convert_daf_cell_type_colors <- function(daf_obj) {
    type_names <- dafr::axis_entries(daf_obj, "type")
    colors <- daf_vec(daf_obj, "type", "color")

    tibble(
        cell_type = type_names,
        color = colors,
        order = seq_along(type_names)
    )
}

# ==============================================================================
# Optional Component Conversions
# ==============================================================================

convert_daf_lateral_genes <- function(daf_obj) {
    convert_daf_flagged_genes(daf_obj, "is_lateral")
}

convert_daf_noisy_genes <- function(daf_obj) {
    convert_daf_flagged_genes(daf_obj, "is_noisy")
}

convert_daf_marker_genes <- function(daf_obj) {
    marker_flags <- daf_vec(daf_obj, "gene", "is_marker", required = FALSE)
    if (is.null(marker_flags) || sum(marker_flags, na.rm = TRUE) == 0) {
        return(tibble(
            gene = character(), metacell = character(),
            rank = integer(), fp = numeric()
        ))
    }

    gene_names <- dafr::axis_entries(daf_obj, "gene")
    marker_genes <- gene_names[marker_flags]

    tibble(
        gene = marker_genes,
        metacell = "all",
        rank = seq_along(marker_genes),
        fp = 1.0
    )
}

convert_daf_gene_modules <- function(daf_obj) {
    modules <- daf_vec(daf_obj, "gene", "module", required = FALSE)
    if (is.null(modules)) {
        return(NULL)
    }

    gene_names <- dafr::axis_entries(daf_obj, "gene")
    result <- tibble::tibble(gene = gene_names, module = modules)
    result <- result[!is.na(result$module), ]
    result
}

convert_daf_inner_fold_mat <- function(daf_obj) {
    daf_mat(daf_obj, "gene", "metacell", "inner_fold", required = FALSE)
}

convert_daf_inner_stdev_mat <- function(daf_obj) {
    daf_mat(daf_obj, "gene", "metacell", "inner_stdev_log", required = FALSE)
}

# ==============================================================================
# Metadata and QC Conversions
# ==============================================================================

convert_daf_metadata <- function(daf_obj) {
    # Enumerate all metacell vectors except core fields
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    # Core fields that are handled elsewhere
    core_fields <- c("type", "x", "y", "u", "v", "total_UMIs", "n_cell")

    # Get list of available metacell properties
    props <- tryCatch(
        {
            daf_obj["/ metacell ?"]
        },
        error = function(e) {
            return(character(0))
        }
    )

    # Filter out core fields and QC fields (handled separately)
    metadata_fields <- setdiff(props, core_fields)

    if (length(metadata_fields) == 0) {
        return(NULL)
    }

    # Batch-load all metadata vectors in a single JuliaCall via get_frame
    result <- tryCatch(
        {
            frame <- dafr::get_frame(daf_obj, "metacell", as.list(metadata_fields))
            # get_frame returns a data.frame; prepend metacell names column
            frame <- tibble::as_tibble(frame)
            frame <- tibble::tibble(metacell = metacell_names, !!!frame)
            frame
        },
        error = function(e) {
            # Fallback: load vectors one-by-one if get_frame fails
            res <- tibble(metacell = metacell_names)
            for (field in metadata_fields) {
                vec <- daf_vec(daf_obj, "metacell", field, required = FALSE)
                if (!is.null(vec)) {
                    res[[field]] <- vec
                }
            }
            res
        }
    )

    if (ncol(result) == 1) {
        return(NULL)
    }

    return(result)
}

convert_daf_cell_metadata <- function(daf_obj) {
    # Check if the DAF has a "cell" axis (cell-level data)
    if (!dafr::has_axis(daf_obj, "cell")) {
        return(NULL)
    }

    cell_names <- dafr::axis_entries(daf_obj, "cell")

    # Read the required cell.metacell vector
    mc_vec <- daf_vec(daf_obj, "cell", "metacell", required = FALSE)
    if (is.null(mc_vec)) {
        return(NULL)
    }

    # Read the required cell.samp_id vector
    samp_id_vec <- daf_vec(daf_obj, "cell", "samp_id", required = FALSE)
    if (is.null(samp_id_vec)) {
        return(NULL)
    }

    # Build the base tibble with cell identifier, metacell, and samp_id
    result <- tibble(
        cell_id = cell_names,
        metacell = as.character(mc_vec),
        samp_id = as.character(samp_id_vec)
    )

    # Discover and add any additional cell-level vectors
    cell_props <- tryCatch(
        {
            daf_obj["/ cell ?"]
        },
        error = function(e) {
            return(character(0))
        }
    )

    # Fields already handled above
    handled_fields <- c("metacell", "samp_id")

    extra_fields <- setdiff(cell_props, handled_fields)

    for (field in extra_fields) {
        vec <- daf_vec(daf_obj, "cell", field, required = FALSE)
        if (!is.null(vec)) {
            result[[field]] <- vec
        }
    }

    return(result)
}

convert_daf_mc_qc_metadata <- function(daf_obj) {
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")
    result <- tibble(metacell = metacell_names)

    # Handle fields with fallback names
    result <- add_optional_vec_with_fallback(result, daf_obj, "metacell", "umis", "total_UMIs")
    result <- add_optional_vec_with_fallback(result, daf_obj, "metacell", "cells", "n_cell")

    # Add standard QC vectors
    result <- add_optional_vecs(result, daf_obj, "metacell", c(
        "max_inner_fold", "max_inner_fold_no_lateral",
        "max_inner_stdev_log", "zero_fold", "rare_gene_module", "is_rare"
    ))

    # Compute max_inner_fold per metacell using DAF aggregation if not present
    if (is.null(result$max_inner_fold) && dafr::has_matrix(daf_obj, "gene", "metacell", "inner_fold")) {
        max_if <- tryCatch(
            daf_obj["/ gene / metacell : inner_fold %> Max"],
            error = function(e) NULL
        )
        if (!is.null(max_if)) {
            result$max_inner_fold <- as.numeric(max_if)
        }
    }

    if (ncol(result) == 1) {
        return(NULL)
    }
    result
}

convert_daf_gene_qc <- function(daf_obj) {
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    result <- tibble(gene = gene_names)

    # Handle max_expr with fallback
    result <- add_optional_vec_with_fallback(result, daf_obj, "gene", "max_expr", "mcview_cache_gene_max_umis")

    # Add standard gene QC vectors
    result <- add_optional_vecs(result, daf_obj, "gene", c(
        "type", "is_marker", "significant_inner_folds_count",
        "is_lateral", "is_noisy", "module", "correction_factor"
    ))

    # Try precomputed vector first, then DAF aggregation, then full matrix fallback
    max_if <- daf_vec(daf_obj, "gene", "max_inner_fold", required = FALSE)
    if (!is.null(max_if)) {
        result$max_inner_fold <- max_if
    } else if (dafr::has_matrix(daf_obj, "gene", "metacell", "inner_fold")) {
        max_if <- tryCatch(
            daf_obj["/ metacell / gene : inner_fold %> Max"],
            error = function(e) NULL
        )
        if (!is.null(max_if)) {
            result$max_inner_fold <- as.numeric(max_if)
        }
    }

    # Include fitted gene metadata if present
    gene_props <- tryCatch(daf_obj["/ gene ?"], error = function(e) character(0))
    fitted_fields <- grep("^fitted_gene_of", gene_props, value = TRUE)
    result <- add_optional_vecs(result, daf_obj, "gene", fitted_fields)

    if (ncol(result) == 1) {
        return(NULL)
    }
    result
}

# ==============================================================================
# Additional Table Conversions
# ==============================================================================

convert_daf_gg_mc_top_cor <- function(daf_obj) {
    if (!dafr::has_axis(daf_obj, "gg_mc_top_cor")) {
        return(NULL)
    }

    gene1 <- daf_vec(daf_obj, "gg_mc_top_cor", "gene1", required = FALSE)
    gene2 <- daf_vec(daf_obj, "gg_mc_top_cor", "gene2", required = FALSE)
    cor <- daf_vec(daf_obj, "gg_mc_top_cor", "cor", required = FALSE)
    type <- daf_vec(daf_obj, "gg_mc_top_cor", "type", required = FALSE)

    if (is.null(gene1) || is.null(gene2) || is.null(cor) || is.null(type)) {
        return(NULL)
    }

    tibble(
        gene1 = as.character(gene1),
        gene2 = as.character(gene2),
        cor = as.numeric(cor),
        type = as.character(type)
    )
}

convert_daf_gene_zero_fold <- function(daf_obj) {
    if (!dafr::has_axis(daf_obj, "gene_zero_fold")) {
        return(NULL)
    }

    required <- c("gene", "metacell", "zero_fold", "avg", "obs", "exp", "type")
    vectors <- purrr::map(required, ~ daf_vec(daf_obj, "gene_zero_fold", .x, required = FALSE))
    if (any(purrr::map_lgl(vectors, is.null))) {
        return(NULL)
    }

    tibble(
        gene = as.character(vectors[[1]]),
        metacell = as.character(vectors[[2]]),
        zero_fold = as.numeric(vectors[[3]]),
        avg = as.numeric(vectors[[4]]),
        obs = as.numeric(vectors[[5]]),
        exp = as.numeric(vectors[[6]]),
        type = as.character(vectors[[7]])
    )
}

convert_daf_metacell_graphs <- function(daf_obj) {
    if (!dafr::has_axis(daf_obj, "metacell_graph")) {
        return(NULL)
    }

    graph_name <- daf_vec(daf_obj, "metacell_graph", "graph_name", required = FALSE)
    from <- daf_vec(daf_obj, "metacell_graph", "from", required = FALSE)
    to <- daf_vec(daf_obj, "metacell_graph", "to", required = FALSE)
    weight <- daf_vec(daf_obj, "metacell_graph", "weight", required = FALSE)

    if (is.null(from) || is.null(to) || is.null(weight)) {
        return(NULL)
    }

    if (is.null(graph_name)) {
        graph_name <- rep("metacell", length(from))
    }

    graph_df <- tibble(
        graph_name = as.character(graph_name),
        from = as.character(from),
        to = as.character(to),
        weight = as.numeric(weight)
    )

    split(graph_df %>% select(from, to, weight), graph_df$graph_name)
}

convert_daf_qc_stats <- function(daf_obj) {
    fields <- c(
        "n_outliers",
        "n_cells",
        "n_umis",
        "median_umis_per_metacell",
        "median_cells_per_metacell"
    )

    stats <- list()
    for (field in fields) {
        scalar_name <- glue("qc_stats_{field}")
        if (dafr::has_scalar(daf_obj, scalar_name)) {
            stats[[field]] <- dafr::get_scalar(daf_obj, scalar_name)
        }
    }

    if (length(stats) == 0) {
        return(NULL)
    }

    stats
}

# ==============================================================================
# Projection/Atlas Conversions
# ==============================================================================

#' Convert a fraction matrix (metacell x gene) to UMI-scale (gene x metacell)
#'
#' Shared helper for the common pattern: load fraction matrix, transpose,
#' multiply by total_UMIs, and label axes. Used by convert_daf_mc_mat_corrected
#' and convert_daf_projected_mat.
#'
#' @param daf_obj DAF object
#' @param property_name Name of the metacell-gene fraction property
#'
#' @return Gene x metacell matrix scaled to UMIs, or NULL if property missing
#' @noRd
convert_daf_fraction_to_umi <- function(daf_obj, property_name) {
    frac_mat <- daf_mat(daf_obj, "metacell", "gene", property_name, required = FALSE)
    if (is.null(frac_mat)) {
        return(NULL)
    }

    mc_sum <- daf_vec(daf_obj, "metacell", "total_UMIs")
    # frac_mat is metacell (rows) x gene (cols).
    # Multiply each row by the corresponding metacell's total_UMIs, then transpose
    # to get gene x metacell. Using frac_mat * mc_sum is correct because R recycles
    # the mc_sum vector (length = nrow) down each column, matching row indices.
    umi_mat <- Matrix::t(frac_mat * mc_sum)
    rownames(umi_mat) <- dafr::axis_entries(daf_obj, "gene")
    colnames(umi_mat) <- dafr::axis_entries(daf_obj, "metacell")

    return(umi_mat)
}

convert_daf_projected_fold <- function(daf_obj) {
    daf_mat(daf_obj, "gene", "metacell", "projected_fold", required = FALSE)
}

convert_daf_mc_mat_corrected <- function(daf_obj) {
    convert_daf_fraction_to_umi(daf_obj, "corrected_fraction")
}

convert_daf_query_metadata <- function(daf_obj) {
    # Query metadata for atlas projection
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    result <- tibble(metacell = metacell_names)

    # Projected type from atlas
    projected_type <- daf_vec(daf_obj, "metacell", "projected_type", required = FALSE)
    if (!is.null(projected_type)) {
        result$projected_type <- projected_type
    }

    # Similarity flag
    similar <- daf_vec(daf_obj, "metacell", "similar", required = FALSE)
    if (!is.null(similar)) {
        result$similar <- similar
    }

    # Projection correlation
    proj_cor <- daf_vec(daf_obj, "metacell", "projected_correlation", required = FALSE)
    if (!is.null(proj_cor)) {
        result$projected_correlation <- proj_cor
    }

    # Most similar atlas metacell
    atlas_metacell <- daf_vec(daf_obj, "metacell", "atlas_metacell", required = FALSE)
    if (!is.null(atlas_metacell)) {
        result$atlas_metacell <- atlas_metacell
    }

    if (ncol(result) == 1) {
        return(NULL)
    }

    return(result)
}

convert_daf_projected_mat <- function(daf_obj) {
    convert_daf_fraction_to_umi(daf_obj, "projected_fraction")
}

convert_daf_proj_weights <- function(daf_obj) {
    # Projection weights stored as JSON
    proj_weights_json <- daf_scalar(daf_obj, "projection_weights_json", default = NULL)
    if (is.null(proj_weights_json)) {
        return(NULL)
    }

    tryCatch(
        {
            proj_weights <- jsonlite::fromJSON(proj_weights_json)
            as_tibble(proj_weights)
        },
        error = function(e) {
            cli_warn("Failed to parse projection weights JSON: {e$message}")
            NULL
        }
    )
}

convert_daf_query_cell_type_fracs <- function(daf_obj) {
    # Query cell type fractions stored as JSON
    fracs_json <- daf_scalar(daf_obj, "query_cell_type_fracs_json", default = NULL)
    if (is.null(fracs_json)) {
        return(NULL)
    }

    tryCatch(
        {
            fracs <- jsonlite::fromJSON(fracs_json)
            as_tibble(fracs)
        },
        error = function(e) {
            cli_warn("Failed to parse query cell type fractions JSON: {e$message}")
            NULL
        }
    )
}
