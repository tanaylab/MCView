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
            query <- glue::glue("@ gene [ {flag_name} ] : name")
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

#' Try multiple DAF names and return the first found
#'
#' Generic helper that tries a sequence of names for a DAF vector and returns
#' the value for the first name that exists. Useful for backward-compatible
#' reading where the canonical name changed across Metacells.jl versions.
#'
#' @param daf_obj DAF object
#' @param axis Axis name (e.g. "metacell", "gene", "type")
#' @param names_vec Character vector of names to try, in priority order
#' @param required If TRUE and no name is found, throws an error
#' @return The vector for the first name found, or NULL if none found
#' @export
try_daf_names <- function(daf_obj, axis, names_vec, required = FALSE) {
    for (nm in names_vec) {
        vec <- daf_vec(daf_obj, axis, nm, required = FALSE)
        if (!is.null(vec)) {
            return(vec)
        }
    }
    if (required) {
        cli_abort("Missing required vector: {axis}.{paste(names_vec, collapse = '|')}")
    }
    NULL
}

#' Try multiple DAF scalar names and return the first found
#'
#' Like try_daf_names but for scalars.
#'
#' @param daf_obj DAF object
#' @param names_vec Character vector of scalar names to try, in priority order
#' @param default Default value if no name is found
#' @return The scalar value for the first name found, or default
#' @export
try_daf_scalar_names <- function(daf_obj, names_vec, default = NULL) {
    for (nm in names_vec) {
        if (dafr::has_scalar(daf_obj, nm)) {
            return(dafr::get_scalar(daf_obj, nm))
        }
    }
    default
}

#' Try multiple DAF matrix names and return the first found
#'
#' Like try_daf_names but for matrices.
#'
#' @param daf_obj DAF object
#' @param axis1 First axis name
#' @param axis2 Second axis name
#' @param names_vec Character vector of matrix names to try, in priority order
#' @param required If TRUE and no name is found, throws an error
#' @return The matrix for the first name found, or NULL if none found
#' @export
try_daf_matrix_names <- function(daf_obj, axis1, axis2, names_vec, required = FALSE) {
    for (nm in names_vec) {
        if (dafr::has_matrix(daf_obj, axis1, axis2, nm)) {
            return(dafr::get_matrix(daf_obj, axis1, axis2, nm))
        }
    }
    if (required) {
        cli_abort("Missing required matrix: {axis1},{axis2}.{paste(names_vec, collapse = '|')}")
    }
    NULL
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
    mc_mat <- daf_query_mc_mat(daf_obj, genes = genes, metacells = metacells)
    mc_sum <- daf_query_mc_sum(daf_obj, metacells = metacells)

    # Filter mc_sum to match metacells in matrix
    if (!is.null(metacells)) {
        mc_sum <- mc_sum[intersect(metacells, names(mc_sum))]
    }

    # Divide each column by the corresponding mc_sum value.
    sweep(mc_mat, 2, mc_sum, "/")
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
    query <- paste0("@ metacell @ gene :: UMIs >> ", agg_op)
    result <- daf_obj[query]
    # dafr may return a named vector; skip redundant names assignment.
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
            # Query each gene individually - avoids loading full gene x metacell matrix
            vecs <- lapply(valid_genes, function(gene) {
                query <- glue::glue("@ metacell @ gene = {escape_daf_value(gene)} :: UMIs")
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
                query <- glue::glue("@ gene @ metacell = {escape_daf_value(mc)} :: UMIs")
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
    # Query: @ metacell @ gene :: UMIs -/ type >- Sum
    # GroupRowsBy type (metacell axis) and reduce columns with Sum
    ct_mat <- daf_obj["@ metacell @ gene :: UMIs -/ type >- Sum"]

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
    # Query: @ metacell : total_UMIs / type >> Sum
    result <- daf_obj["@ metacell : total_UMIs / type >> Sum"]

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
        "marker_genes_projected" = convert_daf_marker_genes_from_mat(daf_obj, "projected_fold"),
        "mc_qc_metadata" = convert_daf_mc_qc_metadata(daf_obj),
        "gene_qc" = convert_daf_gene_qc(daf_obj),
        "gg_mc_top_cor" = convert_daf_gg_mc_top_cor(daf_obj),
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
        {
            cli::cli_warn("Unknown MCView conversion type: {var_name}")
            NULL
        }
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
    # dafr::get_vector() returns a named vector; skip redundant names assignment.
    if (is.null(names(mc_sum))) {
        names(mc_sum) <- dafr::axis_entries(daf_obj, "metacell")
    }
    mc_sum
}

convert_daf_mc_egc <- function(daf_obj) {
    # Try pre-computed fraction matrices in priority order:
    # 1. linear_fraction (Metacells.jl canonical: UMIs / total_UMIs)
    # 2. geomean_fraction (MCView/older pipelines: geometric mean of per-cell fractions)
    # Note: these are DIFFERENT computations, but both serve as EGC representation.
    egc_mat <- try_daf_matrix_names(
        daf_obj, "gene", "metacell",
        c("linear_fraction", "geomean_fraction")
    )
    if (!is.null(egc_mat)) {
        which_name <- if (dafr::has_matrix(daf_obj, "gene", "metacell", "linear_fraction")) {
            "linear_fraction"
        } else {
            "geomean_fraction"
        }
        cli::cli_inform("Using pre-computed {.field {which_name}} matrix from DAF")
        return(egc_mat)
    }

    # Fall back: compute linear fraction from UMIs / total_UMIs
    cli::cli_inform("Computing EGC as UMIs / total_UMIs (no fraction matrix in DAF)")
    mc_mat <- convert_daf_mc_mat(daf_obj)
    mc_sum <- convert_daf_mc_sum(daf_obj)

    # Compute EGC (normalized expression) - gene x metacell
    # Divide each column by its total to get fractions summing to 1.
    sweep(mc_mat, 2, mc_sum, "/")
}

convert_daf_mc2d <- function(daf_obj) {
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    # Try coordinate names in priority order: x/y, then umap_x/umap_y, then u/v
    x_coords <- try_daf_names(daf_obj, "metacell", c("x", "umap_x", "u"))
    y_coords <- try_daf_names(daf_obj, "metacell", c("y", "umap_y", "v"))

    if (is.null(x_coords) || is.null(y_coords)) {
        cli_abort("No 2D coordinates found in DAF (tried x/y, umap_x/umap_y, u/v)")
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

    # Add optional cell count: try canonical "n_cells" first, fall back to legacy "n_cell"
    n_cell <- try_daf_names(daf_obj, "metacell", c("n_cells", "n_cell"))
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
            # Top-2 genes per metacell (mc_egc is gene x metacell).
            gene_names_egc <- rownames(mc_egc)
            tops <- top2_per_col(mc_egc)
            mc_types$top1_gene <- gene_names_egc[tops$top1_idx]
            mc_types$top1_lfp <- log2(tops$top1_val + 1e-5)
            mc_types$top2_gene <- gene_names_egc[tops$top2_idx]
            mc_types$top2_lfp <- log2(tops$top2_val + 1e-5)
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

    # Use ("type", "order") from DAF if available; otherwise default to axis order
    type_order <- daf_vec(daf_obj, "type", "order", required = FALSE)
    if (!is.null(type_order)) {
        ord <- as.integer(type_order)
    } else {
        ord <- seq_along(type_names)
    }

    result <- tibble(
        cell_type = type_names,
        color = colors,
        order = ord
    )

    # Sort by order so that downstream code respects the DAF ordering
    result <- result[order(result$order), ]

    result
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

    # Use marker_rank from DAF if available (Metacells.jl provides this)
    marker_rank <- daf_vec(daf_obj, "gene", "marker_rank", required = FALSE)
    if (!is.null(marker_rank)) {
        # Extract ranks for marker genes only, then sort by rank
        marker_ranks <- marker_rank[marker_flags]
        names(marker_ranks) <- marker_genes
        # Remove any NAs/zeros (unranked genes)
        valid <- !is.na(marker_ranks) & marker_ranks > 0
        if (sum(valid) > 0) {
            marker_ranks <- marker_ranks[valid]
            ord <- order(marker_ranks)
            marker_genes_sorted <- names(marker_ranks)[ord]
            return(tibble(
                gene = marker_genes_sorted,
                metacell = "all",
                rank = seq_along(marker_genes_sorted),
                fp = 1.0
            ))
        }
    }

    # Fallback: no marker_rank, use original order
    tibble(
        gene = marker_genes,
        metacell = "all",
        rank = seq_along(marker_genes),
        fp = 1.0
    )
}

convert_daf_gene_modules <- function(daf_obj) {
    # Try canonical "primary_module" first, fall back to legacy "module"
    modules <- try_daf_names(daf_obj, "gene", c("primary_module", "module"))
    if (is.null(modules)) {
        return(NULL)
    }

    gene_names <- dafr::axis_entries(daf_obj, "gene")
    result <- tibble::tibble(gene = gene_names, module = modules)
    result <- result[!is.na(result$module), ]
    result
}

#' Compute marker genes from a fold-change matrix
#'
#' For Projected-fold heatmaps the marker gene list must be derived from the
#' matrix itself (top variable genes per metacell). This function selects the
#' top-2 genes per metacell by absolute value, using the same
#' select_top_fold_genes_per_metacell() logic used elsewhere in the codebase.
#'
#' @param daf_obj DAF object
#' @param mode Currently only "projected_fold" is supported
#' @return Tibble with columns metacell, gene, rank, fp; or NULL if matrix absent
#' @noRd
convert_daf_marker_genes_from_mat <- function(daf_obj, mode) {
    mat <- switch(mode,
        "projected_fold" = daf_mat(daf_obj, "gene", "metacell", "projected_fold", required = FALSE),
        NULL
    )
    if (is.null(mat)) {
        return(NULL)
    }
    # Keep only genes with non-zero rows
    row_nz <- if (is(mat, "sparseMatrix")) Matrix::rowSums(mat) != 0 else rowSums(mat) != 0
    mat <- mat[row_nz, , drop = FALSE]
    if (nrow(mat) == 0) {
        return(NULL)
    }
    mat_dense <- as.matrix(mat)
    select_top_fold_genes_per_metacell(
        mat_dense,
        genes_per_metacell = 2,
        minimal_relative_log_fraction = 0.5,
        fold_change_reg = 0
    )
}

# ==============================================================================
# Metadata and QC Conversions
# ==============================================================================

convert_daf_metadata <- function(daf_obj) {
    # Enumerate all metacell vectors except core fields
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    # Core fields that are handled elsewhere
    core_fields <- c("type", "x", "y", "u", "v", "umap_x", "umap_y", "total_UMIs", "n_cells", "n_cell")

    # Get list of available metacell properties
    props <- tryCatch(
        {
            daf_obj["@ metacell : ?"]
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

    # Load each metadata vector individually via daf_vec rather than build a
    # full DataFrame up front via get_frame().
    result <- tibble(metacell = metacell_names)
    for (field in metadata_fields) {
        vec <- daf_vec(daf_obj, "metacell", field, required = FALSE)
        if (!is.null(vec)) {
            result[[field]] <- vec
        }
    }

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

    # Determine the sample ID source vector:
    # 1. Check mcview_sample_property scalar (names which cell vector to use)
    # 2. Fall back to literal "samp_id" vector
    # 3. If neither exists, build tibble without samp_id (Samples tab disabled)
    sample_property <- NULL
    if (dafr::has_scalar(daf_obj, "mcview_sample_property")) {
        sample_property <- dafr::get_scalar(daf_obj, "mcview_sample_property")
    }

    samp_id_vec <- NULL
    if (!is.null(sample_property)) {
        samp_id_vec <- daf_vec(daf_obj, "cell", sample_property, required = FALSE)
        if (is.null(samp_id_vec)) {
            cli_warn("mcview_sample_property = '{sample_property}' but cell vector '{sample_property}' not found; falling back to 'samp_id'")
            samp_id_vec <- daf_vec(daf_obj, "cell", "samp_id", required = FALSE)
        }
    } else {
        samp_id_vec <- daf_vec(daf_obj, "cell", "samp_id", required = FALSE)
    }

    # Build the base tibble with cell identifier and metacell
    result <- tibble(
        cell_id = cell_names,
        metacell = as.character(mc_vec)
    )

    # Add samp_id column if we found a source vector
    if (!is.null(samp_id_vec)) {
        result$samp_id <- as.character(samp_id_vec)
    }

    # Discover and add any additional cell-level vectors
    cell_props <- tryCatch(
        {
            daf_obj["@ cell : ?"]
        },
        error = function(e) {
            return(character(0))
        }
    )

    # Fields already handled above (samp_id always listed to avoid duplication,
    # plus the original source property name if different)
    handled_fields <- c("metacell", "samp_id")
    if (!is.null(sample_property) && sample_property != "samp_id") {
        handled_fields <- c(handled_fields, sample_property)
    }

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
    # Try canonical "n_cells" first, then legacy "n_cell"
    n_cells_vec <- try_daf_names(daf_obj, "metacell", c("n_cells", "n_cell"))
    if (!is.null(n_cells_vec)) {
        result[["cells"]] <- n_cells_vec
    }

    # Add standard QC vectors
    result <- add_optional_vecs(result, daf_obj, "metacell", c(
        "rare_gene_module", "is_rare"
    ))

    if (ncol(result) == 1) {
        return(NULL)
    }
    result
}

convert_daf_gene_qc <- function(daf_obj) {
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    result <- tibble(gene = gene_names)

    # Handle max_expr: try precomputed vectors, then compute from UMIs matrix
    result <- add_optional_vec_with_fallback(result, daf_obj, "gene", "max_expr", "mcview_cache_gene_max_umis")
    if (!"max_expr" %in% colnames(result) && dafr::has_matrix(daf_obj, "metacell", "gene", "UMIs")) {
        max_umis <- tryCatch(
            {
                # Compute per-gene max from UMIs matrix (gene x metacell after transpose)
                umat <- dafr::get_matrix(daf_obj, "gene", "metacell", "UMIs")
                matrixStats::rowMaxs(as.matrix(umat))
            },
            error = function(e) NULL
        )
        if (!is.null(max_umis)) {
            result$max_expr <- max_umis
        }
    }

    # Add standard gene QC vectors
    result <- add_optional_vecs(result, daf_obj, "gene", c(
        "type", "is_marker",
        "is_lateral", "is_noisy", "correction_factor"
    ))

    # Add module: try canonical "primary_module" first, fall back to "module"
    mod_vec <- try_daf_names(daf_obj, "gene", c("primary_module", "module"))
    if (!is.null(mod_vec)) {
        result[["module"]] <- mod_vec
    }

    # Synthesize 'type' column from is_lateral/is_noisy if the pre-computed
    # type vector is absent.  Downstream QC scatter plots colour genes by type,
    # so the column must exist even if only as "other".
    if (!"type" %in% colnames(result)) {
        is_lat <- if ("is_lateral" %in% colnames(result)) result$is_lateral else rep(FALSE, nrow(result))
        is_noi <- if ("is_noisy" %in% colnames(result)) result$is_noisy else rep(FALSE, nrow(result))
        result$type <- dplyr::case_when(
            is_lat & is_noi ~ "lateral, noisy",
            is_lat ~ "lateral",
            is_noi ~ "noisy",
            TRUE ~ "other"
        )
    }

    # Include fitted gene metadata if present
    gene_props <- tryCatch(daf_obj["@ gene : ?"], error = function(e) character(0))
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
    # Map from result field name to (canonical, legacy) scalar names.
    # Try canonical names first (no qc_stats_ prefix), fall back to legacy.
    # Special case: "n_cells" uses "n_cells_total" canonically to avoid
    # collision with the ("metacell", "n_cells") vector name.
    field_aliases <- list(
        n_outliers = c("n_outliers", "qc_stats_n_outliers"),
        n_cells    = c("n_cells_total", "qc_stats_n_cells"),
        n_umis     = c("n_umis", "qc_stats_n_umis"),
        median_umis_per_metacell  = c("median_umis_per_metacell", "qc_stats_median_umis_per_metacell"),
        median_cells_per_metacell = c("median_cells_per_metacell", "qc_stats_median_cells_per_metacell")
    )

    stats <- list()
    for (field in names(field_aliases)) {
        val <- try_daf_scalar_names(daf_obj, field_aliases[[field]])
        if (!is.null(val)) {
            stats[[field]] <- val
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
    if (is.null(rownames(umi_mat))) {
        rownames(umi_mat) <- dafr::axis_entries(daf_obj, "gene")
    }
    if (is.null(colnames(umi_mat))) {
        colnames(umi_mat) <- dafr::axis_entries(daf_obj, "metacell")
    }

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

    # Similarity flag: try "similar" first, then "is_similar" as alias
    similar <- try_daf_names(daf_obj, "metacell", c("similar", "is_similar"))
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

#' Read per-type marker genes from DAF
#'
#' Reads type markers from DAF. Tries the new sparse matrix format first:
#' ("type","gene","mcview_marker_rank") and ("type","gene","mcview_marker_fold_change"),
#' then falls back to the legacy mcview_type_markers axis+vectors format.
#'
#' @param daf_obj DAF object (typically a complete_daf chain)
#' @return Tibble with columns: cell_type, gene, rank, fold_change.
#'   Returns NULL if neither format exists.
#' @export
convert_daf_type_markers <- function(daf_obj) {
    # Try new sparse matrix format first
    if (dafr::has_matrix(daf_obj, "type", "gene", "mcview_marker_rank")) {
        result <- tryCatch(
            {
                rank_mat <- dafr::get_matrix(daf_obj, "type", "gene", "mcview_marker_rank")
                fc_mat <- NULL
                if (dafr::has_matrix(daf_obj, "type", "gene", "mcview_marker_fold_change")) {
                    fc_mat <- dafr::get_matrix(daf_obj, "type", "gene", "mcview_marker_fold_change")
                }

                # Convert sparse matrices to tibble
                type_names <- rownames(rank_mat)
                gene_names <- colnames(rank_mat)
                if (is.null(type_names)) type_names <- dafr::axis_entries(daf_obj, "type")
                if (is.null(gene_names)) gene_names <- dafr::axis_entries(daf_obj, "gene")

                # Extract non-zero entries from rank matrix
                if (is(rank_mat, "sparseMatrix")) {
                    sm <- Matrix::summary(rank_mat)
                    if (nrow(sm) == 0) return(NULL)
                    tibble(
                        cell_type = type_names[sm$i],
                        gene = gene_names[sm$j],
                        rank = as.integer(sm$x),
                        fold_change = if (!is.null(fc_mat)) as.numeric(fc_mat[cbind(sm$i, sm$j)]) else NA_real_
                    ) %>% arrange(cell_type, rank)
                } else {
                    # Dense matrix: extract non-zero entries
                    nz <- which(rank_mat != 0, arr.ind = TRUE)
                    if (nrow(nz) == 0) return(NULL)
                    tibble(
                        cell_type = type_names[nz[, 1]],
                        gene = gene_names[nz[, 2]],
                        rank = as.integer(rank_mat[nz]),
                        fold_change = if (!is.null(fc_mat)) as.numeric(fc_mat[nz]) else NA_real_
                    ) %>% arrange(cell_type, rank)
                }
            },
            error = function(e) {
                cli_warn("Failed to read type markers from sparse matrices: {e$message}")
                NULL
            }
        )
        if (!is.null(result)) return(result)
    }

    # Fall back to legacy axis+vectors format
    if (!dafr::has_axis(daf_obj, "mcview_type_markers")) {
        return(NULL)
    }

    tryCatch(
        {
            tibble(
                cell_type = dafr::get_vector(daf_obj, "mcview_type_markers", "cell_type"),
                gene = dafr::get_vector(daf_obj, "mcview_type_markers", "gene"),
                rank = dafr::get_vector(daf_obj, "mcview_type_markers", "rank"),
                fold_change = dafr::get_vector(daf_obj, "mcview_type_markers", "fold_change")
            )
        },
        error = function(e) {
            cli_warn("Failed to read type markers from DAF: {e$message}")
            NULL
        }
    )
}
