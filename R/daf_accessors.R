# daf_accessors.R - Safe DAF accessors and shared optional/try helpers
#
# Split from R/daf_data.R (2026-05-01). Companions:
#   - R/daf_queries.R    - daf_query_* functions
#   - R/daf_conversion.R - convert_daf_* functions


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
#' Returns the gene x metacell expression matrix, preferring whatever
#' pre-computed fraction matrix the DAF ships - in priority order
#' `linear_fraction`, `geomean_fraction`. When neither exists, asks dafr
#' to compute UMIs % Fraction inline (the `% Fraction` eltwise divides
#' each column by its own column sum, matching the historical
#' sweep(UMIs, 2, total_UMIs, "/") whenever total_UMIs == colSums(UMIs) -
#' the standard import-pipeline invariant).
#'
#' Preferring the stored fraction matrix also keeps this cold path in
#' agreement with `convert_daf_mc_egc` (the warm `mc_egc_full` builder)
#' - see project_linear_geomean_egc_mismatch.md for the latent
#' discrepancy this avoids.
#'
#' @param daf_obj DAF object
#' @param genes Optional vector of gene names to filter
#' @param metacells Optional vector of metacell names to filter
#' @return EGC matrix (genes x metacells) with columns summing to 1
#' @export
compute_egc_from_daf <- function(daf_obj, genes = NULL, metacells = NULL) {
    fraction_name <- if (dafr::has_matrix(daf_obj, "gene", "metacell", "linear_fraction")) {
        "linear_fraction"
    } else if (dafr::has_matrix(daf_obj, "gene", "metacell", "geomean_fraction")) {
        "geomean_fraction"
    } else {
        NULL
    }

    if (!is.null(fraction_name)) {
        mat <- dafr::get_matrix(daf_obj, "gene", "metacell", fraction_name)
        if (!is.null(genes)) {
            mat <- mat[intersect(genes, rownames(mat)), , drop = FALSE]
        }
        if (!is.null(metacells)) {
            mat <- mat[, intersect(metacells, colnames(mat)), drop = FALSE]
        }
        return(mat)
    }

    # No fraction matrix in the DAF - let dafr compute it. Full-matrix
    # case is one query; filtered case routes through daf_query_mc_mat's
    # per-gene / per-metacell fast paths to avoid materialising the
    # whole UMI matrix, then divides by total_UMIs in R.
    if (is.null(genes) && is.null(metacells)) {
        return(daf_obj["@ gene @ metacell :: UMIs % Fraction"])
    }

    mc_mat <- daf_query_mc_mat(daf_obj, genes = genes, metacells = metacells)
    mc_sum <- daf_query_mc_sum(daf_obj, metacells = metacells)
    if (!is.null(metacells)) {
        mc_sum <- mc_sum[intersect(metacells, names(mc_sum))]
    }
    sweep(mc_mat, 2, mc_sum, "/")
}

