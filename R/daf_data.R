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
#' @param cache Whether to cache the query result (currently unused, kept for API compatibility)
#'
#' @return Matrix slice with genes as rows and metacells as columns
#' @export
daf_query_mc_mat <- function(daf_obj, genes = NULL, metacells = NULL, cache = FALSE) {
    if (!is.null(genes) && length(genes) == 1 && !is.null(genes[[1]]) && !is.na(genes[[1]])) {
        gene <- as.character(genes[[1]])
        query <- glue::glue("/ metacell / gene = {escape_daf_value(gene)} : UMIs")
        vec <- tryCatch(daf_obj[query], error = function(e) NULL)
        if (!is.null(vec) && length(vec) > 0) {
            if (is.null(names(vec))) {
                names(vec) <- dafr::axis_entries(daf_obj, "metacell")
            }
            if (!is.null(metacells)) {
                metacells <- metacells[metacells %in% names(vec)]
                vec <- vec[metacells]
            }
            mat <- Matrix::Matrix(as.numeric(vec), nrow = 1, sparse = TRUE)
            rownames(mat) <- gene
            colnames(mat) <- names(vec)
            return(mat)
        }
    }

    # Get the full UMIs matrix
    # DAF stores as metacell x gene, we need gene x metacell
    mc_mat <- dafr::get_matrix(daf_obj, "metacell", "gene", "UMIs")
    mc_mat <- Matrix::t(mc_mat)

    # Get axis entries for proper naming
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    rownames(mc_mat) <- gene_names
    colnames(mc_mat) <- metacell_names

    # Filter by genes if specified
    if (!is.null(genes)) {
        valid_genes <- intersect(genes, gene_names)
        if (length(valid_genes) > 0) {
            mc_mat <- mc_mat[valid_genes, , drop = FALSE]
        }
    }

    # Filter by metacells if specified
    if (!is.null(metacells)) {
        valid_metacells <- intersect(metacells, metacell_names)
        if (length(valid_metacells) > 0) {
            mc_mat <- mc_mat[, valid_metacells, drop = FALSE]
        }
    }

    return(mc_mat)
}

#' Retrieve total UMIs per metacell via DAF
#'
#' @param daf_obj DAF object
#' @param metacells Optional vector of metacell names to retrieve
#' @param cache Whether to cache the query result (currently unused, kept for API compatibility)
#'
#' @return Named vector of total UMIs per metacell
#' @export
daf_query_mc_sum <- function(daf_obj, metacells = NULL, cache = FALSE) {
    # Get full total_UMIs vector
    result <- dafr::get_vector(daf_obj, "metacell", "total_UMIs")

    # Get axis entries for naming
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")
    names(result) <- metacell_names

    # Filter by metacells if specified
    if (!is.null(metacells)) {
        valid_metacells <- intersect(metacells, metacell_names)
        if (length(valid_metacells) > 0) {
            result <- result[valid_metacells]
        }
    }

    return(result)
}

#' Retrieve cell type assignments via DAF
#'
#' @param daf_obj DAF object
#' @param metacells Optional vector of metacell names to retrieve
#' @param cache Whether to cache the query result (currently unused, kept for API compatibility)
#'
#' @return Named vector of cell types per metacell
#' @export
daf_query_cell_types <- function(daf_obj, metacells = NULL, cache = FALSE) {
    # Get full type vector
    result <- dafr::get_vector(daf_obj, "metacell", "type")

    # Get axis entries for naming
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")
    names(result) <- metacell_names

    # Filter by metacells if specified
    if (!is.null(metacells)) {
        valid_metacells <- intersect(metacells, metacell_names)
        if (length(valid_metacells) > 0) {
            result <- result[valid_metacells]
        }
    }

    return(result)
}

#' Retrieve 2D coordinates via DAF
#'
#' @param daf_obj DAF object
#' @param metacells Optional vector of metacell names to retrieve
#' @param cache Whether to cache the query result (currently unused, kept for API compatibility)
#'
#' @return Data frame with metacell, x, y columns
#' @export
daf_query_2d_coords <- function(daf_obj, metacells = NULL, cache = FALSE) {
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
    result <- daf_obj["/ metacell / gene : UMIs %> Max"]
    names(result) <- dafr::axis_entries(daf_obj, "gene")
    return(result)
}

#' Get mean UMIs per gene using DAF query
#'
#' @param daf_obj DAF object
#' @return Named vector of mean UMIs per gene
#' @export
daf_query_gene_mean_umis <- function(daf_obj) {
    result <- daf_obj["/ metacell / gene : UMIs %> Mean"]
    names(result) <- dafr::axis_entries(daf_obj, "gene")
    return(result)
}

#' Get sum UMIs per gene using DAF query
#'
#' @param daf_obj DAF object
#' @return Named vector of sum UMIs per gene
#' @export
daf_query_gene_sum_umis <- function(daf_obj) {
    result <- daf_obj["/ metacell / gene : UMIs %> Sum"]
    names(result) <- dafr::axis_entries(daf_obj, "gene")
    return(result)
}

#' Get marker genes using DAF query
#'
#' Returns genes where is_marker is TRUE
#'
#' @param daf_obj DAF object
#' @return Character vector of marker gene names
#' @export
daf_query_marker_genes <- function(daf_obj) {
    if (!dafr::has_vector(daf_obj, "gene", "is_marker")) {
        return(character(0))
    }
    # Use DAF query to filter by is_marker
    tryCatch(
        {
            result <- daf_obj["/ gene & is_marker : name"]
            return(result)
        },
        error = function(e) {
            # Fallback: get all and filter in R
            is_marker <- dafr::get_vector(daf_obj, "gene", "is_marker")
            gene_names <- dafr::axis_entries(daf_obj, "gene")
            return(gene_names[is_marker])
        }
    )
}

#' Get lateral genes using DAF query
#'
#' Returns genes where is_lateral is TRUE
#'
#' @param daf_obj DAF object
#' @return Character vector of lateral gene names
#' @export
daf_query_lateral_genes <- function(daf_obj) {
    if (!dafr::has_vector(daf_obj, "gene", "is_lateral")) {
        return(character(0))
    }
    tryCatch(
        {
            result <- daf_obj["/ gene & is_lateral : name"]
            return(result)
        },
        error = function(e) {
            is_lateral <- dafr::get_vector(daf_obj, "gene", "is_lateral")
            gene_names <- dafr::axis_entries(daf_obj, "gene")
            return(gene_names[is_lateral])
        }
    )
}

#' Get noisy genes using DAF query
#'
#' Returns genes where is_noisy is TRUE
#'
#' @param daf_obj DAF object
#' @return Character vector of noisy gene names
#' @export
daf_query_noisy_genes <- function(daf_obj) {
    if (!dafr::has_vector(daf_obj, "gene", "is_noisy")) {
        return(character(0))
    }
    tryCatch(
        {
            result <- daf_obj["/ gene & is_noisy : name"]
            return(result)
        },
        error = function(e) {
            is_noisy <- dafr::get_vector(daf_obj, "gene", "is_noisy")
            gene_names <- dafr::axis_entries(daf_obj, "gene")
            return(gene_names[is_noisy])
        }
    )
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

    # Transpose to get modules as rows, metacells as columns
    mod_mat <- t(mod_mat)

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
    # Divide each column by its total and multiply by median total
    median_total <- median(mc_sum)
    mc_egc <- sweep(mc_mat, 2, mc_sum, "/") * median_total

    as.matrix(mc_egc)
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
            # Compute top 2 genes per metacell
            top_genes <- lapply(metacell_names, function(mc) {
                if (mc %in% colnames(mc_egc)) {
                    expr <- mc_egc[, mc]
                    top2_idx <- order(expr, decreasing = TRUE)[1:2]
                    list(
                        top1_gene = names(expr)[top2_idx[1]],
                        top2_gene = names(expr)[top2_idx[2]],
                        top1_lfp = log2(expr[top2_idx[1]]),
                        top2_lfp = log2(expr[top2_idx[2]])
                    )
                } else {
                    list(
                        top1_gene = NA_character_, top2_gene = NA_character_,
                        top1_lfp = NA_real_, top2_lfp = NA_real_
                    )
                }
            })
            mc_types$top1_gene <- sapply(top_genes, `[[`, "top1_gene")
            mc_types$top2_gene <- sapply(top_genes, `[[`, "top2_gene")
            mc_types$top1_lfp <- sapply(top_genes, `[[`, "top1_lfp")
            mc_types$top2_lfp <- sapply(top_genes, `[[`, "top2_lfp")
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
    lateral_flags <- daf_vec(daf_obj, "gene", "is_lateral", required = FALSE)
    if (is.null(lateral_flags)) {
        return(character(0))
    }

    gene_names <- dafr::axis_entries(daf_obj, "gene")
    gene_names[lateral_flags]
}

convert_daf_noisy_genes <- function(daf_obj) {
    noisy_flags <- daf_vec(daf_obj, "gene", "is_noisy", required = FALSE)
    if (is.null(noisy_flags)) {
        return(character(0))
    }

    gene_names <- dafr::axis_entries(daf_obj, "gene")
    gene_names[noisy_flags]
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
    daf_vec(daf_obj, "gene", "module", required = FALSE)
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

    # Build metadata tibble
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

convert_daf_mc_qc_metadata <- function(daf_obj) {
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")

    result <- tibble(metacell = metacell_names)

    # Prefer existing QC vectors from DAF if present
    umis <- daf_vec(daf_obj, "metacell", "umis", required = FALSE)
    if (is.null(umis)) {
        umis <- daf_vec(daf_obj, "metacell", "total_UMIs", required = FALSE)
    }
    if (!is.null(umis)) {
        result$umis <- umis
    }

    cells <- daf_vec(daf_obj, "metacell", "cells", required = FALSE)
    if (is.null(cells)) {
        cells <- daf_vec(daf_obj, "metacell", "n_cell", required = FALSE)
    }
    if (!is.null(cells)) {
        result$cells <- cells
    }

    max_inner_fold <- daf_vec(daf_obj, "metacell", "max_inner_fold", required = FALSE)
    if (!is.null(max_inner_fold)) {
        result$max_inner_fold <- max_inner_fold
    }

    max_inner_fold_no_lateral <- daf_vec(daf_obj, "metacell", "max_inner_fold_no_lateral", required = FALSE)
    if (!is.null(max_inner_fold_no_lateral)) {
        result$max_inner_fold_no_lateral <- max_inner_fold_no_lateral
    }

    max_inner_stdev_log <- daf_vec(daf_obj, "metacell", "max_inner_stdev_log", required = FALSE)
    if (!is.null(max_inner_stdev_log)) {
        result$max_inner_stdev_log <- max_inner_stdev_log
    }

    zero_fold <- daf_vec(daf_obj, "metacell", "zero_fold", required = FALSE)
    if (!is.null(zero_fold)) {
        result$zero_fold <- zero_fold
    }

    # Try to compute max_inner_fold from matrix if available
    inner_fold_mat <- daf_mat(daf_obj, "gene", "metacell", "inner_fold", required = FALSE)
    if (!is.null(inner_fold_mat) && is.null(result$max_inner_fold)) {
        result$max_inner_fold <- apply(inner_fold_mat, 2, max, na.rm = TRUE)
    }

    # Add other QC vectors if present
    qc_fields <- c("rare_gene_module", "is_rare")
    for (field in qc_fields) {
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

convert_daf_gene_qc <- function(daf_obj) {
    # Gene QC metadata
    gene_names <- dafr::axis_entries(daf_obj, "gene")

    result <- tibble(gene = gene_names)

    # Add max expression if available
    max_expr <- daf_vec(daf_obj, "gene", "max_expr", required = FALSE)
    if (is.null(max_expr)) {
        max_expr <- daf_vec(daf_obj, "gene", "mcview_cache_gene_max_umis", required = FALSE)
    }
    if (!is.null(max_expr)) {
        result$max_expr <- max_expr
    }

    # Add gene type if available
    gene_type <- daf_vec(daf_obj, "gene", "type", required = FALSE)
    if (!is.null(gene_type)) {
        result$type <- gene_type
    }

    # Add is_marker if available
    is_marker <- daf_vec(daf_obj, "gene", "is_marker", required = FALSE)
    if (!is.null(is_marker)) {
        result$is_marker <- is_marker
    }

    # Add significant_inner_folds_count if available
    sig_count <- daf_vec(daf_obj, "gene", "significant_inner_folds_count", required = FALSE)
    if (!is.null(sig_count)) {
        result$significant_inner_folds_count <- sig_count
    }

    # Compute max expression if inner_fold matrix available
    inner_fold_mat <- daf_mat(daf_obj, "gene", "metacell", "inner_fold", required = FALSE)
    if (!is.null(inner_fold_mat)) {
        max_inner_fold <- apply(inner_fold_mat, 1, max, na.rm = TRUE)
        result$max_inner_fold <- max_inner_fold
    }

    # Add other gene QC vectors if present
    qc_fields <- c("is_lateral", "is_noisy", "module", "correction_factor")
    for (field in qc_fields) {
        vec <- daf_vec(daf_obj, "gene", field, required = FALSE)
        if (!is.null(vec)) {
            result[[field]] <- vec
        }
    }

    # Include fitted gene metadata if present
    gene_props <- tryCatch(
        {
            daf_obj["/ gene ?"]
        },
        error = function(e) {
            character(0)
        }
    )
    fitted_fields <- grep("^fitted_gene_of", gene_props, value = TRUE)
    for (field in fitted_fields) {
        vec <- daf_vec(daf_obj, "gene", field, required = FALSE)
        if (!is.null(vec)) {
            result[[field]] <- vec
        }
    }

    if (ncol(result) == 1) {
        return(NULL)
    }

    return(result)
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

convert_daf_projected_fold <- function(daf_obj) {
    daf_mat(daf_obj, "gene", "metacell", "projected_fold", required = FALSE)
}

convert_daf_mc_mat_corrected <- function(daf_obj) {
    # Corrected fraction matrix (from atlas projection)
    corrected <- daf_mat(daf_obj, "metacell", "gene", "corrected_fraction", required = FALSE)
    if (is.null(corrected)) {
        return(NULL)
    }

    # Transpose to gene x metacell and convert to UMIs
    mc_sum <- daf_vec(daf_obj, "metacell", "total_UMIs")
    mc_mat_corrected <- Matrix::t(corrected) * mc_sum
    rownames(mc_mat_corrected) <- dafr::axis_entries(daf_obj, "gene")
    colnames(mc_mat_corrected) <- dafr::axis_entries(daf_obj, "metacell")

    return(mc_mat_corrected)
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
    # Projected (expected) fraction matrix (from atlas projection)
    projected <- daf_mat(daf_obj, "metacell", "gene", "projected_fraction", required = FALSE)
    if (is.null(projected)) {
        return(NULL)
    }

    # Transpose to gene x metacell and convert to UMIs
    mc_sum <- daf_vec(daf_obj, "metacell", "total_UMIs")
    projected_mat <- Matrix::t(projected) * mc_sum
    rownames(projected_mat) <- dafr::axis_entries(daf_obj, "gene")
    colnames(projected_mat) <- dafr::axis_entries(daf_obj, "metacell")

    return(projected_mat)
}

convert_daf_proj_weights <- function(daf_obj) {
    # Projection weights stored as JSON
    proj_weights_json <- daf_scalar(daf_obj, "projection_weights_json", default = NULL)
    if (is.null(proj_weights_json)) {
        return(NULL)
    }

    tryCatch({
        proj_weights <- jsonlite::fromJSON(proj_weights_json)
        as_tibble(proj_weights)
    }, error = function(e) {
        cli_warn("Failed to parse projection weights JSON: {e$message}")
        NULL
    })
}

convert_daf_query_cell_type_fracs <- function(daf_obj) {
    # Query cell type fractions stored as JSON
    fracs_json <- daf_scalar(daf_obj, "query_cell_type_fracs_json", default = NULL)
    if (is.null(fracs_json)) {
        return(NULL)
    }

    tryCatch({
        fracs <- jsonlite::fromJSON(fracs_json)
        as_tibble(fracs)
    }, error = function(e) {
        cli_warn("Failed to parse query cell type fractions JSON: {e$message}")
        NULL
    })
}
