# daf_conversion.R - DAF -> MCView format converters (convert_daf_*)
#
# Split from R/daf_data.R (2026-05-01). The dispatcher convert_daf_to_mcview
# routes per-data-type converters; each per-type converter reads from the
# DAF (using R/daf_accessors.R + R/daf_queries.R) and reshapes into the
# tibble / matrix / list shapes the rest of MCView expects.

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
        # Known-optional keys that are not stored in the DAF: caller-side
        # fallbacks (calc_* in utils_cache, plot_vein/plot_network for the
        # Flow tab, get_metadata_colors for color overrides). Returning NULL
        # silently lets those fallbacks run without log noise. Anything not
        # listed here is treated as a typo and warns.
        "samp_mc_count" = ,
        "samp_mc_frac" = ,
        "samp_metadata" = ,
        "samp_list" = ,
        "default_markers" = convert_daf_default_markers(daf_obj),
        "default_markers_dist" = convert_daf_default_markers_dist(daf_obj),
        "metadata_colors" = ,
        "umap_anchors" = ,
        "outliers_metadata" = ,
        "project_max_projection_fold_factor" = ,
        "mc_ag" = ,
        "mc_network" = ,
        "mc_rank" = ,
        "mct_probs_trans" = ,
        "time_annot" = ,
        "type_ag" = ,
        "type_flow" = NULL,
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

    # Fall back: ask dafr to compute UMIs % Fraction. The eltwise `% Fraction`
    # divides each column (metacell) by its column sum, matching the historical
    # sweep(UMIs, 2, total_UMIs, "/") whenever total_UMIs == colSums(UMIs) -
    # the standard import-pipeline invariant. No R-side sweep, no extra
    # mc_mat materialisation.
    cli::cli_inform("Computing EGC as UMIs % Fraction (no fraction matrix in DAF)")
    daf_obj["@ gene @ metacell :: UMIs % Fraction"]
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

    # Add optional cell count: accept n_cells (canonical), n_cell (legacy
    # OBK-style), or cells (legacy converted projects like HEP_atlas_v5).
    n_cell <- try_daf_names(daf_obj, "metacell", c("n_cells", "n_cell", "cells"))
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

convert_daf_default_markers <- function(daf_obj) {
    s <- daf_scalar(daf_obj, "mcview_default_markers", default = NULL)
    if (is.null(s) || !nzchar(s)) return(NULL)
    strsplit(as.character(s), ",", fixed = TRUE)[[1L]]
}

convert_daf_default_markers_dist <- function(daf_obj) {
    if (!dafr::has_matrix(daf_obj, "metacell", "metacell",
                          "mcview_default_markers_dist")) {
        return(NULL)
    }
    dafr::get_matrix(daf_obj, "metacell", "metacell",
                     "mcview_default_markers_dist")
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

    # Drop MCView-derived cache vectors (mcview_cache_*, read through their
    # own accessors, not as metadata) and DAF-internal (double-underscore)
    # properties such as __zeros_downsample_UMIs. Surfacing them as
    # user-selectable metadata clutters the field selectors and can break
    # metadata plots (e.g. the character-valued mcview_cache_top1_gene).
    metadata_fields <- metadata_fields[!grepl("^(mcview_cache_|__)", metadata_fields)]

    if (length(metadata_fields) == 0) {
        return(NULL)
    }

    # Single dafr::get_dataframe call instead of N+1 daf_vec loops.
    df <- tryCatch(
        dafr::get_dataframe(daf_obj, "metacell", columns = metadata_fields),
        error = function(e) NULL
    )
    if (is.null(df) || ncol(df) == 0) {
        return(NULL)
    }
    result <- tibble::as_tibble(df) %>%
        tibble::add_column(metacell = metacell_names, .before = 1L)
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

    if (length(extra_fields) > 0) {
        # Batched fetch: one dafr::get_dataframe call instead of N daf_vec.
        # Matters most on cell DAFs with many properties (~60 on OBK cells).
        extra_df <- tryCatch(
            dafr::get_dataframe(daf_obj, "cell", columns = extra_fields),
            error = function(e) NULL
        )
        if (!is.null(extra_df) && ncol(extra_df) > 0) {
            for (field in colnames(extra_df)) {
                result[[field]] <- extra_df[[field]]
            }
        }
    }

    return(result)
}

convert_daf_mc_qc_metadata <- function(daf_obj) {
    metacell_names <- dafr::axis_entries(daf_obj, "metacell")
    result <- tibble(metacell = metacell_names)

    # Handle fields with fallback names
    result <- add_optional_vec_with_fallback(result, daf_obj, "metacell", "umis", "total_UMIs")
    # Cells-per-metacell may be exported as n_cells (canonical), n_cell (legacy
    # OBK-style), or cells (legacy converted projects like HEP_atlas_v5).
    n_cells_vec <- try_daf_names(daf_obj, "metacell", c("n_cells", "n_cell", "cells"))
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
        # Per-gene max via DAF reduction (OpenMP, no R-side densification).
        max_umis <- tryCatch(daf_query_gene_max_umis(daf_obj), error = function(e) NULL)
        if (!is.null(max_umis)) {
            result$max_expr <- as.numeric(max_umis[result$gene])
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

    # Fall back to derived values when scalars weren't pre-stored. Raw DAFs
    # (not produced by convert_project_to_daf) typically lack them, so we
    # compute medians from per-metacell vectors and cell totals from axes
    # rather than letting the QC value boxes silently req() out.
    if (is.null(stats$median_umis_per_metacell)) {
        total_umis <- tryCatch(
            dafr::get_vector(daf_obj, "metacell", "total_UMIs"),
            error = function(e) NULL
        )
        if (!is.null(total_umis) && length(total_umis) > 0) {
            stats$median_umis_per_metacell <- stats::median(total_umis, na.rm = TRUE)
        }
    }
    if (is.null(stats$median_cells_per_metacell)) {
        n_cells_vec <- try_daf_names(daf_obj, "metacell", c("n_cells", "n_cell", "cells"))
        if (!is.null(n_cells_vec) && length(n_cells_vec) > 0) {
            stats$median_cells_per_metacell <- stats::median(n_cells_vec, na.rm = TRUE)
        }
    }
    if (is.null(stats$n_cells)) {
        n_cells_vec <- try_daf_names(daf_obj, "metacell", c("n_cells", "n_cell", "cells"))
        if (!is.null(n_cells_vec) && length(n_cells_vec) > 0) {
            stats$n_cells <- sum(n_cells_vec, na.rm = TRUE)
        } else if (dafr::has_axis(daf_obj, "cell")) {
            stats$n_cells <- dafr::axis_length(daf_obj, "cell")
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
