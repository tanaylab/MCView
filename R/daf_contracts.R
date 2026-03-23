# daf_contracts.R - DAF Contract Definitions and Validation for MCView
#
# This file defines the formal data contract for MCView DAF data.
# Contracts specify required and optional axes, vectors, matrices, and scalars.
#
# Contract types:
#   - RequiredInput: Must exist before loading
#   - OptionalInput: Used if present, but not required
#
# Note: dafr now exposes contract validation. We keep the MCView contract
# definitions here and convert them to dafr contracts for verification.

# ==============================================================================
# Canonical Tab Names
# ==============================================================================

#' Canonical ordered list of MCView tab names
#'
#' Single source of truth for tab names used in contracts, validation,
#' and auto-detection. Tab definitions in app_config.R may include
#' additional UI-only tabs (e.g. "Projected-fold", "Query") that do
#' not have contracts.
#'
#' @noRd
MCVIEW_TAB_NAMES <- c(
    "About", "Manifold", "Genes", "Diff. Expression", "Cell types",
    "QC", "Markers", "Gene modules", "Gene correlation", "Inner-fold",
    "Stdev-fold", "Projection QC", "Atlas", "Annotate", "Samples", "Flow"
)

# ==============================================================================
# Contract Type Definitions
# ==============================================================================

#' Contract expectation types
#' @noRd
CONTRACT_REQUIRED <- "required"
CONTRACT_OPTIONAL <- "optional"

#' Create an axis specification
#' @noRd
axis_spec <- function(expectation, description) {
    list(expectation = expectation, description = description)
}

#' Create a vector specification
#' @noRd
vector_spec <- function(axis, name, expectation, type, description) {
    list(
        axis = axis,
        name = name,
        expectation = expectation,
        type = type,
        description = description
    )
}

#' Create a matrix specification
#' @noRd
matrix_spec <- function(rows_axis, cols_axis, name, expectation, type, description) {
    list(
        rows_axis = rows_axis,
        cols_axis = cols_axis,
        name = name,
        expectation = expectation,
        type = type,
        description = description
    )
}

#' Create a scalar specification
#' @noRd
scalar_spec <- function(name, expectation, type, description) {
    list(
        name = name,
        expectation = expectation,
        type = type,
        description = description
    )
}

# ==============================================================================
# Contract Conversion Helpers
# ==============================================================================

#' Map MCView contract expectation to dafr expectation
#' @noRd
mcview_expectation_to_dafr <- function(expectation) {
    switch(
        expectation,
        "required" = dafr::RequiredInput,
        "optional" = dafr::OptionalInput,
        stop("Unknown contract expectation: ", expectation, call. = FALSE)
    )
}

#' Convert an MCView contract to a dafr contract
#' @noRd
mcview_contract_to_dafr <- function(contract) {
    axes <- list()
    if (!is.null(contract$axes)) {
        for (axis_name in names(contract$axes)) {
            spec <- contract$axes[[axis_name]]
            axes <- c(
                axes,
                list(dafr::axis_contract(
                    axis_name,
                    mcview_expectation_to_dafr(spec$expectation),
                    spec$description
                ))
            )
        }
    }

    data <- list()
    all_vectors <- c(contract$vectors %||% list(), contract$graph_vectors %||% list())
    for (spec in all_vectors) {
        data <- c(
            data,
            list(dafr::vector_contract(
                spec$axis,
                spec$name,
                mcview_expectation_to_dafr(spec$expectation),
                spec$type,
                spec$description
            ))
        )
    }

    if (!is.null(contract$matrices)) {
        for (spec in contract$matrices) {
            data <- c(
                data,
                list(dafr::matrix_contract(
                    spec$rows_axis,
                    spec$cols_axis,
                    spec$name,
                    mcview_expectation_to_dafr(spec$expectation),
                    spec$type,
                    spec$description
                ))
            )
        }
    }

    if (!is.null(contract$scalars)) {
        for (spec in contract$scalars) {
            data <- c(
                data,
                list(dafr::scalar_contract(
                    spec$name,
                    mcview_expectation_to_dafr(spec$expectation),
                    spec$type,
                    spec$description
                ))
            )
        }
    }

    dafr::create_contract(axes = axes, data = data, is_relaxed = TRUE)
}

# ==============================================================================
# Core MCView Contract (Required for all tabs)
# ==============================================================================

#' Get the core MCView contract
#'
#' This contract specifies the minimum required data for MCView to function.
#'
#' @return List containing axes, vectors, matrices, and scalars specifications
#' @export
mcview_core_contract <- function() {
    list(
        name = "MCView Core",
        description = "Minimum required data for MCView",
        axes = list(
            metacell = axis_spec(CONTRACT_REQUIRED, "Metacell identifiers"),
            gene = axis_spec(CONTRACT_REQUIRED, "Gene identifiers"),
            type = axis_spec(CONTRACT_REQUIRED, "Cell type identifiers")
        ),
        vectors = list(
            # Required metacell vectors
            vector_spec(
                "metacell", "total_UMIs", CONTRACT_REQUIRED, "numeric",
                "Total UMIs per metacell"
            ),
            vector_spec(
                "metacell", "type", CONTRACT_REQUIRED, "character",
                "Cell type assignment per metacell"
            ),

            # Required type vectors
            vector_spec(
                "type", "color", CONTRACT_REQUIRED, "character",
                "Hex color per cell type (#RRGGBB format)"
            ),

            # 2D coordinates - special handling: need at least one pair
            vector_spec(
                "metacell", "x", CONTRACT_OPTIONAL, "numeric",
                "X coordinate for 2D projection"
            ),
            vector_spec(
                "metacell", "y", CONTRACT_OPTIONAL, "numeric",
                "Y coordinate for 2D projection"
            ),
            vector_spec(
                "metacell", "u", CONTRACT_OPTIONAL, "numeric",
                "U coordinate (alternative to x)"
            ),
            vector_spec(
                "metacell", "v", CONTRACT_OPTIONAL, "numeric",
                "V coordinate (alternative to y)"
            )
        ),
        matrices = list(
            matrix_spec(
                "metacell", "gene", "UMIs", CONTRACT_REQUIRED, "numeric",
                "UMI counts per metacell-gene"
            )
        ),
        scalars = list(
            # Config scalars (all optional — can come from YAML config instead)
            scalar_spec("mcview_title", CONTRACT_OPTIONAL, "character",
                        "App title displayed in UI header"),
            scalar_spec("mcview_tabs", CONTRACT_OPTIONAL, "character",
                        "Comma-separated list of enabled tab names"),
            scalar_spec("mcview_excluded_tabs", CONTRACT_OPTIONAL, "character",
                        "Comma-separated list of excluded tab names"),
            scalar_spec("mcview_light_version", CONTRACT_OPTIONAL, "logical",
                        "Use light version of UI"),
            scalar_spec("mcview_about_markdown", CONTRACT_OPTIONAL, "character",
                        "Markdown content for About tab"),
            scalar_spec("mcview_cache_in_daf", CONTRACT_OPTIONAL, "logical",
                        "Store runtime cache in writable DAF layer"),
            scalar_spec("mcview_cache_daf_root", CONTRACT_OPTIONAL, "character",
                        "Path to writable cache DAF directory"),

            # Pre-computed metadata (strongly recommended for fast startup)
            scalar_spec("mcview_available_tabs", CONTRACT_OPTIONAL, "character",
                        "Comma-separated auto-detected tabs; eliminates ~1.3s detect_available_tabs on cold start")
        ),

        # Special validation rules
        custom_rules = list(
            coordinates = list(
                description = "Must have either (x,y) or (u,v) coordinates",
                validate = function(daf_obj) {
                    has_xy <- dafr::has_vector(daf_obj, "metacell", "x") &&
                        dafr::has_vector(daf_obj, "metacell", "y")
                    has_uv <- dafr::has_vector(daf_obj, "metacell", "u") &&
                        dafr::has_vector(daf_obj, "metacell", "v")
                    if (!has_xy && !has_uv) {
                        return("Missing 2D coordinates: need either (x,y) or (u,v)")
                    }
                    return(NULL)
                }
            )
        )
    )
}

# ==============================================================================
# Tab-Specific Contracts
# ==============================================================================

#' Contract for Manifold tab
#' @export
mcview_manifold_contract <- function() {
    list(
        name = "Manifold",
        description = "2D projection and cell type visualization",
        extends = "core",
        vectors = list(
            # Optional enhancements
            vector_spec(
                "metacell", "n_cell", CONTRACT_OPTIONAL, "integer",
                "Number of cells per metacell"
            ),
            vector_spec(
                "metacell", "mc_col", CONTRACT_OPTIONAL, "character",
                "Per-metacell color override"
            )
        ),
        axes = list(
            metacell_graph = axis_spec(CONTRACT_OPTIONAL, "Metacell graph edges for visualization")
        ),

        # Graph edge data
        graph_vectors = list(
            vector_spec(
                "metacell_graph", "from", CONTRACT_OPTIONAL, "character",
                "Source metacell of edge"
            ),
            vector_spec(
                "metacell_graph", "to", CONTRACT_OPTIONAL, "character",
                "Target metacell of edge"
            ),
            vector_spec(
                "metacell_graph", "weight", CONTRACT_OPTIONAL, "numeric",
                "Edge weight"
            )
        )
    )
}

#' Contract for Genes tab
#' @export
mcview_genes_contract <- function() {
    list(
        name = "Genes",
        description = "Gene expression visualization",
        extends = "core",
        vectors = list(
            vector_spec(
                "metacell", "top1_gene", CONTRACT_OPTIONAL, "character",
                "Top expressed gene per metacell"
            ),
            vector_spec(
                "metacell", "top2_gene", CONTRACT_OPTIONAL, "character",
                "Second top expressed gene per metacell"
            ),
            vector_spec(
                "metacell", "top1_lfp", CONTRACT_OPTIONAL, "numeric",
                "Log fold-enrichment of top gene"
            ),
            vector_spec(
                "metacell", "top2_lfp", CONTRACT_OPTIONAL, "numeric",
                "Log fold-enrichment of second top gene"
            ),
            vector_spec(
                "gene", "is_lateral", CONTRACT_OPTIONAL, "logical",
                "Whether gene is lateral (excluded from analysis)"
            ),
            vector_spec(
                "gene", "is_noisy", CONTRACT_OPTIONAL, "logical",
                "Whether gene is noisy"
            )
        )
    )
}

#' Contract for Differential Expression tab
#' @export
mcview_diff_expr_contract <- function() {
    list(
        name = "Diff. Expression",
        description = "Differential expression analysis",
        extends = "core",
        vectors = list(
            vector_spec(
                "gene", "is_lateral", CONTRACT_OPTIONAL, "logical",
                "Whether gene is lateral"
            ),
            vector_spec(
                "gene", "is_noisy", CONTRACT_OPTIONAL, "logical",
                "Whether gene is noisy"
            )
        )
    )
}

#' Contract for Cell Types tab
#' @export
mcview_cell_types_contract <- function() {
    list(
        name = "Cell types",
        description = "Cell type annotation and analysis",
        extends = "core",
        vectors = list(
            vector_spec(
                "metacell", "n_cell", CONTRACT_OPTIONAL, "integer",
                "Number of cells per metacell"
            )
        )
    )
}

#' Contract for QC tab
#' @export
mcview_qc_contract <- function() {
    list(
        name = "QC",
        description = "Quality control metrics",
        extends = "core",
        vectors = list(
            vector_spec(
                "metacell", "n_cell", CONTRACT_OPTIONAL, "integer",
                "Number of cells per metacell"
            ),
            vector_spec(
                "metacell", "cells", CONTRACT_OPTIONAL, "integer",
                "Alias for n_cell"
            ),
            vector_spec(
                "metacell", "umis", CONTRACT_OPTIONAL, "numeric",
                "Alias for total_UMIs"
            ),
            vector_spec(
                "metacell", "max_inner_fold", CONTRACT_OPTIONAL, "numeric",
                "Maximum inner fold per metacell"
            ),
            vector_spec(
                "metacell", "max_inner_fold_no_lateral", CONTRACT_OPTIONAL, "numeric",
                "Max inner fold excluding lateral genes"
            ),
            vector_spec(
                "metacell", "max_inner_stdev_log", CONTRACT_OPTIONAL, "numeric",
                "Max inner stdev log per metacell"
            ),
            vector_spec(
                "metacell", "zero_fold", CONTRACT_OPTIONAL, "numeric",
                "Zero fold per metacell"
            ),
            vector_spec(
                "metacell", "rare_metacell", CONTRACT_OPTIONAL, "logical",
                "Whether metacell is rare"
            ),
            vector_spec(
                "metacell", "metacells_rare_gene_module", CONTRACT_OPTIONAL, "integer",
                "Rare gene module assignment"
            )
        ),
        scalars = list(
            scalar_spec(
                "qc_stats_n_outliers", CONTRACT_OPTIONAL, "integer",
                "Number of outlier cells"
            ),
            scalar_spec(
                "qc_stats_n_cells", CONTRACT_OPTIONAL, "integer",
                "Total number of cells"
            ),
            scalar_spec(
                "qc_stats_n_umis", CONTRACT_OPTIONAL, "integer",
                "Total UMI count"
            ),
            scalar_spec(
                "qc_stats_median_umis_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "Median UMIs per metacell"
            ),
            scalar_spec(
                "qc_stats_median_cells_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "Median cells per metacell"
            )
        )
    )
}

#' Contract for Markers tab
#' @export
mcview_markers_contract <- function() {
    list(
        name = "Markers",
        description = "Marker gene visualization",
        extends = "core",
        vectors = list(
            vector_spec(
                "gene", "is_marker", CONTRACT_REQUIRED, "logical",
                "Whether gene is a marker"
            )
        )
    )
}

#' Contract for Gene Modules tab
#' @export
mcview_gene_modules_contract <- function() {
    list(
        name = "Gene modules",
        description = "Gene module analysis",
        extends = "core",
        vectors = list(
            vector_spec(
                "gene", "module", CONTRACT_REQUIRED, "character",
                "Gene module assignment"
            )
        )
    )
}

#' Contract for Inner-fold tab
#' @export
mcview_inner_fold_contract <- function() {
    list(
        name = "Inner-fold",
        description = "Inner fold analysis per gene",
        extends = "core",
        matrices = list(
            matrix_spec(
                "gene", "metacell", "inner_fold", CONTRACT_REQUIRED, "numeric",
                "Inner fold per gene-metacell"
            )
        ),
        vectors = list(
            vector_spec(
                "gene", "is_lateral", CONTRACT_OPTIONAL, "logical",
                "Whether gene is lateral"
            ),
            vector_spec(
                "metacell", "max_inner_fold", CONTRACT_OPTIONAL, "numeric",
                "Maximum inner fold per metacell"
            )
        )
    )
}

#' Contract for Stdev-fold tab
#' @export
mcview_stdev_fold_contract <- function() {
    list(
        name = "Stdev-fold",
        description = "Inner standard deviation analysis",
        extends = "core",
        matrices = list(
            matrix_spec(
                "gene", "metacell", "inner_stdev_log", CONTRACT_REQUIRED, "numeric",
                "Inner stdev log per gene-metacell"
            )
        )
    )
}

#' Contract for Projection QC tab
#' @export
mcview_projection_qc_contract <- function() {
    list(
        name = "Projection QC",
        description = "Atlas projection quality control",
        extends = "core",
        vectors = list(
            vector_spec(
                "metacell", "projected_type", CONTRACT_OPTIONAL, "character",
                "Projected type from atlas"
            ),
            vector_spec(
                "metacell", "projected_correlation", CONTRACT_OPTIONAL, "numeric",
                "Correlation with atlas projection"
            ),
            vector_spec(
                "metacell", "similar", CONTRACT_OPTIONAL, "logical",
                "Whether similar to atlas"
            ),
            vector_spec(
                "metacell", "atlas_metacell", CONTRACT_OPTIONAL, "character",
                "Most similar atlas metacell"
            )
        ),
        matrices = list(
            matrix_spec(
                "gene", "metacell", "projected_fold", CONTRACT_OPTIONAL, "numeric",
                "Projected fold from atlas"
            ),
            matrix_spec(
                "metacell", "gene", "projected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Projected expression fraction from atlas"
            ),
            matrix_spec(
                "metacell", "gene", "corrected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Corrected expression fraction from atlas projection"
            )
        ),
        scalars = list(
            scalar_spec(
                "project_max_projection_fold_factor", CONTRACT_OPTIONAL, "numeric",
                "Maximum projection fold factor"
            ),
            scalar_spec(
                "projection_weights_json", CONTRACT_OPTIONAL, "character",
                "Projection weights as JSON (query-atlas metacell mapping)"
            ),
            scalar_spec(
                "query_cell_type_fracs_json", CONTRACT_OPTIONAL, "character",
                "Query cell type fractions as JSON"
            )
        )
    )
}

#' Contract for Atlas tab
#'
#' The query DAF should contain projection metadata, while the atlas DAF
#' (loaded via init_atlas) contains the reference data.
#'
#' @export
mcview_atlas_contract <- function() {
    list(
        name = "Atlas",
        description = "Atlas visualization (requires separate atlas DAF)",
        extends = "core",

        # Atlas requires projection metadata in query DAF
        vectors = list(
            vector_spec(
                "metacell", "projected_type", CONTRACT_OPTIONAL, "character",
                "Projected type from atlas"
            ),
            vector_spec(
                "metacell", "projected_correlation", CONTRACT_OPTIONAL, "numeric",
                "Correlation with atlas projection"
            ),
            vector_spec(
                "metacell", "similar", CONTRACT_OPTIONAL, "logical",
                "Whether similar to atlas"
            )
        ),
        matrices = list(
            matrix_spec(
                "metacell", "gene", "corrected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Corrected fraction from atlas projection"
            ),
            matrix_spec(
                "metacell", "gene", "projected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Projected expression fraction from atlas"
            )
        ),
        scalars = list(
            scalar_spec(
                "projection_weights_json", CONTRACT_OPTIONAL, "character",
                "Projection weights as JSON (query-atlas metacell mapping)"
            ),
            scalar_spec(
                "query_cell_type_fracs_json", CONTRACT_OPTIONAL, "character",
                "Query cell type fractions as JSON"
            )
        )
    )
}

#' Contract for Annotate tab
#' @export
mcview_annotate_contract <- function() {
    list(
        name = "Annotate",
        description = "Interactive annotation (always available)",
        extends = "core",
        vectors = list() # No additional requirements
    )
}

#' Contract for gene correlations data
#' @export
mcview_gene_correlations_contract <- function() {
    list(
        name = "Gene Correlations",
        description = "Gene-gene correlation data",
        axes = list(
            gg_mc_top_cor = axis_spec(CONTRACT_OPTIONAL, "Gene correlation pairs")
        ),
        vectors = list(
            vector_spec(
                "gg_mc_top_cor", "gene1", CONTRACT_OPTIONAL, "character",
                "First gene in correlation pair"
            ),
            vector_spec(
                "gg_mc_top_cor", "gene2", CONTRACT_OPTIONAL, "character",
                "Second gene in correlation pair"
            ),
            vector_spec(
                "gg_mc_top_cor", "cor", CONTRACT_OPTIONAL, "numeric",
                "Correlation value"
            ),
            vector_spec(
                "gg_mc_top_cor", "type", CONTRACT_OPTIONAL, "character",
                "Cell type for correlation"
            )
        )
    )
}

#' Contract for gene zero-fold data
#' @export
mcview_gene_zero_fold_contract <- function() {
    list(
        name = "Gene Zero Fold",
        description = "Gene zero-fold statistics",
        axes = list(
            gene_zero_fold = axis_spec(CONTRACT_OPTIONAL, "Gene zero-fold statistics")
        ),
        vectors = list(
            vector_spec(
                "gene_zero_fold", "gene", CONTRACT_OPTIONAL, "character",
                "Gene name"
            ),
            vector_spec(
                "gene_zero_fold", "metacell", CONTRACT_OPTIONAL, "character",
                "Metacell name"
            ),
            vector_spec(
                "gene_zero_fold", "zero_fold", CONTRACT_OPTIONAL, "numeric",
                "Zero fold value"
            ),
            vector_spec(
                "gene_zero_fold", "avg", CONTRACT_OPTIONAL, "numeric",
                "Average value"
            ),
            vector_spec(
                "gene_zero_fold", "obs", CONTRACT_OPTIONAL, "numeric",
                "Observed value"
            ),
            vector_spec(
                "gene_zero_fold", "exp", CONTRACT_OPTIONAL, "numeric",
                "Expected value"
            ),
            vector_spec(
                "gene_zero_fold", "type", CONTRACT_OPTIONAL, "character",
                "Cell type"
            )
        )
    )
}

#' Contract for MCView configuration scalars
#' @export
mcview_config_contract <- function() {
    list(
        name = "MCView Configuration",
        description = "MCView-specific configuration stored in DAF",
        scalars = list(
            scalar_spec(
                "mcview_title", CONTRACT_OPTIONAL, "character",
                "Application title"
            ),
            scalar_spec(
                "mcview_tabs", CONTRACT_OPTIONAL, "character",
                "Comma-separated list of tabs to show"
            ),
            scalar_spec(
                "mcview_excluded_tabs", CONTRACT_OPTIONAL, "character",
                "Comma-separated list of tabs to exclude"
            ),
            scalar_spec(
                "mcview_light_version", CONTRACT_OPTIONAL, "logical",
                "Whether to use light version"
            ),
            scalar_spec(
                "mcview_about_markdown", CONTRACT_OPTIONAL, "character",
                "About page markdown content"
            ),
            scalar_spec(
                "mcview_cache_in_daf", CONTRACT_OPTIONAL, "logical",
                "Whether to cache computed data in DAF"
            ),
            scalar_spec(
                "mcview_cache_daf_root", CONTRACT_OPTIONAL, "character",
                "Root directory for cache DAF"
            ),
            scalar_spec(
                "metacells_algorithm", CONTRACT_OPTIONAL, "character",
                "Metacells algorithm version"
            ),
            scalar_spec(
                "mcview_available_tabs", CONTRACT_OPTIONAL, "character",
                "Pre-stored comma-separated list of available tabs (skips tab detection)"
            )
        )
    )
}

#' Contract for Samples tab
#'
#' The Samples tab requires cell-level metadata with sample assignments.
#' This is optional as not all datasets have cell-level data.
#'
#' @export
mcview_samples_contract <- function() {
    list(
        name = "Samples",
        description = "Sample-level analysis requiring cell metadata",
        extends = "core",
        axes = list(
            cell = axis_spec(CONTRACT_REQUIRED, "Cell identifiers")
        ),
        vectors = list(
            # Cell -> metacell mapping
            vector_spec(
                "cell", "metacell", CONTRACT_REQUIRED, "character",
                "Metacell assignment per cell (-1 for outliers)"
            ),

            # Cell -> sample mapping
            vector_spec(
                "cell", "samp_id", CONTRACT_REQUIRED, "character",
                "Sample ID per cell"
            ),

            # Optional cell metadata
            vector_spec(
                "cell", "cell_type", CONTRACT_OPTIONAL, "character",
                "Cell type per cell (if different from metacell type)"
            ),
            vector_spec(
                "cell", "outlier", CONTRACT_OPTIONAL, "logical",
                "Whether cell is an outlier"
            )
        ),

        # Note: Additional cell metadata columns are allowed and will be displayed
        custom_rules = list(
            has_samples = list(
                description = "Cells must have sample assignments",
                validate = function(daf_obj) {
                    if (!dafr::has_axis(daf_obj, "cell")) {
                        return("Missing cell axis")
                    }
                    if (!dafr::has_vector(daf_obj, "cell", "samp_id")) {
                        return("Missing cell.samp_id vector")
                    }
                    return(NULL)
                }
            )
        )
    )
}

#' Contract for Flow tab
#'
#' The Flow tab requires time annotation data.
#'
#' @export
mcview_flow_contract <- function() {
    list(
        name = "Flow",
        description = "Time/flow analysis",
        extends = "core",
        vectors = list(
            vector_spec(
                "metacell", "time", CONTRACT_OPTIONAL, "numeric",
                "Time point per metacell"
            ),
            vector_spec(
                "metacell", "flow_score", CONTRACT_OPTIONAL, "numeric",
                "Flow score per metacell"
            )
        ),

        # Time annotation data stored in metadata
        custom_rules = list(
            has_time = list(
                description = "Must have time annotation data",
                validate = function(daf_obj) {
                    # Check for time-related vectors
                    if (dafr::has_vector(daf_obj, "metacell", "time")) {
                        return(NULL)
                    }
                    return("Missing time annotation data")
                }
            )
        )
    )
}

#' Contract for precomputed cache vectors
#'
#' These vectors are computed by MCView and cached for faster loading.
#' They are optional - MCView will compute them at runtime if missing,
#' but precomputing them improves first-load performance.
#'
#' @export
mcview_cache_contract <- function() {
    list(
        name = "MCView Cache",
        description = paste(
            "Precomputed cache vectors for faster loading.",
            "These are computed at runtime if missing, but precomputing",
            "improves performance on first load."
        ),
        vectors = list(
            # Top genes per metacell (precomputed for hover info)
            vector_spec(
                "metacell", "mcview_cache_top1_gene", CONTRACT_OPTIONAL, "character",
                "Precomputed top expressed gene"
            ),
            vector_spec(
                "metacell", "mcview_cache_top2_gene", CONTRACT_OPTIONAL, "character",
                "Precomputed second top expressed gene"
            ),
            vector_spec(
                "metacell", "mcview_cache_top1_lfp", CONTRACT_OPTIONAL, "numeric",
                "Precomputed log fold-enrichment of top gene"
            ),
            vector_spec(
                "metacell", "mcview_cache_top2_lfp", CONTRACT_OPTIONAL, "numeric",
                "Precomputed log fold-enrichment of second gene"
            ),

            # Gene statistics
            vector_spec(
                "gene", "mcview_cache_max_expr", CONTRACT_OPTIONAL, "numeric",
                "Precomputed maximum expression per gene"
            ),
            vector_spec(
                "gene", "mcview_cache_total_expr", CONTRACT_OPTIONAL, "numeric",
                "Precomputed total expression per gene"
            )
        ),
        axes = list(
            # Precomputed correlation data axis
            mcview_cache_gg_mc_top_cor = axis_spec(
                CONTRACT_OPTIONAL,
                "Precomputed gene correlation pairs"
            )
        )
    )
}

#' Precompute cache vectors for a DAF
#'
#' Computes and stores cache vectors in the DAF for faster loading.
#' This should be run once after data preparation.
#'
#' @param daf_obj DAF object (must be writable)
#' @param verbose If TRUE, print progress
#'
#' @return Modified DAF object with cache vectors
#' @export
precompute_mcview_cache <- function(daf_obj, verbose = TRUE) {
    if (verbose) cli::cli_alert_info("Precomputing MCView cache vectors...")

    # Try Julia-accelerated path for top genes
    jl_top <- julia_precompute_top_genes(daf_obj, egc_epsilon = 1e-5)
    if (!is.null(jl_top)) {
        if (verbose) cli::cli_alert_info("Using Julia for top genes computation")
        daf_obj <- dafr::set_vector(daf_obj, "metacell", "mcview_cache_top1_gene", jl_top$top1_gene)
        daf_obj <- dafr::set_vector(daf_obj, "metacell", "mcview_cache_top2_gene", jl_top$top2_gene)
        daf_obj <- dafr::set_vector(daf_obj, "metacell", "mcview_cache_top1_lfp", jl_top$top1_lfp)
        daf_obj <- dafr::set_vector(daf_obj, "metacell", "mcview_cache_top2_lfp", jl_top$top2_lfp)
    } else {
        # R fallback for top genes
        mc_mat <- dafr::get_matrix(daf_obj, "metacell", "gene", "UMIs")
        mc_sum <- dafr::get_vector(daf_obj, "metacell", "total_UMIs")

        egc <- sweep(mc_mat, 1, mc_sum, "/")
        gene_means <- colMeans(egc)
        gene_means[gene_means == 0] <- 1e-10
        lfp <- log2(sweep(egc, 2, gene_means, "/") + 1e-5)
        lfp <- as.matrix(lfp)

        gene_names <- colnames(mc_mat)

        # Find top-2 genes per metacell without mutation.
        # lfp is metacell x gene; jlview_top2_per_row finds top-2 columns
        # (genes) per row (metacell).
        tops <- jlview::jlview_top2_per_row(lfp)
        top1_gene <- gene_names[tops$top1_idx]
        top2_gene <- gene_names[tops$top2_idx]
        top1_lfp <- tops$top1_val
        top2_lfp <- tops$top2_val

        daf_obj <- dafr::set_vector(daf_obj, "metacell", "mcview_cache_top1_gene", top1_gene)
        daf_obj <- dafr::set_vector(daf_obj, "metacell", "mcview_cache_top2_gene", top2_gene)
        daf_obj <- dafr::set_vector(daf_obj, "metacell", "mcview_cache_top1_lfp", top1_lfp)
        daf_obj <- dafr::set_vector(daf_obj, "metacell", "mcview_cache_top2_lfp", top2_lfp)
    }

    # Try Julia-accelerated path for gene statistics
    jl_stats <- julia_precompute_gene_stats(daf_obj)
    if (!is.null(jl_stats)) {
        if (verbose) cli::cli_alert_info("Using Julia for gene statistics computation")
        daf_obj <- dafr::set_vector(daf_obj, "gene", "mcview_cache_max_expr", jl_stats$max_expr)
        daf_obj <- dafr::set_vector(daf_obj, "gene", "mcview_cache_total_expr", jl_stats$total_expr)
    } else {
        # R fallback for gene statistics
        if (!exists("egc", inherits = FALSE)) {
            mc_mat <- dafr::get_matrix(daf_obj, "metacell", "gene", "UMIs")
            mc_sum <- dafr::get_vector(daf_obj, "metacell", "total_UMIs")
            egc <- sweep(mc_mat, 1, mc_sum, "/")
            gene_names <- colnames(mc_mat)
        }
        max_expr <- matrixStats::colMaxs(egc)
        names(max_expr) <- gene_names
        total_expr <- colSums(egc)
        daf_obj <- dafr::set_vector(daf_obj, "gene", "mcview_cache_max_expr", max_expr)
        daf_obj <- dafr::set_vector(daf_obj, "gene", "mcview_cache_total_expr", total_expr)
    }

    if (verbose) cli::cli_alert_success("Cache vectors precomputed")

    return(daf_obj)
}

# ==============================================================================
# Get All Contracts
# ==============================================================================

#' Get all MCView contracts
#'
#' @return Named list of all contract definitions
#' @export
mcview_all_contracts <- function() {
    list(
        core = mcview_core_contract(),
        manifold = mcview_manifold_contract(),
        genes = mcview_genes_contract(),
        diff_expr = mcview_diff_expr_contract(),
        cell_types = mcview_cell_types_contract(),
        qc = mcview_qc_contract(),
        markers = mcview_markers_contract(),
        gene_modules = mcview_gene_modules_contract(),
        inner_fold = mcview_inner_fold_contract(),
        stdev_fold = mcview_stdev_fold_contract(),
        projection_qc = mcview_projection_qc_contract(),
        atlas = mcview_atlas_contract(),
        annotate = mcview_annotate_contract(),
        samples = mcview_samples_contract(),
        flow = mcview_flow_contract(),
        gene_correlations = mcview_gene_correlations_contract(),
        gene_zero_fold = mcview_gene_zero_fold_contract(),
        config = mcview_config_contract(),
        cache = mcview_cache_contract()
    )
}

#' Get contract for a specific tab
#'
#' @param tab_name Name of the tab
#' @return Contract definition or NULL if not found
#' @export
mcview_tab_contract <- function(tab_name) {
    tab_map <- list(
        "About" = mcview_core_contract(),
        "Manifold" = mcview_manifold_contract(),
        "Genes" = mcview_genes_contract(),
        "Diff. Expression" = mcview_diff_expr_contract(),
        "Cell types" = mcview_cell_types_contract(),
        "QC" = mcview_qc_contract(),
        "Markers" = mcview_markers_contract(),
        "Gene modules" = mcview_gene_modules_contract(),
        "Gene correlation" = mcview_gene_correlations_contract(),
        "Inner-fold" = mcview_inner_fold_contract(),
        "Stdev-fold" = mcview_stdev_fold_contract(),
        "Projection QC" = mcview_projection_qc_contract(),
        "Atlas" = mcview_atlas_contract(),
        "Annotate" = mcview_annotate_contract(),
        "Samples" = mcview_samples_contract(),
        "Flow" = mcview_flow_contract()
    )

    tab_map[[tab_name]]
}

# ==============================================================================
# Contract Validation Functions
# ==============================================================================

#' Validate DAF against MCView contract
#'
#' Validates that a DAF object meets the MCView contract requirements.
#'
#' @param daf_obj DAF object to validate
#' @param contract Contract to validate against (default: core contract)
#' @param strict If TRUE, fail on any missing optional data (default: FALSE)
#' @param verbose If TRUE, print detailed validation info
#'
#' @return List with `valid` (logical), `errors` (character vector),
#'         `warnings` (character vector), and `info` (character vector)
#' @export
validate_mcview_contract <- function(daf_obj,
                                     contract = mcview_core_contract(),
                                     strict = FALSE,
                                     verbose = FALSE) {
    result <- list(
        valid = TRUE,
        errors = character(),
        warnings = character(),
        info = character()
    )

    if (!inherits(daf_obj, "Daf")) {
        result$valid <- FALSE
        result$errors <- c(result$errors, "Object is not a valid DAF object")
        return(result)
    }

    daf_contract <- mcview_contract_to_dafr(contract)
    computation <- if (!is.null(contract$name)) contract$name else "MCView"
    contract_daf <- dafr::contractor(paste0("MCView.", computation), daf_contract, daf_obj)

    verify_result <- dafr::verify_contract(contract_daf$daf, contract_daf$contract)
    if (!verify_result$valid) {
        result$valid <- FALSE
        result$errors <- c(result$errors, verify_result$errors)
    }
    if (length(verify_result$warnings) > 0) {
        result$warnings <- c(result$warnings, verify_result$warnings)
    }

    if (!is.null(contract$axes)) {
        for (axis_name in names(contract$axes)) {
            spec <- contract$axes[[axis_name]]
            has_it <- dafr::has_axis(daf_obj, axis_name)

            if (spec$expectation == CONTRACT_OPTIONAL && !has_it && strict) {
                result$warnings <- c(
                    result$warnings,
                    glue::glue("Missing optional axis: {axis_name}")
                )
            } else if (has_it && verbose) {
                n <- dafr::axis_length(daf_obj, axis_name)
                result$info <- c(
                    result$info,
                    glue::glue("Axis '{axis_name}': {n} entries")
                )
            }
        }
    }

    all_vectors <- c(
        contract$vectors %||% list(),
        contract$graph_vectors %||% list()
    )

    for (spec in all_vectors) {
        has_it <- dafr::has_axis(daf_obj, spec$axis) &&
            dafr::has_vector(daf_obj, spec$axis, spec$name)

        if (spec$expectation == CONTRACT_OPTIONAL && !has_it && strict) {
            result$warnings <- c(
                result$warnings,
                glue::glue("Missing optional vector: {spec$axis}.{spec$name}")
            )
        } else if (has_it && verbose) {
            result$info <- c(
                result$info,
                glue::glue("Vector '{spec$axis}.{spec$name}': present")
            )
        }
    }

    if (!is.null(contract$matrices)) {
        for (spec in contract$matrices) {
            has_it <- dafr::has_axis(daf_obj, spec$rows_axis) &&
                dafr::has_axis(daf_obj, spec$cols_axis) &&
                dafr::has_matrix(daf_obj, spec$rows_axis, spec$cols_axis, spec$name)

            if (spec$expectation == CONTRACT_OPTIONAL && !has_it && strict) {
                result$warnings <- c(
                    result$warnings,
                    glue::glue("Missing optional matrix: {spec$rows_axis},{spec$cols_axis}.{spec$name}")
                )
            } else if (has_it && verbose) {
                result$info <- c(
                    result$info,
                    glue::glue("Matrix '{spec$rows_axis},{spec$cols_axis}.{spec$name}': present")
                )
            }
        }
    }

    if (!is.null(contract$scalars)) {
        for (spec in contract$scalars) {
            has_it <- dafr::has_scalar(daf_obj, spec$name)

            if (spec$expectation == CONTRACT_OPTIONAL && !has_it && strict) {
                result$warnings <- c(
                    result$warnings,
                    glue::glue("Missing optional scalar: {spec$name}")
                )
            } else if (has_it && verbose) {
                result$info <- c(
                    result$info,
                    glue::glue("Scalar '{spec$name}': present")
                )
            }
        }
    }

    # Run custom validation rules
    if (!is.null(contract$custom_rules)) {
        for (rule_name in names(contract$custom_rules)) {
            rule <- contract$custom_rules[[rule_name]]
            error_msg <- rule$validate(daf_obj)
            if (!is.null(error_msg)) {
                result$valid <- FALSE
                result$errors <- c(result$errors, error_msg)
            } else if (verbose) {
                result$info <- c(
                    result$info,
                    glue::glue("Custom rule '{rule_name}': passed")
                )
            }
        }
    }

    return(result)
}

#' Validate DAF for all tabs
#'
#' Validates DAF against all tab-specific contracts and reports which tabs
#' can be enabled.
#'
#' @param daf_obj DAF object to validate
#' @param verbose If TRUE, print detailed validation info
#'
#' @return Named list with tab names as keys and validation results as values
#' @export
validate_mcview_tabs <- function(daf_obj, verbose = FALSE) {
    tabs <- MCVIEW_TAB_NAMES

    results <- list()

    # First validate core contract
    core_result <- validate_mcview_contract(daf_obj, mcview_core_contract(), verbose = verbose)

    if (!core_result$valid) {
        if (verbose) {
            cli::cli_alert_danger("Core contract validation failed:")
            for (err in core_result$errors) {
                cli::cli_alert_warning("  {err}")
            }
        }
        # If core fails, no tabs are available
        for (tab in tabs) {
            results[[tab]] <- list(available = FALSE, reason = "Core contract failed")
        }
        return(results)
    }

    # Validate each tab
    for (tab in tabs) {
        tab_contract <- mcview_tab_contract(tab)
        if (is.null(tab_contract)) {
            results[[tab]] <- list(available = TRUE, reason = "No specific contract")
            next
        }

        tab_result <- validate_mcview_contract(daf_obj, tab_contract, verbose = FALSE)

        if (tab_result$valid) {
            results[[tab]] <- list(available = TRUE, reason = "Contract satisfied")
            if (verbose) {
                cli::cli_alert_success("Tab '{tab}': available")
            }
        } else {
            # Check if it's just missing optional data for extended features
            required_errors <- grep("required", tab_result$errors, value = TRUE, ignore.case = TRUE)
            if (length(required_errors) == 0) {
                results[[tab]] <- list(available = TRUE, reason = "Optional data missing")
            } else {
                results[[tab]] <- list(
                    available = FALSE,
                    reason = paste(required_errors, collapse = "; ")
                )
                if (verbose) {
                    cli::cli_alert_warning("Tab '{tab}': unavailable - {results[[tab]]$reason}")
                }
            }
        }
    }

    return(results)
}

#' Print contract summary
#'
#' Prints a human-readable summary of a contract.
#'
#' @param contract Contract to summarize
#' @export
print_contract <- function(contract) {
    cli::cli_h1(contract$name)
    cli::cli_text(contract$description)

    if (!is.null(contract$extends)) {
        cli::cli_text("Extends: {contract$extends}")
    }

    if (length(contract$axes) > 0) {
        cli::cli_h2("Axes")
        for (name in names(contract$axes)) {
            spec <- contract$axes[[name]]
            cli::cli_li("{name} ({spec$expectation}): {spec$description}")
        }
    }

    vectors <- c(contract$vectors %||% list(), contract$graph_vectors %||% list())
    if (length(vectors) > 0) {
        cli::cli_h2("Vectors")
        for (spec in vectors) {
            cli::cli_li("{spec$axis}.{spec$name} ({spec$expectation}, {spec$type}): {spec$description}")
        }
    }

    if (length(contract$matrices) > 0) {
        cli::cli_h2("Matrices")
        for (spec in contract$matrices) {
            cli::cli_li("{spec$rows_axis},{spec$cols_axis}.{spec$name} ({spec$expectation}): {spec$description}")
        }
    }

    if (length(contract$scalars) > 0) {
        cli::cli_h2("Scalars")
        for (spec in contract$scalars) {
            cli::cli_li("{spec$name} ({spec$expectation}, {spec$type}): {spec$description}")
        }
    }

    if (!is.null(contract$custom_rules) && length(contract$custom_rules) > 0) {
        cli::cli_h2("Custom Rules")
        for (name in names(contract$custom_rules)) {
            rule <- contract$custom_rules[[name]]
            cli::cli_li("{name}: {rule$description}")
        }
    }
}

#' Generate contract documentation as markdown
#'
#' @param contract Contract to document
#' @return Character string with markdown content
#' @export
contract_to_markdown <- function(contract) {
    lines <- character()

    lines <- c(lines, glue::glue("# {contract$name}"))
    lines <- c(lines, "")
    lines <- c(lines, contract$description)
    lines <- c(lines, "")

    if (!is.null(contract$extends)) {
        lines <- c(lines, glue::glue("**Extends:** {contract$extends}"))
        lines <- c(lines, "")
    }

    if (length(contract$axes) > 0) {
        lines <- c(lines, "## Axes")
        lines <- c(lines, "")
        lines <- c(lines, "| Axis | Required | Description |")
        lines <- c(lines, "|------|----------|-------------|")
        for (name in names(contract$axes)) {
            spec <- contract$axes[[name]]
            req <- if (spec$expectation == CONTRACT_REQUIRED) "Yes" else "No"
            lines <- c(lines, glue::glue("| `{name}` | {req} | {spec$description} |"))
        }
        lines <- c(lines, "")
    }

    vectors <- c(contract$vectors %||% list(), contract$graph_vectors %||% list())
    if (length(vectors) > 0) {
        lines <- c(lines, "## Vectors")
        lines <- c(lines, "")
        lines <- c(lines, "| Axis | Name | Required | Type | Description |")
        lines <- c(lines, "|------|------|----------|------|-------------|")
        for (spec in vectors) {
            req <- if (spec$expectation == CONTRACT_REQUIRED) "Yes" else "No"
            lines <- c(lines, glue::glue("| `{spec$axis}` | `{spec$name}` | {req} | {spec$type} | {spec$description} |"))
        }
        lines <- c(lines, "")
    }

    if (length(contract$matrices) > 0) {
        lines <- c(lines, "## Matrices")
        lines <- c(lines, "")
        lines <- c(lines, "| Rows Axis | Cols Axis | Name | Required | Description |")
        lines <- c(lines, "|-----------|-----------|------|----------|-------------|")
        for (spec in contract$matrices) {
            req <- if (spec$expectation == CONTRACT_REQUIRED) "Yes" else "No"
            lines <- c(lines, glue::glue("| `{spec$rows_axis}` | `{spec$cols_axis}` | `{spec$name}` | {req} | {spec$description} |"))
        }
        lines <- c(lines, "")
    }

    if (length(contract$scalars) > 0) {
        lines <- c(lines, "## Scalars")
        lines <- c(lines, "")
        lines <- c(lines, "| Name | Required | Type | Description |")
        lines <- c(lines, "|------|----------|------|-------------|")
        for (spec in contract$scalars) {
            req <- if (spec$expectation == CONTRACT_REQUIRED) "Yes" else "No"
            lines <- c(lines, glue::glue("| `{spec$name}` | {req} | {spec$type} | {spec$description} |"))
        }
        lines <- c(lines, "")
    }

    paste(lines, collapse = "\n")
}

#' Generate full MCView contract documentation
#'
#' @param output_file Path to output markdown file (optional)
#' @return Character string with full documentation
#' @export
generate_contract_docs <- function(output_file = NULL) {
    contracts <- mcview_all_contracts()

    lines <- character()
    lines <- c(lines, "# MCView DAF Contract Specification")
    lines <- c(lines, "")
    lines <- c(lines, "This document specifies the data contract for MCView DAF files.")
    lines <- c(lines, "")
    lines <- c(lines, "## Overview")
    lines <- c(lines, "")
    lines <- c(lines, "MCView requires data to be stored in DAF (Data Axes Format) files.")
    lines <- c(lines, "The contract defines required and optional data for each feature.")
    lines <- c(lines, "")

    for (name in names(contracts)) {
        contract <- contracts[[name]]
        lines <- c(lines, contract_to_markdown(contract))
        lines <- c(lines, "---")
        lines <- c(lines, "")
    }

    result <- paste(lines, collapse = "\n")

    if (!is.null(output_file)) {
        writeLines(result, output_file)
        cli::cli_alert_success("Contract documentation written to {output_file}")
    }

    invisible(result)
}
