# daf_contracts.R - DAF Contract Definitions and Validation for MCView
#
# CHAIN ARCHITECTURE
# ==================
# MCView uses a two-layer DAF chain:
#
#   mcview_derived (writable, MCView-specific precomputations)
#     |
#     v chains on top of
#   input DAF (read-only, from metacells pipeline)
#
# When querying the chain, DAF looks up the derived layer first, then falls
# through to the input layer. This means:
#   - The pipeline produces the input DAF with core biological data.
#   - MCView's derived layer adds precomputed caches (top genes, gene stats,
#     per-type marker tables, fallback marker correlations, etc.).
#   - Queries against the chain see both layers seamlessly.
#
# CONTRACT ORGANIZATION
# =====================
# 1. mcview_core_contract()     - what MCView expects from the input DAF
#    (RequiredInput / OptionalInput)
# 2. mcview_cache_contract()    - what the mcview_derived layer provides
# 3. Tab-specific contracts     - per-tab data requirements
# 4. mcview_config_contract()   - MCView configuration scalars
#
# Contract types:
#   - RequiredInput: Must exist before loading
#   - OptionalInput: Used if present, but not required
#
# Note: dafr now exposes contract validation. We keep the MCView contract
# definitions here and convert them to dafr contracts for verification.
# The Julia-side contracts are in inst/julia/mcview_contracts.jl.
#
# AUDIT CHANGES (2026-03-27)
# ==========================
# - Data classification: mcview_type_markers restructured from axis+vectors
#   to sparse matrices on (type, gene). mcview_marker_correlations kept as
#   axis (pragmatic compromise for 3D gene x gene x type data).
# - Naming: n_cells (Metacells.jl canonical), inner_std_log (canonical),
#   umap_x/umap_y aliases, primary_module, is_similar, etc.
# - QC stats: removed qc_stats_ prefix, added aliases for backward compat.
# - Removed MC x MC matrices (outgoing_weights, obs_balanced_ranks,
#   umap_distances) -- no R code reads these.
# - New entries: type order, marker_rank, significant_inner_folds_count,
#   linear_fraction.
# - gg_mc_top_cor deprecated (kept for backward compat only).

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
    "QC", "Markers", "Gene modules", "Gene correlation",
    "Projection QC", "Atlas", "Annotate", "Samples", "Flow"
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

#' Get the core MCView contract (INPUT layer)
#'
#' This contract specifies what MCView expects from the input DAF (read-only,
#' from the metacells pipeline). Items marked RequiredInput must exist;
#' items marked OptionalInput enrich specific tabs when present.
#'
#' Coordinate note: at least one pair of (x,y) or (u,v) is required for the
#' Manifold tab. The contract marks all four as OptionalInput because the
#' "at least one pair" constraint is enforced by a custom rule.
#'
#' Sample ID note: `mcview_sample_property` is a scalar that names which cell
#' vector to use as the sample identifier (e.g. "embryo", "batch_set_id").
#' There is no hardcoded `samp_id` vector -- the cell vector name is dynamic.
#'
#' Inner stdev note: the pipeline may produce either `inner_stdev_log` or
#' `inner_std_log`. MCView accepts both names; the R layer normalizes to
#' whichever is present. `inner_std_log` is the canonical Metacells.jl name;
#' `inner_stdev_log` is kept as a backward-compatible alias.
#'
#' n_cells note: Metacells.jl uses ("metacell", "n_cells"). MCView also accepts
#' the legacy names ("metacell", "n_cell"). The old ("metacell", "cells") alias
#' has been removed.
#'
#' geomean_fraction vs linear_fraction: these are DIFFERENT computations.
#'   - linear_fraction = UMIs / total_UMIs (Metacells.jl canonical)
#'   - geomean_fraction = geometric mean of per-cell fractions within a metacell
#' Both are accepted as OptionalInput; the R layer uses whichever is present.
#'
#' Gene module note: Metacells.jl has a separate "module" axis for gene modules.
#' The ("gene", "module") vector assigns genes to modules. We add
#' ("gene", "primary_module") as the canonical name to avoid collision with the
#' Metacells.jl module axis. ("gene", "module") is kept as a deprecated alias.
#'
#' @return List containing axes, vectors, matrices, and scalars specifications
#' @export
mcview_core_contract <- function() {
    list(
        name = "MCView Core",
        description = "What MCView expects from the input DAF",
        axes = list(
            # --- Required axes ---
            metacell = axis_spec(CONTRACT_REQUIRED, "Metacell identifiers"),
            gene = axis_spec(CONTRACT_REQUIRED, "Gene identifiers"),
            type = axis_spec(CONTRACT_REQUIRED, "Cell type identifiers"),

            # --- Optional axes ---
            cell = axis_spec(CONTRACT_OPTIONAL,
                "Cell identifiers (enables Samples tab and per-cell metadata)"),

            # Marker correlations: kept as axis (pragmatic compromise).
            # The data is inherently 3D (gene x gene x type) which does not map
            # cleanly to DAF's 2D matrix model. Flattening into an axis with
            # vectors (gene1, gene2, cor, type) is the most practical representation.
            # Pipeline CAN pre-compute these; MCView's derived layer computes
            # them as a fallback.
            mcview_marker_correlations = axis_spec(CONTRACT_OPTIONAL,
                "Pre-computed marker gene correlation pairs (preferred axis name)"),

            # DEPRECATED: gg_mc_top_cor -- kept only for backward compatibility
            # with older pipeline outputs. New pipelines should use
            # mcview_marker_correlations instead.
            gg_mc_top_cor = axis_spec(CONTRACT_OPTIONAL,
                "DEPRECATED: Pre-computed marker gene correlation pairs (legacy axis name)")
        ),
        vectors = list(
            # ==================================================================
            # Required core data
            # ==================================================================
            vector_spec(
                "metacell", "total_UMIs", CONTRACT_REQUIRED, "numeric",
                "Total UMIs per metacell"
            ),
            vector_spec(
                "metacell", "type", CONTRACT_REQUIRED, "character",
                "Cell type assignment per metacell"
            ),
            vector_spec(
                "type", "color", CONTRACT_REQUIRED, "character",
                "Hex color per cell type (#RRGGBB format)"
            ),

            # ==================================================================
            # 2D coordinates (at least one pair required -- see custom rule)
            # ==================================================================
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
                "U coordinate (alternative 2D projection)"
            ),
            vector_spec(
                "metacell", "v", CONTRACT_OPTIONAL, "numeric",
                "V coordinate (alternative 2D projection)"
            ),
            vector_spec(
                "metacell", "w", CONTRACT_OPTIONAL, "numeric",
                "W coordinate (3rd dimension, rarely used)"
            ),
            # Aliases matching Metacells.jl naming convention
            vector_spec(
                "metacell", "umap_x", CONTRACT_OPTIONAL, "numeric",
                "Alias for x (Metacells.jl canonical name)"
            ),
            vector_spec(
                "metacell", "umap_y", CONTRACT_OPTIONAL, "numeric",
                "Alias for y (Metacells.jl canonical name)"
            ),
            vector_spec(
                "metacell", "umap_u", CONTRACT_OPTIONAL, "numeric",
                "Alias for u (Metacells.jl canonical name)"
            ),
            vector_spec(
                "metacell", "umap_v", CONTRACT_OPTIONAL, "numeric",
                "Alias for v (Metacells.jl canonical name)"
            ),
            vector_spec(
                "metacell", "umap_w", CONTRACT_OPTIONAL, "numeric",
                "Alias for w (Metacells.jl canonical name)"
            ),

            # ==================================================================
            # Optional metacell vectors
            # ==================================================================
            # n_cells: Metacells.jl canonical name (preferred)
            vector_spec(
                "metacell", "n_cells", CONTRACT_OPTIONAL, "numeric",
                "Number of cells per metacell (Metacells.jl canonical: n_cells)"
            ),
            # n_cell: backward-compatible alias
            vector_spec(
                "metacell", "n_cell", CONTRACT_OPTIONAL, "numeric",
                "ALIAS for n_cells (legacy MCView name, kept for backward compat)"
            ),
            # NOTE: ("metacell", "cells") removed -- was never in Metacells.jl

            vector_spec(
                "metacell", "is_rare", CONTRACT_OPTIONAL, "logical",
                "Whether metacell is rare"
            ),

            # ==================================================================
            # Optional gene vectors
            # ==================================================================
            vector_spec(
                "gene", "is_marker", CONTRACT_OPTIONAL, "logical",
                "Whether gene is a marker (strongly recommended -- enables Markers tab)"
            ),
            vector_spec(
                "gene", "is_lateral", CONTRACT_OPTIONAL, "logical",
                "Whether gene is lateral (excluded from analysis; colors marker heatmap)"
            ),
            vector_spec(
                "gene", "is_noisy", CONTRACT_OPTIONAL, "logical",
                "Whether gene is noisy (colors marker heatmap)"
            ),

            # Gene module: primary_module is canonical, module is deprecated alias
            # (collision risk with Metacells.jl "module" axis)
            vector_spec(
                "gene", "primary_module", CONTRACT_OPTIONAL, "character",
                "Gene module assignment (enables Gene Modules tab; canonical name)"
            ),
            # DEPRECATED: ("gene", "module") -- kept for backward compat with
            # existing DAF files. New pipelines should use primary_module.
            vector_spec(
                "gene", "module", CONTRACT_OPTIONAL, "character",
                "DEPRECATED alias for primary_module (collision with Metacells.jl module axis)"
            ),

            # Gene boolean flags from cells_clean (informational)
            vector_spec(
                "gene", "is_forbidden", CONTRACT_OPTIONAL, "logical",
                "Whether gene is forbidden"
            ),
            vector_spec(
                "gene", "is_selected", CONTRACT_OPTIONAL, "logical",
                "Whether gene is selected for metacell computation"
            ),
            vector_spec(
                "gene", "is_transcription_factor", CONTRACT_OPTIONAL, "logical",
                "Whether gene is a transcription factor"
            ),
            vector_spec(
                "gene", "is_regulator", CONTRACT_OPTIONAL, "logical",
                "Whether gene is a regulator"
            ),

            # NEW: marker_rank from Metacells.jl
            vector_spec(
                "gene", "marker_rank", CONTRACT_OPTIONAL, "integer",
                "Marker gene ranking from Metacells.jl (1=top marker, typemax=non-marker)"
            ),

            # NEW: significant_inner_folds_count from Metacells.jl anndata import
            vector_spec(
                "gene", "significant_inner_folds_count", CONTRACT_OPTIONAL, "integer",
                "Count of significant inner folds per gene"
            ),

            # ==================================================================
            # Optional per-type vectors
            # ==================================================================
            # NEW: type display order
            vector_spec(
                "type", "order", CONTRACT_OPTIONAL, "integer",
                "Display order of cell types in plots (1=first)"
            ),

            # ==================================================================
            # Cell-level data (enables Samples tab)
            # ==================================================================
            vector_spec(
                "cell", "metacell", CONTRACT_OPTIONAL, "character",
                "Metacell assignment per cell"
            ),
            vector_spec(
                "cell", "type", CONTRACT_OPTIONAL, "character",
                "Cell type per cell"
            ),

            # ==================================================================
            # Projection / Atlas vectors (all optional)
            # ==================================================================
            vector_spec(
                "metacell", "projected_type", CONTRACT_OPTIONAL, "character",
                "Projected type from atlas"
            ),
            vector_spec(
                "metacell", "projected_correlation", CONTRACT_OPTIONAL, "numeric",
                "Correlation with atlas projection"
            ),
            # similar / is_similar: both accepted
            vector_spec(
                "metacell", "similar", CONTRACT_OPTIONAL, "logical",
                "Whether similar to atlas"
            ),
            vector_spec(
                "metacell", "is_similar", CONTRACT_OPTIONAL, "logical",
                "Alias for similar (alternative naming convention)"
            ),
            vector_spec(
                "metacell", "atlas_metacell", CONTRACT_OPTIONAL, "character",
                "Most similar atlas metacell"
            ),

            # ==================================================================
            # Marker correlation vectors (pipeline CAN pre-compute)
            # ==================================================================
            vector_spec(
                "mcview_marker_correlations", "gene1", CONTRACT_OPTIONAL, "character",
                "First gene in correlation pair"
            ),
            vector_spec(
                "mcview_marker_correlations", "gene2", CONTRACT_OPTIONAL, "character",
                "Second gene in correlation pair"
            ),
            vector_spec(
                "mcview_marker_correlations", "cor", CONTRACT_OPTIONAL, "numeric",
                "Correlation value"
            ),
            vector_spec(
                "mcview_marker_correlations", "type", CONTRACT_OPTIONAL, "character",
                "Cell type context for correlation"
            ),

            # DEPRECATED: Legacy axis name for same data
            vector_spec(
                "gg_mc_top_cor", "gene1", CONTRACT_OPTIONAL, "character",
                "DEPRECATED: First gene in correlation pair (legacy axis)"
            ),
            vector_spec(
                "gg_mc_top_cor", "gene2", CONTRACT_OPTIONAL, "character",
                "DEPRECATED: Second gene in correlation pair (legacy axis)"
            ),
            vector_spec(
                "gg_mc_top_cor", "cor", CONTRACT_OPTIONAL, "numeric",
                "DEPRECATED: Correlation value (legacy axis)"
            ),
            vector_spec(
                "gg_mc_top_cor", "type", CONTRACT_OPTIONAL, "character",
                "DEPRECATED: Cell type for correlation (legacy axis)"
            )
        ),
        matrices = list(
            # ==================================================================
            # Required core matrix
            # ==================================================================
            matrix_spec(
                "metacell", "gene", "UMIs", CONTRACT_REQUIRED, "numeric",
                "UMI counts per metacell-gene pair"
            ),

            # ==================================================================
            # Optional gene x metacell matrices
            # ==================================================================
            matrix_spec(
                "gene", "metacell", "inner_fold", CONTRACT_OPTIONAL, "numeric",
                "Inner fold per gene-metacell (enables Inner-fold / QC tabs)"
            ),

            # inner_std_log: Metacells.jl canonical name (preferred)
            matrix_spec(
                "gene", "metacell", "inner_std_log", CONTRACT_OPTIONAL, "numeric",
                "Inner std log per gene-metacell (Metacells.jl canonical name)"
            ),
            # inner_stdev_log: backward-compatible alias
            matrix_spec(
                "gene", "metacell", "inner_stdev_log", CONTRACT_OPTIONAL, "numeric",
                "ALIAS for inner_std_log (legacy name, kept for backward compat)"
            ),

            # geomean_fraction: MCView's own computation (geometric mean of per-cell fractions)
            matrix_spec(
                "gene", "metacell", "geomean_fraction", CONTRACT_OPTIONAL, "numeric",
                "Geometric mean fraction per gene-metacell (MCView computation)"
            ),
            # linear_fraction: Metacells.jl canonical (UMIs / total_UMIs, different from geomean)
            matrix_spec(
                "gene", "metacell", "linear_fraction", CONTRACT_OPTIONAL, "numeric",
                "Linear fraction per gene-metacell (Metacells.jl canonical; UMIs/total_UMIs)"
            ),

            matrix_spec(
                "gene", "metacell", "zeros", CONTRACT_OPTIONAL, "numeric",
                "Zero counts per gene-metacell (enables Zero-fold tab)"
            ),

            # ==================================================================
            # MC x MC matrices -- REMOVED
            # ==================================================================
            # NOTE: The following MC x MC matrices were removed because no R code
            # in MCView reads them:
            #   outgoing_weights, obs_balanced_ranks, umap_distances

            # ==================================================================
            # Projection / Atlas matrices (all optional)
            # ==================================================================
            matrix_spec(
                "gene", "metacell", "projected_fold", CONTRACT_OPTIONAL, "numeric",
                "Projected fold from atlas"
            ),

            # projected_fraction / projected_geomean_fraction: both accepted
            matrix_spec(
                "metacell", "gene", "projected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Projected expression fraction from atlas"
            ),
            matrix_spec(
                "metacell", "gene", "projected_geomean_fraction", CONTRACT_OPTIONAL, "numeric",
                "Alias for projected_fraction"
            ),

            # corrected_fraction / corrected_geomean_fraction: both accepted
            matrix_spec(
                "metacell", "gene", "corrected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Corrected expression fraction from atlas projection"
            ),
            matrix_spec(
                "metacell", "gene", "corrected_geomean_fraction", CONTRACT_OPTIONAL, "numeric",
                "Alias for corrected_fraction"
            )
        ),
        scalars = list(
            # ==================================================================
            # Configuration scalars (all optional -- can come from YAML instead)
            # ==================================================================
            scalar_spec("mcview_title", CONTRACT_OPTIONAL, "character",
                        "Application title displayed in UI header"),
            scalar_spec("mcview_tabs", CONTRACT_OPTIONAL, "character",
                        "Comma-separated list of enabled tab names"),
            scalar_spec("mcview_excluded_tabs", CONTRACT_OPTIONAL, "character",
                        "Comma-separated list of excluded tab names"),
            scalar_spec("mcview_light_version", CONTRACT_OPTIONAL, "logical",
                        "Whether to use the light version of the UI"),
            scalar_spec("mcview_about_markdown", CONTRACT_OPTIONAL, "character",
                        "Markdown content for the About tab"),
            scalar_spec("mcview_sample_property", CONTRACT_OPTIONAL, "character",
                        "Name of the cell vector to use as sample ID (e.g. embryo, batch_set_id, set)"),
            scalar_spec("mcview_available_tabs", CONTRACT_OPTIONAL, "character",
                        "Pre-stored comma-separated list of available tabs (skips runtime detection)"),
            scalar_spec("metacells_algorithm", CONTRACT_OPTIONAL, "character",
                        "Metacells algorithm version (e.g. metacells2)"),

            # Projection configuration
            scalar_spec("project_max_projection_fold_factor", CONTRACT_OPTIONAL, "numeric",
                        "Maximum projection fold factor for atlas tab"),
            scalar_spec("projection_weights_json", CONTRACT_OPTIONAL, "character",
                        "Projection weights as JSON (query-atlas metacell mapping)"),
            scalar_spec("query_cell_type_fracs_json", CONTRACT_OPTIONAL, "character",
                        "Query cell type fractions as JSON"),

            # QC statistics -- canonical names (no qc_stats_ prefix)
            scalar_spec("n_outliers", CONTRACT_OPTIONAL, "numeric",
                        "Number of outlier cells"),
            scalar_spec("n_cells_total", CONTRACT_OPTIONAL, "numeric",
                        "Total number of cells (renamed from qc_stats_n_cells to avoid collision with vector name n_cells)"),
            scalar_spec("n_umis", CONTRACT_OPTIONAL, "numeric",
                        "Total UMI count"),
            scalar_spec("median_umis_per_metacell", CONTRACT_OPTIONAL, "numeric",
                        "Median UMIs per metacell"),
            scalar_spec("median_cells_per_metacell", CONTRACT_OPTIONAL, "numeric",
                        "Median cells per metacell"),

            # QC statistics -- backward-compatible aliases (old qc_stats_ prefix)
            scalar_spec("qc_stats_n_outliers", CONTRACT_OPTIONAL, "numeric",
                        "ALIAS for n_outliers (legacy name)"),
            scalar_spec("qc_stats_n_cells", CONTRACT_OPTIONAL, "numeric",
                        "ALIAS for n_cells_total (legacy name)"),
            scalar_spec("qc_stats_n_umis", CONTRACT_OPTIONAL, "numeric",
                        "ALIAS for n_umis (legacy name)"),
            scalar_spec("qc_stats_median_umis_per_metacell", CONTRACT_OPTIONAL, "numeric",
                        "ALIAS for median_umis_per_metacell (legacy name)"),
            scalar_spec("qc_stats_median_cells_per_metacell", CONTRACT_OPTIONAL, "numeric",
                        "ALIAS for median_cells_per_metacell (legacy name)")
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
#'
#' 2D projection with cell type coloring.
#' - Requires: core + at least one coordinate pair (x,y) or (u,v)
#' - Enriched by: w coordinate, n_cells
#' - Accepts umap_x/umap_y/umap_u/umap_v/umap_w aliases
#' - NOTE: MC x MC outgoing_weights removed -- not read by R code
#'
#' @export
mcview_manifold_contract <- function() {
    list(
        name = "Manifold",
        description = "2D projection and cell type visualization",
        extends = "core",
        vectors = list(
            # Coordinate aliases for Metacells.jl naming
            vector_spec(
                "metacell", "umap_x", CONTRACT_OPTIONAL, "numeric",
                "Alias for x (Metacells.jl canonical name)"
            ),
            vector_spec(
                "metacell", "umap_y", CONTRACT_OPTIONAL, "numeric",
                "Alias for y (Metacells.jl canonical name)"
            ),
            vector_spec(
                "metacell", "umap_u", CONTRACT_OPTIONAL, "numeric",
                "Alias for u (Metacells.jl canonical name)"
            ),
            vector_spec(
                "metacell", "umap_v", CONTRACT_OPTIONAL, "numeric",
                "Alias for v (Metacells.jl canonical name)"
            ),
            vector_spec(
                "metacell", "umap_w", CONTRACT_OPTIONAL, "numeric",
                "Alias for w (Metacells.jl canonical name)"
            ),
            # n_cells / n_cell
            vector_spec(
                "metacell", "n_cells", CONTRACT_OPTIONAL, "numeric",
                "Cells per metacell (Metacells.jl canonical)"
            ),
            vector_spec(
                "metacell", "n_cell", CONTRACT_OPTIONAL, "numeric",
                "Alias for n_cells (legacy)"
            )
        )
    )
}

#' Contract for Genes tab
#'
#' Per-gene expression across metacells.
#' - Requires: core
#' - Enriched by: top-2 genes [derived], is_lateral/is_noisy [input]
#'
#' @export
mcview_genes_contract <- function() {
    list(
        name = "Genes",
        description = "Gene expression visualization",
        extends = "core",
        vectors = list(
            vector_spec(
                "gene", "is_lateral", CONTRACT_OPTIONAL, "logical",
                "Lateral gene [input]"
            ),
            vector_spec(
                "gene", "is_noisy", CONTRACT_OPTIONAL, "logical",
                "Noisy gene [input]"
            ),
            # Top genes come from derived layer
            vector_spec(
                "metacell", "mcview_cache_top1_gene", CONTRACT_OPTIONAL, "character",
                "Top gene [derived]"
            ),
            vector_spec(
                "metacell", "mcview_cache_top2_gene", CONTRACT_OPTIONAL, "character",
                "Second gene [derived]"
            ),
            vector_spec(
                "metacell", "mcview_cache_top1_lfp", CONTRACT_OPTIONAL, "numeric",
                "Top gene LFP [derived]"
            ),
            vector_spec(
                "metacell", "mcview_cache_top2_lfp", CONTRACT_OPTIONAL, "numeric",
                "Second gene LFP [derived]"
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
                "metacell", "n_cells", CONTRACT_OPTIONAL, "numeric",
                "Number of cells per metacell (Metacells.jl canonical)"
            ),
            vector_spec(
                "metacell", "n_cell", CONTRACT_OPTIONAL, "numeric",
                "Alias for n_cells (legacy)"
            )
        )
    )
}

#' Contract for QC tab
#'
#' Quality control metrics.
#' - Requires: core
#' - Enriched by: n_cells, inner_fold, QC scalars [input]
#' - QC scalars accept both new names (no prefix) and old qc_stats_ prefix.
#'
#' @export
mcview_qc_contract <- function() {
    list(
        name = "QC",
        description = "Quality control metrics",
        extends = "core",
        vectors = list(
            vector_spec(
                "metacell", "n_cells", CONTRACT_OPTIONAL, "numeric",
                "Cells per metacell (Metacells.jl canonical)"
            ),
            vector_spec(
                "metacell", "n_cell", CONTRACT_OPTIONAL, "numeric",
                "Alias for n_cells (legacy)"
            ),
            vector_spec(
                "metacell", "is_rare", CONTRACT_OPTIONAL, "logical",
                "Rare metacell flag"
            )
        ),
        matrices = list(
            matrix_spec(
                "gene", "metacell", "inner_fold", CONTRACT_OPTIONAL, "numeric",
                "Inner fold"
            )
        ),
        scalars = list(
            # Canonical QC scalar names
            scalar_spec(
                "n_outliers", CONTRACT_OPTIONAL, "numeric",
                "Outlier count"
            ),
            scalar_spec(
                "n_cells_total", CONTRACT_OPTIONAL, "numeric",
                "Total cells"
            ),
            scalar_spec(
                "n_umis", CONTRACT_OPTIONAL, "numeric",
                "Total UMIs"
            ),
            scalar_spec(
                "median_umis_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "Median UMIs per metacell"
            ),
            scalar_spec(
                "median_cells_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "Median cells per metacell"
            ),
            # Legacy aliases
            scalar_spec(
                "qc_stats_n_outliers", CONTRACT_OPTIONAL, "numeric",
                "ALIAS: Outlier count"
            ),
            scalar_spec(
                "qc_stats_n_cells", CONTRACT_OPTIONAL, "numeric",
                "ALIAS: Total cells"
            ),
            scalar_spec(
                "qc_stats_n_umis", CONTRACT_OPTIONAL, "numeric",
                "ALIAS: Total UMIs"
            ),
            scalar_spec(
                "qc_stats_median_umis_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "ALIAS: Median UMIs/mc"
            ),
            scalar_spec(
                "qc_stats_median_cells_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "ALIAS: Median cells/mc"
            )
        )
    )
}

#' Contract for Markers tab
#'
#' Marker gene heatmap.
#' - Requires: core + is_marker [input]
#' - Enriched by: is_lateral, is_noisy [input], marker_rank [input]
#'
#' @export
mcview_markers_contract <- function() {
    list(
        name = "Markers",
        description = "Marker gene visualization",
        extends = "core",
        vectors = list(
            vector_spec(
                "gene", "is_marker", CONTRACT_REQUIRED, "logical",
                "Marker gene flag [input] -- required for this tab"
            ),
            vector_spec(
                "gene", "is_lateral", CONTRACT_OPTIONAL, "logical",
                "Lateral flag [input]"
            ),
            vector_spec(
                "gene", "is_noisy", CONTRACT_OPTIONAL, "logical",
                "Noisy flag [input]"
            ),
            vector_spec(
                "gene", "marker_rank", CONTRACT_OPTIONAL, "integer",
                "Marker rank from Metacells.jl [input]"
            )
        )
    )
}

#' Contract for Gene Modules tab
#'
#' Gene module heatmap.
#' - Requires: core + gene module assignment [input]
#' - Accepts both primary_module (canonical) and module (deprecated)
#'
#' @export
mcview_gene_modules_contract <- function() {
    list(
        name = "Gene modules",
        description = "Gene module analysis",
        extends = "core",
        vectors = list(
            # Accept both names; R layer checks which is present
            vector_spec(
                "gene", "primary_module", CONTRACT_OPTIONAL, "character",
                "Gene module assignment [input] (canonical name)"
            ),
            # DEPRECATED alias
            vector_spec(
                "gene", "module", CONTRACT_OPTIONAL, "character",
                "DEPRECATED alias for primary_module [input]"
            )
        )
    )
}

#' Contract for Projection QC tab
#'
#' Atlas projection quality control.
#' - Requires: core + projection vectors [input]
#' - Accepts both similar and is_similar
#' - Accepts both projected_fraction and projected_geomean_fraction
#'
#' @export
mcview_projection_qc_contract <- function() {
    list(
        name = "Projection QC",
        description = "Atlas projection quality control",
        extends = "core",
        vectors = list(
            vector_spec(
                "metacell", "projected_type", CONTRACT_OPTIONAL, "character",
                "Projected type"
            ),
            vector_spec(
                "metacell", "projected_correlation", CONTRACT_OPTIONAL, "numeric",
                "Projection corr"
            ),
            vector_spec(
                "metacell", "similar", CONTRACT_OPTIONAL, "logical",
                "Similar to atlas"
            ),
            vector_spec(
                "metacell", "is_similar", CONTRACT_OPTIONAL, "logical",
                "Alias for similar"
            ),
            vector_spec(
                "metacell", "atlas_metacell", CONTRACT_OPTIONAL, "character",
                "Atlas MC"
            )
        ),
        matrices = list(
            matrix_spec(
                "gene", "metacell", "projected_fold", CONTRACT_OPTIONAL, "numeric",
                "Projected fold"
            ),
            matrix_spec(
                "metacell", "gene", "projected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Projected frac"
            ),
            matrix_spec(
                "metacell", "gene", "projected_geomean_fraction", CONTRACT_OPTIONAL, "numeric",
                "Alias for projected_fraction"
            ),
            matrix_spec(
                "metacell", "gene", "corrected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Corrected frac"
            ),
            matrix_spec(
                "metacell", "gene", "corrected_geomean_fraction", CONTRACT_OPTIONAL, "numeric",
                "Alias for corrected_fraction"
            )
        ),
        scalars = list(
            scalar_spec(
                "project_max_projection_fold_factor", CONTRACT_OPTIONAL, "numeric",
                "Max fold factor"
            ),
            scalar_spec(
                "projection_weights_json", CONTRACT_OPTIONAL, "character",
                "Weights JSON"
            ),
            scalar_spec(
                "query_cell_type_fracs_json", CONTRACT_OPTIONAL, "character",
                "Type fracs JSON"
            )
        )
    )
}

#' Contract for Atlas tab
#'
#' Atlas comparison visualization.
#' - Requires: core + projection metadata [input]
#' - Accepts both similar and is_similar
#' - Accepts both projected_fraction and projected_geomean_fraction
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
                "Projected type"
            ),
            vector_spec(
                "metacell", "projected_correlation", CONTRACT_OPTIONAL, "numeric",
                "Projection corr"
            ),
            vector_spec(
                "metacell", "similar", CONTRACT_OPTIONAL, "logical",
                "Similar to atlas"
            ),
            vector_spec(
                "metacell", "is_similar", CONTRACT_OPTIONAL, "logical",
                "Alias for similar"
            )
        ),
        matrices = list(
            matrix_spec(
                "metacell", "gene", "corrected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Corrected frac"
            ),
            matrix_spec(
                "metacell", "gene", "corrected_geomean_fraction", CONTRACT_OPTIONAL, "numeric",
                "Alias for corrected_fraction"
            ),
            matrix_spec(
                "metacell", "gene", "projected_fraction", CONTRACT_OPTIONAL, "numeric",
                "Projected frac"
            ),
            matrix_spec(
                "metacell", "gene", "projected_geomean_fraction", CONTRACT_OPTIONAL, "numeric",
                "Alias for projected_fraction"
            )
        ),
        scalars = list(
            scalar_spec(
                "projection_weights_json", CONTRACT_OPTIONAL, "character",
                "Weights JSON"
            ),
            scalar_spec(
                "query_cell_type_fracs_json", CONTRACT_OPTIONAL, "character",
                "Type fracs JSON"
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
#'
#' Gene-gene correlations per type.
#' - Requires: core + marker correlations [either input or derived]
#' - gg_mc_top_cor is DEPRECATED but still accepted
#'
#' Accepts either the preferred `mcview_marker_correlations` axis name or the
#' legacy `gg_mc_top_cor` axis. MCView's derived layer computes as a fallback
#' if neither is present in the input DAF.
#'
#' @export
mcview_gene_correlations_contract <- function() {
    list(
        name = "Gene Correlations",
        description = "Gene-gene correlation data (from input or derived layer)",
        axes = list(
            mcview_marker_correlations = axis_spec(CONTRACT_OPTIONAL,
                "Correlation pairs [either]"),
            # DEPRECATED
            gg_mc_top_cor = axis_spec(CONTRACT_OPTIONAL,
                "DEPRECATED: Correlation pairs legacy axis [either]")
        ),
        vectors = list(
            # Preferred axis name
            vector_spec(
                "mcview_marker_correlations", "gene1", CONTRACT_OPTIONAL, "character",
                "First gene in correlation pair"
            ),
            vector_spec(
                "mcview_marker_correlations", "gene2", CONTRACT_OPTIONAL, "character",
                "Second gene in correlation pair"
            ),
            vector_spec(
                "mcview_marker_correlations", "cor", CONTRACT_OPTIONAL, "numeric",
                "Correlation value"
            ),
            vector_spec(
                "mcview_marker_correlations", "type", CONTRACT_OPTIONAL, "character",
                "Cell type context for correlation"
            ),

            # DEPRECATED legacy axis vectors
            vector_spec(
                "gg_mc_top_cor", "gene1", CONTRACT_OPTIONAL, "character",
                "DEPRECATED: Gene 1 (legacy)"
            ),
            vector_spec(
                "gg_mc_top_cor", "gene2", CONTRACT_OPTIONAL, "character",
                "DEPRECATED: Gene 2 (legacy)"
            ),
            vector_spec(
                "gg_mc_top_cor", "cor", CONTRACT_OPTIONAL, "numeric",
                "DEPRECATED: Correlation (legacy)"
            ),
            vector_spec(
                "gg_mc_top_cor", "type", CONTRACT_OPTIONAL, "character",
                "DEPRECATED: Cell type (legacy)"
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
                "Application title displayed in UI header"
            ),
            scalar_spec(
                "mcview_tabs", CONTRACT_OPTIONAL, "character",
                "Comma-separated list of enabled tab names"
            ),
            scalar_spec(
                "mcview_excluded_tabs", CONTRACT_OPTIONAL, "character",
                "Comma-separated list of excluded tab names"
            ),
            scalar_spec(
                "mcview_light_version", CONTRACT_OPTIONAL, "logical",
                "Whether to use the light version of the UI"
            ),
            scalar_spec(
                "mcview_about_markdown", CONTRACT_OPTIONAL, "character",
                "Markdown content for the About tab"
            ),
            scalar_spec(
                "mcview_sample_property", CONTRACT_OPTIONAL, "character",
                "Name of the cell vector to use as sample ID (e.g. embryo, batch_set_id, set)"
            ),
            scalar_spec(
                "mcview_available_tabs", CONTRACT_OPTIONAL, "character",
                "Pre-stored comma-separated list of available tabs (skips runtime detection)"
            ),
            scalar_spec(
                "metacells_algorithm", CONTRACT_OPTIONAL, "character",
                "Metacells algorithm version (e.g. metacells2)"
            ),

            # Projection configuration
            scalar_spec(
                "project_max_projection_fold_factor", CONTRACT_OPTIONAL, "numeric",
                "Maximum projection fold factor for atlas tab"
            ),
            scalar_spec(
                "projection_weights_json", CONTRACT_OPTIONAL, "character",
                "Projection weights as JSON (query-atlas metacell mapping)"
            ),
            scalar_spec(
                "query_cell_type_fracs_json", CONTRACT_OPTIONAL, "character",
                "Query cell type fractions as JSON"
            ),

            # QC statistics -- canonical names (no qc_stats_ prefix)
            scalar_spec(
                "n_outliers", CONTRACT_OPTIONAL, "numeric",
                "Number of outlier cells"
            ),
            scalar_spec(
                "n_cells_total", CONTRACT_OPTIONAL, "numeric",
                "Total number of cells (renamed from qc_stats_n_cells to avoid collision with vector name n_cells)"
            ),
            scalar_spec(
                "n_umis", CONTRACT_OPTIONAL, "numeric",
                "Total UMI count"
            ),
            scalar_spec(
                "median_umis_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "Median UMIs per metacell"
            ),
            scalar_spec(
                "median_cells_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "Median cells per metacell"
            ),

            # QC statistics -- backward-compatible aliases (old qc_stats_ prefix)
            scalar_spec(
                "qc_stats_n_outliers", CONTRACT_OPTIONAL, "numeric",
                "ALIAS for n_outliers (legacy name)"
            ),
            scalar_spec(
                "qc_stats_n_cells", CONTRACT_OPTIONAL, "numeric",
                "ALIAS for n_cells_total (legacy name)"
            ),
            scalar_spec(
                "qc_stats_n_umis", CONTRACT_OPTIONAL, "numeric",
                "ALIAS for n_umis (legacy name)"
            ),
            scalar_spec(
                "qc_stats_median_umis_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "ALIAS for median_umis_per_metacell (legacy name)"
            ),
            scalar_spec(
                "qc_stats_median_cells_per_metacell", CONTRACT_OPTIONAL, "numeric",
                "ALIAS for median_cells_per_metacell (legacy name)"
            )
        )
    )
}

#' Contract for Samples tab
#'
#' The Samples tab requires cell-level metadata with sample assignments.
#' The sample ID vector is named by the `mcview_sample_property` scalar
#' (e.g. "embryo", "batch_set_id", "set"). There is no hardcoded `samp_id`
#' vector -- the cell vector name is dynamic. The legacy `samp_id` name
#' is accepted as a fallback.
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
                "Metacell assignment per cell"
            ),

            # Cell -> type (optional)
            vector_spec(
                "cell", "type", CONTRACT_OPTIONAL, "character",
                "Cell type per cell"
            )

            # Note: additional cell vectors (embryo, batch_set_id, set, etc.)
            # are dynamic metadata -- no fixed contract entries needed. The
            # mcview_sample_property scalar names which one to use as sample ID.
        ),
        scalars = list(
            scalar_spec(
                "mcview_sample_property", CONTRACT_OPTIONAL, "character",
                "Name of cell vector to use as sample ID"
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

#' Contract for the mcview_derived layer (DERIVED contract)
#'
#' What the mcview_derived writable layer guarantees to provide on top of the
#' input DAF. These are computed during MCView's cache-building phase.
#'
#' The derived layer is a writable DAF that chains on top of the read-only
#' input DAF. Any data already present in the input (e.g. marker correlations
#' pre-computed by the pipeline) will NOT be overwritten -- the derived layer
#' only fills in what is missing.
#'
#' Type markers note: Per-type marker data is stored as two sparse matrices on
#' (type, gene) axes rather than as an axis with vectors. This is the correct
#' DAF data model: the data is naturally indexed by (type, gene) and sparse
#' matrices handle the sparsity (most gene x type pairs are not markers).
#'
#' @export
mcview_cache_contract <- function() {
    list(
        name = "MCView Derived",
        description = paste(
            "What the mcview_derived layer provides on top of the input DAF.",
            "Computed during MCView's cache-building phase.",
            "Input data takes precedence via the chain lookup order."
        ),
        axes = list(
            # Marker correlations (fallback -- only computed if not in input)
            # Kept as axis: pragmatic compromise because the data is inherently
            # 3D (gene x gene x type). Flattening into an axis with vectors
            # (gene1, gene2, cor, type) is the most practical DAF representation.
            mcview_marker_correlations = axis_spec(
                CONTRACT_OPTIONAL,
                "Marker gene correlation pairs (fallback if not in input DAF)"
            )
        ),
        vectors = list(
            # ==================================================================
            # Top-2 genes per metacell (hover info, gene coloring)
            # ==================================================================
            vector_spec(
                "metacell", "mcview_cache_top1_gene", CONTRACT_OPTIONAL, "character",
                "Top expressed gene per metacell (by log fold-enrichment)"
            ),
            vector_spec(
                "metacell", "mcview_cache_top2_gene", CONTRACT_OPTIONAL, "character",
                "Second top expressed gene per metacell"
            ),
            vector_spec(
                "metacell", "mcview_cache_top1_lfp", CONTRACT_OPTIONAL, "numeric",
                "Log fold-enrichment of top gene per metacell"
            ),
            vector_spec(
                "metacell", "mcview_cache_top2_lfp", CONTRACT_OPTIONAL, "numeric",
                "Log fold-enrichment of second gene per metacell"
            ),

            # ==================================================================
            # Gene statistics (for gene selectors and filtering)
            # ==================================================================
            vector_spec(
                "gene", "mcview_cache_gene_max_umis", CONTRACT_OPTIONAL, "numeric",
                "Maximum UMIs across metacells per gene"
            ),
            vector_spec(
                "gene", "mcview_cache_gene_mean_umis", CONTRACT_OPTIONAL, "numeric",
                "Mean UMIs across metacells per gene"
            ),
            vector_spec(
                "gene", "mcview_cache_gene_sum_umis", CONTRACT_OPTIONAL, "numeric",
                "Sum of UMIs across metacells per gene"
            ),

            # ==================================================================
            # Marker correlations (fallback computation)
            # If the input DAF already provides mcview_marker_correlations or
            # gg_mc_top_cor, the derived layer does NOT recompute -- the input
            # data takes precedence via the chain lookup order.
            # ==================================================================
            vector_spec(
                "mcview_marker_correlations", "gene1", CONTRACT_OPTIONAL, "character",
                "First gene in correlation pair"
            ),
            vector_spec(
                "mcview_marker_correlations", "gene2", CONTRACT_OPTIONAL, "character",
                "Second gene in correlation pair"
            ),
            vector_spec(
                "mcview_marker_correlations", "cor", CONTRACT_OPTIONAL, "numeric",
                "Correlation value"
            ),
            vector_spec(
                "mcview_marker_correlations", "type", CONTRACT_OPTIONAL, "character",
                "Cell type context for correlation"
            ),

            # ==================================================================
            # QC vectors (max inner fold per metacell)
            # ==================================================================
            vector_spec(
                "metacell", "mcview_cache_max_inner_fold", CONTRACT_OPTIONAL, "numeric",
                "Max inner_fold per metacell across all genes"
            ),
            vector_spec(
                "metacell", "mcview_cache_max_inner_fold_no_lateral", CONTRACT_OPTIONAL, "numeric",
                "Max inner_fold per metacell excluding lateral genes"
            ),

            # ==================================================================
            # Per-gene max inner fold (optional derived)
            # ==================================================================
            vector_spec(
                "gene", "mcview_cache_gene_max_inner_fold", CONTRACT_OPTIONAL, "numeric",
                "Max inner_fold per gene across metacells"
            )
        ),
        matrices = list(
            # ==================================================================
            # Per-type marker genes -- sparse matrices on (type, gene)
            # Replaces the old mcview_type_markers axis + vectors approach.
            # Sparse because most (type, gene) pairs are not markers.
            # ==================================================================
            matrix_spec(
                "type", "gene", "mcview_marker_rank", CONTRACT_OPTIONAL, "integer",
                "Per-type marker rank (1=top marker for that type; 0 or missing=not a marker)"
            ),
            matrix_spec(
                "type", "gene", "mcview_marker_fold_change", CONTRACT_OPTIONAL, "numeric",
                "Per-type marker fold change (gene expression in type vs. background)"
            ),

            # ==================================================================
            # Pre-computed EGC (geomean_fraction fallback)
            # ==================================================================
            matrix_spec(
                "gene", "metacell", "mcview_cache_geomean_fraction", CONTRACT_OPTIONAL, "numeric",
                "Pre-computed geomean_fraction when not in input DAF"
            ),

            # ==================================================================
            # Default markers distance matrix
            # ==================================================================
            matrix_spec(
                "metacell", "metacell", "mcview_default_markers_dist", CONTRACT_OPTIONAL, "numeric",
                "Pre-computed hclust distance matrix for default markers"
            )
        ),
        scalars = list(
            # ==================================================================
            # Cache metadata scalars
            # ==================================================================
            scalar_spec(
                "mcview_cache_created", CONTRACT_OPTIONAL, "character",
                "ISO 8601 timestamp when cache was created"
            ),
            scalar_spec(
                "mcview_cache_base_version", CONTRACT_OPTIONAL, "character",
                "MCView version that created this cache"
            ),
            scalar_spec(
                "mcview_cache_base_hash", CONTRACT_OPTIONAL, "character",
                "Hash of the input DAF at cache creation time (for staleness detection)"
            ),

            # ==================================================================
            # Default markers list
            # ==================================================================
            scalar_spec(
                "mcview_default_markers", CONTRACT_OPTIONAL, "character",
                "Comma-separated list of default marker gene names for heatmap ordering"
            ),

            # ==================================================================
            # Cell grouping fields (optional)
            # ==================================================================
            scalar_spec(
                "mcview_cell_grouping_fields", CONTRACT_OPTIONAL, "character",
                "Comma-separated list of categorical cell vectors suitable for grouping in Samples tab"
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

        # Top-2 genes per metacell (lfp is metacell x gene).
        tops <- top2_per_row(lfp)
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
        projection_qc = mcview_projection_qc_contract(),
        atlas = mcview_atlas_contract(),
        annotate = mcview_annotate_contract(),
        samples = mcview_samples_contract(),
        flow = mcview_flow_contract(),
        gene_correlations = mcview_gene_correlations_contract(),
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

    if (!dafr::is_daf(daf_obj)) {
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
