# daf_contracts_definitions.R - MCView contract definitions
#
# Per-tab contract functions, the cache contract, the precompute helper,
# and the all-contracts / by-tab lookup index. Split from R/daf_contracts.R
# (2026-05-01).
#
# Companion files:
#   - R/daf_contracts_primitives.R - axis/vector/matrix/scalar spec helpers
#     and the mcview->dafr adapter.
#   - R/daf_contracts_validators.R - validate_mcview_contract / _tabs.
#   - R/daf_contracts_docs.R - print_contract / contract_to_markdown /
#     generate_contract_docs.

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

    # Compute top genes from lfp
    mc_mat <- dafr::get_matrix(daf_obj, "metacell", "gene", "UMIs")
    mc_sum <- dafr::get_vector(daf_obj, "metacell", "total_UMIs")

    egc <- sweep(mc_mat, 1, mc_sum, "/")
    gene_means <- colMeans(egc)
    gene_means[gene_means == 0] <- 1e-10
    lfp <- dafr::fast_log(as.matrix(sweep(egc, 2, gene_means, "/")),
                          eps = 1e-5, base = 2)

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

    # Compute gene statistics
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
