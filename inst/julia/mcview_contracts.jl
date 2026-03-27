# mcview_contracts.jl - Formal DAF Contracts for MCView
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
# 1. MCVIEW_INPUT_CONTRACT   - what MCView expects from the input DAF
#    (RequiredInput / OptionalInput)
# 2. MCVIEW_DERIVED_CONTRACT - what the mcview_derived layer provides
#    (GuaranteedOutput)
# 3. Tab-specific contracts   - per-tab data requirements (reference only)
# 4. MCVIEW_FULL_CONTRACT     - relaxed union of input + derived
# 5. Validation functions     - validate_mcview_input, validate_mcview_derived,
#    get_available_tabs
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
#
# Usage:
#   using DataAxesFormats.Contracts
#   include("mcview_contracts.jl")
#   verify_input(contractor("MCView.Input", MCVIEW_INPUT_CONTRACT, daf_obj))

using DataAxesFormats.Contracts

# ==============================================================================
# INPUT CONTRACT - What MCView requires from the input DAF
# ==============================================================================

"""
    MCVIEW_INPUT_CONTRACT

What MCView expects to find in the input DAF (read-only, from the metacells
pipeline). Items marked RequiredInput must exist; items marked OptionalInput
enrich specific tabs when present.

Coordinate note: at least one pair of (x,y) or (u,v) is required for the
Manifold tab. The contract marks all four as OptionalInput because the
"at least one pair" constraint is enforced by a custom rule in R.

Sample ID note: `mcview_sample_property` is a scalar that names which cell
vector to use as the sample identifier (e.g. "embryo", "batch_set_id", "set").
There is no hardcoded `samp_id` vector -- the cell vector name is dynamic.

Inner stdev note: the pipeline may produce either `inner_stdev_log` or
`inner_std_log`. MCView accepts both names; the R layer normalizes to
whichever is present. `inner_std_log` is the canonical Metacells.jl name;
`inner_stdev_log` is kept as a backward-compatible alias.

n_cells note: Metacells.jl uses ("metacell", "n_cells"). MCView also accepts
the legacy names ("metacell", "n_cell"). The old ("metacell", "cells") alias
has been removed.

geomean_fraction vs linear_fraction: these are DIFFERENT computations.
  - linear_fraction = UMIs / total_UMIs (Metacells.jl canonical)
  - geomean_fraction = geometric mean of per-cell fractions within a metacell
Both are accepted as OptionalInput; the R layer uses whichever is present.

Gene module note: Metacells.jl has a separate "module" axis for gene modules.
The ("gene", "module") vector assigns genes to modules. We add
("gene", "primary_module") as the canonical name to avoid collision with the
Metacells.jl module axis. ("gene", "module") is kept as a deprecated alias.
"""
MCVIEW_INPUT_CONTRACT = Contract(
    is_relaxed = true,
    axes = [
        # --- Required axes ---
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),

        # --- Optional axes ---
        "cell" => (OptionalInput,
            "Cell identifiers (enables Samples tab and per-cell metadata)"),

        # Marker correlations: kept as axis (pragmatic compromise).
        # The data is inherently 3D (gene x gene x type) which does not map
        # cleanly to DAF's 2D matrix model. Flattening into an axis with
        # vectors (gene1, gene2, cor, type) is the most practical representation.
        # Pipeline CAN pre-compute these; MCView's derived layer computes
        # them as a fallback.
        "mcview_marker_correlations" => (OptionalInput,
            "Pre-computed marker gene correlation pairs (preferred axis name)"),

        # DEPRECATED: gg_mc_top_cor -- kept only for backward compatibility
        # with older pipeline outputs. New pipelines should use
        # mcview_marker_correlations instead.
        "gg_mc_top_cor" => (OptionalInput,
            "DEPRECATED: Pre-computed marker gene correlation pairs (legacy axis name)"),
    ],
    data = [
        # ==================================================================
        # Required core data
        # ==================================================================

        # UMI count matrix (dense UInt32 or numeric)
        ("metacell", "gene", "UMIs") => (RequiredInput, StorageFloat,
            "UMI counts per metacell-gene pair"),

        # Per-metacell vectors
        ("metacell", "total_UMIs") => (RequiredInput, StorageFloat,
            "Total UMIs per metacell"),
        ("metacell", "type") => (RequiredInput, AbstractString,
            "Cell type assignment per metacell"),

        # Per-type vectors
        ("type", "color") => (RequiredInput, AbstractString,
            "Hex color per cell type (#RRGGBB format)"),

        # ==================================================================
        # 2D coordinates (at least one pair required -- see docstring)
        # Metacells.jl uses umap_x, umap_y, umap_u, umap_v, umap_w.
        # MCView also accepts the short names x, y, u, v, w.
        # ==================================================================
        ("metacell", "x") => (OptionalInput, StorageFloat,
            "X coordinate for 2D projection"),
        ("metacell", "y") => (OptionalInput, StorageFloat,
            "Y coordinate for 2D projection"),
        ("metacell", "u") => (OptionalInput, StorageFloat,
            "U coordinate (alternative 2D projection)"),
        ("metacell", "v") => (OptionalInput, StorageFloat,
            "V coordinate (alternative 2D projection)"),
        ("metacell", "w") => (OptionalInput, StorageFloat,
            "W coordinate (3rd dimension, rarely used)"),
        # Aliases matching Metacells.jl naming convention
        ("metacell", "umap_x") => (OptionalInput, StorageFloat,
            "Alias for x (Metacells.jl canonical name)"),
        ("metacell", "umap_y") => (OptionalInput, StorageFloat,
            "Alias for y (Metacells.jl canonical name)"),
        ("metacell", "umap_u") => (OptionalInput, StorageFloat,
            "Alias for u (Metacells.jl canonical name)"),
        ("metacell", "umap_v") => (OptionalInput, StorageFloat,
            "Alias for v (Metacells.jl canonical name)"),
        ("metacell", "umap_w") => (OptionalInput, StorageFloat,
            "Alias for w (Metacells.jl canonical name)"),

        # ==================================================================
        # Optional metacell vectors
        # ==================================================================
        # n_cells: Metacells.jl canonical name (preferred)
        ("metacell", "n_cells") => (OptionalInput, StorageFloat,
            "Number of cells per metacell (Metacells.jl canonical: n_cells)"),
        # n_cell: backward-compatible alias
        ("metacell", "n_cell") => (OptionalInput, StorageFloat,
            "ALIAS for n_cells (legacy MCView name, kept for backward compat)"),
        # NOTE: ("metacell", "cells") removed -- was never in Metacells.jl

        ("metacell", "is_rare") => (OptionalInput, Bool,
            "Whether metacell is rare"),

        # ==================================================================
        # Optional gene vectors
        # ==================================================================
        ("gene", "is_marker") => (OptionalInput, Bool,
            "Whether gene is a marker (strongly recommended -- enables Markers tab)"),
        ("gene", "is_lateral") => (OptionalInput, Bool,
            "Whether gene is lateral (excluded from analysis; colors marker heatmap)"),
        ("gene", "is_noisy") => (OptionalInput, Bool,
            "Whether gene is noisy (colors marker heatmap)"),

        # Gene module: primary_module is canonical, module is deprecated alias
        # (collision risk with Metacells.jl "module" axis)
        ("gene", "primary_module") => (OptionalInput, AbstractString,
            "Gene module assignment (enables Gene Modules tab; canonical name)"),
        # DEPRECATED: ("gene", "module") -- kept for backward compat with
        # existing DAF files. New pipelines should use primary_module.
        ("gene", "module") => (OptionalInput, AbstractString,
            "DEPRECATED alias for primary_module (collision with Metacells.jl module axis)"),

        # Gene boolean flags from cells_clean (informational)
        ("gene", "is_forbidden") => (OptionalInput, Bool,
            "Whether gene is forbidden"),
        ("gene", "is_selected") => (OptionalInput, Bool,
            "Whether gene is selected for metacell computation"),
        ("gene", "is_transcription_factor") => (OptionalInput, Bool,
            "Whether gene is a transcription factor"),
        ("gene", "is_regulator") => (OptionalInput, Bool,
            "Whether gene is a regulator"),

        # NEW: marker_rank from Metacells.jl (("gene", "marker_rank") => StorageUnsigned)
        ("gene", "marker_rank") => (OptionalInput, StorageUnsigned,
            "Marker gene ranking from Metacells.jl (1=top marker, typemax=non-marker)"),

        # NEW: significant_inner_folds_count from Metacells.jl anndata import
        ("gene", "significant_inner_folds_count") => (OptionalInput, StorageUnsigned,
            "Count of significant inner folds per gene"),

        # ==================================================================
        # Optional per-type vectors
        # ==================================================================
        # NEW: type display order
        ("type", "order") => (OptionalInput, StorageUnsigned,
            "Display order of cell types in plots (1=first)"),

        # ==================================================================
        # Optional gene x metacell matrices
        # ==================================================================
        ("gene", "metacell", "inner_fold") => (OptionalInput, StorageFloat,
            "Inner fold per gene-metacell (enables Inner-fold / QC tabs)"),

        # inner_std_log: Metacells.jl canonical name (preferred)
        ("gene", "metacell", "inner_std_log") => (OptionalInput, StorageFloat,
            "Inner std log per gene-metacell (Metacells.jl canonical name)"),
        # inner_stdev_log: backward-compatible alias
        ("gene", "metacell", "inner_stdev_log") => (OptionalInput, StorageFloat,
            "ALIAS for inner_std_log (legacy name, kept for backward compat)"),

        # geomean_fraction: MCView's own computation (geometric mean of per-cell fractions)
        ("gene", "metacell", "geomean_fraction") => (OptionalInput, StorageFloat,
            "Geometric mean fraction per gene-metacell (MCView computation)"),
        # linear_fraction: Metacells.jl canonical (UMIs / total_UMIs, different from geomean)
        ("gene", "metacell", "linear_fraction") => (OptionalInput, StorageFloat,
            "Linear fraction per gene-metacell (Metacells.jl canonical; UMIs/total_UMIs)"),

        ("gene", "metacell", "zeros") => (OptionalInput, StorageFloat,
            "Zero counts per gene-metacell (enables Zero-fold tab)"),

        # ==================================================================
        # MC x MC matrices -- REMOVED
        # ==================================================================
        # NOTE: The following MC x MC matrices were removed because no R code
        # in MCView reads them:
        #   ("metacell", "metacell", "outgoing_weights")
        #   ("metacell", "metacell", "obs_balanced_ranks")
        #   ("metacell", "metacell", "umap_distances")

        # ==================================================================
        # Cell-level data (enables Samples tab)
        # ==================================================================
        ("cell", "metacell") => (OptionalInput, AbstractString,
            "Metacell assignment per cell"),
        ("cell", "type") => (OptionalInput, AbstractString,
            "Cell type per cell"),

        # Note: additional cell vectors (embryo, batch_set_id, set, etc.)
        # are dynamic metadata -- no fixed contract entries needed. The
        # mcview_sample_property scalar names which one to use as sample ID.

        # ==================================================================
        # Projection / Atlas data (all optional)
        # ==================================================================
        ("metacell", "projected_type") => (OptionalInput, AbstractString,
            "Projected type from atlas"),
        ("metacell", "projected_correlation") => (OptionalInput, StorageFloat,
            "Correlation with atlas projection"),

        # similar / is_similar: both accepted
        ("metacell", "similar") => (OptionalInput, Bool,
            "Whether similar to atlas"),
        ("metacell", "is_similar") => (OptionalInput, Bool,
            "Alias for similar (alternative naming convention)"),

        ("metacell", "atlas_metacell") => (OptionalInput, AbstractString,
            "Most similar atlas metacell"),
        ("gene", "metacell", "projected_fold") => (OptionalInput, StorageFloat,
            "Projected fold from atlas"),

        # projected_fraction / projected_geomean_fraction: both accepted
        ("metacell", "gene", "projected_fraction") => (OptionalInput, StorageFloat,
            "Projected expression fraction from atlas"),
        ("metacell", "gene", "projected_geomean_fraction") => (OptionalInput, StorageFloat,
            "Alias for projected_fraction"),

        # corrected_fraction / corrected_geomean_fraction: both accepted
        ("metacell", "gene", "corrected_fraction") => (OptionalInput, StorageFloat,
            "Corrected expression fraction from atlas projection"),
        ("metacell", "gene", "corrected_geomean_fraction") => (OptionalInput, StorageFloat,
            "Alias for corrected_fraction"),

        # ==================================================================
        # Marker correlations (OptionalInput -- pipeline CAN pre-compute)
        # MCView derived layer computes as FALLBACK only.
        # Accept either axis name.
        # ==================================================================
        ("mcview_marker_correlations", "gene1") => (OptionalInput, AbstractString,
            "First gene in correlation pair"),
        ("mcview_marker_correlations", "gene2") => (OptionalInput, AbstractString,
            "Second gene in correlation pair"),
        ("mcview_marker_correlations", "cor") => (OptionalInput, StorageFloat,
            "Correlation value"),
        ("mcview_marker_correlations", "type") => (OptionalInput, AbstractString,
            "Cell type context for correlation"),

        # DEPRECATED: Legacy axis name for same data
        ("gg_mc_top_cor", "gene1") => (OptionalInput, AbstractString,
            "DEPRECATED: First gene in correlation pair (legacy axis)"),
        ("gg_mc_top_cor", "gene2") => (OptionalInput, AbstractString,
            "DEPRECATED: Second gene in correlation pair (legacy axis)"),
        ("gg_mc_top_cor", "cor") => (OptionalInput, StorageFloat,
            "DEPRECATED: Correlation value (legacy axis)"),
        ("gg_mc_top_cor", "type") => (OptionalInput, AbstractString,
            "DEPRECATED: Cell type for correlation (legacy axis)"),

        # ==================================================================
        # Configuration scalars
        # ==================================================================
        "mcview_title" => (OptionalInput, AbstractString,
            "Application title displayed in UI header"),
        "mcview_tabs" => (OptionalInput, AbstractString,
            "Comma-separated list of enabled tab names"),
        "mcview_excluded_tabs" => (OptionalInput, AbstractString,
            "Comma-separated list of excluded tab names"),
        "mcview_light_version" => (OptionalInput, Bool,
            "Whether to use the light version of the UI"),
        "mcview_about_markdown" => (OptionalInput, AbstractString,
            "Markdown content for the About tab"),
        "mcview_sample_property" => (OptionalInput, AbstractString,
            "Name of the cell vector to use as sample ID (e.g. embryo, batch_set_id, set)"),
        "mcview_available_tabs" => (OptionalInput, AbstractString,
            "Pre-stored comma-separated list of available tabs (skips runtime detection)"),
        "metacells_algorithm" => (OptionalInput, AbstractString,
            "Metacells algorithm version (e.g. metacells2)"),

        # Projection configuration
        "project_max_projection_fold_factor" => (OptionalInput, StorageFloat,
            "Maximum projection fold factor for atlas tab"),
        "projection_weights_json" => (OptionalInput, AbstractString,
            "Projection weights as JSON (query-atlas metacell mapping)"),
        "query_cell_type_fracs_json" => (OptionalInput, AbstractString,
            "Query cell type fractions as JSON"),

        # QC statistics -- canonical names (no qc_stats_ prefix)
        "n_outliers" => (OptionalInput, StorageFloat,
            "Number of outlier cells"),
        "n_cells_total" => (OptionalInput, StorageFloat,
            "Total number of cells (renamed from qc_stats_n_cells to avoid collision with vector name n_cells)"),
        "n_umis" => (OptionalInput, StorageFloat,
            "Total UMI count"),
        "median_umis_per_metacell" => (OptionalInput, StorageFloat,
            "Median UMIs per metacell"),
        "median_cells_per_metacell" => (OptionalInput, StorageFloat,
            "Median cells per metacell"),

        # QC statistics -- backward-compatible aliases (old qc_stats_ prefix)
        "qc_stats_n_outliers" => (OptionalInput, StorageFloat,
            "ALIAS for n_outliers (legacy name)"),
        "qc_stats_n_cells" => (OptionalInput, StorageFloat,
            "ALIAS for n_cells_total (legacy name)"),
        "qc_stats_n_umis" => (OptionalInput, StorageFloat,
            "ALIAS for n_umis (legacy name)"),
        "qc_stats_median_umis_per_metacell" => (OptionalInput, StorageFloat,
            "ALIAS for median_umis_per_metacell (legacy name)"),
        "qc_stats_median_cells_per_metacell" => (OptionalInput, StorageFloat,
            "ALIAS for median_cells_per_metacell (legacy name)"),
    ]
)

# ==============================================================================
# DERIVED CONTRACT - What the mcview_derived layer provides
# ==============================================================================

"""
    MCVIEW_DERIVED_CONTRACT

What the mcview_derived writable layer guarantees to provide on top of the
input DAF. These are computed during MCView's cache-building phase.

The derived layer is a writable DAF that chains on top of the read-only
input DAF. Any data already present in the input (e.g. marker correlations
pre-computed by the pipeline) will NOT be overwritten -- the derived layer
only fills in what is missing.

Type markers note: Per-type marker data is stored as two sparse matrices on
(type, gene) axes rather than as an axis with vectors. This is the correct
DAF data model: the data is naturally indexed by (type, gene) and sparse
matrices handle the sparsity (most gene x type pairs are not markers).
"""
MCVIEW_DERIVED_CONTRACT = Contract(
    is_relaxed = true,
    axes = [
        # Marker correlations (fallback -- only computed if not in input)
        # Kept as axis: pragmatic compromise because the data is inherently
        # 3D (gene x gene x type). Flattening into an axis with vectors
        # (gene1, gene2, cor, type) is the most practical DAF representation.
        "mcview_marker_correlations" => (GuaranteedOutput,
            "Marker gene correlation pairs (fallback if not in input DAF)"),
    ],
    data = [
        # ==================================================================
        # Top-2 genes per metacell (hover info, gene coloring)
        # ==================================================================
        ("metacell", "mcview_cache_top1_gene") => (GuaranteedOutput, AbstractString,
            "Top expressed gene per metacell (by log fold-enrichment)"),
        ("metacell", "mcview_cache_top2_gene") => (GuaranteedOutput, AbstractString,
            "Second top expressed gene per metacell"),
        ("metacell", "mcview_cache_top1_lfp") => (GuaranteedOutput, StorageFloat,
            "Log fold-enrichment of top gene per metacell"),
        ("metacell", "mcview_cache_top2_lfp") => (GuaranteedOutput, StorageFloat,
            "Log fold-enrichment of second gene per metacell"),

        # ==================================================================
        # Gene statistics (for gene selectors and filtering)
        # ==================================================================
        ("gene", "mcview_cache_gene_max_umis") => (GuaranteedOutput, StorageFloat,
            "Maximum UMIs across metacells per gene"),
        ("gene", "mcview_cache_gene_mean_umis") => (GuaranteedOutput, StorageFloat,
            "Mean UMIs across metacells per gene"),
        ("gene", "mcview_cache_gene_sum_umis") => (GuaranteedOutput, StorageFloat,
            "Sum of UMIs across metacells per gene"),

        # ==================================================================
        # Per-type marker genes -- sparse matrices on (type, gene)
        # Replaces the old mcview_type_markers axis + vectors approach.
        # Sparse because most (type, gene) pairs are not markers.
        # ==================================================================
        ("type", "gene", "mcview_marker_rank") => (GuaranteedOutput, StorageUnsigned,
            "Per-type marker rank (1=top marker for that type; 0 or missing=not a marker)"),
        ("type", "gene", "mcview_marker_fold_change") => (GuaranteedOutput, StorageFloat,
            "Per-type marker fold change (gene expression in type vs. background)"),

        # ==================================================================
        # Marker correlations (fallback computation)
        # If the input DAF already provides mcview_marker_correlations or
        # gg_mc_top_cor, the derived layer does NOT recompute -- the input
        # data takes precedence via the chain lookup order.
        # ==================================================================
        ("mcview_marker_correlations", "gene1") => (GuaranteedOutput, AbstractString,
            "First gene in correlation pair"),
        ("mcview_marker_correlations", "gene2") => (GuaranteedOutput, AbstractString,
            "Second gene in correlation pair"),
        ("mcview_marker_correlations", "cor") => (GuaranteedOutput, StorageFloat,
            "Correlation value"),
        ("mcview_marker_correlations", "type") => (GuaranteedOutput, AbstractString,
            "Cell type context for correlation"),

        # ==================================================================
        # Cache metadata scalars
        # ==================================================================
        "mcview_cache_created" => (GuaranteedOutput, AbstractString,
            "ISO 8601 timestamp when cache was created"),
        "mcview_cache_base_version" => (GuaranteedOutput, AbstractString,
            "MCView version that created this cache"),
        "mcview_cache_base_hash" => (GuaranteedOutput, AbstractString,
            "Hash of the input DAF at cache creation time (for staleness detection)"),

        # ==================================================================
        # QC vectors (max inner fold per metacell)
        # ==================================================================
        ("metacell", "mcview_cache_max_inner_fold") => (GuaranteedOutput, StorageFloat,
            "Max inner_fold per metacell across all genes"),
        ("metacell", "mcview_cache_max_inner_fold_no_lateral") => (GuaranteedOutput, StorageFloat,
            "Max inner_fold per metacell excluding lateral genes"),

        # ==================================================================
        # Pre-computed EGC (geomean_fraction fallback)
        # ==================================================================
        ("gene", "metacell", "mcview_cache_geomean_fraction") => (GuaranteedOutput, StorageFloat,
            "Pre-computed geomean_fraction when not in input DAF"),

        # ==================================================================
        # Default markers and distance matrix
        # ==================================================================
        # NOTE: mcview_default_markers is a scalar (comma-separated string).
        # Ideally this should be a vector on a dedicated axis, but the current
        # CSV-in-scalar approach is pragmatic and works. Revisit if the number
        # of default markers grows significantly.
        "mcview_default_markers" => (GuaranteedOutput, AbstractString,
            "Comma-separated list of default marker gene names for heatmap ordering"),
        ("metacell", "metacell", "mcview_default_markers_dist") => (GuaranteedOutput, StorageFloat,
            "Pre-computed hclust distance matrix for default markers"),

        # ==================================================================
        # Per-gene max inner fold (optional derived)
        # ==================================================================
        ("gene", "mcview_cache_gene_max_inner_fold") => (OptionalOutput, StorageFloat,
            "Max inner_fold per gene across metacells"),

        # ==================================================================
        # Cell grouping fields (optional derived)
        # ==================================================================
        "mcview_cell_grouping_fields" => (OptionalOutput, AbstractString,
            "Comma-separated list of categorical cell vectors suitable for grouping in Samples tab"),
    ]
)

# ==============================================================================
# Tab-Specific Contracts (Reference)
#
# These document per-tab data requirements. Each tab implicitly requires
# the core data from MCVIEW_INPUT_CONTRACT. The additional items listed
# here indicate what EXTRA data enables or enriches each tab.
#
# Source column indicates where the data typically comes from:
#   [input]   - from the input DAF (pipeline output)
#   [derived] - from the mcview_derived layer
#   [either]  - pipeline can provide, derived layer fills as fallback
# ==============================================================================

"""
Manifold tab: 2D projection with cell type coloring.
- Requires: core + at least one coordinate pair (x,y) or (u,v)
- Enriched by: w coordinate [input]
- NOTE: MC x MC outgoing_weights removed -- not read by R code
"""
MCVIEW_MANIFOLD_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("metacell", "x")             => (OptionalInput, StorageFloat, "X coordinate"),
        ("metacell", "y")             => (OptionalInput, StorageFloat, "Y coordinate"),
        ("metacell", "u")             => (OptionalInput, StorageFloat, "U coordinate"),
        ("metacell", "v")             => (OptionalInput, StorageFloat, "V coordinate"),
        ("metacell", "w")             => (OptionalInput, StorageFloat, "W coordinate"),
        # Aliases for Metacells.jl naming
        ("metacell", "umap_x")        => (OptionalInput, StorageFloat, "Alias for x"),
        ("metacell", "umap_y")        => (OptionalInput, StorageFloat, "Alias for y"),
        ("metacell", "umap_u")        => (OptionalInput, StorageFloat, "Alias for u"),
        ("metacell", "umap_v")        => (OptionalInput, StorageFloat, "Alias for v"),
        ("metacell", "umap_w")        => (OptionalInput, StorageFloat, "Alias for w"),
        ("metacell", "n_cells")       => (OptionalInput, StorageFloat, "Cells per metacell"),
        ("metacell", "n_cell")        => (OptionalInput, StorageFloat, "Alias for n_cells (legacy)"),
    ]
)

"""
Genes tab: per-gene expression across metacells.
- Requires: core
- Enriched by: top-2 genes [derived], is_lateral/is_noisy [input]
"""
MCVIEW_GENES_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("gene", "is_lateral")        => (OptionalInput, Bool, "Lateral gene [input]"),
        ("gene", "is_noisy")          => (OptionalInput, Bool, "Noisy gene [input]"),
        # Top genes come from derived layer
        ("metacell", "mcview_cache_top1_gene") => (OptionalInput, AbstractString,
            "Top gene [derived]"),
        ("metacell", "mcview_cache_top2_gene") => (OptionalInput, AbstractString,
            "Second gene [derived]"),
        ("metacell", "mcview_cache_top1_lfp")  => (OptionalInput, StorageFloat,
            "Top gene LFP [derived]"),
        ("metacell", "mcview_cache_top2_lfp")  => (OptionalInput, StorageFloat,
            "Second gene LFP [derived]"),
    ]
)

"""
Markers tab: marker gene heatmap.
- Requires: core + is_marker [input]
- Enriched by: is_lateral, is_noisy [input], marker_rank [input]
"""
MCVIEW_MARKERS_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("gene", "is_marker")         => (RequiredInput, Bool,
            "Marker gene flag [input] -- required for this tab"),
        ("gene", "is_lateral")        => (OptionalInput, Bool, "Lateral flag [input]"),
        ("gene", "is_noisy")          => (OptionalInput, Bool, "Noisy flag [input]"),
        ("gene", "marker_rank")       => (OptionalInput, StorageUnsigned,
            "Marker rank from Metacells.jl [input]"),
    ]
)

"""
Gene Modules tab: gene module heatmap.
- Requires: core + gene module assignment [input]
- Accepts both primary_module (canonical) and module (deprecated)
"""
MCVIEW_GENE_MODULES_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        # Accept both names; R layer checks which is present
        ("gene", "primary_module")    => (OptionalInput, AbstractString,
            "Gene module assignment [input] (canonical name)"),
        ("gene", "module")            => (OptionalInput, AbstractString,
            "DEPRECATED alias for primary_module [input]"),
    ]
)

"""
Gene Correlation tab: gene-gene correlations per type.
- Requires: core + marker correlations [either input or derived]
- gg_mc_top_cor is DEPRECATED but still accepted
"""
MCVIEW_GENE_CORRELATION_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
        "mcview_marker_correlations" => (OptionalInput,
            "Correlation pairs [either]"),
        # DEPRECATED
        "gg_mc_top_cor" => (OptionalInput,
            "DEPRECATED: Correlation pairs legacy axis [either]"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("mcview_marker_correlations", "gene1") => (OptionalInput, AbstractString, "Gene 1"),
        ("mcview_marker_correlations", "gene2") => (OptionalInput, AbstractString, "Gene 2"),
        ("mcview_marker_correlations", "cor")   => (OptionalInput, StorageFloat, "Correlation"),
        ("mcview_marker_correlations", "type")  => (OptionalInput, AbstractString, "Cell type"),
        # DEPRECATED legacy axis vectors
        ("gg_mc_top_cor", "gene1") => (OptionalInput, AbstractString, "DEPRECATED: Gene 1 (legacy)"),
        ("gg_mc_top_cor", "gene2") => (OptionalInput, AbstractString, "DEPRECATED: Gene 2 (legacy)"),
        ("gg_mc_top_cor", "cor")   => (OptionalInput, StorageFloat, "DEPRECATED: Correlation (legacy)"),
        ("gg_mc_top_cor", "type")  => (OptionalInput, AbstractString, "DEPRECATED: Cell type (legacy)"),
    ]
)

"""
Inner-fold tab: inner fold QC per gene-metacell.
- Requires: core + inner_fold matrix [input]
- Enriched by: is_lateral [input]
"""
MCVIEW_INNER_FOLD_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("gene", "metacell", "inner_fold") => (RequiredInput, StorageFloat,
            "Inner fold matrix [input]"),
        ("gene", "is_lateral")        => (OptionalInput, Bool, "Lateral flag [input]"),
    ]
)

"""
Stdev-fold tab: inner standard deviation per gene-metacell.
- Requires: core + inner_std_log OR inner_stdev_log matrix [input]
  (The contract lists both as OptionalInput; the R layer checks that at
  least one is present.)
- inner_std_log is the Metacells.jl canonical name.
"""
MCVIEW_STDEV_FOLD_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        # inner_std_log: Metacells.jl canonical name
        ("gene", "metacell", "inner_std_log") => (OptionalInput, StorageFloat,
            "Inner std log [input] (Metacells.jl canonical name)"),
        # inner_stdev_log: backward-compatible alias
        ("gene", "metacell", "inner_stdev_log") => (OptionalInput, StorageFloat,
            "ALIAS for inner_std_log [input] (legacy name)"),
    ]
)

"""
QC tab: quality control metrics.
- Requires: core
- Enriched by: n_cells, inner_fold, QC scalars [input]
- QC scalars accept both new names (no prefix) and old qc_stats_ prefix.
"""
MCVIEW_QC_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("metacell", "n_cells")       => (OptionalInput, StorageFloat, "Cells per metacell"),
        ("metacell", "n_cell")        => (OptionalInput, StorageFloat, "Alias for n_cells (legacy)"),
        ("metacell", "is_rare")       => (OptionalInput, Bool, "Rare metacell flag"),
        ("gene", "metacell", "inner_fold") => (OptionalInput, StorageFloat, "Inner fold"),
        # Canonical QC scalar names
        "n_outliers"                      => (OptionalInput, StorageFloat, "Outlier count"),
        "n_cells_total"                   => (OptionalInput, StorageFloat, "Total cells"),
        "n_umis"                          => (OptionalInput, StorageFloat, "Total UMIs"),
        "median_umis_per_metacell"        => (OptionalInput, StorageFloat, "Median UMIs/mc"),
        "median_cells_per_metacell"       => (OptionalInput, StorageFloat, "Median cells/mc"),
        # Legacy aliases
        "qc_stats_n_outliers"             => (OptionalInput, StorageFloat, "ALIAS: Outlier count"),
        "qc_stats_n_cells"               => (OptionalInput, StorageFloat, "ALIAS: Total cells"),
        "qc_stats_n_umis"                => (OptionalInput, StorageFloat, "ALIAS: Total UMIs"),
        "qc_stats_median_umis_per_metacell"  => (OptionalInput, StorageFloat, "ALIAS: Median UMIs/mc"),
        "qc_stats_median_cells_per_metacell" => (OptionalInput, StorageFloat, "ALIAS: Median cells/mc"),
    ]
)

"""
Samples tab: per-sample cell composition.
- Requires: core + cell axis + cell->metacell mapping [input]
- The sample ID vector is named by mcview_sample_property scalar.
  Cell metadata vectors are dynamic -- no fixed contract entries.
"""
MCVIEW_SAMPLES_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
        "cell"     => (RequiredInput, "Cell identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("cell", "metacell")          => (RequiredInput, AbstractString,
            "Metacell assignment per cell [input]"),
        ("cell", "type")              => (OptionalInput, AbstractString,
            "Cell type per cell [input]"),
        "mcview_sample_property"      => (OptionalInput, AbstractString,
            "Name of cell vector to use as sample ID"),
    ]
)

"""
Projection QC tab: atlas projection quality control.
- Requires: core + projection vectors [input]
- Accepts both similar and is_similar
- Accepts both projected_fraction and projected_geomean_fraction
"""
MCVIEW_PROJECTION_QC_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("metacell", "projected_type")        => (OptionalInput, AbstractString, "Projected type"),
        ("metacell", "projected_correlation") => (OptionalInput, StorageFloat, "Projection corr"),
        ("metacell", "similar")               => (OptionalInput, Bool, "Similar to atlas"),
        ("metacell", "is_similar")            => (OptionalInput, Bool, "Alias for similar"),
        ("metacell", "atlas_metacell")        => (OptionalInput, AbstractString, "Atlas MC"),
        ("gene", "metacell", "projected_fold")       => (OptionalInput, StorageFloat, "Projected fold"),
        ("metacell", "gene", "projected_fraction")   => (OptionalInput, StorageFloat, "Projected frac"),
        ("metacell", "gene", "projected_geomean_fraction") => (OptionalInput, StorageFloat, "Alias for projected_fraction"),
        ("metacell", "gene", "corrected_fraction")   => (OptionalInput, StorageFloat, "Corrected frac"),
        ("metacell", "gene", "corrected_geomean_fraction") => (OptionalInput, StorageFloat, "Alias for corrected_fraction"),
        "project_max_projection_fold_factor"         => (OptionalInput, StorageFloat, "Max fold factor"),
        "projection_weights_json"                    => (OptionalInput, AbstractString, "Weights JSON"),
        "query_cell_type_fracs_json"                 => (OptionalInput, AbstractString, "Type fracs JSON"),
    ]
)

"""
Atlas tab: atlas comparison visualization.
- Requires: core + projection metadata [input]
"""
MCVIEW_ATLAS_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("metacell", "projected_type")        => (OptionalInput, AbstractString, "Projected type"),
        ("metacell", "projected_correlation") => (OptionalInput, StorageFloat, "Projection corr"),
        ("metacell", "similar")               => (OptionalInput, Bool, "Similar to atlas"),
        ("metacell", "is_similar")            => (OptionalInput, Bool, "Alias for similar"),
        ("metacell", "gene", "corrected_fraction")   => (OptionalInput, StorageFloat, "Corrected frac"),
        ("metacell", "gene", "corrected_geomean_fraction") => (OptionalInput, StorageFloat, "Alias for corrected_fraction"),
        ("metacell", "gene", "projected_fraction")   => (OptionalInput, StorageFloat, "Projected frac"),
        ("metacell", "gene", "projected_geomean_fraction") => (OptionalInput, StorageFloat, "Alias for projected_fraction"),
        "projection_weights_json"      => (OptionalInput, AbstractString, "Weights JSON"),
        "query_cell_type_fracs_json"   => (OptionalInput, AbstractString, "Type fracs JSON"),
    ]
)

"""
Flow tab: temporal / flow analysis.
- Requires: core + time annotation [input]
"""
MCVIEW_FLOW_TAB = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),
        ("metacell", "time")          => (OptionalInput, StorageFloat, "Time point [input]"),
        ("metacell", "flow_score")    => (OptionalInput, StorageFloat, "Flow score [input]"),
    ]
)

# ==============================================================================
# Combined Full Contract (relaxed union of input + derived)
# ==============================================================================

"""
    MCVIEW_FULL_CONTRACT

Relaxed union of MCVIEW_INPUT_CONTRACT and MCVIEW_DERIVED_CONTRACT.
This represents the full set of data visible through the chain
(mcview_derived on top of input DAF). Used for documentation and
comprehensive validation of a fully populated MCView dataset.
"""
MCVIEW_FULL_CONTRACT = Contract(
    is_relaxed = true,
    axes = [
        # Input axes
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene"     => (RequiredInput, "Gene identifiers"),
        "type"     => (RequiredInput, "Cell type identifiers"),
        "cell"     => (OptionalInput, "Cell identifiers (Samples tab)"),

        # Marker correlations (from input or derived)
        # Kept as axis: pragmatic compromise for 3D gene x gene x type data
        "mcview_marker_correlations" => (OptionalInput,
            "Marker correlations (from input or derived)"),
        # DEPRECATED
        "gg_mc_top_cor" => (OptionalInput,
            "DEPRECATED: Marker correlations (legacy axis name)"),
    ],
    data = [
        # ---- Core required (input) ----
        ("metacell", "gene", "UMIs")  => (RequiredInput, StorageFloat, "UMI counts"),
        ("metacell", "total_UMIs")    => (RequiredInput, StorageFloat, "Total UMIs"),
        ("metacell", "type")          => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color")             => (RequiredInput, AbstractString, "Color"),

        # ---- Coordinates (input) ----
        ("metacell", "x") => (OptionalInput, StorageFloat, "X coordinate"),
        ("metacell", "y") => (OptionalInput, StorageFloat, "Y coordinate"),
        ("metacell", "u") => (OptionalInput, StorageFloat, "U coordinate"),
        ("metacell", "v") => (OptionalInput, StorageFloat, "V coordinate"),
        ("metacell", "w") => (OptionalInput, StorageFloat, "W coordinate"),
        # Metacells.jl aliases
        ("metacell", "umap_x") => (OptionalInput, StorageFloat, "Alias for x"),
        ("metacell", "umap_y") => (OptionalInput, StorageFloat, "Alias for y"),
        ("metacell", "umap_u") => (OptionalInput, StorageFloat, "Alias for u"),
        ("metacell", "umap_v") => (OptionalInput, StorageFloat, "Alias for v"),
        ("metacell", "umap_w") => (OptionalInput, StorageFloat, "Alias for w"),

        # ---- Optional metacell vectors (input) ----
        ("metacell", "n_cells") => (OptionalInput, StorageFloat, "Cells per metacell (canonical)"),
        ("metacell", "n_cell")  => (OptionalInput, StorageFloat, "Alias for n_cells (legacy)"),
        ("metacell", "is_rare") => (OptionalInput, Bool, "Rare metacell"),

        # ---- Optional gene vectors (input) ----
        ("gene", "is_marker")  => (OptionalInput, Bool, "Marker gene"),
        ("gene", "is_lateral") => (OptionalInput, Bool, "Lateral gene"),
        ("gene", "is_noisy")   => (OptionalInput, Bool, "Noisy gene"),
        ("gene", "primary_module") => (OptionalInput, AbstractString, "Gene module (canonical)"),
        ("gene", "module")     => (OptionalInput, AbstractString, "DEPRECATED alias for primary_module"),
        ("gene", "is_forbidden")            => (OptionalInput, Bool, "Forbidden gene"),
        ("gene", "is_selected")             => (OptionalInput, Bool, "Selected gene"),
        ("gene", "is_transcription_factor") => (OptionalInput, Bool, "TF gene"),
        ("gene", "is_regulator")            => (OptionalInput, Bool, "Regulator gene"),
        ("gene", "marker_rank") => (OptionalInput, StorageUnsigned, "Marker rank from Metacells.jl"),
        ("gene", "significant_inner_folds_count") => (OptionalInput, StorageUnsigned, "Significant inner folds count"),

        # ---- Per-type vectors (input) ----
        ("type", "order") => (OptionalInput, StorageUnsigned, "Type display order"),

        # ---- Gene x metacell matrices (input) ----
        ("gene", "metacell", "inner_fold")      => (OptionalInput, StorageFloat, "Inner fold"),
        ("gene", "metacell", "inner_std_log")   => (OptionalInput, StorageFloat, "Inner std log (canonical)"),
        ("gene", "metacell", "inner_stdev_log") => (OptionalInput, StorageFloat, "Alias for inner_std_log (legacy)"),
        ("gene", "metacell", "geomean_fraction") => (OptionalInput, StorageFloat, "Geomean fraction (MCView)"),
        ("gene", "metacell", "linear_fraction") => (OptionalInput, StorageFloat, "Linear fraction (Metacells.jl)"),
        ("gene", "metacell", "zeros")           => (OptionalInput, StorageFloat, "Zeros"),

        # ---- MC x MC matrices -- REMOVED ----
        # outgoing_weights, obs_balanced_ranks, umap_distances removed (no R code reads them)

        # ---- Cell data (input) ----
        ("cell", "metacell") => (OptionalInput, AbstractString, "Cell metacell"),
        ("cell", "type")     => (OptionalInput, AbstractString, "Cell type"),

        # ---- Projection data (input) ----
        ("metacell", "projected_type")        => (OptionalInput, AbstractString, "Projected type"),
        ("metacell", "projected_correlation") => (OptionalInput, StorageFloat, "Projection corr"),
        ("metacell", "similar")               => (OptionalInput, Bool, "Similar to atlas"),
        ("metacell", "is_similar")            => (OptionalInput, Bool, "Alias for similar"),
        ("metacell", "atlas_metacell")        => (OptionalInput, AbstractString, "Atlas MC"),
        ("gene", "metacell", "projected_fold")       => (OptionalInput, StorageFloat, "Projected fold"),
        ("metacell", "gene", "projected_fraction")   => (OptionalInput, StorageFloat, "Projected frac"),
        ("metacell", "gene", "projected_geomean_fraction") => (OptionalInput, StorageFloat, "Alias for projected_fraction"),
        ("metacell", "gene", "corrected_fraction")   => (OptionalInput, StorageFloat, "Corrected frac"),
        ("metacell", "gene", "corrected_geomean_fraction") => (OptionalInput, StorageFloat, "Alias for corrected_fraction"),

        # ---- Marker correlations (input or derived) ----
        ("mcview_marker_correlations", "gene1") => (OptionalInput, AbstractString, "Corr gene 1"),
        ("mcview_marker_correlations", "gene2") => (OptionalInput, AbstractString, "Corr gene 2"),
        ("mcview_marker_correlations", "cor")   => (OptionalInput, StorageFloat, "Correlation"),
        ("mcview_marker_correlations", "type")  => (OptionalInput, AbstractString, "Corr cell type"),
        # DEPRECATED legacy
        ("gg_mc_top_cor", "gene1") => (OptionalInput, AbstractString, "DEPRECATED: Corr gene 1 (legacy)"),
        ("gg_mc_top_cor", "gene2") => (OptionalInput, AbstractString, "DEPRECATED: Corr gene 2 (legacy)"),
        ("gg_mc_top_cor", "cor")   => (OptionalInput, StorageFloat, "DEPRECATED: Correlation (legacy)"),
        ("gg_mc_top_cor", "type")  => (OptionalInput, AbstractString, "DEPRECATED: Corr cell type (legacy)"),

        # ---- Top-2 genes (derived) ----
        ("metacell", "mcview_cache_top1_gene") => (OptionalInput, AbstractString, "Top gene"),
        ("metacell", "mcview_cache_top2_gene") => (OptionalInput, AbstractString, "Second gene"),
        ("metacell", "mcview_cache_top1_lfp")  => (OptionalInput, StorageFloat, "Top gene LFP"),
        ("metacell", "mcview_cache_top2_lfp")  => (OptionalInput, StorageFloat, "Second gene LFP"),

        # ---- Gene stats (derived) ----
        ("gene", "mcview_cache_gene_max_umis")  => (OptionalInput, StorageFloat, "Max UMIs"),
        ("gene", "mcview_cache_gene_mean_umis") => (OptionalInput, StorageFloat, "Mean UMIs"),
        ("gene", "mcview_cache_gene_sum_umis")  => (OptionalInput, StorageFloat, "Sum UMIs"),

        # ---- Per-type markers (derived) -- sparse matrices on (type, gene) ----
        ("type", "gene", "mcview_marker_rank")        => (OptionalInput, StorageUnsigned, "Per-type marker rank"),
        ("type", "gene", "mcview_marker_fold_change") => (OptionalInput, StorageFloat, "Per-type marker fold change"),

        # ---- Config scalars (input) ----
        "mcview_title"          => (OptionalInput, AbstractString, "App title"),
        "mcview_tabs"           => (OptionalInput, AbstractString, "Enabled tabs"),
        "mcview_excluded_tabs"  => (OptionalInput, AbstractString, "Excluded tabs"),
        "mcview_light_version"  => (OptionalInput, Bool, "Light version"),
        "mcview_about_markdown" => (OptionalInput, AbstractString, "About markdown"),
        "mcview_sample_property" => (OptionalInput, AbstractString, "Sample ID property name"),
        "mcview_available_tabs" => (OptionalInput, AbstractString, "Pre-stored available tabs"),
        "metacells_algorithm"   => (OptionalInput, AbstractString, "Algorithm version"),
        "project_max_projection_fold_factor" => (OptionalInput, StorageFloat, "Max fold factor"),
        "projection_weights_json"      => (OptionalInput, AbstractString, "Weights JSON"),
        "query_cell_type_fracs_json"   => (OptionalInput, AbstractString, "Type fracs JSON"),

        # QC scalars -- canonical names
        "n_outliers"                => (OptionalInput, StorageFloat, "Outlier count"),
        "n_cells_total"             => (OptionalInput, StorageFloat, "Total cells"),
        "n_umis"                    => (OptionalInput, StorageFloat, "Total UMIs"),
        "median_umis_per_metacell"  => (OptionalInput, StorageFloat, "Median UMIs/mc"),
        "median_cells_per_metacell" => (OptionalInput, StorageFloat, "Median cells/mc"),
        # QC scalars -- legacy aliases
        "qc_stats_n_outliers"                => (OptionalInput, StorageFloat, "ALIAS: Outlier count"),
        "qc_stats_n_cells"                   => (OptionalInput, StorageFloat, "ALIAS: Total cells"),
        "qc_stats_n_umis"                    => (OptionalInput, StorageFloat, "ALIAS: Total UMIs"),
        "qc_stats_median_umis_per_metacell"  => (OptionalInput, StorageFloat, "ALIAS: Median UMIs/mc"),
        "qc_stats_median_cells_per_metacell" => (OptionalInput, StorageFloat, "ALIAS: Median cells/mc"),

        # ---- Cache metadata (derived) ----
        "mcview_cache_created"      => (OptionalInput, AbstractString, "Cache timestamp"),
        "mcview_cache_base_version" => (OptionalInput, AbstractString, "Cache MCView version"),
        "mcview_cache_base_hash"    => (OptionalInput, AbstractString, "Cache input hash"),

        # ---- QC vectors (derived) ----
        ("metacell", "mcview_cache_max_inner_fold")            => (OptionalInput, StorageFloat, "Max inner fold"),
        ("metacell", "mcview_cache_max_inner_fold_no_lateral") => (OptionalInput, StorageFloat, "Max inner fold (no lateral)"),

        # ---- Pre-computed EGC (derived) ----
        ("gene", "metacell", "mcview_cache_geomean_fraction") => (OptionalInput, StorageFloat, "Cached geomean fraction"),

        # ---- Default markers (derived) ----
        # NOTE: Ideally this should be a vector on a dedicated axis, but the
        # CSV-in-scalar approach is pragmatic and works for now.
        "mcview_default_markers" => (OptionalInput, AbstractString, "Default marker genes"),
        ("metacell", "metacell", "mcview_default_markers_dist") => (OptionalInput, StorageFloat, "Default markers dist matrix"),

        # ---- Per-gene max inner fold (derived) ----
        ("gene", "mcview_cache_gene_max_inner_fold") => (OptionalInput, StorageFloat, "Max inner fold per gene"),

        # ---- Cell grouping fields (derived) ----
        "mcview_cell_grouping_fields" => (OptionalInput, AbstractString, "Cell grouping fields"),
    ]
)

# ==============================================================================
# Validation Helper Functions
# ==============================================================================

"""
    validate_mcview_input(daf::DafReader)

Validate that a DAF object meets the MCView input requirements.
Returns `true` on success, throws on failure.
"""
function validate_mcview_input(daf::DafReader)
    verify_input(contractor("MCView.Input", MCVIEW_INPUT_CONTRACT, daf))
    return true
end

"""
    validate_mcview_derived(daf::DafReader)

Validate that a DAF object (typically the chained view) has the
mcview_derived layer's guaranteed outputs. Call this on the chain,
not on the derived layer alone.
"""
function validate_mcview_derived(daf::DafReader)
    verify_output(contractor("MCView.Derived", MCVIEW_DERIVED_CONTRACT, daf))
    return true
end

"""
    validate_mcview_tab(daf::DafReader, tab::AbstractString)

Validate that a DAF object has sufficient data for a specific tab.
"""
function validate_mcview_tab(daf::DafReader, tab::AbstractString)
    contracts = Dict(
        "Manifold"         => MCVIEW_MANIFOLD_TAB,
        "Genes"            => MCVIEW_GENES_TAB,
        "Markers"          => MCVIEW_MARKERS_TAB,
        "Gene modules"     => MCVIEW_GENE_MODULES_TAB,
        "Gene correlation" => MCVIEW_GENE_CORRELATION_TAB,
        "Inner-fold"       => MCVIEW_INNER_FOLD_TAB,
        "Stdev-fold"       => MCVIEW_STDEV_FOLD_TAB,
        "QC"               => MCVIEW_QC_TAB,
        "Samples"          => MCVIEW_SAMPLES_TAB,
        "Projection QC"    => MCVIEW_PROJECTION_QC_TAB,
        "Atlas"            => MCVIEW_ATLAS_TAB,
        "Flow"             => MCVIEW_FLOW_TAB,
    )

    if !haskey(contracts, tab)
        error("Unknown tab: $tab. Known tabs: $(join(sort(collect(keys(contracts))), ", "))")
    end

    verify_input(contractor("MCView.$tab", contracts[tab], daf))
    return true
end

"""
    get_available_tabs(daf::DafReader)

Determine which tabs have sufficient data in the DAF.
Returns a Vector{String} of available tab names.

Tabs like "About", "Diff. Expression", "Cell types", and "Annotate"
are always available if the core contract passes (they only need core data).
"""
function get_available_tabs(daf::DafReader)
    # Core validation must pass for anything to work
    try
        validate_mcview_input(daf)
    catch
        return String[]
    end

    # Tabs that are always available with core data
    always_available = ["About", "Manifold", "Genes", "Diff. Expression",
                        "Cell types", "Annotate"]

    # Tabs that need extra data -- check each
    conditional_tabs = [
        "QC", "Markers", "Gene modules", "Gene correlation",
        "Inner-fold", "Stdev-fold", "Projection QC", "Atlas",
        "Samples", "Flow"
    ]

    available = copy(always_available)

    for tab in conditional_tabs
        try
            validate_mcview_tab(daf, tab)
            push!(available, tab)
        catch
            # Tab data not available -- skip
        end
    end

    return available
end
