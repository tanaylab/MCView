# mcview_contracts.jl - Formal DAF Contracts for MCView
#
# This file defines the data contracts for MCView using the Daf.jl Contracts module.
# These contracts specify required and optional data for MCView to function.
#
# Usage:
#   using DataAxesFormats.Contracts
#   include("mcview_contracts.jl")
#   verify_input(contractor("MCView", MCVIEW_CORE_CONTRACT, daf_obj))

using DataAxesFormats.Contracts

# ==============================================================================
# Core MCView Contract
# ==============================================================================

"""
Minimum required data for MCView to function.
All tabs require at least this data.
"""
MCVIEW_CORE_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        # Required matrix
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64,
            "UMI counts per metacell-gene"),

        # Required metacell vectors
        ("metacell", "total_UMIs") => (RequiredInput, Float64,
            "Total UMIs per metacell"),
        ("metacell", "type") => (RequiredInput, AbstractString,
            "Cell type assignment per metacell"),

        # Required type vector
        ("type", "color") => (RequiredInput, AbstractString,
            "Hex color per cell type (#RRGGBB format)"),

        # 2D coordinates - at least one pair required (x,y or u,v)
        # Validated by custom rule in R implementation
        ("metacell", "x") => (OptionalInput, Float64,
            "X coordinate for 2D projection"),
        ("metacell", "y") => (OptionalInput, Float64,
            "Y coordinate for 2D projection"),
        ("metacell", "u") => (OptionalInput, Float64,
            "U coordinate (alternative to x)"),
        ("metacell", "v") => (OptionalInput, Float64,
            "V coordinate (alternative to y)"),

        # Optional metacell vectors
        ("metacell", "n_cell") => (OptionalInput, Int64,
            "Number of cells per metacell"),
        ("metacell", "cells") => (OptionalInput, Int64,
            "Alias for n_cell"),
        ("metacell", "mc_col") => (OptionalInput, AbstractString,
            "Per-metacell color override"),
    ]
)

# ==============================================================================
# Manifold Tab Contract
# ==============================================================================

"""
Contract for 2D projection and cell type visualization.
"""
MCVIEW_MANIFOLD_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
        "metacell_graph" => (OptionalInput, "Metacell graph edges for visualization"),
    ],
    data = [
        # Inherit from core
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),
        ("metacell", "x") => (OptionalInput, Float64, "X coordinate"),
        ("metacell", "y") => (OptionalInput, Float64, "Y coordinate"),
        ("metacell", "u") => (OptionalInput, Float64, "U coordinate"),
        ("metacell", "v") => (OptionalInput, Float64, "V coordinate"),

        # Graph edge data
        ("metacell_graph", "from") => (OptionalInput, AbstractString,
            "Source metacell of edge"),
        ("metacell_graph", "to") => (OptionalInput, AbstractString,
            "Target metacell of edge"),
        ("metacell_graph", "weight") => (OptionalInput, Float64, "Edge weight"),
        ("metacell_graph", "graph_name") => (OptionalInput, AbstractString,
            "Name of the graph"),
    ]
)

# ==============================================================================
# Genes Tab Contract
# ==============================================================================

"""
Contract for gene expression visualization.
"""
MCVIEW_GENES_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Optional gene metadata
        ("metacell", "top1_gene") => (OptionalInput, AbstractString,
            "Top expressed gene per metacell"),
        ("metacell", "top2_gene") => (OptionalInput, AbstractString,
            "Second top expressed gene per metacell"),
        ("metacell", "top1_lfp") => (OptionalInput, Float64,
            "Log fold-enrichment of top gene"),
        ("metacell", "top2_lfp") => (OptionalInput, Float64,
            "Log fold-enrichment of second top gene"),
        ("gene", "is_lateral") => (OptionalInput, Bool,
            "Whether gene is lateral (excluded from analysis)"),
        ("gene", "is_noisy") => (OptionalInput, Bool, "Whether gene is noisy"),
    ]
)

# ==============================================================================
# QC Tab Contract
# ==============================================================================

"""
Contract for quality control metrics.
"""
MCVIEW_QC_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # QC vectors
        ("metacell", "n_cell") => (OptionalInput, Int64,
            "Number of cells per metacell"),
        ("metacell", "cells") => (OptionalInput, Int64, "Alias for n_cell"),
        ("metacell", "max_inner_fold") => (OptionalInput, Float64,
            "Maximum inner fold per metacell"),
        ("metacell", "max_inner_fold_no_lateral") => (OptionalInput, Float64,
            "Max inner fold excluding lateral genes"),
        ("metacell", "max_inner_stdev_log") => (OptionalInput, Float64,
            "Max inner stdev log per metacell"),
        ("metacell", "zero_fold") => (OptionalInput, Float64,
            "Zero fold per metacell"),
        ("metacell", "rare_metacell") => (OptionalInput, Bool,
            "Whether metacell is rare"),
        ("metacell", "metacells_rare_gene_module") => (OptionalInput, Int64,
            "Rare gene module assignment"),

        # QC scalars
        "qc_stats_n_outliers" => (OptionalInput, Int64,
            "Number of outlier cells"),
        "qc_stats_n_cells" => (OptionalInput, Int64, "Total number of cells"),
        "qc_stats_n_umis" => (OptionalInput, Int64, "Total UMI count"),
        "qc_stats_median_umis_per_metacell" => (OptionalInput, Float64,
            "Median UMIs per metacell"),
        "qc_stats_median_cells_per_metacell" => (OptionalInput, Float64,
            "Median cells per metacell"),
    ]
)

# ==============================================================================
# Markers Tab Contract
# ==============================================================================

"""
Contract for marker gene visualization.
Requires is_marker vector on gene axis.
"""
MCVIEW_MARKERS_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Required for Markers tab
        ("gene", "is_marker") => (RequiredInput, Bool, "Whether gene is a marker"),
    ]
)

# ==============================================================================
# Gene Modules Tab Contract
# ==============================================================================

"""
Contract for gene module analysis.
Requires module vector on gene axis.
"""
MCVIEW_GENE_MODULES_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Required for Gene modules tab
        ("gene", "module") => (RequiredInput, AbstractString,
            "Gene module assignment"),
    ]
)

# ==============================================================================
# Inner-fold Tab Contract
# ==============================================================================

"""
Contract for inner fold analysis per gene.
Requires inner_fold matrix.
"""
MCVIEW_INNER_FOLD_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Required for Inner-fold tab
        ("gene", "metacell", "inner_fold") => (RequiredInput, Float64,
            "Inner fold per gene-metacell"),

        # Optional
        ("gene", "is_lateral") => (OptionalInput, Bool, "Whether gene is lateral"),
        ("metacell", "max_inner_fold") => (OptionalInput, Float64,
            "Maximum inner fold per metacell"),
    ]
)

# ==============================================================================
# Stdev-fold Tab Contract
# ==============================================================================

"""
Contract for inner standard deviation analysis.
Requires inner_stdev_log matrix.
"""
MCVIEW_STDEV_FOLD_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Required for Stdev-fold tab
        ("gene", "metacell", "inner_stdev_log") => (RequiredInput, Float64,
            "Inner stdev log per gene-metacell"),
    ]
)

# ==============================================================================
# Projection QC Tab Contract
# ==============================================================================

"""
Contract for atlas projection quality control.
"""
MCVIEW_PROJECTION_QC_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Projection vectors
        ("metacell", "projected_type") => (OptionalInput, AbstractString,
            "Projected type from atlas"),
        ("metacell", "projected_correlation") => (OptionalInput, Float64,
            "Correlation with atlas projection"),
        ("metacell", "similar") => (OptionalInput, Bool,
            "Whether similar to atlas"),
        ("metacell", "atlas_metacell") => (OptionalInput, AbstractString,
            "Most similar atlas metacell"),

        # Projection matrices
        ("gene", "metacell", "projected_fold") => (OptionalInput, Float64,
            "Projected fold from atlas"),
        ("metacell", "gene", "projected_fraction") => (OptionalInput, Float64,
            "Projected expression fraction from atlas"),
        ("metacell", "gene", "corrected_fraction") => (OptionalInput, Float64,
            "Corrected expression fraction from atlas projection"),

        # Projection scalars (serialized data)
        "projection_weights_json" => (OptionalInput, AbstractString,
            "Projection weights as JSON (query-atlas metacell mapping)"),
        "query_cell_type_fracs_json" => (OptionalInput, AbstractString,
            "Query cell type fractions as JSON"),

        # Configuration
        "project_max_projection_fold_factor" => (OptionalInput, Float64,
            "Maximum projection fold factor"),
    ]
)

# ==============================================================================
# Atlas Tab Contract
# ==============================================================================

"""
Contract for atlas visualization.
Requires separate atlas DAF for comparison.
The query DAF should contain projection metadata, while the atlas DAF
contains the reference data.
"""
MCVIEW_ATLAS_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Projection vectors
        ("metacell", "projected_type") => (OptionalInput, AbstractString,
            "Projected type from atlas"),
        ("metacell", "projected_correlation") => (OptionalInput, Float64,
            "Correlation with atlas projection"),
        ("metacell", "similar") => (OptionalInput, Bool,
            "Whether similar to atlas"),

        # Projection matrices
        ("metacell", "gene", "corrected_fraction") => (OptionalInput, Float64,
            "Corrected fraction from atlas projection"),
        ("metacell", "gene", "projected_fraction") => (OptionalInput, Float64,
            "Projected expression fraction from atlas"),

        # Projection scalars (serialized data)
        "projection_weights_json" => (OptionalInput, AbstractString,
            "Projection weights as JSON (query-atlas metacell mapping)"),
        "query_cell_type_fracs_json" => (OptionalInput, AbstractString,
            "Query cell type fractions as JSON"),
    ]
)

# ==============================================================================
# Samples Tab Contract
# ==============================================================================

"""
Contract for sample-level analysis.
Requires cell axis with sample assignments.
"""
MCVIEW_SAMPLES_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
        "cell" => (RequiredInput, "Cell identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Cell metadata
        ("cell", "metacell") => (RequiredInput, AbstractString,
            "Metacell assignment per cell (-1 for outliers)"),
        ("cell", "samp_id") => (RequiredInput, AbstractString,
            "Sample ID per cell"),
        ("cell", "cell_type") => (OptionalInput, AbstractString,
            "Cell type per cell"),
        ("cell", "outlier") => (OptionalInput, Bool,
            "Whether cell is an outlier"),
    ]
)

# ==============================================================================
# Flow Tab Contract
# ==============================================================================

"""
Contract for time/flow analysis.
"""
MCVIEW_FLOW_CONTRACT = Contract(
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
    ],
    data = [
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # Time data
        ("metacell", "time") => (OptionalInput, Float64,
            "Time point per metacell"),
        ("metacell", "flow_score") => (OptionalInput, Float64,
            "Flow score per metacell"),
    ]
)

# ==============================================================================
# Gene Correlations Contract
# ==============================================================================

"""
Contract for gene-gene correlation data.
"""
MCVIEW_GENE_CORRELATIONS_CONTRACT = Contract(
    axes = [
        "gg_mc_top_cor" => (OptionalInput, "Gene correlation pairs"),
    ],
    data = [
        ("gg_mc_top_cor", "gene1") => (OptionalInput, AbstractString,
            "First gene in correlation pair"),
        ("gg_mc_top_cor", "gene2") => (OptionalInput, AbstractString,
            "Second gene in correlation pair"),
        ("gg_mc_top_cor", "cor") => (OptionalInput, Float64, "Correlation value"),
        ("gg_mc_top_cor", "type") => (OptionalInput, AbstractString,
            "Cell type for correlation"),
    ]
)

# ==============================================================================
# Gene Zero Fold Contract
# ==============================================================================

"""
Contract for gene zero-fold statistics.
"""
MCVIEW_GENE_ZERO_FOLD_CONTRACT = Contract(
    axes = [
        "gene_zero_fold" => (OptionalInput, "Gene zero-fold statistics"),
    ],
    data = [
        ("gene_zero_fold", "gene") => (OptionalInput, AbstractString, "Gene name"),
        ("gene_zero_fold", "metacell") => (OptionalInput, AbstractString,
            "Metacell name"),
        ("gene_zero_fold", "zero_fold") => (OptionalInput, Float64,
            "Zero fold value"),
        ("gene_zero_fold", "avg") => (OptionalInput, Float64, "Average value"),
        ("gene_zero_fold", "obs") => (OptionalInput, Float64, "Observed value"),
        ("gene_zero_fold", "exp") => (OptionalInput, Float64, "Expected value"),
        ("gene_zero_fold", "type") => (OptionalInput, AbstractString, "Cell type"),
    ]
)

# ==============================================================================
# Configuration Contract
# ==============================================================================

"""
Contract for MCView-specific configuration stored in DAF.
"""
MCVIEW_CONFIG_CONTRACT = Contract(
    data = [
        "mcview_title" => (OptionalInput, AbstractString, "Application title"),
        "mcview_tabs" => (OptionalInput, AbstractString,
            "Comma-separated list of tabs to show"),
        "mcview_excluded_tabs" => (OptionalInput, AbstractString,
            "Comma-separated list of tabs to exclude"),
        "mcview_light_version" => (OptionalInput, Bool,
            "Whether to use light version"),
        "mcview_about_markdown" => (OptionalInput, AbstractString,
            "About page markdown content"),
        "mcview_cache_in_daf" => (OptionalInput, Bool,
            "Whether to cache computed data in DAF"),
        "mcview_cache_daf_root" => (OptionalInput, AbstractString,
            "Root directory for cache DAF"),
        "metacells_algorithm" => (OptionalInput, AbstractString,
            "Metacells algorithm version"),
    ]
)

# ==============================================================================
# Cache Contract
# ==============================================================================

"""
Contract for precomputed cache vectors.
These are optional but improve first-load performance.
"""
MCVIEW_CACHE_CONTRACT = Contract(
    axes = [
        "mcview_cache_gg_mc_top_cor" => (OptionalInput,
            "Precomputed gene correlation pairs"),
    ],
    data = [
        # Top genes per metacell
        ("metacell", "mcview_cache_top1_gene") => (OptionalInput, AbstractString,
            "Precomputed top expressed gene"),
        ("metacell", "mcview_cache_top2_gene") => (OptionalInput, AbstractString,
            "Precomputed second top expressed gene"),
        ("metacell", "mcview_cache_top1_lfp") => (OptionalInput, Float64,
            "Precomputed log fold-enrichment of top gene"),
        ("metacell", "mcview_cache_top2_lfp") => (OptionalInput, Float64,
            "Precomputed log fold-enrichment of second gene"),

        # Gene statistics
        ("gene", "mcview_cache_max_expr") => (OptionalInput, Float64,
            "Precomputed maximum expression per gene"),
        ("gene", "mcview_cache_total_expr") => (OptionalInput, Float64,
            "Precomputed total expression per gene"),
    ]
)

# ==============================================================================
# Combined Full Contract
# ==============================================================================

"""
Full MCView contract combining all tab-specific requirements.
This is a relaxed contract that allows additional data.
"""
MCVIEW_FULL_CONTRACT = Contract(
    is_relaxed = true,
    axes = [
        "metacell" => (RequiredInput, "Metacell identifiers"),
        "gene" => (RequiredInput, "Gene identifiers"),
        "type" => (RequiredInput, "Cell type identifiers"),
        "cell" => (OptionalInput, "Cell identifiers (for Samples tab)"),
        "metacell_graph" => (OptionalInput, "Metacell graph edges"),
        "gg_mc_top_cor" => (OptionalInput, "Gene correlation pairs"),
        "gene_zero_fold" => (OptionalInput, "Gene zero-fold statistics"),
    ],
    data = [
        # Core required
        ("metacell", "gene", "UMIs") => (RequiredInput, Float64, "UMI counts"),
        ("metacell", "total_UMIs") => (RequiredInput, Float64, "Total UMIs"),
        ("metacell", "type") => (RequiredInput, AbstractString, "Cell type"),
        ("type", "color") => (RequiredInput, AbstractString, "Cell type color"),

        # 2D coordinates (at least one pair required)
        ("metacell", "x") => (OptionalInput, Float64, "X coordinate"),
        ("metacell", "y") => (OptionalInput, Float64, "Y coordinate"),
        ("metacell", "u") => (OptionalInput, Float64, "U coordinate"),
        ("metacell", "v") => (OptionalInput, Float64, "V coordinate"),

        # Optional metacell data
        ("metacell", "n_cell") => (OptionalInput, Int64, "Cells per metacell"),
        ("metacell", "mc_col") => (OptionalInput, AbstractString, "Color override"),

        # Optional gene data
        ("gene", "is_marker") => (OptionalInput, Bool, "Is marker gene"),
        ("gene", "is_lateral") => (OptionalInput, Bool, "Is lateral gene"),
        ("gene", "is_noisy") => (OptionalInput, Bool, "Is noisy gene"),
        ("gene", "module") => (OptionalInput, AbstractString, "Gene module"),

        # Optional matrices
        ("gene", "metacell", "inner_fold") => (OptionalInput, Float64, "Inner fold"),
        ("gene", "metacell", "inner_stdev_log") => (OptionalInput, Float64,
            "Inner stdev log"),
        ("gene", "metacell", "projected_fold") => (OptionalInput, Float64,
            "Projected fold"),

        # Cell data (for Samples tab)
        ("cell", "metacell") => (OptionalInput, AbstractString, "Cell metacell"),
        ("cell", "samp_id") => (OptionalInput, AbstractString, "Sample ID"),
    ]
)

# ==============================================================================
# Validation Helper Functions
# ==============================================================================

"""
    validate_mcview_core(daf::DafReader)

Validate that a DAF object meets the core MCView requirements.
"""
function validate_mcview_core(daf::DafReader)
    verify_input(contractor("MCView.Core", MCVIEW_CORE_CONTRACT, daf))
    return true
end

"""
    validate_mcview_tab(daf::DafReader, tab::AbstractString)

Validate that a DAF object has data for a specific tab.
"""
function validate_mcview_tab(daf::DafReader, tab::AbstractString)
    contracts = Dict(
        "Manifold" => MCVIEW_MANIFOLD_CONTRACT,
        "Genes" => MCVIEW_GENES_CONTRACT,
        "QC" => MCVIEW_QC_CONTRACT,
        "Markers" => MCVIEW_MARKERS_CONTRACT,
        "Gene modules" => MCVIEW_GENE_MODULES_CONTRACT,
        "Inner-fold" => MCVIEW_INNER_FOLD_CONTRACT,
        "Stdev-fold" => MCVIEW_STDEV_FOLD_CONTRACT,
        "Projection QC" => MCVIEW_PROJECTION_QC_CONTRACT,
        "Atlas" => MCVIEW_ATLAS_CONTRACT,
        "Samples" => MCVIEW_SAMPLES_CONTRACT,
        "Flow" => MCVIEW_FLOW_CONTRACT,
    )

    if !haskey(contracts, tab)
        error("Unknown tab: $tab")
    end

    verify_input(contractor("MCView.$tab", contracts[tab], daf))
    return true
end

"""
    get_available_tabs(daf::DafReader)

Get list of tabs that have sufficient data in the DAF.
"""
function get_available_tabs(daf::DafReader)
    # Core validation must pass
    try
        validate_mcview_core(daf)
    catch
        return String[]
    end

    # Check each tab
    all_tabs = ["About", "Manifold", "Genes", "Diff. Expression", "Cell types",
                "QC", "Markers", "Gene modules", "Inner-fold", "Stdev-fold",
                "Projection QC", "Atlas", "Annotate", "Samples", "Flow"]

    available = String[]

    for tab in all_tabs
        try
            validate_mcview_tab(daf, tab)
            push!(available, tab)
        catch
            # Tab not available
        end
    end

    return available
end
