# daf_contracts_primitives.R - Contract primitive helpers + dafr adapters
#
# Split from R/daf_contracts.R (2026-05-01).  See R/daf_contracts_definitions.R
# for the contract definitions, R/daf_contracts_validators.R for validation,
# and R/daf_contracts_docs.R for printing / markdown generation.
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

    vectors <- list()
    all_vectors <- c(contract$vectors %||% list(), contract$graph_vectors %||% list())
    for (spec in all_vectors) {
        vectors <- c(
            vectors,
            list(dafr::contract_vector(
                spec$axis,
                spec$name,
                mcview_expectation_to_dafr(spec$expectation),
                spec$type,
                spec$description
            ))
        )
    }

    matrices <- list()
    if (!is.null(contract$matrices)) {
        for (spec in contract$matrices) {
            matrices <- c(
                matrices,
                list(dafr::contract_matrix(
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

    scalars <- list()
    if (!is.null(contract$scalars)) {
        for (spec in contract$scalars) {
            scalars <- c(
                scalars,
                list(dafr::contract_scalar(
                    spec$name,
                    mcview_expectation_to_dafr(spec$expectation),
                    spec$type,
                    spec$description
                ))
            )
        }
    }

    dafr::create_contract(
        scalars = scalars,
        vectors = vectors,
        matrices = matrices,
        axes = axes,
        is_relaxed = TRUE
    )
}
