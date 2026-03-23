# julia_helpers.R - Julia-accelerated computation helpers for MCView
#
# This module defines Julia functions (via JuliaCall) that replace expensive
# R computations. All Julia-backed functions have automatic fallback to the
# existing R implementations when Julia is not available or when the option
# mcview.use_julia_helpers is FALSE.
#
# Usage:
#   After dafr::setup_daf() completes, call init_julia_helpers() to load
#   the Julia helper functions. The R wrapper functions can then be used
#   as drop-in replacements.
#
# Option: getOption("mcview.use_julia_helpers", TRUE)
#   Set to FALSE to disable Julia helpers and always use R fallbacks.

# ==============================================================================
# Internal State
# ==============================================================================

# Use an environment for mutable state so that it survives namespace locking
# (top-level scalar bindings get locked by devtools::load_all / package sealing)
.julia_helpers_state <- new.env(parent = emptyenv())
.julia_helpers_state$initialized <- FALSE

# ==============================================================================
# Initialization
# ==============================================================================

#' Initialize Julia helper functions for MCView
#'
#' Loads the Julia source file containing optimized computation functions.
#' Should be called after dafr::setup_daf() completes.
#'
#' @return TRUE if helpers were loaded successfully, FALSE otherwise
#' @export
init_julia_helpers <- function() {
    if (!julia_helpers_enabled()) {
        return(invisible(FALSE))
    }

    success <- tryCatch(
        {
            # Check if Julia is available via dafr
            if (!requireNamespace("dafr", quietly = TRUE)) {
                return(invisible(FALSE))
            }

            # dafr:::is_julia_initialized is not exported, so check via JuliaCall
            if (!requireNamespace("JuliaCall", quietly = TRUE)) {
                return(invisible(FALSE))
            }
            # Quick probe: if Julia is not initialized, this will error
            JuliaCall::julia_eval("true")

            # Load the Julia helpers source file
            julia_file <- system.file("julia", "mcview_helpers.jl", package = "MCView")
            if (!nzchar(julia_file) || !file.exists(julia_file)) {
                cli::cli_alert_warning("Julia helpers file not found at expected location")
                return(invisible(FALSE))
            }

            JuliaCall::julia_source(julia_file)

            .julia_helpers_state$initialized <- TRUE
            cli::cli_alert_success("Julia helpers initialized")
            TRUE
        },
        error = function(e) {
            cli::cli_alert_warning("Failed to initialize Julia helpers: {e$message}")
            FALSE
        }
    )

    invisible(success)
}

#' Check if Julia helpers are ready
#'
#' Returns TRUE if Julia helpers are enabled and initialized.
#' @return Logical
#' @export
julia_helpers_ready <- function() {
    julia_helpers_enabled() && .julia_helpers_state$initialized
}

#' Check if Julia helpers are enabled by option
#'
#' @return TRUE if the mcview.use_julia_helpers option is not FALSE
#' @noRd
julia_helpers_enabled <- function() {
    isTRUE(getOption("mcview.use_julia_helpers", TRUE))
}

#' Pre-warm Julia correlation cache for a dataset
#'
#' Calls mcview_calc_top_cors once to trigger JIT compilation of all helper
#' functions and build the EGC cache. This shifts ~48s of JIT overhead from
#' the first user-triggered correlation query to app startup.
#'
#' @param daf_obj DAF object
#' @param egc_epsilon Epsilon for log2 computation
#' @noRd
prewarm_julia_cors <- function(daf_obj, egc_epsilon = 1e-5) {
    if (!julia_helpers_ready()) {
        return(invisible(NULL))
    }
    tryCatch(
        {
            # Pick the first gene to force JIT + cache build
            gene <- dafr::axis_entries(daf_obj, "gene")[1]
            JuliaCall::julia_call(
                "mcview_calc_top_cors",
                daf_obj$jl_obj,
                as.character(gene),
                as.numeric(egc_epsilon),
                1L,
                need_return = "R"
            )
            cli::cli_alert_success("Julia correlation cache pre-warmed")
        },
        error = function(e) {
            cli::cli_alert_warning("Julia correlation pre-warm failed: {e$message}")
        }
    )
    invisible(NULL)
}

# ==============================================================================
# dafr dispatch warmup
# ==============================================================================

#' Warm up Julia JIT for dafr FilesDaf operations
#'
#' The first call to dafr::get_vector() on a FilesDaf triggers ~4s of Julia JIT
#' compilation for RCall's sexp()/rcopy() type conversion of NamedVector/ReadOnly
#' return types. The DataAxesFormats functions themselves are precompiled in the
#' sysimage, but RCall/JuliaCall aren't, so the R↔Julia bridge needs JIT.
#'
#' This function creates a tiny temporary FilesDaf and returns one float vector
#' and one string vector through the full julia_call → docall → sexp path.
#' This triggers JIT for the NamedVector{Float64, ReadOnly{...}} and
#' NamedVector{String, ReadOnly{...}} conversion paths.
#'
#' Cost: ~2s. Saves ~4s on every cold-start first data access.
#'
#' @noRd
prewarm_dafr_dispatch <- function() {
    tryCatch(
        {
            start <- proc.time()[["elapsed"]]
            # Create a tiny 1-element FilesDaf to minimize I/O
            JuliaCall::julia_command("using DataAxesFormats")
            JuliaCall::julia_command("_pw_dir = mktempdir()")
            JuliaCall::julia_command("_pw = FilesDaf(_pw_dir, \"w\"; name=\"w\")")
            JuliaCall::julia_command("add_axis!(_pw, \"r\", [\"a\"])")
            JuliaCall::julia_command("add_axis!(_pw, \"c\", [\"x\"])")
            JuliaCall::julia_command("set_vector!(_pw, \"r\", \"nv\", [1.0])")
            JuliaCall::julia_command("set_vector!(_pw, \"r\", \"sv\", [\"x\"])")
            JuliaCall::julia_command("set_matrix!(_pw, \"r\", \"c\", \"m\", reshape([1.0], 1, 1))")
            JuliaCall::julia_command("_pw = nothing; GC.gc()")
            JuliaCall::julia_command("_pw_r = FilesDaf(_pw_dir, \"r\"; name=\"r\")")

            # Return vectors/matrix through the full dafr R pipeline to trigger
            # JIT for: julia_call docall path, sexp for NamedVector types,
            # AND from_julia_array (_strip_wrappers, jlview, name extraction)
            pw_daf <- dafr::Daf(JuliaCall::julia_eval("_pw_r"))
            invisible(dafr::get_vector(pw_daf, "r", "nv"))  # Float64
            invisible(dafr::get_vector(pw_daf, "r", "sv"))  # String
            invisible(dafr::get_matrix(pw_daf, "r", "c", "m"))  # Float64 matrix

            # Cleanup
            JuliaCall::julia_command("_pw_r = nothing; GC.gc()")

            elapsed <- proc.time()[["elapsed"]] - start
            cli::cli_alert_success("dafr dispatch pre-warmed ({round(elapsed, 1)}s)")
        },
        error = function(e) {
            cli::cli_alert_warning("dafr dispatch pre-warm failed: {e$message}")
        }
    )
    invisible(NULL)
}

# ==============================================================================
# R Wrapper: calc_top_cors via Julia
# ==============================================================================

#' Compute top gene correlations using Julia (with R fallback)
#'
#' This wraps the Julia mcview_calc_top_cors function. If Julia helpers are
#' not available, returns NULL to signal the caller to use the R fallback.
#'
#' @param daf_obj DAF object
#' @param gene Gene name
#' @param type "pos", "neg", or "both"
#' @param egc_epsilon Epsilon for log2 computation
#' @param k Number of top correlations to return per direction
#' @param exclude Character vector of gene names to exclude (or NULL)
#'
#' @return Data frame with columns gene2, cor, type; or NULL if Julia unavailable
#' @noRd
julia_calc_top_cors <- function(daf_obj, gene, type = "both", egc_epsilon = 1e-5,
                                k = 30L, exclude = NULL) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            result <- JuliaCall::julia_call(
                "mcview_calc_top_cors",
                daf_obj$jl_obj,
                as.character(gene),
                as.numeric(egc_epsilon),
                as.integer(k),
                need_return = "R"
            )

            if (is.null(result) || length(result$genes) == 0) {
                return(NULL)
            }

            df <- tibble::tibble(
                gene2 = as.character(result$genes),
                cor = as.numeric(result$cors),
                type = as.character(result$types)
            )

            # Apply type filter
            if (type %in% c("pos", "neg")) {
                df <- df[df$type == type, , drop = FALSE]
            }

            # Apply exclusion
            if (!is.null(exclude) && length(exclude) > 0) {
                df <- df[!(df$gene2 %in% exclude), , drop = FALSE]
            }

            return(df)
        },
        error = function(e) {
            cli::cli_alert_warning("Julia calc_top_cors failed: {e$message}")
            NULL
        }
    )
}

#' Compute top correlations with a user-supplied expression vector via Julia
#'
#' @param daf_obj DAF object
#' @param data_vec Named numeric vector of expression values per metacell
#' @param type "pos", "neg", or "both"
#' @param egc_epsilon Epsilon for log2 computation
#' @param k Number of top correlations per direction
#' @param exclude Character vector of gene names to exclude
#'
#' @return Data frame with gene2, cor, type columns; or NULL if Julia unavailable
#' @noRd
julia_calc_top_cors_with_vec <- function(daf_obj, data_vec, type = "both",
                                          egc_epsilon = 1e-5, k = 30L,
                                          exclude = NULL) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            # data_vec should be ordered by metacell
            vec <- as.numeric(data_vec)
            excl <- if (is.null(exclude)) character(0) else as.character(exclude)

            result <- JuliaCall::julia_call(
                "mcview_calc_top_cors_with_vec",
                daf_obj$jl_obj,
                vec,
                as.numeric(egc_epsilon),
                as.integer(k),
                excl,
                need_return = "R"
            )

            if (is.null(result) || length(result$genes) == 0) {
                return(NULL)
            }

            df <- tibble::tibble(
                gene2 = as.character(result$genes),
                cor = as.numeric(result$cors),
                type = as.character(result$types)
            )

            if (type %in% c("pos", "neg")) {
                df <- df[df$type == type, , drop = FALSE]
            }

            return(df)
        },
        error = function(e) {
            cli::cli_alert_warning("Julia calc_top_cors_with_vec failed: {e$message}")
            NULL
        }
    )
}


# ==============================================================================
# R Wrapper: Full EGC matrix retrieval via Julia
# ==============================================================================

#' Retrieve full EGC matrix via Julia (single call, no per-gene loop)
#'
#' @param daf_obj DAF object
#' @param genes Optional character vector of gene names to retrieve
#'
#' @return Matrix (genes x metacells) or NULL if Julia unavailable
#' @noRd
julia_get_egc_matrix <- function(daf_obj, genes = NULL) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            if (!is.null(genes) && length(genes) > 0) {
                result <- JuliaCall::julia_call(
                    "mcview_get_egc_genes",
                    daf_obj$jl_obj,
                    as.character(genes),
                    need_return = "R"
                )
            } else {
                result <- JuliaCall::julia_call(
                    "mcview_get_egc_matrix",
                    daf_obj$jl_obj,
                    need_return = "R"
                )
            }

            if (is.null(result) || is.null(result$matrix)) {
                return(NULL)
            }

            mat <- result$matrix
            rownames(mat) <- as.character(result$gene_names)
            colnames(mat) <- as.character(result$metacell_names)

            return(mat)
        },
        error = function(e) {
            cli::cli_alert_warning("Julia get_egc_matrix failed: {e$message}")
            NULL
        }
    )
}


# ==============================================================================
# R Wrapper: Precompute cache reductions via Julia
# ==============================================================================

#' Precompute top genes per metacell via Julia
#'
#' @param daf_obj DAF object
#' @param egc_epsilon Epsilon for log2
#'
#' @return List with top1_gene, top2_gene, top1_lfp, top2_lfp; or NULL
#' @noRd
julia_precompute_top_genes <- function(daf_obj, egc_epsilon = 1e-5) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            result <- JuliaCall::julia_call(
                "mcview_precompute_top_genes",
                daf_obj$jl_obj,
                as.numeric(egc_epsilon),
                need_return = "R"
            )

            if (is.null(result)) {
                return(NULL)
            }

            list(
                top1_gene = as.character(result$top1_gene),
                top2_gene = as.character(result$top2_gene),
                top1_lfp = as.numeric(result$top1_lfp),
                top2_lfp = as.numeric(result$top2_lfp)
            )
        },
        error = function(e) {
            cli::cli_alert_warning("Julia precompute_top_genes failed: {e$message}")
            NULL
        }
    )
}


#' Precompute gene statistics via Julia
#'
#' @param daf_obj DAF object
#'
#' @return List with max_expr, total_expr, gene_names; or NULL
#' @noRd
julia_precompute_gene_stats <- function(daf_obj) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            result <- JuliaCall::julia_call(
                "mcview_precompute_gene_stats",
                daf_obj$jl_obj,
                need_return = "R"
            )

            if (is.null(result)) {
                return(NULL)
            }

            list(
                max_expr = stats::setNames(as.numeric(result$max_expr),
                    as.character(result$gene_names)),
                total_expr = stats::setNames(as.numeric(result$total_expr),
                    as.character(result$gene_names)),
                gene_names = as.character(result$gene_names)
            )
        },
        error = function(e) {
            cli::cli_alert_warning("Julia precompute_gene_stats failed: {e$message}")
            NULL
        }
    )
}


# ==============================================================================
# R Wrapper: calc_gg_mc_top_cor via Julia
# ==============================================================================

#' Compute top gene-gene correlations for ALL genes via Julia (with R fallback)
#'
#' This wraps the Julia mcview_calc_gg_mc_top_cor function which uses a
#' block-based BLAS approach to avoid materializing the full N×N correlation
#' matrix. If Julia helpers are not available, returns NULL to signal the
#' caller to use the R fallback.
#'
#' @param daf_obj DAF object
#' @param egc_epsilon Epsilon for log2 computation
#' @param k Number of top correlations per direction per gene
#' @param block_size Number of genes to process per BLAS block
#'
#' @return Data frame with columns gene1, gene2, cor, type; or NULL if Julia unavailable
#' @noRd
julia_calc_gg_mc_top_cor <- function(daf_obj, egc_epsilon = 1e-5, k = 30L,
                                      block_size = 1000L) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            result <- JuliaCall::julia_call(
                "mcview_calc_gg_mc_top_cor",
                daf_obj$jl_obj,
                as.numeric(egc_epsilon),
                as.integer(k),
                as.integer(block_size),
                need_return = "R"
            )

            if (is.null(result) || length(result$gene1) == 0) {
                return(NULL)
            }

            tibble::tibble(
                gene1 = as.character(result$gene1),
                gene2 = as.character(result$gene2),
                cor = as.numeric(result$cor),
                type = as.character(result$type)
            )
        },
        error = function(e) {
            cli::cli_alert_warning("Julia calc_gg_mc_top_cor failed: {e$message}")
            NULL
        }
    )
}


# ==============================================================================
# R Wrapper: calc_marker_genes via Julia
# ==============================================================================

#' Compute top marker genes per metacell via Julia (with R fallback)
#'
#' This wraps the Julia mcview_calc_marker_genes function which computes
#' log2-transform, fold-change, and per-metacell top-k extraction entirely
#' in Julia. If Julia helpers are not available, returns NULL to signal the
#' caller to use the R fallback.
#'
#' @param daf_obj DAF object
#' @param genes_per_metacell Number of top genes per metacell
#' @param minimal_max_log_fraction Gene filter threshold
#' @param minimal_relative_log_fraction Fold-change filter threshold
#' @param fold_change_reg Regularization added to fold-change
#' @param use_abs Use absolute fold-change for ranking
#'
#' @return Tibble with columns metacell, gene, rank, fp; or NULL if Julia unavailable
#' @noRd
julia_calc_marker_genes <- function(daf_obj,
                                     genes_per_metacell = 2L,
                                     minimal_max_log_fraction = -13,
                                     minimal_relative_log_fraction = 2,
                                     fold_change_reg = 0.1,
                                     use_abs = TRUE) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            result <- JuliaCall::julia_call(
                "mcview_calc_marker_genes",
                daf_obj$jl_obj,
                as.integer(genes_per_metacell),
                as.numeric(minimal_max_log_fraction),
                as.numeric(minimal_relative_log_fraction),
                as.numeric(fold_change_reg),
                as.logical(use_abs),
                need_return = "R"
            )

            if (is.null(result) || length(result$metacell) == 0) {
                return(NULL)
            }

            tibble::tibble(
                metacell = as.character(result$metacell),
                gene = as.character(result$gene),
                rank = as.integer(result$rank),
                fp = as.numeric(result$fp)
            )
        },
        error = function(e) {
            cli::cli_alert_warning("Julia calc_marker_genes failed: {e$message}")
            NULL
        }
    )
}


# ==============================================================================
# R Wrapper: Batched group QC stats via Julia
# ==============================================================================

#' Compute per-group QC stats in a single Julia call
#'
#' Batches Count, Sum, and Median of total_UMIs per group into one Julia call,
#' eliminating 3 sequential R<->Julia round-trips.
#'
#' @param cells_daf DAF object with cell axis and group_field vector
#' @param group_field Cell-level categorical field name
#'
#' @return Tibble with group_id, n_cells, total_umis, median_umis_per_cell;
#'   or NULL if Julia unavailable
#' @noRd
julia_group_qc_stats <- function(cells_daf, group_field) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            result <- JuliaCall::julia_call(
                "mcview_group_qc_stats",
                cells_daf$jl_obj,
                as.character(group_field),
                need_return = "R"
            )

            if (is.null(result) || length(result$group_id) == 0) {
                return(NULL)
            }

            tibble::tibble(
                group_id = as.character(result$group_id),
                n_cells = as.integer(result$n_cells),
                total_umis = as.numeric(result$total_umis),
                median_umis_per_cell = as.numeric(result$median_umis_per_cell)
            ) %>% dplyr::arrange(group_id)
        },
        error = function(e) {
            cli::cli_alert_warning("Julia group_qc_stats failed: {e$message}")
            NULL
        }
    )
}


# ==============================================================================
# R Wrapper: Check string vector cardinalities via Julia
# ==============================================================================

#' Filter string vectors by cardinality using Julia (with R fallback)
#'
#' Given a set of candidate vector names known to be string-typed, checks their
#' cardinality (unique value count) on the Julia side without transferring full
#' vectors to R. Returns only those with 1 < n_unique <= max_cardinality.
#'
#' @param daf_obj DAF object (must have a cell axis)
#' @param axis Axis name (e.g., "cell")
#' @param candidate_names Character vector of string-typed vector names to check
#' @param max_cardinality Maximum number of unique values (default: 1000)
#'
#' @return Character vector of field names passing the cardinality filter,
#'   or NULL if Julia unavailable (caller should use R fallback)
#' @noRd
julia_get_grouping_fields <- function(daf_obj, axis, candidate_names,
                                       max_cardinality = 1000L) {
    if (!julia_helpers_ready()) {
        return(NULL)
    }

    tryCatch(
        {
            result <- JuliaCall::julia_call(
                "mcview_get_string_vector_cardinalities",
                daf_obj$jl_obj,
                as.character(axis),
                as.character(candidate_names),
                as.integer(max_cardinality),
                need_return = "R"
            )

            if (is.null(result)) {
                return(character(0))
            }
            as.character(result)
        },
        error = function(e) {
            cli::cli_alert_warning(
                "Julia get_grouping_fields failed: {e$message}"
            )
            NULL
        }
    )
}
