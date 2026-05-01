# daf_cache_populate.R - Cache precomputation pipeline
#
# Split from R/daf_cache.R (2026-05-01). populate_dataset_cache orchestrates
# the precompute_daf_* family below. See R/daf_cache_init.R for the cache
# DAF layer this populates and R/daf_contracts_definitions.R for the cache
# contract that constrains which keys live here.


cache_correlations_daf <- function(daf_obj, gene1, df) {
    if (is.null(daf_obj) || !cache_in_daf_enabled()) {
        return(invisible(FALSE))
    }
    if (!daf_is_writable(daf_obj)) {
        return(invisible(FALSE))
    }
    if (is.null(df) || nrow(df) == 0) {
        return(invisible(FALSE))
    }
    if (!has_name(df, "gene2") || !has_name(df, "cor")) {
        return(invisible(FALSE))
    }

    if (!has_name(df, "type")) {
        df$type <- ifelse(df$cor >= 0, "pos", "neg")
    }

    new_df <- tibble(
        gene1 = as.character(gene1),
        gene2 = as.character(df$gene2),
        cor = as.numeric(df$cor),
        type = as.character(df$type)
    ) %>%
        filter(!is.na(gene1), !is.na(gene2), !is.na(cor), !is.na(type))

    if (nrow(new_df) == 0) {
        return(invisible(FALSE))
    }

    axis_name <- "mcview_cache_gg_mc_top_cor"
    existing <- NULL
    if (dafr::has_axis(daf_obj, axis_name)) {
        existing_gene1 <- daf_vec(daf_obj, axis_name, "gene1", required = FALSE)
        existing_gene2 <- daf_vec(daf_obj, axis_name, "gene2", required = FALSE)
        existing_cor <- daf_vec(daf_obj, axis_name, "cor", required = FALSE)
        existing_type <- daf_vec(daf_obj, axis_name, "type", required = FALSE)
        if (!is.null(existing_gene1) && !is.null(existing_gene2) &&
            !is.null(existing_cor) && !is.null(existing_type)) {
            existing <- tibble(
                gene1 = as.character(existing_gene1),
                gene2 = as.character(existing_gene2),
                cor = as.numeric(existing_cor),
                type = as.character(existing_type)
            )
        }
    }

    combined <- if (is.null(existing) || nrow(existing) == 0) {
        new_df
    } else {
        new_df <- anti_join(new_df, existing, by = c("gene1", "gene2", "type"))
        if (nrow(new_df) == 0) {
            return(invisible(FALSE))
        }
        bind_rows(existing, new_df)
    }

    success <- tryCatch(
        {
            dafr::add_axis(daf_obj, axis_name, as.character(seq_len(nrow(combined))), overwrite = TRUE)
            dafr::set_vector(daf_obj, axis_name, "gene1", as.character(combined$gene1))
            dafr::set_vector(daf_obj, axis_name, "gene2", as.character(combined$gene2))
            dafr::set_vector(daf_obj, axis_name, "cor", as.numeric(combined$cor))
            dafr::set_vector(daf_obj, axis_name, "type", as.character(combined$type))
            TRUE
        },
        error = function(e) {
            cli::cli_alert_warning("Cache operation failed: {e$message}")
            FALSE
        }
    )

    invisible(success)
}

precompute_daf_correlations <- function(daf_obj, k = 30, egc_epsilon = 1e-5, force = FALSE) {
    if (is.null(daf_obj) || !daf_is_writable(daf_obj)) {
        return(invisible(FALSE))
    }

    cache_axis <- "mcview_cache_gg_mc_top_cor"
    if (!force && dafr::has_axis(daf_obj, cache_axis)) {
        return(invisible(FALSE))
    }

    gg_df <- NULL
    if (dafr::has_axis(daf_obj, "gg_mc_top_cor") && !force) {
        gg_df <- convert_daf_gg_mc_top_cor(daf_obj)
    }

    if (is.null(gg_df) || nrow(gg_df) == 0) {
        mc_egc <- convert_daf_mc_egc(daf_obj)
        gg_df <- calc_gg_mc_top_cor(mc_egc, k = k, egc_epsilon = egc_epsilon, daf_obj = daf_obj)
    }

    gg_df <- gg_df %>%
        filter(!is.na(gene1), !is.na(gene2), !is.na(cor), !is.na(type))
    if (nrow(gg_df) == 0) {
        return(invisible(FALSE))
    }

    success <- tryCatch(
        {
            dafr::add_axis(daf_obj, cache_axis, as.character(seq_len(nrow(gg_df))), overwrite = TRUE)
            dafr::set_vector(daf_obj, cache_axis, "gene1", as.character(gg_df$gene1))
            dafr::set_vector(daf_obj, cache_axis, "gene2", as.character(gg_df$gene2))
            dafr::set_vector(daf_obj, cache_axis, "cor", as.numeric(gg_df$cor))
            dafr::set_vector(daf_obj, cache_axis, "type", as.character(gg_df$type))
            TRUE
        },
        error = function(e) {
            cli::cli_alert_warning("Cache operation failed: {e$message}")
            FALSE
        }
    )

    invisible(success)
}

precompute_daf_metacell_top_genes <- function(daf_obj, egc_epsilon = 1e-5, force = FALSE) {
    if (is.null(daf_obj) || !daf_is_writable(daf_obj)) {
        return(invisible(FALSE))
    }

    if (!force && dafr::has_vector(daf_obj, "metacell", "top1_gene") &&
        dafr::has_vector(daf_obj, "metacell", "top2_gene") &&
        dafr::has_vector(daf_obj, "metacell", "top1_lfp") &&
        dafr::has_vector(daf_obj, "metacell", "top2_lfp")) {
        return(invisible(FALSE))
    }

    if (!force && dafr::has_vector(daf_obj, "metacell", "mcview_cache_top1_gene") &&
        dafr::has_vector(daf_obj, "metacell", "mcview_cache_top2_gene") &&
        dafr::has_vector(daf_obj, "metacell", "mcview_cache_top1_lfp") &&
        dafr::has_vector(daf_obj, "metacell", "mcview_cache_top2_lfp")) {
        return(invisible(FALSE))
    }

    mc_egc <- convert_daf_mc_egc(daf_obj)
    if (is.null(mc_egc) || ncol(mc_egc) == 0) {
        return(invisible(FALSE))
    }

    metacell_names <- dafr::axis_entries(daf_obj, "metacell")
    gene_names <- rownames(mc_egc)

    # Top-2 genes per metacell via dafr::top_k_per_col (OpenMP, dense +
    # sparse paths). mc_egc is gene × metacell, so per-column top-K picks
    # the top genes per metacell.
    tops <- dafr::top_k_per_col(mc_egc, k = 2L)
    top1_gene <- gene_names[tops$indices[1L, ]]
    top2_gene <- gene_names[tops$indices[2L, ]]
    top1_lfp <- log2(tops$values[1L, ] + egc_epsilon)
    top2_lfp <- log2(tops$values[2L, ] + egc_epsilon)

    if (any(is.na(top1_gene)) || any(is.na(top2_gene)) ||
        any(is.na(top1_lfp)) || any(is.na(top2_lfp))) {
        return(invisible(FALSE))
    }

    success <- tryCatch(
        {
            dafr::set_vector(daf_obj, "metacell", "mcview_cache_top1_gene", as.character(top1_gene), overwrite = TRUE)
            dafr::set_vector(daf_obj, "metacell", "mcview_cache_top2_gene", as.character(top2_gene), overwrite = TRUE)
            dafr::set_vector(daf_obj, "metacell", "mcview_cache_top1_lfp", as.numeric(top1_lfp), overwrite = TRUE)
            dafr::set_vector(daf_obj, "metacell", "mcview_cache_top2_lfp", as.numeric(top2_lfp), overwrite = TRUE)
            TRUE
        },
        error = function(e) {
            cli::cli_alert_warning("Cache operation failed: {e$message}")
            FALSE
        }
    )

    invisible(success)
}

precompute_daf_gene_stats <- function(daf_obj, force = FALSE) {
    if (is.null(daf_obj) || !daf_is_writable(daf_obj)) {
        return(invisible(FALSE))
    }

    if (!force && dafr::has_vector(daf_obj, "gene", "mcview_cache_gene_max_umis") &&
        dafr::has_vector(daf_obj, "gene", "mcview_cache_gene_mean_umis") &&
        dafr::has_vector(daf_obj, "gene", "mcview_cache_gene_sum_umis")) {
        return(invisible(FALSE))
    }

    max_umis <- daf_query_gene_max_umis(daf_obj)
    mean_umis <- daf_query_gene_mean_umis(daf_obj)
    sum_umis <- daf_query_gene_sum_umis(daf_obj)

    if (any(is.na(max_umis)) || any(is.na(mean_umis)) || any(is.na(sum_umis))) {
        return(invisible(FALSE))
    }

    success <- tryCatch(
        {
            dafr::set_vector(daf_obj, "gene", "mcview_cache_gene_max_umis", as.numeric(max_umis), overwrite = TRUE)
            dafr::set_vector(daf_obj, "gene", "mcview_cache_gene_mean_umis", as.numeric(mean_umis), overwrite = TRUE)
            dafr::set_vector(daf_obj, "gene", "mcview_cache_gene_sum_umis", as.numeric(sum_umis), overwrite = TRUE)
            TRUE
        },
        error = function(e) {
            cli::cli_alert_warning("Cache operation failed: {e$message}")
            FALSE
        }
    )

    invisible(success)
}

precompute_daf_type_markers <- function(daf_obj, genes_per_type = 50,
                                       min_log_fraction = -13, force = FALSE) {
    if (is.null(daf_obj) || !daf_is_writable(daf_obj)) {
        return(invisible(FALSE))
    }

    # Check both new sparse matrix format and legacy axis format
    if (!force &&
        (dafr::has_matrix(daf_obj, "type", "gene", "mcview_marker_rank") ||
         dafr::has_axis(daf_obj, "mcview_type_markers"))) {
        return(invisible(FALSE))
    }

    # Need both EGC matrix and type annotations
    if (!dafr::has_vector(daf_obj, "metacell", "type")) {
        return(invisible(FALSE))
    }

    mc_egc <- convert_daf_mc_egc(daf_obj)
    if (is.null(mc_egc) || ncol(mc_egc) == 0) {
        return(invisible(FALSE))
    }

    mc_types <- dafr::get_vector(daf_obj, "metacell", "type")

    markers <- calc_type_marker_genes(mc_egc, mc_types,
        genes_per_type = genes_per_type,
        min_log_fraction = min_log_fraction
    )

    if (is.null(markers) || nrow(markers) == 0) {
        return(invisible(FALSE))
    }

    # Write as sparse matrices ("type","gene","mcview_marker_rank") and
    # ("type","gene","mcview_marker_fold_change")
    success <- tryCatch(
        {
            cache_type_markers_daf(daf_obj, markers)
        },
        error = function(e) {
            # Fall back to legacy axis+vectors format
            cli::cli_alert_warning("Sparse matrix write failed ({e$message}), using legacy axis format")
            cache_type_markers_daf_legacy(daf_obj, markers)
        }
    )

    invisible(success)
}

#' Write type markers as sparse matrices to DAF
#'
#' Writes ("type","gene","mcview_marker_rank") and
#' ("type","gene","mcview_marker_fold_change") sparse matrices.
#'
#' @param daf_obj Writable DAF object
#' @param markers Tibble with columns: cell_type, gene, rank, fold_change
#' @return TRUE on success, FALSE on failure
#' @noRd
cache_type_markers_daf <- function(daf_obj, markers) {
    type_names <- dafr::axis_entries(daf_obj, "type")
    gene_names <- dafr::axis_entries(daf_obj, "gene")

    # Build sparse matrices from markers tibble
    type_idx <- match(markers$cell_type, type_names)
    gene_idx <- match(markers$gene, gene_names)

    # Filter out any unmatched entries
    valid <- !is.na(type_idx) & !is.na(gene_idx)
    if (sum(valid) == 0) return(FALSE)

    type_idx <- type_idx[valid]
    gene_idx <- gene_idx[valid]
    ranks <- as.numeric(markers$rank[valid])
    folds <- as.numeric(markers$fold_change[valid])

    n_types <- length(type_names)
    n_genes <- length(gene_names)

    rank_mat <- Matrix::sparseMatrix(
        i = type_idx, j = gene_idx, x = ranks,
        dims = c(n_types, n_genes),
        dimnames = list(type_names, gene_names)
    )

    fc_mat <- Matrix::sparseMatrix(
        i = type_idx, j = gene_idx, x = folds,
        dims = c(n_types, n_genes),
        dimnames = list(type_names, gene_names)
    )

    dafr::set_matrix(daf_obj, "type", "gene", "mcview_marker_rank", rank_mat, overwrite = TRUE)
    dafr::set_matrix(daf_obj, "type", "gene", "mcview_marker_fold_change", fc_mat, overwrite = TRUE)

    TRUE
}

#' Write type markers using legacy axis+vectors format
#'
#' Fallback for DAF backends that don't support sparse matrices between
#' existing axes.
#'
#' @param daf_obj Writable DAF object
#' @param markers Tibble with columns: cell_type, gene, rank, fold_change
#' @return TRUE on success, FALSE on failure
#' @noRd
cache_type_markers_daf_legacy <- function(daf_obj, markers) {
    cache_axis <- "mcview_type_markers"
    tryCatch(
        {
            dafr::add_axis(daf_obj, cache_axis,
                as.character(seq_len(nrow(markers))),
                overwrite = TRUE
            )
            dafr::set_vector(daf_obj, cache_axis, "cell_type",
                as.character(markers$cell_type))
            dafr::set_vector(daf_obj, cache_axis, "gene",
                as.character(markers$gene))
            dafr::set_vector(daf_obj, cache_axis, "rank",
                as.numeric(markers$rank))
            dafr::set_vector(daf_obj, cache_axis, "fold_change",
                as.numeric(markers$fold_change))
            TRUE
        },
        error = function(e) {
            cli::cli_alert_warning("Cache operation failed: {e$message}")
            FALSE
        }
    )
}

precompute_daf_derived <- function(daf_obj,
                                   correlations = TRUE,
                                   metacell_top_genes = TRUE,
                                   gene_stats = TRUE,
                                   type_markers = TRUE,
                                   qc_vectors = TRUE,
                                   egc_cache = TRUE,
                                   default_markers = TRUE,
                                   correlations_k = 30,
                                   egc_epsilon = 1e-5,
                                   force = FALSE) {
    if (is.null(daf_obj) || !daf_is_writable(daf_obj)) {
        return(invisible(FALSE))
    }

    any_done <- FALSE
    if (isTRUE(correlations)) {
        any_done <- precompute_daf_correlations(daf_obj, k = correlations_k, egc_epsilon = egc_epsilon, force = force) || any_done
    }
    if (isTRUE(metacell_top_genes)) {
        any_done <- precompute_daf_metacell_top_genes(daf_obj, egc_epsilon = egc_epsilon, force = force) || any_done
    }
    if (isTRUE(gene_stats)) {
        any_done <- precompute_daf_gene_stats(daf_obj, force = force) || any_done
    }
    if (isTRUE(type_markers)) {
        any_done <- precompute_daf_type_markers(daf_obj, force = force) || any_done
    }
    if (isTRUE(qc_vectors)) {
        precompute_daf_qc_vectors(daf_obj, daf_obj, force = force)
    }
    if (isTRUE(egc_cache)) {
        precompute_daf_egc(daf_obj, daf_obj, force = force)
    }
    if (isTRUE(default_markers)) {
        precompute_daf_default_markers(daf_obj, daf_obj, force = force)
    }

    invisible(any_done)
}

# Backward-compatible alias
precompute_daf_cache <- precompute_daf_derived

# ==============================================================================
# Stub Pre-computation Functions (not yet implemented)
# ==============================================================================

#' Pre-compute max inner_fold QC vectors per metacell
#'
#' Computes mcview_cache_max_inner_fold (max inner_fold per metacell across
#' all genes) and mcview_cache_max_inner_fold_no_lateral (same but excluding
#' lateral genes).
#'
#' @param complete_daf The chained DAF for reads (derived has priority)
#' @param derived_daf The writable derived DAF layer
#' @param force Recompute even if already present
#'
#' @return Invisibly returns TRUE if vectors were written, FALSE otherwise
#' @noRd
precompute_daf_qc_vectors <- function(complete_daf, derived_daf, force = FALSE) {
    # TODO: Pre-compute max_inner_fold and max_inner_fold_no_lateral per metacell
    # Uses: daf["@ gene @ metacell :: inner_fold >> Max"] for max_inner_fold
    # Uses: filter lateral genes, then same query for max_inner_fold_no_lateral
    cli::cli_inform("Skipping QC vector pre-computation (not yet implemented)")
    invisible(FALSE)
}

#' Pre-compute geomean_fraction (EGC) when not in input DAF
#'
#' Stores the gene x metacell geomean_fraction matrix in the derived layer
#' so downstream computations can use zero-copy reads instead of computing
#' EGC on the fly.
#'
#' @param complete_daf The chained DAF for reads (derived has priority)
#' @param derived_daf The writable derived DAF layer
#' @param force Recompute even if already present
#'
#' @return Invisibly returns TRUE if matrix was written, FALSE otherwise
#' @noRd
precompute_daf_egc <- function(complete_daf, derived_daf, force = FALSE) {
    # TODO: Pre-compute geomean_fraction when not in input DAF
    # Store as gene x metacell :: mcview_cache_geomean_fraction in derived layer
    cli::cli_inform("Skipping EGC pre-computation (not yet implemented)")
    invisible(FALSE)
}

#' Pre-compute default markers list and distance matrix
#'
#' Computes the default marker gene set and a metacell x metacell distance
#' matrix based on hclust of marker correlations. Saves ~2-5s on first
#' Markers tab visit.
#'
#' @param complete_daf The chained DAF for reads (derived has priority)
#' @param derived_daf The writable derived DAF layer
#' @param force Recompute even if already present
#'
#' @return Invisibly returns TRUE if data was written, FALSE otherwise
#' @noRd
precompute_daf_default_markers <- function(complete_daf, derived_daf, force = FALSE) {
    # TODO: Pre-compute default markers list and distance matrix
    # Uses: choose_markers() to select default set, then tgs_cor + tgs_dist + hclust
    # Stores: mcview_default_markers scalar and mcview_default_markers_dist matrix
    cli::cli_inform("Skipping default markers pre-computation (not yet implemented)")
    invisible(FALSE)
}

# ==============================================================================
# External Populate API
# ==============================================================================

#' Populate MCView derived data for a DAF dataset
#'
#' This function can be called standalone (outside MCView) to precompute
#' derived data. Useful for CI/CD pipelines, batch processing, or scheduled jobs.
#'
#' @param daf_path Path to base DAF (directory for files_daf, or H5 file)
#' @param cache_dir Derived data directory path (absolute or relative to daf_path)
#' @param cache_type Storage type: "files" (persistent) or "memory" (session only)
#' @param what What to compute: NULL (all), or character vector subset of:
#'   "correlations", "metacell_top_genes", "gene_stats", "type_markers"
#' @param force Recompute even if derived data exists and is valid
#' @param verbose Print progress messages
#'
#' @return Invisibly returns the derived path (for files) or TRUE (for memory)
#'
#' @examples
#' \dontrun{
#' # Populate all derived data
#' populate_mcview_derived("/path/to/daf")
#'
#' # Populate specific items
#' populate_mcview_derived(
#'     "/path/to/daf",
#'     what = c("correlations", "gene_stats"),
#'     verbose = TRUE
#' )
#'
#' # Use custom derived directory
#' populate_mcview_derived(
#'     "/path/to/daf",
#'     cache_dir = "/shared/derived/my_dataset"
#' )
#' }
#' @export
populate_mcview_derived <- function(
    daf_path,
    cache_dir = ".mcview_derived",
    cache_type = "files",
    what = NULL,
    force = FALSE,
    verbose = TRUE) {
    # Initialize dafr if needed
    if (!requireNamespace("dafr", quietly = TRUE)) {
        cli::cli_abort("dafr package is required for cache population")
    }

    if (verbose) cli::cli_alert_info("Opening DAF: {daf_path}")

    # Open base DAF
    base_daf <- tryCatch(
        {
            if (dir.exists(daf_path)) {
                dafr::files_daf(daf_path, mode = "r")
            } else if (file.exists(daf_path) &&
                       grepl("\\.h5ad$|\\.h5$", daf_path, ignore.case = TRUE)) {
                dafr::h5ad_as_daf(daf_path)
            } else {
                dafr::open_daf(daf_path, mode = "r")
            }
        },
        error = function(e) {
            cli::cli_abort("Failed to open DAF at {daf_path}: {e$message}")
        }
    )

    # Create cache configuration
    cache_config <- create_cache_config(
        enabled = TRUE,
        type = cache_type,
        cache_dir = cache_dir,
        precompute_on_startup = FALSE
    )

    # Derive dataset name from path
    dataset_name <- basename(normalizePath(daf_path, mustWork = FALSE))
    if (dataset_name == "" || dataset_name == ".") {
        dataset_name <- "dataset"
    }

    if (verbose) cli::cli_alert_info("Initializing derived DAF (type: {cache_type})")

    # Initialize derived DAF
    derived_result <- init_derived_daf(
        base_daf = base_daf,
        dataset_name = dataset_name,
        cache_config = cache_config,
        base_path = daf_path
    )

    if (is.null(derived_result$cache_daf)) {
        cli::cli_abort("Failed to create derived DAF")
    }

    derived_daf <- derived_result$cache_daf
    # Precompute functions read source matrices (UMIs, total_UMIs) and write
    # derived vectors back. The complete chain satisfies both — reads fall
    # through to base, writes land in derived via chain_writer. Passing the
    # bare derived_daf aborts on missing UMIs.
    chain_daf <- derived_result$complete_daf %||% derived_daf

    # Determine what to compute
    all_items <- c("correlations", "metacell_top_genes", "gene_stats", "type_markers")
    if (is.null(what)) {
        what <- all_items
    } else {
        what <- intersect(what, all_items)
        if (length(what) == 0) {
            cli::cli_abort("No valid items specified. Choose from: {paste(all_items, collapse=', ')}")
        }
    }

    if (verbose) cli::cli_alert_info("Populating derived data: {paste(what, collapse=', ')}")

    # Populate derived data
    results <- list()

    if ("correlations" %in% what) {
        if (verbose) cli::cli_alert("Computing gene correlations...")
        results$correlations <- precompute_daf_correlations(chain_daf, force = force)
        if (verbose && results$correlations) cli::cli_alert_success("Gene correlations computed")
    }

    if ("metacell_top_genes" %in% what) {
        if (verbose) cli::cli_alert("Computing metacell top genes...")
        results$metacell_top_genes <- precompute_daf_metacell_top_genes(chain_daf, force = force)
        if (verbose && results$metacell_top_genes) cli::cli_alert_success("Metacell top genes computed")
    }

    if ("gene_stats" %in% what) {
        if (verbose) cli::cli_alert("Computing gene statistics...")
        results$gene_stats <- precompute_daf_gene_stats(chain_daf, force = force)
        if (verbose && results$gene_stats) cli::cli_alert_success("Gene statistics computed")
    }

    if ("type_markers" %in% what) {
        if (verbose) cli::cli_alert("Computing per-type marker genes...")
        results$type_markers <- precompute_daf_type_markers(chain_daf, force = force)
        if (verbose && results$type_markers) cli::cli_alert_success("Per-type marker genes computed")
    }

    # Update derived metadata
    update_cache_metadata(derived_daf, base_daf)

    if (verbose) {
        n_done <- sum(unlist(results))
        if (n_done > 0) {
            cli::cli_alert_success("Derived data population complete: {n_done} items computed")
        } else {
            cli::cli_alert_info("Derived data is already up to date (use force=TRUE to recompute)")
        }
    }

    invisible(derived_result$cache_path %||% TRUE)
}

#' @rdname populate_mcview_derived
#' @export
populate_mcview_cache <- populate_mcview_derived

#' Populate derived data for a dataset in a running MCView session
#'
#' @param dataset Dataset name
#' @param what What to compute (NULL = all)
#' @param force Recompute even if present
#' @param verbose Print progress
#'
#' @export
populate_dataset_cache <- function(dataset, what = NULL, force = FALSE, verbose = TRUE) {
    derived_daf <- get_cache_daf(dataset)
    base_daf <- get_base_daf(dataset)
    # `get_dataset_daf` returns the complete chain (base + derived) — reads
    # hit base (UMIs, total_UMIs), writes land in the derived layer.
    chain_daf <- get_dataset_daf(dataset) %||% derived_daf

    if (is.null(derived_daf)) {
        if (verbose) cli::cli_alert_warning("No derived DAF available for dataset: {dataset}")
        return(invisible(FALSE))
    }

    # Determine what to compute
    all_items <- c("correlations", "metacell_top_genes", "gene_stats", "type_markers")
    if (is.null(what)) {
        what <- all_items
    }

    results <- list()

    if ("correlations" %in% what) {
        if (verbose) cli::cli_alert("Computing gene correlations for {dataset}...")
        results$correlations <- precompute_daf_correlations(chain_daf, force = force)
    }

    if ("metacell_top_genes" %in% what) {
        if (verbose) cli::cli_alert("Computing metacell top genes for {dataset}...")
        results$metacell_top_genes <- precompute_daf_metacell_top_genes(chain_daf, force = force)
    }

    if ("gene_stats" %in% what) {
        if (verbose) cli::cli_alert("Computing gene statistics for {dataset}...")
        results$gene_stats <- precompute_daf_gene_stats(chain_daf, force = force)
    }

    if ("type_markers" %in% what) {
        if (verbose) cli::cli_alert("Computing per-type marker genes for {dataset}...")
        results$type_markers <- precompute_daf_type_markers(chain_daf, force = force)
    }

    # Update derived metadata
    if (!is.null(base_daf)) {
        update_cache_metadata(derived_daf, base_daf)
    }

    if (verbose) {
        n_done <- sum(unlist(results))
        cli::cli_alert_success("Derived data population for {dataset}: {n_done} items computed")
    }

    invisible(any(unlist(results)))
}
