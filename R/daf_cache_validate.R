# daf_cache_validate.R - Validate the cache against the cache contract
#
# Split from R/daf_cache.R (2026-05-01). Calls into the validators in
# R/daf_contracts_validators.R using the cache contract from
# R/daf_contracts_definitions.R.


#' Validate MCView cache
#'
#' Check if cache is valid and complete for a DAF dataset.
#'
#' @param daf_path Path to base DAF
#' @param cache_dir Cache directory path (or NULL to auto-detect)
#'
#' @return List with validation results:
#'   - valid: TRUE/FALSE
#'   - complete: TRUE/FALSE (has all expected data)
#'   - stale_reason: character (if invalid)
#'   - missing: character vector of missing cache items
#'   - cache_path: resolved cache path
#'
#' @export
validate_mcview_cache <- function(daf_path, cache_dir = NULL) {
    # Resolve derived path (try new name first, fall back to legacy)
    if (is.null(cache_dir)) {
        cache_dir <- ".mcview_derived"
    }

    dataset_name <- basename(normalizePath(daf_path, mustWork = FALSE))
    cache_path <- resolve_cache_path(cache_dir, daf_path, dataset_name)

    result <- list(
        valid = FALSE,
        complete = FALSE,
        stale_reason = NULL,
        missing = character(),
        cache_path = cache_path
    )

    # Check if cache exists
    if (!dir.exists(cache_path) || !file.exists(fs::path(cache_path, "daf.json"))) {
        result$stale_reason <- "Cache directory does not exist"
        return(result)
    }

    # Open cache and base DAFs
    cache_daf <- tryCatch(dafr::files_daf(cache_path, mode = "r"), error = function(e) NULL)
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
        error = function(e) NULL
    )

    if (is.null(cache_daf)) {
        result$stale_reason <- "Failed to open cache DAF"
        return(result)
    }

    if (is.null(base_daf)) {
        result$stale_reason <- "Failed to open base DAF"
        return(result)
    }

    # Check validity
    if (is_cache_valid(cache_daf, base_daf, list(strategy = "hash"))) {
        result$valid <- TRUE
    } else {
        result$stale_reason <- "Cache hash does not match base DAF"
    }

    # Check completeness
    expected_items <- list(
        correlations = dafr::has_axis(cache_daf, "mcview_cache_gg_mc_top_cor"),
        metacell_top_genes = dafr::has_vector(cache_daf, "metacell", "mcview_cache_top1_gene"),
        gene_stats = dafr::has_vector(cache_daf, "gene", "mcview_cache_gene_max_umis"),
        type_markers = dafr::has_axis(cache_daf, "mcview_type_markers")
    )

    result$missing <- names(expected_items)[!unlist(expected_items)]
    result$complete <- length(result$missing) == 0

    result
}
