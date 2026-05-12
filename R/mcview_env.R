#' MCView Global State Management
#'
#' Centralized environment for all MCView global state to improve
#' isolation, testing, and maintainability.
#'
#' MCView uses DAF (Data Access Framework) as its data backend.
#'
#' Dataset structure (mc_data[[dataset_name]]):
#'   - base_daf: Original DAF object (read-only or writable)
#'   - cache_daf: Cache DAF object (writable, for computed data)
#'   - daf_obj: Complete DAF (cache + base chained via complete_daf)
#'   - top_cor_genes: In-memory correlation cache
#'   - cache_path: Path to cache directory (for files cache)
#'   - needs_population: Whether cache needs to be populated

# Create the global environment
mcview_env <- new.env(parent = emptyenv())

#' Initialize MCView environment with default values
#'
#' @return TRUE invisibly
#' @export
init_mcview_env <- function() {
    mcview_env$config <- NULL
    mcview_env$mc_data <- NULL # Named list of datasets (see structure above)
    mcview_env$atlas <- NULL # Atlas DAF object (optional)
    mcview_env$tab_defs <- NULL
    mcview_env$about_file <- NULL
    mcview_env$about_markdown <- NULL

    invisible(TRUE)
}

#' Get value from MCView environment
#'
#' @param var_name Name of the variable
#' @return The value stored in the environment
#' @export
mcv_get <- function(var_name) {
    if (!exists(var_name, envir = mcview_env, inherits = FALSE)) {
        return(NULL)
    }
    get(var_name, envir = mcview_env, inherits = FALSE)
}

#' Set value in MCView environment
#'
#' @param var_name Name of the variable
#' @param value Value to set
#' @return The value, invisibly
#' @export
mcv_set <- function(var_name, value) {
    assign(var_name, value, envir = mcview_env)
    invisible(value)
}

#' Check if variable exists in MCView environment
#'
#' @param var_name Name of the variable
#' @return TRUE if variable exists
#' @export
mcv_exists <- function(var_name) {
    exists(var_name, envir = mcview_env, inherits = FALSE)
}

#' Get value from mc_data cache for a dataset
#'
#' Reads directly from the mcview_env$mc_data environment slot,
#' avoiding the copy-on-modify overhead of mcv_get/mcv_set round-trips.
#'
#' @param dataset Dataset name
#' @param key Cache key within the dataset
#' @param default Value to return if not found (default NULL)
#' @return Cached value or default
#' @export
mcv_cache_get <- function(dataset, key, default = NULL) {
    mc_data <- mcview_env$mc_data
    if (is.null(mc_data[[dataset]])) {
        return(default)
    }
    mc_data[[dataset]][[key]] %||% default
}

#' Set value in mc_data cache for a dataset
#'
#' Writes directly into the mcview_env$mc_data environment slot,
#' avoiding the copy-on-modify overhead of mcv_get/mcv_set round-trips.
#'
#' @param dataset Dataset name
#' @param key Cache key within the dataset
#' @param value Value to store
#' @return The value, invisibly
#' @export
mcv_cache_set <- function(dataset, key, value) {
    if (is.null(mcview_env$mc_data[[dataset]])) {
        mcview_env$mc_data[[dataset]] <- list()
    }
    mcview_env$mc_data[[dataset]][[key]] <- value
    invisible(value)
}

#' Clean up MCView environment for orderly shutdown
#'
#' Iterates over all datasets, clears in-memory caches (top_cor_genes),
#' and resets the mcview_env to a clean state. All errors are suppressed
#' since this runs during shutdown when partial teardown is expected.
#'
#' @return TRUE invisibly
#' @export
cleanup_mcview_env <- function() {
    # Clear per-dataset in-memory caches and DAF internal caches
    mc_data <- mcv_get("mc_data")
    if (!is.null(mc_data)) {
        for (ds_name in names(mc_data)) {
            tryCatch(
                {
                    # Clear in-memory correlation cache
                    if (!is.null(mc_data[[ds_name]][["top_cor_genes"]])) {
                        mc_data[[ds_name]][["top_cor_genes"]] <- list()
                    }
                    # Empty DAF internal caches
                    for (daf_field in c("daf_obj", "cache_daf", "base_daf", "cells_daf", "chained_cells_daf")) {
                        daf_ref <- mc_data[[ds_name]][[daf_field]]
                        if (!is.null(daf_ref) && dafr::is_daf(daf_ref)) {
                            tryCatch(
                                dafr::empty_cache(daf_ref),
                                error = function(e) NULL
                            )
                        }
                    }
                    # Null out DAF references to allow GC
                    mc_data[[ds_name]][["daf_obj"]] <- NULL
                    mc_data[[ds_name]][["cache_daf"]] <- NULL
                    mc_data[[ds_name]][["base_daf"]] <- NULL
                    mc_data[[ds_name]][["cells_daf"]] <- NULL
                    mc_data[[ds_name]][["chained_cells_daf"]] <- NULL
                },
                error = function(e) {
                    # Suppress errors during cleanup
                }
            )
        }
    }

    # Clear atlas DAF cache
    atlas <- mcv_get("atlas")
    if (!is.null(atlas) && dafr::is_daf(atlas)) {
        tryCatch(
            dafr::empty_cache(atlas),
            error = function(e) NULL
        )
    }

    # Reset the environment to clean state
    tryCatch(
        init_mcview_env(),
        error = function(e) {
            # Suppress errors during cleanup
        }
    )

    invisible(TRUE)
}

#' Get DAF object for a dataset
#'
#' Returns the complete DAF (cache + base chained), which provides
#' unified read access with cache data taking priority.
#'
#' For more granular access, see:
#' - get_cache_daf(): Get writable cache DAF only
#' - get_base_daf(): Get original base DAF only
#' - get_complete_daf(): Alias for this function
#'
#' @param dataset Name of the dataset
#' @return Complete DAF object or NULL if not found
#' @export
get_dataset_daf <- function(dataset) {
    mc_data <- mcv_get("mc_data")
    if (is.null(mc_data) || is.null(mc_data[[dataset]])) {
        return(NULL)
    }
    mc_data[[dataset]]$daf_obj
}

#' Get atlas DAF object
#'
#' @return Atlas DAF object or NULL if not set
#' @export
get_atlas_daf <- function() {
    mcv_get("atlas")
}

#' List available datasets
#'
#' @return Character vector of dataset names
#' @export
dataset_names <- function() {
    mc_data <- mcv_get("mc_data")
    if (is.null(mc_data)) {
        return(character(0))
    }
    names(mc_data)
}

#' Validate that the MCView environment is fully populated post-init
#'
#' Called by `run_app()` after the init chain (init_single_daf_mode /
#' init_multi_daf_mode -> init_atlas -> init_defs). Treats partial
#' initialization as a hard error so modules don't have to defensively
#' check for NULL on every read.
#'
#' Required keys: `config`, `mc_data`, `tab_defs`, `egc_epsilon`,
#' `expr_breaks`, `default_gene1`, `default_gene2`. Each entry in
#' `mc_data` must be a list with at least `base_daf` and `daf_obj` as
#' valid DAF objects. `atlas` is optional; if present it must be a DAF
#' object.
#'
#' @return TRUE invisibly on success. Calls `cli::cli_abort()` on failure.
#' @export
validate_mcview_env <- function() {
    errors <- character(0)

    required_keys <- c(
        "config", "mc_data", "tab_defs",
        "egc_epsilon", "expr_breaks",
        "default_gene1", "default_gene2"
    )
    for (key in required_keys) {
        if (is.null(mcv_get(key))) {
            errors <- c(errors, sprintf("required key {.val %s} is missing or NULL", key))
        }
    }

    mc_data <- mcv_get("mc_data")
    if (!is.null(mc_data)) {
        if (!is.list(mc_data) || length(mc_data) == 0) {
            errors <- c(errors, "{.var mc_data} must be a non-empty list")
        } else if (is.null(names(mc_data)) || any(!nzchar(names(mc_data)))) {
            errors <- c(errors, "{.var mc_data} entries must all be named")
        } else {
            for (ds_name in names(mc_data)) {
                ds <- mc_data[[ds_name]]
                if (!is.list(ds)) {
                    errors <- c(errors, sprintf("{.code mc_data[['%s']]} must be a list", ds_name))
                    next
                }
                for (slot in c("base_daf", "daf_obj")) {
                    val <- ds[[slot]]
                    if (is.null(val)) {
                        errors <- c(errors, sprintf("{.code mc_data[['%s']][['%s']]} is missing", ds_name, slot))
                    } else if (!dafr::is_daf(val)) {
                        errors <- c(errors, sprintf("{.code mc_data[['%s']][['%s']]} is not a DAF object", ds_name, slot))
                    }
                }
                if (!is.null(ds$cache_daf) && !dafr::is_daf(ds$cache_daf)) {
                    errors <- c(errors, sprintf("{.code mc_data[['%s']]$cache_daf} is set but not a DAF object", ds_name))
                }
            }
        }
    }

    tab_defs <- mcv_get("tab_defs")
    if (!is.null(tab_defs) && (!is.list(tab_defs) || length(tab_defs) == 0)) {
        errors <- c(errors, "{.var tab_defs} must be a non-empty list")
    }

    atlas <- mcv_get("atlas")
    if (!is.null(atlas) && !dafr::is_daf(atlas)) {
        errors <- c(errors, "{.var atlas} is set but not a DAF object")
    }

    if (length(errors) > 0) {
        names(errors) <- rep_len("x", length(errors))
        cli::cli_abort(c("MCView environment validation failed:", errors))
    }

    invisible(TRUE)
}

# Initialize on package load
init_mcview_env()
