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
                    # Empty DAF internal caches (releases Julia-side memory)
                    for (daf_field in c("daf_obj", "cache_daf", "base_daf", "cells_daf", "chained_cells_daf")) {
                        daf_ref <- mc_data[[ds_name]][[daf_field]]
                        if (!is.null(daf_ref) && inherits(daf_ref, "Daf")) {
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
    if (!is.null(atlas) && inherits(atlas, "Daf")) {
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

# Initialize on package load
init_mcview_env()
