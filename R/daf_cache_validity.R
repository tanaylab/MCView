# daf_cache_validity.R - Cache validity, invalidation, and accessors
#
# Split from R/daf_cache.R (2026-05-01). See R/daf_cache_init.R for cache
# construction.

# ==============================================================================
# Cache Validation
# ==============================================================================

#' Check if cache is valid (not stale)
#'
#' @param cache_daf Cache DAF object
#' @param base_daf Base DAF object
#' @param invalidation Invalidation config (list with strategy)
#'
#' @return TRUE if cache is valid, FALSE if stale or missing
#' @export
is_cache_valid <- function(cache_daf, base_daf, invalidation = NULL) {
    if (is.null(cache_daf)) {
        return(FALSE)
    }

    strategy <- invalidation$strategy %||% "version"

    switch(strategy,
        "version" = {
            base_version <- daf_scalar(base_daf, "mcview_data_version", default = NULL)
            cache_base_version <- daf_scalar(cache_daf, "mcview_cache_base_version", default = NULL)

            # If base has no version, fall back to hash
            if (is.null(base_version)) {
                return(is_cache_valid(cache_daf, base_daf, list(strategy = "hash")))
            }

            !is.null(cache_base_version) && base_version == cache_base_version
        },
        "hash" = {
            current_hash <- compute_base_hash(base_daf)
            cached_hash <- daf_scalar(cache_daf, "mcview_cache_base_hash", default = "")
            nzchar(cached_hash) && current_hash == cached_hash
        },
        "mtime" = {
            base_mtime <- daf_scalar(base_daf, "mcview_data_mtime", default = NULL)
            cache_created <- daf_scalar(cache_daf, "mcview_cache_created", default = 0)
            !is.null(base_mtime) && cache_created >= base_mtime
        },
        "manual" = {
            # Always valid unless explicitly invalidated
            daf_scalar(cache_daf, "mcview_cache_valid", default = TRUE)
        },
        FALSE # Unknown strategy = invalid
    )
}

#' Invalidate cache for a dataset
#'
#' Clears cache metadata to force repopulation on next access.
#'
#' @param dataset Dataset name
#' @param clear_data If TRUE, also clear cached vectors/axes
#'
#' @export
invalidate_cache <- function(dataset, clear_data = FALSE) {
    cache_daf <- get_cache_daf(dataset)
    if (is.null(cache_daf)) {
        return(invisible(FALSE))
    }

    # Clear validation scalars
    tryCatch(
        {
            if (dafr::has_scalar(cache_daf, "mcview_cache_base_version")) {
                dafr::delete_scalar(cache_daf, "mcview_cache_base_version")
            }
            if (dafr::has_scalar(cache_daf, "mcview_cache_base_hash")) {
                dafr::delete_scalar(cache_daf, "mcview_cache_base_hash")
            }
            if (dafr::has_scalar(cache_daf, "mcview_cache_created")) {
                dafr::delete_scalar(cache_daf, "mcview_cache_created")
            }
        },
        error = function(e) NULL
    )

    if (clear_data) {
        # Clear all mcview_cache_* vectors
        for (axis_name in c("metacell", "gene")) {
            if (dafr::has_axis(cache_daf, axis_name)) {
                vectors <- tryCatch(dafr::vectors_set(cache_daf, axis_name), error = function(e) character())
                for (vec_name in vectors[startsWith(vectors, "mcview_cache_")]) {
                    tryCatch(dafr::delete_vector(cache_daf, axis_name, vec_name), error = function(e) NULL)
                }
            }
        }

        # Clear correlation axis
        if (dafr::has_axis(cache_daf, "mcview_cache_gg_mc_top_cor")) {
            tryCatch(dafr::delete_axis(cache_daf, "mcview_cache_gg_mc_top_cor"), error = function(e) NULL)
        }
    }

    cli::cli_alert_info("Derived data invalidated for dataset: {dataset}")
    invisible(TRUE)
}

#' Update cache metadata after population
#'
#' @param cache_daf Cache DAF
#' @param base_daf Base DAF
update_cache_metadata <- function(cache_daf, base_daf) {
    if (is.null(cache_daf)) {
        return(invisible(FALSE))
    }

    # Update creation timestamp
    tryCatch(
        dafr::set_scalar(cache_daf, "mcview_cache_created", as.numeric(Sys.time()), overwrite = TRUE),
        error = function(e) {
            cli::cli_alert_warning("Cache operation failed: {e$message}")
            NULL
        }
    )

    # Update base version
    base_version <- daf_scalar(base_daf, "mcview_data_version", default = NULL)
    if (!is.null(base_version)) {
        tryCatch(
            dafr::set_scalar(cache_daf, "mcview_cache_base_version", base_version, overwrite = TRUE),
            error = function(e) {
                cli::cli_alert_warning("Cache operation failed: {e$message}")
                NULL
            }
        )
    }

    # Update base hash
    base_hash <- compute_base_hash(base_daf)
    tryCatch(
        dafr::set_scalar(cache_daf, "mcview_cache_base_hash", base_hash, overwrite = TRUE),
        error = function(e) {
            cli::cli_alert_warning("Cache operation failed: {e$message}")
            NULL
        }
    )

    invisible(TRUE)
}

# ==============================================================================
# Cache DAF Accessors
# ==============================================================================

#' Get cache DAF for a dataset
#'
#' @param dataset Dataset name
#'
#' @return Cache DAF object (writable) or NULL
#' @export
get_cache_daf <- function(dataset) {
    mc_data <- mcv_get("mc_data")
    if (is.null(mc_data) || is.null(mc_data[[dataset]])) {
        return(NULL)
    }
    mc_data[[dataset]]$cache_daf
}

#' Get complete DAF for a dataset (cache + base chained)
#'
#' @param dataset Dataset name
#'
#' @return Complete DAF object for reading (cache has priority)
#' @export
get_complete_daf <- function(dataset) {
    mc_data <- mcv_get("mc_data")
    if (is.null(mc_data) || is.null(mc_data[[dataset]])) {
        return(NULL)
    }
    # complete_daf is stored as daf_obj for backward compatibility
    mc_data[[dataset]]$daf_obj
}

#' Get base DAF for a dataset (without cache)
#'
#' @param dataset Dataset name
#'
#' @return Base DAF object or NULL
#' @export
get_base_daf <- function(dataset) {
    mc_data <- mcv_get("mc_data")
    if (is.null(mc_data) || is.null(mc_data[[dataset]])) {
        return(NULL)
    }
    mc_data[[dataset]]$base_daf
}

# ==============================================================================
# Legacy Helper Functions
# ==============================================================================

normalize_cache_flag <- function(value) {
    if (is.null(value)) {
        return(NULL)
    }
    if (is.logical(value) && length(value) == 1 && !is.na(value)) {
        return(isTRUE(value))
    }
    if (is.numeric(value) && length(value) == 1 && !is.na(value)) {
        return(value != 0)
    }
    if (is.character(value) && length(value) == 1) {
        norm <- tolower(trimws(value))
        if (norm %in% c("true", "t", "yes", "y", "1")) {
            return(TRUE)
        }
        if (norm %in% c("false", "f", "no", "n", "0")) {
            return(FALSE)
        }
    }
    NULL
}

daf_description_fields <- function(daf_obj, deep = TRUE) {
    desc <- tryCatch(dafr::description(daf_obj, deep = deep), error = function(e) NULL)
    if (is.null(desc)) {
        return(list())
    }
    lines <- strsplit(desc, "\n", fixed = TRUE)[[1]]
    fields <- list()
    for (line in lines) {
        matches <- regexec("^\\s*([A-Za-z0-9_]+):\\s*(.*)$", line)
        parts <- regmatches(line, matches)[[1]]
        if (length(parts) == 3) {
            fields[[parts[2]]] <- parts[3]
        }
    }
    fields
}

daf_is_writable <- function(daf_obj) {
    # Native dafr tags every writable DAF with the `dafr::DafWriter` S7 class;
    # read-only wrappers (`dafr::read_only(...)`) strip it. Inherit check is
    # both cheaper and more reliable than parsing the description text, which
    # omits a `mode:` line for FilesDaf and would fall through to FALSE.
    if (inherits(daf_obj, "dafr::DafWriter")) {
        return(TRUE)
    }
    # Fall back to description parsing for any non-native / legacy DAF object
    # that still exposes `mode` / `type` fields.
    fields <- daf_description_fields(daf_obj, deep = TRUE)
    if (!is.null(fields$mode)) {
        return(!identical(fields$mode, "r"))
    }
    if (!is.null(fields$type)) {
        if (grepl("ReadOnly", fields$type, fixed = TRUE)) {
            return(FALSE)
        }
        if (grepl("Write", fields$type, fixed = TRUE) || grepl("Writer", fields$type, fixed = TRUE)) {
            return(TRUE)
        }
    }
    FALSE
}

daf_storage_path <- function(daf_obj) {
    path <- tryCatch(dafr::complete_path(daf_obj), error = function(e) NULL)
    if (!is.null(path) && nzchar(path)) {
        return(path)
    }
    NULL
}

cache_in_daf_enabled <- function(config = mcv_get("config")) {
    isTRUE(config$cache_in_daf)
}
