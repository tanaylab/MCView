# daf_cache_init.R - Construction of the mcview_derived cache layer
#
# Split from R/daf_cache.R (2026-05-01). Companions:
#   - R/daf_cache_validity.R   - is_cache_valid, invalidate_cache, getters
#   - R/daf_cache_populate.R   - precompute_daf_*, populate_dataset_cache
#   - R/daf_cache_validate.R   - validate_mcview_cache (cache vs. contract)


# ==============================================================================
# Derived DAF Initialization (complete_daf based)
# ==============================================================================

#' Initialize derived DAF for a dataset
#'
#' Creates or opens a derived DAF and chains it with the base DAF using
#' complete_daf. The derived DAF stores computed data that can be persisted
#' across sessions.
#'
#' @param base_daf Base DAF object (read-only or writable)
#' @param dataset_name Dataset name (used for derived directory)
#' @param cache_config Cache configuration from extract_cache_config()
#' @param base_path Path to base DAF (for relative path calculation)
#'
#' @return List with:
#'   - cache_daf: The writable derived DAF (or NULL if disabled)
#'   - complete_daf: The chained DAF for reads (derived has priority)
#'   - needs_population: Whether derived data needs to be populated
#'   - cache_path: Path to derived directory (for files mode)
#'
#' @export
init_derived_daf <- function(base_daf, dataset_name, cache_config, base_path = NULL) {
    # If disabled, return base DAF as-is
    if (!isTRUE(cache_config$enabled)) {
        return(list(
            cache_daf = NULL,
            complete_daf = base_daf,
            needs_population = FALSE,
            cache_path = NULL
        ))
    }

    # Get base path if not provided
    if (is.null(base_path)) {
        base_path <- daf_storage_path(base_daf)
    }

    # Determine derived path (accept legacy ".mcview_cache" for backward compat)
    derived_dir <- cache_config$cache_dir %||% ".mcview_derived"
    derived_path <- resolve_cache_path(derived_dir, base_path, dataset_name)

    # Create derived DAF based on type
    derived_daf <- NULL
    if (cache_config$type == "memory") {
        derived_daf <- create_memory_cache_daf(base_daf, dataset_name)
    } else {
        derived_daf <- create_files_cache_daf(base_daf, derived_path, base_path)
    }

    if (is.null(derived_daf)) {
        cli::cli_alert_warning("Failed to create derived DAF for {dataset_name}, using base DAF only")
        return(list(
            cache_daf = NULL,
            complete_daf = base_daf,
            needs_population = FALSE,
            cache_path = NULL
        ))
    }

    # Check if derived data is valid
    needs_population <- !is_cache_valid(derived_daf, base_daf, cache_config$invalidation)

    # Create complete DAF chain (derived as leaf, base as parent)
    complete_daf <- NULL
    if (cache_config$type == "files" && !is.null(derived_path)) {
        # For files mode, use complete_daf which follows base_daf_repository
        complete_daf <- tryCatch(
            dafr::complete_daf(derived_path, mode = "r+"),
            error = function(e) {
                cli::cli_alert_warning("Failed to create complete DAF: {e$message}")
                NULL
            }
        )
    }

    # Fallback to chain_writer if complete_daf failed or memory mode
    if (is.null(complete_daf)) {
        complete_daf <- tryCatch(
            dafr::chain_writer(list(dafr::read_only(base_daf), derived_daf)),
            error = function(e) {
                cli::cli_alert_warning("Failed to chain DAFs: {e$message}")
                base_daf
            }
        )
    }

    label <- if (cache_config$type == "memory") "in-memory" else derived_path
    cli::cli_alert_info("Derived DAF initialized for {dataset_name}: {label}")

    list(
        cache_daf = derived_daf,
        complete_daf = complete_daf,
        needs_population = needs_population,
        cache_path = if (cache_config$type == "files") derived_path else NULL
    )
}

#' @rdname init_derived_daf
#' @export
init_cache_daf <- init_derived_daf

#' Resolve cache path from configuration
#'
#' @param cache_dir Cache directory from config (relative or absolute)
#' @param base_path Path to base DAF
#' @param dataset_name Dataset name
#'
#' @return Absolute path to cache directory
resolve_cache_path <- function(cache_dir, base_path, dataset_name) {
    if (is.null(cache_dir) || !nzchar(cache_dir)) {
        cache_dir <- ".mcview_derived"
    }

    # If cache_dir is absolute, use it directly
    if (fs::is_absolute_path(cache_dir)) {
        return(fs::path(cache_dir, dataset_name))
    }

    # Otherwise, resolve relative to base DAF path
    if (!is.null(base_path) && nzchar(base_path)) {
        base_dir <- if (dir.exists(base_path)) base_path else dirname(base_path)
        return(normalizePath(fs::path(base_dir, cache_dir, dataset_name), mustWork = FALSE))
    }

    # Fallback to working directory
    normalizePath(fs::path(cache_dir, dataset_name), mustWork = FALSE)
}

#' Create in-memory cache DAF
#'
#' @param base_daf Base DAF to copy axes from
#' @param name Name for the memory DAF
#'
#' @return Memory DAF object with axes copied from base
create_memory_cache_daf <- function(base_daf, name = "cache") {
    cache_daf <- tryCatch(
        dafr::memory_daf(name = paste0(name, "_cache")),
        error = function(e) {
            cli::cli_alert_warning("Cache operation failed: {e$message}")
            NULL
        }
    )

    if (is.null(cache_daf)) {
        return(NULL)
    }

    # Copy axes from base DAF
    copy_axes_to_cache(base_daf, cache_daf)

    cache_daf
}

#' Create files-based cache DAF
#'
#' @param base_daf Base DAF to reference
#' @param cache_path Path to cache directory
#' @param base_path Path to base DAF (for relative path)
#'
#' @return Files DAF object with base_daf_repository set
create_files_cache_daf <- function(base_daf, cache_path, base_path) {
    if (is.null(cache_path) || !nzchar(cache_path)) {
        return(NULL)
    }

    cache_path <- normalizePath(cache_path, mustWork = FALSE)

    # Check if cache already exists
    cache_daf <- NULL
    if (dir.exists(cache_path) && file.exists(fs::path(cache_path, "daf.json"))) {
        cache_daf <- tryCatch(
            dafr::files_daf(cache_path, mode = "r+"),
            error = function(e) {
                cli::cli_alert_warning("Cache operation failed: {e$message}")
                NULL
            }
        )
    } else {
        # Create new cache directory
        fs::dir_create(cache_path, recurse = TRUE)
        cache_daf <- tryCatch(
            dafr::files_daf(cache_path, mode = "w+"),
            error = function(e) {
                cli::cli_alert_warning("Cache operation failed: {e$message}")
                NULL
            }
        )
    }

    if (is.null(cache_daf)) {
        return(NULL)
    }

    # Set base_daf_repository scalar with relative path
    if (!dafr::has_scalar(cache_daf, "base_daf_repository")) {
        if (!is.null(base_path) && nzchar(base_path)) {
            rel_path <- tryCatch(
                fs::path_rel(base_path, start = cache_path),
                error = function(e) base_path
            )
            tryCatch(
                dafr::set_scalar(cache_daf, "base_daf_repository", as.character(rel_path)),
                error = function(e) {
                    cli::cli_alert_warning("Cache operation failed: {e$message}")
                    NULL
                }
            )
        }
    }

    # Copy axes from base DAF if not present
    copy_axes_to_cache(base_daf, cache_daf)

    # Set cache metadata
    set_cache_metadata(cache_daf, base_daf)

    cache_daf
}

#' Copy axes from base DAF to cache DAF
#'
#' @param base_daf Source DAF
#' @param cache_daf Target cache DAF
copy_axes_to_cache <- function(base_daf, cache_daf) {
    required_axes <- c("metacell", "gene", "type")

    for (axis_name in required_axes) {
        if (!dafr::has_axis(cache_daf, axis_name) && dafr::has_axis(base_daf, axis_name)) {
            entries <- tryCatch(
                dafr::axis_entries(base_daf, axis_name),
                error = function(e) {
                    cli::cli_alert_warning("Cache operation failed: {e$message}")
                    NULL
                }
            )
            if (!is.null(entries)) {
                tryCatch(
                    dafr::add_axis(cache_daf, axis_name, entries),
                    error = function(e) {
                        cli::cli_alert_warning("Cache operation failed: {e$message}")
                        NULL
                    }
                )
            }
        }
    }
}

#' Set cache metadata scalars
#'
#' @param cache_daf Cache DAF to update
#' @param base_daf Base DAF to read version from
set_cache_metadata <- function(cache_daf, base_daf) {
    # Set creation timestamp
    if (!dafr::has_scalar(cache_daf, "mcview_cache_created")) {
        tryCatch(
            dafr::set_scalar(cache_daf, "mcview_cache_created", as.numeric(Sys.time())),
            error = function(e) {
                cli::cli_alert_warning("Cache operation failed: {e$message}")
                NULL
            }
        )
    }

    # Copy base version if available
    base_version <- daf_scalar(base_daf, "mcview_data_version", default = NULL)
    if (!is.null(base_version) && !dafr::has_scalar(cache_daf, "mcview_cache_base_version")) {
        tryCatch(
            dafr::set_scalar(cache_daf, "mcview_cache_base_version", base_version),
            error = function(e) {
                cli::cli_alert_warning("Cache operation failed: {e$message}")
                NULL
            }
        )
    }

    # Store base hash for invalidation
    base_hash <- compute_base_hash(base_daf)
    if (!dafr::has_scalar(cache_daf, "mcview_cache_base_hash")) {
        tryCatch(
            dafr::set_scalar(cache_daf, "mcview_cache_base_hash", base_hash),
            error = function(e) {
                cli::cli_alert_warning("Cache operation failed: {e$message}")
                NULL
            }
        )
    }
}

#' Fingerprint a DAF storage directory
#'
#' Digests each file (excluding the derived/cache subtree) by
#' (relative path, size, mtime). mtime is included so a same-size content
#' rewrite is detected - e.g. re-annotating metacell types without changing
#' the metacell count rewrites a fixed-width vector `.data` at identical byte
#' size, which a size-only fingerprint would miss (serving a stale cache).
#' `rsync -a` and docker COPY preserve mtime, so ordinary redeploys do not
#' falsely invalidate; a non-preserving copy triggers one harmless rebuild.
#'
#' @param base_path DAF storage directory
#' @return A digest string ("" if the path is absent/empty)
#' @noRd
daf_dir_fingerprint <- function(base_path) {
    if (is.null(base_path) || !dir.exists(base_path)) {
        return("")
    }
    files <- list.files(base_path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
    # Exclude derived/cache subtrees if they live inside the DAF dir.
    # The current default is .mcview_derived; .mcview_cache is the legacy
    # name. Without this exclusion every precompute write flips the hash
    # and invalidates the cache that was just written.
    files <- files[!grepl("[/\\\\]\\.mcview_(derived|cache)(/|\\\\|$)", files)]
    if (length(files) == 0) {
        return("")
    }
    # Sort for deterministic order; pair each file with size + mtime, keyed by
    # its path relative to base_path (distinguishes same-basename files).
    files <- sort(files)
    rel <- fs::path_rel(files, base_path)
    sizes <- file.size(files)
    mtimes <- as.numeric(file.mtime(files))
    paste(rel, sizes, mtimes, sep = ":", collapse = "|")
}

#' Compute hash of base DAF for cache invalidation
#'
#' Combines axis-length and version fingerprints with a per-file
#' (relative path, size, mtime) digest of the DAF's storage directory.
#'
#' @param base_daf Base DAF object
#'
#' @return Hash string
compute_base_hash <- function(base_daf) {
    base_path <- tryCatch(daf_storage_path(base_daf), error = function(e) NULL)
    file_fingerprint <- daf_dir_fingerprint(base_path)

    components <- c(
        as.character(tryCatch(dafr::axis_length(base_daf, "metacell"), error = function(e) 0)),
        as.character(tryCatch(dafr::axis_length(base_daf, "gene"), error = function(e) 0)),
        as.character(tryCatch(dafr::axis_length(base_daf, "type"), error = function(e) 0)),
        daf_scalar(base_daf, "mcview_data_version", default = ""),
        file_fingerprint
    )
    digest::digest(paste(components, collapse = "|"), algo = "md5")
}

