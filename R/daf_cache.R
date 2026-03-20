# Helpers for DAF-backed cache settings and write access.
# This module implements a cache DAF layer that chains with the base DAF
# using dafr::complete_daf() for unified read access.

# ==============================================================================
# Cache DAF Initialization (complete_daf based)
# ==============================================================================

#' Initialize cache DAF for a dataset
#'
#' Creates or opens a cache DAF and chains it with the base DAF using
#' complete_daf. The cache DAF stores computed data that can be persisted
#' across sessions.
#'
#' @param base_daf Base DAF object (read-only or writable)
#' @param dataset_name Dataset name (used for cache directory)
#' @param cache_config Cache configuration from extract_cache_config()
#' @param base_path Path to base DAF (for relative path calculation)
#'
#' @return List with:
#'   - cache_daf: The writable cache DAF (or NULL if caching disabled)
#'   - complete_daf: The chained DAF for reads (cache has priority)
#'   - needs_population: Whether cache needs to be populated
#'   - cache_path: Path to cache directory (for files cache)
#'
#' @export
init_cache_daf <- function(base_daf, dataset_name, cache_config, base_path = NULL) {
    # If caching disabled, return base DAF as-is
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

    # Determine cache path
    cache_dir <- cache_config$cache_dir %||% ".mcview_cache"
    cache_path <- resolve_cache_path(cache_dir, base_path, dataset_name)

    # Create cache DAF based on type
    cache_daf <- NULL
    if (cache_config$type == "memory") {
        cache_daf <- create_memory_cache_daf(base_daf, dataset_name)
    } else {
        cache_daf <- create_files_cache_daf(base_daf, cache_path, base_path)
    }

    if (is.null(cache_daf)) {
        cli::cli_alert_warning("Failed to create cache DAF for {dataset_name}, using base DAF only")
        return(list(
            cache_daf = NULL,
            complete_daf = base_daf,
            needs_population = FALSE,
            cache_path = NULL
        ))
    }

    # Check if cache is valid
    needs_population <- !is_cache_valid(cache_daf, base_daf, cache_config$invalidation)

    # Create complete DAF chain (cache as leaf, base as parent)
    complete_daf <- NULL
    if (cache_config$type == "files" && !is.null(cache_path)) {
        # For files cache, use complete_daf which follows base_daf_repository
        complete_daf <- tryCatch(
            dafr::complete_daf(cache_path, mode = "r+"),
            error = function(e) {
                cli::cli_alert_warning("Failed to create complete DAF: {e$message}")
                NULL
            }
        )
    }

    # Fallback to chain_writer if complete_daf failed or memory cache
    if (is.null(complete_daf)) {
        complete_daf <- tryCatch(
            dafr::chain_writer(list(dafr::read_only(base_daf), cache_daf)),
            error = function(e) {
                cli::cli_alert_warning("Failed to chain DAFs: {e$message}")
                base_daf
            }
        )
    }

    label <- if (cache_config$type == "memory") "in-memory" else cache_path
    cli::cli_alert_info("Cache DAF initialized for {dataset_name}: {label}")

    list(
        cache_daf = cache_daf,
        complete_daf = complete_daf,
        needs_population = needs_population,
        cache_path = if (cache_config$type == "files") cache_path else NULL
    )
}

#' Resolve cache path from configuration
#'
#' @param cache_dir Cache directory from config (relative or absolute)
#' @param base_path Path to base DAF
#' @param dataset_name Dataset name
#'
#' @return Absolute path to cache directory
resolve_cache_path <- function(cache_dir, base_path, dataset_name) {
    if (is.null(cache_dir) || !nzchar(cache_dir)) {
        cache_dir <- ".mcview_cache"
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
        error = function(e) NULL
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
            error = function(e) NULL
        )
    } else {
        # Create new cache directory
        fs::dir_create(cache_path, recurse = TRUE)
        cache_daf <- tryCatch(
            dafr::files_daf(cache_path, mode = "w+"),
            error = function(e) NULL
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
                error = function(e) NULL
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
                error = function(e) NULL
            )
            if (!is.null(entries)) {
                tryCatch(
                    dafr::add_axis(cache_daf, axis_name, entries),
                    error = function(e) NULL
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
            error = function(e) NULL
        )
    }

    # Copy base version if available
    base_version <- daf_scalar(base_daf, "mcview_data_version", default = NULL)
    if (!is.null(base_version) && !dafr::has_scalar(cache_daf, "mcview_cache_base_version")) {
        tryCatch(
            dafr::set_scalar(cache_daf, "mcview_cache_base_version", base_version),
            error = function(e) NULL
        )
    }

    # Store base hash for invalidation
    base_hash <- compute_base_hash(base_daf)
    if (!dafr::has_scalar(cache_daf, "mcview_cache_base_hash")) {
        tryCatch(
            dafr::set_scalar(cache_daf, "mcview_cache_base_hash", base_hash),
            error = function(e) NULL
        )
    }
}

#' Compute hash of base DAF for cache invalidation
#'
#' @param base_daf Base DAF object
#'
#' @return Hash string based on axis lengths and version
compute_base_hash <- function(base_daf) {
    components <- c(
        as.character(tryCatch(dafr::axis_length(base_daf, "metacell"), error = function(e) 0)),
        as.character(tryCatch(dafr::axis_length(base_daf, "gene"), error = function(e) 0)),
        as.character(tryCatch(dafr::axis_length(base_daf, "type"), error = function(e) 0)),
        daf_scalar(base_daf, "mcview_data_version", default = "")
    )
    digest::digest(paste(components, collapse = "|"), algo = "md5")
}

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

    cli::cli_alert_info("Cache invalidated for dataset: {dataset}")
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
        error = function(e) NULL
    )

    # Update base version
    base_version <- daf_scalar(base_daf, "mcview_data_version", default = NULL)
    if (!is.null(base_version)) {
        tryCatch(
            dafr::set_scalar(cache_daf, "mcview_cache_base_version", base_version, overwrite = TRUE),
            error = function(e) NULL
        )
    }

    # Update base hash
    base_hash <- compute_base_hash(base_daf)
    tryCatch(
        dafr::set_scalar(cache_daf, "mcview_cache_base_hash", base_hash, overwrite = TRUE),
        error = function(e) NULL
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
    fields <- daf_description_fields(daf_obj, deep = TRUE)
    if (!is.null(fields$path) && nzchar(fields$path)) {
        return(fields$path)
    }
    if (!is.null(fields$base) && nzchar(fields$base)) {
        base_path <- sub("^\\S+\\s+", "", fields$base)
        if (nzchar(base_path)) {
            return(base_path)
        }
    }
    if (!is.null(fields$root) && nzchar(fields$root)) {
        match <- regexpr("/[^ ]+\\.h5[^ ]*", fields$root)
        if (match > 0) {
            return(regmatches(fields$root, match))
        }
    }
    NULL
}

cache_in_daf_enabled <- function(config = mcv_get("config")) {
    isTRUE(config$cache_in_daf)
}

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
        error = function(e) FALSE
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

    # Vectorized top gene computation using max.col on transposed matrix
    # mc_egc is gene x metacell, transpose to metacell x gene for max.col
    # jlview_t returns a lazy ALTREP view; as.matrix() materializes for max.col
    egc_t <- as.matrix(jlview::jlview_t(mc_egc))

    # Top 1: find column index of max value per metacell
    top1_idx <- max.col(egc_t, ties.method = "first")
    top1_gene <- gene_names[top1_idx]
    top1_lfp <- log2(egc_t[cbind(seq_len(nrow(egc_t)), top1_idx)] + egc_epsilon)

    # Top 2: mask top1 values, find next max
    egc_t[cbind(seq_len(nrow(egc_t)), top1_idx)] <- -Inf
    top2_idx <- max.col(egc_t, ties.method = "first")
    top2_gene <- gene_names[top2_idx]
    top2_lfp <- log2(egc_t[cbind(seq_len(nrow(egc_t)), top2_idx)] + egc_epsilon)

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
        error = function(e) FALSE
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
        error = function(e) FALSE
    )

    invisible(success)
}

precompute_daf_cache <- function(daf_obj,
                                 correlations = TRUE,
                                 metacell_top_genes = TRUE,
                                 gene_stats = TRUE,
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

    invisible(any_done)
}

# ==============================================================================
# External Populate API
# ==============================================================================

#' Populate MCView cache for a DAF dataset
#'
#' This function can be called standalone (outside MCView) to precompute
#' cached data. Useful for CI/CD pipelines, batch processing, or scheduled jobs.
#'
#' @param daf_path Path to base DAF (directory for files_daf, or H5 file)
#' @param cache_dir Cache directory path (absolute or relative to daf_path)
#' @param cache_type Cache type: "files" (persistent) or "memory" (session only)
#' @param what What to cache: NULL (all), or character vector subset of:
#'   "correlations", "metacell_top_genes", "gene_stats"
#' @param force Recompute even if cache exists and is valid
#' @param verbose Print progress messages
#'
#' @return Invisibly returns the cache path (for files) or TRUE (for memory)
#'
#' @examples
#' \dontrun{
#' # Populate all cache data
#' populate_mcview_cache("/path/to/daf")
#'
#' # Populate specific cache items
#' populate_mcview_cache(
#'     "/path/to/daf",
#'     what = c("correlations", "gene_stats"),
#'     verbose = TRUE
#' )
#'
#' # Use custom cache directory
#' populate_mcview_cache(
#'     "/path/to/daf",
#'     cache_dir = "/shared/cache/my_dataset"
#' )
#' }
#' @export
populate_mcview_cache <- function(
    daf_path,
    cache_dir = ".mcview_cache",
    cache_type = "files",
    what = NULL,
    force = FALSE,
    verbose = TRUE) {
    # Initialize dafr if needed
    if (!requireNamespace("dafr", quietly = TRUE)) {
        cli::cli_abort("dafr package is required for cache population")
    }

    # Ensure Julia is set up (use sysimage if available for faster startup)
    tryCatch(
        {
            setup_args <- list()
            sysimage <- Sys.getenv("JULIA_SYSIMAGE", "")
            if (!nzchar(sysimage)) {
                conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
                if (nzchar(conda_prefix)) {
                    candidate <- file.path(conda_prefix, "share", "julia", "sysimage_daf.so")
                    if (file.exists(candidate)) sysimage <- candidate
                }
            }
            if (nzchar(sysimage) && file.exists(sysimage)) {
                setup_args$sysimage_path <- sysimage
            }
            do.call(dafr::setup_daf, setup_args)
        },
        error = function(e) {
            cli::cli_abort("Failed to initialize dafr: {e$message}")
        }
    )

    if (verbose) cli::cli_alert_info("Opening DAF: {daf_path}")

    # Open base DAF
    base_daf <- tryCatch(
        {
            if (dir.exists(daf_path)) {
                dafr::files_daf(daf_path, mode = "r")
            } else if (file.exists(daf_path) && grepl("\\.h5$", daf_path)) {
                dafr::h5df(daf_path, mode = "r")
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

    if (verbose) cli::cli_alert_info("Initializing cache DAF (type: {cache_type})")

    # Initialize cache DAF
    cache_result <- init_cache_daf(
        base_daf = base_daf,
        dataset_name = dataset_name,
        cache_config = cache_config,
        base_path = daf_path
    )

    if (is.null(cache_result$cache_daf)) {
        cli::cli_abort("Failed to create cache DAF")
    }

    cache_daf <- cache_result$cache_daf

    # Determine what to cache
    all_items <- c("correlations", "metacell_top_genes", "gene_stats")
    if (is.null(what)) {
        what <- all_items
    } else {
        what <- intersect(what, all_items)
        if (length(what) == 0) {
            cli::cli_abort("No valid cache items specified. Choose from: {paste(all_items, collapse=', ')}")
        }
    }

    if (verbose) cli::cli_alert_info("Populating cache: {paste(what, collapse=', ')}")

    # Populate cache
    results <- list()

    if ("correlations" %in% what) {
        if (verbose) cli::cli_alert("Computing gene correlations...")
        results$correlations <- precompute_daf_correlations(cache_daf, force = force)
        if (verbose && results$correlations) cli::cli_alert_success("Gene correlations cached")
    }

    if ("metacell_top_genes" %in% what) {
        if (verbose) cli::cli_alert("Computing metacell top genes...")
        results$metacell_top_genes <- precompute_daf_metacell_top_genes(cache_daf, force = force)
        if (verbose && results$metacell_top_genes) cli::cli_alert_success("Metacell top genes cached")
    }

    if ("gene_stats" %in% what) {
        if (verbose) cli::cli_alert("Computing gene statistics...")
        results$gene_stats <- precompute_daf_gene_stats(cache_daf, force = force)
        if (verbose && results$gene_stats) cli::cli_alert_success("Gene statistics cached")
    }

    # Update cache metadata
    update_cache_metadata(cache_daf, base_daf)

    if (verbose) {
        n_done <- sum(unlist(results))
        if (n_done > 0) {
            cli::cli_alert_success("Cache population complete: {n_done} items computed")
        } else {
            cli::cli_alert_info("Cache is already up to date (use force=TRUE to recompute)")
        }
    }

    invisible(cache_result$cache_path %||% TRUE)
}

#' Populate cache for a dataset in a running MCView session
#'
#' @param dataset Dataset name
#' @param what What to cache (NULL = all)
#' @param force Recompute even if present
#' @param verbose Print progress
#'
#' @export
populate_dataset_cache <- function(dataset, what = NULL, force = FALSE, verbose = TRUE) {
    cache_daf <- get_cache_daf(dataset)
    base_daf <- get_base_daf(dataset)

    if (is.null(cache_daf)) {
        if (verbose) cli::cli_alert_warning("No cache DAF available for dataset: {dataset}")
        return(invisible(FALSE))
    }

    # Determine what to cache
    all_items <- c("correlations", "metacell_top_genes", "gene_stats")
    if (is.null(what)) {
        what <- all_items
    }

    results <- list()

    if ("correlations" %in% what) {
        if (verbose) cli::cli_alert("Computing gene correlations for {dataset}...")
        results$correlations <- precompute_daf_correlations(cache_daf, force = force)
    }

    if ("metacell_top_genes" %in% what) {
        if (verbose) cli::cli_alert("Computing metacell top genes for {dataset}...")
        results$metacell_top_genes <- precompute_daf_metacell_top_genes(cache_daf, force = force)
    }

    if ("gene_stats" %in% what) {
        if (verbose) cli::cli_alert("Computing gene statistics for {dataset}...")
        results$gene_stats <- precompute_daf_gene_stats(cache_daf, force = force)
    }

    # Update cache metadata
    if (!is.null(base_daf)) {
        update_cache_metadata(cache_daf, base_daf)
    }

    if (verbose) {
        n_done <- sum(unlist(results))
        cli::cli_alert_success("Cache population for {dataset}: {n_done} items computed")
    }

    invisible(any(unlist(results)))
}

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
    # Resolve cache path
    if (is.null(cache_dir)) {
        cache_dir <- ".mcview_cache"
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
        gene_stats = dafr::has_vector(cache_daf, "gene", "mcview_cache_gene_max_umis")
    )

    result$missing <- names(expected_items)[!unlist(expected_items)]
    result$complete <- length(result$missing) == 0

    result
}
