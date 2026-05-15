# daf_core.R - Core DAF functionality for MCView
# Consolidates: daf_mode.R, daf_validation.R, daf_helpers.R

# ==============================================================================
# DAF Initialization
# ==============================================================================

#' Initialize MCView in single DAF mode
#'
#' Sets up environment variables to work with a single DAF object
#'
#' @param daf_obj DAF object
#' @param dataset_name Name for the dataset
#' @param config_file Optional YAML config file
#' @param profile Enable profiling
#' @param cache_in_daf Legacy cache flag (use cache_config instead)
#' @param cache_config Cache configuration from create_cache_config()
#'
#' @export
init_single_daf_mode <- function(daf_obj, dataset_name, config_file = NULL, profile = FALSE,
                                 cache_in_daf = NULL, cache_config = NULL) {
    # Validate DAF object
    validate_daf_for_mcview(daf_obj)

    # Load config from DAF scalars, YAML file, or auto-detect
    config <- extract_config_for_daf(daf_obj, dataset_name, config_file)
    config$profile <- profile

    # Handle legacy cache_in_daf parameter
    if (!is.null(cache_in_daf)) {
        config$cache_in_daf <- normalize_cache_flag(cache_in_daf)
    }

    # Extract cache configuration (new system)
    if (is.null(cache_config)) {
        cache_config <- extract_cache_config(config, daf_obj, dataset_name)
    }

    # Legacy support: if cache_in_daf was explicitly set, update cache_config
    if (!is.null(config$cache_in_daf)) {
        cache_config$enabled <- isTRUE(config$cache_in_daf)
    }

    # Get base path for derived DAF
    base_path <- daf_storage_path(daf_obj)

    # Initialize derived DAF
    cache_result <- init_derived_daf(
        base_daf = daf_obj,
        dataset_name = dataset_name,
        cache_config = cache_config,
        base_path = base_path
    )

    # Update config with derived/cache status
    config$cache_in_daf <- cache_config$enabled
    config$cache_config <- cache_config
    mcv_set("config", config)

    # Initialize mc_data with DAF dataset
    # Store base_daf, cache_daf, and complete_daf (as daf_obj for compatibility)
    mc_data <- list()
    mc_data[[dataset_name]] <- list(
        base_daf = daf_obj,
        cache_daf = cache_result$cache_daf,
        daf_obj = cache_result$complete_daf,
        top_cor_genes = list(),
        cache_path = cache_result$cache_path,
        needs_population = cache_result$needs_population
    )
    mcv_set("mc_data", mc_data)

    # Auto-detect or explicitly set cells DAF
    cells_daf_path <- config$cells_daf
    if (is.null(cells_daf_path)) {
        # Auto-detect: look for cells_clean sibling of metacells_clean
        metacells_path <- daf_storage_path(daf_obj)
        if (!is.null(metacells_path)) {
            candidate <- gsub("metacells", "cells", metacells_path)
            if (candidate != metacells_path && fs::dir_exists(candidate)) {
                cells_daf_path <- candidate
            }
        }
    }
    if (!is.null(cells_daf_path) && fs::dir_exists(cells_daf_path)) {
        tryCatch(
            {
                set_cells_daf(dataset_name, cells_daf_path)
                cli_alert_success("Cells DAF auto-detected for '{dataset_name}': {cells_daf_path}")
            },
            error = function(e) {
                cli_warn("Failed to set cells DAF for '{dataset_name}': {conditionMessage(e)}")
            }
        )
    }

    cli_alert_success("DAF dataset '{dataset_name}' loaded successfully")

    # Precompute top genes into cache DAF if vectors are missing. Reads go
    # through the complete chain (UMIs live in base); writes land in the
    # derived layer via chain_writer. Saves ~6 s per cold `metacell_types`
    # access on subsequent sessions.
    if (!is.null(cache_result$cache_daf) && !is.null(cache_result$complete_daf)) {
        precompute_daf_metacell_top_genes(cache_result$complete_daf)
    }

    # Return cache population status for caller to handle
    invisible(cache_result$needs_population)
}

#' Initialize MCView in multi-DAF mode
#'
#' Sets up environment variables to work with multiple DAF objects
#'
#' @param daf_list Named list of DAF objects
#' @param config_file Optional YAML config file
#' @param profile Enable profiling
#' @param cache_in_daf Legacy cache flag (use cache_config instead)
#' @param cache_config Cache configuration from create_cache_config()
#'
#' @export
init_multi_daf_mode <- function(daf_list, config_file = NULL, profile = FALSE,
                                cache_in_daf = NULL, cache_config = NULL) {
    # Validate all DAF objects
    purrr::walk(daf_list, validate_daf_for_mcview)

    # Use first DAF's config as base, allow YAML override
    config <- extract_config_for_daf(daf_list[[1]], names(daf_list)[1], config_file)
    config$profile <- profile

    # Handle legacy cache_in_daf parameter
    if (!is.null(cache_in_daf)) {
        config$cache_in_daf <- normalize_cache_flag(cache_in_daf)
    }

    # Initialize mc_data with cache for each dataset
    mc_data <- list()
    needs_population <- logical(length(daf_list))
    names(needs_population) <- names(daf_list)

    for (dataset_name in names(daf_list)) {
        daf_obj <- daf_list[[dataset_name]]

        # Extract per-dataset cache configuration
        ds_cache_config <- if (!is.null(cache_config)) {
            cache_config
        } else {
            extract_cache_config(config, daf_obj, dataset_name)
        }

        # Legacy support
        if (!is.null(config$cache_in_daf)) {
            ds_cache_config$enabled <- isTRUE(config$cache_in_daf)
        }

        # Get base path for derived DAF
        base_path <- daf_storage_path(daf_obj)

        # Initialize derived DAF
        cache_result <- init_derived_daf(
            base_daf = daf_obj,
            dataset_name = dataset_name,
            cache_config = ds_cache_config,
            base_path = base_path
        )

        mc_data[[dataset_name]] <- list(
            base_daf = daf_obj,
            cache_daf = cache_result$cache_daf,
            daf_obj = cache_result$complete_daf,
            top_cor_genes = list(),
            cache_path = cache_result$cache_path,
            needs_population = cache_result$needs_population
        )

        needs_population[[dataset_name]] <- cache_result$needs_population
    }

    mcv_set("mc_data", mc_data)

    # Update config
    config$cache_in_daf <- any(sapply(mc_data, function(x) !is.null(x$cache_daf)))
    mcv_set("config", config)

    # Auto-detect or explicitly set cells DAF for each dataset
    for (dataset_name in names(daf_list)) {
        daf_obj <- daf_list[[dataset_name]]
        cells_daf_path <- config$cells_daf
        if (is.null(cells_daf_path)) {
            # Auto-detect: look for cells_clean sibling of metacells_clean
            metacells_path <- daf_storage_path(daf_obj)
            if (!is.null(metacells_path)) {
                candidate <- gsub("metacells", "cells", metacells_path)
                if (candidate != metacells_path && fs::dir_exists(candidate)) {
                    cells_daf_path <- candidate
                }
            }
        }
        if (!is.null(cells_daf_path) && fs::dir_exists(cells_daf_path)) {
            tryCatch(
                {
                    set_cells_daf(dataset_name, cells_daf_path)
                    cli_alert_success("Cells DAF auto-detected for '{dataset_name}': {cells_daf_path}")
                },
                error = function(e) {
                    cli_warn("Failed to set cells DAF for '{dataset_name}': {conditionMessage(e)}")
                }
            )
        }
    }

    cli_alert_success("Loaded {length(daf_list)} DAF datasets: {paste(names(daf_list), collapse=', ')}")

    # Return which datasets need cache population
    invisible(needs_population)
}

#' Initialize atlas for projection comparison
#'
#' @param atlas_daf Atlas DAF object
#'
#' @export
init_atlas <- function(atlas_daf) {
    # Validate atlas DAF
    validate_daf_for_mcview(atlas_daf)

    # Store atlas
    mcv_set("atlas", atlas_daf)

    cli_alert_success("Atlas loaded successfully")
}

# ==============================================================================
# DAF Validation
# ==============================================================================

#' Validate DAF object for MCView compatibility
#'
#' Validates that a DAF object meets the MCView contract requirements.
#' Uses the formal contract system defined in daf_contracts.R.
#'
#' @param daf_obj DAF object to validate
#' @param strict If TRUE, performs additional content checks (negative values, etc.)
#' @param verbose If TRUE, prints detailed validation information
#'
#' @return TRUE if validation passes, throws error otherwise
#' @export
validate_daf_for_mcview <- function(daf_obj, strict = FALSE, verbose = FALSE) {
    if (verbose) cli_alert_info("Validating DAF object for MCView...")

    # Use contract-based validation
    result <- validate_mcview_contract(
        daf_obj,
        contract = mcview_core_contract(),
        strict = strict,
        verbose = verbose
    )

    if (!result$valid) {
        cli_abort(c("DAF validation failed:", result$errors))
    }

    # Print info if verbose
    if (verbose) {
        for (info in result$info) {
            cli_alert_success(info)
        }
    }

    # Content validation (data quality checks)
    if (strict) {
        validate_daf_content(daf_obj, verbose)
    }

    if (verbose) cli_alert_success("DAF object validation passed")
    invisible(TRUE)
}

validate_daf_content <- function(daf_obj, verbose = FALSE) {
    if (verbose) cli_alert_info("Validating data content...")

    # Validate UMI matrix content via two scalar reductions instead of
    # materialising the gene x metacell UMI matrix in R just to spot-check.
    umi_min <- tryCatch(daf_obj["@ metacell @ gene :: UMIs >> Min"],
                        error = function(e) NULL)
    umi_max <- tryCatch(daf_obj["@ metacell @ gene :: UMIs >> Max"],
                        error = function(e) NULL)

    if (is.null(umi_min) || is.null(umi_max)) {
        cli_abort("UMIs matrix is missing or unreadable")
    }
    if (!is.na(umi_min) && umi_min < 0) {
        cli_abort("UMIs matrix contains negative values")
    }
    if (!is.na(umi_max) && umi_max == 0) {
        cli_abort("UMIs matrix contains only zeros")
    }

    # Validate type consistency
    mc_types <- dafr::get_vector(daf_obj, "metacell", "type")
    type_axis <- dafr::axis_entries(daf_obj, "type")
    unknown_types <- setdiff(unique(mc_types), type_axis)

    if (length(unknown_types) > 0) {
        cli_abort("Unknown cell types: {paste(unknown_types, collapse=', ')}")
    }

    if (verbose) cli_alert_success("\u2713 Data content validation passed")
    invisible(TRUE)
}

# ==============================================================================
# Configuration Extraction
# ==============================================================================

#' Extract configuration for DAF mode
#'
#' @param daf_obj DAF object
#' @param dataset_name Default dataset name
#' @param config_file Optional YAML config file
#'
#' @export
extract_config_for_daf <- function(daf_obj, dataset_name, config_file = NULL) {
    config <- list()

    # Priority 1: YAML file (if provided)
    if (!is.null(config_file)) {
        if (!fs::file_exists(config_file)) {
            cli_abort("Config file {config_file} not found")
        }
        yaml_config <- yaml::read_yaml(config_file)
        config <- utils::modifyList(config, yaml_config)
        if (!is.null(config$cache_in_daf)) {
            config$cache_in_daf <- normalize_cache_flag(config$cache_in_daf)
        }
        cli_alert_info("Loaded configuration from {config_file}")
    }

    # Priority 2: DAF scalars (override auto-detection but not YAML)
    # Batch: get the set of all available scalars in one call, then only
    # fetch the ones that exist. This avoids ~16 individual has_scalar +
    # get_scalar roundtrips (each going through Julia).
    available_scalars <- dafr::scalars_set(daf_obj)

    # Helper: fetch a scalar only if it exists in the pre-queried set
    get_if_exists <- function(name) {
        if (name %in% available_scalars) {
            dafr::get_scalar(daf_obj, name)
        } else {
            NULL
        }
    }

    # Fetch all config scalars that exist (only those present, skip the rest)
    daf_title <- get_if_exists("mcview_title")
    daf_tabs <- get_if_exists("mcview_tabs")
    daf_light <- get_if_exists("mcview_light_version")
    daf_excluded <- get_if_exists("mcview_excluded_tabs")
    daf_cache <- get_if_exists("mcview_cache_in_daf")
    daf_cache_root <- get_if_exists("mcview_cache_daf_root")
    daf_sample_property <- get_if_exists("mcview_sample_property")

    # Apply DAF scalars if not overridden by YAML
    if (is.null(config$title) && !is.null(daf_title)) {
        config$title <- daf_title
    }

    if (is.null(config$tabs) && !is.null(daf_tabs)) {
        # Handle special "all" keyword - show all possible tabs
        if (trimws(daf_tabs) == "all") {
            config$tabs <- get_all_tab_names()
        } else {
            config$tabs <- trimws(strsplit(daf_tabs, ",")[[1]])
        }
    }

    if (is.null(config$light_version) && !is.null(daf_light)) {
        config$light_version <- daf_light
    }

    if (is.null(config$excluded_tabs) && !is.null(daf_excluded)) {
        config$excluded_tabs <- trimws(strsplit(daf_excluded, ",")[[1]])
    }
    if (is.null(config$cache_in_daf) && !is.null(daf_cache)) {
        config$cache_in_daf <- normalize_cache_flag(daf_cache)
    }
    if (is.null(config$cache_daf_root) && !is.null(daf_cache_root)) {
        config$cache_daf_root <- daf_cache_root
    }
    if (is.null(config$sample_property) && !is.null(daf_sample_property)) {
        config$sample_property <- daf_sample_property
    }

    # Priority 3: Auto-detection (lowest priority)
    if (is.null(config$title)) {
        config$title <- paste("MCView -", dataset_name)
    }

    if (is.null(config$tabs)) {
        config$tabs <- detect_available_tabs(daf_obj)
        cli_alert_info("Auto-detected tabs: {paste(config$tabs, collapse=', ')}")
    }

    if (is.null(config$light_version)) {
        config$light_version <- FALSE
    }
    if (is.null(config$cache_in_daf)) {
        # Default to enabled: see extract_cache_config() for rationale. Users
        # can opt out via mcview_cache_in_daf DAF scalar, YAML config, or the
        # cache_in_daf= argument of init_*_daf_mode().
        config$cache_in_daf <- TRUE
    }

    # Add metacells version if available (uses the already-fetched scalar set)
    metacells_version <- get_if_exists("metacells_algorithm")
    if (!is.null(metacells_version)) {
        config$metacells_version <- metacells_version
    }

    # Store about markdown if provided in DAF
    about_markdown <- get_if_exists("mcview_about_markdown")
    if (!is.null(about_markdown)) {
        mcv_set("about_markdown", about_markdown)
    }

    return(config)
}

#' Auto-detect available tabs based on DAF content
#'
#' Uses a fast batch-query approach: queries all available axes, vectors,
#' matrices, and scalars from the DAF once, then checks tab requirements
#' in pure R without per-tab Julia roundtrips.
#'
#' If the DAF contains a pre-stored `mcview_available_tabs` scalar (a
#' comma-separated list of tab names), that value is used directly and
#' all validation is skipped.
#'
#' @param daf_obj DAF object
#'
#' @return Character vector of available tab names
#' @export
detect_available_tabs <- function(daf_obj) {
    # Fast path: check for pre-stored available tabs scalar
    if (dafr::has_scalar(daf_obj, "mcview_available_tabs")) {
        stored_tabs <- dafr::get_scalar(daf_obj, "mcview_available_tabs")
        tabs <- trimws(strsplit(stored_tabs, ",")[[1]])
        # Filter to canonical names only
        tabs <- tabs[tabs %in% MCVIEW_TAB_NAMES]
        if (length(tabs) > 0) {
            return(tabs)
        }
        # If stored value was empty/invalid, fall through to detection
    }

    # Batch-query DAF metadata with targeted checks (not N² axis-pair queries)
    available_axes <- dafr::axes_set(daf_obj)

    # Query vectors only for axes relevant to tab detection
    tab_relevant_axes <- intersect(available_axes, c("metacell", "gene", "type", "cell", "gg_mc_top_cor"))
    available_vectors <- list()
    for (axis in tab_relevant_axes) {
        vecs <- tryCatch(dafr::vectors_set(daf_obj, axis), error = function(e) character(0))
        available_vectors[[axis]] <- vecs
    }

    # Only query the specific matrix pairs needed for tab detection
    # (3 targeted queries instead of N² queries for all axis pairs)
    available_matrices <- list()
    matrix_checks <- list(
        c("metacell", "gene"),
        c("gene", "metacell")
    )
    for (pair in matrix_checks) {
        if (pair[1] %in% available_axes && pair[2] %in% available_axes) {
            mats <- tryCatch(
                dafr::matrices_set(daf_obj, pair[1], pair[2]),
                error = function(e) character(0)
            )
            if (length(mats) > 0) {
                key <- paste0(pair[1], ",", pair[2])
                available_matrices[[key]] <- mats
            }
        }
    }

    # Helper: check if a vector exists in the cached metadata
    has_vec <- function(axis, name) {
        axis %in% available_axes && name %in% available_vectors[[axis]]
    }

    # Helper: check if a matrix exists in the cached metadata
    has_mat <- function(rows_axis, cols_axis, name) {
        key <- paste0(rows_axis, ",", cols_axis)
        key %in% names(available_matrices) && name %in% available_matrices[[key]]
    }

    # Helper: check if an axis exists
    has_ax <- function(axis) {
        axis %in% available_axes
    }

    # Validate core contract first (required for all tabs)
    core_ok <- has_ax("metacell") && has_ax("gene") && has_ax("type") &&
        has_vec("metacell", "total_UMIs") &&
        has_vec("metacell", "type") &&
        has_vec("type", "color") &&
        (has_mat("metacell", "gene", "UMIs") || has_mat("gene", "metacell", "UMIs")) &&
        ((has_vec("metacell", "x") && has_vec("metacell", "y")) ||
            (has_vec("metacell", "u") && has_vec("metacell", "v")))

    if (!core_ok) {
        return(character(0))
    }

    # Check each tab's requirements against cached metadata.
    # A tab is available if all its REQUIRED fields exist.
    # Optional fields do not affect availability.
    tabs <- character(0)

    # About - always available if core passes
    tabs <- c(tabs, "About")

    # Manifold - core is sufficient (graph is optional)
    tabs <- c(tabs, "Manifold")

    # Genes - core is sufficient (top genes, lateral, noisy are optional)
    tabs <- c(tabs, "Genes")

    # Diff. Expression - core is sufficient (lateral, noisy are optional)
    tabs <- c(tabs, "Diff. Expression")

    # Cell types - core is sufficient (n_cell is optional)
    tabs <- c(tabs, "Cell types")

    # QC - core is sufficient (all QC-specific fields are optional)
    tabs <- c(tabs, "QC")

    # Markers - requires gene.is_marker
    if (has_vec("gene", "is_marker")) {
        tabs <- c(tabs, "Markers")
    }

    # Gene modules - requires gene.module
    if (has_vec("gene", "module")) {
        tabs <- c(tabs, "Gene modules")
    }

    # Projection QC / Atlas - all individual fields are marked optional in
    # their contracts, but the tabs are meaningful only when the DAF actually
    # carries projection metadata. Without any of these markers both tabs
    # render empty panels (no fitted-genes table, no projected scatter, etc.).
    has_projection <- has_vec("metacell", "projected_correlation") ||
        has_vec("metacell", "projected_type") ||
        has_vec("metacell", "similar") ||
        has_vec("metacell", "is_similar") ||
        has_vec("metacell", "atlas_metacell") ||
        has_mat("gene", "metacell", "projected_fold") ||
        has_mat("metacell", "gene", "projected_fraction") ||
        has_mat("metacell", "gene", "projected_geomean_fraction") ||
        has_mat("metacell", "gene", "corrected_fraction") ||
        has_mat("metacell", "gene", "corrected_geomean_fraction")

    if (has_projection) {
        tabs <- c(tabs, "Projection QC", "Atlas")
    }

    # Annotate - always available if core passes
    tabs <- c(tabs, "Annotate")

    # Samples - always include optimistically. Cell data may live in the main
    # DAF or in a separate cells DAF loaded after tab detection. The runtime
    # has_samples() check in app_server.R removes the tab if no grouping
    # fields are actually available.
    tabs <- c(tabs, "Samples")

    # Gene correlation - always available with core data (UMIs + genes).
    # Pre-computed correlations (gg_mc_top_cor or mcview_marker_correlations axes)
    # are optional; correlations can be computed on-the-fly.
    tabs <- c(tabs, "Gene correlation")

    # Flow - requires metacell.time (custom rule)
    if (has_vec("metacell", "time")) {
        tabs <- c(tabs, "Flow")
    }

    # Return in canonical order
    tabs <- MCVIEW_TAB_NAMES[MCVIEW_TAB_NAMES %in% tabs]

    return(tabs)
}

#' Pre-store available tabs in a DAF
#'
#' Runs detect_available_tabs() and stores the result as a scalar
#' in the DAF. On subsequent startups, detect_available_tabs() will
#' find this scalar and skip all tab validation entirely.
#'
#' @param daf_obj DAF object (must be writable)
#' @return Character vector of available tabs (invisibly)
#' @export
store_available_tabs <- function(daf_obj) {
    tabs <- detect_available_tabs(daf_obj)
    tabs_str <- paste(tabs, collapse = ",")
    dafr::set_scalar(daf_obj, "mcview_available_tabs", tabs_str)
    cli_alert_success("Stored {length(tabs)} available tabs in DAF: {tabs_str}")
    invisible(tabs)
}

# ==============================================================================
# Configuration Helpers
# ==============================================================================

#' Add MCView configuration to DAF object
#'
#' @param daf_obj DAF object
#' @param title Application title
#' @param tabs Character vector of tabs or comma-separated string
#' @param light_version Whether to use light version
#' @param excluded_tabs Character vector of excluded tabs
#' @param cache_in_daf Whether to store cache in the DAF
#' @param cache_daf_root Root path for cache DAF storage
#' @export
add_mcview_config <- function(daf_obj,
                              title = NULL,
                              tabs = NULL,
                              light_version = FALSE,
                              excluded_tabs = NULL,
                              cache_in_daf = NULL,
                              cache_daf_root = NULL) {
    if (!is.null(title)) {
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_title", title)
    }

    if (!is.null(tabs)) {
        if (is.character(tabs) && length(tabs) > 1) {
            tabs <- paste(tabs, collapse = ",")
        }
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_tabs", tabs)
    }

    daf_obj <- dafr::set_scalar(daf_obj, "mcview_light_version", light_version)

    if (!is.null(excluded_tabs)) {
        if (is.character(excluded_tabs) && length(excluded_tabs) > 1) {
            excluded_tabs <- paste(excluded_tabs, collapse = ",")
        }
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_excluded_tabs", excluded_tabs)
    }

    if (!is.null(cache_in_daf)) {
        cache_in_daf <- normalize_cache_flag(cache_in_daf)
        if (!is.null(cache_in_daf)) {
            daf_obj <- dafr::set_scalar(daf_obj, "mcview_cache_in_daf", cache_in_daf)
        }
    }
    if (!is.null(cache_daf_root)) {
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_cache_daf_root", cache_daf_root)
    }

    daf_obj
}

#' Create template YAML config file
#'
#' @param output_file Path to output YAML file
#' @export
create_mcview_config_template <- function(output_file) {
    template_config <- list(
        title = "MCView Analysis",
        tabs = c("About", "Manifold", "Genes", "Diff. Expression", "Markers"),
        light_version = FALSE,
        excluded_tabs = c("QC", "Samples"),
        cache_in_daf = TRUE,
        cache_daf_root = "cache_daf"
    )

    yaml::write_yaml(template_config, output_file)
    cli_alert_success("Created config template at {output_file}")
    cli_alert_info("Use with: run_app(daf_obj, config_file = '{output_file}')")
}

#' Load YAML configuration into DAF object
#'
#' @param daf_obj DAF object
#' @param config_file Path to YAML config file
#' @export
load_yaml_config_to_daf <- function(daf_obj, config_file) {
    if (!fs::file_exists(config_file)) {
        cli_abort("Config file {config_file} not found")
    }

    yaml_config <- yaml::read_yaml(config_file)

    if (!is.null(yaml_config$title)) {
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_title", yaml_config$title)
    }

    if (!is.null(yaml_config$tabs)) {
        tabs_str <- if (length(yaml_config$tabs) > 1) {
            paste(yaml_config$tabs, collapse = ",")
        } else {
            yaml_config$tabs
        }
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_tabs", tabs_str)
    }

    if (!is.null(yaml_config$light_version)) {
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_light_version", yaml_config$light_version)
    }

    if (!is.null(yaml_config$excluded_tabs)) {
        excluded_str <- if (length(yaml_config$excluded_tabs) > 1) {
            paste(yaml_config$excluded_tabs, collapse = ",")
        } else {
            yaml_config$excluded_tabs
        }
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_excluded_tabs", excluded_str)
    }

    if (!is.null(yaml_config$cache_in_daf)) {
        cache_in_daf <- normalize_cache_flag(yaml_config$cache_in_daf)
        if (!is.null(cache_in_daf)) {
            daf_obj <- dafr::set_scalar(daf_obj, "mcview_cache_in_daf", cache_in_daf)
        }
    }
    if (!is.null(yaml_config$cache_daf_root)) {
        daf_obj <- dafr::set_scalar(daf_obj, "mcview_cache_daf_root", yaml_config$cache_daf_root)
    }

    cli_alert_success("Loaded configuration from {config_file} into DAF object")
    daf_obj
}
