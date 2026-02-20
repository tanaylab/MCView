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

    # Get base path for cache
    base_path <- daf_storage_path(daf_obj)

    # Initialize cache DAF using new system
    cache_result <- init_cache_daf(
        base_daf = daf_obj,
        dataset_name = dataset_name,
        cache_config = cache_config,
        base_path = base_path
    )

    # Update config with cache status
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

    cli_alert_success("DAF dataset '{dataset_name}' loaded successfully")

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

        # Get base path for cache
        base_path <- daf_storage_path(daf_obj)

        # Initialize cache DAF
        cache_result <- init_cache_daf(
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
        strict = FALSE,
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

#' Legacy validation function (deprecated, use validate_mcview_contract)
#' @noRd
validate_daf_structure <- function(daf_obj, verbose = FALSE) {
    # Delegate to contract-based validation
    result <- validate_mcview_contract(
        daf_obj,
        contract = mcview_core_contract(),
        strict = FALSE,
        verbose = verbose
    )

    if (!result$valid) {
        cli_abort(c("DAF validation failed:", result$errors))
    }

    invisible(TRUE)
}

validate_daf_content <- function(daf_obj, verbose = FALSE) {
    if (verbose) cli_alert_info("Validating data content...")

    # Validate UMI matrix content
    umat <- dafr::get_matrix(daf_obj, "metacell", "gene", "UMIs")

    if (any(umat < 0, na.rm = TRUE)) {
        cli_abort("UMIs matrix contains negative values")
    }

    if (all(umat == 0, na.rm = TRUE)) {
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
        config <- modifyList(config, yaml_config)
        if (!is.null(config$cache_in_daf)) {
            config$cache_in_daf <- normalize_cache_flag(config$cache_in_daf)
        }
        cli_alert_info("Loaded configuration from {config_file}")
    }

    # Priority 2: DAF scalars (override auto-detection but not YAML)
    daf_title <- daf_scalar(daf_obj, "mcview_title", default = NULL)
    daf_tabs <- daf_scalar(daf_obj, "mcview_tabs", default = NULL)
    daf_light <- daf_scalar(daf_obj, "mcview_light_version", default = NULL)
    daf_excluded <- daf_scalar(daf_obj, "mcview_excluded_tabs", default = NULL)
    daf_cache <- daf_scalar(daf_obj, "mcview_cache_in_daf", default = NULL)
    daf_cache_root <- daf_scalar(daf_obj, "mcview_cache_daf_root", default = NULL)

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
        config$cache_in_daf <- FALSE
    }

    # Add metacells version if available
    metacells_version <- daf_scalar(daf_obj, "metacells_algorithm", default = NULL)
    if (!is.null(metacells_version)) {
        config$metacells_version <- metacells_version
    }

    # Store about markdown if provided in DAF
    about_markdown <- daf_scalar(daf_obj, "mcview_about_markdown", default = NULL)
    if (!is.null(about_markdown)) {
        mcv_set("about_markdown", about_markdown)
    }

    return(config)
}

#' Auto-detect available tabs based on DAF content
#'
#' Uses the contract-based validation system to determine which tabs
#' have sufficient data to be displayed.
#'
#' @param daf_obj DAF object
#'
#' @return Character vector of available tab names
#' @export
detect_available_tabs <- function(daf_obj) {
    # Use contract-based tab validation
    tab_results <- validate_mcview_tabs(daf_obj, verbose = FALSE)

    # Get tabs that are available
    available_tabs <- names(tab_results)[sapply(tab_results, function(x) x$available)]

    # Ensure proper ordering
    all_tabs <- c(
        "About", "Manifold", "Genes", "Diff. Expression", "Cell types",
        "QC", "Markers", "Gene modules", "Inner-fold", "Stdev-fold",
        "Projection QC", "Atlas", "Annotate"
    )

    # Return in order, only those that are available
    tabs <- all_tabs[all_tabs %in% available_tabs]

    return(tabs)
}

#' Get list of all possible tab names
#'
#' @return Character vector of all tab names
#' @export
get_all_tab_names <- function() {
    c(
        "About", "Manifold", "Genes", "Diff. Expression", "Cell types",
        "QC", "Markers", "Gene modules", "Inner-fold", "Stdev-fold",
        "Projection QC", "Atlas", "Annotate"
    )
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
