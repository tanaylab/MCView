# ==============================================================================
# Lazy future::plan initialization
# ==============================================================================

# Track whether future::plan(multisession) has been activated
.future_plan_state <- new.env(parent = emptyenv())
.future_plan_state$initialized <- FALSE

#' Ensure ggplotly's measurement device can open
#'
#' `plotly::ggplotly()` opens a temporary PNG device to measure panel sizes.
#' On R builds without X11 capability, the default `bitmapType = "Xlib"`
#' makes that device fail. Switch to cairo when it's available.
#' @return TRUE invisibly
#' @noRd
.ensure_safe_bitmap_type <- function() {
    if (identical(getOption("bitmapType"), "Xlib") &&
        isFALSE(capabilities("X11")) && isTRUE(capabilities("cairo"))) {
        options(bitmapType = "cairo")
    }
    invisible(TRUE)
}

#' Assert that dafr memory-mapping is enabled
#'
#' MCView requires `dafr.mmap = TRUE` to avoid materializing full matrices
#' through R copies. Aborts with a helpful message if the option is not set.
#'
#' @return TRUE invisibly
#' @noRd
.assert_dafr_mmap <- function() {
    # dafr's .onLoad sets getOption("dafr.mmap") = TRUE (see dafr-native
    # R/zzz.R + R/options.R). We read via getOption rather than
    # dafr::dafr_opt — the latter is internal (not exported) so calling it
    # from a subprocess that hasn't pulled MCView's full import surface
    # raises a namespace error. requireNamespace forces .onLoad to fire even
    # if dafr hasn't been attached yet.
    requireNamespace("dafr", quietly = TRUE)
    if (!isTRUE(getOption("dafr.mmap"))) {
        cli::cli_abort(c(
            "MCView requires {.code dafr.mmap = TRUE}.",
            "i" = "Set via {.code options(dafr.mmap = TRUE)} before {.code library(MCView)}."
        ))
    }
    invisible(TRUE)
}

#' Ensure future plan is set to multisession (lazy, one-time)
#'
#' Called on first use of future_promise (e.g., Flow tab's plot_vein).
#' Avoids the ~5.7s cost of spawning worker processes at app startup
#' when the Flow tab may never be visited.
#'
#' @return TRUE invisibly
#' @noRd
ensure_future_plan <- function() {
    if (!.future_plan_state$initialized) {
        .future_plan_state$initialized <- TRUE
        # Use multisession instead of multicore to avoid forking a process
        # with embedded Julia (fork is unsupported and causes segfaults)
        future::plan(future::multisession)
    }
    invisible(TRUE)
}

#' Run the MCView Application
#'
#' Run the MCView application with DAF data.
#'
#' @param daf A DAF object, a named list of DAF objects, or a path to a DAF directory/H5 file
#' @param name Name for the dataset when using a single DAF object (default: "data").
#'   Ignored when `daf` is a named list.
#' @param cells_daf Optional path to a cells-level DAF directory. When provided,
#'   enables cell-level features (grouping, pseudobulk DE, QC stats). If NULL,
#'   auto-detection looks for a sibling directory with "metacells" replaced by "cells".
#' @param atlas Optional atlas DAF object for projection comparison
#' @param config_file Optional path to YAML configuration file
#' @param port App port
#' @param host App host
#' @param launch.browser Launch web browser after app start
#' @param profile Enable profiling for debugging
#' @param cache_in_daf Whether to store runtime cache in the DAF object
#' @param ... Additional options passed to shinyApp
#'
#' @examples
#' \dontrun{
#' # Single DAF object
#' library(dafr)
#' daf_obj <- dafr::files_daf("path/to/daf")
#' run_app(daf_obj)
#' run_app(daf_obj, name = "my_experiment")
#'
#' # Multiple DAF objects
#' embryo_daf <- dafr::files_daf("path/to/embryo")
#' adult_daf <- dafr::files_daf("path/to/adult")
#' run_app(list(
#'     embryo = embryo_daf,
#'     adult = adult_daf
#' ))
#'
#' # With atlas for projection
#' query_daf <- dafr::files_daf("path/to/query")
#' atlas_daf <- dafr::files_daf("path/to/atlas")
#' run_app(query_daf, atlas = atlas_daf)
#'
#' # With cell-level data
#' run_app(daf_obj, cells_daf = "path/to/cells_clean")
#'
#' # From path (auto-loads DAF)
#' run_app("path/to/daf")
#' run_app("path/to/data.h5")
#'
#' # With external YAML configuration
#' run_app(daf_obj, config_file = "config.yaml")
#' }
#'
#' @inheritDotParams shiny::shinyApp
#'
#' @export
run_app <- function(daf,
                    name = "data",
                    cells_daf = NULL,
                    atlas = NULL,
                    config_file = NULL,
                    port = NULL,
                    host = NULL,
                    launch.browser = FALSE,
                    profile = FALSE,
                    cache_in_daf = NULL,
                    ...) {
    timed <- function(label, expr) {
        if (!isTRUE(profile)) {
            return(force(expr))
        }
        start <- proc.time()[["elapsed"]]
        on.exit(
            {
                elapsed <- proc.time()[["elapsed"]] - start
                cli::cli_alert_info("Timing {label}: {round(elapsed, digits = 2)}s")
            },
            add = TRUE
        )
        force(expr)
    }

    # Initialize environment
    init_mcview_env()

    # Force a usable bitmap device for ggplotly() panel-size measurement.
    # callr/Docker/plain Rscript don't load ~/.Rprofile, so the child can
    # inherit `bitmapType = "Xlib"` while the R build lacks X11 — every
    # renderPlotly() then errors with "unable to start device PNG".
    .ensure_safe_bitmap_type()

    # Enforce dafr memory-mapping before opening any DAF
    .assert_dafr_mmap()

    # Handle path input - convert to DAF object
    if (is.character(daf)) {
        daf <- load_daf_from_path(daf)
    }

    # Initialize based on input type
    if (dafr::is_daf(daf)) {
        # Single DAF object
        timed("init_single_daf_mode", init_single_daf_mode(daf, name, config_file, profile, cache_in_daf = cache_in_daf))
    } else if (is.list(daf) && length(daf) > 0 && all(sapply(daf, dafr::is_daf))) {
        # Multiple DAF objects (must be named)
        if (is.null(names(daf)) || any(names(daf) == "")) {
            cli_abort("When providing multiple DAF objects, the list must be named")
        }
        timed("init_multi_daf_mode", init_multi_daf_mode(daf, config_file, profile, cache_in_daf = cache_in_daf))
    } else {
        cli_abort("daf must be a DAF object, a named list of DAF objects, or a path to DAF data")
    }

    # Set cells DAF if explicitly provided (overrides auto-detection)
    if (!is.null(cells_daf)) {
        cells_daf_path <- if (is.character(cells_daf)) cells_daf else NULL
        if (!is.null(cells_daf_path) && dir.exists(cells_daf_path)) {
            dataset_name <- if (dafr::is_daf(daf)) name else names(daf)[1]
            set_cells_daf(dataset_name, cells_daf_path)
            cli::cli_alert_success("Cells DAF set for '{dataset_name}': {cells_daf_path}")
        }
    }

    # Initialize atlas if provided
    if (!is.null(atlas)) {
        if (is.character(atlas)) {
            atlas <- load_daf_from_path(atlas)
        }
        init_atlas(atlas)
    }

    timed("init_defs", init_defs())
    validate_mcview_env()
    timed("config_shiny_cache", config_shiny_cache())
    # Defer future::plan(multisession) until first use (Flow tab).
    # This saves ~5.7s at startup since spawning worker processes is expensive
    # and only the Flow tab (plot_vein) uses future_promise.
    # See ensure_future_plan() below.

    # Register onStop hook for orderly shutdown
    shiny::onStop(function() {
        mcview_on_stop()
    })

    with_golem_options(
        app = shinyApp(
            ui = app_ui,
            server = app_server,
            options = list(port = port, host = host, launch.browser = launch.browser)
        ),
        golem_opts = list(...)
    )
}

#' Load DAF from file path
#'
#' @param path Path to DAF directory or H5 file
#' @return DAF object
#' @noRd
load_daf_from_path <- function(path) {
    if (!file.exists(path) && !dir.exists(path)) {
        cli_abort("Path does not exist: {.path {path}}")
    }

    # h5ad (AnnData-shaped HDF5)
    if (grepl("\\.h5ad$|\\.h5$", path, ignore.case = TRUE)) {
        cli_alert_info("Loading DAF from h5ad file: {.path {path}}")
        return(dafr::h5ad_as_daf(path))
    }

    # files_daf directory
    if (dir.exists(path)) {
        if (file.exists(file.path(path, "daf.json"))) {
            cli_alert_info("Loading DAF from directory: {.path {path}}")
            return(dafr::files_daf(path, mode = "r"))
        }
        cli_abort("Directory {.path {path}} does not appear to be a valid DAF (missing daf.json)")
    }

    cli_abort("Could not determine DAF format for: {.path {path}}")
}

#' Orderly shutdown callback for MCView
#'
#' Called by shiny::onStop() when the app is shutting down.
#' Resets the future plan and cleans up MCView state.
#'
#' @noRd
mcview_on_stop <- function() {
    cli::cli_alert_info("MCView shutting down...")

    # Reset future plan to sequential to avoid orphaned worker processes
    # (only if the multisession plan was actually activated)
    if (.future_plan_state$initialized) {
        tryCatch(
            {
                future::plan(future::sequential)
                .future_plan_state$initialized <- FALSE
            },
            error = function(e) {
                # Suppress errors during shutdown
            }
        )
    }

    # Clean up MCView environment (DAF caches, in-memory state)
    tryCatch(
        cleanup_mcview_env(),
        error = function(e) {
            # Suppress errors during shutdown
        }
    )

    cli::cli_alert_success("MCView shutdown complete")
}

config_shiny_cache <- function() {
    config <- mcv_get("config")
    max_size <- config$shiny_cache_max_size %||% 200e6

    if (!is.null(config$shiny_cache_dir)) {
        # Generate cache dir name based on config title or random
        cache_name <- gsub("[^a-zA-Z0-9]", "_", config$title %||% "mcview")
        if (is.logical(config$shiny_cache_dir) && config$shiny_cache_dir) {
            shiny_cache_dir <- tempfile(paste0("shiny_cache_", cache_name), tmpdir = tempdir())
        } else {
            shiny_cache_dir <- tempfile(paste0("shiny_cache_", cache_name), tmpdir = config$shiny_cache_dir)
        }
        shinyOptions(cache = cachem::cache_disk(shiny_cache_dir, max_size = max_size))
    } else {
        shinyOptions(cache = cachem::cache_mem(max_size = max_size))
    }
}
