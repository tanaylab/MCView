#!/usr/bin/env Rscript
# convert_mcview_app.R
#
# Converts an old-format MCView app directory (cache/ + config/) to a new
# DAF-backed app directory suitable for MCView::run_app().
#
# Usage:
#   Rscript inst/scripts/convert_mcview_app.R \
#     --input /path/to/old/app \
#     --output /path/to/new/app \
#     [--format files|h5] \
#     [--dataset NAME] \
#     [--dry-run] \
#     [--verbose]

# ==============================================================================
# Messaging helpers (cli if available, message() fallback)
# ==============================================================================
has_cli <- requireNamespace("cli", quietly = TRUE)

msg_info <- if (has_cli) cli::cli_alert_info else function(...) message("[INFO] ", ...)
msg_ok <- if (has_cli) cli::cli_alert_success else function(...) message("[OK]   ", ...)
msg_warn <- if (has_cli) cli::cli_alert_warning else function(...) message("[WARN] ", ...)
msg_err <- if (has_cli) cli::cli_alert_danger else function(...) message("[ERR]  ", ...)
msg_bullet <- if (has_cli) cli::cli_li else function(...) message("  - ", ...)

stop_msg <- function(...) {
    msg <- paste0(...)
    if (has_cli) {
        cli::cli_abort(msg)
    } else {
        stop(msg, call. = FALSE)
    }
}

# ==============================================================================
# CLI argument parsing
# ==============================================================================
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
    opts <- list(
        input = NULL,
        output = NULL,
        format = "files",
        dataset = NULL,
        dry_run = FALSE,
        verbose = FALSE
    )

    i <- 1
    while (i <= length(args)) {
        arg <- args[i]
        if (arg == "--input" && i < length(args)) {
            i <- i + 1
            opts$input <- args[i]
        } else if (arg == "--output" && i < length(args)) {
            i <- i + 1
            opts$output <- args[i]
        } else if (arg == "--format" && i < length(args)) {
            i <- i + 1
            opts$format <- args[i]
        } else if (arg == "--dataset" && i < length(args)) {
            i <- i + 1
            opts$dataset <- args[i]
        } else if (arg == "--dry-run") {
            opts$dry_run <- TRUE
        } else if (arg == "--verbose") {
            opts$verbose <- TRUE
        } else if (arg == "--help" || arg == "-h") {
            cat(
                "Usage: Rscript convert_mcview_app.R [options]\n",
                "\n",
                "Options:\n",
                "  --input PATH      Path to existing MCView app directory (required)\n",
                "  --output PATH     Path to output DAF app directory (required)\n",
                "  --format FORMAT   \"files\" (default) or \"h5\"\n",
                "  --dataset NAME    Dataset name override (default: auto-detect)\n",
                "  --dry-run         Only analyze, don't convert\n",
                "  --verbose         Detailed progress\n",
                "  --help, -h        Show this help message\n",
                sep = ""
            )
            quit(status = 0)
        } else {
            stop_msg("Unknown argument: ", arg, "\nUse --help for usage information.")
        }
        i <- i + 1
    }

    # Validate required arguments
    if (is.null(opts$input)) {
        stop_msg("--input is required. Use --help for usage information.")
    }
    if (is.null(opts$output)) {
        stop_msg("--output is required. Use --help for usage information.")
    }
    if (!opts$format %in% c("files", "h5")) {
        stop_msg("--format must be \"files\" or \"h5\", got: ", opts$format)
    }

    opts
}

# ==============================================================================
# Detect old app structure
# ==============================================================================
detect_project_dir <- function(input_path) {
    if (!dir.exists(input_path)) {
        stop_msg("Input directory does not exist: ", input_path)
    }

    # Check for nested project/ subdirectory layout first
    nested_dirs <- list.dirs(input_path, recursive = FALSE, full.names = FALSE)
    for (subdir in nested_dirs) {
        candidate <- file.path(input_path, subdir)
        config_dir <- file.path(candidate, "config")
        cache_dir <- file.path(candidate, "cache")
        if (dir.exists(config_dir) && dir.exists(cache_dir)) {
            msg_info("Detected nested project layout: {.path {candidate}}")
            return(candidate)
        }
    }

    # Check for flat layout (config/ and cache/ directly under input)
    config_dir <- file.path(input_path, "config")
    cache_dir <- file.path(input_path, "cache")
    if (dir.exists(config_dir) && dir.exists(cache_dir)) {
        msg_info("Detected flat project layout: {.path {input_path}}")
        return(input_path)
    }

    stop_msg(
        "Could not detect MCView project structure in: ", input_path, "\n",
        "Expected config/ and cache/ subdirectories (either directly or inside a project/ subdirectory)."
    )
}

# ==============================================================================
# Read old config
# ==============================================================================
read_old_config <- function(project_dir) {
    config_file <- file.path(project_dir, "config", "config.yaml")
    if (!file.exists(config_file)) {
        msg_warn("No config.yaml found at {.path {config_file}}")
        return(list())
    }

    if (!requireNamespace("yaml", quietly = TRUE)) {
        stop_msg("Package 'yaml' is required to read config.yaml")
    }

    config <- yaml::read_yaml(config_file)
    msg_info("Read configuration from {.path {config_file}}")
    config
}

# ==============================================================================
# Detect datasets in cache
# ==============================================================================
detect_datasets <- function(project_dir, dataset_override = NULL) {
    cache_dir <- file.path(project_dir, "cache")
    if (!dir.exists(cache_dir)) {
        stop_msg("Cache directory not found: ", cache_dir)
    }

    if (!is.null(dataset_override)) {
        dataset_path <- file.path(cache_dir, dataset_override)
        if (!dir.exists(dataset_path)) {
            stop_msg("Specified dataset not found in cache: ", dataset_override)
        }
        return(dataset_override)
    }

    # Auto-detect datasets (subdirectories, excluding "atlas")
    all_dirs <- list.dirs(cache_dir, recursive = FALSE, full.names = FALSE)
    datasets <- all_dirs[all_dirs != "atlas"]

    if (length(datasets) == 0) {
        stop_msg("No dataset directories found in cache: ", cache_dir)
    }

    datasets
}

# ==============================================================================
# Detect atlases
# ==============================================================================
detect_atlases <- function(project_dir, datasets) {
    cache_dir <- file.path(project_dir, "cache")
    atlases <- character(0)
    for (ds in datasets) {
        atlas_dir <- file.path(cache_dir, ds, "atlas")
        if (dir.exists(atlas_dir)) {
            atlases <- c(atlases, ds)
        }
    }
    atlases
}

# ==============================================================================
# Load MCView package
# ==============================================================================
load_mcview <- function(verbose = FALSE) {
    # Try loading as installed package
    if (requireNamespace("MCView", quietly = TRUE)) {
        if (verbose) msg_info("Using installed MCView package")
        return(TRUE)
    }

    # Fall back to devtools::load_all for development
    # Look for DESCRIPTION file to find package root
    script_dir <- getwd()
    candidates <- c(
        script_dir,
        dirname(dirname(script_dir)),  # inst/scripts -> package root
        Sys.getenv("MCVIEW_ROOT", unset = "")
    )
    candidates <- candidates[nzchar(candidates)]

    for (candidate in candidates) {
        if (file.exists(file.path(candidate, "DESCRIPTION"))) {
            if (!requireNamespace("devtools", quietly = TRUE)) {
                stop_msg(
                    "MCView is not installed and devtools is not available.\n",
                    "Install MCView or run from within the package directory."
                )
            }
            if (verbose) msg_info("Loading MCView via devtools::load_all from {.path {candidate}}")
            devtools::load_all(candidate, quiet = !verbose)
            return(TRUE)
        }
    }

    stop_msg(
        "Could not load MCView package.\n",
        "Either install it (install.packages/remotes::install_github) or run from the package source directory."
    )
}

# ==============================================================================
# Build transformed config for the new app
# ==============================================================================
build_new_config <- function(old_config) {
    # Keep only the settings relevant to the new DAF app
    keep_keys <- c(
        "title", "tabs", "selected_gene1", "selected_gene2",
        "selected_mc1", "selected_mc2", "light_version",
        "excluded_tabs", "cache_in_daf", "cache_daf_root",
        "shiny_cache_dir", "shiny_cache_max_size",
        "cells_daf"
    )

    new_config <- old_config[intersect(names(old_config), keep_keys)]

    # Remove old-format keys that don't apply
    new_config$datasets <- NULL
    new_config$help <- NULL

    new_config
}

# ==============================================================================
# Generate app.R content
# ==============================================================================
generate_app_r <- function(datasets, atlases, format) {
    lines <- c(
        "# MCView DAF App - converted from old format",
        paste0("# Generated by convert_mcview_app.R on ", Sys.time()),
        "",
        'options("golem.app.prod" = TRUE)',
        "library(MCView)",
        ""
    )

    daf_ext <- if (format == "h5") ".h5" else ""
    open_fn <- if (format == "h5") "dafr::h5df" else "dafr::files_daf"

    if (length(datasets) == 1) {
        ds <- datasets[1]
        daf_name <- paste0(ds, daf_ext)
        has_atlas <- ds %in% atlases

        lines <- c(lines,
            paste0('daf_path <- file.path(getwd(), "', daf_name, '")'),
            paste0("daf <- ", open_fn, "(daf_path)")
        )

        if (has_atlas) {
            atlas_name <- paste0(ds, "_atlas", daf_ext)
            lines <- c(lines,
                "",
                paste0('atlas_path <- file.path(getwd(), "', atlas_name, '")'),
                paste0("atlas <- if (", if (format == "h5") "file.exists" else "dir.exists", "(atlas_path)) ", open_fn, "(atlas_path) else NULL")
            )
        }

        lines <- c(lines,
            "",
            paste0(
                'MCView::run_app(daf, name = "', ds, '"',
                if (has_atlas) ", atlas = atlas" else "",
                ', config_file = "config/config.yaml")'
            )
        )
    } else {
        # Multi-dataset
        lines <- c(lines,
            "daf_list <- list()"
        )
        for (ds in datasets) {
            daf_name <- paste0(ds, daf_ext)
            lines <- c(lines,
                paste0('daf_list[["', ds, '"]] <- ', open_fn, '(file.path(getwd(), "', daf_name, '"))')
            )
        }

        # Check for atlases in multi-dataset mode
        atlas_datasets <- intersect(datasets, atlases)
        if (length(atlas_datasets) > 0) {
            lines <- c(lines, "")
            # Use the first atlas found
            ds <- atlas_datasets[1]
            atlas_name <- paste0(ds, "_atlas", daf_ext)
            lines <- c(lines,
                paste0('atlas_path <- file.path(getwd(), "', atlas_name, '")'),
                paste0("atlas <- if (", if (format == "h5") "file.exists" else "dir.exists", "(atlas_path)) ", open_fn, "(atlas_path) else NULL")
            )
            lines <- c(lines,
                "",
                'MCView::run_app(daf_list, atlas = atlas, config_file = "config/config.yaml")'
            )
        } else {
            lines <- c(lines,
                "",
                'MCView::run_app(daf_list, config_file = "config/config.yaml")'
            )
        }
    }

    paste(lines, collapse = "\n")
}

# ==============================================================================
# Embed config into DAF after conversion
# ==============================================================================
embed_config_in_daf <- function(output_dir, datasets, format, old_config, verbose = FALSE) {
    daf_ext <- if (format == "h5") ".h5" else ""

    for (ds in datasets) {
        daf_path <- file.path(output_dir, paste0(ds, daf_ext))

        if (!file.exists(daf_path) && !dir.exists(daf_path)) {
            msg_warn("DAF not found at {.path {daf_path}}, skipping config embedding")
            next
        }

        if (verbose) msg_info("Embedding config in DAF: {.path {daf_path}}")

        daf <- if (format == "h5") {
            dafr::h5df(daf_path, mode = "w+")
        } else {
            dafr::files_daf(daf_path, mode = "w+")
        }

        # Embed title
        title <- old_config$title
        if (!is.null(title) && nzchar(title)) {
            tryCatch(
                dafr::set_scalar(daf, "mcview_title", title),
                error = function(e) {
                    if (verbose) msg_warn("Could not set mcview_title: {conditionMessage(e)}")
                }
            )
        }

        # Embed tabs
        tabs <- old_config$tabs
        if (!is.null(tabs) && length(tabs) > 0) {
            tabs_str <- if (length(tabs) > 1) paste(tabs, collapse = ",") else tabs
            tryCatch(
                dafr::set_scalar(daf, "mcview_tabs", tabs_str),
                error = function(e) {
                    if (verbose) msg_warn("Could not set mcview_tabs: {conditionMessage(e)}")
                }
            )
        }

        # Embed excluded_tabs
        excluded_tabs <- old_config$excluded_tabs
        if (!is.null(excluded_tabs) && length(excluded_tabs) > 0) {
            excluded_str <- if (length(excluded_tabs) > 1) paste(excluded_tabs, collapse = ",") else excluded_tabs
            tryCatch(
                dafr::set_scalar(daf, "mcview_excluded_tabs", excluded_str),
                error = function(e) {
                    if (verbose) msg_warn("Could not set mcview_excluded_tabs: {conditionMessage(e)}")
                }
            )
        }

        # Embed light_version
        if (!is.null(old_config$light_version)) {
            tryCatch(
                dafr::set_scalar(daf, "mcview_light_version", old_config$light_version),
                error = function(e) {
                    if (verbose) msg_warn("Could not set mcview_light_version: {conditionMessage(e)}")
                }
            )
        }

        if (verbose) msg_ok("Config embedded in {.path {daf_path}}")
    }
}

# ==============================================================================
# Copy about.Rmd if present
# ==============================================================================
copy_about_rmd <- function(project_dir, output_dir, verbose = FALSE) {
    about_src <- file.path(project_dir, "config", "about.Rmd")
    if (file.exists(about_src)) {
        about_dst <- file.path(output_dir, "config", "about.Rmd")
        file.copy(about_src, about_dst, overwrite = TRUE)
        if (verbose) msg_info("Copied about.Rmd")
        return(TRUE)
    }
    FALSE
}

# ==============================================================================
# Print conversion summary
# ==============================================================================
print_summary <- function(opts, project_dir, datasets, atlases, old_config, dry_run = FALSE) {
    cat("\n")
    if (dry_run) {
        msg_info("=== DRY RUN - Analysis Only ===")
    } else {
        msg_ok("=== Conversion Complete ===")
    }
    cat("\n")
    msg_info("Input:      {.path {opts$input}}")
    msg_info("Project:    {.path {project_dir}}")
    if (!dry_run) {
        msg_info("Output:     {.path {opts$output}}")
    }
    msg_info("Format:     {opts$format}")
    msg_info("Datasets:   {paste(datasets, collapse = ', ')}")
    if (length(atlases) > 0) {
        msg_info("Atlases:    {paste(atlases, collapse = ', ')}")
    }
    if (!is.null(old_config$title)) {
        msg_info("Title:      {old_config$title}")
    }
    if (!is.null(old_config$tabs)) {
        msg_info("Tabs:       {paste(old_config$tabs, collapse = ', ')}")
    }
    cat("\n")

    if (!dry_run) {
        msg_ok("To run the converted app:")
        cat("  cd ", opts$output, "\n")
        cat("  Rscript app.R\n\n")
    }
}

# ==============================================================================
# Main conversion logic
# ==============================================================================
main <- function() {
    opts <- parse_args()

    msg_info("MCView App Converter")
    cat("\n")

    # Step 1: Detect project structure
    if (opts$verbose) msg_info("Detecting project structure...")
    project_dir <- detect_project_dir(opts$input)

    # Step 2: Read old config
    if (opts$verbose) msg_info("Reading configuration...")
    old_config <- read_old_config(project_dir)

    # Step 3: Detect datasets
    if (opts$verbose) msg_info("Detecting datasets...")
    datasets <- detect_datasets(project_dir, opts$dataset)
    msg_info("Found {length(datasets)} dataset(s): {paste(datasets, collapse = ', ')}")

    # Step 4: Detect atlases
    atlases <- detect_atlases(project_dir, datasets)
    if (length(atlases) > 0) {
        msg_info("Found atlas data for: {paste(atlases, collapse = ', ')}")
    }

    # Dry-run: print summary and exit
    if (opts$dry_run) {
        print_summary(opts, project_dir, datasets, atlases, old_config, dry_run = TRUE)
        return(invisible(NULL))
    }

    # Step 5: Load MCView and initialize Julia/dafr
    if (opts$verbose) msg_info("Loading MCView package...")
    load_mcview(verbose = opts$verbose)

    if (opts$verbose) msg_info("Initializing Julia/dafr...")
    julia_home <- file.path(Sys.getenv("CONDA_PREFIX", unset = ""), "bin")
    if (nzchar(julia_home) && dir.exists(julia_home)) {
        options(dafr.JULIA_HOME = julia_home)
    }
    sysimage <- Sys.getenv("JULIA_SYSIMAGE", unset = "")
    if (nzchar(sysimage) && file.exists(sysimage)) {
        options(dafr.JULIA_SYSIMAGE = sysimage)
    }
    dafr::setup_daf(pkg_check = FALSE, julia_environment = "custom")

    # Step 6: Create output directory
    if (!dir.exists(opts$output)) {
        dir.create(opts$output, recursive = TRUE)
    }

    # Step 7: Run conversion
    msg_info("Converting datasets to DAF format...")
    daf_output_dir <- opts$output

    MCView::convert_project_to_daf(
        project = project_dir,
        output = daf_output_dir,
        format = opts$format,
        datasets = datasets
    )

    # Step 8: Embed config into DAF scalars
    if (opts$verbose) msg_info("Embedding configuration into DAF...")
    embed_config_in_daf(daf_output_dir, datasets, opts$format, old_config, verbose = opts$verbose)

    # Step 9: Create output app structure
    msg_info("Creating app structure...")
    config_out_dir <- file.path(opts$output, "config")
    if (!dir.exists(config_out_dir)) {
        dir.create(config_out_dir, recursive = TRUE)
    }

    # Write transformed config.yaml
    new_config <- build_new_config(old_config)
    config_out_file <- file.path(config_out_dir, "config.yaml")
    yaml::write_yaml(new_config, config_out_file)
    msg_ok("Wrote {.path {config_out_file}}")

    # Copy about.Rmd if present
    copy_about_rmd(project_dir, opts$output, verbose = opts$verbose)

    # Step 10: Generate app.R
    app_r_content <- generate_app_r(datasets, atlases, opts$format)
    app_r_file <- file.path(opts$output, "app.R")
    writeLines(app_r_content, app_r_file)
    msg_ok("Wrote {.path {app_r_file}}")

    # Summary
    print_summary(opts, project_dir, datasets, atlases, old_config)
}

# ==============================================================================
# Entry point
# ==============================================================================
main()
