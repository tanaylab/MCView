verify_config_file <- function(config) {
    required_fields <- c(
        "title",
        "tabs"
    )

    for (field in required_fields) {
        if (!(field %in% names(config))) {
            cli_abort("The field {.field {field}} does not exist in the config file.")
        }
    }

    # TODO: verify datasets config

    invisible(TRUE)
}


#' Initialize application config file
#'
#' @param project path of the project to initialize
#'
#' @noRd
init_config <- function(project, profile = FALSE) {
    config_file <- project_config_file(project)
    mcv_set("project", project)
    mcv_set("about_file", fs::path_abs(project_about_file(project)))
    mcv_set("cache_dir", project_cache_dir(project))

    config <- yaml::read_yaml(config_file)
    verify_config_file(config)
    if (is.null(config$light_version)) {
        config$light_version <- FALSE
    }

    if (file.exists(project_metacells_algorithm_file(project))) {
        metacells_version <- readLines(project_metacells_algorithm_file(project))
        config$metacells_version <- metacells_version
    }
    config$profile <- profile
    mcv_set("config", config)
}

init_defs <- function() {
    mcv_set("egc_epsilon", 1e-5)
    options(spinner.type = 6)

    theme_set(theme_classic())

    init_tab_defs()

    if (isTRUE(app_config("profile"))) {
        start <- proc.time()[["elapsed"]]
        init_selected_genes()
        elapsed <- proc.time()[["elapsed"]] - start
        cli::cli_alert_info("Timing init_selected_genes: {round(elapsed, digits = 2)}s")
    } else {
        init_selected_genes()
    }

    mcv_set("expr_breaks", c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1))
}

#' Get all possible tab names from tab definitions
#'
#' @return Character vector of all tab names
get_all_tab_names <- function() {
    tab_defs <- mcv_get("tab_defs")
    if (is.null(tab_defs)) {
        # Fallback if tab_defs not yet initialized
        return(c(
            "About", "QC", "Projection QC", "Manifold", "Genes", "Markers",
            "Diff. Expression", "Gene modules", "Cell types", "Annotate",
            "Inner-fold", "Stdev-fold", "Outliers", "Projected-fold",
            "Flow", "Samples", "Query", "Atlas"
        ))
    }
    names(tab_defs)
}

init_selected_genes <- function() {
    config <- mcv_get("config")
    # if selected genes are not set - choose them from the first dataset
    if (is.null(config$selected_gene1) || is.null(config$selected_gene2)) {
        cache_enabled <- isTRUE(config$cache_in_daf)
        if (cache_enabled) {
            daf_obj <- get_dataset_daf(dataset_ls()[1])
            cached_gene1 <- daf_scalar(daf_obj, "mcview_cache_default_gene1", default = NULL)
            cached_gene2 <- daf_scalar(daf_obj, "mcview_cache_default_gene2", default = NULL)
            if (!is.null(cached_gene1) && !is.null(cached_gene2)) {
                genes <- c(cached_gene1, cached_gene2)
            }
        }

        if (!exists("genes", inherits = FALSE)) {
            mc_egc <- get_mc_egc(dataset_ls()[1])
            if (has_atlas(dataset_ls()[1])) {
                mc_egc_atlas <- get_mc_egc(dataset_ls()[1], atlas = TRUE)
                mc_egc <- mc_egc[intersect(rownames(mc_egc), rownames(mc_egc_atlas)), ]
            }

            f <- rep(TRUE, nrow(mc_egc))
            lateral_genes <- get_mc_data(dataset_ls()[1], "lateral_genes")
            if (!is.null(lateral_genes)) {
                f <- f & !(rownames(mc_egc) %in% lateral_genes)
            }
            noisy_genes <- get_mc_data(dataset_ls()[1], "noisy_genes")
            if (!is.null(noisy_genes)) {
                f <- f & !(rownames(mc_egc) %in% noisy_genes)
            }
            mc_egc <- mc_egc[f, ]

            minmax <- matrixStats::rowMaxs(mc_egc, na.rm = TRUE) - matrixStats::rowMins(mc_egc, na.rm = TRUE)
            names(minmax) <- rownames(mc_egc)
            genes <- names(head(sort(minmax, decreasing = TRUE), n = 2))

            if (cache_enabled) {
                tryCatch(
                    {
                        daf_obj <- get_dataset_daf(dataset_ls()[1])
                        dafr::set_scalar(daf_obj, "mcview_cache_default_gene1", genes[1])
                        dafr::set_scalar(daf_obj, "mcview_cache_default_gene2", genes[2])
                    },
                    error = function(e) NULL
                )
            }
        }
    }

    mcv_set("default_gene1", config$selected_gene1 %||% genes[1])
    mcv_set("default_gene2", config$selected_gene2 %||% genes[2])
}


init_tab_defs <- function() {
    tab_defs <- list(
        "About" = list(
            title = "About",
            module_name = "about",
            icon = "info"
        ),
        "QC" = list(
            title = "QC",
            module_name = "qc",
            icon = "check"
        ),
        "Projection QC" = list(
            title = "Projection QC",
            module_name = "projection_qc",
            icon = "clipboard-check"
        ),
        "Manifold" = list(
            title = "Manifold",
            module_name = "manifold",
            icon = "project-diagram"
        ),
        "Genes" = list(
            title = "Genes",
            module_name = "gene_mc",
            icon = "wind"
        ),
        "Markers" = list(
            title = "Markers",
            module_name = "markers",
            icon = "map-marker"
        ),
        "Diff. Expression" = list(
            title = "Diff. Expression",
            module_name = "mc_mc",
            icon = "chart-bar"
        ),
        "Gene modules" = list(
            title = "Gene modules",
            module_name = "gene_modules",
            icon = "layer-group"
        ),
        "Cell types" = list(
            title = "Cell types",
            module_name = "cell_type",
            icon = "bacteria"
        ),
        "Annotate" = list(
            title = "Annotate",
            module_name = "annotate",
            icon = "pen"
        ),
        "Inner-fold" = list(
            title = "Inner-fold",
            module_name = "inner_fold",
            icon = "cloud-sun-rain"
        ),
        "Stdev-fold" = list(
            title = "Stdev-fold",
            module_name = "stdev_fold",
            icon = "cloud-sun-rain"
        ),
        "Outliers" = list(
            title = "Outliers",
            module_name = "outliers",
            icon = "user-astronaut"
        ),
        "Projected-fold" = list(
            title = "Projected-fold",
            module_name = "proj_fold",
            icon = "cloud-moon-rain"
        ),
        "Flow" = list(
            title = "Flow",
            module_name = "flow",
            icon = "water"
        ),
        "Samples" = list(
            title = "Samples",
            module_name = "samples",
            icon = "vials"
        ),
        "Query" = list(
            title = "Query",
            module_name = "query",
            icon = "video"
        ),
        "Atlas" = list(
            title = "Atlas",
            module_name = "atlas",
            icon = "atlas"
        )
    )
    mcv_set("tab_defs", tab_defs)

    default_tabs <- c("About", "Genes", "Diff. Expression")

    cur_config <- mcv_get("config")

    if (!is.null(cur_config$tabs)) {
        cur_config$original_tabs <- cur_config$tabs
        cur_config$tabs[cur_config$tabs == "Metacells"] <- "Diff. Expression" # here for backward compatibility
        cur_config$tabs <- cur_config$tabs[cur_config$tabs != "Metadata"] # ignore "Metadata" for backward compatibility
        for (.x in cur_config$tabs) {
            if (!(.x %in% names(mcv_get("tab_defs")))) {
                cli_warn("{.x} is not a valid tab name. Update `tabs` in your configuration file.")
                cur_config$tabs <- cur_config$tabs[cur_config$tabs != .x]
            }
        }
    } else {
        cur_config$tabs <- default_tabs
    }

    if (!is.null(cur_config$excluded_tabs)) {
        tab_defs <- mcv_get("tab_defs")
        mcv_set("tab_defs", tab_defs[!(names(tab_defs) %in% cur_config$excluded_tabs)])
        cur_config$tabs <- cur_config$tabs[!(cur_config$tabs %in% cur_config$excluded_tabs)]
    }

    if (!rmarkdown::pandoc_available() && "About" %in% cur_config$tabs) {
        warning("pandoc is not available, removing 'About' tab'")
        cur_config$tabs <- cur_config$tabs[cur_config$tabs != "About"]
    }

    if (!is.null(cur_config$light_version) && cur_config$light_version) {
        # make the About tab first if exists
        if ("About" %in% cur_config$tabs) {
            cur_config$tabs <- c("About", cur_config$tabs[cur_config$tabs != "About"])
        }
    }

    mcv_set("config", cur_config)
}

order_tabs <- function(tabs) {
    tabs_order <- c("QC", "Projection QC", "Manifold", "Genes", "Query", "Atlas", "Markers", "Gene modules", "Projected-fold", "Diff. Expression", "Samples", "Cell types", "Annotate")

    # order the tabs according to the order in tabs_order. Tabs that are not in tabs_order will be added at the end
    new_tabs <- c(tabs_order[tabs_order %in% tabs], tabs[!(tabs %in% tabs_order)])

    if ("About" %in% tabs) {
        new_tabs <- c(setdiff(new_tabs, "About"), "About")
    }

    return(new_tabs)
}

app_config <- function(key = NULL) {
    config <- mcv_get("config")
    if (is.null(key)) {
        return(config)
    }
    config[[key]]
}

#' Get gene names from DAF
#'
#' @param dataset Dataset name
#' @param atlas Whether to use atlas data
#' @return Character vector of gene names
gene_names <- function(dataset, atlas = FALSE) {
    # Get DAF object
    if (atlas) {
        daf_obj <- get_atlas_daf()
    } else {
        daf_obj <- get_dataset_daf(dataset)
    }

    if (is.null(daf_obj)) {
        return(character(0))
    }

    # Use DAF axis entries instead of downloading full matrix
    dafr::axis_entries(daf_obj, "gene")
}

gene_names_label <- function(dataset, atlas = FALSE, gene_modules = NULL) {
    gene_names <- gene_names(dataset, atlas)
    if (is.null(gene_names)) {
        return(gene_names)
    }
    names(gene_names) <- gene_label(gene_names, dataset, gene_modules)
    return(gene_names)
}


# ==============================================================================
# Cache Configuration
# ==============================================================================

#' Extract cache configuration
#'
#' Parses cache configuration from multiple sources with priority:
#' 1. YAML config file (highest priority)
#' 2. DAF scalars
#' 3. Defaults (lowest priority)
#'
#' @param config Full configuration list (from YAML or existing config)
#' @param daf_obj Optional DAF object to read scalars from
#' @param dataset_name Optional dataset name for per-dataset overrides
#'
#' @return Cache configuration list with fields:
#'   - enabled: Whether caching is enabled
#'   - type: "files" or "memory"
#'   - cache_dir: Directory for files cache (relative path)
#'   - precompute_on_startup: Whether to populate cache at startup
#'   - invalidation: List with strategy and version_scalar
#'
#' @export
extract_cache_config <- function(config = NULL, daf_obj = NULL, dataset_name = NULL) {
    # Default cache configuration
    defaults <- list(
        enabled = FALSE,
        type = "files",
        cache_dir = ".mcview_cache",
        precompute_on_startup = TRUE,
        invalidation = list(
            strategy = "version",
            version_scalar = "mcview_data_version"
        )
    )

    cache_config <- defaults

    # Priority 3: DAF scalars (lowest priority for cache-specific settings)
    if (!is.null(daf_obj)) {
        daf_cache_enabled <- daf_scalar(daf_obj, "mcview_cache_enabled", default = NULL)
        daf_cache_type <- daf_scalar(daf_obj, "mcview_cache_type", default = NULL)
        daf_cache_dir <- daf_scalar(daf_obj, "mcview_cache_dir", default = NULL)
        daf_cache_precompute <- daf_scalar(daf_obj, "mcview_cache_precompute", default = NULL)

        # Also support legacy scalars
        daf_cache_in_daf <- daf_scalar(daf_obj, "mcview_cache_in_daf", default = NULL)
        daf_cache_daf_root <- daf_scalar(daf_obj, "mcview_cache_daf_root", default = NULL)

        if (!is.null(daf_cache_enabled)) {
            cache_config$enabled <- isTRUE(normalize_cache_flag(daf_cache_enabled))
        } else if (!is.null(daf_cache_in_daf)) {
            # Legacy support
            cache_config$enabled <- isTRUE(normalize_cache_flag(daf_cache_in_daf))
        }

        if (!is.null(daf_cache_type) && daf_cache_type %in% c("files", "memory")) {
            cache_config$type <- daf_cache_type
        }

        if (!is.null(daf_cache_dir) && nzchar(daf_cache_dir)) {
            cache_config$cache_dir <- daf_cache_dir
        } else if (!is.null(daf_cache_daf_root) && nzchar(daf_cache_daf_root)) {
            # Legacy support
            cache_config$cache_dir <- daf_cache_daf_root
        }

        if (!is.null(daf_cache_precompute)) {
            cache_config$precompute_on_startup <- isTRUE(normalize_cache_flag(daf_cache_precompute))
        }
    }

    # Priority 2 & 1: Config from YAML (overrides DAF scalars)
    if (!is.null(config)) {
        # Support new nested cache config structure
        if (!is.null(config$cache) && is.list(config$cache)) {
            yaml_cache <- config$cache

            if (!is.null(yaml_cache$enabled)) {
                cache_config$enabled <- isTRUE(normalize_cache_flag(yaml_cache$enabled))
            }
            if (!is.null(yaml_cache$type) && yaml_cache$type %in% c("files", "memory")) {
                cache_config$type <- yaml_cache$type
            }
            if (!is.null(yaml_cache$cache_dir) && nzchar(yaml_cache$cache_dir)) {
                cache_config$cache_dir <- yaml_cache$cache_dir
            }
            if (!is.null(yaml_cache$precompute_on_startup)) {
                cache_config$precompute_on_startup <- isTRUE(normalize_cache_flag(yaml_cache$precompute_on_startup))
            }
            if (!is.null(yaml_cache$invalidation) && is.list(yaml_cache$invalidation)) {
                cache_config$invalidation <- modifyList(cache_config$invalidation, yaml_cache$invalidation)
            }

            # Per-dataset overrides
            if (!is.null(dataset_name) && !is.null(config$datasets[[dataset_name]]$cache)) {
                ds_cache <- config$datasets[[dataset_name]]$cache
                if (!is.null(ds_cache$enabled)) {
                    cache_config$enabled <- isTRUE(normalize_cache_flag(ds_cache$enabled))
                }
                if (!is.null(ds_cache$type) && ds_cache$type %in% c("files", "memory")) {
                    cache_config$type <- ds_cache$type
                }
                if (!is.null(ds_cache$cache_dir) && nzchar(ds_cache$cache_dir)) {
                    cache_config$cache_dir <- ds_cache$cache_dir
                }
                if (!is.null(ds_cache$precompute_on_startup)) {
                    cache_config$precompute_on_startup <- isTRUE(normalize_cache_flag(ds_cache$precompute_on_startup))
                }
            }
        }

        # Also support legacy flat config structure for backward compatibility
        if (!is.null(config$cache_in_daf)) {
            cache_config$enabled <- isTRUE(normalize_cache_flag(config$cache_in_daf))
        }
        if (!is.null(config$cache_daf_root) && nzchar(config$cache_daf_root)) {
            cache_config$cache_dir <- config$cache_daf_root
        }
    }

    cache_config
}

#' Create default cache configuration
#'
#' @param enabled Whether caching is enabled
#' @param type Cache type: "files" or "memory"
#' @param cache_dir Directory for files cache
#' @param precompute_on_startup Whether to populate cache at startup
#'
#' @return Cache configuration list
#' @export
create_cache_config <- function(
  enabled = TRUE,
  type = "files",
  cache_dir = ".mcview_cache",
  precompute_on_startup = TRUE
) {
    list(
        enabled = enabled,
        type = type,
        cache_dir = cache_dir,
        precompute_on_startup = precompute_on_startup,
        invalidation = list(
            strategy = "version",
            version_scalar = "mcview_data_version"
        )
    )
}

#' Access files in the current app
#'
#' @param ... Character vector specifying directory and or file to
#'     point to inside the current package.
#'
#' @noRd
app_sys <- function(...) {
    system.file(..., package = "MCView")
}
