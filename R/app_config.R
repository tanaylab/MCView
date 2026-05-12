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

#' Compute optimal default genes (expensive)
#'
#' Loads the full EGC matrix and finds the two genes with highest
#' expression range (max - min across metacells), excluding lateral
#' and noisy genes. This is the original init_selected_genes() logic
#' but extracted so it can be called lazily.
#'
#' @return Character vector of length 2 with gene names
#' @noRd
compute_optimal_default_genes <- function() {
    mc_egc <- get_mc_egc(dataset_names()[1])
    if (has_atlas(dataset_names()[1])) {
        mc_egc_atlas <- get_mc_egc(dataset_names()[1], atlas = TRUE)
        mc_egc <- mc_egc[intersect(rownames(mc_egc), rownames(mc_egc_atlas)), ]
    }

    f <- rep(TRUE, nrow(mc_egc))
    lateral_genes <- get_mc_data(dataset_names()[1], "lateral_genes")
    if (!is.null(lateral_genes)) {
        f <- f & !(rownames(mc_egc) %in% lateral_genes)
    }
    noisy_genes <- get_mc_data(dataset_names()[1], "noisy_genes")
    if (!is.null(noisy_genes)) {
        f <- f & !(rownames(mc_egc) %in% noisy_genes)
    }
    mc_egc <- mc_egc[f, ]

    minmax <- matrixStats::rowMaxs(mc_egc, na.rm = TRUE) - matrixStats::rowMins(mc_egc, na.rm = TRUE)
    names(minmax) <- rownames(mc_egc)
    genes <- names(head(sort(minmax, decreasing = TRUE), n = 2))

    # Cache for future use
    config <- mcv_get("config")
    cache_enabled <- isTRUE(config$cache_in_daf)
    if (cache_enabled) {
        tryCatch(
            {
                daf_obj <- get_dataset_daf(dataset_names()[1])
                dafr::set_scalar(daf_obj, "mcview_cache_default_gene1", genes[1])
                dafr::set_scalar(daf_obj, "mcview_cache_default_gene2", genes[2])
            },
            error = function(e) NULL
        )
    }

    genes
}

#' Get all possible tab names from tab definitions
#'
#' @return Character vector of all tab names
get_all_tab_names <- function() {
    tab_defs <- mcv_get("tab_defs")
    if (is.null(tab_defs)) {
        # Fallback if tab_defs not yet initialized.
        # MCVIEW_TAB_NAMES covers contract-backed tabs; add UI-only tabs too.
        return(unique(c(MCVIEW_TAB_NAMES, "Projected-fold", "Query")))
    }
    names(tab_defs)
}

init_selected_genes <- function() {
    config <- mcv_get("config")

    # If config explicitly provides selected genes, use them directly (no computation)
    if (!is.null(config$selected_gene1) && !is.null(config$selected_gene2)) {
        mcv_set("default_gene1", config$selected_gene1)
        mcv_set("default_gene2", config$selected_gene2)
        return(invisible(NULL))
    }

    # Try fast sources: cached scalars in DAF, then gene axis fallback
    genes <- NULL

    # Source 1: Check DAF cache scalars (single roundtrip each)
    daf_obj <- get_dataset_daf(dataset_names()[1])
    if (!is.null(daf_obj)) {
        cached_gene1 <- daf_scalar(daf_obj, "mcview_cache_default_gene1", default = NULL)
        cached_gene2 <- daf_scalar(daf_obj, "mcview_cache_default_gene2", default = NULL)
        if (!is.null(cached_gene1) && !is.null(cached_gene2)) {
            genes <- c(cached_gene1, cached_gene2)
        }
    }

    # Source 2: Use first 2 genes from the gene axis (very fast, no matrix load)
    if (is.null(genes) && !is.null(daf_obj)) {
        gene_axis <- dafr::axis_entries(daf_obj, "gene")
        if (length(gene_axis) >= 2) {
            genes <- gene_axis[1:2]
            # Mark that we should compute optimal genes lazily
            mcv_set("default_genes_need_optimization", TRUE)
        }
    }

    mcv_set("default_gene1", config$selected_gene1 %||% genes[1])
    mcv_set("default_gene2", config$selected_gene2 %||% genes[2])
}

#' Get default gene with lazy optimization
#'
#' Returns the current default gene. On first call where optimization
#' is pending, triggers async computation of optimal genes via
#' compute_optimal_default_genes() and updates the defaults.
#'
#' @param which Which gene (1 or 2)
#' @return Gene name string
#' @noRd
get_default_gene <- function(which = 1) {
    var_name <- paste0("default_gene", which)
    gene <- mcv_get(var_name)

    # If genes need optimization, compute them now (one-time cost on first access)
    if (isTRUE(mcv_get("default_genes_need_optimization"))) {
        mcv_set("default_genes_need_optimization", FALSE)
        cli::cli_alert_info("Computing optimal default genes (one-time)...")
        tryCatch(
            {
                optimal_genes <- compute_optimal_default_genes()
                if (length(optimal_genes) >= 2) {
                    mcv_set("default_gene1", optimal_genes[1])
                    mcv_set("default_gene2", optimal_genes[2])
                    gene <- mcv_get(var_name)
                }
            },
            error = function(e) {
                cli::cli_alert_warning("Could not compute optimal default genes: {e$message}")
            }
        )
    }

    gene
}


init_tab_defs <- function() {
    tab_defs <- list(
        "About" = list(
            title = "About",
            module_name = "about",
            icon = "info",
            ui_fn = mod_about_ui,
            sidebar_ui_fn = mod_about_sidebar_ui,
            server_fn = mod_about_server
        ),
        "QC" = list(
            title = "QC",
            module_name = "qc",
            icon = "check",
            ui_fn = mod_qc_ui,
            sidebar_ui_fn = mod_qc_sidebar_ui,
            server_fn = mod_qc_server
        ),
        "Projection QC" = list(
            title = "Projection QC",
            module_name = "projection_qc",
            icon = "clipboard-check",
            ui_fn = mod_projection_qc_ui,
            sidebar_ui_fn = mod_projection_qc_sidebar_ui,
            server_fn = mod_projection_qc_server
        ),
        "Manifold" = list(
            title = "Manifold",
            module_name = "manifold",
            icon = "project-diagram",
            ui_fn = mod_manifold_ui,
            sidebar_ui_fn = mod_manifold_sidebar_ui,
            server_fn = mod_manifold_server
        ),
        "Genes" = list(
            title = "Genes",
            module_name = "gene_mc",
            icon = "wind",
            ui_fn = mod_gene_mc_ui,
            sidebar_ui_fn = mod_gene_mc_sidebar_ui,
            server_fn = mod_gene_mc_server
        ),
        "Markers" = list(
            title = "Markers",
            module_name = "markers",
            icon = "map-marker",
            ui_fn = mod_markers_ui,
            sidebar_ui_fn = mod_markers_sidebar_ui,
            server_fn = mod_markers_server
        ),
        "Diff. Expression" = list(
            title = "Diff. Expression",
            module_name = "mc_mc",
            icon = "chart-bar",
            ui_fn = mod_mc_mc_ui,
            sidebar_ui_fn = mod_mc_mc_sidebar_ui,
            server_fn = mod_mc_mc_server
        ),
        "Gene modules" = list(
            title = "Gene modules",
            module_name = "gene_modules",
            icon = "layer-group",
            ui_fn = mod_gene_modules_ui,
            sidebar_ui_fn = mod_gene_modules_sidebar_ui,
            server_fn = mod_gene_modules_server
        ),
        "Cell types" = list(
            title = "Cell types",
            module_name = "cell_type",
            icon = "bacteria",
            ui_fn = mod_cell_type_ui,
            sidebar_ui_fn = mod_cell_type_sidebar_ui,
            server_fn = mod_cell_type_server
        ),
        "Annotate" = list(
            title = "Annotate",
            module_name = "annotate",
            icon = "pen",
            ui_fn = mod_annotate_ui,
            sidebar_ui_fn = mod_annotate_sidebar_ui,
            server_fn = mod_annotate_server
        ),
        "Projected-fold" = list(
            title = "Projected-fold",
            module_name = "proj_fold",
            icon = "cloud-moon-rain",
            ui_fn = mod_proj_fold_ui,
            sidebar_ui_fn = mod_proj_fold_sidebar_ui,
            server_fn = mod_proj_fold_server
        ),
        "Flow" = list(
            title = "Flow",
            module_name = "flow",
            icon = "water",
            ui_fn = mod_flow_ui,
            sidebar_ui_fn = mod_flow_sidebar_ui,
            server_fn = mod_flow_server
        ),
        "Samples" = list(
            title = "Samples",
            module_name = "samples",
            icon = "vials",
            ui_fn = mod_samples_ui,
            sidebar_ui_fn = mod_samples_sidebar_ui,
            server_fn = mod_samples_server
        ),
        "Query" = list(
            title = "Query",
            module_name = "query",
            icon = "video",
            ui_fn = mod_query_ui,
            sidebar_ui_fn = mod_query_sidebar_ui,
            server_fn = mod_query_server
        ),
        "Atlas" = list(
            title = "Atlas",
            module_name = "atlas",
            icon = "atlas",
            ui_fn = mod_atlas_ui,
            sidebar_ui_fn = mod_atlas_sidebar_ui,
            server_fn = mod_atlas_server
        ),
        "Gene correlation" = list(
            title = "Gene correlation",
            module_name = "gene_correlation",
            icon = "exchange-alt",
            ui_fn = mod_gene_correlation_ui,
            sidebar_ui_fn = NULL,
            server_fn = mod_gene_correlation_server
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
        # exclude Gene correlation tab in light version
        if (is.null(cur_config$excluded_tabs)) {
            cur_config$excluded_tabs <- c("Gene correlation")
        } else {
            cur_config$excluded_tabs <- unique(c(cur_config$excluded_tabs, "Gene correlation"))
        }

        # make the About tab first if exists
        if ("About" %in% cur_config$tabs) {
            cur_config$tabs <- c("About", cur_config$tabs[cur_config$tabs != "About"])
        }
    }

    mcv_set("config", cur_config)
}

order_tabs <- function(tabs) {
    tabs_order <- c("QC", "Projection QC", "Manifold", "Genes", "Query", "Atlas", "Markers", "Gene modules", "Gene correlation", "Projected-fold", "Diff. Expression", "Samples", "Cell types", "Annotate")

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
    # Default cache configuration.
    # `enabled = TRUE`: the derived chain is additive, not destructive — it
    # creates a `.mcview_cache` sibling directory next to the base DAF. Without
    # it, every cold session pays ~6 s to recompute top1/top2 gene vectors in
    # `convert_daf_metacell_types` because the fast-path vectors are never
    # persisted. Users who want the pre-April-2026 behavior can set
    # `mcview_cache_enabled = FALSE` as a DAF scalar or in YAML config.
    defaults <- list(
        enabled = TRUE,
        type = "files",
        cache_dir = ".mcview_cache",
        precompute_on_startup = TRUE,
        invalidation = list(
            strategy = "hash",
            version_scalar = "mcview_data_version"
        )
    )

    cache_config <- defaults

    # Priority 3: DAF scalars (lowest priority for cache-specific settings)
    # Batch: get the set of all available scalars in one call to avoid
    # individual has_scalar roundtrips for each cache config scalar.
    if (!is.null(daf_obj)) {
        avail_scalars <- dafr::scalars_set(daf_obj)
        get_cache_scalar <- function(name) {
            if (name %in% avail_scalars) dafr::get_scalar(daf_obj, name) else NULL
        }

        daf_cache_enabled <- get_cache_scalar("mcview_cache_enabled")
        daf_cache_type <- get_cache_scalar("mcview_cache_type")
        daf_cache_dir <- get_cache_scalar("mcview_cache_dir")
        daf_cache_precompute <- get_cache_scalar("mcview_cache_precompute")

        # Also support legacy scalars
        daf_cache_in_daf <- get_cache_scalar("mcview_cache_in_daf")
        daf_cache_daf_root <- get_cache_scalar("mcview_cache_daf_root")

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
                cache_config$invalidation <- utils::modifyList(cache_config$invalidation, yaml_cache$invalidation)
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
    precompute_on_startup = TRUE) {
    list(
        enabled = enabled,
        type = type,
        cache_dir = cache_dir,
        precompute_on_startup = precompute_on_startup,
        invalidation = list(
            strategy = "hash",
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
