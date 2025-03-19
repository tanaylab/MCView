verify_config_file <- function(config) {
    required_fields <- c(
        "title",
        "tabs",
        "help"
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
    help_file <- project_help_file(project)
    project <<- project
    about_file <<- fs::path_abs(project_about_file(project))
    cache_dir <<- project_cache_dir(project)

    config <<- yaml::read_yaml(config_file)
    verify_config_file(config)
    if (is.null(config$light_version)) {
        config$light_version <<- FALSE
    }

    verify_app_cache(project)

    help_config <<- yaml::read_yaml(help_file)

    if (file.exists(project_metacells_algorithm_file(project))) {
        metacells_version <- readLines(project_metacells_algorithm_file(project))
        config$metacells_version <<- metacells_version
    }
    config$profile <<- profile
}

init_defs <- function() {
    egc_epsilon <<- 1e-5
    options(spinner.type = 6)

    theme_set(theme_classic())

    init_tab_defs()

    init_selected_genes()

    expr_breaks <<- c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1)
}

init_selected_genes <- function() {
    # if selected genes are not set - choose them from the first dataset
    if (is.null(config$selected_gene1) || is.null(config$selected_gene2)) {
        mc_egc <- get_mc_egc(dataset_ls(project)[1])
        if (has_atlas(dataset_ls(project)[1])) {
            mc_egc_atlas <- get_mc_egc(dataset_ls(project)[1], atlas = TRUE)
            mc_egc <- mc_egc[intersect(rownames(mc_egc), rownames(mc_egc_atlas)), ]
        }

        f <- rep(TRUE, nrow(mc_egc))
        lateral_genes <- get_mc_data(dataset_ls(project)[1], "lateral_genes")
        if (!is.null(lateral_genes)) {
            f <- f & !(rownames(mc_egc) %in% lateral_genes)
        }
        noisy_genes <- get_mc_data(dataset_ls(project)[1], "noisy_genes")
        if (!is.null(noisy_genes)) {
            f <- f & !(rownames(mc_egc) %in% noisy_genes)
        }
        mc_egc <- mc_egc[f, ]

        minmax <- matrixStats::rowMaxs(mc_egc, na.rm = TRUE) - matrixStats::rowMins(mc_egc, na.rm = TRUE)
        names(minmax) <- rownames(mc_egc)
        genes <- names(head(sort(minmax, decreasing = TRUE), n = 2))
    }

    default_gene1 <<- config$selected_gene1 %||% genes[1]
    default_gene2 <<- config$selected_gene2 %||% genes[2]
}


init_tab_defs <- function() {
    tab_defs <<- list(
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

    default_tabs <- c("About", "Genes", "Diff. Expression")

    cur_config <- config

    if (!is.null(cur_config$tabs)) {
        cur_config$original_tabs <- cur_config$tabs
        cur_config$tabs[cur_config$tabs == "Metacells"] <- "Diff. Expression" # here for backward compatibility
        cur_config$tabs <- cur_config$tabs[config$tabs != "Metadata"] # ignore "Metadata" for backward compatibility
        for (.x in cur_config$tabs) {
            if (!(.x %in% names(tab_defs))) {
                cli_warn("{.x} is not a valid tab name. Update `tabs` in your configuration file.")
                cur_config$tabs <- cur_config$tabs[cur_config$tabs != .x]
            }
        }
    } else {
        cur_config$tabs <- default_tabs
    }

    if (!is.null(cur_config$excluded_tabs)) {
        tab_defs <<- tab_defs[!(names(tab_defs) %in% cur_config$excluded_tabs)]
        cur_config$tabs <- cur_config$tabs[!(cur_config$tabs %in% cur_config$excluded_tabs)]
    }

    if (!rmarkdown::pandoc_available() && "About" %in% cur_config$tabs) {
        warning("pandoc is not available, removing 'About' tab'")
        cur_config$tabs <- cur_config$tabs[cur_config$tabs != "About"]
    }

    if (!is.null(config$light_version) && config$light_version) {
        # make the About tab first if exists
        if ("About" %in% cur_config$tabs) {
            cur_config$tabs <- c("About", cur_config$tabs[cur_config$tabs != "About"])
        }
    }

    config <<- cur_config
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

gene_names <- function(dataset, atlas = FALSE) {
    gene_names <- rownames(get_mc_data(dataset, "mc_mat", atlas = atlas))
    return(gene_names)
}

gene_names_label <- function(dataset, atlas = FALSE, gene_modules = NULL) {
    gene_names <- gene_names(dataset, atlas)
    if (is.null(gene_names)) {
        return(gene_names)
    }
    names(gene_names) <- gene_label(gene_names, dataset, gene_modules)
    return(gene_names)
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
