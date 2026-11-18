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
    project <<- project
    about_file <<- fs::path_abs(project_about_file(project))
    cache_dir <<- project_cache_dir(project)

    config <<- yaml::read_yaml(config_file)
    verify_config_file(config)
    if (is.null(config$light_version)) {
        config$light_version <<- FALSE
    }

    verify_app_cache(project)

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
        dataset <- dataset_ls(project)[1]
        mc_mat <- get_mc_data(dataset, "mc_mat")

        # filter rows with `rows =` below instead of copying the (genes x metacells) matrix
        f <- rep(TRUE, nrow(mc_mat))
        if (has_atlas(dataset)) {
            f <- rownames(mc_mat) %in% gene_names(dataset, atlas = TRUE)
        }
        lateral_genes <- get_mc_data(dataset, "lateral_genes")
        if (!is.null(lateral_genes)) {
            f <- f & !(rownames(mc_mat) %in% lateral_genes)
        }
        noisy_genes <- get_mc_data(dataset, "noisy_genes")
        if (!is.null(noisy_genes)) {
            f <- f & !(rownames(mc_mat) %in% noisy_genes)
        }

        # choose the 2 genes with the largest range (max - min) of expression over the metacells.
        # Computing the expression of every gene in every metacell takes 20 s on 30K metacells.
        # Instead, bound each gene's range using its range of UMIs and the smallest and largest
        # metacell sizes, and compute the expression only of the genes whose bound reaches the
        # second largest range among the 20 genes with the highest bounds (10-1000 genes).
        rows <- which(f)
        mc_sum <- as.numeric(get_mc_data(dataset, "mc_sum"))
        # both steps below run in 8 forked processes, which also free their memory when they exit
        cores <- if (.Platform$OS.type == "windows") 1 else 8
        # mclapply gives NULL (with a warning) for a process that was killed, e.g. out of memory, and a
        # try-error for one that failed: compute those in this process. n: the length of each result.
        mclapply_or_here <- function(X, FUN, n) {
            res <- parallel::mclapply(X, FUN, mc.cores = cores)
            failed <- !purrr::map2_lgl(res, n, ~ is.numeric(.x) && NROW(.x) == .y)
            res[failed] <- lapply(X[failed], FUN)
            res
        }

        # range of UMIs of each gene, over 8 blocks of metacells (2.7 s -> 0.45 s on 30K metacells)
        blocks <- split(seq_len(ncol(mc_mat)), cut(seq_len(ncol(mc_mat)), 8, labels = FALSE))
        umis <- mclapply_or_here(blocks, function(cols) {
            matrixStats::rowRanges(mc_mat, rows = rows, cols = cols, na.rm = TRUE)
        }, length(rows))
        umis <- cbind(do.call(pmin, lapply(umis, function(x) x[, 1])), do.call(pmax, lapply(umis, function(x) x[, 2])))
        sizes <- range(mc_sum)
        bound <- pmax(umis[, 2] / sizes[1], umis[, 2] / sizes[2]) - pmin(umis[, 1] / sizes[1], umis[, 1] / sizes[2])

        expr_range <- function(i) {
            # 100 genes at a time: the expression of all the genes takes as much memory as mc_mat
            chunks <- split(rows[i], ceiling(seq_along(i) / 100))
            unlist(mclapply_or_here(chunks, function(genes) {
                mc_egc <- mc_mat[genes, , drop = FALSE] / rep(mc_sum, each = length(genes))
                matrixStats::rowMaxs(mc_egc, na.rm = TRUE) - matrixStats::rowMins(mc_egc, na.rm = TRUE)
            }, lengths(chunks)), use.names = FALSE)
        }
        top <- head(order(bound, decreasing = TRUE), n = 20)
        candidates <- which(bound >= min(head(sort(expr_range(top), decreasing = TRUE), n = 2)))

        minmax <- expr_range(candidates)
        names(minmax) <- rownames(mc_mat)[rows[candidates]]
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
        ),
        "Gene correlation" = list(
            title = "Gene correlation",
            module_name = "gene_correlation",
            icon = "exchange-alt"
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

    config <<- cur_config
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
