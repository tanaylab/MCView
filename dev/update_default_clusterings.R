# Store the clusterings of the default heatmaps (update_default_clusterings()) in MCView projects that were
# imported before import_dataset() stored them, then check in a new R process that the first render of each
# heatmap tab uses them.
#
#   Rscript dev/update_default_clusterings.R [--overwrite] <project or bundle dir> ...
#
# Runs the MCView of this repository, which should be the version the app runs: the app finds a stored
# clustering by a key that its code computes. Writes only <project>/cache/<dataset>/default_clusterings.qs2,
# and skips datasets that already have it unless --overwrite is given.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
suppressMessages(pkgload::load_all(dirname(dirname(script)), export_all = TRUE, helpers = FALSE, quiet = TRUE))
library(shiny)

# Open each heatmap tab once, as in a new app process, and report whether its first render clustered the
# metacells (calls of fastcluster::hclust) or used the stored clustering.
check <- function(path) {
    project <- init_project_dir(path)
    init_config(project)
    load_all_data(cache_dir = project_cache_dir(project))
    init_defs()

    counts <- new.env()
    suppressMessages({
        trace("hclust", where = asNamespace("fastcluster"), print = FALSE, tracer = bquote(assign("hclust", .(counts)$hclust + 1, envir = .(counts))))
        trace("cluster_metacells", where = asNamespace("MCView"), print = FALSE, tracer = bquote(assign("clusterings", .(counts)$clusterings + 1, envir = .(counts))))
    })

    # the data each tab needs (the Markers tab saves marker_genes when they are missing, so we don't open it then)
    tabs <- c(markers = "marker_genes", gene_modules = "gene_modules", inner_fold = "inner_fold_mat", stdev_fold = "inner_stdev_mat", proj_fold = "projected_fold", outliers = "deviant_fold_mat")
    res <- list()
    for (dataset in dataset_ls(project)) {
        for (tab in names(tabs)[purrr::map_lgl(tabs, ~ has_mc_data(dataset, .x))]) {
            shinyOptions(cache = cachem::cache_mem(max_size = 200e6)) # the app cache of a new process
            counts$hclust <- 0
            counts$clusterings <- 0
            heatmap <- paste0(tab, "-", tab, "_heatmap-")
            # the defaults of the heatmap inputs
            inputs <- list(
                lfp_range = c(-3, 3), midpoint = 0, low_color = "blue", high_color = "red", mid_color = "white",
                plot_legend = TRUE, plot_cell_type_legend = TRUE, plot_genes_legend = TRUE, highlight_color = "#fba236",
                filter_by_clipboard = FALSE, mat_value = "Global", metadata_order_var = "Hierarchical-Clustering",
                metadata_order_cell_type_var = "Hierarchical-Clustering", force_cell_type = TRUE, max_gene_num = 100,
                brush_action = "Zoom", enable_categorical_filter = FALSE, legend_width = 2
            )
            names(inputs) <- paste0(heatmap, names(inputs))
            # the browser sends the selected dataset before the server starts
            server <- function(input, output, session) {
                session$setInputs(dataset = dataset)
                app_server(input, output, session)
            }
            status <- NULL
            tryCatch(
                suppressWarnings(testServer(server, {
                    session$setInputs(tab_sidebar = tab)
                    if (tab == "gene_modules") { # the modules the tab's selector starts with
                        html <- as.character(output[["gene_modules-shown_gene_modules_ui"]]$html)
                        selected <- regmatches(html, gregexpr("<option value=\"[^\"]*\" selected>", html))[[1]]
                        inputs[["gene_modules-shown_gene_modules"]] <- gsub("^<option value=\"|\" selected>$", "", selected)
                    }
                    seconds <- system.time({
                        do.call(session$setInputs, inputs)
                        output[[paste0(heatmap, "heatmap")]]
                    })[["elapsed"]]
                    status <<- sprintf("%.1fs", seconds)
                })),
                error = function(e) status <<- paste("error:", conditionMessage(e))
            )
            hit <- counts$clusterings > 0 && counts$hclust == 0 && !startsWith(status, "error")
            tab_name <- names(tab_defs)[purrr::map_chr(tab_defs, "module_name") == tab]
            res[[length(res) + 1]] <- data.frame(dataset = dataset, tab = tab_name, first_render = status, clusterings = counts$clusterings, hclust_calls = counts$hclust, stored_clustering = if (hit) "used" else "NOT used")
        }
    }
    res <- do.call(rbind, res)
    cat("\n", project, "\n", sep = "")
    options(width = 200)
    print(res, row.names = FALSE, right = FALSE)
    all(res$stored_clustering == "used")
}

args <- commandArgs(TRUE)
if (identical(args[1], "--check")) {
    quit(save = "no", status = if (check(args[2])) 0 else 1)
}
ok <- TRUE
for (path in setdiff(args, "--overwrite")) {
    updated <- tryCatch(
        {
            update_default_clusterings(path, overwrite = "--overwrite" %in% args)
            TRUE
        },
        error = function(e) {
            cli::cli_alert_danger("{path}: {conditionMessage(e)}")
            FALSE
        }
    )
    ok <- updated && system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), "--check", shQuote(path))) == 0 && ok
}
quit(save = "no", status = if (ok) 0 else 1)
