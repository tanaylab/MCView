#' Run the MCView Application
#'
#' Load project cache and run the MCView application.
#'
#'
#' @param project path of the project to run.
#' @param port app port
#' @param host app host
#' @param launch.browser launch web browser after app start
#' @param profile enable profiling for debugging
#' @param ... A series of options to be used inside the app.
#'
#' @examples
#' \dontrun{
#' MCView::run_app("PBMC")
#' MCView::run_app(project = "PBMC", port = 5555, host = "127.0.0.1")
#' }
#'
#' @inheritDotParams shiny::shinyApp
#'
#' @export
run_app <- function(project,
                    port = NULL,
                    host = NULL,
                    launch.browser = FALSE,
                    profile = FALSE,
                    ...) {
    project <- init_project_dir(project)
    verify_version(project)
    init_config(project = project, profile = profile)
    load_all_data(cache_dir = project_cache_dir(project))
    init_defs()

    config_shiny_cache()

    future::plan(future::multicore)

    # tgstat guesses whether R uses an optimized BLAS from the library's file name, which misses
    # Debian/Ubuntu's OpenBLAS (.../openblas-pthread/libblas.so.3), and then takes a path that is
    # ~10x slower. Check the full path instead, unless the user set the option. The option is
    # global: it applies to every tgstat call in this R process, not only MCView's.
    # Drop this when MCView requires tgstat >= 2.5.0 (tanaylab/tgstat#26 fixed its BLAS detection).
    if (is.null(getOption("tgs_use.blas")) && grepl("openblas|mkl|atlas", extSoftVersion()[["BLAS"]], ignore.case = TRUE)) {
        options(tgs_use.blas = TRUE)
    }

    # The UI depends only on the project, so build it once per process instead of on every
    # page load. Building it here also loads shinyWidgets and shinyjqui, which register
    # their input handlers on load: a browser that reconnects to a restarted app sends
    # those inputs before requesting any page, and without the handlers the app crashed.
    ui <- app_ui()

    # shiny also renders the UI to HTML on every page request, which takes ~1 s for this UI.
    # The HTML is the same for every request, so keep the first one. shiny's renderPage() is
    # not exported: if it fails, shiny renders the page as usual, and it isn't tried again.
    page <- NULL # FALSE once renderPage() failed
    render_page <- function() {
        tryCatch(
            httpResponse(200, content = utils::getFromNamespace("renderPage", "shiny")(ui, testMode = isTRUE(getShinyOption("testmode")))),
            error = function(e) FALSE
        )
    }

    # Render it now, in a forked process, while this one reads the expression matrix: the first
    # session reads it anyway (to check the genes of the gene modules), and both take ~1.5 s for a
    # large project. Rendering registers the resource paths of the UI's dependencies (its JS and
    # CSS files): register them here too. Without fork (Windows), or when the forked process fails
    # or doesn't finish within a minute, the page is rendered on the first request, as before.
    # The forked process must not read cache files (qs2 can't start its threads after a fork, and
    # this process wouldn't get the data), so it runs without them: UI code that read data would
    # fail there, and the page would be rendered here instead.
    if (.Platform$OS.type == "unix") {
        job <- tryCatch(
            parallel::mcparallel(
                {
                    rm(list = "mc_data", envir = globalenv())
                    list(page = render_page(), paths = resourcePaths())
                },
                silent = TRUE
            ),
            error = function(e) NULL
        )
        if (has_mc_data(dataset_ls(project)[1], "gene_modules")) {
            invisible(gene_names(dataset_ls(project)[1]))
        }
        if (!is.null(job)) {
            res <- parallel::mccollect(job, wait = FALSE, timeout = 60)[[1]]
            if (is.null(res)) {
                tools::pskill(job$pid)
                suppressWarnings(parallel::mccollect(job, wait = FALSE, timeout = 5))
            }
            if (is.list(res) && inherits(res$page, "httpResponse")) {
                for (prefix in setdiff(names(res$paths), names(resourcePaths()))) {
                    addResourcePath(prefix, res$paths[[prefix]])
                }
                page <- res$page
            }
        }
    }

    with_golem_options(
        app = shinyApp(
            ui = function(request) {
                page <<- page %||% render_page()
                if (isFALSE(page)) ui else page
            },
            server = app_server,
            options = list(port = port, host = host, launch.browser = launch.browser)
        ),
        golem_opts = list(...)
    )
}

config_shiny_cache <- function() {
    if (!is.null(config$shiny_cache_max_size)) {
        max_size <- config$shiny_cache_max_size
    } else {
        max_size <- 200e6
    }

    if (!is.null(config$shiny_cache_dir)) {
        if (is.logical(config$shiny_cache_dir) && config$shiny_cache_dir) {
            shiny_cache_dir <- tempfile(paste0("shiny_cache_", basename(project)), tmpdir = tempdir())
        } else {
            shiny_cache_dir <- tempfile(paste0("shiny_cache_", basename(project)), tmpdir = config$shiny_cache_dir)
        }

        shinyOptions(cache = cachem::cache_disk(shiny_cache_dir, max_size = max_size))
    } else {
        shinyOptions(cache = cachem::cache_mem(max_size = max_size))
    }
}
