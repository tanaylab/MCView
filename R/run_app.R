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
#' @details The app caches the plots and tables it computes (the shiny cache), by default in memory, so a restarted
#' app computes them again. When the project's config sets \code{shiny_cache_dir} to a path, or when it doesn't
#' set it and the environment variable \code{MCVIEW_SHINY_CACHE_DIR} is set (e.g. for all the apps of a shiny
#' server), the cache is also kept on disk, in a directory under that path that the next processes of the app
#' reuse until the project changes. Directories under that path that weren't used for 30 days are removed.
#' Apps that set \code{MCVIEW_CLUSTERING_CACHE_DIR} share the clusterings of their heatmaps in that directory
#' (see \code{\link{update_default_clusterings}}).
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

    root <- config$shiny_cache_dir
    if (is.null(root) && nzchar(Sys.getenv("MCVIEW_SHINY_CACHE_DIR"))) {
        root <- Sys.getenv("MCVIEW_SHINY_CACHE_DIR")
    }

    cache <- NULL
    if (isTRUE(root)) {
        cache <- cachem::cache_disk(tempfile(paste0("shiny_cache_", basename(project)), tmpdir = tempdir()), max_size = max_size)
    } else if (is.character(root)) {
        # A directory under root that the next processes of the app reuse, until the project changes: its name
        # changes with the project's path, the MCView version, and the names, sizes and times of the project's
        # cache and config files. Outputs are looked up in memory, then on disk.
        cache <- tryCatch(
            {
                # Remove the directories that weren't used for 30 days (e.g. of projects that changed since). Only
                # by age: a process started before its project changed may still use its directory.
                if (dir.exists(root)) {
                    dirs <- list.dirs(root, recursive = FALSE)
                    dirs <- dirs[startsWith(basename(dirs), "shiny_cache_")]
                    unlink(dirs[difftime(Sys.time(), file.info(dirs)$mtime, units = "days") > 30], recursive = TRUE)
                }

                files <- list.files(c(project_cache_dir(project), dirname(project_config_file(project))), recursive = TRUE, full.names = TRUE)
                info <- file.info(files)
                key <- rlang::hash(list(as.character(fs::path_abs(project)), as.character(utils::packageVersion("MCView")), files, info$size, as.numeric(info$mtime)))
                dir <- fs::path(root, paste0("shiny_cache_", basename(project), "_", key))
                if (!dir.exists(dir)) {
                    dir.create(dir, recursive = TRUE)
                }
                # writable by the group, like the root of a server's apps (e.g. setgid to their group), so that
                # any of the apps can remove it when it's old; and marked as used now
                Sys.chmod(dir, file.mode(dir) | as.octmode("070"), use_umask = FALSE)
                Sys.setFileTime(dir, Sys.time())

                # qs2 instead of saveRDS()'s gzip, which took ~14 s to write a 45 MB cached heatmap (qs2 built
                # without TBB uses one thread and warns)
                no_tbb <- function(w) if (grepl("TBB not available", conditionMessage(w))) invokeRestart("muffleWarning")
                disk <- cachem::cache_disk(dir,
                    max_size = max_size,
                    read_fn = function(file) withCallingHandlers(qs2::qs_read(file, nthreads = 4), warning = no_tbb),
                    write_fn = function(value, file) withCallingHandlers(qs2::qs_save(value, file, nthreads = 4), warning = no_tbb),
                    extension = ".qs2"
                )
                mem <- cachem::cache_mem(max_size = max_size)
                # Errors are kept in memory only: shiny caches an output's error with its stack trace, whose
                # calls can hold the data (one took 3 s to write and 6 s to read). A failure on disk (e.g. the
                # directory was removed) is a miss.
                list(
                    get = function(key, missing = cachem::key_missing()) {
                        value <- mem$get(key)
                        if (cachem::is.key_missing(value)) {
                            value <- tryCatch(disk$get(key), error = function(e) cachem::key_missing())
                            if (cachem::is.key_missing(value)) {
                                return(missing)
                            }
                            mem$set(key, value)
                        }
                        value
                    },
                    set = function(key, value) {
                        mem$set(key, value)
                        if (!(is.list(value) && isTRUE(value$error))) {
                            try(disk$set(key, value), silent = TRUE)
                        }
                        invisible(TRUE)
                    }
                )
            },
            error = function(e) {
                cli_alert_warning("Can't keep the shiny cache in {.file {root}}, keeping it in memory: {conditionMessage(e)}")
                NULL
            }
        )
    }

    shinyOptions(cache = cache %||% cachem::cache_mem(max_size = max_size))
}
