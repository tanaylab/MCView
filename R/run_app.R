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

    with_golem_options(
        app = shinyApp(
            ui = app_ui,
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
