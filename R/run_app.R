#' Run the MCView Application
#'
#' Load project cache and run the MCView application.
#'
#'
#' @param project Either a project path (existing behavior) or a DAF object
#' @param dataset_name Name for the dataset when using DAF object (default: "daf_data")
#' @param config_file Optional YAML configuration file (for DAF mode)
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
                    dataset_name = "daf_data",
                    config_file = NULL,
                    port = NULL,
                    host = NULL,
                    launch.browser = FALSE,
                    profile = FALSE,
                    ...) {
    # Initialize environment
    init_mcview_env()

    if (is.character(project)) {
        set_backend("project", path = project)
        mcv_set("project", init_project_dir(project))
        verify_version(mcv_get("project"))
        init_config(project = mcv_get("project"), profile = profile)
        load_all_data(cache_dir = project_cache_dir(mcv_get("project")))
    } else if (inherits(project, "Daf")) {
        stop("DAF mode not yet implemented")
    } else {
        cli_abort("project must be either a project path or a DAF object")
    }

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
    config <- mcv_get("config")
    project <- mcv_get("project")
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
