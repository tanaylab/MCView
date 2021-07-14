#' Run the MCView Application
#'
#' Loads project data and runs the MCView application.
#'
#'
#' @param project path of the project to run.
#' @param port app port
#' @param host app host
#' @param launch.browser launch web browser after app start
#' @param ... A series of options to be used inside the app.
#'
#' @examples 
#'\dontrun{
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
                    ...) {
    init_config(project = project)
    load_all_data(cache_dir = project_cache_dir(project))
    init_defs()

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
