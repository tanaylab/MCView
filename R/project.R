project_config_file <- function(path) {
    fs::path(path, "config", "config.yaml")
}

project_help_file <- function(path) {
    fn <- fs::path(path, "config", "help.yaml")
    if (!file.exists(fn)) {
        fn <- app_sys("default-config", "help.yaml")
    }
    return(fn)
}

project_about_file <- function(path) {
    fn <- fs::path(path, "config", "about.Rmd")
    if (!file.exists(fn)) {
        fn <- app_sys("default-config", "about.Rmd")
    }
    return(fn)
}

project_cache_dir <- function(path) {
    fs::path(path, "cache")
}

verify_project_dir <- function(path, create = FALSE, atlas = FALSE) {
    if (!dir.exists(path)) {
        if (create) {
            create_project(project = basename(path), edit_config = FALSE, atlas = atlas)
        } else {
            cli_abort("{.path path} does not exist. Maybe there is a typo? You can start a new project by running {.code MCView::create_project}.")
        }
    } else {
        config_file <- project_config_file(path)
        if (!file.exists(config_file)) {
            cli_abort("No config file found in {.file {config_file}}. Are you sure this is an MCView project dir? You can start a new project by running {.code MCView::create_project}.")
        }
    }
}

create_project_dirs <- function(project_dir, atlas = FALSE) {
    fs::dir_create(project_dir)
    cli_alert_info("creating {project_dir}")

    fs::dir_create(project_cache_dir(project_dir))
    fs::dir_create(fs::path(project_dir, "config"))

    if (atlas) {
        defaults_dir <- app_sys("atlas-proj-config")
    } else {
        defaults_dir <- app_sys("default-config")
    }

    files <- c("config.yaml", "help.yaml", "about.Rmd")
    for (file in files) {
        if (!fs::file_exists(fs::path(project_dir, "config", file))) {
            fs::file_copy(fs::path(defaults_dir, file), fs::path(project_dir, "config", file))
        }
    }

    return(project_dir)
}


#' Create a configuration folder for a project
#'
#' Create a project directory with the default configuration files and directory structure.
#' A text editor would be opened in order to edit the project \code{config.yaml} file.
#'
#' @param project path of the project
#' @param edit_config open file editor for config file editing
#' @param atlas use default configuration for atlas projections
#'
#' @examples
#' \dontrun{
#' dir.create("raw")
#' download.file("http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz", "raw/#' PBMC_processed.tar.gz")
#' untar("raw/PBMC_processed.tar.gz", exdir = "raw")
#' create_project("PBMC")
#' }
#'
#' @export
create_project <- function(project, edit_config = TRUE, atlas = FALSE) {
    project_dir <- create_project_dirs(project, atlas = atlas)
    project_dir <- fs::path(project_dir, "config")

    if (rlang::is_interactive() && edit_config) {
        utils::file.edit(fs::path(project_dir, "config.yaml"))
    }
    cli_alert("please edit {.file {fs::path(project_dir, 'config.yaml')}}")
}

#' Generate a 'deployment ready' bundle of the a project app
#'
#' Generate a 'deployment ready' bundle of the a project app
#'
#' Create a minimal shiny app in \code{path}/\code{name} directory which would contain:
#' \itemize{
#' \item{}{app.R file. }
#' \item{}{project config and cache. }
#' }
#'
#' The bundle can then be deployed in shiny-server, shinyapps.io or any other environment that supports serving shiny apps.
#'
#' Note: when deploying to these services - make sure you have the MCView package installed.
#'
#' @param project path to the project directory
#' @param path path in which to create the bundle.
#' @param name name of the folder in which to create the bundle. The bundle would be created at \code{path}/\code{name}
#' @param overwrite overwrite bundle if already exists
#'
#' @examples
#' \dontrun{
#' MCView::create_bundle(project = "PBMC", path = getwd(), name = "PBMC")
#' }
#'
#' @export
create_bundle <- function(project, path = getwd(), name = "MCView_bundle", overwrite = FALSE) {
    bundle_dir <- fs::path(path, name)
    if (!(fs::dir_exists(project))) {
        cli::cli_abort("{.path {project}} does not exists.")
    }
    if (fs::dir_exists(bundle_dir)) {
        if (overwrite) {
            fs::dir_delete(bundle_dir)
            fs::dir_create(bundle_dir)
        } else {
            cli::cli_abort("{.path {bundle_dir}} already exists. Run with {.code overwrite=TRUE} to force overwriting it.")
        }
    } else {
        fs::dir_create(bundle_dir)
    }

    fs::file_copy(app_sys("app.R"), fs::path(bundle_dir, "app.R"))
    fs::dir_copy(project, fs::path(bundle_dir, "project"))

    cli::cat_line("Bundle files:")
    fs::dir_tree(bundle_dir)
    cli::cat_line("")
    cli::cli_alert_success("created a bundle at {bundle_dir}")
    cli::cli_bullets("To deploy to shinyapps.io, run: {.code rsconnect::deployApp(appDir = '{bundle_dir}')}\n")
    cli::cli_bullets("To deploy to another shiny-server service, upload {.path {bundle_dir}} to the service.")
}
