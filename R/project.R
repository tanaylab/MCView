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

project_data_dir <- function(path) {
    fs::path(path, "data")
}

verify_project_dir <- function(path) {
    if (!dir.exists(path)) {
        cli_abort("{path} does not exist. Maybe there is a typo? You can start a new project by running MCView::create_project.")
    }
    config_file <- project_config_file(path)
    if (!file.exists(config_file)) {
        cli_abort("No config file found in {config_file}. Are you sure this is an MCView project dir? You can start a new project by running MCView::create_project.")
    }
}

create_project_dirs <- function(project_dir, project) {
    if (is.null(project_dir)) {
        project_dir <- fs::path(getwd(), project)
        cli_alert_info("creating {project} under current directory")
    }

    fs::dir_create(project_dir)
    cli_alert_info("creating {project_dir}")

    fs::dir_create(project_dir, "data")

    defaults_dir <- app_sys("default-config")
    fs::dir_copy(path = defaults_dir, new_path = fs::path(project_dir, "config"), overwrite = FALSE)

    return(project_dir)
}


#' Create a configuration folder for a project
#'
#' This would create a project directory with the default configuration files and directory structure.
#' A text editor would be opened in order to edit the project \code{config.yaml} file.
#'
#' @param project name of the project
#' @param project_dir directory of the project. Would be created if not existing. If \code{NULL} - the project would be created under current directory.
#'
#' @examples
#'\dontrun{
#' dir.create("raw")
#' download.file("http://www.wisdom.weizmann.ac.il/~atanay/metac_data/PBMC_processed.tar.gz", "raw/#' PBMC_processed.tar.gz")
#' untar("raw/PBMC_processed.tar.gz", exdir = "raw")
#' create_project("PBMC")
#' }
#'
#' @export
create_project <- function(project, project_dir = NULL) {
    project_dir <- create_project_dirs(project_dir, project)
    project_dir <- fs::path(project_dir, "config")

    utils::file.edit(fs::path(project_dir, "config.yaml"))
    cli_alert_info("please edit {fs::path(project_dir, 'config.yaml')}")
}


get_default_project <- function() {
    default_dir <- fs::path(app_sys("data-raw"), "default")
    verify_project_dir(default_dir)
    return(default_dir)
}

#' Set default project
#'
#' This creates a symbolic link to the \code{project} in order to run MCView without
#' specifying a project directory. Note that in order for this to work you must have
#' write permissions in the installation directory of MCView. Used mainly for debugging purposes.
#'
#' @param project path to the project directory
#'
#'
#' @export
set_default_project <- function(project) {
    default_dir <- fs::path(app_sys("data-raw"), "default")
    if (fs::link_exists(default_dir)) {
        fs::link_delete(default_dir)
    }
    fs::link_create(fs::path_abs(project), default_dir, symbolic = TRUE)

    cli_alert_info("default project was set to {project}")
}


#' Generate a 'deployment ready' bundle of the a project app
#'
#' This would create a minimal shiny app in \code{path}/\code{name} directory which would contain:
#' \itemize{
#' \item{}{app.R file. }
#' \item{}{project config and data. }
#' }
#' The bundle can then be deployed in shiny-server, shinyapps.io or any other environment that supports serving shiny apps.
#' Note: when deploying to these services - make sure you have the MCView package installed.
#'
#' @param project path to the project directory
#' @param path path in which to create the bundle.
#'
#' @examples
#'\dontrun{
#' MCView::create_bundle(project = "PBMC", path = getwd(), name = "PBMC")
#' }
#'
#' @export
create_bundle <- function(project, path = getwd(), name = "MCView_bundle") {
    bundle_dir <- fs::path(path, name)
    if (!fs::dir_exists(bundle_dir)) {
        fs::dir_create(bundle_dir)
    }

    fs::file_copy(app_sys("app.R"), fs::path(bundle_dir, "app.R"))
    fs::dir_copy(project, fs::path(bundle_dir, "project"))

    fs::dir_tree(bundle_dir)
    cli::cli_alert_success("created a bundle at {bundle_dir}")
}
