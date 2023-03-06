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

project_version_file <- function(path) {
    fs::path(path, "config", "MCVIEW_VERSION")
}

verify_version <- function(path) {
    version_file <- project_version_file(path)
    if (!file.exists(version_file)) {
        cli_alert_warning("No version file found. This is probably an old project. Please re-import the project with the latest version of MCView.")
    } else {
        version <- readLines(version_file)
        if (utils::compareVersion(version, as.character(packageVersion("MCView"))) < 0) {
            cli_alert_warning("Project was created with MCView version {.field {version}} while the current version is {.field {packageVersion('MCView')}}. Please re-import the project with the latest version of MCView.")
        }
    }
}

verify_project_dir <- function(path, create = FALSE, atlas = FALSE, ...) {
    if (!dir.exists(path)) {
        if (create) {
            create_project(project = path, edit_config = FALSE, atlas = atlas, ...)
        } else {
            cli_abort("{.path path} does not exist. Maybe there is a typo? You can start a new project by running {.code MCView::create_project}.")
        }
    } else {
        config_file <- project_config_file(path)
        if (!file.exists(config_file)) {
            if (file.exists(project_config_file(file.path(path, "project")))) {
                path <- file.path(path, "project")
            } else {
                cli_abort("No config file found at {.file {config_file}}. Are you sure this is an MCView project dir? You can start a new project by running {.code MCView::create_project}.")
            }
        }
    }
    invisible(path)
}

create_project_dirs <- function(project_dir, atlas = FALSE) {
    fs::dir_create(project_dir)
    cli_alert_info("creating {project_dir}")

    fs::dir_create(project_cache_dir(project_dir))
    fs::dir_create(fs::path(project_dir, "config"))

    defaults_dir <- app_sys("default-config")

    files <- c("help.yaml", "about.Rmd")
    for (file in files) {
        if (!fs::file_exists(fs::path(project_dir, "config", file))) {
            fs::file_copy(fs::path(defaults_dir, file), fs::path(project_dir, "config", file))
        }
    }

    return(project_dir)
}

create_project_config_file <- function(project_dir,
                                       title = "MCView",
                                       tabs = NULL,
                                       help = FALSE,
                                       selected_gene1 = NULL,
                                       selected_gene2 = NULL,
                                       selected_mc1 = NULL,
                                       selected_mc2 = NULL,
                                       datasets = NULL,
                                       other_params = NULL,
                                       atlas = FALSE) {
    config <- list()
    config$title <- title
    if (is.null(tabs)) {
        if (atlas) {
            tabs <- c("Manifold", "Genes", "Query", "Atlas", "Markers", "Gene modules", "Projected-fold", "Diff. Expression", "Cell types", "Annotate", "About")
        } else {
            tabs <- c("Manifold", "Genes", "Markers", "Gene modules", "Diff. Expression", "Cell types", "Annotate", "About")
        }
    }
    config$tabs <- tabs
    config$help <- help
    config$selected_gene1 <- selected_gene1
    config$selected_gene2 <- selected_gene2
    config$selected_mc1 <- selected_mc1
    config$selected_mc2 <- selected_mc2
    if (!is.null(other_params)) {
        config <- c(config, other_params)
    }
    config$datasets <- datasets
    yaml::write_yaml(config, fs::path(project_dir, "config", "config.yaml"))
}


#' Create a configuration folder for a project
#'
#' Create a project directory with the default configuration files and directory structure.
#' If \code{edit_config == TRUE}, a text editor would be opened in order to edit the project \code{config.yaml} file.
#'
#' @param project path of the project
#' @param title The title of the app. This would be shown on the top left of the screen.
#' @param tabs Controls which tabs to show in the left sidebar and their order. Options are: "Manifold", "Genes", "Query", "Atlas", "Markers", "Gene modules", "Projected-fold", "Diff. Expression", "Cell types", "Flow", "Annotate", "About". When NULL - default tabs would be set. For projects with atlas projections, please set \code{atlas} to TRUE.
#' @param help Controls wether to start the app with a help modal (from introjs). Help messages can be edited in help.yaml file (see 'Architecture' vignette).
#' @param selected_gene1,selected_gene2 The default genes that would be selected (in any screen with gene selection). If this parameter is missing, the 2 genes with highest max(expr)-min(expr) in the first dataset would be chosen.
#' @param selected_mc1,selected_mc2 The default metacells that would be selected in the Diff. Expression tab.
#' @param datasets A named list with additional per-dataset parameters. Current parameters include default visualization properties of projection and scatter plots.
#' @param other_params Named list of additional parameters such as projection_point_size, projection_point_stroke, scatters_point_size and scatters_stroke_size
#' @param edit_config open file editor for config file editing
#' @param atlas use default configuration for atlas projections (relevant only when \code{tabs} is NULL)
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
create_project <- function(project,
                           title = "MCView",
                           tabs = NULL,
                           help = FALSE,
                           selected_gene1 = NULL,
                           selected_gene2 = NULL,
                           selected_mc1 = NULL,
                           selected_mc2 = NULL,
                           datasets = NULL,
                           other_params = NULL,
                           edit_config = TRUE,
                           atlas = FALSE) {
    project_dir <- create_project_dirs(project, atlas = atlas)
    config <- create_project_config_file(
        project = project,
        title = title,
        tabs = tabs,
        help = help,
        selected_gene1 = selected_gene1,
        selected_gene2 = selected_gene2,
        selected_mc1 = selected_mc1,
        selected_mc2 = selected_mc2,
        datasets = datasets,
        other_params = other_params,
        atlas = atlas
    )
    project_dir <- fs::path(project_dir, "config")

    if (rlang::is_interactive() && edit_config) {
        utils::file.edit(fs::path(project_dir, "config.yaml"))
    }
    cli_alert("You can edit the app configuration at {.file {fs::path(project_dir, 'config.yaml')}}")
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
#' @param self_contained include the source code of \code{MCView} in the bundle
#' and use it to run the app. Use this in order to ensure that the package would always
#' run the same way, regardless of MCView changes. When this option is FALSE,
#' the version of \code{MCView} which is installed on the server would be loaded, which can be occasionally
#' be different than the one used when creating the app. By default, the code uses the latest \code{MCView} release would be used, see \code{branch} for
#' other options.
#' @param branch name of the \code{MCView} branch to include when \code{self_contained=TRUE}. By default, the latest release would be used. You can set this
#' parameter to NULL in order to include the current development version
#' ('master' branch), or set it to any other branch in the 'tanaylab/MCView' github
#' repository.
#' @param restart add a file named 'restart.txt' to the bundle. This would force shiny-server to restart the app when updated.
#' @param permissions change the file permissions of the bundle after creation, e.g. "777". When NULL -
#' permissions would not be changed.
#' @param light_version create a light version of the bundle, which would not include features that require heavy computation (e.g. changing Marker genes, Gene modules etc.)
#'
#' @inheritDotParams gert::git_clone
#' @examples
#' \dontrun{
#' MCView::create_bundle(project = "PBMC", path = getwd(), name = "PBMC")
#'
#' # latest release
#' MCView::create_bundle(project = "PBMC", path = getwd(), name = "PBMC", self_contained = TRUE)
#'
#' # development version
#' MCView::create_bundle(project = "PBMC", path = getwd(), name = "PBMC", self_contained = TRUE, branch = NULL)
#'
#' # specific branch
#' MCView::create_bundle(project = "PBMC", path = getwd(), name = "PBMC", self_contained = TRUE, branch = "feat@atlas-projection")
#' }
#'
#' @export
create_bundle <- function(project, path = getwd(), name = "MCView_bundle", overwrite = FALSE, self_contained = FALSE, branch = "latest_release", restart = overwrite, permissions = NULL, light_version = FALSE, ...) {
    bundle_dir <- fs::path(path, name)
    if (!(fs::dir_exists(project))) {
        cli::cli_abort("{.path {project}} does not exists.")
    }
    if (fs::dir_exists(bundle_dir)) {
        if (overwrite) {
            fs::dir_delete(bundle_dir)
            fs::dir_create(bundle_dir)
            cli::cli_li("Removing previous bundle ({.field overwrite = TRUE})")
        } else {
            cli::cli_abort("{.path {bundle_dir}} already exists. Run with {.code overwrite=TRUE} to force overwriting it.")
        }
    } else {
        fs::dir_create(bundle_dir)
    }

    if (self_contained) {
        cli::cli_alert("Creating a self-contained bundle")
        code_dir <- fs::path(bundle_dir, "code")
        if (!is.null(branch) && branch == "latest_release") {
            gert::git_clone("git@github.com:tanaylab/MCView", path = code_dir, ...)
            tag_list <- gert::git_tag_list(repo = code_dir)
            latest_tag <- tail(tag_list, n = 1)
            gert::git_branch_create(
                branch = latest_tag$name,
                ref = latest_tag$commit,
                repo = code_dir,
                checkout = TRUE
            )
            cli::cli_alert_info("Using latest release: {.file {latest_tag$name}}")
        } else {
            gert::git_clone("git@github.com:tanaylab/MCView", path = code_dir, branch = branch, ...)
        }
    }

    fs::file_copy(app_sys("app.R"), fs::path(bundle_dir, "app.R"))
    fs::dir_copy(project, fs::path(bundle_dir, "project"))

    if (light_version) {
        bundle_config_file <- project_config_file(fs::path(bundle_dir, "project"))
        bundle_config <- yaml::read_yaml(bundle_config_file)
        bundle_config$light_version <- TRUE
        bundle_config$excluded_tabs <- c("Gene modules", "Annotate", "Inner-fold")
        yaml::write_yaml(bundle_config, bundle_config_file)
        cli::cli_alert("Creating a light version of the bundle. Excluded tabs: {.field Gene modules, Annotate, Inner-fold}. To change this, edit the {.file project/config.yaml} file.")
    }

    if (restart) {
        fs::file_touch(fs::path(bundle_dir, "restart.txt"))
        cli::cli_li("Adding a file called {.field restart.txt}")
    }

    if (!is.null(permissions)) {
        if (!is.character(permissions)) {
            cli::cli_abort("{.field permissions} must be a character string.")
        }
        fs::file_chmod(c(bundle_dir, fs::dir_ls(bundle_dir, recurse = TRUE)), mode = permissions)
        cli::cli_li("Changing permissions to {.field {permissions}}")
    }

    cli::cli_li("Bundle files:")
    fs::dir_tree(bundle_dir)
    cli::cat_line("")
    cli::cli_alert_success("created a bundle at {bundle_dir}")
    cli::cli_li("To deploy to shinyapps.io, run: {.field rsconnect::deployApp(appDir = \"{as.character(bundle_dir)}\")}")
    cli::cli_li("To deploy to another shiny-server service, upload {.path {bundle_dir}} to the service.")
}
