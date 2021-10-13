#' create a bundle for current project
#'
#' @return path to the new bundle
#'
#' @noRd
download_project <- function(project, dir) {
    create_bundle(
        project,
        path = dir,
        name = project,
        overwrite = TRUE
    )
    return(fs::path(dir, project))
}
