#' create a bundle for current project
#'
#' @return path to the new bundle
#'
#' @noRd
download_project <- function(file, project) {
    dir <- tempdir()
    bundle_dir <- fs::path(dir, project)

    fs::dir_copy(project, bundle_dir)
    on.exit(fs::dir_delete(bundle_dir))
    zip::zip(file,
        bundle_dir,
        include_directories = TRUE,
        recurse = TRUE,
        mode = "cherry-pick",
        compression_level = 1
    )

    invisible(file)
}
