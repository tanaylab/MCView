#' MCView Global State Management
#'
#' Centralized environment for all MCView global state to improve
#' isolation, testing, and maintainability.

# Create the global environment
mcview_env <- new.env(parent = emptyenv())

#' Initialize MCView environment with default values
init_mcview_env <- function() {
    mcview_env$backend <- NULL
    mcview_env$config <- NULL
    mcview_env$mc_data <- NULL
    mcview_env$cache_dir <- NULL
    mcview_env$about_file <- NULL
    mcview_env$project <- NULL
    mcview_env$dataset <- NULL
    mcview_env$tab_defs <- NULL

    invisible(TRUE)
}

#' Get value from MCView environment
#' @param var_name Name of the variable
mcv_get <- function(var_name) {
    get(var_name, envir = mcview_env, inherits = FALSE)
}

#' Set value in MCView environment
#' @param var_name Name of the variable
#' @param value Value to set
mcv_set <- function(var_name, value) {
    assign(var_name, value, envir = mcview_env)
    invisible(value)
}

#' Check if variable exists in MCView environment
#' @param var_name Name of the variable
mcv_exists <- function(var_name) {
    exists(var_name, envir = mcview_env, inherits = FALSE)
}

#' Get current backend information
current_backend <- function() {
    backend <- mcv_get("backend")
    if (is.null(backend)) {
        # Default to project mode for backward compatibility
        return(list(kind = "project", path = mcv_get("project")))
    }
    backend
}

#' Set backend mode
#' @param kind Backend type: "project" or "daf"
#' @param path Project path (for project mode)
#' @param daf_obj DAF object (for DAF mode)
set_backend <- function(kind, path = NULL, daf_obj = NULL) {
    backend <- list(kind = kind)

    if (kind == "project") {
        if (is.null(path)) stop("path required for project backend")
        backend$path <- path
    } else if (kind == "daf") {
        if (is.null(daf_obj)) stop("daf_obj required for daf backend")
        backend$daf_obj <- daf_obj
    } else {
        stop("kind must be 'project' or 'daf'")
    }

    mcv_set("backend", backend)
    invisible(backend)
}

# Initialize on package load
init_mcview_env()
