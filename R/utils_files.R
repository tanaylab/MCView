command_file_path <- function(path) {
    fs::path(path, "config", "IMPORT_COMMAND.R")
}

#' Save function call to a file
#'
#' This helper function captures the calling function's arguments and
#' saves them to a specified file.
#'
#' @param command_file name or path of the file to save the command
#' @param call_depth integer specifying how many levels up the call stack to capture (default is 1, which captures the direct caller)
#' @param add_details logical specifying whether to add details about the user and time to the saved command
#'
#' @return invisibly returns the path to the saved file
#'
#' @examples
#' \dontrun{
#' # Inside another function:
#' save_function_call()
#' save_function_call("import.R")
#' }
#' @noRd
save_function_call <- function(command_file = "IMPORT_COMMAND", call_depth = 1, add_details = TRUE, project = NULL) {
    # Get the call from the specified depth in the call stack
    parent_call <- sys.call(-call_depth)

    # Convert to a string with a generous width to avoid unnecessary line breaks
    call_string <- deparse(parent_call, width.cutoff = 500)

    # Combine lines if it's a multi-line string
    if (length(call_string) > 1) {
        call_string <- paste(call_string, collapse = "\n")
    }

    call_string <- paste0("setwd(\"", getwd(), "\")\n", call_string)
    call_string <- paste0("library(MCView)\n\n", call_string)

    if (add_details) {
        call_string <- paste0(call_string, "\n\n# Command was run by ", Sys.info()["user"], " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n# MCView version: ", utils::packageVersion("MCView"))
        if (!is.null(project) && fs::file_exists(project_metacells_algorithm_file(project))) {
            mc_version <- readLines(project_metacells_algorithm_file(project))
            call_string <- paste0(call_string, "\n# Metacells algorithm version: ", mc_version)
        }
    }

    if (fs::file_exists(command_file)) {
        fs::file_delete(command_file)
    }

    writeLines(call_string, command_file)
    return(invisible(command_file))
}
