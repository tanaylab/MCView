#' Unified notification helper for MCView
#'
#' Routes a single notification to the console (via `cli::cli_alert_*`)
#' and/or the Shiny UI (via `shiny::showNotification`). One call site,
#' two channels - so a user running from the R console sees the same
#' message as a user in the browser.
#'
#' Severity mapping:
#' \itemize{
#'   \item `"info"`    -> `cli_alert_info`    + `showNotification(type = "default")`
#'   \item `"success"` -> `cli_alert_success` + `showNotification(type = "message")`
#'   \item `"warning"` -> `cli_alert_warning` + `showNotification(type = "warning")`
#'   \item `"error"`   -> `cli_alert_danger`  + `showNotification(type = "error")`
#' }
#'
#' UI emission is silently skipped when no Shiny reactive domain is
#' active (e.g. when called from package load, scripts, or tests). Use
#' [cli::cli_abort()] directly for fatal errors that must halt execution.
#'
#' @param severity One of `"info"`, `"success"`, `"warning"`, `"error"`.
#' @param message  Message string. Passed through `cli` formatters - so
#'   inline markup like `{.val foo}` and `{.code bar}` is honored.
#' @param where    `"both"` (default), `"ui"`, or `"console"`.
#' @param duration UI toast duration in seconds. Defaults to 10 for
#'   warning/error, 5 for info/success. Pass `NA` for sticky toasts.
#' @return Invisible `NULL`.
#' @export
mcview_notify <- function(severity = c("info", "success", "warning", "error"),
                          message,
                          where = c("both", "ui", "console"),
                          duration = NULL) {
    severity <- match.arg(severity)
    where <- match.arg(where)

    if (is.null(duration)) {
        duration <- if (severity %in% c("warning", "error")) 10 else 5
    }

    if (where %in% c("both", "console")) {
        switch(severity,
            info    = cli::cli_alert_info(message),
            success = cli::cli_alert_success(message),
            warning = cli::cli_alert_warning(message),
            error   = cli::cli_alert_danger(message)
        )
    }

    if (where %in% c("both", "ui") && !is.null(shiny::getDefaultReactiveDomain())) {
        ui_type <- switch(severity,
            info    = "default",
            success = "message",
            warning = "warning",
            error   = "error"
        )
        shiny::showNotification(message, type = ui_type, duration = duration)
    }

    invisible(NULL)
}
