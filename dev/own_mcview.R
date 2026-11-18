# Sourced by the scripts in dev/ that write to existing bundles.
#
# Why the bundle isn't served by the installed MCView, or NULL. A bundle whose app.R loads its own MCView from a
# `code` directory may run an MCView that fails at startup on files newer versions write (MCView < 0.2.54 can't
# read .qs2, and lists every file whose name contains ".qs", so it also tries to read default_clusterings.qs2).
# `code` counts even when it is a symlink that doesn't resolve here (e.g. to /MCView/0.2.5, which exists in the
# server's container).
own_mcview <- function(bundle) {
    app <- file.path(bundle, "app.R")
    if (!file.exists(app)) {
        return(NULL)
    }
    lines <- readLines(app, warn = FALSE)
    if (!any(grepl("MCView::run_app", lines))) {
        return("its app.R doesn't call MCView::run_app")
    }
    code <- file.path(bundle, "code")
    target <- Sys.readlink(code) # "" for a directory, NA when there is no code
    link <- !is.na(target) && nzchar(target)
    if (any(grepl("load_all", lines)) && (dir.exists(code) || link)) {
        return(paste0("its app.R loads its own MCView from code/", if (link) paste0(" (-> ", target, ")")))
    }
    NULL
}
