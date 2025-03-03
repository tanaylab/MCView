# Launch the ShinyApp

options("golem.app.prod" = TRUE)
if (dir.exists("code")) {
    pkgload::load_all("code", export_all = FALSE)
}

if (dir.exists("project")) { # here for backward compatibility
    MCView::run_app("project")
} else {
    MCView::run_app(getwd())
}


