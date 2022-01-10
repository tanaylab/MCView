# Launch the ShinyApp

options("golem.app.prod" = TRUE)
if (dir.exists("code")){
    pkgload::load_all("code", export_all = FALSE)
}
MCView::run_app(project = "project")
