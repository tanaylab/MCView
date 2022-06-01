download_modal_reactives <- function(input, output, session, globals) {
    observeEvent(
        input$download_modal,
        showModal(modalDialog(
            title = "Run MCView locally",
            "To run the app locally (on a unix or mac machine*), download the bundle by pressing the download button below:",
            br(),
            br(),
            downloadButton("download_bundle", "Download", style = "align-items: center;"),
            br(),
            br(),
            "Then, run the following lines in R (make sure that you are at the download directory):",
            br(),
            br(),
            glue("# install dependencies if needed"),
            br(),
            glue("if (!require('remotes')) install.packages('remotes')"),
            br(),
            glue("if (!require('MCView')) remotes::install_github('tanaylab/MCView', ref = remotes::github_release())"),
            br(),
            glue("zip::unzip('MCView-{basename(project)}.zip')"),
            br(),
            br(),
            glue("# run the app"),
            br(),
            glue("MCView::run_app('{basename(project)}', launch.browser = TRUE)"),
            br(),
            br(),
            "* It is possible to run on a windows machine using WSL",
            br(),
            easyClose = TRUE
        ))
    )

    output$download_bundle <- downloadHandler(
        filename = function() {
            glue("MCView-{project}.zip")
        },
        content = function(file) {
            download_project(file, project)
        }
    )
}
