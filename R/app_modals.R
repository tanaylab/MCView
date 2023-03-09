download_modal_reactives <- function(input, output, session, globals) {
    dl_text <- glue("
    Then, run the following lines in R (make sure that you are at the download directory):

    ```r
    # install dependencies if needed
    if (!require('remotes')) install.packages('remotes')
    if (!require('MCView')) remotes::install_github('tanaylab/MCView', ref = remotes::github_release())

    # unzip the project
    zip::unzip('MCView-{basename(project)}.zip')

    # run the app
    MCView::run_app('{basename(project)}', launch.browser = TRUE)
    ```
    ")
    dl_text <- markdown::markdownToHTML(text = dl_text, fragment.only = TRUE)
    Encoding(dl_text) <- "UTF-8"
    dl_text <- HTML(dl_text)

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
            dl_text,
            "* It is possible to run on a windows machine using WSL",
            easyClose = TRUE
        ))
    )

    output$download_bundle <- downloadHandler(
        filename = function() {
            glue("MCView-{basename(project)}.zip")
        },
        content = function(file) {
            download_project(file, project)
        }
    )
}
