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

download_data_modal_reactives <- function(input, output, session, globals) {
    creation_text <- ""
    if (fs::file_exists(command_file_path(project)) && !(!is.null(config$show_command) && !config$show_command) && !config$light_version) {
        cmd_text <- readLines(command_file_path(project))
        cmd_text <- paste0(cmd_text, collapse = "\n")
        cmd_text <- paste0("```r\n", cmd_text, "\n```")
        creation_text <- glue("The MCView app was created using the following commands:\n\n{cmd_text}")
        creation_text <- markdown::markdownToHTML(text = creation_text, fragment.only = TRUE)
        Encoding(creation_text) <- "UTF-8"
        creation_text <- HTML(creation_text)
    }

    if (file.exists(source_metacells_file_path(project))) {
        file_size <- fs::file_info(source_metacells_file_path(project), follow = TRUE)$size
        if (!is.na(file_size)) {
            dialog <- modalDialog(
                title = "Download data",
                "To download the metacells anndata file, press the download button below:",
                br(),
                br(),
                downloadButton("download_mc_data", glue("Download ({file_size})"), style = "align-items: center;"),
                br(),
                br(),
                creation_text,
                easyClose = TRUE
            )
        }
    } else {
        dialog <- modalDialog(
            title = "Download data",
            "The metacells anndata file is not available. Run the import commands with copy_source_file=TRUE to make it available.",
            br(),
            br(),
            creation_text,
            easyClose = TRUE
        )
    }

    observeEvent(
        input$download_data_modal,
        showModal(dialog)
    )

    output$download_mc_data <- downloadHandler(
        filename = function() {
            glue("{basename(project)}_metacells.h5ad")
        },
        content = function(file) {
            fs::file_copy(source_metacells_file_path(project), file)
            invisible(file)
        }
    )
}
