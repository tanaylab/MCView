download_modal_reactives <- function(input, output, session, globals, state) {
    project <- mcv_get("project")

    # Handle DAF mode where project might be a dummy path
    if (is.null(project) || !fs::file_exists(project)) {
        dl_text <- glue("
        Then, run the following lines in R (make sure that you are at the download directory):

        ```r
        # install dependencies if needed
        if (!require('remotes')) install.packages('remotes')
        if (!require('MCView')) remotes::install_github('tanaylab/MCView', ref = remotes::github_release())

        # run the app with your DAF object
        library(dafr)
        daf_obj <- complete_daf('your_daf_path')
        MCView::run_app(daf_obj, launch.browser = TRUE)
        ```
        ")
    } else {
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
    }

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
            if (is.null(project) || !fs::file_exists(project)) {
                "MCView-daf-project.zip"
            } else {
                glue("MCView-{basename(project)}.zip")
            }
        },
        content = function(file) {
            if (is.null(project) || !fs::file_exists(project)) {
                # Create a dummy zip file for DAF mode
                writeLines("DAF mode - no project bundle available", file)
            } else {
                download_project(file, project)
            }
        }
    )
}

download_data_modal_reactives <- function(input, output, session, globals, state) {
    project <- mcv_get("project")
    creation_text <- ""

    # Handle DAF mode where project might be a dummy path
    if (is.null(project) || !fs::file_exists(project)) {
        dialog <- modalDialog(
            title = "Download data",
            "Data download is not available in DAF mode. The data is already loaded in memory.",
            br(),
            br(),
            creation_text,
            easyClose = TRUE
        )
    } else if (fs::file_exists(command_file_path(project)) && !(!is.null(app_config("show_command")) && !app_config("show_command")) && !app_config("light_version")) {
        cmd_text <- readLines(command_file_path(project))
        cmd_text <- paste0(cmd_text, collapse = "\n")
        cmd_text <- paste0("```r\n", cmd_text, "\n```")
        creation_text <- glue("The MCView app was created using the following commands:\n\n{cmd_text}")
        creation_text <- markdown::markdownToHTML(text = creation_text, fragment.only = TRUE)
        Encoding(creation_text) <- "UTF-8"
        creation_text <- HTML(creation_text)

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
            if (is.null(project) || !fs::file_exists(project)) {
                "daf_metacells.h5ad"
            } else {
                glue("{basename(project)}_metacells.h5ad")
            }
        },
        content = function(file) {
            if (is.null(project) || !fs::file_exists(project)) {
                # Create a dummy file for DAF mode
                writeLines("DAF mode - no data file available", file)
            } else {
                fs::file_copy(source_metacells_file_path(project), file)
            }
            invisible(file)
        }
    )
}
