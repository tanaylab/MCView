#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
    if (length(dataset_ls(project)) > 1) {
        dataset <- reactive(input$dataset)
    } else {
        dataset <- function() {
            dataset_ls(project)[1]
        }
    }

    # annotation reactives
    metacell_types <- reactiveVal()
    cell_type_colors <- reactiveVal()

    callModule(mod_manifold_server, "manifold_ui_1", dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors)
    callModule(mod_gene_mc_server, "gene_mc_ui_1", dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors)
    callModule(mod_mc_mc_server, "mc_mc_ui_1", dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors)
    callModule(mod_annotate_mc_server, "annotate_ui_1", dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors)

    show_help <- function(input, output, session) {
        if (input$tab_sidebar == "about") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = FALSE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("about")
                )
            )
        } else if (input$tab_sidebar == "gene_mc") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("genes", "gene_mc_ui_1")
                )
            )
        } else if (input$tab_sidebar == "mc_mc") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("metacells", "mc_mc_ui_1")
                )
            )
        } else if (input$tab_sidebar == "annotate") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("metacells", "annotate_ui_1")
                )
            )
        }
    }

    observe({
        req(config$help)
        rintrojs::introjs(session,
            options = list(
                "showProgress" = FALSE,
                "showBullets" = FALSE,
                "nextLabel" = "next",
                "prevLabel" = "back",
                "scrollToElement" = TRUE,
                steps =
                    rbind(
                        data.frame(element = NA, intro = help_config$introduction),
                        get_tab_steps("about")
                    )
            )
        )
    })

    observeEvent(
        input$help,
        show_help(input, output, session)
    )

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
            glue("zip::unzip('MCView-{project}.zip')"),
            br(),
            br(),
            glue("# run the app"),
            br(),
            glue("MCView::run_app('{project}', launch.browser = TRUE)"),
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
