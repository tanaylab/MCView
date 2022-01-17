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

    globals <- reactiveValues()

    observe({
        globals$screen_width <- input$screen_width
        globals$screen_height <- input$screen_height
    })

    # annotation reactives
    metacell_types <- reactiveVal()
    cell_type_colors <- reactiveVal()

    observe({
        initial_cell_type_colors <- get_mc_data(dataset(), "cell_type_colors")
        initial_metacell_types <- get_mc_data(dataset(), "metacell_types")

        # remove metacell color column if exists
        initial_metacell_types$mc_col <- NULL

        # add cell type color from initial cell type annotation
        initial_metacell_types <- initial_metacell_types %>%
            left_join(initial_cell_type_colors %>% select(cell_type, mc_col = color), by = "cell_type")

        metacell_types(initial_metacell_types)
        cell_type_colors(initial_cell_type_colors)
    })

    load_tab <- function(tab_name, module) {
        callModule(module, glue("{tab_name}_ui_1"), dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors, globals = globals)
    }

    load_tab("manifold", mod_manifold_server)
    load_tab("gene_mc", mod_gene_mc_server)
    load_tab("markers", mod_markers_server)
    load_tab("inner_fold", mod_inner_fold_server)
    load_tab("samples", mod_samples_server)
    
    if (any_has_atlas(project)) {
        load_tab("query", mod_query_server)
        load_tab("atlas", mod_atlas_server)
        load_tab("proj_fold", mod_proj_fold_server)
    }

    load_tab("mc_mc", mod_mc_mc_server)
    load_tab("annotate", mod_annotate_server)
    load_tab("about", mod_about_server)

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
        } else if (input$tab_sidebar == "manifold") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("manifold", "manifold_ui_1")
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
        } else if (input$tab_sidebar == "samples") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("samples", "samples_ui_1")
                )
            )
        } else if (input$tab_sidebar == "Query") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("query", "query_ui_1")
                )
            )
        } else if (input$tab_sidebar == "Atlas") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("atlas", "atlas_ui_1")
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

    # Rprof(strftime(Sys.time(), "%Y-%m-%d-%H-%M-%S.Rprof"),
    #     interval = 0.01, line.profiling = TRUE,
    #     gc.profiling = FALSE, memory.profiling = FALSE
    # )

    # onStop(function() {
    #     Rprof(NULL)
    # })
}
