#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
    if (length(dataset_ls()) > 1) {
        dataset <- reactive(input$dataset)
    } else {
        dataset <- function() {
            dataset_ls()[1]
        }
    }

    # annotation reactives
    mc_annot <- reactiveVal()
    cell_type_annot <- reactiveVal()

    callModule(mod_manifold_server, "manifold_ui_1", dataset = dataset, mc_annot = mc_annot, cell_type_annot = cell_type_annot)
    callModule(mod_gene_mc_server, "gene_mc_ui_1", dataset = dataset, mc_annot = mc_annot, cell_type_annot = cell_type_annot)
    callModule(mod_mc_mc_server, "mc_mc_ui_1", dataset = dataset, mc_annot = mc_annot, cell_type_annot = cell_type_annot)
    callModule(mod_annotate_mc_server, "annotate_ui_1", dataset = dataset, mc_annot = mc_annot, cell_type_annot = cell_type_annot)

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
                        # get_tab_steps("genes", "gene_mc_ui_1"),
                        # get_tab_steps("metacells", "mc_mc_ui_1")
                    )
            )
        )
    })

    observeEvent(
        input$help,
        show_help(input, output, session)
    )
}
