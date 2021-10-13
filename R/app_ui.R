#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
    items_list <- purrr::map(tab_defs, ~ {
        shinydashboard::menuItem(.x$title, tabName = .x$module_name, icon = icon(.x$icon))
    })

    sidebar <- shinydashboard::dashboardSidebar(
        div(
            id = "tab_sidebar_help",
            shinydashboard::sidebarMenu(
                id = "tab_sidebar",
                items_list
            )
        ),
        tags$hr(),
        conditionalPanel(
            condition = "input.tab_sidebar == 'about'",
            mod_about_sidebar_ui("about_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'manifold'",
            mod_manifold_sidebar_ui("manifold_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'gene_mc'",
            mod_gene_mc_sidebar_ui("gene_mc_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'mc_mc'",
            mod_mc_mc_sidebar_ui("mc_mc_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'annotate'",
            mod_annotate_mc_sidebar_ui("annotate_ui_1")
        )
    )

    if (length(dataset_ls(project)) > 1) {
        right_sidebar <- shinydashboardPlus::dashboardControlbar(
            selectInput("dataset", label = "Dataset", choices = dataset_ls(project), selected = dataset_ls(project)[1], multiple = FALSE, selectize = FALSE)
        )
    } else {
        right_sidebar <- NULL
    }


    body <- shinydashboard::dashboardBody(
        shinydashboard::tabItems(
            shinydashboard::tabItem(tabName = "about", mod_about_ui("about_ui_1")),
            shinydashboard::tabItem(tabName = "manifold", mod_manifold_ui("manifold_ui_1")),
            shinydashboard::tabItem(tabName = "gene_mc", mod_gene_mc_ui("gene_mc_ui_1")),
            shinydashboard::tabItem(tabName = "mc_mc", mod_mc_mc_ui("mc_mc_ui_1")),
            shinydashboard::tabItem(tabName = "annotate", mod_annotate_mc_ui("annotate_ui_1"))
        )
    )

    app_title <- config$title
    if (is.null(app_title) || app_title == "MCView") {
        app_title <- glue("MCView {version}", version = packageVersion("MCView"))
    }

    dashboard_page <- shinydashboardPlus::dashboardPage(
        shinydashboardPlus::dashboardHeader(
            title = app_title,
            tags$li(
                title = "",
                class = "dropdown",
                shinyWidgets::actionBttn(
                    inputId = "download_modal",
                    label = "Run locally",
                    icon = NULL,
                    # style = "fill",
                    color = "default",
                    size = "sm",
                    block = FALSE,
                    no_outline = FALSE
                ),
                shinyWidgets::actionBttn(
                    inputId = "help",
                    label = "Help",
                    icon = NULL,
                    # style = "fill",
                    color = "default",
                    size = "sm",
                    block = FALSE,
                    no_outline = FALSE
                )
            )
        ),
        sidebar = sidebar,
        body = body,
        controlbar = right_sidebar
    )

    tagList(
        rintrojs::introjsUI(),
        dashboard_page
    )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
    add_resource_path(
        "www", app_sys("app/www")
    )

    tags$head(
        favicon(),
        bundle_resources(
            path = app_sys("app/www"),
            app_title = config$title
        )
        # Add here other external resources
        # for example, you can add shinyalert::useShinyalert()
    )
}
