#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
    items_list <- purrr::map(tab_defs, ~ {
        shinydashboard::menuSubItem(.x$title, tabName = .x$module_name, icon = icon(.x$icon))
    })

    sidebar <- shinydashboard::dashboardSidebar(
        div(
            id = "tab_sidebar_help",
            shinydashboard::sidebarMenu(
                id = "tab_sidebar",
                shinydashboard::menuItem("Tabs",
                    tabname = "tabs",
                    startExpanded = TRUE,
                    items_list
                )
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
            condition = "input.tab_sidebar == 'flow'",
            mod_flow_sidebar_ui("flow_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'markers'",
            mod_markers_sidebar_ui("markers_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'inner_fold'",
            mod_inner_fold_sidebar_ui("inner_fold_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'proj_fold'",
            mod_proj_fold_sidebar_ui("proj_fold_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'consistency_fold'",
            mod_consistency_fold_sidebar_ui("consistency_fold_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'query'",
            mod_query_sidebar_ui("query_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'atlas'",
            mod_atlas_sidebar_ui("atlas_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'cell_type'",
            mod_cell_type_sidebar_ui("cell_type_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'samples'",
            mod_samples_sidebar_ui("samples_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'mc_mc'",
            mod_mc_mc_sidebar_ui("mc_mc_ui_1")
        ),
        conditionalPanel(
            condition = "input.tab_sidebar == 'annotate'",
            mod_annotate_sidebar_ui("annotate_ui_1")
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
            shinydashboard::tabItem(tabName = "flow", mod_flow_ui("flow_ui_1")),
            shinydashboard::tabItem(tabName = "markers", mod_markers_ui("markers_ui_1")),
            shinydashboard::tabItem(tabName = "inner_fold", mod_inner_fold_ui("inner_fold_ui_1")),
            shinydashboard::tabItem(tabName = "proj_fold", mod_proj_fold_ui("proj_fold_ui_1")),
            shinydashboard::tabItem(tabName = "consistency_fold", mod_consistency_fold_ui("consistency_fold_ui_1")),
            shinydashboard::tabItem(tabName = "query", mod_query_ui("query_ui_1")),
            shinydashboard::tabItem(tabName = "atlas", mod_atlas_ui("atlas_ui_1")),
            shinydashboard::tabItem(tabName = "mc_mc", mod_mc_mc_ui("mc_mc_ui_1")),
            shinydashboard::tabItem(tabName = "samples", mod_samples_ui("samples_ui_1")),
            shinydashboard::tabItem(tabName = "cell_type", mod_cell_type_ui("cell_type_ui_1")),
            shinydashboard::tabItem(tabName = "annotate", mod_annotate_ui("annotate_ui_1"))
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

    screen_size_jscode <-
        '$(document).on("shiny:connected", function(e) {
         var screenWidth = screen.width;
         var screenHeight = screen.height;
        Shiny.onInputChange("screen_width", screenWidth)
        Shiny.onInputChange("screen_height", screenHeight)
        });
    '

    tagList(
        tags$script(screen_size_jscode),
        shinyjs::useShinyjs(), # Set up shinyjs
        rintrojs::introjsUI(),
        shinybusy::add_busy_spinner(spin = "breeding-rhombus", position = "bottom-right", timeout = 100),
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
