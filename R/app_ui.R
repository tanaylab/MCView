#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
    modules <- purrr::map_chr(tab_defs, "module_name")

    ui_sidebar_funcs <- purrr::map(modules, ~ {
        func_name <- glue("mod_{.x}_sidebar_ui")
        if (!exists(func_name)) {
            return(NULL)
        } else {
            return(get(func_name))
        }
    })

    cond_panels <- purrr::map2(
        modules,
        ui_sidebar_funcs,
        function(mod, func) {
            if (is.null(func)) {
                return(NULL)
            } else {
                conditionalPanel(
                    condition = glue("input.tab_sidebar == '{mod}'"),
                    func(mod)
                )
            }
        }
    )

    sidebar <- do.call(
        shinydashboard::dashboardSidebar,
        c(
            list(div(
                id = "tab_sidebar_help",
                shinydashboard::sidebarMenuOutput("menu")
            )),
            list(tags$hr()),
            cond_panels
        )
    )

    right_sidebar <- shinydashboardPlus::dashboardControlbar(
        width = 400,
        shinydashboardPlus::controlbarMenu(
            shinydashboardPlus::controlbarItem(
                "Tabs",
                checkboxGroupInput(
                    "selected_tabs",
                    label = "Tabs",
                    choices = names(tab_defs),
                    selected = intersect(config$tabs, names(tab_defs)),
                ),
                actionButton("update_tabs", "Update Tabs")
            ),
            shinydashboardPlus::controlbarItem(
                "Datasets",
                selectInput("dataset", label = "Dataset", choices = dataset_ls(project), selected = dataset_ls(project)[1], multiple = FALSE, selectize = FALSE)
            ),
            shinydashboardPlus::controlbarItem(
                "Theme",
                shinydashboardPlus::skinSelector()
            )
        )
    )

    ui_funcs <- purrr::map(modules, ~ {
        func_name <- glue("mod_{.x}_ui")
        if (!exists(func_name)) {
            return(NULL)
        } else {
            return(get(func_name))
        }
    })

    tab_items <- div( # shinydashboard::tabItems
        class = "tab-content",
        purrr::map2(
            modules,
            ui_funcs,
            function(mod, func) {
                if (!is.null(func)) {
                    shinydashboard::tabItem(tabName = mod, func(mod))
                }
            }
        )
    )

    body <- shinydashboard::dashboardBody(
        tab_items
    )

    app_title <- config$title
    if (is.null(app_title) || app_title == "MCView") {
        if (length(dataset_ls(project)) == 1) {
            app_title <- dataset_ls(project)[1]
            config$title <- app_title
        } else {
            app_title <- glue("MCView {version}", version = utils::packageVersion("MCView"))
        }
    }

    app_footer <- glue("MCView {version}", version = utils::packageVersion("MCView"))
    if (!is.null(config$metacells_version)) {
        app_footer <- glue("{app_footer} | {config$metacells_version}")
    }

    dashboard_page <- shinydashboardPlus::dashboardPage(
        title = app_title,
        shinydashboardPlus::dashboardHeader(
            title = app_title,
            tags$li(
                title = "",
                class = "dropdown",
                style = "margin-top: 7.5px; margin-bottom: 7.5px;",
                shinyWidgets::actionBttn(
                    inputId = "clipboard_modal",
                    label = "Clipboard",
                    icon = NULL,
                    color = "default",
                    size = "sm",
                    block = FALSE,
                    no_outline = FALSE
                ),
                shinyWidgets::actionBttn(
                    inputId = "download_modal",
                    label = "Run locally",
                    icon = NULL,
                    color = "default",
                    size = "sm",
                    block = FALSE,
                    no_outline = FALSE
                ),
                shinyWidgets::actionBttn(
                    inputId = "help",
                    label = "Help",
                    icon = NULL,
                    color = "default",
                    size = "sm",
                    block = FALSE,
                    no_outline = FALSE
                )
            )
        ),
        sidebar = sidebar,
        body = body,
        controlbar = right_sidebar,
        preloader = list(html = tagList(waiter::spin_1(), h4("Loading MCView..."), h5(config$title)), color = "#333e48"),
        footer = shinydashboardPlus::dashboardFooter(
            left = app_footer,
            right = glue("(C) Weizmann Institute of Science, 2020-{year}", year = format(Sys.time(), "%Y"))
        ),
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
        # profvis::profvis_ui("profiler"),
        tags$script(screen_size_jscode),
        shinyjs::useShinyjs(), # Set up shinyjs
        rintrojs::introjsUI(),
        golem_add_external_resources(),
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
