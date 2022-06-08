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
        globals$clipboard <- character(0)
        globals$active_tabs <- config$tabs
    })

    output$menu <- shinydashboard::renderMenu({
        items_list <- purrr::map(tab_defs[globals$active_tabs], ~ {
            shinydashboard::menuSubItem(.x$title, tabName = .x$module_name, icon = icon(.x$icon))
        })

        shinydashboard::sidebarMenu(
            id = "tab_sidebar",
            shinydashboard::menuItem("Tabs",
                tabname = "tabs",
                startExpanded = TRUE,
                items_list
            )
        )
    })

    observeEvent(input$update_tabs, {
        globals$active_tabs <- c(config$tabs, setdiff(input$selected_tabs, config$tabs))
        globals$active_tabs <- globals$active_tabs[globals$active_tabs %in% input$selected_tabs]
    })

    observe({
        available_tabs <- names(tab_defs)
        if (!has_atlas(dataset())) {
            available_tabs <- available_tabs[!(available_tabs %in% c("Atlas", "Query", "Projected-fold"))]
        }
        if (!has_samples(dataset())) {
            available_tabs <- available_tabs[available_tabs != "Samples"]
        }
        if (is.null(get_mc_data(dataset(), "inner_fold_mat"))) {
            available_tabs <- available_tabs[available_tabs != "Inner-fold"]
        }
        if (is.null(get_mc_data(dataset(), "deviant_fold_mat"))) {
            available_tabs <- available_tabs[available_tabs != "Outliers"]
        }
        if (is.null(get_mc_data(dataset(), "type_flow"))) {
            available_tabs <- available_tabs[available_tabs != "Flow"]
        }
        updateCheckboxGroupInput(
            inputId = "selected_tabs",
            selected = globals$active_tabs,
            choices = available_tabs
        )
    })

    # annotation reactives
    metacell_types <- reactiveVal()
    cell_type_colors <- reactiveVal()
    gene_modules <- reactiveVal()

    observe({
        initial_cell_type_colors <- get_cell_type_data(dataset())
        initial_metacell_types <- get_metacell_types_data(dataset())
        initial_gene_modules <- get_mc_data(dataset(), "gene_modules")

        # remove metacell color column if exists
        initial_metacell_types$mc_col <- NULL

        # add cell type color from initial cell type annotation
        initial_metacell_types <- initial_metacell_types %>%
            left_join(initial_cell_type_colors %>% select(cell_type, mc_col = color), by = "cell_type")

        metacell_types(initial_metacell_types)
        cell_type_colors(initial_cell_type_colors)
        gene_modules(initial_gene_modules)
    })

    observe({
        req(gene_modules())
        if (!is.factor(gene_modules()$module)) {
            gene_modules(gene_modules() %>%
                mutate(module = factor(module)))
        }
    })

    load_tab <- function(tab_name) {
        func_name <- glue("mod_{tab_name}_server")
        if (exists(func_name)) {
            module <- get(func_name)
            module(tab_name, dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors, gene_modules = gene_modules, globals = globals)
        } else {
            warning(paste0("Tab ", tab_name, " not found"))
        }
    }

    purrr::map(tab_defs, ~ load_tab(.x$module_name))

    clipboard_reactives(dataset, input, output, session, metacell_types, cell_type_colors, gene_modules, globals)

    download_modal_reactives(input, output, session, globals)

    help_reactives(input, output, session, globals)

    # callModule(profvis::profvis_server, "profiler")
    # Rprof(strftime(Sys.time(), "%Y-%m-%d-%H-%M-%S.Rprof"),
    #     interval = 0.01, line.profiling = TRUE,
    #     gc.profiling = FALSE, memory.profiling = FALSE
    # )

    # onStop(function() {
    #     Rprof(NULL)
    # })
}
