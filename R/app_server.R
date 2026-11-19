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
        globals$mc2d <- get_mc_data(dataset(), "mc2d")
        globals$anchor_genes <- get_mc_data(dataset(), "umap_anchors")
        globals$plotly_scale <- 1
        globals$plotly_format <- "svg"
        globals$plotly_width <- NULL
        globals$plotly_height <- NULL
    })

    observe({
        globals$plotly_format <- input$plotly_format
        globals$plotly_width <- input$plotly_width
        globals$plotly_height <- input$plotly_height
        globals$plotly_scale <- input$plotly_scale
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
        globals$active_tabs <- order_tabs(globals$active_tabs)
    })

    observe({
        updateCheckboxGroupInput(
            inputId = "selected_tabs",
            selected = globals$active_tabs,
            choices = available_tabs(dataset())
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

        if (!is.null(initial_gene_modules)) {
            initial_gene_modules <- initial_gene_modules %>%
                filter(gene %in% gene_names(dataset())) %>%
                mutate(gene = as.character(gene))
        }

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

    # start a tab's server once: when the tab is first opened, so tabs the user never visits don't compute
    # anything or send data to the browser, or when another tab needs its results (globals$load_tab, e.g. the
    # heatmaps' "Use DE genes" reads the Diff. Expression comparison)
    loaded_tabs <- character(0)
    load_tab <- function(tab_name) {
        if (tab_name %in% loaded_tabs) {
            return(invisible(NULL))
        }
        loaded_tabs <<- c(loaded_tabs, tab_name)
        func_name <- glue("mod_{tab_name}_server")
        if (exists(func_name)) {
            module <- get(func_name)
            # in the app's session, also when called from inside another module
            withReactiveDomain(session, module(tab_name, dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors, gene_modules = gene_modules, globals = globals))
        } else {
            warning(paste0("Tab ", tab_name, " not found"))
        }
    }

    globals$load_tab <- load_tab
    observeEvent(input$tab_sidebar, load_tab(input$tab_sidebar))

    clipboard_reactives(dataset, input, output, session, metacell_types, cell_type_colors, gene_modules, globals)

    download_modal_reactives(input, output, session, globals)
    download_data_modal_reactives(input, output, session, globals)

    if (!is.null(config$profile) && config$profile) {
        if (!requireNamespace("profvis", quietly = TRUE)) {
            stop("Please install profvis R package in order to use profiling")
        }
        callModule(profvis::profvis_server, "profiler")
        # Rprof(strftime(Sys.time(), "%Y-%m-%d-%H-%M-%S.Rprof"),
        #     interval = 0.01, line.profiling = TRUE,
        #     gc.profiling = FALSE, memory.profiling = FALSE
        # )

        # onStop(function() {
        #     Rprof(NULL)
        # })
    }
}

# The tabs that dataset has the data for
available_tabs <- function(dataset) {
    tabs <- names(tab_defs)
    if (!has_atlas(dataset)) {
        tabs <- tabs[!(tabs %in% c("Atlas", "Query", "Projected-fold"))]
    }
    if (!has_samples(dataset)) {
        tabs <- tabs[tabs != "Samples"]
    }
    if (!has_mc_data(dataset, "inner_fold_mat")) {
        tabs <- tabs[tabs != "Inner-fold"]
    }
    if (!has_mc_data(dataset, "inner_stdev_mat")) {
        tabs <- tabs[tabs != "Stdev-fold"]
    }
    if (!has_mc_data(dataset, "deviant_fold_mat")) {
        tabs <- tabs[tabs != "Outliers"]
    }
    if (!has_mc_data(dataset, "type_flow")) {
        tabs <- tabs[tabs != "Flow"]
    }
    tabs
}
