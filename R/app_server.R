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
        globals$active_tabs <- input$selected_tabs
    })

    observe({
        updateCheckboxGroupInput(
            inputId = "selected_tabs",
            selected = globals$active_tabs
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
        module <- get(glue("mod_{tab_name}_server"))
        module(tab_name, dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors, gene_modules = gene_modules, globals = globals)
    }

    load_tab("manifold")
    load_tab("gene_mc")
    load_tab("flow")
    load_tab("markers")
    load_tab("gene_modules")
    load_tab("inner_fold")
    load_tab("outliers")
    load_tab("samples")
    load_tab("cell_type")

    if (any_has_atlas(project)) {
        load_tab("query")
        load_tab("atlas")
        load_tab("proj_fold")
    }

    load_tab("mc_mc")
    load_tab("annotate")
    load_tab("about")


    clipboard_reactives(dataset, input, output, session, metacell_types, cell_type_colors, gene_modules, globals)

    download_modal_reactives(input, output, session, globals)

    help_reactives(input, output, session, globals)
}
