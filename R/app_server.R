#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
    if (length(dataset_names()) > 1) {
        dataset <- reactive(input$dataset)
    } else {
        dataset <- function() {
            dataset_names()[1]
        }
    }

    globals <- reactiveValues()

    observe({
        globals$screen_width <- input$screen_width
        globals$screen_height <- input$screen_height
    })

    observe({
        globals$mc2d <- get_mc_data(dataset(), "mc2d")
        globals$anchor_genes <- get_mc_data(dataset(), "umap_anchors")
    })

    globals$clipboard <- character(0)
    globals$active_tabs <- app_config("tabs")
    globals$plotly_scale <- 1
    globals$plotly_format <- "svg"
    globals$plotly_width <- NULL
    globals$plotly_height <- NULL

    observe({
        globals$plotly_format <- input$plotly_format
        globals$plotly_width <- input$plotly_width
        globals$plotly_height <- input$plotly_height
        globals$plotly_scale <- input$plotly_scale
    })

    # Track current sidebar tab so modules can defer computation until visited
    observe({
        globals$current_tab <- input$tab_sidebar
    })

    output$menu <- shinydashboard::renderMenu({
        items_list <- purrr::map(mcv_get("tab_defs")[globals$active_tabs], ~ {
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
        globals$active_tabs <- c(app_config("tabs"), setdiff(input$selected_tabs, app_config("tabs")))
        globals$active_tabs <- globals$active_tabs[globals$active_tabs %in% input$selected_tabs]
        globals$active_tabs <- order_tabs(globals$active_tabs)
    })

    observe({
        available_tabs <- names(mcv_get("tab_defs"))
        if (!has_atlas(dataset())) {
            available_tabs <- available_tabs[!(available_tabs %in% c("Atlas", "Query", "Projected-fold"))]
        }
        if (!has_samples(dataset()) && !("Samples" %in% mcv_get("config")$tabs)) {
            available_tabs <- available_tabs[available_tabs != "Samples"]
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

        if (!is.null(initial_gene_modules)) {
            initial_gene_modules <- initial_gene_modules %>%
                filter(gene %in% gene_names(dataset())) %>%
                mutate(gene = as.character(gene))
        }

        # Only proceed if we have valid data
        if (!is.null(initial_metacell_types) && !is.null(initial_cell_type_colors)) {
            # remove metacell color column if exists
            if ("mc_col" %in% colnames(initial_metacell_types)) {
                initial_metacell_types$mc_col <- NULL
            }

            # add cell type color from initial cell type annotation
            initial_metacell_types <- initial_metacell_types %>%
                left_join(initial_cell_type_colors %>% select(cell_type, mc_col = color), by = "cell_type")
        }

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

    load_tab <- function(tab_def) {
        server_fn <- tab_def$server_fn
        if (!is.null(server_fn)) {
            server_fn(tab_def$module_name, dataset = dataset, metacell_types = metacell_types, cell_type_colors = cell_type_colors, gene_modules = gene_modules, globals = globals)
        } else {
            warning(paste0("Tab ", tab_def$module_name, " has no server_fn"))
        }
    }

    purrr::walk(mcv_get("tab_defs"), load_tab)

    clipboard_reactives(dataset, input, output, session, metacell_types, cell_type_colors, gene_modules, globals)

    download_modal_reactives(input, output, session, globals)
    download_data_modal_reactives(input, output, session, globals)

    if (!is.null(app_config("profile")) && app_config("profile")) {
        # Profiling UI disabled; keep profile flag for timing logs only.
    }

    # Per-session cleanup: release resources when a browser tab closes
    session$onSessionEnded(function() {
        tryCatch(
            {
                cli::cli_alert_info("Shiny session ended, cleaning up session resources")
                # Trigger R garbage collection to release any session-held
                # DAF object references promptly rather than deferring to exit
                gc(verbose = FALSE)
            },
            error = function(e) {
                # Suppress errors during session cleanup
            }
        )
    })
}
