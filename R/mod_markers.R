#' gene_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_markers_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 12,
                shinydashboardPlus::box(
                    id = ns("markers_heatmap_box"),
                    title = "Markers Heatmap",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "80vh",
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("markers_heatmap_sidebar"),
                        checkboxInput(ns("force_cell_type"), "Force cell type", value = TRUE),
                        shinyWidgets::numericRangeInput(ns("lfp_range"), "Fold change range", c(-3, 3), width = "80%", separator = " to "),
                        checkboxInput(ns("plot_legend"), "Plot legend", value = TRUE)
                    ),
                    shinycssloaders::withSpinner(
                        plotOutput(ns("markers_heatmap"), height = "80vh")
                    )
                )
            )
        )
    )
}


#' markers sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_markers_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            shinyWidgets::actionGroupButtons(ns("apply"), labels = "Draw Heatmap"),
            uiOutput(ns("cell_type_list")),
            uiOutput(ns("metadata_list")),
            shinyWidgets::actionGroupButtons(ns("update_markers"), labels = "Update markers", size = "sm"),
            shinyWidgets::radioGroupButtons(
                inputId = ns("mode"),
                label = "Mode:",
                choices = c(
                    "Markers",
                    "Inner-folds"
                ),
                justified = TRUE
            ),
            uiOutput(ns("marker_genes_list")),
            uiOutput(ns("add_genes_ui"))
        )
    )
}

#' gene_mc Server Function
#'
#' @noRd
mod_markers_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    markers <- reactiveVal()
    markers_matrix <- reactiveVal()
    lfp_range <- reactiveVal()

    output$cell_type_list <- renderUI({
        cell_types_hex <- col2hex(get_cell_type_colors(dataset()))
        cell_types <- names(get_cell_type_colors(dataset()))
        shinyWidgets::pickerInput(ns("selected_cell_types"), "Cell types",
            choices = cell_types,
            selected = cell_types,
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `dropup-auto` = FALSE),
            choicesOpt = list(
                style = paste0("color: ", cell_types_hex, ";")
            )
        )
    })

    output$metadata_list <- renderUI({
        metadata <- get_mc_data(dataset(), "metadata")
        req(metadata)
        metadata_fields <- colnames(metadata)[-1]
        shinyWidgets::pickerInput(ns("selected_md"), "Metadata",
            choices = metadata_fields,
            selected = NULL,
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `dropup-auto` = FALSE)
        )
    })

    output$marker_genes_list <- renderUI({
        tagList(
            selectInput(
                ns("selected_marker_genes"),
                "Marker genes",
                choices = markers(),
                selected = NULL,
                multiple = TRUE,
                size = 30,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(ns("remove_genes"), labels = "Remove selected genes", size = "sm")
        )
    })

    observe({
        initial_markers <- choose_markers(get_marker_genes(dataset(), mode = input$mode), 100)
        markers(initial_markers)
        req(metacell_types())
        req(all(input$selected_cell_types %in% metacell_types()$cell_type))

        mat <- get_marker_matrix(
            dataset(),
            initial_markers,
            input$selected_cell_types,
            metacell_types(),
            force_cell_type = input$force_cell_type,
            mode = input$mode
        )
        markers_matrix(mat)

        if (input$mode == "Inner-folds") {
            lfp_range(c(0, 4))
        } else {
            lfp_range(c(-3, 3))
        }
    })

    observeEvent(input$update_markers, {
        req(metacell_types())
        req(input$selected_cell_types)
        req(all(input$selected_cell_types %in% metacell_types()$cell_type))

        markers_df <- metacell_types() %>%
            filter(cell_type %in% input$selected_cell_types) %>%
            select(metacell) %>%
            inner_join(get_marker_genes(dataset(), mode = input$mode), by = "metacell")
        new_markers <- choose_markers(markers_df, 100)

        markers(new_markers)
    })


    observeEvent(input$remove_genes, {
        new_markers <- markers()[!(markers() %in% input$selected_marker_genes)]
        markers(new_markers)
        shinyWidgets::updatePickerInput(session, ns("genes_to_add"), selected = c())
    })

    output$add_genes_ui <- renderUI({
        tagList(
            shinyWidgets::pickerInput(ns("genes_to_add"), "",
                choices = gene_names,
                selected = c(),
                multiple = TRUE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
            ),
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm")
        )
    })

    observeEvent(input$add_genes, {
        new_markers <- sort(unique(c(markers(), input$genes_to_add)))
        markers(new_markers)
        shinyWidgets::updatePickerInput(session = session, inputId = "genes_to_add", selected = character(0))
    })

    observeEvent(input$apply, {
        shinyjs::hide("markers_heatmap")
        selected_cell_types <- input$selected_cell_types %||% cell_type_colors()$cell_type
        force_cell_type <- input$force_cell_type %||% TRUE
        req(markers())
        req(metacell_types())
        req(input$mode)

        mat <- get_marker_matrix(
            dataset(),
            markers(),
            input$selected_cell_types,
            metacell_types(),
            force_cell_type = input$force_cell_type,
            mode = input$mode
        )

        markers_matrix(mat)

        lfp_range(input$lfp_range)
        shinyjs::show("markers_heatmap")
    })

    output$markers_heatmap <- renderPlot({
        req(dataset())
        req(markers_matrix())
        req(input$mode)
        req(lfp_range())
        req(metacell_types())
        req(cell_type_colors())
        req(input$plot_legend)

        if (input$mode == "Markers") {
            colors <- c("darkblue", "blue", "lightblue", "white", "red", "darkred")
            mid_color <- 4
        } else {
            colors <- c("white", "red", "darkred", "black")
            mid_color <- 1
        }

        if (!is.null(input$selected_md)) {
            metadata <- get_mc_data(dataset(), "metadata") %>%
                select(metacell, one_of(input$selected_md))
        } else {
            metadata <- NULL
        }

        plot_markers_mat(
            markers_matrix(),
            metacell_types(),
            cell_type_colors(),
            dataset(),
            min_lfp = lfp_range()[1],
            max_lfp = lfp_range()[2],
            plot_legend = input$plot_legend %||% TRUE,
            colors = colors,
            mid_color = mid_color,
            metadata = metadata
        )
    }) %>% bindCache(dataset(), markers_matrix(), metacell_types(), cell_type_colors(), lfp_range(), input$plot_legend, input$selected_md)
}

get_marker_matrix <- function(dataset, markers, cell_types = NULL, metacell_types = NULL, force_cell_type = TRUE, mode = "Markers") {
    if (mode == "Inner-folds") {
        mc_fp <- get_mc_data(dataset, "inner_fold_mat")
        mc_fp <- as.matrix(mc_fp[Matrix::rowSums(mc_fp) > 0, ])
        mc_fp <- mc_fp[intersect(markers, rownames(mc_fp)), ]
        epsilon <- 1e-5
    } else {
        mc_fp <- get_mc_fp(dataset, markers)
        epsilon <- 0
    }

    if (!is.null(cell_types)) {
        mat <- filter_mat_by_cell_types(mc_fp, cell_types, metacell_types)
    } else {
        mat <- mc_fp
    }

    if (ncol(mat) > 1) {
        mc_order <- order_mc_by_most_var_genes(mat, force_cell_type = force_cell_type, metacell_types = metacell_types, epsilon = epsilon)
        mat <- mat[, mc_order]
    }

    return(mat)
}

get_marker_genes <- function(dataset, mode = "Markers") {
    if (mode == "Inner-folds") {
        return(get_mc_data(dataset, "marker_genes_inner_fold"))
    } else {
        return(get_markers(dataset))
    }
}
