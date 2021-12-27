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
            shinyWidgets::radioGroupButtons(
                inputId = ns("mode"),
                label = "Mode:",
                choices = c(
                    "Markers",
                    "Inner",
                    "Proj"
                ),
                justified = TRUE
            ),
            uiOutput(ns("cell_type_list")),
            uiOutput(ns("metadata_list")),
            shinyWidgets::actionGroupButtons(ns("update_markers"), labels = "Update markers", size = "sm"),
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

    markers <- reactiveValues()
    lfp_range <- reactiveValues()

    output$cell_type_list <- cell_type_selector(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = "all", cell_type_colors = cell_type_colors())

    output$metadata_list <- metadata_selector(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata")

    output$marker_genes_list <- renderUI({
        tagList(
            selectInput(
                ns("selected_marker_genes"),
                "Marker genes",
                choices = markers[[input$mode]],
                selected = NULL,
                multiple = TRUE,
                size = 30,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(ns("remove_genes"), labels = "Remove selected genes", size = "sm")
        )
    })

    observe({
        req(input$mode)
        if (is.null(markers[[input$mode]])) {
            initial_markers <- choose_markers(get_marker_genes(dataset(), mode = input$mode), 100)
            markers[[input$mode]] <- initial_markers
        }

        lfp_range[[input$mode]] <- input$lfp_range
    })

    markers_matrix <- reactive({
        req(input$mode)
        req(markers[[input$mode]])
        req(metacell_types())
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

        get_marker_matrix(
            dataset(),
            markers[[input$mode]],
            input$selected_cell_types,
            metacell_types(),
            force_cell_type = input$force_cell_type,
            mode = input$mode,
            notify_var_genes = TRUE
        )
    })


    observeEvent(input$update_markers, {
        req(metacell_types())
        req(input$mode)
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

        if (!is.null(input$selected_cell_types)) {
            markers_df <- metacell_types() %>%
                filter(cell_type %in% input$selected_cell_types)
        } else {
            markers_df <- metacell_types()
        }

        markers_df <- markers_df %>%
            select(metacell) %>%
            inner_join(get_marker_genes(dataset(), mode = input$mode), by = "metacell")
        new_markers <- choose_markers(markers_df, 100)

        markers[[input$mode]] <- new_markers
    })


    observeEvent(input$remove_genes, {
        req(markers[[input$mode]])
        new_markers <- markers[[input$mode]][!(markers[[input$mode]] %in% input$selected_marker_genes)]
        markers[[input$mode]] <- new_markers
        shinyWidgets::updatePickerInput(session, ns("genes_to_add"), selected = c())
    })

    output$add_genes_ui <- renderUI({
        tagList(
            shinyWidgets::pickerInput(ns("genes_to_add"), "",
                choices = gene_names(dataset()),
                selected = c(),
                multiple = TRUE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
            ),
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm")
        )
    })

    observeEvent(input$add_genes, {
        markers[[input$mode]]
        new_markers <- sort(unique(c(markers[[input$mode]], input$genes_to_add)))
        markers[[input$mode]] <- new_markers
        shinyWidgets::updatePickerInput(session = session, inputId = "genes_to_add", selected = character(0))
    })

    output$markers_heatmap <- renderPlot({
        req(dataset())
        req(input$mode)
        req(lfp_range[[input$mode]])
        req(metacell_types())
        req(cell_type_colors())

        mat <- markers_matrix()
        req(mat)

        if (input$mode == "Markers") {
            mat <- log2(mat)
        }

        # if (input$mode == "Markers") {
        #     colors <- c("darkblue", "blue", "lightblue", "white", "red", "darkred")
        #     mid_color <- 4
        # } else {
        #     colors <- c("white", "red", "darkred", "black")
        #     mid_color <- 1
        # }

        if (!is.null(input$selected_md)) {
            metadata <- get_mc_data(dataset(), "metadata") %>%
                select(metacell, one_of(input$selected_md))
        } else {
            metadata <- NULL
        }


        plot_markers_mat(
            mat,
            metacell_types(),
            cell_type_colors(),
            dataset(),
            min_lfp = lfp_range[[input$mode]][1],
            max_lfp = lfp_range[[input$mode]][2],
            plot_legend = input$plot_legend %||% TRUE,
            # colors = colors,
            # mid_color = mid_color,
            metadata = metadata
        )
    }) %>% bindCache(dataset(), metacell_types(), cell_type_colors(), lfp_range[[input$mode]], input$plot_legend, input$selected_md, input$mode, markers[[input$mode]], input$selected_cell_types, input$force_cell_type)
}

get_marker_matrix <- function(dataset, markers, cell_types = NULL, metacell_types = NULL, force_cell_type = TRUE, mode = "Markers", notify_var_genes = FALSE) {
    if (mode == "Inner") {
        mc_fp <- get_mc_data(dataset, "inner_fold_mat")
        req(mc_fp)
        mc_fp <- as.matrix(mc_fp[Matrix::rowSums(mc_fp) > 0, ])
        epsilon <- 1e-5
    } else if (mode == "Proj") {
        mc_fp <- get_mc_data(dataset, "projected_fold")
        req(mc_fp)
        mc_fp <- as.matrix(mc_fp[Matrix::rowSums(mc_fp) > 0, ])
        mc_fp <- mc_fp[intersect(markers, rownames(mc_fp)), ]
        epsilon <- 1e-5
    } else {
        mc_fp <- get_mc_fp(dataset, markers)
        epsilon <- 0
    }

    req(dim(mc_fp))

    if (!is.null(cell_types)) {
        mat <- filter_mat_by_cell_types(mc_fp, cell_types, metacell_types)
    } else {
        mat <- mc_fp
    }

    if (ncol(mat) > 1) {
        mc_order <- order_mc_by_most_var_genes(mat, force_cell_type = force_cell_type, metacell_types = metacell_types, epsilon = epsilon, notify_var_genes = notify_var_genes)
        mat <- mat[, mc_order]
    }

    return(mat)
}

get_marker_genes <- function(dataset, mode = "Markers") {
    if (mode == "Inner") {
        if (is.null(get_mc_data(dataset, "inner_fold_mat"))) {
            showNotification(glue("Inner-fold matrix was not computed. Please compute it in python using the metacells package and rerun the import"), type = "error")
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_inner_fold"))
    } else if (mode == "Proj") {
        if (is.null(get_mc_data(dataset, "projected_fold"))) {
            showNotification(glue("Proj-fold matrix was not computed. Please compute it in python using the metacells package and rerun the import"), type = "error")
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_projected"))
    } else {
        return(get_markers(dataset))
    }

    return(get_markers(dataset))
}
