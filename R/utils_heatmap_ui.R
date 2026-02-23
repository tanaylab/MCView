heatmap_box <- function(id,
                        title = "Heatmap",
                        fold_change_range = c(-3, 3),
                        midpoint = 0,
                        low_color = "blue",
                        mid_color = "white",
                        high_color = "red",
                        highlight_color = "#fba236",
                        gene_select_label = "Select on double-click (gene)",
                        gene_select_choices = c("X axis", "Y axis"),
                        legend_width = 2,
                        height = "80vh") {
    ns <- NS(id)
    tagList(
    rclipboard::rclipboardSetup(),
    div(
        generic_box(
            id = ns("heatmap_box"),
            title = span(
                title,
                actionButton(
                    ns("show_help"),
                    icon = icon("question-circle"),
                    label = "",
                    class = "btn-link",
                    style = "padding: 0 0 0 0; color: white; font-size: 16px; background: transparent; border: none;"
                )
            ),
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            height = height,
            sidebar = shinydashboardPlus::boxSidebar(
                startOpen = FALSE,
                width = 25,
                id = ns("heatmap_sidebar"),
                checkboxInput(ns("filter_by_clipboard"), "Filter by clipboard", value = FALSE),
                shinyWidgets::numericRangeInput(ns("lfp_range"), "Fold change range", fold_change_range, width = "80%", separator = " to "),
                numericInput(ns("midpoint"), "Midpoint", midpoint),
                colourpicker::colourInput(ns("low_color"), "Low color", low_color),
                colourpicker::colourInput(ns("high_color"), "High color", high_color),
                colourpicker::colourInput(ns("mid_color"), "Mid color", mid_color),
                checkboxInput(ns("plot_legend"), "Show legend", value = TRUE),
                checkboxInput(ns("plot_cell_type_legend"), "Show cell type legend", value = TRUE),
                checkboxInput(ns("plot_genes_legend"), "Show genes legend", value = mcv_get("config")$show_heatmap_genes_legend %||% TRUE),
                colourpicker::colourInput(ns("highlight_color"), "Gene highlight color", highlight_color),
                numericInput(ns("legend_width"), "Legend width", min = 1, max = 11, step = 1, value = legend_width),
                shinyWidgets::prettyRadioButtons(
                    inputId = ns("gene_select"),
                    label = gene_select_label,
                    choices = gene_select_choices,
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                ),
                shinyWidgets::prettyRadioButtons(
                    inputId = ns("metacell_select"),
                    label = "Select on double-click (metacell)",
                    choices = c("Metacell A", "Metacell B"),
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                )
            ),
            uiOutput(ns("plotting_area"))
        ),
        style = "position:relative",
        uiOutput(ns("hover_info"), style = "pointer-events: none")
    )
    )
}

heatmap_sidebar <- function(id, ..., show_fitted_filter = FALSE) {
    ns <- NS(id)
    show_only_fitted_ui <- NULL
    config <- mcv_get("config")
    if (config$light_version) {
        max_gene_num_ui <- NULL
        remove_genes_ui <- NULL
        add_genes_ui <- NULL
        update_genes_ui <- NULL
        load_genes_ui <- NULL
        use_de_genes_ui <- NULL
        include_lateral_ui <- NULL
        include_noisy_ui <- NULL
        highlight_genes_ui <- NULL
        include_metadata_ui <- NULL
    } else {
        max_gene_num_ui <- numericInput(ns("max_gene_num"), "Maximal number of genes", value = 100)
        highlight_genes_ui <- shinyWidgets::actionGroupButtons(
            inputIds = c(ns("highlight_genes"), ns("clear_highlights")),
            labels = c("Highlight selected genes", "Clear highlights"),
            status = c("primary", "default"),
            size = "sm"
        )
        remove_genes_ui <- shinyWidgets::actionGroupButtons(ns("remove_genes"), labels = "Remove selected genes", size = "sm")
        add_genes_ui <- uiOutput(ns("add_genes_ui"))
        update_genes_ui <- shinyWidgets::actionGroupButtons(ns("update_genes"), labels = "Update genes", size = "sm")
        use_de_genes_ui <- shinyWidgets::actionGroupButtons(ns("use_de_genes"), labels = "Use DE genes", size = "sm")
        include_lateral_ui <- shinyWidgets::awesomeCheckbox(
            inputId = ns("include_lateral"),
            label = "Include lateral",
            value = TRUE
        )
        include_noisy_ui <- shinyWidgets::awesomeCheckbox(
            inputId = ns("include_noisy"),
            label = "Include noisy",
            value = TRUE
        )
        if (show_fitted_filter) {
            show_only_fitted_ui <- shinyWidgets::awesomeCheckbox(
                inputId = ns("show_only_fitted"),
                label = "Show only fitted",
                value = FALSE
            )
        }
        load_genes_ui <- fileInput(ns("load_genes"), label = NULL, buttonLabel = "Load genes", multiple = FALSE, accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            "text/tab-separated-values",
            ".csv",
            ".tsv"
        ))
        include_metadata_ui <- shinyWidgets::awesomeCheckbox(
            inputId = ns("include_metadata"),
            label = "Include metadata",
            value = FALSE
        )
    }

    list(
        uiOutput(ns("apply_heatmap_changes_ui")),
        tags$hr(),
        uiOutput(ns("reset_zoom_ui")),
        shinyWidgets::radioGroupButtons(
            inputId = ns("brush_action"),
            label = "Brush action:",
            choices = c("Zoom", "Select"),
            selected = "Zoom",
            size = "sm",
            justified = TRUE
        ),
        uiOutput(ns("mat_value_ui")),
        uiOutput(ns("copy_metacells_ui")),
        tags$hr(),
        uiOutput(ns("cell_type_list")),
        uiOutput(ns("metadata_list")),
        checkboxInput(ns("force_cell_type"), "Force cell type", value = TRUE),
        shinyWidgets::virtualSelectInput(ns("metadata_order_cell_type_var"), "Order cell types by", choices = NULL, selected = NULL, multiple = FALSE, search = TRUE),
        shinyWidgets::virtualSelectInput(ns("metadata_order_var"), "Order by", choices = NULL, selected = NULL, multiple = FALSE, search = TRUE),
        checkboxInput(ns("enable_categorical_filter"), "Enable categorical filtering", value = FALSE),
        uiOutput(ns("categorical_filter_ui")),
        tags$hr(),
        ...,
        highlight_genes_ui,
        selectInput(
            ns("selected_marker_genes"),
            "Genes",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            size = 30,
            selectize = FALSE
        ),
        remove_genes_ui,
        update_genes_ui,
        max_gene_num_ui,
        add_genes_ui,
        use_de_genes_ui,
        show_only_fitted_ui,
        include_lateral_ui,
        include_noisy_ui,
        tags$hr(),
        load_genes_ui,
        downloadButton(ns("download_genes"), "Save genes", align = "center", style = "margin: 5px 5px 5px 15px; "),
        uiOutput(ns("copy_genes_button")),
        tags$hr(),
        include_metadata_ui,
        downloadButton(ns("download_matrix"), "Download matrix", align = "center", style = "margin: 5px 5px 5px 15px; ")
    )
}
