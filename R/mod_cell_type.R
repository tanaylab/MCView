#' cell type UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_cell_type_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_column(
                width = 12,
                generic_box(
                    id = ns("boxplot_box"),
                    title = "Cell types",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    height = "70vh",
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 50,
                        id = ns("cell_type_sidebar"),
                        shinyWidgets::switchInput(ns("custom_ylim"), "Custom Y limits", value = FALSE, onLabel = "Yes", offLabel = "No", onStatus = "success", offStatus = "danger", size = "mini"),
                        conditionalPanel(
                            condition = paste0("input['", ns("custom_ylim"), "'] == true"),
                            fluidRow(
                                column(6, numericInput(ns("ylim_min"), "Y min:", value = NULL, step = 0.1)),
                                column(6, numericInput(ns("ylim_max"), "Y max:", value = NULL, step = 0.1))
                            ),
                            uiOutput(ns("ylim_warning"))
                        )
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("cell_type_boxplot"), height = "70vh")
                    )
                )
            )
        ),
        fluidRow(
            generic_column(
                width = 12,
                generic_box(
                    title = "Advanced Options",
                    status = "info",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    collapsed = TRUE,
                    closable = FALSE,
                    width = 12,
                    fluidRow(
                        column(
                            3,
                            shinyWidgets::prettyRadioButtons(
                                ns("plot_type"),
                                label = "Plot type:",
                                choices = c("Boxplot" = "boxplot", "Violin plot" = "violin", "Sina plot" = "sina"),
                                selected = "sina",
                                inline = FALSE,
                                status = "primary",
                                fill = TRUE
                            )
                        ),
                        column(
                            3,
                            shinyWidgets::prettyRadioButtons(
                                ns("facet_by"),
                                label = "Facet by:",
                                choices = c("None" = "none", "Cell type" = "cell_type", "Categorical metadata" = "metadata"),
                                selected = "none",
                                inline = FALSE,
                                status = "warning",
                                fill = TRUE
                            )
                        ),
                        column(
                            2,
                            uiOutput(ns("facet_metadata_selector"))
                        ),
                        column(
                            2,
                            shinyWidgets::switchInput(
                                ns("coord_flip"),
                                "Flip coordinates",
                                value = FALSE,
                                onLabel = "Yes",
                                offLabel = "No",
                                onStatus = "info",
                                offStatus = "secondary",
                                size = "mini"
                            ),
                            br(),
                            shinyWidgets::switchInput(
                                ns("log_scale"),
                                "Log scale Y-axis",
                                value = FALSE,
                                onLabel = "Yes",
                                offLabel = "No",
                                onStatus = "success",
                                offStatus = "secondary",
                                size = "mini"
                            )
                        ),
                        column(
                            2,
                            div(
                                style = "margin-top: 20px; font-size: 12px; color: #666;",
                                strong("Plot Controls:"), br(),
                                "• Plot type affects visualization", br(),
                                "• Faceting splits into panels", br(),
                                "• Flip coords for long labels", br(),
                                "• Y limits in plot sidebar"
                            )
                        )
                    )
                )
            )
        )
    )
}


#' cell type sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_cell_type_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            axis_selector("boxplot_axis", "Gene", ns, choices = c("Metadata", "Gene", "Gene module"), orientation = "vertical", wrap_in_box = FALSE),
            shinyWidgets::prettyRadioButtons(
                ns("x_axis_type"),
                label = "X axis:",
                choices = c("Cell types" = "cell_types", "Categorical metadata" = "metadata"),
                selected = "cell_types",
                inline = TRUE,
                status = "info",
                fill = TRUE
            ),
            uiOutput(ns("x_axis_metadata_selector")),
            uiOutput(ns("x_axis_categories_selector")),
            uiOutput(ns("cell_type_selector_ui")),
            uiOutput(ns("confusion_color_by_selector")),
            shinyWidgets::switchInput(ns("show_correlations"), "Show correlations", value = FALSE, onLabel = "Yes", offLabel = "No", onStatus = "success", offStatus = "danger", size = "mini"),
            uiOutput(ns("top_correlated_select_boxplot_axis"))
        )
    )
}

#' cell type Server Function
#'
#' @noRd
mod_cell_type_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            top_correlated_selector("boxplot_axis_var", "boxplot_axis", "boxplot_axis_type", input, output, session, dataset, ns, button_labels = c("Select"), ids = c("boxplot"), gene_modules = gene_modules, metacell_types = metacell_types, selected_cell_types = "boxplot_cell_types")

            render_axis_select_ui("boxplot_axis", "Data", "boxplot_axis_select", md_choices = dataset_metadata_fields(dataset()), md_selected = dataset_metadata_fields(dataset())[1], selected_gene = default_gene1, input = input, output = output, ns = ns, dataset = dataset, gene_modules = gene_modules, session = session)

            # Dynamic cell type selector based on x-axis type
            output$cell_type_selector_ui <- renderUI({
                req(input$x_axis_type)

                if (input$x_axis_type == "cell_types") {
                    # For cell types on x-axis, we need to call the cell_type_selector function
                    # But we can't nest renderUI calls, so we'll create the UI directly here
                    req(metacell_types())
                    req(cell_type_colors())

                    cell_types <- unique(metacell_types()$cell_type)
                    cell_types <- cell_types[!is.na(cell_types)]

                    # Create a simple version similar to the original cell_type_selector
                    tagList(
                        div(
                            title = "Select which cell types to include in the analysis. Deselecting cell types will exclude them from the visualization.",
                            style = "cursor: help;",
                            checkboxGroupInput(ns("boxplot_cell_types"),
                                "Cell types:",
                                choices = as.character(cell_types),
                                selected = as.character(cell_types)
                            )
                        ),
                        div(
                            style = "margin-top: 5px;",
                            div(
                                title = "Select all available cell types",
                                style = "cursor: help; display: inline-block;",
                                actionButton(ns("select_all_cell_types"), "Select All", size = "xs", style = "margin-right: 5px;")
                            ),
                            div(
                                title = "Deselect all cell types",
                                style = "cursor: help; display: inline-block;",
                                actionButton(ns("clear_all_cell_types"), "Clear All", size = "xs")
                            )
                        )
                    )
                } else {
                    # Cell type filter for metadata x-axis (similar to category selector)
                    req(metacell_types())
                    cell_types <- unique(metacell_types()$cell_type)
                    cell_types <- cell_types[!is.na(cell_types)]

                    div(
                        title = "Filter the analysis to include only specific cell types when using metadata on X-axis. This controls which cell types contribute to each metadata category.",
                        style = "cursor: help;",
                        checkboxGroupInput(ns("cell_type_filter"),
                            "Filter by cell types:",
                            choices = as.character(cell_types),
                            selected = as.character(cell_types)
                        )
                    )
                }
            })

            output$confusion_color_by_selector <- renderUI({
                div(
                    title = "Choose normalization for the confusion matrix: 'X axis' shows proportions within each X category, 'Y axis' shows proportions within each Y category.",
                    style = "cursor: help;",
                    shinyWidgets::prettyRadioButtons(
                        ns("confusion_color_by"),
                        label = "Normalize by:",
                        choices = c("X axis", "Y Axis"),
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    )
                )
            })

            # X-axis metadata selector
            output$x_axis_metadata_selector <- renderUI({
                req(input$x_axis_type == "metadata")
                metadata <- get_mc_data(dataset(), "metadata")
                req(metadata)

                # Get categorical metadata columns
                categorical_cols <- names(metadata)[sapply(metadata, function(x) !is.numeric(x) && !is.logical(x))]
                categorical_cols <- categorical_cols[categorical_cols != "metacell"]

                if (length(categorical_cols) > 0) {
                    div(
                        title = "Choose which categorical metadata variable to use for grouping on the X-axis (e.g., treatment, condition, batch, tissue type).",
                        style = "cursor: help;",
                        selectInput(ns("x_axis_metadata_var"),
                            "Select metadata variable:",
                            choices = categorical_cols,
                            selected = categorical_cols[1]
                        )
                    )
                } else {
                    tags$div(style = "color: orange;", "No categorical metadata variables found")
                }
            })

            # X-axis categories selector
            output$x_axis_categories_selector <- renderUI({
                req(input$x_axis_type == "metadata")
                req(input$x_axis_metadata_var)
                metadata <- get_mc_data(dataset(), "metadata")
                req(metadata)
                req(input$x_axis_metadata_var %in% colnames(metadata))

                categories <- unique(metadata[[input$x_axis_metadata_var]])
                categories <- categories[!is.na(categories)]

                if (length(categories) > 0) {
                    tagList(
                        div(
                            title = "Select which specific categories from the metadata variable to include in the plot. Useful for comparing specific subsets.",
                            style = "cursor: help;",
                            checkboxGroupInput(ns("x_axis_categories"),
                                "Select categories to plot:",
                                choices = as.character(categories),
                                selected = as.character(categories)
                            )
                        ),
                        div(
                            style = "margin-bottom: 5px;",
                            div(
                                title = "Select all available categories",
                                style = "cursor: help; display: inline-block;",
                                actionButton(ns("select_all_categories"), "Select All", size = "xs", style = "margin-right: 5px;")
                            ),
                            div(
                                title = "Deselect all categories",
                                style = "cursor: help; display: inline-block;",
                                actionButton(ns("clear_all_categories"), "Clear All", size = "xs")
                            )
                        )
                    )
                }
            })

            # Facet metadata selector
            output$facet_metadata_selector <- renderUI({
                req(input$facet_by == "metadata")
                metadata <- get_mc_data(dataset(), "metadata")
                req(metadata)

                # Get categorical metadata columns, excluding the x-axis variable if it's already metadata
                categorical_cols <- names(metadata)[sapply(metadata, function(x) !is.numeric(x) && !is.logical(x))]
                categorical_cols <- categorical_cols[categorical_cols != "metacell"]

                # Remove x-axis metadata variable if it exists
                if (input$x_axis_type == "metadata" && !is.null(input$x_axis_metadata_var)) {
                    categorical_cols <- categorical_cols[categorical_cols != input$x_axis_metadata_var]
                }

                if (length(categorical_cols) > 0) {
                    div(
                        title = "Choose a different categorical metadata variable to create separate panels. Each unique value will get its own panel in the plot.",
                        style = "cursor: help;",
                        selectInput(ns("facet_metadata_var"),
                            "Select metadata variable for faceting:",
                            choices = categorical_cols,
                            selected = categorical_cols[1]
                        )
                    )
                } else {
                    tags$div(style = "color: orange;", "No additional categorical metadata variables available for faceting")
                }
            })

            observe({
                req(input$boxplot_axis_type)
                req(input$boxplot_axis_var)
                metadata <- get_mc_data(dataset(), "metadata")
                req(metadata)

                shinyjs::toggle(id = "confusion_color_by_selector", condition = input$boxplot_axis_type == "Metadata" && input$boxplot_axis_var %in% colnames(metadata) && !is_numeric_field(metadata, input$boxplot_axis_var))
            })

            # Set default ylim values based on current data
            observe({
                req(input$boxplot_axis_type)
                req(input$boxplot_axis_var)
                req(dataset())
                req(metacell_types())
                # Handle cell type requirements based on x-axis type
                if (input$x_axis_type == "cell_types") {
                    req(input$boxplot_cell_types)
                    cell_types_to_use <- input$boxplot_cell_types
                } else {
                    req(input$cell_type_filter)
                    cell_types_to_use <- input$cell_type_filter
                }

                tryCatch(
                    {
                        if (input$boxplot_axis_type %in% c("Gene", "Gene module")) {
                            if (input$boxplot_axis_type == "Gene module") {
                                req(input$boxplot_axis_var %in% levels(gene_modules()$module))
                                genes <- get_module_genes(input$boxplot_axis_var, gene_modules())
                                egc_gene <- colSums(get_mc_egc(dataset(), genes = genes), na.rm = TRUE) + egc_epsilon
                            } else {
                                req(input$boxplot_axis_var %in% gene_names(dataset()))
                                egc_gene <- get_gene_egc(input$boxplot_axis_var, dataset()) + egc_epsilon
                            }

                            df <- metacell_types() %>% filter(cell_type %in% cell_types_to_use)
                            if (nrow(df) > 0) {
                                egc_subset <- egc_gene[df$metacell]
                                egc_subset <- egc_subset[!is.na(egc_subset)]
                                if (length(egc_subset) > 0) {
                                    y_min <- floor(min(egc_subset, na.rm = TRUE) * 10) / 10
                                    y_max <- ceiling(max(egc_subset, na.rm = TRUE) * 10) / 10
                                    updateNumericInput(session, "ylim_min", value = y_min)
                                    updateNumericInput(session, "ylim_max", value = y_max)
                                }
                            }
                        } else {
                            metadata <- get_mc_data(dataset(), "metadata")
                            req(!is.null(metadata))
                            req(input$boxplot_axis_var %in% colnames(metadata))
                            if (is_numeric_field(metadata, input$boxplot_axis_var)) {
                                df <- metacell_types() %>%
                                    filter(cell_type %in% cell_types_to_use) %>%
                                    left_join(metadata %>% select(metacell, !!input$boxplot_axis_var), by = "metacell")

                                if (nrow(df) > 0) {
                                    var_values <- df[[input$boxplot_axis_var]]
                                    var_values <- var_values[!is.na(var_values)]
                                    if (length(var_values) > 0) {
                                        y_min <- floor(min(var_values, na.rm = TRUE) * 10) / 10
                                        y_max <- ceiling(max(var_values, na.rm = TRUE) * 10) / 10
                                        updateNumericInput(session, "ylim_min", value = y_min)
                                        updateNumericInput(session, "ylim_max", value = y_max)
                                    }
                                }
                            }
                        }
                    },
                    error = function(e) {
                        # Silently ignore errors in default value setting
                    }
                )
            })

            # Validate ylim inputs and show warnings
            output$ylim_warning <- renderUI({
                if (isTRUE(input$custom_ylim)) {
                    ylim_min <- input$ylim_min
                    ylim_max <- input$ylim_max

                    if (is.null(ylim_min) || is.null(ylim_max)) {
                        return(tags$div(
                            style = "color: orange; font-size: 12px; margin-top: 5px;",
                            icon("exclamation-triangle"), " Please enter both min and max values"
                        ))
                    }

                    tryCatch(
                        {
                            ylim_min <- as.numeric(ylim_min)
                            ylim_max <- as.numeric(ylim_max)

                            if (is.na(ylim_min) || is.na(ylim_max)) {
                                return(tags$div(
                                    style = "color: red; font-size: 12px; margin-top: 5px;",
                                    icon("times-circle"), " Invalid numeric values"
                                ))
                            }

                            if (ylim_min >= ylim_max) {
                                return(tags$div(
                                    style = "color: red; font-size: 12px; margin-top: 5px;",
                                    icon("times-circle"), " Min must be less than max"
                                ))
                            }

                            return(tags$div(
                                style = "color: green; font-size: 12px; margin-top: 5px;",
                                icon("check-circle"), " Valid Y limits"
                            ))
                        },
                        error = function(e) {
                            return(tags$div(
                                style = "color: red; font-size: 12px; margin-top: 5px;",
                                icon("times-circle"), " Invalid input"
                            ))
                        }
                    )
                }
            })

            # Select All / Clear All observers for cell types
            observeEvent(input$select_all_cell_types, {
                req(metacell_types())
                cell_types <- unique(metacell_types()$cell_type)
                cell_types <- cell_types[!is.na(cell_types)]
                updateCheckboxGroupInput(session, "boxplot_cell_types", selected = as.character(cell_types))
            })

            observeEvent(input$clear_all_cell_types, {
                updateCheckboxGroupInput(session, "boxplot_cell_types", selected = character(0))
            })

            # Select All / Clear All observers for x-axis categories
            observeEvent(input$select_all_categories, {
                req(input$x_axis_metadata_var)
                metadata <- get_mc_data(dataset(), "metadata")
                req(metadata)
                req(input$x_axis_metadata_var %in% colnames(metadata))

                categories <- unique(metadata[[input$x_axis_metadata_var]])
                categories <- categories[!is.na(categories)]
                updateCheckboxGroupInput(session, "x_axis_categories", selected = as.character(categories))
            })

            observeEvent(input$clear_all_categories, {
                updateCheckboxGroupInput(session, "x_axis_categories", selected = character(0))
            })

            output$cell_type_boxplot <- plotly::renderPlotly({
                req(input$boxplot_axis_type)
                req(dataset())
                req(metacell_types())
                req(cell_type_colors())
                # Handle cell type requirements based on x-axis type
                if (input$x_axis_type == "cell_types") {
                    req(input$boxplot_cell_types)
                } else {
                    req(input$cell_type_filter)
                }
                req(input$boxplot_axis_var)
                req(input$plot_type)
                req(input$x_axis_type)
                req(input$facet_by)
                # coord_flip has a default value, so no req needed

                # Handle x-axis requirements
                if (input$x_axis_type == "metadata") {
                    req(input$x_axis_metadata_var)
                    req(input$x_axis_categories)
                }

                # Handle facet requirements
                if (input$facet_by == "metadata") {
                    req(input$facet_metadata_var)
                }

                if (input$boxplot_axis_type %in% c("Gene", "Gene module")) {
                    if (input$boxplot_axis_type == "Gene module") {
                        req(input$boxplot_axis_var %in% levels(gene_modules()$module))
                        genes <- get_module_genes(input$boxplot_axis_var, gene_modules())
                        egc_gene <- colSums(get_mc_egc(dataset(), genes = genes), na.rm = TRUE) + egc_epsilon
                    } else {
                        req(input$boxplot_axis_var %in% gene_names(dataset()))
                        egc_gene <- NULL
                    }

                    custom_ylim <- if (isTRUE(input$custom_ylim) && !is.null(input$ylim_min) && !is.null(input$ylim_max)) {
                        tryCatch(
                            {
                                ylim_min <- as.numeric(input$ylim_min)
                                ylim_max <- as.numeric(input$ylim_max)
                                if (!is.na(ylim_min) && !is.na(ylim_max) && ylim_min < ylim_max) {
                                    c(ylim_min, ylim_max)
                                } else {
                                    NULL
                                }
                            },
                            error = function(e) {
                                NULL
                            }
                        )
                    } else {
                        NULL
                    }

                    # Determine x-axis parameters
                    if (input$x_axis_type == "cell_types") {
                        x_axis_var <- "cell_type"
                        x_axis_categories <- input$boxplot_cell_types
                        cell_type_filter <- input$boxplot_cell_types
                    } else {
                        x_axis_var <- input$x_axis_metadata_var
                        x_axis_categories <- input$x_axis_categories
                        cell_type_filter <- input$cell_type_filter
                    }

                    # Determine facet parameters
                    if (input$facet_by == "metadata") {
                        facet_var <- input$facet_metadata_var
                    } else if (input$facet_by == "cell_type") {
                        facet_var <- "cell_type"
                    } else {
                        facet_var <- NULL
                    }

                    p <- cell_type_gene_boxplot(
                        input$boxplot_axis_var,
                        dataset(),
                        cell_types = cell_type_filter,
                        metacell_types = metacell_types(),
                        cell_type_colors = cell_type_colors(),
                        egc_gene = egc_gene,
                        plot_type = input$plot_type,
                        custom_ylim = custom_ylim,
                        x_axis_var = x_axis_var,
                        x_axis_categories = x_axis_categories,
                        facet_var = facet_var,
                        coord_flip = isTRUE(input$coord_flip),
                        log_scale = isTRUE(input$log_scale)
                    )
                } else {
                    metadata <- get_mc_data(dataset(), "metadata")
                    req(!is.null(metadata))
                    req(input$boxplot_axis_var %in% colnames(metadata))
                    if (is_numeric_field(metadata, input$boxplot_axis_var)) {
                        custom_ylim <- if (isTRUE(input$custom_ylim) && !is.null(input$ylim_min) && !is.null(input$ylim_max)) {
                            tryCatch(
                                {
                                    ylim_min <- as.numeric(input$ylim_min)
                                    ylim_max <- as.numeric(input$ylim_max)
                                    if (!is.na(ylim_min) && !is.na(ylim_max) && ylim_min < ylim_max) {
                                        c(ylim_min, ylim_max)
                                    } else {
                                        NULL
                                    }
                                },
                                error = function(e) {
                                    NULL
                                }
                            )
                        } else {
                            NULL
                        }

                        # Determine x-axis parameters
                        if (input$x_axis_type == "cell_types") {
                            x_axis_var <- "cell_type"
                            x_axis_categories <- input$boxplot_cell_types
                            cell_type_filter <- input$boxplot_cell_types
                        } else {
                            x_axis_var <- input$x_axis_metadata_var
                            x_axis_categories <- input$x_axis_categories
                            cell_type_filter <- input$cell_type_filter
                        }

                        # Determine facet parameters
                        if (input$facet_by == "metadata") {
                            facet_var <- input$facet_metadata_var
                        } else if (input$facet_by == "cell_type") {
                            facet_var <- "cell_type"
                        } else {
                            facet_var <- NULL
                        }

                        p <- cell_type_metadata_boxplot(
                            input$boxplot_axis_var,
                            dataset(),
                            cell_types = cell_type_filter,
                            metacell_types = metacell_types(),
                            cell_type_colors = cell_type_colors(),
                            plot_type = input$plot_type,
                            custom_ylim = custom_ylim,
                            x_axis_var = x_axis_var,
                            x_axis_categories = x_axis_categories,
                            facet_var = facet_var,
                            coord_flip = isTRUE(input$coord_flip),
                            log_scale = isTRUE(input$log_scale)
                        )
                    } else {
                        p <- cell_type_metadata_confusion(
                            input$boxplot_axis_var,
                            dataset(),
                            color_by = input$confusion_color_by,
                            cell_types = input$boxplot_cell_types,
                            metacell_types = metacell_types()
                        )
                    }
                }

                req(p)

                fig <- plotly::ggplotly(p, source = "cell_type_boxplot") %>%
                    sanitize_plotly_buttons() %>%
                    sanitize_plotly_download(globals)

                return(fig)
            })
        }
    )
}
