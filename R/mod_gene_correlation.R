#' gene_correlation UI Function
#'
#' @description A shiny Module for gene correlations.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @importFrom rclipboard rclipboardSetup rclipButton
mod_gene_correlation_ui <- function(id) {
    ns <- NS(id)

    tagList(
        # Setup clipboard functionality
        if (requireNamespace("rclipboard", quietly = TRUE)) {
            rclipboard::rclipboardSetup()
        } else {
            tags$head(tags$script(
                "function copyToClipboard(text) {
                    if (navigator.clipboard) {
                        navigator.clipboard.writeText(text).then(function() {
                            console.log('Copied to clipboard');
                        }).catch(function(err) {
                            console.error('Failed to copy: ', err);
                        });
                    } else {
                        // Fallback for older browsers
                        var textArea = document.createElement('textarea');
                        textArea.value = text;
                        document.body.appendChild(textArea);
                        textArea.select();
                        document.execCommand('copy');
                        document.body.removeChild(textArea);
                    }
                }"
            ))
        },
        fluidRow(
            # Left panel: Input and controls
            generic_column(
                width = 4,
                generic_box(
                    id = ns("gene_input_box"),
                    title = "Gene Input",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,

                    # Gene list input methods
                    textAreaInput(ns("gene_list"),
                        "Gene List (one per line or comma-separated)",
                        rows = 8,
                        value = "",
                        placeholder = "Enter gene names here...\nExample:\nACTB\nGAPDH\nTP53"
                    ),
                    fluidRow(
                        generic_column(
                            width = 12,
                            actionButton(ns("clear_genes"), "Clear",
                                class = "btn-outline-secondary", style = "width: 100%;"
                            )
                        )
                    ),
                    br(),
                    fileInput(ns("gene_file"), "Or upload file",
                        accept = c(".csv", ".txt", ".tsv"),
                        buttonLabel = "Browse...",
                        placeholder = "No file selected"
                    ),
                    tags$hr(),

                    # Mode selection
                    shinyWidgets::radioGroupButtons(
                        inputId = ns("correlation_mode"),
                        label = "Correlation Mode:",
                        choices = list(
                            "Individual genes" = "individual",
                            "Gene module/anchor" = "module"
                        ),
                        selected = "individual",
                        justified = TRUE
                    ),

                    # Toggle for gene search vs correlation calculation
                    conditionalPanel(
                        condition = paste0("input['", ns("correlation_mode"), "'] == 'individual'"),
                        shinyWidgets::radioGroupButtons(
                            inputId = ns("analysis_type"),
                            label = "Analysis Type:",
                            choices = list(
                                "Find correlated genes" = "find_genes",
                                "Calculate gene-gene correlation" = "gene_gene_cor"
                            ),
                            selected = "find_genes",
                            justified = TRUE
                        )
                    ),

                    # Help text for modes
                    conditionalPanel(
                        condition = paste0("input['", ns("correlation_mode"), "'] == 'individual' && input['", ns("analysis_type"), "'] == 'find_genes'"),
                        div(
                            class = "alert alert-info", style = "font-size: 12px;",
                            icon("info-circle"), " Find genes correlated with each input gene separately"
                        )
                    ),
                    conditionalPanel(
                        condition = paste0("input['", ns("correlation_mode"), "'] == 'individual' && input['", ns("analysis_type"), "'] == 'gene_gene_cor'"),
                        div(
                            class = "alert alert-info", style = "font-size: 12px;",
                            icon("info-circle"), " Calculate correlations between the input genes (all correlations shown, filters disabled)"
                        )
                    ),
                    conditionalPanel(
                        condition = paste0("input['", ns("correlation_mode"), "'] == 'module'"),
                        div(
                            class = "alert alert-info", style = "font-size: 12px;",
                            icon("info-circle"), " Treat genes as a module and find genes correlated with the combined expression"
                        )
                    ),
                    tags$hr(),

                    # Parameters
                    conditionalPanel(
                        condition = paste0("!(input['", ns("correlation_mode"), "'] == 'individual' && input['", ns("analysis_type"), "'] == 'gene_gene_cor')"),
                        numericInput(ns("n_correlations"), "Top correlations per gene",
                            value = 30, min = 5, max = 100, step = 5
                        ),
                        sliderInput(ns("cor_threshold"), "Correlation threshold",
                            min = 0, max = 1, value = 0, step = 0.05
                        ),

                        # Correlation direction filter
                        shinyWidgets::radioGroupButtons(
                            inputId = ns("correlation_direction"),
                            label = "Correlation Direction:",
                            choices = list(
                                "Positive only" = "positive",
                                "Negative only" = "negative",
                                "Both" = "both"
                            ),
                            selected = "positive",
                            justified = TRUE
                        )
                    ),

                    # Cell type filtering
                    uiOutput(ns("cell_type_filter_ui")),
                    br(),

                    # Calculate button
                    actionButton(ns("calculate_correlations"), "Calculate Correlations",
                        class = "btn-primary btn-lg",
                        style = "width: 100%; font-weight: bold;",
                        icon = icon("calculator")
                    ),
                    br(), br(),

                    # Status display
                    uiOutput(ns("status_display"))
                )
            ),

            # Right panel: Visualizations
            generic_column(
                width = 8,
                # Conditional visualization based on input
                conditionalPanel(
                    condition = paste0("output['", ns("show_heatmap"), "']"),
                    generic_box(
                        title = "Correlation Heatmap",
                        status = "success",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        closable = FALSE,
                        width = 12,
                        sidebar = shinydashboardPlus::boxSidebar(
                            startOpen = FALSE,
                            width = 50,
                            id = ns("heatmap_sidebar"),
                            checkboxInput(ns("cluster_heatmap"), "Cluster genes", value = TRUE),
                            checkboxInput(ns("mask_low_correlations"), "Mask low correlations", value = FALSE),
                            numericInput(ns("heatmap_top_genes"), "Max genes to show",
                                value = 50, min = 10, max = 200, step = 10
                            )
                        ),
                        shinycssloaders::withSpinner(
                            plotOutput(ns("correlation_heatmap"), height = "500px")
                        )
                    )
                ),
                conditionalPanel(
                    condition = paste0("output['", ns("show_barplot"), "']"),
                    generic_box(
                        title = "Top Correlations",
                        status = "success",
                        solidHeader = TRUE,
                        collapsible = TRUE,
                        closable = FALSE,
                        width = 12,
                        shinycssloaders::withSpinner(
                            plotly::plotlyOutput(ns("correlation_barplot"), height = "500px")
                        )
                    )
                ),

                # Instructions when no results
                conditionalPanel(
                    condition = paste0("!output['", ns("show_heatmap"), "'] && !output['", ns("show_barplot"), "']"),
                    generic_box(
                        title = "Instructions",
                        status = "info",
                        solidHeader = TRUE,
                        width = 12,
                        div(
                            class = "text-center", style = "padding: 40px;",
                            icon("lightbulb", class = "fa-3x text-muted"), br(), br(),
                            h4("How to use:", class = "text-muted"),
                            p("1. Enter gene names in the text area (one per line)", class = "text-left"),
                            p("2. Choose correlation mode (individual vs module)", class = "text-left"),
                            p("3. Adjust parameters as needed", class = "text-left"),
                            p("4. Click 'Calculate Correlations' to generate results", class = "text-left"),
                            br(),
                            p("Visualizations will appear here after calculation.", class = "text-muted")
                        )
                    )
                )
            )
        ),

        # Results table row
        fluidRow(
            generic_column(
                width = 12,
                generic_box(
                    id = ns("results_box"),
                    title = "Correlation Results",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    conditionalPanel(
                        condition = paste0("output['", ns("show_results_table"), "']"),
                        shinycssloaders::withSpinner(
                            DT::DTOutput(ns("correlation_table"))
                        ),
                        br(),
                        fluidRow(
                            generic_column(
                                width = 6,
                                downloadButton(ns("download_correlations"), "Export Results",
                                    class = "btn-success",
                                    icon = icon("download"),
                                    style = "width: 100%;"
                                )
                            ),
                            generic_column(
                                width = 6,
                                uiOutput(ns("copy_genes_button"))
                            )
                        )
                    ),
                    conditionalPanel(
                        condition = paste0("!output['", ns("show_results_table"), "']"),
                        div(
                            class = "text-center text-muted", style = "padding: 20px;",
                            "Results will appear here after calculation."
                        )
                    )
                )
            )
        )
    )
}

#' gene_correlation Server Function
#'
#' @noRd
mod_gene_correlation_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            # Reactive values
            parsed_genes <- reactiveVal(character(0))
            correlation_results <- reactiveVal(NULL)
            calculation_status <- reactiveVal("")

            # Parse genes from input
            observe({
                genes <- c()

                # From text area
                if (!is.null(input$gene_list) && nzchar(input$gene_list)) {
                    text_genes <- strsplit(input$gene_list, "[\n\r,;\\s]+")[[1]]
                    genes <- c(genes, text_genes)
                }

                # Clean and deduplicate
                genes <- unique(trimws(genes[nzchar(genes)]))
                parsed_genes(genes)
            })

            # Handle file upload
            observeEvent(input$gene_file, {
                req(input$gene_file)

                tryCatch(
                    {
                        ext <- tools::file_ext(input$gene_file$name)

                        if (ext %in% c("csv", "tsv")) {
                            sep <- ifelse(ext == "csv", ",", "\t")
                            file_data <- read.table(input$gene_file$datapath,
                                sep = sep,
                                header = FALSE, stringsAsFactors = FALSE
                            )
                            file_genes <- as.character(unlist(file_data))
                        } else {
                            file_genes <- readLines(input$gene_file$datapath)
                        }

                        # Combine with existing genes
                        current_genes <- if (nzchar(input$gene_list)) input$gene_list else ""
                        new_genes <- paste(file_genes, collapse = "\n")
                        combined <- if (nzchar(current_genes)) {
                            paste(current_genes, new_genes, sep = "\n")
                        } else {
                            new_genes
                        }

                        updateTextAreaInput(session, "gene_list", value = combined)
                    },
                    error = function(e) {
                        showNotification(paste("Error reading file:", e$message), type = "error")
                    }
                )
            })


            # Clear genes
            observeEvent(input$clear_genes, {
                updateTextAreaInput(session, "gene_list", value = "")
            })

            # Cell type filter UI
            output$cell_type_filter_ui <- renderUI({
                req(cell_type_colors())
                req(metacell_types())

                cell_types <- unique(cell_type_colors()$cell_type)
                cell_types <- cell_types[cell_types %in% metacell_types()$cell_type]

                if (length(cell_types) > 1) {
                    cell_types_hex <- col2hex(cell_type_colors()$color)
                    names(cell_types_hex) <- cell_type_colors()$cell_type

                    tagList(
                        checkboxInput(ns("enable_cell_type_filter"),
                            "Filter by cell types",
                            value = FALSE
                        ),
                        conditionalPanel(
                            condition = paste0("input['", ns("enable_cell_type_filter"), "']"),
                            shinyWidgets::pickerInput(
                                ns("selected_cell_types"),
                                "Cell types:",
                                choices = cell_types,
                                selected = cell_types,
                                multiple = TRUE,
                                options = list(
                                    `actions-box` = TRUE,
                                    `dropup-auto` = FALSE,
                                    `live-search` = TRUE
                                ),
                                choicesOpt = list(
                                    style = paste0("color: ", cell_types_hex[cell_types], ";")
                                )
                            )
                        )
                    )
                } else {
                    div()
                }
            })

            # Main calculation triggered by button
            observeEvent(input$calculate_correlations, {
                req(length(parsed_genes()) > 0)

                # Validate genes
                available_genes <- gene_names(dataset())
                valid_genes <- parsed_genes()[parsed_genes() %in% available_genes]
                invalid_genes <- setdiff(parsed_genes(), valid_genes)

                if (length(invalid_genes) > 0) {
                    showNotification(
                        paste(
                            "Invalid genes (ignored):", paste(head(invalid_genes, 5), collapse = ", "),
                            if (length(invalid_genes) > 5) paste("and", length(invalid_genes) - 5, "more")
                        ),
                        type = "warning",
                        duration = 5
                    )
                }

                if (length(valid_genes) == 0) {
                    showNotification("No valid genes found", type = "error")
                    return()
                }

                # Apply gene limits based on mode
                if (input$correlation_mode == "individual") {
                    if (!is.null(input$analysis_type) && input$analysis_type == "gene_gene_cor") {
                        # Limit to 250 genes for gene-gene correlation
                        if (length(valid_genes) > 250) {
                            valid_genes <- valid_genes[1:250]
                            showNotification(
                                paste("Limited to first 250 genes for gene-gene correlation analysis. Using:", paste(head(valid_genes, 5), collapse = ", "), "..."),
                                type = "warning",
                                duration = 5
                            )
                        }
                    } else {
                        # Limit to 300 genes for finding neighbors
                        if (length(valid_genes) > 300) {
                            valid_genes <- valid_genes[1:300]
                            showNotification(
                                paste("Limited to first 300 genes for correlation analysis. Using:", paste(head(valid_genes, 5), collapse = ", "), "..."),
                                type = "warning",
                                duration = 5
                            )
                        }
                    }
                } else {
                    # Module mode - no additional limit beyond existing logic
                }

                # Show progress
                calculation_status("Calculating correlations...")

                # Get cell type filter
                cell_type_filter <- NULL
                if (!is.null(input$enable_cell_type_filter) && input$enable_cell_type_filter &&
                    !is.null(input$selected_cell_types) && length(input$selected_cell_types) > 0) {
                    all_cell_types <- unique(metacell_types()$cell_type)
                    if (length(input$selected_cell_types) < length(all_cell_types)) {
                        cell_type_filter <- input$selected_cell_types
                    }
                }

                withProgress(message = "Calculating correlations...", value = 0, {
                    tryCatch(
                        {
                            if (input$correlation_mode == "individual") {
                                if (!is.null(input$analysis_type) && input$analysis_type == "gene_gene_cor") {
                                    incProgress(0.2, detail = "Gene-gene correlations")
                                    # For gene-gene correlations, don't apply threshold - show all correlations between input genes
                                    results <- calc_gene_gene_correlations(dataset(), valid_genes, cell_type_filter, threshold = 0)
                                } else {
                                    incProgress(0.2, detail = "Individual gene correlations")
                                    results <- calc_individual_correlations(dataset(), valid_genes, input$n_correlations, cell_type_filter, input$cor_threshold)
                                }
                            } else {
                                incProgress(0.2, detail = "Module correlation")
                                results <- calc_module_correlations(dataset(), valid_genes, input$n_correlations, cell_type_filter, input$cor_threshold)
                            }

                            incProgress(0.6, detail = "Processing results")

                            # Add additional columns and filter by direction
                            results <- results %>%
                                mutate(
                                    mode = input$correlation_mode,
                                    correlation_type = ifelse(cor > 0, "positive", "negative"),
                                    abs_cor = abs(cor)
                                )

                            # Filter by correlation direction (except for gene-gene correlations)
                            if (input$correlation_mode == "individual" &&
                                !is.null(input$analysis_type) && input$analysis_type == "gene_gene_cor") {
                                # For gene-gene correlations, don't apply direction filter
                            } else {
                                if (input$correlation_direction == "positive") {
                                    results <- results %>% filter(cor > 0)
                                } else if (input$correlation_direction == "negative") {
                                    results <- results %>% filter(cor < 0)
                                }
                            }

                            results <- results %>% arrange(input_gene, desc(abs_cor))

                            incProgress(0.2, detail = "Finalizing")

                            correlation_results(results)
                            calculation_status(paste("Calculated correlations for", length(valid_genes), "genes"))

                            showNotification(
                                paste("Successfully calculated correlations for", length(valid_genes), "genes"),
                                type = "default"
                            )
                        },
                        error = function(e) {
                            calculation_status(paste("Error:", e$message))
                            showNotification(paste("Calculation failed:", e$message), type = "error")
                            correlation_results(NULL)
                        }
                    )
                })
            })

            # Status display
            output$status_display <- renderUI({
                if (nzchar(calculation_status())) {
                    if (startsWith(calculation_status(), "Error:")) {
                        div(
                            class = "alert alert-danger",
                            icon("exclamation-triangle"), " ", calculation_status()
                        )
                    } else if (startsWith(calculation_status(), "Calculated")) {
                        div(
                            class = "alert alert-success",
                            icon("check"), " ", calculation_status()
                        )
                    } else {
                        div(
                            class = "alert alert-info",
                            icon("info"), " ", calculation_status()
                        )
                    }
                }
            })

            # Control visualization display
            output$show_heatmap <- reactive({
                !is.null(correlation_results()) && length(unique(correlation_results()$input_gene)) > 1
            })
            outputOptions(output, "show_heatmap", suspendWhenHidden = FALSE)

            output$show_barplot <- reactive({
                !is.null(correlation_results()) && length(unique(correlation_results()$input_gene)) == 1
            })
            outputOptions(output, "show_barplot", suspendWhenHidden = FALSE)

            output$show_results_table <- reactive({
                !is.null(correlation_results())
            })
            outputOptions(output, "show_results_table", suspendWhenHidden = FALSE)

            # Heatmap visualization with caching
            output$correlation_heatmap <- renderPlot({
                req(correlation_results())
                req(length(unique(correlation_results()$input_gene)) > 1)

                input_genes <- unique(correlation_results()$input_gene)
                plot_correlation_heatmap(
                    correlation_results(),
                    input_genes,
                    dataset(),
                    cluster = input$cluster_heatmap,
                    max_genes = input$heatmap_top_genes,
                    mask_low_correlations = input$mask_low_correlations,
                    correlation_mode = input$correlation_mode
                )
            }) %>% bindCache(correlation_results(), dataset(), input$heatmap_top_genes, input$cluster_heatmap, input$mask_low_correlations, input$correlation_mode)

            # Barplot visualization
            output$correlation_barplot <- plotly::renderPlotly({
                req(correlation_results())
                req(length(unique(correlation_results()$input_gene)) == 1)

                input_gene <- unique(correlation_results()$input_gene)[1]
                plot_correlation_barplot(
                    correlation_results(),
                    gene = input_gene
                )
            }) %>% bindCache(correlation_results())

            # Results table
            output$correlation_table <- DT::renderDataTable({
                req(correlation_results())

                # Filter by direction (except for gene-gene correlations)
                filtered_data <- correlation_results()

                # Apply direction filter (except for gene-gene correlations)
                if (input$correlation_mode == "individual" &&
                    !is.null(input$analysis_type) && input$analysis_type == "gene_gene_cor") {
                    # For gene-gene correlations, don't apply direction filter
                } else {
                    if (input$correlation_direction == "positive") {
                        filtered_data <- filtered_data %>% filter(cor > 0)
                    } else if (input$correlation_direction == "negative") {
                        filtered_data <- filtered_data %>% filter(cor < 0)
                    }
                }

                filtered_data <- filtered_data %>%
                    select(
                        `Input Gene` = input_gene,
                        `Correlated Gene` = gene2,
                        `Correlation` = cor,
                        `Type` = correlation_type,
                        `Rank` = rank
                    ) %>%
                    arrange(`Input Gene`, desc(abs(`Correlation`)))

                DT::datatable(
                    filtered_data,
                    escape = FALSE,
                    rownames = FALSE,
                    options = list(
                        dom = "Bfrtip",
                        pageLength = 25,
                        scrollX = TRUE,
                        columnDefs = list(
                            list(className = "dt-center", targets = c(2, 3, 4))
                        )
                    )
                ) %>%
                    DT::formatRound("Correlation", digits = 3) %>%
                    DT::formatStyle(
                        "Correlation",
                        backgroundColor = DT::styleInterval(
                            c(-0.7, -0.3, 0, 0.3, 0.7),
                            c("#d73027", "#fc8d59", "#fee08b", "#e0f3f8", "#91bfdb", "#4575b4")
                        )
                    )
            })

            # Reactive for genes to copy
            copy_genes_data <- reactive({
                req(correlation_results())

                all_genes <- correlation_results()

                # Apply direction filter (except for gene-gene correlations)
                if (input$correlation_mode == "individual" &&
                    !is.null(input$analysis_type) && input$analysis_type == "gene_gene_cor") {
                    # For gene-gene correlations, don't apply direction filter
                } else {
                    if (input$correlation_direction == "positive") {
                        all_genes <- all_genes %>% filter(cor > 0)
                    } else if (input$correlation_direction == "negative") {
                        all_genes <- all_genes %>% filter(cor < 0)
                    }
                }

                unique(all_genes$gene2)
            })

            # Render copy button using utility function
            output$copy_genes_button <- clipboard_copy_button_ui(
                ns, "copy_all_genes", copy_genes_data,
                label = "Copy All Genes",
                style = "width: 100%; background-color: #17a2b8; color: white; border: none;",
                tooltip = "Copy all found genes to clipboard"
            )

            # Copy all genes functionality using utility function
            clipboard_copy_button_server(
                input, "copy_all_genes", copy_genes_data, globals,
                message_template = "Copied {count} genes to clipboard"
            )

            # Export functionality
            output$download_correlations <- downloadHandler(
                filename = function() {
                    mode_suffix <- ifelse(input$correlation_mode == "individual", "individual", "module")
                    paste0("gene_correlations_", mode_suffix, "_", Sys.Date(), ".csv")
                },
                content = function(file) {
                    req(correlation_results())

                    data <- correlation_results()

                    # Apply direction filter (except for gene-gene correlations)
                    if (input$correlation_mode == "individual" &&
                        !is.null(input$analysis_type) && input$analysis_type == "gene_gene_cor") {
                        # For gene-gene correlations, don't apply direction filter
                    } else {
                        if (input$correlation_direction == "positive") {
                            data <- data %>% filter(cor > 0)
                        } else if (input$correlation_direction == "negative") {
                            data <- data %>% filter(cor < 0)
                        }
                    }

                    data <- data %>%
                        select(
                            input_gene,
                            correlated_gene = gene2,
                            correlation = cor,
                            correlation_type,
                            rank,
                            mode
                        ) %>%
                        arrange(input_gene, desc(abs(correlation)))

                    # Add metadata header
                    writeLines(c(
                        "# MCView Gene Correlation Export",
                        paste("# Date:", Sys.time()),
                        paste("# Mode:", input$correlation_mode),
                        paste("# Analysis type:", if (!is.null(input$analysis_type)) input$analysis_type else "find_genes"),
                        paste("# Direction:", input$correlation_direction),
                        paste("# Threshold:", input$cor_threshold),
                        paste("# Input genes:", paste(parsed_genes(), collapse = ", ")),
                        paste("# Total correlations:", nrow(data)),
                        ""
                    ), file)

                    # Append data
                    write.csv(data, file, row.names = FALSE, append = TRUE)
                }
            )
        }
    )
}
