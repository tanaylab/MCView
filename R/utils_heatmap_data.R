# utils_heatmap_data.R - Heatmap matrix reactives (data layer)
#
# Split from R/utils_heatmap.R (2026-05-01).  Hosts the
# heatmap_matrix_reactives() function that builds the heatmap matrix and
# its filter / categorical-controls reactives. Companions:
#   - R/utils_heatmap_server.R   - heatmap_reactives() (module body + plot)
#   - R/utils_heatmap_handlers.R - tooltip + download handlers
#   - R/utils_heatmap_ui.R       - UI builders (existing, unchanged)
#   - R/utils_heatmap_help.R     - help-modal content (existing, unchanged)

heatmap_matrix_reactives <- function(ns, input, output, session, dataset, metacell_types, cell_type_colors, state, markers, lfp_range, mode, metacell_filter, mat) {
    observe({
        choices <- markers()
        if (!is.null(choices)) {
            names(choices) <- gene_label(choices, dataset())
            updateSelectInput(session, "selected_marker_genes", choices = choices)
        }
        shinyWidgets::updateVirtualSelect(session = session, inputId = "metadata_order_var", choices = c("Hierarchical-Clustering", dataset_metadata_fields_numeric(dataset())), selected = "Hierarchical-Clustering")
        shinyWidgets::updateVirtualSelect(session = session, inputId = "metadata_order_cell_type_var", choices = c("Hierarchical-Clustering", "Colors table", dataset_metadata_fields_numeric(dataset())), selected = "Hierarchical-Clustering")

        shinyjs::toggle(id = "metadata_order_cell_type_var", condition = input$force_cell_type)
        shinyjs::toggle(id = "categorical_filter_ui", condition = input$enable_categorical_filter)
    })

    output$categorical_filter_ui <- renderUI({
        categorical_fields <- dataset_metadata_fields_categorical(dataset())
        if (length(categorical_fields) == 0) {
            return(p("No filterable categorical metadata fields available",
                style = "color: #6c757d; font-style: italic;"
            ))
        }

        tagList(
            shinyWidgets::virtualSelectInput(
                ns("categorical_filter_field"),
                "Filter by categorical variable:",
                choices = categorical_fields,
                selected = categorical_fields[1],
                multiple = FALSE,
                search = TRUE
            ),
            uiOutput(ns("categorical_filter_values_ui"))
        )
    })

    output$categorical_filter_values_ui <- renderUI({
        req(input$categorical_filter_field)

        metadata <- get_mc_data(dataset(), "metadata")
        req(metadata)
        req(input$categorical_filter_field %in% colnames(metadata))

        unique_values <- unique(metadata[[input$categorical_filter_field]])
        # Remove NA and empty strings
        unique_values <- unique_values[!is.na(unique_values) & unique_values != ""]
        unique_values <- sort(unique_values)

        shinyWidgets::virtualSelectInput(
            ns("categorical_filter_values"),
            "Select values to include:",
            choices = unique_values,
            selected = unique_values,
            multiple = TRUE,
            search = TRUE
        )
    })

    observe({
        if (is.null(markers())) {
            config <- mcv_get("config")
            if (config$light_version) {
                max_gene_num <- 100
            } else {
                req(input$max_gene_num)
                max_gene_num <- input$max_gene_num
            }

            initial_markers <- choose_markers(get_marker_genes(dataset(), mode = mode), max_gene_num, dataset = dataset())
            markers(initial_markers)
        }

        lfp_range(input$lfp_range)
    })

    observeEvent(input$load_genes, {
        req(input$load_genes)
        req(input$load_genes$datapath)
        req(input$load_genes$datapath != "")
        new_markers <- data.table::fread(input$load_genes$datapath, header = FALSE, data.table = FALSE)[, 1]
        new_markers <- new_markers[new_markers %in% gene_names(dataset())]
        if (length(new_markers) == 0) {
            showNotification("No valid genes were loaded", type = "warning")
        } else {
            markers(new_markers)
            showNotification(glue("Loaded {length(new_markers)} genes"), type = "message")
        }
    })

    observeEvent(input$use_de_genes, {
        if (!is.null(state$selection$significant_genes) && length(state$selection$significant_genes) > 0) {
            markers(state$selection$significant_genes)
            showNotification(glue("Updated markers with {length(state$selection$significant_genes)} significant genes from current Diff. Expression comparison"), type = "message")
        } else {
            showNotification("No significant genes available in the current Diff. Expression comparison.", type = "warning")
        }
    })

    observeEvent(input$update_genes, {
        req(metacell_types())
        req(cell_type_colors())
        req(is.null(input$selected_cell_types) || all(input$selected_cell_types %in% c(cell_type_colors()$cell_type, "(Missing)")))

        if (input$mat_value == "Local") {
            mc_egc <- get_mc_egc(dataset(), metacells = colnames(mat()))
            markers_df <- calc_marker_genes(mc_egc, genes_per_metacell = 20, daf_obj = get_dataset_daf(dataset()))
        } else {
            if (!is.null(input$selected_cell_types)) {
                markers_df <- metacell_types() %>%
                    filter(cell_type %in% input$selected_cell_types)
            } else {
                markers_df <- metacell_types()
            }

            if (!is.null(input$use_markers) && input$use_markers) {
                new_markers_df <- get_marker_genes(dataset(), mode = "Markers")
            } else {
                new_markers_df <- get_marker_genes(dataset(), mode = mode)
            }

            if (has_name(new_markers_df, "metacell")) {
                markers_df <- markers_df %>%
                    select(metacell) %>%
                    inner_join(new_markers_df, by = "metacell")
            } else {
                markers_df <- new_markers_df
            }

            if (!is.null(input$filter_by_clipboard) && input$filter_by_clipboard) {
                if (!is.null(state$selection$clipboard) && length(state$selection$clipboard) > 0) {
                    markers_df <- markers_df %>%
                        filter(metacell %in% state$selection$clipboard)
                }
            }
        }

        req(input$max_gene_num)
        markers_df <- filter_genes_by_flags(
            markers_df,
            lateral_genes = get_mc_data(dataset(), "lateral_genes"),
            noisy_genes = get_mc_data(dataset(), "noisy_genes"),
            include_lateral = is.null(input$include_lateral) || input$include_lateral,
            include_noisy = is.null(input$include_noisy) || input$include_noisy
        )
        if (!is.null(input$show_only_fitted) && input$show_only_fitted) {
            gene_metadata <- get_mc_data(dataset(), "gene_metadata")

            if (!is.null(gene_metadata)) {
                fitted_genes <- gene_metadata %>%
                    filter(fitted_gene) %>%
                    pull(gene) %>%
                    unique()
                markers_df <- markers_df %>%
                    filter(gene %in% fitted_genes)
            }
        }

        new_markers <- choose_markers(markers_df, input$max_gene_num, dataset = dataset())

        # If we did not choose all the cell types
        if (!is.null(input$selected_cell_types) && length(input$selected_cell_types) != nrow(cell_type_colors())) {
            # TODO: add also genes that are distictive for the specific cell type (although not variable within it)
        }

        markers(new_markers)
    })

    observeEvent(input$remove_genes, {
        req(markers())
        new_markers <- markers()[!(markers() %in% input$selected_marker_genes)]
        markers(new_markers)
        shinyWidgets::updateVirtualSelect(session = session, inputId = "genes_to_add", selected = c())
    })

    output$add_genes_ui <- renderUI({
        gene_choices <- gene_names(dataset())


        tagList(
            shinyWidgets::virtualSelectInput(ns("genes_to_add"),
                label = "Add genes",
                choices = gene_choices,
                selected = c(),
                multiple = TRUE,
                showSelectedOptionsFirst = TRUE,
                search = TRUE,
                markSearchResults = TRUE,
                searchByStartsWith = TRUE,
                disableSelectAll = TRUE
            ),
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm")
        )
    })

    observeEvent(input$add_genes, {
        new_markers <- sort(unique(c(markers(), input$genes_to_add)))
        markers(new_markers)
        shinyWidgets::updateVirtualSelect(session = session, inputId = "genes_to_add", selected = character(0))
    })
}


