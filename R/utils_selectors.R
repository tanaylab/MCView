cell_type_selector <- function(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = NULL, cell_type_colors = NULL, metacell_types = NULL, apply_button = FALSE) {
    renderUI({
        if (is.null(cell_type_colors)) {
            colors <- NULL
        } else {
            colors <- cell_type_colors()
        }

        if (!is.null(metacell_types)) {
            colors <- colors %>% filter(cell_type %in% metacell_types()$cell_type)
        }

        cell_types_hex <- col2hex(get_cell_type_colors(dataset(), colors))
        cell_types <- names(get_cell_type_colors(dataset(), colors))

        if (!is.null(selected) && selected == "all") {
            selected <- cell_types
        }

        if ("(Missing)" %in% cell_types) {
            cell_types_hex <- c(cell_types_hex[cell_types != "(Missing)"], cell_types_hex[cell_types == "(Missing)"])
            cell_types <- c(cell_types[cell_types != "(Missing)"], "(Missing)")
        }

        picker <- shinyWidgets::pickerInput(ns(id), label,
            choices = cell_types,
            selected = selected,
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `dropup-auto` = FALSE),
            choicesOpt = list(
                style = paste0("color: ", cell_types_hex, ";")
            )
        )

        if (apply_button) {
            tagList(
                picker,
                shiny::actionButton(ns("apply_cell_types"), "Change cell types", size = "sm")
            )
        } else {
            picker
        }
    })
}

gene_modules_selector <- function(dataset, gene_modules, ns, id, label = "Gene modules", selected = NULL, with_genes = TRUE) {
    renderUI({
        modules <- levels(gene_modules()$module)
        if (with_genes) {
            modules <- modules[modules %in% gene_modules()$module]
        }
        if (!is.null(selected) && selected == "all") {
            selected <- modules
        }

        shinyWidgets::pickerInput(ns(id), label,
            choices = modules,
            selected = selected,
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `dropup-auto` = FALSE)
        )
    })
}

metacell_selector <- function(dataset, ns, id, label, selected = NULL, atlas = FALSE, ...) {
    renderUI({
        metacell_names <- colnames(get_mc_data(dataset(), "mc_mat", atlas = atlas))
        shinyWidgets::pickerInput(ns(id), label,
            choices = metacell_names,
            selected = selected %||% config$selected_mc1 %||% metacell_names[1],
            multiple = FALSE,
            options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE, ...)
        )
    })
}

colored_metacell_selector <- function(dataset, ns, id, label, metacell_colors, metacell_names, selected = NULL) {
    renderUI({
        req(metacell_colors())
        req(metacell_names())

        # add 'similar' annotation
        md <- get_mc_data(dataset(), "metadata")
        if (!is.null(md) && has_name(md, "similar")) {
            choices <- metacell_names()
            names(choices) <- ifelse(md$similar == "dissimilar", paste0(metacell_names(), " (dissimilar)"), metacell_names())
        } else {
            choices <- metacell_names()
        }

        cell_types_hex <- col2hex(metacell_colors())
        shinyWidgets::pickerInput(ns(id), label,
            choices = choices,
            selected = selected %||% config$selected_mc1 %||% metacell_names()[1],
            multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
            choicesOpt = list(
                style = paste0("color: ", cell_types_hex, ";")
            )
        )
    })
}

metadata_selector <- function(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata", selected = NULL, multiple = TRUE, additional_fields = c()) {
    renderUI({
        metadata <- get_mc_data(dataset(), metadata_id)
        req(metadata)
        metadata_fields <- c(additional_fields, colnames(metadata)[-1])
        shinyWidgets::virtualSelectInput(ns(id),
            label = label,
            choices = metadata_fields,
            selected = selected,
            multiple = multiple,
            showSelectedOptionsFirst = TRUE,
            search = TRUE,
            markSearchResults = TRUE
        )
    })
}


top_correlated_selector_multiple_genes <- function(input, output, session, dataset, ns, id, label, gene, action_id, action_label = "Add") {
    req(has_gg_mc_top_cor(project, dataset()))
    tagList(
        selectInput(
            ns(id),
            label,
            choices = c(get_top_cor_gene(dataset(), gene, type = "pos"), get_top_cor_gene(dataset(), gene, type = "neg")),
            selected = NULL,
            size = 10,
            selectize = FALSE,
            multiple = TRUE
        ),
        shinyWidgets::actionGroupButtons(ns(action_id), labels = action_label, size = "sm"),
        shiny::actionButton(
            inputId = ns(glue("genecards_{id}")), label = glue("GeneCards: {gene}"),
            size = "sm", onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}')")
        )
    )
}

top_correlated_selector <- function(gene_id, id, type_id, input, output, session, dataset, ns, button_labels = c("X", "Y", "Color", "2D"), ids = c("x", "y", "color", "proj2d"), gene_modules = NULL, metacell_types = NULL, selected_cell_types = NULL) {
    output[[glue("top_correlated_select_{id}")]] <- renderUI({
        if (!is.null(input$show_correlations) && !input$show_correlations) {
            if (input[[type_id]] == "Gene") {
                gene <- input[[gene_id]]
                req(gene %in% gene_names(dataset()))
                return(
                    shiny::actionButton(
                        inputId = ns(glue("genecards_{id}")), label = glue("GeneCards: {gene}"),
                        size = "sm", onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}')")
                    )
                )
            }
            req(input$show_correlations)
        }
        req(has_gg_mc_top_cor(project, dataset()))
        req(input[[type_id]] == "Gene" || input[[type_id]] == "Gene module" || input[[type_id]] == "Metadata")
        req(input[[gene_id]])

        metacell_filter <- NULL
        gene <- input[[gene_id]]
        data_vec <- NULL
        exclude <- NULL
        if (input[[type_id]] == "Gene") {
            req(gene %in% gene_names(dataset()))
            if (!is.null(selected_cell_types) &&
                length(selected_cell_types()) > 0 &&
                !is.null(metacell_types) &&
                length(setdiff(metacell_types()$cell_type, selected_cell_types())) > 0) {
                metacell_filter <- metacell_types()$metacell[metacell_types()$cell_type %in% selected_cell_types()]
            }
        } else if (input[[type_id]] == "Gene module") {
            req(!is.null(gene_modules))
            req(gene %in% gene_modules()$module)
            data_vec <- get_gene_module_egc(gene, dataset(), gene_modules()) + egc_epsilon
            exclude <- gene_modules()$gene[gene_modules()$module == gene]
        } else {
            metadata <- get_mc_data(dataset(), "metadata")
            req(gene %in% colnames(metadata))
            req(is_numeric_field(metadata, gene))
            data_vec <- metadata[[gene]]
            names(data_vec) <- metadata$metacell
        }

        input_ids <- purrr::map_chr(
            ids,
            ~ {
                ns(glue("select_top_cor_{id}_{.x}"))
            }
        )

        si <- selectInput(
            ns(glue("selected_top_{id}")),
            glue("Top correlated to {gene}:"),
            choices = get_top_cor_gene(dataset(), gene, type = "both", data_vec = data_vec, exclude = exclude, metacell_filter = metacell_filter),
            selected = NULL,
            size = 10,
            selectize = FALSE
        )
        btn <- shinyWidgets::actionGroupButtons(
            input_ids[1:length(button_labels)],
            labels = button_labels, size = "sm", fullwidth = FALSE
        )

        if (input[[type_id]] == "Gene") {
            tagList(
                si,
                btn,
                shiny::actionButton(
                    inputId = ns(glue("genecards_{id}")), label = glue("GeneCards: {gene}"),
                    size = "sm", onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}')")
                )
            )
        } else {
            tagList(
                si,
                btn
            )
        }
    })

    observeEvent(input[[glue("select_top_cor_{id}_axis")]], {
        req(input[["axis_type"]] == "Gene")
        shinyWidgets::updateVirtualSelect(session = session, inputId = "axis_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_x")]], {
        req(input[["x_axis_type"]] == "Gene")
        shinyWidgets::updateVirtualSelect(session = session, inputId = "x_axis_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_y")]], {
        req(input[["y_axis_type"]] == "Gene")
        shinyWidgets::updateVirtualSelect(session = session, inputId = "y_axis_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_color")]], {
        req(input[["color_by_type"]] == "Gene")
        shinyWidgets::updateVirtualSelect(session = session, inputId = "color_by_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_proj2d")]], {
        req(input[["color_proj"]] == "Gene")
        shinyWidgets::updateVirtualSelect(session = session, inputId = "color_proj_gene", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_boxplot")]], {
        req(input[["boxplot_axis_type"]] == "Gene")
        shinyWidgets::updateVirtualSelect(session = session, inputId = "boxplot_axis_var", selected = input[[glue("selected_top_{id}")]])
    })
}

top_correlated_selectors <- function(input, output, session, dataset, metacell_types, ns, button_labels = c("X", "Y", "Color", "2D"), gene_modules = NULL, selected_cell_types = NULL) {
    top_correlated_selector("x_axis_var", "x_axis", "x_axis_type", input, output, session, dataset, ns, button_labels = button_labels, gene_modules = gene_modules, metacell_types = metacell_types, selected_cell_types = selected_cell_types)
    top_correlated_selector("y_axis_var", "y_axis", "y_axis_type", input, output, session, dataset, ns, button_labels = button_labels, gene_modules = gene_modules, metacell_types = metacell_types, selected_cell_types = selected_cell_types)
    top_correlated_selector("color_by_var", "color_by", "color_by_type", input, output, session, dataset, ns, button_labels = button_labels, gene_modules = gene_modules, metacell_types = metacell_types, selected_cell_types = selected_cell_types)
    top_correlated_selector("color_proj_gene", "color_proj", "color_proj", input, output, session, dataset, ns, button_labels = button_labels, gene_modules = gene_modules, metacell_types = metacell_types, selected_cell_types = selected_cell_types)
}
