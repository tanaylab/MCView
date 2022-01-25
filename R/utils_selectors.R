cell_type_selector <- function(dataset, ns, id = "selected_cell_types", label = "Cell types", selected = NULL, cell_type_colors = NULL) {
    renderUI({
        if (is.null(cell_type_colors)) {
            colors <- NULL
        } else {
            colors <- cell_type_colors()
        }
        cell_types_hex <- col2hex(get_cell_type_colors(dataset(), colors))
        cell_types <- names(get_cell_type_colors(dataset(), colors))

        if (!is.null(selected) && selected == "all") {
            selected <- cell_types
        }
        shinyWidgets::pickerInput(ns(id), label,
            choices = cell_types,
            selected = selected,
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `dropup-auto` = FALSE),
            choicesOpt = list(
                style = paste0("color: ", cell_types_hex, ";")
            )
        )
    })
}

metacell_selector <- function(dataset, ns, id, label, selected = NULL, atlas = FALSE, ...) {
    renderUI({
        metacell_names <- colnames(get_mc_data(dataset(), "mc_mat", atlas = atlas))
        shinyWidgets::pickerInput(ns(id), label,
            choices = metacell_names,
            selected = selected %||% config$selected_mc1,
            multiple = FALSE,
            options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", ...)
        )
    })
}

colored_metacell_selector <- function(dataset, ns, id, label, metacell_colors, metacell_names, selected = NULL) {
    renderUI({
        req(metacell_colors())
        req(metacell_names())
        cell_types_hex <- col2hex(metacell_colors())
        shinyWidgets::pickerInput(ns(id), label,
            choices = metacell_names(),
            selected = selected %||% config$selected_mc1,
            multiple = FALSE, options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith"),
            choicesOpt = list(
                style = paste0("color: ", cell_types_hex, ";")
            )
        )
    })
}

metadata_selector <- function(dataset, ns, id = "selected_md", label = "Metadata", metadata_id = "metadata", selected = NULL, multiple = TRUE) {
    renderUI({
        metadata <- get_mc_data(dataset(), metadata_id)
        req(metadata)
        metadata_fields <- colnames(metadata)[-1]
        shinyWidgets::pickerInput(ns(id), label,
            choices = metadata_fields,
            selected = selected,
            multiple = multiple,
            options = list(`actions-box` = TRUE, `dropup-auto` = FALSE)
        )
    })
}


axis_selector <- function(axis, selected, ns, choices = c("Metadata", "Gene")) {
    fluidRow(
        column(
            width = 9,
            uiOutput(ns(glue("{axis}_select")))
        ),
        column(
            width = 3,
            shinyWidgets::prettyRadioButtons(
                ns(glue("{axis}_type")),
                label = "",
                choices = choices,
                inline = TRUE,
                status = "danger",
                fill = TRUE,
                selected = selected
            )
        )
    )
}

render_axis_select_ui <- function(axis, title, md_choices, md_selected, selected_gene, ns, input, dataset, cell_types = NULL, cell_type_selected = NULL) {
    picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE)

    renderUI({
        req(dataset())
        req(input[[glue("{axis}_type")]])
        if (input[[glue("{axis}_type")]] == "Metadata") {
            shinyWidgets::pickerInput(
                ns(glue("{axis}_var")),
                title,
                choices = md_choices,
                selected = md_selected,
                multiple = FALSE,
                options = picker_options
            )
        } else if (input[[glue("{axis}_type")]] == "Gene") {
            shinyWidgets::pickerInput(
                ns(glue("{axis}_var")),
                title,
                choices = gene_names(dataset()),
                selected = selected_gene,
                multiple = FALSE,
                options = picker_options
            )
        } else if (input[[glue("{axis}_type")]] == "Cell type") {
            req(cell_types)
            shinyWidgets::pickerInput(
                ns(glue("{axis}_var")),
                title,
                choices = cell_types,
                selected = cell_type_selected %||% cell_types[1],
                multiple = FALSE,
                options = picker_options
            )
        }
    })
}

top_correlated_selector <- function(gene_id, id, type_id, input, output, session, dataset, ns, button_labels = c("X", "Y", "Color", "2D"), ids = c("x", "y", "color", "proj2d")) {
    output[[glue("top_correlated_select_{id}")]] <- renderUI({
        req(has_gg_mc_top_cor(project, dataset()))
        req(input[[type_id]] == "Gene")
        req(input[[gene_id]])
        gene <- input[[gene_id]]
        req(gene %in% gene_names(dataset()))
        input_ids <- purrr::map_chr(
            ids,
            ~ {
                ns(glue("select_top_cor_{id}_{.x}"))
            }
        )
        # input_ids <- c(ns(glue("select_top_cor_{id}_x")), ns(glue("select_top_cor_{id}_y")), ns(glue("select_top_cor_{id}_color")), ns(glue("select_top_cor_{id}_proj2d")))
        tagList(
            selectInput(
                ns(glue("selected_top_{id}")),
                glue("Top correlated to {gene}:"),
                choices = c(get_top_cor_gene(dataset(), gene, type = "pos"), rev(get_top_cor_gene(dataset(), gene, type = "neg"))),
                selected = NULL,
                size = 10,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(
                input_ids[1:length(button_labels)],
                labels = button_labels, size = "sm", fullwidth = FALSE
            ),
            shiny::actionButton(
                inputId = ns(glue("genecards_{id}")), label = glue("GeneCards: {gene}"),
                size = "sm", onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}')")
            )
        )
    })
    observeEvent(input[[glue("select_top_cor_{id}_axis")]], {
        req(input[["axis_type"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "axis_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_x")]], {
        req(input[["x_axis_type"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "x_axis_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_y")]], {
        req(input[["y_axis_type"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "y_axis_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_color")]], {
        req(input[["color_by_type"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "color_by_var", selected = input[[glue("selected_top_{id}")]])
    })
    observeEvent(input[[glue("select_top_cor_{id}_proj2d")]], {
        req(input[["color_proj"]] == "Gene")
        shinyWidgets::updatePickerInput(session, "color_proj_gene", selected = input[[glue("selected_top_{id}")]])
    })
}

top_correlated_selectors <- function(input, output, session, dataset, ns, button_labels = c("X", "Y", "Color", "2D")) {
    top_correlated_selector("x_axis_var", "x_axis", "x_axis_type", input, output, session, dataset, ns, button_labels = button_labels)
    top_correlated_selector("y_axis_var", "y_axis", "y_axis_type", input, output, session, dataset, ns, button_labels = button_labels)
    top_correlated_selector("color_by_var", "color_by", "color_by_type", input, output, session, dataset, ns, button_labels = button_labels)
    top_correlated_selector("color_proj_gene", "color_proj", "color_proj", input, output, session, dataset, ns, button_labels = button_labels)
}

axis_vars_ok <- function(dataset, input, md_id, axes = c("x_axis", "y_axis", "color_by"), atlas = FALSE) {
    metadata <- get_mc_data(dataset, md_id, atlas = atlas)
    vars_ok <- purrr::map_lgl(axes, function(v) {
        type <- input[[glue("{v}_type")]]
        var <- input[[glue("{v}_var")]]
        if (type == "Metadata" && (var %in% c(colnames(metadata), "None", "Cell type"))) {
            return(TRUE)
        } else if (type == "Metadata" && atlas && var %in% paste0(colnames(metadata), "_atlas")) {
            return(TRUE)
        } else if (type == "Gene" && var %in% gene_names(dataset)) {
            return(TRUE)
        } else if (type == "Cell type" && (var %in% names(get_cell_type_colors(dataset)))) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    })
    return(all(vars_ok))
}

scatter_selectors <- function(ns, dataset, output, globals, prefix = "gene_gene") {
    output[[glue("{prefix}_point_size_ui")]] <- renderUI({
        req(globals$screen_width)
        req(globals$screen_height)
        numericInput(ns(glue("{prefix}_point_size")), label = "Point size", value = initial_scatters_point_size(dataset(), globals$screen_width, globals$screen_height), min = 0.05, max = 3, step = 0.1)
    })

    output[[glue("{prefix}_stroke_ui")]] <- renderUI({
        numericInput(ns(glue("{prefix}_stroke")), label = "Stroke width", value = initial_scatters_stroke(dataset()), min = 0, max = 3, step = 0.01)
    })
}

projection_selectors <- function(ns, dataset, output, input, globals, weight = 1, atlas = FALSE) {
    output$proj_stat_ui <- renderUI({
        req(input$color_proj == "Gene" || input$color_proj == "Gene A" || input$color_proj == "Gene B")
        selectInput(ns("proj_stat"), label = "Statistic", choices = c("Expression" = "expression", "Enrichment" = "enrichment"), selected = "Expression", multiple = FALSE, selectize = FALSE)
    })

    # Expression range
    output$set_range_ui <- renderUI({
        req(input$color_proj == "Gene" || input$color_proj == "Gene A" || input$color_proj == "Gene B")
        req(input$proj_stat == "expression")
        checkboxInput(ns("set_range"), "Manual range", value = FALSE)
    })

    output$expr_range_ui <- renderUI({
        req(input$color_proj == "Gene" || input$color_proj == "Gene A" || input$color_proj == "Gene B")
        req(input$proj_stat == "expression")
        req(input$set_range)
        shinyWidgets::numericRangeInput(ns("expr_range"), "Expression range", c(-18, -5), width = "80%", separator = " to ")
    })

    # Enrichment range
    output$enrich_range_ui <- renderUI({
        req(input$color_proj == "Gene" || input$color_proj == "Gene A" || input$color_proj == "Gene B")
        req(input$proj_stat == "enrichment")
        shinyWidgets::numericRangeInput(ns("lfp"), "Enrichment range", c(-3, 3), width = "80%", separator = " to ")
    })

    # Point size selectors
    output$point_size_ui <- renderUI({
        req(globals$screen_height)
        req(globals$screen_width)
        req(dataset())
        numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset(), globals$screen_width, globals$screen_height, weight = weight, atlas = atlas), min = 0.1, max = 3, step = 0.1)
    })

    output$stroke_ui <- renderUI({
        numericInput(ns("stroke"), label = "Stroke width", value = initial_proj_stroke(dataset()), min = 0, max = 3, step = 0.01)
    })

    # Minimal edge length selector
    output$edge_distance_ui <- renderUI({
        sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
    })
}
