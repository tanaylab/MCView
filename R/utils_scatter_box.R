scatter_box <- function(ns, id, title = "Gene/Gene", x_selected = "Gene", y_selected = "Gene", color_selected = "Metadata") {
    shinydashboardPlus::box(
        id = ns(id),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        sidebar = shinydashboardPlus::boxSidebar(
            startOpen = FALSE,
            width = 100,
            id = ns("gene_gene_sidebar"),
            axis_selector("x_axis", x_selected, ns),
            axis_selector("y_axis", y_selected, ns),
            axis_selector("color_by", color_selected, ns),
            uiOutput(ns("gene_gene_xyline_ui")),
            uiOutput(ns("gene_gene_fixed_limits_ui")),
            uiOutput(ns("use_atlas_limits_ui")),
            uiOutput(ns("gene_gene_point_size_ui")),
            uiOutput(ns("gene_gene_stroke_ui"))
        ),
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("plot_gene_gene_mc"))
        )
    )
}

scatter_box_outputs <- function(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, ns, plotly_source, plotly_buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = NULL) {
    output$x_axis_select <- render_axis_select_ui("x_axis", "X axis", md_choices = dataset_metadata_fields_numeric(dataset()), md_selected = dataset_metadata_fields_numeric(dataset())[1], selected_gene = default_gene1, input = input, ns = ns, dataset = dataset, gene_modules = gene_modules) %>% bindCache(dataset(), ns, ns("x_axis"), input$x_axis_type, gene_modules())

    output$y_axis_select <- render_axis_select_ui("y_axis", "Y axis", md_choices = dataset_metadata_fields_numeric(dataset()), md_selected = dataset_metadata_fields_numeric(dataset())[2], selected_gene = default_gene2, input = input, ns = ns, dataset = dataset, gene_modules = gene_modules) %>% bindCache(dataset(), ns, ns("y_axis"), input$y_axis_type, gene_modules())

    output$color_by_select <- render_axis_select_ui("color_by", "Color", md_choices = c("Cell type", dataset_metadata_fields(dataset())), md_selected = "Cell type", selected_gene = default_gene1, input = input, ns = ns, dataset = dataset, gene_modules = gene_modules) %>% bindCache(dataset(), ns, ns("color_by"), input$color_by_type, gene_modules())

    output$use_atlas_limits_ui <- renderUI({
        req(has_atlas(dataset()))
        checkboxInput(ns("use_atlas_limits"), label = "Use atlas limits", value = FALSE)
    })

    output$plot_gene_gene_mc <- plotly::renderPlotly({
        req(input$x_axis_var)
        req(input$y_axis_var)
        req(input$color_by_var)
        req(input$x_axis_type)
        req(input$y_axis_type)
        req(input$color_by_type)
        req(input$gene_gene_point_size)
        req(input$gene_gene_stroke)
        req(!is.null(input$gene_gene_fixed_limits))
        req(axis_vars_ok(dataset(), input, "metadata", gene_modules))

        color_var <- input$color_by_var
        if (input$color_by_var == "Cell type") {
            color_var <- NULL
        }

        x_limits <- NULL
        y_limits <- NULL
        if (!is.null(input$use_atlas_limits) && input$use_atlas_limits) {
            req(has_atlas(dataset()))
            if (input$x_axis_type == "Gene") {
                egc_x <- get_gene_egc(input$x_axis_var, dataset(), atlas = TRUE) + egc_epsilon
                x_limits <- c(min(egc_x), max(egc_x))
            }

            if (input$y_axis_type == "Gene") {
                egc_y <- get_gene_egc(input$y_axis_var, dataset(), atlas = TRUE) + egc_epsilon
                y_limits <- c(min(egc_y), max(egc_y))
            }
        }

        fig <- plot_mc_scatter(
            dataset(),
            input$x_axis_var,
            input$y_axis_var,
            color_var,
            gene_modules = gene_modules(),
            x_type = input$x_axis_type,
            y_type = input$y_axis_type,
            color_type = input$color_by_type,
            metacell_types = metacell_types(),
            cell_type_colors = cell_type_colors(),
            point_size = input$gene_gene_point_size,
            stroke = input$gene_gene_stroke,
            plot_text = FALSE,
            x_limits = x_limits,
            y_limits = y_limits,
            fixed_limits = input$gene_gene_fixed_limits,
            xyline = input$gene_gene_xyline %||% FALSE
        ) %>%
            plotly::ggplotly(tooltip = "tooltip_text", source = plotly_source) %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons(buttons = plotly_buttons)

        if (!is.null(dragmode)) {
            fig <- fig %>%
                plotly::layout(dragmode = dragmode)
        }

        if (input$color_by_var == "Cell type") {
            fig <- plotly::hide_legend(fig)
        } else {
            # This ugly hack is due to https://github.com/ropensci/plotly/issues/1234
            # We need to remove the legend generated by scale_color_identity
            fig$x$data <- fig$x$data %>% purrr::map(~ {
                .x$showlegend <- FALSE
                .x
            })
        }

        return(fig)
    }) %>% bindCache(dataset(), input$x_axis_var, input$x_axis_type, input$y_axis_var, input$y_axis_type, input$color_by_type, input$color_by_var, metacell_types(), cell_type_colors(), gene_modules(), input$gene_gene_point_size, input$gene_gene_stroke, input$use_atlas_limits, input$gene_gene_fixed_limits, input$gene_gene_xyline, dragmode, plotly_buttons)
}

axis_selector <- function(axis, selected, ns, choices = c("Metadata", "Gene", "Gene module")) {
    fluidRow(
        column(
            width = 7,
            uiOutput(ns(glue("{axis}_select")))
        ),
        column(
            width = 5,
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

axis_vars_ok <- function(dataset, input, md_id, gene_modules, axes = c("x_axis", "y_axis", "color_by"), atlas = FALSE) {
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
        } else if (type == "Gene module" && var %in% levels(gene_modules()$module)) {
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

    output[[glue("{prefix}_fixed_limits_ui")]] <- renderUI({
        checkboxInput(ns(glue("{prefix}_fixed_limits")), label = "X axis limits = Y axis limits", value = FALSE)
    })

    output[[glue("{prefix}_xyline_ui")]] <- renderUI({
        checkboxInput(ns(glue("{prefix}_xyline")), label = "Show X = Y line", value = FALSE)
    })
}

render_axis_select_ui <- function(axis, title, md_choices, md_selected, selected_gene, ns, input, dataset, gene_modules, cell_types = NULL, cell_type_selected = NULL) {
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
        } else if (input[[glue("{axis}_type")]] == "Gene module") {
            shinyWidgets::pickerInput(
                ns(glue("{axis}_var")),
                title,
                choices = levels(gene_modules()$module),
                selected = levels(gene_modules()$module)[1],
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