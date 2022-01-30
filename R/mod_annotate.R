#' annotate UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_annotate_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 8,
                style = "padding-right:0px;",
                shinydashboardPlus::box(
                    id = ns("gene_projection"),
                    title = "2D Projection",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 80,
                        id = ns("gene_projection_sidebar"),
                        shinyWidgets::prettyRadioButtons(
                            ns("color_proj"),
                            label = "Color by:",
                            choices = c("Cell type", "Gene", "Metadata", "Selected"),
                            inline = TRUE,
                            status = "danger",
                            fill = TRUE
                        ),
                        checkboxInput(ns("show_selected_metacells"), "Show selected metacells", value = FALSE),
                        uiOutput(ns("gene_selector")),
                        uiOutput(ns("metadata_selector")),
                        uiOutput(ns("proj_stat_ui")),
                        uiOutput(ns("set_range_ui")),
                        uiOutput(ns("expr_range_ui")),
                        uiOutput(ns("enrich_range_ui")),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("stroke_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_proj_2d"))
                    )
                ),
                fluidRow(
                    uiOutput(ns("time_box_ui_column")),
                    column(
                        width = 6,
                        offset = 0,
                        style = "padding-right:0px;",
                        shinydashboardPlus::box(
                            id = ns("gene_gene_box"),
                            title = "Gene/Gene",
                            status = "primary",
                            solidHeader = TRUE,
                            collapsible = TRUE,
                            closable = FALSE,
                            width = 12,
                            sidebar = shinydashboardPlus::boxSidebar(
                                startOpen = FALSE,
                                width = 100,
                                id = ns("gene_gene_sidebar"),
                                axis_selector("x_axis", "Gene", ns),
                                axis_selector("y_axis", "Gene", ns),
                                axis_selector("color_by", "Metadata", ns),
                                uiOutput(ns("gene_gene_xyline_ui")),
                                uiOutput(ns("gene_gene_fixed_limits_ui")),
                                uiOutput(ns("gene_gene_point_size_ui")),
                                uiOutput(ns("gene_gene_stroke_ui"))
                            ),
                            shinycssloaders::withSpinner(
                                plotly::plotlyOutput(ns("plot_gene_gene_mc"))
                            )
                        )
                    ),
                    column(
                        width = 6,
                        offset = 0,
                        style = "padding-left:0px;",
                        shinydashboardPlus::box(
                            id = ns("mc_mc_box"),
                            title = "Differential expression",
                            status = "primary",
                            solidHeader = TRUE,
                            collapsible = TRUE,
                            closable = FALSE,
                            width = 12,
                            sidebar = shinydashboardPlus::boxSidebar(
                                startOpen = FALSE,
                                width = 80,
                                id = ns("mc_mc_sidebar"),
                                shinyWidgets::radioGroupButtons(
                                    inputId = ns("mode"),
                                    label = "Compare:",
                                    choices = c(
                                        "MCs",
                                        "Types"
                                    ),
                                    justified = TRUE
                                ),
                                uiOutput(ns("metacell1_select")),
                                uiOutput(ns("metacell2_select")),
                                shinyWidgets::actionGroupButtons(ns("switch_metacells"), labels = c("Switch"), size = "sm")
                            ),
                            shinycssloaders::withSpinner(
                                plotly::plotlyOutput(ns("plot_mc_mc_gene_scatter"))
                            ),
                            shinyWidgets::prettySwitch(inputId = ns("show_diff_expr_table"), value = FALSE, label = "Show table"),
                            DT::DTOutput(ns("diff_expr_table"))
                        )
                    )
                )
            ),
            column(
                width = 4,
                style = "padding-right:0px; padding-left:0px;",
                shinydashboardPlus::box(
                    id = ns("metacell_types_box"),
                    title = "Metacell annotation",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    splitLayout(
                        fileInput(ns("metacell_types_fn"),
                            label = NULL,
                            buttonLabel = "Load",
                            multiple = FALSE,
                            accept =
                                c(
                                    "text/csv",
                                    "text/comma-separated-values,text/plain",
                                    "text/tab-separated-values",
                                    ".csv",
                                    ".tsv"
                                )
                        ),
                        actionButton(ns("reset_metacell_types"), "Reset", style = "align-items: center;"),
                        downloadButton(ns("metacell_types_download"), "Export", style = "align-items: center;")
                    ),
                    uiOutput(ns("annotation_box")),
                    uiOutput(ns("update_all_selectors")),
                    shinycssloaders::withSpinner(
                        DT::dataTableOutput(ns("mc_type_table"))
                    )
                )
            ),
            column(
                width = 4,
                style = "padding-right:0px; padding-left:0px;",
                shinydashboardPlus::box(
                    id = ns("cell_type_colors"),
                    title = "Cell Types",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    splitLayout(
                        actionButton(ns("reset_cell_type_colors"), "Reset", style = "align-items: center;"),
                        downloadButton(ns("cell_type_colors_download"), "Export", style = "align-items: center;"),
                        actionButton(ns("delete_cell_type_colors"), "Delete"),
                        actionButton(ns("add_cell_type_colors"), "Add")
                    ),
                    uiOutput(ns("annot_color_picker")),
                    shinycssloaders::withSpinner(
                        DT::dataTableOutput(ns("cell_type_table"))
                    )
                )
            )
        )
    )
}


#' annotate sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_annotate_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("top_correlated_select_x_axis")),
        uiOutput(ns("top_correlated_select_y_axis")),
        uiOutput(ns("top_correlated_select_color_by")),
        uiOutput(ns("top_correlated_select_color_proj"))
    )
}


#' annotate Server Function
#'
#' @noRd
mod_annotate_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals) {
    ns <- session$ns

    # gene selectors
    values <- reactiveValues(file_status = NULL)
    top_correlated_selectors(input, output, session, dataset, ns)

    picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith", dropupAuto = FALSE)

    output$gene_selector <- renderUI({
        shinyWidgets::pickerInput(
            ns("color_proj_gene"),
            label = "Gene:",
            choices = gene_names(dataset()),
            selected = default_gene1,
            width = "70%",
            multiple = FALSE,
            options = picker_options
        )
    })

    output$metadata_selector <- renderUI({
        if (!has_metadata(dataset())) {
            print(glue("Dataset doesn't have any metadata."))
        } else {
            shinyWidgets::pickerInput(
                ns("color_proj_metadata"),
                label = "Metadata:",
                choices = dataset_metadata_fields(dataset()),
                selected = dataset_metadata_fields(dataset())[1],
                width = "70%",
                multiple = FALSE,
                options = picker_options
            )
        }
    })

    observe({
        req(input$color_proj)
        shinyjs::toggle(id = "gene_selector", condition = input$color_proj == "Gene")
        shinyjs::toggle(id = "metadata_selector", condition = input$color_proj == "Metadata")
    })


    scatter_selectors(ns, dataset, output, globals)
    projection_selectors(ns, dataset, output, input, globals, weight = 0.6)

    observeEvent(input$metacell_types_fn, {
        values$file_status <- "uploaded"
    })

    observe({
        req(input$metacell_types_fn)
        req(values$file_status)
        new_metacell_types <- tgutil::fread(input$metacell_types_fn$datapath, colClasses = c("cell_type" = "character", "metacell" = "character")) %>% as_tibble()

        input_ok <- TRUE
        required_fields <- c("cell_type", "metacell")
        if (!all(required_fields %in% colnames(new_metacell_types))) {
            showNotification(glue("Please provide a file with the following fields: cell_type, metacell"), type = "error")
            input_ok <- FALSE
        }

        metacells <- get_metacell_ids(project, dataset())

        unknown_metacells <- new_metacell_types$metacell[!(new_metacell_types$metacell %in% metacells)]
        if (length(unknown_metacells) > 0) {
            mcs <- paste(unknown_metacells, collapse = ", ")
            showNotification(glue("Metacell types contains metacells that are missing from the data: {mcs}"), type = "error")
            input_ok <- FALSE
        }

        missing_metacells <- metacells[!(metacells %in% new_metacell_types$metacell)]
        if (length(missing_metacells) > 0) {
            mcs <- paste(missing_metacells, collapse = ", ")
            showNotification(glue("Some metacells are missing from metacell types: {mcs}"), type = "warning")
        }

        if (input_ok) {
            if (has_name(new_metacell_types, "color")) {
                new_cell_type_colors <- new_metacell_types %>%
                    distinct(cell_type, color) %>%
                    select(cell_type, color) %>%
                    filter(cell_type != "(Missing)") %>%
                    arrange(cell_type) %>%
                    mutate(order = 1:n())

                cell_type_colors(new_cell_type_colors)
            }

            cur_metacell_types <- metacell_types()
            new_metacell_types <- cur_metacell_types %>%
                select(-any_of(c("cell_type"))) %>%
                left_join(new_metacell_types %>% select(metacell, cell_type), by = "metacell") %>%
                mutate(cell_type = ifelse(cell_type == "(Missing)", NA, cell_type)) %>%
                mutate(cell_type = as.character(forcats::fct_explicit_na(factor(cell_type))))

            new_metacell_types <- sanitize_metacell_types(new_metacell_types, cell_type_colors(), dataset())
            metacell_types(new_metacell_types)
            values$file_status <- NULL
        }
    })

    # export metacell types file
    output$metacell_types_download <- downloadHandler(
        filename = function() {
            paste("metacell_types-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            fwrite(
                metacell_types() %>%
                    select(metacell, cell_type, top1_gene, top1_lfp, top2_gene, top2_lfp) %>%
                    left_join(cell_type_colors() %>% select(cell_type, color), by = "cell_type"),
                file
            )
        }
    )

    # export cell type colors file
    output$cell_type_colors_download <- downloadHandler(
        filename = function() {
            paste("cell_type_colors-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            fwrite(
                cell_type_colors() %>%
                    select(cell_type, color),
                file
            )
        }
    )

    # set reactive values
    selected_metacell_types <- reactiveVal(tibble(metacell = character(), cell_type = character()))
    to_show <- reactiveVal()

    # keep the last cell type that was chosen in order for it to be defaultly selected
    last_chosen_cell_type <- reactiveVal("(Missing)")

    observeEvent(input$reset_metacell_types, {
        metacell_types(get_mc_data(dataset(), "metacell_types"))
        selected_metacell_types(tibble(metacell = character(), cell_type = character()))
        to_show(NULL)
        last_chosen_cell_type("(Missing)")
        values$file_status <- NULL
    })

    observeEvent(input$reset_cell_type_colors, {
        cell_type_colors(get_mc_data(dataset(), "cell_type_colors"))
    })

    output$annotation_box <- renderUI({
        if (nrow(selected_metacell_types()) == 0) {
            textOutput(ns("please_select_metacells"))
        } else {
            list(
                textOutput(ns("number_of_selected_metacells")),
                actionButton(ns("update_annotation"), "Apply"),
                actionButton(ns("reset_annotation"), "Reset Selection"),
                shinyWidgets::radioGroupButtons(
                    inputId = ns("update_option"),
                    label = "",
                    choices = c(
                        "Change all",
                        "Change selected"
                    ),
                    justified = TRUE
                )
            )
        }
    })

    output$number_of_selected_metacells <- renderPrint(glue("Selected {nrow(selected_metacell_types())} metacells"))
    output$please_select_metacells <- renderPrint(glue("Please select metacells"))

    output$update_all_selectors <- renderUI({
        req(nrow(selected_metacell_types()) > 0)
        shinyWidgets::pickerInput(ns("selected_cell_type_update_all"), "Cell type", choices = c("(Missing)", cell_type_colors() %>% pull(cell_type) %>% as.character() %>% unique() %>% sort()), multiple = FALSE, selected = last_chosen_cell_type())
    })


    observeEvent(input$update_annotation, {
        new_metacell_types <- metacell_types()
        changed <- FALSE

        req(input$update_option)
        req(input$selected_cell_type_update_all)
        req(selected_metacell_types())

        if (input$update_option == "Change all") {
            req(input$mc_type_table_rows_all)
            metacells <- selected_metacell_types()[input$mc_type_table_rows_all, ] %>% pull(metacell)
        } else {
            req(input$mc_type_table_rows_selected)
            metacells <- selected_metacell_types()[input$mc_type_table_rows_selected, ] %>% pull(metacell)
        }

        new_metacell_types <- new_metacell_types %>% mutate(
            cell_type = ifelse(metacell %in% metacells, input$selected_cell_type_update_all, cell_type)
        )
        new_selected_annot <- selected_metacell_types()
        new_selected_annot <- new_selected_annot %>% mutate(
            cell_type = ifelse(metacell %in% metacells, input$selected_cell_type_update_all, cell_type),
        )
        selected_metacell_types(new_selected_annot)
        last_chosen_cell_type(input$selected_cell_type_update_all)
        changed <- TRUE

        if (changed) {
            metacell_types(new_metacell_types)
        }
    })

    observeEvent(input$reset_annotation, {
        selected_metacell_types(tibble(metacell = character(), cell_type = character()))
        to_show(NULL)
    })


    observe({
        req(metacell_types)
        if (nrow(selected_metacell_types()) == 0) {
            to_show(NULL)
        }

        req(nrow(selected_metacell_types()) != 0)
        to_show_new <- metacell_types() %>%
            select(metacell, cell_type) %>%
            filter(metacell %in% selected_metacell_types()$metacell)

        to_show(to_show_new)
    })

    output$mc_type_table <- DT::renderDataTable(
        to_show(),
        escape = FALSE,
        server = FALSE,
        rownames = FALSE,
        filter = "top",
        options = list(
            dom = "t",
            paging = FALSE,
            language = list(emptyTable = "Please select metacells")
        )
    )

    output$cell_type_table <- DT::renderDataTable(
        DT::datatable(cell_type_colors() %>% select(cell_type, color),
            editable = "cell",
            rownames = FALSE,
            options = list(
                paging = FALSE
            )
        ) %>%
            DT::formatStyle(
                "color", "cell_type",
                backgroundColor = DT::styleEqual(
                    cell_type_colors()$cell_type,
                    col2hex(cell_type_colors()$color)
                )
            ),
        server = TRUE # see https://github.com/rstudio/DT/issues/598
    )

    cell_type_table_proxy <- DT::dataTableProxy("cell_type_table")

    observeEvent(input$cell_type_table_cell_edit, {
        # fix column number to be 1 based
        new_input <- input$cell_type_table_cell_edit %>% mutate(col = col + 1)
        edited_data <- DT::editData(cell_type_colors() %>% select(cell_type, color), new_input, "cell_type_table")

        # the data table was given only cell_type and color columns so we add the rest
        edited_data <- bind_cols(
            cell_type_colors() %>%
                select(-cell_type, -color),
            edited_data
        ) %>%
            select(cell_type, color, order)

        # change corresponding metacell_type entries
        if (new_input$col == 1) {
            old_cell_type <- as.character(cell_type_colors()$cell_type[new_input$row])
            new_cell_type <- as.character(edited_data$cell_type[new_input$row])
            new_metacell_types <- metacell_types() %>% mutate(cell_type = ifelse(cell_type == old_cell_type, new_cell_type, cell_type))
            metacell_types(new_metacell_types)
        }


        cell_type_colors(edited_data)

        DT::replaceData(cell_type_table_proxy, cell_type_colors() %>% select(cell_type, color), resetPaging = FALSE)
    })

    observeEvent(input$delete_cell_type_colors, {
        rows <- input$cell_type_table_rows_selected

        if (!is.null(rows) && length(rows) > 0) {
            to_delete <- cell_type_colors()[rows, ]
            cell_type_colors(cell_type_colors()[-rows, ])
            metacell_types(
                metacell_types() %>%
                    mutate(
                        cell_type = ifelse(cell_type %in% to_delete$cell_type, NA, cell_type),
                    )
            )
        }
    })

    observeEvent(input$add_cell_type_colors, {
        rows <- input$cell_type_table_rows_selected
        if (!is.null(rows) && length(rows) > 0) {
            place <- rows[1] + 1
        } else {
            place <- 1
        }

        # TODO: allow creation of more than one new row without editing (e.g. by adding a suffix to cell_type)
        new_data <- cell_type_colors() %>% arrange(order)

        new_name <- vctrs::vec_as_names(c("Cell type", new_data$cell_type), repair = "unique") %>%
            setdiff(new_data$cell_type) %>%
            sort() %>%
            tail(1)

        new_row <- tibble(cell_type = new_name, color = "red", order = place)
        new_data <- bind_rows(
            new_data %>% filter(order < place),
            new_row,
            new_data %>% filter(order >= place) %>% mutate(order = order + 1)
        )

        new_data <- new_data %>%
            arrange(order) %>%
            distinct(cell_type, .keep_all = TRUE) %>%
            mutate(order = 1:n())
        cell_type_colors(new_data)
    })

    output$annot_color_picker <- renderUI({
        fluidRow(
            column(6, actionButton(ns("submit_new_color"), "Change color")),
            column(6, colourpicker::colourInput(ns("selected_new_color"), NULL, "black"))
        )
    })

    observe({
        req(input$cell_type_table_rows_selected)
        row <- tail(input$cell_type_table_rows_selected, n = 1)
        colourpicker::updateColourInput(session, "selected_new_color", value = cell_type_colors()$color[row])
    })

    observeEvent(input$submit_new_color, {
        rows <- input$cell_type_table_rows_selected
        new_data <- cell_type_colors()
        new_data$color[rows] <- input$selected_new_color
        cell_type_colors(new_data)
    })


    # Select metacell when clicking on it
    observe_mc_click_event("proj_annot_plot", input, session, cell_type_colors, metacell_types, selected_metacell_types)
    observe_mc_click_event("gene_gene_plot_annot", input, session, cell_type_colors, metacell_types, selected_metacell_types)
    observe_mc_click_event("gene_time_mc_plot1_annot", input, session, cell_type_colors, metacell_types, selected_metacell_types)
    observe_mc_click_event("gene_time_mc_plot2_annot", input, session, cell_type_colors, metacell_types, selected_metacell_types)

    # Select multiple metacells
    observer_mc_select_event("proj_annot_plot", input, cell_type_colors, metacell_types, selected_metacell_types)
    observer_mc_select_event("gene_gene_plot_annot", input, cell_type_colors, metacell_types, selected_metacell_types)
    observer_mc_select_event("gene_time_mc_plot1_annot", input, cell_type_colors, metacell_types, selected_metacell_types)
    observer_mc_select_event("gene_time_mc_plot2_annot", input, cell_type_colors, metacell_types, selected_metacell_types)

    projection_selectors(ns, dataset, output, input, globals, weight = 0.6)
    scatter_selectors(ns, dataset, output, globals)

    # Projection plots
    output$plot_gene_proj_2d <- render_2d_plotly(
        input,
        output,
        session,
        dataset,
        metacell_types,
        cell_type_colors,
        source = "proj_annot_plot",
        buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"),
        dragmode = "select",
        selected_metacell_types = selected_metacell_types
    )

    # Metadata/Metadata plots
    output$x_axis_select <- render_axis_select_ui("x_axis", "X axis", md_choices = dataset_metadata_fields_numeric(dataset()), md_selected = dataset_metadata_fields_numeric(dataset())[1], selected_gene = default_gene1, input = input, ns = ns, dataset = dataset) %>% bindCache(dataset(), ns, ns("x_axis"), input$x_axis_type)

    output$y_axis_select <- render_axis_select_ui("y_axis", "Y axis", md_choices = dataset_metadata_fields_numeric(dataset()), md_selected = dataset_metadata_fields_numeric(dataset())[2], selected_gene = default_gene2, input = input, ns = ns, dataset = dataset) %>% bindCache(dataset(), ns, ns("y_axis"), input$y_axis_type)

    output$color_by_select <- render_axis_select_ui("color_by", "Color", md_choices = c("Cell type", dataset_metadata_fields_numeric(dataset())), md_selected = "Cell type", selected_gene = default_gene1, input = input, ns = ns, dataset = dataset) %>% bindCache(dataset(), ns, ns("color_by"), input$color_by_type)

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
        req(axis_vars_ok(dataset(), input, "metadata"))

        color_var <- input$color_by_var
        if (input$color_by_var == "Cell type") {
            color_var <- NULL
        }

        fig <- plot_mc_scatter(
            dataset(),
            input$x_axis_var,
            input$y_axis_var,
            color_var,
            x_type = input$x_axis_type,
            y_type = input$y_axis_type,
            color_type = input$color_by_type,
            metacell_types = metacell_types(),
            cell_type_colors = cell_type_colors(),
            point_size = input$gene_gene_point_size,
            stroke = input$gene_gene_stroke,
            plot_text = FALSE,
            fixed_limits = input$gene_gene_fixed_limits,
            xyline = input$gene_gene_xyline %||% FALSE
        ) %>%
            plotly::ggplotly(tooltip = "tooltip_text", source = "gene_gene_plot_annot") %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons(buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines")) %>%
            plotly::layout(dragmode = "select")

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
    }) %>% bindCache(dataset(), input$x_axis_var, input$x_axis_type, input$y_axis_var, input$y_axis_type, input$color_by_type, input$color_by_var, metacell_types(), cell_type_colors(), input$gene_gene_point_size, input$gene_gene_stroke, input$gene_gene_fixed_limits, input$gene_gene_xyline)



    connect_gene_plots(input, output, session, ns, source = "proj_annot_plot")

    # MC/MC diff gene expression plots
    metacell_names <- metacell_names_reactive(dataset)
    metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)

    group_selectors(input, output, session, dataset, ns)
    metacell_selectors(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors, )

    mc_mc_gene_scatter_df <- mc_mc_gene_scatter_df_reactive(dataset, input, output, session, metacell_types, cell_type_colors)

    diff_expr_switch_metacells(dataset, input, output, session)

    output$plot_mc_mc_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, mc_mc_gene_scatter_df, metacell_names(), cell_type_colors(), source_suffix = "_annot")

    mod_gene_mc_plotly_observers(input, session, source = "mc_mc_plot_annot", notification_suffix = "")

    output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, mc_mc_gene_scatter_df)
}


observe_mc_click_event <- function(source, input, session, cell_type_colors, metacell_types, selected_metacell_types) {
    observeEvent(plotly::event_data("plotly_click", source = source), {
        el <- plotly::event_data("plotly_click", source = source)

        selected_metacell <- el$customdata

        new_selected_annot <- metacell_types() %>% filter(metacell == selected_metacell)

        selected_metacell_types(
            bind_rows(
                selected_metacell_types(),
                new_selected_annot
            ) %>% distinct(metacell, cell_type)
        )

        shinyWidgets::updatePickerInput(session, "metacell1", selected = selected_metacell)
    })
}

observer_mc_select_event <- function(source, input, cell_type_colors, metacell_types, selected_metacell_types) {
    observeEvent(plotly::event_data("plotly_selected", source = source), {
        el <- plotly::event_data("plotly_selected", source = source)

        selected_metacells <- unique(el$customdata)

        new_selected_annot <- metacell_types() %>% filter(metacell %in% selected_metacells)
        selected_metacell_types(
            bind_rows(
                selected_metacell_types(),
                new_selected_annot
            ) %>% distinct(metacell, cell_type)
        )
    })
}
