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
                projection_box(ns, "gene_projection", checkboxInput(ns("show_selected_metacells"), "Show selected metacells", value = FALSE), color_choices = c("Cell type", "Gene", "Gene module", "Metadata", "Selected")),
                fluidRow(
                    uiOutput(ns("time_box_ui_column")),
                    column(
                        width = 6,
                        offset = 0,
                        style = "padding-right:0px;",
                        scatter_box(ns, "gene_gene_box")
                    ),
                    column(
                        width = 6,
                        offset = 0,
                        style = "padding-left:0px;",
                        diff_expr_box(ns, "mc_mc_box")
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
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 50,
                        id = ns("metacell_types_box_sidebar"),
                        checkboxInput(ns("add_to_selection"), label = "Add to\ncurrent selection", value = TRUE),
                        checkboxInput(ns("reset_on_apply"), label = "Reset selection\non apply", value = FALSE)
                    ),
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
                        actionButton(ns("paste_metacells"), "Paste", style = "align-items: center;"),
                        actionButton(ns("copy_metacells"), "Copy", style = "align-items: center;"),
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
                        actionButton(ns("reset_cell_type_colors"), "Reset", style = "align-items: left;"),
                        downloadButton(ns("cell_type_colors_download"), "Export", style = "align-items: left;"),
                        actionButton(ns("add_cell_type_modal"), "Add", style = "align-items: left;")
                    ),
                    br(),
                    br(),
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
mod_annotate_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            # gene selectors
            values <- reactiveValues(file_status = NULL)
            top_correlated_selectors(input, output, session, dataset, ns)
            scatter_selectors(ns, dataset, output, globals)

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
                metacell_types(get_metacell_types_data(dataset()))
                selected_metacell_types(tibble(metacell = character(), cell_type = character()))
                to_show(NULL)
                last_chosen_cell_type("(Missing)")
                values$file_status <- NULL
            })

            observeEvent(input$paste_metacells, {
                selected_metacells <- unique(globals$clipboard)

                new_selected_annot <- metacell_types() %>% filter(metacell %in% selected_metacells)
                if (!is.null(input$add_to_selection) && input$add_to_selection) {
                    selected_metacell_types(
                        bind_rows(
                            selected_metacell_types(),
                            new_selected_annot
                        ) %>% distinct(metacell, cell_type)
                    )
                } else {
                    selected_metacell_types(new_selected_annot %>% distinct(metacell, cell_type))
                }
            })

            observeEvent(input$copy_metacells, {
                selected_metacells <- selected_metacell_types()$metacell
                globals$clipboard <- selected_metacells
                showNotification(glue("Copied {length(selected_metacells)} metacells to clipboard"))
            })

            observeEvent(input$reset_cell_type_colors, {
                cell_type_colors(get_cell_type_data(dataset()))
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

                req(input$reset_on_apply)
                if (input$reset_on_apply) {
                    selected_metacell_types(tibble(metacell = character(), cell_type = character()))
                    to_show(NULL)
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
                edited_data <- DT::editData(cell_type_colors() %>% select(cell_type, color), new_input)
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

                if (!any(duplicated(edited_data$cell_type))) {
                    cell_type_colors(edited_data)
                } else {
                    dups <- paste(edited_data$cell_type[duplicated(edited_data$cell_type)], collapse = ", ")
                    showNotification(glue("Cell types cannot be duplicated: {dups}"))
                }

                DT::replaceData(cell_type_table_proxy, cell_type_colors() %>% select(cell_type, color), resetPaging = FALSE, rownames = FALSE)
            })

            observeEvent(input$merge_cell_types_modal, {
                rows <- input$cell_type_table_rows_selected
                req(rows)
                cell_types <- cell_type_colors()$cell_type[input$cell_type_table_rows_selected]
                default_color <- cell_type_colors()$color[input$cell_type_table_rows_selected[1]]
                showModal({
                    modalDialog(
                        title = "Merge cell types",
                        textInput(ns("new_merged_cell_type_name"), "Cell type name"),
                        colourpicker::colourInput(ns("new_merged_cell_type_color"), NULL, default_color),
                        glue("Are you sure you want to merge the following cell types: {paste(cell_types, collapse = ',')}?"),
                        footer = tagList(
                            modalButton("Cancel"),
                            actionButton(ns("merge_cell_types"), "OK")
                        )
                    )
                })
            })

            observeEvent(input$merge_cell_types, {
                rows <- input$cell_type_table_rows_selected
                req(rows)
                cell_types <- cell_type_colors()$cell_type[input$cell_type_table_rows_selected]
                req(input$new_merged_cell_type_name)
                req(input$new_merged_cell_type_color)
                if (input$new_merged_cell_type_name %in% cell_type_colors()$cell_type[-rows]) {
                    showNotification(glue("Cell type {input$new_merged_cell_type_name} already exists"), type = "error")
                    removeModal()
                    req(FALSE)
                }

                new_cell_type_colors <- cell_type_colors() %>%
                    filter(!(cell_type %in% cell_types)) %>%
                    tibble::add_row(cell_type = input$new_merged_cell_type_name, color = input$new_merged_cell_type_color, order = rows[1], .before = rows[1]) %>%
                    arrange(order) %>%
                    distinct(cell_type, .keep_all = TRUE) %>%
                    mutate(order = 1:n())

                cell_type_colors(new_cell_type_colors)

                new_metacell_types <- metacell_types() %>%
                    mutate(cell_type = ifelse(cell_type %in% cell_types, input$new_merged_cell_type_name, cell_type))
                metacell_types(new_metacell_types)

                removeModal()
            })

            observeEvent(input$delete_cell_type_colors_modal, {
                req(input$cell_type_table_rows_selected)

                cell_types <- paste(cell_type_colors()$cell_type[input$cell_type_table_rows_selected], collapse = ", ")
                showModal({
                    modalDialog(
                        title = "Remove cell type(s)",
                        glue("Are you sure you want to delete the following cell types: {cell_types}?"),
                        footer = tagList(
                            modalButton("Cancel"),
                            actionButton(ns("delete_cell_type_colors"), "OK")
                        )
                    )
                })
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
                removeModal()
            })

            observeEvent(input$add_cell_type_modal, {
                showModal({
                    modalDialog(
                        title = "Add a new cell type",
                        textInput(ns("new_cell_type_name"), "Cell type name"),
                        colourpicker::colourInput(ns("new_cell_type_color"), NULL, "red"),
                        footer = tagList(
                            modalButton("Cancel"),
                            actionButton(ns("add_cell_type"), "OK")
                        )
                    )
                })
            })

            observeEvent(input$add_cell_type, {
                req(input$new_cell_type_name)
                req(input$new_cell_type_color)
                if (input$new_cell_type_name %in% cell_type_colors()$cell_type) {
                    showNotification(glue("Cell type {input$new_cell_type_name} already exists"), type = "error")
                    removeModal()
                    req(FALSE)
                }

                rows <- input$cell_type_table_rows_selected
                if (!is.null(rows) && length(rows) > 0) {
                    place <- rows[1] + 1
                } else {
                    place <- 1
                }

                new_data <- cell_type_colors() %>% arrange(order)

                new_row <- tibble(cell_type = input$new_cell_type_name, color = input$new_cell_type_color, order = place)
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
                removeModal()
            })

            output$annot_color_picker <- renderUI({
                fluidRow(
                    column(3, actionButton(ns("submit_new_color"), "Change color")),
                    column(3, colourpicker::colourInput(ns("selected_new_color"), NULL, "black")),
                    column(3, actionButton(ns("delete_cell_type_colors_modal"), "Delete")),
                    column(3, actionButton(ns("merge_cell_types_modal"), "Merge"))
                )
            })

            observe({
                shinyjs::toggle(id = "submit_new_color", condition = !is.null(input$cell_type_table_rows_selected))
                shinyjs::toggle(id = "selected_new_color", condition = !is.null(input$cell_type_table_rows_selected))
                shinyjs::toggle(id = "delete_cell_type_colors_modal", condition = !is.null(input$cell_type_table_rows_selected))
                shinyjs::toggle(id = "merge_cell_types_modal", condition = !is.null(input$cell_type_table_rows_selected) && length(input$cell_type_table_rows_selected) > 1)
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

            projection_selectors(ns, dataset, output, input, gene_modules, globals, weight = 0.6)
            scatter_selectors(ns, dataset, output, globals)

            # Projection plots
            output$plot_gene_proj_2d <- render_2d_plotly(
                input,
                output,
                session,
                dataset,
                metacell_types,
                cell_type_colors,
                gene_modules,
                globals,
                source = "proj_annot_plot",
                buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"),
                dragmode = "select",
                selected_metacell_types = selected_metacell_types
            )

            scatter_box_outputs(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, ns, plotly_source = "gene_gene_plot_annot", plotly_buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = "select")

            connect_gene_plots(input, output, session, ns, source = "proj_annot_plot")

            # MC/MC diff gene expression plots
            diff_expr_outputs(input, output, session, dataset, metacell_types, cell_type_colors, globals, ns, source_suffix = "_annot")

            mod_gene_mc_plotly_observers(input, session, source = "mc_mc_plot_annot", notification_suffix = "")

            observe({
                req(input$color_proj)
                shinyjs::toggle(id = "show_selected_metacells", condition = input$color_proj == "Cell type")
            })
        }
    )
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
        if (!is.null(input$add_to_selection) && input$add_to_selection) {
            selected_metacell_types(
                bind_rows(
                    selected_metacell_types(),
                    new_selected_annot
                ) %>% distinct(metacell, cell_type)
            )
        } else {
            selected_metacell_types(new_selected_annot %>% distinct(metacell, cell_type))
        }
    })
}
