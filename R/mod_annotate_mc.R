#' annotate_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_annotate_mc_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                width = 8,
                shinydashboardPlus::box(
                    id = ns("gene_projection"),
                    title = "Gene projections",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 25,
                        id = ns("gene_projection_sidebar"),
                        selectInput(ns("proj_stat"), label = "Statistic", choices = c("Expression" = "expression", "Enrichment" = "enrichment"), selected = "Expression", multiple = FALSE, selectize = FALSE),
                        uiOutput(ns("set_range_ui")),
                        uiOutput(ns("expr_range_ui")),
                        uiOutput(ns("enrich_range_ui")),
                        uiOutput(ns("point_size_ui")),
                        uiOutput(ns("edge_distance_ui"))
                    ),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_gene_proj_2d"))
                    ),
                    shinyWidgets::prettyRadioButtons(
                        ns("color_proj"),
                        label = "Color by:",
                        choices = c("Cell type", "Gene A", "Gene B"),
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    )
                ),
                fluidRow(
                    uiOutput(ns("time_box_ui_column")),
                    column(
                        width = 6,
                        offset = 0,
                        shinydashboardPlus::box(
                            id = ns("gene_gene_box"),
                            title = "Gene/Gene",
                            status = "primary",
                            solidHeader = TRUE,
                            collapsible = TRUE,
                            closable = FALSE,
                            width = 12,
                            shinycssloaders::withSpinner(
                                plotly::plotlyOutput(ns("plot_gene_gene_mc"))
                            )
                        )
                    ),
                    column(
                        width = 6,
                        offset = 0,
                        shinydashboardPlus::box(
                            title = "Metacell/Metacell",
                            status = "primary",
                            solidHeader = TRUE,
                            collapsible = TRUE,
                            closable = FALSE,
                            width = 12,
                            uiOutput(ns("metacell1_select")),
                            uiOutput(ns("metacell2_select")),
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
                shinydashboardPlus::box(
                    id = ns("metacell_typesation"),
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
                    shinyWidgets::prettySwitch(inputId = ns("show_all_annotation"), value = FALSE, label = "All metacells"),
                    uiOutput(ns("annotation_box")),
                    uiOutput(ns("update_all_selectors")),
                    shinycssloaders::withSpinner(
                        DT::dataTableOutput(ns("mc_type_table"))
                    )
                )
            ),
            column(
                width = 4,
                shinydashboardPlus::box(
                    id = ns("cell_type_colorsation"),
                    title = "Cell Types",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    splitLayout(
                        fileInput(ns("cell_type_colors_fn"),
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
                        actionButton(ns("reset_cell_type_colors"), "Reset", style = "align-items: center;"),
                        downloadButton(ns("cell_type_colors_download"), "Export", style = "align-items: center;"),
                        actionButton(ns("delete_cell_type_colorsation"), "Delete"),
                        actionButton(ns("add_cell_type_colorsation"), "Add")
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


#' annotate_mc sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_annotate_mc_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("gene_selectors")),
        tags$hr(),
        uiOutput(ns("top_correlated_select_gene1")),
        uiOutput(ns("top_correlated_select_gene2")),
        tags$hr(),
        uiOutput(ns("genecards_buttons"))
    )
}


#' annotate_mc Server Function
#'
#' @noRd
mod_annotate_mc_server <- function(input, output, session, dataset, metacell_types, cell_type_colors) {
    ns <- session$ns

    # gene selectors
    values <- reactiveValues(gene1 = default_gene1, gene2 = default_gene2)
    server_gene_selectors(input, output, session, values, dataset, ns)

    observe({
        initial_cell_type_colors <- get_mc_data(dataset(), "cell_type_colors")
        initial_metacell_types <- get_mc_data(dataset(), "metacell_types")

        # remove metacell color column if exists
        initial_metacell_types$mc_col <- NULL

        # add cell_type_id and cell type color from initial cell type annotation
        initial_metacell_types <- initial_metacell_types %>%
            left_join(initial_cell_type_colors %>% select(cell_type, cell_type_id, mc_col = color), by = "cell_type")

        metacell_types(initial_metacell_types)
        cell_type_colors(initial_cell_type_colors)
    })

    observe({
        req(input$metacell_types_fn)
        new_metacell_types <- tgutil::fread(input$metacell_types_fn$datapath, colClasses = c("cell_type_id" = "character", "cell_type" = "character", "metacell" = "character")) %>% as_tibble()

        cur_metacell_types <- metacell_types()
        new_metacell_types <- cur_metacell_types %>%
            select(-any_of(c("cell_type", "cell_type_id"))) %>%
            left_join(new_metacell_types %>% select(metacell, cell_type, cell_type_id), by = "metacell") %>%
            mutate(cell_type = as.character(forcats::fct_explicit_na(factor(cell_type))))

        new_metacell_types <- sanitize_metacell_types(new_metacell_types, cell_type_colors(), dataset())

        metacell_types(new_metacell_types)
    })

    observe({
        req(input$cell_type_colors_fn)
        new_cell_type_colors <- tgutil::fread(input$cell_type_colors_fn$datapath, colClasses = c("cell_type_id" = "character", "cell_type" = "character", "color" = "character")) %>% as_tibble()
        if ("order" %in% colnames(new_cell_type_colors)) {
            new_cell_type_colors <- new_cell_type_colors %>% arrange(order)
        }
        if (!rlang::has_name(new_cell_type_colors, "cell_type_id")) {
            new_cell_type_colors <- new_cell_type_colors %>% mutate(cell_type_id = as.character(1:n()))
        }

        # Metacell annotations that are now invalid would get cell_type of NA.
        cur_metacell_types <- metacell_types()
        new_metacell_types <- cur_metacell_types %>%
            select(-cell_type, -mc_col) %>%
            left_join(
                new_cell_type_colors %>% select(cell_type_id, cell_type, mc_col = color),
                by = "cell_type_id"
            )


        cell_type_colors(new_cell_type_colors)
        metacell_types(new_metacell_types)
    })

    output$metacell_types_download <- downloadHandler(
        filename = function() {
            paste("metacell_types-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            fwrite(
                metacell_types() %>%
                    select(metacell, cell_type_id, cell_type, top1_gene, top1_lfp, top2_gene, top2_lfp),
                file
            )
        }
    )

    output$cell_type_colors_download <- downloadHandler(
        filename = function() {
            paste("cell_type_colors-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            fwrite(
                cell_type_colors() %>%
                    select(cell_type_id, cell_type, color),
                file
            )
        }
    )

    selected_metacell_types <- reactiveVal(tibble(metacell = character(), cell_type_id = character(), cell_type = character()))
    to_show <- reactiveVal()

    observeEvent(input$reset_metacell_types, {
        metacell_types(get_mc_data(dataset(), "metacell_types"))
        selected_metacell_types(tibble(metacell = character(), cell_type_id = character(), cell_type = character()))
        to_show(NULL)
    })

    observeEvent(input$reset_cell_type_colors, {
        cell_type_colors(get_mc_data(dataset(), "cell_type_colors"))
    })


    output$annotation_box <- renderUI({
        if (!input$show_all_annotation) {
            if (nrow(selected_metacell_types()) == 0) {
                print("Please select metacells")
            } else {
                list(
                    actionButton(ns("update_annotation"), "Apply"),
                    actionButton(ns("reset_annotation"), "Reset Selection"),
                    shinyWidgets::radioGroupButtons(
                        inputId = ns("update_option"),
                        label = "",
                        choices = c(
                            "Change all",
                            "Change one by one"
                        ),
                        justified = TRUE
                    )
                )
            }
        } else {
            list(
                uiOutput(ns("cell_type_select")),
                actionButton(ns("update_annotation"), "Apply")
            )
        }
    })

    output$update_all_selectors <- renderUI({
        req(input$update_option)
        if (!input$show_all_annotation && input$update_option == "Change all") {
            selectizeInput(ns("selected_cell_type_update_all"), "Cell type", choices = c("(Missing)", cell_type_colors() %>% pull(cell_type) %>% as.character() %>% unique() %>% sort()), multiple = FALSE, selected = "(Missing)")
        }
    })

    output$cell_type_select <- renderUI({
        selectizeInput(ns("selected_cell_type"), "Show cell type", choices = c("All", "(Missing)", cell_type_colors() %>% pull(cell_type) %>% as.character() %>% unique() %>% sort()), multiple = FALSE, selected = "(Missing)")
    })

    observeEvent(input$update_annotation, {
        new_metacell_types <- metacell_types()
        changed <- FALSE

        if (!input$show_all_annotation && !is.null(input$update_option) && input$update_option == "Change all") {
            req(input$selected_cell_type_update_all)
            metacells <- selected_metacell_types() %>% pull(metacell)
            new_metacell_types <- new_metacell_types %>% mutate(
                cell_type = ifelse(metacell %in% metacells, input$selected_cell_type_update_all, cell_type),
                cell_type_id = ifelse(metacell %in% metacells, cell_type_to_cell_type_id(input$selected_cell_type_update_all, cell_type_colors()), cell_type_id)
            )
            new_selected_annot <- selected_metacell_types()
            new_selected_annot <- new_selected_annot %>% mutate(
                cell_type = ifelse(metacell %in% metacells, input$selected_cell_type_update_all, cell_type),
                cell_type_id = ifelse(metacell %in% metacells, cell_type_to_cell_type_id(input$selected_cell_type_update_all, cell_type_colors()), cell_type_id)
            )
            selected_metacell_types(new_selected_annot)
            changed <- TRUE
        } else {
            for (i in 1:nrow(new_metacell_types)) {
                cur_input <- input[[glue("select_type_{metacell}", metacell = new_metacell_types$metacell[i])]]
                if (!is.null(cur_input)) {
                    if (cur_input != new_metacell_types[i, ]$metacell) {
                        new_metacell_types[i, ]$cell_type <- cur_input
                        new_metacell_types[i, ]$cell_type_id <- cell_type_to_cell_type_id(cur_input, cell_type_colors())
                        changed <- TRUE
                    }
                }
            }
        }
        if (changed) {
            metacell_types(new_metacell_types)
        }
    })

    observeEvent(input$reset_annotation, {
        selected_metacell_types(tibble(metacell = character(), cell_type_id = character(), cell_type = character()))
    })


    observe({
        if (input$show_all_annotation) {
            req(metacell_types)
            req(input$selected_cell_type)
            to_show_new <- metacell_types() %>% select(metacell, cell_type)
            if (input$selected_cell_type != "All") {
                to_show_new <- to_show_new %>% filter(cell_type == input$selected_cell_type)
            }
        } else {
            to_show_new <- selected_metacell_types() %>% select(metacell, cell_type)
        }

        if (input$show_all_annotation || (!is.null(input$update_option) && input$update_option != "Change all")) {
            to_show_new <- dt_selector_column(
                to_show_new,
                "select_type",
                ns,
                c("(Missing)", sort(as.character(unique(cell_type_colors()$cell_type))))
            )
        }

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
            paging = FALSE
        ),
        callback = DT::JS("table.rows().every(function(i, tab, row) {
                        var $this = $(this.node());
                        $this.attr('id', this.data()[0]);
                        $this.addClass('shiny-input-container');
                    });
                    Shiny.unbindAll(table.table().node());
                    Shiny.bindAll(table.table().node());")
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
                    cell_type_colors()$color
                )
            ),
        server = TRUE # see https://github.com/rstudio/DT/issues/598
    )

    cell_type_table_proxy <- DT::dataTableProxy("cell_type_table")

    observeEvent(input$cell_type_table_cell_edit, {
        # fix column number to be 1 based
        new_input <- input$cell_type_table_cell_edit %>% mutate(col = col + 1)
        edited_data <- DT::editData(cell_type_colors(), new_input, "cell_type_table")

        # Should we reorder or is it annoying?
        # if (!is.numeric(edited_data$order)){
        #     edited_data <- edited_data %>% mutate(order = 1:n())
        # }

        # edited_data <- edited_data %>%
        #     arrange(order) %>%
        #     distinct(cell_type, .keep_all=TRUE) %>%
        #     mutate(order = 1:n())

        # change corresponding metacell_typeentries
        if (new_input$col == 1) {
            old_cell_type <- as.character(cell_type_colors()$cell_type[new_input$row])
            new_cell_type <- as.character(edited_data$cell_type[new_input$row])
            new_metacell_types <- metacell_types() %>% mutate(cell_type = ifelse(cell_type == old_cell_type, new_cell_type, cell_type))
            metacell_types(new_metacell_types)
        }

        cell_type_colors(edited_data)

        DT::replaceData(cell_type_table_proxy, cell_type_colors(), resetPaging = FALSE)
    })

    observeEvent(input$delete_cell_type_colorsation, {
        rows <- input$cell_type_table_rows_selected

        if (!is.null(rows) && length(rows) > 0) {
            to_delete <- cell_type_colors()[rows, ]
            cell_type_colors(cell_type_colors()[-rows, ])
            metacell_types(
                metacell_types() %>%
                    mutate(
                        cell_type = ifelse(cell_type_id %in% to_delete$cell_type_id, NA, cell_type),
                        cell_type_id = ifelse(cell_type_id %in% to_delete$cell_type_id, NA, cell_type_id)
                    )
            )
        }
    })

    observeEvent(input$add_cell_type_colorsation, {
        rows <- input$cell_type_table_rows_selected
        if (!is.null(rows) && length(rows) > 0) {
            place <- rows[1] + 1
        } else {
            place <- 1
        }

        # TODO: allow creation of more than one new row without editing (e.g. by adding a suffix to cell_type)
        new_data <- cell_type_colors() %>% arrange(order)
        new_id <- as.character(max(as.numeric(cell_type_colors()$cell_type_id)) + 1)
        new_row <- tibble(cell_type = "Cell Type", cell_type_id = new_id, color = "red", order = place)
        new_data <- bind_rows(
            new_data %>% filter(order < place),
            new_row,
            new_data %>% filter(order >= place) %>% mutate(order = order + 1)
        )
        new_data <- new_data %>%
            arrange(order) %>%
            distinct(cell_type, cell_type_id, .keep_all = TRUE) %>%
            mutate(order = 1:n())
        cell_type_colors(new_data)
    })

    output$annot_color_picker <- renderUI({
        req(input$cell_type_table_rows_selected)
        fluidRow(
            column(6, actionButton(ns("submit_new_color"), "Change color")),
            column(6, colourpicker::colourInput(ns("selected_new_color"), NULL, "black"))
        )
    })

    observeEvent(input$submit_new_color, {
        rows <- input$cell_type_table_rows_selected
        new_data <- cell_type_colors()
        new_data$color[rows] <- input$selected_new_color
        cell_type_colors(new_data)
    })


    # Select metacell when clicking on it
    observe_mc_click_event("proj_annot_plot", input, cell_type_colors, metacell_types)
    observe_mc_click_event("gene_gene_plot_annot", input, cell_type_colors, metacell_types)
    observe_mc_click_event("gene_time_mc_plot1_annot", input, cell_type_colors, metacell_types)
    observe_mc_click_event("gene_time_mc_plot2_annot", input, cell_type_colors, metacell_types)

    # Select multiple metacells
    observer_mc_select_event("proj_annot_plot", input, cell_type_colors, metacell_types, selected_metacell_types)
    observer_mc_select_event("gene_gene_plot_annot", input, cell_type_colors, metacell_types, selected_metacell_types)
    observer_mc_select_event("gene_time_mc_plot1_annot", input, cell_type_colors, metacell_types, selected_metacell_types)
    observer_mc_select_event("gene_time_mc_plot2_annot", input, cell_type_colors, metacell_types, selected_metacell_types)

    # De-select multiple metacells
    observer_mc_deselect_event("proj_annot_plot", selected_metacell_types)
    observer_mc_deselect_event("gene_gene_plot_annot", selected_metacell_types)
    observer_mc_deselect_event("gene_time_mc_plot1_annot", selected_metacell_types)
    observer_mc_deselect_event("gene_time_mc_plot2_annot", selected_metacell_types)

    # Expression range
    output$set_range_ui <- renderUI({
        req(input$proj_stat == "expression")
        checkboxInput(ns("set_range"), "Manual range", value = FALSE)
    })

    output$expr_range_ui <- renderUI({
        req(input$proj_stat == "expression")
        req(input$set_range)
        shinyWidgets::numericRangeInput(ns("expr_range"), "Expression range", c(-18, -5), width = "80%", separator = " to ")
    })

    # Enrichment range
    output$enrich_range_ui <- renderUI({
        req(input$proj_stat == "enrichment")
        shinyWidgets::numericRangeInput(ns("lfp"), "Enrichment range", c(-3, 3), width = "80%", separator = " to ")
    })

    # Point size selector
    output$point_size_ui <- renderUI({
        numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset()), min = 0.1, max = 3, step = 0.1)
    })

    # Minimal edge length selector
    output$edge_distance_ui <- renderUI({
        sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
    })
    # Projection plots
    output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, values, metacell_types, cell_type_colors, source = "proj_annot_plot", buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = "select")


    output$plot_gene_gene_mc <- plotly::renderPlotly({
        req(values$gene1)
        req(values$gene2)
        req(metacell_types())
        req(dataset())

        p_gg <- plotly::ggplotly(plot_gg_over_mc(dataset(), values$gene1, values$gene2, metacell_types = metacell_types(), cell_type_colors = cell_type_colors(), plot_text = FALSE), tooltip = "tooltip_text", source = "gene_gene_plot_annot") %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons(buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines")) %>%
            plotly::hide_legend() %>%
            plotly::layout(dragmode = "select")

        return(p_gg)
    })


    output$time_box_ui_column <- renderUI({
        req(has_time(dataset()))

        column(
            width = 6,
            offset = 0,
            uiOutput(ns("gene_time_box_ui")),
            uiOutput(ns("mc_time_box_ui"))
        )
    })

    output$gene_time_box_ui <- renderUI({
        req(has_time(dataset()))

        shinydashboardPlus::box(
            id = ns("gene_time_box"),
            title = "Gene Expression/Time",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            shinycssloaders::withSpinner(
                plotly::plotlyOutput(ns("plot_gene_age_mc1"))
            ),
            shinycssloaders::withSpinner(
                plotly::plotlyOutput(ns("plot_gene_age_mc2"))
            )
        )
    })

    output$mc_time_box_ui <- renderUI({
        req(has_time(dataset()))

        shinydashboardPlus::box(
            id = ns("mc_time_box"),
            title = "Metacell Age",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            uiOutput(ns("time_dist_metacell_select")),
            shinycssloaders::withSpinner(
                plotly::plotlyOutput(ns("mc_time_dist_plot"))
            )
        )
    })

    # Expression/Time plots
    output$plot_gene_age_mc1 <- plotly::renderPlotly({
        req(values$gene1)

        plotly::ggplotly(plot_gene_time_over_mc(dataset(), values$gene1, metacell_types = metacell_types(), cell_type_colors = cell_type_colors()), source = "gene_time_mc_plot1_annot", tooltip = "tooltip_text") %>%
            plotly::hide_legend() %>%
            sanitize_plotly_buttons(buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines")) %>%
            plotly::layout(dragmode = "select")
    })

    output$plot_gene_age_mc2 <- plotly::renderPlotly({
        req(values$gene2)

        plotly::ggplotly(plot_gene_time_over_mc(dataset(), values$gene2, metacell_types = metacell_types(), cell_type_colors = cell_type_colors()), source = "gene_time_mc_plot2_annot", tooltip = "tooltip_text") %>%
            plotly::hide_legend() %>%
            sanitize_plotly_buttons(buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines")) %>%
            plotly::layout(dragmode = "select")
    })

    connect_gene_plots(input, output, session, ns, source = "proj_annot_plot")


    # MC/MC diff gene expression plots
    metacell_names <- reactive({
        req(dataset())
        colnames(get_mc_data(dataset(), "mc_mat"))
    })

    output$metacell1_select <- renderUI({
        req(dataset())
        selectizeInput(ns("metacell1"), "Metacell A", choices = metacell_names(), multiple = FALSE, selected = config$selected_mc1, options = list(maxOptions = 1e4))
    })

    output$metacell2_select <- renderUI({
        req(dataset())
        selectizeInput(ns("metacell2"), "Metacell B", choices = metacell_names(), multiple = FALSE, selected = config$selected_mc2, options = list(maxOptions = 1e4))
    })

    mc_mc_gene_scatter_df <- reactive({
        calc_mc_mc_gene_df(dataset(), input$metacell1, input$metacell2)
    })

    output$plot_mc_mc_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, mc_mc_gene_scatter_df)

    output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, mc_mc_gene_scatter_df)

    output$time_dist_metacell_select <- renderUI({
        selectizeInput(ns("time_dist_metacell"), "Metacell", choices = metacell_names(), multiple = FALSE, selected = config$selected_mc1, options = list(maxOptions = 1e4))
    })

    output$mc_time_dist_plot <- plotly::renderPlotly({
        req(input$time_dist_metacell)
        plotly::ggplotly(plot_mc_time_dist(dataset(), input$time_dist_metacell, ylab = glue("# of cells #{input$time_dist_metacell}")), tooltip = "tooltip_text", height = 200) %>%
            sanitize_plotly_buttons()
    })
}


observe_mc_click_event <- function(source, input, cell_type_colors, metacell_types) {
    observeEvent(plotly::event_data("plotly_click", source = source), {
        el <- plotly::event_data("plotly_click", source = source)

        selected_metacell <- el$customdata

        if (input$show_all_annotation && input$selected_cell_type %in% cell_type_colors()$cell_type) {
            new_metacell_types <- metacell_types() %>% mutate(cell_type = ifelse(metacell == selected_metacell, input$selected_cell_type, cell_type))
            metacell_types(new_metacell_types)
            showNotification(glue("Added metacell #{selected_metacell} to {input$selected_cell_type}"))
        }
    })
}

observer_mc_select_event <- function(source, input, cell_type_colors, metacell_types, selected_metacell_types) {
    observeEvent(plotly::event_data("plotly_selected", source = source), {
        el <- plotly::event_data("plotly_selected", source = source)

        selected_metacells <- el$customdata

        new_selected_annot <- metacell_types() %>% filter(metacell %in% selected_metacells)
        selected_metacell_types(new_selected_annot)
    })
}

observer_mc_deselect_event <- function(source, selected_metacell_types) {
    observeEvent(plotly::event_data("plotly_deselect", source = source), {
        selected_metacell_types(tibble(metacell = character(), cell_type = character()))
    })
}
