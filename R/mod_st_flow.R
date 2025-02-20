#' st_flow UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_st_flow_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_column(
                width = 3,
                projection_box_id(ns, id="flow_proj8", title = "tb8", collapsed_accordion = FALSE, show_legend = FALSE, 
                            color_choices = c("Cell type", "Flow in", "Flow out"), plotly_height = "30vh", height = "30vh")
            ),
            generic_column(
                width = 3,
                projection_box_id(ns, id="flow_proj9", title = "tb9", collapsed_accordion = FALSE, show_legend = FALSE, 
                            color_choices = c("Cell type", "Flow in", "Flow out"), plotly_height = "30vh", height = "30vh")
            ),
            generic_column(
                width = 3,
                projection_box_id(ns, id="flow_proj10", title = "tb10", collapsed_accordion = FALSE, show_legend = FALSE, 
                            color_choices = c("Cell type", "Flow in", "Flow out"), plotly_height = "30vh", height = "30vh")
            ),
            generic_column(
                width = 3,
                projection_box_id(ns, id="flow_proj11", title = "tb11", collapsed_accordion = FALSE, show_legend = FALSE, 
                            color_choices = c("Cell type", "Flow in", "Flow out"), plotly_height = "30vh", height = "30vh")
            )
        ),
        fluidRow(
            generic_column(
                width = 3,
                projection_box_id(ns, id="flow_proj12", title = "tb12", collapsed_accordion = FALSE, show_legend = FALSE, 
                            color_choices = c("Cell type", "Flow in", "Flow out"), plotly_height = "30vh", height = "30vh")
            ),
            generic_column(
                width = 3,
                projection_box_id(ns, id="flow_proj13", title = "tb13", collapsed_accordion = FALSE, show_legend = FALSE, 
                            color_choices = c("Cell type", "Flow in", "Flow out"), plotly_height = "30vh", height = "30vh")
            ),
            generic_column(
                width = 3,
                projection_box_id(ns, id="flow_proj14", title = "tb14", collapsed_accordion = FALSE, show_legend = FALSE, 
                            color_choices = c("Cell type", "Flow in", "Flow out"), plotly_height = "30vh", height = "30vh")
            ),
            generic_column(
                width = 3,
                projection_box_id(ns, id="flow_proj15", title = "tb15", collapsed_accordion = FALSE, show_legend = FALSE, 
                            color_choices = c("Cell type", "Flow in", "Flow out"), plotly_height = "30vh", height = "30vh")
            )
        )
    )
}


#' st_flow sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_st_flow_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            shinyWidgets::radioGroupButtons(
                    inputId = ns("mode"),
                    label = "Select display:",
                    choices = c(
                        "SMCs",
                        "Types"
                    ),
                    selected = "SMCs",
                    justified = TRUE
                ),

            uiOutput(ns("display_select")),

            fileInput(ns("load_projection"),
                label = NULL,
                buttonLabel = "Load 2D layout",
                multiple = FALSE,
                accept =
                    c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        "text/tab-separated-values",
                        ".csv",
                        ".tsv"
                    )
            )
        )
    )
}


#' st_flow Server Function
#'
#' @noRd
mod_st_flow_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            metacell_names <- metacell_names_reactive(dataset)
            metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)
            display_selectors(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors)

            projection_selectors_id(ns, id="flow_proj8", dataset, output, input, gene_modules, globals, session)
            projection_selectors_id(ns, id="flow_proj9", dataset, output, input, gene_modules, globals, session)
            projection_selectors_id(ns, id="flow_proj10", dataset, output, input, gene_modules, globals, session)
            projection_selectors_id(ns, id="flow_proj11", dataset, output, input, gene_modules, globals, session)
            projection_selectors_id(ns, id="flow_proj12", dataset, output, input, gene_modules, globals, session)
            projection_selectors_id(ns, id="flow_proj13", dataset, output, input, gene_modules, globals, session)
            projection_selectors_id(ns, id="flow_proj14", dataset, output, input, gene_modules, globals, session)
            projection_selectors_id(ns, id="flow_proj15", dataset, output, input, gene_modules, globals, session)

            data <- get_mc_data(dataset(), "spatial_flow_data")

            observeEvent(input$load_graph, {
                req(input$load_graph)
                req(input$load_graph$datapath)
                mc2d <- layout_and_graph_to_mc2d(mc2d_to_df(globals$mc2d), input$load_graph$datapath, metacells = get_metacell_ids(project, dataset()), warn_function = function(msg) {
                    showNotification(msg, type = "error")
                }, error_function = function(msg) {
                    showNotification(msg, type = "error")
                })
                req(mc2d)
                globals$mc2d <- mc2d
            })
            # Projection plots
            output$flow_proj8_plot_gene_proj_2d <- render_2d_plotly_temporal_id('flow_proj8', input, output, session, dataset, data, time_bin = '8', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
            output$flow_proj9_plot_gene_proj_2d <- render_2d_plotly_temporal_id('flow_proj9', input, output, session, dataset, data, time_bin = '9', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
            output$flow_proj10_plot_gene_proj_2d <- render_2d_plotly_temporal_id('flow_proj10', input, output, session, dataset, data, time_bin = '10', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
            output$flow_proj11_plot_gene_proj_2d <- render_2d_plotly_temporal_id('flow_proj11', input, output, session, dataset, data, time_bin = '11', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
            output$flow_proj12_plot_gene_proj_2d <- render_2d_plotly_temporal_id('flow_proj12', input, output, session, dataset, data, time_bin = '12', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
            output$flow_proj13_plot_gene_proj_2d <- render_2d_plotly_temporal_id('flow_proj13', input, output, session, dataset, data, time_bin = '13', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
            output$flow_proj14_plot_gene_proj_2d <- render_2d_plotly_temporal_id('flow_proj14', input, output, session, dataset, data, time_bin = '14', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
            output$flow_proj15_plot_gene_proj_2d <- render_2d_plotly_temporal_id('flow_proj15', input, output, session, dataset, data, time_bin = '15', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
        }
    )
}

render_2d_plotly_temporal_id <- function(id="flow_proj", input, output, session, dataset, data, time_bin, metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source,  buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = NULL, selected_metacell_types = NULL, selected_cell_types = NULL) {
    plotly::renderPlotly({

        req(input[[glue("{id}_point_size")]])
        req(input[[glue("{id}_min_edge_size")]])
        req(input[[glue("{id}_color_proj")]])
        req(metacell_types())
        req(cell_type_colors())

        proj_opt <- input[[glue("{id}_color_proj")]]

        if(proj_opt == 'Cell type'){
            fig <- mc2d_plot_metadata_ggp(
                dataset(),
                "Cell type",
                point_size = input[[glue("{id}_point_size")]],
                min_d = input[[glue("{id}_min_edge_size")]],
                metacell_types = metacell_types(),
                metadata = metacell_types() %>% rename(`Cell type` = cell_type),
                colors = get_cell_type_colors(dataset, cell_type_colors = cell_type_colors()),
                graph_name = input[[glue("{id}_graph_name")]],
                selected_cell_types = selected_cell_types)

        } else if(proj_opt == 'Flow in' | proj_opt == 'Flow out') {

            req(metacell_names())
            smcs = metacell_names()

            req(input$display_select %in% metacell_names())
            smc = input$display_select

            flow_stroke = rep(initial_scatters_stroke(dataset), length(smcs))
            names(flow_stroke) = smcs
            flow_stroke[smc] = initial_scatters_stroke(dataset)*3

            flow_color = rep('#ffffff', length(smcs))
            names(flow_color) = smcs 
            flow_color[smc] = '#000000'

            if(proj_opt == 'Flow in'){
                colgrad <- colorRampPalette(c("white", "Red"))
                flow_to = data$flow_to
                source_smcs = flow_to[flow_to$smc2 == smc & flow_to$time_bin == time_bin,]
                flow_color[source_smcs$smc1] = colgrad(100)[(round(source_smcs$f_norm,2)*100) + 1]
            } else {
                colgrad <- colorRampPalette(c("white", "Blue"))
                flow_from = data$flow_from
                target_smcs = flow_from[flow_from$smc1 == smc & flow_from$time_bin == time_bin,]
                flow_color[target_smcs$smc2] = colgrad(100)[(round(target_smcs$f_norm,2)*100) + 1]
            }
            
            fig <- mc2d_plot_metadata_ggp(
                dataset(),
                "metacell",
                point_size = input[[glue("{id}_point_size")]],
                min_d = input[[glue("{id}_min_edge_size")]],
                metacell_types = metacell_types(),
                colors = flow_color,
                stroke = flow_stroke,
                graph_name = input[[glue("{id}_graph_name")]],
                selected_cell_types = selected_cell_types
            )
        }
        
        fig <- fig %>% plotly::event_register("plotly_restyle")
        fig$x$source <- source

        if (!is.null(dragmode)) {
            fig <- fig %>% plotly::layout(dragmode = dragmode)
        } else if (!is.null(input$mode) && input$mode %in% c("Groups", "Group")) {
            fig <- fig %>% plotly::layout(dragmode = "select")
            buttons <- buttons[!(buttons %in% c("select2d", "lasso2d"))]
        }

        fig <- fig %>% sanitize_plotly_buttons(buttons = buttons)

        if (!is.null(input[[glue("{id}_legend_orientation")]])) {
            if (input[[glue("{id}_legend_orientation")]] == "Horizontal") {
                orientation <- "h"
            } else if (input[[glue("{id}_legend_orientation")]] == "Vertical") {
                orientation <- "v"
            }
            fig <- fig %>% plotly::layout(legend = list(orientation = orientation))
        }

        if (!is.null(input[[glue("{id}_show_legend_projection")]]) && !input[[glue("{id}_show_legend_projection")]]) {
            fig <- plotly::hide_legend(fig)
        }

        fig <- fig %>%
            sanitize_for_WebGL()
        fig <- fig %>%
            plotly::toWebGL()
        fig <- fig %>%
            rm_plotly_grid()

        return(fig)
    })
}

display_selectors <- function(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors) {
    output$display_select <- renderUI({
        req(dataset())
        req(input$mode)
        if (input$mode == "SMCs") {
            req(metacell_colors())
            req(metacell_names())
            cell_types_hex <- col2hex(metacell_colors())
            # add 'similar' annotation
            md <- get_mc_data(dataset(), "metadata")
            if (!is.null(md) && has_name(md, "similar")) {
                choices <- metacell_names()
                names(choices) <- ifelse(md$similar == "dissimilar", paste0(metacell_names(), " (dissimilar)"), metacell_names())
            } else {
                choices <- metacell_names()
            }
            shinyWidgets::pickerInput(ns("display_select"), "Smc",
                choices = choices,
                selected = 'M267.82',
                multiple = FALSE,
                options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                choicesOpt = list(
                    style = paste0("color: ", cell_types_hex, ";")
                )
            )
        }else if(input$mode == "Types"){
            req(cell_type_colors())
            req(metacell_types())
            types_df <- cell_type_colors() %>% filter(cell_type %in% metacell_types()$cell_type)
            cell_types_hex <- col2hex(types_df$color)
            cell_types <- types_df$cell_type
            tagList(
                shinyWidgets::pickerInput(ns("display_select"), "Cell Type",
                    choices = cell_types,
                    selected = 'Epiblast',
                    multiple = FALSE,
                    options = shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
                    choicesOpt = list(
                        style = paste0("color: ", cell_types_hex, ";")
                    )
                ),
                shinyWidgets::switchInput(
                inputId = ns("Expand_Smcs"),
                label = "Show SMCs",
                size = "sm",
                value = FALSE
                )
            )
        }
    })}

projection_box_id <- function(ns,
                           id,
                           ...,
                           color_choices = c("Cell type", "Gene", "Gene module", "Metadata"),
                           title = "2D Projection",
                           height = NULL,
                           plotly_height = "400px",
                           additional_elements = NULL,
                           collapsed_accordion = TRUE,
                           legend_orientation = "Vertical",
                           show_legend = TRUE) {
    generic_box(
        id = ns(id),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        height = height,
        sidebar = shinydashboardPlus::boxSidebar(
            startOpen = FALSE,
            width = 80,
            id = ns(glue("{id}_sidebar")),
            uiOutput(ns(glue("{id}_graph_select"))),
            uiOutput(ns(glue("{id}_proj_stat"))),
            uiOutput(ns(glue("{id}_set_range"))),
            uiOutput(ns(glue("{id}_expr_range"))),
            uiOutput(ns(glue("{id}_enrich_range"))),
            uiOutput(ns(glue("{id}_point_size"))),
            uiOutput(ns(glue("{id}_stroke"))),
            uiOutput(ns(glue("{id}_edge_distance"))),
            shinyWidgets::prettyRadioButtons(
                ns(glue("{id}_legend_orientation")),
                label = "Legend orientation:",
                choices = c("Vertical", "Horizontal"),
                selected = legend_orientation,
                inline = TRUE,
                status = "danger",
                fill = TRUE
            )
        ),
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns(glue("{id}_plot_gene_proj_2d")), height = plotly_height)
        ),
        shinydashboardPlus::accordion(
            id = ns(glue("{id}_proj_accordion")),
            shinydashboardPlus::accordionItem(
                title = "Modify colors",
                collapsed = collapsed_accordion,
                shinyWidgets::prettyRadioButtons(
                    ns(glue("{id}_color_proj")),
                    label = "Color by:",
                    choices = color_choices,
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                ),
                ...,
                shinyWidgets::virtualSelectInput(
                    ns(glue("{id}_color_proj_gene")),
                    "Gene:",
                    choices = c(),
                    multiple = FALSE,
                    search = TRUE,
                    dropboxWrapper = "body",
                    markSearchResults = TRUE,
                    searchByStartsWith = TRUE
                ),
                shinyWidgets::virtualSelectInput(
                    ns(glue("{id}_color_proj_metadata")),
                    "Metadata:",
                    choices = c(),
                    multiple = FALSE,
                    search = TRUE,
                    dropboxWrapper = "body"
                ),
                shinyWidgets::virtualSelectInput(
                    ns(glue("{id}_color_proj_gene_module")),
                    "Gene module:",
                    choices = c(),
                    multiple = FALSE,
                    search = TRUE,
                    dropboxWrapper = "body"
                ),
                shinyWidgets::prettyRadioButtons(
                    ns(glue("{id}_scatter_axis_proj")),
                    label = "Axis:",
                    choices = c("x", "y"),
                    inline = TRUE,
                    status = "danger",
                    fill = TRUE
                ),
                checkboxInput(ns(glue("{id}_show_legend_projection")), "Show legend", value = show_legend)
            )
        )
    )
}

projection_selectors_id <- function(ns, id, dataset, output, input, gene_modules, globals, session, weight = 1, atlas = FALSE) {
    observe({
        shinyWidgets::updateVirtualSelect(
            session = session,
            inputId = glue("{id}_color_proj_gene"),
            choices = gene_names_label(dataset(), atlas = atlas),
            selected = default_gene1
        )

        shinyWidgets::updateVirtualSelect(
            session = session,
            inputId = glue("{id}_color_proj_metadata"),
            choices = c("Clipboard", dataset_metadata_fields(dataset(), atlas = atlas)),
            selected = dataset_metadata_fields(dataset(), atlas = atlas)[1]
        )

        shinyWidgets::updateVirtualSelect(
            session = session,
            inputId = glue("{id}_color_proj_gene_module"),
            choices = levels(gene_modules()$module),
            selected = NULL
        )
    })

    picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE)

    observe({
        req(input[[glue("{id}_color_proj")]])
        shinyjs::toggle(id = glue("{id}_color_proj_gene"), condition = input[[glue("{id}_color_proj")]] == "Gene")
        shinyjs::toggle(id = glue("{id}_color_proj_metadata"), condition = input[[glue("{id}_color_proj")]] == "Metadata")
        shinyjs::toggle(id = glue("{id}_color_proj_gene_module"), condition = input[[glue("{id}_color_proj")]] == "Gene module")
        shinyjs::toggle(id = glue("{id}_scatter_axis_proj"), condition = input[[glue("{id}_color_proj")]] == "Scatter Axis")
    })

    output[[glue("{id}_proj_stat")]] <- renderUI({
        req(input[[glue("{id}_color_proj")]] == "Gene" ||
            input[[glue("{id}_color_proj")]] == "Gene A" ||
            input[[glue("{id}_color_proj")]] == "Gene B" ||
            input[[glue("{id}_color_proj")]] == "Gene module")
        selectInput(ns(glue("{id}_proj_stat")),
            label = "Statistic",
            choices = c("Expression" = "expression", "Enrichment" = "enrichment"),
            selected = "Expression",
            multiple = FALSE,
            selectize = FALSE
        )
    })

    output[[glue("{id}_graph_select")]] <- renderUI({
        choices <- c("metacell")
        graphs <- get_mc_data(dataset(), "metacell_graphs")
        if (!is.null(graphs)) {
            choices <- c(choices, names(graphs))
        }
        selectInput(ns(glue("{id}_graph_name")),
            label = "Graph",
            choices = choices,
            selected = "metacell",
            multiple = FALSE,
            selectize = FALSE
        )
    })

    output[[glue("{id}_set_range")]] <- renderUI({
        req(input[[glue("{id}_color_proj")]] == "Gene" ||
            input[[glue("{id}_color_proj")]] == "Gene A" ||
            input[[glue("{id}_color_proj")]] == "Gene B" ||
            input[[glue("{id}_color_proj")]] == "Gene module")
        req(input[[glue("{id}_proj_stat")]] == "expression")
        checkboxInput(ns(glue("{id}_set_range")), "Manual range", value = FALSE)
    })

    output[[glue("{id}_expr_range")]] <- renderUI({
        req(input[[glue("{id}_color_proj")]] == "Gene" ||
            input[[glue("{id}_color_proj")]] == "Gene A" ||
            input[[glue("{id}_color_proj")]] == "Gene B" ||
            input[[glue("{id}_color_proj")]] == "Gene module")
        req(input[[glue("{id}_proj_stat")]] == "expression")
        req(input[[glue("{id}_set_range")]])
        shinyWidgets::numericRangeInput(ns(glue("{id}_expr_range")),
            "Expression range",
            c(-18, -5),
            width = "80%",
            separator = " to "
        )
    })

    output[[glue("{id}_enrich_range")]] <- renderUI({
        req(input[[glue("{id}_color_proj")]] == "Gene" ||
            input[[glue("{id}_color_proj")]] == "Gene A" ||
            input[[glue("{id}_color_proj")]] == "Gene B" ||
            input[[glue("{id}_color_proj")]] == "Gene module")
        req(input[[glue("{id}_proj_stat")]] == "enrichment")
        shinyWidgets::numericRangeInput(ns(glue("{id}_lfp")),
            "Enrichment range",
            c(-3, 3),
            width = "80%",
            separator = " to "
        )
    })

    output[[glue("{id}_point_size")]] <- renderUI({
        req(globals$screen_height)
        req(globals$screen_width)
        req(dataset())
        numericInput(ns(glue("{id}_point_size")),
            label = "Point size",
            value = initial_proj_point_size(dataset(), globals$screen_width, globals$screen_height, weight = weight, atlas = atlas),
            min = 0.1,
            max = 3,
            step = 0.1
        )
    })

    output[[glue("{id}_stroke")]] <- renderUI({
        numericInput(ns(glue("{id}_stroke")),
            label = "Stroke width",
            value = initial_proj_stroke(dataset()),
            min = 0,
            max = 3,
            step = 0.01
        )
    })

    output[[glue("{id}_edge_distance")]] <- renderUI({
        graph <- input[[glue("{id}_graph_name")]]
        if (is.null(graph) || graph == "metacell") {
            sliderInput(ns(glue("{id}_min_edge_size")),
                label = "Min edge length",
                min = 0,
                max = 0.3,
                value = min_edge_length(dataset()),
                step = 0.001
            )
        } else {
            graph <- get_mc_data(dataset(), "metacell_graphs")[[graph]]
            sliderInput(ns(glue("{id}_min_edge_size")),
                label = "Min weight",
                min = min(graph$weight, na.rm = TRUE),
                max = max(graph$weight, na.rm = TRUE),
                value = median(graph$weight),
                step = (max(graph$weight, na.rm = TRUE) - min(graph$weight, na.rm = TRUE)) / 50
            )
        }
    })
}