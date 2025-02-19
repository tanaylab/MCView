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
        generic_column(
            width = 3,
            projection_box(ns, "gene_projection1", title = "2D Projection", collapsed_accordion = FALSE, show_legend = FALSE, 
                           color_choices = c("Flow in", "Flow out", "Cell type"), plotly_height = "30vh", height = "30vh")
        ),
        generic_column(
            width = 3,
            projection_box(ns, "gene_projection2", title = "2D Projection", collapsed_accordion = FALSE, show_legend = FALSE, 
                           color_choices = c("Flow in", "Flow out", "Cell type"), plotly_height = "30vh", height = "30vh")
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

            projection_selectors(ns, dataset, output, input, gene_modules, globals, session)
            projection_selectors(ns, dataset, output, input, gene_modules, globals, session)

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
            output$plot_gene_proj_2d <- render_2d_plotly_temporal(input, output, session, dataset, data, time_bin = '8', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")
            output$plot_gene_proj_2d <- render_2d_plotly_temporal(input, output, session, dataset, data, time_bin = '9', metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot")

        }
    )
}


render_2d_plotly_temporal <- function(input, output, session, dataset, data, time_bin, metacell_types, metacell_names, cell_type_colors, gene_modules, globals, source, buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"), dragmode = NULL, refresh_on_gene_change = FALSE, atlas = FALSE, query_types = NULL, group = NULL, groupA = NULL, groupB = NULL, selected_metacell_types = NULL, selected_cell_types = NULL) {
    plotly::renderPlotly({

        req(input$point_size)
        req(input$min_edge_size)
        req(input$color_proj)
        req(metacell_types())
        req(cell_type_colors())

        if (atlas) {
            mc2d <- NULL
        } else {
            mc2d <- globals$mc2d %||% get_mc_data(dataset, "mc2d")
        }

        proj_opt <- input$color_proj

        if(proj_opt == 'Cell type'){
            fig <- mc2d_plot_metadata_ggp(
                dataset(),
                "Cell type",
                point_size = input$point_size,
                min_d = input$min_edge_size,
                metacell_types = metacell_types(),
                atlas = atlas,
                metadata = metacell_types() %>% rename(`Cell type` = cell_type),
                colors = get_cell_type_colors(dataset, cell_type_colors = cell_type_colors(), atlas = atlas),
                graph_name = input$graph_name,
                mc2d = mc2d,
                selected_cell_types = selected_cell_types)

        }else if(proj_opt == 'Flow in'){

            req(metacell_names())
            smcs = metacell_names()

            req(input$display_select %in% metacell_names())
            smc = input$display_select
            flow_to = data$flow_to
            source_smcs = flow_to[flow_to$smc2 == smc & flow_to$time_bin == time_bin,]

            colfunc <- colorRampPalette(c("white", "Red"))
            colfunc(100)

            flow_color = rep('#ffffff', length(smcs))
            names(flow_color) = smcs 
            flow_color[smc] = '#000000'
            flow_color[source_smcs$smc1] = colfunc(100)[(round(source_smcs$f_norm,2)*100) + 1]

            flow_stroke = rep(initial_scatters_stroke(dataset), length(smcs))
            names(flow_stroke) = smcs
            flow_stroke[smc] = initial_scatters_stroke(dataset)*3

            fig <- mc2d_plot_metadata_ggp(
                dataset(),
                "metacell",
                point_size = input$point_size,
                min_d = input$min_edge_size,
                metacell_types = metacell_types(),
                atlas = atlas,
                metadata = metacell_types() %>% rename(`Cell type` = cell_type),
                colors = flow_color,
                stroke = flow_stroke,
                graph_name = input$graph_name,
                mc2d = mc2d,
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

        if (!is.null(input$legend_orientation)) {
            if (input$legend_orientation == "Horizontal") {
                orientation <- "h"
            } else if (input$legend_orientation == "Vertical") {
                orientation <- "v"
            }
            fig <- fig %>% plotly::layout(legend = list(orientation = orientation))
        }

        if (!is.null(input$show_legend_projection) && !input$show_legend_projection) {
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
