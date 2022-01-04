# #' gene_mc UI Function
# #'
# #' @description A shiny Module.
# #'
# #' @param id,input,output,session Internal parameters for {shiny}.
# #'
# #' @noRd
# #'
# #' @importFrom shiny NS tagList
# mod_flow_mc_ui <- function(id) {
#     ns <- NS(id)
#     tagList(
#         fluidRow(
#             column(
#                 width = 7,
#                 shinydashboardPlus::box(
#                     id = ns("gene_projection"),
#                     title = "Gene projections",
#                     status = "primary",
#                     solidHeader = TRUE,
#                     collapsible = TRUE,
#                     closable = FALSE,
#                     width = 12,
#                     sidebar = shinydashboardPlus::boxSidebar(
#                         startOpen = FALSE,
#                         width = 25,
#                         id = ns("gene_projection_sidebar"),
#                         selectInput(ns("proj_stat"), label = "Statistic", choices = c("Expression" = "expression", "Enrichment" = "enrichment"), selected = "Expression", multiple = FALSE, selectize = FALSE),
#                         uiOutput(ns("set_range_ui")),
#                         uiOutput(ns("expr_range_ui")),
#                         uiOutput(ns("enrich_range_ui")),
#                         uiOutput(ns("point_size_ui")),
#                         uiOutput(ns("stroke_ui")),
#                         uiOutput(ns("edge_distance_ui"))
#                     ),
#                     shinycssloaders::withSpinner(
#                         plotly::plotlyOutput(ns("plot_gene_proj_2d"))
#                     ),
#                     shinyWidgets::prettyRadioButtons(
#                         ns("color_proj"),
#                         label = "Color by:",
#                         choices = c("Cell type", "Gene A", "Gene B"),
#                         inline = TRUE,
#                         status = "danger",
#                         fill = TRUE
#                     )
#                 ),
#                 uiOutput(ns("vein_box"))
#             ),
#             column(
#                 width = 5,
#                 shinydashboardPlus::box(
#                     id = ns("gene_gene_box"),
#                     title = "Gene/Gene",
#                     status = "primary",
#                     solidHeader = TRUE,
#                     collapsible = TRUE,
#                     closable = FALSE,
#                     width = 12,
#                     sidebar = shinydashboardPlus::boxSidebar(
#                         startOpen = FALSE,
#                         width = 25,
#                         id = ns("gene_gene_sidebar"),
#                         uiOutput(ns("gene_gene_point_size_ui")),
#                         uiOutput(ns("gene_gene_stroke_ui"))
#                     ),
#                     shinycssloaders::withSpinner(
#                         plotly::plotlyOutput(ns("plot_gene_gene_mc"))
#                     )
#                 ),
#                 uiOutput(ns("gene_time_box_ui"))
#             )
#         )
#     )
# }


# #' gene_mc sidebar UI Function
# #'
# #' @description A shiny Module.
# #'
# #' @param id,input,output,session Internal parameters for {shiny}.
# #'
# #' @noRd
# #'
# #' @importFrom shiny NS tagList
# mod_flow_mc_sidebar_ui <- function(id) {
#     ns <- NS(id)
#     tagList(
#         list(
#             uiOutput(ns("gene_selectors")),
#             uiOutput(ns("top_correlated_select_gene1")),
#             uiOutput(ns("top_correlated_select_gene2")),
#             uiOutput(ns("genecards_buttons"))
#         )
#     )
# }

# #' gene_mc Server Function
# #'
# #' @noRd
# mod_flow_mc_server <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals) {
#     ns <- session$ns

#     # gene selectors
#     values <- reactiveValues(gene1 = default_gene1, gene2 = default_gene2)
#     server_gene_selectors(input, output, session, values, dataset, ns)

#     observeEvent(plotly::event_data("plotly_click", source = "traj_plot"), {
#         el <- plotly::event_data("plotly_click", source = "traj_plot")

#         gene <- el$customdata

#         updateSelectInput(session, "gene1", selected = gene)
#         showNotification(glue("Selected {gene} in \"Genes\" tab"))
#     })

#     observeEvent(plotly::event_data("plotly_click", source = "mc_mc_plot"), {
#         el <- plotly::event_data("plotly_click", source = "mc_mc_plot")

#         gene <- el$customdata

#         updateSelectInput(session, "gene1", selected = gene)
#         showNotification(glue("Selected {gene} in \"Genes\" tab"))
#     })

#     # Expression range
#     output$set_range_ui <- renderUI({
#         req(input$proj_stat == "expression")
#         checkboxInput(ns("set_range"), "Manual range", value = FALSE)
#     })

#     output$expr_range_ui <- renderUI({
#         req(input$proj_stat == "expression")
#         req(input$set_range)
#         shinyWidgets::numericRangeInput(ns("expr_range"), "Expression range", c(-18, -5), width = "80%", separator = " to ")
#     })

#     # Enrichment range
#     output$enrich_range_ui <- renderUI({
#         req(input$proj_stat == "enrichment")
#         shinyWidgets::numericRangeInput(ns("lfp"), "Enrichment range", c(-3, 3), width = "80%", separator = " to ")
#     })

#     # Point size selectors
#     output$point_size_ui <- renderUI({
#         numericInput(ns("point_size"), label = "Point size", value = initial_proj_point_size(dataset()), min = 0.1, max = 3, step = 0.1)
#     })

#     output$gene_gene_point_size_ui <- renderUI({
#         numericInput(ns("gene_gene_point_size"), label = "Point size", value = initial_scatters_point_size(dataset()), min = 0.05, max = 3, step = 0.1)
#     })

#     output$gene_gene_stroke_ui <- renderUI({
#         numericInput(ns("gene_gene_stroke"), label = "Stroke width", value = initial_scatters_stroke(dataset()), min = 0, max = 3, step = 0.01)
#     })

#     output$stroke_ui <- renderUI({
#         numericInput(ns("stroke"), label = "Stroke width", value = initial_proj_stroke(dataset()), min = 0, max = 3, step = 0.01)
#     })

#     # Minimal edge length selector
#     output$edge_distance_ui <- renderUI({
#         sliderInput(ns("min_edge_size"), label = "Min edge length", min = 0, max = 0.3, value = min_edge_length(dataset()), step = 0.001)
#     })

#     # Projection plots
#     output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, values, metacell_types, cell_type_colors, source = "proj_mc_plot_gene_tab", refresh_on_gene_change = TRUE)

#     output$plot_gene_gene_mc <- plotly::renderPlotly({
#         req(values$gene1)
#         req(values$gene2)
#         req(input$gene_gene_point_size)

#         p_gg <- plotly::ggplotly(plot_gg_over_mc(dataset(), values$gene1, values$gene2, metacell_types = metacell_types(), cell_type_colors = cell_type_colors(), point_size = input$gene_gene_point_size, stroke = input$gene_gene_stroke, plot_text = FALSE), tooltip = "tooltip_text", source = "gene_gene_plot") %>%
#             sanitize_for_WebGL() %>%
#             plotly::toWebGL() %>%
#             sanitize_plotly_buttons() %>%
#             plotly::hide_legend()

#         return(p_gg)
#     })

#     # Expression/Time plots
#     output$gene_time_box_ui <- renderUI({
#         req(has_time(dataset()))

#         shinydashboardPlus::box(
#             id = ns("gene_time_box"),
#             title = "Gene Expression/Time",
#             status = "primary",
#             solidHeader = TRUE,
#             collapsible = TRUE,
#             closable = FALSE,
#             width = 12,
#             shinycssloaders::withSpinner(
#                 plotly::plotlyOutput(ns("plot_gene_age_mc1"))
#             ),
#             shinycssloaders::withSpinner(
#                 plotly::plotlyOutput(ns("plot_gene_age_mc2"))
#             )
#         )
#     })


#     output$plot_gene_age_mc1 <- plotly::renderPlotly({
#         req(values$gene1)
#         req(has_time(dataset()))

#         plotly::ggplotly(plot_gene_time_over_mc(dataset(), values$gene1, metacell_types = metacell_types(), cell_type_colors = cell_type_colors()), source = "gene_time_mc_plot1", tooltip = "tooltip_text") %>%
#             plotly::hide_legend() %>%
#             sanitize_plotly_buttons()
#     })

#     output$plot_gene_age_mc2 <- plotly::renderPlotly({
#         req(values$gene2)
#         req(has_time(dataset()))

#         plotly::ggplotly(plot_gene_time_over_mc(dataset(), values$gene2, metacell_types = metacell_types(), cell_type_colors = cell_type_colors()), source = "gene_time_mc_plot2", tooltip = "tooltip_text") %>%
#             plotly::hide_legend() %>%
#             sanitize_plotly_buttons()
#     })

#     connect_gene_plots(input, output, session, ns, source = "proj_mc_plot_gene_tab")

#     output$vein_gene_foc_type_select <- renderUI({
#         req(dataset())
#         selectInput(
#             ns("vein_gene_foc_type"),
#             "Focus on type:",
#             choices = c("All", sort(levels(get_mc_data(dataset(), "cell_type_colors")$cell_type))),
#             selected = "All"
#         )
#     })

#     output$vein_box <- renderUI({
#         req(has_network(dataset()))

#         shinydashboardPlus::box(
#             title = "Vein plot",
#             status = "primary",
#             solidHeader = TRUE,
#             collapsible = TRUE,
#             closable = FALSE,
#             width = 12,
#             shinyWidgets::prettyRadioButtons(
#                 ns("color_gene_vein"),
#                 label = "Color by:",
#                 choices = c("Gene A", "Gene B", "Cell type"),
#                 selected = "Cell type",
#                 inline = TRUE,
#                 status = "danger",
#                 fill = TRUE
#             ),
#             div(
#                 id = ns("vein_gene_foc_type_help"),
#                 uiOutput(ns("vein_gene_foc_type_select"))
#             ),
#             div(
#                 id = ns("vein_box"),
#                 shinycssloaders::withSpinner(
#                     plotOutput(ns("gene_vein_plot"), height = 600)
#                 )
#             )
#         )
#     })

#     output$gene_vein_plot <- renderPlot({
#         req(input$color_gene_vein)
#         req(input$vein_gene_foc_type)

#         if (input$color_gene_vein == "Gene A") {
#             req(values$gene1)
#             gene <- values$gene1
#         } else if (input$color_gene_vein == "Gene B") {
#             req(values$gene2)
#             gene <- values$gene2
#         } else {
#             gene <- NULL
#         }

#         cell_type_colors <- get_mc_data(dataset(), "cell_type_colors")

#         color_order <- cell_type_colors$color
#         if (input$vein_gene_foc_type == "All") {
#             foc_type <- NULL
#             vein_rm_cell_types <- get_mc_config(dataset(), "vein_rm_cell_types")
#             if (!is.null(vein_rm_cell_types)) {
#                 color_order <- cell_type_colors %>%
#                     filter(!(cell_type %in% vein_rm_cell_types)) %>%
#                     pull(color)
#             }
#         } else {
#             foc_type <- input$vein_gene_foc_type
#         }

#         plot_vein(dataset(), gene = gene, foc_type = foc_type, color_order = color_order, metacell_types = metacell_types())
#     })

#     # Output priorities
#     outputOptions(output, "plot_gene_proj_2d", priority = 6)
#     outputOptions(output, "plot_gene_age_mc1", priority = 5)
#     outputOptions(output, "plot_gene_age_mc2", priority = 4)
#     outputOptions(output, "top_correlated_select_gene1", priority = 4)
#     outputOptions(output, "top_correlated_select_gene2", priority = 3)
#     outputOptions(output, "gene_vein_plot", priority = 0)
# }
