#' gene_mc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_flow_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_column(
                width = 6,
                generic_box(
                    title = "Vein plot",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    shinycssloaders::withSpinner(
                        plotOutput(ns("gene_vein_plot"), height = 600)
                    ),
                    shinyWidgets::prettyRadioButtons(
                        ns("color_gene_vein"),
                        label = "Color by:",
                        choices = c("Gene", "Cell type"),
                        selected = "Cell type",
                        inline = TRUE,
                        status = "danger",
                        fill = TRUE
                    ),
                    shinyWidgets::virtualSelectInput(
                        ns("gene"),
                        "",
                        choices = c(),
                        multiple = FALSE,
                        search = TRUE,
                        dropboxWrapper = "body"
                    ),
                    uiOutput(ns("vein_gene_foc_type_select"))
                ),
                generic_box(
                    title = "Genes Trajectory",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    uiOutput(ns("traj_gene_selector")),
                    shinycssloaders::withSpinner(
                        plotly::plotlyOutput(ns("plot_mc_traj"))
                    ),
                    uiOutput(ns("metacell_selector_traj"))
                )
            ),
            generic_column(
                width = 6,
                generic_box(
                    title = "Metacell flow",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    # sidebar = shinydashboardPlus::boxSidebar(
                    #     startOpen = FALSE,
                    #     width = 100,
                    #     id = ns("flow_box_sidebar")
                    # ),
                    shinycssloaders::withSpinner(
                        plotOutput(ns("plot_metacell_flow"), height = 600)
                    ),
                    uiOutput(ns("metacell_selector"))
                )
            )
        )
    )
}


#' gene_mc sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_flow_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list()
    )
}

#' gene_mc Server Function
#'
#' @noRd
mod_flow_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            metacell_names <- metacell_names_reactive(dataset)
            metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)

            output$metacell_selector <- colored_metacell_selector(dataset, ns, "selected_metacell", "Metacell", metacell_colors, metacell_names)
            output$metacell_selector_traj <- colored_metacell_selector(dataset, ns, "selected_metacell_traj", "Metacell", metacell_colors, metacell_names)

            picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE)

            observe({
                shinyWidgets::updateVirtualSelect(
                    session = session,
                    inputId = "gene",
                    choices = gene_names(dataset()),
                    selected = default_gene1,
                    label = "Gene"
                )
            })

            output$vein_gene_foc_type_select <- renderUI({
                req(dataset())
                colors <- cell_type_colors()
                cell_types_hex <- col2hex(get_cell_type_colors(dataset(), colors))
                cell_types <- names(get_cell_type_colors(dataset(), colors))

                shinyWidgets::pickerInput(
                    ns("vein_gene_foc_type"),
                    "Focus on type:",
                    choices = c("All", cell_types),
                    selected = "All",
                    choicesOpt = list(
                        style = paste0("color: ", cell_types_hex, ";")
                    )
                )
            })

            observe({
                req(input$color_gene_vein)
                shinyjs::toggle(id = "gene", condition = input$color_gene_vein == "Gene")
                shinyjs::toggle(id = "vein_gene_foc_type_select", condition = input$color_gene_vein == "Cell type")
            })

            output$gene_vein_plot <- renderPlot(
                {
                    req(input$color_gene_vein)
                    req(input$vein_gene_foc_type)

                    if (input$color_gene_vein == "Gene") {
                        req(input$gene)
                        gene <- input$gene
                    } else {
                        gene <- NULL
                    }

                    cell_type_colors <- get_cell_type_data(dataset())

                    color_order <- cell_type_colors$color
                    if (input$vein_gene_foc_type == "All") {
                        foc_type <- NULL
                        vein_rm_cell_types <- get_mc_config(dataset(), "vein_rm_cell_types")
                        if (!is.null(vein_rm_cell_types)) {
                            color_order <- cell_type_colors %>%
                                filter(!(cell_type %in% vein_rm_cell_types)) %>%
                                pull(color)
                        }
                    } else {
                        foc_type <- input$vein_gene_foc_type
                    }

                    plot_vein(dataset(), gene = gene, foc_type = foc_type, color_order = color_order, metacell_types = metacell_types())
                },
                res = 96
            ) %>% bindCache(dataset(), input$color_gene_vein, input$gene, input$vein_gene_foc_type, metacell_types(), cell_type_colors())

            # Metacell flow
            output$plot_metacell_flow <- renderPlot(
                {
                    req(has_network(dataset()))
                    req(input$selected_metacell)

                    mc_data <- metacell_types() %>% filter(metacell == input$selected_metacell)
                    plot_propagation_net_metacell(dataset(), input$selected_metacell, metacell_types = metacell_types()) +
                        ggtitle(glue("metacell #{input$selected_metacell} ({mc_data$cell_type[1]})")) +
                        theme(plot.title = element_text(colour = "darkblue"))
                } # ,
                # res = 96
            ) %>% bindCache(dataset(), input$selected_metacell, metacell_types())


            top_var_genes <- reactive({
                req(input$selected_metacell)
                req(has_network(dataset()))

                # Initially select the genes with highest variance along the trajectory
                get_top_var_genes(dataset(), input$selected_metacell)
            })

            output$traj_gene_selector <- renderUI({
                selectizeInput(ns("traj_genes"), "Select genes:", choices = gene_names(dataset()), selected = top_var_genes(), multiple = TRUE)
            })

            observe({
                req(input$selected_metacell)
                req(has_network(dataset()))
                updateSelectizeInput(session, "traj_genes", choices = gene_names(dataset()), server = TRUE, selected = top_var_genes(), options = list(maxItems = 5))
            })

            output$plot_mc_traj <- plotly::renderPlotly({
                req(input$traj_genes)
                req(has_network(dataset()))
                req(input$selected_metacell_traj)

                plotly::ggplotly(plot_gene_trajectory(dataset(), input$traj_genes, input$selected_metacell_traj, anchor_gene = NULL) + theme(axis.title.y = element_text(color = "darkblue")), source = "traj_plot", tooltip = "tooltip_text") %>% sanitize_plotly_buttons()
            }) %>% bindCache(dataset(), input$selected_metacell_traj, input$traj_genes)
        }
    )
}
