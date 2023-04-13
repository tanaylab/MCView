#' manifold UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_manifold_ui <- function(id) {
    ns <- NS(id)
    tagList(
        generic_column(
            width = 12,
            generic_box(
                id = ns("gene_projection"),
                title = "2D Projection",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                closable = FALSE,
                width = 12,
                height = "80vh",
                sidebar = shinydashboardPlus::boxSidebar(
                    startOpen = FALSE,
                    width = 25,
                    id = ns("gene_projection_sidebar"),
                    uiOutput(ns("graph_select_ui")),
                    selectInput(ns("proj_stat"), label = "Statistic", choices = c("Expression" = "expression", "Enrichment" = "enrichment"), selected = "Expression", multiple = FALSE, selectize = FALSE),
                    uiOutput(ns("set_range_ui")),
                    uiOutput(ns("expr_range_ui")),
                    uiOutput(ns("enrich_range_ui")),
                    uiOutput(ns("point_size_ui")),
                    uiOutput(ns("stroke_ui")),
                    uiOutput(ns("edge_distance_ui"))
                ),
                shinycssloaders::withSpinner(
                    plotly::plotlyOutput(ns("plot_manifold_proj_2d"), height = "80vh")
                )
            )
        )
    )
}


#' manifold sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_manifold_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            uiOutput(ns("color_by_selector")),
            uiOutput(ns("metadata_selector")),
            uiOutput(ns("gene_selectors")),
            uiOutput(ns("proj_gene_module_selector")),
            uiOutput(ns("top_correlated_select_gene1")),
            uiOutput(ns("top_correlated_select_gene2")),
            uiOutput(ns("genecards_buttons"))
        )
    )
}

#' manifold Server Function
#'
#' @noRd
mod_manifold_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            # gene selectors
            manifold_tab_gene_selectors(input, output, session, dataset, ns)

            picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE)

            output$color_by_selector <- renderUI({
                shinyWidgets::prettyRadioButtons(
                    ns("color_proj"),
                    label = "Color by:",
                    choices = c("Cell type", "Gene A", "Gene B", "Gene module", "Metadata"),
                    inline = FALSE,
                    status = "danger",
                    fill = TRUE
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

            output$proj_gene_module_selector <- renderUI({
                req(gene_modules())
                req(levels(gene_modules()$module))
                shinyWidgets::pickerInput(
                    ns("color_proj_gene_module"),
                    "Gene module:",
                    choices = levels(gene_modules()$module),
                    selected = NULL,
                    multiple = FALSE,
                    options = picker_options
                )
            })

            observe({
                req(input$color_proj)
                shinyjs::toggle(id = "metadata_selector", condition = input$color_proj == "Metadata")
                shinyjs::toggle(id = "proj_gene_module_selector", condition = input$color_proj == "Gene module")
            })

            projection_selectors(ns, dataset, output, input, gene_modules, globals, session, weight = 1)

            clipboard_changed <- clipboard_changed_2d_reactive(input, globals)

            # Projection plots
            output$plot_manifold_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot") %>%
                bindCache(dataset(), input$gene1, input$gene2, input$color_proj, metacell_types(), cell_type_colors(), input$point_size, input$stroke, input$min_edge_size, input$metacell1, input$metacell2, input$proj_stat, input$expr_range, input$lfp, input$set_range, input$color_proj_metadata, input$color_proj_gene_module, clipboard_changed(), input$graph_name)
        }
    )
}

manifold_tab_gene_selectors <- function(input, output, session, dataset, ns) {
    output$gene_selectors <- renderUI({
        div(
            id = ns("sidebar_select"),
            shinyWidgets::virtualSelectInput(ns("gene1"), "Gene A",
                choices = gene_names(dataset()),
                selected = default_gene1,
                multiple = FALSE,
                search = TRUE
            ),
            shinyWidgets::virtualSelectInput(ns("gene2"), "Gene B",
                choices = gene_names(dataset()),
                selected = default_gene2,
                multiple = FALSE,
                search = TRUE
            ),
            shinyWidgets::actionGroupButtons(ns("switch_genes"), labels = c("Switch"), size = "sm"),
            tags$hr()
        )
    })


    output$top_correlated_select_gene1 <- renderUI({
        req(input$gene1)
        req(has_gg_mc_top_cor(project, dataset()))
        tagList(
            selectInput(
                ns("selected_top_gene1"),
                glue("Top correlated to {input$gene1}:"),
                choices = c(get_top_cor_gene(dataset(), input$gene1, type = "pos"), get_top_cor_gene(dataset(), input$gene1, type = "neg")),
                selected = NULL,
                size = 10,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(c(ns("select_top_cor1_gene1"), ns("select_top_cor1_gene2")), labels = c("Select as Gene A", "Select as Gene B"), size = "sm")
        )
    })

    output$top_correlated_select_gene2 <- renderUI({
        req(input$gene2)
        req(has_gg_mc_top_cor(project, dataset()))
        tagList(
            selectInput(
                ns("selected_top_gene2"),
                glue("Top correlated to {input$gene2}:"),
                choices = c(get_top_cor_gene(dataset(), input$gene2, type = "pos"), get_top_cor_gene(dataset(), input$gene2, type = "neg")),
                selected = NULL,
                size = 10,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(c(ns("select_top_cor2_gene1"), ns("select_top_cor2_gene2")), labels = c("Select as Gene A", "Select as Gene B"), size = "sm"),
            tags$hr()
        )
    })

    output$genecards_buttons <- renderUI({
        req(input$gene1)
        req(input$gene2)
        tagList(
            shiny::actionButton(inputId = ns("genecards1"), label = glue("GeneCards: {input$gene1}"), onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={input$gene1}')")),
            shiny::actionButton(inputId = ns("genecards2"), label = glue("GeneCards: {input$gene2}"), onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={input$gene2}')"))
        )
    })

    observeEvent(input$select_top_cor1_gene1, {
        updateSelectInput(session, "gene1", selected = input$selected_top_gene1)
    })

    observeEvent(input$select_top_cor1_gene2, {
        updateSelectInput(session, "gene2", selected = input$selected_top_gene1)
    })

    observeEvent(input$select_top_cor2_gene1, {
        updateSelectInput(session, "gene1", selected = input$selected_top_gene2)
    })

    observeEvent(input$select_top_cor2_gene2, {
        updateSelectInput(session, "gene2", selected = input$selected_top_gene2)
    })

    observeEvent(input$switch_genes, {
        gene1 <- input$gene1
        gene2 <- input$gene2
        updateSelectInput(session, "gene1", selected = gene2)
        updateSelectInput(session, "gene2", selected = gene1)
    })
}
