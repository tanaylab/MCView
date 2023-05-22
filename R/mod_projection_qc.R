#' projection QC UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_projection_qc_ui <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            width = 12,
            fluidRow(
                shinydashboard::valueBoxOutput(ns("num_metacells_atlas"), width = 3),
                shinydashboard::valueBoxOutput(ns("num_metacells_query"), width = 3),
                shinydashboard::valueBoxOutput(ns("num_metacells_similar"), width = 3),
                shinydashboard::valueBoxOutput(ns("num_disjoined_genes"), width = 3),
            )
        ),
        generic_column(
            width = 6,
            qc_stat_box(ns, id, "Projected correlation per metacell", "plot_projected_correlation"),
        ),
        generic_column(
            width = 6,
            gene_correction_factor_stat_box(ns, id, "Correction factor per gene", "plot_correction_factor_scatter")
        )
    )
}


#' projection QC sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_projection_qc_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            div()
        )
    )
}

#' QC Server Function
#'
#' @noRd
mod_projection_qc_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            output$num_metacells_atlas <- shinydashboard::renderValueBox({
                mc_mat <- get_mc_data(dataset(), "mc_mat", atlas = TRUE)
                req(mc_mat)
                num_metacells <- ncol(mc_mat)
                shinydashboard::valueBox(
                    scales::comma(num_metacells),
                    "Number of atlas metacells",
                    color = "black"
                )
            })
            output$num_metacells_query <- qc_value_box("n_metacells", "Number of query metacells", dataset, color = "purple")
            output$num_metacells_similar <- shinydashboard::renderValueBox({
                mc_mat <- get_mc_data(dataset(), "mc_mat")
                req(mc_mat)

                num_metacells <- ncol(mc_mat)
                req(num_metacells)

                md <- get_mc_data(dataset(), "metadata")
                req(md)
                req(md$similar)
                num_similar <- sum(md$similar == "similar", na.rm = TRUE)

                p_similar <- num_similar / num_metacells
                if (p_similar <= 0.8) {
                    color <- "red"
                } else {
                    color <- "green"
                }
                p_similar <- scales::percent(p_similar)

                shinydashboard::valueBox(
                    glue("{scales::comma(num_similar)} ({p_similar})"),
                    "Number of similar metacells",
                    color = color
                )
            })

            output$num_disjoined_genes <- shinydashboard::renderValueBox({
                disjoined_genes_no_atlas <- get_mc_data(dataset(), "disjoined_genes_no_atlas")
                disjoined_genes_no_query <- get_mc_data(dataset(), "disjoined_genes_no_query")

                disjoined <- intersect(disjoined_genes_no_atlas, disjoined_genes_no_query)

                req(!is.null(disjoined))

                num_disjoined <- length(disjoined)

                shinydashboard::valueBox(
                    scales::comma(num_disjoined),
                    "Number of disjoined genes",
                    color = "maroon"
                )
            })

            output$plot_projected_correlation <- qc_stat_plot("projected_correlation", "Projected correlation per metacell", dataset, input, "plot_projected_correlation_type")

            output$plot_correction_factor_scatter <- gene_correction_factor_scatter_plot(dataset, input)
            output$gene_correction_factor_table <- gene_correction_factor_table(dataset, input)
        }
    )
}



gene_correction_factor_stat_box <- function(ns, id, title, output_id, width = 12, height = "35vh") {
    generic_box(
        id = ns(id),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = width,
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns(output_id), height = height)
        ),
        shinyWidgets::prettySwitch(inputId = ns("show_correction_factor_table"), value = FALSE, label = "Show table"),
        DT::DTOutput(ns("gene_correction_factor_table"))
    )
}

gene_correction_factor_scatter_plot <- function(dataset, input) {
    plotly::renderPlotly({
        gene_qc <- get_mc_data(dataset(), "gene_inner_fold")
        if (is.null(gene_qc)) {
            return(plotly_text_plot("Please recompute the metacells\nusing the latest version\nin order to see this plot."))
        }
        req(gene_qc)

        req(gene_qc$correction_factor)


        p <- gene_qc %>%
            mutate(max_expr = log2(max_expr + 1e-5)) %>%
            rename(Gene = gene, `Correction factor` = correction_factor, `Max expression` = max_expr, Type = type) %>%
            ggplot(aes(x = `Max expression`, y = `Correction factor`, label = Gene, color = Type)) +
            scale_color_manual(values = c("other" = "gray", "lateral" = "red", "noisy" = "purple")) +
            geom_point(size = 0.5) +
            xlab("log2(gene expression)") +
            ylab("Correction factor")

        plotly::ggplotly(p) %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons()
    }) %>% bindCache(dataset())
}

gene_correction_factor_table <- function(dataset, input) {
    DT::renderDT(
        if (input$show_correction_factor_table) {
            gene_qc <- get_mc_data(dataset(), "gene_inner_fold")
            req(gene_qc)
            req(gene_qc$correction_factor)
            gene_qc %>%
                arrange(desc(correction_factor)) %>%
                mutate(max_expr = log2(max_expr + 1e-5)) %>%
                select(Gene = gene, `Correction factor` = correction_factor, `Max expression` = max_expr, Type = type) %>%
                mutate(`Max expression` = round(`Max expression`, digits = 2)) %>%
                DT::datatable(
                    rownames = FALSE,
                    options = list(
                        pageLength = 20,
                        scrollX = TRUE,
                        scrollY = "300px",
                        scrollCollapse = TRUE,
                        dom = "ftp",
                        columnDefs = list(
                            list(
                                targets = 0,
                                width = "100px"
                            )
                        )
                    )
                )
        }
    )
}
