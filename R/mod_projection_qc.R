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
                shinydashboard::valueBoxOutput(ns("num_metacells_atlas"), width = 2),
                shinydashboard::valueBoxOutput(ns("num_metacells_query"), width = 2),
                shinydashboard::valueBoxOutput(ns("num_metacells_similar"), width = 2),
                shinydashboard::valueBoxOutput(ns("avg_projection_cor"), width = 2),
                shinydashboard::valueBoxOutput(ns("common_genes"), width = 2),
                shinydashboard::valueBoxOutput(ns("fitted_genes"), width = 2)
            )
        ),
        generic_column(
            width = 7,
            generic_box(
                id = ns("metacell_projection"),
                title = "Type Projections",
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                closable = FALSE,
                width = 12,
                shinycssloaders::withSpinner(
                    plotOutput(ns("plot_mc_stacked_type"))
                )
            ),
            qc_stat_box(ns, id, "Projected correlation per metacell", "plot_projected_correlation")
        ),
        generic_column(
            width = 5,
            uiOutput(ns("correction_factor_box")),
            fitted_genes_per_cell_type_stat_box(ns, id, "Fitted genes per cell type", "plot_fitted_genes_per_cell_type"),
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
                    p_similar,
                    "Percentage of 'similar' metacells",
                    color = color
                )
            })

            output$common_genes <- shinydashboard::renderValueBox({
                query_mat <- get_mc_data(dataset(), "mc_mat")
                atlas_mat <- get_mc_data(dataset(), "mc_mat", atlas = TRUE)
                req(query_mat)
                req(atlas_mat)

                common_genes <- intersect(rownames(query_mat), rownames(atlas_mat))

                p_atlas <- scales::percent(length(common_genes) / nrow(atlas_mat))
                p_query <- scales::percent(length(common_genes) / nrow(query_mat))

                shinydashboard::valueBox(
                    scales::comma(length(common_genes)),
                    glue("Common genes ({p_atlas} of atlas, {p_query} of query)"),
                    color = "blue"
                )
            })

            output$fitted_genes <- shinydashboard::renderValueBox({
                gene_md <- get_mc_data(dataset(), "gene_metadata")
                req(gene_md)
                req(has_name(gene_md, "fitted_gene_any"))

                gene_md <- gene_md %>%
                    filter(marker_gene) %>%
                    distinct(gene, marker_gene, fitted_gene_any)

                num_fitted_genes <- sum(gene_md$fitted_gene_any, na.rm = TRUE)
                num_markers <- sum(gene_md$marker_gene, na.rm = TRUE)
                p_fitted <- num_fitted_genes / num_markers
            
                if (p_fitted <= 0.33) {
                    color <- "red"
                } else if (p_fitted <= 0.66) {
                    color <- "yellow"
                } else {
                    color <- "green"
                }

                shinydashboard::valueBox(
                    scales::percent(p_fitted, accuracy = 0.1),
                    "% of 'fitted' genes out of markers",
                    color = color
                )
            })

            output$avg_projection_cor <- shinydashboard::renderValueBox({
                qc_df <- as_tibble(get_mc_data(dataset(), "mc_qc_metadata"))
                req(qc_df)
                req(qc_df$projected_correlation)

                avg_projection_cor <- mean(qc_df$projected_correlation, na.rm = TRUE)

                shinydashboard::valueBox(
                    round(avg_projection_cor, 2),
                    "Average projection correlation",
                    color = "maroon"
                )
            })

            output$plot_projected_correlation <- qc_stat_plot("projected_correlation", "Projected correlation per metacell", dataset, input, "plot_projected_correlation_type")

            output$correction_factor_box <- gene_correction_factor_stat_box(ns, id, dataset, "Correction factor per gene", "plot_correction_factor_scatter")

            output$plot_correction_factor_scatter <- gene_correction_factor_scatter_plot(dataset, input)
            output$gene_correction_factor_table <- gene_correction_factor_table(dataset, input)

            output$plot_mc_stacked_type <- plot_type_predictions_bar(dataset, metacell_types, cell_type_colors)

            output$plot_fitted_genes_per_cell_type <- fitted_genes_per_cell_type_plot(dataset, input)
            output$fitted_genes_per_cell_type_table <- fitted_genes_per_cell_type_table(dataset, input)

            output$fitted_gene_per_cell_type_selector <- fitted_gene_per_cell_type_selector(ns, dataset, input)
        }
    )
}



gene_correction_factor_stat_box <- function(ns, id, dataset, title, output_id, width = 12, height = "35vh") {
    renderUI({
        gene_qc <- get_gene_qc(dataset())
        req(gene_qc)
        req(has_name(gene_qc, "correction_factor"))
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
    })
}

gene_correction_factor_scatter_plot <- function(dataset, input) {
    plotly::renderPlotly({
        gene_qc <- get_gene_qc(dataset())
        if (is.null(gene_qc) || is.null(gene_qc$correction_factor)) {
            return(plotly_text_plot("Please recompute the metacells\nusing the latest version\nin order to see this plot."))
        }
        req(gene_qc)

        req(gene_qc$correction_factor)

        p <- gene_qc %>%
            filter(correction_factor != 0, correction_factor != 1) %>%
            mutate(max_expr = log2(max_expr + 1e-5)) %>%
            mutate(correction_factor = log2(correction_factor)) %>%
            rename(Gene = gene, `Correction factor` = correction_factor, `Max expression` = max_expr, Type = type) %>%
            ggplot(aes(x = `Max expression`, y = `Correction factor`, label = Gene, color = Type)) +
            scale_color_manual(values = c("other" = "gray", "lateral" = "blue", "noisy" = "red", "lateral, noisy" = "purple")) +
            geom_point(size = 0.5) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            xlab("log2(gene expression)") +
            ylab("log2(correction factor)")

        plotly::ggplotly(p) %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons()
    }) %>% bindCache(dataset())
}

gene_correction_factor_table <- function(dataset, input) {
    DT::renderDT(
        if (input$show_correction_factor_table) {
            gene_qc <- get_gene_qc(dataset())
            req(gene_qc)
            req(gene_qc$correction_factor)
            gene_qc %>%
                filter(correction_factor != 0, correction_factor != 1) %>%
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

fitted_genes_per_cell_type_stat_box <- function(ns, id, title, output_id, width = 12, height = "35vh") {
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
        shinyWidgets::prettySwitch(inputId = ns("show_genes_per_cell_type_table"), value = FALSE, label = "Show table"),
        uiOutput(ns("fitted_gene_per_cell_type_selector")),
        DT::DTOutput(ns("fitted_genes_per_cell_type_table"))
    )
}

fitted_gene_per_cell_type_selector <- function(ns, dataset, input) {
    renderUI({
        req(has_atlas(dataset()))
        cell_types <- get_mc_data(dataset(), "cell_type_colors", atlas = TRUE)$cell_type
        req(cell_types)
        req(input$show_genes_per_cell_type_table)
        shinyWidgets::pickerInput(
            inputId = ns("fitted_gene_per_cell_type_selector"),
            label = "Select cell type",
            choices = c("all", cell_types),
            selected = "all",
            multiple = FALSE
        )
    })
}


fitted_genes_per_cell_type_table <- function(dataset, input) {
    DT::renderDT(
        if (input$show_genes_per_cell_type_table) {
            gene_qc <- get_gene_qc(dataset())
            req(gene_qc)
            req(any(grepl("fitted_gene_of", colnames(gene_qc))))

            m <- gene_qc %>%
                select(starts_with("fitted_gene_of")) %>%
                as.matrix()

            sums <- rowSums(m)
            common_set <- gene_qc$gene[sums == ncol(m)]
            fitted_set <- sums > 0

            df <- gene_qc[fitted_set, ] %>%
                select(gene, starts_with("fitted_gene_of")) %>%
                gather(key = "type", value = "fitted_gene_of", -gene) %>%
                filter(fitted_gene_of) %>%
                mutate(type = gsub("fitted_gene_of_", "", type)) %>%
                mutate(common = gene %in% common_set)

            if (!is.null(input$fitted_gene_per_cell_type_selector) && input$fitted_gene_per_cell_type_selector != "all") {
                df <- df %>%
                    filter(type == input$fitted_gene_per_cell_type_selector)
            }
            df %>%
                select(gene, type, common) %>%
                arrange(common, type, gene) %>%
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

fitted_genes_per_cell_type_plot <- function(dataset, input) {
    plotly::renderPlotly({
        gene_qc <- get_gene_qc(dataset())
        if (is.null(gene_qc)) {
            return(plotly_text_plot("Please recompute the metacells\nusing the latest version\nin order to see this plot."))
        }

        req(gene_qc)

        req(any(grepl("fitted_gene_of", colnames(gene_qc))))

        m <- gene_qc %>%
            select(starts_with("fitted_gene_of")) %>%
            as.matrix()

        common_set <- rowSums(m) == ncol(m)

        m_f <- m[!common_set, , drop = FALSE]

        fitted_per_type <- tibble::enframe(colSums(m_f), "type", "n") %>%
            mutate(type = gsub("fitted_gene_of_", "", type))


        req(has_atlas(dataset()))
        atlas_colors <- get_mc_data(dataset(), "cell_type_colors", atlas = TRUE) %>%
            select(cell_type, color) %>%
            tibble::deframe()

        fitted_per_type <- fitted_per_type %>%
            mutate(type = factor(type, levels = names(atlas_colors)))

        p <- fitted_per_type %>%
            ggplot(aes(x = type, y = n, fill = type)) +
            geom_col() +
            xlab("Cell type") +
            ylab("Number of non-common fitted genes") +
            scale_fill_manual(values = atlas_colors) +
            guides(fill = "none") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            ggtitle(glue("Common set: {sum(common_set)} genes"))

        plotly::ggplotly(p) %>% sanitize_plotly_buttons()
    })
}
