#' QC UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_qc_ui <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            width = 12,
            fluidRow(
                shinydashboard::valueBoxOutput(ns("num_umis"), width = 3),
                shinydashboard::valueBoxOutput(ns("num_cells"), width = 2),
                shinydashboard::valueBoxOutput(ns("num_outliers"), width = 3),
                shinydashboard::valueBoxOutput(ns("median_umis_per_metacell"), width = 2),
                shinydashboard::valueBoxOutput(ns("median_cells_per_metacell"), width = 2)
            )
        ),
        resizable_column(
            width = 6,
            qc_stat_box(ns, id, "# of UMIs per metacell", "plot_qc_umis"),
            qc_stat_box(ns, id, "Max inner-fold per metacell", "plot_qc_inner_fold")
        ),
        resizable_column(
            width = 6,
            qc_stat_box(ns, id, "# of cells per metacell", "plot_qc_cell_num")
            # qc_stat_box(ns, id, "Max zero-fold per metacell", "plot_qc_zero_fold")
        )
    )
}


#' QC sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_qc_sidebar_ui <- function(id) {
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
mod_qc_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            # Value boxes
            output$num_umis <- qc_value_box("n_umis", "Total number of UMIs", dataset, color = "black")
            output$num_cells <- qc_value_box("n_cells", "Number of cells", dataset, color = "purple")
            output$num_outliers <- shinydashboard::renderValueBox({
                num_cells <- get_mc_data(dataset(), "qc_stats")$n_cells
                num_outliers <- get_mc_data(dataset(), "qc_stats")$n_outliers
                req(num_cells)
                req(num_outliers)
                p_outliers <- num_outliers / num_cells
                if (p_outliers >= 0.2) {
                    color <- "red"
                } else {
                    color <- "green"
                }
                shinydashboard::valueBox(
                    glue("{scales::comma(num_outliers)} ({scales::percent(p_outliers)})"),
                    "Number of outlier cells",
                    color = color
                )
            })
            output$median_umis_per_metacell <- qc_value_box("median_umis_per_metacell", "Median UMIs / MC", dataset, color = "blue")
            output$median_cells_per_metacell <- qc_value_box("median_cells_per_metacell", "Median cells / MC", dataset, color = "maroon")

            output$plot_qc_umis <- qc_stat_plot("umis", "Number of UMIs per metacell", dataset)
            output$plot_qc_cell_num <- qc_stat_plot("cells", "Number of cells per metacell", dataset)
            output$plot_qc_inner_fold <- qc_stat_plot("max_inner_fold", "Max inner-fold per metacell", dataset)
        }
    )
}

qc_value_box <- function(field, title, dataset, color = "black") {
    shinydashboard::renderValueBox({
        qc_stats <- get_mc_data(dataset(), "qc_stats")
        req(qc_stats)
        stat <- qc_stats[[field]]
        req(stat)
        shinydashboard::valueBox(
            scales::comma(stat),
            title,
            color = color
        )
    })
}

qc_stat_box <- function(ns, id, title, output_id, width = 12, height = "30vh") {
    shinydashboardPlus::box(
        id = ns(id),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = width,
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns(output_id), height = height)
        )
    )
}

# qc_stat_plot <- function(field, xlab, dataset, ylab = "Density"){
#     plotly::renderPlotly({
#         qc_df <- as_tibble(get_mc_data(dataset(), "mc_qc_metadata"))
#         req(qc_df[[field]])

#         med <- median(qc_df[[field]])

#         p <- qc_df %>%
#             ggplot(aes(x = !!sym(field))) +
#             geom_density(fill = "black", alpha = 0.5) +
#             xlab(xlab) +
#             ylab(ylab) +
#             geom_vline(xintercept = med, color = "red", linetype = "dashed")
#         plotly::ggplotly(p) %>%
#             sanitize_plotly_buttons()
#     })
# }

qc_stat_plot <- function(field, xlab, dataset, ylab = "% of metacells <= x") {
    plotly::renderPlotly({
        qc_df <- as_tibble(get_mc_data(dataset(), "mc_qc_metadata"))
        req(qc_df[[field]])

        p <- tibble(
            x = qc_df[[field]],
            y = 1 - ecdf(x)(x)
        ) %>%
            ggplot(aes(x = x, y = y)) +
            geom_line() +
            xlab(xlab) +
            ylab(ylab) +
            scale_y_continuous(labels = scales::percent) +
            theme_bw()

        plotly::ggplotly(p) %>%
            sanitize_plotly_buttons()
    })
}
