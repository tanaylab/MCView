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
            qc_stat_box(ns, id, "# of cells per metacell", "plot_qc_cell_num"),
            qc_stat_box(ns, id, "Max zero-fold per metacell", "plot_mc_zero_fold"),
        ),
        resizable_column(
            width = 6,
            zero_fold_stat_box(ns, id, "# of cells with zero UMIs per gene", "plot_zero_fold")
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

            output$plot_qc_umis <- qc_stat_plot("umis", "Number of UMIs per metacell", dataset, input, "plot_qc_umis_type", log_scale = TRUE)
            output$plot_qc_cell_num <- qc_stat_plot("cells", "Number of cells per metacell", dataset, input, "plot_qc_cell_num_type")
            output$plot_qc_inner_fold <- qc_stat_plot("max_inner_fold", "Max inner-fold per metacell", dataset, input, "plot_qc_inner_fold_type")
            output$plot_mc_zero_fold <- qc_stat_plot("zero_fold", "Max zero-fold per metacell", dataset, input, "plot_mc_zero_fold_type")
            output$plot_zero_fold <- zero_fold_gene_plot(dataset, input)
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

qc_stat_box <- function(ns, id, title, output_id, width = 12, height = "27vh") {
    shinydashboardPlus::box(
        id = ns(id),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = width,
        sidebar = shinydashboardPlus::boxSidebar(
            startOpen = FALSE,
            id = ns(glue("stat_selector_{output_id}")),
            shinyWidgets::prettyRadioButtons(
                ns(glue("{output_id}_type")),
                label = "Plot type",
                choices = c("ECDF", "Density"),
                inline = TRUE,
                status = "danger",
                fill = TRUE,
                selected = "Density"
            )
        ),
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns(output_id), height = height)
        )
    )
}

zero_fold_stat_box <- function(ns, id, title, output_id, width = 12, height = "35vh") {
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

zero_fold_gene_plot <- function(dataset, input) {
    plotly::renderPlotly({
        zero_fold_df <- get_mc_data(dataset(), "gene_zero_fold")
        req(zero_fold_df)

        limits <- c(0, max(zero_fold_df$obs, zero_fold_df$exp) + 1)

        p <- zero_fold_df %>%
            rename(Observed = obs, Expected = exp) %>%
            ggplot(aes(x = Observed, y = Expected, label = gene, metacell = metacell, FC = zero_fold)) +
            geom_point(size = 0.5) +
            geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
            xlab("# of cells with 0 UMIs (observed)") +
            ylab("# of cells with 0 UMIs (expected)") +
            ylim(limits) +
            xlim(limits)

        plotly::ggplotly(p) %>%
            sanitize_plotly_buttons()
    }) %>% bindCache(dataset())
}

qc_stat_plot <- function(field, xlab, dataset, input, plot_type_id, ylab = NULL, log_scale = FALSE) {
    plotly::renderPlotly({
        qc_df <- as_tibble(get_mc_data(dataset(), "mc_qc_metadata"))
        req(qc_df[[field]])

        req(input[[plot_type_id]])

        if (input[[plot_type_id]] == "ECDF") {
            p <- qc_ecdf(qc_df, field, xlab, ylab, log_scale = log_scale)
        } else {
            p <- qc_density(qc_df, field, xlab, ylab, log_scale = log_scale)
        }

        return(p)
    }) %>% bindCache(dataset(), input[[plot_type_id]])
}

qc_density <- function(qc_df, field, xlab, ylab, log_scale = FALSE) {
    if (is.null(ylab)) {
        ylab <- "Density"
    }
    quants <- quantile(qc_df[[field]], probs = c(0.1, 0.5, 0.9))

    p <- qc_df %>%
        ggplot(aes(x = !!sym(field))) +
        geom_density(fill = "black", alpha = 0.2, color = "black") +
        xlab(xlab) +
        ylab(ylab) +
        geom_vline(xintercept = quants[2], color = "red", linetype = "dashed") +
        geom_vline(xintercept = quants[c(1, 3)], color = "black", linetype = "dashed")

    quants_text <- c(
        glue("bottom 10%: {round(quants[1], digits = 2)}"),
        glue("median: {round(quants[2], digits = 2)}"),
        glue("top 10%: {round(quants[3], digits = 2)}")
    )

    if (log_scale) {
        p <- p + scale_x_log10(labels = scales::scientific)
        quants <- log10(quants)
    }

    max_y <- layer_scales(p)$y$get_limits()[2]

    p <- plotly::ggplotly(p) %>%
        plotly::add_annotations(
            x = quants,
            y = c(0.75, 0.75, 0.75) * max_y,
            text = quants_text,
            showarrow = FALSE,
            font = list(color = "black", size = 10),
            yshift = 0,
            xshift = -10,
            textangle = 90,
            ax = 0,
            ay = 0
        ) %>%
        sanitize_plotly_buttons()

    return(p)
}

qc_ecdf <- function(qc_df, field, xlab, ylab, log_scale = FALSE) {
    if (is.null(ylab)) {
        ylab <- "% of metacells <= x"
    }

    quants <- quantile(qc_df[[field]], probs = c(0.1, 0.5, 0.9))

    p <- tibble(
        x = qc_df[[field]],
        y = 1 - ecdf(x)(x)
    ) %>%
        ggplot(aes(x = x, y = y)) +
        geom_line() +
        xlab(xlab) +
        ylab(ylab) +
        scale_y_continuous(labels = scales::percent) +
        geom_vline(xintercept = quants[2], color = "red", linetype = "dashed") +
        geom_vline(xintercept = quants[c(1, 3)], color = "gray", linetype = "dashed") +
        theme_bw()

    quants_text <- c(
        glue("bottom 10%: {round(quants[1], digits = 2)}"),
        glue("median: {round(quants[2], digits = 2)}"),
        glue("top 10%: {round(quants[3], digits = 2)}")
    )

    if (log_scale) {
        p <- p + scale_x_log10(labels = scales::scientific)
        quants <- log10(quants)
    }

    p <- plotly::ggplotly(p) %>%
        plotly::add_annotations(
            x = quants,
            y = c(0.3, 0.3, 0.75),
            text = quants_text,
            showarrow = FALSE,
            font = list(color = "darkgray", size = 10),
            yshift = 0,
            xshift = -10,
            textangle = 90,
            ax = 0,
            ay = 0
        ) %>%
        sanitize_plotly_buttons()


    return(p)
}
