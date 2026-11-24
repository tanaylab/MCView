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
        value_box_row(
            value_box_output(ns("num_metacells")),
            value_box_output(ns("median_umis_per_metacell")),
            value_box_output(ns("num_cells")),
            value_box_output(ns("median_cells_per_metacell")),
            value_box_output(ns("num_outliers"))
        ),
        fluidRow(
            generic_column(
                width = 6,
                qc_stat_box(ns, id, "# of UMIs per metacell", "plot_qc_umis"),
                qc_stat_box(ns, id, "Max inner-fold per metacell", "plot_qc_inner_fold",
                    width = 12, height = "27vh",
                    checkboxInput(ns("include_lateral"), label = "Include lateral", value = FALSE)
                ),
                gene_inner_fold_stat_box(ns, id, "# of metacells with significant inner-fold", "plot_gene_inner_fold_scatter"),
                qc_stat_box(ns, id, "Max inner-stdev per metacell", "plot_qc_std")
            ),
            generic_column(
                width = 6,
                qc_stat_box(ns, id, "# of cells per metacell", "plot_qc_cell_num"),
                qc_stat_box(ns, id, "Max zero-fold per metacell", "plot_mc_zero_fold"),
                zero_fold_stat_box(ns, id, "# of cells with zero UMIs per gene", "plot_zero_fold")
            )
        )
    )
}

# A row of value boxes (of value_box_output()s), aligned with the boxes below it: as many equal columns
# as fit (at least `min_width` wide each, see mcview.css), with numbers that shrink to fit their box.
# Bootstrap's columns (valueBoxOutput(width = 2)) clipped the numbers at the box's edge, or ran them into
# the next box, on narrower screens.
value_box_row <- function(..., min_width = NULL) {
    div(
        class = "mcv-value-boxes",
        style = if (!is.null(min_width)) glue("--mcv-value-box-min-width: {min_width};"),
        ...
    )
}

# shinydashboard::valueBoxOutput() without bootstrap's column width (see value_box_row())
value_box_output <- function(output_id) {
    shiny::uiOutput(output_id, class = "mcv-value-box")
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
            output$num_metacells <- qc_value_box("n_metacells", "Number of metacells", dataset, color = "black")
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
                if (p_outliers < 0.01) {
                    p_outliers <- scales::percent(p_outliers, accuracy = 0.01)
                } else {
                    p_outliers <- scales::percent(p_outliers)
                }
                shinydashboard::valueBox(
                    # (the percentage smaller, so that the value fits in the box)
                    tagList(scales::comma(num_outliers), " ", tags$small(glue("({p_outliers})"))),
                    "Number of outlier cells",
                    color = color
                )
            })
            output$median_umis_per_metacell <- qc_value_box("median_umis_per_metacell", "Median UMIs / MC", dataset, color = "blue")
            output$median_cells_per_metacell <- qc_value_box("median_cells_per_metacell", "Median cells / MC", dataset, color = "maroon")

            output$plot_qc_umis <- qc_stat_plot("umis", "Number of UMIs per metacell", dataset, input, "plot_qc_umis_type", globals, log_scale = TRUE)
            output$plot_qc_cell_num <- qc_stat_plot("cells", "Number of cells per metacell", dataset, input, "plot_qc_cell_num_type", globals)
            output$plot_qc_inner_fold <- qc_stat_plot(c("max_inner_fold", "max_inner_fold_no_lateral"), "Max inner-fold per metacell", dataset, input, "plot_qc_inner_fold_type", globals, field_input = "include_lateral")
            output$plot_qc_std <- qc_stat_plot("max_inner_stdev_log", "Max stdev(log(fractions)) per metacell", dataset, input, "plot_qc_std_type", globals)
            output$plot_mc_zero_fold <- qc_stat_plot("zero_fold", "Max log2(# of zero cells / expected) per metacell", dataset, input, "plot_mc_zero_fold_type", globals)

            output$plot_zero_fold <- zero_fold_gene_plot(dataset, input, globals)
            output$zero_fold_table <- zero_fold_table(dataset, input)

            output$plot_gene_inner_fold_scatter <- gene_inner_fold_scatter_plot(dataset, input, globals)
            output$gene_inner_fold_table <- gene_inner_fold_table(dataset, input)
        }
    )
}

qc_value_box <- function(field, title, dataset, color = "black") {
    shinydashboard::renderValueBox({
        if (field == "n_metacells") {
            # (mc_sum has a value per metacell: reading mc_mat here took seconds for large projects)
            stat <- length(get_mc_data(dataset(), "mc_sum"))
        } else {
            qc_stats <- get_mc_data(dataset(), "qc_stats")
            req(qc_stats)
            stat <- qc_stats[[field]]
        }
        req(stat)

        shinydashboard::valueBox(
            scales::comma(stat),
            title,
            color = color
        )
    })
}

qc_stat_box <- function(ns, id, title, output_id, width = 12, height = "27vh", ...) {
    generic_box(
        id = ns(paste0(output_id, "_box")),
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
            ),
            ...
        ),
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns(output_id), height = height)
        )
    )
}

zero_fold_stat_box <- function(ns, id, title, output_id, width = 12, height = "35vh") {
    generic_box(
        id = ns(paste0(output_id, "_box")),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = width,
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns(output_id), height = height)
        ),
        shinyWidgets::prettySwitch(inputId = ns("show_zero_fold_table"), value = FALSE, label = "Show table"),
        DT::DTOutput(ns("zero_fold_table"))
    )
}

zero_fold_table <- function(dataset, input) {
    DT::renderDT(
        if (input$show_zero_fold_table) {
            zero_fold_df <- get_mc_data(dataset(), "gene_zero_fold")
            req(zero_fold_df)
            zero_fold_df %>%
                filter(avg >= -10) %>%
                slice(1:100) %>%
                select(Gene = gene, Observed = obs, Expected = exp, FC = zero_fold, Type = type, Expression = avg, Metacell = metacell) %>%
                mutate(Expected = round(Expected, digits = 1), FC = round(FC, digits = 2), Expression = round(Expression, digits = 2)) %>%
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

zero_fold_gene_plot <- function(dataset, input, globals) {
    plotly::renderPlotly({
        zero_fold_df <- get_mc_data(dataset(), "gene_zero_fold")
        req(zero_fold_df)

        df <- zero_fold_df %>%
            rename(Expression = avg, FC = zero_fold, Gene = gene, Type = type, Metacell = metacell) %>%
            mutate(Observed = obs, Expected = formatC(exp, format = "f", digits = 1))

        # group = Type: by default ggplot makes a group of every Type and Gene pair (a group per point)
        p <- df %>%
            ggplot(aes(x = Expression, y = FC, color = Type, customdata = Gene, group = Type)) +
            scale_color_manual(values = c("other" = "gray", "lateral" = "blue", "noisy" = "red", "lateral, noisy" = "purple")) +
            geom_point(size = 0.5) +
            geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
            xlab("log2(gene expression)") +
            ylab("log2(# of zero cells / expected)")

        # set_hovertemplate() replaces the hover texts
        plotly::ggplotly(p, tooltip = "none") %>%
            set_hovertemplate(df, "Gene", c("Observed", "Expected", "Metacell"), "Expression: %{x:.3~r}<br />FC: %{y:.3~r}<br />Gene: %{customdata[0]}<br />Type: %{fullData.name}<br />Observed: %{customdata[1]}<br />Expected: %{customdata[2]}<br />Metacell: %{customdata[3]}<extra></extra>") %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons() %>%
            sanitize_plotly_download(globals)
    }) %>% bindCache(dataset(), globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)
}

gene_inner_fold_stat_box <- function(ns, id, title, output_id, width = 12, height = "35vh") {
    generic_box(
        id = ns(paste0(output_id, "_box")),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = width,
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns(output_id), height = height)
        ),
        shinyWidgets::prettySwitch(inputId = ns("show_gene_inner_fold_table"), value = FALSE, label = "Show table"),
        DT::DTOutput(ns("gene_inner_fold_table"))
    )
}

gene_inner_fold_table <- function(dataset, input) {
    DT::renderDT(
        if (input$show_gene_inner_fold_table) {
            gene_inner_fold_df <- get_gene_qc(dataset())
            req(gene_inner_fold_df)
            req(gene_inner_fold_df$significant_inner_folds_count)
            gene_inner_fold_df %>%
                filter(significant_inner_folds_count > 0) %>%
                mutate(max_expr = log2(max_expr + 1e-5)) %>%
                select(Gene = gene, `# of metacells` = significant_inner_folds_count, `Max expression` = max_expr, Type = type) %>%
                mutate(`# of metacells` = round(`# of metacells`, digits = 1), `Max expression` = round(`Max expression`, digits = 2)) %>%
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

gene_inner_fold_scatter_plot <- function(dataset, input, globals) {
    plotly::renderPlotly({
        gene_inner_fold_df <- get_gene_qc(dataset())
        if (is.null(gene_inner_fold_df) || is.null(gene_inner_fold_df$significant_inner_folds_count)) {
            return(plotly_text_plot("Please recompute the metacells\nusing the latest version\nin order to see this plot."))
        }
        req(gene_inner_fold_df)
        req(gene_inner_fold_df$significant_inner_folds_count)

        p <- gene_inner_fold_df %>%
            mutate(max_expr = log2(max_expr + 1e-5)) %>%
            rename(Gene = gene, `# of metacells` = significant_inner_folds_count, `Max expression` = max_expr, Type = type) %>%
            ggplot(aes(x = `Max expression`, y = `# of metacells`, color = Type, customdata = Gene, group = Type)) +
            scale_color_manual(values = c("other" = "gray", "lateral" = "blue", "noisy" = "red", "lateral, noisy" = "purple")) +
            geom_point(size = 0.5) +
            xlab("log2(gene expression)") +
            ylab("# of metacells with significant inner-fold")

        # set_hovertemplate() replaces the hover texts
        plotly::ggplotly(p, tooltip = "none") %>%
            set_hovertemplate(gene_inner_fold_df, "gene", c(), "Max expression: %{x:.3~r}<br /># of metacells: %{y:~r}<br />Gene: %{customdata[0]}<br />Type: %{fullData.name}<extra></extra>") %>%
            sanitize_for_WebGL() %>%
            plotly::toWebGL() %>%
            sanitize_plotly_buttons() %>%
            sanitize_plotly_download(globals)
    }) %>% bindCache(dataset(), globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)
}

qc_stat_plot <- function(field, xlab, dataset, input, plot_type_id, globals, ylab = NULL, log_scale = FALSE, field_input = plot_type_id) {
    plotly::renderPlotly({
        qc_df <- as_tibble(get_mc_data(dataset(), "mc_qc_metadata"))

        if (length(field) > 1) {
            req(length(field) == 2)
            req(!is.null(field_input))
            req(field_input != plot_type_id)
            req(length(field_input) == 1)
            req(!is.null(input[[field_input]]))
            if (input[[field_input]]) {
                field <- field[1]
            } else {
                field <- field[2]
            }
        }

        req(qc_df[[field]])

        req(input[[plot_type_id]])

        if (input[[plot_type_id]] == "ECDF") {
            p <- qc_ecdf(qc_df, field, xlab, ylab, globals, log_scale = log_scale)
        } else {
            p <- qc_density(qc_df, field, xlab, ylab, globals, log_scale = log_scale)
        }

        p <- sanitize_plotly_download(p, globals)

        return(p)
    }) %>% bindCache(dataset(), input[[plot_type_id]], field, plot_type_id, input[[field_input]], globals$plotly_format, globals$plotly_width, globals$plotly_height, globals$plotly_scale)
}

# The QC plots are drawn with plotly directly, as ggplotly() drew the ggplot of geom_density() with
# geom_vline()s and theme_classic(): converting the ggplot took ~0.3 s per plot, and the QC tab, often the
# first tab users see, has five of them.
qc_density <- function(qc_df, field, xlab, ylab, globals, log_scale = FALSE) {
    if (is.null(ylab)) {
        ylab <- "Density"
    }
    values_quants <- quantile(qc_df[[field]], probs = c(0.1, 0.5, 0.9))

    x <- qc_df[[field]]
    quants <- values_quants
    if (log_scale) {
        x <- log10(x)
        quants <- log10(quants)
    }
    x <- x[is.finite(x)]
    # as geom_density(): bw.nrd0(), 512 points over the range of the values
    d <- stats::density(x, n = 512, from = min(x), to = max(x))
    max_y <- max(d$y)
    # ggplot's axes: 5% wider than the data, its breaks and labels
    y_range <- c(-0.05, 1.05) * max_y
    y_breaks <- scales::breaks_extended()(c(0, max_y))
    # the hover texts in the units of the axis (e.g. UMIs, not log10(UMIs))
    d_values <- if (log_scale) 10^d$x else d$x

    p <- plotly::plot_ly(
        x = d$x, y = d$y, text = paste0(xlab, ": ", qc_format_number(d_values), "<br />Density: ", qc_format_number(d$y)),
        type = "scatter", mode = "lines", fill = "tozeroy", fillcolor = "rgba(0,0,0,0.2)",
        line = list(width = 1.8898, color = "rgba(0,0,0,0.2)"), hoverinfo = "text", showlegend = FALSE
    ) %>%
        qc_add_quantiles(values_quants, quants, y_range) %>%
        plotly::layout(
            margin = list(t = 23.3059, r = 7.3059, b = 37.2603, l = 43.105),
            font = list(color = "rgba(0,0,0,1)", size = 14.6119),
            xaxis = qc_x_axis(x, xlab, log_scale),
            yaxis = qc_plot_axis(ylab, y_range, list(tickvals = y_breaks, ticktext = format(y_breaks, trim = TRUE))),
            showlegend = FALSE, hovermode = "closest"
        ) %>%
        collapse_hoverinfo() %>%
        sanitize_plotly_buttons() %>%
        sanitize_plotly_download(globals)

    return(p)
}

# The share of metacells with a value above x (1 - the ECDF of the values)
qc_ecdf <- function(qc_df, field, xlab, ylab, globals, log_scale = FALSE) {
    if (is.null(ylab)) {
        ylab <- "% of metacells > x"
    }
    values <- qc_df[[field]]
    values <- values[!is.na(values)]
    values_quants <- quantile(values, probs = c(0.1, 0.5, 0.9))

    # at the values (at most ~1000 of them)
    v <- sort(unique(values))
    if (length(v) > 1000) {
        v <- unique(stats::quantile(values, seq(0, 1, length.out = 1000), type = 1, names = FALSE))
    }
    y <- 1 - stats::ecdf(values)(v)

    x <- v
    quants <- values_quants
    if (log_scale) {
        x <- log10(x)
        quants <- log10(quants)
    }
    finite <- is.finite(x)
    x <- x[finite]
    y <- y[finite]
    v <- v[finite]
    y_range <- c(-0.05, 1.05)
    y_breaks <- c(0, 0.25, 0.5, 0.75, 1)

    p <- plotly::plot_ly(
        x = x, y = y, text = paste0(xlab, ": ", qc_format_number(v), "<br />Metacells above: ", scales::percent(y, accuracy = 0.1)),
        # (the share at a value holds until the next value)
        type = "scatter", mode = "lines", line = list(width = 1.8898, color = "rgba(0,0,0,1)", shape = "hv"),
        hoverinfo = "text", showlegend = FALSE
    ) %>%
        qc_add_quantiles(values_quants, quants, y_range) %>%
        plotly::layout(
            margin = list(t = 23.3059, r = 7.3059, b = 37.2603, l = 43.105),
            font = list(color = "rgba(0,0,0,1)", size = 14.6119),
            xaxis = qc_x_axis(x, xlab, log_scale),
            yaxis = qc_plot_axis(ylab, y_range, list(tickvals = y_breaks, ticktext = scales::percent(y_breaks))),
            showlegend = FALSE, hovermode = "closest"
        ) %>%
        collapse_hoverinfo() %>%
        sanitize_plotly_buttons() %>%
        sanitize_plotly_download(globals)

    return(p)
}

# A number in the labels and hover texts of the QC plots: big numbers with thousands separators (178,755),
# the others with 3 significant digits (2.77, 0.0153)
qc_format_number <- function(x) {
    ifelse(
        abs(x) >= 1000,
        scales::comma(x, accuracy = 1),
        trimws(formatC(signif(x, 3), digits = 3, format = "fg"))
    )
}

# An axis of the QC plots, styled as ggplot's theme_classic()
qc_plot_axis <- function(title, range, ticks) {
    c(list(
        title = list(text = title, font = list(color = "rgba(0,0,0,1)", size = 14.6119)),
        type = "linear", autorange = FALSE, range = range, tickmode = "array", ticks = "outside",
        tickcolor = "rgba(0,0,0,1)", ticklen = 3.653, tickwidth = 0.6642, tickfont = list(color = "rgba(0,0,0,1)", size = 11.6895),
        showline = TRUE, linecolor = "rgba(0,0,0,1)", linewidth = 0.6642, showgrid = FALSE, zeroline = FALSE
    ), ticks)
}

# The x axis of a QC plot of `x` (the log10 of the values when log_scale) as ggplot's: 5% wider than the
# data, with its breaks and labels
qc_x_axis <- function(x, xlab, log_scale) {
    x_range <- range(x) + c(-0.05, 0.05) * diff(range(x))
    if (log_scale) {
        x_breaks <- scales::breaks_log(base = 10)(10^range(x))
        x_breaks <- x_breaks[log10(x_breaks) >= x_range[1] & log10(x_breaks) <= x_range[2]]
        ticks <- list(tickvals = log10(x_breaks), ticktext = scales::scientific(x_breaks))
    } else {
        x_breaks <- scales::breaks_extended()(range(x))
        ticks <- list(tickvals = x_breaks, ticktext = format(x_breaks, trim = TRUE))
    }
    qc_plot_axis(xlab, x_range, ticks)
}

# Dashed lines at the bottom 10%, the median and the top 10% of the values (values_quants), at `quants` (the
# same in the units of the x axis, e.g. their log10), with their labels. Each label is beside its line, so that
# the line doesn't strike through it, and away from the other lines, so that the labels don't overlap: the
# bottom 10% to the left of its line, the top 10% to the right, and the median on the side with more room.
qc_add_quantiles <- function(p, values_quants, quants, y_range) {
    labels <- paste0(c("bottom 10%", "median", "top 10%"), ": ", qc_format_number(values_quants))
    line_colors <- c("rgba(0,0,0,1)", "rgba(255,0,0,1)", "rgba(0,0,0,1)")
    label_colors <- c("rgba(0,0,0,1)", "rgba(200,0,0,1)", "rgba(0,0,0,1)")
    sides <- c("left", if (quants[3] - quants[2] >= quants[2] - quants[1]) "right" else "left", "right")

    # (points along the lines, so that a line shows its label when hovered anywhere and not only at its ends)
    line_y <- seq(y_range[1], y_range[2], length.out = 20)
    for (i in seq_along(quants)) {
        p <- plotly::add_trace(p,
            x = rep(quants[[i]], length(line_y)), y = line_y, text = labels[i],
            type = "scatter", mode = "lines", line = list(width = 1.8898, color = line_colors[i], dash = "dash"),
            hoverinfo = "text", showlegend = FALSE, fill = "none", inherit = FALSE
        )
    }

    annotations <- purrr::map(seq_along(quants), function(i) {
        list(
            x = quants[[i]], xref = "x", xanchor = if (sides[i] == "left") "right" else "left", xshift = if (sides[i] == "left") -2 else 2,
            y = 1, yref = "paper", yanchor = "top",
            text = labels[i], textangle = -90, showarrow = FALSE,
            font = list(color = label_colors[i], size = 10), bgcolor = "rgba(255,255,255,0.85)"
        )
    })
    plotly::layout(p, annotations = annotations)
}
