# This ugly code defines the titles for plotly subplots
get_plotly_subplot_title <- function(.x) {
    list(
        text = .x,
        font = list(size = 14),
        xref = "paper",
        yref = "paper",
        yanchor = "bottom",
        xanchor = "center",
        align = "center",
        x = 0.5,
        y = 1,
        showarrow = FALSE
    )
}


# Remove unnecessary atributes from plots that were created with ggplotly before converting to WebGL.
# See: see: https://stackoverflow.com/questions/62911241/ggplot2-to-plotly-webgl-warning-suppression-in-r-shiny
sanitize_for_WebGL <- function(p) {
    for (i in 1:length(p$x$data)) {
        if ("hoveron" %in% names(p$x$data[[i]])) {
            p$x$data[[i]]$hoveron <- NULL
        }
    }

    return(p)
}

# See https://github.com/plotly/plotly.js/issues/4999
arrange_2d_proj_tooltip <- function(fig) {
    if (fig$x$data[[1]]$mode == "lines") {
        fig$x$data[[1]]$hoverinfo <- "none"
        point_data <- fig$x$data[2:length(fig$x$data)]
        point_data_no_legend <- purrr::map(point_data, ~ {
            .x$showlegend <- FALSE
            .x
        })
        fig$x$data <- c(point_data, fig$x$data[1], point_data_no_legend)
    }
    return(fig)
}

rm_plotly_grid <- function(fig) {
    fig %>%
        plotly::layout(xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE))
}


sanitize_plotly_buttons <- function(p,
                                    buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"),
                                    ...) {
    p %>% plotly::config(displaylogo = FALSE, modeBarButtonsToRemove = buttons, ...)
}

plotly_text_plot <- function(text) {
    p <- ggplot() +
        annotate("text",
            x = 1,
            y = 1,
            label = text
        ) +
        theme_void()
    return(plotly:::ggplotly(p))
}
