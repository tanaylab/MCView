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

# plotly_build() repeats a trace's hoverinfo for every point (with NA for NA points). Send it
# once per trace. Needed for hoverinfo = "skip" to work, and smaller.
collapse_hoverinfo <- function(p) {
    p <- plotly::plotly_build(p)
    p$x$data <- purrr::map(p$x$data, function(tr) {
        h <- unique(tr$hoverinfo[!is.na(tr$hoverinfo)])
        if (length(tr$hoverinfo) > 1 && length(h) == 1) {
            tr$hoverinfo <- h
        }
        tr
    })
    p
}

# ggplotly writes the hover text of every point as one string ("label: value<br />label: value..."),
# so the labels are sent once per point. Send them once per trace instead, as a hovertemplate.
# The points' customdata must be df[[id]] (unique values). It is replaced by the point's values of
# c(id, fields) as strings, which `template` shows with %{customdata[0]}, %{customdata[1]}, ...
# (%{x:.3~r}, %{y} and %{fullData.name} for the trace name can be used too).
set_hovertemplate <- function(fig, df, id, fields, template) {
    fig <- plotly::plotly_build(fig)
    values <- do.call(cbind, lapply(df[c(id, fields)], as.character))
    fig$x$data <- purrr::map(fig$x$data, function(tr) {
        if (is.null(tr$customdata)) {
            return(tr)
        }
        tr$customdata <- values[match(unlist(tr$customdata), df[[id]]), , drop = FALSE]
        tr$hovertemplate <- template
        tr$text <- NULL
        tr$hoverinfo <- NULL
        tr
    })
    fig
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
    p %>% plotly::config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = buttons,
        ...
    )
}

sanitize_plotly_download <- function(p, globals) {
    to_image_options <- list(format = globals$plotly_format, width = globals$plotly_width, height = globals$plotly_height, scale = globals$plotly_scale)
    # (done here because all our plotly outputs pass through this function)
    attr(p$x, "TOJSON_FUNC") <- plotly_to_json
    # By default plotly sends every hover, relayout and redraw of a plot to the server as an input,
    # and shiny goes over all the outputs of the session for each input it receives. Send only
    # the events that MCView observes (click, selection and restyle).
    for (event in c("plotly_hover", "plotly_relayout", "plotly_brushed", "plotly_brushing", "plotly_clickannotation", "plotly_doubleclick", "plotly_afterplot", "plotly_sunburstclick")) {
        p <- plotly::event_unregister(p, event)
    }
    p %>% plotly::config(
        toImageButtonOptions = to_image_options
    )
}

# plotly's own serializer (plotly:::to_JSON) writes numbers with 15-17 significant digits. 6 are
# plenty for drawing (hover texts are formatted in R), make the plots' JSON ~10% smaller and are
# faster to write.
plotly_to_json <- function(x, ...) {
    jsonlite::toJSON(x, digits = I(6), auto_unbox = TRUE, force = TRUE, null = "null", na = "null", time_format = "%Y-%m-%d %H:%M:%OS6", ...)
}

plotly_text_plot <- function(text) {
    p <- ggplot() +
        annotate("text",
            x = 1,
            y = 1,
            label = text
        ) +
        theme_void()
    return(plotly::ggplotly(p))
}
