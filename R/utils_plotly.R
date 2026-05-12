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
    p %>% plotly::config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = buttons,
        ...
    )
}

sanitize_plotly_download <- function(p, state) {
    to_image_options <- list(format = state$session_ui$plotly_format, width = state$session_ui$plotly_width, height = state$session_ui$plotly_height, scale = state$session_ui$plotly_scale)
    p %>% plotly::config(
        toImageButtonOptions = to_image_options
    )
}

#' Prepare a ggplot for plotly with WebGL and standard sanitization
#'
#' Consolidates the repeated pipeline of:
#'   plotly::ggplotly() %>% sanitize_for_WebGL() %>% toWebGL() %>%
#'   sanitize_plotly_buttons() %>% sanitize_plotly_download()
#'
#' @param p ggplot object
#' @param tooltip Tooltip text column name (default: "tooltip_text"), or NULL
#' @param source Plotly source name for event handling
#' @param buttons Buttons to remove from mode bar
#' @param state Domain-scoped reactiveValues (session_ui / tab_state / selection / manifold_state).
#' @param hide_legend Whether to hide the legend (default: FALSE)
#' @return Sanitized plotly object
#' @export
prepare_plotly_scatter <- function(p, tooltip = "tooltip_text", source = NULL,
                                   buttons = c(
                                       "select2d", "lasso2d", "hoverClosestCartesian",
                                       "hoverCompareCartesian", "toggleSpikelines"
                                   ),
                                   state = NULL, hide_legend = FALSE) {
    if (is.null(tooltip)) {
        fig <- plotly::ggplotly(p, source = source)
    } else {
        fig <- plotly::ggplotly(p, tooltip = tooltip, source = source)
    }

    if (hide_legend) {
        fig <- fig %>% plotly::hide_legend()
    }

    fig <- fig %>%
        sanitize_for_WebGL() %>%
        plotly::toWebGL() %>%
        plotly::partial_bundle() %>%
        sanitize_plotly_buttons(buttons = buttons)

    if (!is.null(state)) {
        fig <- fig %>% sanitize_plotly_download(state)
    }

    fig %>%
        plotly::event_register("plotly_click")
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

#' Pre-warm plotly's internal partial_bundle cache.
#'
#' First call to `plotly::partial_bundle()` in an R process parses the
#' ~3 MB plotly.js JSON; subsequent calls reuse the cached result. We pay
#' this cost once at session start so the first user-visible plot does
#' not.
#'
#' Idempotent at process level: the `mcv_get(".partial_bundle_warmed")`
#' guard avoids the function-call overhead on subsequent Shiny sessions
#' that share the same R process.
#'
#' @return Invisibly NULL.
#' @noRd
prewarm_plotly_bundle <- function() {
    if (isTRUE(mcv_get(".partial_bundle_warmed"))) {
        return(invisible(NULL))
    }
    tryCatch(
        invisible(plotly::partial_bundle(plotly::plot_ly(type = "scatter", mode = "markers"))),
        error = function(e) NULL  # opportunistic; failure is acceptable
    )
    # Set the flag unconditionally so a transient failure doesn't make
    # every subsequent session retry the (potentially slow) pre-warm.
    mcv_set(".partial_bundle_warmed", TRUE)
    invisible(NULL)
}
