#' Plot a continuous color bar with values
#'
#' @param vals to plot
#' @param labels labels of vals
#' @param cols colors of the values (same length as vals)
#' @param title of legend
#'
#' @noRd
plot_color_bar <- function(vals, labels, cols, title = "") {
    plot.new()
    plot.window(xlim = c(0, 100), ylim = c(-30, length(cols) + 30))
    rect(45, 1:length(cols), 17, 1:length(cols) + 1, border = NA, col = cols)
    rect(45, 1, 17, length(cols) + 1, col = NA, border = "black")

    text(50, vals, labels = labels, pos = 4)
    text(7, length(cols) / 2 + 1, labels = title, srt = 90, cex = 1)
}

#' Wait until a plot output has a size
#'
#' Use as the first key of bindCache() on renderPlot(). shiny opens the PNG device of a plot with
#' the output's size before drawing, and when the width is zero (or not reported yet by the
#' browser) it fails and leaves that device open (see shiny:::startPNG). Every ggplotly() in the
#' process then fails with "invalid 'width' argument", for all users, until the app restarts.
#'
#' @noRd
req_plot_size <- function(session = getDefaultReactiveDomain()) {
    output <- paste0("output_", getCurrentOutputInfo()$name)
    req(session$clientData[[paste0(output, "_width")]] > 0, session$clientData[[paste0(output, "_height")]] > 0)
}

#' The sizePolicy of bindCache() on renderPlot(): shiny's default, but never zero
#'
#' bindCache() takes the plot's size from an observer of the output's size, so right after a
#' zero width becomes positive, req_plot_size() can pass while the size is still zero.
#'
#' @noRd
plot_size_policy <- function(dims) {
    sizeGrowthRatio(width = 400, height = 400, growthRate = 1.2)(pmax(dims, 1))
}

#' Convert R color to hex (for css usage)
#'
#' @param x R color name
#'
#' @noRd
col2hex <- function(x) {
    grDevices::rgb(t(grDevices::col2rgb(x) / 255))
}
