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
