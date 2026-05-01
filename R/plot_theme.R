#' MCView ggplot2 theme
#'
#' A clean, minimal ggplot theme inheriting from `theme_minimal()`. Used by
#' every `plot_*()` function. For HTML rendering (via plotly), typography
#' aligns with `inst/app/www/mcview.css` because plotly inherits browser
#' fonts; for bitmap export the device default font is used.
#'
#' @param base_size base font size in pt. Default 12.
#' @param base_family base font family. Default `""` (device default).
#' @return A ggplot2 theme.
#' @export
theme_mcview <- function(base_size = 12, base_family = "") {
    half_line <- base_size / 2
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = ggplot2::rel(1.15), face = "bold", hjust = 0,
                margin = ggplot2::margin(b = half_line)
            ),
            plot.subtitle = ggplot2::element_text(
                size = ggplot2::rel(0.95), color = "#5b6770", hjust = 0,
                margin = ggplot2::margin(b = half_line)
            ),
            axis.title = ggplot2::element_text(size = ggplot2::rel(0.95), color = "#1f2933"),
            axis.text = ggplot2::element_text(size = ggplot2::rel(0.85), color = "#5b6770"),
            axis.line = ggplot2::element_line(color = "#1f2933", linewidth = 0.3),
            axis.ticks = ggplot2::element_line(color = "#1f2933", linewidth = 0.3),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(color = "#e5e9ee", linewidth = 0.3),
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(fill = "#f7f9fb", color = NA),
            strip.text = ggplot2::element_text(face = "bold", color = "#1f2933"),
            legend.position = "bottom",
            legend.key.size = grid::unit(0.8, "lines"),
            legend.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
            legend.title = ggplot2::element_text(size = ggplot2::rel(0.9), face = "bold"),
            plot.background = ggplot2::element_rect(fill = "white", color = NA),
            panel.background = ggplot2::element_rect(fill = "white", color = NA)
        )
}
