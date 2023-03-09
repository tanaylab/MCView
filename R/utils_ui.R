generic_column <- function(width, ..., offset = 0) {
    column(
        width = width,
        offset = offset,
        ...
    )
}

generic_box <- function(...) {
    shinyjqui::jqui_resizable(
        shinydashboardPlus::box(
            ...
        )
    )
}
