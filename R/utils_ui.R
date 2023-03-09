generic_column <- function(width, ..., offset = 0) {
    column(
        width = width,
        offset = offset,
        ...
    )
}

movable_box <- function(...) {
    shinyjqui::jqui_draggable(
        shinyjqui::jqui_resizable(
            shinydashboardPlus::box(
                ...
            )
        )
    )
}
