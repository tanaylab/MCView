generic_column <- function(width, ..., offset = 0) {
    column(
        width = width,
        offset = offset,
        ...
    )
}

# A resizable box. shinyjqui::jqui_resizable() would add an iframe next to each box to run its
# init script, which made the page load slower: app_ui() makes all the boxes with the
# "jqui-resizable" class resizable at once instead.
generic_box <- function(...) {
    htmltools::tagAppendAttributes(
        shinydashboardPlus::box(
            ...
        ),
        class = "jqui-resizable"
    )
}
