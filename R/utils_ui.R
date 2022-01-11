resizable_column <- function(width, ..., offset = 0, operation = c("enable", "disable", "destroy", "save", "load"), options = NULL) {
    shinyjqui::jqui_resizable(
        column(
            width = width,
            offset = offset,
            ...
        ),
        operation = operation,
        options = options
    )
}
