clipboard_changed_2d_reactive <- function(input, globals) {
    # We use this reactive in order to invalidate the cache only when needed
    reactive({
        if (is.null(input$color_proj) || is.null(input$color_proj_metadata) || input$color_proj != "metadata" || input$color_proj_metadata != "Clipboard") {
            return(FALSE)
        } else {
            return(globals$clipboard)
        }
    })
}

clipboard_changed_scatter_reactive <- function(input, globals) {
    reactive({
        if (is.null(input$color_by_type) || is.null(input$color_by_var) || input$color_by_type != "Metadata" || input$color_by_var != "Clipboard") {
            return(FALSE)
        } else {
            return(globals$clipboard)
        }
    })
}
