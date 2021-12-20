metacell_click_observer <- function(source_id, session, input_id = "metacell1") {
    plotly_click_observer(source_id, session, notification_prefix = "Selected Metacell #", input_id = input_id, update_function = updateSelectInput)
}

sample_click_observer <- function(source_id, session, input_id = "samp1") {
    plotly_click_observer(source_id, session, notification_prefix = "Selected Sample #", input_id = input_id, update_function = shinyWidgets::updatePickerInput)
}

plotly_click_observer <- function(source_id, session, input_id, notification_prefix = NULL, update_function = updateSelectInput) {
    observeEvent(plotly::event_data("plotly_click", source = source_id), {
        el <- plotly::event_data("plotly_click", source = source_id)
        selected <- el$customdata
        update_function(session, input_id, selected = selected)
        if (!is.null(notification_prefix)) {
            showNotification(glue("{notification_prefix}{selected}"))
        }
    })
}
