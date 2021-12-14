metacell_click_observer <- function(source_id, session, input_id = "metacell1", update_function = updateSelectInput) {
    observeEvent(plotly::event_data("plotly_click", source = source_id), {
        el <- plotly::event_data("plotly_click", source = source_id)
        metacell <- el$customdata
        update_function(session, input_id, selected = metacell)
        showNotification(glue("Selected Metacell #{metacell}"))
    })
}

sample_click_observer <- function(source_id, session, input_id = "metacell1", update_function = shinyWidgets::updatePickerInput) {
    observeEvent(plotly::event_data("plotly_click", source = source_id), {
        el <- plotly::event_data("plotly_click", source = source_id)
        sample <- el$customdata
        update_function(session, input_id, selected = sample)
        showNotification(glue("Selected Sample #{sample}"))
    })
}
