metacell_click_observer <- function(source_id, session, input_id = "metacell1") {
    plotly_click_observer(source_id, session, notification_prefix = "Selected Metacell #", input_id = input_id, update_function = updateSelectInput)
}

sample_click_observer <- function(source_id, session, input_id = "samp1") {
    plotly_click_observer(source_id, session, notification_prefix = "Selected Sample #", input_id = input_id, update_function = shinyWidgets::updatePickerInput)
}

plotly_click_observer <- function(source_id, session, input_id, notification_prefix = NULL, update_function = updateSelectInput) {
    # ignoreInit: a tab starts when it is first opened, and must not apply a click made in another tab before that
    # (e.g. a metacell id into a picker that shows cell types)
    observeEvent(plotly::event_data("plotly_click", source = source_id), ignoreInit = TRUE, {
        el <- plotly::event_data("plotly_click", source = source_id)
        selected <- plotly_event_ids(el)
        update_function(session, input_id, selected = selected)
        if (!is.null(notification_prefix)) {
            showNotification(glue("{notification_prefix}{selected}"))
        }
    })
}

# The ids (metacells, genes...) of the points of a plotly click or selection event. Plots with a
# hovertemplate (set_hovertemplate()) keep the id first in each point's customdata, followed by
# its hover values.
plotly_event_ids <- function(el) {
    ids <- el$customdata
    if (is.list(ids)) {
        ids <- purrr::map_chr(ids, 1, .default = NA)
    }
    ids
}
