# utils_annotate.R - Helpers used by the Annotate module
#
# Split from R/mod_annotate.R (2026-05-01). Pure helpers that operate on
# the cell-type-colors reactiveVal + plotly click/selection events. Kept
# at top-level scope so mod_annotate_server stays focused on the moduleServer
# body.



#' Move selected cell types up or down in the ordering
#'
#' @param cell_type_colors reactive value holding cell type colors data frame
#' @param selected_rows sorted vector of selected row indices
#' @param direction -1 for up, +1 for down
#' @noRd
move_cell_type <- function(cell_type_colors, selected_rows, direction) {
    new_data <- cell_type_colors()
    total_rows <- nrow(new_data)

    # Check boundary: can't move up past top or down past bottom
    if (direction == -1 && min(selected_rows) <= 1) {
        return(invisible(NULL))
    }
    if (direction == 1 && max(selected_rows) >= total_rows) {
        return(invisible(NULL))
    }

    # Find the edge order value and the adjacent row to swap with
    if (direction == -1) {
        edge_order <- min(new_data$order[selected_rows])
    } else {
        edge_order <- max(new_data$order[selected_rows])
    }
    neighbor_row <- which(new_data$order == (edge_order + direction))

    if (length(neighbor_row) != 1) {
        return(invisible(NULL))
    }

    selected_cell_types_names <- new_data$cell_type[selected_rows]

    new_data <- new_data %>%
        mutate(temp_order = order)

    # Move the neighbor to the opposite side of the selection
    new_data$temp_order[neighbor_row] <- edge_order - direction * (length(selected_rows) - 1)

    # Shift all selected rows in the direction
    new_data$temp_order[selected_rows] <- new_data$temp_order[selected_rows] + direction

    new_data <- new_data %>%
        mutate(order = rank(temp_order, ties.method = "first")) %>%
        select(-temp_order) %>%
        arrange(order)

    cell_type_colors(new_data)

    new_indices <- which(new_data$cell_type %in% selected_cell_types_names)
    shinyjs::delay(100, {
        DT::selectRows(DT::dataTableProxy("cell_type_table"), new_indices)
    })
}

observe_mc_click_event <- function(source, input, session, cell_type_colors, metacell_types, selected_metacell_types) {
    observeEvent(plotly::event_data("plotly_click", source = source), {
        el <- plotly::event_data("plotly_click", source = source)

        selected_metacell <- el$customdata

        new_selected_annot <- metacell_types() %>% filter(metacell == selected_metacell)

        selected_metacell_types(
            bind_rows(
                selected_metacell_types(),
                new_selected_annot
            ) %>% distinct(metacell, cell_type)
        )

        shinyWidgets::updatePickerInput(session, "metacell1", selected = selected_metacell)
    })
}

observer_mc_select_event <- function(source, input, cell_type_colors, metacell_types, selected_metacell_types) {
    observeEvent(plotly::event_data("plotly_selected", source = source), {
        el <- plotly::event_data("plotly_selected", source = source)

        selected_metacells <- unique(el$customdata)

        new_selected_annot <- metacell_types() %>% filter(metacell %in% selected_metacells)
        if (!is.null(input$add_to_selection) && input$add_to_selection) {
            selected_metacell_types(
                bind_rows(
                    selected_metacell_types(),
                    new_selected_annot
                ) %>% distinct(metacell, cell_type)
            )
        } else {
            selected_metacell_types(new_selected_annot %>% distinct(metacell, cell_type))
        }
    })
}
