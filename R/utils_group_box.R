#' Render a DT::datatable for a metacell group
#'
#' Creates a DT::renderDataTable output that displays metacells with cell-type-colored
#' backgrounds. Used by both mod_mc_mc (groupA/groupB) and mod_query (single group).
#'
#' @param output shiny output object
#' @param output_id character id for the DT output (e.g. "groupA_table")
#' @param group_reactive reactiveVal holding the current group vector
#' @param metacell_types reactive returning metacell-to-cell_type mapping
#' @param cell_type_colors reactive returning cell_type-to-color mapping
#'
#' @noRd
render_group_table <- function(output, output_id, group_reactive, metacell_types, cell_type_colors) {
    output[[output_id]] <- DT::renderDataTable(
        {
            req(metacell_types())
            req(cell_type_colors())
            req(group_reactive())
            DT::datatable(
                tibble(metacell = group_reactive()) %>%
                    left_join(metacell_types() %>% select(metacell, cell_type), by = "metacell"),
                escape = FALSE,
                rownames = FALSE,
                colnames = "",
                filter = "none",
                options = list(
                    dom = "t",
                    paging = FALSE,
                    language = list(emptyTable = "Please select metacells"),
                    columnDefs = list(list(visible = FALSE, targets = c(1)))
                )
            ) %>%
                DT::formatStyle(
                    "metacell", "cell_type",
                    backgroundColor = DT::styleEqual(
                        cell_type_colors()$cell_type,
                        col2hex(cell_type_colors()$color)
                    )
                )
        },
        server = FALSE
    )
}

#' Render the UI box for a metacell group
#'
#' Creates a renderUI output containing a generic_box with Reset/Remove/Paste buttons
#' and a DT::dataTableOutput.
#'
#' @param output shiny output object
#' @param output_id character id for the UI output (e.g. "groupA_box")
#' @param ns namespace function from the module session
#' @param group_name character display name for the group (e.g. "Group A", "Group")
#' @param table_output_id character id for the DT table inside the box (e.g. "groupA_table")
#' @param reset_id character input id for the reset button (e.g. "reset_groupA")
#' @param remove_id character input id for the remove button (e.g. "remove_groupA_metacells")
#' @param paste_id character input id for the paste button (e.g. "paste_groupA_metacells")
#' @param mode_condition expression string for the mode check (e.g. "Groups" or "Group")
#' @param input shiny input object
#'
#' @noRd
render_group_box <- function(output, output_id, ns, group_name, table_output_id,
                             reset_id, remove_id, paste_id, mode_condition, input) {
    output[[output_id]] <- renderUI({
        req(input$mode == mode_condition)
        generic_box(
            id = ns(paste0(gsub("_box$", "", output_id), "_box_1")),
            title = paste(group_name, "metacells"),
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            closable = FALSE,
            width = 12,
            actionButton(ns(reset_id), "Reset"),
            actionButton(ns(remove_id), "Remove"),
            actionButton(ns(paste_id), "Paste"),
            shinycssloaders::withSpinner(
                DT::dataTableOutput(ns(table_output_id))
            )
        )
    })
}

#' Set up observers for group management (add, remove, paste, reset)
#'
#' Creates observeEvent handlers for the four standard group management actions:
#' adding a metacell from a picker, removing selected rows, pasting from clipboard,
#' and resetting the group.
#'
#' @param input shiny input object
#' @param group_reactive reactiveVal holding the current group vector
#' @param add_id character input id for the add button (e.g. "add_metacell_to_groupA")
#' @param remove_id character input id for the remove button (e.g. "remove_groupA_metacells")
#' @param paste_id character input id for the paste button (e.g. "paste_groupA_metacells")
#' @param reset_id character input id for the reset button (e.g. "reset_groupA")
#' @param table_id character input id for the DT table (e.g. "groupA_table")
#' @param metacell_input_id character input id for the metacell picker (e.g. "metacell")
#' @param state Domain-scoped reactiveValues (session_ui / tab_state / selection / manifold_state).
#'
#' @noRd
setup_group_observers <- function(input, group_reactive, add_id, remove_id, paste_id,
                                  reset_id, table_id, metacell_input_id, state) {
    observeEvent(input[[add_id]], {
        add_to_group(group_reactive, input[[metacell_input_id]])
    })

    observeEvent(input[[remove_id]], {
        rows <- input[[paste0(table_id, "_rows_selected")]]
        req(rows)
        req(length(rows) > 0)
        group_reactive(group_reactive()[-rows])
    })

    observeEvent(input[[paste_id]], {
        metacells <- state$selection$clipboard
        group_reactive(unique(c(group_reactive(), metacells)))
    })

    observeEvent(input[[reset_id]], {
        group_reactive(NULL)
    })
}

#' Add metacells to a group reactiveVal
#'
#' Appends metacells to a group, handling the NULL initial case. Ensures uniqueness.
#'
#' @param group_reactive reactiveVal holding the current group vector
#' @param metacells character vector of metacell ids to add
#'
#' @noRd
add_to_group <- function(group_reactive, metacells) {
    if (is.null(group_reactive())) {
        group_reactive(metacells)
    } else {
        group_reactive(unique(c(group_reactive(), metacells)))
    }
}
