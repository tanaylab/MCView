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
        if (
            (!is.null(input$color_by_type) && !is.null(input$color_by_var) && input$color_by_type == "Metadata" && input$color_by_var == "Clipboard") ||
                (!is.null(input$filter_by_clipboard_scatter) && input$filter_by_clipboard_scatter)) {
            return(globals$clipboard)
        } else {
            return(FALSE)
        }
    })
}

clipboard_reactives <- function(dataset, input, output, session, metacell_types, cell_type_colors, gene_modules, globals) {
    observeEvent(
        input$clipboard_modal,
        showModal(modalDialog(
            title = "Clipboard",
            actionButton("clear_clipboard", "Clear clipboard"),
            actionButton("delete_clipboard_row", "Remove selected"),
            downloadButton("download_clipboard", "Download"),
            br(),
            br(),
            shinycssloaders::withSpinner(
                DT::dataTableOutput("clipboard_table")
            ),
            br(),
            br(),
            uiOutput("metacell_selector"),
            actionButton("add_to_clipboard", "Add to clipboard"),
            easyClose = TRUE
        ))
    )

    metacell_names <- metacell_names_reactive(dataset)
    metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)

    output$metacell_selector <- renderUI({
        cell_types_hex <- col2hex(metacell_colors())
        shinyWidgets::pickerInput("metacells_to_add", "",
            choices = metacell_names(),
            multiple = TRUE,
            options = shinyWidgets::pickerOptions(`actions-box` = TRUE, liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "contains", dropupAuto = FALSE),
            choicesOpt = list(
                style = paste0("color: ", cell_types_hex, ";")
            )
        )
    })

    observeEvent(input$add_to_clipboard, {
        globals$clipboard <- unique(c(globals$clipboard, input$metacells_to_add))
    })

    observeEvent(input$clear_clipboard, {
        globals$clipboard <- character(0)
    })

    observeEvent(input$delete_clipboard_row, {
        globals$clipboard <- globals$clipboard[-input$clipboard_table_rows_selected]
    })

    output$clipboard_table <- DT::renderDataTable(
        metacell_types() %>% filter(metacell %in% globals$clipboard) %>% select(metacell, cell_type),
        escape = FALSE,
        server = FALSE,
        rownames = FALSE,
        filter = "top",
        options = list(
            dom = "t",
            paging = FALSE,
            language = list(emptyTable = "No metacells in clipboard")
        )
    )

    output$download_clipboard <- downloadHandler(
        filename = function() {
            paste("metacell_clipboard-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            fwrite(
                metacell_types() %>% filter(metacell %in% globals$clipboard) %>% select(metacell, cell_type),
                file
            )
        }
    )
}

#' Create a reactive copy button using rclipboard
#'
#' @param ns Namespace function
#' @param id Button ID
#' @param data_reactive Reactive expression that returns the data to copy
#' @param label Button label
#' @param style Button styling
#' @param tooltip Tooltip text
#' @param disabled_label Label when button is disabled
#' @return UI output for the copy button
#' @importFrom rclipboard rclipButton
clipboard_copy_button_ui <- function(ns, id, data_reactive, label = "Copy to Clipboard",
                                     style = "background-color: #17a2b8; color: white; border: none;",
                                     tooltip = "Copy to system clipboard",
                                     disabled_label = NULL) {
    renderUI({
        data <- data_reactive()
        if (!is.null(data) && length(data) > 0) {
            clip_text <- if (is.character(data)) {
                paste(data, collapse = "\n")
            } else {
                paste(as.character(data), collapse = "\n")
            }

            rclipboard::rclipButton(
                inputId = ns(id),
                label = label,
                clipText = clip_text,
                icon = icon("copy"),
                tooltip = tooltip,
                placement = "top",
                style = style
            )
        } else {
            actionButton(ns(paste0(id, "_disabled")),
                disabled_label %||% label,
                icon = icon("copy"),
                style = style,
                disabled = TRUE
            )
        }
    })
}

#' Create server logic for copy button notifications
#'
#' @param input Shiny input
#' @param id Button ID
#' @param data_reactive Reactive expression that returns the data
#' @param globals Global reactive values (optional, for internal clipboard)
#' @param message_template Template for success message with {count} placeholder
#' @return Observer for button clicks
clipboard_copy_button_server <- function(input, id, data_reactive, globals = NULL,
                                         message_template = "Copied {count} items to clipboard") {
    observeEvent(input[[id]], {
        data <- data_reactive()
        if (!is.null(data) && length(data) > 0) {
            # Update internal clipboard if globals provided
            if (!is.null(globals) && !is.null(globals$clipboard)) {
                globals$clipboard <- as.character(data)
            }

            # Show notification
            message <- glue::glue(message_template, count = length(data))
            showNotification(message, type = "default")
        } else {
            showNotification("No data to copy", type = "warning")
        }
    })
}
