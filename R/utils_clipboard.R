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
