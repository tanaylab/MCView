#' gene_modules UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_gene_modules_ui <- function(id) {
    ns <- NS(id)

    tagList(
        fluidRow(
            generic_column(
                width = 9,
                style = "padding-right:0px;",
                heatmap_box(ns("gene_modules_heatmap"), "Gene modules Heatmap", legend_width = 2),
                fluidRow(
                    generic_column(
                        width = 6,
                        scatter_box(ns, "gene_gene_box", x_selected = "Gene module", y_selected = "Gene")
                    ),
                    generic_column(
                        width = 6,
                        diff_expr_box(
                            ns,
                            "mc_mc_box",
                            "Differential expression",
                            c("MCs", "Types"),
                            "Types",
                            uiOutput(ns("add_selected_genes_button"))
                        )
                    )
                )
            ),
            generic_column(
                width = 3,
                style = "padding-right:0px; padding-left:0px;",
                generic_box(
                    id = ns("gene_modules_box"),
                    title = "Gene modules",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,
                    sidebar = shinydashboardPlus::boxSidebar(
                        startOpen = FALSE,
                        width = 80,
                        id = ns("gene_modules_sidebar"),
                        span(
                            actionButton(ns("remove_gene_module_modal"), "Remove module", style = "margin-bottom: 15px;"),
                            actionButton(ns("add_gene_module_modal"), "New module", style = "margin-bottom: 15px; margin-left: 5px;"),
                            downloadButton(ns("gene_modules_download"), "Export", style = "margin-bottom: 15px;"),
                            fileInput(ns("gene_modules_fn"),
                                label = NULL,
                                buttonLabel = "Load",
                                multiple = FALSE,
                                accept =
                                    c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        "text/tab-separated-values",
                                        ".csv",
                                        ".tsv"
                                    )
                            )
                        )
                    ),
                    uiOutput(ns("gene_module_selector")),
                    span(
                        actionButton(ns("update_module"), "Update", style = "align-items: center;"),
                        actionButton(ns("clear_genes"), "Clear", style = "align-items: center;"),
                        actionButton(ns("reset_genes"), "Reset", style = "align-items: center;"),
                        actionButton(ns("reset_gene_modules"), "Reset all", style = "align-items: left;"),
                    ),
                    br(),
                    br(),
                    actionButton(ns("remove_genes"), "Remove", style = "align-items: center;"),
                    actionButton(ns("move_selected_genes_modal"), "Move", style = "align-items: left;"),
                    br(),
                    br(),
                    shinycssloaders::withSpinner(
                        DT::dataTableOutput(ns("genes_table"))
                    )
                )
            )
        )
    )
}


#' gene_modules sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_gene_modules_sidebar_ui <- function(id) {
    ns <- NS(id)
    ns_heatmap <- NS(ns("gene_modules_heatmap"))
    ns_genes <- NS(ns("genes_heatmap"))
    tagList(
        list(
            uiOutput(ns("add_genes_sidebar")),
            tags$hr(),
            uiOutput(ns_heatmap("reset_zoom_ui")),
            shinyWidgets::switchInput(
                inputId = ns_heatmap("show_genes"),
                label = "Show genes",
                size = "small",
                labelWidth = "100px",
                value = TRUE
            ),
            uiOutput(ns("shown_gene_modules_ui")),
            uiOutput(ns_heatmap("cell_type_list")),
            uiOutput(ns_heatmap("metadata_list")),
            tags$hr(),
            uiOutput(ns("top_correlated_gene_selector")),
            uiOutput(ns("top_correlated_ui"))
        )
    )
}

#' gene_modules Server Function
#'
#' @noRd
mod_gene_modules_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            shown_gene_modules <- reactiveVal()
            genes <- reactiveVal() # genes to show below the gene modules
            lfp_range <- reactiveVal()
            selected_module <- reactiveVal()
            selected_genes <- reactiveVal() # selected genes at the diff. expr plot

            scatter_selectors(ns, dataset, output, globals)

            output$shown_gene_modules_ui <- gene_modules_selector(
                dataset,
                gene_modules,
                ns,
                label = "Shown gene modules",
                id = "shown_gene_modules",
                selected = "all"
            )

            observe({
                req(input$shown_gene_modules)
                shown_gene_modules(as.character(input$shown_gene_modules[input$shown_gene_modules %in% gene_modules()$module]))
            })

            output$top_correlated_gene_selector <- renderUI({
                # gene_choices <- gene_names(dataset())
                tagList(
                    shinyWidgets::virtualSelectInput(ns("top_correlated_gene"),
                        label = "Top correlated to:",
                        choices = genes(),
                        selected = c(),
                        multiple = FALSE,
                        search = TRUE,
                        markSearchResults = TRUE,
                        searchByStartsWith = TRUE
                    )
                )
            })

            output$top_correlated_ui <- renderUI({
                req(input$top_correlated_gene)
                top_correlated_selector_multiple_genes(input, output, session, dataset, ns, "selected_top_genes", "", gene = input$top_correlated_gene, action_id = "add_genes_from_top_cor_list", action_label = "Add to gene module")
            })

            heatmap_reactives("gene_modules_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, shown_gene_modules, lfp_range, "Gene modules", genes = genes, highlighted_genes = selected_module)

            mod_gene_module_controllers(ns, dataset, input, output, session, gene_modules, genes, selected_module, selected_genes, globals)

            # Scatter plot
            selected_cell_types <- reactiveVal(NULL)
            scatter_box_outputs(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, ns, selected_cell_types = selected_cell_types, plotly_source = "gene_modules_md_md_plot")

            # Diff. expression
            diff_expr_outputs(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, ns, source_suffix = "_gene_modules", dragmode = "select", plotly_buttons = c("hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines"))

            output$add_selected_genes_button <- renderUI({
                actionButton(ns("add_selected_genes"), "Add to gene module")
            })

            selected_genes_event_observer("mc_mc_plot_gene_modules", selected_genes)
            selected_genes_event_observer("ct_ct_plot_gene_modules", selected_genes)
            observe({
                req(selected_genes)
                shinyjs::toggle(id = "add_selected_genes_button", condition = length(selected_genes()) > 0)
            })
        }
    )
}

selected_genes_event_observer <- function(source, selected_genes) {
    observeEvent(plotly::event_data("plotly_selected", source = source), {
        el <- plotly::event_data("plotly_selected", source = source)
        selected_genes(unique(el$customdata))
    })

    observeEvent(plotly::event_data("plotly_deselect", source = source), {
        selected_genes(c())
    })
}

mod_gene_module_controllers <- function(ns, dataset, input, output, session, gene_modules, genes, selected_module, selected_genes, globals) {
    values <- reactiveValues(file_status = NULL, module_list = NULL)
    observe({
        req(gene_modules())
        values$module_list <- levels(gene_modules()$module)
    })

    # gene module selector
    output$gene_module_selector <- renderUI({
        req(values$module_list)
        shinyWidgets::pickerInput(
            ns("selected_gene_module"),
            "Gene module",
            inline = TRUE,
            choices = values$module_list,
            selected = NULL,
            multiple = FALSE
        )
    })

    observe({
        selected_module(input$selected_gene_module)
    })

    # update the 'genes' reactive value when gene module changes
    observe({
        if (is.null(input$selected_gene_module) || length(input$selected_gene_module) == 0) {
            genes(character(0))
        } else {
            req(gene_modules())
            new_genes <- gene_modules() %>%
                filter(module %in% input$selected_gene_module) %>%
                pull(gene)
            genes(new_genes)
        }
    })

    observe({
        req(gene_modules)
        req(input$selected_gene_module)
        prev_genes <- gene_modules() %>%
            filter(module == input$selected_gene_module) %>%
            pull(gene)
        if (!setequal(genes(), prev_genes)) {
            shinydashboardPlus::updateBox(
                "gene_modules_box",
                action = "update",
                options = list(
                    title = h2("Gene modules (modified)"),
                    status = "warning"
                )
            )
        } else {
            shinydashboardPlus::updateBox(
                "gene_modules_box",
                action = "update",
                options = list(
                    title = h2("Gene modules"),
                    status = "primary"
                )
            )
        }
    })

    # 'Update' button
    observeEvent(input$update_module, {
        req(input$selected_gene_module)
        modules <- levels(gene_modules()$module)
        new_gene_modules <- bind_rows(
            gene_modules() %>%
                filter(module != input$selected_gene_module),
            tibble(gene = genes(), module = input$selected_gene_module)
        ) %>%
            mutate(module = factor(module, levels = modules))
        gene_modules(new_gene_modules)
    })

    # 'Remove' button
    observeEvent(input$remove_gene_module_modal, {
        req(input$selected_gene_module)
        showModal({
            modalDialog(
                title = "Remove gene module",
                glue("Are you sure you want to delete the following module: {input$selected_gene_module}?"),
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("remove_gene_module"), "OK")
                )
            )
        })
    })

    observeEvent(input$remove_gene_module, {
        req(input$selected_gene_module)
        new_gene_modules <- gene_modules() %>%
            filter(module != input$selected_gene_module) %>%
            mutate(module = forcats::fct_drop(module, only = input$selected_gene_module))
        gene_modules(new_gene_modules)
        removeModal()
    })

    # 'Add' button
    observeEvent(input$add_gene_module_modal, {
        showModal({
            modalDialog(
                title = "Add a new gene module",
                textInput(ns("new_gene_module_name"), "Gene module name"),
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("add_gene_module"), "OK")
                )
            )
        })
    })

    observeEvent(input$add_gene_module, {
        req(input$selected_gene_module)
        if (input$new_gene_module_name %in% levels(gene_modules()$module)) {
            showNotification(glue("Module {input$new_gene_module_name} already exists"), type = "error")
        } else {
            new_gene_modules <- gene_modules() %>%
                mutate(module = forcats::fct_expand(module, input$new_gene_module_name))
            gene_modules(new_gene_modules)
        }

        shinyWidgets::updatePickerInput(session, "selected_gene_module", selected = input$new_gene_module_name)

        removeModal()
    })

    # Load gene modules file
    observeEvent(input$gene_modules_fn, {
        values$file_status <- "uploaded"
    })

    observe({
        req(input$gene_modules_fn)
        req(values$file_status)

        new_gene_modules <- fread(input$gene_modules_fn$datapath, colClasses = c("gene" = "character", "module" = "character")) %>% as_tibble()

        values$file_status <- NULL

        for (field in c("gene", "module")) {
            if (!has_name(new_gene_modules, field)) {
                showNotification(glue("File should have a column named '{field}'"), type = "error")
                req(FALSE)
            }
        }

        bad_genes <- new_gene_modules %>%
            count(gene) %>%
            filter(n > 1) %>%
            pull(gene)

        if (length(bad_genes) > 0) {
            showNotification(glue("The following genes appear in more than one module: {genes}", genes = paste(bad_genes, collapse = ", ")), type = "error")
            req(FALSE)
        }

        gene_modules(new_gene_modules)
    })

    # download button
    output$gene_modules_download <- downloadHandler(
        filename = function() {
            paste("gene_modules-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            fwrite(
                gene_modules() %>%
                    select(gene, module),
                file
            )
        }
    )

    # 'Reset' button
    observeEvent(input$reset_gene_modules, {
        req(dataset())
        gene_modules(get_mc_data(dataset(), "gene_modules"))
        values$file_status <- NULL
    })

    # GENES

    # Clear genes
    observeEvent(input$clear_genes, {
        genes(character(0))
    })

    # Reset genes
    observeEvent(input$reset_genes, {
        req(gene_modules())
        req(input$selected_gene_module)
        genes(gene_modules() %>% filter(module == input$selected_gene_module) %>% pull(gene))
    })

    # remove genes
    observeEvent(input$remove_genes, {
        rows <- input$genes_table_rows_selected
        req(!is.null(rows) && length(rows) > 0)
        req(genes())
        genes_to_remove <- genes()[rows]
        new_genes <- genes()[!(genes() %in% genes_to_remove)]
        genes(new_genes)
    })

    # Add genes
    output$add_genes_sidebar <- renderUI({
        gene_choices <- gene_names(dataset())
        tagList(
            shinyWidgets::virtualSelectInput(ns("genes_to_add"),
                label = "Add genes:",
                choices = gene_choices,
                selected = c(),
                multiple = TRUE,
                showSelectedOptionsFirst = TRUE,
                search = TRUE,
                markSearchResults = TRUE,
                searchByStartsWith = TRUE,
                disableSelectAll = TRUE
            ),
            shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add selected genes", size = "sm")
        )
    })

    observeEvent(input$add_genes, {
        add_genes_to_gene_module(input$genes_to_add, gene_modules, genes)
        shinyWidgets::updateVirtualSelect(session = session, inputId = "genes_to_add", selected = character(0))
    })

    # Add from top-correlated list
    observeEvent(input$add_genes_from_top_cor_list, {
        add_genes_to_gene_module(input$selected_top_genes, gene_modules, genes)
    })

    # Add from Diff. Expr selection
    observeEvent(input$add_selected_genes, {
        add_genes_to_gene_module(selected_genes(), gene_modules, genes)
    })

    # Move genes
    observeEvent(input$move_selected_genes_modal, {
        req(input$selected_gene_module)
        showModal({
            req(values$module_list)
            modalDialog(
                title = "Move gene module",
                shinyWidgets::pickerInput(
                    ns("selected_gene_module_move"),
                    "Selected a gene module to move the gene(s) to:",
                    inline = TRUE,
                    choices = values$module_list,
                    selected = NULL,
                    multiple = FALSE
                ),
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("move_genes"), "OK")
                )
            )
        })
    })

    observeEvent(input$move_genes, {
        req(input$selected_gene_module_move)
        rows <- input$genes_table_rows_selected
        req(!is.null(rows) && length(rows) > 0)
        req(genes())

        modules <- levels(gene_modules()$module)

        # get the genes to move
        genes_to_move <- genes()[rows]
        gene_module_genes <- gene_modules() %>%
            filter(module == input$selected_gene_module_move) %>%
            pull(gene)
        genes_to_move <- genes_to_move[!(genes_to_move %in% gene_module_genes)]

        # update the other gene module
        new_gene_modules <- bind_rows(
            gene_modules(),
            tibble(gene = genes_to_move, module = input$selected_gene_module_move)
        )
        # update the current gene module
        genes(genes()[!(genes() %in% genes_to_move)])
        new_gene_modules <- bind_rows(
            new_gene_modules %>%
                filter(module != input$selected_gene_module),
            tibble(gene = genes(), module = input$selected_gene_module)
        ) %>%
            mutate(module = factor(module, levels = modules))

        # update the gene modules
        gene_modules(new_gene_modules)
        removeModal()
    })


    # Hide and show selected genes buttons
    observe({
        shinyjs::toggle(id = "remove_genes", condition = !is.null(input$genes_table_rows_selected))
        shinyjs::toggle(id = "move_selected_genes_modal", condition = !is.null(input$genes_table_rows_selected))
    })

    # Genes table
    output$genes_table <- DT::renderDataTable(
        tibble(gene = genes()),
        escape = FALSE,
        server = FALSE,
        rownames = FALSE,
        colnames = "",
        options = list(
            dom = "Bfrtip",
            paging = FALSE,
            language = list(emptyTable = "Please select a gene module.")
        )
    )
}

add_genes_to_gene_module <- function(new_genes_input, gene_modules, genes) {
    req(new_genes_input)
    to_add <- new_genes_input[!(new_genes_input %in% gene_modules()$gene)]
    if (any(new_genes_input %in% gene_modules()$gene)) {
        existing_genes <- gene_modules() %>%
            filter(gene %in% new_genes_input) %>%
            mutate(str = glue("{gene} ({module})")) %>%
            pull(str) %>%
            paste(collapse = ",")

        showNotification(glue("The following genes already exist in other gene modules: {existing_genes}"), type = "error")
    }

    if (length(to_add) > 0) {
        genes(c(to_add, genes()))
    }
}
