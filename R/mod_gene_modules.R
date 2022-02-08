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
            resizable_column(
                width = 9,
                style = "padding-right:0px;",
                heatmap_box(ns("gene_modules_heatmap"), "Gene modules Heatmap", legend_width = 2)
            ),
            resizable_column(
                width = 3, 
                style = "padding-right:0px; padding-left:0px;",
                shinydashboardPlus::box(
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
                        ),
                        downloadButton(ns("gene_modules_download"), "Export", style = "align-items: left;")
                        
                    ),                        
                    actionButton(ns("remove_gene_module_modal"), "Remove", style = "align-items: left;"),
                    actionButton(ns("merge_gene_module_modal"), "Merge", style = "align-items: left;"),
                    actionButton(ns("reset_gene_modules"), "Reset", style = "align-items: left;"),                
                    br(),
                    shinycssloaders::withSpinner(
                        DT::dataTableOutput(ns("gene_modules_table"))
                    )                    
                )

            )
        ),
        fluidRow(
            resizable_column(
                width = 9,
                style = "padding-right:0px;",
                heatmap_box(ns("genes_heatmap"), "Genes Heatmap", legend_width = 2)
            ),
            resizable_column(
                width = 3,
                style = "padding-right:0px; padding-left:0px;",
                shinydashboardPlus::box(
                    id = ns("genes_box"),
                    title = "Gene modules",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    closable = FALSE,
                    width = 12,                                        
                    actionButton(ns("show_module"), "Show", style = "align-items: center;"),
                    actionButton(ns("update_module"), "Update", style = "align-items: center;"),                    
                    actionButton(ns("add_genes_modal"), "Add", style = "align-items: center;"),
                    actionButton(ns("clear_genes"), "Clear", style = "align-items: center;"),
                    br(),                    
                    br(),
                    shinyWidgets::prettyRadioButtons(
                            ns("gene_module_type"),  
                            label = NULL,                           
                            choices = c("Existing module", "New module"),
                            inline = TRUE,                            
                            fill = TRUE
                    ),
                    uiOutput(ns("gene_module_selector")),                                        
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
            HTML("<h5><b><center>Genes modules heatmap</center></b></h5>"),
            uiOutput(ns_heatmap("reset_zoom_ui")),
            uiOutput(ns("shown_gene_modules_ui")),
            uiOutput(ns_heatmap("cell_type_list")),
            uiOutput(ns_heatmap("metadata_list")),
            tags$hr(),
            HTML("<h5><b><center>Genes heatmap</center></b></h5>"),
            uiOutput(ns_genes("reset_zoom_ui")),
            uiOutput(ns("selected_gene_modules_ui")),
            uiOutput(ns_genes("cell_type_list")),
            uiOutput(ns_genes("metadata_list"))
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
            genes <- reactiveVal()
            lfp_range <- reactiveVal()

            output$shown_gene_modules_ui <- gene_modules_selector(
                dataset,
                gene_modules,
                ns,
                label = "Shown gene modules",
                id = "shown_gene_modules",
                selected = "all"
            )

            observe({
                shown_gene_modules(input$shown_gene_modules)
            })

            heatmap_reactives("gene_modules_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, shown_gene_modules, lfp_range, "Gene modules")

            gene_modules_table_reactives(dataset, input, output, session, gene_modules, globals)

            # Genes of a specific gene module
            output$selected_gene_modules_ui <- renderUI({
                modules <- unique(gene_modules()$module)

                shinyWidgets::pickerInput(
                    ns("selected_gene_modules"),
                    "Selected gene modules",
                    choices = modules,
                    selected = NULL,
                    multiple = TRUE,
                    options = shinyWidgets::pickerOptions(
                        liveSearch = TRUE,
                        liveSearchNormalize = TRUE,
                        liveSearchStyle = "startsWith",
                        `dropup-auto` = FALSE,
                        `max-options` = 3,
                        `max-options-text` = "Cannot choose more than 3 gene modules"
                    )
                )
            })

            observe({
                if (is.null(input$selected_gene_modules) || length(input$selected_gene_modules) == 0) {
                    genes(character(0))
                } else {
                    req(gene_modules())
                    new_genes <- gene_modules() %>%
                        filter(module %in% input$selected_gene_modules) %>%
                        pull(gene)
                    genes(new_genes)
                }
            })

            heatmap_reactives("genes_heatmap", dataset, metacell_types, gene_modules, cell_type_colors, globals, genes, lfp_range, "Markers")

            gene_module_reactives(dataset, input, output, session, gene_modules, globals)
        }
    )
}

gene_modules_table_reactives <- function(dataset, input, output, session, gene_modules, globals){
    ns <- session$ns

    values <- reactiveValues(file_status = NULL)

    gene_modules_tab <- reactive({
        gene_modules() %>% distinct(module)
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

    observeEvent(input$reset_gene_modules, {
        req(dataset())
        gene_modules(get_mc_data(dataset(), "gene_modules"))
        values$file_status <- NULL
    })

    observeEvent(input$remove_gene_module_modal, {
        rows <- input$gene_modules_table_rows_selected
        req(!is.null(rows) && length(rows) > 0)
        modules_to_remove <- gene_modules_tab() %>%
            slice(rows) %>%
            pull(module) %>%
            paste(collapse = ", ")
        
        showModal({            
            modalDialog(
                title = "Remove gene module(s)",
                glue("Are you sure you want to delete the following modules: {modules_to_remove}?"),
                footer = tagList(
                    modalButton("Cancel"),
                    actionButton(ns("remove_gene_module"), "OK")
                )
            )
        })
    })

    observeEvent(input$remove_gene_module, {
        rows <- input$gene_modules_table_rows_selected
        req(!is.null(rows) && length(rows) > 0)
        modules_to_remove <- gene_modules_tab() %>%
            slice(rows) %>%
            pull(module)
        gene_modules(gene_modules() %>% filter(!(module %in% modules_to_remove)))   
        removeModal()
    })

    observeEvent(input$merge_gene_module_modal, {
        rows <- input$gene_modules_table_rows_selected
        req(!is.null(rows) && length(rows) >= 2)     
        showModal(gene_module_name_modal(ns, dataset, gene_modules, "merge_gene_module", "Merge gene module"))
    })

    observeEvent(input$merge_gene_module, {
        req(input$new_gene_module_name)
        rows <- input$gene_modules_table_rows_selected
        req(!is.null(rows) && length(rows) >= 2)
        modules_to_merge <- gene_modules_tab() %>%
            slice(rows) %>%
            pull(module)
        if (input$new_gene_module_name %in% setdiff(gene_modules_tab()$module, modules_to_merge)){
            showModal(gene_module_name_modal(ns, dataset, gene_modules, "merge_gene_module", "Merge gene module", name_exists = TRUE))
        } else {        
            new_module <- gene_modules() %>%
                filter(module %in% modules_to_merge) %>%
                mutate(module = input$new_gene_module_name)
            
            new_modules <- bind_rows(
                gene_modules() %>% filter(!(module %in% modules_to_merge)),
                new_module
            )                
            
            gene_modules(new_modules)
            removeModal()
        }
    })

    gene_modules_table_proxy <- DT::dataTableProxy("gene_modules_table")

    observeEvent(input$gene_modules_table_cell_edit, {
        new_value <- input$gene_modules_table_cell_edit$value

        if (new_value %in% gene_modules_tab()$module){
            new_value <- paste0(new_value, "_1")
            showNotification(glue("Module with the name {input$gene_modules_table_cell_edit$value} already exists. Changing name to {new_value}"), type = "warning")            
        }

        old_value <- gene_modules_tab()$module[input$gene_modules_table_cell_edit$row]        

        new_gene_modules <- gene_modules() %>%
            mutate(module = ifelse(module == old_value, new_value, module))

        gene_modules(new_gene_modules)        
    })

    output$gene_modules_table <- DT::renderDataTable(
        DT::datatable(
            gene_modules_tab(),
            escape = FALSE,            
            rownames = FALSE,
            editable = "cell",
            colnames = "",        
            options = list(
                dom = "Bfrtip",
                paging = TRUE,
                language = list(emptyTable = "No gene modules. Load or create to get started"),
                pageLength = 14
            )
        ),
        server = FALSE        
    )

    # add_gene_module
}

gene_module_reactives <- function(dataset, input, output, session, gene_modules, globals){
    ns <- session$ns
    cur_module <- reactiveVal()

    observe({
        cur_module <- tibble(gene = character(0))
    })

    output$gene_module_selector <- renderUI({
        req(gene_modules())
        shinyWidgets::pickerInput(
            ns("existing_gene_module"),
            "Gene module",
            inline = TRUE,
            choices = unique(gene_modules()$module),
            selected = NULL,
            multiple = FALSE
        )
    })

    observe({
        req(input$gene_module_type)
        shinyjs::toggle(id = "gene_module_selector", condition = input$gene_module_type == "Existing module")
        shinyjs::toggle(id = "show_module", condition = input$gene_module_type == "Existing module")        
    })

    observe({        
        shinyjs::toggle(id = "remove_genes", condition = !is.null(input$genes_table_rows_selected))
        shinyjs::toggle(id = "move_selected_genes_modal", condition = !is.null(input$genes_table_rows_selected))
    })

    observe({
        req(input$existing_gene_module)
        req(input$gene_module_type == "Existing module")
        cur_module(
            gene_modules() %>%
                filter(module == input$existing_gene_module) %>% 
                select(gene)
        )
    })

    observeEvent(input$clear_genes, {
        cur_module(tibble(gene = character(0)))
    })

    observeEvent(input$show_module, {
        req(input$existing_gene_module)        
        shinyWidgets::updatePickerInput(session, "selected_gene_modules", selected = input$existing_gene_module)
    })

    observeEvent(input$update_module, {
        req(input$gene_module_type)
        if (input$gene_module_type == "New module"){
            # modal etc. 
            browser()
        } else {
            req(input$existing_gene_module)
            new_gene_modules <- gene_modules() %>% filter(module != input$existing_gene_module)
            gene_modules(
                bind_rows(
                    gene_modules() %>% filter(module != input$existing_gene_module),
                    cur_module() %>%
                        mutate(module = input$existing_gene_module) %>%
                        select(gene, module)
                )
            )
        }
    })

    observeEvent(input$remove_genes, {
        rows <- input$genes_table_rows_selected
        req(!is.null(rows) && length(rows) > 0)
        genes_to_remove <- cur_module() %>%
            slice(rows) %>%
            pull(gene)

        cur_module(
            cur_module() %>% filter(!(gene %in% genes_to_remove))
        )
    })

    # add_genes_modal
    # add_gene_module_modal        
    
    # add_selected_genes
    # move_selected_genes_modal


    output$genes_table <- DT::renderDataTable(
        cur_module(),
        escape = FALSE,
        server = FALSE,
        rownames = FALSE,
        colnames = "",
        options = list(
            dom = "Bfrtip",
            paging = TRUE,
            language = list(emptyTable = "Please select a gene module."),
            pageLength = 14
        )
    )



}

gene_module_name_modal <- function(ns, dataset, gene_modules, ok_id, title, name_exists = FALSE) {
    modalDialog(
        title = title,
        textInput(ns("new_gene_module_name"), "Gene module name"),
        if (name_exists) {
            span(tags$b("Gene module already exists", style = "color: red;"))
        },
        footer = tagList(
            modalButton("Cancel"),
            actionButton(ns(ok_id), "OK")
        )
    )
}