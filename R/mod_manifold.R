#' manifold UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_manifold_ui <- function(id) {
    ns <- NS(id)
    tagList(
        generic_column(
            width = 12,
            projection_box(ns, "gene_projection", title = "2D Projection", collapsed_accordion = FALSE, show_legend = TRUE, color_choices = c("Cell type", "Gene", "Gene module", "Metadata"), plotly_height = "70vh", height = "70vh")
        )
    )
}


#' manifold sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_manifold_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list(
            uiOutput(ns("top_correlated_select_color_proj")),
            shinyWidgets::actionGroupButtons(ns("recompute"), labels = "Recompute 2D projection", size = "sm"),
            shinyWidgets::actionGroupButtons(ns("reset"), labels = "Restore default", size = "sm"),
            tags$hr(),
            uiOutput(ns("add_genes_ui")),
            selectInput(
                ns("selected_anchor_genes"),
                "Anchor Genes",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                size = 15,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(ns("remove_genes"), labels = "Remove selected genes", size = "sm"),
            uiOutput(ns("add_gene_modules_ui")),
            tags$hr(),
            downloadButton(ns("download_genes"), "Save genes", align = "center", style = "margin: 5px 5px 5px 15px; "),
            fileInput(ns("load_genes"),
                label = NULL,
                buttonLabel = "Load genes",
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
            tags$hr(),
            numericInput(ns("genes_per_anchor"), "Genes per anchor", value = 30, min = 1, max = 100, step = 1),
            numericInput(ns("n_neighbors"), "Number of neighbors", value = 10, min = 1, max = 100, step = 1),
            numericInput(ns("min_dist"), "Minimum distance", value = 0.96, min = 0, max = 1, step = 0.01),
            numericInput(ns("n_epoch"), "Number of epochs", value = 500, min = 1, max = 10000, step = 1),
            numericInput(ns("min_log_expr"), "Minimum log expression", value = -14, min = -50, max = 0, step = 0.1),
            tags$hr(),
            downloadButton(ns("download_projection"), "Download 2D layout", align = "center", style = "margin: 5px 5px 5px 15px; "),
            fileInput(ns("load_projection"),
                label = NULL,
                buttonLabel = "Load 2D layout",
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
            tags$hr(),
            downloadButton(ns("download_graph"), "Download graph", align = "center", style = "margin: 5px 5px 5px 15px; "),
            fileInput(ns("load_graph"),
                label = NULL,
                buttonLabel = "Load graph",
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
    )
}


#' manifold Server Function
#'
#' @noRd
mod_manifold_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns

            top_correlated_selectors(input, output, session, dataset, ns, gene_modules = gene_modules)
            projection_selectors(ns, dataset, output, input, gene_modules, globals, session)

            clipboard_changed <- clipboard_changed_2d_reactive(input, globals)

            # 2D projection re-calculations
            observeEvent(input$recompute, {
                req(globals$anchor_genes)
                showNotification("Recomputing 2D projection")
                mc_egc <- get_mc_egc(dataset())

                mc2d <- compute_umap(mc_egc, globals$anchor_genes, n_neighbors = input$n_neighbors, min_dist = input$min_dist, n_epoch = input$n_epoch, min_log_expr = input$min_log_expr, genes_per_anchor = input$genes_per_anchor)
                if (is.null(mc2d)) {
                    showNotification("Recomputing 2D projection failed", type = "error")
                }
                req(mc2d)
                globals$mc2d <- mc2d
                shinyjs::show("reset")
            })

            observeEvent(input$reset, {
                globals$mc2d <- get_mc_data(dataset(), "mc2d")
                shinyjs::hide("reset")
            })

            observe({
                shinyjs::hide("reset")
                updateSelectInput(session, "selected_anchor_genes", choices = globals$anchor_genes, label = glue("Anchor Genes ({length(globals$anchor_genes)})"))
            })

            output$add_genes_ui <- renderUI({
                gene_choices <- gene_names(dataset())
                tagList(
                    shinyWidgets::virtualSelectInput(ns("genes_to_add"),
                        label = "Add genes",
                        choices = gene_choices,
                        selected = c(),
                        multiple = TRUE,
                        showSelectedOptionsFirst = TRUE,
                        search = TRUE,
                        markSearchResults = TRUE,
                        searchByStartsWith = TRUE
                    ),
                    shinyWidgets::actionGroupButtons(ns("add_genes"), labels = "Add genes", size = "sm")
                )
            })

            observeEvent(input$add_genes, {
                new_anchors <- sort(unique(c(globals$anchor_genes, input$genes_to_add)))
                names(new_anchors) <- add_gene_modules(new_anchors, dataset(), gene_modules())
                globals$anchor_genes <- new_anchors
                shinyWidgets::updateVirtualSelect(session = session, inputId = "genes_to_add", selected = character(0))
            })

            observeEvent(input$remove_genes, {
                new_anchors <- setdiff(globals$anchor_genes, input$selected_anchor_genes)
                globals$anchor_genes <- new_anchors
            })

            output$add_gene_modules_ui <- renderUI({
                modules <- unique(gene_modules()$module)
                tagList(
                    shinyWidgets::virtualSelectInput(ns("gene_modules_to_add"),
                        label = "Add gene modules genes",
                        choices = modules,
                        selected = c(),
                        multiple = TRUE,
                        showSelectedOptionsFirst = TRUE,
                        search = TRUE,
                        markSearchResults = TRUE
                    ),
                    shinyWidgets::actionGroupButtons(c(ns("add_gene_modules"), ns("remove_gene_modules")), labels = c("Add modules genes", "Remove modules genes"), size = "sm")
                )
            })

            observeEvent(input$add_gene_modules, {
                req(input$gene_modules_to_add)
                new_anchors <- sort(unique(c(globals$anchor_genes, gene_modules()$gene[gene_modules()$module %in% input$gene_modules_to_add])))
                names(new_anchors) <- add_gene_modules(new_anchors, dataset(), gene_modules())
                globals$anchor_genes <- new_anchors
                shinyWidgets::updateVirtualSelect(session = session, inputId = "gene_modules_to_add", selected = character(0))
            })

            observeEvent(input$remove_gene_modules, {
                req(input$gene_modules_to_add)
                new_anchors <- setdiff(globals$anchor_genes, gene_modules()$gene[gene_modules()$module %in% input$gene_modules_to_add])
                globals$anchor_genes <- new_anchors
                shinyWidgets::updateVirtualSelect(session = session, inputId = "gene_modules_to_add", selected = character(0))
            })

            output$download_genes <- downloadHandler(
                filename = function() {
                    paste0("anchor_genes_", dataset(), "_", Sys.Date(), ".csv")
                },
                content = function(file) {
                    fwrite(tibble(gene = globals$anchor_genes), file, row.names = FALSE, col.names = FALSE)
                }
            )

            observeEvent(input$load_genes, {
                req(input$load_genes)
                req(input$load_genes$datapath)
                req(file.exists(input$load_genes$datapath))
                new_anchors <- fread(input$load_genes$datapath, header = FALSE)[, 1]

                # check if all genes are present in the dataset
                unknown_genes <- setdiff(new_anchors, gene_names(dataset()))
                if (length(unknown_genes) > 0) {
                    showNotification(paste0("Unknown genes: ", paste(unknown_genes, collapse = ", ")), type = "error")
                    new_anchors <- setdiff(new_anchors, unknown_genes)
                }

                req(length(new_anchors) > 0)
                globals$anchor_genes <- new_anchors
            })

            output$download_projection <- downloadHandler(
                filename = function() {
                    paste0("2d_layout_", dataset(), "_", Sys.Date(), ".csv")
                },
                content = function(file) {
                    mc2d <- globals$mc2d
                    req(mc2d)
                    fwrite(mc2d_to_df(mc2d), file, row.names = FALSE)
                }
            )

            observeEvent(input$load_projection, {
                req(input$load_projection)
                req(input$load_projection$datapath)
                mc2d <- layout_and_graph_to_mc2d(input$load_projection$datapath, globals$mc2d$graph %>% rename(from = mc1, to = mc2), metacells = get_metacell_ids(project, dataset()), warn_function = function(msg) {
                    showNotification(msg, type = "error")
                }, error_function = function(msg) {
                    showNotification(msg, type = "error")
                })
                req(mc2d)
                globals$mc2d <- mc2d
            })


            output$download_graph <- downloadHandler(
                filename = function() {
                    paste0("graph_", dataset(), "_", Sys.Date(), ".csv")
                },
                content = function(file) {
                    mc2d <- globals$mc2d
                    req(mc2d)
                    fwrite(mc2d$graph %>% rename(from = mc1, to = mc2), file, row.names = FALSE)
                }
            )

            observeEvent(input$load_graph, {
                req(input$load_graph)
                req(input$load_graph$datapath)
                mc2d <- layout_and_graph_to_mc2d(mc2d_to_df(globals$mc2d), input$load_graph$datapath, metacells = get_metacell_ids(project, dataset()), warn_function = function(msg) {
                    showNotification(msg, type = "error")
                }, error_function = function(msg) {
                    showNotification(msg, type = "error")
                })
                req(mc2d)
                globals$mc2d <- mc2d
            })

            # Projection plots
            output$plot_gene_proj_2d <- render_2d_plotly(input, output, session, dataset, metacell_types, cell_type_colors, gene_modules, globals, source = "proj_manifold_plot") %>%
                bindCache(
                    dataset(),
                    input$color_proj,
                    metacell_types(),
                    cell_type_colors(),
                    input$point_size,
                    input$stroke,
                    input$min_edge_size,
                    input$set_range,
                    input$metacell1,
                    input$metacell2,
                    input$proj_stat,
                    input$expr_range,
                    input$lfp,
                    input$color_proj_gene,
                    input$color_proj_metadata,
                    input$color_proj_gene_module,
                    clipboard_changed(),
                    input$graph_name,
                    input$legend_orientation,
                    input$show_legend_projection,
                    globals$mc2d
                )
        }
    )
}
