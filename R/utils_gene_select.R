
server_gene_selectors <- function(input, output, session, values, dataset, ns, show_button = FALSE) {
    # Gene selectors
    # When a gene changes via the picker input the query string is changed => values are changed.
    # When a gene changes via the query string the gene selectors are not rendered until values are updated, and then they are initialized with the query string genes.

    observe({
        query <- getQueryString()
        values$gene1 <- query$gene1 %||% default_gene1
        values$gene2 <- query$gene2 %||% default_gene2
    })

    # We render the gene selectors on demand in order to allow faster initialization
    if (show_button) {
        output$toggle_gene_selectors_button <- renderUI({
            div(
                actionButton(
                    ns("show_gene_selectors"),
                    "Choose genes",
                    onclick = "var $btn=$(this); setTimeout(function(){$btn.remove();},0);"
                )
            )
        })
    }

    output$gene_selectors <- renderUI({
        if (show_button) {
            req(input$show_gene_selectors)
        }

        req(values$gene1)
        req(values$gene2)
        picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
        div(
            id = ns("sidebar_select"),
            shinyWidgets::pickerInput(ns("gene1"), "Gene A",
                choices = gene_names,
                selected = values$gene1,
                multiple = FALSE,
                options = picker_options
            ),
            shinyWidgets::pickerInput(ns("gene2"), "Gene B",
                choices = gene_names,
                selected = values$gene2,
                multiple = FALSE,
                options = picker_options
            ),
            shinyWidgets::actionGroupButtons(ns("switch_genes"), labels = c("Switch"), size = "sm"),
            tags$hr()
        )
    })


    gene_changed <- reactive({
        list(input$gene1, input$gene2)
    })

    observeEvent(gene_changed(), {
        query_string <- glue("?gene1={input$gene1}&gene2={input$gene2}")
        updateQueryString(query_string, mode = "push")
    })

    output$top_correlated_select_gene1 <- renderUI({
        if (show_button) {
            req(input$show_gene_selectors)
        }
        req(values$gene1)
        req(input$gene1)
        req(has_gg_mc_top_cor(project, dataset()))
        tagList(
            selectInput(
                ns("selected_top_gene1"),
                glue("Top correlated to {values$gene1}:"),
                choices = c(get_top_cor_gene(dataset(), values$gene1, type = "pos"), rev(get_top_cor_gene(dataset(), values$gene1, type = "neg"))),
                selected = NULL,
                size = 10,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(c(ns("select_top_cor1_gene1"), ns("select_top_cor1_gene2")), labels = c("Select as Gene A", "Select as Gene B"), size = "sm")
        )
    })

    output$top_correlated_select_gene2 <- renderUI({
        if (show_button) {
            req(input$show_gene_selectors)
        }
        req(values$gene2)
        req(input$gene2)
        req(has_gg_mc_top_cor(project, dataset()))
        tagList(
            selectInput(
                ns("selected_top_gene2"),
                glue("Top correlated to {values$gene2}:"),
                choices = c(get_top_cor_gene(dataset(), values$gene2, type = "pos"), rev(get_top_cor_gene(dataset(), values$gene2, type = "neg"))),
                selected = NULL,
                size = 10,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(c(ns("select_top_cor2_gene1"), ns("select_top_cor2_gene2")), labels = c("Select as Gene A", "Select as Gene B"), size = "sm"),
            tags$hr()
        )
    })

    output$genecards_buttons <- renderUI({
        if (show_button) {
            req(input$show_gene_selectors)
        }
        req(values$gene1)
        req(values$gene2)
        tagList(
            shiny::actionButton(inputId = ns("genecards1"), label = glue("GeneCards: {values$gene1}"), onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={values$gene1}')")),
            shiny::actionButton(inputId = ns("genecards2"), label = glue("GeneCards: {values$gene2}"), onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={values$gene2}')"))
        )
    })

    observeEvent(input$select_top_cor1_gene1, {
        updateSelectInput(session, "gene1", selected = input$selected_top_gene1)
    })

    observeEvent(input$select_top_cor1_gene2, {
        updateSelectInput(session, "gene2", selected = input$selected_top_gene1)
    })

    observeEvent(input$select_top_cor2_gene1, {
        updateSelectInput(session, "gene1", selected = input$selected_top_gene2)
    })

    observeEvent(input$select_top_cor2_gene2, {
        updateSelectInput(session, "gene2", selected = input$selected_top_gene2)
    })

    observeEvent(input$switch_genes, {
        gene1 <- input$gene1
        gene2 <- input$gene2
        updateSelectInput(session, "gene1", selected = gene2)
        updateSelectInput(session, "gene2", selected = gene1)
    })
}
