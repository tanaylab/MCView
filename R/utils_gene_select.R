
server_gene_selectors <- function(input, output, session, dataset, ns) {
    output$gene_selectors <- renderUI({
        picker_options <- shinyWidgets::pickerOptions(liveSearch = TRUE, liveSearchNormalize = TRUE, liveSearchStyle = "startsWith")
        div(
            id = ns("sidebar_select"),
            shinyWidgets::pickerInput(ns("gene1"), "Gene A",
                choices = gene_names(dataset()),
                selected = default_gene1,
                multiple = FALSE,
                options = picker_options
            ),
            shinyWidgets::pickerInput(ns("gene2"), "Gene B",
                choices = gene_names(dataset()),
                selected = default_gene2,
                multiple = FALSE,
                options = picker_options
            ),
            shinyWidgets::actionGroupButtons(ns("switch_genes"), labels = c("Switch"), size = "sm"),
            tags$hr()
        )
    })


    output$top_correlated_select_gene1 <- renderUI({
        req(input$gene1)
        req(has_gg_mc_top_cor(project, dataset()))
        tagList(
            selectInput(
                ns("selected_top_gene1"),
                glue("Top correlated to {input$gene1}:"),
                choices = c(get_top_cor_gene(dataset(), input$gene1, type = "pos"), rev(get_top_cor_gene(dataset(), input$gene1, type = "neg"))),
                selected = NULL,
                size = 10,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(c(ns("select_top_cor1_gene1"), ns("select_top_cor1_gene2")), labels = c("Select as Gene A", "Select as Gene B"), size = "sm")
        )
    })

    output$top_correlated_select_gene2 <- renderUI({
        req(input$gene2)
        req(has_gg_mc_top_cor(project, dataset()))
        tagList(
            selectInput(
                ns("selected_top_gene2"),
                glue("Top correlated to {input$gene2}:"),
                choices = c(get_top_cor_gene(dataset(), input$gene2, type = "pos"), rev(get_top_cor_gene(dataset(), input$gene2, type = "neg"))),
                selected = NULL,
                size = 10,
                selectize = FALSE
            ),
            shinyWidgets::actionGroupButtons(c(ns("select_top_cor2_gene1"), ns("select_top_cor2_gene2")), labels = c("Select as Gene A", "Select as Gene B"), size = "sm"),
            tags$hr()
        )
    })

    output$genecards_buttons <- renderUI({
        req(input$gene1)
        req(input$gene2)
        tagList(
            shiny::actionButton(inputId = ns("genecards1"), label = glue("GeneCards: {input$gene1}"), onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={input$gene1}')")),
            shiny::actionButton(inputId = ns("genecards2"), label = glue("GeneCards: {input$gene2}"), onclick = glue("window.open('https://www.genecards.org/cgi-bin/carddisp.pl?gene={input$gene2}')"))
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
