get_tab_steps <- function(tab, module_id = NULL) {
    steps <- tibble(element = names(help_config[[tab]]), intro = help_config[[tab]])

    if (!is.null(module_id)) {
        steps <- steps %>% mutate(element = NS(module_id)(element))
    }

    steps <- steps %>%
        mutate(element = ifelse(element == "NA", NA, paste0("#", element)))

    steps <- as.data.frame(steps)

    return(steps)
}


help_reactives <- function(input, output, session, globals) {
    show_help <- function(input, output, session) {
        if (input$tab_sidebar == "about") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = FALSE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("about")
                )
            )
        } else if (input$tab_sidebar == "manifold") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("manifold", "manifold")
                )
            )
        } else if (input$tab_sidebar == "gene_mc") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("genes", "gene_mc")
                )
            )
        } else if (input$tab_sidebar == "mc_mc") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("metacells", "mc_mc")
                )
            )
        } else if (input$tab_sidebar == "cell_type") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("cell_type", "cell_type")
                )
            )
        } else if (input$tab_sidebar == "samples") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("samples", "samples")
                )
            )
        } else if (input$tab_sidebar == "Query") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("query", "query")
                )
            )
        } else if (input$tab_sidebar == "Atlas") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("atlas", "atlas")
                )
            )
        } else if (input$tab_sidebar == "annotate") {
            rintrojs::introjs(session,
                options = list(
                    "showProgress" = TRUE,
                    "showBullets" = FALSE,
                    "nextLabel" = "next",
                    "prevLabel" = "back",
                    steps = get_tab_steps("metacells", "annotate")
                )
            )
        }
    }

    observe({
        req(config$help)
        rintrojs::introjs(session,
            options = list(
                "showProgress" = FALSE,
                "showBullets" = FALSE,
                "nextLabel" = "next",
                "prevLabel" = "back",
                "scrollToElement" = TRUE,
                steps =
                    rbind(
                        data.frame(element = NA, intro = help_config$introduction),
                        get_tab_steps("about")
                    )
            )
        )
    })

    observeEvent(
        input$help,
        show_help(input, output, session)
    )
}
