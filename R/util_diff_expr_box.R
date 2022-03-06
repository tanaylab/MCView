diff_expr_box <- function(ns,
                          id,
                          title = "Differential expression",
                          choices = c("MCs", "Types"),
                          selected = "Types",
                          ...) {
    shinydashboardPlus::box(
        id = ns(id),
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        closable = FALSE,
        width = 12,
        sidebar = shinydashboardPlus::boxSidebar(
            startOpen = FALSE,
            width = 80,
            id = ns("mc_mc_sidebar"),
            shinyWidgets::radioGroupButtons(
                inputId = ns("mode"),
                label = "Compare:",
                choices = choices,
                selected = selected,
                justified = TRUE
            ),
            uiOutput(ns("metacell1_select")),
            uiOutput(ns("metacell2_select")),
            shinyWidgets::actionGroupButtons(ns("switch_metacells"), labels = c("Switch"), size = "sm")
        ),
        ...,
        shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("plot_mc_mc_gene_scatter"))
        ),
        shinyWidgets::prettySwitch(inputId = ns("show_diff_expr_table"), value = FALSE, label = "Show table"),
        DT::DTOutput(ns("diff_expr_table"))
    )
}

diff_expr_outputs <- function(input, output, session, dataset, metacell_types, cell_type_colors, globals, ns, source_suffix, dragmode = NULL, plotly_buttons = c("select2d", "lasso2d", "hoverClosestCartesian", "hoverCompareCartesian", "toggleSpikelines")) {
    metacell_names <- metacell_names_reactive(dataset)
    metacell_colors <- metacell_colors_reactive(dataset, metacell_names, metacell_types)

    metacell_selectors(input, output, session, dataset, ns, metacell_names, metacell_colors, metacell_types, cell_type_colors)
    group_selectors(input, output, session, dataset, ns)

    mc_mc_gene_scatter_df <- mc_mc_gene_scatter_df_reactive(dataset, input, output, session, metacell_types, cell_type_colors)

    diff_expr_switch_metacells(dataset, input, output, session)

    output$plot_mc_mc_gene_scatter <- render_mc_mc_gene_plotly(input, output, session, ns, dataset, mc_mc_gene_scatter_df, metacell_names(), cell_type_colors(), source_suffix = source_suffix, dragmode = dragmode, plotly_buttons = plotly_buttons)

    output$diff_expr_table <- render_mc_mc_gene_diff_table(input, output, session, ns, dataset, mc_mc_gene_scatter_df)
}
