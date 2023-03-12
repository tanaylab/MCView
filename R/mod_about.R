#' about UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_about_ui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            generic_box(
                id = ns("about"),
                title = "About",
                collapsible = FALSE,
                closable = FALSE,
                width = 12,
                height = "80vh",
                includeRMarkdown(about_file),
                plotOutput(ns("about_2d_plot"))
            )
        )
    )
}

#' about sidebar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_about_sidebar_ui <- function(id) {
    ns <- NS(id)
    tagList(
        list()
    )
}

#' about Server Function
#'
#' @noRd
mod_about_server <- function(id, dataset, metacell_types, cell_type_colors, gene_modules, globals) {
    moduleServer(
        id,
        function(input, output, session) {
            ns <- session$ns
            output$about_2d_plot <- renderPlot(
                {
                    mc2d <- get_mc_data(dataset(), "mc2d")
                    ct <- tibble::deframe(cell_type_colors() %>% select(cell_type, color))
                    graph <- mc2d_to_graph_df(mc2d, min_d = min_edge_length(dataset()))
                    p_proj <- mc2d_to_df(mc2d) %>%
                        left_join(metacell_types(), by = "metacell") %>%
                        ggplot(aes(x = x, y = y, fill = cell_type)) +
                        geom_segment(data = graph, inherit.aes = FALSE, aes(x = x_mc1, y = y_mc1, xend = x_mc2, yend = y_mc2), color = "black", size = 0.1) +
                        geom_point(size = initial_proj_point_size(dataset()), shape = 21, color = "black", stroke = 0.1) +
                        scale_fill_manual(values = ct) +
                        theme_void() +
                        theme(legend.position = "none") +
                        theme(aspect.ratio = 1) +
                        guides(color = "none")
                    p_gg <- plot_gg_over_mc(dataset(), default_gene1, default_gene2, metacell_types = metacell_types(), cell_type_colors = cell_type_colors(), plot_text = FALSE) +
                        theme(legend.position = "none") +
                        theme(aspect.ratio = 1) +
                        guides(color = "none")
                    cowplot::plot_grid(p_proj, p_gg, ncol = 2)
                },
                res = 96
            ) %>% bindCache(dataset(), metacell_types(), cell_type_colors())
        }
    )
}
