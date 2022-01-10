metacell_names_reactive <- function(dataset) {
    reactive({
        req(dataset())
        colnames(get_mc_data(dataset(), "mc_mat"))
    })
}

metacell_colors_reactive <- function(dataset, metacell_names, metacell_types) {
    reactive({
        req(metacell_names())
        req(metacell_types())
        tibble(metacell = metacell_names()) %>%
            left_join(metacell_types() %>% select(metacell, mc_col), by = "metacell") %>%
            pull(mc_col)
    })
}
