mc2d_to_df <- function(mc2d) {
    tibble(metacell = names(mc2d$mc_x), x = mc2d$mc_x, y = mc2d$mc_y)
}

mc2d_to_graph_df <- function(mc2d, min_d = 0.05) {
    mc2d_df <- mc2d_to_df(mc2d)

    max_range <- mc2d_df %>%
        filter(is.finite(x), is.finite(y)) %>%
        summarise(range_x = max(x) - min(x), range_y = max(y) - min(y)) %>%
        summarise(range = sqrt(range_x^2 + range_y^2)) %>%
        pull(range)

    graph <- mc2d$graph %>%
        mutate(mc1 = as.character(mc1), mc2 = as.character(mc2)) %>%
        left_join(mc2d_df %>% rename(mc1 = metacell, x_mc1 = x, y_mc1 = y), by = "mc1") %>%
        left_join(mc2d_df %>% rename(mc2 = metacell, x_mc2 = x, y_mc2 = y), by = "mc2") %>%
        mutate(dx = x_mc1 - x_mc2, dy = y_mc1 - y_mc2, d = sqrt(dx^2 + dy^2), d_norm = (d - min(d)) / (max(d) - min(d))) %>%
        as_tibble()

    # we show only edges that are at least min_d% of the maximum possible width
    graph <- graph %>%
        filter(d / max_range >= min_d)

    return(graph)
}
