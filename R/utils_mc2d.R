mc2d_to_df <- function(mc2d) {
    if (is.null(mc2d)) {
        return(tibble(metacell = character(), x = numeric(), y = numeric()))
    }

    if (is.data.frame(mc2d)) {
        if (all(c("metacell", "x", "y") %in% colnames(mc2d))) {
            return(mc2d %>%
                transmute(
                    metacell = as.character(metacell),
                    x = as.numeric(x),
                    y = as.numeric(y)
                ))
        }
        if (all(c("mc_id", "mc_x", "mc_y") %in% colnames(mc2d))) {
            return(tibble(
                metacell = as.character(mc2d$mc_id),
                x = as.numeric(mc2d$mc_x),
                y = as.numeric(mc2d$mc_y)
            ))
        }
    }

    if (is.list(mc2d) && !is.null(mc2d$mc_x) && !is.null(mc2d$mc_y)) {
        metacells <- names(mc2d$mc_x)
        if (is.null(metacells) || any(metacells == "")) {
            metacells <- mc2d$mc_id %||% seq_along(mc2d$mc_x)
        }
        return(tibble(
            metacell = as.character(metacells),
            x = as.numeric(mc2d$mc_x),
            y = as.numeric(mc2d$mc_y)
        ))
    }

    if (is.list(mc2d) && !is.null(mc2d$x) && !is.null(mc2d$y)) {
        metacells <- mc2d$metacell %||% names(mc2d$x) %||% seq_along(mc2d$x)
        return(tibble(
            metacell = as.character(metacells),
            x = as.numeric(mc2d$x),
            y = as.numeric(mc2d$y)
        ))
    }

    tibble(metacell = character(), x = numeric(), y = numeric())
}

mc2d_to_graph_df <- function(mc2d, min_d = 0.05, graph = NULL) {
    mc2d_df <- mc2d_to_df(mc2d)
    empty_graph <- tibble(
        mc1 = character(),
        mc2 = character(),
        x_mc1 = numeric(),
        y_mc1 = numeric(),
        x_mc2 = numeric(),
        y_mc2 = numeric(),
        dx = numeric(),
        dy = numeric(),
        d = numeric(),
        d_norm = numeric()
    )

    if (nrow(mc2d_df) == 0) {
        return(empty_graph)
    }

    if (is.null(graph)) {
        graph <- if (is.list(mc2d)) {
            mc2d$graph
        } else {
            NULL
        }
    }

    if (is.null(graph) || (is.data.frame(graph) && nrow(graph) == 0)) {
        return(empty_graph)
    }

    max_range <- mc2d_df %>%
        filter(is.finite(x), is.finite(y)) %>%
        summarise(range_x = max(x) - min(x), range_y = max(y) - min(y)) %>%
        summarise(range = sqrt(range_x^2 + range_y^2)) %>%
        pull(range)

    if (is.null(max_range) || !is.finite(max_range) || max_range <= 0) {
        max_range <- NA_real_
    }

    if (all(c("from", "to") %in% colnames(graph))) {
        graph <- graph %>%
            mutate(mc1 = as.character(from), mc2 = as.character(to), d = weight) %>%
            left_join(mc2d_df %>% rename(mc1 = metacell, x_mc1 = x, y_mc1 = y), by = "mc1") %>%
            left_join(mc2d_df %>% rename(mc2 = metacell, x_mc2 = x, y_mc2 = y), by = "mc2") %>%
            filter(d >= min_d)
        return(graph)
    }

    graph <- graph %>%
        mutate(mc1 = as.character(mc1), mc2 = as.character(mc2)) %>%
        left_join(mc2d_df %>% rename(mc1 = metacell, x_mc1 = x, y_mc1 = y), by = "mc1") %>%
        left_join(mc2d_df %>% rename(mc2 = metacell, x_mc2 = x, y_mc2 = y), by = "mc2") %>%
        mutate(dx = x_mc1 - x_mc2, dy = y_mc1 - y_mc2, d = sqrt(dx^2 + dy^2), d_norm = (d - min(d)) / (max(d) - min(d))) %>%
        as_tibble()

    if (is.finite(max_range)) {
        graph <- graph %>%
            filter(d / max_range >= min_d)
    }

    return(graph)
}
