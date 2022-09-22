read_metacell_graphs <- function(graph_list, metacells) {
    graph_names <- names(graph_list)
    if (is.null(graph_names)) {
        graph_names <- paste0("graph", 1:length(graph_list))
    }
    if ("metacell" %in% graph_names) {
        cli_abort("The graph name 'metacell' is reserved.")
    }
    cli_alret_info("graph names: {.val {graph_names}}")
    cl <- map_chr(graph_list, class)
    if (all(cl == "character")) {
        cli_alert_info("Reading graph files: {.val {paste0(graph_list, collapse = ', ')}}")
        graph_list <- map(graph_list, fread)
    }

    graphs <- map(graph_list, function(g) {
        if (!all(c("from", "to") %in% names(g))) {
            cli_abort("Graphs must have columns 'from' and 'to'")
        }

        if (!has_name(g, "weight")) {
            g$weight <- 1
        }

        if (!all(g$from %in% metacells)) {
            cli_abort("Some metacells in the 'from' column are not in the metacell list")
        }

        if (!all(g$to %in% metacells)) {
            cli_abort("Some metacells in the 'to' column are not in the metacell list")
        }

        g <- g %>%
            filter(
                from != to,
                weight != 0
            ) %>%
            mutate(
                from = as.character(from),
                to = as.character(to)
            )

        return(g)
    })

    names(graphs) <- graph_names

    return(graphs)
}
