plot_mc_time_dist <- function(dataset, metacell, ylab = "# of cells") {
    metacell_type<- get_mc_data(dataset, "metacell_types")
    mc_ag <- get_mc_data(dataset, "mc_ag")
    time_annot <- get_mc_data(dataset, "time_annot")

    metacell <- as.character(metacell)

    cell_type <- metacell_type%>%
        filter(metacell == !!metacell) %>%
        pull(cell_type)

    df <- mc_ag[metacell, ] %>%
        tibble::enframe("time_bin", "n_cells") %>%
        mutate(time_bin = as.numeric(time_bin))
    if (!is.null(time_annot)) {
        df <- df %>%
            left_join(time_annot, by = "time_bin") %>%
            rename(time = time_desc) %>%
            mutate(time = as.numeric(time))
        xlim <- c(min(time_annot$time_desc), max(time_annot$time_desc))
    } else {
        df <- df %>% rename(time = time_bin)
        xlim <- c(min(df$time, na.rm = TRUE), max(df$time, na.rm = TRUE))
    }

    df <- df %>%
        mutate(
            `Metacell age (E[t])` = paste(
                glue("{round(time, digits=2)}"),
                glue("MC #{metacell} ({cell_type})"),
                glue("# of cells: {n_cells}"),
                sep = "\n"
            )
        )

    p <- df %>%
        ggplot(aes(x = time, y = n_cells, tooltip_text = `Metacell age (E[t])`)) +
        # geom_col(fill = metacell_type%>% filter(metacell == !!metacell) %>% pull(mc_col)) +
        geom_col(fill = "black") +
        xlab("Metacell age (E[t])") +
        ylab(ylab) +
        coord_cartesian(ylim = c(0, max(mc_ag)), xlim = xlim)

    return(p)
}
