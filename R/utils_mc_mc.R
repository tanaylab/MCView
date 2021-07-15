metacell_from_coords_proj <- function(dataset, x, y) {
    mc2d <- get_mc_data(dataset, "mc2d")
    metacell_types <- get_mc_data(dataset, "metacell_types")

    df <- mc2d_to_df(mc2d) %>% left_join(metacell_types, by = "metacell")

    metacell <- df %>%
        mutate(diff1 = x - !!x, diff2 = y - !!y) %>%
        arrange(abs(diff1), abs(diff2)) %>%
        slice(1) %>%
        pull(metacell)

    return(metacell)
}


get_top_var_genes <- function(dataset, metacell) {
    metacell <- as.character(metacell)
    mc_egc_t <- calc_mc_egc_t(dataset, metacell)
    vec <- log2(mc_egc_t + 3e-4)
    diff <- matrixStats::rowMaxs(vec, na.rm = TRUE) - matrixStats::rowMins(vec, na.rm = TRUE)
    names(diff) <- rownames(mc_egc_t)
    top_var_genes <- names(sort(diff, decreasing = TRUE)[1:4])
    return(top_var_genes)
}
