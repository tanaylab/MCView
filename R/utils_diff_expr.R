calc_diff_expr <- function(mat, egc, columns, diff_thresh = 1.5, pval_thresh = 0.01) {
    df <- egc %>%
        as.data.frame()

    df$diff <- log2(df[, 1]) - log2(df[, 2])

    df$pval <- NA

    f <- rownames(df)[abs(df$diff) >= diff_thresh]

    m <- mat[f, columns]

    if (length(f) == 1) {
        m <- t(as.matrix(m))
        rownames(m) <- f
    }

    if (nrow(m) > 0) {
        tots <- colSums(m)

        pvals <- apply(m, 1, function(x) suppressWarnings(chisq.test(matrix(c(x, tots), nrow = 2))$p.value))

        df[f, ]$pval <- pvals
    }

    df <- df %>%
        rownames_to_column("gene") %>%
        as_tibble()

    df <- df %>%
        mutate(col = case_when(
            diff >= diff_thresh & pval <= pval_thresh ~ "darkred",
            diff <= -diff_thresh & pval <= pval_thresh ~ "darkblue",
            TRUE ~ "gray"
        ))

    # arrange in order for the significant genes to apear above the gray
    df <- df %>%
        mutate(col = factor(col, levels = c("gray", "darkred", "darkblue"))) %>%
        arrange(col)

    return(df)
}


#' calculate mc mc gene expression dataframe
#'
#' @param dataset name of metacell object
#' @param metacell1 id of the first metacell
#' @param metacell2 id of the second metacell
#'
#' @noRd
calc_mc_mc_gene_df <- function(dataset, metacell1, metacell2, diff_thresh = 1.5, pval_thresh = 0.01) {
    mat <- get_mc_data(dataset, "mc_mat")

    egc <- get_metacells_egc(c(metacell1, metacell2), dataset) + egc_epsilon

    df <- calc_diff_expr(mat, egc, c(metacell1, metacell2), diff_thresh, pval_thresh)


    return(df)
}

#' calculate cell type / cell type gene expression dataframe
#'
#' @param dataset name of metacell object
#' @param cell_type1 id of the first cell type
#' @param cell_type2 id of the second cell type
#'
#' @noRd
calc_ct_ct_gene_df <- function(dataset, cell_type1, cell_type2, metacell_types, diff_thresh = 1.5, pval_thresh = 0.01) {
    mat <- get_cell_types_mat(c(cell_type1, cell_type2), metacell_types, dataset)
    egc <- get_cell_types_egc(c(cell_type1, cell_type2), metacell_types, dataset) + egc_epsilon

    df <- calc_diff_expr(mat, egc, c(cell_type1, cell_type2), diff_thresh, pval_thresh)

    return(df)
}

#' calculate sample / sample gene expression dataframe
#'
#' @param dataset name of metacell object
#' @param samp1 id of the first cell type
#' @param samp2 id of the second cell type
#'
#' @noRd
calc_samp_samp_gene_df <- function(dataset, samp1, samp2, metacell_types, cell_types, diff_thresh = 1.5, pval_thresh = 0.01) {
    mat <- get_samples_mat(cell_types, metacell_types, dataset)
    egc <- get_samples_egc(cell_types, metacell_types, dataset) + egc_epsilon

    df <- calc_diff_expr(mat, egc, c(samp1, samp2), diff_thresh, pval_thresh)

    return(df)
}

calc_obs_exp_mc_df <- function(dataset, metacell, diff_thresh = 1.2, pval_thresh = 0.05) {
    obs_mat <- get_mc_data(dataset, "mc_mat")
    exp_mat <- get_mc_data(dataset, "projected_mat")
    obs_egc <- get_metacells_egc(metacell, dataset) + egc_epsilon
    exp_egc <- get_metacells_egc(metacell, dataset, projected = TRUE) + egc_epsilon

    genes <- intersect(rownames(obs_mat), rownames(exp_mat))
    mat <- as.matrix(data.frame(Observed = obs_mat[genes, 1], Projected = exp_mat[, 1]))
    egc <- as.matrix(data.frame(Observed = obs_egc[genes, 1], Projected = exp_egc[, 1]))

    df <- calc_diff_expr(mat, egc, c("Observed", "Projected"), diff_thresh, pval_thresh)

    return(df)
}

calc_obs_exp_type_df <- function(dataset, cell_type, metacell_types, diff_thresh = 1.2, pval_thresh = 0.05) {
    obs_mat <- get_cell_types_mat(cell_type, metacell_types, dataset)
    exp_mat <- get_cell_types_mat(cell_type, metacell_types, dataset, projected = TRUE)

    obs_egc <- get_cell_types_egc(cell_type, metacell_types, dataset) + egc_epsilon
    exp_egc <- get_cell_types_egc(cell_type, metacell_types, dataset, projected = TRUE) + egc_epsilon

    genes <- intersect(rownames(obs_mat), rownames(exp_mat))
    mat <- as.matrix(data.frame(Observed = obs_mat[genes, 1], Projected = exp_mat[, 1]))
    egc <- as.matrix(data.frame(Observed = obs_egc[genes, 1], Projected = exp_egc[, 1]))

    df <- calc_diff_expr(mat, egc, c("Observed", "Projected"), diff_thresh, pval_thresh)

    return(df)
}
