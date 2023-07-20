calc_diff_expr <- function(mat, egc, columns, diff_thresh = 1.5, pval_thresh = 0.01) {
    df <- egc %>%
        as.data.frame()

    df$diff <- log2(df[, columns[1]]) - log2(df[, columns[2]])

    df$pval <- NA

    f <- rownames(df)[abs(df$diff) >= diff_thresh]

    m <- mat[f, columns, drop = FALSE]

    df$pval <- NA
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

calc_obs_exp_mc_df <- function(dataset, metacell, diff_thresh = 1.5, pval_thresh = 0.01, corrected = TRUE) {
    obs_mat <- get_mc_data(dataset, "mc_mat_corrected")
    if (is.null(obs_mat)) {
        obs_mat <- get_mc_data(dataset, "mc_mat")
    }
    exp_mat <- get_mc_data(dataset, "projected_mat")
    obs_egc <- get_metacells_egc(metacell, dataset, corrected = TRUE) + egc_epsilon
    exp_egc <- get_metacells_egc(metacell, dataset, projected = TRUE) + egc_epsilon

    genes <- intersect(rownames(obs_mat), rownames(exp_mat))
    mat <- as.matrix(data.frame(Observed = obs_mat[genes, 1], Projected = exp_mat[, 1]))
    egc <- as.matrix(data.frame(Observed = obs_egc[genes, 1], Projected = exp_egc[, 1]))

    df <- calc_diff_expr(mat, egc, c("Observed", "Projected"), diff_thresh, pval_thresh)

    return(df)
}

calc_obs_exp_type_df <- function(dataset, cell_type, metacell_types, diff_thresh = 1.5, pval_thresh = 0.01) {
    obs_mat <- get_cell_types_mat(cell_type, metacell_types, dataset, corrected = TRUE)
    exp_mat <- get_cell_types_mat(cell_type, metacell_types, dataset, projected = TRUE)

    obs_egc <- get_cell_types_egc(cell_type, metacell_types, dataset, corrected = TRUE) + egc_epsilon
    exp_egc <- get_cell_types_egc(cell_type, metacell_types, dataset, projected = TRUE) + egc_epsilon

    genes <- intersect(rownames(obs_mat), rownames(exp_mat))
    mat <- as.matrix(data.frame(Observed = obs_mat[genes, 1], Projected = exp_mat[genes, 1]))
    egc <- as.matrix(data.frame(Observed = obs_egc[genes, 1], Projected = exp_egc[genes, 1]))

    df <- calc_diff_expr(mat, egc, c("Observed", "Projected"), diff_thresh, pval_thresh)

    return(df)
}


mc_mc_gene_scatter_df_reactive <- function(dataset, input, output, session, metacell_types, cell_type_colors, globals, groupA = NULL, groupB = NULL) {
    reactive({
        req(input$mode)
        if (!is.null(input$filter_by_clipboard) && input$filter_by_clipboard && length(globals$clipboard) > 0) {
            metacell_filter <- globals$clipboard
        } else {
            metacell_filter <- NULL
        }
        if (input$mode == "MCs") {
            req(input$metacell1)
            req(input$metacell2)
            calc_mc_mc_gene_df(dataset(), input$metacell1, input$metacell2)
        } else if (input$mode == "Types") {
            req(metacell_types())
            # both types should be present as an annotation for at least one metacell
            req(input$metacell1 %in% cell_type_colors()$cell_type)
            req(input$metacell1 %in% metacell_types()$cell_type)
            req(input$metacell2 %in% cell_type_colors()$cell_type)
            req(input$metacell2 %in% metacell_types()$cell_type)
            types_df <- metacell_types()
            if (!is.null(metacell_filter)) {
                types_df <- types_df %>%
                    filter(cell_type %in% c(input$metacell1, input$metacell2)) %>%
                    filter(metacell %in% metacell_filter)
                req(nrow(types_df) > 0)
                req(all(c(input$metacell1, input$metacell2) %in% types_df$cell_type))
            }
            calc_ct_ct_gene_df(dataset(), input$metacell1, input$metacell2, types_df)
        } else if (input$mode == "Groups") {
            req(groupA)
            req(groupB)
            req(groupA())
            req(groupB())
            group_types_df <- bind_rows(
                tibble(metacell = groupA(), cell_type = "Group A"),
                tibble(metacell = groupB(), cell_type = "Group B")
            )
            if (!is.null(metacell_filter)) {
                group_types_df <- group_types_df %>% filter(metacell %in% metacell_filter)
                req(all(c("Group A", "Group B") %in% group_types_df$cell_type))
                req(nrow(group_types_df) > 0)
            }

            calc_ct_ct_gene_df(dataset(), "Group A", "Group B", group_types_df)
        }
    })
}

diff_expr_switch_metacells <- function(dataset, input, output, session, groupA = NULL, groupB = NULL) {
    observeEvent(input$switch_metacells, {
        if (input$mode == "Groups") {
            req(groupA)
            req(groupB)
            temp <- groupA()
            groupA(groupB())
            groupB(temp)
        } else {
            mc1 <- input$metacell1
            mc2 <- input$metacell2
            updateSelectInput(session, "metacell1", selected = mc2)
            updateSelectInput(session, "metacell2", selected = mc1)
        }
    })
}
