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

        # chisq.test(matrix(c(m[i, ], tots), nrow = 2))$p.value of every row i (with Yates' continuity
        # correction), computed for all rows at once: calling chisq.test per gene took ~0.3s for 3K
        # genes. The sums are rowSums of the cells in the order chisq.test sums them, so the p-values
        # are the same to the last bit. Unlike chisq.test, this doesn't stop on the tiny negative
        # values (down to -1e-5) of projected counts.
        x <- cbind(m[, 1], m[, 2], tots[1], tots[2]) # the cells of each table, column by column
        n <- rowSums(x)
        sr <- cbind(rowSums(x[, c(1, 3), drop = FALSE]), rowSums(x[, c(2, 4), drop = FALSE]))
        sc <- cbind(rowSums(x[, c(1, 2), drop = FALSE]), rowSums(x[, c(3, 4), drop = FALSE]))
        e <- cbind(sr[, 1] * sc[, 1], sr[, 2] * sc[, 1], sr[, 1] * sc[, 2], sr[, 2] * sc[, 2]) / n
        d <- abs(x - e)
        yates <- pmin(0.5, d[, 1], d[, 2], d[, 3], d[, 4])
        pvals <- stats::pchisq(rowSums((d - yates)^2 / e), 1, lower.tail = FALSE)

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
    egc <- get_cell_types_egc(c(cell_type1, cell_type2), metacell_types, dataset, mat = mat) + egc_epsilon

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
    samples <- unique(c(samp1, samp2))
    mat <- get_samples_mat(cell_types, metacell_types, dataset, samples = samples)
    egc <- get_samples_egc(cell_types, metacell_types, dataset, samples = samples) + egc_epsilon

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
    # the counts of the chosen metacell (this used the first metacell of the dataset)
    mat <- as.matrix(data.frame(Observed = obs_mat[genes, metacell], Projected = exp_mat[genes, metacell]))
    egc <- as.matrix(data.frame(Observed = obs_egc[genes, 1], Projected = exp_egc[genes, 1]))

    df <- calc_diff_expr(mat, egc, c("Observed", "Projected"), diff_thresh, pval_thresh)

    return(df)
}

calc_obs_exp_type_df <- function(dataset, cell_type, metacell_types, diff_thresh = 1.5, pval_thresh = 0.01) {
    obs_mat <- get_cell_types_mat(cell_type, metacell_types, dataset, corrected = TRUE)
    exp_mat <- get_cell_types_mat(cell_type, metacell_types, dataset, projected = TRUE)

    obs_egc <- get_cell_types_egc(cell_type, metacell_types, dataset, mat = obs_mat, corrected = TRUE) + egc_epsilon
    exp_egc <- get_cell_types_egc(cell_type, metacell_types, dataset, mat = exp_mat, projected = TRUE) + egc_epsilon

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
            req(input$metacell1 %in% metacell_types()$metacell)
            req(input$metacell2 %in% metacell_types()$metacell)
            # (comparing a metacell with itself failed with an error about duplicated column names)
            validate(need(input$metacell1 != input$metacell2, "Please choose two different metacells."))
            calc_mc_mc_gene_df(dataset(), input$metacell1, input$metacell2)
        } else if (input$mode == "Types") {
            req(metacell_types())
            # both types should be present as an annotation for at least one metacell
            req(input$metacell1 %in% cell_type_colors()$cell_type)
            req(input$metacell1 %in% metacell_types()$cell_type)
            req(input$metacell2 %in% cell_type_colors()$cell_type)
            req(input$metacell2 %in% metacell_types()$cell_type)
            validate(need(input$metacell1 != input$metacell2, "Please choose two different cell types."))
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
            # only the metacells of the dataset: req(groupA() %in% ...) passed when any of them was, and the
            # others (e.g. pasted from a clipboard of another dataset) broke the comparison
            group_a <- intersect(groupA(), metacell_types()$metacell)
            group_b <- intersect(groupB(), metacell_types()$metacell)
            validate(need(
                length(group_a) > 0 && length(group_b) > 0,
                "Please add metacells to group A and to group B: select them in the 2D projection (\"Select: Group A / Group B\"), add them with the buttons on the left, or paste them from the clipboard."
            ))
            group_types_df <- bind_rows(
                tibble(metacell = group_a, cell_type = "Group A"),
                tibble(metacell = group_b, cell_type = "Group B")
            )
            if (!is.null(metacell_filter)) {
                group_types_df <- group_types_df %>% filter(metacell %in% metacell_filter)
                validate(need(
                    all(c("Group A", "Group B") %in% group_types_df$cell_type),
                    "\"Filter by clipboard\" is on, and group A or group B has no metacells in the clipboard."
                ))
            }

            calc_ct_ct_gene_df(dataset(), "Group A", "Group B", group_types_df)
        }
    }) %>% bindCache(
        # the comparison of two cell types or groups of a large dataset takes seconds
        dataset(),
        input$mode,
        if (isTRUE(input$mode == "MCs")) list(input$metacell1, input$metacell2),
        if (isTRUE(input$mode == "Types")) list(input$metacell1, input$metacell2, cell_type_colors()$cell_type, metacell_types()),
        if (isTRUE(input$mode == "Groups")) list(if (!is.null(groupA)) groupA(), if (!is.null(groupB)) groupB()),
        # the metacells (all the modes check that the compared ones exist)
        metacell_types()$metacell,
        if (isTRUE(input$filter_by_clipboard)) list(TRUE, globals$clipboard)
    )
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

diff_expr_auto_update_globals <- function(mc_mc_gene_scatter_df, globals) {
    observe({
        req(mc_mc_gene_scatter_df())
        globals$significant_genes <- mc_mc_gene_scatter_df() %>%
            filter(col %in% c("darkred", "darkblue")) %>%
            pull(gene) %>%
            unique()
    })
}
