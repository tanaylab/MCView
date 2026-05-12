calc_diff_expr <- function(mat, egc, columns, diff_thresh = 1.5, pval_thresh = 0.01) {
    df <- as.data.frame(egc)

    df$diff <- log2(df[, columns[1]]) - log2(df[, columns[2]])

    f <- rownames(df)[abs(df$diff) >= diff_thresh]

    m <- mat[f, columns, drop = FALSE]

    df$pval <- NA
    if (nrow(m) > 0) {
        tots <- colSums(m)

        # Vectorized chi-squared test for 2x2 contingency tables with
        # Yates' continuity correction (matching stats::chisq.test default
        # for 2x2 tables).
        # Each gene's table is:  matrix(c(x[1], x[2], tots[1], tots[2]), nrow=2)
        # i.e.  a=x[1], c=x[2], b=tots[1], d=tots[2]
        # chi2_yates = (max(|a*d - b*c| - N/2, 0))^2 * N / (r1 * r2 * c1 * c2)
        a <- m[, 1]
        c_val <- m[, 2]
        b <- rep(tots[1], length(a))
        d <- rep(tots[2], length(a))
        N <- a + b + c_val + d
        r1 <- a + b
        r2 <- c_val + d
        c1 <- a + c_val
        c2 <- b + d
        denom <- r1 * r2 * c1 * c2
        ad_bc <- abs(a * d - b * c_val)
        chi2 <- ifelse(denom == 0, 0, pmax(ad_bc - N / 2, 0)^2 * N / denom)
        pvals <- ifelse(denom == 0, 1, stats::pchisq(chi2, df = 1, lower.tail = FALSE))

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
    # Prefer the session-cached full UMIs matrix when available; otherwise
    # query only the two requested metacell columns from DAF. Loading the
    # full matrix here would cost ~1.7 s on OBK; the targeted query is
    # ~30 ms warm / ~700 ms cold.
    mc_mat <- mcv_cache_get(dataset, "mc_mat")
    if (!is.null(mc_mat)) {
        mat <- mc_mat[, c(metacell1, metacell2), drop = FALSE]
    } else {
        daf_obj <- get_daf_for_query(dataset, atlas = FALSE)
        mat <- daf_query_mc_mat(daf_obj, metacells = c(metacell1, metacell2))
    }

    egc <- get_metacells_egc(c(metacell1, metacell2), dataset) + mcv_get("egc_epsilon")

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
    egc <- get_cell_types_egc(c(cell_type1, cell_type2), metacell_types, dataset) + mcv_get("egc_epsilon")

    df <- calc_diff_expr(mat, egc, c(cell_type1, cell_type2), diff_thresh, pval_thresh)

    return(df)
}

#' calculate sample / sample gene expression dataframe
#'
#' @param dataset name of metacell object
#' @param samp1 id of the first sample/group
#' @param samp2 id of the second sample/group
#' @param metacell_types metacell_types tibble
#' @param cell_types cell types to filter by
#' @param group_field optional grouping field for cell-level pseudobulk
#' @param diff_thresh log2 fold-change threshold
#' @param pval_thresh p-value threshold
#'
#' @noRd
calc_samp_samp_gene_df <- function(dataset, samp1, samp2, metacell_types, cell_types,
                                   group_field = NULL, diff_thresh = 1.5, pval_thresh = 0.01) {
    # Cell-level pseudobulk path
    if (!is.null(group_field) && has_cell_gene_umis(dataset)) {
        df <- calc_group_diff_expr(
            dataset, group_field,
            group1_values = samp1,
            group2_values = samp2,
            cell_types = cell_types,
            metacell_types = metacell_types,
            diff_thresh = diff_thresh,
            pval_thresh = pval_thresh
        )
        return(df)
    }

    # Original metacell-weighted approach
    mat <- get_samples_mat(cell_types, metacell_types, dataset)
    egc <- get_samples_egc(cell_types, metacell_types, dataset) + mcv_get("egc_epsilon")

    df <- calc_diff_expr(mat, egc, c(samp1, samp2), diff_thresh, pval_thresh)

    return(df)
}

calc_obs_exp_mc_df <- function(dataset, metacell, diff_thresh = 1.5, pval_thresh = 0.01, corrected = TRUE) {
    daf_obj <- get_dataset_daf(dataset)
    mc_sum_val <- daf_query_mc_sum(daf_obj, metacells = metacell)

    # Query single-metacell UMIs from corrected/projected via DAF instead of loading full matrices
    obs_umis <- get_single_mc_fraction_umis(daf_obj, metacell, mc_sum_val, "corrected_fraction")
    if (is.null(obs_umis)) {
        # Fallback to regular UMIs if no corrected fraction
        obs_umis <- daf_query_mc_mat(daf_obj, metacells = metacell)
        obs_umis <- as.numeric(obs_umis[, 1])
        names(obs_umis) <- dafr::axis_entries(daf_obj, "gene")
    }
    exp_umis <- get_single_mc_fraction_umis(daf_obj, metacell, mc_sum_val, "projected_fraction")
    if (is.null(exp_umis)) {
        return(NULL)
    }

    genes <- intersect(names(obs_umis), names(exp_umis))
    mat <- as.matrix(data.frame(Observed = obs_umis[genes], Projected = exp_umis[genes]))
    rownames(mat) <- genes

    obs_egc <- obs_umis[genes] / sum(obs_umis[genes])
    exp_egc <- exp_umis[genes] / sum(exp_umis[genes])
    egc <- as.matrix(data.frame(
        Observed = obs_egc + mcv_get("egc_epsilon"),
        Projected = exp_egc + mcv_get("egc_epsilon")
    ))
    rownames(egc) <- genes

    df <- calc_diff_expr(mat, egc, c("Observed", "Projected"), diff_thresh, pval_thresh)

    return(df)
}

# Helper: extract a single metacell's UMIs from a fraction matrix via DAF
# Note: get_matrix returns a zero-copy ALTREP view, and column extraction from the
# column-major (gene x metacell) layout is O(n_genes) with contiguous memory access,
# so this is already efficient without needing a DAF single-column query.
get_single_mc_fraction_umis <- function(daf_obj, metacell, mc_sum_val, fraction_name) {
    if (!dafr::has_matrix(daf_obj, "metacell", "gene", fraction_name)) {
        return(NULL)
    }
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    frac_vec <- tryCatch(
        {
            mat <- dafr::get_matrix(daf_obj, "gene", "metacell", fraction_name)
            mat[, metacell]
        },
        error = function(e) NULL
    )
    if (is.null(frac_vec)) {
        return(NULL)
    }
    umis <- as.numeric(frac_vec) * as.numeric(mc_sum_val)
    if (is.null(names(umis))) {
        names(umis) <- gene_names
    }
    umis
}

calc_obs_exp_type_df <- function(dataset, cell_type, metacell_types, diff_thresh = 1.5, pval_thresh = 0.01) {
    obs_mat <- get_cell_types_mat(cell_type, metacell_types, dataset, corrected = TRUE)
    exp_mat <- get_cell_types_mat(cell_type, metacell_types, dataset, projected = TRUE)

    obs_egc <- get_cell_types_egc(cell_type, metacell_types, dataset, corrected = TRUE) + mcv_get("egc_epsilon")
    exp_egc <- get_cell_types_egc(cell_type, metacell_types, dataset, projected = TRUE) + mcv_get("egc_epsilon")

    genes <- intersect(rownames(obs_mat), rownames(exp_mat))
    mat <- as.matrix(data.frame(Observed = obs_mat[genes, 1], Projected = exp_mat[genes, 1]))
    egc <- as.matrix(data.frame(Observed = obs_egc[genes, 1], Projected = exp_egc[genes, 1]))

    df <- calc_diff_expr(mat, egc, c("Observed", "Projected"), diff_thresh, pval_thresh)

    return(df)
}


mc_mc_gene_scatter_df_reactive <- function(dataset, input, output, session, metacell_types, cell_type_colors, state, groupA = NULL, groupB = NULL) {
    reactive({
        req(input$mode)
        if (!is.null(input$filter_by_clipboard) && input$filter_by_clipboard && length(state$selection$clipboard) > 0) {
            metacell_filter <- state$selection$clipboard
        } else {
            metacell_filter <- NULL
        }
        if (input$mode == "MCs") {
            req(input$metacell1)
            req(input$metacell2)
            req(input$metacell1 %in% metacell_types()$metacell)
            req(input$metacell2 %in% metacell_types()$metacell)
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
            req(groupA() %in% metacell_types()$metacell)
            req(groupB() %in% metacell_types()$metacell)
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

diff_expr_auto_update_selection <- function(mc_mc_gene_scatter_df, state) {
    observe({
        req(mc_mc_gene_scatter_df())
        state$selection$significant_genes <- mc_mc_gene_scatter_df() %>%
            filter(col %in% c("darkred", "darkblue")) %>%
            pull(gene) %>%
            unique()
    })
}
