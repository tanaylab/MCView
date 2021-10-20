#' calculate the top k marker genes for each metacell
#'
#' @param mc_egc egc matrix (normalized metacell counts per gene)
#' @param minimal_max_log_fraction take only genes with at least one value
#' (in log fraction units - normalized egc) above this threshold
#' @param minimal_relative_log_fraction take only genes with relative
#' log fraction (mc_fp) above this this value
#'
#' @noRd
calc_marker_genes <- function(mc_egc,
                              genes_per_metacell = 2,
                              minimal_max_log_fraction = -10,
                              minimal_relative_log_fraction = 2,
                              fold_change_reg = 0.1) {
    max_log_fractions_of_genes <- apply(mc_egc, 1, max)

    interesting_genes_mask <- (max_log_fractions_of_genes
    >= minimal_max_log_fraction)

    mc_egc_norm <- mc_egc + 1e-5
    mc_fp <- mc_egc_norm / apply(mc_egc_norm, 1, median, na.rm = TRUE)

    mc_fp_f <- mc_fp[interesting_genes_mask, ] + fold_change_reg
    mc_fp_f[mc_fp_f < minimal_relative_log_fraction] <- NA

    mc_top_genes <- apply(mc_fp_f, 2, function(fp) {
        top_ind <- order(-fp)[1:genes_per_metacell]
        return(tibble(gene = rownames(mc_fp_f)[top_ind], rank = 1:length(top_ind), fp = fp[top_ind]))
    }) %>%
        purrr::imap_dfr(~ .x %>% mutate(metacell = .y)) %>%
        arrange(metacell, rank) %>%
        select(metacell, gene, rank, fp)

    return(mc_top_genes)
}

choose_markers <- function(marker_genes, max_markers) {
    markers <- marker_genes %>%
        group_by(metacell) %>%
        slice(2) %>%
        pull(gene) %>%
        unique()

    # make sure we do not have more than max markers
    if (length(markers) > max_markers) {
        markers <- markers[1:max_markers]
    } else if (length(markers) < max_markers) {
        more_genes <- marker_genes %>%
            filter(!(gene %in% markers)) %>%
            arrange(desc(fp)) %>%
            slice(1:(max_markers - length(markers))) %>%
            pull(gene) %>%
            unique()
        markers <- c(markers, more_genes)
    }

    markers <- sort(unique(markers))

    return(markers)
}

#' Order metacells based on 2 most variable genes
#'
#'
#' @return named vector with metacell order
order_mc_by_most_var_genes <- function(gene_folds, marks = NULL, filter_markers = FALSE, force_cell_type = FALSE, metacell_types = NULL) {
    if (filter_markers) {
        good_marks <- unique(as.vector(unlist(
            apply(
                gene_folds,
                2,
                function(x) {
                    names(head(sort(-x[x > 0.5]), n = 10))
                }
            )
        )))

        if (is.null(good_marks) | length(good_marks) < 4) {
            good_marks <- rownames(gene_folds)
        }
    } else {
        good_marks <- rownames(gene_folds)
    }

    feat <- log2(gene_folds[good_marks, ])

    if (is.null(marks)) {
        g_ncover <- apply(feat > 1, 1, sum)
        main_mark <- names(g_ncover)[which.max(g_ncover)]
        f <- feat[main_mark, ] < 0.25
        if (sum(f) > 0.2 * ncol(feat)) {
            g_score <- apply(feat[, f] > 1, 1, sum)
        } else {
            g_score <- -apply(feat, 1, cor, feat[main_mark, ])
        }
        second_mark <- names(g_score)[which.max(g_score)]
        cli_alert_info("Ordering metacells based on {.file {main_mark}} vs {.file {second_mark}}")
    } else {
        main_mark <- marks[1]
        second_mark <- marks[2]
    }

    if (ncol(feat) == 1) {
        return(1)
    }

    hc <- hclust(tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE)), method = "ward.D2")

    d <- reorder(
        as.dendrogram(hc),
        feat[main_mark, ] - feat[second_mark, ],
        agglo.FUN = mean
    )
    ord <- as.hclust(d)$order
    names(ord) <- colnames(feat)

    if (force_cell_type) {
        ord_df <- enframe(ord, "metacell", "glob_ord") %>%
            left_join(
                metacell_types %>% select(metacell, cell_type),
                by = "metacell"
            )

        ord_inside_cell_type <- function(x) {
            if (nrow(x) == 1) {
                ct_ord <- 1
            } else {
                ct_ord <- suppressWarnings(order_mc_by_most_var_genes(gene_folds[, x$metacell]))
            }

            tibble(
                metacell = x$metacell,
                orig_ord = x$orig_ord,
                glob_ord = x$glob_ord,
                ct_ord = ct_ord
            )
        }

        ord <- ord_df %>%
            mutate(orig_ord = 1:n()) %>%
            group_by(cell_type) %>%
            do(ord_inside_cell_type(.)) %>%
            mutate(glob_ord = mean(glob_ord)) %>%
            ungroup() %>%
            arrange(glob_ord, ct_ord) %>%
            select(metacell, orig_ord) %>%
            deframe()

        # The commented strategy is faster - only order according to the
        # global markers
        # ord <- ord_df %>%
        #     mutate(orig_ord = 1:n()) %>%
        #     group_by(cell_type) %>%
        #     mutate(ct_ord = mean(glob_ord)) %>%
        #     ungroup() %>%
        #     arrange(ct_ord, glob_ord) %>%
        #     select(metacell, orig_ord) %>%
        #     deframe()
    }

    return(ord)
}
