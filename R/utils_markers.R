choose_markers <- function(marker_genes, max_markers = 80) {
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
order_mc_by_most_var_genes <- function(gene_folds, marks = NULL, filter_markers = FALSE) {
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
    return(ord)
}
