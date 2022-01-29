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

    mc_top_genes <- select_top_fold_genes(
        mc_fp[interesting_genes_mask, ],
        genes_per_metacell = genes_per_metacell,
        minimal_relative_log_fraction = minimal_relative_log_fraction,
        fold_change_reg = fold_change_reg
    )

    return(mc_top_genes)
}

select_top_fold_genes <- function(fold_matrix, genes_per_metacell = 2, minimal_relative_log_fraction = 2, fold_change_reg = 0.1, use_abs = FALSE) {
    fold_matrix <- fold_matrix + fold_change_reg
    fold_matrix[fold_matrix < minimal_relative_log_fraction] <- NA

    mc_top_genes <- apply(fold_matrix, 2, function(fp) {
        if (use_abs) {
            top_ind <- order(-abs(fp))[1:genes_per_metacell]
        } else {
            top_ind <- order(-fp)[1:genes_per_metacell]
        }

        return(tibble(gene = rownames(fold_matrix)[top_ind], rank = 1:length(top_ind), fp = fp[top_ind]))
    }) %>%
        purrr::imap_dfr(~ .x %>% mutate(metacell = .y)) %>%
        arrange(metacell, rank) %>%
        select(metacell, gene, rank, fp)

    return(mc_top_genes)
}

choose_markers <- function(marker_genes, max_markers, dataset = NULL, add_systematic = FALSE) {
    if (has_name(marker_genes, "metacell")) {
        n_metacells <- length(unique(marker_genes$metacell))
        genes_per_metacell <- min(10, max(1, round(max_markers / n_metacells)))

        markers <- marker_genes %>%
            group_by(metacell) %>%
            slice(1:genes_per_metacell) %>%
            pull(gene) %>%
            unique()
    } else {
        markers <- unique(marker_genes$gene)
    }

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

    if (add_systematic && !is.null(dataset)) {
        markers <- add_systematic_markers(dataset, markers)
        markers <- markers[1:min(length(markers), max_markers)]
    }

    return(markers)
}

add_systematic_markers <- function(dataset, markers) {
    systematic_genes <- get_mc_data(dataset, "systematic_genes")
    if (!is.null(systematic_genes)) {
        m <- get_mc_data(dataset, "projected_fold")
        f <- Matrix::rowSums(m[intersect(rownames(m), systematic_genes), ]) > 0
        return(unique(c(systematic_genes[f], markers)))
    }
    return(markers)
}

filter_markers_mat <- function(gene_folds) {
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

    return(gene_folds[good_marks, ])
}

get_top_marks <- function(feat, notify_var_genes = TRUE) {
    g_ncover <- apply(feat > 1, 1, sum)
    main_mark <- names(g_ncover)[which.max(g_ncover)]
    f <- feat[main_mark, ] < 0.25
    if (sum(f) > 0.2 * ncol(feat)) {
        g_score <- apply(feat[, f] > 1, 1, sum)
    } else {
        g_score <- -apply(feat, 1, cor, feat[main_mark, ])
    }
    second_mark <- names(g_score)[which.max(g_score)]
    if (notify_var_genes) {
        cli_alert_info("Ordering metacells based on {.file {main_mark}} vs {.file {second_mark}}")
    }

    return(c(main_mark, second_mark))
}

#' Order metacells based on 2 most variable genes
#'
#' @noRd
#' @return named vector with metacell order
order_mc_by_most_var_genes <- function(gene_folds, marks = NULL, filter_markers = FALSE, force_cell_type = FALSE, metacell_types = NULL, order_each_cell_type = FALSE, epsilon = 0, notify_var_genes = FALSE, log_transform = TRUE) {
    if (filter_markers) {
        gene_folds <- filter_markers_mat(gene_folds)
    }

    if (log_transform) {
        feat <- log2(gene_folds + epsilon)
    } else {
        feat <- gene_folds
    }


    if (is.null(marks)) {
        marks <- get_top_marks(feat, notify_var_genes = notify_var_genes)
        main_mark <- marks[1]
        second_mark <- marks[2]
    } else {
        main_mark <- marks[1]
        second_mark <- marks[2]
    }

    if (notify_var_genes) {
        showNotification(glue("Ordering metacells based on {main_mark} vs {second_mark}"))
    }

    if (ncol(feat) == 1) {
        return(1)
    }

    zero_mcs <- colSums(abs(feat) > 0) < 2
    if (any(zero_mcs)) {
        feat_all <- feat
        feat <- feat_all[, !zero_mcs]
        feat_zero <- feat_all[, zero_mcs]
    }

    if (ncol(feat) == 0) { # all the metacells do not have enough non-zero values
        return(1:ncol(feat_all))
    }

    hc <- hclust(tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE)), method = "ward.D2")

    d <- reorder(
        as.dendrogram(hc),
        feat[main_mark, ] - feat[second_mark, ],
        agglo.FUN = mean
    )
    ord <- as.hclust(d)$order

    if (any(zero_mcs)) {
        mc_order <- c(colnames(feat)[ord], colnames(feat_zero))
        feat <- feat_all
    } else {
        mc_order <- colnames(feat)[ord]
    }

    if (force_cell_type) {
        ord_df <- tibble(metacell = colnames(feat)) %>%
            mutate(orig_ord = 1:n()) %>%
            left_join(
                tibble(metacell = mc_order) %>% mutate(glob_ord = 1:n()),
                by = "metacell"
            ) %>%
            left_join(
                metacell_types %>% select(metacell, cell_type),
                by = "metacell"
            )

        ord_inside_cell_type <- function(x) {
            if (nrow(x) == 1) {
                ct_ord <- 1
            } else {
                ct_ord <- suppressWarnings(order_mc_by_most_var_genes(gene_folds[, x$metacell], notify_var_genes = FALSE))
            }

            tibble(
                metacell = x$metacell,
                orig_ord = x$orig_ord,
                glob_ord = x$glob_ord,
                ct_ord = ct_ord
            )
        }

        if (order_each_cell_type) {
            ord <- ord_df %>%
                group_by(cell_type) %>%
                do(ord_inside_cell_type(.)) %>%
                mutate(glob_ord = mean(glob_ord)) %>%
                ungroup() %>%
                arrange(glob_ord, ct_ord) %>%
                pull(orig_ord)
        } else { # only order according to the global markers
            ord <- ord_df %>%
                mutate(orig_ord = 1:n()) %>%
                group_by(cell_type) %>%
                mutate(ct_ord = mean(glob_ord)) %>%
                ungroup() %>%
                arrange(ct_ord, glob_ord) %>%
                pull(orig_ord)
        }
    }

    return(ord)
}

get_markers <- function(dataset) {
    marker_genes <- get_mc_data(dataset, "marker_genes")
    if (!is.null(marker_genes)) {
        return(marker_genes)
    }

    marker_genes <- calc_marker_genes(get_mc_egc(dataset), 20)
    serialize_shiny_data(marker_genes, "marker_genes", dataset = dataset, cache_dir = cache_dir)

    return(marker_genes)
}

get_marker_matrix <- function(dataset, markers, cell_types = NULL, metacell_types = NULL, force_cell_type = TRUE, mode = "Markers", notify_var_genes = FALSE) {
    if (mode == "Inner") {
        mc_fp <- get_mc_data(dataset, "inner_fold_mat")
        req(mc_fp)
        mc_fp <- as.matrix(mc_fp[Matrix::rowSums(mc_fp) > 0, ])
        mc_fp <- mc_fp[intersect(markers, rownames(mc_fp)), ]
        epsilon <- 1e-5
        log_transform <- FALSE
    } else if (mode == "Proj") {
        mc_fp <- get_mc_data(dataset, "projected_fold")
        req(mc_fp)
        mc_fp <- as.matrix(mc_fp[abs(Matrix::rowSums(mc_fp)) > 0, ])
        mc_fp <- mc_fp[intersect(markers, rownames(mc_fp)), ]
        epsilon <- 1e-5
        log_transform <- FALSE
    } else if (mode == "Markers") {
        mc_fp <- get_mc_fp(dataset, markers)
        epsilon <- 0
        log_transform <- TRUE
    } else {
        stop("unknown mode")
    }

    req(dim(mc_fp))
    req(nrow(mc_fp) > 0)

    if (!is.null(cell_types)) {
        mat <- filter_mat_by_cell_types(mc_fp, cell_types, metacell_types)
    } else {
        mat <- mc_fp
    }

    if (ncol(mat) > 1) {
        mc_order <- order_mc_by_most_var_genes(mat, force_cell_type = force_cell_type, metacell_types = metacell_types, epsilon = epsilon, notify_var_genes = notify_var_genes, log_transform = log_transform)
        mat <- mat[, mc_order]
    }

    if (mode == "Markers") {
        mat <- log2(mat)
    }

    metacells <- colnames(mat)

    gene_ord <- order(apply(mat, 1, which.max))
    mat <- mat[gene_ord, ]

    # matrix has only a single metacell
    if (!is.matrix(mat)) {
        mat <- as.matrix(mat)
        colnames(mat) <- metacells
    }

    return(mat)
}

get_marker_genes <- function(dataset, mode = "Markers") {
    if (mode == "Inner") {
        if (is.null(get_mc_data(dataset, "inner_fold_mat"))) {
            if ("Inner-fold" %in% config$original_tabs) {
                showNotification(glue("Inner-fold matrix was not computed. Please compute it in python using the metacells package and rerun the import"), type = "error")
            }
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_inner_fold"))
    } else if (mode == "Proj") {
        if (is.null(get_mc_data(dataset, "projected_fold"))) {
            if ("Projected-fold" %in% config$original_tabs) {
                showNotification(glue("Projected-fold matrix was not computed. Please compute it in python using the metacells package and rerun the import"), type = "error")
            }
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_projected"))
    } else {
        return(get_markers(dataset))
    }

    return(get_markers(dataset))
}
