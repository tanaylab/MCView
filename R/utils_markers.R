#' calculate the top k marker genes for each metacell
#'
#' @param mc_egc egc matrix (normalized metacell counts per gene)
#' @param minimal_max_log_fraction take only genes with at least one value
#' (in log fraction units - normalized egc) above this threshold
#' @param minimal_relative_log_fraction take only genes with relative
#' log fraction (mc_fp) above this this value
#' @param use_abs take top k that are both enriched and anti-enriched
#'
#' @noRd
calc_marker_genes <- function(mc_egc,
                              genes_per_metacell = 2,
                              minimal_max_log_fraction = -13,
                              minimal_relative_log_fraction = 2,
                              fold_change_reg = 0.1,
                              use_abs = TRUE,
                              daf_obj = NULL) {
    mc_egc <- dafr::fast_log(mc_egc, eps = 1e-5, base = 2)

    max_log_fractions_of_genes <- matrixStats::rowMaxs(mc_egc)

    interesting_genes_mask <- (max_log_fractions_of_genes
    >= minimal_max_log_fraction)

    if (length(interesting_genes_mask) == 0) {
        cli_abort("No genes with at least one value above the {.field minimal_max_log_fraction} threshold ({.val {.value minimal_max_log_fraction}})")
    }

    mc_fp <- sweep(mc_egc, 1, matrixStats::rowMedians(mc_egc, na.rm = TRUE))

    mc_top_genes <- top_fold_per_col(
        fp = mc_fp[interesting_genes_mask, , drop = FALSE],
        k = genes_per_metacell,
        min_val = minimal_relative_log_fraction,
        fold_reg = fold_change_reg,
        use_abs = use_abs
    )
    # Drop NA-padded rows (happens when n_rows < k for some metacell)
    mc_top_genes <- mc_top_genes[!is.na(mc_top_genes$gene), , drop = FALSE]

    return(mc_top_genes)
}

select_top_fold_genes_per_metacell <- function(fold_matrix, genes_per_metacell = 2, minimal_relative_log_fraction = 2, fold_change_reg = 0.1, use_abs = TRUE) {
    fold_matrix <- fold_matrix + fold_change_reg
    fold_matrix[fold_matrix < minimal_relative_log_fraction] <- NA

    if (nrow(fold_matrix) == 0) {
        cli_abort("No genes with at least one value above the {.field minimal_relative_log_fraction} threshold ({.val {.value minimal_relative_log_fraction}})")
    }

    fm <- fold_matrix
    gene_nms <- rownames(fold_matrix)
    mc_nms <- colnames(fold_matrix)
    n_mc <- ncol(fm)

    sort_mat <- if (use_abs) abs(fm) else fm
    results <- vector("list", n_mc)
    for (j in seq_len(n_mc)) {
        col_vals <- sort_mat[, j]
        top_ind <- order(-col_vals, na.last = TRUE)[seq_len(genes_per_metacell)]
        results[[j]] <- tibble(
            metacell = mc_nms[j],
            gene = gene_nms[top_ind],
            rank = seq_len(genes_per_metacell),
            fp = fm[top_ind, j]
        )
    }

    mc_top_genes <- do.call(rbind, results)
    return(mc_top_genes)
}

#' Calculate top marker genes per cell type
#'
#' For each cell type, finds genes most enriched compared to the global mean.
#' Uses log2 fold-change of type mean EGC vs global mean EGC.
#'
#' @param mc_egc EGC matrix (genes x metacells), e.g. from convert_daf_mc_egc()
#' @param mc_types Named character vector: names = metacell IDs, values = cell type labels
#' @param genes_per_type Number of top marker genes to keep per type (default 50)
#' @param min_log_fraction Minimum log2 expression threshold for a gene to be
#'   considered (default -13); genes whose type-mean log2 EGC is below this
#'   are excluded from that type's ranking
#'
#' @return Tibble with columns: cell_type, gene, rank, fold_change
#' @noRd
calc_type_marker_genes <- function(mc_egc,
                                   mc_types,
                                   genes_per_type = 50,
                                   min_log_fraction = -13) {
    # Subset to metacells present in both EGC matrix and type annotation
    common_mcs <- intersect(colnames(mc_egc), names(mc_types))
    if (length(common_mcs) == 0) {
        cli_abort("No metacells in common between EGC matrix and mc_types")
    }
    mc_egc <- mc_egc[, common_mcs, drop = FALSE]
    mc_types <- mc_types[common_mcs]

    # Global mean EGC per gene (across all metacells)
    global_mean <- rowMeans(mc_egc)

    # Unique cell types
    types <- sort(unique(mc_types))

    results <- vector("list", length(types))
    for (i in seq_along(types)) {
        ct <- types[i]
        ct_mcs <- names(mc_types)[mc_types == ct]

        # Type mean EGC per gene
        if (length(ct_mcs) == 1) {
            type_mean <- mc_egc[, ct_mcs, drop = TRUE]
        } else {
            type_mean <- rowMeans(mc_egc[, ct_mcs, drop = FALSE])
        }

        # Log2 fold-change vs global (with epsilon to avoid log(0))
        eps <- 1e-5
        log2_fc <- log2(type_mean + eps) - log2(global_mean + eps)

        # Filter: gene must have sufficient expression in this type
        log2_expr <- log2(type_mean + eps)
        keep <- log2_expr >= min_log_fraction
        if (sum(keep) == 0) next

        log2_fc_filtered <- log2_fc[keep]

        # Take top genes_per_type by absolute fold change
        k <- min(genes_per_type, length(log2_fc_filtered))
        top_idx <- order(-abs(log2_fc_filtered))[seq_len(k)]

        results[[i]] <- tibble(
            cell_type = ct,
            gene = names(log2_fc_filtered)[top_idx],
            rank = seq_len(k),
            fold_change = as.numeric(log2_fc_filtered[top_idx])
        )
    }

    bind_rows(results)
}

choose_markers <- function(marker_genes, max_markers, dataset = NULL) {
    if (is.null(marker_genes)) {
        return(character(0))
    }

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

    return(gene_folds[good_marks, , drop = FALSE])
}

get_top_marks <- function(feat, notify_var_genes = TRUE) {
    g_ncover <- rowSums(feat > 1)
    main_mark <- names(g_ncover)[which.max(g_ncover)]
    f <- feat[main_mark, , drop = FALSE] < 0.25
    if (sum(f) > 0.2 * ncol(feat)) {
        g_score <- rowSums(feat[, f, drop = FALSE] > 1)
    } else {
        g_score <- -apply(feat, 1, cor, feat[main_mark, ])
    }
    second_mark <- names(g_score)[which.max(g_score)]
    if (notify_var_genes) {
        cli_alert_info("Ordering metacells based on {.file {main_mark}} vs {.file {second_mark}}")
    }

    return(c(main_mark, second_mark))
}

# Environment for caching per-cell-type clustering results.
# Key = digest of (cell_type, gene_names, metacell_ids), value = ordering vector.
# This avoids re-computing fastcluster::hclust when toggling force_cell_type
# back and forth, since the per-cell-type subsets do not change.
.ct_cluster_cache <- new.env(parent = emptyenv())

#' Invalidate the per-cell-type clustering cache.
#'
#' Called when markers or data change, ensuring stale results are not reused.
#' @noRd
clear_ct_cluster_cache <- function() {
    rm(list = ls(.ct_cluster_cache, all.names = TRUE), envir = .ct_cluster_cache)
}

# Maximum number of cached entries to prevent unbounded memory growth.
# Each entry is a small integer vector (ordering), so 500 entries is negligible.
.ct_cluster_cache_max <- 500L

#' Order metacells based on 2 most variable genes
#'
#' @noRd
#' @return named vector with metacell order
order_mc_by_most_var_genes <- function(gene_folds, marks = NULL, filter_markers = FALSE, force_cell_type = FALSE, metacell_types = NULL, order_each_cell_type = FALSE, epsilon = 0, notify_var_genes = FALSE, log_transform = TRUE, cached_dist = NULL, secondary_order = NULL) {
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
        feat <- feat_all[, !zero_mcs, drop = FALSE]
        feat_zero <- feat_all[, zero_mcs, drop = FALSE]
    }

    if (ncol(feat) == 0) { # all the metacells do not have enough non-zero values
        return(1:ncol(feat_all))
    }

    if (!is.null(cached_dist)) {
        metacells <- colnames(feat)
        hc <- fastcluster::hclust(as.dist(as.matrix(cached_dist)[metacells, metacells]), method = "ward.D2")
    } else {
        hc <- fastcluster::hclust(tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE)), method = "ward.D2")
    }

    reorder.hclust <- utils::getFromNamespace("reorder.hclust", "vegan")
    ord <- reorder.hclust(hc, as.numeric(feat[main_mark, ] - feat[second_mark, ]), agglo.FUN = "uwmean")$order

    if (any(zero_mcs)) {
        mc_order <- c(colnames(feat)[ord], colnames(feat_zero))
        feat <- feat_all
    } else {
        mc_order <- colnames(feat)[ord]
    }


    if (force_cell_type) {
        ord_df <- tibble(metacell = colnames(feat)) %>%
            mutate(orig_ord = seq_len(n())) %>%
            left_join(
                tibble(metacell = mc_order) %>% mutate(glob_ord = seq_len(n())),
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
                # Build a cache key from (cell_type, sorted gene names, sorted metacell ids)
                ct_name <- x$cell_type[1]
                ct_mcs <- sort(x$metacell)
                ct_genes <- sort(rownames(gene_folds))
                cache_key <- rlang::hash(list(ct_name, ct_genes, ct_mcs))

                if (exists(cache_key, envir = .ct_cluster_cache, inherits = FALSE)) {
                    ct_ord <- get(cache_key, envir = .ct_cluster_cache, inherits = FALSE)
                } else {
                    ct_ord <- suppressWarnings(order_mc_by_most_var_genes(gene_folds[, x$metacell], notify_var_genes = FALSE))
                    # Evict all entries if cache exceeds size limit
                    if (length(ls(.ct_cluster_cache, all.names = TRUE)) >= .ct_cluster_cache_max) {
                        clear_ct_cluster_cache()
                    }
                    assign(cache_key, ct_ord, envir = .ct_cluster_cache)
                }
            }

            tibble(
                metacell = x$metacell,
                orig_ord = x$orig_ord,
                glob_ord = x$glob_ord,
                ct_ord = ct_ord,
                rand = x$rand
            )
        }

        # Deterministic tiebreaker for cell types with identical mean(glob_ord).
        # Using a hash-derived value instead of runif ensures stable ordering
        # across repeated calls (important for caching and toggle consistency).
        ord_df <- ord_df %>% left_join(ord_df %>%
            distinct(cell_type) %>%
            mutate(rand = vapply(cell_type, function(ct) {
                strtoi(substr(rlang::hash(ct), 1, 7), 16L) * 1e-12
            }, numeric(1))), by = "cell_type")
        if (order_each_cell_type) {
            ord <- ord_df %>%
                group_by(cell_type) %>%
                do(ord_inside_cell_type(.)) %>%
                mutate(glob_ord = mean(glob_ord) + rand) %>%
                ungroup() %>%
                arrange(glob_ord, ct_ord) %>%
                pull(orig_ord)
        } else { # only order according to the global markers
            ord <- ord_df %>%
                mutate(orig_ord = 1:n()) %>%
                group_by(cell_type) %>%
                mutate(ct_ord = mean(glob_ord) + rand) %>%
                ungroup()
            if (!is.null(secondary_order)) {
                ord <- ord %>%
                    left_join(tibble(metacell = secondary_order, sec_ord = seq_along(secondary_order)), by = "metacell") %>%
                    arrange(ct_ord, sec_ord) %>%
                    pull(orig_ord)
            } else {
                ord <- ord %>%
                    arrange(ct_ord, glob_ord) %>%
                    pull(orig_ord)
            }
        }
    }

    return(ord)
}

get_markers <- function(dataset) {
    marker_genes <- get_mc_data(dataset, "marker_genes")
    if (!is.null(marker_genes)) {
        return(marker_genes)
    }

    # Calculate marker genes on the fly if not in DAF
    daf_obj <- get_dataset_daf(dataset)
    marker_genes <- calc_marker_genes(get_mc_egc(dataset), 20, daf_obj = daf_obj)

    # Store in session memory cache
    mc_data <- mcv_get("mc_data")
    mc_data[[dataset]][["marker_genes"]] <- marker_genes
    mcv_set("mc_data", mc_data)

    return(marker_genes)
}

get_marker_matrix <- function(dataset, markers, cell_types = NULL, metacell_types = NULL, cell_type_colors = NULL, gene_modules = NULL, force_cell_type = TRUE, mode = "Markers", notify_var_genes = FALSE, metadata_order = NULL, cell_type_metadata_order = NULL, recalc = FALSE, metacells = NULL) {
    cached_dist <- NULL
    if (mode == "Proj") {
        mc_fp <- get_mc_data(dataset, "projected_fold")
        req(mc_fp)
        # Subset to markers first, then drop zero-only rows. Avoids a full-matrix
        # rowSums pass on the (typically sparse) projected_fold and keeps the
        # densification scoped to the marker subset.
        candidate_rows <- intersect(markers, rownames(mc_fp))
        sub <- mc_fp[candidate_rows, , drop = FALSE]
        nz <- if (methods::is(sub, "sparseMatrix")) {
            Matrix::rowSums(abs(sub)) > 0
        } else {
            matrixStats::rowSums2(abs(sub)) > 0
        }
        mc_fp <- as.matrix(sub[nz, , drop = FALSE])
        epsilon <- 1e-5
        log_transform <- FALSE
    } else if (mode == "Markers") {
        if (recalc) {
            if (!is.null(cell_types)) {
                ct_metacells <- metacell_types %>%
                    filter(cell_type %in% cell_types) %>%
                    pull(metacell)
                if (!is.null(metacells)) {
                    metacells <- intersect(metacells, ct_metacells)
                } else {
                    metacells <- ct_metacells
                }
            }
            mc_fp <- get_mc_fp(dataset, markers, metacells = metacells)
        } else {
            mc_fp <- get_mc_fp(dataset, markers)
        }

        epsilon <- 0
        log_transform <- TRUE
        default_markers <- get_mc_data(dataset, "default_markers")
        if (!is.null(default_markers) && all(markers %in% default_markers) && is.null(metacells)) {
            cached_dist <- get_mc_data(dataset, "default_markers_dist")
        }
    } else if (mode == "Gene modules") {
        mc_fp <- get_mc_gene_modules_fp(dataset, markers, gene_modules)
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
        # order cell types
        cell_type_order <- NULL
        if (force_cell_type && !is.null(cell_type_metadata_order) && cell_type_metadata_order != "" && cell_type_metadata_order != "Hierarchical-Clustering" && cell_type_metadata_order %in% c(dataset_metadata_fields_numeric(dataset), "Colors table")) {
            req(!is.null(metacell_types))
            req(!is.null(cell_type_colors))
            if (cell_type_metadata_order == "Colors table") {
                cell_type_order <- cell_type_colors$cell_type
            } else {
                cell_type_order <- get_mc_data(dataset, "metadata") %>%
                    left_join(metacell_types, by = "metacell") %>%
                    group_by(cell_type) %>%
                    summarise(ord_mean = mean(!!sym(cell_type_metadata_order))) %>%
                    arrange(ord_mean) %>%
                    pull(cell_type)
            }
        }

        # order within cell types (if force_cell_type) / global (if !force_cell_type)
        secondary_order <- NULL
        if (!is.null(metadata_order) && metadata_order != "" && metadata_order != "Hierarchical-Clustering" && metadata_order %in% dataset_metadata_fields_numeric(dataset)) {
            secondary_order <- get_mc_data(dataset, "metadata") %>%
                arrange(!!sym(metadata_order)) %>%
                pull(metacell)
        }

        if (!is.null(secondary_order) && !force_cell_type) {
            mc_order <- intersect(secondary_order, colnames(mat))
            # transform to match the order of mat
            mc_order <- match(mc_order, colnames(mat))
        } else {
            mc_order <- order_mc_by_most_var_genes(mat, force_cell_type = force_cell_type, metacell_types = metacell_types, epsilon = epsilon, notify_var_genes = notify_var_genes, log_transform = log_transform, cached_dist = cached_dist, secondary_order = secondary_order)

            if (!is.null(cell_type_order)) {
                mc_order_names <- tibble(metacell = colnames(mat)[mc_order]) %>%
                    left_join(metacell_types, by = "metacell") %>%
                    mutate(cell_type = factor(cell_type, levels = cell_type_order)) %>%
                    arrange(cell_type) %>%
                    pull(metacell)
                mc_order <- match(mc_order_names, colnames(mat))
            }
        }

        mat <- mat[, mc_order, drop = FALSE]
    }

    if (mode %in% c("Markers", "Gene modules")) {
        mat <- log2(mat)
    }

    gene_ord <- order(max.col(mat, ties.method = "first"))
    mat <- mat[gene_ord, , drop = FALSE]

    return(mat)
}

add_genes_to_marker_matrix <- function(mat, genes, dataset) {
    genes_fp <- log2(get_mc_fp(dataset, genes))[, colnames(mat), drop = FALSE]
    genes_fp <- genes_fp[rev(rownames(genes_fp)), , drop = FALSE]
    m <- rbind(genes_fp, mat)
    return(m)
}

get_marker_genes <- function(dataset, mode = "Markers") {
    daf_obj <- get_dataset_daf(dataset)
    if (mode == "Proj") {
        if (is.null(get_mc_data(dataset, "projected_fold"))) {
            if ("Projected-fold" %in% app_config("original_tabs")) {
                showNotification(glue("Projected-fold matrix was not computed. Please compute it in python using the metacells package and rerun the import"), type = "error")
            }
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_projected"))
    } else {
        return(get_markers(dataset))
    }
}
