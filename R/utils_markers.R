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
                              use_abs = TRUE) {
    mc_egc <- log2(mc_egc + 1e-5)

    max_log_fractions_of_genes <- matrixStats::rowMaxs(mc_egc)

    interesting_genes_mask <- (max_log_fractions_of_genes
    >= minimal_max_log_fraction)

    if (length(interesting_genes_mask) == 0) {
        cli_abort("No genes with at least one value above the {.field minimal_max_log_fraction} threshold ({.val {.value minimal_max_log_fraction}})")
    }

    mc_fp <- mc_egc - matrixStats::rowMedians(mc_egc, na.rm = TRUE)

    mc_top_genes <- select_top_fold_genes_per_metacell(
        mc_fp[interesting_genes_mask, , drop = FALSE],
        genes_per_metacell = genes_per_metacell,
        minimal_relative_log_fraction = minimal_relative_log_fraction,
        fold_change_reg = fold_change_reg,
        use_abs = use_abs
    )

    return(mc_top_genes)
}

select_top_fold_genes_per_metacell <- function(fold_matrix, genes_per_metacell = 2, minimal_relative_log_fraction = 2, fold_change_reg = 0.1, use_abs = TRUE) {
    fold_matrix <- fold_matrix + fold_change_reg
    fold_matrix[fold_matrix < minimal_relative_log_fraction] <- NA

    if (nrow(fold_matrix) == 0) {
        cli_abort("No genes with at least one value above the {.field minimal_relative_log_fraction} threshold ({.val {.value minimal_relative_log_fraction}})")
    }

    top_ind <- apply(fold_matrix, 2, function(fp) {
        if (use_abs) {
            order(-abs(fp))[1:genes_per_metacell]
        } else {
            order(-fp)[1:genes_per_metacell]
        }
    })
    top_ind <- matrix(top_ind, nrow = genes_per_metacell)

    genes <- rownames(fold_matrix)[top_ind]
    mc_top_genes <- tibble(
        metacell = rep(colnames(fold_matrix), each = genes_per_metacell),
        gene = genes,
        rank = rep(seq_len(genes_per_metacell), ncol(fold_matrix)),
        fp = stats::setNames(fold_matrix[cbind(as.vector(top_ind), as.vector(col(top_ind)))], genes)
    ) %>%
        arrange(metacell, rank)

    return(mc_top_genes)
}

choose_markers <- function(marker_genes, max_markers, dataset = NULL) {
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
    g_ncover <- apply(feat > 1, 1, sum)
    main_mark <- names(g_ncover)[which.max(g_ncover)]
    f <- feat[main_mark, , drop = FALSE] < 0.25
    if (sum(f) > 0.2 * ncol(feat)) {
        g_score <- apply(feat[, f, drop = FALSE] > 1, 1, sum)
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
order_mc_by_most_var_genes <- function(gene_folds, marks = NULL, filter_markers = FALSE, force_cell_type = FALSE, metacell_types = NULL, order_each_cell_type = FALSE, epsilon = 0, notify_var_genes = FALSE, log_transform = TRUE, cached_dist = NULL, secondary_order = NULL, precomputed_clusterings = NULL) {
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

    # The key of the clustering in the app cache and in the clusterings stored with the dataset (calc_default_clusterings).
    # The app can run on another machine than the one that stored them, so the key must not depend on the versions of
    # R and rlang (rlang::hash() changed in rlang 1.3.0), the math library or the locale: it is a hash of the serialized
    # features without the header (which has the R version), taken before the log and with the genes in a fixed order
    # (choose_markers sorts them by the locale).
    key <- list("cluster_metacells", gene_folds[order(rownames(gene_folds), method = "radix"), , drop = FALSE], epsilon, log_transform, !is.null(cached_dist))
    key <- qs2::xxhash_raw(serialize(key, NULL, xdr = FALSE, version = 2)[-(1:14)])
    hc <- cluster_metacells(feat, key, cached_dist, precomputed_clusterings)

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
                ct_ord <- suppressWarnings(order_mc_by_most_var_genes(gene_folds[, x$metacell], notify_var_genes = FALSE))
            }

            tibble(
                metacell = x$metacell,
                orig_ord = x$orig_ord,
                glob_ord = x$glob_ord,
                ct_ord = ct_ord
            )
        }

        # sample a random small number for each cell type and left_join it to ord_df
        ord_df <- ord_df %>% left_join(ord_df %>%
            distinct(cell_type) %>%
            mutate(rand = runif(n(), 0, 1e-5)), by = "cell_type")
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

# The clustering is the slow part of the heatmaps and depends only on the features, while the
# heatmap is recomputed on every zoom / annotation / ordering change. Keep it in the app cache
# (shared by all sessions) when running inside shiny.
# key: identifies the features (see order_mc_by_most_var_genes).
# cached_dist: NULL, or a function returning the distances between the metacells to cluster by instead of the
# features. It is called only when the clustering is neither cached nor stored, since reading them is slow.
# precomputed_clusterings: clusterings stored with the dataset (calc_default_clusterings), by key.
cluster_metacells <- function(feat, key, cached_dist = NULL, precomputed_clusterings = NULL) {
    cache <- shiny::getShinyOption("cache")
    if (!is.null(cache)) {
        hc <- cache$get(key)
        if (!cachem::is.key_missing(hc)) {
            return(hc)
        }
    }

    if (!is.null(precomputed_clusterings[[key]])) {
        return(precomputed_clusterings[[key]])
    }

    # Clusterings shared by all the apps that set MCVIEW_CLUSTERING_CACHE_DIR (e.g. on a server), across processes
    # and restarts, for projects that don't store their own (the key identifies the features, not the project).
    shared <- NULL
    if (nzchar(Sys.getenv("MCVIEW_CLUSTERING_CACHE_DIR"))) {
        shared <- tryCatch(cachem::cache_disk(Sys.getenv("MCVIEW_CLUSTERING_CACHE_DIR"), max_size = 10 * 1024^3), error = function(e) NULL)
        hc <- tryCatch(shared$get(key), error = function(e) cachem::key_missing())
        if (!cachem::is.key_missing(hc)) {
            if (!is.null(cache)) {
                cache$set(key, hc)
            }
            return(hc)
        }
    }

    if (!is.null(cached_dist)) {
        hc <- fastcluster::hclust(subset_dist(cached_dist(), colnames(feat)), method = "ward.D2")
    } else {
        hc <- fastcluster::hclust(cor_dist(feat), method = "ward.D2")
    }

    if (!is.null(cache)) {
        cache$set(key, hc)
    }
    if (!is.null(shared)) {
        try(shared$set(key, hc), silent = TRUE)
    }
    return(hc)
}

# Cluster the metacells of the default view of each heatmap tab whose data exist, since clustering is the slowest
# part of their first render. The result is keyed like the app cache (order_mc_by_most_var_genes), so the app uses
# a clustering only when the heatmap features are unchanged (e.g. not after editing gene modules).
# Reads the dataset from cache_dir as the app does, and returns the clusterings by key.
calc_default_clusterings <- function(dataset, cache_dir) {
    # the heatmap code reads the data through the global mc_data, and stores the clusterings in the shiny cache
    prev_mc_data <- get0("mc_data", envir = globalenv())
    assign("mc_data", list(), envir = globalenv())
    prev_cache <- shiny::getShinyOption("cache")
    cache <- cachem::cache_mem(max_size = Inf)
    shiny::shinyOptions(cache = cache)
    on.exit({
        shiny::shinyOptions(cache = prev_cache)
        if (is.null(prev_mc_data)) {
            rm("mc_data", envir = globalenv())
        } else {
            assign("mc_data", prev_mc_data, envir = globalenv())
        }
    })

    load_all_mc_data(dataset, cache_dir)
    if (has_mc_data(dataset, "default_clusterings")) { # recompute them
        rm("default_clusterings", envir = mc_data[[dataset]])
    }

    # the inputs each tab starts with (see app_server, heatmap_reactives, mod_gene_modules and mod_outliers)
    metacell_types <- get_metacell_types_data(dataset)
    cell_type_colors <- get_cell_type_data(dataset)
    gene_modules <- get_mc_data(dataset, "gene_modules")
    if (!is.null(gene_modules)) {
        gene_modules <- gene_modules %>%
            filter(gene %in% gene_names(dataset)) %>%
            mutate(gene = as.character(gene), module = factor(module))
    }

    # the data each heatmap needs (without marker_genes the Markers tab computes and saves them, which we leave to the app)
    needs <- c("Markers" = "marker_genes", "Gene modules" = "gene_modules", "Inner" = "inner_fold_mat", "Stdev" = "inner_stdev_mat", "Proj" = "projected_fold", "Outliers" = "deviant_fold_mat")
    for (mode in names(needs)[purrr::map_lgl(needs, ~ has_mc_data(dataset, .x))]) {
        tryCatch(
            {
                mode_metacell_types <- metacell_types
                if (mode == "Outliers") {
                    mode_metacell_types <- get_mc_data(dataset, "outliers_metadata") %>%
                        select(cell_id, metacell = most_similar) %>%
                        mutate(metacell = as.character(metacell)) %>%
                        left_join(metacell_types, by = "metacell") %>%
                        mutate(most_similar_metacell = metacell) %>%
                        select(-metacell) %>%
                        rename(metacell = cell_id)
                }
                cell_types <- cell_type_colors %>%
                    filter(cell_type %in% mode_metacell_types$cell_type) %>%
                    pull(cell_type)
                if (mode == "Gene modules") {
                    modules <- levels(gene_modules$module)
                    markers <- modules[modules %in% gene_modules$module]
                } else {
                    markers <- choose_markers(get_marker_genes(dataset, mode = mode), 100, dataset = dataset)
                }
                cli_alert_info("Clustering the metacells of the {.field {mode}} heatmap")
                get_marker_matrix(dataset, markers, cell_types, mode_metacell_types, cell_type_colors, gene_modules, mode = mode)
            },
            error = function(e) {
                cli_alert_warning("Could not cluster the {.field {mode}} heatmap (the app will cluster it on first use): {conditionMessage(e)}")
            }
        )
    }

    keys <- cache$keys()
    return(stats::setNames(lapply(keys, cache$get), keys))
}

#' Store the clusterings of the default heatmaps
#'
#' Cluster the metacells of the heatmap that each heatmap tab (Markers, Gene modules, Inner-fold, Stdev-fold,
#' Projected-fold and Outliers) starts with, and store the clusterings in the project, so that the app doesn't
#' cluster them the first time a tab is opened, which takes minutes for large datasets.
#' \code{import_dataset} already does this. Use this function for projects imported with an older version of
#' MCView, or after changing the gene modules with \code{update_gene_modules}.
#'
#' Under the hood - MCView writes a file named "default_clusterings.qs2" under \code{project/cache/dataset}, and
#' skips datasets whose directory is a symlink to outside the project.
#' When that directory isn't writable and the environment variable \code{MCVIEW_CLUSTERING_CACHE_DIR} is set, the
#' clusterings are stored in that directory instead. Apps that set the variable (e.g. all the apps of a shiny
#' server) also store there every clustering they compute and look for clusterings there, so each clustering is
#' computed once per server rather than once per app process.
#'
#' @param project path to the project directory
#' @param overwrite recompute the clusterings of datasets that already have them
#'
#' @examples
#' \dontrun{
#' update_default_clusterings("PBMC")
#' }
#'
#' @export
update_default_clusterings <- function(project, overwrite = FALSE) {
    project <- init_project_dir(project)
    cache_dir <- project_cache_dir(project)
    project_real <- fs::path_real(project)
    for (dataset in dataset_ls(project)) {
        dataset_dir <- fs::path(cache_dir, dataset)
        # a dataset linked from another project: that project may be served by an older MCView
        dataset_real <- fs::path_real(dataset_dir)
        if (!fs::path_has_parent(dataset_real, project_real)) {
            cli_alert_info("{.field {dataset}} is a symlink to {.file {dataset_real}}, outside the project, skipping it (update that project instead)")
            next
        }
        if (!overwrite && any(fs::file_exists(fs::path(dataset_dir, c("default_clusterings.qs2", "default_clusterings.qs"))))) {
            cli_alert_info("{.field {dataset}} already has default_clusterings, skipping it (use {.code overwrite = TRUE} to recompute them)")
            next
        }
        # a project we can't write to (e.g. another user's) can use the directory shared by the apps instead
        shared_dir <- Sys.getenv("MCVIEW_CLUSTERING_CACHE_DIR")
        writable <- file.access(dataset_dir, 2) == 0
        if (!writable && !nzchar(shared_dir)) {
            cli_alert_warning("Can't write to {.file {dataset_dir}}, skipping {.field {dataset}} (set MCVIEW_CLUSTERING_CACHE_DIR to store its clusterings in a directory shared by the apps)")
            next
        }
        cli_alert_info("Clustering the metacells of the default heatmaps of {.field {dataset}}")
        clusterings <- tryCatch(calc_default_clusterings(dataset, cache_dir), error = function(e) {
            cli_alert_warning(conditionMessage(e))
            list()
        })
        if (length(clusterings) == 0) {
            cli_alert_warning("Could not cluster any heatmap of {.field {dataset}}")
            next
        }
        if (!writable) {
            shared <- cachem::cache_disk(shared_dir, max_size = 10 * 1024^3)
            purrr::iwalk(clusterings, function(hc, key) shared$set(key, hc))
            cli_alert_success("Saved the clusterings of {length(clusterings)} heatmap{?s} of {.field {dataset}} in {.file {shared_dir}}")
            next
        }
        # write to a temporary file and rename it, so that an app that starts meanwhile doesn't read a partial file,
        # and read it back first (some file systems, e.g. FUSE mounts, drop qs2's final write)
        tmp <- tempfile(".default_clusterings.qs2-", tmpdir = dataset_dir) # hidden: MCView doesn't list it
        written <- tryCatch(
            {
                qs2::qs_save(clusterings, tmp)
                identical(qs2::qs_read(tmp, validate_checksum = TRUE), clusterings)
            },
            error = function(e) FALSE
        )
        if (!written) {
            unlink(tmp)
            cli_alert_warning("Could not write the clusterings of {.field {dataset}} to {.file {dataset_dir}}, skipping it")
            next
        }
        fs::file_move(tmp, fs::path(dataset_dir, "default_clusterings.qs2"))
        cli_alert_success("Saved the clusterings of {length(clusterings)} heatmap{?s} of {.field {dataset}}")
    }
}

# Same as tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE)): the euclidean distances between the rows of
# the metacell correlation matrix, without computing it. With z the standardized columns of feat, the correlation
# matrix is t(z) %*% z, so the dot products of its rows are t(z) %*% w with w = z %*% t(z) %*% z. That costs
# n^2 * nrow(feat) instead of n^3 for n metacells.
cor_dist <- function(feat) {
    z <- scale(feat) / sqrt(nrow(feat) - 1)
    if (!all(is.finite(z))) { # missing values or constant metacells
        return(tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE)))
    }
    if (nrow(z) < ncol(z) && isTRUE(getOption("tgs_use.blas"))) {
        # The distance between rows i and j of the correlation matrix is |t(z) %*% (z_i - z_j)| = |M^(1/2) %*% (z_i - z_j)|
        # with M = z %*% t(z), i.e. the distance between columns i and j of y = M^(1/2) %*% z, which tgs_dist computes
        # faster than the loop below when it uses BLAS (and slower when it doesn't).
        e <- eigen(tcrossprod(z), symmetric = TRUE)
        y <- e$vectors %*% (sqrt(pmax(e$values, 0)) * crossprod(e$vectors, z))
        d <- tgs_dist(t(y))
        attributes(d) <- list(Size = ncol(z), Labels = colnames(feat), Diag = FALSE, Upper = FALSE, method = "euclidean", class = "dist")
        return(d)
    }
    w <- tcrossprod(z) %*% z
    sq <- colSums(z * w)
    # the squared distance sq_i + sq_j - 2 * z_i . w_j as a single matrix product
    z <- rbind(z, sq, 1)
    w <- rbind(-2 * w, 1, sq)
    n <- ncol(z)
    d <- numeric(n * (n - 1) / 2)
    pos <- 0
    # by blocks of columns, filling d in place
    for (js in split(seq_len(n - 1), ceiling(seq_len(n - 1) / 200))) {
        d2 <- crossprod(z[, js[1]:n, drop = FALSE], w[, js, drop = FALSE])
        # d(i, j) for i > j, in the column-major order of a dist object
        d2 <- d2[lower.tri(d2)]
        d2[d2 < 0] <- 0 # rounding errors
        d[(pos + 1):(pos + length(d2))] <- sqrt(d2)
        pos <- pos + length(d2)
    }
    attributes(d) <- list(Size = n, Labels = colnames(feat), Diag = FALSE, Upper = FALSE, method = "euclidean", class = "dist")

    return(d)
}

# Same as as.dist(as.matrix(d)[labels, labels]) without expanding d to a full n x n matrix
subset_dist <- function(d, labels) {
    idx <- match(labels, attr(d, "Labels"))
    n <- attr(d, "Size")
    if (length(idx) == n && identical(idx, seq_len(n))) {
        return(d)
    }
    k <- length(idx)
    if (k < 2 || anyNA(idx) || is.unsorted(idx, strictly = TRUE)) {
        return(as.dist(as.matrix(d)[labels, labels]))
    }
    # d(i, j) for i < j is stored at n * (i - 1) - i * (i - 1) / 2 + j - i
    i <- as.numeric(idx[-k])
    offset <- n * (i - 1) - i * (i - 1) / 2 - i
    pos <- unlist(lapply(seq_len(k - 1), function(a) offset[a] + idx[(a + 1):k]))
    structure(d[pos], Size = k, Labels = labels, Diag = FALSE, Upper = FALSE, class = "dist")
}

get_markers <- function(dataset) {
    marker_genes <- get_mc_data(dataset, "marker_genes")
    if (!is.null(marker_genes)) {
        return(marker_genes)
    }

    marker_genes <- calc_marker_genes(get_mc_egc(dataset), 20)
    # keep the markers for this process even when they can't be saved (e.g. a read-only cache dir)
    mc_data[[dataset]][["marker_genes"]] <<- marker_genes
    tryCatch(
        serialize_shiny_data(marker_genes, "marker_genes", dataset = dataset, cache_dir = cache_dir),
        error = function(e) {
            try(unlink(fs::path(cache_dir, dataset, "marker_genes.qs2")), silent = TRUE) # don't leave a partial file
            cli_alert_warning("Could not save {.field marker_genes}: {conditionMessage(e)}")
        }
    )

    return(marker_genes)
}

get_marker_matrix <- function(dataset, markers, cell_types = NULL, metacell_types = NULL, cell_type_colors = NULL, gene_modules = NULL, force_cell_type = TRUE, mode = "Markers", notify_var_genes = FALSE, metadata_order = NULL, cell_type_metadata_order = NULL, recalc = FALSE, metacells = NULL) {
    cached_dist <- NULL
    if (mode == "Inner") {
        mc_fp <- get_mc_data(dataset, "inner_fold_mat")
        req(mc_fp)
        mc_fp <- as.matrix(mc_fp[intersect(markers, rownames(mc_fp)[Matrix::rowSums(mc_fp) > 0]), , drop = FALSE])
        epsilon <- 1e-5
        log_transform <- FALSE
    } else if (mode == "Stdev") {
        mc_fp <- get_mc_data(dataset, "inner_stdev_mat")
        req(mc_fp)
        mc_fp <- as.matrix(mc_fp[intersect(markers, rownames(mc_fp)[Matrix::rowSums(mc_fp) > 0]), , drop = FALSE])
        epsilon <- 1e-5
        log_transform <- FALSE
    } else if (mode == "Proj") {
        mc_fp <- get_mc_data(dataset, "projected_fold")
        req(mc_fp)
        mc_fp <- as.matrix(mc_fp[intersect(markers, rownames(mc_fp)[abs(Matrix::rowSums(mc_fp)) > 0]), , drop = FALSE])
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
        if (!is.null(default_markers) && all(markers %in% default_markers) && is.null(metacells) && has_mc_data(dataset, "default_markers_dist")) {
            cached_dist <- function() get_mc_data(dataset, "default_markers_dist")
        }
    } else if (mode == "Gene modules") {
        mc_fp <- get_mc_gene_modules_fp(dataset, markers, gene_modules)
        epsilon <- 0
        log_transform <- TRUE
    } else if (mode == "Outliers") {
        mc_fp <- get_mc_data(dataset, "deviant_fold_mat")
        req(mc_fp)
        mc_fp <- as.matrix(mc_fp[intersect(markers, rownames(mc_fp)[Matrix::rowSums(mc_fp) > 0]), , drop = FALSE])
        epsilon <- 1e-5
        log_transform <- FALSE
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
            mc_order <- order_mc_by_most_var_genes(mat, force_cell_type = force_cell_type, metacell_types = metacell_types, epsilon = epsilon, notify_var_genes = notify_var_genes, log_transform = log_transform, cached_dist = cached_dist, secondary_order = secondary_order, precomputed_clusterings = tryCatch(get_mc_data(dataset, "default_clusterings"), error = function(e) NULL)) # stored clusterings only save time: if they can't be read, cluster

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

    gene_ord <- order(apply(mat, 1, which.max))
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
    if (mode == "Inner") {
        if (!has_mc_data(dataset, "inner_fold_mat")) {
            if ("Inner-fold" %in% config$original_tabs) {
                showNotification(glue("Inner-fold matrix was not computed. Please compute it in python using the metacells package and rerun the import"), type = "error")
            }
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_inner_fold"))
    } else if (mode == "Stdev") {
        if (!has_mc_data(dataset, "inner_stdev_mat")) {
            if ("Stdev-fold" %in% config$original_tabs) {
                showNotification(glue("Inner-stdev matrix was not computed. Please compute it in python using the metacells package and rerun the import"), type = "error")
            }
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_inner_stdev"))
    } else if (mode == "Proj") {
        if (!has_mc_data(dataset, "projected_fold")) {
            if ("Projected-fold" %in% config$original_tabs) {
                showNotification(glue("Projected-fold matrix was not computed. Please compute it in python using the metacells package and rerun the import"), type = "error")
            }
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_projected"))
    } else if (mode == "Outliers") {
        if (!has_mc_data(dataset, "deviant_fold_mat")) {
            if ("Outliers" %in% config$original_tabs) {
                showNotification(glue("Outliers matrix was not imported. Please use outliers_anndata_file paramter to import it and re-run the app"), type = "error")
            }
            req(FALSE)
        }
        return(get_mc_data(dataset, "marker_genes_deviant_fold"))
    } else {
        return(get_markers(dataset))
    }

    return(get_markers(dataset))
}
