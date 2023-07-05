#' Compute a umap 2D projection based on gene anchors
#'
#' @noRd
compute_umap <- function(mc_egc, anchors, min_dist = 0.96, n_neighbors = 10, n_epoch = 500, min_log_expr = -14, config = NULL) {
    if (!all(anchors %in% rownames(mc_egc))) {
        cli::cli_abort("Umap gene{?s} {.val {anchors[!(anchors %in% rownames(mc_egc))]}} not found in metacell gene expression data")
    }

    legc <- log2(mc_egc + 1e-5)

    f_gcov <- matrixStats::rowMaxs(legc) >= min_log_expr
    legc <- legc[f_gcov, ]

    if (!all(anchors %in% rownames(legc))) {
        cli::cli_warn("Umap genes: {.val {anchors[!(anchors %in% rownames(legc))]}} do not have a minimum expression of {.val {min_log_expr}} in the metacell gene expression data. Please consider using a lower value for {.val {min_log_expr}}")
    }
    anchors <- anchors[anchors %in% rownames(legc)]

    if (length(anchors) < 2) {
        cli::cli_warn("At least 2 anchors are required to compute a umap projection")
        return(NULL)
    }

    knn <- tgs_cor_knn(t(legc[anchors, ]), t(legc), knn = 30)
    gmods <- knn$col1
    names(gmods) <- knn$col2
    legc_anchors <- legc[names(gmods), ]

    feats <- t(tgs_matrix_tapply(t(legc_anchors), gmods, mean, na.rm = TRUE))

    if (is.null(config)) {
        config <- umap::umap.defaults
        config$min_dist <- min_dist
        config$n_neighbors <- n_neighbors
        config$n_epoch <- n_epoch
    }


    uf <- umap::umap(feats, config = config)

    # graph
    idx_df <- uf$knn$indexes %>%
        as.data.frame() %>%
        rownames_to_column("mc1") %>%
        gather("k", "idx", -mc1) %>%
        as_tibble()

    distances_df <- uf$knn$distances %>%
        as.data.frame() %>%
        rownames_to_column("mc1") %>%
        gather("k", "weight", -mc1) %>%
        as_tibble()

    graph <- idx_df %>%
        left_join(distances_df, by = join_by(mc1, k)) %>%
        mutate(mc2 = rownames(feats)[idx]) %>%
        filter(mc1 != mc2) %>%
        select(mc1, mc2, weight)

    layout_df <- as.data.frame(uf$layout) %>%
        rownames_to_column("mc") %>%
        as_tibble()

    colnames(layout_df) <- c("mc", "umap_x", "umap_y")

    mc2d_list <- list(
        graph = graph,
        mc_id = layout_df$mc,
        mc_x = tibble::deframe(layout_df[, c("mc", "umap_x")]),
        mc_y = tibble::deframe(layout_df[, c("mc", "umap_y")])
    )

    return(mc2d_list)
}
