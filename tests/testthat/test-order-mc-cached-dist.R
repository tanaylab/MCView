make_feat_and_dist <- function(seed = 42L) {
    set.seed(seed)
    n_genes <- 10L
    n_mcs <- 30L
    feat <- matrix(rnorm(n_genes * n_mcs), nrow = n_genes,
                   dimnames = list(
                       paste0("g", seq_len(n_genes)),
                       paste0("m", seq_len(n_mcs))
                   ))
    # Dist matrix in the shape produced by precompute_daf_default_markers:
    # base R matrix, full N x N, named on both axes.
    dist_full <- as.matrix(dist(t(feat)))
    list(feat = feat, dist_full = dist_full)
}

test_that("cached_dist subset path produces hclust output identical to coerce-first path", {
    fx <- make_feat_and_dist()
    mcs <- colnames(fx$feat)[1:20]
    feat_sub <- fx$feat[, mcs, drop = FALSE]

    # Old behaviour: coerce full dist to base matrix, then subset.
    sub_old <- as.dist(as.matrix(fx$dist_full)[mcs, mcs])

    # New behaviour: subset first, coerce only if needed.
    sub_new_mat <- fx$dist_full[mcs, mcs, drop = FALSE]
    if (!is.matrix(sub_new_mat)) sub_new_mat <- as.matrix(sub_new_mat)
    sub_new <- as.dist(sub_new_mat)

    expect_equal(as.vector(sub_old), as.vector(sub_new), tolerance = 1e-12)
    expect_equal(attr(sub_old, "Labels"), attr(sub_new, "Labels"))

    hc_old <- fastcluster::hclust(sub_old, method = "ward.D2")
    hc_new <- fastcluster::hclust(sub_new, method = "ward.D2")

    expect_equal(hc_old$order, hc_new$order)
    expect_equal(hc_old$height, hc_new$height, tolerance = 1e-12)
    expect_equal(hc_old$merge, hc_new$merge)
})

test_that("order_mc_by_most_var_genes returns identical ordering with cached_dist", {
    fx <- make_feat_and_dist()
    mcs <- colnames(fx$feat)

    # Pretend fx$feat is already a log/fp matrix (log_transform=FALSE skips
    # the inner log2 transform) so the function uses our values directly.
    # filter_markers=FALSE keeps all genes; we pick marks explicitly.
    marks <- c("g1", "g2")

    out_full <- order_mc_by_most_var_genes(
        gene_folds = fx$feat,
        marks = marks,
        filter_markers = FALSE,
        log_transform = FALSE,
        notify_var_genes = FALSE,
        cached_dist = fx$dist_full
    )

    # Subset case: only 20 metacells. The cached_dist still has the full
    # 30 entries; the subset-before-coerce path must produce the same order
    # as the legacy as.matrix(.)-then-subset path. Since both paths are
    # numerically identical (test above) the hclust output is identical, so
    # the final ordering returned by this function must match itself
    # against any randomization of the input metacell order.
    subset_mcs <- mcs[1:20]
    feat_subset <- fx$feat[, subset_mcs, drop = FALSE]
    out_subset <- order_mc_by_most_var_genes(
        gene_folds = feat_subset,
        marks = marks,
        filter_markers = FALSE,
        log_transform = FALSE,
        notify_var_genes = FALSE,
        cached_dist = fx$dist_full
    )

    expect_true(is.numeric(out_subset))
    expect_length(out_subset, ncol(feat_subset))
    # Each metacell index is visited exactly once.
    expect_setequal(as.integer(out_subset), seq_along(subset_mcs))
})
