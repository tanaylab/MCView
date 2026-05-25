# Regression test: module-mode gene correlation must run in LOG space (module
# log-score vs gene lfp), matching single-gene correlation (lfp vs lfp). It
# used to correlate a LINEAR module fraction against log gene expression.

setup_module_corr_dataset <- function(n_genes = 40L, n_mcs = 60L, seed = 17L) {
    set.seed(seed)
    egc <- matrix(abs(rnorm(n_genes * n_mcs, mean = 1, sd = 0.5)),
        nrow = n_genes,
        dimnames = list(paste0("g", seq_len(n_genes)), paste0("m", seq_len(n_mcs)))
    )
    egc <- sweep(egc, 2, colSums(egc), "/")
    d <- dafr::memory_daf("module_corr_fixture")
    dafr::add_axis(d, "gene", rownames(egc))
    dafr::add_axis(d, "metacell", colnames(egc))
    dafr::set_matrix(d, "gene", "metacell", "linear_fraction", egc)
    dafr::set_vector(d, "metacell", "total_UMIs", rep(1.0, n_mcs))
    dafr::set_matrix(d, "gene", "metacell", "UMIs", round(egc * 1000))
    init_mcview_env()
    prev_mc_data <- mcv_get("mc_data")
    prev_eps <- mcv_get("egc_epsilon")
    mcv_set("egc_epsilon", 1e-5)
    mcv_set("mc_data", list(test_ds = list(daf_obj = d)))
    list(restore = function() {
        mcv_set("mc_data", prev_mc_data)
        mcv_set("egc_epsilon", prev_eps)
    })
}

test_that("calc_module_correlations correlates the log module score against gene lfp", {
    ctx <- setup_module_corr_dataset()
    on.exit(ctx$restore(), add = TRUE)

    lfp <- get_mc_lfp("test_ds")
    module_genes <- rownames(lfp)[1:3]

    res <- calc_module_correlations("test_ds", module_genes, n_top = 10, threshold = 0)

    # Reference correlation VALUES: log2(mean module fraction + eps) correlated
    # against each gene's lfp (log-vs-log). The old linear code correlated
    # colMeans(egc) (linear) against lfp, giving different cor values - so this
    # value comparison (not just ranking) discriminates the fix.
    egc <- get_mc_egc("test_ds", genes = module_genes)
    mod_score <- log2(colMeans(egc) + 1e-5)
    ref <- apply(lfp, 1, function(r) stats::cor(r, mod_score))

    res_cor <- stats::setNames(res$cor, res$gene2)
    common <- intersect(names(res_cor), names(ref))
    expect_gt(length(common), 0)
    expect_equal(res_cor[common], ref[common], tolerance = 1e-6)
})
