# Regression test: the masked fast path in daf_query_mc_mat() and
# compute_egc_from_daf() must return rows/cols in the REQUESTED order, like
# the full-matrix fallback does. DAF name-masks preserve the native axis
# order (they do NOT reorder to the mask's listing), so the two paths used to
# disagree. Positional consumers - e.g. `mc_egc %*% t(samp_mc_frac)` in
# get_samples_egc(), where samp_mc_frac is explicitly column-indexed to the
# requested metacell order - silently misalign when the masked path runs and
# the requested order differs from native order.

make_order_daf <- function() {
    genes <- c("gA", "gB", "gC", "gD", "gE")
    mcs <- c("m1", "m2", "m3", "m4", "m5")
    umis <- matrix(seq_len(length(genes) * length(mcs)),
        nrow = length(genes), ncol = length(mcs),
        dimnames = list(genes, mcs)
    )
    frac <- sweep(umis, 2, colSums(umis), "/")

    d <- dafr::memory_daf("order_fixture")
    dafr::add_axis(d, "gene", genes)
    dafr::add_axis(d, "metacell", mcs)
    dafr::set_matrix(d, "gene", "metacell", "UMIs", umis)
    dafr::set_matrix(d, "gene", "metacell", "linear_fraction", frac)
    dafr::set_vector(d, "metacell", "total_UMIs", colSums(umis))
    list(daf = d, umis = umis, frac = frac)
}

test_that("daf_query_mc_mat returns requested gene/metacell order (mask path)", {
    fx <- make_order_daf()
    req_genes <- c("gD", "gA", "gC") # scrambled vs native gA..gE
    req_mcs <- c("m4", "m1") # scrambled vs native m1..m5

    res <- daf_query_mc_mat(fx$daf, genes = req_genes, metacells = req_mcs)

    expect_equal(rownames(res), req_genes)
    expect_equal(colnames(res), req_mcs)
    # Values must match the same slice of the source matrix in requested order.
    expect_equal(as.matrix(res), fx$umis[req_genes, req_mcs, drop = FALSE])
})

test_that("compute_egc_from_daf returns requested metacell order (mask path)", {
    fx <- make_order_daf()
    req_mcs <- c("m5", "m2", "m1")

    res <- compute_egc_from_daf(fx$daf, metacells = req_mcs)

    expect_equal(colnames(res), req_mcs)
    expect_equal(as.matrix(res), fx$frac[, req_mcs, drop = FALSE])
})
