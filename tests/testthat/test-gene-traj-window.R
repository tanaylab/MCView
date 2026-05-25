# Regression test for the non-zero time-window selection in calc_mc_egc_t().
#
# The original code used `which.max(c(1:ncol)[non_zero])` /
# `which.min(...)`, which return the POSITION WITHIN the subset of
# non-zero columns rather than the original column index. That collapses
# to "take the first sum(non_zero) columns", which is only correct when
# the non-zero flow window starts at column 1 and is contiguous. This
# test uses an offset non-zero window (cols 3:5 of 7) to lock the fix.

test_that("calc_mc_egc_t selects the correct offset non-zero time window", {
    # flow_prob: metacells (rows) x time bins (cols). Non-zero only in cols 3:5.
    flow_prob <- matrix(0, nrow = 3, ncol = 7)
    flow_prob[, 3] <- c(0.2, 0.3, 0.5)
    flow_prob[, 4] <- c(0.1, 0.6, 0.3)
    flow_prob[, 5] <- c(0.4, 0.4, 0.2)
    rownames(flow_prob) <- paste0("mc", 1:3)

    mc_egc <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
    rownames(mc_egc) <- c("geneA", "geneB")
    colnames(mc_egc) <- paste0("mc", 1:3)

    testthat::local_mocked_bindings(
        get_mc_data = function(dataset, var_name, atlas = FALSE) {
            if (var_name == "mct_probs_trans") {
                list(mc1 = list(probs = flow_prob))
            } else {
                NULL
            }
        },
        get_mc_egc = function(dataset, ...) mc_egc
    )

    res <- calc_mc_egc_t("dummy", "mc1")

    # Correct window is cols 3:5, each time bin column normalized to sum 1.
    win <- flow_prob[, 3:5, drop = FALSE]
    win_n <- t(t(win) / colSums(win))
    expected <- mc_egc %*% win_n

    expect_equal(dim(res), c(2, 3))
    expect_false(any(is.nan(res)))
    expect_equal(unname(res), unname(expected))
})
