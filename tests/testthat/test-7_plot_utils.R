test_that("thin_log_breaks keeps at most max_n breaks, at round values", {
    # log2 labels: every log2 unit, thinned to multiples of a step
    expect_equal(log2(thin_log_breaks(2^(-17:0))), c(-15, -12, -9, -6, -3, 0))
    expect_equal(log2(thin_log_breaks(2^(-17:-4))), c(-15, -12, -9, -6))
    expect_equal(thin_log_breaks(2^(-8:-4)), 2^(-8:-4))
    # 1, 2, 4 x 10^k: the powers of 10
    breaks <- c(1e-5, 2e-5, 4e-5, 1e-4, 2e-4, 4e-4, 1e-3, 2e-3, 4e-3, 1e-2, 2e-2, 4e-2, 1e-1, 2e-1, 4e-1, 1)
    expect_equal(thin_log_breaks(breaks), 10^(-5:0))
    expect_equal(thin_log_breaks(breaks[2:10]), c(1e-4, 1e-3, 1e-2))
    expect_equal(thin_log_breaks(breaks[4:8]), breaks[4:8])
    expect_equal(thin_log_breaks(breaks, max_n = 3), c(1e-4, 1e-2, 1))
    for (n in 2:16) {
        expect_lte(length(thin_log_breaks(breaks[1:n])), 6)
        expect_gte(length(thin_log_breaks(breaks[1:n])), min(n, 2))
    }
})

test_that("scatter_legend_categories counts the legend entries of the metacells shown", {
    mt <- tibble::tibble(metacell = paste0("m", 1:30), cell_type = rep(paste0("type", 1:15), 2))
    md <- tibble::tibble(metacell = mt$metacell, batch = rep(c("a", "b", "c"), 10), score = seq(0, 1, length.out = 30), Clipboard = "not selected")
    expect_equal(scatter_legend_categories("Metadata", "Cell type", mt, md), 15)
    expect_equal(scatter_legend_categories("Metadata", "Cell type", mt, md, metacell_filter = paste0("m", 1:4)), 4)
    expect_equal(scatter_legend_categories("Metadata", "batch", mt, md), 3)
    expect_equal(scatter_legend_categories("Metadata", "batch", mt, md, metacell_filter = c("m1", "m4")), 1)
    expect_equal(scatter_legend_categories("Metadata", "Clipboard", mt, md), 1)
    expect_equal(scatter_legend_categories("Metadata", "score", mt, md), 0)
    expect_equal(scatter_legend_categories("Gene", "Actb", mt, md), 0)
})
