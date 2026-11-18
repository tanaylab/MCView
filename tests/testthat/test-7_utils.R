test_that("subset_dist matches subsetting the full distance matrix", {
    x <- matrix(rnorm(60), nrow = 20, dimnames = list(paste0("m", 1:20), NULL))
    d <- dist(x)
    expect_identical(subset_dist(d, labels(d)), d)
    for (labs in list(paste0("m", c(2, 5, 6, 11, 20)), paste0("m", c(7, 3, 9)), paste0("m", 1:19))) {
        expect_equal(as.matrix(subset_dist(d, labs)), as.matrix(d)[labs, labs])
    }
})

test_that("cor_dist matches tgs_dist(tgs_cor())", {
    set.seed(60427)
    feat <- matrix(rnorm(30 * 700), nrow = 30, dimnames = list(NULL, paste0("m", 1:700)))
    expected <- tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE))
    d <- cor_dist(feat)
    expect_equal(as.vector(d), as.vector(expected))
    expect_identical(attr(d, "Labels"), attr(expected, "Labels"))
    expect_identical(fastcluster::hclust(d, method = "ward.D2")$merge, fastcluster::hclust(expected, method = "ward.D2")$merge)

    feat[, 5] <- 1 # constant metacell: falls back to tgs_dist(tgs_cor())
    expect_identical(cor_dist(feat), tgs_dist(tgs_cor(feat, pairwise.complete.obs = TRUE)))
})

test_that("cluster_metacells returns the same clustering with and without the app cache", {
    feat <- matrix(rnorm(200), nrow = 10, dimnames = list(NULL, paste0("m", 1:20)))
    hc <- cluster_metacells(feat)
    shiny::shinyOptions(cache = cachem::cache_mem())
    withr::defer(shiny::shinyOptions(cache = NULL))
    expect_identical(cluster_metacells(feat)$order, hc$order)
    expect_identical(cluster_metacells(feat)$order, hc$order) # cache hit
})

test_that("calc_diff_expr p-values match chisq.test", {
    set.seed(60427)
    # counts and projected (non-integer) counts, with many genes far from the diagonal
    mat <- cbind(a = rpois(3000, rexp(3000, 0.01)), b = rexp(3000, 0.01) * runif(3000, 0, 5))
    rownames(mat) <- paste0("g", 1:3000)
    egc <- t(t(mat) / colSums(mat)) + 1e-5
    df <- calc_diff_expr(mat, egc, c("a", "b"))
    f <- abs(df$diff) >= 1.5
    expect_gt(sum(f), 500)
    m <- mat[df$gene[f], ]
    tots <- colSums(m)
    expected <- apply(m, 1, function(x) suppressWarnings(stats::chisq.test(matrix(c(x, tots), nrow = 2))$p.value))
    expect_identical(df$pval[f], unname(expected))
    expect_true(all(is.na(df$pval[!f])))

    # a single gene far from the diagonal
    mat1 <- rbind(g1 = c(a = 500, b = 3), g2 = c(a = 5000, b = 300))
    df1 <- calc_diff_expr(mat1, t(t(mat1) / colSums(mat1)) + 1e-5, c("a", "b"))
    expect_identical(df1$pval[df1$gene == "g1"], suppressWarnings(stats::chisq.test(matrix(c(500, 3, 500, 3), nrow = 2))$p.value))
})
