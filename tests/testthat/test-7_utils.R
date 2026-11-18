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
    hc <- cluster_metacells(feat, "key")
    shiny::shinyOptions(cache = cachem::cache_mem())
    withr::defer(shiny::shinyOptions(cache = NULL))
    expect_identical(cluster_metacells(feat, "key")$order, hc$order)
    expect_identical(cluster_metacells(feat, "key")$order, hc$order) # cache hit
})

test_that("cluster_metacells uses a stored clustering without reading the distances", {
    feat <- matrix(rnorm(200), nrow = 10, dimnames = list(NULL, paste0("m", 1:20)))
    stored <- list(key = fastcluster::hclust(dist(t(feat))))
    expect_identical(cluster_metacells(feat, "key", cached_dist = function() stop("read the distances"), precomputed_clusterings = stored), stored$key)
    d <- dist(t(feat))
    expect_identical(cluster_metacells(feat, "other", cached_dist = function() d, precomputed_clusterings = stored)$merge, fastcluster::hclust(d, method = "ward.D2")$merge)
})
