test_that("subset_dist matches subsetting the full distance matrix", {
    x <- matrix(rnorm(60), nrow = 20, dimnames = list(paste0("m", 1:20), NULL))
    d <- dist(x)
    expect_identical(subset_dist(d, labels(d)), d)
    for (labs in list(paste0("m", c(2, 5, 6, 11, 20)), paste0("m", c(7, 3, 9)), paste0("m", 1:19))) {
        expect_equal(as.matrix(subset_dist(d, labs)), as.matrix(d)[labs, labs])
    }
})

test_that("cluster_metacells returns the same clustering with and without the app cache", {
    feat <- matrix(rnorm(200), nrow = 10, dimnames = list(NULL, paste0("m", 1:20)))
    hc <- cluster_metacells(feat)
    shiny::shinyOptions(cache = cachem::cache_mem())
    withr::defer(shiny::shinyOptions(cache = NULL))
    expect_identical(cluster_metacells(feat)$order, hc$order)
    expect_identical(cluster_metacells(feat)$order, hc$order) # cache hit
})
