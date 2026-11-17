test_that("subset_dist matches subsetting the full distance matrix", {
    x <- matrix(rnorm(60), nrow = 20, dimnames = list(paste0("m", 1:20), NULL))
    d <- dist(x)
    expect_identical(subset_dist(d, labels(d)), d)
    for (labs in list(paste0("m", c(2, 5, 6, 11, 20)), paste0("m", c(7, 3, 9)), paste0("m", 1:19))) {
        expect_equal(as.matrix(subset_dist(d, labs)), as.matrix(d)[labs, labs])
    }
})
