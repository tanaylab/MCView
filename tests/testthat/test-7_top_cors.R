test_that("calc_top_cors matches stats::cor on log2 egc", {
    set.seed(60427)
    mat <- matrix(rpois(300 * 40, 2), 300, 40, dimnames = list(paste0("g", 1:300), paste0("mc", 1:40)))
    mat[7, ] <- 0 # constant gene, should be dropped (NA correlation)

    prev <- mget(c("mc_data", "egc_epsilon"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    assign("mc_data", list(toy = list(mc_mat = mat, mc_sum = colSums(mat))), envir = globalenv())
    assign("egc_epsilon", 1e-5, envir = globalenv())

    lfp <- log2(t(t(mat) / colSums(mat)) + 1e-5)
    top_both <- function(cors) c(head(sort(cors, decreasing = TRUE), n = 30), head(sort(cors), n = 30))

    # gene vs. genes on a subset of metacells
    mcs <- paste0("mc", 1:25)
    expected <- top_both(suppressWarnings(cor(lfp["g1", mcs], t(lfp[, mcs])))[1, ])
    res <- calc_top_cors("toy", "g1", "both", NULL, mcs, NULL)
    expect_equal(res$gene2, names(expected))
    expect_equal(res$cor, unname(expected))

    # data vector with missing values, excluding genes
    vec <- setNames(rnorm(40), colnames(mat))
    vec[c(3, 9)] <- NA
    ok <- !is.na(vec)
    genes <- setdiff(rownames(mat), c("g2", "g3"))
    expected <- top_both(suppressWarnings(cor(vec[ok], t(lfp[genes, ok])))[1, ])
    res <- calc_top_cors("toy", "x", "both", vec, NULL, c("g2", "g3"))
    expect_equal(res$gene2, names(expected))
    expect_equal(res$cor, unname(expected))

    # several genes at once
    genes <- c("g4", "g1", "g10")
    expected <- suppressWarnings(cor(t(lfp[, mcs]), t(lfp[genes, mcs])))
    expect_equal(cor_log_egc("toy", t(lfp[genes, mcs]), rownames(mat), mcs), expected)
})

test_that("calc_individual_correlations skips a gene without counts", {
    set.seed(60427)
    mat <- matrix(rpois(300 * 40, 2), 300, 40, dimnames = list(paste0("g", 1:300), paste0("mc", 1:40)))
    mat[7, ] <- 0 # no correlations

    prev <- mget(c("mc_data", "egc_epsilon"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    assign("mc_data", list(toy = list(mc_mat = mat, mc_sum = colSums(mat))), envir = globalenv())
    assign("egc_epsilon", 1e-5, envir = globalenv())

    expect_identical(select_top_cors(c(g1 = NaN, g2 = NaN), "both"), tibble(gene2 = character(0), cor = numeric(0)))
    res <- calc_individual_correlations("toy", c("g7", "g1"))
    expect_identical(res, calc_individual_correlations("toy", "g1"))
    expect_gt(nrow(res), 0)
})
