test_that("calc_top_cors keeps its results with the data of the dataset", {
    set.seed(60427)
    mat <- matrix(rpois(300 * 40, 2), 300, 40, dimnames = list(paste0("g", 1:300), paste0("mc", 1:40)))

    prev <- mget(c("mc_data", "egc_epsilon"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    assign("mc_data", list(toy = list2env(list(mc_mat = mat, mc_sum = colSums(mat)))), envir = globalenv()) # an environment per dataset, as in the app
    assign("egc_epsilon", 1e-5, envir = globalenv())

    mcs <- paste0("mc", 1:25)
    vec <- setNames(rnorm(40), colnames(mat))
    direct <- calc_top_cors_uncached("toy", "g1", "both", NULL, mcs, NULL)
    expect_identical(calc_top_cors("toy", "g1", "both", NULL, mcs, NULL), direct)
    expect_identical(top_cors_cache("toy")$size(), 1L)
    # from the cache
    expect_identical(calc_top_cors("toy", "g1", "both", NULL, mcs, NULL), direct)
    expect_identical(top_cors_cache("toy")$size(), 1L)

    direct_vec <- calc_top_cors_uncached("toy", "x", "pos", vec, NULL, c("g2", "g3"))
    expect_identical(calc_top_cors("toy", "x", "pos", vec, NULL, c("g2", "g3")), direct_vec)
    expect_identical(calc_top_cors("toy", "x", "pos", vec, NULL, c("g2", "g3")), direct_vec)
    expect_identical(top_cors_cache("toy")$size(), 2L)

    # every input the result depends on changes the key
    key <- function(...) {
        args <- utils::modifyList(list(dataset = "toy", gene = "g1", type = "both", data_vec = NULL, metacell_filter = mcs, exclude = NULL, atlas = FALSE, n = 30), list(...), keep.null = TRUE)
        do.call(top_cors_key, args)
    }
    keys <- c(
        key(), key(dataset = "other"), key(gene = "g2"), key(type = "pos"), key(data_vec = vec),
        key(data_vec = vec * 2), key(metacell_filter = mcs[-1]), key(metacell_filter = NULL), key(exclude = "g5"),
        key(atlas = TRUE), key(n = 10)
    )
    expect_identical(anyDuplicated(keys), 0L)
    expect_identical(key(), key())
    # the gene isn't used with a data vector
    expect_identical(key(data_vec = vec, gene = "a"), key(data_vec = vec, gene = "b"))
    assign("egc_epsilon", 1e-4, envir = globalenv())
    expect_false(identical(key(), keys[1]))
    assign("egc_epsilon", 1e-5, envir = globalenv())

    # new data of the dataset start a new cache
    assign("mc_data", list(toy = list2env(list(mc_mat = mat[, 40:1], mc_sum = colSums(mat)[40:1]))), envir = globalenv())
    expect_identical(top_cors_cache("toy")$size(), 0L)
    # without an environment per dataset nothing is kept
    assign("mc_data", list(toy = list(mc_mat = mat, mc_sum = colSums(mat))), envir = globalenv())
    expect_null(top_cors_cache("toy"))
    expect_identical(calc_top_cors("toy", "g1", "both", NULL, mcs, NULL), direct)
})

test_that("natural_sort_cell_types orders the numbers in the names as numbers", {
    expect_identical(natural_sort_cell_types(c("10", "2", "1", "20", "3", "11")), c("1", "2", "3", "10", "11", "20"))
    expect_identical(natural_sort_cell_types(c("C10", "C2", "C1")), c("C1", "C2", "C10"))
    expect_identical(natural_sort_cell_types(c("T cells 10", "B cells", "T cells 9")), c("B cells", "T cells 9", "T cells 10"))
    expect_identical(natural_sort_cell_types(c("12345678901234567890", "9")), c("9", "12345678901234567890"))
    # names without numbers as sort() orders them, missing names last
    plain <- c("NK", "Mono", "B", "T", "DC")
    expect_identical(natural_sort_cell_types(plain), sort(plain))
    expect_identical(natural_sort_cell_types(c("2", NA, "1")), c("1", "2", NA))
    expect_identical(natural_sort_cell_types(character(0)), character(0))
    expect_identical(natural_order_cell_types(factor(c("10", "9"))), c(2L, 1L))

    # the order of sort() becomes natural, any other order is kept; "(Missing)" keeps its place
    expect_identical(cell_types_display_order(c("(Missing)", "1", "10", "2")), c(1L, 2L, 4L, 3L))
    expect_identical(cell_types_display_order(c("1", "10", "2", "(Missing)")), c(1L, 3L, 2L, 4L))
    expect_identical(cell_types_display_order(c("10", "1", "2")), 1:3)
    expect_identical(cell_types_display_order(c("Mono", "B", "T")), 1:3)
    expect_identical(cell_types_display_order(c("B", "Mono", "T")), 1:3)
    expect_identical(cell_types_display_order("1"), 1L)
})
