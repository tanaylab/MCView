test_that("top2_per_col returns correct indices and values for dense matrix", {
    # Fill column-major so each column is [1,9,4], [5,3,8], [2,7,6]
    m <- matrix(
        c(
            1, 9, 4,
            5, 3, 8,
            2, 7, 6
        ),
        nrow = 3, byrow = FALSE
    )
    # columns: [1,9,4], [5,3,8], [2,7,6]
    # top-2 by value: col1 -> (9@2, 4@3), col2 -> (8@3, 5@1), col3 -> (7@2, 6@3)
    r <- top2_per_col(m)
    expect_equal(r$top1_idx, c(2L, 3L, 2L))
    expect_equal(r$top1_val, c(9, 8, 7))
    expect_equal(r$top2_idx, c(3L, 1L, 3L))
    expect_equal(r$top2_val, c(4, 5, 6))
})

test_that("top2_per_row returns correct indices and values for dense matrix", {
    m <- matrix(
        c(
            1, 9, 4,
            5, 3, 8,
            2, 7, 6
        ),
        nrow = 3, byrow = TRUE
    )
    # rows: [1,9,4], [5,3,8], [2,7,6]
    r <- top2_per_row(m)
    expect_equal(r$top1_idx, c(2L, 3L, 2L))
    expect_equal(r$top1_val, c(9, 8, 7))
    expect_equal(r$top2_idx, c(3L, 1L, 3L))
    expect_equal(r$top2_val, c(4, 5, 6))
})

test_that("top2_per_col returns NA when a column has fewer than 2 finite values", {
    m <- matrix(c(1, NA, NA, 2, 3, NA), nrow = 3)
    r <- top2_per_col(m)
    # col1: only 1 is finite -> top1=1@1, top2=NA
    # col2: 2 and 3 are finite -> top1=3@2, top2=2@1
    expect_equal(r$top1_idx, c(1L, 2L))
    expect_equal(r$top1_val, c(1, 3))
    expect_equal(r$top2_idx, c(NA_integer_, 1L))
    expect_equal(r$top2_val, c(NA_real_, 2))
})
