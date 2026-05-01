test_that("mcview_palette has the expected named entries", {
    expect_true(is.list(mcview_palette))
    expect_named(mcview_palette,
        c("expression", "expression_with_zero", "cell_type_score", "edge"),
        ignore.order = TRUE)
})

test_that("mcview_palette$expression is the canonical 11-color diverging RdBu", {
    expect_length(mcview_palette$expression, 11)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", mcview_palette$expression)))
    expect_equal(mcview_palette$expression[1], "#053061")
    expect_equal(mcview_palette$expression[11], "#67001F")
    expect_equal(mcview_palette$expression[6], "#F7F7F7")
})

test_that("mcview_palette$cell_type_score matches the historical white->Reds->black ramp", {
    expect_equal(mcview_palette$cell_type_score,
        c("white", "#F4A582", "#D6604D", "#B2182B", "#67001F", "black"))
})

test_that("mcview_palette$expression_with_zero is white followed by viridis(6)", {
    expect_length(mcview_palette$expression_with_zero, 7)
    expect_equal(mcview_palette$expression_with_zero[1], "white")
    expect_equal(mcview_palette$expression_with_zero[-1], viridis::viridis_pal()(6))
})

test_that("mcview_palette$edge is a single near-black hex", {
    expect_length(mcview_palette$edge, 1)
    expect_match(mcview_palette$edge, "^#[0-9A-Fa-f]{6}$")
})
