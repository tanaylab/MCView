test_that("theme_mcview returns a ggplot2 theme object", {
    th <- theme_mcview()
    expect_s3_class(th, "theme")
    expect_s3_class(th, "gg")
})

test_that("theme_mcview accepts base_size override", {
    expect_equal(theme_mcview(base_size = 12)$text$size, 12)
    expect_equal(theme_mcview(base_size = 18)$text$size, 18)
})

test_that("theme_mcview composes with ggplot without warnings", {
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, hp)) +
        ggplot2::geom_point() +
        theme_mcview()
    expect_silent(invisible(ggplot2::ggplot_build(p)))
})

test_that("theme_mcview overrides survive subsequent + theme(...) layering", {
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, hp)) +
        ggplot2::geom_point() +
        theme_mcview() +
        ggplot2::theme(legend.position = "top")
    th <- p$theme
    expect_equal(th$legend.position, "top")
    # base theme retained: panel.grid.major.x should still be element_blank
    expect_s3_class(th$panel.grid.major.x, "element_blank")
})
