# Regression test: parse_cell_type_colors must reject a cell type mapped to
# more than one DISTINCT color. The original guard ran after
# distinct(cell_type), which had already collapsed the conflict away, so it
# never fired (a conflicting color was silently resolved to the first).

test_that("parse_cell_type_colors rejects conflicting colors for one cell type", {
    bad <- tibble::tibble(
        cell_type = c("A", "A", "B"),
        color = c("red", "blue", "green")
    )
    expect_error(parse_cell_type_colors(bad), "more than one color")
})

test_that("parse_cell_type_colors allows duplicate rows with the same color", {
    ok <- tibble::tibble(
        cell_type = c("A", "A", "B"),
        color = c("red", "red", "green")
    )
    res <- parse_cell_type_colors(ok)
    expect_equal(nrow(res), 2)
    expect_setequal(res$cell_type, c("A", "B"))
})

test_that("parse_cell_type_colors returns one row per cell type for clean input", {
    clean <- tibble::tibble(
        cell_type = c("A", "B", "C"),
        color = c("red", "green", "blue")
    )
    res <- parse_cell_type_colors(clean)
    expect_equal(nrow(res), 3)
    expect_true(all(c("cell_type", "color", "order") %in% colnames(res)))
})
