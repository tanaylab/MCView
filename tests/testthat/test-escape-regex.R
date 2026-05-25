# Regression test for escape_regex(), used by make_cell_type_view() to build
# literal cell-type match patterns for DAF `~` regex queries. The previous
# inline gsub was a no-op (it placed [, ], \ inside a bracket class), so type
# names with regex metacharacters matched the wrong cells.

test_that("escape_regex escapes metacharacters so patterns match literally", {
    types <- c("B.cell", "type(1)", "a+b", "x|y", "T*", "g$", "c^d", "q?", "k[1]", "m{2}")

    escaped <- escape_regex(types)

    # Each escaped pattern, anchored, matches its own literal...
    for (i in seq_along(types)) {
        pat <- paste0("^(", escaped[i], ")$")
        expect_true(grepl(pat, types[i]), info = types[i])
        # ...and does not match the literal plus a trailing character.
        expect_false(grepl(pat, paste0(types[i], "X")), info = types[i])
    }
})

test_that("escape_regex makes '.' literal (no wildcard match)", {
    pat <- paste0("^(", escape_regex("B.cell"), ")$")
    expect_true(grepl(pat, "B.cell"))
    expect_false(grepl(pat, "BXcell"))
})

test_that("escape_regex leaves non-metacharacters untouched", {
    expect_equal(escape_regex("NK-T"), "NK-T")
    expect_equal(escape_regex("CD4 T"), "CD4 T")
    expect_equal(escape_regex("plain"), "plain")
})
