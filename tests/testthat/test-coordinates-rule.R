# Regression test: the core-contract "coordinates" custom rule must accept the
# same coordinate names the converter (convert_daf_mc2d) resolves - x/y, u/v,
# AND umap_x/umap_y. It used to require (x,y) or (u,v) only, so a Metacells.jl
# DAF that stores coordinates as umap_x/umap_y failed validate_daf_for_mcview()
# and the dataset could not load, even though it renders fine.

coord_rule <- function() {
    mcview_core_contract()$custom_rules$coordinates$validate
}

make_coord_daf <- function(xname, yname) {
    d <- dafr::memory_daf("coord_fixture")
    dafr::add_axis(d, "metacell", paste0("m", 1:3))
    dafr::set_vector(d, "metacell", xname, c(0.1, 0.2, 0.3))
    dafr::set_vector(d, "metacell", yname, c(0.4, 0.5, 0.6))
    d
}

test_that("coordinates rule accepts x/y, u/v, and umap_x/umap_y", {
    rule <- coord_rule()
    expect_null(rule(make_coord_daf("x", "y")))
    expect_null(rule(make_coord_daf("u", "v")))
    expect_null(rule(make_coord_daf("umap_x", "umap_y")))
})

test_that("coordinates rule rejects a DAF with no coordinate pair", {
    rule <- coord_rule()
    d <- dafr::memory_daf("coord_none")
    dafr::add_axis(d, "metacell", paste0("m", 1:3))
    dafr::set_vector(d, "metacell", "type", rep("T", 3))
    expect_false(is.null(rule(d)))
})
