make_module_fixture <- function(prop = "module") {
    d <- dafr::memory_daf("module_fixture")
    dafr::add_axis(d, "gene", paste0("g", 1:6))
    dafr::add_axis(d, "metacell", paste0("m", 1:3))
    dafr::set_vector(d, "gene", prop, c("M1", "M1", "M2", "M2", "M1", "M2"))
    mat <- matrix(
        c(
            1, 2, 3, 4, 5, 6,
            7, 8, 9, 10, 11, 12,
            13, 14, 15, 16, 17, 18
        ),
        nrow = 6,
        dimnames = list(paste0("g", 1:6), paste0("m", 1:3))
    )
    dafr::set_matrix(d, "gene", "metacell", "UMIs", mat)
    dafr::set_vector(d, "metacell", "total_UMIs", colSums(mat))
    d
}

test_that("daf_query_module_fraction matches the legacy UMIs / mc_sum normalization", {
    d <- make_module_fixture()

    # Legacy path: two queries + double-transpose / divide.
    umis <- daf_query_module_umis(d)
    mc_sum <- as.numeric(dafr::get_vector(d, "metacell", "total_UMIs"))
    ref <- t(t(as.matrix(umis)) / mc_sum)

    # New path: one query with % Fraction.
    out <- daf_query_module_fraction(d)
    expect_false(is.null(out))
    expect_equal(dim(out), c(2L, 3L))
    expect_equal(as.matrix(out), ref, tolerance = 1e-12)
    # Per-metacell module fractions sum to 1 in this fixture (all genes
    # carry a module assignment, no orphan UMIs).
    expect_equal(colSums(as.matrix(out)), setNames(rep(1, 3), paste0("m", 1:3)),
                 tolerance = 1e-12)
})

test_that("daf_query_module_fraction prefers primary_module over module", {
    d <- make_module_fixture(prop = "primary_module")
    # Also set a legacy `module` with different labels so the wrong choice
    # would produce different output.
    dafr::set_vector(d, "gene", "module", c("X", "X", "X", "X", "X", "X"))

    out <- daf_query_module_fraction(d)
    expect_setequal(rownames(out), c("M1", "M2"))
})

test_that("daf_query_module_fraction returns NULL when no module property exists", {
    d <- dafr::memory_daf("no_module")
    dafr::add_axis(d, "gene", paste0("g", 1:2))
    dafr::add_axis(d, "metacell", paste0("m", 1:2))
    dafr::set_matrix(d, "gene", "metacell", "UMIs",
                     matrix(c(1, 2, 3, 4), nrow = 2))
    dafr::set_vector(d, "metacell", "total_UMIs", c(3, 7))

    expect_null(daf_query_module_fraction(d))
})

test_that("daf_query_module_fraction modules filter subsets rows", {
    d <- make_module_fixture()

    only_m1 <- daf_query_module_fraction(d, modules = "M1")
    expect_equal(rownames(only_m1), "M1")
    expect_equal(dim(only_m1), c(1L, 3L))

    # Filter with no intersection: matches existing daf_query_module_umis
    # semantic - returns unfiltered (defensive against typos).
    bogus <- daf_query_module_fraction(d, modules = "DOES_NOT_EXIST")
    expect_equal(rownames(bogus), c("M1", "M2"))
})

test_that("get_mc_gene_modules_egc uses daf_query_module_fraction when module prop exists", {
    d <- make_module_fixture()
    # Wire a minimal dataset entry through the MCView env so
    # get_daf_for_query("test_ds") resolves to our fixture.
    init_mcview_env()
    prev_mc_data <- mcv_get("mc_data")
    mcv_set("mc_data", list(test_ds = list(daf_obj = d)))
    on.exit(mcv_set("mc_data", prev_mc_data), add = TRUE)

    out <- get_mc_gene_modules_egc("test_ds")
    expect_false(is.null(out))
    # Same shape as the legacy path: module x metacell, per-metacell sums = 1.
    expect_equal(dim(out), c(2L, 3L))
    expect_equal(unname(colSums(as.matrix(out))), rep(1, 3), tolerance = 1e-12)
})
