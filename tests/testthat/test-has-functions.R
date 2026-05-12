# Parametric tests for the has_*() family.
#
# Architecture review section 8.3 flagged that the has_time() bug shipped because
# the sole real test DAF had no time annotations, so the FALSE-return regression
# went undetected. These tests build small in-memory DAFs that carry exactly the
# feature under test (or omit it) and pin both branches.

build_minimal_daf <- function(features = character(0)) {
    d <- dafr::memory_daf("test")
    dafr::add_axis(d, "metacell", c("mc1", "mc2"))
    dafr::add_axis(d, "gene", c("geneA", "geneB"))
    dafr::add_axis(d, "type", "typeA")
    dafr::set_vector(d, "metacell", "type", c("typeA", "typeA"))

    if ("metacell_graph" %in% features) {
        dafr::add_axis(d, "metacell_graph", "e1")
    }
    if ("time" %in% features) {
        dafr::set_vector(d, "metacell", "time", c(0.1, 0.5))
    }
    if ("projected_type" %in% features) {
        dafr::set_vector(d, "metacell", "projected_type", c("typeA", "typeA"))
    }
    if ("corrected_fraction" %in% features) {
        m <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2)
        rownames(m) <- c("mc1", "mc2")
        colnames(m) <- c("geneA", "geneB")
        dafr::set_matrix(d, "metacell", "gene", "corrected_fraction", m, relayout = TRUE)
    }
    if ("cell" %in% features) {
        dafr::add_axis(d, "cell", c("c1", "c2", "c3"))
    }
    if ("samp_id" %in% features) {
        if (!("cell" %in% features)) {
            dafr::add_axis(d, "cell", c("c1", "c2", "c3"))
        }
        dafr::set_vector(d, "cell", "samp_id", c("s1", "s1", "s2"))
    }

    d
}

# Register a freshly-built DAF directly into mcview_env, bypassing the contract
# validator (which requires UMI matrices etc. these unit tests don't need).
# Saves the prior env on entry and restores it on test exit so test ordering
# stays independent.
with_test_daf <- function(features, name = "test_fixture", code) {
    saved <- as.list.environment(mcview_env, all.names = TRUE)
    on.exit({
        rm(list = ls(envir = mcview_env, all.names = TRUE), envir = mcview_env)
        list2env(saved, envir = mcview_env)
    })

    init_mcview_env()
    daf <- build_minimal_daf(features)
    mc_data <- list()
    mc_data[[name]] <- list(
        base_daf = daf,
        daf_obj = daf,
        top_cor_genes = list()
    )
    mcv_set("mc_data", mc_data)

    force(code)
}

# ---- Individual has_*() functions --------------------------------------------

test_that("has_time returns TRUE iff metacell.time vector exists", {
    skip_if_no_daf()
    with_test_daf(features = character(0), code = {
        expect_false(has_time("test_fixture"))
    })
    with_test_daf(features = "time", code = {
        expect_true(has_time("test_fixture"))
    })
})

test_that("has_network returns TRUE iff metacell_graph axis exists", {
    skip_if_no_daf()
    with_test_daf(features = character(0), code = {
        expect_false(has_network("test_fixture"))
    })
    with_test_daf(features = "metacell_graph", code = {
        expect_true(has_network("test_fixture"))
    })
})

test_that("has_projection returns TRUE iff metacell.projected_type vector exists", {
    skip_if_no_daf()
    with_test_daf(features = character(0), code = {
        expect_false(has_projection("test_fixture"))
    })
    with_test_daf(features = "projected_type", code = {
        expect_true(has_projection("test_fixture"))
    })
})

test_that("has_corrected returns TRUE iff metacell.gene corrected_fraction matrix exists", {
    skip_if_no_daf()
    with_test_daf(features = character(0), code = {
        expect_false(has_corrected("test_fixture"))
    })
    with_test_daf(features = "corrected_fraction", code = {
        expect_true(has_corrected("test_fixture"))
    })
})

test_that("has_cell_metadata returns TRUE iff cell axis exists", {
    skip_if_no_daf()
    with_test_daf(features = character(0), code = {
        expect_false(has_cell_metadata("test_fixture"))
    })
    with_test_daf(features = "cell", code = {
        expect_true(has_cell_metadata("test_fixture"))
    })
})

test_that("has_atlas returns TRUE iff an atlas DAF is set in mcview_env", {
    skip_if_no_daf()
    with_test_daf(features = character(0), code = {
        # No atlas registered for the fixture
        expect_false(has_atlas("test_fixture"))

        # Register an atlas and re-check
        mcv_set("atlas", build_minimal_daf())
        expect_true(has_atlas("test_fixture"))
    })
})

# ---- NULL / missing-dataset path ---------------------------------------------

test_that("has_*() return FALSE for an unknown dataset name", {
    skip_if_no_daf()
    with_test_daf(features = c("time", "metacell_graph", "projected_type", "corrected_fraction", "cell"), code = {
        # All features present for the real fixture, but a different name => NULL daf
        expect_false(has_time("nope"))
        expect_false(has_network("nope"))
        expect_false(has_projection("nope"))
        expect_false(has_corrected("nope"))
        # has_cell_metadata falls through to a get_mc_data check on NULL → FALSE
        expect_false(has_cell_metadata("nope"))
    })
})
