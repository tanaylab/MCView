test_that("MCView environment initialization works", {
    init_mcview_env()

    expect_true(mcv_exists("mc_data"))
    # mc_data is initially NULL, becomes a list when datasets are added
    expect_null(mcv_get("mc_data"))
})

test_that("Environment get/set operations work", {
    mcv_set("test_var", "test_value")
    expect_equal(mcv_get("test_var"), "test_value")
    expect_true(mcv_exists("test_var"))
})

# ---- validate_mcview_env() ---------------------------------------------------

# Helper: build a minimally-valid mcview_env state for the validator. Saves the
# real env, replaces it with the fixture, runs the body, then restores.
with_valid_env <- function(daf_obj, code) {
    saved <- as.list.environment(mcview_env, all.names = TRUE)
    on.exit({
        rm(list = ls(envir = mcview_env, all.names = TRUE), envir = mcview_env)
        list2env(saved, envir = mcview_env)
    })

    init_mcview_env()
    mcv_set("config", list(tabs = c("Manifold", "QC")))
    mcv_set("mc_data", list(test = list(base_daf = daf_obj, daf_obj = daf_obj)))
    mcv_set("tab_defs", list(qc = list(module_name = "qc")))
    mcv_set("egc_epsilon", 1e-5)
    mcv_set("expr_breaks", c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1))
    mcv_set("default_gene1", "GeneA")
    mcv_set("default_gene2", "GeneB")

    force(code)
}

test_that("validate_mcview_env passes on a fully populated env", {
    skip_if_no_daf()
    daf <- dafr::files_daf(get_test_daf_path(), mode = "r")
    on.exit(try(dafr::close_daf(daf), silent = TRUE), add = TRUE)

    with_valid_env(daf, {
        expect_true(validate_mcview_env())
    })
})

test_that("validate_mcview_env errors on missing required key", {
    skip_if_no_daf()
    daf <- dafr::files_daf(get_test_daf_path(), mode = "r")
    on.exit(try(dafr::close_daf(daf), silent = TRUE), add = TRUE)

    with_valid_env(daf, {
        mcv_set("egc_epsilon", NULL)
        expect_error(validate_mcview_env(), regexp = "egc_epsilon")
    })
})

test_that("validate_mcview_env errors on missing default_gene fields", {
    skip_if_no_daf()
    daf <- dafr::files_daf(get_test_daf_path(), mode = "r")
    on.exit(try(dafr::close_daf(daf), silent = TRUE), add = TRUE)

    with_valid_env(daf, {
        mcv_set("default_gene1", NULL)
        mcv_set("default_gene2", NULL)
        expect_error(validate_mcview_env(), regexp = "default_gene1")
    })
})

test_that("validate_mcview_env errors when mc_data is empty list", {
    skip_if_no_daf()
    daf <- dafr::files_daf(get_test_daf_path(), mode = "r")
    on.exit(try(dafr::close_daf(daf), silent = TRUE), add = TRUE)

    with_valid_env(daf, {
        mcv_set("mc_data", list())
        expect_error(validate_mcview_env(), regexp = "mc_data")
    })
})

test_that("validate_mcview_env errors when mc_data entry has unnamed slot", {
    skip_if_no_daf()
    daf <- dafr::files_daf(get_test_daf_path(), mode = "r")
    on.exit(try(dafr::close_daf(daf), silent = TRUE), add = TRUE)

    with_valid_env(daf, {
        mcv_set("mc_data", list(test = list(daf_obj = daf))) # missing base_daf
        expect_error(validate_mcview_env(), regexp = "base_daf")
    })
})

test_that("validate_mcview_env errors when mc_data slot is not a DAF", {
    skip_if_no_daf()
    daf <- dafr::files_daf(get_test_daf_path(), mode = "r")
    on.exit(try(dafr::close_daf(daf), silent = TRUE), add = TRUE)

    with_valid_env(daf, {
        mcv_set("mc_data", list(test = list(base_daf = "not a daf", daf_obj = daf)))
        expect_error(validate_mcview_env(), regexp = "not a DAF")
    })
})

test_that("validate_mcview_env accepts NULL atlas (optional)", {
    skip_if_no_daf()
    daf <- dafr::files_daf(get_test_daf_path(), mode = "r")
    on.exit(try(dafr::close_daf(daf), silent = TRUE), add = TRUE)

    with_valid_env(daf, {
        expect_null(mcv_get("atlas"))
        expect_true(validate_mcview_env())
    })
})

test_that("validate_mcview_env errors when atlas is set to non-DAF", {
    skip_if_no_daf()
    daf <- dafr::files_daf(get_test_daf_path(), mode = "r")
    on.exit(try(dafr::close_daf(daf), silent = TRUE), add = TRUE)

    with_valid_env(daf, {
        mcv_set("atlas", "not a daf either")
        expect_error(validate_mcview_env(), regexp = "atlas")
    })
})
