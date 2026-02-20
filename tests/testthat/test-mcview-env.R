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
