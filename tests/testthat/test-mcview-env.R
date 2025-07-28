test_that("MCView environment initialization works", {
    init_mcview_env()

    expect_true(mcv_exists("backend"))
    expect_null(mcv_get("backend"))
})

test_that("Environment get/set operations work", {
    mcv_set("test_var", "test_value")
    expect_equal(mcv_get("test_var"), "test_value")
    expect_true(mcv_exists("test_var"))
})

test_that("Backend system works", {
    # Test project backend
    set_backend("project", path = "/test/path")
    backend <- current_backend()
    expect_equal(backend$kind, "project")
    expect_equal(backend$path, "/test/path")

    # Test DAF backend (when implemented)
    mock_daf <- structure(list(), class = "Daf")
    set_backend("daf", daf_obj = mock_daf)
    backend <- current_backend()
    expect_equal(backend$kind, "daf")
})
