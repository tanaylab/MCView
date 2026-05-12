test_that("mcview_notify routes each severity to the matching cli channel", {
    out_info <- testthat::capture_messages(mcview_notify("info", "info msg", where = "console"))
    expect_match(paste(out_info, collapse = ""), "info msg")

    out_success <- testthat::capture_messages(mcview_notify("success", "success msg", where = "console"))
    expect_match(paste(out_success, collapse = ""), "success msg")

    out_warning <- testthat::capture_messages(mcview_notify("warning", "warn msg", where = "console"))
    expect_match(paste(out_warning, collapse = ""), "warn msg")

    out_error <- testthat::capture_messages(mcview_notify("error", "err msg", where = "console"))
    expect_match(paste(out_error, collapse = ""), "err msg")
})

test_that("mcview_notify with where='ui' skips the console emission outside a Shiny session", {
    # No reactive domain active here, so where='ui' should produce no output and not error.
    out <- testthat::capture_messages(mcview_notify("info", "should not appear on console", where = "ui"))
    expect_equal(out, character(0))
})

test_that("mcview_notify with where='console' emits exactly once on the cli channel", {
    out <- testthat::capture_messages(mcview_notify("warning", "console only", where = "console"))
    expect_equal(sum(grepl("console only", out)), 1)
})

test_that("mcview_notify rejects unknown severity", {
    expect_error(mcview_notify("fatal", "boom"), regexp = "should be one of")
})

test_that("mcview_notify accepts cli inline markup in message", {
    out <- testthat::capture_messages(mcview_notify("info", "loaded {.val foo}", where = "console"))
    expect_match(paste(out, collapse = ""), "foo")
})

test_that("mcview_notify is a no-op for UI when no Shiny session is active", {
    # We can't directly observe a non-existent toast, but we can confirm
    # the function returns NULL invisibly without erroring out.
    expect_null(mcview_notify("error", "no session here", where = "both"))
})
