# helper-daf.R - Shared DAF gate for test files.
#
# Auto-sourced by testthat before any test file runs. Under native dafr there
# is no per-session setup to cache; we just gate on data availability.

#' Get the path to the test DAF dataset
get_test_daf_path <- function() {
    "/home/obk/data/mcview/metacells_clean"
}

#' Skip test if DAF test data or the dafr package is not available.
skip_if_no_daf <- function() {
    if (!dir.exists(get_test_daf_path())) {
        skip("Test DAF data not available")
    }
    if (!requireNamespace("dafr", quietly = TRUE)) {
        skip("dafr package not installed")
    }
}
