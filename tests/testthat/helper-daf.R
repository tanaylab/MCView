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

#' Stand up a fresh MCView session for cold-path perf tests.
#'
#' Opens the OBK DAF, initialises the MCView environment, and runs the
#' init steps the app server normally performs. Caches like `mc_mat`
#' and `mc_egc_full` are left empty so tests can exercise the cold
#' targeted-DAF-query paths from Phase 4.
#'
#' @return The opened DAF object (`daf_obj`), invisibly. Callers usually
#'   bind it to introspect axes via `dafr::axis_entries()`.
setup_cold_session <- function() {
    daf_obj <- dafr::open_daf(get_test_daf_path())
    init_mcview_env()
    init_single_daf_mode(daf_obj, "data", NULL, FALSE)
    init_defs()
    config_shiny_cache()
    invisible(daf_obj)
}
