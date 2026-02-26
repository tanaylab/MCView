# helper-daf.R - Shared DAF setup for all test files
#
# This file is auto-sourced by testthat before any test file runs.
# It provides a single session-global DAF setup so that Julia/DAF
# initialization (which takes ~43s for precompilation) only happens ONCE
# per test session, not once per test file.

# ==============================================================================
# Shared State
# ==============================================================================

.daf_test_state <- new.env(parent = emptyenv())
.daf_test_state$attempted <- FALSE
.daf_test_state$ok <- FALSE
.daf_test_state$error_msg <- NULL

# ==============================================================================
# Shared Helpers
# ==============================================================================

#' Get the path to the test DAF dataset
get_test_daf_path <- function() {
    "/home/obk/data/mcview/metacells_clean"
}

#' Skip test if DAF is not available or Julia setup fails
#'
#' Calls dafr::setup_daf() at most once per session, caching the result.
skip_if_no_daf <- function() {
    test_daf_path <- get_test_daf_path()

    if (!dir.exists(test_daf_path)) {
        skip("Test DAF data not available")
    }
    if (!requireNamespace("dafr", quietly = TRUE)) {
        skip("dafr package not installed")
    }

    # Setup Julia/DAF once per session, caching the result
    if (!.daf_test_state$attempted) {
        .daf_test_state$attempted <- TRUE

        # Point dafr at the conda environment's Julia binary (if available)
        # Without this, dafr picks up whatever Julia is on PATH (e.g. juliaup's
        # 1.11.8) instead of the conda env's 1.12.5, causing silent fallback.
        conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
        if (nzchar(conda_prefix)) {
            options(dafr.JULIA_HOME = file.path(conda_prefix, "bin"))
            message(
                "[helper-daf] Using Julia from conda: ",
                file.path(conda_prefix, "bin", "julia")
            )
        } else {
            message("[helper-daf] CONDA_PREFIX not set; using Julia from PATH")
        }

        tryCatch(
            {
                dafr::setup_daf(pkg_check = FALSE, julia_environment = "custom")
                .daf_test_state$ok <- TRUE

                # Initialize Julia helpers so tests can exercise Julia-accelerated paths
                julia_ok <- tryCatch(
                    {
                        init_julia_helpers()
                    },
                    error = function(e) {
                        message("[helper-daf] Julia helpers init failed: ", e$message)
                        FALSE
                    }
                )
                if (isTRUE(julia_ok)) {
                    message("[helper-daf] Julia helpers initialized successfully")
                } else {
                    message("[helper-daf] Julia helpers not available; tests will use R fallbacks")
                }
            },
            error = function(e) {
                .daf_test_state$error_msg <- conditionMessage(e)
            }
        )
    }

    if (!.daf_test_state$ok) {
        skip(paste(
            "dafr::setup_daf() failed:",
            if (is.null(.daf_test_state$error_msg)) "unknown error" else .daf_test_state$error_msg
        ))
    }
}
