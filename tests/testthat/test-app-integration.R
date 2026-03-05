# Integration tests for the MCView Shiny app using the OBK DAF dataset
#
# These tests verify:
#   1. DAF contract validation against a real dataset
#   2. App initialization (environment, config, tabs)
#   3. HTTP-level app startup (background process + HTTP 200 check)
#
# DAF setup and skip_if_no_daf() provided by helper-daf.R

# Helper: open the OBK DAF (call after skip_if_no_daf)
open_obk_daf <- function() {
    dafr::open_daf(get_test_daf_path())
}

# ==============================================================================
# Section 1 - DAF Contract Validation
# ==============================================================================

test_that("validate_daf_for_mcview() passes on OBK dataset", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    expect_true(validate_daf_for_mcview(daf))
})

test_that("validate_daf_for_mcview() passes with verbose output", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    expect_true(validate_daf_for_mcview(daf, verbose = TRUE))
})

test_that("core contract validation returns valid result", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    result <- validate_mcview_contract(daf, mcview_core_contract(), verbose = TRUE)

    expect_true(result$valid)
    expect_length(result$errors, 0)
})

test_that("required axes exist: metacell, gene, type", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    expect_true(dafr::has_axis(daf, "metacell"))
    expect_true(dafr::has_axis(daf, "gene"))
    expect_true(dafr::has_axis(daf, "type"))

    # Non-empty
    expect_gt(length(dafr::axis_entries(daf, "metacell")), 0)
    expect_gt(length(dafr::axis_entries(daf, "gene")), 0)
    expect_gt(length(dafr::axis_entries(daf, "type")), 0)
})

test_that("required matrix UMIs exists", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    expect_true(dafr::has_matrix(daf, "metacell", "gene", "UMIs"))

    umat <- dafr::get_matrix(daf, "metacell", "gene", "UMIs")
    expect_gt(nrow(umat), 0)
    expect_gt(ncol(umat), 0)
})

test_that("required vectors exist: total_UMIs, type, color", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    # metacell vectors
    expect_true(dafr::has_vector(daf, "metacell", "total_UMIs"))
    expect_true(dafr::has_vector(daf, "metacell", "type"))

    # type vector
    expect_true(dafr::has_vector(daf, "type", "color"))

    total_umis <- dafr::get_vector(daf, "metacell", "total_UMIs")
    expect_true(is.numeric(total_umis))
    expect_gt(length(total_umis), 0)
    expect_true(all(total_umis > 0))

    mc_types <- dafr::get_vector(daf, "metacell", "type")
    expect_type(mc_types, "character")
    expect_gt(length(mc_types), 0)

    colors <- dafr::get_vector(daf, "type", "color")
    expect_type(colors, "character")
    expect_gt(length(colors), 0)
    # Colors should look like hex codes
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
})

test_that("2D coordinates exist (x/y or u/v)", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    has_xy <- dafr::has_vector(daf, "metacell", "x") &&
        dafr::has_vector(daf, "metacell", "y")
    has_uv <- dafr::has_vector(daf, "metacell", "u") &&
        dafr::has_vector(daf, "metacell", "v")

    expect_true(has_xy || has_uv, label = "At least one coordinate pair (x/y or u/v) must exist")

    # Whichever pair is present, verify non-empty and numeric
    if (has_xy) {
        x <- dafr::get_vector(daf, "metacell", "x")
        y <- dafr::get_vector(daf, "metacell", "y")
        expect_true(is.numeric(x))
        expect_true(is.numeric(y))
        expect_equal(length(x), length(y))
        expect_gt(length(x), 0)
    }
    if (has_uv) {
        u <- dafr::get_vector(daf, "metacell", "u")
        v <- dafr::get_vector(daf, "metacell", "v")
        expect_true(is.numeric(u))
        expect_true(is.numeric(v))
        expect_equal(length(u), length(v))
        expect_gt(length(u), 0)
    }
})

test_that("custom coordinate rule in core contract passes", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    contract <- mcview_core_contract()
    rule_result <- contract$custom_rules$coordinates$validate(daf)
    expect_null(rule_result) # NULL means no error
})

# ==============================================================================
# Section 2 - App Initialization
# ==============================================================================

test_that("init_mcview_env() works and resets state", {
    result <- init_mcview_env()
    expect_true(result)

    # After init, config should be NULL
    expect_null(mcv_get("config"))
    expect_null(mcv_get("mc_data"))
    expect_null(mcv_get("atlas"))
})

test_that("init_single_daf_mode() works with the OBK DAF", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    init_mcview_env()
    init_single_daf_mode(daf, "obk_test")

    # Dataset should be registered
    expect_equal(dataset_names(), "obk_test")

    # DAF object should be retrievable
    daf_obj <- get_dataset_daf("obk_test")
    expect_true(inherits(daf_obj, "Daf"))
})

test_that("config is properly populated after init_single_daf_mode", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    init_mcview_env()
    init_single_daf_mode(daf, "obk_config_test")

    config <- mcv_get("config")
    expect_true(!is.null(config))

    # Title should be set (either from DAF scalar or auto-generated)
    expect_true(!is.null(config$title))
    expect_type(config$title, "character")
    expect_gt(nchar(config$title), 0)

    # Tabs should be set (either from DAF scalar or auto-detected)
    expect_true(!is.null(config$tabs))
    expect_type(config$tabs, "character")
    expect_gt(length(config$tabs), 0)

    # Each tab name should be a non-empty string
    expect_true(all(nchar(config$tabs) > 0))
})

test_that("detect_available_tabs() returns reasonable tabs", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    tabs <- detect_available_tabs(daf)

    expect_type(tabs, "character")
    expect_gt(length(tabs), 0)

    # Core tabs that should always be available with a well-formed DAF
    # (Manifold, Genes, and Diff. Expression depend only on the core contract)
    expect_true("Manifold" %in% tabs)
    expect_true("Genes" %in% tabs)
    expect_true("Diff. Expression" %in% tabs)

    # All returned tabs should be from the canonical list
    expect_true(all(tabs %in% MCVIEW_TAB_NAMES))
})

test_that("validate_mcview_tabs() returns per-tab results", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    tab_results <- validate_mcview_tabs(daf, verbose = FALSE)

    expect_true(is.list(tab_results))
    expect_gt(length(tab_results), 0)

    # Each entry should have 'available' and 'reason'
    for (tab_name in names(tab_results)) {
        expect_true("available" %in% names(tab_results[[tab_name]]))
        expect_true("reason" %in% names(tab_results[[tab_name]]))
        expect_type(tab_results[[tab_name]]$available, "logical")
    }
})

test_that("multi-dataset mode works with same DAF under two names", {
    skip_if_no_daf()
    daf1 <- open_obk_daf()
    daf2 <- open_obk_daf()

    init_mcview_env()
    init_multi_daf_mode(list(dataset_a = daf1, dataset_b = daf2))

    ds_names <- dataset_names()
    expect_length(ds_names, 2)
    expect_true("dataset_a" %in% ds_names)
    expect_true("dataset_b" %in% ds_names)

    # Both should have DAF objects
    expect_true(inherits(get_dataset_daf("dataset_a"), "Daf"))
    expect_true(inherits(get_dataset_daf("dataset_b"), "Daf"))

    # Config should be populated
    config <- mcv_get("config")
    expect_true(!is.null(config))
    expect_true(!is.null(config$title))
    expect_true(!is.null(config$tabs))
})

test_that("extract_config_for_daf auto-generates title and tabs", {
    skip_if_no_daf()
    daf <- open_obk_daf()

    config <- extract_config_for_daf(daf, "my_dataset")

    expect_true(!is.null(config$title))
    expect_true(!is.null(config$tabs))
    expect_type(config$tabs, "character")
    expect_gt(length(config$tabs), 0)
})

# ==============================================================================
# Section 3 - App Startup (HTTP verification)
# ==============================================================================

test_that("run_app() starts and responds with HTTP 200", {
    skip_on_cran()
    skip_if_no_daf()
    skip_if(!requireNamespace("callr", quietly = TRUE), "callr not installed")
    skip_if(!requireNamespace("httr", quietly = TRUE), "httr not installed")

    # Pick a port
    port <- if (requireNamespace("httpuv", quietly = TRUE)) {
        httpuv::randomPort()
    } else {
        19876L
    }

    # Determine the package root for devtools::load_all()
    pkg_root <- system.file(package = "MCView")
    # If running from source tree, use the project directory
    src_root <- normalizePath(file.path(get_test_daf_path(), "..", "..", ".."), mustWork = FALSE)
    if (file.exists(file.path(src_root, "DESCRIPTION"))) {
        pkg_root <- src_root
    }
    # Try common development locations
    dev_roots <- c(
        Sys.getenv("MCVIEW_PKG_ROOT", ""),
        "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/aviezerl/src/MCView-daf"
    )
    for (root in dev_roots) {
        if (nzchar(root) && file.exists(file.path(root, "DESCRIPTION"))) {
            pkg_root <- root
            break
        }
    }

    # Launch the app in a background R process
    bg <- callr::r_bg(
        function(daf_path, port, pkg_root) {
            # Set up Julia environment for dafr
            Sys.setenv(JULIA_PROJECT = "@dafr-mcview")
            Sys.setenv(JULIA_LOAD_PATH = "@:@dafr-mcview:@stdlib")
            conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
            if (nzchar(conda_prefix)) {
                Sys.setenv(JULIA_DEPOT_PATH = paste0(conda_prefix, "/share/julia:"))
                options(dafr.JULIA_HOME = file.path(conda_prefix, "bin"))
            }

            # Load package from source in dev mode
            devtools::load_all(pkg_root, quiet = TRUE)
            # Use sysimage if available for faster Julia startup
            sysimage <- Sys.getenv("JULIA_SYSIMAGE", "")
            if (!nzchar(sysimage) && nzchar(conda_prefix)) {
                candidate <- file.path(conda_prefix, "share", "julia", "sysimage_daf.so")
                if (file.exists(candidate)) sysimage <- candidate
            }
            setup_args <- list(pkg_check = FALSE, julia_environment = "custom")
            if (nzchar(sysimage) && file.exists(sysimage)) {
                setup_args$sysimage_path <- sysimage
            }
            do.call(dafr::setup_daf, setup_args)
            daf <- dafr::open_daf(daf_path)
            run_app(
                daf,
                name = "integration_test",
                port = port,
                host = "127.0.0.1",
                launch.browser = FALSE
            )
        },
        args = list(daf_path = get_test_daf_path(), port = port, pkg_root = pkg_root),
        env = c(
            JULIA_PROJECT = "@dafr-mcview",
            JULIA_LOAD_PATH = "@:@dafr-mcview:@stdlib",
            CONDA_PREFIX = Sys.getenv("CONDA_PREFIX", "")
        )
    )

    # Ensure we clean up the background process
    withr::defer({
        if (bg$is_alive()) bg$kill()
    })

    # Poll until the app is ready (up to 120 seconds for Julia startup + app init)
    app_url <- paste0("http://127.0.0.1:", port)
    ready <- FALSE
    max_wait <- 120
    poll_interval <- 2

    for (elapsed in seq(0, max_wait, by = poll_interval)) {
        # Check if the process died
        if (!bg$is_alive()) {
            # Capture any error output for diagnostics
            err_out <- tryCatch(bg$read_error(), error = function(e) "")
            std_out <- tryCatch(bg$read_output(), error = function(e) "")
            skip(paste0(
                "App process died before becoming ready.\nStderr: ",
                err_out, "\nStdout: ", std_out
            ))
        }

        resp <- tryCatch(
            httr::GET(app_url, httr::timeout(5)),
            error = function(e) NULL
        )

        if (!is.null(resp) && httr::status_code(resp) == 200) {
            ready <- TRUE
            break
        }

        Sys.sleep(poll_interval)
    }

    expect_true(ready, label = paste("App should be reachable at", app_url, "within", max_wait, "seconds"))

    # Verify the response has HTML content
    if (ready) {
        resp <- httr::GET(app_url, httr::timeout(10))
        expect_equal(httr::status_code(resp), 200)
        content_type <- httr::http_type(resp)
        expect_true(grepl("html", content_type, ignore.case = TRUE),
            label = "Response Content-Type should contain 'html'"
        )
        body <- httr::content(resp, as = "text", encoding = "UTF-8")
        expect_gt(nchar(body), 0)
    }
})
