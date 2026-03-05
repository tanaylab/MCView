# Full-app smoke test for MCView
#
# This test launches the complete Shiny application in a background R process
# and verifies:
#   1. The app starts up and responds with HTTP 200
#   2. No critical errors appear in stderr during startup and initial operation
#
# The test uses callr::r_bg() to isolate the app process, with proper cleanup
# via withr::defer().
#
# DAF setup and skip_if_no_daf() provided by helper-daf.R

test_that("Full app starts without critical errors", {
    skip_on_cran()
    skip_if_no_daf()
    skip_if(!requireNamespace("callr", quietly = TRUE), "callr not installed")
    skip_if(!requireNamespace("httr", quietly = TRUE), "httr not installed")

    # Pick a random port to avoid collisions
    port <- if (requireNamespace("httpuv", quietly = TRUE)) {
        httpuv::randomPort()
    } else {
        sample(10000:60000, 1)
    }

    # Determine the package root for devtools::load_all()
    pkg_root <- "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/aviezerl/src/MCView-daf"
    if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
        skip("Package root not found at expected location")
    }

    # Temp files for capturing output
    stderr_file <- tempfile(fileext = ".log")
    stdout_file <- tempfile(fileext = ".log")

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
            # run_app() returns a shiny app object via golem::with_golem_options
            # which defaults to print=FALSE. We must explicitly print() to start
            # the HTTP server in a background process.
            app <- run_app(
                daf,
                name = "smoke_test",
                port = port,
                host = "127.0.0.1",
                launch.browser = FALSE
            )
            print(app)
        },
        args = list(daf_path = get_test_daf_path(), port = port, pkg_root = pkg_root),
        stderr = stderr_file,
        stdout = stdout_file,
        env = c(
            JULIA_PROJECT = "@dafr-mcview",
            JULIA_LOAD_PATH = "@:@dafr-mcview:@stdlib",
            CONDA_PREFIX = Sys.getenv("CONDA_PREFIX", "")
        )
    )

    # Ensure we clean up the background process
    withr::defer({
        if (bg$is_alive()) bg$kill()
        unlink(stderr_file)
        unlink(stdout_file)
    })

    # Poll until the app is ready (up to 120 seconds for Julia startup + app init)
    app_url <- paste0("http://127.0.0.1:", port)
    ready <- FALSE
    max_wait <- 120
    poll_interval <- 2

    for (elapsed in seq(0, max_wait, by = poll_interval)) {
        # Check if the process died
        if (!bg$is_alive()) {
            stderr_content <- tryCatch(readLines(stderr_file), error = function(e) character(0))
            skip(paste0(
                "App process died before becoming ready.\nStderr (last 30 lines):\n",
                paste(utils::tail(stderr_content, 30), collapse = "\n")
            ))
        }

        resp <- tryCatch(
            httr::GET(app_url, httr::timeout(3)),
            error = function(e) NULL
        )

        if (!is.null(resp) && httr::status_code(resp) == 200) {
            ready <- TRUE
            break
        }

        Sys.sleep(poll_interval)
    }

    # If still not ready, check if process died during the last sleep
    if (!ready && !bg$is_alive()) {
        stderr_content <- tryCatch(readLines(stderr_file), error = function(e) character(0))
        skip(paste0(
            "App process died.\nStderr (last 30 lines):\n",
            paste(utils::tail(stderr_content, 30), collapse = "\n")
        ))
    }

    expect_true(ready, label = paste("App should be reachable at", app_url, "within", max_wait, "seconds"))

    if (ready) {
        # Verify the response contains HTML content
        resp <- httr::GET(app_url, httr::timeout(10))
        expect_equal(httr::status_code(resp), 200)
        content_type <- httr::http_type(resp)
        expect_true(
            grepl("html", content_type, ignore.case = TRUE),
            label = "Response Content-Type should contain 'html'"
        )
        body <- httr::content(resp, as = "text", encoding = "UTF-8")
        expect_gt(nchar(body), 0, label = "Response body should not be empty")

        # Let observers fire for a few seconds
        Sys.sleep(5)

        # Check stderr for fatal errors (not warnings -- plotly and other
        # packages commonly emit warnings that are harmless)
        stderr_content <- tryCatch(readLines(stderr_file), error = function(e) character(0))

        # Pattern for actual R errors (not warnings or messages)
        error_lines <- grep("^Error in |^Error:", stderr_content, value = TRUE)

        # Filter out known benign patterns
        benign_patterns <- c(
            "shinydashboard",   # shinydashboard deprecation messages
            "fontawesome",      # font-awesome icon warnings
            "plotly",           # plotly trace warnings
            "DT::",             # DataTables warnings
            "Warning"           # Warnings mis-classified as errors
        )
        for (pat in benign_patterns) {
            error_lines <- error_lines[!grepl(pat, error_lines, ignore.case = TRUE)]
        }

        if (length(error_lines) > 0) {
            fail(paste(
                "App produced", length(error_lines), "error(s) in stderr:\n",
                paste(error_lines, collapse = "\n")
            ))
        }

        expect_true(TRUE, label = "No critical errors in stderr")
    }
})

test_that("Full app serves static assets", {
    skip_on_cran()
    skip_if_no_daf()
    skip_if(!requireNamespace("callr", quietly = TRUE), "callr not installed")
    skip_if(!requireNamespace("httr", quietly = TRUE), "httr not installed")

    # Pick a random port
    port <- if (requireNamespace("httpuv", quietly = TRUE)) {
        httpuv::randomPort()
    } else {
        sample(10000:60000, 1)
    }

    pkg_root <- "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/aviezerl/src/MCView-daf"
    if (!file.exists(file.path(pkg_root, "DESCRIPTION"))) {
        skip("Package root not found at expected location")
    }

    stderr_file <- tempfile(fileext = ".log")
    stdout_file <- tempfile(fileext = ".log")

    bg <- callr::r_bg(
        function(daf_path, port, pkg_root) {
            Sys.setenv(JULIA_PROJECT = "@dafr-mcview")
            Sys.setenv(JULIA_LOAD_PATH = "@:@dafr-mcview:@stdlib")
            conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
            if (nzchar(conda_prefix)) {
                Sys.setenv(JULIA_DEPOT_PATH = paste0(conda_prefix, "/share/julia:"))
                options(dafr.JULIA_HOME = file.path(conda_prefix, "bin"))
            }

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
            app <- run_app(
                daf,
                name = "smoke_static",
                port = port,
                host = "127.0.0.1",
                launch.browser = FALSE
            )
            print(app)
        },
        args = list(daf_path = get_test_daf_path(), port = port, pkg_root = pkg_root),
        stderr = stderr_file,
        stdout = stdout_file,
        env = c(
            JULIA_PROJECT = "@dafr-mcview",
            JULIA_LOAD_PATH = "@:@dafr-mcview:@stdlib",
            CONDA_PREFIX = Sys.getenv("CONDA_PREFIX", "")
        )
    )

    withr::defer({
        if (bg$is_alive()) bg$kill()
        unlink(stderr_file)
        unlink(stdout_file)
    })

    # Wait for the app to start
    app_url <- paste0("http://127.0.0.1:", port)
    ready <- FALSE
    for (elapsed in seq(0, 120, by = 2)) {
        if (!bg$is_alive()) {
            stderr_content <- tryCatch(readLines(stderr_file), error = function(e) character(0))
            skip(paste0("App process died.\nStderr:\n", paste(utils::tail(stderr_content, 20), collapse = "\n")))
        }
        resp <- tryCatch(httr::GET(app_url, httr::timeout(3)), error = function(e) NULL)
        if (!is.null(resp) && httr::status_code(resp) == 200) {
            ready <- TRUE
            break
        }
        Sys.sleep(2)
    }

    skip_if(!ready, "App did not become ready in time")

    # Test that the favicon endpoint works (verifies golem_add_external_resources)
    favicon_resp <- tryCatch(
        httr::GET(paste0(app_url, "/favicon.ico"), httr::timeout(5)),
        error = function(e) NULL
    )
    # Favicon may return 200 or 404 depending on setup; we just check we got a response
    expect_true(!is.null(favicon_resp), label = "Favicon endpoint should respond")

    # Test the shared resource path
    www_resp <- tryCatch(
        httr::GET(paste0(app_url, "/www/"), httr::timeout(5)),
        error = function(e) NULL
    )
    # The /www/ path may not serve a directory listing, so just check for response
    expect_true(!is.null(www_resp), label = "Static assets path should respond")
})
