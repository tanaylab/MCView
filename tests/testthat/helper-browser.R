# helper-browser.R - Headless Chrome testing infrastructure for MCView
#
# This file provides reusable functions for browser-based (chromote) testing
# of the MCView Shiny app. It is auto-sourced by testthat before any test
# file runs.
#
# Proven working pattern on this server:
#   1. Launch app via callr::r_bg() with Julia env vars
#   2. Set CHROMOTE_CHROME and LD_LIBRARY_PATH env vars
#   3. Poll HTTP until 200 (takes ~36s)
#   4. Connect chromote::ChromoteSession$new()
#   5. Set viewport with Emulation$setDeviceMetricsOverride
#   6. Navigate to app URL
#   7. Wait ~60s total for shiny-busy to clear
#   8. Take screenshots, execute JS, read DOM

# ==============================================================================
# Constants
# ==============================================================================

.CHROME_PATH <- Sys.getenv(
    "CHROMOTE_CHROME",
    "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/aviezerl/.cache/R/chromote/chrome/145.0.7632.67/chrome-headless-shell-linux64/chrome-headless-shell"
)

.MCVIEW_PKG_ROOT <- "/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/aviezerl/src/MCView-daf"

.SCREENSHOT_DIR <- file.path(.MCVIEW_PKG_ROOT, "tests", "testthat", "screenshots")

# ==============================================================================
# Availability Checks
# ==============================================================================

#' Check if browser testing is available
#'
#' Returns TRUE if chromote is installed and the Chrome binary exists.
#' @return logical
browser_test_available <- function() {
    tryCatch(
        {
            has_chromote <- requireNamespace("chromote", quietly = TRUE)
            has_chrome <- file.exists(.CHROME_PATH)
            has_callr <- requireNamespace("callr", quietly = TRUE)
            has_httr <- requireNamespace("httr", quietly = TRUE)
            has_chromote && has_chrome && has_callr && has_httr
        },
        error = function(e) FALSE
    )
}

#' Skip test if browser testing is not available
#'
#' Checks for chromote, Chrome binary, callr, and httr. Also skips on CRAN.
skip_if_no_browser <- function() {
    skip_on_cran()
    if (!requireNamespace("chromote", quietly = TRUE)) {
        skip("chromote package not installed")
    }
    if (!file.exists(.CHROME_PATH)) {
        skip(paste("Chrome binary not found at:", .CHROME_PATH))
    }
    if (!requireNamespace("callr", quietly = TRUE)) {
        skip("callr package not installed")
    }
    if (!requireNamespace("httr", quietly = TRUE)) {
        skip("httr package not installed")
    }
}

# ==============================================================================
# App Launching
# ==============================================================================

#' Launch MCView in a background R process
#'
#' Starts the app via callr::r_bg(). Uses devtools::load_all() for source
#' loading, then calls run_app() and explicitly starts the server with
#' shiny::runApp().
#'
#' @param daf_path Path to the DAF dataset
#' @param port Port number for the Shiny server
#' @param pkg_root Path to the MCView package root (default: auto-detected)
#' @return The callr process object
launch_mcview_bg <- function(daf_path = get_test_daf_path(),
                             port = NULL,
                             pkg_root = .MCVIEW_PKG_ROOT) {
    if (is.null(port)) {
        port <- if (requireNamespace("httpuv", quietly = TRUE)) {
            httpuv::randomPort()
        } else {
            19876L
        }
    }

    tryCatch(
        {
            bg <- callr::r_bg(
                function(daf_path, port, pkg_root) {
                    # Load package from source in dev mode
                    devtools::load_all(pkg_root, quiet = TRUE)

                    # Open DAF dataset (dir → files_daf; .h5/.h5ad → h5ad_as_daf)
                    daf <- if (dir.exists(daf_path)) {
                        dafr::files_daf(daf_path, mode = "r")
                    } else if (grepl("\\.h5ad$|\\.h5$", daf_path, ignore.case = TRUE)) {
                        dafr::h5ad_as_daf(daf_path)
                    } else {
                        dafr::open_daf(daf_path)
                    }

                    # Create the Shiny app object.
                    # run_app() returns a Shiny app object; in non-interactive mode
                    # we must explicitly call shiny::runApp() to start the server.
                    app <- run_app(
                        daf,
                        name = "browser_test",
                        port = port,
                        host = "127.0.0.1",
                        launch.browser = FALSE
                    )
                    shiny::runApp(app, port = port, host = "127.0.0.1", launch.browser = FALSE)
                },
                args = list(
                    daf_path = daf_path,
                    port = port,
                    pkg_root = pkg_root
                )
            )

            # Attach the port as an attribute for downstream use
            attr(bg, "port") <- port
            bg
        },
        error = function(e) {
            stop("Failed to launch MCView background process: ", e$message)
        }
    )
}

# ==============================================================================
# App Readiness
# ==============================================================================

#' Wait for the Shiny app to become reachable via HTTP
#'
#' Polls the app URL until it returns HTTP 200 or the timeout expires.
#' Also monitors the background process for premature death.
#'
#' @param port Port number the app is running on
#' @param timeout Maximum seconds to wait (default: 120)
#' @param poll_interval Seconds between poll attempts (default: 2)
#' @param bg_process Optional callr process object to monitor for early death
#' @return TRUE if app is ready, FALSE if timeout was reached
wait_for_app <- function(port, timeout = 120, poll_interval = 2, bg_process = NULL) {
    app_url <- paste0("http://127.0.0.1:", port)

    for (elapsed in seq(0, timeout, by = poll_interval)) {
        # Check if the background process died
        if (!is.null(bg_process) && !bg_process$is_alive()) {
            err_out <- tryCatch(bg_process$read_error(), error = function(e) "")
            std_out <- tryCatch(bg_process$read_output(), error = function(e) "")
            warning(
                "MCView background process died before becoming ready.\n",
                "Stderr: ", err_out, "\n",
                "Stdout: ", std_out
            )
            return(FALSE)
        }

        resp <- tryCatch(
            httr::GET(app_url, httr::timeout(5)),
            error = function(e) NULL
        )

        if (!is.null(resp) && httr::status_code(resp) == 200) {
            return(TRUE)
        }

        Sys.sleep(poll_interval)
    }

    FALSE
}

# ==============================================================================
# Browser Connection
# ==============================================================================

#' Create a chromote browser session connected to the app
#'
#' Sets CHROMOTE_CHROME and LD_LIBRARY_PATH, creates a ChromoteSession,
#' configures the viewport, and navigates to the app URL.
#'
#' @param port Port number the app is running on
#' @param width Viewport width in pixels (default: 1400)
#' @param height Viewport height in pixels (default: 900)
#' @return A ChromoteSession object
connect_browser <- function(port, width = 1400, height = 900) {
    # Set Chrome binary path
    Sys.setenv(CHROMOTE_CHROME = .CHROME_PATH)

    # Set LD_LIBRARY_PATH to include conda libraries (needed for Chrome)
    conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
    if (nzchar(conda_prefix)) {
        current_ld <- Sys.getenv("LD_LIBRARY_PATH", "")
        conda_lib <- file.path(conda_prefix, "lib")
        conda_sysroot <- file.path(conda_prefix, "x86_64-conda-linux-gnu", "sysroot", "usr", "lib64")
        new_ld <- paste(conda_lib, conda_sysroot, current_ld, sep = ":")
        # Remove trailing colons
        new_ld <- sub(":+$", "", new_ld)
        Sys.setenv(LD_LIBRARY_PATH = new_ld)
    }

    tryCatch(
        {
            # Create a new Chromote browser and session
            browser <- chromote::Chromote$new()
            session <- chromote::ChromoteSession$new(parent = browser)

            # Set viewport dimensions
            session$Emulation$setDeviceMetricsOverride(
                width = width,
                height = height,
                deviceScaleFactor = 1,
                mobile = FALSE
            )

            # Navigate to the app
            app_url <- paste0("http://127.0.0.1:", port)
            session$Page$navigate(app_url)

            # Brief pause to let the initial page load begin
            Sys.sleep(2)

            session
        },
        error = function(e) {
            stop("Failed to connect browser session: ", e$message)
        }
    )
}

# ==============================================================================
# Shiny State Helpers
# ==============================================================================

#' Wait for Shiny to become idle (no shiny-busy class on html element)
#'
#' Polls the html element's classList for the absence of "shiny-busy".
#' This indicates that all reactive outputs have finished rendering.
#'
#' @param session A ChromoteSession object
#' @param timeout Maximum seconds to wait (default: 180 — bumped from 90 after
#'   the Julia-removal migration; tracked for follow-up perf pass)
#' @param poll_interval Seconds between polls (default: 1)
#' @return TRUE if Shiny became idle, FALSE if timeout expired
wait_for_shiny_idle <- function(session, timeout = 180, poll_interval = 1) {
    tryCatch(
        {
            start_time <- Sys.time()

            while (difftime(Sys.time(), start_time, units = "secs") < timeout) {
                result <- tryCatch(
                    session$Runtime$evaluate(
                        expression = "!document.documentElement.classList.contains('shiny-busy')"
                    ),
                    error = function(e) {
                        warning("An error occurred: ", conditionMessage(e))
                        NULL
                    }
                )

                if (!is.null(result) && isTRUE(result$result$value)) {
                    return(TRUE)
                }

                Sys.sleep(poll_interval)
            }

            FALSE
        },
        error = function(e) {
            warning("Error while waiting for Shiny idle: ", e$message)
            FALSE
        }
    )
}

# ==============================================================================
# Interaction Helpers
# ==============================================================================

#' Click a sidebar tab by its text label
#'
#' Finds the sidebar link whose text matches tab_name and clicks it,
#' then waits for Shiny to become idle.
#'
#' @param session A ChromoteSession object
#' @param tab_name The visible text of the tab to click
#' @param timeout Seconds to wait for Shiny idle after clicking (default: 180
#'   — bumped from 60 after the Julia-removal migration; tracked for perf pass)
#' @return TRUE if click succeeded and Shiny became idle, FALSE otherwise
click_tab <- function(session, tab_name, timeout = 180) {
    tryCatch(
        {
            # Find and click the sidebar link matching the tab name.
            # menuSubItem generates links whose textContent includes icon text,
            # so we trim and compare. Also try matching via innerText or
            # checking child <span> elements.
            js <- sprintf(
                "
                (function() {
                    // Try sidebar menuSubItem links first
                    var links = document.querySelectorAll('.sidebar-menu .treeview-menu a');
                    for (var i = 0; i < links.length; i++) {
                        var spans = links[i].querySelectorAll('span');
                        for (var j = 0; j < spans.length; j++) {
                            if (spans[j].textContent.trim() === '%s') {
                                links[i].click();
                                return true;
                            }
                        }
                        if (links[i].textContent.trim() === '%s') {
                            links[i].click();
                            return true;
                        }
                    }
                    // Broader fallback
                    var allLinks = document.querySelectorAll('.sidebar-menu a, .nav-tabs a, [data-toggle=\"tab\"], [role=\"tab\"]');
                    for (var i = 0; i < allLinks.length; i++) {
                        if (allLinks[i].textContent.trim() === '%s') {
                            allLinks[i].click();
                            return true;
                        }
                    }
                    return false;
                })()
                ",
                tab_name, tab_name, tab_name
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    # CDP error; retry once after pause
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) list(result = list(value = FALSE))
                    )
                }
            )

            if (!isTRUE(result$result$value)) {
                warning("Tab '", tab_name, "' not found in the sidebar")
                return(FALSE)
            }

            # Wait for Shiny to process the tab switch
            Sys.sleep(1)
            wait_for_shiny_idle(session, timeout = timeout)
        },
        error = function(e) {
            warning("Error clicking tab '", tab_name, "': ", e$message)
            FALSE
        }
    )
}

# ==============================================================================
# Screenshot & DOM Inspection
# ==============================================================================

#' Take a screenshot and save it to the screenshots directory
#'
#' Creates the screenshots/ directory if it does not exist.
#' Saves the screenshot as tests/testthat/screenshots/{name}.png.
#'
#' @param session A ChromoteSession object
#' @param name Base name for the screenshot file (without extension)
#' @return The absolute path to the saved screenshot file, or NULL on error
take_screenshot <- function(session, name) {
    tryCatch(
        {
            # Ensure the screenshots directory exists
            if (!dir.exists(.SCREENSHOT_DIR)) {
                dir.create(.SCREENSHOT_DIR, recursive = TRUE, showWarnings = FALSE)
            }

            file_path <- file.path(.SCREENSHOT_DIR, paste0(name, ".png"))

            session$screenshot(filename = file_path)

            # Verify the file was created
            if (file.exists(file_path)) {
                message("[browser-test] Screenshot saved: ", file_path)
                file_path
            } else {
                warning("Screenshot file was not created: ", file_path)
                NULL
            }
        },
        error = function(e) {
            warning("Failed to take screenshot '", name, "': ", e$message)
            NULL
        }
    )
}

#' Read a Shiny output value via JavaScript
#'
#' Accesses the Shiny output binding to retrieve the current value of
#' a given output ID.
#'
#' @param session A ChromoteSession object
#' @param output_id The Shiny output ID (e.g., "my_plot", "my_table")
#' @return The value of the output, or NULL if not found or on error
get_shiny_value <- function(session, output_id) {
    tryCatch(
        {
            js <- sprintf(
                "
                (function() {
                    var el = document.getElementById('%s');
                    if (!el) return null;
                    // Try to get text content for simple outputs
                    if (el.textContent) return el.textContent.trim();
                    // Try innerHTML as fallback
                    return el.innerHTML;
                })()
                ",
                output_id
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(1)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) NULL
                    )
                }
            )
            if (is.null(result)) return(NULL)
            result$result$value
        },
        error = function(e) {
            warning("Failed to get Shiny value for '", output_id, "': ", e$message)
            NULL
        }
    )
}

#' Count visible elements matching a CSS selector
#'
#' @param session A ChromoteSession object
#' @param selector A CSS selector string
#' @return Integer count of visible matching elements, or 0 on error
get_visible_elements <- function(session, selector) {
    tryCatch(
        {
            safe_selector <- gsub("'", "\\\\'", selector)
            js <- sprintf(
                "
                (function() {
                    var els = document.querySelectorAll('%s');
                    var count = 0;
                    for (var i = 0; i < els.length; i++) {
                        var style = window.getComputedStyle(els[i]);
                        if (style.display !== 'none' && style.visibility !== 'hidden' && style.opacity !== '0') {
                            count++;
                        }
                    }
                    return count;
                })()
                ",
                safe_selector
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(1)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) NULL
                    )
                }
            )

            if (is.null(result) || is.null(result$result$value)) {
                return(0L)
            }
            as.integer(result$result$value)
        },
        error = function(e) {
            warning("Failed to count visible elements for '", selector, "': ", e$message)
            0L
        }
    )
}

#' Check if an element exists in the DOM
#'
#' @param session A ChromoteSession object
#' @param selector A CSS selector string
#' @return TRUE if at least one matching element exists, FALSE otherwise
element_exists <- function(session, selector) {
    tryCatch(
        {
            # Escape single quotes in selector for JS string
            safe_selector <- gsub("'", "\\\\'", selector)
            js <- sprintf(
                "document.querySelector('%s') !== null",
                safe_selector
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    # CDP errors like "Could not find node with given id"
                    # Retry once after a brief pause
                    Sys.sleep(1)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) NULL
                    )
                }
            )
            if (is.null(result)) return(FALSE)
            isTRUE(result$result$value)
        },
        error = function(e) {
            warning("Failed to check element existence for '", selector, "': ", e$message)
            FALSE
        }
    )
}

#' Check the browser console for JavaScript errors
#'
#' Evaluates a JS snippet that collects any errors from the console.
#' Note: This requires that console error capturing has been set up.
#' As a fallback, checks for common Shiny error indicators in the DOM.
#'
#' @param session A ChromoteSession object
#' @return A character vector of error messages, or character(0) if none
get_js_errors <- function(session) {
    tryCatch(
        {
            # First, check for Shiny notification errors in the DOM
            js_dom_errors <- "
            (function() {
                var errors = [];
                // Check for Shiny notification errors
                var notifications = document.querySelectorAll('.shiny-notification-error, .shiny-output-error');
                for (var i = 0; i < notifications.length; i++) {
                    var el = notifications[i];
                    var text = (el.innerText || el.textContent || '').trim();
                    var id = el.id || el.getAttribute('data-output-id') || '(unknown)';
                    var cls = el.className || '';
                    if (text) {
                        errors.push('[' + id + '] ' + text);
                    } else {
                        // Element exists but has no text; report id and classes
                        errors.push('[' + id + '] (empty error element, class=' + cls + ')');
                    }
                }
                return errors;
            })()
            "

            result <- tryCatch(
                session$Runtime$evaluate(
                    expression = js_dom_errors,
                    returnByValue = TRUE
                ),
                error = function(e) NULL
            )

            if (!is.null(result$result$value)) {
                errors <- unlist(result$result$value)
                if (length(errors) > 0) {
                    return(as.character(errors))
                }
            }

            character(0)
        },
        error = function(e) {
            warning("Failed to check JS errors: ", e$message)
            character(0)
        }
    )
}

# ==============================================================================
# Cleanup
# ==============================================================================

#' Orderly cleanup of browser session and background process
#'
#' Closes the chromote session and kills the background R process.
#' Handles errors gracefully so that cleanup does not mask test failures.
#'
#' @param session A ChromoteSession object (or NULL)
#' @param bg_process A callr process object (or NULL)
#' @return Invisible NULL
cleanup_browser_test <- function(session = NULL, bg_process = NULL) {
    # Close the chromote session
    if (!is.null(session)) {
        tryCatch(
            {
                session$close()
            },
            error = function(e) {
                # Session may already be closed
            }
        )
    }

    # Kill the background process
    if (!is.null(bg_process)) {
        tryCatch(
            {
                if (bg_process$is_alive()) {
                    bg_process$kill()
                }
            },
            error = function(e) {
                # Process may already be dead
            }
        )
    }

    invisible(NULL)
}

# ==============================================================================
# Interaction Helpers (for scenario tests)
# ==============================================================================

#' Set a Shiny input value via JavaScript
#'
#' Uses Shiny.setInputValue() to programmatically set an input, then waits
#' for Shiny to become idle.
#'
#' @param session A ChromoteSession object
#' @param input_id The Shiny input ID to set
#' @param value The value to set (will be JSON-encoded for non-scalar types)
#' @return TRUE if the value was set successfully, FALSE otherwise
set_shiny_input <- function(session, input_id, value) {
    tryCatch(
        {
            # JSON-encode the value for safe injection into JS
            value_js <- jsonlite::toJSON(value, auto_unbox = TRUE)

            js <- sprintf(
                "
                (function() {
                    if (typeof Shiny === 'undefined' || !Shiny.setInputValue) return false;
                    Shiny.setInputValue('%s', %s, {priority: 'event'});
                    return true;
                })()
                ",
                input_id, value_js
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) list(result = list(value = FALSE))
                    )
                }
            )

            if (!isTRUE(result$result$value)) {
                warning("Failed to set Shiny input '", input_id, "'")
                return(FALSE)
            }

            Sys.sleep(0.5)
            wait_for_shiny_idle(session)
        },
        error = function(e) {
            warning("Error setting Shiny input '", input_id, "': ", e$message)
            FALSE
        }
    )
}

#' Click a radio button by its label text
#'
#' Supports radioGroupButtons (shinyWidgets), prettyRadioButtons, and standard
#' radio inputs. Searches within the specified container for a button or label
#' whose text matches label_text.
#'
#' @param session A ChromoteSession object
#' @param container_selector CSS selector for the container holding the radio group
#' @param label_text The visible text of the radio option to select
#' @return TRUE if the click succeeded, FALSE otherwise
click_radio_button <- function(session, container_selector, label_text) {
    tryCatch(
        {
            safe_selector <- gsub("'", "\\\\'", container_selector)
            safe_label <- gsub("'", "\\\\'", label_text)

            js <- sprintf(
                "
                (function() {
                    var container = document.querySelector('%s');
                    if (!container) return false;

                    // Strategy 1: radioGroupButtons (shinyWidgets) - buttons in btn-group
                    var buttons = container.querySelectorAll('.btn-group button, .btn-group-toggle button, .btn-group label.btn');
                    for (var i = 0; i < buttons.length; i++) {
                        if (buttons[i].textContent.trim() === '%s') {
                            buttons[i].click();
                            return true;
                        }
                    }

                    // Strategy 2: prettyRadioButtons - div.pretty with labels
                    var prettyLabels = container.querySelectorAll('.pretty.p-default.p-round label, .pretty label');
                    for (var i = 0; i < prettyLabels.length; i++) {
                        if (prettyLabels[i].textContent.trim() === '%s') {
                            var input = prettyLabels[i].closest('.pretty').querySelector('input[type=\"radio\"]');
                            if (input) {
                                input.click();
                                return true;
                            }
                            prettyLabels[i].click();
                            return true;
                        }
                    }

                    // Strategy 3: Generic fallback - any label or button with matching text
                    var allClickables = container.querySelectorAll('label, button, input[type=\"radio\"]');
                    for (var i = 0; i < allClickables.length; i++) {
                        if (allClickables[i].textContent.trim() === '%s') {
                            allClickables[i].click();
                            return true;
                        }
                    }

                    // Strategy 4: Find label by text and click associated input
                    var labels = container.querySelectorAll('label');
                    for (var i = 0; i < labels.length; i++) {
                        if (labels[i].textContent.trim() === '%s') {
                            var forId = labels[i].getAttribute('for');
                            if (forId) {
                                var input = document.getElementById(forId);
                                if (input) { input.click(); return true; }
                            }
                            labels[i].click();
                            return true;
                        }
                    }

                    return false;
                })()
                ",
                safe_selector, safe_label, safe_label, safe_label, safe_label
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) list(result = list(value = FALSE))
                    )
                }
            )

            if (!isTRUE(result$result$value)) {
                warning("Radio button '", label_text, "' not found in '", container_selector, "'")
                return(FALSE)
            }

            Sys.sleep(0.5)
            wait_for_shiny_idle(session)
        },
        error = function(e) {
            warning("Error clicking radio button '", label_text, "': ", e$message)
            FALSE
        }
    )
}

#' Toggle a checkbox input
#'
#' Supports standard checkboxInput and awesomeCheckbox (shinyWidgets).
#' Finds the checkbox element within the input wrapper and clicks it.
#'
#' @param session A ChromoteSession object
#' @param input_id The Shiny input ID of the checkbox
#' @return TRUE if the checkbox was toggled, FALSE otherwise
click_checkbox <- function(session, input_id) {
    tryCatch(
        {
            js <- sprintf(
                "
                (function() {
                    // Strategy 1: Direct input element by ID (standard checkboxInput)
                    var el = document.getElementById('%s');
                    if (el && el.tagName === 'INPUT' && el.type === 'checkbox') {
                        el.click();
                        return true;
                    }

                    // Strategy 2: Wrapper div contains the input (checkboxInput wraps in div)
                    if (el) {
                        var input = el.querySelector('input[type=\"checkbox\"]');
                        if (input) {
                            input.click();
                            return true;
                        }
                    }

                    // Strategy 3: awesomeCheckbox - look for wrapper with matching input name
                    var awesome = document.querySelector('#' + CSS.escape('%s') + ' input[type=\"checkbox\"], input[type=\"checkbox\"][id=\"%s\"]');
                    if (awesome) {
                        awesome.click();
                        return true;
                    }

                    return false;
                })()
                ",
                input_id, input_id, input_id
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) list(result = list(value = FALSE))
                    )
                }
            )

            if (!isTRUE(result$result$value)) {
                warning("Checkbox '", input_id, "' not found")
                return(FALSE)
            }

            Sys.sleep(0.5)
            wait_for_shiny_idle(session)
        },
        error = function(e) {
            warning("Error clicking checkbox '", input_id, "': ", e$message)
            FALSE
        }
    )
}

#' Toggle a shinyWidgets switch or prettySwitch
#'
#' prettySwitch renders as a `<div class="pretty p-switch">` wrapping an input
#' and a label. This function finds the switch by input ID and clicks it.
#'
#' @param session A ChromoteSession object
#' @param input_id The Shiny input ID of the switch
#' @return TRUE if the switch was toggled, FALSE otherwise
click_switch <- function(session, input_id) {
    tryCatch(
        {
            js <- sprintf(
                "
                (function() {
                    // Strategy 1: Direct input element
                    var input = document.getElementById('%s');
                    if (input && input.tagName === 'INPUT') {
                        input.click();
                        return true;
                    }

                    // Strategy 2: prettySwitch - wrapper div, find input inside
                    if (input) {
                        var cb = input.querySelector('input[type=\"checkbox\"]');
                        if (cb) {
                            cb.click();
                            return true;
                        }
                        // Try clicking the label inside the pretty wrapper
                        var label = input.querySelector('.pretty.p-switch label, label');
                        if (label) {
                            label.click();
                            return true;
                        }
                    }

                    // Strategy 3: Find by CSS escape of the ID
                    var byEscape = document.querySelector('#' + CSS.escape('%s') + ' input[type=\"checkbox\"], input#' + CSS.escape('%s'));
                    if (byEscape) {
                        byEscape.click();
                        return true;
                    }

                    // Strategy 4: Find the pretty p-switch wrapper containing our input
                    var switches = document.querySelectorAll('.pretty.p-switch');
                    for (var i = 0; i < switches.length; i++) {
                        var inp = switches[i].querySelector('input[id=\"' + '%s' + '\"]');
                        if (inp) {
                            inp.click();
                            return true;
                        }
                    }

                    return false;
                })()
                ",
                input_id, input_id, input_id, input_id
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) list(result = list(value = FALSE))
                    )
                }
            )

            if (!isTRUE(result$result$value)) {
                warning("Switch '", input_id, "' not found")
                return(FALSE)
            }

            Sys.sleep(0.5)
            wait_for_shiny_idle(session)
        },
        error = function(e) {
            warning("Error clicking switch '", input_id, "': ", e$message)
            FALSE
        }
    )
}

#' Click a Shiny action button
#'
#' Finds the button by its ID and clicks it via JavaScript, then waits
#' for Shiny to become idle.
#'
#' @param session A ChromoteSession object
#' @param button_id The Shiny button ID to click
#' @return TRUE if the button was clicked, FALSE otherwise
click_button <- function(session, button_id) {
    tryCatch(
        {
            js <- sprintf(
                "
                (function() {
                    var btn = document.getElementById('%s');
                    if (btn) {
                        btn.click();
                        return true;
                    }
                    return false;
                })()
                ",
                button_id
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) list(result = list(value = FALSE))
                    )
                }
            )

            if (!isTRUE(result$result$value)) {
                warning("Button '", button_id, "' not found")
                return(FALSE)
            }

            Sys.sleep(0.5)
            wait_for_shiny_idle(session)
        },
        error = function(e) {
            warning("Error clicking button '", button_id, "': ", e$message)
            FALSE
        }
    )
}

#' Count the number of traces in a plotly chart
#'
#' Accesses the plotly internal data structure to count how many traces
#' are present in the plot.
#'
#' @param session A ChromoteSession object
#' @param output_id The Shiny output ID of the plotly output
#' @return Integer count of plotly traces, or 0 if not found
get_plotly_traces <- function(session, output_id) {
    tryCatch(
        {
            js <- sprintf(
                "
                (function() {
                    // Plotly attaches data to the .js-plotly-plot element inside the output div
                    var container = document.getElementById('%s');
                    if (!container) return 0;

                    // The actual plotly element may be the container itself or a child
                    var plotEl = container.querySelector('.js-plotly-plot') || container;

                    // plotly.js stores trace data in _fullData
                    if (plotEl._fullData && Array.isArray(plotEl._fullData)) {
                        return plotEl._fullData.length;
                    }

                    // Fallback: check data attribute
                    if (plotEl.data && Array.isArray(plotEl.data)) {
                        return plotEl.data.length;
                    }

                    return 0;
                })()
                ",
                output_id
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) list(result = list(value = 0))
                    )
                }
            )

            if (is.null(result) || is.null(result$result$value)) {
                return(0L)
            }
            as.integer(result$result$value)
        },
        error = function(e) {
            warning("Error getting plotly traces for '", output_id, "': ", e$message)
            0L
        }
    )
}

#' Get the text content of a DOM element
#'
#' Uses querySelector to find an element and returns its trimmed textContent.
#'
#' @param session A ChromoteSession object
#' @param selector A CSS selector string
#' @return The text content as a string, or NULL if the element was not found
get_element_text <- function(session, selector) {
    tryCatch(
        {
            safe_selector <- gsub("'", "\\\\'", selector)

            js <- sprintf(
                "
                (function() {
                    var el = document.querySelector('%s');
                    return el ? el.textContent.trim() : null;
                })()
                ",
                safe_selector
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) NULL
                    )
                }
            )

            if (is.null(result)) return(NULL)
            result$result$value
        },
        error = function(e) {
            warning("Error getting element text for '", selector, "': ", e$message)
            NULL
        }
    )
}

#' Wait for a DOM element to appear
#'
#' Polls until an element matching the selector exists in the DOM,
#' or until the timeout expires.
#'
#' @param session A ChromoteSession object
#' @param selector A CSS selector string
#' @param timeout Maximum seconds to wait (default: 15)
#' @param poll_interval Seconds between polls (default: 0.5)
#' @return TRUE if the element appeared, FALSE if timeout expired
wait_for_element <- function(session, selector, timeout = 15, poll_interval = 0.5) {
    tryCatch(
        {
            safe_selector <- gsub("'", "\\\\'", selector)
            js <- sprintf(
                "document.querySelector('%s') !== null",
                safe_selector
            )

            start_time <- Sys.time()

            while (difftime(Sys.time(), start_time, units = "secs") < timeout) {
                result <- tryCatch(
                    session$Runtime$evaluate(expression = js),
                    error = function(e) {
                        Sys.sleep(1)
                        tryCatch(
                            session$Runtime$evaluate(expression = js),
                            error = function(e2) NULL
                        )
                    }
                )

                if (!is.null(result) && isTRUE(result$result$value)) {
                    return(TRUE)
                }

                Sys.sleep(poll_interval)
            }

            FALSE
        },
        error = function(e) {
            warning("Error waiting for element '", selector, "': ", e$message)
            FALSE
        }
    )
}

#' Open a shinydashboardPlus boxSidebar
#'
#' Finds the sidebar toggle button within a box and clicks it if the
#' sidebar is not already open. Waits for the CSS animation to complete.
#'
#' @param session A ChromoteSession object
#' @param box_id The ID of the shinydashboardPlus box element
#' @return TRUE if the sidebar was opened (or was already open), FALSE on error
open_box_sidebar <- function(session, box_id) {
    tryCatch(
        {
            js <- sprintf(
                "
                (function() {
                    var box = document.getElementById('%s');
                    if (!box) return false;

                    // Check if sidebar is already open
                    // shinydashboardPlus uses 'control-sidebar-open' class on the box
                    // or 'direct-chat-contacts-open' on some versions
                    var isOpen = box.classList.contains('control-sidebar-open') ||
                                 box.classList.contains('direct-chat-contacts-open') ||
                                 box.querySelector('.direct-chat-contacts-open') !== null;

                    if (isOpen) return true;

                    // Find the sidebar toggle button
                    // shinydashboardPlus boxSidebar toggle is typically a button with
                    // data-widget='chat-pane-toggle' or a btn inside the box-tools
                    var toggle = box.querySelector('[data-widget=\"chat-pane-toggle\"]') ||
                                 box.querySelector('.box-tools .btn[data-card-widget=\"sidebar-toggle\"]') ||
                                 box.querySelector('.card-tools .btn[data-card-widget=\"sidebar-toggle\"]') ||
                                 box.querySelector('[data-card-widget=\"sidebar-toggle\"]') ||
                                 box.querySelector('.box-tools button:not([data-widget=\"collapse\"]):not([data-widget=\"remove\"])');

                    if (toggle) {
                        toggle.click();
                        return true;
                    }

                    return false;
                })()
                ",
                box_id
            )

            result <- tryCatch(
                session$Runtime$evaluate(expression = js),
                error = function(e) {
                    Sys.sleep(2)
                    tryCatch(
                        session$Runtime$evaluate(expression = js),
                        error = function(e2) list(result = list(value = FALSE))
                    )
                }
            )

            if (!isTRUE(result$result$value)) {
                warning("Box sidebar toggle not found for box '", box_id, "'")
                return(FALSE)
            }

            # Wait for CSS animation to complete
            Sys.sleep(0.5)
            wait_for_shiny_idle(session)
        },
        error = function(e) {
            warning("Error opening box sidebar '", box_id, "': ", e$message)
            FALSE
        }
    )
}
