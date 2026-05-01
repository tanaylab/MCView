# daf_contracts_validators.R - Contract validation against a DAF object
#
# Split from R/daf_contracts.R (2026-05-01). See R/daf_contracts_definitions.R
# for the contract definitions being validated.

# ==============================================================================
# Contract Validation Functions
# ==============================================================================

#' Validate DAF against MCView contract
#'
#' Validates that a DAF object meets the MCView contract requirements.
#'
#' @param daf_obj DAF object to validate
#' @param contract Contract to validate against (default: core contract)
#' @param strict If TRUE, fail on any missing optional data (default: FALSE)
#' @param verbose If TRUE, print detailed validation info
#'
#' @return List with `valid` (logical), `errors` (character vector),
#'         `warnings` (character vector), and `info` (character vector)
#' @export
validate_mcview_contract <- function(daf_obj,
                                     contract = mcview_core_contract(),
                                     strict = FALSE,
                                     verbose = FALSE) {
    result <- list(
        valid = TRUE,
        errors = character(),
        warnings = character(),
        info = character()
    )

    if (!dafr::is_daf(daf_obj)) {
        result$valid <- FALSE
        result$errors <- c(result$errors, "Object is not a valid DAF object")
        return(result)
    }

    daf_contract <- mcview_contract_to_dafr(contract)

    # native dafr's verify_contract(contract, daf) stops on violation and
    # returns `daf` invisibly on success. Capture any error into result$errors
    # instead of letting it propagate.
    tryCatch(
        dafr::verify_contract(daf_contract, daf_obj),
        error = function(e) {
            result$valid <<- FALSE
            result$errors <<- c(result$errors, conditionMessage(e))
        }
    )

    if (!is.null(contract$axes)) {
        for (axis_name in names(contract$axes)) {
            spec <- contract$axes[[axis_name]]
            has_it <- dafr::has_axis(daf_obj, axis_name)

            if (spec$expectation == CONTRACT_OPTIONAL && !has_it && strict) {
                result$warnings <- c(
                    result$warnings,
                    glue::glue("Missing optional axis: {axis_name}")
                )
            } else if (has_it && verbose) {
                n <- dafr::axis_length(daf_obj, axis_name)
                result$info <- c(
                    result$info,
                    glue::glue("Axis '{axis_name}': {n} entries")
                )
            }
        }
    }

    all_vectors <- c(
        contract$vectors %||% list(),
        contract$graph_vectors %||% list()
    )

    for (spec in all_vectors) {
        has_it <- dafr::has_axis(daf_obj, spec$axis) &&
            dafr::has_vector(daf_obj, spec$axis, spec$name)

        if (spec$expectation == CONTRACT_OPTIONAL && !has_it && strict) {
            result$warnings <- c(
                result$warnings,
                glue::glue("Missing optional vector: {spec$axis}.{spec$name}")
            )
        } else if (has_it && verbose) {
            result$info <- c(
                result$info,
                glue::glue("Vector '{spec$axis}.{spec$name}': present")
            )
        }
    }

    if (!is.null(contract$matrices)) {
        for (spec in contract$matrices) {
            has_it <- dafr::has_axis(daf_obj, spec$rows_axis) &&
                dafr::has_axis(daf_obj, spec$cols_axis) &&
                dafr::has_matrix(daf_obj, spec$rows_axis, spec$cols_axis, spec$name)

            if (spec$expectation == CONTRACT_OPTIONAL && !has_it && strict) {
                result$warnings <- c(
                    result$warnings,
                    glue::glue("Missing optional matrix: {spec$rows_axis},{spec$cols_axis}.{spec$name}")
                )
            } else if (has_it && verbose) {
                result$info <- c(
                    result$info,
                    glue::glue("Matrix '{spec$rows_axis},{spec$cols_axis}.{spec$name}': present")
                )
            }
        }
    }

    if (!is.null(contract$scalars)) {
        for (spec in contract$scalars) {
            has_it <- dafr::has_scalar(daf_obj, spec$name)

            if (spec$expectation == CONTRACT_OPTIONAL && !has_it && strict) {
                result$warnings <- c(
                    result$warnings,
                    glue::glue("Missing optional scalar: {spec$name}")
                )
            } else if (has_it && verbose) {
                result$info <- c(
                    result$info,
                    glue::glue("Scalar '{spec$name}': present")
                )
            }
        }
    }

    # Run custom validation rules
    if (!is.null(contract$custom_rules)) {
        for (rule_name in names(contract$custom_rules)) {
            rule <- contract$custom_rules[[rule_name]]
            error_msg <- rule$validate(daf_obj)
            if (!is.null(error_msg)) {
                result$valid <- FALSE
                result$errors <- c(result$errors, error_msg)
            } else if (verbose) {
                result$info <- c(
                    result$info,
                    glue::glue("Custom rule '{rule_name}': passed")
                )
            }
        }
    }

    return(result)
}

#' Validate DAF for all tabs
#'
#' Validates DAF against all tab-specific contracts and reports which tabs
#' can be enabled.
#'
#' @param daf_obj DAF object to validate
#' @param verbose If TRUE, print detailed validation info
#'
#' @return Named list with tab names as keys and validation results as values
#' @export
validate_mcview_tabs <- function(daf_obj, verbose = FALSE) {
    tabs <- MCVIEW_TAB_NAMES

    results <- list()

    # First validate core contract
    core_result <- validate_mcview_contract(daf_obj, mcview_core_contract(), verbose = verbose)

    if (!core_result$valid) {
        if (verbose) {
            cli::cli_alert_danger("Core contract validation failed:")
            for (err in core_result$errors) {
                cli::cli_alert_warning("  {err}")
            }
        }
        # If core fails, no tabs are available
        for (tab in tabs) {
            results[[tab]] <- list(available = FALSE, reason = "Core contract failed")
        }
        return(results)
    }

    # Validate each tab
    for (tab in tabs) {
        tab_contract <- mcview_tab_contract(tab)
        if (is.null(tab_contract)) {
            results[[tab]] <- list(available = TRUE, reason = "No specific contract")
            next
        }

        tab_result <- validate_mcview_contract(daf_obj, tab_contract, verbose = FALSE)

        if (tab_result$valid) {
            results[[tab]] <- list(available = TRUE, reason = "Contract satisfied")
            if (verbose) {
                cli::cli_alert_success("Tab '{tab}': available")
            }
        } else {
            # Check if it's just missing optional data for extended features
            required_errors <- grep("required", tab_result$errors, value = TRUE, ignore.case = TRUE)
            if (length(required_errors) == 0) {
                results[[tab]] <- list(available = TRUE, reason = "Optional data missing")
            } else {
                results[[tab]] <- list(
                    available = FALSE,
                    reason = paste(required_errors, collapse = "; ")
                )
                if (verbose) {
                    cli::cli_alert_warning("Tab '{tab}': unavailable - {results[[tab]]$reason}")
                }
            }
        }
    }

    return(results)
}
