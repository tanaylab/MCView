# daf_contracts_docs.R - Pretty-print + markdown generation for contracts
#
# Split from R/daf_contracts.R (2026-05-01). The companion contract
# definition / validator files live in R/daf_contracts_definitions.R and
# R/daf_contracts_validators.R.

#' Print contract summary
#'
#' Prints a human-readable summary of a contract.
#'
#' @param contract Contract to summarize
#' @export
print_contract <- function(contract) {
    cli::cli_h1(contract$name)
    cli::cli_text(contract$description)

    if (!is.null(contract$extends)) {
        cli::cli_text("Extends: {contract$extends}")
    }

    if (length(contract$axes) > 0) {
        cli::cli_h2("Axes")
        for (name in names(contract$axes)) {
            spec <- contract$axes[[name]]
            cli::cli_li("{name} ({spec$expectation}): {spec$description}")
        }
    }

    vectors <- c(contract$vectors %||% list(), contract$graph_vectors %||% list())
    if (length(vectors) > 0) {
        cli::cli_h2("Vectors")
        for (spec in vectors) {
            cli::cli_li("{spec$axis}.{spec$name} ({spec$expectation}, {spec$type}): {spec$description}")
        }
    }

    if (length(contract$matrices) > 0) {
        cli::cli_h2("Matrices")
        for (spec in contract$matrices) {
            cli::cli_li("{spec$rows_axis},{spec$cols_axis}.{spec$name} ({spec$expectation}): {spec$description}")
        }
    }

    if (length(contract$scalars) > 0) {
        cli::cli_h2("Scalars")
        for (spec in contract$scalars) {
            cli::cli_li("{spec$name} ({spec$expectation}, {spec$type}): {spec$description}")
        }
    }

    if (!is.null(contract$custom_rules) && length(contract$custom_rules) > 0) {
        cli::cli_h2("Custom Rules")
        for (name in names(contract$custom_rules)) {
            rule <- contract$custom_rules[[name]]
            cli::cli_li("{name}: {rule$description}")
        }
    }
}

#' Generate contract documentation as markdown
#'
#' @param contract Contract to document
#' @return Character string with markdown content
#' @export
contract_to_markdown <- function(contract) {
    lines <- character()

    lines <- c(lines, glue::glue("# {contract$name}"))
    lines <- c(lines, "")
    lines <- c(lines, contract$description)
    lines <- c(lines, "")

    if (!is.null(contract$extends)) {
        lines <- c(lines, glue::glue("**Extends:** {contract$extends}"))
        lines <- c(lines, "")
    }

    if (length(contract$axes) > 0) {
        lines <- c(lines, "## Axes")
        lines <- c(lines, "")
        lines <- c(lines, "| Axis | Required | Description |")
        lines <- c(lines, "|------|----------|-------------|")
        for (name in names(contract$axes)) {
            spec <- contract$axes[[name]]
            req <- if (spec$expectation == CONTRACT_REQUIRED) "Yes" else "No"
            lines <- c(lines, glue::glue("| `{name}` | {req} | {spec$description} |"))
        }
        lines <- c(lines, "")
    }

    vectors <- c(contract$vectors %||% list(), contract$graph_vectors %||% list())
    if (length(vectors) > 0) {
        lines <- c(lines, "## Vectors")
        lines <- c(lines, "")
        lines <- c(lines, "| Axis | Name | Required | Type | Description |")
        lines <- c(lines, "|------|------|----------|------|-------------|")
        for (spec in vectors) {
            req <- if (spec$expectation == CONTRACT_REQUIRED) "Yes" else "No"
            lines <- c(lines, glue::glue("| `{spec$axis}` | `{spec$name}` | {req} | {spec$type} | {spec$description} |"))
        }
        lines <- c(lines, "")
    }

    if (length(contract$matrices) > 0) {
        lines <- c(lines, "## Matrices")
        lines <- c(lines, "")
        lines <- c(lines, "| Rows Axis | Cols Axis | Name | Required | Description |")
        lines <- c(lines, "|-----------|-----------|------|----------|-------------|")
        for (spec in contract$matrices) {
            req <- if (spec$expectation == CONTRACT_REQUIRED) "Yes" else "No"
            lines <- c(lines, glue::glue("| `{spec$rows_axis}` | `{spec$cols_axis}` | `{spec$name}` | {req} | {spec$description} |"))
        }
        lines <- c(lines, "")
    }

    if (length(contract$scalars) > 0) {
        lines <- c(lines, "## Scalars")
        lines <- c(lines, "")
        lines <- c(lines, "| Name | Required | Type | Description |")
        lines <- c(lines, "|------|----------|------|-------------|")
        for (spec in contract$scalars) {
            req <- if (spec$expectation == CONTRACT_REQUIRED) "Yes" else "No"
            lines <- c(lines, glue::glue("| `{spec$name}` | {req} | {spec$type} | {spec$description} |"))
        }
        lines <- c(lines, "")
    }

    paste(lines, collapse = "\n")
}

#' Generate full MCView contract documentation
#'
#' @param output_file Path to output markdown file (optional)
#' @return Character string with full documentation
#' @export
generate_contract_docs <- function(output_file = NULL) {
    contracts <- mcview_all_contracts()

    lines <- character()
    lines <- c(lines, "# MCView DAF Contract Specification")
    lines <- c(lines, "")
    lines <- c(lines, "This document specifies the data contract for MCView DAF files.")
    lines <- c(lines, "")
    lines <- c(lines, "## Overview")
    lines <- c(lines, "")
    lines <- c(lines, "MCView requires data to be stored in DAF (Data Axes Format) files.")
    lines <- c(lines, "The contract defines required and optional data for each feature.")
    lines <- c(lines, "")

    for (name in names(contracts)) {
        contract <- contracts[[name]]
        lines <- c(lines, contract_to_markdown(contract))
        lines <- c(lines, "---")
        lines <- c(lines, "")
    }

    result <- paste(lines, collapse = "\n")

    if (!is.null(output_file)) {
        writeLines(result, output_file)
        cli::cli_alert_success("Contract documentation written to {output_file}")
    }

    invisible(result)
}
