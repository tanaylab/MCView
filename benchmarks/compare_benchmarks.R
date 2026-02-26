#!/usr/bin/env Rscript
# ==============================================================================
# MCView DAF Benchmark Comparison
# ==============================================================================
#
# Reads two benchmark result JSON files and shows a before/after comparison
# table with speedup calculations.
#
# Usage:
#   Rscript benchmarks/compare_benchmarks.R <before.json> <after.json>
#   Rscript benchmarks/compare_benchmarks.R benchmarks/results/bench_*.json  # latest two
#
# Output:
#   Formatted table printed to stdout showing:
#   - Operation name
#   - Before median (s)
#   - After median (s)
#   - Speedup factor (before/after)
#   - Absolute improvement (seconds saved)
#   - Relative improvement (percentage)
# ==============================================================================

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) < 2) {
        # Try to find the two most recent result files
        result_dir <- file.path(dirname(dirname(
            tryCatch({
                cmd_args <- commandArgs(trailingOnly = FALSE)
                file_arg <- grep("^--file=", cmd_args, value = TRUE)
                if (length(file_arg) > 0) {
                    normalizePath(sub("^--file=", "", file_arg[1]))
                } else {
                    file.path(getwd(), "benchmarks", "compare_benchmarks.R")
                }
            }, error = function(e) file.path(getwd(), "benchmarks", "compare_benchmarks.R"))
        )), "results")

        json_files <- sort(list.files(result_dir, pattern = "\\.json$", full.names = TRUE))

        if (length(json_files) >= 2) {
            args <- tail(json_files, 2)
            cat(sprintf("Auto-detected result files:\n  Before: %s\n  After:  %s\n\n", args[1], args[2]))
        } else {
            cat("Usage: Rscript compare_benchmarks.R <before.json> <after.json>\n")
            cat("\nAlternatively, provide no arguments and the two most recent result\n")
            cat("files in benchmarks/results/ will be used.\n")
            cat(sprintf("\nFound %d result file(s) in %s\n", length(json_files), result_dir))
            quit(status = 1)
        }
    }

    before_path <- args[1]
    after_path <- args[2]

    if (!file.exists(before_path)) stop("File not found: ", before_path)
    if (!file.exists(after_path)) stop("File not found: ", after_path)

    before <- jsonlite::read_json(before_path, simplifyVector = FALSE)
    after <- jsonlite::read_json(after_path, simplifyVector = FALSE)

    # -- Header ----------------------------------------------------------------
    cat("==============================================================================\n")
    cat("MCView DAF Benchmark Comparison\n")
    cat("==============================================================================\n\n")

    cat("Before:\n")
    print_metadata(before$metadata, before_path)
    cat("\nAfter:\n")
    print_metadata(after$metadata, after_path)

    # -- Init timings comparison -----------------------------------------------
    cat("\n--- Initialization Timings ---\n\n")
    cat(sprintf("%-25s %10s %10s %10s %10s\n",
                "Phase", "Before(s)", "After(s)", "Speedup", "Saved(s)"))
    cat(paste(rep("-", 70), collapse = ""), "\n")

    init_fields <- c("pkg_load_s", "julia_init_s", "daf_open_s", "mcview_init_s", "total_init_s")
    init_labels <- c("Package load", "Julia init", "DAF open", "MCView init", "Total init")

    for (i in seq_along(init_fields)) {
        field <- init_fields[i]
        b_val <- before$init_timings[[field]] %||% NA_real_
        a_val <- after$init_timings[[field]] %||% NA_real_

        if (!is.na(b_val) && !is.na(a_val) && a_val > 0) {
            speedup <- b_val / a_val
            saved <- b_val - a_val
            cat(sprintf("%-25s %10.2f %10.2f %9.2fx %10.2f\n",
                        init_labels[i], b_val, a_val, speedup, saved))
        } else {
            cat(sprintf("%-25s %10s %10s %10s %10s\n",
                        init_labels[i],
                        if (is.na(b_val)) "N/A" else sprintf("%.2f", b_val),
                        if (is.na(a_val)) "N/A" else sprintf("%.2f", a_val),
                        "N/A", "N/A"))
        }
    }

    # -- Benchmark comparison --------------------------------------------------
    cat("\n--- Benchmark Timings ---\n\n")
    cat(sprintf("%-35s %10s %10s %10s %10s %8s\n",
                "Operation", "Before(s)", "After(s)", "Speedup", "Saved(s)", "Saved%"))
    cat(paste(rep("-", 90), collapse = ""), "\n")

    # Get union of all benchmark names
    all_names <- union(names(before$benchmarks), names(after$benchmarks))

    # Sort by before median time (descending) so slowest ops are at top
    before_medians <- sapply(all_names, function(n) {
        b <- before$benchmarks[[n]]
        if (is.null(b)) return(0)
        b$median_s %||% 0
    })
    all_names <- all_names[order(-before_medians)]

    total_before <- 0
    total_after <- 0
    n_compared <- 0

    for (bname in all_names) {
        b <- before$benchmarks[[bname]]
        a <- after$benchmarks[[bname]]

        b_median <- if (!is.null(b)) b$median_s %||% NA_real_ else NA_real_
        a_median <- if (!is.null(a)) a$median_s %||% NA_real_ else NA_real_

        if (!is.na(b_median) && !is.na(a_median) && a_median > 0) {
            speedup <- b_median / a_median
            saved <- b_median - a_median
            saved_pct <- (1 - a_median / b_median) * 100

            color_prefix <- ""
            color_suffix <- ""

            cat(sprintf("%-35s %10.4f %10.4f %9.2fx %10.4f %7.1f%%\n",
                        bname, b_median, a_median, speedup, saved, saved_pct))

            total_before <- total_before + b_median
            total_after <- total_after + a_median
            n_compared <- n_compared + 1
        } else {
            cat(sprintf("%-35s %10s %10s %10s %10s %8s\n",
                        bname,
                        if (is.na(b_median)) "N/A" else sprintf("%.4f", b_median),
                        if (is.na(a_median)) "N/A" else sprintf("%.4f", a_median),
                        "N/A", "N/A", "N/A"))
        }
    }

    # Summary row
    if (n_compared > 0 && total_after > 0) {
        cat(paste(rep("-", 90), collapse = ""), "\n")
        overall_speedup <- total_before / total_after
        overall_saved <- total_before - total_after
        overall_saved_pct <- (1 - total_after / total_before) * 100
        cat(sprintf("%-35s %10.4f %10.4f %9.2fx %10.4f %7.1f%%\n",
                    "TOTAL (sum of medians)", total_before, total_after,
                    overall_speedup, overall_saved, overall_saved_pct))
    }

    # -- Dataset dimension comparison ------------------------------------------
    cat("\n--- Dataset Dimensions ---\n")
    cat(sprintf("  Before: %d genes x %d metacells\n",
                before$metadata$n_genes %||% 0,
                before$metadata$n_metacells %||% 0))
    cat(sprintf("  After:  %d genes x %d metacells\n",
                after$metadata$n_genes %||% 0,
                after$metadata$n_metacells %||% 0))

    if ((before$metadata$n_genes %||% 0) != (after$metadata$n_genes %||% 0) ||
        (before$metadata$n_metacells %||% 0) != (after$metadata$n_metacells %||% 0)) {
        cat("  WARNING: Dataset dimensions differ between runs. Comparison may not be meaningful.\n")
    }

    cat("\n")
}

print_metadata <- function(md, path) {
    cat(sprintf("  File:       %s\n", basename(path)))
    cat(sprintf("  Timestamp:  %s\n", md$timestamp %||% "unknown"))
    cat(sprintf("  Host:       %s\n", md$hostname %||% "unknown"))
    cat(sprintf("  Dimensions: %d genes x %d metacells\n",
                md$n_genes %||% 0, md$n_metacells %||% 0))
    cat(sprintf("  Iterations: %d\n", md$iterations %||% 0))
    cat(sprintf("  DAF path:   %s\n", md$daf_path %||% "unknown"))
}

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1 && is.na(x))) y else x

# Run
main()
