#!/usr/bin/env Rscript
# ==============================================================================
# MCView DAF Benchmark Suite
# ==============================================================================
#
# Measures timing of key DAF operations and MCView workflows.
#
# Usage:
#   Rscript benchmarks/benchmark_daf.R --daf /path/to/daf
#   Rscript benchmarks/benchmark_daf.R --daf /path/to/daf --iterations 5
#   Rscript benchmarks/benchmark_daf.R --daf /path/to/daf --output benchmarks/results/my_run.json
#   Rscript benchmarks/benchmark_daf.R --daf /path/to/daf --benchmarks "single_gene_egc,full_egc_matrix"
#
# Requires:
#   - conda environment dafr-mcview activated
#   - Julia env vars set (JULIA_PROJECT, JULIA_LOAD_PATH, JULIA_DEPOT_PATH)
#   - dafr.JULIA_HOME option set
#
# Output:
#   JSON file in benchmarks/results/ with timing data for each operation.
# ==============================================================================

# -- Parse command-line arguments ----------------------------------------------

parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    opts <- list(
        daf = NULL,
        output = NULL,
        iterations = 5L,
        dataset_name = "bench",
        pkg_root = NULL,
        benchmarks = NULL  # NULL means run all; otherwise comma-separated names
    )

    i <- 1L
    while (i <= length(args)) {
        if (args[i] == "--daf" && i < length(args)) {
            opts$daf <- args[i + 1L]
            i <- i + 2L
        } else if (args[i] == "--output" && i < length(args)) {
            opts$output <- args[i + 1L]
            i <- i + 2L
        } else if (args[i] == "--iterations" && i < length(args)) {
            opts$iterations <- as.integer(args[i + 1L])
            i <- i + 2L
        } else if (args[i] == "--dataset-name" && i < length(args)) {
            opts$dataset_name <- args[i + 1L]
            i <- i + 2L
        } else if (args[i] == "--pkg-root" && i < length(args)) {
            opts$pkg_root <- args[i + 1L]
            i <- i + 2L
        } else if (args[i] == "--benchmarks" && i < length(args)) {
            opts$benchmarks <- trimws(strsplit(args[i + 1L], ",")[[1]])
            i <- i + 2L
        } else if (args[i] == "--help") {
            cat("Usage: Rscript benchmark_daf.R --daf <path> [options]\n")
            cat("Options:\n")
            cat("  --daf <path>           Path to DAF directory or H5 file (required)\n")
            cat("  --output <path>        Output JSON file (default: auto-generated in benchmarks/results/)\n")
            cat("  --iterations <n>       Number of timed iterations per operation (default: 5)\n")
            cat("  --dataset-name <name>  Dataset name for MCView init (default: bench)\n")
            cat("  --pkg-root <path>      Package root for devtools::load_all (default: auto-detect)\n")
            cat("  --benchmarks <list>    Comma-separated list of benchmarks to run (default: all)\n")
            cat("  --help                 Show this help message\n")
            cat("\nAvailable benchmarks:\n")
            cat("  single_gene_egc, two_gene_scatter, full_egc_matrix, calc_top_cors,\n")
            cat("  gene_gene_correlations, gg_mc_top_cor, marker_genes, mc_ordering,\n")
            cat("  gene_modules, diff_expr, umap, daf_single_gene_vec, daf_full_matrix,\n")
            cat("  daf_gene_aggregations\n")
            quit(status = 0)
        } else {
            stop("Unknown argument: ", args[i], "\nUse --help for usage information.")
        }
    }

    if (is.null(opts$daf)) {
        stop("--daf argument is required. Use --help for usage information.")
    }

    opts
}

# -- Timing utilities ----------------------------------------------------------

#' Run a single benchmark with warmup and multiple iterations
#'
#' @param name Human-readable name for the benchmark
#' @param expr_fn A zero-argument function to benchmark
#' @param iterations Number of timed iterations
#' @param warmup Number of warmup iterations (untimed)
#' @return List with name, timings, and summary statistics
run_benchmark <- function(name, expr_fn, iterations = 5L, warmup = 1L) {
    cat(sprintf("  [%s] warming up (%d)...", name, warmup))
    for (w in seq_len(warmup)) {
        tryCatch(expr_fn(), error = function(e) {
            cat(sprintf(" WARMUP ERROR: %s", conditionMessage(e)))
        })
    }
    cat(" done.\n")

    cat(sprintf("  [%s] timing (%d iterations)...", name, iterations))
    timings <- numeric(iterations)
    errors <- character(0)

    for (i in seq_len(iterations)) {
        gc(verbose = FALSE)  # Clean GC before each iteration
        t0 <- proc.time()[["elapsed"]]
        tryCatch(
            expr_fn(),
            error = function(e) {
                errors <<- c(errors, conditionMessage(e))
            }
        )
        timings[i] <- proc.time()[["elapsed"]] - t0
    }

    if (length(errors) > 0) {
        cat(sprintf(" ERRORS: %d/%d\n", length(errors), iterations))
    } else {
        cat(sprintf(" median=%.4fs\n", median(timings)))
    }

    list(
        name = name,
        timings_s = timings,
        median_s = median(timings),
        mean_s = mean(timings),
        min_s = min(timings),
        max_s = max(timings),
        sd_s = sd(timings),
        iterations = iterations,
        errors = if (length(errors) > 0) errors else NULL
    )
}

# -- Pick representative genes from the dataset --------------------------------

pick_test_genes <- function(daf_obj) {
    gene_names <- dafr::axis_entries(daf_obj, "gene")
    n_genes <- length(gene_names)

    # Try to pick biologically interesting genes if they exist
    candidates <- c(
        # Common housekeeping / high-expression genes
        "Actb", "Gapdh", "Rpl13a", "Hbb-bs", "Tmsb4x",
        # Common marker genes
        "Cd79a", "Cd3e", "Lyz2", "Nkg7", "Mki67",
        # Typical developmental genes
        "Sox2", "Pax6", "Epcam", "Col1a1", "Hba-a1"
    )
    available <- intersect(candidates, gene_names)

    if (length(available) >= 4) {
        gene1 <- available[1]
        gene2 <- available[2]
        gene_set <- available[1:min(4, length(available))]
        gene_large_set <- available
    } else {
        # Fallback: pick genes by position (spread across the alphabet)
        idx <- round(seq(1, n_genes, length.out = 15))
        available <- gene_names[idx]
        gene1 <- available[1]
        gene2 <- available[2]
        gene_set <- available[1:4]
        gene_large_set <- available
    }

    list(
        gene1 = gene1,
        gene2 = gene2,
        gene_set_small = gene_set,
        gene_set_large = gene_large_set,
        all_gene_names = gene_names,
        n_genes = n_genes
    )
}

# -- Pick representative metacells from the dataset ----------------------------

pick_test_metacells <- function(daf_obj) {
    mc_names <- dafr::axis_entries(daf_obj, "metacell")
    n_mc <- length(mc_names)

    list(
        mc1 = mc_names[1],
        mc2 = mc_names[min(2, n_mc)],
        mc_pair = mc_names[1:min(2, n_mc)],
        all_metacell_names = mc_names,
        n_metacells = n_mc
    )
}

# -- Define all benchmark operations ------------------------------------------

define_benchmarks <- function(dataset_name, daf_obj, test_genes, test_mcs, iterations) {

    benchmarks <- list()

    # ---- 1. Single gene EGC extraction (get_gene_egc) ----
    benchmarks$single_gene_egc <- function() {
        run_benchmark(
            "single_gene_egc",
            function() get_gene_egc(test_genes$gene1, dataset_name),
            iterations = iterations
        )
    }

    # ---- 2. Two gene scatter data (plot_gg_over_mc data extraction) ----
    benchmarks$two_gene_scatter <- function() {
        run_benchmark(
            "two_gene_scatter",
            function() {
                egc1 <- get_gene_egc(test_genes$gene1, dataset_name)
                egc2 <- get_gene_egc(test_genes$gene2, dataset_name)
                # Simulate what plot_gg_over_mc does: get 2 genes + combine
                data.frame(g1 = egc1, g2 = egc2)
            },
            iterations = iterations
        )
    }

    # ---- 3. Full EGC matrix retrieval (get_mc_egc) ----
    benchmarks$full_egc_matrix <- function() {
        run_benchmark(
            "full_egc_matrix",
            function() get_mc_egc(dataset_name),
            iterations = iterations
        )
    }

    # ---- 4. Gene-gene correlation for one gene vs all (calc_top_cors) ----
    benchmarks$calc_top_cors <- function() {
        run_benchmark(
            "calc_top_cors",
            function() calc_top_cors(dataset_name, test_genes$gene1, "both",
                                     NULL, NULL, NULL),
            iterations = iterations
        )
    }

    # ---- 5. Gene-gene correlation for gene subset (calc_gene_gene_correlations) ----
    benchmarks$gene_gene_correlations <- function() {
        run_benchmark(
            "gene_gene_correlations",
            function() calc_gene_gene_correlations(dataset_name, test_genes$gene_set_small),
            iterations = iterations
        )
    }

    # ---- 6. Full gene-gene top correlation matrix (calc_gg_mc_top_cor) ----
    benchmarks$gg_mc_top_cor <- function() {
        run_benchmark(
            "gg_mc_top_cor",
            function() {
                egc <- get_mc_egc(dataset_name)
                calc_gg_mc_top_cor(egc, k = 30, daf_obj = daf_obj)
            },
            iterations = max(1L, iterations %/% 2L)  # Fewer iterations - this is expensive
        )
    }

    # ---- 7. Marker gene computation (calc_marker_genes) ----
    benchmarks$marker_genes <- function() {
        run_benchmark(
            "marker_genes",
            function() {
                egc <- get_mc_egc(dataset_name)
                calc_marker_genes(egc, genes_per_metacell = 20, daf_obj = daf_obj)
            },
            iterations = iterations
        )
    }

    # ---- 8. Metacell ordering for heatmap (order_mc_by_most_var_genes) ----
    benchmarks$mc_ordering <- function() {
        run_benchmark(
            "mc_ordering",
            function() {
                mc_fp <- get_mc_fp(dataset_name)
                order_mc_by_most_var_genes(mc_fp, notify_var_genes = FALSE)
            },
            iterations = iterations
        )
    }

    # ---- 9. Gene module computation (calc_gene_modules) ----
    benchmarks$gene_modules <- function() {
        run_benchmark(
            "gene_modules",
            function() {
                mc_mat <- daf_query_mc_mat(daf_obj)
                calc_gene_modules(mc_mat, verbose = FALSE)
            },
            iterations = max(1L, iterations %/% 2L)  # Fewer iterations - this is expensive
        )
    }

    # ---- 10. Differential expression (calc_diff_expr) ----
    benchmarks$diff_expr <- function() {
        run_benchmark(
            "diff_expr",
            function() {
                calc_mc_mc_gene_df(dataset_name, test_mcs$mc1, test_mcs$mc2)
            },
            iterations = iterations
        )
    }

    # ---- 11. UMAP computation (compute_umap) ----
    benchmarks$umap <- function() {
        # Pick 2-3 anchor genes for UMAP (needs >= 2)
        umap_anchors <- test_genes$gene_set_small[1:min(3, length(test_genes$gene_set_small))]
        run_benchmark(
            "umap",
            function() {
                egc <- get_mc_egc(dataset_name)
                compute_umap(egc, anchors = umap_anchors, n_epoch = 50)
            },
            iterations = max(1L, iterations %/% 2L)  # Fewer iterations - UMAP is expensive
        )
    }

    # ---- 12. DAF query for single gene vector ----
    benchmarks$daf_single_gene_vec <- function() {
        gene <- test_genes$gene1
        escaped_gene <- gsub("([^A-Za-z0-9_+\\-.])", "\\\\\\1", gene, perl = TRUE)
        run_benchmark(
            "daf_single_gene_vec",
            function() {
                query <- paste0("/ metacell / gene = ", escaped_gene, " : UMIs")
                daf_obj[query]
            },
            iterations = iterations * 2L  # More iterations for fast ops
        )
    }

    # ---- 13. DAF query for full matrix ----
    benchmarks$daf_full_matrix <- function() {
        run_benchmark(
            "daf_full_matrix",
            function() {
                dafr::get_matrix(daf_obj, "gene", "metacell", "UMIs")
            },
            iterations = iterations
        )
    }

    # ---- 14. DAF query for aggregated values (gene max/mean/sum UMIs) ----
    benchmarks$daf_gene_aggregations <- function() {
        run_benchmark(
            "daf_gene_aggregations",
            function() {
                max_umis <- daf_query_gene_max_umis(daf_obj)
                mean_umis <- daf_query_gene_mean_umis(daf_obj)
                sum_umis <- daf_query_gene_sum_umis(daf_obj)
                list(max = max_umis, mean = mean_umis, sum = sum_umis)
            },
            iterations = iterations
        )
    }

    benchmarks
}

# -- Main entry point ----------------------------------------------------------

main <- function() {
    opts <- parse_args()
    cat("==============================================================================\n")
    cat("MCView DAF Benchmark Suite\n")
    cat("==============================================================================\n")
    cat(sprintf("DAF path:    %s\n", opts$daf))
    cat(sprintf("Iterations:  %d\n", opts$iterations))
    cat(sprintf("Dataset:     %s\n", opts$dataset_name))
    cat(sprintf("Time:        %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    cat("==============================================================================\n\n")

    # -- Phase 0: Load MCView package ------------------------------------------
    cat("Phase 0: Loading MCView package...\n")
    pkg_root <- opts$pkg_root
    if (is.null(pkg_root)) {
        # Auto-detect: look for DESCRIPTION relative to this script
        script_dir <- tryCatch({
            # When run via Rscript, find the script's directory
            cmd_args <- commandArgs(trailingOnly = FALSE)
            file_arg <- grep("^--file=", cmd_args, value = TRUE)
            if (length(file_arg) > 0) {
                dirname(normalizePath(sub("^--file=", "", file_arg[1])))
            } else {
                getwd()
            }
        }, error = function(e) getwd())

        # Try parent of benchmarks/
        candidate <- dirname(script_dir)
        if (file.exists(file.path(candidate, "DESCRIPTION"))) {
            pkg_root <- candidate
        } else if (file.exists(file.path(getwd(), "DESCRIPTION"))) {
            pkg_root <- getwd()
        } else {
            stop("Cannot auto-detect package root. Use --pkg-root or run from package directory.")
        }
    }

    t0_load <- proc.time()[["elapsed"]]
    suppressMessages(devtools::load_all(pkg_root, quiet = TRUE))
    t_pkg_load <- proc.time()[["elapsed"]] - t0_load
    cat(sprintf("  Package loaded in %.2fs\n\n", t_pkg_load))

    # -- Phase 1: Julia/DAF initialization -------------------------------------
    cat("Phase 1: Julia/DAF initialization...\n")
    conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
    if (nzchar(conda_prefix)) {
        options(dafr.JULIA_HOME = file.path(conda_prefix, "bin"))
    }
    t0_julia <- proc.time()[["elapsed"]]
    dafr::setup_daf(pkg_check = FALSE, julia_environment = "custom")
    t_julia_init <- proc.time()[["elapsed"]] - t0_julia
    cat(sprintf("  Julia initialized in %.2fs\n", t_julia_init))

    # -- Phase 2: Open DAF dataset ---------------------------------------------
    cat("\nPhase 2: Opening DAF dataset...\n")
    t0_daf <- proc.time()[["elapsed"]]
    daf_obj <- dafr::open_daf(opts$daf)
    t_daf_open <- proc.time()[["elapsed"]] - t0_daf
    cat(sprintf("  DAF opened in %.2fs\n", t_daf_open))

    # -- Phase 3: MCView initialization ----------------------------------------
    cat("\nPhase 3: MCView initialization...\n")
    t0_mcview <- proc.time()[["elapsed"]]
    init_mcview_env()
    init_single_daf_mode(daf_obj, opts$dataset_name, cache_in_daf = FALSE)
    init_defs()
    t_mcview_init <- proc.time()[["elapsed"]] - t0_mcview
    cat(sprintf("  MCView initialized in %.2fs\n", t_mcview_init))

    # -- Phase 4: Discover dataset dimensions and pick test data ---------------
    cat("\nPhase 4: Dataset discovery...\n")
    test_genes <- pick_test_genes(daf_obj)
    test_mcs <- pick_test_metacells(daf_obj)
    cat(sprintf("  Genes:      %d total, test gene1=%s, gene2=%s\n",
                test_genes$n_genes, test_genes$gene1, test_genes$gene2))
    cat(sprintf("  Metacells:  %d total, test mc1=%s, mc2=%s\n",
                test_mcs$n_metacells, test_mcs$mc1, test_mcs$mc2))
    cat(sprintf("  Gene set:   [%s]\n", paste(test_genes$gene_set_small, collapse = ", ")))

    # -- Phase 5: Run benchmarks -----------------------------------------------
    cat("\n==============================================================================\n")
    cat("Phase 5: Running benchmarks\n")
    cat("==============================================================================\n\n")

    all_benchmarks <- define_benchmarks(
        opts$dataset_name, daf_obj, test_genes, test_mcs, opts$iterations
    )

    # Filter benchmarks if --benchmarks was specified
    if (!is.null(opts$benchmarks)) {
        unknown <- setdiff(opts$benchmarks, names(all_benchmarks))
        if (length(unknown) > 0) {
            warning("Unknown benchmark names: ", paste(unknown, collapse = ", "))
        }
        all_benchmarks <- all_benchmarks[intersect(opts$benchmarks, names(all_benchmarks))]
    }

    results <- list()
    for (bname in names(all_benchmarks)) {
        tryCatch({
            results[[bname]] <- all_benchmarks[[bname]]()
        }, error = function(e) {
            cat(sprintf("  [%s] SKIPPED due to error: %s\n", bname, conditionMessage(e)))
            results[[bname]] <<- list(
                name = bname,
                timings_s = numeric(0),
                median_s = NA_real_,
                mean_s = NA_real_,
                min_s = NA_real_,
                max_s = NA_real_,
                sd_s = NA_real_,
                iterations = 0L,
                errors = list(conditionMessage(e))
            )
        })
    }

    # -- Phase 6: Assemble and save results ------------------------------------
    cat("\n==============================================================================\n")
    cat("Phase 6: Results summary\n")
    cat("==============================================================================\n\n")

    # Print summary table
    cat(sprintf("%-35s %10s %10s %10s %10s\n",
                "Operation", "Median(s)", "Mean(s)", "Min(s)", "Max(s)"))
    cat(paste(rep("-", 80), collapse = ""), "\n")
    for (r in results) {
        if (!is.na(r$median_s)) {
            cat(sprintf("%-35s %10.4f %10.4f %10.4f %10.4f\n",
                        r$name, r$median_s, r$mean_s, r$min_s, r$max_s))
        } else {
            cat(sprintf("%-35s %10s %10s %10s %10s\n",
                        r$name, "ERROR", "ERROR", "ERROR", "ERROR"))
        }
    }

    # Build output JSON
    output <- list(
        metadata = list(
            timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
            daf_path = opts$daf,
            dataset_name = opts$dataset_name,
            iterations = opts$iterations,
            n_genes = test_genes$n_genes,
            n_metacells = test_mcs$n_metacells,
            test_gene1 = test_genes$gene1,
            test_gene2 = test_genes$gene2,
            test_gene_set = test_genes$gene_set_small,
            test_mc1 = test_mcs$mc1,
            test_mc2 = test_mcs$mc2,
            r_version = paste(R.version$major, R.version$minor, sep = "."),
            platform = Sys.info()[["sysname"]],
            hostname = Sys.info()[["nodename"]]
        ),
        init_timings = list(
            pkg_load_s = t_pkg_load,
            julia_init_s = t_julia_init,
            daf_open_s = t_daf_open,
            mcview_init_s = t_mcview_init,
            total_init_s = t_pkg_load + t_julia_init + t_daf_open + t_mcview_init
        ),
        benchmarks = results
    )

    # Determine output path
    output_path <- opts$output
    if (is.null(output_path)) {
        timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
        hostname <- gsub("[^a-zA-Z0-9_-]", "", Sys.info()[["nodename"]])
        output_dir <- file.path(pkg_root, "benchmarks", "results")
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        output_path <- file.path(output_dir, sprintf("bench_%s_%s.json", hostname, timestamp_str))
    }

    jsonlite::write_json(output, output_path, pretty = TRUE, auto_unbox = TRUE)
    cat(sprintf("\nResults saved to: %s\n", output_path))

    # Print init timing summary
    cat("\n--- Initialization timings ---\n")
    cat(sprintf("  Package load:    %6.2fs\n", t_pkg_load))
    cat(sprintf("  Julia init:      %6.2fs\n", t_julia_init))
    cat(sprintf("  DAF open:        %6.2fs\n", t_daf_open))
    cat(sprintf("  MCView init:     %6.2fs\n", t_mcview_init))
    cat(sprintf("  Total init:      %6.2fs\n", t_pkg_load + t_julia_init + t_daf_open + t_mcview_init))

    invisible(output)
}

# Run if executed as a script
if (!interactive()) {
    main()
}
