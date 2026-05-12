# Phase 4 perf budgets - server-side cold-render time for the gene-gene
# scatter. The plot is the marquee user-facing cost on the Genes tab and
# was the primary motivator for Phase 4.
#
# Budgets are set conservatively (~30% above measured post-fix numbers
# on OBK) to absorb run-to-run page-cache variance. If this test fails
# locally, re-run it: a single bad run often passes on retry. CI should
# run it at least twice with the slower of the two used for the budget
# check.

# Helper: run the cold render once and return milliseconds for the
# data + ggplot construction step (the part Phase 4 optimised).
measure_genes_cold_ms <- function(daf_path = get_test_daf_path()) {
    daf_obj <- dafr::open_daf(daf_path)
    init_mcview_env()
    init_single_daf_mode(daf_obj, "data", NULL, FALSE)
    init_defs()
    config_shiny_cache()
    init_selected_genes()

    metacell_types <- get_mc_data("data", "metacell_types")
    cell_type_colors <- get_mc_data("data", "cell_type_colors")
    metadata <- get_mc_data("data", "metadata")
    gene_modules <- get_mc_data("data", "gene_modules")

    g1 <- mcv_get("default_gene1")
    g2 <- mcv_get("default_gene2")
    if (is.na(g1) || is.na(g2)) {
        gene_axis <- dafr::axis_entries(daf_obj, "gene")
        g1 <- gene_axis[1]; g2 <- gene_axis[2]
    }

    md <- metadata
    if (is.null(md)) md <- metacell_types %>% dplyr::select(metacell)
    md <- md %>% dplyr::mutate(Clipboard = "not selected")

    gc(verbose = FALSE)
    start <- proc.time()[["elapsed"]]
    plot_mc_scatter(
        "data", g1, g2,
        color_var = NULL, gene_modules = gene_modules, metadata = md,
        x_type = "Gene", y_type = "Gene", color_type = "Cell type",
        metacell_types = metacell_types, cell_type_colors = cell_type_colors,
        point_size = 2, stroke = 0.1, plot_text = FALSE
    )
    round((proc.time()[["elapsed"]] - start) * 1000)
}

test_that("genes-tab cold scatter build stays under budget", {
    skip_if_not(dir.exists(get_test_daf_path()), "OBK fixture not available")
    # Phase 4 takes plot_mc_scatter cold from ~1850 ms to ~680 ms on the
    # OBK fixture by skipping the eager mc_mat load (Task 2) and using a
    # targeted DAF query. The remaining ~680 ms is dominated by the two
    # `compute_egc_from_daf` calls with cold DAF / cold OS page cache.
    # Budget = 1000 ms - ~50% headroom over cold numbers, still catches a
    # regression back to the eager-load baseline.
    elapsed <- measure_genes_cold_ms()
    expect_lt(elapsed, 1000)
})
