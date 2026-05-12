# Tests that Phase 4 cold-path optimisations keep their invariants:
# - partial_bundle pre-warm is silent and idempotent
# - get_gene_egc / calc_mc_mc_gene_df do not load full mc_mat on cold path

test_that("get_gene_egc does not eagerly load full mc_mat on cold path", {
    skip_if_not(dir.exists(get_test_daf_path()), "OBK fixture not available")

    # Fresh session: open DAF, init env. Do NOT touch mc_mat / mc_egc_full.
    daf_obj <- dafr::open_daf(get_test_daf_path())
    init_mcview_env()
    init_single_daf_mode(daf_obj, "data", NULL, FALSE)
    init_defs()
    config_shiny_cache()

    expect_null(mcv_cache_get("data", "mc_mat"))
    expect_null(mcv_cache_get("data", "mc_egc_full"))

    g1 <- dafr::axis_entries(daf_obj, "gene")[1]
    vec <- get_gene_egc(g1, "data")
    expect_true(is.numeric(vec))
    expect_true(length(vec) > 0)

    # The cold call must NOT have populated mc_mat (the regression we are
    # guarding against).
    expect_null(mcv_cache_get("data", "mc_mat"))
})

test_that("partial_bundle pre-warm is silent and idempotent", {
    mcv_set(".partial_bundle_warmed", NULL)
    on.exit(mcv_set(".partial_bundle_warmed", NULL), add = TRUE)

    expect_null(mcv_get(".partial_bundle_warmed"))

    # First call goes through the live path: no warning, no error.
    # This is the regression guard against a future change that drops
    # the typed plot_ly() and re-introduces the "No trace type" warning.
    expect_silent(prewarm_plotly_bundle())
    expect_true(isTRUE(mcv_get(".partial_bundle_warmed")))

    # Second call short-circuits via the guard: also silent, flag stays set.
    expect_silent(prewarm_plotly_bundle())
    expect_true(isTRUE(mcv_get(".partial_bundle_warmed")))
})

test_that("get_gene_egc cold path matches eager-load value (no silent drift)", {
    skip_if_not(dir.exists(get_test_daf_path()), "OBK fixture not available")

    daf_obj <- dafr::open_daf(get_test_daf_path())
    init_mcview_env()
    init_single_daf_mode(daf_obj, "data", NULL, FALSE)
    init_defs()
    config_shiny_cache()

    gene_axis <- dafr::axis_entries(daf_obj, "gene")
    # Sample a few genes: first, middle, last. Catches positional / ordering
    # bugs in compute_egc_from_daf's row extraction.
    sample_genes <- gene_axis[c(1, floor(length(gene_axis) / 2), length(gene_axis))]

    # Reference: eager-load path. Triggers full mc_mat into cache.
    mat <- get_mc_data("data", "mc_mat")
    mc_sum <- get_mc_sum("data")
    ref <- lapply(sample_genes, function(g) {
        umis <- mat[g, ]
        umis / mc_sum[names(umis)]
    })
    names(ref) <- sample_genes

    # Clear caches so the next call goes through the cold targeted-query path.
    mcv_cache_set("data", "mc_mat", NULL)
    mcv_cache_set("data", "mc_egc_full", NULL)
    expect_null(mcv_cache_get("data", "mc_mat"))

    for (g in sample_genes) {
        v <- get_gene_egc(g, "data")
        # Align names: ref keeps metacell axis order; v should too.
        expect_equal(v[names(ref[[g]])], ref[[g]], tolerance = 1e-10,
                     info = paste("gene =", g))
    }

    # And confirm cold path did NOT re-populate mc_mat.
    expect_null(mcv_cache_get("data", "mc_mat"))
})
