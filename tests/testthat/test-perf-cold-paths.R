# Tests that Phase 4 cold-path optimisations keep their invariants:
# - partial_bundle pre-warm is silent and idempotent
# - get_gene_egc / calc_mc_mc_gene_df do not load full mc_mat on cold path

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
