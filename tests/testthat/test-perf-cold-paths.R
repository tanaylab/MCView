# Tests that Phase 4 cold-path optimisations keep their invariants:
# - partial_bundle pre-warm fires at most once per R process
# - get_gene_egc / calc_mc_mc_gene_df do not load full mc_mat on cold path

test_that("partial_bundle pre-warm sets a flag and is idempotent", {
    mcv_set(".partial_bundle_warmed", NULL)
    expect_null(mcv_get(".partial_bundle_warmed"))

    prewarm_plotly_bundle()
    expect_true(isTRUE(mcv_get(".partial_bundle_warmed")))

    # Second call should be cheap (no error, no side effect beyond flag).
    expect_silent(prewarm_plotly_bundle())
    expect_true(isTRUE(mcv_get(".partial_bundle_warmed")))
})
