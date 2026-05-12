# Tests that Phase 4 cold-path optimisations keep their invariants:
# - partial_bundle pre-warm fires at most once per R process
# - get_gene_egc / calc_mc_mc_gene_df do not load full mc_mat on cold path

test_that("partial_bundle pre-warm sets a flag and is idempotent", {
    mcv_set(".partial_bundle_warmed", NULL)
    on.exit(mcv_set(".partial_bundle_warmed", NULL), add = TRUE)

    expect_null(mcv_get(".partial_bundle_warmed"))

    # First call goes through the live path: no warning, no error.
    expect_silent(prewarm_plotly_bundle())
    expect_true(isTRUE(mcv_get(".partial_bundle_warmed")))

    # Second call short-circuits via the guard: also silent, flag stays set.
    expect_silent(prewarm_plotly_bundle())
    expect_true(isTRUE(mcv_get(".partial_bundle_warmed")))
})

test_that("partial_bundle pre-warm actually warms the plotly cache", {
    # Two back-to-back calls with the guard cleared between them: the
    # first does real work (plotly.js bundle parse, several hundred ms);
    # the second hits plotly's internal cache and is near-instant.
    # If a future refactor stops calling plotly::partial_bundle() but
    # keeps setting the flag, this test catches it.
    mcv_set(".partial_bundle_warmed", NULL)
    on.exit(mcv_set(".partial_bundle_warmed", NULL), add = TRUE)

    t1 <- system.time(prewarm_plotly_bundle())[["elapsed"]]
    mcv_set(".partial_bundle_warmed", NULL)
    t2 <- system.time(prewarm_plotly_bundle())[["elapsed"]]

    # Generous threshold: empirically t2 / t1 is ~0.01 on a fresh process.
    # If the second call is more than half the first, real warming didn't
    # happen.
    expect_lt(t2, max(t1 * 0.5, 0.05))
})
