# Regression tests for two broken validation guards.

test_that("calc_marker_genes aborts when no gene clears the max-log-fraction threshold", {
    # All EGC values tiny => log2(egc + 1e-5) ~ -16.6, below the -13 default,
    # so the interesting-genes mask is all FALSE. The guard used length() of the
    # per-gene logical vector (never 0), so it never fired and silently returned
    # an empty tibble. It must abort instead.
    mc_egc <- matrix(1e-10,
        nrow = 5, ncol = 4,
        dimnames = list(paste0("g", 1:5), paste0("m", 1:4))
    )
    expect_error(
        calc_marker_genes(mc_egc, minimal_max_log_fraction = -13),
        "No genes"
    )
})

test_that("calc_marker_genes returns markers when some genes clear the threshold", {
    # Baseline-low matrix with explicit per-metacell spikes so each spiked gene
    # has a large relative (median-subtracted) log fold -> selected as a marker.
    mc_egc <- matrix(0.001,
        nrow = 6, ncol = 4,
        dimnames = list(paste0("g", 1:6), paste0("m", 1:4))
    )
    mc_egc["g1", "m1"] <- 0.5
    mc_egc["g2", "m2"] <- 0.5
    mc_egc["g3", "m3"] <- 0.5
    mc_egc["g4", "m4"] <- 0.5
    res <- calc_marker_genes(mc_egc,
        minimal_max_log_fraction = -13,
        minimal_relative_log_fraction = 2
    )
    expect_true(nrow(res) > 0)
})

test_that("parse_metacell_types gives a column-name error (not 'object file not found')", {
    bad <- tibble::tibble(not_metacell = 1:3, cell_type = c("A", "B", "C"))
    err <- tryCatch(parse_metacell_types(bad), error = function(e) conditionMessage(e))
    expect_match(err, "metacell")
    expect_false(grepl("object 'file' not found", err))
})
