# Regression test: get_cell_types_mat() must aggregate UMIs by the PASSED
# metacell_types mapping (so DE / obs-exp reflect session re-annotations),
# not the DAF's stored `type` vector. The DAF-query rewrite (75dfe01) had
# dropped both the metacell_types argument and the projected/corrected flags.

setup_ct_mat_dataset <- function() {
    genes <- c("g1", "g2", "g3")
    mcs <- c("m1", "m2", "m3", "m4")
    umis <- matrix(
        c(
            10, 20, 30, 40, # g1
            1, 2, 3, 4, # g2
            5, 6, 7, 8 # g3
        ),
        nrow = 3, byrow = TRUE, dimnames = list(genes, mcs)
    )

    d <- dafr::memory_daf("ct_mat_fixture")
    dafr::add_axis(d, "gene", genes)
    dafr::add_axis(d, "metacell", mcs)
    dafr::set_matrix(d, "gene", "metacell", "UMIs", umis)
    dafr::set_vector(d, "metacell", "total_UMIs", colSums(umis))
    # DAF's stored type: m1,m2 -> A ; m3,m4 -> B
    dafr::set_vector(d, "metacell", "type", c("A", "A", "B", "B"))

    init_mcview_env()
    prev_mc_data <- mcv_get("mc_data")
    prev_eps <- mcv_get("egc_epsilon")
    mcv_set("egc_epsilon", 1e-5)
    mcv_set("mc_data", list(test_ds = list(daf_obj = d)))
    list(
        umis = umis,
        restore = function() {
            mcv_set("mc_data", prev_mc_data)
            mcv_set("egc_epsilon", prev_eps)
        }
    )
}

test_that("get_cell_types_mat aggregates by the passed metacell_types, not the DAF type", {
    ctx <- setup_ct_mat_dataset()
    on.exit(ctx$restore(), add = TRUE)

    # Re-annotation: move m2 from A to B. A = {m1}, B = {m2, m3, m4}.
    reannotated <- tibble::tibble(
        metacell = c("m1", "m2", "m3", "m4"),
        cell_type = c("A", "B", "B", "B")
    )

    mat <- get_cell_types_mat(c("A", "B"), reannotated, "test_ds")

    expect_equal(sort(colnames(mat)), c("A", "B"))
    expect_equal(rownames(mat), c("g1", "g2", "g3"))
    # A = m1 only; B = m2 + m3 + m4 (per the PASSED mapping).
    expect_equal(unname(mat[, "A"]), unname(ctx$umis[, "m1"]))
    expect_equal(unname(mat[, "B"]), unname(ctx$umis[, "m2"] + ctx$umis[, "m3"] + ctx$umis[, "m4"]))
})

test_that("get_cell_types_egc normalizes each cell type by its own total", {
    ctx <- setup_ct_mat_dataset()
    on.exit(ctx$restore(), add = TRUE)

    types <- tibble::tibble(metacell = c("m1", "m2", "m3", "m4"), cell_type = c("A", "A", "B", "B"))
    egc <- get_cell_types_egc(c("A", "B"), types, "test_ds")
    # Each column should sum to 1 (fractions per cell type).
    expect_equal(unname(colSums(egc)), c(1, 1), tolerance = 1e-12)
})
