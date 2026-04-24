test_that("get_group_qc_stats new path matches legacy 3-query output", {
    skip_if_not(dir.exists("/home/obk/data/mcview/metacells_clean"),
                "OBK dataset unavailable")
    skip_if_not(requireNamespace("dafr", quietly = TRUE))

    daf_path <- "/home/obk/data/mcview/metacells_clean"
    dataset <- "_test_qc_stats"
    MCView:::init_mcview_env()
    daf_obj <- dafr::files_daf(daf_path, mode = "r")
    MCView:::init_single_daf_mode(daf_obj, dataset, NULL, FALSE)

    cells_daf <- MCView:::get_cells_daf(dataset)
    skip_if(is.null(cells_daf), "cells DAF not present")

    fields <- MCView:::get_cell_grouping_fields(dataset)
    skip_if(length(fields) == 0L, "no grouping fields available")
    gf <- fields[[1L]]

    result <- MCView::get_group_qc_stats(dataset, gf)

    count_q <- sprintf("@ cell : %s / %s >> Count", gf, gf)
    sum_q   <- sprintf("@ cell : total_UMIs / %s >> Sum", gf)
    med_q   <- sprintf("@ cell : total_UMIs / %s >> Median", gf)
    legacy_counts <- cells_daf[count_q]
    legacy_sums   <- tryCatch(cells_daf[sum_q], error = function(e) NULL)
    legacy_meds   <- tryCatch(cells_daf[med_q], error = function(e) NULL)

    expect_setequal(result$group_id, names(legacy_counts))
    expect_equal(
        result$n_cells[match(names(legacy_counts), result$group_id)],
        as.integer(legacy_counts)
    )
    if (!is.null(legacy_sums)) {
        expect_equal(
            result$total_umis[match(names(legacy_sums), result$group_id)],
            as.numeric(legacy_sums),
            tolerance = 1e-9
        )
    }
    if (!is.null(legacy_meds)) {
        expect_equal(
            result$median_umis_per_cell[match(names(legacy_meds), result$group_id)],
            as.numeric(legacy_meds),
            tolerance = 1e-9
        )
    }
})
