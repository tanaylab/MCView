# Test cell-level DAF data access (daf_cells.R)
# Uses two DAF datasets:
#   - metacells_clean at /home/obk/data/mcview/metacells_clean (metacell axis, cell.metacell)
#   - cells_clean at /home/obk/data/mcview/cells_clean (cell axis, cell metadata, cell x gene UMIs)
# DAF setup and skip_if_no_daf() provided by helper-daf.R

get_cells_daf_path <- function() {
    "/home/obk/data/mcview/cells_clean"
}

# Helper: initialize metacells DAF then attach cells DAF
setup_cells_daf <- function() {
    skip_if_no_daf()

    cells_path <- get_cells_daf_path()
    if (!dir.exists(cells_path)) {
        skip("Cells DAF data not available at /home/obk/data/mcview/cells_clean")
    }

    daf <- dafr::open_daf(get_test_daf_path())
    init_mcview_env()
    init_single_daf_mode(daf, "test_data")

    set_cells_daf("test_data", cells_path)
    invisible(daf)
}

# ==============================================================================
# set_cells_daf / get_cells_daf
# ==============================================================================

test_that("set_cells_daf configures cells DAF path", {
    setup_cells_daf()

    mc_data <- mcv_get("mc_data")
    expect_equal(mc_data[["test_data"]][["cells_daf_path"]], get_cells_daf_path())
})

test_that("get_cells_daf returns a chained Daf object", {
    setup_cells_daf()

    chained <- get_cells_daf("test_data")
    expect_true(!is.null(chained))
    expect_true(dafr::is_daf(chained))
})

test_that("get_cells_daf auto-detects sibling cells DAF", {
    skip_if_no_daf()
    init_mcview_env()
    daf <- dafr::open_daf(get_test_daf_path())
    init_single_daf_mode(daf, "test_data")

    # With auto-detection, cells DAF should be found automatically
    # because metacells_clean has a sibling cells_clean directory
    result <- get_cells_daf("test_data")
    cells_path <- gsub("metacells", "cells", get_test_daf_path())
    if (dir.exists(cells_path)) {
        expect_true(!is.null(result))
        expect_true(dafr::is_daf(result))
    } else {
        expect_null(result)
    }
})

test_that("get_cells_daf returns NULL for dataset with no sibling cells dir", {
    skip_if_no_daf()
    init_mcview_env()
    daf <- dafr::open_daf(get_test_daf_path())
    init_single_daf_mode(daf, "test_data")

    # Manually clear the auto-detected cells DAF to test the NULL path
    mc_data <- mcv_get("mc_data")
    mc_data[["test_data"]][["cells_daf_path"]] <- NULL
    mc_data[["test_data"]][["cells_daf"]] <- NULL
    mc_data[["test_data"]][["chained_cells_daf"]] <- NULL
    mcv_set("mc_data", mc_data)

    result <- get_cells_daf("test_data")
    expect_null(result)
})

test_that("get_cells_daf caches the chained DAF on second call", {
    setup_cells_daf()

    daf1 <- get_cells_daf("test_data")
    daf2 <- get_cells_daf("test_data")
    # Both should be the exact same cached object
    expect_identical(daf1, daf2)
})

# ==============================================================================
# has_cell_gene_umis
# ==============================================================================

test_that("has_cell_gene_umis returns TRUE when cells_clean is configured", {
    setup_cells_daf()

    expect_true(has_cell_gene_umis("test_data"))
})

test_that("has_cell_gene_umis returns FALSE when cells DAF cleared", {
    skip_if_no_daf()
    init_mcview_env()
    daf <- dafr::open_daf(get_test_daf_path())
    init_single_daf_mode(daf, "test_data")

    # Clear the auto-detected cells DAF to test the FALSE path
    mc_data <- mcv_get("mc_data")
    mc_data[["test_data"]][["cells_daf_path"]] <- NULL
    mc_data[["test_data"]][["cells_daf"]] <- NULL
    mc_data[["test_data"]][["chained_cells_daf"]] <- NULL
    mcv_set("mc_data", mc_data)

    expect_false(has_cell_gene_umis("test_data"))
})

# ==============================================================================
# get_cell_metacell_map
# ==============================================================================

test_that("get_cell_metacell_map returns named character vector", {
    setup_cells_daf()

    mc_map <- get_cell_metacell_map("test_data")

    expect_type(mc_map, "character")
    expect_gt(length(mc_map), 0)
    expect_true(!is.null(names(mc_map)))
})

test_that("get_cell_metacell_map values are valid metacell IDs", {
    setup_cells_daf()

    mc_map <- get_cell_metacell_map("test_data")
    daf_obj <- get_dataset_daf("test_data")
    valid_metacells <- dafr::axis_entries(daf_obj, "metacell")

    # Filter out empty strings (outlier/unassigned cells)
    non_empty <- mc_map[nchar(mc_map) > 0]
    expect_gt(length(non_empty), 0)
    # All non-empty mapped metacell IDs should be in the metacell axis
    expect_true(all(unique(non_empty) %in% valid_metacells))
})

# ==============================================================================
# get_cell_field_map
# ==============================================================================

test_that("get_cell_field_map returns named character vector for batch_set_id", {
    setup_cells_daf()

    field_map <- get_cell_field_map("test_data", "batch_set_id")

    expect_type(field_map, "character")
    expect_gt(length(field_map), 0)
    expect_true(!is.null(names(field_map)))
    # batch_set_id values are strings
    expect_true(all(nchar(field_map) > 0))
})

test_that("get_cell_field_map returns NULL for nonexistent field", {
    setup_cells_daf()

    result <- get_cell_field_map("test_data", "nonexistent_field_xyz")
    expect_null(result)
})

# ==============================================================================
# get_cell_grouping_fields
# ==============================================================================

test_that("get_cell_grouping_fields returns character vector", {
    setup_cells_daf()

    fields <- get_cell_grouping_fields("test_data")

    expect_type(fields, "character")
    expect_gt(length(fields), 0)
})

test_that("get_cell_grouping_fields includes batch_set_id", {
    setup_cells_daf()

    fields <- get_cell_grouping_fields("test_data")

    expect_true("batch_set_id" %in% fields)
})

test_that("get_cell_grouping_fields returns empty when cells DAF cleared", {
    skip_if_no_daf()
    init_mcview_env()
    daf <- dafr::open_daf(get_test_daf_path())
    init_single_daf_mode(daf, "test_data")

    # Clear the auto-detected cells DAF to test the empty path
    mc_data <- mcv_get("mc_data")
    mc_data[["test_data"]][["cells_daf_path"]] <- NULL
    mc_data[["test_data"]][["cells_daf"]] <- NULL
    mc_data[["test_data"]][["chained_cells_daf"]] <- NULL
    mcv_set("mc_data", mc_data)

    fields <- get_cell_grouping_fields("test_data")
    expect_equal(length(fields), 0)
})

# ==============================================================================
# get_group_cell_type_composition
# ==============================================================================

test_that("get_group_cell_type_composition returns tibble with expected columns", {
    setup_cells_daf()

    comp <- get_group_cell_type_composition("test_data", "batch_set_id")

    expect_s3_class(comp, "tbl_df")
    expect_true(all(c("group_id", "cell_type", "n_cells", "fraction") %in% names(comp)))
    expect_gt(nrow(comp), 0)
})

test_that("get_group_cell_type_composition fractions sum to ~1 per group", {
    setup_cells_daf()

    comp <- get_group_cell_type_composition("test_data", "batch_set_id")

    fraction_sums <- comp %>%
        dplyr::group_by(group_id) %>%
        dplyr::summarise(total = sum(fraction), .groups = "drop")

    expect_true(all(abs(fraction_sums$total - 1) < 1e-6))
})

test_that("get_group_cell_type_composition n_cells are positive integers", {
    setup_cells_daf()

    comp <- get_group_cell_type_composition("test_data", "batch_set_id")

    expect_true(all(comp$n_cells > 0))
    expect_true(all(comp$n_cells == as.integer(comp$n_cells)))
})

# ==============================================================================
# get_group_gene_expression
# ==============================================================================

test_that("get_group_gene_expression returns named numeric vector", {
    setup_cells_daf()

    # Get a known gene from the gene axis
    cells_daf <- get_cells_daf("test_data")
    test_gene <- dafr::axis_entries(cells_daf, "gene")[1]

    expr <- get_group_gene_expression("test_data", test_gene, "batch_set_id")

    expect_type(expr, "double")
    expect_gt(length(expr), 0)
    expect_true(!is.null(names(expr)))
})

test_that("get_group_gene_expression values are non-negative", {
    setup_cells_daf()

    cells_daf <- get_cells_daf("test_data")
    test_gene <- dafr::axis_entries(cells_daf, "gene")[1]

    expr <- get_group_gene_expression("test_data", test_gene, "batch_set_id")

    expect_true(all(expr >= 0))
})

# ==============================================================================
# get_group_qc_stats
# ==============================================================================

test_that("get_group_qc_stats returns tibble with expected columns", {
    setup_cells_daf()

    qc <- get_group_qc_stats("test_data", "batch_set_id")

    expect_s3_class(qc, "tbl_df")
    expect_true(all(c("group_id", "n_cells", "total_umis", "median_umis_per_cell") %in% names(qc)))
    expect_gt(nrow(qc), 0)
})

test_that("get_group_qc_stats n_cells are positive", {
    setup_cells_daf()

    qc <- get_group_qc_stats("test_data", "batch_set_id")

    expect_true(all(qc$n_cells > 0))
})

test_that("get_group_qc_stats total_umis are positive when available", {
    setup_cells_daf()

    qc <- get_group_qc_stats("test_data", "batch_set_id")

    # cells_clean has total_UMIs, so these should be populated
    expect_true(all(!is.na(qc$total_umis)))
    expect_true(all(qc$total_umis > 0))
})

test_that("get_group_qc_stats median_umis_per_cell are positive when available", {
    setup_cells_daf()

    qc <- get_group_qc_stats("test_data", "batch_set_id")

    expect_true(all(!is.na(qc$median_umis_per_cell)))
    expect_true(all(qc$median_umis_per_cell > 0))
})

# ==============================================================================
# get_default_sample_field
# ==============================================================================

test_that("get_default_sample_field returns a string", {
    setup_cells_daf()

    field <- get_default_sample_field("test_data")

    expect_type(field, "character")
    expect_equal(length(field), 1)
    expect_true(nchar(field) > 0)
})

test_that("get_default_sample_field returns NULL for unconfigured dataset", {
    skip_if_no_daf()
    init_mcview_env()

    # Create a minimal dataset with no cell axis at all
    daf <- dafr::open_daf(get_test_daf_path())
    init_single_daf_mode(daf, "test_data")

    # Clear the auto-detected cells DAF to test the NULL path
    mc_data <- mcv_get("mc_data")
    mc_data[["test_data"]][["cells_daf_path"]] <- NULL
    mc_data[["test_data"]][["cells_daf"]] <- NULL
    mc_data[["test_data"]][["chained_cells_daf"]] <- NULL
    mcv_set("mc_data", mc_data)

    # metacells_clean does have a cell axis with no samp_id,
    # but no cells DAF is set; get_default_sample_field will check
    # the metacells DAF first for samp_id, then fall back to cells DAF.
    # Without cells DAF, and without samp_id in metacells, this should
    # return NULL (or possibly a field from the metacells DAF cell axis).
    field <- get_default_sample_field("test_data")

    # The metacells_clean DAF has a cell axis but only "metacell" as a vector,
    # no samp_id. Without cells DAF set, cells_daf will be NULL, so no
    # preferred fields found -> grouping_fields will be empty -> NULL.
    expect_null(field)
})
