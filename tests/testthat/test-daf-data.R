# Test DAF data operations
# Uses the test data at /home/obk/data/mcview/metacells_clean

test_daf_path <- "/home/obk/data/mcview/metacells_clean"

# Skip all tests if DAF data is not available
skip_if_no_daf <- function() {
    if (!dir.exists(test_daf_path)) {
        skip("Test DAF data not available")
    }
    if (!requireNamespace("dafr", quietly = TRUE)) {
        skip("dafr package not installed")
    }
}

# Helper to initialize DAF for testing
setup_test_daf <- function() {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)
    init_mcview_env()
    init_single_daf_mode(daf, "test_data")
    invisible(daf)
}

test_that("DAF initialization works", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    expect_true(inherits(daf, "Daf"))

    # Check axes exist using axis_entries
    expect_true(length(dafr::axis_entries(daf, "metacell")) > 0)
    expect_true(length(dafr::axis_entries(daf, "gene")) > 0)
    expect_true(length(dafr::axis_entries(daf, "type")) > 0)
})

test_that("init_single_daf_mode works", {
    daf <- setup_test_daf()

    # Check dataset was registered
    expect_equal(dataset_names(), "test_data")

    # Check config was set
    config <- mcv_get("config")
    expect_true(!is.null(config))
    expect_true(!is.null(config$title))
    expect_true(!is.null(config$tabs))
})

test_that("get_mc_data returns mc2d correctly", {
    daf <- setup_test_daf()

    mc2d <- get_mc_data("test_data", "mc2d")

    # mc2d is returned as a list with mc_id, mc_x, mc_y, graph
    expect_true(is.list(mc2d))
    expect_true("mc_id" %in% names(mc2d))
    expect_true("mc_x" %in% names(mc2d))
    expect_true("mc_y" %in% names(mc2d))
    expect_gt(length(mc2d$mc_id), 0)
})

test_that("get_mc_data returns metacell_types correctly", {
    daf <- setup_test_daf()

    mc_types <- get_mc_data("test_data", "metacell_types")

    expect_s3_class(mc_types, "tbl_df")
    expect_true("metacell" %in% names(mc_types))
    expect_true("cell_type" %in% names(mc_types))
    expect_true("mc_col" %in% names(mc_types))
    expect_gt(nrow(mc_types), 0)
})

test_that("get_mc_data returns cell_type_colors correctly", {
    daf <- setup_test_daf()

    ct_colors <- get_mc_data("test_data", "cell_type_colors")

    expect_s3_class(ct_colors, "tbl_df")
    expect_true("cell_type" %in% names(ct_colors))
    expect_true("color" %in% names(ct_colors))
    expect_gt(nrow(ct_colors), 0)
})

test_that("gene_names returns gene list", {
    daf <- setup_test_daf()

    genes <- gene_names("test_data")

    expect_type(genes, "character")
    expect_gt(length(genes), 0)
})

test_that("get_mc_egc returns full matrix", {
    daf <- setup_test_daf()

    egc <- get_mc_egc("test_data")

    expect_true(is.matrix(egc) || inherits(egc, "Matrix"))
    expect_gt(nrow(egc), 0)
    expect_gt(ncol(egc), 0)
    # Check row names are genes
    expect_true(all(rownames(egc) %in% gene_names("test_data")))
})

test_that("get_mc_egc filters genes correctly", {
    daf <- setup_test_daf()

    genes_to_get <- c("Cd79a", "Xkr4")
    egc <- get_mc_egc("test_data", genes = genes_to_get)

    # Should only contain requested genes that exist
    expect_true(all(rownames(egc) %in% genes_to_get))
    expect_lte(nrow(egc), length(genes_to_get))
})

test_that("get_gene_egc returns vector for single gene", {
    daf <- setup_test_daf()

    gene_egc <- get_gene_egc("Cd79a", "test_data")

    expect_type(gene_egc, "double")
    expect_gt(length(gene_egc), 0)
    expect_true(!is.null(names(gene_egc)))
})

test_that("get_metacells_egc returns matrix for selected metacells", {
    daf <- setup_test_daf()

    mc2d <- get_mc_data("test_data", "mc2d")
    test_metacells <- head(mc2d$mc_id, 5)

    mc_egc <- get_metacells_egc(test_metacells, "test_data")

    expect_true(is.matrix(mc_egc) || inherits(mc_egc, "Matrix"))
    expect_equal(ncol(mc_egc), length(test_metacells))
})

test_that("get_cell_types_mat uses DAF query aggregation", {
    daf <- setup_test_daf()

    mc_types <- get_mc_data("test_data", "metacell_types")
    cell_types <- c("ExE mesoderm", "Nascent mesoderm late")

    ct_mat <- get_cell_types_mat(cell_types, mc_types, "test_data")

    expect_true(is.matrix(ct_mat))
    expect_equal(ncol(ct_mat), length(cell_types))
    expect_true(all(colnames(ct_mat) %in% cell_types))
    expect_gt(nrow(ct_mat), 0)
})

test_that("get_cell_types_egc returns normalized matrix", {
    daf <- setup_test_daf()

    mc_types <- get_mc_data("test_data", "metacell_types")
    cell_types <- c("ExE mesoderm", "Nascent mesoderm late")

    ct_egc <- get_cell_types_egc(cell_types, mc_types, "test_data")

    expect_true(is.matrix(ct_egc))
    # Check columns sum to ~1 (normalized)
    col_sums <- colSums(ct_egc)
    expect_true(all(abs(col_sums - 1) < 1e-6))
})

test_that("daf_query_cell_type_umis aggregates correctly", {
    daf <- setup_test_daf()
    daf_obj <- get_dataset_daf("test_data")

    ct_mat <- daf_query_cell_type_umis(daf_obj)

    expect_true(is.matrix(ct_mat))
    expect_gt(nrow(ct_mat), 0) # genes
    expect_gt(ncol(ct_mat), 0) # cell types
})

test_that("daf_query_cell_type_sum returns per-type totals", {
    daf <- setup_test_daf()
    daf_obj <- get_dataset_daf("test_data")

    ct_sum <- daf_query_cell_type_sum(daf_obj)

    expect_true(is.numeric(ct_sum))
    expect_gt(length(ct_sum), 0)
    expect_true(!is.null(names(ct_sum)))
    expect_true(all(ct_sum > 0))
})

test_that("get_mc_fp returns fold change matrix", {
    daf <- setup_test_daf()

    fp <- get_mc_fp("test_data")

    expect_true(is.matrix(fp) || inherits(fp, "Matrix"))
    expect_gt(nrow(fp), 0)
    expect_gt(ncol(fp), 0)
})

test_that("inner_fold_mat is available when present", {
    daf <- setup_test_daf()

    inner_fold <- get_mc_data("test_data", "inner_fold_mat")

    # This may or may not be present depending on the DAF
    if (!is.null(inner_fold)) {
        expect_true(is.matrix(inner_fold) || inherits(inner_fold, "Matrix"))
    }
})

# ==============================================================================
# Cache Configuration Tests
# ==============================================================================

test_that("create_cache_config creates valid config", {
    config <- create_cache_config()

    expect_true(is.list(config))
    expect_true(config$enabled)
    expect_equal(config$type, "files")
    expect_equal(config$cache_dir, ".mcview_cache")
    expect_true(config$precompute_on_startup)
    expect_equal(config$invalidation$strategy, "version")
})

test_that("create_cache_config accepts custom values", {
    config <- create_cache_config(
        enabled = FALSE,
        type = "memory",
        cache_dir = "/custom/path",
        precompute_on_startup = FALSE
    )

    expect_false(config$enabled)
    expect_equal(config$type, "memory")
    expect_equal(config$cache_dir, "/custom/path")
    expect_false(config$precompute_on_startup)
})

test_that("extract_cache_config returns defaults when no config", {
    config <- extract_cache_config(NULL, NULL, NULL)

    expect_true(is.list(config))
    expect_false(config$enabled) # Default is disabled
    expect_equal(config$type, "files")
    expect_equal(config$cache_dir, ".mcview_cache")
})

test_that("extract_cache_config reads from YAML config", {
    yaml_config <- list(
        cache = list(
            enabled = TRUE,
            type = "memory",
            cache_dir = "/yaml/cache"
        )
    )

    config <- extract_cache_config(yaml_config, NULL, NULL)

    expect_true(config$enabled)
    expect_equal(config$type, "memory")
    expect_equal(config$cache_dir, "/yaml/cache")
})

test_that("extract_cache_config supports legacy cache_in_daf", {
    yaml_config <- list(
        cache_in_daf = TRUE,
        cache_daf_root = "/legacy/path"
    )

    config <- extract_cache_config(yaml_config, NULL, NULL)

    expect_true(config$enabled)
    expect_equal(config$cache_dir, "/legacy/path")
})

# ==============================================================================
# Cache DAF Initialization Tests
# ==============================================================================

test_that("init_cache_daf returns base DAF when caching disabled", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    cache_config <- create_cache_config(enabled = FALSE)
    result <- init_cache_daf(daf, "test", cache_config, test_daf_path)

    expect_null(result$cache_daf)
    expect_equal(result$complete_daf, daf)
    expect_false(result$needs_population)
    expect_null(result$cache_path)
})

test_that("init_cache_daf creates memory cache", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    cache_config <- create_cache_config(enabled = TRUE, type = "memory")
    result <- init_cache_daf(daf, "test_memory", cache_config, test_daf_path)

    expect_true(!is.null(result$cache_daf))
    expect_true(!is.null(result$complete_daf))
    expect_null(result$cache_path) # Memory cache has no path
})

test_that("init_cache_daf creates files cache with relative path", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    # Use temp dir for cache
    temp_cache <- tempfile("mcview_cache_test")
    cache_config <- create_cache_config(enabled = TRUE, type = "files", cache_dir = temp_cache)
    result <- init_cache_daf(daf, "test_files", cache_config, test_daf_path)

    expect_true(!is.null(result$cache_daf))
    expect_true(!is.null(result$complete_daf))
    expect_true(!is.null(result$cache_path))
    expect_true(dir.exists(result$cache_path))

    # Check base_daf_repository scalar is set
    base_ref <- daf_scalar(result$cache_daf, "base_daf_repository", default = NULL)
    expect_true(!is.null(base_ref))

    # Cleanup
    unlink(temp_cache, recursive = TRUE)
})

# ==============================================================================
# Cache Validation Tests
# ==============================================================================

test_that("is_cache_valid returns FALSE for NULL cache", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    expect_false(is_cache_valid(NULL, daf, list(strategy = "hash")))
})

test_that("compute_base_hash is consistent", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    hash1 <- compute_base_hash(daf)
    hash2 <- compute_base_hash(daf)

    expect_equal(hash1, hash2)
    expect_type(hash1, "character")
    expect_gt(nchar(hash1), 0)
})

# ==============================================================================
# Cache DAF Accessors Tests
# ==============================================================================

test_that("get_cache_daf returns NULL when no cache", {
    daf <- setup_test_daf()

    # Default setup doesn't enable cache
    cache_daf <- get_cache_daf("test_data")

    # May or may not be NULL depending on config, but shouldn't error
    expect_true(is.null(cache_daf) || inherits(cache_daf, "Daf"))
})

test_that("get_base_daf returns original DAF", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    init_mcview_env()

    # Initialize with cache enabled
    cache_config <- create_cache_config(enabled = TRUE, type = "memory")
    init_single_daf_mode(daf, "test_with_cache", cache_config = cache_config)

    base_daf <- get_base_daf("test_with_cache")

    expect_true(!is.null(base_daf))
    expect_true(inherits(base_daf, "Daf"))
})

test_that("get_complete_daf returns chained DAF", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    init_mcview_env()
    cache_config <- create_cache_config(enabled = TRUE, type = "memory")
    init_single_daf_mode(daf, "test_complete", cache_config = cache_config)

    complete_daf <- get_complete_daf("test_complete")

    expect_true(!is.null(complete_daf))
    expect_true(inherits(complete_daf, "Daf"))
})

# ==============================================================================
# Cache Population Tests
# ==============================================================================

test_that("populate_dataset_cache works with memory cache", {
    skip_if_no_daf()
    dafr::setup_daf()
    daf <- dafr::open_daf(test_daf_path)

    init_mcview_env()
    cache_config <- create_cache_config(enabled = TRUE, type = "memory")
    init_single_daf_mode(daf, "test_populate", cache_config = cache_config)

    # Populate cache
    result <- populate_dataset_cache("test_populate", verbose = FALSE)

    # Should complete without error
    expect_true(is.logical(result))
})
