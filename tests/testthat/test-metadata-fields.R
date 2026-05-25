# Regression test: convert_daf_metadata() must not surface MCView-derived
# cache vectors (mcview_cache_*) or DAF-internal (double-underscore) properties
# as user-selectable metadata fields. They cluttered the metadata selectors and
# could break metadata plots (e.g. the character mcview_cache_top1_gene).

test_that("convert_daf_metadata drops derived cache and internal properties", {
    d <- dafr::memory_daf("md_fixture")
    mcs <- paste0("m", 1:4)
    dafr::add_axis(d, "metacell", mcs)

    # A genuine metadata field, a core field, plus the artifacts to exclude.
    dafr::set_vector(d, "metacell", "age", c(1.0, 2.0, 3.0, 4.0))
    dafr::set_vector(d, "metacell", "type", rep("T", 4))
    dafr::set_vector(d, "metacell", "mcview_cache_top1_gene", rep("GeneA", 4))
    dafr::set_vector(d, "metacell", "mcview_cache_top1_lfp", c(0.1, 0.2, 0.3, 0.4))
    dafr::set_vector(d, "metacell", "__zeros_downsample_UMIs", c(0.0, 0.0, 0.0, 0.0))

    md <- convert_daf_metadata(d)

    expect_true("age" %in% colnames(md))
    expect_false(any(grepl("^mcview_cache_", colnames(md))))
    expect_false(any(grepl("^__", colnames(md))))
    # "type" is a core field handled elsewhere - also excluded.
    expect_false("type" %in% colnames(md))
})
