# Regression test: daf_dir_fingerprint() (used by compute_base_hash for cache
# invalidation) must change when file CONTENT changes, including same-size
# rewrites, so a stale derived cache is not served after the base DAF is
# re-annotated. The previous fingerprint was size-only and missed same-size
# rewrites.

test_that("daf_dir_fingerprint is deterministic and excludes the derived subtree", {
    dir <- tempfile()
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    writeLines("hello", file.path(dir, "a.txt"))
    dir.create(file.path(dir, ".mcview_derived"))
    writeLines("cache", file.path(dir, ".mcview_derived", "c.txt"))

    fp1 <- daf_dir_fingerprint(dir)
    fp2 <- daf_dir_fingerprint(dir)
    expect_identical(fp1, fp2)
    expect_true(nzchar(fp1))
    # Writing into the excluded derived subtree must not change the fingerprint.
    writeLines("more cache", file.path(dir, ".mcview_derived", "d.txt"))
    expect_identical(daf_dir_fingerprint(dir), fp1)
})

test_that("daf_dir_fingerprint changes on a same-size content rewrite (via mtime)", {
    dir <- tempfile()
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    f <- file.path(dir, "vec.data")
    writeBin(as.raw(c(1, 2, 3, 4)), f)

    fp_before <- daf_dir_fingerprint(dir)

    # Same byte length, different content, later mtime (simulates re-annotation
    # rewriting a fixed-width vector at identical size).
    writeBin(as.raw(c(9, 9, 9, 9)), f)
    Sys.setFileTime(f, Sys.time() + 5)

    expect_equal(file.size(f), 4) # size unchanged - a size-only hash would miss this
    expect_false(identical(daf_dir_fingerprint(dir), fp_before))
})

test_that("daf_dir_fingerprint returns empty for missing/empty dirs", {
    expect_identical(daf_dir_fingerprint(NULL), "")
    expect_identical(daf_dir_fingerprint(file.path(tempdir(), "does-not-exist-xyz")), "")
})
