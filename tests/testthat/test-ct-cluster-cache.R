test_that("ct_cluster_cache put/touch/evict behaves as LRU", {
    on.exit(clear_ct_cluster_cache(), add = TRUE)
    clear_ct_cluster_cache()

    # Shrink the cap so we can exercise eviction in O(1) calls.
    cap <- .ct_cluster_cache_max
    assignInNamespace(".ct_cluster_cache_max", 3L, ns = "MCView")
    on.exit(assignInNamespace(".ct_cluster_cache_max", cap, ns = "MCView"),
            add = TRUE)

    .ct_cluster_cache_put("a", 1L)
    .ct_cluster_cache_put("b", 2L)
    .ct_cluster_cache_put("c", 3L)
    expect_setequal(ls(.ct_cluster_cache), c("a", "b", "c"))
    expect_identical(.ct_cluster_cache_meta$keys, c("a", "b", "c"))

    # Insert past the cap -> oldest ("a") is evicted, others stay.
    .ct_cluster_cache_put("d", 4L)
    expect_setequal(ls(.ct_cluster_cache), c("b", "c", "d"))
    expect_identical(.ct_cluster_cache_meta$keys, c("b", "c", "d"))
    expect_false(exists("a", envir = .ct_cluster_cache, inherits = FALSE))

    # Touch "b" -> it moves to the tail; next eviction takes "c", not "b".
    .ct_cluster_cache_touch("b")
    expect_identical(.ct_cluster_cache_meta$keys, c("c", "d", "b"))

    .ct_cluster_cache_put("e", 5L)
    expect_setequal(ls(.ct_cluster_cache), c("d", "b", "e"))
    expect_identical(.ct_cluster_cache_meta$keys, c("d", "b", "e"))
    expect_false(exists("c", envir = .ct_cluster_cache, inherits = FALSE))

    # Re-put existing key updates value AND moves to tail (no eviction).
    .ct_cluster_cache_put("d", 40L)
    expect_identical(get("d", envir = .ct_cluster_cache), 40L)
    expect_identical(.ct_cluster_cache_meta$keys, c("b", "e", "d"))
    expect_setequal(ls(.ct_cluster_cache), c("b", "d", "e"))
})

test_that("clear_ct_cluster_cache wipes both values and key order", {
    on.exit(clear_ct_cluster_cache(), add = TRUE)

    .ct_cluster_cache_put("x", 1L)
    .ct_cluster_cache_put("y", 2L)
    expect_gt(length(.ct_cluster_cache_meta$keys), 0L)

    clear_ct_cluster_cache()
    expect_identical(.ct_cluster_cache_meta$keys, character(0))
    expect_identical(ls(.ct_cluster_cache), character(0))
})
