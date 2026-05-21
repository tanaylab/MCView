make_lfp_fixture <- function(n_genes = 40L, n_mcs = 60L, seed = 17L) {
    set.seed(seed)
    egc <- matrix(abs(rnorm(n_genes * n_mcs, mean = 1, sd = 0.5)),
                  nrow = n_genes,
                  dimnames = list(paste0("g", seq_len(n_genes)),
                                  paste0("m", seq_len(n_mcs))))
    # Normalise columns so they sum to ~1 (looks like a fraction matrix).
    egc <- sweep(egc, 2, colSums(egc), "/")

    d <- dafr::memory_daf("calc_top_cors_fixture")
    dafr::add_axis(d, "gene", rownames(egc))
    dafr::add_axis(d, "metacell", colnames(egc))
    dafr::set_matrix(d, "gene", "metacell", "linear_fraction", egc)
    # mc_sum doesn't need to be exact for calc_top_cors; provide it for
    # the get_mc_lfp path that goes via fast_log.
    dafr::set_vector(d, "metacell", "total_UMIs", rep(1.0, n_mcs))
    # UMIs is required by some converters (defensive).
    dafr::set_matrix(d, "gene", "metacell", "UMIs", round(egc * 1000))
    d
}

setup_calc_top_cors_dataset <- function() {
    d <- make_lfp_fixture()
    init_mcview_env()
    prev_mc_data <- mcv_get("mc_data")
    prev_eps <- mcv_get("egc_epsilon")
    mcv_set("egc_epsilon", 1e-5)
    mcv_set("mc_data", list(test_ds = list(daf_obj = d)))
    list(daf = d, restore = function() {
        mcv_set("mc_data", prev_mc_data)
        mcv_set("egc_epsilon", prev_eps)
    })
}

test_that("calc_top_cors: data_vec without filter matches gene-self path", {
    ctx <- setup_calc_top_cors_dataset()
    on.exit(ctx$restore(), add = TRUE)

    lfp <- get_mc_lfp("test_ds")
    g1 <- rownames(lfp)[1]

    self <- calc_top_cors("test_ds", g1, "pos", NULL, NULL, NULL, FALSE)
    via_vec <- calc_top_cors("test_ds", "module", "pos",
                             lfp[g1, ], NULL, NULL, FALSE)

    # Same correlation set in same order.
    expect_equal(self$gene2, via_vec$gene2)
    expect_equal(self$cor, via_vec$cor, tolerance = 1e-12)
})

test_that("calc_top_cors: data_vec overrides metacell_filter", {
    ctx <- setup_calc_top_cors_dataset()
    on.exit(ctx$restore(), add = TRUE)

    lfp <- get_mc_lfp("test_ds")
    g1 <- rownames(lfp)[1]
    mcf <- colnames(lfp)[1:20]    # arbitrary subset

    via_vec_only <- calc_top_cors("test_ds", "module", "pos",
                                  lfp[g1, ], NULL, NULL, FALSE)
    via_vec_with_filter <- calc_top_cors("test_ds", "module", "pos",
                                         lfp[g1, ], mcf, NULL, FALSE)

    # data_vec semantically overrides metacell_filter; both calls run on
    # the full metacell axis. Outputs must match.
    expect_equal(via_vec_only$gene2, via_vec_with_filter$gene2)
    expect_equal(via_vec_only$cor, via_vec_with_filter$cor,
                 tolerance = 1e-12)
})

test_that("calc_top_cors: exclude drops excluded genes from the result", {
    ctx <- setup_calc_top_cors_dataset()
    on.exit(ctx$restore(), add = TRUE)

    lfp <- get_mc_lfp("test_ds")
    g1 <- rownames(lfp)[1]
    exc <- rownames(lfp)[2:6]    # exclude these from the cor output

    # Run with exclude=NULL to discover what the unfiltered top-30 looks like
    # AFTER post-filtering excluded genes ourselves; calc_top_cors does this
    # as a post-filter step in get_top_cor_gene, but inside calc_top_cors
    # itself the exclude is applied BEFORE the cor (row-subsets lfp).
    df_excl <- calc_top_cors("test_ds", g1, "pos", NULL, NULL, exc, FALSE)
    expect_false(any(df_excl$gene2 %in% exc))
})

test_that("calc_top_cors: metacell_filter restricts cor to filtered metacells", {
    ctx <- setup_calc_top_cors_dataset()
    on.exit(ctx$restore(), add = TRUE)

    lfp <- get_mc_lfp("test_ds")
    g1 <- rownames(lfp)[1]
    mcf <- colnames(lfp)[1:30]

    # Reference: manually compute the same cor on the filtered subset.
    egc_sub <- get_mc_egc("test_ds", metacells = mcf)
    lfp_sub <- dafr::fast_log(egc_sub, eps = 1e-5, base = 2)
    g1_vec <- lfp_sub[g1, ]
    ref_cors <- apply(lfp_sub, 1, function(r) cor(r, g1_vec))
    ref_top_pos <- names(sort(ref_cors, decreasing = TRUE)[1:30])

    out <- calc_top_cors("test_ds", g1, "pos", NULL, mcf, NULL, FALSE)
    expect_equal(out$gene2[1:5], ref_top_pos[1:5])
})

test_that("calc_top_cors: lfp_full_t is invalidated when lfp_full is cleared", {
    ctx <- setup_calc_top_cors_dataset()
    on.exit(ctx$restore(), add = TRUE)

    lfp <- get_mc_lfp("test_ds")
    lfp_t <- get_mc_lfp_t("test_ds")
    expect_equal(dim(lfp_t), rev(dim(lfp)))
    expect_equal(rownames(lfp_t), colnames(lfp))
    expect_equal(colnames(lfp_t), rownames(lfp))

    # Externally evict lfp_full. The cached lfp_full_t is now stale - it
    # belongs to a matrix that no longer exists in the cache. The next
    # read MUST drop the stale transpose and rebuild from scratch (which
    # implies rebuilding lfp_full too, via the get_mc_lfp() chain).
    mcv_cache_set("test_ds", "lfp_full", NULL)
    # Bug-detection: if get_mc_lfp_t doesn't notice lfp_full is gone, it
    # short-circuits on the cached lfp_full_t and never repopulates
    # lfp_full. The cleaner-than-digest signal that the staleness guard
    # fired is "lfp_full was rebuilt".
    expect_false(is.null(get_mc_lfp_t("test_ds")))
    expect_false(is.null(mcv_cache_get("test_ds", "lfp_full")),
                 info = "lfp_full must be repopulated by the staleness rebuild path")

    # Sanity: the rebuilt transpose still has the right shape/labels
    # (data is deterministic so values match the pre-eviction version).
    rebuilt_t <- mcv_cache_get("test_ds", "lfp_full_t")
    expect_equal(dim(rebuilt_t), dim(lfp_t))
    expect_equal(rownames(rebuilt_t), rownames(lfp_t))
})
