test_that("top_fold_per_col matches legacy select_top_fold_genes_per_metacell", {
    set.seed(42)
    n_genes <- 200
    n_mcs <- 30
    fm <- matrix(rnorm(n_genes * n_mcs), nrow = n_genes, ncol = n_mcs)
    fm[sample(length(fm), 200L)] <- NA_real_
    rownames(fm) <- paste0("g", seq_len(n_genes))
    colnames(fm) <- paste0("mc", seq_len(n_mcs))

    k <- 5L
    min_val <- 2
    fold_reg <- 0.1
    use_abs <- TRUE

    legacy <- MCView:::select_top_fold_genes_per_metacell(
        fm, genes_per_metacell = k,
        minimal_relative_log_fraction = min_val,
        fold_change_reg = fold_reg,
        use_abs = use_abs
    )

    new <- MCView::top_fold_per_col(
        fp = fm, k = k, min_val = min_val,
        fold_reg = fold_reg, use_abs = use_abs
    )

    # Drop NA-padded rows from both sides. Legacy pads with arbitrary gene
    # names when fewer than k values pass the min_val mask (the gene slot is
    # filled from row_names in index order, but fp is NA -- so gene is junk).
    # The kernel emits NA gene for those same padding positions. Drop fp==NA
    # from both so we compare only the rows with real top-K matches.
    legacy <- legacy[!is.na(legacy$fp), , drop = FALSE]
    new <- new[!is.na(new$fp), , drop = FALSE]

    # Sort both by (metacell, rank) for a stable compare.
    legacy <- legacy[order(legacy$metacell, legacy$rank), ]
    new <- new[order(new$metacell, new$rank), ]

    expect_equal(nrow(new), nrow(legacy))
    expect_equal(as.character(new$metacell), as.character(legacy$metacell))
    expect_equal(as.character(new$gene), as.character(legacy$gene))
    expect_equal(as.integer(new$rank), as.integer(legacy$rank))
    expect_equal(unname(new$fp), unname(legacy$fp), tolerance = 1e-10)
})

test_that("top_fold_per_col handles k > n_rows (pads with NA gene/fp)", {
    fm <- matrix(runif(6, 2, 5), nrow = 3, ncol = 2)
    rownames(fm) <- c("g1", "g2", "g3")
    colnames(fm) <- c("m1", "m2")
    out <- MCView::top_fold_per_col(fm, k = 10L, min_val = 0, fold_reg = 0, use_abs = FALSE)
    expect_equal(nrow(out), 2L * 10L)
    # 3 genes qualify per metacell (all values > 0), rest are NA
    expect_equal(sum(!is.na(out$gene)), 2L * 3L)
    # rank column is always 1..10
    expect_equal(sort(unique(out$rank)), 1:10)
})

test_that("top_fold_per_col is deterministic across calls", {
    set.seed(7)
    fm <- matrix(rnorm(100), nrow = 20, ncol = 5)
    rownames(fm) <- paste0("g", 1:20)
    colnames(fm) <- paste0("m", 1:5)
    r1 <- MCView::top_fold_per_col(fm, k = 3L, min_val = -100, fold_reg = 0, use_abs = FALSE)
    r2 <- MCView::top_fold_per_col(fm, k = 3L, min_val = -100, fold_reg = 0, use_abs = FALSE)
    expect_equal(r1, r2)
})
