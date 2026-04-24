test_that("streaming_top_k_cor matches tgs_cross_cor_blas + tgs_knn", {
    skip_if_not(requireNamespace("tgstat", quietly = TRUE))
    set.seed(1)
    G <- 500L; M <- 50L; k <- 10L
    X <- matrix(rnorm(G * M), nrow = G, ncol = M)
    X <- X + matrix(rnorm(G, sd = 0.01), nrow = G, ncol = M)
    rownames(X) <- paste0("g", seq_len(G))
    colnames(X) <- paste0("m", seq_len(M))

    # Legacy path: full cor + tgs_knn (pos) + tgs_knn(-cor) (neg)
    cm <- tgstat::tgs_cor(t(X), y = t(X),
                         pairwise.complete.obs = TRUE, spearman = FALSE)
    options(tgs_max.processes = 1L)
    pos <- tgstat::tgs_knn(cm, knn = k, diag = FALSE)
    names(pos)[names(pos) == "col1"] <- "gene1"
    names(pos)[names(pos) == "col2"] <- "gene2"
    names(pos)[names(pos) == "val"] <- "cor"
    pos$type <- "pos"

    neg <- tgstat::tgs_knn(-cm, knn = k, diag = FALSE)
    names(neg)[names(neg) == "col1"] <- "gene1"
    names(neg)[names(neg) == "col2"] <- "gene2"
    names(neg)[names(neg) == "val"] <- "cor"
    neg$cor <- -neg$cor
    neg$type <- "neg"

    legacy <- rbind(pos, neg)
    legacy <- legacy[!is.na(legacy$cor), ]

    # New kernel
    new <- MCView::streaming_top_k_cor(
        X = X, k = k, diag = FALSE, min_cor = 0
    )

    # Same number of (type, gene1, gene2) triples
    key_new <- paste(new$type, new$gene1, new$gene2)
    key_leg <- paste(legacy$type, legacy$gene1, legacy$gene2)
    expect_setequal(key_new, key_leg)

    # Values match per key (allow tolerance)
    merged <- merge(
        new, legacy,
        by = c("type", "gene1", "gene2"),
        suffixes = c(".new", ".leg")
    )
    expect_equal(nrow(merged), nrow(new))
    expect_equal(merged$cor.new, merged$cor.leg, tolerance = 1e-8)
})

test_that("streaming_top_k_cor respects min_cor threshold", {
    set.seed(2)
    X <- matrix(rnorm(200 * 30), nrow = 200, ncol = 30)
    rownames(X) <- paste0("g", 1:200)
    colnames(X) <- paste0("m", 1:30)
    out <- MCView::streaming_top_k_cor(X = X, k = 50L, diag = FALSE, min_cor = 0.5)
    # All returned cors should have |cor| >= 0.5
    expect_true(all(abs(out$cor) >= 0.5))
})

test_that("streaming_top_k_cor is deterministic", {
    set.seed(3)
    X <- matrix(rnorm(300 * 40), nrow = 300, ncol = 40)
    rownames(X) <- paste0("g", 1:300)
    colnames(X) <- paste0("m", 1:40)
    r1 <- MCView::streaming_top_k_cor(X, k = 5L, diag = FALSE)
    r2 <- MCView::streaming_top_k_cor(X, k = 5L, diag = FALSE)
    # Same (type, gene1, gene2) set
    expect_setequal(paste(r1$type, r1$gene1, r1$gene2),
                    paste(r2$type, r2$gene1, r2$gene2))
})
