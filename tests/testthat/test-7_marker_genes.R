# calc_marker_genes / select_top_fold_genes_per_metacell against copies of the implementation they replaced

old_calc_marker_genes <- function(mc_egc,
                                  genes_per_metacell = 2,
                                  minimal_max_log_fraction = -13,
                                  minimal_relative_log_fraction = 2,
                                  fold_change_reg = 0.1,
                                  use_abs = TRUE) {
    mc_egc <- log2(mc_egc + 1e-5)

    max_log_fractions_of_genes <- matrixStats::rowMaxs(mc_egc)

    interesting_genes_mask <- (max_log_fractions_of_genes
    >= minimal_max_log_fraction)

    if (length(interesting_genes_mask) == 0) {
        cli_abort("No genes with at least one value above the {.field minimal_max_log_fraction} threshold ({.val {.value minimal_max_log_fraction}})")
    }

    mc_fp <- mc_egc - matrixStats::rowMedians(mc_egc, na.rm = TRUE)

    mc_top_genes <- old_select_top_fold_genes_per_metacell(
        mc_fp[interesting_genes_mask, , drop = FALSE],
        genes_per_metacell = genes_per_metacell,
        minimal_relative_log_fraction = minimal_relative_log_fraction,
        fold_change_reg = fold_change_reg,
        use_abs = use_abs
    )

    return(mc_top_genes)
}

old_select_top_fold_genes_per_metacell <- function(fold_matrix, genes_per_metacell = 2, minimal_relative_log_fraction = 2, fold_change_reg = 0.1, use_abs = TRUE) {
    fold_matrix <- fold_matrix + fold_change_reg
    fold_matrix[fold_matrix < minimal_relative_log_fraction] <- NA

    if (nrow(fold_matrix) == 0) {
        cli_abort("No genes with at least one value above the {.field minimal_relative_log_fraction} threshold ({.val {.value minimal_relative_log_fraction}})")
    }

    top_ind <- apply(fold_matrix, 2, function(fp) {
        if (use_abs) {
            order(-abs(fp))[1:genes_per_metacell]
        } else {
            order(-fp)[1:genes_per_metacell]
        }
    })
    top_ind <- matrix(top_ind, nrow = genes_per_metacell)

    genes <- rownames(fold_matrix)[top_ind]
    mc_top_genes <- tibble(
        metacell = rep(colnames(fold_matrix), each = genes_per_metacell),
        gene = genes,
        rank = rep(seq_len(genes_per_metacell), ncol(fold_matrix)),
        fp = stats::setNames(fold_matrix[cbind(as.vector(top_ind), as.vector(col(top_ind)))], genes)
    ) %>%
        arrange(metacell, rank)

    return(mc_top_genes)
}

# metacell fractions (egc) from sparse poisson counts: many all-zero genes, many genes below the thresholds, a few
# strongly enriched genes per metacell, and exact ties (the same counts in genes and in metacells)
synthetic_egc <- function(n_genes, n_metacells, seed, lambda = 0.3, umis = 2e4) {
    set.seed(seed)
    rates <- rexp(n_genes, 1 / lambda) * rbinom(n_genes, 1, 0.6)
    counts <- matrix(rpois(n_genes * n_metacells, rep(rates, n_metacells)), n_genes, n_metacells)
    # enriched genes in groups of metacells
    groups <- sample(1:5, n_metacells, replace = TRUE)
    for (g in 1:5) {
        genes <- sample(n_genes, min(n_genes, 8))
        counts[genes, groups == g] <- counts[genes, groups == g] + rpois(length(genes) * sum(groups == g), 30)
    }
    # exact ties: duplicated genes and metacells
    if (n_genes > 10) {
        counts[n_genes - 1, ] <- counts[1, ]
        counts[n_genes, ] <- 0
    }
    if (n_metacells > 3) {
        counts[, n_metacells] <- counts[, 1]
    }
    sizes <- umis / sample(c(1, 2), n_metacells, replace = TRUE)
    egc <- t(t(counts) / pmax(colSums(counts), 1)) * (colSums(counts) / sizes)
    dimnames(egc) <- list(paste0("g", seq_len(n_genes)), paste0("mc", seq_len(n_metacells)))
    egc
}

expect_same_markers <- function(new, old) {
    new <- tryCatch(new, error = function(e) conditionMessage(e))
    old <- tryCatch(old, error = function(e) conditionMessage(e))
    expect_identical(new, old)
    # also the names of fp (identical() does not compare the attributes of the columns of a tibble in all versions)
    if (is.data.frame(new) && is.data.frame(old)) {
        expect_identical(names(new$fp), names(old$fp))
    }
}

test_that("calc_marker_genes returns the same markers as before", {
    shapes <- list(c(200, 30), c(1000, 7), c(50, 120), c(12, 3), c(3, 40), c(500, 1), c(1, 9))
    for (shape in shapes) {
        egc <- synthetic_egc(shape[1], shape[2], seed = shape[1] * 7 + shape[2])
        for (k in c(1, 2, 5, 20)) {
            for (use_abs in c(TRUE, FALSE)) {
                expect_same_markers(
                    calc_marker_genes(egc, k, use_abs = use_abs),
                    old_calc_marker_genes(egc, k, use_abs = use_abs)
                )
            }
            # thresholds that leave metacells with fewer than k genes above them, or with none
            for (thresholds in list(c(-10, 4), c(-16, 0.5), c(-13, 9), c(-17, -Inf), c(-6, 1))) {
                expect_same_markers(
                    calc_marker_genes(egc, k, minimal_max_log_fraction = thresholds[1], minimal_relative_log_fraction = thresholds[2]),
                    old_calc_marker_genes(egc, k, minimal_max_log_fraction = thresholds[1], minimal_relative_log_fraction = thresholds[2])
                )
            }
        }
    }
})

test_that("calc_marker_genes returns the same markers as before on the import and heatmap inputs", {
    egc <- synthetic_egc(3000, 60, seed = 11)
    lateral <- seq_len(nrow(egc)) %% 7 == 0
    # import: the lateral genes are excluded before
    expect_same_markers(
        calc_marker_genes(egc[!lateral, ], 20, minimal_max_log_fraction = -13, minimal_relative_log_fraction = 2),
        old_calc_marker_genes(egc[!lateral, ], 20, minimal_max_log_fraction = -13, minimal_relative_log_fraction = 2)
    )
    # heatmap "Local": the visible metacells
    expect_same_markers(
        calc_marker_genes(egc[, c(5, 1, 30:40)], genes_per_metacell = 20),
        old_calc_marker_genes(egc[, c(5, 1, 30:40)], genes_per_metacell = 20)
    )
})

test_that("calc_marker_genes returns the same markers as before across blocks of genes", {
    egc <- synthetic_egc(600, 25, seed = 5)
    old <- old_calc_marker_genes(egc, 20)
    for (block_size in c(1, 2, 7, 50)) {
        local_mocked_bindings(marker_block_size = function(n_cols) block_size)
        expect_same_markers(calc_marker_genes(egc, 20), old)
        expect_same_markers(calc_marker_genes(egc, 3, use_abs = FALSE, minimal_relative_log_fraction = 6), old_calc_marker_genes(egc, 3, use_abs = FALSE, minimal_relative_log_fraction = 6))
    }
})

test_that("the heatmap's Local markers skip genes below the threshold without changing them", {
    prev <- mget(c("mc_data"), envir = globalenv(), ifnotfound = list(NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    set.seed(9)
    counts <- matrix(rpois(4000 * 50, rep(rexp(4000, 1 / 0.5) * rbinom(4000, 1, 0.5), 50)), 4000, 50, dimnames = list(paste0("g", 1:4000), paste0("mc", 1:50)))
    counts[, 7] <- counts[, 7] * 20 # metacells of different sizes
    assign("mc_data", list(toy = list2env(list(mc_mat = counts, mc_sum = colSums(counts)))), envir = globalenv())
    for (metacells in list(paste0("mc", 1:50), paste0("mc", c(30, 7, 2, 11)), "mc7", paste0("mc", 40:50))) {
        egc <- get_marker_candidates_egc("toy", metacells)
        expect_lt(nrow(egc), 4000)
        expect_same_markers(
            calc_marker_genes(egc, genes_per_metacell = 20),
            old_calc_marker_genes(get_mc_egc("toy", metacells = metacells), genes_per_metacell = 20)
        )
    }
})

test_that("calc_marker_genes handles missing values, integers and unnamed genes as before", {
    egc <- synthetic_egc(80, 10, seed = 3)
    egc_na <- egc
    egc_na[5, 3] <- NA
    egc_na[1, ] <- NA
    egc_na[60, 2] <- NaN
    expect_same_markers(suppressWarnings(calc_marker_genes(egc_na, 20)), suppressWarnings(old_calc_marker_genes(egc_na, 20)))
    expect_same_markers(suppressWarnings(calc_marker_genes(egc_na, 3)), suppressWarnings(old_calc_marker_genes(egc_na, 3)))

    counts <- round(egc * 1e4)
    storage.mode(counts) <- "integer"
    expect_same_markers(calc_marker_genes(counts, 4, minimal_max_log_fraction = 3), old_calc_marker_genes(counts, 4, minimal_max_log_fraction = 3))

    unnamed <- unname(egc)
    colnames(unnamed) <- colnames(egc)
    expect_identical(calc_marker_genes(unnamed, 5), old_calc_marker_genes(unnamed, 5))
})

test_that("calc_marker_genes fails as before without genes above the thresholds", {
    egc <- synthetic_egc(50, 10, seed = 2)
    expect_error(calc_marker_genes(egc, 5, minimal_max_log_fraction = 10), "minimal_relative_log_fraction")
    expect_error(old_calc_marker_genes(egc, 5, minimal_max_log_fraction = 10), "minimal_relative_log_fraction")
    expect_error(calc_marker_genes(egc[0, ], 5), "minimal_max_log_fraction")
    expect_error(old_calc_marker_genes(egc[0, ], 5), "minimal_max_log_fraction")
})

test_that("select_top_fold_genes_per_metacell returns the same genes as before", {
    set.seed(17)
    # projected folds, as in atlas projection (all the folds, by absolute value), with ties and zeros
    fold <- matrix(round(rnorm(400 * 30), 1), 400, 30, dimnames = list(paste0("g", 1:400), paste0("q", 1:30)))
    fold[1:50, ] <- 0
    fold[sample(length(fold), 100)] <- -fold[sample(length(fold), 100)]
    for (k in c(1, 2, 50, 450)) {
        for (use_abs in c(TRUE, FALSE)) {
            for (threshold in c(-Inf, 0, 1.5, 3)) {
                expect_same_markers(
                    select_top_fold_genes_per_metacell(fold, minimal_relative_log_fraction = threshold, use_abs = use_abs, genes_per_metacell = k),
                    old_select_top_fold_genes_per_metacell(fold, minimal_relative_log_fraction = threshold, use_abs = use_abs, genes_per_metacell = k)
                )
            }
        }
    }
    local_mocked_bindings(marker_block_size = function(n_cols) 13)
    fold[3, 4] <- NA
    fold[300, ] <- NaN
    expect_same_markers(
        select_top_fold_genes_per_metacell(fold, minimal_relative_log_fraction = -Inf, use_abs = TRUE, genes_per_metacell = 50),
        old_select_top_fold_genes_per_metacell(fold, minimal_relative_log_fraction = -Inf, use_abs = TRUE, genes_per_metacell = 50)
    )
    expect_same_markers(
        select_top_fold_genes_per_metacell(fold, genes_per_metacell = 20),
        old_select_top_fold_genes_per_metacell(fold, genes_per_metacell = 20)
    )
    # blocks with more than k rows and many folds above the threshold, with missing values in them
    for (use_abs in c(TRUE, FALSE)) {
        for (threshold in c(-Inf, -1, 0.5)) {
            expect_same_markers(
                select_top_fold_genes_per_metacell(fold, minimal_relative_log_fraction = threshold, use_abs = use_abs, genes_per_metacell = 3),
                old_select_top_fold_genes_per_metacell(fold, minimal_relative_log_fraction = threshold, use_abs = use_abs, genes_per_metacell = 3)
            )
        }
    }
    fold[10, ] <- -Inf
    fold[11, 1:5] <- Inf
    expect_same_markers(
        select_top_fold_genes_per_metacell(fold, minimal_relative_log_fraction = -Inf, use_abs = FALSE, genes_per_metacell = 4),
        old_select_top_fold_genes_per_metacell(fold, minimal_relative_log_fraction = -Inf, use_abs = FALSE, genes_per_metacell = 4)
    )
})
