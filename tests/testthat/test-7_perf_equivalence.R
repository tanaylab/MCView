# Faster implementations against copies of the code they replaced, on synthetic data

# A toy dataset in the global mc_data, as the app keeps it (an environment per dataset), until the calling test ends
local_toy_dataset <- function(mat, ..., env = parent.frame()) {
    prev <- mget(c("mc_data", "egc_epsilon"), envir = globalenv(), ifnotfound = list(NULL, NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())), envir = env)
    assign("mc_data", list(toy = list2env(list(mc_mat = mat, mc_sum = colSums(mat), ...))), envir = globalenv())
    assign("egc_epsilon", 1e-5, envir = globalenv())
}

toy_mc_mat <- function(n_genes = 300, n_metacells = 40, seed = 60427) {
    set.seed(seed)
    mat <- matrix(rpois(n_genes * n_metacells, 2), n_genes, n_metacells, dimnames = list(paste0("g", seq_len(n_genes)), paste0("mc", seq_len(n_metacells))))
    # genes correlated with g1 and g2, so that there are strong correlations in both directions
    mat[3:12, ] <- mat[3:12, ] + 3 * mat[rep(1, 10), ]
    mat[13:22, ] <- pmax(0, mat[13:22, ] + 5 - 2 * mat[rep(2, 10), ])
    mat[7, ] <- 0 # no counts: no correlations
    mat
}

test_that("calc_gene_gene_correlations matches the implementation with mapply", {
    # the implementation before: a lookup of the correlation matrix per gene pair
    old_calc_gene_gene_correlations <- function(dataset, genes, cell_type_filter = NULL, threshold = 0, atlas = FALSE, metacell_types = NULL) {
        mc_egc <- get_mc_egc(dataset, genes = genes, atlas = atlas)
        available_genes <- rownames(mc_egc)
        valid_genes <- genes[genes %in% available_genes]
        if (length(valid_genes) < 2) {
            stop("Need at least 2 valid genes for gene-gene correlation")
        }
        metacell_filter <- NULL
        if (!is.null(cell_type_filter)) {
            metacell_types <- metacell_types %||% get_mc_data(dataset, "metacell_types")
            if (!is.null(metacell_types)) {
                metacell_filter <- metacell_types %>%
                    filter(cell_type %in% cell_type_filter) %>%
                    pull(metacell)
            }
        }
        expr_data <- mc_egc[valid_genes, , drop = FALSE]
        if (!is.null(metacell_filter)) {
            metacell_filter <- intersect(metacell_filter, colnames(expr_data))
            expr_data <- expr_data[, metacell_filter, drop = FALSE]
        }
        cor_matrix <- cor(t(expr_data), use = "pairwise.complete.obs")
        expand.grid(input_gene = valid_genes, gene2 = valid_genes, stringsAsFactors = FALSE) %>%
            filter(input_gene != gene2) %>%
            mutate(
                cor = mapply(function(g1, g2) cor_matrix[g1, g2], input_gene, gene2),
                correlation_type = ifelse(cor > 0, "positive", "negative")
            ) %>%
            filter(abs(cor) >= threshold) %>%
            mutate(rank = ave(abs(cor), input_gene, FUN = function(x) rank(-x))) %>%
            arrange(input_gene, desc(abs(cor)))
    }
    # mapply named its result by input_gene, and cor, correlation_type and rank (computed from cor) inherited the
    # names. Nothing reads them, and the new columns have none.
    unname_columns <- function(df) {
        df[] <- lapply(df, unname)
        df
    }

    mat <- toy_mc_mat()
    local_toy_dataset(mat)
    metacell_types <- tibble(metacell = colnames(mat), cell_type = rep(c("A", "B", "C"), length.out = ncol(mat)))
    genes <- c("g5", "g1", "g13", "g2", "g7", "g40", "g3", "g300", "not_a_gene")
    suppressWarnings({ # the gene without counts has no standard deviation
        cases <- list(
            list(genes),
            list(genes, threshold = 0.3),
            list(genes, cell_type_filter = c("A", "C"), metacell_types = metacell_types),
            list(c("g1", "g2"))
        )
        for (args in cases) {
            old <- do.call(old_calc_gene_gene_correlations, c(list("toy"), args))
            new <- do.call(calc_gene_gene_correlations, c(list("toy"), args))
            expect_gt(nrow(new), 0)
            expect_identical(new, unname_columns(old))
        }
    })
})
