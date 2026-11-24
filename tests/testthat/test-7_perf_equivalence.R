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

# the implementations before: a dplyr pipeline per input gene, on the top 30 correlations in each direction
old_select_top_cors <- function(cors, type) {
    if (type == "pos") {
        cors <- head(sort(cors, decreasing = TRUE), n = 30)
    } else if (type == "neg") {
        cors <- head(sort(cors), n = 30)
    } else {
        cors <- c(head(sort(cors, decreasing = TRUE), n = 30), head(sort(cors), n = 30))
    }
    tibble(gene2 = names(cors) %||% character(0), cor = unname(cors))
}

old_calc_individual_correlations <- function(dataset, genes, n_top = 30, cell_type_filter = NULL, threshold = 0, atlas = FALSE, metacell_types = NULL) {
    available_genes <- gene_names(dataset, atlas = atlas)
    valid_genes <- genes[genes %in% available_genes]
    if (length(valid_genes) == 0) {
        stop("No valid genes found in dataset")
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
    mc_mat <- get_mc_data(dataset, "mc_mat", atlas = atlas)
    metacells <- colnames(mc_mat)
    if (!is.null(metacell_filter)) {
        req(metacell_filter %in% metacells)
        metacells <- metacell_filter
    }
    x <- log2(t(mc_mat[valid_genes, metacells, drop = FALSE]) / get_mc_sum(dataset, atlas = atlas)[metacells] + egc_epsilon)
    all_cors <- cor_log_egc(dataset, x, rownames(mc_mat), metacells, atlas = atlas)

    purrr::map_dfr(valid_genes, function(gene) {
        cors <- old_select_top_cors(all_cors[, gene], "both") %>%
            filter(abs(cor) >= threshold)
        pos_cors <- cors %>%
            filter(cor > 0) %>%
            arrange(desc(cor)) %>%
            head(n_top)
        neg_cors <- cors %>%
            filter(cor < 0) %>%
            arrange(cor) %>%
            head(n_top)
        bind_rows(pos_cors, neg_cors) %>%
            arrange(desc(abs(cor))) %>%
            head(n_top * 2) %>%
            mutate(
                input_gene = gene,
                rank = row_number(),
                correlation_type = if_else(cor > 0, "positive", "negative")
            )
    })
}

test_that("calc_individual_correlations matches the dplyr implementation, and n_top bounds what it keeps", {
    mat <- toy_mc_mat()
    # ties, which keep the order of the genes
    mat["g250", ] <- mat["g251", ] <- mat["g252", ] <- mat["g4", ]
    local_toy_dataset(mat)
    metacell_types <- tibble(metacell = colnames(mat), cell_type = rep(c("A", "B", "C"), length.out = ncol(mat)))
    genes <- c("g1", "g7", "g2", "g4", "g100", "not_a_gene", "g1")

    cases <- list(
        list(genes),
        list(genes, n_top = 5),
        list(genes, n_top = 10, threshold = 0.3),
        list(genes, n_top = 30, threshold = 0.5),
        list(genes, n_top = 12, cell_type_filter = c("A", "C"), metacell_types = metacell_types),
        list("g2", n_top = 7.5),
        list("g7"), # only a gene without counts: no rows
        list(genes, threshold = 0.999) # no rows
    )
    for (args in cases) {
        old <- do.call(old_calc_individual_correlations, c(list("toy"), args))
        new <- do.call(calc_individual_correlations, c(list("toy"), args))
        expect_identical(new, old)
    }
    expect_identical(nrow(calc_individual_correlations("toy", "g7")), 0L)
    expect_gt(nrow(calc_individual_correlations("toy", genes, threshold = 0.5)), 0)

    # more than 30: n_top per direction (before, the 30 of select_top_cors at most)
    old <- old_calc_individual_correlations("toy", c("g1", "g2"), n_top = 30)
    new <- calc_individual_correlations("toy", c("g1", "g2"), n_top = 50)
    expect_equal(as.vector(table(new$input_gene, new$correlation_type)), rep(50L, 4))
    expect_identical(dplyr::semi_join(new, old, by = c("input_gene", "gene2")) %>% select(-rank), old %>% select(-rank))
    # the 50 highest and lowest of all the correlations
    all_cors <- suppressWarnings(cor(log2(t(mat) / colSums(mat) + 1e-5)))[, "g1"]
    g1 <- new[new$input_gene == "g1", ]
    expect_equal(sort(g1$cor[g1$cor > 0], decreasing = TRUE), unname(head(sort(all_cors, decreasing = TRUE), 50)))
    expect_equal(sort(g1$cor[g1$cor < 0]), unname(head(sort(all_cors), 50)))
})

test_that("calc_module_correlations is unchanged for n_top <= 30, and returns n_top per direction above", {
    # the implementation before: always the top 30 of each direction (calc_top_cors' default)
    old_calc_module_correlations <- function(dataset, genes, n_top = 30, threshold = 0) {
        valid_genes <- genes[genes %in% rownames(get_mc_data(dataset, "mc_mat"))]
        module_expr <- colMeans(get_mc_egc(dataset, genes = valid_genes), na.rm = TRUE)
        cors <- calc_top_cors(dataset, "module", "both", module_expr, NULL, valid_genes)
        gene_display <- paste(valid_genes, collapse = ",")
        cors %>%
            filter(abs(cor) >= threshold) %>%
            arrange(desc(abs(cor))) %>%
            head(n_top * 2) %>%
            mutate(input_gene = paste0("Module(", gene_display, ")"), rank = row_number(), correlation_type = ifelse(cor > 0, "positive", "negative"))
    }
    mat <- toy_mc_mat()
    local_toy_dataset(mat)
    genes <- c("g1", "g3", "g4")
    for (n_top in c(5, 20, 30)) {
        for (threshold in c(0, 0.4)) {
            expect_identical(calc_module_correlations("toy", genes, n_top = n_top, threshold = threshold), old_calc_module_correlations("toy", genes, n_top = n_top, threshold = threshold))
        }
    }
    new <- calc_module_correlations("toy", genes, n_top = 50)
    expect_equal(as.vector(table(new$correlation_type)), c(50L, 50L))
})
