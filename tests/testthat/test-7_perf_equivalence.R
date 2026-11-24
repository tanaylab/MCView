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

test_that("the stored top correlations of a gene are the same through the index of gg_mc_top_cor", {
    # the implementation before: filter the whole table for every lookup
    old_get_cached_cor <- function(gg_mc_top_cor, gene, type, exclude = NULL) {
        df <- gg_mc_top_cor %>%
            filter(gene1 == gene)
        if (type %in% c("pos", "neg")) {
            df <- df %>%
                filter(type == !!type)
        }
        df <- df %>%
            filter(gene1 != gene2)
        if (!is.null(exclude)) {
            df <- df %>%
                filter(!(gene2 %in% exclude))
        }
        return(df)
    }
    old_get_top_cor_gene <- function(dataset, gene, type = "pos", atlas = FALSE, data_vec = NULL, exclude = NULL, metacell_filter = NULL) {
        gg_mc_top_cor <- get_mc_data(dataset, "gg_mc_top_cor", atlas = atlas)
        if (is.null(data_vec) && is.null(metacell_filter) && gene %in% gg_mc_top_cor$gene1) {
            df <- old_get_cached_cor(gg_mc_top_cor, gene, type, exclude = exclude)
        } else {
            df <- calc_top_cors(dataset, gene, type, data_vec, metacell_filter, exclude, atlas = atlas)
        }
        df %>%
            mutate(
                gene2 = as.character(gene2),
                gene_name = gene_label(gene2, dataset),
                label = glue("{gene_name} ({cr})", cr = round(cor, digits = 2))
            ) %>%
            select(label, gene2) %>%
            tibble::deframe()
    }

    mat <- toy_mc_mat()
    genes <- rownames(mat)
    set.seed(5)
    # in no particular order (calc_gg_mc_top_cor doesn't sort it), with self-correlations, and without g300 as gene1
    gg <- tibble(gene1 = rep(genes[-300], each = 8), gene2 = sample(genes, 299 * 8, replace = TRUE), cor = runif(299 * 8, -1, 1)) %>%
        mutate(type = ifelse(cor > 0, "pos", "neg"))
    gg$gene2[c(1, 50, 51)] <- gg$gene1[c(1, 50, 51)]
    gg <- gg[sample(nrow(gg)), ]
    lookups <- c("g1", "g6", "g7", "g120", "g299", "g300", "not_a_gene")
    for (gg_mc_top_cor in list(gg, as.data.frame(gg), mutate(gg, gene1 = factor(gene1, levels = genes), gene2 = factor(gene2, levels = genes)))) {
        local_toy_dataset(mat, gg_mc_top_cor = gg_mc_top_cor, lateral_genes = c("g3", "g10"), noisy_genes = "g10")
        for (gene in lookups) {
            for (type in c("pos", "neg", "both")) {
                for (exclude in list(NULL, c("g2", "g5", gg$gene2[gg$gene1 == gene][1:2]))) {
                    rows <- get_gg_mc_top_cor_index("toy")[[gene]]
                    expect_identical(get_cached_cor(gg_mc_top_cor, gene, type, exclude = exclude, rows = rows), old_get_cached_cor(gg_mc_top_cor, gene, type, exclude = exclude))
                    expect_identical(get_cached_cor(gg_mc_top_cor, gene, type, exclude = exclude), old_get_cached_cor(gg_mc_top_cor, gene, type, exclude = exclude))
                    if (gene != "not_a_gene") { # calc_top_cors() req()s a gene of the dataset
                        expect_identical(get_top_cor_gene("toy", gene, type = type, exclude = exclude), old_get_top_cor_gene("toy", gene, type = type, exclude = exclude))
                    }
                }
            }
        }
        # the index is kept with the dataset
        expect_true(exists("gg_mc_top_cor_index", envir = mc_data$toy, inherits = FALSE))
    }

    # without stored correlations
    local_toy_dataset(mat)
    expect_null(get_gg_mc_top_cor_index("toy"))
    expect_identical(get_top_cor_gene("toy", "g1", type = "both"), old_get_top_cor_gene("toy", "g1", type = "both"))
})

test_that("choose_markers matches the implementation with group_by() and slice()", {
    old_choose_markers <- function(marker_genes, max_markers, dataset = NULL) {
        if (has_name(marker_genes, "metacell")) {
            n_metacells <- length(unique(marker_genes$metacell))
            genes_per_metacell <- min(10, max(1, round(max_markers / n_metacells)))
            markers <- marker_genes %>%
                group_by(metacell) %>%
                slice(1:genes_per_metacell) %>%
                pull(gene) %>%
                unique()
        } else {
            markers <- unique(marker_genes$gene)
        }
        if (length(markers) > max_markers) {
            markers <- markers[1:max_markers]
        } else if (length(markers) < max_markers) {
            more_genes <- marker_genes %>%
                filter(!(gene %in% markers)) %>%
                arrange(desc(fp)) %>%
                slice(1:(max_markers - length(markers))) %>%
                pull(gene) %>%
                unique()
            markers <- c(markers, more_genes)
        }
        markers <- sort(unique(markers))
        return(markers)
    }

    set.seed(17)
    fold <- matrix(rnorm(500 * 60, sd = 2), 500, 60, dimnames = list(paste0("g", 1:500), paste0("M", 1:60)))
    marker_genes <- select_top_fold_genes_per_metacell(fold, genes_per_metacell = 20, minimal_relative_log_fraction = 0)
    # metacells whose order differs between the C locale (group_by) and others, and fewer genes for some metacells
    marker_genes$metacell[marker_genes$metacell == "M1"] <- "m1"
    marker_genes$metacell[marker_genes$metacell == "M2"] <- "a2"
    marker_genes <- marker_genes[!(marker_genes$metacell == "M3" & marker_genes$rank > 2), ]
    variants <- list(
        sorted = marker_genes,
        shuffled = marker_genes[sample(nrow(marker_genes)), ], # the rows of each metacell in no particular order
        factor = mutate(marker_genes[sample(nrow(marker_genes)), ], metacell = factor(metacell, levels = sample(unique(metacell)))),
        no_rank = select(marker_genes[sample(nrow(marker_genes)), ], -rank),
        no_metacell = tibble(gene = sample(rownames(fold), 300, replace = TRUE), fp = runif(300))
    )
    for (variant in names(variants)) {
        for (max_markers in c(5, 30, 60, 100, 300, 1000)) {
            expect_identical(choose_markers(variants[[variant]], max_markers), old_choose_markers(variants[[variant]], max_markers), label = paste(variant, max_markers))
        }
    }
})

test_that("sparse matrices stored as triplets are read as compressed columns, with the same values", {
    prev <- mget("mc_data", envir = globalenv(), ifnotfound = list(NULL))
    withr::defer(purrr::iwalk(prev, ~ assign(.y, .x, envir = globalenv())))
    cache_dir <- withr::local_tempdir()
    fs::dir_create(fs::path(cache_dir, "ds"))
    set.seed(3)
    dims <- list(paste0("g", 1:50), paste0("M", 1:20))
    sparse <- Matrix::rsparsematrix(50, 20, density = 0.1, dimnames = dims) # dgCMatrix
    triplets <- methods::as(sparse, "TsparseMatrix") # as projects imported with older versions store them
    dense_triplets <- methods::as(methods::as(Matrix::Matrix(matrix(runif(1000), 50, 20, dimnames = dims), sparse = TRUE), "generalMatrix"), "TsparseMatrix")
    rows <- methods::as(sparse, "RsparseMatrix")
    expect_s4_class(triplets, "dgTMatrix")
    expect_s4_class(dense_triplets, "dgTMatrix")
    qs2::qs_save(triplets, fs::path(cache_dir, "ds", "inner_fold_mat.qs2"))
    qs2::qs_save(rows, fs::path(cache_dir, "ds", "inner_stdev_mat.qs2"))
    qs2::qs_save(dense_triplets, fs::path(cache_dir, "ds", "projected_mat.qs2"))
    qs2::qs_save(triplets, fs::path(cache_dir, "ds", "projected_fold.qs2"))

    load_all_data(cache_dir, datasets = "ds")
    inner_fold <- get_mc_data("ds", "inner_fold_mat")
    expect_s4_class(inner_fold, "dgCMatrix")
    expect_identical(inner_fold, sparse)
    # what the app reads from it
    expect_identical(inner_fold["g3", ], triplets["g3", ])
    expect_identical(inner_fold[c("g9", "g1"), , drop = FALSE], methods::as(triplets[c("g9", "g1"), , drop = FALSE], "CsparseMatrix"))
    expect_identical(as.matrix(inner_fold[, 4:9]), as.matrix(triplets[, 4:9]))
    expect_identical(Matrix::rowSums(inner_fold), Matrix::rowSums(triplets))
    # other sparse classes are kept, and a mostly non-zero projected matrix is still made dense
    expect_identical(get_mc_data("ds", "inner_stdev_mat"), rows)
    expect_identical(get_mc_data("ds", "projected_mat"), as.matrix(dense_triplets))
    expect_identical(get_mc_data("ds", "projected_fold"), sparse)
})
