# ==============================================================================
# Gene Module Utilities
# ==============================================================================

#' Parse gene modules file
#'
#' @param file Path to gene modules file
#' @return Tibble with gene and module columns
#' @noRd
parse_gene_modules_file <- function(file) {
    gene_modules <- fread(file, colClasses = c("gene" = "character", "module" = "character")) %>% as_tibble()

    for (field in c("gene", "module")) {
        if (!has_name(gene_modules, field)) {
            cli_abort("File should have a column named {.field {field}}")
        }
    }

    verify_gene_modules(gene_modules)

    return(gene_modules)
}

verify_gene_modules <- function(gene_modules) {
    # verify that every gene exists only in a single module

    bad_genes <- gene_modules %>%
        count(gene) %>%
        filter(n > 1) %>%
        pull(gene)

    if (length(bad_genes) > 0) {
        cli_abort("The following genes appear in more than one module: {genes}", genes = paste(bad_genes, collapse = ", "))
    }
}

#' Calculate gene modules by clustering the egc matrix (metacell gene fractions)
#'
#' @param mat a matrix where genes are rows and metacells are columns
#' @param k number of clusters. If NULL - the number of clusters would be determined such that an gene module would contain 16 genes on average.
#' @param n_genes target number of genes to consider.
#' @param minimal_max_log_fraction take only genes with at least one value
#' (in log2 fraction units - normalized egc) above this threshold
#' @param minimal_relative_log_fraction take only genes with relative
#' log2 fraction (mc_fp) above this this value
#' @param n_downsamp number of UMIs to downsample from each gene. Genes with less umis would be
#' inflated to that number.
#' @param egc_epsilon egc regularization
#'
#' @inheritParams tglkmeans::TGL_kmeans_tidy
#'
#' @noRd
calc_gene_modules <- function(mat, k = NULL, n_genes = 1000, minimal_max_log_fraction = -13, minimal_relative_log_fraction = 3, n_downsamp = 5000, egc_epsilon = 1e-5, verbose = FALSE, seed = 60427) {
    cand_genes <- select_top_fold_genes(mat, n_genes, minimal_max_log_fraction, minimal_relative_log_fraction, egc_epsilon)

    k <- k %||% max(1, round(length(cand_genes) / (16 * 2)))

    m <- mat[cand_genes, , drop = FALSE]

    egc <- t(t(m) / colSums(m))
    legc <- dafr::fast_log(egc, eps = egc_epsilon, base = 2)
    legc_norm <- sweep(legc, 1, matrixStats::rowMedians(legc, na.rm = TRUE))
    # legc_norm <- sweep(legc, 1, matrixStats::rowSds(legc, na.rm = TRUE))

    cli_alert_info("Clustering in order to get gene modules. k = {.val {k}}")
    cli_alert_info("Number of genes considered = {.val {nrow(m)}}")

    # m_ds <- downsample_mat_rows(m, n_downsamp = n_downsamp, inflate = TRUE)
    km <- tglkmeans::TGL_kmeans_tidy(legc_norm, k = k, id_column = FALSE, verbose = verbose, seed = seed)

    gene_modules <- km$cluster %>%
        select(gene = id, clust) %>%
        arrange(clust, gene)

    # Name every gene module by its top gene (according to the order given by 'select_top_fold_genes').
    # Uses "M<cluster>.<top_gene>" convention (Metacells.jl style) instead of "mod_<gene>".
    gene_modules <- gene_modules %>%
        left_join(tibble(gene = cand_genes, rank = 1:length(cand_genes)), by = "gene") %>%
        group_by(clust) %>%
        mutate(module = paste0("M", clust, ".", gene[which.min(rank)])) %>%
        ungroup() %>%
        add_count(module) %>%
        arrange(desc(n), rank) %>%
        select(gene, module)

    return(gene_modules)
}

select_top_fold_genes <- function(mat, n_genes = 1000,
                                  minimal_max_log_fraction = -13, minimal_relative_log_fraction = 3,
                                  egc_epsilon = 1e-5) {
    egc <- t(t(mat) / colSums(mat))

    legc <- dafr::fast_log(egc, eps = egc_epsilon, base = 2)
    fold_matrix <- sweep(legc, 1, matrixStats::rowMedians(legc, na.rm = TRUE))
    fold_matrix <- abs(fold_matrix)

    max_log_fractions_of_genes <- matrixStats::rowMaxs(legc, na.rm = TRUE, useNames = TRUE)
    max_relative_fraction_of_genes <- matrixStats::rowMaxs(fold_matrix, na.rm = TRUE, useNames = TRUE)

    interesting_genes_mask <- (max_log_fractions_of_genes >= minimal_max_log_fraction) &
        (max_relative_fraction_of_genes >= minimal_relative_log_fraction)

    fold_matrix <- fold_matrix[interesting_genes_mask, ]

    if (nrow(fold_matrix) > n_genes) {
        gene_ranks <- matrixStats::colRanks(fold_matrix, useNames = TRUE, preserveShape = TRUE, ties.method = "first")
        gene_ranks[fold_matrix < 2] <- NA
        max_gene_ranks <- matrixStats::rowMaxs(gene_ranks, na.rm = TRUE, useNames = TRUE)

        cand_genes <- sort(max_gene_ranks, decreasing = TRUE) %>%
            head(n = n_genes) %>%
            names()
    } else {
        cand_genes <- rownames(fold_matrix)
        # sort cand genes by maximum relative fraction
        cand_genes <- names(sort(max_relative_fraction_of_genes[cand_genes], decreasing = TRUE))
    }

    return(cand_genes)
}

get_module_genes <- function(module, gene_modules) {
    gene_modules %>%
        filter(module == !!module) %>%
        pull(gene)
}
