#' Update gene modules
#'
#' Change the gene modules of an MCView app.
#'
#' This is usually done after a first iteration of annotation using the "Gene modules" tab in the MCView annotation, which can export a valid \code{gene_modules_file}.
#' The file should have a column named "gene" with the gene names and another column named "module" with the id of the gene module.
#' Note that the exported file from the __MCView__ app might contain additional fields which will be ignored in this function.
#'
#' Under the hood - MCView updates a file named "gene_modules.tsv" under \code{project/cache/dataset}, which can also be edited manually.
#'
#'
#' @inheritParams import_dataset
#'
#' @export
#'
#' @examples
#' \dontrun{
#' update_gene_modules("PBMC", "PBMC163k", "raw/gene-modules.csv")
#' }
#'
#' @export
update_gene_modules <- function(project, dataset, gene_modules_file) {
    verify_project_dir(project)
    verify_app_cache(project)

    gene_modules <- parse_geneparse_gene_modules_file(gene_modules_file)

    serialize_shiny_data(gene_modules, "gene_modules", dataset = dataset, cache_dir = project_cache_dir(project), flat = TRUE)

    cli_alert_success("Succesfully changed gene modules")
}


parse_gene_modules_file <- function(file) {
    gene_modules <- fread(file) %>% as_tibble()

    for (field in c("gene", "module")) {
        if (!has_name(gene_modules, field)) {
            cli_abort("{.field {file}} should have a column named {.field {field}}")
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

    cli_alert_info("Clustering in order to get gene modules. k = {.val {k}}")
    cli_alert_info("Number of genes considered = {.val {nrow(m)}}")

    m_ds <- downsample_mat_rows(m, n_downsamp = n_downsamp, inflate = TRUE)
    km <- tglkmeans::TGL_kmeans_tidy(m_ds, k = k, id_column = FALSE, verbose = verbose, seed = seed)

    gene_modules <- km$cluster %>%
        select(gene = id, clust) %>%
        arrange(clust, gene)

    # name every gene module by it's top gene (according to the order given by 'select_top_fold_genes')
    gene_modules <- gene_modules %>%
        left_join(tibble(gene = cand_genes, rank = 1:length(cand_genes)), by = "gene") %>%
        group_by(clust) %>%
        mutate(module = paste0("mod_", gene[which.min(rank)])) %>%
        ungroup() %>%
        arrange(module, rank) %>%
        select(gene, module)

    return(gene_modules)
}

select_top_fold_genes <- function(mat, n_genes = 1000,
                                  minimal_max_log_fraction = -13, minimal_relative_log_fraction = 3,
                                  egc_epsilon = 1e-5) {
    egc <- t(t(mat) / colSums(mat))

    legc <- log2(egc + egc_epsilon)
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

downsample_mat_rows <- function(m, n_downsamp = NULL, inflate = FALSE) {
    row_sums <- rowSums(m, na.rm = TRUE)
    n_downsamp <- n_downsamp %||% min(row_sums)

    f <- row_sums > n_downsamp
    n_obs <- ncol(m)

    m_ds <- apply(m[f, , drop = FALSE], 1, function(x) {
        tabulate(sample(rep(1:length(x), times = round(x)), replace = FALSE, size = n_downsamp), n_obs)
    }) %>% t()

    if (inflate) {
        m_ds <- rbind(
            m_ds,
            m[!f, , drop = FALSE] * (n_downsamp / row_sums[!f])
        )
    }

    return(m_ds)
}

##' @param n_boot number of bootstrap iterations
# bootstrap_clusters <- function(m, n_downsamp, k, n_boot){
#     N_genes <- nrow(m)
#     tot_coclust <- matrix(0, nrow = N_genes, ncol = N_genes, dimnames = list(cand_genes, cand_genes))
#     doMC::registerDoMC(parallel::detectCores() / 2)

#     bootstrap_res <- plyr::llply(1:n_boot, function(i) {
#         cli_alert_info("Bootstrap: {.val {i}}")
#         m_ds <- downsample_mat_rows(m, n_downsamp = n_downsamp, inflate = TRUE)
#         km <- tglkmeans::TGL_kmeans(m_ds, k = k, id_column = FALSE, verbose = verbose, seed = seed)
#         isclust_ci <- diag(max(km$cluster))[, km$cluster]
#         coclust_ij <- t(isclust_ci) %*% isclust_ci
#         tot_coclust <- tot_coclust + coclust_ij
#     }, .parallel = TRUE)


# }
