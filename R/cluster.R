
cluster_egc <- function(egc, k = NULL, filter_genes = FALSE, min_max_thresh = 3, min_gene_expr = -12, verbose = FALSE) {
    legc <- log2(egc + 1e-5)
    if (filter_genes) {
        gmax <- matrixStats::rowMaxs(legc, na.rm = TRUE)
        gmin <- matrixStats::rowMins(legc, na.rm = TRUE)
        f <- (gmax - gmin) >= min_max_thresh & gmax >= min_gene_expr
        m <- legc[f, ]
    } else {
        m <- legc
    }

    cli_alert_info("using {nrow(m)} genes")

    if (is.null(k)) {
        k <- round(max(5, min(64, ncol(legc) / 30)))
    }
    
    cli_alert_info("clustering k = {.val {k}}")
    cli_alert_info("number of features = {.val {nrow(egc)}}")
    km <- tglkmeans::TGL_kmeans(t(m), k = k, id_column = FALSE, verbose = verbose)

    res <- list(
        clusters = tibble(metacell = colnames(legc), cluster = factor(km$cluster)),
        centers = km$centers
    )

    return(res)
}
