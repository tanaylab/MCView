#' Top-K per column of a fold-point matrix with marker thresholds
#'
#' Replaces the R `select_top_fold_genes_per_metacell` loop with a
#' single-pass cpp11 + OpenMP kernel.
#'
#' @param fp Numeric matrix (genes x metacells) with rownames and colnames.
#' @param k Positive integer -- rows returned per column.
#' @param min_val Mask values `< min_val` after adding `fold_reg`.
#' @param fold_reg Constant added to every non-NA value before masking.
#' @param use_abs Rank by `abs(x)` when TRUE, else by raw `x`.
#'
#' @return tibble with columns `(metacell, gene, rank, fp)` of length
#'   `k * ncol(fp)`. Rows where fewer than `k` values qualified receive
#'   `NA` for `gene` and `fp`.
#' @export
top_fold_per_col <- function(fp, k, min_val, fold_reg, use_abs) {
    if (!is.matrix(fp)) {
        fp <- as.matrix(fp)
    }
    if (!is.numeric(fp)) stop("fp must be numeric")
    if (is.null(rownames(fp)) || is.null(colnames(fp))) {
        stop("fp must have rownames and colnames")
    }
    k <- as.integer(k)
    if (is.na(k) || k <= 0L) stop("k must be a positive integer")
    stopifnot(length(min_val) == 1L, is.numeric(min_val))
    stopifnot(length(fold_reg) == 1L, is.numeric(fold_reg))
    stopifnot(length(use_abs) == 1L, is.logical(use_abs))

    out <- top_fold_per_col_cpp(
        fp = fp,
        k = k,
        min_val = as.numeric(min_val),
        fold_reg = as.numeric(fold_reg),
        use_abs = use_abs,
        row_names = rownames(fp),
        col_names = colnames(fp)
    )
    tibble::as_tibble(out)
}

#' Streaming top-K Pearson correlation (block gemm + per-row heap)
#'
#' Replaces the full G×G correlation materialization in
#' `calc_gg_mc_top_cor`. Tiles row-blocks of size `block_rows` and
#' maintains per-row bounded heaps for positive and negative
#' correlations. Never holds the full cor matrix in memory.
#'
#' @param X Numeric matrix (genes × metacells).
#' @param k Positive integer — top-K retained per sign per gene1.
#' @param diag Include self-correlations? Default FALSE.
#' @param min_cor Drop any `|cor| < min_cor`. Default 0.
#' @param block_rows Tile size. Default `getOption("mcview.streaming_block_rows", 512L)`.
#'
#' @return tibble with columns `(gene1, gene2, cor, type)` of length up
#'   to `2 * k * nrow(X)`.
#' @export
streaming_top_k_cor <- function(X, k, diag = FALSE, min_cor = 0,
                                block_rows = getOption("mcview.streaming_block_rows", 512L)) {
    if (!is.matrix(X)) X <- as.matrix(X)
    if (!is.numeric(X)) stop("X must be numeric")
    if (is.null(rownames(X))) stop("X must have rownames")
    k <- as.integer(k)
    if (is.na(k) || k <= 0L) stop("k must be a positive integer")
    block_rows <- as.integer(block_rows)
    if (is.na(block_rows) || block_rows <= 0L) block_rows <- 512L

    # dgemm runs between OpenMP parallel regions (not inside them), so BLAS
    # threading and outer OpenMP do not overlap. Leave BLAS at its configured
    # thread count so the gemm — the dominant cost at G=O(10k), M=O(100) — can
    # use all cores.

    out <- streaming_top_k_cor_cpp(
        X = X, k = k, diag = isTRUE(diag),
        min_cor = as.numeric(min_cor),
        row_names = rownames(X),
        block_rows = block_rows
    )
    out <- tibble::as_tibble(out)
    out <- out[!is.na(out$cor), , drop = FALSE]
    out
}
