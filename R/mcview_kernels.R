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
